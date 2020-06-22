
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_cuda_error.cuh"
#include "pdm_cuda.cuh"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_handles.h"
#include "pdm_handles.cuh"
#include "pdm_para_octree.h"
#include "pdm_para_octree.cuh"
#include "pdm_timer.h"
#include "pdm_timer.cuh"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define NTIMER 11

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \struct _l_octant_t
 * \brief  Define a list of octants
 *
 */

typedef struct  {

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */

  PDM_morton_code_t *codes;        /*!< Morton codes */

  int  *n_points;          /*!< Number of points in octant*/
  int  *range;             /*!< Start index of point list for each octant */

  int   *neighbour_idx;
  int   *neighbours;               /*!< rank + id_node size = 2 * n_nodes */
  int   dim;

} _l_octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct  {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /*!< Translation for the normalization */
  double      d[3];           /*!< Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t        t_n_points;     /*!< total number of points */
  int                n_points;       /*!< Number of points in each cloud */
  double            *points;         /*!< Point coordinates */
  int               *points_icloud;  /*!< Point cloud */
  PDM_g_num_t       *points_gnum;    /*!< Point global number */
  PDM_morton_code_t *points_code;    /*!< Morton codes */

  PDM_morton_code_t *rank_octants_index;
  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

  int  n_part_boundary_elt;    /*!< Number of partitioning boundary element */
  int *part_boundary_elt_idx; /*!< Index for part_boundary_elt (size=\ref n_part_boundary_elt + 1 */
  int *part_boundary_elt;     /*!< Partitioning boundary elements description (proc number + element number) */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

  int neighboursToBuild;

  int  n_connected;
  int *connected_idx;

} _octree_t;



/*============================================================================
 * Global variable
 *============================================================================*/

__managed__ static PDM_Handles_t *_octrees    = NULL;

//static const double _eps_default  = 1.e-12;

/*============================================================================
 * Private function definitions
 *============================================================================*/

 /**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

__device__
static _octree_t *
_get_from_id
(
 int  id
 )
{
  _octree_t *octree = (_octree_t *) PDM_Handles_get_GPU (_octrees, id);
  if (octree == NULL) {
    PDM_error_GPU(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}

//remove duplicate if cc >= 3.5, and use __host__ __device__ keyword before function
__host__
static _octree_t *
_get_from_id2
(
 int  id
 )
{
  _octree_t *octree = (_octree_t *) PDM_Handles_get (_octrees, id);

  if (octree == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create an octree structure on a GPU
 *
 * \param [in]   n_point_cloud          Number of point cloud
 * \param [in]   depth_max              Maximum depth
 * \param [in]   points_in_leaf_max     Maximum points in a leaf
 * \param [in]   build_leaf_neighbours  Build leaf nieghbours (1 = true)
 * \param [in]   comm                   MPI communicator
 *
 * \return     Identifier
 */


int
PDM_para_octree_create_GPU
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const int build_leaf_neighbours,
 const PDM_MPI_Comm comm
 )
{

  if (_octrees == NULL) {
    _octrees = PDM_Handles_create_GPU (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store_GPU (_octrees, octree);

  octree->dim = 3;

  for (int i = 0; i < octree->dim; i++) {
    octree->global_extents[i]   = -HUGE_VAL;
    octree->global_extents[octree->dim+i] =  HUGE_VAL;
    octree->s[i]         = 0.;
    octree->d[i]         = 0.;
  }

  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;

  octree->n_point_clouds = n_point_cloud;
  octree->t_n_points = 0;
  octree->n_points = 0;
  octree->points = NULL;
  octree->points_icloud = NULL;
  octree->points_gnum = NULL;
  octree->points_code = NULL;

  octree->rank_octants_index = NULL;
  octree->octants = NULL;

  octree->n_part_boundary_elt = 0;
  octree->part_boundary_elt_idx = NULL;
  octree->part_boundary_elt = NULL;

  octree->neighboursToBuild = build_leaf_neighbours;

  octree->comm = comm;

  octree->n_connected = 0;
  octree->connected_idx = NULL;

  octree->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    octree->times_elapsed[i] = 0.;
    octree->times_cpu[i] = 0.;
    octree->times_cpu_u[i] = 0.;
    octree->times_cpu_s[i] = 0.;
  }

  PDM_printf("Depth max from cpu : %d\n", octree->depth_max);

  return id;

}

/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

__host__
void
PDM_para_octree_free_GPU
(
  const int          id
  )
{
  _octree_t *octree = _get_from_id2 (id);

  if (octree->points != NULL) {
    gpuErrchk(cudaFree (octree->points));
  }

  if (octree->points_icloud != NULL) {
    gpuErrchk(cudaFree (octree->points_icloud));
  }

  if (octree->points_gnum != NULL) {
    gpuErrchk(cudaFree (octree->points_gnum));
  }

  if (octree->points_code != NULL) {
    gpuErrchk(cudaFree (octree->points_code));
  }

  if (octree->part_boundary_elt_idx != NULL) {
    gpuErrchk(cudaFree (octree->part_boundary_elt_idx));
  }

  if (octree->part_boundary_elt != NULL) {
    gpuErrchk(cudaFree (octree->part_boundary_elt));
  }

  if (octree->rank_octants_index != NULL) {
    gpuErrchk(cudaFree (octree->rank_octants_index));
  }

  if (octree->octants != NULL) {

    if (octree->octants->codes != NULL) {
      gpuErrchk(cudaFree (octree->octants->codes));
    }

    if (octree->octants->range != NULL) {
      gpuErrchk(cudaFree (octree->octants->range));
    }

    if (octree->octants->n_points != NULL) {
      gpuErrchk(cudaFree (octree->octants->n_points));
    }

    if (octree->octants->neighbour_idx != NULL) {
      gpuErrchk(cudaFree (octree->octants->neighbour_idx));
    }

    if (octree->octants->neighbours != NULL) {
      gpuErrchk(cudaFree (octree->octants->neighbours));
    }

    gpuErrchk(cudaFree (octree->octants));
  }

  if (octree->connected_idx != NULL) {
    gpuErrchk(cudaFree (octree->connected_idx));
  }

  PDM_timer_free_GPU (octree->timer);

  gpuErrchk(cudaFree (octree));

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }
}

//test
__global__
void
print_from_gpu
(
 int id
 )
{
  printf("Hello World! from thread [%d,%d] From device\n", threadIdx.x,blockIdx.x);
  _octree_t *octree = _get_from_id (id);
  printf("Depth max from GPU : %d\n", octree->depth_max);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
