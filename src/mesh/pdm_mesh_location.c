/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
/* #include "pdm_mesh_dist.h" */
#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
/* #include "pdm_octree.h" */
#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_timer.h"
#include "pdm_hash_tab.h"
#include "pdm_mesh_location.h"

#include "pdm_binary_search.h"
#include "pdm_para_octree.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 6

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                            = 0,
  BUILD_BOUNDING_BOXES             = 1,
  SEARCH_CANDIDATES                = 2,
  DISTRIBUTE_ELEMENTARY_OPERATIONS = 3,
  COMPUTE_ELEMENTARY_LOCATIONS     = 4,
  END                              = 5,

} _ol_timer_step_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  PDM_g_num_t **location;
  double      **uvw;
  int         **weights_idx;
  double      **weights; /*!< Barycentric coordinates */

} _point_cloud_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  int  n_point_cloud; /*!< Number of point clouds */
  PDM_MPI_Comm comm;  /*!< MPI communicator */

  PDM_mesh_nature_t mesh_nature;  /*!< Nature of the mesh */

  int  shared_nodal;   /*!< 1 if mesh nodal is shared, 0 otherwise */
  int  mesh_nodal_id;  /*!< Mesh identifier */
  int _mesh_nodal_id;

  _point_cloud_t *point_clouds; /*!< Point clouds */

  double tolerance;

  PDM_mesh_location_method_t method;

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


} _PDM_location_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_locations   = NULL;

static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_location_t *
_get_from_id
(
 int  id
 )
{
  _PDM_location_t *location = (_PDM_location_t *) PDM_Handles_get (_locations, id);

  if (location == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_location error : Bad identifier\n");
  }

  return location;
}




static void
_get_candidate_elements_from_octree
(
 const PDM_MPI_Comm   comm,
 _point_cloud_t      *pcloud,
 const int            n_boxes,
 const double        *box_extents,
 const PDM_g_num_t   *box_g_num,
 PDM_g_num_t        **candidates_g_num,
 int                **candidates_idx
 )
{
  int my_rank;
  PDM_MPI_Comm_rank (comm, &my_rank);

  int n_ranks;
  PDM_MPI_Comm_rank (comm, &n_ranks);

  /**************************************
   * Build octree from point cloud
   ***************************************/
  const int depth_max = 15;//?
  const int points_in_leaf_max = 1;
  const int build_leaf_neighbours = 0;


  /* Concatenate partitions */
  int n_points = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    n_points += pcloud->n_points[ipart];
  }

  double      *pts_coord = malloc (sizeof(double) *      n_points*3);
  PDM_g_num_t *pts_g_num = malloc (sizeof(PDM_g_num_t) * n_points);
  int idx = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
      for (int idim = 0; idim < 3; idim++) {
        pts_coord[3*idx + idim] = pcloud->coords[ipart][3*ipt + idim];
      }
      pts_g_num[idx] = pcloud->gnum[ipart][ipt];
      idx++;
    }
  }

  /* Create empty parallel octree structure */
  int octree_id = PDM_para_octree_create (1,
                                          depth_max,
                                          points_in_leaf_max,
                                          build_leaf_neighbours,
                                          comm);

  /* Set octree point cloud */
  PDM_para_octree_point_cloud_set (octree_id,
                                   0,
                                   n_points,
                                   pts_coord,
                                   pts_g_num);

  /* Build parallel octree */
  PDM_para_octree_build (octree_id);
  //PDM_para_octree_dump (octree_id);
  PDM_para_octree_dump_times (octree_id);


  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_MERGE,
                                                       1.,
                                                       &pts_g_num,
                                                       NULL,
                                                       &n_points,
                                                       1,
                                                       comm);

  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  int *block_n_candidates = NULL;
  PDM_g_num_t *block_candidates_g_num = NULL;

  PDM_para_octree_location_boxes_get (octree_id,
                                      n_boxes,
                                      box_extents,
                                      box_g_num,
                                      &block_candidates_g_num,
                                      &block_n_candidates);

  PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                       (const PDM_g_num_t **) &pts_g_num,
                                                       &n_points,
                                                       1,
                                                       comm);
  int *part_stride = malloc (sizeof(int) * n_points);
  int one = 1;
  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_CST,
                          &one,
                          (void *) block_n_candidates,
                          NULL,
                          (void **) &part_stride);

  *candidates_idx = malloc (sizeof(int) * (n_points+1));
  int *_candidates_idx = *candidates_idx;
  _candidates_idx[0] = 0;
  for (int i = 0; i < n_points; i++) {
    _candidates_idx[i+1] = _candidates_idx[i] + part_stride[i];
  }

  *candidates_g_num = malloc (sizeof(PDM_g_num_t) * _candidates_idx[n_points]);

  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          block_n_candidates,
                          (void *) block_candidates_g_num,
                          &part_stride,
                          (void **) candidates_g_num);

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);

  free (part_stride);
  free (block_n_candidates);
  free (block_candidates_g_num);
  free (pts_coord);
  free (pts_g_num);

  /***************************************
   * Free octree
   ***************************************/
  PDM_para_octree_free (octree_id);
}





static void
_location_points_in_boxes_octree
(
 const PDM_MPI_Comm   comm,
 _point_cloud_t      *pcloud,
 const int            n_boxes,
 const double        *box_extents,
 const PDM_g_num_t   *box_g_num,
 int                **pts_in_box_idx,
 PDM_g_num_t        **pts_in_box_g_num,
 double             **pts_in_box_coord
 )
{
  int my_rank;
  PDM_MPI_Comm_rank (comm, &my_rank);

  int n_ranks;
  PDM_MPI_Comm_rank (comm, &n_ranks);

  /**************************************
   * Build octree from point cloud
   ***************************************/
  const int depth_max = 15;//?
  const int points_in_leaf_max = 1;
  const int build_leaf_neighbours = 0;


  /* Concatenate partitions */
  int n_points = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    n_points += pcloud->n_points[ipart];
  }

  double      *pts_coord = malloc (sizeof(double) *      n_points*3);
  PDM_g_num_t *pts_g_num = malloc (sizeof(PDM_g_num_t) * n_points);
  int idx = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
      for (int idim = 0; idim < 3; idim++) {
        pts_coord[3*idx + idim] = pcloud->coords[ipart][3*ipt + idim];
      }
      pts_g_num[idx] = pcloud->gnum[ipart][ipt];
      idx++;
    }
  }

  /* Create empty parallel octree structure */
  int octree_id = PDM_para_octree_create (1,
                                          depth_max,
                                          points_in_leaf_max,
                                          build_leaf_neighbours,
                                          comm);

  /* Set octree point cloud */
  PDM_para_octree_point_cloud_set (octree_id,
                                   0,
                                   n_points,
                                   pts_coord,
                                   pts_g_num);

  /* Build parallel octree */
  PDM_para_octree_build (octree_id);
  //PDM_para_octree_dump (octree_id);
  PDM_para_octree_dump_times (octree_id);

  /***************************************
   * Locate points inside boxes
   ***************************************/
  PDM_para_octree_points_inside_boxes (octree_id,
                                       n_boxes,
                                       box_extents,
                                       box_g_num,
                                       &pts_in_box_idx,
                                       &pts_in_box_g_num,
                                       &pts_in_box_coord);

  /***************************************
   * Free octree
   ***************************************/
  PDM_para_octree_free (octree_id);
}











static void
_get_candidate_elements_from_dbbtree
(
 const PDM_MPI_Comm   comm,
 _point_cloud_t      *pcloud,
 const int            n_boxes,
 PDM_dbbtree_t       *dbbt,
 PDM_g_num_t        **candidates_g_num,
 int                **candidates_idx
 )
{
  int my_rank;
  PDM_MPI_Comm_rank (comm, &my_rank);

  int n_ranks;
  PDM_MPI_Comm_rank (comm, &n_ranks);


  /* Concatenate partitions */
  int n_points = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    n_points += pcloud->n_points[ipart];
  }

  double      *pts_coord = malloc (sizeof(double)      * n_points*3);
  PDM_g_num_t *pts_g_num = malloc (sizeof(PDM_g_num_t) * n_points);
  int idx = 0;
  for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
    for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
      for (int idim = 0; idim < 3; idim++) {
        pts_coord[3*idx + idim] = pcloud->coords[ipart][3*ipt + idim];
      }
      pts_g_num[idx] = pcloud->gnum[ipart][ipt];
      idx++;
    }
  }

  /* Find candidates in dbbtree */
  PDM_dbbtree_location_boxes_get (dbbt,
                                  n_points,
                                  pts_coord,
                                  pts_g_num,
                                  candidates_idx,
                                  candidates_g_num);

  free (pts_coord);
  free (pts_g_num);
}


/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute the location of point clouds inta a mesh
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 *
 */

int
PDM_mesh_location_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
 )
{
  if (_locations == NULL) {
    _locations = PDM_Handles_create (4);
  }

  _PDM_location_t *location = (_PDM_location_t *) malloc(sizeof(_PDM_location_t));

  int id = PDM_Handles_store (_locations, location);

  location->n_point_cloud = n_point_cloud;
  location->comm = comm;
  location->mesh_nature = mesh_nature;

  location->shared_nodal   = 0;
  location->mesh_nodal_id  = -1;
  location->_mesh_nodal_id = -1;

  location->point_clouds =
    (_point_cloud_t*) malloc (sizeof(_point_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    location->point_clouds[i].n_part = -1;
    location->point_clouds[i].n_points = NULL;
    location->point_clouds[i].coords = NULL;
    location->point_clouds[i].gnum = NULL;
    location->point_clouds[i].location = NULL;
    location->point_clouds[i].uvw = NULL;
    location->point_clouds[i].weights = NULL;
    location->point_clouds[i].weights_idx = NULL;
  }

  location->tolerance = 0.;

  location->method = PDM_MESH_LOCATION_OCTREE;

  location->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    location->times_elapsed[i] = 0.;
    location->times_cpu[i]     = 0.;
    location->times_cpu_u[i]   = 0.;
    location->times_cpu_s[i]   = 0.;
  }

  return id;

}

void
PDM_mesh_location_create_cf
(
 const PDM_mesh_nature_t mesh_nature,
 const int n_point_cloud,
 const PDM_MPI_Fint comm,
 int *id
 )
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_mesh_location_create (mesh_nature, n_point_cloud, _comm);

}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_location_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_part = n_part;
  location->point_clouds[i_point_cloud].n_points =
    realloc(location->point_clouds[i_point_cloud].n_points, n_part * sizeof(int));
  location->point_clouds[i_point_cloud].coords =
    realloc(location->point_clouds[i_point_cloud].coords,
            n_part * sizeof(double *));
  location->point_clouds[i_point_cloud].gnum =
    realloc(location->point_clouds[i_point_cloud].gnum,
            n_part * sizeof(PDM_g_num_t *));

  for (int i = 0; i < n_part; i++) {
    location->point_clouds[i_point_cloud].n_points[i] = -1;
    location->point_clouds[i_point_cloud].coords[i] = NULL;
    location->point_clouds[i_point_cloud].gnum[i] = NULL;
  }

}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_mesh_location_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
 double      *coords,
 PDM_g_num_t *gnum
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->point_clouds[i_point_cloud].n_points[i_part] = n_points;
  location->point_clouds[i_point_cloud].coords[i_part] = coords;
  location->point_clouds[i_point_cloud].gnum[i_part] = gnum;

}


/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   id             Identifier
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 *
 */

void
PDM_mesh_location_shared_nodal_mesh_set
(
 const int  id,
 const int  mesh_nodal_id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->mesh_nodal_id = mesh_nodal_id;
  location->shared_nodal = 1;
}


/**
 *
 * \brief Set global data of a mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_mesh_location_mesh_global_data_set
(
 const int  id,
 const int  n_part
 )
{
  _PDM_location_t *location = _get_from_id (id);

  if ((location->shared_nodal == 0) && (location->mesh_nodal_id != -1)) {
    PDM_Mesh_nodal_free (location->mesh_nodal_id);
  }

  location->mesh_nodal_id = PDM_Mesh_nodal_create (n_part, location->comm);
}


/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_mesh_location_part_set
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
 )
{
  _PDM_location_t *location = _get_from_id (id);

  /*
   * Creation de mesh nodal
   */

  PDM_Mesh_nodal_coord_set (location->mesh_nodal_id,
                            i_part,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);



  PDM_l_num_t *face_vtx_nb  = malloc (sizeof(PDM_l_num_t) * n_face);
  PDM_l_num_t *cell_face_nb = malloc (sizeof(PDM_l_num_t) * n_cell);

  for (int i = 0; i < n_face; i++) {
    face_vtx_nb[i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  }

  for (int i = 0; i < n_cell; i++) {
    cell_face_nb[i] = cell_face_idx[i+1] - cell_face_idx[i];
  }

  PDM_Mesh_nodal_cell3d_cellface_add (location->mesh_nodal_id,
                                      i_part,
                                      n_cell,
                                      n_face,
                                      face_vtx_idx,
                                      face_vtx_nb,
                                      face_vtx,
                                      cell_face_idx,
                                      cell_face_nb,
                                      cell_face,
                                      cell_ln_to_gn);
  free (face_vtx_nb);//?
  free (cell_face_nb);//?
  /*PDM_Mesh_nodal_cell2d_celledge_add
    (
    const int          idx,
    const int          id_part,
    const int          n_elt,
    const int          n_edge,
    const PDM_l_num_t *edge_vtx_idx,
    const PDM_l_num_t *edge_vtx_nb,
    const PDM_l_num_t *edge_vtx,
    const PDM_l_num_t *cell_edge_idx,
    const PDM_l_num_t *cell_edge_nb,
    const PDM_l_num_t *cell_edge,
    const PDM_g_num_t *numabs
    );*/
}


/**
 *
 * \brief Set the tolerance for bounding boxes
 *
 * \param [in]   id              Identifier
 * \param [in]   tol             Tolerance
 *
 */

void
PDM_mesh_location_tolerance_set
(
 const int    id,
 const double tol
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->tolerance = tol;
}


/**
 *
 * \brief Set the method for computing location
 *
 * \param [in]   id              Identifier
 * \param [in]   method          Method
 *
 */

void
PDM_mesh_location_method_set
(
 const int                        id,
 const PDM_mesh_location_method_t method
 )
{
  _PDM_location_t *location = _get_from_id (id);

  location->method = method;
}




/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_location_compute
(
 const int id
 )
{
  const int dim = 3;

  _PDM_location_t *location = _get_from_id (id);

  int my_rank;
  PDM_MPI_Comm_rank (location->comm, &my_rank);

  int n_procs;
  PDM_MPI_Comm_size (location->comm, &n_procs);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  location->times_elapsed[BEGIN] = PDM_timer_elapsed(location->timer);
  location->times_cpu[BEGIN]     = PDM_timer_cpu(location->timer);
  location->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(location->timer);
  location->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(location->timer);

  b_t_elapsed = location->times_elapsed[BEGIN];
  b_t_cpu     = location->times_cpu[BEGIN];
  b_t_cpu_u   = location->times_cpu_u[BEGIN];
  b_t_cpu_s   = location->times_cpu_s[BEGIN];
  PDM_timer_resume(location->timer);

  /*
   * Build the bounding boxes of mesh elements
   */
  int n_blocks = PDM_Mesh_nodal_n_blocks_get (location->mesh_nodal_id);
  int n_parts  = PDM_Mesh_nodal_n_part_get (location->mesh_nodal_id);
  int *blocks_id = PDM_Mesh_nodal_blocks_id_get (location->mesh_nodal_id);

  int n_boxes = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    n_boxes += PDM_Mesh_nodal_n_cell_get (location->mesh_nodal_id,
                                          ipart);
  }

  const int n_origin = 4; // rank, part, block, lnum
  int *box_origin = malloc (sizeof(int) * n_origin * n_boxes);

  PDM_g_num_t *box_g_num   = malloc (sizeof(PDM_g_num_t) * n_boxes);
  double      *box_extents = malloc (sizeof(double)      * n_boxes * 6);
  int ibox = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {

    int id_block = blocks_id[iblock];

    for (int ipart = 0; ipart < n_parts; ipart++) {
      /* get element extents */
      PDM_Mesh_nodal_compute_cell_extents (location->mesh_nodal_id,
                                           id_block,
                                           ipart,
                                           location->tolerance,
                                           (&box_extents + ibox));

      /* get elements gnum */
      PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (location->mesh_nodal_id,
                                                     id_block,
                                                     ipart);

      int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                  id_block,
                                                  ipart);

      for (int ielt = 0; ielt < n_elt; ielt++) {
        box_g_num[ibox] = _gnum[ielt];

        box_origin[n_origin*ibox]   = my_rank;
        box_origin[n_origin*ibox+1] = ipart;
        box_origin[n_origin*ibox+2] = iblock;
        box_origin[n_origin*ibox+3] = ielt;

        ibox++;
      }
    }
  }

  /* Store boxes origin (i.e. rank, ipart, iblock, ielt) in block data */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                       &box_g_num,
                                                       NULL,
                                                       &n_boxes,
                                                       1,
                                                       location->comm);

  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  int *block_stride = NULL;
  int *block_box_origin = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(int),
                          PDM_STRIDE_CST,
                          n_origin,
                          NULL,
                          (void **) &box_origin,
                          &block_stride,
                          (void **) &block_box_origin);
  free (block_stride);
  free (box_origin);


  PDM_timer_hang_on(location->timer);
  e_t_elapsed = PDM_timer_elapsed(location->timer);
  e_t_cpu     = PDM_timer_cpu(location->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

  location->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  location->times_cpu[BUILD_BOUNDING_BOXES]     += e_t_cpu - b_t_cpu;
  location->times_cpu_u[BUILD_BOUNDING_BOXES]   += e_t_cpu_u - b_t_cpu_u;
  location->times_cpu_s[BUILD_BOUNDING_BOXES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(location->timer);


  PDM_dbbtree_t *dbbt = NULL;
  if (location->method == PDM_MESH_LOCATION_DBBTREE) {
    dbbt = PDM_dbbtree_create (location->comm, dim);

    PDM_dbbtree_boxes_set (dbbt,
                           1,//const int n_part,
                           &n_boxes,
                           &box_extents,
                           &box_g_num);
  }


  /*
   * Search candidates
   *   - coder localisation dans ddbtree
   */
  PDM_g_num_t *candidates_g_num = NULL;
  int         *candidates_idx   = NULL;

  for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {
    _point_cloud_t *pcloud = location->point_clouds + icloud;

    switch (location->method) {
    case PDM_MESH_LOCATION_OCTREE:
      _get_candidate_elements_from_octree (location->comm,
                                           pcloud,
                                           n_boxes,
                                           box_extents,
                                           box_g_num,
                                           &candidates_g_num,
                                           &candidates_idx);
      break;

    case PDM_MESH_LOCATION_DBBTREE:
      _get_candidate_elements_from_dbbtree (location->comm,
                                            pcloud,
                                            n_boxes,
                                            dbbt,
                                            &candidates_g_num,
                                            &candidates_idx);
      break;

    default:
      printf("Error: unknown location method %d\n", location->method);
      assert (1 == 0);
    }

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

    //----------->>>> Remove that as soon as location in ELEMENTS is implemented
#if 1
    assert (pcloud->location == NULL);

    int idx = 0;
    pcloud->location = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * pcloud->n_points[ipart]);
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {

        assert (candidates_idx[idx+1] > candidates_idx[idx]);

        if (candidates_idx[idx+1] != candidates_idx[idx] + 1) {
          printf("[%d] cloud %d part %d: pt %d (%ld) has %d candidate boxes\n",
                 my_rank,
                 icloud,
                 ipart,
                 ipt,
                 pcloud->gnum[ipart][ipt],
                 candidates_idx[idx+1] - candidates_idx[idx]);
        }
        pcloud->location[ipart][ipt] = candidates_g_num[candidates_idx[idx]];
        idx++;
      }
    }
#endif
    //<<<<-----------

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);

    int n_points = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      n_points += pcloud->n_points[ipart];
    }


    /* Get origin (i.e. rank, ipart, iblock, ielt) of candidate boxes from their gnum */
    PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx,
                                                         (const PDM_g_num_t **) &candidates_g_num,
                                                         &(candidates_idx[n_points]),
                                                         1,
                                                         location->comm);

    int *candidates_origin = malloc (sizeof(int) * n_origin * candidates_idx[n_points]);
    PDM_block_to_part_exch (btp,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            &n_origin,
                            (void *) block_box_origin,
                            NULL,
                            (void **) &candidates_origin);
    free (block_box_origin);

    ptb = PDM_part_to_block_free (ptb);
    btp = PDM_block_to_part_free (btp);

    /*
      idx = 0;
      for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      printf("[%d] part %d\n", my_rank, ipart);
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
      printf("[%d]  point %d\n", my_rank, ipt);
      for (int i = candidates_idx[idx]; i < candidates_idx[idx+1]; i++) {
      printf("[%d]    (%ld)\t{%d, %d, %d, %d}\n",
      my_rank,
      candidates_g_num[i],
      candidates_origin[n_origin*i],
      candidates_origin[n_origin*i+1],
      candidates_origin[n_origin*i+2],
      candidates_origin[n_origin*i+3]);
      }
      idx++;
      }
      }
    */
    free (candidates_g_num);


    /*
     * Distribution of elementary operations
     */
    int *send_count = malloc (sizeof(int) * n_procs);
    for (int i = 0; i < n_procs; i++) {
      send_count[i] = 0;
    }

    idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        for (int i = candidates_idx[idx]; i < candidates_idx[idx+1]; i++) {
          send_count[candidates_origin[n_origin*i]] += 3;
        }
        idx++;
      }
    }

    int *recv_count = malloc (sizeof(int) * n_procs);
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      location->comm);

    int *send_shift = malloc (sizeof(int) * (n_procs + 1));
    int *recv_shift = malloc (sizeof(int) * (n_procs + 1));
    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (int i = 0; i < n_procs; i++) {
      send_shift[i+1] = send_shift[i] + send_count[i];
      recv_shift[i+1] = recv_shift[i] + recv_count[i];
      send_count[i] = 0;
    }

    //-->>
    int *send_origin = malloc (sizeof(int) * send_shift[n_procs]);
    int *recv_origin = malloc (sizeof(int) * recv_shift[n_procs]);

    double *send_coord = malloc (sizeof(double) * send_shift[n_procs]);
    double *recv_coord = malloc (sizeof(double) * recv_shift[n_procs]);

    idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        for (int i = candidates_idx[idx]; i < candidates_idx[idx+1]; i++) {
          int rank = candidates_origin[n_origin*i];

          for (int j = 0; j < dim; j++) {
            send_origin[send_shift[rank] + send_count[rank]] = candidates_origin[n_origin*i+1+j];
            send_coord[send_shift[rank] + send_count[rank]] = pcloud->coords[ipart][dim*ipt+j];
            send_count[rank]++;
          }
        }
        idx++;
      }
    }

    PDM_MPI_Alltoallv (send_origin, send_count, send_shift, PDM_MPI_INT,
                       recv_origin, recv_count, recv_shift, PDM_MPI_INT,
                       location->comm);

    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       location->comm);

    // CHECK
    int n_recv_pts = recv_shift[n_procs] / 3;
    for (int ipt = 0; ipt < n_recv_pts; ipt++) {

      int ipart  = recv_origin[3*ipt];
      int iblock = recv_origin[3*ipt+1];
      int ielt   = recv_origin[3*ipt+2];

      int _ibox = 0;
      for (int jblock = 0; jblock < iblock; jblock++) {
        int id_block = blocks_id[jblock];
        for (int jpart = 0; jpart < ipart; jpart++) {
          int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                      id_block,
                                                      jpart);
          _ibox += n_elt;
        }
      }
      _ibox += ielt;

      int inside = 1;
      double *_pt = recv_coord + dim*ipt;
      double *box_min = box_extents + 2*dim*_ibox;
      double *box_max = box_min + dim;
      for (int i = 0; i < dim; i++) {
        if (_pt[i] < box_min[i] || _pt[i] > box_max[i]) {
          inside = 0;
          break;
        }
      }

      assert (inside);
    }
    //<<--

    //...

    free (send_count);
    free (recv_count);
    free (send_shift);
    free (recv_shift);

    free (send_coord);
    free (send_origin);



    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[DISTRIBUTE_ELEMENTARY_OPERATIONS] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[DISTRIBUTE_ELEMENTARY_OPERATIONS]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[DISTRIBUTE_ELEMENTARY_OPERATIONS]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[DISTRIBUTE_ELEMENTARY_OPERATIONS]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(location->timer);

    /*
     * Elementary location computations
     */
    //..

    free (candidates_idx);
    free (candidates_origin);
    free (recv_coord);
    free (recv_origin);

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[COMPUTE_ELEMENTARY_LOCATIONS] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[COMPUTE_ELEMENTARY_LOCATIONS]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[COMPUTE_ELEMENTARY_LOCATIONS]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;

    PDM_timer_resume(location->timer);

  }
  free (box_g_num);
  free (box_extents);

  if (dbbt != NULL) {
    PDM_dbbtree_free (dbbt);
  }


  PDM_timer_hang_on(location->timer);

  location->times_elapsed[END] = PDM_timer_elapsed(location->timer);
  location->times_cpu[END]     = PDM_timer_cpu(location->timer);
  location->times_cpu_u[END]   = PDM_timer_cpu_user(location->timer);
  location->times_cpu_s[END]   = PDM_timer_cpu_sys(location->timer);

  b_t_elapsed = location->times_elapsed[END];
  b_t_cpu     = location->times_cpu[END];
  b_t_cpu_u   = location->times_cpu_u[END];
  b_t_cpu_s   = location->times_cpu_s[END];
  PDM_timer_resume(location->timer);
}


/**
 *
 * \brief Get mesh location
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  location_elt_g_num    The global number of the closest element if the point is located,
 *                                     -1 otherwise
 *
 */

void
PDM_mesh_location_get
(
 const int           id,
 const int           i_point_cloud,
 const int           i_part,
 PDM_g_num_t **location_elt_gnum
 )
{
  _PDM_location_t *location = _get_from_id (id);

  assert (location->point_clouds != NULL);
  assert (i_point_cloud < location->n_point_cloud);

  _point_cloud_t *pcloud = location->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *location_elt_gnum = pcloud->location[i_part];
}


/**
 *
 * \brief Free a locationd mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_mesh_location_free
(
 const int id,
 const int partial
 )
{
  _PDM_location_t *location = _get_from_id (id);

  /* Free point clouds */
  if (location->point_clouds != NULL) {
    for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {
      _point_cloud_t *pcloud = location->point_clouds + icloud;

      if (pcloud->n_points != NULL) {
        free (pcloud->n_points);
      }

      if (pcloud->coords != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->coords[ipart] != NULL) {
            free (pcloud->coords[ipart]);
          }
        }
        free (pcloud->coords);
      }

      if (pcloud->gnum != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->gnum[ipart] != NULL) {
            free (pcloud->gnum[ipart]);
          }
        }
        free (pcloud->gnum);
      }

      if (pcloud->location != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->location[ipart] != NULL) {
            free (pcloud->location[ipart]);
          }
        }
        free (pcloud->location);
      }

      if (pcloud->uvw != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->uvw[ipart] != NULL) {
            free (pcloud->uvw[ipart]);
          }
        }
        free (pcloud->uvw);
      }

      if (pcloud->weights_idx != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights_idx[ipart] != NULL) {
            free (pcloud->weights_idx[ipart]);
          }
        }
        free (pcloud->weights_idx);
      }

      if (pcloud->weights != NULL) {
        for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
          if (pcloud->weights[ipart] != NULL) {
            free (pcloud->weights[ipart]);
          }
        }
        free (pcloud->weights);
      }

    }
    free (location->point_clouds);
  }

  /* Free mesh nodal */
  // ! if shared (partial ?)
  //PDM_Mesh_nodal_partial_free (location->mesh_nodal_id);
  PDM_Mesh_nodal_free (location->mesh_nodal_id);
}

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_mesh_location_dump_times
(
 const int id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  double t1 = location->times_elapsed[END] - location->times_elapsed[BEGIN];
  double t2 = location->times_cpu[END] - location->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, location->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, location->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (location->times_elapsed,
                     t_elaps_max,
                     NTIMER,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     location->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (location->times_cpu,
                     t_cpu_max, NTIMER,
                     PDM_MPI_DOUBLE,
                     PDM_MPI_MAX,
                     location->comm);

  int rank;
  PDM_MPI_Comm_rank (location->comm, &rank);

  if (rank == 0) {

    PDM_printf( "mesh_location timer : all (elapsed and cpu) :                                      "
                " %12.5es %12.5es\n",
                t1max, t2max);

    PDM_printf( "mesh_location timer : build bounding boxes (elapsed and cpu) :                     "
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BOUNDING_BOXES],
                t_cpu_max[BUILD_BOUNDING_BOXES]);

    PDM_printf( "mesh_location timer : build aux. struct + search candidates (elapsed and cpu) :    "
                " %12.5es %12.5es\n",
                t_elaps_max[SEARCH_CANDIDATES],
                t_cpu_max[SEARCH_CANDIDATES]);

    PDM_printf( "mesh_location timer : distribute elementary operations (+CHECK) (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[DISTRIBUTE_ELEMENTARY_OPERATIONS],
                t_cpu_max[DISTRIBUTE_ELEMENTARY_OPERATIONS]);

    PDM_printf( "mesh_location timer : compute elementary locations (elapsed and cpu) :             "
                " %12.5es %12.5es\n",
                t_elaps_max[COMPUTE_ELEMENTARY_LOCATIONS],
                t_cpu_max[COMPUTE_ELEMENTARY_LOCATIONS]);
  }
}







/**
 *
 * \brief Compute point location
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_location_compute2
(
 const int id
 )
{
  const int dim = 3;

  _PDM_location_t *location = _get_from_id (id);

  int my_rank;
  PDM_MPI_Comm_rank (location->comm, &my_rank);

  int n_procs;
  PDM_MPI_Comm_size (location->comm, &n_procs);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  location->times_elapsed[BEGIN] = PDM_timer_elapsed(location->timer);
  location->times_cpu[BEGIN]     = PDM_timer_cpu(location->timer);
  location->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(location->timer);
  location->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(location->timer);

  b_t_elapsed = location->times_elapsed[BEGIN];
  b_t_cpu     = location->times_cpu[BEGIN];
  b_t_cpu_u   = location->times_cpu_u[BEGIN];
  b_t_cpu_s   = location->times_cpu_s[BEGIN];
  PDM_timer_resume(location->timer);

  /*
   * Build the bounding boxes of mesh elements
   */
  int n_blocks = PDM_Mesh_nodal_n_blocks_get (location->mesh_nodal_id);
  int n_parts  = PDM_Mesh_nodal_n_part_get (location->mesh_nodal_id);
  int *blocks_id = PDM_Mesh_nodal_blocks_id_get (location->mesh_nodal_id);

  int n_boxes = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    n_boxes += PDM_Mesh_nodal_n_cell_get (location->mesh_nodal_id,
                                          ipart);
  }

  /* const int n_origin = 4; // rank, part, block, lnum */
  /* int *box_origin = malloc (sizeof(int) * n_origin * n_boxes); */

  int *nodal_block_idx = malloc (sizeof(int) * (n_blocks+1));
  nodal_block_idx[0] = 0;

  PDM_g_num_t *box_g_num   = malloc (sizeof(PDM_g_num_t) * n_boxes);
  double      *box_extents = malloc (sizeof(double)      * n_boxes * 6);
  int ibox = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {

    nodal_block_idx[iblock+1] = nodal_block_idx[iblock];
    int id_block = blocks_id[iblock];

    for (int ipart = 0; ipart < n_parts; ipart++) {
      /* get element extents */
      PDM_Mesh_nodal_compute_cell_extents (location->mesh_nodal_id,
                                           id_block,
                                           ipart,
                                           location->tolerance,
                                           (&box_extents + ibox));

      /* get elements gnum */
      PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (location->mesh_nodal_id,
                                                     id_block,
                                                     ipart);

      int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                  id_block,
                                                  ipart);
      nodal_block_idx[iblock+1] += n_elt;

      for (int ielt = 0; ielt < n_elt; ielt++) {
        box_g_num[ibox] = _gnum[ielt];

        /* box_origin[n_origin*ibox]   = my_rank; */
        /* box_origin[n_origin*ibox+1] = ipart; */
        /* box_origin[n_origin*ibox+2] = iblock; */
        /* box_origin[n_origin*ibox+3] = ielt; */

        ibox++;
      }
    }
  }

  PDM_timer_hang_on(location->timer);
  e_t_elapsed = PDM_timer_elapsed(location->timer);
  e_t_cpu     = PDM_timer_cpu(location->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

  location->times_elapsed[BUILD_BOUNDING_BOXES] += e_t_elapsed - b_t_elapsed;
  location->times_cpu[BUILD_BOUNDING_BOXES]     += e_t_cpu - b_t_cpu;
  location->times_cpu_u[BUILD_BOUNDING_BOXES]   += e_t_cpu_u - b_t_cpu_u;
  location->times_cpu_s[BUILD_BOUNDING_BOXES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;
  PDM_timer_resume(location->timer);


#if 1
  assert (location->method == PDM_MESH_LOCATION_OCTREE);
#else
  PDM_dbbtree_t *dbbt = NULL;
  if (location->method == PDM_MESH_LOCATION_DBBTREE) {
    dbbt = PDM_dbbtree_create (location->comm, dim);

    PDM_dbbtree_boxes_set (dbbt,
                           1,//const int n_part,
                           &n_boxes,
                           &box_extents,
                           &box_g_num);
  }
#endif

  /*
   * Locate points
   */
  int         *pts_in_box_idx   = NULL;
  PDM_g_num_t *pts_in_box_g_num = NULL;
  double      *pts_in_box_coord = NULL;

  for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {
    _point_cloud_t *pcloud = location->point_clouds + icloud;

    /*
     * Get points inside bounding boxes of elements
     */

    switch (location->method) {
    case PDM_MESH_LOCATION_OCTREE:
      _location_points_in_boxes_octree (location->comm,
                                        pcloud,
                                        n_boxes,
                                        box_extents,
                                        box_g_num,
                                        &pts_in_box_idx,
                                        &pts_in_box_g_num,
                                        &pts_in_box_coord);

      break;

    case PDM_MESH_LOCATION_DBBTREE:
      /* ... */
      /* break; */

    default:
      printf("Error: unknown location method %d\n", location->method);
      assert (1 == 0);
    }



    #if 1
    assert (pcloud->location == NULL);

    int idx = 0;
    pcloud->location = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * pcloud->n_points[ipart]);


    }
    #endif



    free (pts_in_box_idx);
    free (pts_in_box_g_num);
    free (pts_in_box_coord);
  }

}

#ifdef	__cplusplus
}
#endif
