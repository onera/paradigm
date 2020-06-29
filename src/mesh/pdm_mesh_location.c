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
#include "pdm_point_location.h"

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
_location_points_in_boxes_octree
(
 const PDM_MPI_Comm   comm,
 const int            n_points,
 const PDM_g_num_t   *pts_g_num,
 const double        *pts_coord,
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
  //PDM_para_octree_dump_times (octree_id);


  /***************************************
   * Locate points inside boxes
   ***************************************/
  PDM_para_octree_points_inside_boxes (octree_id,
                                       n_boxes,
                                       box_extents,
                                       box_g_num,
                                       pts_in_box_idx,
                                       pts_in_box_g_num,
                                       pts_in_box_coord);


  /***************************************
   * Free octree
   ***************************************/
  PDM_para_octree_free (octree_id);
}




#if 0
static void
_get_candidate_elements_from_dbbtree
(
 const PDM_MPI_Comm   comm,
 _point_cloud_t      *pcloud,
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
#endif

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


void
PDM_mesh_location_get2
(
 const int     id,
 const int     i_point_cloud,
 const int     i_part,
 int          *n_points,
 double      **coords,
 PDM_g_num_t **gnum,
 PDM_g_num_t **location,
 int         **weights_idx,
 double      **weights
 )
{
  _PDM_location_t *mesh_location = _get_from_id (id);

  assert (mesh_location->point_clouds != NULL);
  assert (i_point_cloud < mesh_location->n_point_cloud);

  _point_cloud_t *pcloud = mesh_location->point_clouds + i_point_cloud;

  assert (i_part < pcloud->n_part);

  *n_points    = pcloud->n_points[i_part];
  *coords      = pcloud->coords[i_part];
  *gnum        = pcloud->gnum[i_part];
  *location    = pcloud->location[i_part];
  *weights_idx = pcloud->weights_idx[i_part];
  *weights     = pcloud->weights[i_part];
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
PDM_mesh_location_compute
(
 const int id
 )
{
  const float eps_dist = 1.e-6;

  const int DEBUG = 1;
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
                                           box_extents + 6*ibox);

      /* get elements gnum */
      PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (location->mesh_nodal_id,
                                                     id_block,
                                                     ipart);

      int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                  id_block,
                                                  ipart);

      for (int ielt = 0; ielt < n_elt; ielt++) {
        box_g_num[ibox++] = _gnum[ielt];
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




  PDM_dbbtree_t *dbbt = NULL;
  if (location->method == PDM_MESH_LOCATION_DBBTREE) {

    dbbt = PDM_dbbtree_create (location->comm, dim);

    PDM_dbbtree_boxes_set (dbbt,
                           1,//const int n_part,
                           &n_boxes,
                           &box_extents,
                           &box_g_num);

    PDM_timer_hang_on(location->timer);
    e_t_elapsed = PDM_timer_elapsed(location->timer);
    e_t_cpu     = PDM_timer_cpu(location->timer);
    e_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    e_t_cpu_s   = PDM_timer_cpu_sys(location->timer);

    location->times_elapsed[SEARCH_CANDIDATES] += e_t_elapsed - b_t_elapsed;
    location->times_cpu[SEARCH_CANDIDATES]     += e_t_cpu - b_t_cpu;
    location->times_cpu_u[SEARCH_CANDIDATES]   += e_t_cpu_u - b_t_cpu_u;
    location->times_cpu_s[SEARCH_CANDIDATES]   += e_t_cpu_s - b_t_cpu_s;

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);
  }


  /*
   * Locate points
   */
  int         *pts_idx   = NULL;
  PDM_g_num_t *pts_g_num = NULL;
  double      *pts_coord = NULL;

  for (int icloud = 0; icloud < location->n_point_cloud; icloud++) {

    PDM_timer_hang_on(location->timer);
    b_t_elapsed = PDM_timer_elapsed(location->timer);
    b_t_cpu     = PDM_timer_cpu(location->timer);
    b_t_cpu_u   = PDM_timer_cpu_user(location->timer);
    b_t_cpu_s   = PDM_timer_cpu_sys(location->timer);
    PDM_timer_resume(location->timer);


    _point_cloud_t *pcloud = location->point_clouds + icloud;


    /*
     * Concatenate point cloud partitions
     */
    int n_pts_pcloud = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      n_pts_pcloud += pcloud->n_points[ipart];
    }
    PDM_g_num_t *pcloud_g_num = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
    double      *pcloud_coord = malloc (sizeof(double)      * n_pts_pcloud * dim);
    int idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {

        pcloud_g_num[idx] = pcloud->gnum[ipart][ipt];

        for (int idim = 0; idim < dim; idim++) {
          pcloud_coord[dim*idx + idim] = pcloud->coords[ipart][dim*ipt + idim];
        }

        idx++;
      }
    }


    /*
     * Get points inside bounding boxes of elements
     */

    switch (location->method) {
    case PDM_MESH_LOCATION_OCTREE:
      _location_points_in_boxes_octree (location->comm,
                                        n_pts_pcloud,
                                        pcloud_g_num,
                                        pcloud_coord,
                                        n_boxes,
                                        box_extents,
                                        box_g_num,
                                        &pts_idx,
                                        &pts_g_num,
                                        &pts_coord);
      break;

    case PDM_MESH_LOCATION_DBBTREE:
      /*_location_points_in_boxes_dbbtree (location->comm,
                                         n_pts_pcloud,
                                         pcloud_g_num,
                                         pcloud_coord,
                                         n_boxes,
                                         box_g_num,
                                         dbbt,
                                         &pts_idx,
                                         &pts_g_num,
                                         &pts_coord);*/
      PDM_dbbtree_points_inside_boxes (dbbt,
                                       n_pts_pcloud,
                                       pcloud_g_num,
                                       pcloud_coord,
                                       n_boxes,
                                       box_g_num,
                                       &pts_idx,
                                       &pts_g_num,
                                       &pts_coord);
      break;

    default:
      printf("Error: unknown location method %d\n", location->method);
      assert (1 == 0);
    }
    free (pcloud_coord);


    if (DEBUG) {
      printf("\n[%d] --- Pts in box ---\n", my_rank);
      for (ibox = 0; ibox < n_boxes; ibox++) {

        if (pts_idx[ibox+1] <= pts_idx[ibox]) {
          continue;
        }

        printf("[%d] %d (%ld): ", my_rank, ibox, box_g_num[ibox]);
        for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
          /*printf("((%ld); %f %f %f) ",
            pts_g_num[i], pts_coord[dim*i], pts_coord[dim*i+1], pts_coord[dim*i+2]);*/
          printf("(%ld) ", pts_g_num[i]);
        }
        printf("\n");
      }
      printf("[%d] ------------------\n\n\n", my_rank);
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

    b_t_elapsed = e_t_elapsed;
    b_t_cpu     = e_t_cpu;
    b_t_cpu_u   = e_t_cpu_u;
    b_t_cpu_s   = e_t_cpu_s;
    PDM_timer_resume(location->timer);



    /*
     * Location in elements
     */
    int n_pts = pts_idx[n_boxes];
    PDM_g_num_t *pts_location = malloc (sizeof(PDM_g_num_t) * n_pts);
    for (ibox = 0; ibox < n_boxes; ibox++) {
      for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
        pts_location[i] = box_g_num[ibox];
      }
    }

    //-->>
    PDM_Mesh_nodal_elt_t *pts_elt_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_pts);
    ibox = 0;
    for (int iblock = 0; iblock < n_blocks; iblock++) {
      int id_block = blocks_id[iblock];

      PDM_Mesh_nodal_elt_t elt_type = PDM_Mesh_nodal_block_type_get (location->mesh_nodal_id,
                                                                     id_block);

      for (int ipart = 0; ipart < n_parts; ipart++) {
        int n_elt = PDM_Mesh_nodal_block_n_elt_get (location->mesh_nodal_id,
                                                    id_block,
                                                    ipart);
        for (int ielt = 0; ielt < n_elt; ielt++) {
          for (int i = pts_idx[ibox]; i < pts_idx[ibox+1]; i++) {
            pts_elt_type[i] = elt_type;
          }
          ibox++;
        }
      }
    }
    //<<--

    float  *distance         = NULL;
    double *projected_coords = NULL;
    int    *bar_coords_idx   = NULL;
    double *bar_coords       = NULL;

    const double tolerance = 1e-6;
    int base_element_num = 0;//???

    PDM_point_location_nodal (location->mesh_nodal_id,
                              n_pts,
                              pts_idx,
                              pts_coord,
                              pts_g_num,//debug only
                              tolerance,
                              base_element_num,
                              &distance,
                              &projected_coords,
                              &bar_coords_idx,
                              &bar_coords);

    if (DEBUG) {
      for (int i = 0; i < n_pts; i++) {
        printf("Point gnum = (%ld)\n", pts_g_num[i]);
        printf("\t  coords = (%f, %f, %f)\n", pts_coord[dim*i], pts_coord[dim*i+1], pts_coord[dim*i+2]);
        printf("\tlocation = (%ld)\n", pts_location[i]);
        printf("\tdistance = %f\n", distance[i]);
        printf("\t  bar_co =");
        double sum = 0;
        for (int j = bar_coords_idx[i]; j < bar_coords_idx[i+1]; j++) {
          printf(" %f", (float) bar_coords[j]);
          sum += bar_coords[j];
        }
        printf("  (sum = %f)\n\n", sum);
      }
    }
    free (pts_coord);
    free (pts_idx);



    PDM_g_num_t *pcloud_location       = NULL;
    int         *pcloud_weights_stride = NULL;
    int         *pcloud_weights_idx    = NULL;
    double      *pcloud_weights        = NULL;

    /* Multiple ranks */
    if (n_procs > 1) {

      /*
       * Merge location data of each point
       *   part_to_block_exchange : pts_location, pts_distance, bar_coords...
       *   (for each point, keep data linked to element with minimal distance,
       *    if two elements are at the same distance, keep the one with smallest type)
       *
       */

      PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           &pts_g_num,
                                                           NULL,
                                                           &n_pts,
                                                           1,
                                                           location->comm);

      const int n_pts_block = PDM_part_to_block_n_elt_block_get (ptb);

      int *block_stride = NULL;
      PDM_g_num_t *block_location1 = NULL;
      PDM_Mesh_nodal_elt_t *block_elt_type1 = NULL;//
      float       *block_distance1 = NULL;
      double      *block_bar_coords1 = NULL;

      int *stride = malloc (sizeof(int) * n_pts);
      for (int i = 0; i < n_pts; i++) {
        stride[i] = 1;
      }

      PDM_part_to_block_exch (ptb,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_VAR,
                              1,
                              &stride,
                              (void **) &pts_location,
                              &block_stride,
                              (void **) &block_location1);
      free (block_stride);
      free (pts_location);

      PDM_part_to_block_exch (ptb,
                              sizeof(PDM_Mesh_nodal_elt_t),
                              PDM_STRIDE_VAR,
                              1,
                              &stride,
                              (void **) &pts_elt_type,
                              &block_stride,
                              (void **) &block_elt_type1);
      free (block_stride);
      free (pts_elt_type);

      PDM_part_to_block_exch (ptb,
                              sizeof(float),
                              PDM_STRIDE_VAR,
                              1,
                              &stride,
                              (void **) &distance,
                              &block_stride,
                              (void **) &block_distance1);
      free (block_stride);
      free (distance);


      int *bar_coords_stride = malloc (sizeof(int) * n_pts);
      for (int i = 0; i < n_pts; i++) {
        bar_coords_stride[i] = bar_coords_idx[i+1] - bar_coords_idx[i];
      }
      free (bar_coords_idx);

      PDM_part_to_block_exch (ptb,
                              sizeof(double),
                              PDM_STRIDE_VAR,
                              1,
                              &bar_coords_stride,
                              (void **) &bar_coords,
                              &block_stride,
                              (void **) &block_bar_coords1);
      free (block_stride);
      free (bar_coords);

      int *block_n_vtx_elt1 = NULL;
      PDM_part_to_block_exch (ptb,
                              sizeof(int),
                              PDM_STRIDE_VAR,
                              1,
                              &stride,
                              (void **) &bar_coords_stride,
                              &block_stride,
                              (void **) &block_n_vtx_elt1);
      free (bar_coords_stride);
      free (stride);


      /*
       * Keep closest elements
       */
      PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get (ptb);
      PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);
      PDM_g_num_t n_pts_block2 = block_distrib_idx[my_rank+1] - block_distrib_idx[my_rank];
      PDM_g_num_t *block_location2 = malloc (sizeof(PDM_g_num_t) * n_pts_block2);
      int *block_bar_coords_stride2 = malloc (sizeof(int) * n_pts_block2);
      for (int i = 0; i < n_pts_block2; i++) {
        block_location2[i] = -1;
      }

      int *block_idx1 = malloc (sizeof(int) * (n_pts_block + 1));
      block_idx1[0] = 0;
      for (int i = 0; i < n_pts_block; i++) {
        block_idx1[i+1] = block_idx1[i] + block_stride[i];
      }

      int *idx_min = malloc (sizeof(int) * n_pts_block2);
      idx = 0;
      for (int i = 0; i < n_pts_block; i++) {
        int n_elt = block_stride[i];

        idx_min[i] = idx;

        if (n_elt > 1) {
          float min_dist = HUGE_VAL;
          PDM_Mesh_nodal_elt_t type_min = -1;
          idx_min[i] = -1;

          for (int j = idx; j < idx + n_elt; j++) {
            if (min_dist > block_distance1[j] ||
                (min_dist < block_distance1[j] + eps_dist &&
                 type_min > block_elt_type1[j])) {
              min_dist = block_distance1[j];
              type_min = block_elt_type1[j];
              idx_min[i] = j;
            }
          }
        }

        int _ipt = block_g_num[i] - 1 - block_distrib_idx[my_rank];
        block_location2[_ipt] = block_location1[idx_min[i]];
        block_bar_coords_stride2[_ipt] = block_n_vtx_elt1[idx_min[i]];

        idx += n_elt;
      }



      int *block_bar_coords_idx1 = malloc (sizeof(int) * (block_idx1[n_pts_block] + 1));
      block_bar_coords_idx1[0] = 0;
      for (int i = 0; i < block_idx1[n_pts_block]; i++) {
        block_bar_coords_idx1[i+1] = block_bar_coords_idx1[i] + block_n_vtx_elt1[i];
      }

      int *block_bar_coords_idx2 = malloc (sizeof(int) * (n_pts_block2 + 1));
      block_bar_coords_idx2[0] = 0;
      for (int i = 0; i < n_pts_block2; i++) {
        if (block_location2[i] < 0) {
          block_bar_coords_stride2[i] = 0;
        }

        block_bar_coords_idx2[i+1] = block_bar_coords_idx2[i] + block_bar_coords_stride2[i];
      }

      double *block_bar_coords2 = malloc (sizeof(double) * block_bar_coords_idx2[n_pts_block2]);
      for (int i = 0; i < n_pts_block; i++) {
        int _ipt = block_g_num[i] - 1 - block_distrib_idx[my_rank];

        for (int j = 0; j < block_bar_coords_stride2[_ipt]; j++) {
          block_bar_coords2[block_bar_coords_idx2[_ipt] + j] =
            block_bar_coords1[block_bar_coords_idx1[idx_min[i]] + j];
        }
      }
      free (idx_min);
      free (block_stride);
      free (block_n_vtx_elt1);
      free (block_bar_coords_idx1);
      free (block_distance1);
      free (block_elt_type1);
      free (block_location1);
      free (block_bar_coords1);
      free (block_bar_coords_idx2);

      // block_to_part_exch...
      PDM_part_to_block_t *ptb2 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_MERGE,
                                                            1.,
                                                            &pcloud_g_num,
                                                            NULL,
                                                            &n_pts_pcloud,
                                                            1,
                                                            location->comm);

      PDM_g_num_t *block_distrib_idx2 = PDM_part_to_block_distrib_index_get (ptb2);

      PDM_block_to_part_t *btp = PDM_block_to_part_create (block_distrib_idx2,
                                                           (const PDM_g_num_t **) &pcloud_g_num,
                                                           &n_pts_pcloud,
                                                           1,
                                                           location->comm);
      free (pcloud_g_num);

      int one = 1;

      pcloud_location = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
      if (1) {
        for (int i = 0; i < n_pts_pcloud; i++) {
          pcloud_location[i] = -1;
        }
      }
      PDM_block_to_part_exch (btp,
                              sizeof(PDM_g_num_t),
                              PDM_STRIDE_CST,
                              &one,
                              block_location2,
                              NULL,
                              (void **) &pcloud_location);
      free (block_location2);

      /*int *block_bar_coords_stride2 = malloc (sizeof(int) * n_pts_block);
      for (int i = 0; i < n_pts_block; i++) {
        block_bar_coords_stride2[i] = block_bar_coords_idx2[i+1] - block_bar_coords_idx2[i];
      }
      free (block_bar_coords_idx2);*/

      pcloud_weights_stride = malloc (sizeof(int) * n_pts_pcloud);
      PDM_block_to_part_exch (btp,
                              sizeof(int),
                              PDM_STRIDE_CST,
                              &one,
                              (void *) block_bar_coords_stride2,
                              NULL,
                              (void **) &pcloud_weights_stride);

      pcloud_weights_idx = malloc (sizeof(int) * (n_pts_pcloud + 1));
      pcloud_weights_idx[0] = 0;
      for (int i = 0; i < n_pts_pcloud; i++) {
        pcloud_weights_idx[i+1] = pcloud_weights_idx[i] + pcloud_weights_stride[i];
      }

      pcloud_weights = malloc (sizeof(double) * pcloud_weights_idx[n_pts_pcloud]);
      PDM_block_to_part_exch (btp,
                              sizeof(double),
                              PDM_STRIDE_VAR,
                              block_bar_coords_stride2,
                              (void *) block_bar_coords2,
                              &pcloud_weights_stride,
                              (void **) &pcloud_weights);
      free (block_bar_coords2);
      free (block_bar_coords_stride2);

      PDM_part_to_block_free (ptb);
      PDM_part_to_block_free (ptb2);
      PDM_block_to_part_free (btp);
    }

    /* Single rank */
    else {
      /*
       * Keep closest elements
       * (if two elements are at the same distance, keep the one with smallest type)
       * !!! PB si plusieurs nuages de points (max(gnum) > n_pts_pcloud...)
       */
      int   *idx_min  = malloc (sizeof(int)   * n_pts_pcloud);
      float *dist_min = malloc (sizeof(float) * n_pts_pcloud);
      PDM_Mesh_nodal_elt_t *type_min = malloc (sizeof(PDM_Mesh_nodal_elt_t) * n_pts_pcloud);
      for (int ipt = 0; ipt < n_pts_pcloud; ipt++) {
        dist_min[ipt] = HUGE_VAL;
        type_min[ipt] = -1;
        idx_min[ipt] = 0;
      }

      for (int i = 0; i < n_pts; i++) {
        PDM_g_num_t id_pt = pts_g_num[i] - 1;

        if (dist_min[id_pt] > distance[i] ||
            (dist_min[id_pt] < distance[i] + eps_dist &&
             type_min[id_pt] > pts_elt_type[i])) {
          idx_min[id_pt] = i;
          dist_min[id_pt] = distance[i];
          type_min[id_pt] = pts_elt_type[i];
        }
      }
      free (dist_min);
      free (distance);
      free (type_min);
      free (pts_elt_type);

      pcloud_location       = malloc (sizeof(PDM_g_num_t) * n_pts_pcloud);
      pcloud_weights_stride = malloc (sizeof(int)         * n_pts_pcloud);
      pcloud_weights_idx    = malloc (sizeof(int)         * (n_pts_pcloud + 1));
      pcloud_weights_idx[0] = 0;

      for (int ipt = 0; ipt < n_pts_pcloud; ipt++) {
        PDM_g_num_t id_pt = pcloud_g_num[ipt] - 1;
        int i = idx_min[id_pt];

        pcloud_location[ipt] = pts_location[i];
        pcloud_weights_stride[ipt] = bar_coords_idx[i+1] - bar_coords_idx[i];
        pcloud_weights_idx[ipt+1] = pcloud_weights_idx[ipt] + pcloud_weights_stride[ipt];
      }
      free (pts_location);

      pcloud_weights = malloc (sizeof(double) * pcloud_weights_idx[n_pts_pcloud]);
      for (int ipt = 0; ipt < n_pts_pcloud; ipt++) {
        PDM_g_num_t id_pt = pcloud_g_num[ipt] - 1;
        int i = idx_min[id_pt];

        for (int j = 0; j < pcloud_weights_stride[ipt]; j++) {
          pcloud_weights[pcloud_weights_idx[ipt] + j] =
            bar_coords[bar_coords_idx[i] + j];
        }
      }

#if 0
      printf("* * * * * *\n");
      for (int ipt = 0; ipt < n_pts_pcloud; ipt++) {
        printf("Pt (%ld), weights =", pcloud_g_num[ipt]);
        for (int i = pcloud_weights_idx[ipt]; i < pcloud_weights_idx[ipt+1]; i++) {
          printf(" %f", pcloud_weights[i]);
        }
        printf("\n");
      }
      printf("* * * * * *\n");
#endif

      free (idx_min);
      free (bar_coords_idx);
      free (bar_coords);
      free (pcloud_g_num);
    }



    /*
     * Conform to original partitioning of current point cloud
     */
    assert (pcloud->location == NULL);

    pcloud->location    = malloc (sizeof(PDM_g_num_t *) * pcloud->n_part);
    pcloud->weights_idx = malloc (sizeof(int *)         * pcloud->n_part);
    pcloud->weights     = malloc (sizeof(double *)      * pcloud->n_part);

    idx = 0;
    for (int ipart = 0; ipart < pcloud->n_part; ipart++) {
      pcloud->location[ipart]    = malloc (sizeof(PDM_g_num_t) * pcloud->n_points[ipart]);
      pcloud->weights_idx[ipart] = malloc (sizeof(int)         * (pcloud->n_points[ipart] + 1));

      pcloud->weights_idx[ipart][0] = 0;
      int idx_tmp = idx;
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        pcloud->location[ipart][ipt] = pcloud_location[idx];
        pcloud->weights_idx[ipart][ipt+1] = pcloud->weights_idx[ipart][ipt] + pcloud_weights_stride[idx];

        idx++;
      }

      pcloud->weights[ipart] = malloc (sizeof(double) * pcloud->weights_idx[ipart][pcloud->n_points[ipart]]);

      idx = idx_tmp;
      for (int ipt = 0; ipt < pcloud->n_points[ipart]; ipt++) {
        for (int j = 0; j < pcloud_weights_stride[idx]; j++) {
          pcloud->weights[ipart][pcloud->weights_idx[ipart][ipt] + j] =
            pcloud_weights[pcloud_weights_idx[idx] + j];
        }
        idx++;
      }
    }
    free (pcloud_location);
    free (pcloud_weights_stride);
    free (pcloud_weights_idx);
    free (pcloud_weights);

    free (pts_g_num);


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

  } // Loop over point clouds

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


int
PDM_mesh_location_mesh_nodal_id_get
(
 const int id
 )
{
  _PDM_location_t *location = _get_from_id (id);

  return location->mesh_nodal_id;
}

#ifdef	__cplusplus
}
#endif
