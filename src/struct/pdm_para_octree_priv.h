#ifndef __PDM_PARA_OCTREE_PRIV_H__
#define __PDM_PARA_OCTREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_error.h"
#include "pdm_printf.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

#define PARA_OCTREE_NTIMER 12

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  BUILD_ORDER_POINTS            = 1,
  BUILD_BLOCK_PARTITION         = 2,
  BUILD_LOCAL_NODES             = 3,
  BUILD_LOCAL_NEIGHBOURS_STEP1  = 4,
  BUILD_LOCAL_NEIGHBOURS_STEP2  = 5,
  BUILD_LOCAL_NEIGHBOURS_STEP3  = 6,
  BUILD_LOCAL_NEIGHBOURS        = 7,
  BUILD_DISTANT_NEIGHBOURS      = 8,
  BUILD_EXPLICIT_NODES          = 9,
  BUILD_TOTAL                   = 10,
  END                           = 11,

} _ol_timer_step_t;


/**
 * \struct _heap_t
 * \brief  Heap used to recursively subdivide nodes
 *
 */

typedef struct  {

  int   top;                  /*!< Top of head  */
  int   size;                 /*!< Size of heap */
  PDM_morton_code_t *codes;   /*!< Morton codes */
  int *range;                 /*!< Points range */
  int *n_points;              /*!< Points number */
  int   max_top;

} _heap_t;


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


typedef struct {
  int                n_nodes;
  PDM_morton_code_t *codes;
  int               *n_points;
  int               *range;
  int               *ancestor_id;
  int               *children_id;
  int               *leaf_id; //-1 if internal, >=0 if leaf
  double            *pts_extents;
} _l_explicit_node_t;


typedef struct {

  PDM_mpi_win_shared_t *w_codes;
  PDM_mpi_win_shared_t *w_n_points;
  PDM_mpi_win_shared_t *w_range;
  PDM_mpi_win_shared_t *w_ancestor_id;
  PDM_mpi_win_shared_t *w_children_id;
  PDM_mpi_win_shared_t *w_leaf_id;
  PDM_mpi_win_shared_t *w_pts_extents;

  int                n_nodes;
  PDM_morton_code_t *pt_w_codes;
  int               *pt_w_n_points;
  int               *pt_w_range;
  int               *pt_w_ancestor_id;
  int               *pt_w_children_id;
  int               *pt_w_leaf_id; //-1 if internal, >=0 if leaf
  double            *pt_w_pts_extents;

} _w_l_explicit_node_t;

typedef struct {

  PDM_mpi_win_shared_t *w_codes;
  PDM_mpi_win_shared_t *w_n_points;
  PDM_mpi_win_shared_t *w_range;

  int                n_nodes;
  PDM_morton_code_t *pt_w_codes;
  int               *pt_w_n_points;
  int               *pt_w_range;

  /*PDM_mpi_win_shared_t *w_neighbour_idx;
    PDM_mpi_win_shared_t *w_neighbours;
    int   *pt_w_neighbour_idx;
    int   *pt_w_neighbours;
  */

} _w_l_octant_t;


typedef struct  {

  PDM_mpi_win_shared_t *w_points;
  //PDM_mpi_win_shared_t *w_points_icloud;
  PDM_mpi_win_shared_t *w_points_gnum;
  PDM_mpi_win_shared_t *w_points_code;

  int                n_points;
  double            *pt_w_points;
  int               *pt_w_points_icloud;
  PDM_g_num_t       *pt_w_points_gnum;
  PDM_morton_code_t *pt_w_points_code;

} _w_points_t;


typedef struct {

  PDM_MPI_Request *req_oct;
  PDM_MPI_Request *req_pts;
  PDM_MPI_Request *req_exp;

} _copy_requests_t;


/**
 * \struct _pdm_para_octree_t
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

  double times_elapsed[PARA_OCTREE_NTIMER]; /*!< Elapsed time */

  double times_cpu[PARA_OCTREE_NTIMER];     /*!< CPU time */

  double times_cpu_u[PARA_OCTREE_NTIMER];  /*!< User CPU time */

  double times_cpu_s[PARA_OCTREE_NTIMER];  /*!< System CPU time */

  int neighboursToBuild;

  int  n_connected;
  int *connected_idx;



  PDM_box_set_t  *rank_boxes;            /*!< Rank Boxes */
  int             n_used_rank;           /*!< Number of used ranks */
  int            *used_rank;             /*!< used ranks */
  double         *used_rank_extents;     /*!< Extents of processes */
  PDM_box_tree_t *bt_shared;             /*!< Shared Boundary box tree */
  PDM_MPI_Comm    rank_comm;             /*!< MPI communicator */



  int                 n_copied_ranks;      /*!< Number of copies from other ranks */
  int                *copied_ranks;        /*!< Copied ranks */
  _l_octant_t       **copied_octants;      /*!< Octants from copied ranks */
  int                *n_copied_points;     /*!< Number of points copied from other ranks */
  double            **copied_points;       /*!< Coordinates of copied points */
  PDM_g_num_t       **copied_points_gnum;  /*!< Global numbers of copied points  */
  PDM_morton_code_t **copied_points_code;  /*!< Morton codes of copied points */


  /*
   *  Shared 'coarse' octree
   */
  int               *shared_rank_idx;
  PDM_morton_code_t *shared_codes;
  int               *shared_pts_n;
  double            *shared_pts_extents;


  int explicit_nodes_to_build;
  int use_win_shared;

  _l_explicit_node_t   *explicit_nodes;
  _l_explicit_node_t  **copied_explicit_nodes;
  _l_explicit_node_t  **shm_explicit_nodes;

  _w_l_octant_t        **w_copied_octants;
  _w_points_t          **w_copied_points;
  _w_l_explicit_node_t **w_copied_explicit_nodes;

  _copy_requests_t       copy_requests;



  /* Shared */
  int                 shared_among_nodes;
  PDM_MPI_Comm        comm_shared;

  int                 n_shm_ranks;      /*!< Number of copies from other ranks */
  int                *shm_ranks;        /*!< shm ranks */
  _l_octant_t       **shm_octants;      /*!< Octants from shm ranks */
  int                *n_shm_points;     /*!< Number of points shm from other ranks */
  double            **shm_points;       /*!< Coordinates of shm points */
  PDM_g_num_t       **shm_points_gnum;  /*!< Global numbers of shm points  */
  PDM_morton_code_t **shm_points_code;  /*!< Morton codes of shm points */

  _w_l_octant_t        *w_shm_octants;
  _w_points_t          *w_shm_points;
  _w_l_explicit_node_t *w_shm_explicit_nodes;

  PDM_mpi_win_shared_t* wshared_all_rank_idx;
  PDM_mpi_win_shared_t* wshared_all_node_idx;
  PDM_mpi_win_shared_t* wshared_all_pts_n;
  PDM_mpi_win_shared_t* wshared_all_pts_extents;
  PDM_mpi_win_shared_t* wshared_all_codes;

} _pdm_para_octree_t;



/**
 * \struct _neighbours_tmp_t
 * \brief  Define a temporary neighbour structure
 *
 */


typedef struct  {

  int n_neighbour[6];     /*!< Number of neighbours in the arrays  */
  int s_neighbour[6];     /*!< Size of arrays */
  int *neighbours[6];     /*!< Arrays */

} _neighbours_tmp_t;



/**
 * \struct _min_heap_t
 * \brief  Binary heap used as (min-)priority queue
 *
 */


typedef struct {

  int                size;
  int                count;
  PDM_morton_code_t *code;
  int               *start;
  int               *end;
  double            *dist2;

} _min_heap_t;



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PARA_OCTREE_PRIV_H__ */
