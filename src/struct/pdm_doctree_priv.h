#ifndef __PDM_DOCTREE_PRIV_H__
#define __PDM_DOCTREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_point_tree_seq.h"
#include "pdm_part_to_block.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
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

struct _pdm_doctree_t {

  PDM_MPI_Comm               comm;                       /*!< MPI communicator                          */
  int                        dim;                        /*!< Dimension                                 */
  int                        n_part_cloud;               /*!< Dimension                                 */
  PDM_doctree_local_tree_t   local_tree_kind;


  int                        coarse_depth_max;           /*!< global_octree depth_max                   */
  int                        coarse_points_in_leaf_max;  /*!< global_octree max pts in leaf             */

  int                        local_depth_max;
  int                        local_points_in_leaf_max;
  double                     local_tolerance;

  PDM_point_tree_seq_t      *coarse_tree;              /*! coarse tree to orient among procs        */
  PDM_point_tree_seq_t      *local_tree;               /*! Local tree                               */
  PDM_point_tree_seq_shm_t  *shmem_tree;               /*! Shared tree among cores in current nodes */

  PDM_MPI_Comm               comm_dist_graph;
  int                        n_degree_in;
  int*                       neighbor_in;

  PDM_MPI_Comm               comm_shared;


  /* Cloud - Just reference */
  int          *n_point_cloud;
  PDM_g_num_t **pts_g_num;
  double      **pts_coords;
  int         **pts_init_location;

  /* Solicitation */
  PDM_tree_solicitation_t    solicitation_kind;
  int                        n_part;
  int                       *n_entity;
  int                      **init_location_entity;
  PDM_g_num_t              **entity_gnum;
  double                   **entity_coords;

  /* Equilibrate results */
  PDM_part_to_block_t       *ptb_unit_op_equi;
  int                       *block_pts_in_box_n;
  PDM_g_num_t               *block_pts_in_box_g_num;
  double                    *block_pts_in_box_coord;


};

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*=============================================================================
 * Static function definitions
 *============================================================================*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOCTREE_PRIV_H__ */
