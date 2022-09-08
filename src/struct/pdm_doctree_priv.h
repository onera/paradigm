#ifndef __PDM_DOCTREE_PRIV_H__
#define __PDM_DOCTREE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_octree_seq.h"
#include "pdm_kdtree_seq.h"

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


  int                        global_depth_max;           /*!< global_octree depth_max                   */
  int                        global_points_in_leaf_max;  /*!< global_octree max pts in leaf             */

  int                        local_depth_max;
  int                        local_points_in_leaf_max;
  double                     local_tolerance;

  PDM_octree_seq_t          *global_octree;              /*! Global octree to orient among procs        */
  PDM_octree_seq_t          *local_octree;               /*! Local octree                               */
  PDM_octree_seq_t          *shmem_octree;               /*! Shared octree among cores in current nodes */

  PDM_kdtree_seq_t          *global_kdtree;              /*! Global octree to orient among procs        */
  PDM_kdtree_seq_t          *local_kdtree;               /*! Local octree                               */
  PDM_kdtree_seq_t          *shmem_kdtree;               /*! Shared octree among cores in current nodes */


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
