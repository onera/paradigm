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

  int                        global_depth_max;           /*!< global_octree depth_max                   */
  int                        global_points_in_leaf_max;  /*!< global_octree max pts in leaf             */

  PDM_octree_seq_t          *global_octree;              /*! Global octree to orient among procs        */
  PDM_octree_seq_t          *local_octree;               /*! Local octree                               */
  PDM_octree_seq_t          *shmem_octree;               /*! Shared octree among cores in current nodes */

  PDM_MPI_Comm               comm_shared;


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
