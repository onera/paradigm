#ifndef __PDM_DMESH_NODAL_TO_DMESH_PRIV_H__
#define __PDM_DMESH_NODAL_TO_DMESH_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  PDM_MPI_Comm         comm;                    /*!< MPI communicator */
  PDM_ownership_t      owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t           results_is_getted;       /*!< Flags to indicate if result is getted      */

  int                  n_mesh;                  /*!< Number of meshes to manages                */

  int                 *dmesh_nodal_ids;         /*!< Identifier list of dmesh_nodal struture
                                                 *   (size = n_mesh )                           */

  int                 *dmesh_ids;               /*!< Identifier list of dmesh structure
                                                 *   (size = n_mesh)                            */

  _pdm_dmesh_nodal_t **dmesh_nodal;
  _pdm_dmesh_t       **dmesh;

} _pdm_dmesh_nodal_to_dmesh_t;

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_DMESH_NODAL_TO_DMESH_PRIV_H__ */
