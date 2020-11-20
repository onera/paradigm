#ifndef __pdm_dmesh_nodal_tO_DMESH_PRIV_H__
#define __pdm_dmesh_nodal_tO_DMESH_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh.h"
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

struct _pdm_dmesh_nodal_to_dmesh_t {

  PDM_MPI_Comm         comm;                    /*!< MPI communicator */
  PDM_ownership_t      owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t           results_is_getted;       /*!< Flags to indicate if result is getted      */

  int                  n_mesh;                  /*!< Number of meshes to manages                */

  PDM_dmesh_nodal_t **dmesh_nodal;
  PDM_dmesh_t       **dmesh;

};

#ifdef  __cplusplus
}
#endif

#endif  /* __pdm_dmesh_nodal_tO_DMESH_PRIV_H__ */
