#ifndef __PDM_DMESH_NODAL_TO_DMESH_H__
#define __PDM_DMESH_NODAL_TO_DMESH_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */



/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_dmesh_nodal_to_dmesh_t PDM_dmesh_nodal_to_dmesh_t;


typedef enum {
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE = 0,
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE = 1,
} PDM_dmesh_nodal_to_dmesh_transform_t;


typedef enum {
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_NONE    = 0,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE = 1,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE = 2,
  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_VTX  = 3,
} PDM_dmesh_nodal_to_dmesh_translate_group_t;


/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
);


void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         i_mesh,
        PDM_dmesh_nodal_t          *dmn
);

void
PDM_dmesh_nodal_to_dmesh_compute
(
        PDM_dmesh_nodal_to_dmesh_t                 *dmn_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t        transform_kind,
  const PDM_dmesh_nodal_to_dmesh_translate_group_t  transform_group_kind
);

void
PDM_dmesh_nodal_to_dmesh_get_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t  *dmn_to_dm,
  const int                          i_mesh,
        PDM_dmesh_t                **dm
);

void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
);

void
PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         extract_dim
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
