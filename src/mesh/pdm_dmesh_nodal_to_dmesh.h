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

// typedef struct _PDM_DMesh_nodal_t PDM_DMesh_nodal_t;


typedef enum {
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE = 1,
  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE = 2,
} PDM_dmesh_nodal_to_dmesh_transform_t;


typedef enum {
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
int
PDM_dmesh_nodal_to_dmesh_create
(
      int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
);


void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
  const int hdl,
  const int i_mesh,
  const int dmesh_nodal_id
);

void
PDM_dmesh_nodal_to_dmesh_compute
(
  const int                                  hdl,
  const PDM_dmesh_nodal_to_dmesh_transform_t transform_kind
);

void
PDM_dmesh_nodal_to_dmesh_free
(
const int hdl
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_H__ */
