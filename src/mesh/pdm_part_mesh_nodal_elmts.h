#ifndef __PDM_PART_MESH_NODAL_ELMTS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"

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

typedef struct _pdm_part_mesh_nodal_elmts_t PDM_part_mesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);


void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t*
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_H__ */
