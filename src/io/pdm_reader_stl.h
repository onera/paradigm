#ifndef __PDM_READER_STL_H__
#define __PDM_READER_STL_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh_nodal.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/**
 *
 * \brief Create a dmesh nodal from a file in ASCII STL mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */
PDM_dmesh_nodal_t *
PDM_reader_stl_dmesh_nodal
(
 PDM_MPI_Comm   comm,
 const char    *filename
 );


#ifdef __cplusplus
}
#endif
#endif  /* __PDM_READER_STL_H__ */
