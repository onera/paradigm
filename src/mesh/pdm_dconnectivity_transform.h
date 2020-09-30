#ifndef PDM_DCONNECTIVITY_TRANSFORM_H_
#define PDM_DCONNECTIVITY_TRANSFORM_H_

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

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *cell_distrib,
 const PDM_g_num_t     *face_distrib,
 const int             *dcell_face_idx,
 const PDM_g_num_t     *dcell_face,
 const int             *dface_vtx_idx,
 const PDM_g_num_t     *dface_vtx,
       int            **dcell_vtx_idx,
       PDM_g_num_t    **dcell_vtx
);


void
PDM_deduce_dual_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
       int              is_signed,
       int            **dentity2_entity1_idx,
       PDM_g_num_t    **dentity2_entity1
);


#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DCONNECTIVITY_TRANSFORM_H_ */
