#ifndef __PDM_DCUBE_NODAL_GEN2_H__
#define __PDM_DCUBE_NODAL_GEN2_H__

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dcube_nodal2_t PDM_dcube_nodal2_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a distributed cube
 *
 * \param [out]  id             dcube_nodal identifier
 * \param [in]   comm           Communicator
 * \param [in]   n_vtx_seg        Number of vertices in segments
 * \param [in]   length         Segment length
 * \param [in]   zero_x         Coordinates of the origin
 * \param [in]   zero_y         Coordinates of the origin
 * \param [in]   zero_z         Coordinates of the origin
 *
 */

PDM_dcube_nodal2_t*
PDM_dcube_nodal_gen2_init
(
 PDM_MPI_Comm          comm,
 const PDM_g_num_t     nx,
 const PDM_g_num_t     ny,
 const PDM_g_num_t     nz,
 const double          length,
 const double          zero_x,
 const double          zero_y,
 const double          zero_z,
 PDM_Mesh_nodal_elt_t  t_elt,
 const int             order,
 PDM_ownership_t       owner
 );

PDM_dmesh_nodal_t*
PDM_dcube_nodal_gen2_dmesh_nodal_get
(
 PDM_dcube_nodal2_t  *dcube
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN2_H__ */
