#ifndef __PDM_MESH_CHECK_H__
#define __PDM_MESH_CHECK_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Remove unconnected vertices in a mesh connectivity
 *
 * If a vertex is not cited in the face->vtx connectivity, the function
 * removes it from the mesh to ensure contiguity
 *
 * /TODO : PDM_mesh_check_unconnected_vertex complexity is n^2. This function must be optimized
 *
 * \param [in, out] n_vtx       Number of vertices
 * \param [in, out] l_face_vtx  Size of face->vtx connectivity
 * \param [in, out] face_vtx    Face->vtx connectivity
 * \param [in, out] coords      Vertices coordinates
 * \param [in, out] n_holes     Number of holes
 *
 */

void PDM_mesh_check_unconnected_vertex
(
PDM_g_num_t* nb_vtx,
PDM_g_num_t* l_face_vtx,
PDM_g_num_t* face_vtx,
double*      coords,
int*         nb_holes
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_GNUM_H__ */
