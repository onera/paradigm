/*
 * \file
 */

#ifndef __PDM_PART_MESH_GEOM_H__
#define __PDM_PART_MESH_GEOM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/
#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/**
 *
 * \brief Compute dual volumes
 *
 * \param [in]  pm          Pointer to \ref PDM_part_mesh_t instance
 * \param [in]  dual_vol    For each part, dual volume for each vertices (synchronise at partition interface)
 *
 */
void
PDM_part_mesh_dual_volume_compute
(
  PDM_part_mesh_t   *pm,
  double          ***out_dual_vol
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_GEOM_H__ */
