/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_surf_part.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_surf_part_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_surf_part_t *
PDM_surf_part_create
(
  const int          n_face,
  const int         *face_vtx_idx,
  const int         *face_vtx,
  const PDM_g_num_t *face_ln_to_gn,
  const int          n_vtx,
  const double      *coords,
  const PDM_g_num_t *vtx_ln_to_gn
)
{
  PDM_surf_part_t *_part;
  PDM_malloc(_part, 1 ,PDM_surf_part_t);

  _part->n_face        = n_face;
  _part->face_vtx_idx  = face_vtx_idx;
  _part->face_vtx      = face_vtx;
  _part->face_ln_to_gn = face_ln_to_gn;
  _part->n_vtx         = n_vtx;
  _part->coords        = coords;
  _part->vtx_ln_to_gn  = vtx_ln_to_gn;
  _part->extents       = NULL;

  return _part;
}


PDM_surf_part_t *
PDM_surf_part_free
(
 PDM_surf_part_t * part
)
{
  assert (part != NULL);

  if (part != NULL) {
    part->face_vtx_idx = NULL;
    part->face_vtx = NULL;

    part->face_ln_to_gn = NULL;
    part->coords = NULL;
    part->vtx_ln_to_gn = NULL;

    if (part->extents != NULL)
      PDM_free(part->extents);

    PDM_free(part);
  }

  return NULL;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
