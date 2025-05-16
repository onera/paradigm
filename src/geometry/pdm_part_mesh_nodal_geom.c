/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_geom.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
int
_check_if_all_simplices
(
  PDM_part_mesh_nodal_t  *pmn
)
{
  int all_simplices = 1;

  // Nodal
  PDM_geometry_kind_t  geom_kind    = (pmn->mesh_dimension == 2) ? PDM_GEOMETRY_KIND_SURFACIC : PDM_GEOMETRY_KIND_VOLUMIC;
  PDM_Mesh_nodal_elt_t simplex_type = (pmn->mesh_dimension == 2) ? PDM_MESH_NODAL_TRIA3       : PDM_MESH_NODAL_TETRA4;

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_Mesh_nodal_elt_t elt_type = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                               sections_id[i_section]);

    if (elt_type != simplex_type) {
      all_simplices = 0;
      break;
    }
  }

  return all_simplices;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_part_mesh_nodal_dual_volume_compute
(
  PDM_part_mesh_nodal_t  *pmn,
  double                **dual_vol
)
{

  int all_simplices = _check_if_all_simplices(pmn);

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
