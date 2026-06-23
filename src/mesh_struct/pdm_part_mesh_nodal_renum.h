#ifndef __PDM_PART_MESH_NODAL_RENUM_H__
#define __PDM_PART_MESH_NODAL_RENUM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

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
PDM_part_mesh_nodal_renum
(
  PDM_part_mesh_nodal_t *pmn,
  const char            *renum_entity_method    [4],
  const int             *renum_entity_properties[4]
);


void
PDM_part_mesh_nodal_reorder_vtx
(
  PDM_part_mesh_nodal_t  *pmn,
  int                   **old_to_new_order_vtx,
  int                   **new_to_old_order_vtx
);


/**
 *
 * \brief Compute new order for a \ref PDM_part_mesh_nodal_elmts_t instance
 *
 * \param [in, out] pmn              Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in,out]  new_to_old_order New order for each part (size=n_part)
 *
 */
void
PDM_part_mesh_nodal_elmts_renum_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old_order
);


/**
 *
 * \brief Perform elmts renumbering from a new order
 *        Actualise all elmts array according to the new numbering
 *        Connectivities
 *
 * \param [in, out] pmn              Pointer to \ref PDM_part_mesh_nodal_elmts_t instance
 * \param [in]      new_to_old_order New order for each part (size=n_part)
 *
 */
void
PDM_part_mesh_nodal_elmts_reorder_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old_order
);


void
PDM_part_mesh_nodal_elmts_reorder_vtx
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                   **new_to_old_order_vtx
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_RENUM_H__ */
