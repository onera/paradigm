#ifndef __PDM_PART_MESH_NODAL_RENUM_H__
#define __PDM_PART_MESH_NODAL_RENUM_H__

/*----------------------------------------------------------------------------
 *  Standard headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * Numbering conventions used by all the functions below
 * -----------------------------------------------------
 *
 *  - "Part numbering" of the elements of a \ref PDM_part_mesh_nodal_elmts_t :
 *    index in [0, n_elmt[ used by \ref PDM_part_mesh_nodal_elmts_cell_vtx_connect_get,
 *    the groups (group_elmt, 1-based), the per-part global numbering cache and
 *    the element \ref PDM_part_comm_graph_t.
 *      * if parent_num == NULL : concatenation of the sections (in sections_id order)
 *      * if parent_num != NULL : given by parent_num (0-based)
 *
 *  - Storage is always section by section : sections are never mixed.
 *
 *  - Permutations are 0-based. new_to_old[i_new] = i_old.
 *    A field is updated in place with PDM_order_array(n, size, new_to_old, field).
 *
 *  - Per-element data that point *outward* (connec, numabs, parent_entity_g_num,
 *    elt_to_entity, cell_centers) are moved with the element.
 *    Data that point *to* the elements (groups, comm graph, parent_num values,
 *    num_elmt_parent_to_local) are remapped. The entries of each group are
 *    re-sorted by increasing element id (see the group permutation).
 *
 *  Only standard (and high-order) sections are supported : POLY_2D / POLY_3D raise an error.
 *
 */


/**
 * \brief Permute the elements inside each section of one partition according to the new_to_old numbering
 *
 * All the arrays attached to the elements of each section are moved
 * (connec, numabs, internal numabs, parent_num, parent_entity_g_num,
 *  elt_to_entity, cell_centers).
 * The *values* of parent_num are not modified and nothing outside the
 * sections is updated : prefer \ref PDM_part_mesh_nodal_elmts_reorder_elmts.
 *
 * \param [in] pmne                   Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in] i_part                 Partition identifier
 * \param [in] new_to_old_by_section  Permutation of each section, indexed by the
 *                                    position of the section in sections_id
 *                                    (size = n_section, each of size n_elt of the section)
 *
 */
void
PDM_part_mesh_nodal_elmts_reorder_sections_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                     i_part,
  int                         **new_to_old_by_section
);


/**
 * \brief Apply a permutation of the part numbering to the elements
 *
 * Sections are re-sorted by increasing new part number and :
 *  - if parent_num != NULL : parent_num values are remapped (any permutation is allowed)
 *  - if parent_num == NULL : the permutation must keep the section blocks
 *                            (error otherwise)
 * Groups are remapped and their entries re-sorted, the per-part global numbering
 * cache and num_elmt_parent_to_local are updated.
 * The comm graph is not known here : see \ref PDM_part_mesh_nodal_reorder_elmts.
 *
 * \param [in]  pmne              Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in]  new_to_old        Permutation of the part numbering, for each partition
 * \param [out] group_new_to_old  Permutation of the entries of each group
 *                                ([i_part][i_group], size n_group_elmt), or NULL if not wanted. The caller owns the arrays.
 *
 */
void
PDM_part_mesh_nodal_elmts_reorder_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old,
  int                       ****group_new_to_old
);


/**
 * \brief Remap the vertices in the element connectivities
 *
 * \param [in] pmne                  Pointer to \ref PDM_part_mesh_nodal_elmts_t object
 * \param [in] old_to_new_order_vtx  Vertex permutation (0-based), for each partition
 *
 */
void
PDM_part_mesh_nodal_elmts_reorder_vtx
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                   **old_to_new_order_vtx
);


/**
 * \brief Apply a permutation to the elements on selected geometry kind
 *
 * \param [in]  pmn               Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in]  geom_kind         Geometry kind
 * \param [in]  new_to_old        Permutation of the part numbering, for each partition
 * \param [out] group_new_to_old  Permutation of the entries of each group ([i_part][i_group]), or NULL if not wanted
 *
 */
void
PDM_part_mesh_nodal_reorder_elmts
(
  PDM_part_mesh_nodal_t    *pmn,
  PDM_geometry_kind_t       geom_kind,
  int                     **new_to_old,
  int                   ****group_new_to_old
);


/**
 * \brief Apply a permutation to the vertices
 *
 * Elements connectivities for all dimensions, coordinates, global numbering, parent numbering of the vertices are remapped
 * \note Vertice graph \ref PDM_part_comm_graph_t is updated if any.
 * \note The permutation is checked (error if it is not a permutation).
 *
 * \param [in] pmn                  Pointer to \ref PDM_part_mesh_nodal_t object
 * \param [in] new_to_old_order_vtx Vertex permutation (0-based), for each partition
 *
 */
void
PDM_part_mesh_nodal_reorder_vtx
(
  PDM_part_mesh_nodal_t  *pmn,
  int                   **new_to_old_order_vtx
);


/**
 * \brief Renumber partitioned nodal mesh entities
 *
 *  1. Elements of the principal dimension : reverse Cuthill-McKee (RCM) on the
 *     element graph (elements sharing a vertex). Sections are kept :
 *     inside each section the elements follow the RCM order.
 *  2. Vertices : order of first appearance while traversing the principal
 *     elements in RCM order, then the vertices only reached by the lower
 *     dimensions.
 *  3. Lower dimensions : elements sorted lexicographically by their sorted
 *     new vertex ids, sections kept.
 *
 * What this function guarantees :
 *    - global numbers are permuted, not regenerated (renumbering them is a
 *      separate service)
 *    - the vertex permutation covers all the vertices of the partition : a
 *      vertex referenced by no element at all (which a well-formed mesh does
 *      not have) is left at the end, in its initial order
 *    - the \ref PDM_part_comm_graph_t attached to \p pmn are updated
 *    - the permutations are stored in \p pmn and can be retrieved with
 *      \ref PDM_part_mesh_nodal_vtx_new_to_old_get,
 *      \ref PDM_part_mesh_nodal_elmts_new_to_old_get and
 *      \ref PDM_part_mesh_nodal_group_new_to_old_get (pdm_part_mesh_nodal.h)
 *
 * \note The stored permutations describe the last call only : a new call to
 *       this function replaces them.
 *
 * \note This is the only renumbering method available so far. A general entry
 *       point taking a method name per entity (as
 *       \ref PDM_multipart_set_reordering_options does for \ref PDM_multipart_t)
 *       is left for when a second method exists.
 *
 * \param [in]     pmn  Pointer to \ref PDM_part_mesh_nodal_t object
 *
 */

void
PDM_part_mesh_nodal_renum_cuthill
(
  PDM_part_mesh_nodal_t *pmn
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PART_MESH_NODAL_RENUM_H__ */
