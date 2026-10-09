/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_cuthill.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_order.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_renum.h"
#include "pdm_part_renum.h"
#include "pdm_priv.h"
#include "pdm_sort.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
PDM_ownership_t *
_init_array_at_ownership_keep
(
  const int n
)
{
  PDM_ownership_t *ownership = NULL;
  PDM_malloc(ownership, n, PDM_ownership_t);
  for (int i = 0; i < n; i++) {
    ownership[i] = PDM_OWNERSHIP_KEEP;
  }
  return ownership;
}


/*
 * Keep in pmn the permutations of one geometry kind (elements and group entries)
 */
static
void
_elmt_new_to_old_store
(
  PDM_part_mesh_nodal_t        *pmn,
  PDM_geometry_kind_t           geom_kind,
  int                         **elmt_new_to_old,
  int                        ***group_new_to_old
)
{
  int dim     = PDM_geometry_kind_to_dimension(geom_kind);
  int n_group = PDM_part_mesh_nodal_n_group_get(pmn, geom_kind);

  pmn->elmt_new_to_old          [dim] = elmt_new_to_old;
  pmn->elmt_new_to_old_ownership[dim] = _init_array_at_ownership_keep (pmn->n_part);

  pmn->n_group_new_to_old        [dim] = n_group;
  pmn->group_new_to_old          [dim] = group_new_to_old;
  PDM_malloc(pmn->group_new_to_old_ownership[dim], pmn->n_part, PDM_ownership_t *);
  for (int i_part = 0; i_part < pmn->n_part; i_part++) {
    pmn->group_new_to_old_ownership[dim][i_part] = _init_array_at_ownership_keep (n_group);
  }
}

static
void
_check_pmne_have_no_poly
(
  PDM_part_mesh_nodal_elmts_t *pmne
)
{
  if (pmne == NULL) {
    return;
  }

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            sections_id[i_section]);
    if (t_elt == PDM_MESH_NODAL_POLY_2D ||
        t_elt == PDM_MESH_NODAL_POLY_3D) {
      PDM_error("POLY_2D and POLY_3D sections are not supported yet\n");
    }
  }
}


/*
 * For one partition, locate each element of the part numbering :
 *   elmt_section[i_elmt] = position of its section in sections_id
 *   elmt_local  [i_elmt] = index of the element inside its section
 * section_elmt_idx (size n_section+1, allocated by the caller) : offsets of
 * the sections in the concatenated storage.
 *
 * Return 1 if the sections carry a parent_num, 0 otherwise.
 */
static
int
_elmt_location_compute
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                     i_part,
  int                          *section_elmt_idx,
  int                         **elmt_section,
  int                         **elmt_local
)
{
  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);
  int  n_elmt      = PDM_part_mesh_nodal_elmts_n_elmts_get    (pmne, i_part);

  /* parent_num must be defined on all the (non-empty) sections or on none */
  int has_parent_num = -1;
  section_elmt_idx[0] = 0;
  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = sections_id[i_section];
    int n_elt      = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);
    section_elmt_idx[i_section+1] = section_elmt_idx[i_section] + n_elt;

    if (n_elt == 0) {
      continue;
    }
    int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                               id_section,
                                                               i_part,
                                                               PDM_OWNERSHIP_BAD_VALUE);
    int lhas_parent_num = (parent_num != NULL);
    if (has_parent_num == -1) {
      has_parent_num = lhas_parent_num;
    }
    else if (has_parent_num != lhas_parent_num) {
      PDM_error("_elmt_location_compute : mix of sections with and without parent_num (i_part = %d)\n",
                   i_part);
    }
  }
  if (has_parent_num == -1) {
    has_parent_num = 0;
  }

  if (section_elmt_idx[n_section] != n_elmt) {
    PDM_error("_elmt_location_compute : inconsistent number of elements (%d / %d)\n",
                 section_elmt_idx[n_section], n_elmt);
  }

  int *_elmt_section = PDM_array_const_int(n_elmt, -1);
  int *_elmt_local   = NULL;
  PDM_malloc(_elmt_local, n_elmt, int);

  for (int i_section = 0; i_section < n_section; i_section++) {
    int id_section = sections_id[i_section];
    int n_elt      = section_elmt_idx[i_section+1] - section_elmt_idx[i_section];
    int *parent_num = NULL;
    if (has_parent_num && n_elt > 0) {
      parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                            id_section,
                                                            i_part,
                                                            PDM_OWNERSHIP_BAD_VALUE);
    }
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_elmt = (parent_num != NULL) ? parent_num[i_elt] : section_elmt_idx[i_section] + i_elt;
      if (i_elmt < 0 || i_elmt >= n_elmt || _elmt_section[i_elmt] != -1) {
        PDM_error("_elmt_location_compute : parent_num is not a permutation (i_part = %d)\n",
                     i_part);
      }
      _elmt_section[i_elmt] = i_section;
      _elmt_local  [i_elmt] = i_elt;
    }
  }

  *elmt_section = _elmt_section;
  *elmt_local   = _elmt_local;

  return has_parent_num;
}


/*
 * Order of the part numbering respecting the sections :
 *   - parent_num != NULL : the order is kept as is (storage stays per section)
 *   - parent_num == NULL : stable partition by section (the part numbering is
 *                          the concatenation of the sections)
 */
static
void
_section_stable_new_to_old
(
  PDM_part_mesh_nodal_elmts_t *pmne,
  const int                    i_part,
  const int                   *order,
  int                         *new_to_old
)
{
  int n_section = PDM_part_mesh_nodal_elmts_n_section_get(pmne);
  int n_elmt    = PDM_part_mesh_nodal_elmts_n_elmts_get  (pmne, i_part);

  int *section_elmt_idx = NULL;
  PDM_malloc(section_elmt_idx, n_section+1, int);
  int *elmt_section = NULL;
  int *elmt_local   = NULL;
  int has_parent_num = _elmt_location_compute(pmne,
                                              i_part,
                                              section_elmt_idx,
                                              &elmt_section,
                                              &elmt_local);

  if (has_parent_num) {
    memcpy(new_to_old, order, n_elmt * sizeof(int));
  }
  else {
    int *count = PDM_array_zeros_int(n_section);
    for (int i = 0; i < n_elmt; i++) {
      int i_old     = order[i];
      int i_section = elmt_section[i_old];
      new_to_old[section_elmt_idx[i_section] + count[i_section]++] = i_old;
    }
    PDM_free(count);
  }

  PDM_free(section_elmt_idx);
  PDM_free(elmt_section);
  PDM_free(elmt_local);
}


/*
 * Reverse Cuthill-McKee on the element graph (two elements are neighbours
 * if they share a vertex). elmt_vtx is 1-based.
 */
static
void
_elmts_rcm_compute
(
  const int   n_elmt,
  const int   n_vtx,
  int        *elmt_vtx_idx,
  int        *elmt_vtx,
  int        *new_to_old
)
{
  if (n_elmt == 0) {
    return;
  }

  int *vtx_elmt_idx = NULL;
  int *vtx_elmt     = NULL;
  PDM_connectivity_transpose(n_elmt,
                             n_vtx,
                             elmt_vtx_idx,
                             elmt_vtx,
                             &vtx_elmt_idx,
                             &vtx_elmt);

  int *elmt_elmt_idx = NULL;
  int *elmt_elmt     = NULL;
  PDM_combine_connectivity(n_elmt,
                           elmt_vtx_idx,
                           elmt_vtx,
                           vtx_elmt_idx,
                           vtx_elmt,
                           &elmt_elmt_idx,
                           &elmt_elmt);
  PDM_free(vtx_elmt_idx);
  PDM_free(vtx_elmt);

  /*
   * PDM_combine_connectivity keeps the element itself (it shares all its vertices) :
   * remove the diagonal and switch to 0-based (PDM_cuthill_generate expects a 0-based graph)
   */
  int idx_write = 0;
  int idx_beg   = elmt_elmt_idx[0];
  for (int i_elmt = 0; i_elmt < n_elmt; i_elmt++) {
    int idx_end = elmt_elmt_idx[i_elmt+1];
    for (int idx = idx_beg; idx < idx_end; idx++) {
      int j_elmt = PDM_ABS(elmt_elmt[idx]) - 1;
      if (j_elmt != i_elmt) {
        elmt_elmt[idx_write++] = j_elmt;
      }
    }
    idx_beg = idx_end;
    elmt_elmt_idx[i_elmt+1] = idx_write;
  }
  elmt_elmt_idx[0] = 0;

  PDM_cuthill_generate(n_elmt, elmt_elmt_idx, elmt_elmt, new_to_old);

  PDM_free(elmt_elmt_idx);
  PDM_free(elmt_elmt);

  /* Safety : the result must be a permutation */
  if (!PDM_array_is_permutation_int(new_to_old, n_elmt)) {
    PDM_error("_elmts_rcm_compute : PDM_cuthill_generate did not return a permutation\n");
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_part_mesh_nodal_elmts_reorder_sections_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                     i_part,
  int                         **new_to_old_by_section
)
{
  if (pmne == NULL) {
    return;
  }

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = sections_id[i_section];
    int n_elt      = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);
    if (n_elt == 0) {
      continue;
    }

    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

    /* NB : std_ho_get does not reset parent_num / parent_entity_g_num when they are NULL */
    int         *connec              = NULL;
    PDM_g_num_t *numabs              = NULL;
    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    int          order               = 1;
    const char  *ho_ordering         = NULL;
    PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                 id_section,
                                                 i_part,
                                                 &connec,
                                                 &numabs,
                                                 &parent_num,
                                                 &parent_entity_g_num,
                                                 &order,
                                                 &ho_ordering,
                                                 PDM_OWNERSHIP_BAD_VALUE);
    int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

    const int *new_to_old = new_to_old_by_section[i_section];

    PDM_order_array(n_elt, n_vtx_elt * sizeof(int), new_to_old, connec);

    if (numabs != NULL) {
      PDM_order_array(n_elt, sizeof(PDM_g_num_t), new_to_old, numabs);
    }
    if (parent_num != NULL) {
      PDM_order_array(n_elt, sizeof(int), new_to_old, parent_num);
    }
    if (parent_entity_g_num != NULL) {
      PDM_order_array(n_elt, sizeof(PDM_g_num_t), new_to_old, parent_entity_g_num);
    }

    int *elt_to_entity = PDM_part_mesh_nodal_elmts_section_elmt_to_entity_get(pmne,
                                                                              id_section,
                                                                              i_part,
                                                                              PDM_OWNERSHIP_BAD_VALUE);
    if (elt_to_entity != NULL) {
      PDM_order_array(n_elt, sizeof(int), new_to_old, elt_to_entity);
    }

    PDM_Mesh_nodal_block_std_t *block = pmne->sections_std[id_section - PDM_BLOCK_ID_BLOCK_STD];
    if (block->numabs_int != NULL && block->numabs_int[i_part] != NULL) {
      PDM_order_array(n_elt, sizeof(PDM_g_num_t), new_to_old, block->numabs_int[i_part]);
    }
    if (block->cell_centers != NULL && block->cell_centers[i_part] != NULL) {
      PDM_order_array(n_elt, 3 * sizeof(double), new_to_old, block->cell_centers[i_part]);
    }
  }
}


void
PDM_part_mesh_nodal_elmts_reorder_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old,
  int                       ****group_new_to_old
)
{
  if (pmne == NULL) {
    if (group_new_to_old != NULL) {
      *group_new_to_old = NULL;
    }
    return;
  }

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);
  int  n_group     = PDM_part_mesh_nodal_elmts_n_group_get    (pmne);

  /* Permutation of the entries of each group : identity unless the group is re-sorted */
  int ***_group_new_to_old = NULL;
  PDM_malloc(_group_new_to_old, pmne->n_part, int **);
  for (int i_part = 0; i_part < pmne->n_part; i_part++) {
    PDM_malloc(_group_new_to_old[i_part], n_group, int *);
    for (int i_group = 0; i_group < n_group; i_group++) {
      int          n_group_elmt   = 0;
      int         *group_elmt     = NULL;
      PDM_g_num_t *group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                          &n_group_elmt,
                                          &group_elmt,
                                          &group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);
      PDM_malloc(_group_new_to_old[i_part][i_group], n_group_elmt, int);
      for (int j = 0; j < n_group_elmt; j++) {
        _group_new_to_old[i_part][i_group][j] = j;
      }
    }
  }

  int *section_elmt_idx = NULL;
  PDM_malloc(section_elmt_idx, n_section+1, int);

  int **new_to_old_by_section = NULL;
  int **new_parent_by_section = NULL;
  PDM_malloc(new_to_old_by_section, n_section, int *);
  PDM_malloc(new_parent_by_section, n_section, int *);

  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    int n_elmt = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
    if (n_elmt == 0) {
      continue;
    }

    int *elmt_section = NULL;
    int *elmt_local   = NULL;
    int has_parent_num = _elmt_location_compute(pmne,
                                                i_part,
                                                section_elmt_idx,
                                                &elmt_section,
                                                &elmt_local);

    const int *_new_to_old = new_to_old[i_part];
    if (!PDM_array_is_permutation_int(_new_to_old, n_elmt)) {
      PDM_error("PDM_part_mesh_nodal_elmts_reorder_elmts : new_to_old is not a permutation (i_part = %d)\n", i_part);
    }
    int *old_to_new = PDM_array_new_inverse_permutation_int(_new_to_old, n_elmt);

    /*
     * Walk the new part numbering in increasing order : inside each section the
     * elements end up sorted by new part number
     */
    int *count = PDM_array_zeros_int(n_section);
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt = section_elmt_idx[i_section+1] - section_elmt_idx[i_section];
      PDM_malloc(new_to_old_by_section[i_section], n_elt, int);
      PDM_malloc(new_parent_by_section[i_section], n_elt, int);
    }

    for (int i_new = 0; i_new < n_elmt; i_new++) {
      int i_old     = _new_to_old[i_new];
      int i_section = elmt_section[i_old];
      int i_elt     = count[i_section]++;
      new_to_old_by_section[i_section][i_elt] = elmt_local[i_old];
      new_parent_by_section[i_section][i_elt] = i_new;
    }

    if (!has_parent_num) {
      /* The part numbering is the concatenation of the sections : blocks must be kept */
      for (int i_section = 0; i_section < n_section; i_section++) {
        int n_elt = section_elmt_idx[i_section+1] - section_elmt_idx[i_section];
        for (int i_elt = 0; i_elt < n_elt; i_elt++) {
          if (new_parent_by_section[i_section][i_elt] != section_elmt_idx[i_section] + i_elt) {
            PDM_error("PDM_part_mesh_nodal_elmts_reorder_elmts : without parent_num, "
                         "the permutation must keep the section blocks (i_part = %d, i_section = %d)\n",
                         i_part, i_section);
          }
        }
      }
    }

    /* Move the per-element data */
    PDM_part_mesh_nodal_elmts_reorder_sections_elmts(pmne,
                                                     i_part,
                                                     new_to_old_by_section);

    /* Remap parent_num values (after the move : position k holds the element whose new number is new_parent[k]) */
    if (has_parent_num) {
      for (int i_section = 0; i_section < n_section; i_section++) {
        int n_elt = section_elmt_idx[i_section+1] - section_elmt_idx[i_section];
        if (n_elt == 0) {
          continue;
        }
        int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                   sections_id[i_section],
                                                                   i_part,
                                                                   PDM_OWNERSHIP_BAD_VALUE);
        memcpy(parent_num, new_parent_by_section[i_section], n_elt * sizeof(int));
      }
    }

    /* Groups : 1-based element ids in the part numbering, kept sorted */
    for (int i_group = 0; i_group < n_group; i_group++) {
      int          n_group_elmt   = 0;
      int         *group_elmt     = NULL;
      PDM_g_num_t *group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                          &n_group_elmt,
                                          &group_elmt,
                                          &group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);
      if (n_group_elmt == 0 || group_elmt == NULL) {
        continue;
      }

      PDM_part_renum_array(n_group_elmt, old_to_new, group_elmt);

      int *order = _group_new_to_old[i_part][i_group];
      PDM_order_lnum_s(group_elmt, 1, order, n_group_elmt);
      PDM_order_array(n_group_elmt, sizeof(int), order, group_elmt);
      if (group_ln_to_gn != NULL) {
        PDM_order_array(n_group_elmt, sizeof(PDM_g_num_t), order, group_ln_to_gn);
      }
    }

    /* Cached global numbering by partition (part numbering) */
    if (pmne->numabs != NULL && pmne->numabs[i_part] != NULL) {
      PDM_order_array(n_elmt, sizeof(PDM_g_num_t), _new_to_old, pmne->numabs[i_part]);
    }

    /* Part number -> index inside the section */
    if (pmne->num_elmt_parent_to_local != NULL && pmne->num_elmt_parent_to_local[i_part] != NULL) {
      for (int i_section = 0; i_section < n_section; i_section++) {
        int n_elt = section_elmt_idx[i_section+1] - section_elmt_idx[i_section];
        for (int i_elt = 0; i_elt < n_elt; i_elt++) {
          pmne->num_elmt_parent_to_local[i_part][new_parent_by_section[i_section][i_elt]] = i_elt;
        }
      }
    }

    for (int i_section = 0; i_section < n_section; i_section++) {
      PDM_free(new_to_old_by_section[i_section]);
      PDM_free(new_parent_by_section[i_section]);
    }
    PDM_free(count);
    PDM_free(old_to_new);
    PDM_free(elmt_section);
    PDM_free(elmt_local);
  }

  PDM_free(new_to_old_by_section);
  PDM_free(new_parent_by_section);
  PDM_free(section_elmt_idx);

  if (group_new_to_old != NULL) {
    *group_new_to_old = _group_new_to_old;
  }
  else {
    for (int i_part = 0; i_part < pmne->n_part; i_part++) {
      for (int i_group = 0; i_group < n_group; i_group++) {
        PDM_free(_group_new_to_old[i_part][i_group]);
      }
      PDM_free(_group_new_to_old[i_part]);
    }
    PDM_free(_group_new_to_old);
  }
}


void
PDM_part_mesh_nodal_elmts_reorder_vtx
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                   **old_to_new_order_vtx
)
{
  if (pmne == NULL) {
    return;
  }

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_part = 0; i_part < pmne->n_part; i_part++) {
    for (int i_section = 0; i_section < n_section; i_section++) {

      int id_section = sections_id[i_section];
      int n_elt      = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);
      if (n_elt == 0) {
        continue;
      }

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int         *connec              = NULL;
      PDM_g_num_t *numabs              = NULL;
      int         *parent_num          = NULL;
      PDM_g_num_t *parent_entity_g_num = NULL;
      int          order               = 1;
      const char  *ho_ordering         = NULL;
      PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec,
                                                   &numabs,
                                                   &parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_BAD_VALUE);
      int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      PDM_part_renum_array(n_elt * n_vtx_elt, old_to_new_order_vtx[i_part], connec);
    }
  }
}


void
PDM_part_mesh_nodal_reorder_elmts
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  int                   **new_to_old,
  int                 ****group_new_to_old
)
{
  if (group_new_to_old != NULL) {
    *group_new_to_old = NULL;
  }

  if (pmn == NULL) {
    return;
  }

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);
  if (pmne == NULL) {
    return;
  }

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);

  PDM_part_mesh_nodal_elmts_reorder_elmts(pmne, new_to_old, group_new_to_old);

  PDM_part_comm_graph_t *pcg = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn, geom_kind, &pcg, PDM_OWNERSHIP_BAD_VALUE);

  /* Collective : the existence of the graph must be the same on all ranks */
  if (pcg != NULL) {
    int **old_to_new = NULL;
    PDM_malloc(old_to_new, n_part, int *);
    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_elmt = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
      /* new_to_old already checked by PDM_part_mesh_nodal_elmts_reorder_elmts */
      old_to_new[i_part] = PDM_array_new_inverse_permutation_int(new_to_old[i_part], n_elmt);
    }

    PDM_part_comm_graph_reorder(pcg, old_to_new);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(old_to_new[i_part]);
    }
    PDM_free(old_to_new);
  }
}


void
PDM_part_mesh_nodal_reorder_vtx
(
  PDM_part_mesh_nodal_t  *pmn,
  int                   **new_to_old_order_vtx
)
{
  if (pmn == NULL) {
    return;
  }

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);

  /* Check and invert the permutation */
  int **old_to_new_order_vtx = NULL;
  PDM_malloc(old_to_new_order_vtx, n_part, int *);
  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    if (!PDM_array_is_permutation_int(new_to_old_order_vtx[i_part], n_vtx)) {
      PDM_error("PDM_part_mesh_nodal_reorder_vtx : new_to_old is not a permutation (i_part = %d)\n", i_part);
    }
    old_to_new_order_vtx[i_part] = PDM_array_new_inverse_permutation_int(new_to_old_order_vtx[i_part], n_vtx);
  }

  /* Connectivities of all geometry kinds (corners included) */
  PDM_geometry_kind_t principal_geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);
  for (PDM_geometry_kind_t geom_kind = principal_geom_kind; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);
    PDM_part_mesh_nodal_elmts_reorder_vtx(pmne, (const int **) old_to_new_order_vtx);
  }

  /* Data attached to the vertices */
  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);

    double *coords = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    if (coords != NULL) {
      PDM_order_array(n_vtx, 3 * sizeof(double), new_to_old_order_vtx[i_part], coords);
    }

    PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    if (vtx_ln_to_gn != NULL) {
      PDM_order_array(n_vtx, sizeof(PDM_g_num_t), new_to_old_order_vtx[i_part], vtx_ln_to_gn);
    }

    int *vtx_parent = (int *) PDM_part_mesh_nodal_vertices_parent_get(pmn, i_part);
    if (vtx_parent != NULL) {
      PDM_order_array(n_vtx, sizeof(int), new_to_old_order_vtx[i_part], vtx_parent);
    }
  }

  /* Collective : the existence of the graph must be the same on all ranks */
  PDM_part_comm_graph_t *pcg_vtx = NULL;
  PDM_part_mesh_nodal_part_comm_graph_vtx_get(pmn, &pcg_vtx, PDM_OWNERSHIP_BAD_VALUE);
  if (pcg_vtx != NULL) {
    PDM_part_comm_graph_reorder(pcg_vtx, old_to_new_order_vtx);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(old_to_new_order_vtx[i_part]);
  }
  PDM_free(old_to_new_order_vtx);
}


void
PDM_part_mesh_nodal_renum
(
  PDM_part_mesh_nodal_t *pmn,
  const char            *renum_entity_method    [4],
  const int             *renum_entity_properties[4]
)
{
  PDM_UNUSED(renum_entity_method);
  PDM_UNUSED(renum_entity_properties);

  if (pmn == NULL) {
    return;
  }

  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);

  PDM_geometry_kind_t principal_geom_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);

  /* Permutations of a previous call are replaced */
  PDM_part_mesh_nodal_new_to_old_free(pmn);

  /* Check everything before modifying anything */
  for (PDM_geometry_kind_t geom_kind = principal_geom_kind; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);
    _check_pmne_have_no_poly(pmne);
  }

  PDM_part_mesh_nodal_elmts_t *pmne_principal = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                              principal_geom_kind);

  int **elmt_new_to_old = NULL;
  int **vtx_old_to_new  = NULL;
  int **vtx_new_to_old  = NULL;
  PDM_malloc(elmt_new_to_old, n_part, int *);
  PDM_malloc(vtx_old_to_new , n_part, int *);
  PDM_malloc(vtx_new_to_old , n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    vtx_old_to_new[i_part] = PDM_array_const_int(n_vtx, -1);
    PDM_malloc(vtx_new_to_old[i_part], n_vtx, int);
    int n_vtx_new = 0;

    elmt_new_to_old[i_part] = NULL;

    /*
     * Step 1 : RCM on the elements of the principal dimension (part numbering)
     * Step 2 : vertices by order of first appearance in the RCM order
     */
    if (pmne_principal != NULL) {
      int *elmt_vtx_idx = NULL;
      int *elmt_vtx     = NULL;
      int n_elmt = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne_principal,
                                                                  i_part,
                                                                  &elmt_vtx_idx,
                                                                  &elmt_vtx);

      int *rcm_new_to_old = NULL;
      PDM_malloc(rcm_new_to_old, n_elmt, int);
      _elmts_rcm_compute(n_elmt,
                         n_vtx,
                         elmt_vtx_idx,
                         elmt_vtx,
                         rcm_new_to_old);

      PDM_part_renum_by_first_appearance(n_elmt,
                                         rcm_new_to_old,
                                         elmt_vtx_idx,
                                         elmt_vtx,
                                         vtx_old_to_new[i_part],
                                         vtx_new_to_old[i_part],
                                         &n_vtx_new);

      /* Sections are never mixed */
      PDM_malloc(elmt_new_to_old[i_part], n_elmt, int);
      _section_stable_new_to_old(pmne_principal,
                                 i_part,
                                 rcm_new_to_old,
                                 elmt_new_to_old[i_part]);

      PDM_free(rcm_new_to_old);
      PDM_free(elmt_vtx_idx);
      PDM_free(elmt_vtx);
    }

    /* Vertices only reached by the lower dimensions */
    for (PDM_geometry_kind_t geom_kind = principal_geom_kind+1; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
      PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);
      if (pmne == NULL) {
        continue;
      }
      int *elmt_vtx_idx = NULL;
      int *elmt_vtx     = NULL;
      int n_elmt = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                                  i_part,
                                                                  &elmt_vtx_idx,
                                                                  &elmt_vtx);
      PDM_part_renum_by_first_appearance(n_elmt,
                                         NULL,
                                         elmt_vtx_idx,
                                         elmt_vtx,
                                         vtx_old_to_new[i_part],
                                         vtx_new_to_old[i_part],
                                         &n_vtx_new);
      PDM_free(elmt_vtx_idx);
      PDM_free(elmt_vtx);
    }

    /* Unreferenced vertices, in their initial order */
    for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
      if (vtx_old_to_new[i_part][i_vtx] == -1) {
        vtx_old_to_new[i_part][i_vtx]     = n_vtx_new;
        vtx_new_to_old[i_part][n_vtx_new] = i_vtx;
        n_vtx_new++;
      }
    }
    assert(n_vtx_new == n_vtx);
  }

  /* Apply and keep the permutations in pmn */
  if (pmne_principal != NULL) {
    int ***group_new_to_old = NULL;
    PDM_part_mesh_nodal_reorder_elmts(pmn,
                                      principal_geom_kind,
                                      elmt_new_to_old,
                                      &group_new_to_old);
    _elmt_new_to_old_store(pmn, principal_geom_kind, elmt_new_to_old, group_new_to_old);
  }
  else {
    PDM_free(elmt_new_to_old);
  }

  PDM_part_mesh_nodal_reorder_vtx(pmn, vtx_new_to_old);
  pmn->vtx_new_to_old           = vtx_new_to_old;
  pmn->vtx_new_to_old_ownership = _init_array_at_ownership_keep (n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(vtx_old_to_new[i_part]);
  }
  PDM_free(vtx_old_to_new);

  /*
   * Step 3 : lower dimensions, lexicographic order on their (new) sorted vertices
   */
  for (PDM_geometry_kind_t geom_kind = principal_geom_kind+1; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);
    if (pmne == NULL) {
      continue;
    }

    int **lower_new_to_old = NULL;
    PDM_malloc(lower_new_to_old, n_part, int *);

    for (int i_part = 0; i_part < n_part; i_part++) {
      int *elmt_vtx_idx = NULL;
      int *elmt_vtx     = NULL;
      int n_elmt = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne,
                                                                  i_part,
                                                                  &elmt_vtx_idx,
                                                                  &elmt_vtx);
      int *lex_new_to_old = NULL;
      PDM_malloc(lex_new_to_old, n_elmt, int);
      PDM_part_renum_by_lexicographic_order(n_elmt,
                                            elmt_vtx_idx,
                                            elmt_vtx,
                                            lex_new_to_old);

      PDM_malloc(lower_new_to_old[i_part], n_elmt, int);
      _section_stable_new_to_old(pmne,
                                 i_part,
                                 lex_new_to_old,
                                 lower_new_to_old[i_part]);

      PDM_free(lex_new_to_old);
      PDM_free(elmt_vtx_idx);
      PDM_free(elmt_vtx);
    }

    int ***group_new_to_old = NULL;
    PDM_part_mesh_nodal_reorder_elmts(pmn,
                                      geom_kind,
                                      lower_new_to_old,
                                      &group_new_to_old);
    _elmt_new_to_old_store(pmn, geom_kind, lower_new_to_old, group_new_to_old);
  }
}
