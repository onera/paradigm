/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_compare_operator.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"

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

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/* Tables for decomposing standard elments in edges and faces */
// https://cgns.github.io/CGNS_docs_current/sids/conv.html

static const int std_elt_n_edge[] = {
   0, // PDM_MESH_NODAL_POINT
   1, // PDM_MESH_NODAL_BAR2
   3, // PDM_MESH_NODAL_TRIA3
   4, // PDM_MESH_NODAL_QUAD4
  -1, // PDM_MESH_NODAL_POLY_2D
   6, // PDM_MESH_NODAL_TETRA4
   8, // PDM_MESH_NODAL_PYRAMID5
   9, // PDM_MESH_NODAL_PRISM6
  12, // PDM_MESH_NODAL_HEXA8
  -1, // PDM_MESH_NODAL_POLY_3D
   2, // PDM_MESH_NODAL_BARHO
   3, // PDM_MESH_NODAL_TRIAHO
   4, // PDM_MESH_NODAL_QUADHO
   6, // PDM_MESH_NODAL_TETRAHO
   8, // PDM_MESH_NODAL_PYRAMIDHO
   9, // PDM_MESH_NODAL_PRISMHO
  12, // PDM_MESH_NODAL_HEXAHO
   1, // PDM_MESH_NODAL_BARHO_BEZIER
   3  // PDM_MESH_NODAL_TRIAHO_BEZIER
};

static const int bar_edge_vtx[] = {
  0, 1
};

static const int tria_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0
};

static const int quad_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0
};

static const int tetra_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0,
  0, 3,
  1, 3,
  2, 3
};

static const int pyramid_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  0, 4,
  1, 4,
  2, 4,
  3, 4
};

static const int prism_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 0,
  0, 3,
  1, 4,
  2, 5,
  3, 4,
  4, 5,
  5, 3
};

static const int hexa_edge_vtx[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  0, 4,
  1, 5,
  2, 6,
  3, 7,
  4, 5,
  5, 6,
  6, 7,
  7, 4
};



static const int tria_face_vtx_idx[] = {
  0, 3
};

static const int tria_face_vtx[] = {
  0, 1, 2
};

static const int quad_face_vtx_idx[] = {
  0, 4
};

static const int quad_face_vtx[] = {
  0, 1, 2, 3
};

static const int tetra_face_vtx_idx[] = {
  0, 3, 6, 9, 12
};

static const int tetra_face_vtx[] = {
  0, 2, 1,
  0, 1, 3,
  0, 3, 2,
  1, 2, 3
};

static const int pyramid_face_vtx_idx[] = {
  0, 4, 7, 10, 13, 16
};

static const int pyramid_face_vtx[] = {
  3, 2, 1, 0,
  4, 0, 1,
  4, 1, 2,
  4, 2, 3,
  4, 3, 0
};

static const int prism_face_vtx_idx[] = {
  0, 3, 6, 10, 14, 18
};

static const int prism_face_vtx[] = {
  2, 1, 0,
  4, 5, 3,
  5, 4, 1, 2,
  4, 3, 0, 1,
  3, 5, 2, 0
};

static const int hexa_face_vtx_idx[] = {
  0, 4, 8, 12, 16, 20, 24
};

static const int hexa_face_vtx[] = {
  3, 2, 1, 0,
  6, 7, 4, 5,
  4, 7, 3, 0,
  7, 6, 2, 3,
  2, 6, 5, 1,
  1, 5, 4, 0
};



static
void
_compute_keys
(
  int   n_entity1,
  int  *entity1_vtx_idx,
  int  *entity1_vtx,
  int **keys_out
)
{
  int *keys = NULL;
  PDM_malloc(keys, n_entity1, int);
  for (int i = 0; i < n_entity1; ++i) {
    int key = 0;
    for (int j = entity1_vtx_idx[i]; j < entity1_vtx_idx[i+1]; ++j) {
      key += entity1_vtx[j];
    }
    keys[i] = key;
  }
  *keys_out = keys;
}


static
void
_generate_entity_connectivity
(
  int          n_parent_decompose_entity2,
  int         *parent_decompose_entity2_vtx_idx,
  int         *parent_decompose_entity2_vtx,
  int         *parent_decompose_parent_entity2,
  int          n_child_decompose_entity2,
  int         *child_decompose_entity2_vtx_idx,
  int         *child_decompose_entity2_vtx,
  int        **out_child_to_parent_idx,
  int        **out_child_to_parent,
  int          pn_entity1,
  int         *pentity1_entity2_idx,
  int         *pentity1_entity2,
  int         *out_pn_entity2,
  int         *pentity2_vtx_idx,
  int         *pentity2_vtx,
  PDM_bool_t   compute_parent_child
)
{
  /* Hash table from key = sum of vtx */
  int *parent_keys = NULL;
  int *child_keys  = NULL;

  _compute_keys(n_parent_decompose_entity2,
                parent_decompose_entity2_vtx_idx,
                parent_decompose_entity2_vtx,
                &parent_keys);

  _compute_keys(n_child_decompose_entity2,
                child_decompose_entity2_vtx_idx,
                child_decompose_entity2_vtx,
                &child_keys);

  if (0 == 1) {
    PDM_log_trace_array_int(parent_keys, n_parent_decompose_entity2, "parent_keys ::");
    PDM_log_trace_array_int(child_keys , n_child_decompose_entity2 , "child_keys  ::");
  }

  int *order     = NULL;
  int *is_solved = NULL;
  PDM_malloc(order    , n_parent_decompose_entity2, int);
  PDM_malloc(is_solved, n_parent_decompose_entity2, int);

  for (int i = 0; i < n_parent_decompose_entity2; ++i) {
    order    [i] = i;
    is_solved[i] = 0;
  }
  PDM_sort_int(parent_keys, order, n_parent_decompose_entity2);

  /* Identify keys in conflict */
  int n_conflit_to_solve = 0;
  int last_key = -1;

  int *key_conflict_idx = NULL;
  int *key_to_conflict  = NULL;
  PDM_malloc(key_conflict_idx, n_parent_decompose_entity2+1, int);
  PDM_malloc(key_to_conflict , n_parent_decompose_entity2  , int);

  key_conflict_idx[0] = 0;
  for (int i = 0; i < n_parent_decompose_entity2; ++i) {
    if (parent_keys[i] != last_key){
      key_conflict_idx[n_conflit_to_solve+1] = key_conflict_idx[n_conflit_to_solve]+1;
      n_conflit_to_solve++;
      last_key = parent_keys[i];
    } else {
      key_conflict_idx[n_conflit_to_solve]++;
    }
    key_to_conflict[i] = n_conflit_to_solve-1;
  }



  int n_max_entity_per_key = 0;
  for (int i = 0; i < n_conflit_to_solve; ++i) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, key_conflict_idx[i+1]-key_conflict_idx[i]);
  }
  for (int i = 0; i < n_child_decompose_entity2; ++i) {
    int n_vtx_child = child_decompose_entity2_vtx_idx[i+1] - child_decompose_entity2_vtx_idx[i];
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, n_vtx_child);
  }

  int *pentity1_entity2_n = NULL;
  if (compute_parent_child == PDM_TRUE) {
    pentity1_entity2_n = PDM_array_zeros_int(pn_entity1);
  }


  int  pn_entity2 = 0;
  int *child_to_parent_idx = NULL;
  int *child_to_parent     = NULL;
  PDM_malloc(child_to_parent_idx, n_child_decompose_entity2+1                     , int);
  PDM_malloc(child_to_parent    , n_child_decompose_entity2 * n_max_entity_per_key, int);


  /*
   * Solve conflicts
   */
  if (0 == 1) {
    PDM_log_trace_array_int(key_conflict_idx, n_conflit_to_solve, "key_conflict_idx ::  ");
    PDM_log_trace_array_int(key_to_conflict , n_parent_decompose_entity2, "key_to_conflict ::  ");
    for (int i = 0; i < n_conflit_to_solve; ++i) {
      log_trace(" ------ i = %i \n", i);
      for (int i_key = key_conflict_idx[i]; i_key < key_conflict_idx[i+1]; ++i_key) {
        int i_conflict = order[i_key];
        int beg = parent_decompose_entity2_vtx_idx[i_conflict];
        int n_vtx_in_entity = parent_decompose_entity2_vtx_idx[i_conflict+1] - beg;
        log_trace(" \t i_key = %i | beg = %i / n = %i : ", parent_keys[i_key], beg, n_vtx_in_entity);
        PDM_log_trace_array_int(parent_decompose_entity2_vtx + beg, n_vtx_in_entity, "");
      }
    }
    log_trace("-----------------------------------\n\n\n");
  }

  int *sens_entity      = NULL;
  int *same_entity_num  = NULL;
  PDM_malloc(sens_entity    , n_max_entity_per_key, int);
  PDM_malloc(same_entity_num, n_max_entity_per_key, int);

  /* If we are only interested in the link child->parent */
  child_to_parent_idx[0] = 0;
  for (int i = 0; i < n_child_decompose_entity2; ++i) {

    int child_key = child_keys[i];

    child_to_parent_idx[i+1] = child_to_parent_idx[i];

    int pos_key = PDM_binary_search_int(child_key, parent_keys, n_parent_decompose_entity2);
    assert(pos_key >= 0);
    int i_conflict = key_to_conflict[pos_key];

    int beg_1       = child_decompose_entity2_vtx_idx[i  ];
    int n_vtx_child = child_decompose_entity2_vtx_idx[i+1] - beg_1;


    /* Brute force */
    int n_conflict_entitys = key_conflict_idx[i_conflict+1] - key_conflict_idx[i_conflict];

    int idx_next_same_entity = 0;
    sens_entity    [idx_next_same_entity] = 1;
    same_entity_num[idx_next_same_entity++] = i_conflict;
    // PDM_log_trace_array_int(&child_decompose_entity2_vtx[beg_1], n_vtx_child, "child   : ");

    for (int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {

      int idx_decompose_entity1 = order[key_conflict_idx[i_conflict]+idx_entity];

      int beg_2        = parent_decompose_entity2_vtx_idx[idx_decompose_entity1  ];
      int n_vtx_parent = parent_decompose_entity2_vtx_idx[idx_decompose_entity1+1] - beg_2;

      // PDM_log_trace_array_int(&parent_decompose_entity2_vtx[beg_2], n_vtx_parent, "  parent : ");

      int is_same_entity = PDM_compare_unsigned_ordered_nuplets_int(n_vtx_child,
                                                                    &child_decompose_entity2_vtx[beg_1],
                                                                    n_vtx_parent,
                                                                    &parent_decompose_entity2_vtx[beg_2]);
      // log_trace("  is_same_entity = %d\n", is_same_entity);
      if (is_same_entity != 0) {
        sens_entity    [idx_next_same_entity  ] = is_same_entity;
        same_entity_num[idx_next_same_entity++] = idx_decompose_entity1;

        int i_entity1 = parent_decompose_parent_entity2[idx_decompose_entity1]-1;

        // Keep link between both
        child_to_parent[child_to_parent_idx[i+1]++] = i_entity1+1; // Attention il peut avoir le parent_num
        is_solved[idx_decompose_entity1] = 1;
      }
    } /* End conflict */

    /* Normalement on a trouvé mais on check quand même :p */
    assert(idx_next_same_entity > 1);

    if (compute_parent_child == PDM_TRUE) {
      /* Generate new entity2 and append in connectivity entity2_vtx*/
      pentity2_vtx_idx[pn_entity2+1] = pentity2_vtx_idx[pn_entity2];
      for (int idx_vtx = 0; idx_vtx < n_vtx_child; ++idx_vtx) {
        pentity2_vtx[pentity2_vtx_idx[pn_entity2+1]++] = child_decompose_entity2_vtx[beg_1+idx_vtx];
      }

      /* Apply sens for all entity found - First is excluded because it's BC !!!! */
      for (int k = 1; k < idx_next_same_entity; ++k) {
        int idx_decompose_entity = same_entity_num[k];
        int t_entity             = parent_decompose_parent_entity2[idx_decompose_entity]-1;
        int sign                 = sens_entity[k];

        int idx_write = pentity1_entity2_idx[t_entity] + pentity1_entity2_n[t_entity]++;
        pentity1_entity2[idx_write] = sign * (pn_entity2+1); // Vient de l'exterieur

        is_solved[idx_decompose_entity] = 1;
      }

      // Go to next
      pn_entity2++;
    }
  }

  PDM_realloc(child_to_parent, child_to_parent, child_to_parent_idx[n_child_decompose_entity2], int);

  if (0 == 1) {
    PDM_log_trace_connectivity_int(child_to_parent_idx, child_to_parent, n_child_decompose_entity2, "child_to_parent ::");
    PDM_log_trace_connectivity_int(pentity2_vtx_idx, pentity2_vtx, pn_entity2, "pentity2_vtx ::");
  }

  *out_child_to_parent_idx = child_to_parent_idx;
  *out_child_to_parent     = child_to_parent;

  if (compute_parent_child == PDM_FALSE) {
    *out_pn_entity2 = 0;
    PDM_free(key_conflict_idx);
    PDM_free(key_to_conflict );
    PDM_free(parent_keys     );
    PDM_free(child_keys      );
    PDM_free(order           );
    PDM_free(sens_entity     );
    PDM_free(same_entity_num );
    PDM_free(is_solved       );
    return;
  }


  /*
   * Si on veut la nouvelle connectivité
   *   - On utilise le child_to_parent et on marque "as solved" le conflict associé
   *   - On remplis directement le parent_child
   *   - Attention au parent_num
   */
  for (int i = 0; i < n_conflit_to_solve; ++i) {

    int n_conflict_entitys = key_conflict_idx[i+1] - key_conflict_idx[i];

    for (int j = 0; j < n_conflict_entitys; ++j ) {
      int idx_decompose_entity1 = order[key_conflict_idx[i]+j];
      if (is_solved[idx_decompose_entity1] == 1) {
        continue;
      }

      int beg_1  = parent_decompose_entity2_vtx_idx[idx_decompose_entity1  ];
      int n_vtx1 = parent_decompose_entity2_vtx_idx[idx_decompose_entity1+1] - beg_1;

      int idx_next_same_entity = 0;
      sens_entity    [idx_next_same_entity] = 1;
      same_entity_num[idx_next_same_entity++] = idx_decompose_entity1;

      for (int k = j+1; k < n_conflict_entitys; ++k) {

        int idx_decompose_entity2 = order[key_conflict_idx[i]+k];
        if (is_solved[idx_decompose_entity2] == 1) {
          continue;
        }

        int beg_2  = parent_decompose_entity2_vtx_idx[idx_decompose_entity2  ];
        int n_vtx2 = parent_decompose_entity2_vtx_idx[idx_decompose_entity2+1] - beg_2;

        int is_same_entity = PDM_compare_unsigned_ordered_nuplets_int(n_vtx1,
                                                                      &parent_decompose_entity2_vtx[beg_1],
                                                                      n_vtx2,
                                                                      &parent_decompose_entity2_vtx[beg_2]);

        if (is_same_entity != 0) {
          sens_entity    [idx_next_same_entity  ] = is_same_entity;
          same_entity_num[idx_next_same_entity++] = idx_decompose_entity2;
        }
      } /* End conflict */

      /* Generate new entity2 and append in connectivity entity2_vtx*/
      pentity2_vtx_idx[pn_entity2+1] = pentity2_vtx_idx[pn_entity2];
      for (int idx_vtx = 0; idx_vtx < n_vtx1; ++idx_vtx) {
        pentity2_vtx[pentity2_vtx_idx[pn_entity2+1]++] = parent_decompose_entity2_vtx[beg_1+idx_vtx];
      }

      /* Apply sens for all entity found */
      for (int k = 0; k < idx_next_same_entity; ++k) {
        int idx_decompose_entity = same_entity_num[k];
        int t_entity             = parent_decompose_parent_entity2[idx_decompose_entity]-1;
        int sign                 = sens_entity[k];

        int idx_write = pentity1_entity2_idx[t_entity] + pentity1_entity2_n[t_entity]++;
        pentity1_entity2[idx_write] = sign * (pn_entity2+1); // Vient de l'exterieur

        is_solved[idx_decompose_entity] = 1;
      }

      // Go to next
      pn_entity2++;

    }
  }

  if (compute_parent_child == PDM_TRUE) {
    PDM_free(pentity1_entity2_n);
  }

  *out_pn_entity2 = pn_entity2;

  PDM_free(key_conflict_idx);
  PDM_free(key_to_conflict );
  PDM_free(parent_keys     );
  PDM_free(child_keys      );
  PDM_free(order           );
  PDM_free(sens_entity     );
  PDM_free(same_entity_num );
  PDM_free(is_solved       );
}



/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_part_mesh_nodal_poly2d_decomposes_edges
(
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const PDM_g_num_t          *vtx_ln_to_gn,
 const int                  *connectivity_elmt_vtx,
 const int                  *connectivity_elmt_vtx_idx,
 const PDM_g_num_t          *elmt_ln_to_gn,
       int                  *elmt_edge_vtx_idx,
       PDM_g_num_t          *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       PDM_g_num_t          *elmt_edge_cell,
       int                  *parent_elmt_position
)
{

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int         *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  PDM_g_num_t *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int         *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int         *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;
  PDM_g_num_t *_elmt_edge_cell            = elmt_edge_cell       + _n_edge_current;

  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    // Reminder for poly2d -> Number of vertex = Number of edge
    int n_edge_elt = connectivity_elmt_vtx_idx[ielt+1] - connectivity_elmt_vtx_idx[ielt];
    *n_edge_current += n_edge_elt;
    int idx2 = connectivity_elmt_vtx_idx[ielt];
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[idx + 1] = _current_elmt_edge_vtx_idx[idx] + 2;
      _parent_elmt_position     [idx    ] = i_edge;

      _elmt_edge_cell           [ielt * n_edge_elt + i_edge    ] = elmt_ln_to_gn[ielt];

      _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

      int inext = (i_edge + 1) % n_edge_elt;
      _current_elmt_edge_vtx[2 * idx    ]  = vtx_ln_to_gn[connectivity_elmt_vtx[idx2 + i_edge]-1];
      _current_elmt_edge_vtx[2 * idx + 1]  = vtx_ln_to_gn[connectivity_elmt_vtx[idx2 + inext ]-1];

      idx += 1;
    }
  }


  *n_elt_current  += n_elt;

}


void
PDM_part_mesh_nodal_elmts_sections_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_face_vtx_idx,
  PDM_g_num_t                  *elmt_face_vtx,
  int                          *elmt_cell_face_idx,
  PDM_g_num_t                  *elmt_face_cell,
  int                          *parent_elmt_position
)
{
  PDM_UNUSED(elmt_cell_face_idx);

  /* Perform local decomposition */
  int  *local_n_decompose_elmt_face = NULL;
  int **local_elmt_face_idx         = NULL;
  int **local_elmt_face_vtx_idx     = NULL;
  int **local_elmt_face_vtx         = NULL;
  int **local_parent_elmt           = NULL;
  int **local_parent_elmt_position  = NULL;
  PDM_part_mesh_nodal_elmts_sections_local_decompose_faces(pmne,
                                                           &local_n_decompose_elmt_face,
                                                           &local_elmt_face_idx,
                                                           &local_elmt_face_vtx_idx,
                                                           &local_elmt_face_vtx,
                                                           &local_parent_elmt,
                                                           &local_parent_elmt_position);



  /* Concatenate partitions and translate to global IDs */
  int idx_face = 0;
  // int idx_cell = 0;
  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    PDM_g_num_t *cell_ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne,
                                                                               i_part,
                                                                               PDM_OWNERSHIP_BAD_VALUE);

    for (int i_face = 0; i_face < local_n_decompose_elmt_face[i_part]; i_face++) {
      elmt_face_vtx_idx[idx_face+1] = elmt_face_vtx_idx[idx_face];
      for (int idx_vtx = local_elmt_face_vtx_idx[i_part][i_face]; idx_vtx < local_elmt_face_vtx_idx[i_part][i_face+1]; idx_vtx++) {
        int i_vtx = local_elmt_face_vtx[i_part][idx_vtx] - 1;
        elmt_face_vtx[elmt_face_vtx_idx[idx_face+1]++] = vtx_ln_to_gn[i_part][i_vtx];
      }

      int i_cell = PDM_ABS (local_parent_elmt[i_part][i_face]) - 1;
      int sign   = PDM_SIGN(local_parent_elmt[i_part][i_face]);
      elmt_face_cell[idx_face] = sign * cell_ln_to_gn[i_cell];

      parent_elmt_position[idx_face] = local_parent_elmt_position[i_part][i_face];

      idx_face++;
    }

    PDM_free(local_elmt_face_idx       [i_part]);
    PDM_free(local_elmt_face_vtx_idx   [i_part]);
    PDM_free(local_elmt_face_vtx       [i_part]);
    PDM_free(local_parent_elmt         [i_part]);
    PDM_free(local_parent_elmt_position[i_part]);
  }
  PDM_free(local_n_decompose_elmt_face);
  PDM_free(local_elmt_face_idx        );
  PDM_free(local_elmt_face_vtx_idx    );
  PDM_free(local_elmt_face_vtx        );
  PDM_free(local_parent_elmt          );
  PDM_free(local_parent_elmt_position );
}


void
PDM_part_mesh_nodal_elmts_sections_decompose_edges
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  PDM_g_num_t                 **vtx_ln_to_gn,
  int                          *elmt_edge_vtx_idx,
  PDM_g_num_t                  *elmt_edge_vtx,
  int                          *elmt_cell_edge_idx,
  PDM_g_num_t                  *elmt_edge_cell,
  int                          *parent_elmt_position
)
{
  PDM_UNUSED(elmt_cell_edge_idx);

  /* Perform local decomposition */
  int  *local_n_decompose_elmt_edge = NULL;
  int **local_elmt_edge_idx         = NULL;
  int **local_elmt_edge_vtx_idx     = NULL;
  int **local_elmt_edge_vtx         = NULL;
  int **local_parent_elmt           = NULL;
  int **local_parent_elmt_position  = NULL;
  PDM_part_mesh_nodal_elmts_sections_local_decompose_edges(pmne,
                                                           &local_n_decompose_elmt_edge,
                                                           &local_elmt_edge_idx,
                                                           &local_elmt_edge_vtx_idx,
                                                           &local_elmt_edge_vtx,
                                                           &local_parent_elmt,
                                                           &local_parent_elmt_position);

  /* Concatenate partitions and translate to global IDs */
  int idx_edge = 0;
  // int idx_cell = 0;
  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    PDM_g_num_t *cell_ln_to_gn = PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne,
                                                                               i_part,
                                                                               PDM_OWNERSHIP_BAD_VALUE);

    for (int i_edge = 0; i_edge < local_n_decompose_elmt_edge[i_part]; i_edge++) {
      elmt_edge_vtx_idx[idx_edge+1] = elmt_edge_vtx_idx[idx_edge];
      for (int idx_vtx = local_elmt_edge_vtx_idx[i_part][i_edge]; idx_vtx < local_elmt_edge_vtx_idx[i_part][i_edge+1]; idx_vtx++) {
        int i_vtx = local_elmt_edge_vtx[i_part][idx_vtx] - 1;
        elmt_edge_vtx[elmt_edge_vtx_idx[idx_edge+1]++] = vtx_ln_to_gn[i_part][i_vtx];
      }

      int i_cell = PDM_ABS (local_parent_elmt[i_part][i_edge]) - 1;
      int sign   = PDM_SIGN(local_parent_elmt[i_part][i_edge]);
      elmt_edge_cell[idx_edge] = sign * cell_ln_to_gn[i_cell];

      parent_elmt_position[idx_edge] = local_parent_elmt_position[i_part][i_edge];

      idx_edge++;
    }

    PDM_free(local_elmt_edge_idx       [i_part]);
    PDM_free(local_elmt_edge_vtx_idx   [i_part]);
    PDM_free(local_elmt_edge_vtx       [i_part]);
    PDM_free(local_parent_elmt         [i_part]);
    PDM_free(local_parent_elmt_position[i_part]);
  }
  PDM_free(local_n_decompose_elmt_edge);
  PDM_free(local_elmt_edge_idx        );
  PDM_free(local_elmt_edge_vtx_idx    );
  PDM_free(local_elmt_edge_vtx        );
  PDM_free(local_parent_elmt          );
  PDM_free(local_parent_elmt_position );
}



void
PDM_part_mesh_nodal_elmts_decompose_faces_get_size
(
 PDM_part_mesh_nodal_elmts_t  *pmne,
 int                          *n_elt_tot,
 int                          *n_face_elt_tot,
 int                          *n_sum_vtx_face_tot,
 int                         **elmt_face_vtx_idx,
 int                         **elmt_cell_face_idx
)
{
  /* Get current structure to treat */
  *n_face_elt_tot     = 0;
  *n_sum_vtx_face_tot = 0;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {

      int n_face_elt     = PDM_n_face_elt_per_elmt    (pmne->sections_std[i_section]->t_elt);
      int n_sum_vtx_face = PDM_n_sum_vtx_face_per_elmt(pmne->sections_std[i_section]->t_elt);

      *n_elt_tot          += pmne->sections_std[i_section]->n_elt[i_part];
      *n_face_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_face_elt;
      *n_sum_vtx_face_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_face;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      // int _n_face      = pmne->sections_poly3d[i_section]->n_face[i_part];
      // *n_face_elt_tot += _n_face;
      // int n_face_vtx   = pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      // *n_sum_vtx_face_tot += n_face_vtx;
      int n_elt = pmne->sections_poly3d[i_section]->n_elt[i_part];
      int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
      int *cell_face     = pmne->sections_poly3d[i_section]->_cellfac    [i_part];
      int *face_vtx_idx  = pmne->sections_poly3d[i_section]->_facvtx_idx [i_part];
      for (int i = 0; i < n_elt; i++) {
        for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
          (*n_face_elt_tot)++;
          int face_id = PDM_ABS(cell_face[j]) - 1;
          *n_sum_vtx_face_tot += face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
        }
      }

      *n_elt_tot += n_elt;
    }

    // Not so sure about this...
    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face      = pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_face_elt_tot += _n_face;
      int n_edge_vtx   = pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_face_tot += n_edge_vtx;
      *n_elt_tot          += pmne->sections_poly2d[i_section]->n_elt[i_part];
    }

  }


  *elmt_face_vtx_idx  = PDM_array_zeros_int(*n_face_elt_tot + 1);
  *elmt_cell_face_idx = PDM_array_zeros_int(*n_elt_tot      + 1);


  int n_section    = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    int i_cell = 0;

    for (int isection = 0; isection < n_section; isection++) {

      int id_section = sections_id[isection];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              id_section,
                                                              i_part);


      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_error(__FILE__, __LINE__, 0, "TODO\n");
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int i_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
        int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
        for (int i = 0; i < n_elt; i++) {
          int icell = i;
          if (parent_num != NULL) {
            icell = parent_num[i];
          }
          (*elmt_cell_face_idx)[icell+1] = cell_face_idx[i+1] - cell_face_idx[i];
        }
      }
      else {
        int cell_face_n = PDM_n_face_elt_per_elmt(pmne->sections_std[isection]->t_elt);

        for (int i = 0; i < n_elt; i++) {
          (*elmt_cell_face_idx)[i_cell+1] = cell_face_n;
          i_cell++;
        }
      }

    }

  }

  PDM_array_accumulate_int(*elmt_cell_face_idx + 1, *n_elt_tot);


  int cell_face_vtx_n[6];


  for (int i_part = 0; i_part < pmne->n_part; i_part++) {

    int i_cell = 0;

    for (int isection = 0; isection < n_section; isection++) {


      int id_section = sections_id[isection];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_KEEP);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              id_section,
                                                              i_part);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_error(__FILE__, __LINE__, 0, "TODO\n");
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        int i_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;
        int *cell_face_idx = pmne->sections_poly3d[i_section]->_cellfac_idx[i_part];
        int *cell_face     = pmne->sections_poly3d[i_section]->_cellfac    [i_part];
        int *face_vtx_idx  = pmne->sections_poly3d[i_section]->_facvtx_idx [i_part];
        for (int i = 0; i < n_elt; i++) {
          int icell = i;
          if (parent_num != NULL) {
            icell = parent_num[i];
          }
          int cell_face_n = cell_face_idx[i+1] - cell_face_idx[i];
          for (int j = 0; j < cell_face_n; j++) {
            int face_id = PDM_ABS(cell_face[cell_face_idx[i]+j]) - 1;
            (*elmt_face_vtx_idx)[(*elmt_cell_face_idx)[icell]+j+1] = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
          }
        }
      }
      else {
        int cell_face_n = PDM_n_face_elt_per_elmt(pmne->sections_std[isection]->t_elt);

        switch (t_elt) {
          case PDM_MESH_NODAL_TRIA3:
          case PDM_MESH_NODAL_TRIAHO: {
            cell_face_vtx_n[0] = 3;
            break;
          }
          case PDM_MESH_NODAL_QUAD4:
          case PDM_MESH_NODAL_QUADHO: {
            cell_face_vtx_n[0] = 4;
            break;
          }
          case PDM_MESH_NODAL_TETRA4:
          case PDM_MESH_NODAL_TETRAHO: {
            for (int i = 0; i < 4; i++) {
              cell_face_vtx_n[i] = 3;
            }
            break;
          }
          case PDM_MESH_NODAL_PYRAMID5:
          case PDM_MESH_NODAL_PYRAMIDHO: {
            cell_face_vtx_n[0] = 4;
            cell_face_vtx_n[1] = 3;
            cell_face_vtx_n[2] = 3;
            cell_face_vtx_n[3] = 3;
            cell_face_vtx_n[4] = 3;
            break;
          }
          case PDM_MESH_NODAL_PRISM6:
          case PDM_MESH_NODAL_PRISMHO: {
            cell_face_vtx_n[0] = 3;
            cell_face_vtx_n[1] = 3;
            cell_face_vtx_n[2] = 4;
            cell_face_vtx_n[3] = 4;
            cell_face_vtx_n[4] = 4;
            break;
          }
          case PDM_MESH_NODAL_HEXA8:
          case PDM_MESH_NODAL_HEXAHO: {
            for (int i = 0; i < 6; i++) {
              cell_face_vtx_n[i] = 4;
            }
            break;
          }
          default:
            PDM_error(__FILE__, __LINE__, 0, "Invalid elt type %d\n", (int) t_elt);
        }

        for (int i = 0; i < n_elt; i++) {

          for (int j = 0; j < cell_face_n; j++) {
            (*elmt_face_vtx_idx)[(*elmt_cell_face_idx)[i_cell]+j+1] = cell_face_vtx_n[j];
          }

          i_cell++;
        }
      }

    }

  }


  PDM_array_accumulate_int(*elmt_face_vtx_idx + 1, *n_face_elt_tot);
  // printf("n_face_elt_tot     ::%i\n", *n_face_elt_tot   );
  // printf("n_sum_vtx_face_tot::%i\n" , *n_sum_vtx_face_tot);
}



void
PDM_part_mesh_nodal_elmts_decompose_edges_get_size
(
 PDM_part_mesh_nodal_elmts_t *pmne,
 int                         *n_elt_tot,
 int                         *n_edge_elt_tot,
 int                         *n_sum_vtx_edge_tot
)
{
  /* Get current structure to treat */
  *n_edge_elt_tot     = 0;
  *n_sum_vtx_edge_tot = 0;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {

      int n_edge_elt     = PDM_n_nedge_elt_per_elmt   (pmne->sections_std[i_section]->t_elt);
      int n_sum_vtx_edge = PDM_n_sum_vtx_edge_per_elmt(pmne->sections_std[i_section]->t_elt);

      *n_elt_tot          += pmne->sections_std[i_section]->n_elt[i_part];
      *n_edge_elt_tot     += pmne->sections_std[i_section]->n_elt[i_part] * n_edge_elt;
      *n_sum_vtx_edge_tot += pmne->sections_std[i_section]->n_elt[i_part] * n_sum_vtx_edge;

    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face          = pmne->sections_poly3d[i_section]->n_face[i_part];
      *n_edge_elt_tot     +=     pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_sum_vtx_edge_tot += 2 * pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      *n_elt_tot          += pmne->sections_poly3d[i_section]->n_elt[i_part];
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face          =     pmne->sections_poly2d[i_section]->n_elt[i_part];
      *n_edge_elt_tot     +=     pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_sum_vtx_edge_tot += 2 * pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      *n_elt_tot          +=     pmne->sections_poly2d[i_section]->n_elt[i_part];
    }
  }

  // printf("n_edge_elt_tot     ::%i\n", *n_edge_elt_tot   );
  // printf("n_sum_vtx_edge_tot::%i\n" , *n_sum_vtx_edge_tot);
}




/* Local decomposition functions */



void
PDM_part_mesh_nodal_std_decompose_local_edges
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *parent_num,
       int                  *elmt_edge_vtx_idx,
       int                  *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
)
{
  int parent_node_std[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = parent_node_std;
  } else {
    _parent_node = parent_node;
  }

  int n_edge_elt     = std_elt_n_edge[t_elt];
  int n_sum_vtx_edge = n_edge_elt * 2;
  int n_vtx_elt      = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);


  const int *elt_edge_vtx = NULL;
  switch (t_elt) {
    case PDM_MESH_NODAL_POINT: {
      return;
    }
    case PDM_MESH_NODAL_BAR2:
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      elt_edge_vtx = bar_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_TRIA3:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      elt_edge_vtx = tria_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_QUAD4:
    case PDM_MESH_NODAL_QUADHO: {
      elt_edge_vtx = quad_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_TETRA4:
    case PDM_MESH_NODAL_TETRAHO: {
      elt_edge_vtx = tetra_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_PYRAMID5:
    case PDM_MESH_NODAL_PYRAMIDHO: {
      elt_edge_vtx = pyramid_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_PRISM6:
    case PDM_MESH_NODAL_PRISMHO: {
      elt_edge_vtx = prism_edge_vtx;
      break;
    }
    case PDM_MESH_NODAL_HEXA8:
    case PDM_MESH_NODAL_HEXAHO: {
      elt_edge_vtx = hexa_edge_vtx;
      break;
    }
    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", t_elt);
    }
  }


  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  int *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int *_parent_elmt               = parent_elmt          + _n_edge_current;
  int *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;

  if (parent_num != NULL) {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = parent_num[i_elt] + 1;
      for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
        _parent_elmt[i_elt * n_edge_elt + i_edge] = i_parent;
      }
    }
  }
  else {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = _n_elt_current + i_elt + 1;
      for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
        _parent_elmt[i_elt * n_edge_elt + i_edge] = i_parent;
      }
    }
  }

  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _parent_elmt_position[i_elt * n_edge_elt + i_edge] = i_edge;
    }
    _elmt_cell_edge_idx[i_elt+1] = _elmt_cell_edge_idx[i_elt] + n_edge_elt;

    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[i_elt * n_edge_elt + i_edge + 1] = _current_elmt_edge_vtx_idx[i_elt * n_edge_elt + i_edge] + 2;
      for (int i = 2*i_edge; i < 2*(i_edge+1); i++) {
        int i_vtx = elt_edge_vtx[i];
        _current_elmt_edge_vtx[n_sum_vtx_edge*i_elt + i] = connectivity_elmt_vtx[n_vtx_elt * i_elt + _parent_node[i_vtx]];
      }
    }
  }

  *n_elt_current  += n_elt;
  *n_edge_current += n_elt * n_edge_elt;
}


void
PDM_part_mesh_nodal_poly2d_decompose_local_edges
(
       int                   n_elt,
       int                  *n_elt_current,
       int                  *n_edge_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *connectivity_elmt_vtx_idx,
 const int                  *parent_num,
       int                  *elmt_edge_vtx_idx,
       int                  *elmt_edge_vtx,
       int                  *elmt_cell_edge_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
)
{

  int _n_edge_current = *n_edge_current;
  int _n_elt_current  = *n_elt_current;
  int *_current_elmt_edge_vtx_idx = elmt_edge_vtx_idx    + _n_edge_current;
  int *_current_elmt_edge_vtx     = elmt_edge_vtx        + elmt_edge_vtx_idx[_n_edge_current];
  int *_parent_elmt_position      = parent_elmt_position + _n_edge_current;
  int *_parent_elmt               = parent_elmt          + _n_edge_current;
  int *_elmt_cell_edge_idx        = elmt_cell_edge_idx   + _n_elt_current;

  if (parent_num != NULL) {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = parent_num[i_elt] + 1;
      int n_edge_elt = connectivity_elmt_vtx_idx[i_elt+1] - connectivity_elmt_vtx_idx[i_elt];
      for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
        _parent_elmt[i_elt * n_edge_elt + i_edge] = i_parent;
      }
    }
  }
  else {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = _n_elt_current + i_elt + 1;
      int n_edge_elt = connectivity_elmt_vtx_idx[i_elt+1] - connectivity_elmt_vtx_idx[i_elt];
      for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
        _parent_elmt[i_elt * n_edge_elt + i_edge] = i_parent;
      }
    }
  }

  int idx = 0;
  for (int ielt = 0; ielt < n_elt; ielt++) {
    // Reminder for poly2d -> Number of vertex = Number of edge
    int n_edge_elt = connectivity_elmt_vtx_idx[ielt+1] - connectivity_elmt_vtx_idx[ielt];
    *n_edge_current += n_edge_elt;
    int idx2 = connectivity_elmt_vtx_idx[ielt];
    for (int i_edge = 0; i_edge < n_edge_elt; i_edge++) {
      _current_elmt_edge_vtx_idx[idx + 1] = _current_elmt_edge_vtx_idx[idx] + 2;
      _parent_elmt_position     [idx    ] = i_edge;

      _elmt_cell_edge_idx[ielt+1] = _elmt_cell_edge_idx[ielt] + n_edge_elt;

      int inext = (i_edge + 1) % n_edge_elt;
      _current_elmt_edge_vtx[2 * idx    ]  = connectivity_elmt_vtx[idx2 + i_edge];
      _current_elmt_edge_vtx[2 * idx + 1]  = connectivity_elmt_vtx[idx2 + inext ];

      idx += 1;
    }
  }

  *n_elt_current  += n_elt;
}


void
PDM_part_mesh_nodal_elmts_sections_local_decompose_edges
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                          **out_n_decompose_elmt_edge,
  int                         ***out_elmt_edge_idx,
  int                         ***out_elmt_edge_vtx_idx,
  int                         ***out_elmt_edge_vtx,
  int                         ***out_parent_elmt,
  int                         ***out_parent_elmt_position
)
{

  int  *n_decompose_elmt_edge = NULL;
  int **elmt_edge_idx         = NULL;
  int **elmt_edge_vtx_idx     = NULL;
  int **elmt_edge_vtx         = NULL;
  int **parent_elmt           = NULL;
  int **parent_elmt_position  = NULL;
  PDM_malloc(n_decompose_elmt_edge, pmne->n_part, int  );
  PDM_malloc(elmt_edge_idx        , pmne->n_part, int *);
  PDM_malloc(elmt_edge_vtx_idx    , pmne->n_part, int *);
  PDM_malloc(elmt_edge_vtx        , pmne->n_part, int *);
  PDM_malloc(parent_elmt          , pmne->n_part, int *);
  PDM_malloc(parent_elmt_position , pmne->n_part, int *);

  int  n_section  = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *section_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int parent_node[8];
  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    int n_elmt = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);

    int n_elmt_edge     = 0;
    int n_elmt_edge_vtx = 0;

    /* Count to evaluate size of elmt_edge_vtx and elmt_edge_vtx */
    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {
      int n_elt_section = pmne->sections_std[i_section]->n_elt[i_part];
      n_elmt_edge     += n_elt_section * PDM_n_nedge_elt_per_elmt   (pmne->sections_std[i_section]->t_elt);
      n_elmt_edge_vtx += n_elt_section * PDM_n_sum_vtx_edge_per_elmt(pmne->sections_std[i_section]->t_elt);
    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face = pmne->sections_poly3d[i_section]->n_face[i_part];
      n_elmt_edge     +=     pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      n_elmt_edge_vtx += 2 * pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face = pmne->sections_poly2d[i_section]->n_elt[i_part];
      n_elmt_edge     +=     pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
      n_elmt_edge_vtx += 2 * pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
    }

    int n_elt_current  = 0;
    int n_edge_current = 0;

    n_decompose_elmt_edge    [i_part] = n_elmt_edge;

    // printf("n_elmt          = %i \n", n_elmt);
    // printf("n_elmt_edge     = %i \n", n_elmt_edge);
    // printf("n_elmt_edge_vtx = %i \n", n_elmt_edge_vtx);

    elmt_edge_idx       [i_part] = NULL;
    elmt_edge_vtx_idx   [i_part] = NULL;
    elmt_edge_vtx       [i_part] = NULL;
    parent_elmt         [i_part] = NULL;
    parent_elmt_position[i_part] = NULL;

    PDM_malloc(elmt_edge_idx       [i_part], n_elmt + 1     , int);
    PDM_malloc(elmt_edge_vtx_idx   [i_part], n_elmt_edge + 1, int);
    PDM_malloc(elmt_edge_vtx       [i_part], n_elmt_edge_vtx, int);
    PDM_malloc(parent_elmt         [i_part], n_elmt_edge    , int);
    PDM_malloc(parent_elmt_position[i_part], n_elmt_edge    , int);

    elmt_edge_idx       [i_part][0] = 0;
    elmt_edge_vtx_idx   [i_part][0] = 0;
    for (int i_section = 0; i_section < n_section; i_section++) {
      int id_section = section_id[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_BAD_VALUE);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *elt_vtx_idx = NULL;
        int *elt_vtx     = NULL;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                     id_section,
                                                     i_part,
                                                     &elt_vtx_idx,
                                                     &elt_vtx,
                                                     PDM_OWNERSHIP_BAD_VALUE);

        PDM_part_mesh_nodal_poly2d_decompose_local_edges(n_elt,
                                                          &n_elt_current,
                                                          &n_edge_current,
                                                          elt_vtx,
                                                          elt_vtx_idx,
                                                          parent_num,
                                                          elmt_edge_vtx_idx   [i_part],
                                                          elmt_edge_vtx       [i_part],
                                                          elmt_edge_idx       [i_part],
                                                          parent_elmt         [i_part],
                                                          parent_elmt_position[i_part]);
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        PDM_error(__FILE__, __LINE__, 0, "Poly3d not handled yet\n");
      }
      else {

        int order = -1;

        int         *elt_vtx             = NULL;
        PDM_g_num_t *elt_ln_to_gn        = NULL;
        int         *_parent_num         = NULL;
        PDM_g_num_t *parent_entity_g_num = NULL;

        int *_parent_node = NULL;

        if (PDM_Mesh_nodal_elmt_is_ho(t_elt)) {
          const char *ho_ordering = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &elt_vtx,
                                                       &elt_ln_to_gn,
                                                       &_parent_num,
                                                       &parent_entity_g_num,
                                                       &order,
                                                       &ho_ordering,
                                                       PDM_OWNERSHIP_BAD_VALUE);
          PDM_Mesh_nodal_ho_parent_node(t_elt,
                                        order,
                                        ho_ordering,
                                        parent_node);
          _parent_node = parent_node;
        }
        else {
          order = 1;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    id_section,
                                                    i_part,
                                                    &elt_vtx,
                                                    &elt_ln_to_gn,
                                                    &_parent_num,
                                                    &parent_entity_g_num,
                                                    PDM_OWNERSHIP_BAD_VALUE);
        }

        PDM_part_mesh_nodal_std_decompose_local_edges(t_elt,
                                                      n_elt,
                                                      order,
                                                      _parent_node,
                                                      &n_elt_current,
                                                      &n_edge_current,
                                                      elt_vtx,
                                                      parent_num,
                                                      elmt_edge_vtx_idx   [i_part],
                                                      elmt_edge_vtx       [i_part],
                                                      elmt_edge_idx       [i_part],
                                                      parent_elmt         [i_part],
                                                      parent_elmt_position[i_part]);
      }

    } // End loop on sections

  } // End loop on parts

  *out_n_decompose_elmt_edge     = n_decompose_elmt_edge;
  *out_elmt_edge_idx             = elmt_edge_idx;
  *out_elmt_edge_vtx_idx         = elmt_edge_vtx_idx;
  *out_elmt_edge_vtx             = elmt_edge_vtx;
  *out_parent_elmt               = parent_elmt;
  *out_parent_elmt_position      = parent_elmt_position;

}




void
PDM_part_mesh_nodal_std_decompose_local_faces
(
       PDM_Mesh_nodal_elt_t  t_elt,
       int                   n_elt,
       int                   order,
       int                  *parent_node,
       int                  *n_elt_current,
       int                  *n_face_current,
 const int                  *connectivity_elmt_vtx,
 const int                  *parent_num,
       int                  *elmt_face_vtx_idx,
       int                  *elmt_face_vtx,
       int                  *elmt_cell_face_idx,
       int                  *parent_elmt,
       int                  *parent_elmt_position
)
{
  int parent_node_std[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  int *_parent_node;
  if (parent_node == NULL) {
    _parent_node = parent_node_std;
  } else {
    _parent_node = parent_node;
  }

  int n_face_elt     = 0;
  int n_sum_vtx_face = 0;
  int n_vtx_elt      = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

  const int *elt_face_vtx_idx = NULL;
  const int *elt_face_vtx     = NULL;

  switch (t_elt) {
    case PDM_MESH_NODAL_POINT:
    case PDM_MESH_NODAL_BAR2:
    case PDM_MESH_NODAL_BARHO:
    case PDM_MESH_NODAL_BARHO_BEZIER: {
      return;
    }
    case PDM_MESH_NODAL_TRIA3:
    case PDM_MESH_NODAL_TRIAHO:
    case PDM_MESH_NODAL_TRIAHO_BEZIER: {
      n_face_elt       = 1;
      n_sum_vtx_face   = 3;
      elt_face_vtx_idx = tria_face_vtx_idx;
      elt_face_vtx     = tria_face_vtx    ;
      break;
    }
    case PDM_MESH_NODAL_QUAD4:
    case PDM_MESH_NODAL_QUADHO: {
      n_face_elt       = 1;
      n_sum_vtx_face   = 4;
      elt_face_vtx_idx = quad_face_vtx_idx;
      elt_face_vtx     = quad_face_vtx    ;
      break;
    }
    case PDM_MESH_NODAL_TETRA4:
    case PDM_MESH_NODAL_TETRAHO: {
      n_face_elt       = 4;
      n_sum_vtx_face   = 4*3;
      elt_face_vtx_idx = tetra_face_vtx_idx;
      elt_face_vtx     = tetra_face_vtx    ;
      break;
    }
    case PDM_MESH_NODAL_PYRAMID5:
    case PDM_MESH_NODAL_PYRAMIDHO: {
      n_face_elt       = 5;
      n_sum_vtx_face   = 4*3 + 4;
      elt_face_vtx_idx = pyramid_face_vtx_idx;
      elt_face_vtx     = pyramid_face_vtx    ;
      break;
    }
    case PDM_MESH_NODAL_PRISM6:
    case PDM_MESH_NODAL_PRISMHO: {
      n_face_elt       = 5;
      n_sum_vtx_face   = 2*3 + 3*4;
      elt_face_vtx_idx = prism_face_vtx_idx;
      elt_face_vtx     = prism_face_vtx    ;
      break;
    }
    case PDM_MESH_NODAL_HEXA8:
    case PDM_MESH_NODAL_HEXAHO: {
      n_face_elt       = 6;
      n_sum_vtx_face   = 6*4;
      elt_face_vtx_idx = hexa_face_vtx_idx;
      elt_face_vtx     = hexa_face_vtx    ;
      break;
    }
    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", t_elt);
    }
  }


  int _n_face_current = *n_face_current;
  int _n_elt_current  = *n_elt_current;
  int *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  int *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  int *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  int *_parent_elmt               = parent_elmt          + _n_face_current;
  int *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;


  if (parent_num != NULL) {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = parent_num[i_elt] + 1;
      for (int i_face = 0; i_face < n_face_elt; i_face++) {
        _parent_elmt_position[i_elt * n_face_elt + i_face] = i_parent;
      }
    }
  }
  else {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = *n_elt_current + i_elt + 1;
      for (int i_face = 0; i_face < n_face_elt; i_face++) {
        _parent_elmt_position[i_elt * n_face_elt + i_face] = i_parent;
      }
    }
  }

  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _parent_elmt_position[i_elt * n_face_elt + i_face] = i_face;
      _parent_elmt         [i_elt * n_face_elt + i_face] = i_elt+1;
    }
    _elmt_cell_face_idx[i_elt+1] = _elmt_cell_face_idx[i_elt] + n_face_elt;

    for (int i_face = 0; i_face < n_face_elt; i_face++) {
      _current_elmt_face_vtx_idx[i_elt * n_face_elt + i_face + 1] = _current_elmt_face_vtx_idx[i_elt * n_face_elt + i_face];
      for (int i = elt_face_vtx_idx[i_face]; i < elt_face_vtx_idx[i_face+1]; i++) {
        int i_vtx = elt_face_vtx[i];
        _current_elmt_face_vtx_idx[i_elt * n_face_elt + i_face + 1]++;
        _current_elmt_face_vtx[n_sum_vtx_face*i_elt + i] = connectivity_elmt_vtx[n_vtx_elt * i_elt + _parent_node[i_vtx]];
      }
    }
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt * n_face_elt;
}


void
PDM_part_mesh_nodal_poly2d_decompose_local_faces
(
       int  n_elt,
       int *n_elt_current,
       int *n_face_current,
 const int *elt_vtx_idx,
 const int *elt_vtx,
 const int *parent_num,
       int *elmt_face_vtx_idx,
       int *elmt_face_vtx,
       int *elmt_cell_face_idx,
       int *parent_elmt,
       int *parent_elmt_position
)
{
  int _n_face_current = *n_face_current;
  int _n_elt_current  = *n_elt_current;
  int *_current_elmt_face_vtx_idx = elmt_face_vtx_idx    + _n_face_current;
  int *_current_elmt_face_vtx     = elmt_face_vtx        + elmt_face_vtx_idx[_n_face_current];
  int *_parent_elmt_position      = parent_elmt_position + _n_face_current;
  int *_parent_elmt               = parent_elmt          + _n_face_current;
  int *_elmt_cell_face_idx        = elmt_cell_face_idx   + _n_elt_current;

  if (parent_num != NULL) {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = parent_num[i_elt] + 1;
      _parent_elmt[i_elt] = i_parent;
    }
  }
  else {
    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int i_parent = *n_elt_current + i_elt + 1;
      _parent_elmt[i_elt] = i_parent;
    }
  }

  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    _parent_elmt_position[i_elt] = 0;

    _elmt_cell_face_idx[i_elt+1] = _elmt_cell_face_idx[i_elt] + 1;
    _current_elmt_face_vtx_idx[i_elt+1] = _current_elmt_face_vtx_idx[i_elt];
    for (int i = elt_vtx_idx[i_elt]; i < elt_vtx_idx[i_elt+1]; i++) {
      _current_elmt_face_vtx[_current_elmt_face_vtx_idx[i_elt+1]++] = elt_vtx[i];
    }
  }

  *n_elt_current  += n_elt;
  *n_face_current += n_elt;
}


void
PDM_part_mesh_nodal_elmts_sections_local_decompose_faces
(
  PDM_part_mesh_nodal_elmts_t   *pmne,
  int                          **out_n_decompose_elmt_face,
  int                         ***out_elmt_face_idx,
  int                         ***out_elmt_face_vtx_idx,
  int                         ***out_elmt_face_vtx,
  int                         ***out_parent_elmt,
  int                         ***out_parent_elmt_position
)
{
  int  *n_decompose_elmt_face = NULL;
  int **elmt_face_idx         = NULL;
  int **elmt_face_vtx_idx     = NULL;
  int **elmt_face_vtx         = NULL;
  int **parent_elmt           = NULL;
  int **parent_elmt_position  = NULL;
  PDM_malloc(n_decompose_elmt_face, pmne->n_part, int  );
  PDM_malloc(elmt_face_idx        , pmne->n_part, int *);
  PDM_malloc(elmt_face_vtx_idx    , pmne->n_part, int *);
  PDM_malloc(elmt_face_vtx        , pmne->n_part, int *);
  PDM_malloc(parent_elmt          , pmne->n_part, int *);
  PDM_malloc(parent_elmt_position , pmne->n_part, int *);

  int  n_section  = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *section_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int parent_node[8];

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    int n_elmt = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);

    int n_elmt_face     = 0;
    int n_elmt_face_vtx = 0;

    /* Count to evaluate size of elmt_face_vtx and elmt_face_vtx */
    for (int i_section = 0; i_section < pmne->n_section_std; i_section++) {
      int n_elt_section = pmne->sections_std[i_section]->n_elt[i_part];
      n_elmt_face     += n_elt_section * PDM_n_face_elt_per_elmt    (pmne->sections_std[i_section]->t_elt);
      n_elmt_face_vtx += n_elt_section * PDM_n_sum_vtx_face_per_elmt(pmne->sections_std[i_section]->t_elt);
    }

    for (int i_section = 0; i_section < pmne->n_section_poly3d; i_section++) {
      int _n_face = pmne->sections_poly3d[i_section]->n_face[i_part];
      n_elmt_face     += pmne->sections_poly3d[i_section]->_cellfac_idx[i_part][_n_face];
      // n_elmt_face_vtx += pmne->sections_poly3d[i_section]->_facvtx_idx[i_part][_n_face];
      //???
    }

    for (int i_section = 0; i_section < pmne->n_section_poly2d; i_section++) {
      int _n_face = pmne->sections_poly2d[i_section]->n_elt[i_part];
      n_elmt_face     += _n_face;
      n_elmt_face_vtx += pmne->sections_poly2d[i_section]->_connec_idx[i_part][_n_face];
    }

    int n_elt_current  = 0;
    int n_face_current = 0;

    n_decompose_elmt_face[i_part] = n_elmt_face;

    elmt_face_idx       [i_part] = NULL;
    elmt_face_vtx_idx   [i_part] = NULL;
    elmt_face_vtx       [i_part] = NULL;
    parent_elmt         [i_part] = NULL;
    parent_elmt_position[i_part] = NULL;

    PDM_malloc(elmt_face_idx       [i_part], n_elmt + 1     , int);
    PDM_malloc(elmt_face_vtx_idx   [i_part], n_elmt_face + 1, int);
    PDM_malloc(elmt_face_vtx       [i_part], n_elmt_face_vtx, int);
    PDM_malloc(parent_elmt         [i_part], n_elmt_face    , int);
    PDM_malloc(parent_elmt_position[i_part], n_elmt_face    , int);


    elmt_face_idx    [i_part][0] = 0;
    elmt_face_vtx_idx[i_part][0] = 0;

    for (int i_section = 0; i_section < n_section; i_section++) {
      int id_section = section_id[i_section];
      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne, id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, id_section, i_part);

      int *parent_num = PDM_part_mesh_nodal_elmts_parent_num_get(pmne,
                                                                 id_section,
                                                                 i_part,
                                                                 PDM_OWNERSHIP_BAD_VALUE);

      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        int *elt_vtx_idx = NULL;
        int *elt_vtx     = NULL;
        PDM_part_mesh_nodal_elmts_section_poly2d_get(pmne,
                                                     id_section,
                                                     i_part,
                                                     &elt_vtx_idx,
                                                     &elt_vtx,
                                                     PDM_OWNERSHIP_KEEP);

        PDM_part_mesh_nodal_poly2d_decompose_local_faces(n_elt,
                                                         &n_elt_current,
                                                         &n_face_current,
                                                         elt_vtx,
                                                         elt_vtx_idx,
                                                         parent_num,
                                                         elmt_face_vtx_idx   [i_part],
                                                         elmt_face_vtx       [i_part],
                                                         elmt_face_idx       [i_part],
                                                         parent_elmt         [i_part],
                                                         parent_elmt_position[i_part]);
      }
      else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        PDM_error(__FILE__, __LINE__, 0, "Poly3d not handled yet\n");
      }
      else {

        int order = -1;

        int         *elt_vtx             = NULL;
        PDM_g_num_t *elt_ln_to_gn        = NULL;
        int         *_parent_num         = NULL;
        PDM_g_num_t *parent_entity_g_num = NULL;

        int *_parent_node = NULL;

        if (PDM_Mesh_nodal_elmt_is_ho(t_elt)) {
          const char *ho_ordering = NULL;
          PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                       id_section,
                                                       i_part,
                                                       &elt_vtx,
                                                       &elt_ln_to_gn,
                                                       &_parent_num,
                                                       &parent_entity_g_num,
                                                       &order,
                                                       &ho_ordering,
                                                       PDM_OWNERSHIP_KEEP);
          PDM_Mesh_nodal_ho_parent_node(t_elt,
                                        order,
                                        ho_ordering,
                                        parent_node);
          _parent_node = parent_node;
        }
        else {
          order = 1;
          PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                    id_section,
                                                    i_part,
                                                    &elt_vtx,
                                                    &elt_ln_to_gn,
                                                    &_parent_num,
                                                    &parent_entity_g_num,
                                                    PDM_OWNERSHIP_KEEP);
        }

        PDM_part_mesh_nodal_std_decompose_local_faces(t_elt,
                                                      n_elt,
                                                      order,
                                                      _parent_node,
                                                      &n_elt_current,
                                                      &n_face_current,
                                                      elt_vtx,
                                                      parent_num,
                                                      elmt_face_vtx_idx   [i_part],
                                                      elmt_face_vtx       [i_part],
                                                      elmt_face_idx       [i_part],
                                                      parent_elmt         [i_part],
                                                      parent_elmt_position[i_part]);

      }

    } // End loop on sections

  } // End loop on parts

  *out_n_decompose_elmt_face = n_decompose_elmt_face;
  *out_elmt_face_idx         = elmt_face_idx;
  *out_elmt_face_vtx_idx     = elmt_face_vtx_idx;
  *out_elmt_face_vtx         = elmt_face_vtx;
  *out_parent_elmt           = parent_elmt;
  *out_parent_elmt_position  = parent_elmt_position;

}



void
PDM_part_mesh_nodal_elmts_compute_child_parent
(
  PDM_part_mesh_nodal_elmts_t   *pmne_parent,
  PDM_part_mesh_nodal_elmts_t   *pmne_child,
  PDM_bool_t                     compute_parent_child,
  int                         ***out_child_to_parent_idx,
  int                         ***out_child_to_parent,
  int                          **out_n_entity,
  int                         ***out_entity_to_vtx_idx,
  int                         ***out_entity_to_vtx,
  int                         ***out_parent_to_entity_idx,
  int                         ***out_parent_to_entity
)
{
  int dim_parent = pmne_parent->mesh_dimension;
  int dim_child  = pmne_child ->mesh_dimension;

  if (dim_parent <= dim_child) {
    PDM_error(__FILE__, __LINE__, 0, "Parent dimension must be greater than child dimension\n");
  }

  int n_part = pmne_parent->n_part;
  if (pmne_child->n_part != n_part) {
    PDM_error(__FILE__, __LINE__, 0, "Parent and child pmne must have same n_part\n");
  }


  int  *n_decompose_parent_entity             = NULL;
  int **decompose_parent_entity_idx           = NULL; // Size = n_elemt+1
  int **decompose_parent_entity_vtx_idx       = NULL; // Size = n_decompose_parent_entity_vtx+1
  int **decompose_parent_entity_vtx           = NULL;
  int **decompose_parent_parent_elmt          = NULL;
  int **decompose_parent_parent_elmt_position = NULL;

  int  *n_decompose_child_entity              = NULL;
  int **decompose_child_entity_idx            = NULL; // Size = n_elemt+1
  int **decompose_child_entity_vtx_idx        = NULL; // Size = n_decompose_child_entity_vtx+1
  int **decompose_child_entity_vtx            = NULL;
  int **decompose_child_parent_elmt           = NULL;
  int **decompose_child_parent_elmt_position  = NULL;

  void (*_decompose) (PDM_part_mesh_nodal_elmts_t   *,
                      int                          **,
                      int                         ***,
                      int                         ***,
                      int                         ***,
                      int                         ***,
                      int                         ***);

  switch (dim_child) {
    case 2: {
      _decompose = &PDM_part_mesh_nodal_elmts_sections_local_decompose_faces;
      break;
    }
    case 1: {
      _decompose = &PDM_part_mesh_nodal_elmts_sections_local_decompose_edges;
      break;
    }
    default : {
      PDM_error(__FILE__, __LINE__, 0, "Invalid dim_child %d\n", dim_child);
    }
  }

  _decompose(pmne_parent,
             &n_decompose_parent_entity,
             &decompose_parent_entity_idx,
             &decompose_parent_entity_vtx_idx,
             &decompose_parent_entity_vtx,
             &decompose_parent_parent_elmt,
             &decompose_parent_parent_elmt_position);

  _decompose(pmne_child,
             &n_decompose_child_entity,
             &decompose_child_entity_idx,
             &decompose_child_entity_vtx_idx,
             &decompose_child_entity_vtx,
             &decompose_child_parent_elmt,
             &decompose_child_parent_elmt_position);




  int **decompose_parent_entity = NULL;
  int  *pn_entity               = NULL;
  int **pentity_vtx_idx         = NULL;
  int **pentity_vtx             = NULL;
  int **child_to_parent_idx     = NULL;
  int **child_to_parent         = NULL;


  if (compute_parent_child == PDM_TRUE) {
    PDM_malloc(decompose_parent_entity, n_part, int *);
  }
  PDM_malloc(pn_entity          , n_part, int  );
  PDM_malloc(pentity_vtx_idx    , n_part, int *);
  PDM_malloc(pentity_vtx        , n_part, int *);
  PDM_malloc(child_to_parent_idx, n_part, int *);
  PDM_malloc(child_to_parent    , n_part, int *);

  for (int i_part = 0; i_part < n_part; i_part++) {

    int pn_parent = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_parent, i_part);

    int *_decompose_parent_entity = NULL;
    if (compute_parent_child == PDM_TRUE) {
      PDM_malloc(decompose_parent_entity[i_part], decompose_parent_entity_idx[i_part][pn_parent], int);
      _decompose_parent_entity = decompose_parent_entity[i_part];
    }
    PDM_malloc(pentity_vtx_idx[i_part], n_decompose_parent_entity      [i_part] + 1,                                int);
    PDM_malloc(pentity_vtx    [i_part], decompose_parent_entity_vtx_idx[i_part][n_decompose_parent_entity[i_part]], int);

    pentity_vtx_idx[i_part][0] = 0;
    _generate_entity_connectivity(n_decompose_parent_entity      [i_part],
                                  decompose_parent_entity_vtx_idx[i_part],
                                  decompose_parent_entity_vtx    [i_part],
                                  decompose_parent_parent_elmt   [i_part],
                                  n_decompose_child_entity       [i_part],
                                  decompose_child_entity_vtx_idx [i_part],
                                  decompose_child_entity_vtx     [i_part],
                                 &child_to_parent_idx            [i_part],
                                 &child_to_parent                [i_part],
                                  pn_parent,
                                  decompose_parent_entity_idx    [i_part],
                                  _decompose_parent_entity,
                                 &pn_entity                      [i_part],
                                  pentity_vtx_idx                [i_part],
                                  pentity_vtx                    [i_part],
                                  compute_parent_child);

    PDM_free(decompose_parent_entity_vtx_idx      [i_part]);
    PDM_free(decompose_parent_entity_vtx          [i_part]);
    PDM_free(decompose_parent_parent_elmt         [i_part]);
    PDM_free(decompose_parent_parent_elmt_position[i_part]);
    PDM_free(decompose_child_entity_idx           [i_part]);
    PDM_free(decompose_child_entity_vtx_idx       [i_part]);
    PDM_free(decompose_child_entity_vtx           [i_part]);
    PDM_free(decompose_child_parent_elmt          [i_part]); // if needed, pass as [out] arguments
    PDM_free(decompose_child_parent_elmt_position [i_part]); // if needed, pass as [out] arguments

    if (compute_parent_child == PDM_FALSE) {
      PDM_free(pentity_vtx_idx            [i_part]);
      PDM_free(pentity_vtx                [i_part]);
      PDM_free(decompose_parent_entity_idx[i_part]);
    }
    else {
      PDM_realloc(pentity_vtx_idx[i_part], pentity_vtx_idx[i_part], pn_entity[i_part] + 1                     , int);
      PDM_realloc(pentity_vtx    [i_part], pentity_vtx    [i_part], pentity_vtx_idx[i_part][pn_entity[i_part]], int);
    }
  }

  PDM_free(n_decompose_parent_entity            );
  PDM_free(decompose_parent_entity_vtx_idx      );
  PDM_free(decompose_parent_entity_vtx          );
  PDM_free(decompose_parent_parent_elmt         );
  PDM_free(decompose_parent_parent_elmt_position);
  PDM_free(n_decompose_child_entity             );
  PDM_free(decompose_child_entity_idx           );
  PDM_free(decompose_child_entity_vtx_idx       );
  PDM_free(decompose_child_entity_vtx           );
  PDM_free(decompose_child_parent_elmt          );
  PDM_free(decompose_child_parent_elmt_position );

  if (compute_parent_child == PDM_FALSE) {
    PDM_free(pn_entity                  );
    PDM_free(pentity_vtx_idx            );
    PDM_free(pentity_vtx                );
    PDM_free(decompose_parent_entity_idx);
    PDM_free(decompose_parent_entity    );
  }

  *out_n_entity             = pn_entity;
  *out_child_to_parent_idx  = child_to_parent_idx;
  *out_child_to_parent      = child_to_parent;
  *out_entity_to_vtx_idx    = pentity_vtx_idx;
  *out_entity_to_vtx        = pentity_vtx;
  *out_parent_to_entity_idx = decompose_parent_entity_idx;
  *out_parent_to_entity     = decompose_parent_entity;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
