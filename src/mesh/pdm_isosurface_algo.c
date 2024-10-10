/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_to_part.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"
#include "pdm_unique.h"
#include "pdm_vtk.h"
#include "pdm_order.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_partitioning_algorithm.h"

#include "pdm_writer_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static const int tria_pairs[] = {
  1, 2,
  2, 0,
  0, 1
};

static const int tetra_pairs[] = {
  0, 1,
  0, 2,
  0, 3,
  1, 2,
  1, 3,
  2, 3
};

static const int *simplex_pairs[] = {
  tria_pairs,
  tetra_pairs
};

/* From split operator  >>> */
static int _pattern_permutation_2d[8] = {
  0, 0, 1, 0, 2, 2, 1, 0
};


static int _pattern_cell_vtx_permutation[64][4] = {
 {0, 1, 2, 3}, // 0
 {0, 1, 2, 3}, // 1
 {0, 2, 3, 1}, // 2
 {0, 1, 2, 3}, // 3
 {0, 3, 1, 2}, // 4
 {0, 3, 1, 2}, // 5
 {0, 2, 3, 1}, // 6
 {0, 1, 2, 3}, // 7
 {1, 2, 0, 3}, // 8
 {1, 2, 0, 3}, // 9
 {2, 0, 1, 3}, // 10
 {3, 0, 2, 1}, // 11
 {0, 1, 2, 3}, // 12
 {0, 1, 2, 3}, // 13
 {2, 0, 3, 1}, // 14 (! vol < 0)
 {0, 1, 2, 3}, // 15
 {1, 3, 2, 0}, // 16
 {1, 0, 3, 2}, // 17
 {0, 3, 1, 2}, // 18
 {0, 1, 3, 2}, // 19 (! vol < 0)
 {3, 1, 0, 2}, // 20
 {2, 0, 1, 3}, // 21
 {0, 3, 1, 2}, // 22
 {0, 3, 1, 2}, // 23
 {1, 3, 2, 0}, // 24
 {1, 2, 0, 3}, // 25
 {1, 2, 0, 3}, // 26
 {1, 2, 0, 3}, // 27
 {1, 3, 0, 2}, // 28 (! vol < 0)
 {1, 0, 3, 2}, // 29
 {0, 1, 2, 3}, // 30
 {0, 1, 2, 3}, // 31
 {2, 3, 0, 1}, // 32
 {0, 2, 3, 1}, // 33
 {2, 3, 0, 1}, // 34
 {0, 2, 3, 1}, // 35
 {3, 0, 2, 1}, // 36
 {0, 3, 2, 1}, // 37 (! vol < 0)
 {1, 0, 3, 2}, // 38
 {0, 2, 3, 1}, // 39
 {2, 1, 3, 0}, // 40
 {1, 2, 3, 0}, // 41 (! vol < 0)
 {2, 0, 1, 3}, // 42
 {2, 0, 1, 3}, // 43
 {2, 3, 0, 1}, // 44
 {1, 3, 2, 0}, // 45
 {2, 3, 0, 1}, // 46
 {2, 0, 1, 3}, // 47
 {3, 2, 1, 0}, // 48
 {1, 3, 2, 0}, // 49
 {2, 3, 1, 0}, // 50 (! vol < 0)
 {1, 2, 0, 3}, // 51
 {3, 0, 2, 1}, // 52
 {3, 1, 0, 2}, // 53
 {3, 0, 2, 1}, // 54
 {3, 0, 2, 1}, // 55
 {0, 1, 2, 3}, // 56
 {1, 3, 2, 0}, // 57
 {2, 1, 3, 0}, // 58
 {1, 2, 0, 3}, // 59
 {3, 2, 1, 0}, // 60
 {3, 1, 0, 2}, // 61
 {2, 3, 0, 1}, // 62
 {0, 1, 2, 3}  // 63
};

static int _pattern_cell_edge_permutation[64][6] = { // from split_operator
 { 1,  2,  3,  4,  5,  6}, // 0
 { 1,  2,  3,  4,  5,  6}, // 1
 { 2,  3,  1,  6, -4, -5}, // 2
 { 1,  2,  3,  4,  5,  6}, // 3
 { 3,  1,  2, -5, -6,  4}, // 4
 { 3,  1,  2, -5, -6,  4}, // 5
 { 2,  3,  1,  6, -4, -5}, // 6
 { 1,  2,  3,  4,  5,  6}, // 7
 { 4, -1,  5, -2,  6,  3}, // 8
 { 4, -1,  5, -2,  6,  3}, // 9
 {-2, -4,  6,  1,  3,  5}, // 10
 {-3, -6, -5,  2,  1, -4}, // 11
 { 1,  2,  3,  4,  5,  6}, // 12
 { 1,  2,  3,  4,  5,  6}, // 13
 {-2,  6, -4,  3,  1, -5}, // 14 (! vol < 0)
 { 1,  2,  3,  4,  5,  6}, // 15
 { 5,  4, -1, -6, -3, -2}, // 16
 {-1,  5,  4,  3,  2, -6}, // 17
 { 3,  1,  2, -5, -6,  4}, // 18
 { 1,  3,  2,  5,  4, -6}, // 19 (! vol < 0)
 {-5, -3, -6, -1,  4,  2}, // 20
 {-2, -4,  6,  1,  3,  5}, // 21
 { 3,  1,  2, -5, -6,  4}, // 22
 { 3,  1,  2, -5, -6,  4}, // 23
 { 5,  4, -1, -6, -3, -2}, // 24
 { 4, -1,  5, -2,  6,  3}, // 25
 { 4, -1,  5, -2,  6,  3}, // 26
 { 4, -1,  5, -2,  6,  3}, // 27
 { 5, -1,  4, -3, -6,  2}, // 28 (! vol < 0)
 {-1,  5,  4,  3,  2, -6}, // 29
 { 1,  2,  3,  4,  5,  6}, // 30
 { 1,  2,  3,  4,  5,  6}, // 31
 { 6, -2, -4, -3, -5,  1}, // 32
 { 2,  3,  1,  6, -4, -5}, // 33
 { 6, -2, -4, -3, -5,  1}, // 34
 { 2,  3,  1,  6, -4, -5}, // 35
 {-3, -6, -5,  2,  1, -4}, // 36
 { 3,  2,  1, -6, -5, -4}, // 37 (! vol < 0)
 {-1,  5,  4,  3,  2, -6}, // 38
 { 2,  3,  1,  6, -4, -5}, // 39
 {-4,  6, -2,  5, -1, -3}, // 40
 { 4,  5, -1,  6, -2, -3}, // 41 (! vol < 0)
 {-2, -4,  6,  1,  3,  5}, // 42
 {-2, -4,  6,  1,  3,  5}, // 43
 { 6, -2, -4, -3, -5,  1}, // 44
 { 5,  4, -1, -6, -3, -2}, // 45
 { 6, -2, -4, -3, -5,  1}, // 46
 {-2, -4,  6,  1,  3,  5}, // 47
 {-6, -5, -3, -4, -2, -1}, // 48
 { 5,  4, -1, -6, -3, -2}, // 49
 { 6, -4, -2, -5, -3, -1}, // 50 (! vol < 0)
 { 4, -1,  5, -2,  6,  3}, // 51
 {-3, -6, -5,  2,  1, -4}, // 52
 {-5, -3, -6, -1,  4,  2}, // 53
 {-3, -6, -5,  2,  1, -4}, // 54
 {-3, -6, -5,  2,  1, -4}, // 55
 { 1,  2,  3,  4,  5,  6}, // 56
 { 5,  4, -1, -6, -3, -2}, // 57
 {-4,  6, -2,  5, -1, -3}, // 58
 { 4, -1,  5, -2,  6,  3}, // 59
 {-6, -5, -3, -4, -2, -1}, // 60
 {-5, -3, -6, -1,  4,  2}, // 61
 { 6, -2, -4, -3, -5,  1}, // 62
 { 1,  2,  3,  4,  5,  6}  // 63
};

static const double ISOSURFACE_EPS_T = 1e-6; // Epsilon for snapping isovtx to vtx (ngon algo)

static inline int
_is_at_0_level(
  const double v,
  const double tol
)
{
  return (PDM_ABS(v) <= tol);
}


static inline int
_is_at_any_level
(
  const double v,
  const int    n_isovalues,
  const double isovalues[],
  const double tol
)
{
  int n_crossings = 0;
  for (int i = 0; i < n_isovalues; i++) {
    n_crossings += _is_at_0_level(v - isovalues[i], tol);
  }

  return n_crossings;
}




/*
 * \brief Compare two nuplets 
 * \return 1 if the nuplets are reversed up to a cyclic permutation, 0 otherwise
 */
static inline int
_compare_nuplets
(
  const int n1,
  const int i1[],
  const int n2,
  const int i2[]
)
{
  if (n1 != n2) {
    return 0;
  }

  for (int i = 0; i < n1; i++) {
    if (i2[i] == i1[0]) {
      for (int j = 1; j < n1; j++) {
        if (i2[(i+n1-j)%n1] != i1[j]) {
          return 0;
        }
      }
      return 1;
    }
  }

  return 0;
}


/**
 * \brief
 * Go though each extracted section of entry pmn, rebuild active edge (cross by iso or on iso) of elements
 * making edge shared by multiple elements unique. For each element of each section, link is kept
 * between the element edge and the iso edge lnum. For edges on iso, parent elements are kept.
 */
static void
_build_active_edges
(
  int                     n_section,
  PDM_Mesh_nodal_elt_t   *sections_elt_t,
  int                    *sections_n_elt,
  int                   **sections_elt_vtx,
  int                   **sections_elt_tag,
  const int               n_isovalues,
  const double           *isovalues,
  const double            tol,
  double                 *vtx_field,
  int                    *n_edge,
  int                   **edge_vtx,
  int                   **edge_bnd_tag_idx,
  int                   **edge_bnd_tag,
  int                   **edge_parent_idx,
  int                   **edge_parent,
  int                   **elt_edge,
  int                    *n_crossings
)
{
  int debug = 0;

  *n_crossings = 0;


  // > Prepare edge hash table
  const int max_key = 1024;
  int key_count[max_key];
  int key_idx  [max_key+1];
  PDM_array_reset_int(key_count, max_key, 0);


  /** 
   * First loop to count
   */
  for (int i_section = 0; i_section < n_section; i_section++) {

    // Get section type
    PDM_Mesh_nodal_elt_t t_elt = sections_elt_t[i_section];

    // Get table of vertex pairs for current element type
    const int *pairs = simplex_pairs[i_section];
    int n_pair = PDM_n_edge_elt_per_elmt(t_elt);

    // Get section information : number of elts and elt->vtx connectivity
    int  n_elt  = sections_n_elt  [i_section];
    int *connec = sections_elt_vtx[i_section];

    if (connec == NULL) {
      // Skip empty sections
      continue;
    }

    // Increment key occurrences
    int elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int *_connec = connec + elt_n_vtx*i_elt;
      for (int i_pair = 0; i_pair < n_pair; i_pair++) {
        int i_vtx0 = _connec[pairs[2*i_pair  ]];
        int i_vtx1 = _connec[pairs[2*i_pair+1]];

        if (_isosurface_cross_any_level(vtx_field[i_vtx0-1],
                                        vtx_field[i_vtx1-1],
                                        n_isovalues,
                                        isovalues,
                                        tol)) {
          // This is an active edge
          int key = (i_vtx0 + i_vtx1 - 2) % max_key;
          key_count[key]++;
        }
        else {
          for (int i = 0; i < n_isovalues; i++) {
            if (_is_at_0_level(vtx_field[i_vtx0-1] - isovalues[i], tol) &&
                _is_at_0_level(vtx_field[i_vtx1-1] - isovalues[i], tol)) {
              // This is an active edge with vertices on entry mesh vertices
              int key = (i_vtx0 + i_vtx1 - 2) % max_key;
              key_count[key]++;
            }
          }
        }
      } // End of loop on pairs
    } // End of loop on elements
  } // End of loop on sections


  /* Set up index and reset counter */
  // TODO: use array function
  key_idx[0] = 0;
  for (int i = 0; i < max_key; i++) {
    key_idx[i+1] = key_idx[i] + key_count[i];
    key_count[i] = 0;
  }


  /**
   * Second loop to fill
   */
  if (debug==1) log_trace("n_active_edges = %d\n", key_idx[max_key]);
  int *key_edge  = NULL;
  int *_edge_vtx = NULL;
  PDM_malloc( key_edge, key_idx[max_key]  , int);
  PDM_malloc(_edge_vtx, key_idx[max_key]*2, int);
  *n_edge = 0;

  int *_edge_count_parent  = PDM_array_zeros_int(key_idx[max_key]);
  int *_edge_count_bnd_tag = PDM_array_zeros_int(key_idx[max_key]);

  for (int i_section = 0; i_section < n_section; i_section++) {

    // Get section type
    PDM_Mesh_nodal_elt_t t_elt = sections_elt_t[i_section];

    // Get table of vertex pairs for current element type
    const int *pairs = simplex_pairs[i_section];
    int n_pair = PDM_n_edge_elt_per_elmt(t_elt);

    // Get section information : number of elts and elt->vtx connectivity
    int  n_elt  = sections_n_elt  [i_section];
    int *connec = sections_elt_vtx[i_section];

    int elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    elt_edge[i_section] = PDM_array_zeros_int(n_pair*n_elt);
    int has_bnd = sections_elt_tag[i_section]!=NULL;

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int is_bnd = 0;
      if (has_bnd) {
        is_bnd = sections_elt_tag[i_section][i_elt]!=0;
      }

      int *_connec = connec + elt_n_vtx*i_elt;
      for (int i_pair = 0; i_pair < n_pair; i_pair++) {
        int i_vtx0 = _connec[pairs[2*i_pair  ]];
        int i_vtx1 = _connec[pairs[2*i_pair+1]];

        int edge_id = 0;
        int key = (i_vtx0 + i_vtx1 - 2) % max_key;

        int _n_crossings = _isosurface_cross_any_level(vtx_field[i_vtx0-1],
                                                       vtx_field[i_vtx1-1],
                                                       n_isovalues,
                                                       isovalues,
                                                       tol);
        
        /* Look for possible collision */
        for (int i = 0; i < key_count[key]; i++) {
          int i_edge = key_edge[key_idx[key] + i];
          if (_edge_vtx[2*i_edge] == i_vtx0 && _edge_vtx[2*i_edge+1] == i_vtx1) {
            edge_id = i_edge+1;
            break;
          }
          else if (_edge_vtx[2*i_edge] == i_vtx1 && _edge_vtx[2*i_edge+1] == i_vtx0) {
            edge_id = -(i_edge+1);
            break;
          }
        }
        
        if (_n_crossings > 0) {

          if (edge_id == 0) {
            // Create new edge
            *n_crossings += _n_crossings;
            _edge_vtx[2*(*n_edge)  ] = i_vtx0;
            _edge_vtx[2*(*n_edge)+1] = i_vtx1;
            key_edge[key_idx[key] + key_count[key]++] = (*n_edge);
            edge_id = ++(*n_edge);
          }

        } else {
          for (int i = 0; i < n_isovalues; i++) {
            if (_is_at_0_level(vtx_field[i_vtx0-1] - isovalues[i], tol) &&
                _is_at_0_level(vtx_field[i_vtx1-1] - isovalues[i], tol)) {
              
              if (edge_id == 0) {
                _edge_vtx[2*(*n_edge)  ] = i_vtx0;
                _edge_vtx[2*(*n_edge)+1] = i_vtx1;
                key_edge[key_idx[key] + key_count[key]++] = (*n_edge);
                edge_id = ++(*n_edge);
              }
              if (t_elt==PDM_MESH_NODAL_TRIA3) {
                _edge_count_parent[PDM_ABS(edge_id)-1]++;
              }
              if (is_bnd==1) {
                _edge_count_bnd_tag[PDM_ABS(edge_id)-1]++;
              }
            }
          }
        } // End if active edge

        elt_edge[i_section][n_pair*i_elt+i_pair] = edge_id;
      } // End of loop on pairs
    } // End of loop on elements
    if (debug==1) {
      log_trace("i_section = %d\n", i_section);
      log_trace("\t (*n_edge) = %d\n", (*n_edge));
      PDM_log_trace_array_int(elt_edge[i_section], n_pair*n_elt      , "\t elt_edge[i_section] ::");
      PDM_log_trace_array_int(_edge_count_bnd_tag, key_idx[max_key]  , "\t edge_count_bnd_tag  ::");
      PDM_log_trace_array_int(_edge_count_parent , key_idx[max_key]  , "\t edge_count_parent   ::");
      PDM_log_trace_array_int(_edge_vtx          , (*n_edge)*2       , "\t edge_vtx            ::");
    }
  } // End of loop on sections
  PDM_free(key_edge);




  /**
   * Third loop to set edge parent
   */
  int *_edge_bnd_tag_idx = PDM_array_new_idx_from_sizes_int(_edge_count_bnd_tag, *n_edge);
  int n_bnd_tag_tot  = _edge_bnd_tag_idx[*n_edge];
  int *_edge_bnd_tag = PDM_array_zeros_int(n_bnd_tag_tot);
  
  int *_edge_parent_idx  = PDM_array_new_idx_from_sizes_int(_edge_count_parent , *n_edge);
  int n_parent_tot       = _edge_parent_idx [*n_edge];
  int *_edge_parent      = PDM_array_zeros_int(n_parent_tot);
  
  PDM_array_reset_int(_edge_count_bnd_tag, *n_edge, 0);
  PDM_array_reset_int(_edge_count_parent , *n_edge, 0);
  
  if (debug==1) {
    log_trace("n_bnd_tag_tot = %d\n", n_bnd_tag_tot);
    PDM_log_trace_array_int(_edge_count_bnd_tag,   key_idx[max_key], "_edge_count_bnd_tag ::");
    log_trace("n_parent_tot = %d\n", n_parent_tot);
    PDM_log_trace_array_int(_edge_count_parent,   key_idx[max_key], "_edge_count_parent ::");
    PDM_log_trace_array_int(_edge_parent_idx  , (*n_edge)+1       , "_edge_parent_idx   ::");
  }

  for (int i_section = 0; i_section < n_section; i_section++) {

    // Get section type
    PDM_Mesh_nodal_elt_t t_elt = sections_elt_t[i_section];

    // Get table of vertex pairs for current element type
    int n_pair = PDM_n_edge_elt_per_elmt(t_elt);

    // Get section information : number of elts and elt->vtx connectivity
    int n_elt = sections_n_elt[i_section];

    int has_bnd = sections_elt_tag[i_section]!=NULL;

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {

      int elt_tag = 0;
      if (has_bnd==1) {
        elt_tag = sections_elt_tag[i_section][i_elt];
      }

      for (int i_pair = 0; i_pair < n_pair; i_pair++) {
        int edge_id = PDM_ABS(elt_edge[i_section][n_pair*i_elt+i_pair]);

        // > Fill parent
        if (t_elt==PDM_MESH_NODAL_TRIA3) {
          if (edge_id!=0 && _edge_parent_idx[edge_id]-_edge_parent_idx[edge_id-1]!=0) {
            int i_write = _edge_parent_idx[edge_id-1] + _edge_count_parent[edge_id-1];
            // TODO: check if its ok after extract
            // if (parent_num) {
            //   _edge_parent[i_write] = parent_num[i_elt]+1;
            // }
            // else {
              _edge_parent[i_write] = i_elt+1;
              // _edge_parent_section[i_write] = i_section+1;
            // }
            _edge_count_parent[edge_id-1]++;
          } // End if has parent
        } // End TRIA_3
          
        // > Fill bnd tag
        if (edge_id!=0 && elt_tag!=0 && _edge_bnd_tag_idx[edge_id]-_edge_bnd_tag_idx[edge_id-1]!=0) {
          int i_write = _edge_bnd_tag_idx[edge_id-1] + _edge_count_bnd_tag[edge_id-1];
          _edge_bnd_tag[i_write] = elt_tag;
          _edge_count_bnd_tag[edge_id-1]++;
        } // End if has bnd tag
          
      } // End of loop on pairs
    } // End of loop on elements
   
    if (debug==1) {
      PDM_log_trace_array_int(_edge_count_bnd_tag , *n_edge     , "_edge_count_bnd_tag ::");
      PDM_log_trace_array_int(_edge_count_parent  , *n_edge     , "_edge_count_parent  ::");
      PDM_log_trace_array_int(_edge_parent        , n_parent_tot, "_edge_parent        ::");
      // PDM_log_trace_array_int(_edge_parent_section, n_parent_tot, "_edge_parent_section ::");
    }
   
  } // End of loop on sections
  

  

  // > Output
  *edge_bnd_tag_idx = _edge_bnd_tag_idx;
  *edge_bnd_tag     = _edge_bnd_tag;
  *edge_parent_idx  = _edge_parent_idx;
  *edge_parent      = _edge_parent;
  PDM_realloc(_edge_vtx, *edge_vtx, (*n_edge) * 2, int);

  // > Free tmp arrays
  PDM_free(_edge_count_parent);
  PDM_free(_edge_count_bnd_tag);
}


/**
 * \brief
 * Go though each extracted section of entry pmn, rebuild active face (on iso) of elements
 * making faces shared by multiple elements unique. For each element of each section, link is kept
 * between the element face and the iso face lnum. For faces on iso, parent elements are kept,
 * as well as the face connectivity in mesh numbering.
 */
static void
_build_active_faces
(
  int                     n_tetra,
  int                    *tetra_vtx,
  // PDM_g_num_t           **sections_elt_gnum,
  // int                   **sections_elt_tag,
  const int               n_isovalues,
  const double           *isovalues,
  const double            tol,
  double                 *vtx_field,
  int                    *vtx_to_iso_vtx,
  int                    *n_face,
  int                   **face_vtx_idx,
  int                   **face_vtx,
  int                   **face_parent_idx,
  int                   **face_parent,
  int                   **elt_face
)
{
  // TODO: simplify this routine since we decided that there will be only tetra ??

  int debug      = 0;
  // int debug_loop = 0;

  // > Prepare edge hash table
  const int max_key = 1024;
  int key_count[max_key];
  int key_idx  [max_key+1];
  PDM_array_reset_int(key_count, max_key, 0);

  int s_face_vtx = 0;
  int tetra_n_face_tot = 0;


  /**
   * First loop to count
   */

  // Get table of faces for current element type
  const int *_tetra_face_vtx_idx = NULL;
  const int *_tetra_face_vtx     = NULL;
  int _n_face = PDM_face_vtx_per_elmt(PDM_MESH_NODAL_TETRA4,
                                      &_tetra_face_vtx_idx,
                                      &_tetra_face_vtx);
  int elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TETRA4, 1);

  // Increment key occurrences
  for (int i_elt = 0; i_elt < n_tetra; i_elt++) {
    int *_connec = tetra_vtx + elt_n_vtx*i_elt;
    for (int i_face = 0; i_face < _n_face; i_face++) {
      for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
        int all_zero = 1;
        for (int i = _tetra_face_vtx_idx[i_face]; i < _tetra_face_vtx_idx[i_face+1]; i++) {
          int i_vtx = _connec[_tetra_face_vtx[i]] - 1;
          if (PDM_ABS(vtx_field[i_vtx] - isovalues[i_isovalue]) > tol) {
            all_zero = 0;
            break;
          }
        }

        if (all_zero) {
          // This is an active face
          int key = 0;
          for (int i = _tetra_face_vtx_idx[i_face]; i < _tetra_face_vtx_idx[i_face+1]; i++) {
            int id_vtx     = _connec[_tetra_face_vtx[i]];
            int id_vtx_iso = vtx_to_iso_vtx[id_vtx-1] - 1;
            assert(id_vtx_iso >= 0);
            key += id_vtx;
          }
          key = key % max_key;
          key_count[key]++;
          s_face_vtx += _tetra_face_vtx_idx[i_face+1] - _tetra_face_vtx_idx[i_face];
          break;
        } // End if all_zero
      } // End of loop on isovalues
      tetra_n_face_tot++;
    } // End of loop on faces
  } // End of loop on elements
  if (debug==1) log_trace("\tetra_n_face_tot = %d\n", tetra_n_face_tot);


  /* Set up index and reset counter */
  key_idx[0] = 0;
  for (int i = 0; i < max_key; i++) {
    key_idx[i+1] = key_idx[i] + key_count[i];
    key_count[i] = 0;
  }

  int *key_face = NULL;
  PDM_malloc(key_face, key_idx[max_key], int);

  if (debug==1) log_trace("n_active_faces = %d\n", key_idx[max_key]);
  if (debug==1) log_trace("face_vtx_size  = %d\n", s_face_vtx);
  int *_face_vtx_idx = NULL;
  int *_face_vtx     = NULL;
  PDM_malloc(_face_vtx_idx, key_idx[max_key] + 1, int);
  PDM_malloc(_face_vtx    , s_face_vtx          , int);
  int *_face_parent_count = PDM_array_zeros_int(key_idx[max_key]);

  _face_vtx_idx[0] = 0;



  /**
   * Second loop to fill
   */
  int tmp_face_vtx[4];

  int *_elt_face = PDM_array_zeros_int(tetra_n_face_tot);
  tetra_n_face_tot = 0;
  for (int i_elt = 0; i_elt < n_tetra; i_elt++) {
    // if (debug_loop) log_trace("i_elt = %d/%d\n", i_elt, n_elt);
    int *_connec = tetra_vtx + elt_n_vtx*i_elt;
    for (int i_face = 0; i_face < _n_face; i_face++) {
      int face_id = 0;
      for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
        int all_zero = 1;
        for (int i = _tetra_face_vtx_idx[i_face]; i < _tetra_face_vtx_idx[i_face+1]; i++) {
          int i_vtx = _connec[_tetra_face_vtx[i]] - 1;
          if (PDM_ABS(vtx_field[i_vtx] - isovalues[i_isovalue]) > tol) {
            all_zero = 0;
            break;
          }
        }

        if (all_zero) {
          int key = 0;
          int _n_vtx = _tetra_face_vtx_idx[i_face+1] - _tetra_face_vtx_idx[i_face];
          // if (debug_loop) log_trace("\t _n_vtx = %d\n", _n_vtx);
          for (int i = 0; i < _n_vtx; i++) {
            tmp_face_vtx[i] = _connec[_tetra_face_vtx[_tetra_face_vtx_idx[i_face]+i]];
            key += tmp_face_vtx[i];
          }
          key = key % max_key;

          /* Look for possible collision */
          // int face_id = 0;
          face_id = 0;
          for (int i = 0; i < key_count[key]; i++) {
            int j_face = key_face[key_idx[key] + i];

            // if (debug_loop) log_trace("\t _face_vtx_idx[j_face  ] = %d\n", _face_vtx_idx[j_face  ]);
            // if (debug_loop) log_trace("\t _face_vtx_idx[j_face+1] = %d\n", _face_vtx_idx[j_face+1]);
            if (_compare_nuplets(_n_vtx,
                                 tmp_face_vtx,
                                 _face_vtx_idx[j_face+1] - _face_vtx_idx[j_face],
                                 _face_vtx + _face_vtx_idx[j_face])) {
              face_id = -(j_face+1);
              break;
            }
          }

          // if (debug_loop) log_trace("\t (*n_face) = %d\n", (*n_face));
          if (face_id == 0) {
            // Create new face
            // TODO: make sure this face is properly oriented wrt the isosurface!! (but how??)
            _face_vtx_idx[(*n_face)+1] = _face_vtx_idx[(*n_face)] + _n_vtx;
            for (int i = 0; i < _n_vtx; i++) {
              // if (debug_loop) log_trace("\t               (*n_face)]   = %d\n",               (*n_face)   );
              // if (debug_loop) log_trace("\t _face_vtx_idx[(*n_face)]   = %d\n", _face_vtx_idx[(*n_face)]  );
              // if (debug_loop) log_trace("\t _face_vtx_idx[(*n_face)]+i = %d\n", _face_vtx_idx[(*n_face)]+i);
              _face_vtx[_face_vtx_idx[(*n_face)]+i] = tmp_face_vtx[i];
            }
            key_face[key_idx[key] + key_count[key]++] = (*n_face);
            face_id = ++(*n_face);

          }
          // if (debug_loop) log_trace("\t face_id = %d\n", face_id);

          _face_parent_count[PDM_ABS(face_id)-1]++;

          break;
        } // End if all_zero
      } // End of loop on isovalues

      _elt_face[tetra_n_face_tot++] = face_id;

    } // End of loop on faces
  } // End of loop on elements
  if (debug==1) {
    PDM_log_trace_array_int(_elt_face, tetra_n_face_tot, "elt_face ::");
  }
  PDM_free(key_face);

  if (debug==1) {
    int _face_vtx_size = _face_vtx_idx[key_idx[max_key]];
    PDM_log_trace_array_int(_face_parent_count, key_idx[max_key]  , "_face_parent_count ::");
    PDM_log_trace_array_int(_face_vtx_idx     , key_idx[max_key]+1, "_face_vtx_idx      ::");
    PDM_log_trace_array_int(_face_vtx         , _face_vtx_size    , "_face_vtx          ::");
  }


  /**
   * Third loop to set face parent
   */
  int *_face_parent_idx   = PDM_array_new_idx_from_sizes_int(_face_parent_count, *n_face);
  int n_parent_tot  = _face_parent_idx[*n_face];
  int *_face_parent = PDM_array_zeros_int(n_parent_tot);
  PDM_array_reset_int(_face_parent_count, *n_face, 0);
  if (debug==1) {
    log_trace("n_parent_tot = %d\n", n_parent_tot);
    PDM_log_trace_array_int(_face_parent_count, key_idx[max_key], "_face_parent_count ::");
    PDM_log_trace_array_int(_face_parent_idx  , *n_face+1       , "_face_parent_idx   ::");
  }

  tetra_n_face_tot = 0;
  for (int i_elt = 0; i_elt < n_tetra; i_elt++) {
    for (int i_face = 0; i_face < _n_face; i_face++) {
      int i_read_face = tetra_n_face_tot;
      int face_id = PDM_ABS(_elt_face[i_read_face]);
      if (face_id!=0 && _face_parent_idx[face_id]-_face_parent_idx[face_id-1] != 0) {
        int i_write_data = _face_parent_idx[face_id-1] + _face_parent_count[face_id-1];
        // TODO: check if its ok after extract
        // if (parent_num) {
        //   _face_parent[i_write_data] = parent_num[i_elt]+1;
        // }
        // else {
          _face_parent[i_write_data] = i_elt+1;
        // }
        _face_parent_count[face_id-1]++;
      } // End if has parent
      tetra_n_face_tot++;

    } // End of loop on pairs
  } // End of loop on elements


  if (debug==1) {
    PDM_log_trace_array_int(_face_parent_count, *n_face     , "_face_parent_count ::");
    PDM_log_trace_array_int(_face_parent      , n_parent_tot, "_face_parent       ::");
  }
  

  // > Output
  int _face_vtx_size = _face_vtx_idx[*n_face];
  *face_parent_idx = _face_parent_idx;
  *face_parent     = _face_parent;
  *elt_face        = _elt_face;
  PDM_realloc(_face_vtx_idx, *face_vtx_idx, *n_face+1     , int);
  PDM_realloc(_face_vtx    , *face_vtx    , _face_vtx_size, int);


  // > Free
  PDM_free(_face_parent_count);
}

static void
_build_iso_vtx
(
  const int           n_isovalues,
  const double       *isovalues,
  const double        tol,
  const int           n_vtx,
  const PDM_g_num_t  *vtx_gnum,
  const double       *vtx_coord,
  const double       *vtx_field,
  const int           n_edge,
  const int          *edge_vtx,
  const int           n_crossings,
        int         **out_vtx_to_iso_vtx,
        int          *out_iso_n_vtx,
        double      **out_iso_vtx_coord,
        int         **out_iso_vtx_parent_idx,
        int         **out_iso_vtx_parent,
        double      **out_iso_vtx_parent_weight,
        PDM_g_num_t **out_iso_vtx_parent_gnum,
        int         **out_iso_vtx_to_edge,
        int         **out_isovalue_vtx_idx
)
{
  int debug = 0;

  /**
   * Count number of iso vertices on entry mesh vertices
   * and allocate.
   */
  int n_vtx_on_vtx = 0;
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
      if (_is_at_0_level(vtx_field[i_vtx] - isovalues[i_isovalue], tol)) {
        n_vtx_on_vtx++;
      }
    }
  }

  int          iso_n_vtx             = n_vtx_on_vtx + n_crossings;
  double      *iso_vtx_coord         = NULL;
  int         *iso_vtx_parent_idx    = NULL;
  int         *iso_vtx_parent        = NULL;
  double      *iso_vtx_parent_weight = NULL;
  PDM_g_num_t *iso_vtx_parent_gnum   = NULL;
  PDM_malloc(iso_vtx_coord        , iso_n_vtx * 3               , double     );
  PDM_malloc(iso_vtx_parent_idx   , iso_n_vtx + 1               , int        );
  PDM_malloc(iso_vtx_parent       , n_vtx_on_vtx + 2*n_crossings, int        );
  PDM_malloc(iso_vtx_parent_weight, n_vtx_on_vtx + 2*n_crossings, double     );
  PDM_malloc(iso_vtx_parent_gnum  , iso_n_vtx * 3               , PDM_g_num_t);
  int *vtx_to_iso_vtx   = PDM_array_zeros_int(n_vtx);
  int *iso_vtx_to_edge  = PDM_array_zeros_int(iso_n_vtx);
  int *isovalue_vtx_idx = PDM_array_zeros_int(n_isovalues+1);
  iso_vtx_parent_idx[0] = 0;


  /**
   * Because all isovalues are on same mesh,
   * some entities will have multiple descendent
   * so we need to count them to generate unique gnum...
   */
  int *edge_n_child = NULL;
  PDM_calloc(edge_n_child, n_edge, int);


  /*
   * Fill iso vertices infos
   */
  iso_n_vtx = 0;

  for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
    double isovalue = isovalues[i_isovalue];

    // > Iso vertices on mesh vertices (one array for all isovalue is okay thx to isos->ISOSURFACE_EPS)
    for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
      if (_is_at_0_level(vtx_field[i_vtx] - isovalue, tol)) {
        vtx_to_iso_vtx[i_vtx] = iso_n_vtx+1;
        memcpy(iso_vtx_coord + 3*iso_n_vtx,
               vtx_coord     + 3*i_vtx,
               sizeof(double) * 3);
        iso_vtx_parent_idx   [                   iso_n_vtx+1] = iso_vtx_parent_idx[iso_n_vtx] + 1;
        iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx ]] = i_vtx + 1;
        iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx ]] = 1.;
        iso_vtx_parent_gnum  [                 3*iso_n_vtx  ] = vtx_gnum[i_vtx];
        iso_vtx_parent_gnum  [                 3*iso_n_vtx+1] = 0;
        iso_vtx_parent_gnum  [                 3*iso_n_vtx+2] = 0;
        iso_n_vtx++;
      }
    }
    
    // > Iso vertices on mesh edges
    for (int i_edge = 0; i_edge < n_edge; i_edge++) {
      int i_vtx0 = edge_vtx[2*i_edge  ] - 1;
      int i_vtx1 = edge_vtx[2*i_edge+1] - 1;

      double val0 = vtx_field[i_vtx0] - isovalue;
      double val1 = vtx_field[i_vtx1] - isovalue;

      if (_isosurface_cross_0_level(val0, val1, tol)) {
        double t = val0 / (val0 - val1);
        for (int i = 0; i < 3; i++) {
          iso_vtx_coord[3*iso_n_vtx+i] = (1-t)*vtx_coord[3*i_vtx0+i] + t*vtx_coord[3*i_vtx1+i];
        }
        iso_vtx_parent_idx[iso_n_vtx+1] = iso_vtx_parent_idx[iso_n_vtx] + 2;

        if (vtx_gnum[i_vtx0]>vtx_gnum[i_vtx1]) {
          iso_vtx_parent_gnum[3*iso_n_vtx  ] =   vtx_gnum[i_vtx1];
          iso_vtx_parent_gnum[3*iso_n_vtx+1] =   vtx_gnum[i_vtx0];
          iso_vtx_parent_gnum[3*iso_n_vtx+2] = --edge_n_child[i_edge];
          
          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]  ] = i_vtx1 + 1;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]  ] = t; 

          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]+1] = i_vtx0 + 1;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]+1] = 1. - t;
        }
        else {
          iso_vtx_parent_gnum[3*iso_n_vtx  ] =   vtx_gnum[i_vtx0];
          iso_vtx_parent_gnum[3*iso_n_vtx+1] =   vtx_gnum[i_vtx1];
          iso_vtx_parent_gnum[3*iso_n_vtx+2] = --edge_n_child[i_edge];
          
          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]  ] = i_vtx0 + 1;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]  ] = 1. - t;

          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]+1] = i_vtx1 + 1;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]+1] = t;
        }
        
        iso_vtx_to_edge[iso_n_vtx] = i_edge+1;

        iso_n_vtx++;
      }

    } // End of loop on edges
    isovalue_vtx_idx[i_isovalue+1] = iso_n_vtx;
  } // End of loop on isovalues

  free(edge_n_child);

  if (debug==1) {
    int n_parent_tot = iso_vtx_parent_idx[iso_n_vtx];
    log_trace("\n");
    PDM_log_trace_array_int (iso_vtx_parent_idx ,   iso_n_vtx    , "iso_vtx_parent_idx ::");
    PDM_log_trace_array_int (iso_vtx_parent     ,   n_parent_tot , "iso_vtx_parent     ::");
    PDM_log_trace_array_long(iso_vtx_parent_gnum, 3*iso_n_vtx    , "iso_vtx_parent_gnum::");
    PDM_log_trace_array_int (iso_vtx_to_edge    ,   iso_n_vtx    , "iso_vtx_to_edge    ::");
    PDM_log_trace_array_int (vtx_to_iso_vtx     ,   n_vtx        , "vtx_to_iso_vtx     ::");
    PDM_log_trace_array_int (isovalue_vtx_idx   ,   n_isovalues+1, "isovalue_vtx_idx   ::");
  }


  /*
   * Output
   */
  *out_vtx_to_iso_vtx        = vtx_to_iso_vtx;
  *out_iso_n_vtx             = iso_n_vtx;
  *out_iso_vtx_coord         = iso_vtx_coord;
  *out_iso_vtx_parent_idx    = iso_vtx_parent_idx;
  *out_iso_vtx_parent        = iso_vtx_parent;
  *out_iso_vtx_parent_weight = iso_vtx_parent_weight;
  *out_iso_vtx_parent_gnum   = iso_vtx_parent_gnum;
  *out_iso_vtx_to_edge       = iso_vtx_to_edge;
  *out_isovalue_vtx_idx      = isovalue_vtx_idx;
}



static void
_contouring_triangles
(
  int           n_elt,
  int          *elt_vtx,
  PDM_g_num_t  *elt_gnum,
  int          *elt_bnd_tag,
  int           n_edge,
  int          *edge_bnd_tag_idx,
  int          *edge_bnd_tag,
  int          *edge_parent_idx,
  int          *edge_parent,
  int          *elt_edge,
  PDM_g_num_t  *vtx_gnum,
  double       *vtx_field,
  int          *vtx_to_iso_vtx,
  int          *edge_to_iso_vtx,
  int          *out_iso_n_edge,
  int         **out_iso_edge_vtx,
  PDM_g_num_t **out_iso_edge_parent_gnum,
  int          *out_iso_n_edge_bnd_tag,
  int         **out_iso_edge_bnd_tag_idx,
  int         **out_iso_edge_bnd_tag,
  int          *out_iso_n_edge_parent,
  int         **out_iso_edge_parent_idx,
  int         **out_iso_edge_parent
)
{
  int debug      = 0;
  int debug_loop = 0;


  /* First loop to count */
  int *iso_edge_def = PDM_array_zeros_int(n_edge);
  int  iso_n_edge         = 0;
  int  iso_n_edge_parent  = 0;
  int  iso_n_edge_bnd_tag = 0;

  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    if (debug_loop==1) log_trace("\ti_elt = %d\n", i_elt);

    int is_bnd = elt_bnd_tag[i_elt]!=0;

    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 3; i_edge++) {
      char bit = 0;
      if (elt_edge[3*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[3*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges
    if (debug_loop==1) log_trace("\t\tpattern = %d\n", pattern);

    if (pattern > 0 && pattern < 7) {
      iso_n_edge++;
      iso_n_edge_parent++;
      if (is_bnd) {
        iso_n_edge_bnd_tag++;
      }
    }
    else if (pattern == 0) {

      int offset = _pattern_permutation_2d[pattern];
      int perm_elt_vtx [3];
      int perm_elt_edge[3];
      for (int i = 0; i < 3; i++) {
        perm_elt_vtx [i] = elt_vtx [3*i_elt + (offset+i)%3] - 1;
        perm_elt_edge[i] = elt_edge[3*i_elt + (offset+i)%3];
      }


      // > For iso edge on entry mesh, count only once edge
      int vtx_on_vtx0 = vtx_to_iso_vtx[perm_elt_vtx[0]];
      int vtx_on_vtx1 = vtx_to_iso_vtx[perm_elt_vtx[1]];
      int vtx_on_vtx2 = vtx_to_iso_vtx[perm_elt_vtx[2]];
      
      if (vtx_on_vtx0!=0 && vtx_on_vtx1!=0) {
        int edge_id2 = PDM_ABS(perm_elt_edge[2]) - 1;
        if (debug_loop==1) log_trace("\t\tedge_id2 = %d\n", edge_id2);
        if (iso_edge_def[edge_id2]==0) { // edge not build yet
          iso_edge_def[edge_id2]=1;
          iso_n_edge++;
        }
        iso_n_edge_parent += edge_parent_idx[edge_id2+1]-edge_parent_idx[edge_id2];
        if (is_bnd) {
          iso_n_edge_bnd_tag += edge_bnd_tag_idx[edge_id2+1]-edge_bnd_tag_idx[edge_id2];
        }
      }

      else if (vtx_on_vtx0!=0 && vtx_on_vtx2!=0) {
        int edge_id1 = PDM_ABS(perm_elt_edge[1]) - 1;
        if (debug_loop==1) log_trace("\t\tedge_id1 = %d\n", edge_id1);
        if (iso_edge_def[edge_id1]==0) { // edge not build yet
          iso_edge_def[edge_id1]=1;
          iso_n_edge++;
        }
        iso_n_edge_parent += edge_parent_idx[edge_id1+1]-edge_parent_idx[edge_id1];
        if (is_bnd) {
          iso_n_edge_bnd_tag += edge_bnd_tag_idx[edge_id1+1]-edge_bnd_tag_idx[edge_id1];
        }
      }

      else if (vtx_on_vtx1!=0 && vtx_on_vtx2!=0) {
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        if (debug_loop==1) log_trace("\t\tedge_id0 = %d\n", edge_id0);
        if (iso_edge_def[edge_id0]==0) { // edge not build yet
          iso_edge_def[edge_id0]=1;
          iso_n_edge++;
        }
        iso_n_edge_parent += edge_parent_idx[edge_id0+1]-edge_parent_idx[edge_id0];
        if (is_bnd) {
          iso_n_edge_bnd_tag += edge_bnd_tag_idx[edge_id0+1]-edge_bnd_tag_idx[edge_id0];
        }
      }
    }
  } // End of loop on elements
  if (debug==1) log_trace("iso_n_edge         = %d\n", iso_n_edge);
  if (debug==1) log_trace("iso_n_edge_bnd_tag = %d\n", iso_n_edge_bnd_tag);
  if (debug==1) log_trace("iso_n_edge_parent  = %d\n", iso_n_edge_parent);


  /* Allocate */

  iso_n_edge         += *out_iso_n_edge;
  iso_n_edge_parent  += *out_iso_n_edge_parent;
  iso_n_edge_bnd_tag += *out_iso_n_edge_bnd_tag;
  int         *iso_edge_vtx          = NULL;
  PDM_g_num_t *iso_edge_parent_gnum  = NULL;
  int         *iso_edge_bnd_tag_idx  = NULL;
  int         *iso_edge_bnd_tag      = NULL;
  int         *iso_edge_parent_idx   = NULL;
  int         *iso_edge_parent       = NULL;
  PDM_realloc(*out_iso_edge_vtx        , iso_edge_vtx        , 2 * iso_n_edge        , int        );
  PDM_realloc(*out_iso_edge_parent_gnum, iso_edge_parent_gnum, 2 * iso_n_edge        , PDM_g_num_t);
  PDM_realloc(*out_iso_edge_bnd_tag_idx, iso_edge_bnd_tag_idx,     iso_n_edge+1      , int        );
  PDM_realloc(*out_iso_edge_bnd_tag    , iso_edge_bnd_tag    ,     iso_n_edge_bnd_tag, int        );
  PDM_realloc(*out_iso_edge_parent_idx , iso_edge_parent_idx ,     iso_n_edge+1      , int        );
  PDM_realloc(*out_iso_edge_parent     , iso_edge_parent     ,     iso_n_edge_parent , int        );
  PDM_array_reset_int(iso_edge_def, n_edge, 0);
  iso_edge_parent_idx [0] = 0;
  iso_edge_bnd_tag_idx[0] = 0;

  iso_n_edge         = *out_iso_n_edge;
  iso_n_edge_bnd_tag = *out_iso_n_edge_bnd_tag;
  iso_n_edge_parent  = *out_iso_n_edge_parent;

  if (debug==1) log_trace("iso_n_edge         = %d\n", iso_n_edge       );
  if (debug==1) log_trace("iso_n_edge_bnd_tag = %d\n", iso_n_edge_bnd_tag);
  if (debug==1) log_trace("iso_n_edge_parent  = %d\n", iso_n_edge_parent);


  /**
   * Because all isovalues are on same mesh,
   * some entities will have multiple descendent
   * so we need to count them to generate unique gnum...
   */
  int *elt_n_child = NULL;
  PDM_calloc(elt_n_child, n_elt, int);


  /* Second loop to fill */
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {

    int elt_tag = elt_bnd_tag[i_elt];
    if (debug_loop==1) log_trace("\ti_elt = %d\n", i_elt);
    if (debug_loop==1) PDM_log_trace_array_int (iso_edge_def, n_edge , "iso_edge_def ::");
    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 3; i_edge++) {
      char bit = 0;
      if (elt_edge[3*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[3*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges

    if (debug_loop==1) log_trace("\t\tpattern = %d\n", pattern);
    if (pattern == 7) { // Three edge are active -> impossible
      continue;
    }

    int offset = _pattern_permutation_2d[pattern];
    int perm_elt_vtx [3];
    int perm_elt_edge[3];
    for (int i = 0; i < 3; i++) {
      perm_elt_vtx [i] = elt_vtx [3*i_elt + (offset+i)%3] - 1;
      perm_elt_edge[i] = elt_edge[3*i_elt + (offset+i)%3];
    }



    switch (pattern) {
      case 0:{
        int vtx_on_vtx0 = vtx_to_iso_vtx[perm_elt_vtx[0]];
        int vtx_on_vtx1 = vtx_to_iso_vtx[perm_elt_vtx[1]];
        int vtx_on_vtx2 = vtx_to_iso_vtx[perm_elt_vtx[2]];
        if (vtx_on_vtx0!=0 && vtx_on_vtx1!=0) {

          int edge_id2 = PDM_ABS(perm_elt_edge[2]) - 1;
          if (iso_edge_def[edge_id2]==0) { // edge not build yet
            iso_edge_def[edge_id2]=1;
          
            iso_edge_vtx[2*iso_n_edge  ] = vtx_on_vtx0;
            iso_edge_vtx[2*iso_n_edge+1] = vtx_on_vtx1;
            if (debug_loop==1) log_trace("i_elt = %d (pattern %d 1) : ivtx0 = %d ; ivtx1 = %d\n", i_elt, pattern, vtx_on_vtx0, vtx_on_vtx2);

            if (vtx_gnum[perm_elt_vtx[0]] > vtx_gnum[perm_elt_vtx[1]]) {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[1]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[0]];
            }
            else {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[0]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[1]];
            }
        
            // > Fill parent
            int i_beg_parent = edge_parent_idx[edge_id2  ];
            int i_end_parent = edge_parent_idx[edge_id2+1];
            int n_edge_parent = i_end_parent - i_beg_parent;
            for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
              iso_edge_parent[iso_n_edge_parent++] = edge_parent[i_parent];
            }
            iso_edge_parent_idx[iso_n_edge+1] = iso_edge_parent_idx[iso_n_edge]+n_edge_parent;

            // > Fill bnd tag
            if (elt_tag!=0) {
              int i_beg_bnd_tag = edge_bnd_tag_idx[edge_id2  ];
              int i_end_bnd_tag = edge_bnd_tag_idx[edge_id2+1];
              int n_edge_bnd_tag = i_end_bnd_tag - i_beg_bnd_tag;
              for (int i_bnd_tag=i_beg_bnd_tag; i_bnd_tag<i_end_bnd_tag; ++i_bnd_tag) {
                iso_edge_bnd_tag[iso_n_edge_bnd_tag++] = edge_bnd_tag[i_bnd_tag];
              }
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge]+n_edge_bnd_tag;
            }
            else {
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge];
            }

            iso_n_edge++;
          }
        }
        else if (vtx_on_vtx0!=0 && vtx_on_vtx2!=0) {

          int edge_id1 = PDM_ABS(perm_elt_edge[1]) - 1;
          if (iso_edge_def[edge_id1]==0) { // edge not build yet
            iso_edge_def[edge_id1]=1;
            
            iso_edge_vtx[2*iso_n_edge  ] = vtx_on_vtx0;
            iso_edge_vtx[2*iso_n_edge+1] = vtx_on_vtx2;

            if (debug_loop==1) log_trace("i_elt = %d (pattern %d 2) : ivtx0 = %d ; ivtx1 = %d\n", i_elt, pattern, vtx_on_vtx0, vtx_on_vtx2);
            if (vtx_gnum[perm_elt_vtx[0]] > vtx_gnum[perm_elt_vtx[2]]) {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[2]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[0]];
            }
            else {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[0]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[2]];
            }
        
            // > Fill parent
            int i_beg_parent = edge_parent_idx[edge_id1  ];
            int i_end_parent = edge_parent_idx[edge_id1+1];
            int n_edge_parent = i_end_parent - i_beg_parent;
            for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
              iso_edge_parent[iso_n_edge_parent++] = edge_parent[i_parent];
            }
            iso_edge_parent_idx[iso_n_edge+1] = iso_edge_parent_idx[iso_n_edge]+n_edge_parent;

            // > Fill bnd tag
            if (elt_tag!=0) {
              int i_beg_bnd_tag = edge_bnd_tag_idx[edge_id1  ];
              int i_end_bnd_tag = edge_bnd_tag_idx[edge_id1+1];
              int n_edge_bnd_tag = i_end_bnd_tag - i_beg_bnd_tag;
              for (int i_bnd_tag=i_beg_bnd_tag; i_bnd_tag<i_end_bnd_tag; ++i_bnd_tag) {
                iso_edge_bnd_tag[iso_n_edge_bnd_tag++] = edge_bnd_tag[i_bnd_tag];
              }
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge]+n_edge_bnd_tag;
            }
            else {
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge];
            }

            iso_n_edge++;
          }
        }
        else if (vtx_on_vtx1!=0 && vtx_on_vtx2!=0) {

          int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
          if (iso_edge_def[edge_id0]==0) { // edge not build yet
            iso_edge_def[edge_id0]=1;
            
            iso_edge_vtx[2*iso_n_edge  ] = vtx_on_vtx1;
            iso_edge_vtx[2*iso_n_edge+1] = vtx_on_vtx2;
            if (debug_loop==1) log_trace("i_elt = %d (pattern %d 3) : ivtx0 = %d ; ivtx1 = %d\n", i_elt, pattern, vtx_on_vtx0, vtx_on_vtx2);

            if (vtx_gnum[perm_elt_vtx[1]] > vtx_gnum[perm_elt_vtx[2]]) {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[2]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[1]];
            }
            else {
              iso_edge_parent_gnum[2*iso_n_edge  ] = vtx_gnum[perm_elt_vtx[1]];
              iso_edge_parent_gnum[2*iso_n_edge+1] = vtx_gnum[perm_elt_vtx[2]];
            }
        
            // > Fill parent
            int i_beg_parent = edge_parent_idx[edge_id0  ];
            int i_end_parent = edge_parent_idx[edge_id0+1];
            int n_edge_parent = i_end_parent - i_beg_parent;
            for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
              iso_edge_parent[iso_n_edge_parent++] = edge_parent[i_parent];
            }
            iso_edge_parent_idx[iso_n_edge+1] = iso_edge_parent_idx[iso_n_edge]+n_edge_parent;

            // > Fill bnd tag
            if (elt_tag!=0) {
              int i_beg_bnd_tag = edge_bnd_tag_idx[edge_id0  ];
              int i_end_bnd_tag = edge_bnd_tag_idx[edge_id0+1];
              int n_edge_bnd_tag = i_end_bnd_tag - i_beg_bnd_tag;
              for (int i_bnd_tag=i_beg_bnd_tag; i_bnd_tag<i_end_bnd_tag; ++i_bnd_tag) {
                iso_edge_bnd_tag[iso_n_edge_bnd_tag++] = edge_bnd_tag[i_bnd_tag];
              }
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge]+n_edge_bnd_tag;
            }
            else {
              iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge];
            }

            iso_n_edge++;

          }
        }
        else {
          // 3 points on element -> quite ambiguous case  
          if (vtx_on_vtx0!=0 && vtx_on_vtx1!=0 && vtx_on_vtx2!=0) {
            printf("WARNING: element "PDM_FMT_G_NUM" is fully on isosurface.\n", elt_gnum[i_elt]);
            // TODO: keep warning ? if on BC of volumic mesh, can be anoying
          }
          // 1 point on isosurface: may be managed/generated by neighbourhood
          // 0 point on isosurface: nothing to do
          continue;
        }

        break;
      }
      case 1: case 2: case 4: {
        int vtx_on_vtx0 = vtx_to_iso_vtx[perm_elt_vtx[0]];
        assert(vtx_on_vtx0);
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        int vtx_on_edge0 = edge_to_iso_vtx[edge_id0];
        assert(vtx_on_edge0);

        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[2]]) {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_vtx0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_edge0;

          iso_edge_parent_gnum[2*iso_n_edge  ] =  elt_gnum   [i_elt];
          iso_edge_parent_gnum[2*iso_n_edge+1] =--elt_n_child[i_elt];
        }
        else {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_edge0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_vtx0;

          iso_edge_parent_gnum[2*iso_n_edge  ] =  elt_gnum   [i_elt];
          iso_edge_parent_gnum[2*iso_n_edge+1] =--elt_n_child[i_elt];
        }

        iso_edge_parent_idx[iso_n_edge+1] = iso_edge_parent_idx[iso_n_edge]+1;
        iso_edge_parent[iso_n_edge_parent++] = i_elt+1;

        if (elt_tag!=0) {
          iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge]+1;
          iso_edge_bnd_tag[iso_n_edge_bnd_tag++] = elt_tag;
        }
        else {
          iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge];
        }

        iso_n_edge++;

        break;
      }

      case 3: case 5: case 6: {
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        int vtx_on_edge0 = edge_to_iso_vtx[edge_id0];
        assert(vtx_on_edge0);

        int edge_id1 = PDM_ABS(perm_elt_edge[1]) - 1;
        int vtx_on_edge1 = edge_to_iso_vtx[edge_id1];
        assert(vtx_on_edge1);

        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[2]]) {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_edge1;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_edge0;

          iso_edge_parent_gnum[2*iso_n_edge  ] =  elt_gnum   [i_elt];
          iso_edge_parent_gnum[2*iso_n_edge+1] =--elt_n_child[i_elt];
        }
        else {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_edge0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_edge1;

          iso_edge_parent_gnum[2*iso_n_edge  ] =  elt_gnum   [i_elt];
          iso_edge_parent_gnum[2*iso_n_edge+1] =--elt_n_child[i_elt];
        }

        iso_edge_parent_idx[iso_n_edge+1] = iso_edge_parent_idx[iso_n_edge]+1;
        iso_edge_parent[iso_n_edge_parent++] = i_elt+1;

        if (elt_tag!=0) {
          iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge]+1;
          iso_edge_bnd_tag[iso_n_edge_bnd_tag++] = elt_tag;
        }
        else {
          iso_edge_bnd_tag_idx[iso_n_edge+1] = iso_edge_bnd_tag_idx[iso_n_edge];
        }

        iso_n_edge++;
        
        break;
      }

      default: {
        break;
      }
    }

  } // End of loop on elements


  *out_iso_n_edge           = iso_n_edge;
  *out_iso_edge_vtx         = iso_edge_vtx;
  *out_iso_edge_parent_gnum = iso_edge_parent_gnum;
  *out_iso_n_edge_bnd_tag   = iso_n_edge_bnd_tag;
  *out_iso_edge_bnd_tag_idx = iso_edge_bnd_tag_idx;
  *out_iso_edge_bnd_tag     = iso_edge_bnd_tag;
  *out_iso_n_edge_parent    = iso_n_edge_parent;
  *out_iso_edge_parent_idx  = iso_edge_parent_idx;
  *out_iso_edge_parent      = iso_edge_parent;

  if (debug==1) {
    log_trace("\n");
    log_trace("iso_n_edge = %d\n", iso_n_edge);
    PDM_log_trace_array_int(iso_edge_vtx, 2*iso_n_edge, "iso_edge_vtx ::");
    log_trace("\n");
    int edge_bnd_tag_size = iso_edge_bnd_tag_idx[iso_n_edge];
    PDM_log_trace_array_int(*out_iso_edge_bnd_tag_idx, iso_n_edge+1     , "out_iso_edge_bnd_tag_idx ::");
    PDM_log_trace_array_int(*out_iso_edge_bnd_tag    , edge_bnd_tag_size, "out_iso_edge_bnd_tag     ::");
    log_trace("\n");
    int edge_parent_size = iso_edge_parent_idx[iso_n_edge];
    PDM_log_trace_array_int(*out_iso_edge_parent_idx, iso_n_edge+1    , "out_iso_edge_parent_idx ::");
    PDM_log_trace_array_int(*out_iso_edge_parent    , edge_parent_size, "out_iso_edge_parent     ::");
    log_trace("\n");
    PDM_log_trace_array_long(iso_edge_parent_gnum, 2*iso_n_edge, "iso_edge_parent_gnum ::");
  }

  PDM_free(iso_edge_def);
  PDM_free(elt_n_child);
}




static void
_contouring_tetrahedra
(
  int           n_elt,
  int          *elt_vtx,
  PDM_g_num_t  *elt_gnum,
  int           n_face,
  int          *face_vtx_idx,
  int          *face_vtx,
  int          *face_parent_idx,
  int          *face_parent,
  int          *elt_edge,
  int          *elt_face,
  PDM_g_num_t  *vtx_gnum,
  double       *vtx_field,
  int          *vtx_to_iso_vtx,
  int          *edge_to_iso_vtx,
  int          *out_iso_n_face,
  int         **out_iso_face_vtx_idx,
  int         **out_iso_face_vtx,
  PDM_g_num_t **out_iso_face_parent_gnum,
  int          *out_iso_n_face_parent,
  int         **out_iso_face_parent_idx,
  int         **out_iso_face_parent
)
{
  int debug      = 1;
  int debug_loop = 0;
  // PDM_UNUSED(elt_gnum);
  // PDM_UNUSED(face_vtx_idx);
  // PDM_UNUSED(face_vtx);
  // PDM_UNUSED(face_parent_idx);
  // PDM_UNUSED(face_parent);
  // PDM_UNUSED(vtx_gnum);
  // PDM_UNUSED(out_iso_face_parent_gnum);
  // PDM_UNUSED(out_iso_n_face_parent);
  // PDM_UNUSED(out_iso_face_parent_idx);
  // PDM_UNUSED(out_iso_face_parent);

  /* First loop to count */
  int s_face_vtx = 0;
  int *iso_face_def     = PDM_array_zeros_int(n_face);
  int iso_n_face        = *out_iso_n_face;
  int iso_n_face_parent = 0;
  
  if (debug) log_trace("iso_n_face = %d\n", iso_n_face);
  
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    if (debug_loop==1) log_trace("i_elt = %d\n", i_elt);

    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 6; i_edge++) {
      char bit = 0;
      if (elt_edge[6*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[6*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges


    switch (pattern) {
      case 0: {
        // Configurations 0: not edge cross by isosurface, but on face can be on isosurface
        
        // > Check that vertices on face are on current isosurface
        int n_active_face_on_elt = 0;
        int n_active_vtx = 0;
        for (int i_vtx=0; i_vtx<4; i_vtx++) {
          int vtx_id = elt_vtx[4*i_elt + i_vtx] - 1;
          if (vtx_to_iso_vtx[vtx_id]!=0) {
            n_active_vtx++;
          }
        }
        if (debug_loop==1) log_trace("\tn_active_vtx = %d\n", n_active_vtx);
        
        if (n_active_vtx>=3) {
          for (int i_face = 0; i_face < 4; i_face++) {
            int face_id = PDM_ABS(elt_face[4*i_elt+i_face]);
            if (face_id!=0 && iso_face_def[face_id-1]==0){            
              if (debug_loop==1) log_trace("\tface_id = %d\n", face_id);
              iso_n_face += 1;
              s_face_vtx += 3;
              n_active_face_on_elt++;
              iso_face_def[face_id-1] = 1;
              iso_n_face_parent += face_parent_idx[face_id]-face_parent_idx[face_id-1];
            }
          }
        }
        if (n_active_face_on_elt>1) {
          printf("WARNING: element "PDM_FMT_G_NUM" tetra has more than one face on isosurface.\n", elt_gnum[i_elt]);
        }

        break;
      }

      case 1: case 2: case 4: case 8: case 16: case 32:
      case 3: case 5: case 6: case 9: case 10: case 17: case 20: case 24: case 34: case 36: case 40: case 48:
      case 7: case 25: case 42: case 52: {
        // Configurations 1, 2a and 3b : one triangle
        iso_n_face += 1;
        s_face_vtx += 3;
        iso_n_face_parent++;
        break;
      }
      case 30: case 45: case 51: {
        // Configuration 4b : one quadrangle (or two triangles)
        iso_n_face += 1;
        s_face_vtx += 4;
        iso_n_face_parent++;
        break;
      }
      default: {
        break;
      }
    }

  } // End of loop on elements
  if (debug==1) log_trace("iso_n_face        = %d\n", iso_n_face);
  if (debug==1) log_trace("iso_n_face_parent = %d\n", iso_n_face_parent);




  /* Allocate */
  iso_n_face_parent += *out_iso_n_face_parent;
  int prev_size = 0;
  if (*out_iso_face_vtx_idx != NULL) {
    prev_size = (*out_iso_face_vtx_idx)[*out_iso_n_face];
  }
  if (debug==1) log_trace("\n");
  if (debug==1) log_trace("iso_n_face        = %d\n", iso_n_face);
  if (debug==1) log_trace("iso_n_face_parent = %d\n", iso_n_face_parent);
  if (debug==1) log_trace("s_face_vtx        = %d\n", s_face_vtx);

  PDM_g_num_t *iso_face_parent_gnum = NULL;
  int         *iso_face_vtx_idx     = NULL;
  int         *iso_face_vtx         = NULL;
  int         *iso_face_parent_idx  = NULL;
  int         *iso_face_parent      = NULL;
  PDM_realloc(*out_iso_face_parent_gnum, iso_face_parent_gnum, 3*iso_n_face,              PDM_g_num_t);
  PDM_realloc(*out_iso_face_vtx_idx,     iso_face_vtx_idx,       iso_n_face + 1,          int        );
  PDM_realloc(*out_iso_face_vtx,         iso_face_vtx,           prev_size  + s_face_vtx, int        );
  PDM_realloc(*out_iso_face_parent_idx,  iso_face_parent_idx,    iso_n_face + 1,          int        );
  PDM_realloc(*out_iso_face_parent    ,  iso_face_parent,        iso_n_face_parent,       int        );
  iso_face_vtx_idx   [0] = 0;
  PDM_array_reset_int(iso_face_def, n_face, 0);
  iso_face_parent_idx[0] = 0.;
  
  iso_n_face        = *out_iso_n_face;
  iso_n_face_parent = *out_iso_n_face_parent;


  /**
   * Because all isovalues are on same mesh,
   * some entities will have multiple descendent
   * so we need to count them to generate unique gnum...
   */
  int *elt_n_child = NULL;
  PDM_calloc(elt_n_child, n_elt, int);


  /* Second loop to fill */
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {
    if (debug_loop==1) log_trace("\n");
    if (debug_loop==1) log_trace("i_elt = %d\n", i_elt);

    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 6; i_edge++) {
      char bit = 0;
      if (elt_edge[6*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[6*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges


    int perm_elt_vtx [4];
    int perm_elt_edge[6];
    for (int i = 0; i < 4; i++) {
      int j = _pattern_cell_vtx_permutation[pattern][i];
      perm_elt_vtx[i] = elt_vtx[4*i_elt + j] - 1;
    }
    for (int i = 0; i < 6; i++) {
      int j = _pattern_cell_edge_permutation[pattern][i];
      perm_elt_edge[i] = PDM_SIGN(j) * elt_edge[6*i_elt + PDM_ABS(j)-1];
    }

    switch (pattern) {
      case 0: {
        // > Check that vertices on face are on current isosurface
        double face_val = 0.;
        double vtx_val  = 0.;
        int n_active_vtx = 0;
        for (int i_vtx=0; i_vtx<4; i_vtx++) {
          int vtx_id = elt_vtx[4*i_elt + i_vtx] - 1;
          if (vtx_to_iso_vtx[vtx_id]!=0) {
            n_active_vtx++;
            face_val = vtx_field[vtx_id];
          }
          else {
            vtx_val = vtx_field[vtx_id];
          }
        }
        if (n_active_vtx>=3) {
          for (int i_face = 0; i_face < 4; i_face++) {
            int face_id = PDM_ABS(elt_face[4*i_elt+i_face]);
            if (face_id!=0 && iso_face_def[face_id-1]==0){            

              if (debug_loop==1) log_trace("\tface_id = %d\n", face_id);

              iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face] + 3;
              int     idx =     face_vtx_idx[face_id-1];
              int iso_idx = iso_face_vtx_idx[iso_n_face];
              if (debug_loop==1) log_trace("\tidx = %d\n", idx);

              if (face_val<vtx_val) {
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx+2]-1];
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx+1]-1];
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx  ]-1];
              }
              else if (vtx_val<=face_val) {
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx  ]-1];
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx+1]-1];
                iso_face_vtx[iso_idx++] = vtx_to_iso_vtx[face_vtx[idx+2]-1];
              }

              iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face]+3;
          
              int i_beg_parent  = face_parent_idx[face_id-1];
              int i_end_parent  = face_parent_idx[face_id  ];
              int n_face_parent = i_end_parent - i_beg_parent;
              for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
                iso_face_parent[iso_n_face_parent++] = face_parent[i_parent];
              }
              iso_face_parent_idx[iso_n_face+1] = iso_face_parent_idx[iso_n_face]+n_face_parent;

              iso_face_parent_gnum[3*iso_n_face  ] = vtx_gnum[face_vtx[idx  ]-1];
              iso_face_parent_gnum[3*iso_n_face+1] = vtx_gnum[face_vtx[idx+1]-1];
              iso_face_parent_gnum[3*iso_n_face+2] = vtx_gnum[face_vtx[idx+2]-1];

              iso_n_face++;

              iso_face_def[face_id-1] = 1;
            }
          }
        }

        break;
      }

      case 1: case 2: case 4: case 8: case 16: case 32: {
        // Configuration 1
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        int vtx_on_edge0 = edge_to_iso_vtx[edge_id0];
        assert(vtx_on_edge0);

        int vtx_on_vtx2 = vtx_to_iso_vtx[perm_elt_vtx[2]];
        assert(vtx_on_vtx2);
        int vtx_on_vtx3 = vtx_to_iso_vtx[perm_elt_vtx[3]];
        assert(vtx_on_vtx3);

        iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face] + 3;
        int idx = iso_face_vtx_idx[iso_n_face];
        iso_face_vtx[idx++] = vtx_on_edge0;
        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[0]]) {
          iso_face_vtx[idx++] = vtx_on_vtx2;
          iso_face_vtx[idx++] = vtx_on_vtx3;
        }
        else {
          iso_face_vtx[idx++] = vtx_on_vtx3;
          iso_face_vtx[idx++] = vtx_on_vtx2;
        }

        iso_face_parent_idx[iso_n_face+1] = iso_face_parent_idx[iso_n_face] + 1;
        iso_face_parent[iso_n_face_parent++] = i_elt + 1;

        iso_face_parent_gnum[3*iso_n_face  ] =  elt_gnum   [i_elt];
        iso_face_parent_gnum[3*iso_n_face+1] =--elt_n_child[i_elt];
        iso_face_parent_gnum[3*iso_n_face+2] =  0;

        iso_n_face++;
        
        break;
      }

      case 3: case 5: case 6: case 9: case 10: case 17: case 20: case 24: case 34: case 36: case 40: case 48: {
        // Configuration 2a
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        int vtx_on_edge0 = edge_to_iso_vtx[edge_id0];
        assert(vtx_on_edge0);
        int edge_id1 = PDM_ABS(perm_elt_edge[1]) - 1;
        int vtx_on_edge1 = edge_to_iso_vtx[edge_id1];
        assert(vtx_on_edge1);

        int vtx_on_vtx3 = vtx_to_iso_vtx[perm_elt_vtx[3]];
        assert(vtx_on_vtx3);

        iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face] + 3;
        int idx = iso_face_vtx_idx[iso_n_face];
        iso_face_vtx[idx++] = vtx_on_vtx3;
        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[0]]) {
          iso_face_vtx[idx++] = vtx_on_edge0;
          iso_face_vtx[idx++] = vtx_on_edge1;
        }
        else {
          iso_face_vtx[idx++] = vtx_on_edge1;
          iso_face_vtx[idx++] = vtx_on_edge0;
        }

        iso_face_parent_idx[iso_n_face+1] = iso_face_parent_idx[iso_n_face] + 1;
        iso_face_parent[iso_n_face_parent++] = i_elt + 1;

        iso_face_parent_gnum[3*iso_n_face  ] =  elt_gnum   [i_elt];
        iso_face_parent_gnum[3*iso_n_face+1] =--elt_n_child[i_elt];
        iso_face_parent_gnum[3*iso_n_face+2] =  0;

        iso_n_face++;

        break;
      }

      case 7: case 25: case 42: case 52: {
        // Configuration 3b
        int vtx_on_edge[3];
        for (int i_edge = 0; i_edge < 3; i_edge++) {
          int edge_id = PDM_ABS(perm_elt_edge[i_edge]) - 1;
          vtx_on_edge[i_edge] = edge_to_iso_vtx[edge_id];
          assert(vtx_on_edge[i_edge]);
        }

        iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face] + 3;
        int idx = iso_face_vtx_idx[iso_n_face];
        iso_face_vtx[idx++] = vtx_on_edge[0];
        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[0]]) {
          iso_face_vtx[idx++] = vtx_on_edge[1];
          iso_face_vtx[idx++] = vtx_on_edge[2];
        }
        else {
          iso_face_vtx[idx++] = vtx_on_edge[2];
          iso_face_vtx[idx++] = vtx_on_edge[1];
        }

        iso_face_parent_idx[iso_n_face+1] = iso_face_parent_idx[iso_n_face] + 1;
        iso_face_parent[iso_n_face_parent++] = i_elt + 1;

        iso_face_parent_gnum[3*iso_n_face  ] =  elt_gnum   [i_elt];
        iso_face_parent_gnum[3*iso_n_face+1] =--elt_n_child[i_elt];
        iso_face_parent_gnum[3*iso_n_face+2] =  0;

        iso_n_face++;

        break;
      }

      case 30: case 45: case 51: {
        // Configuration 4b
        int vtx_on_edge[4];
        for (int i_edge = 0; i_edge < 4; i_edge++) {
          int edge_id = PDM_ABS(perm_elt_edge[i_edge+1]) - 1;
          vtx_on_edge[i_edge] = edge_to_iso_vtx[edge_id];
          assert(vtx_on_edge[i_edge]);
        }

        iso_face_vtx_idx[iso_n_face+1] = iso_face_vtx_idx[iso_n_face] + 4;
        int idx = iso_face_vtx_idx[iso_n_face];
        iso_face_vtx[idx++] = vtx_on_edge[0];
        if (vtx_field[perm_elt_vtx[3]] > vtx_field[perm_elt_vtx[0]]) {
          iso_face_vtx[idx++] = vtx_on_edge[1];
          iso_face_vtx[idx++] = vtx_on_edge[3];
          iso_face_vtx[idx++] = vtx_on_edge[2];
        }
        else {
          iso_face_vtx[idx++] = vtx_on_edge[2];
          iso_face_vtx[idx++] = vtx_on_edge[3];
          iso_face_vtx[idx++] = vtx_on_edge[1];
        }

        iso_face_parent_idx[iso_n_face+1] = iso_face_parent_idx[iso_n_face] + 1;
        iso_face_parent[iso_n_face_parent++] = i_elt + 1;

        iso_face_parent_gnum[3*iso_n_face  ] =  elt_gnum   [i_elt];
        iso_face_parent_gnum[3*iso_n_face+1] =--elt_n_child[i_elt];
        iso_face_parent_gnum[3*iso_n_face+2] =  0;

        iso_n_face++;
        
        break;
      }

      default: {
        continue;
        break;
      }
    }


  } // End of loop on elements

  if (debug==1) {
    log_trace("\n");
    log_trace("iso_n_face   = %d\n", iso_n_face);
    int face_vtx_size = iso_face_vtx_idx[iso_n_face];
    PDM_log_trace_array_int(iso_face_vtx_idx, iso_n_face+1 , "iso_face_vtx_idx ::");
    PDM_log_trace_array_int(iso_face_vtx    , face_vtx_size, "iso_face_vtx     ::");
    log_trace("\n");
    int face_parent_size = iso_face_parent_idx[iso_n_face];
    PDM_log_trace_array_int(iso_face_parent_idx, iso_n_face+1    , "iso_face_parent_idx ::");
    PDM_log_trace_array_int(iso_face_parent    , face_parent_size, "iso_face_parent     ::");
    log_trace("\n");
    PDM_log_trace_array_long(iso_face_parent_gnum, 3*iso_n_face, "iso_face_parent_gnum ::");
  }

  *out_iso_n_face           = iso_n_face;
  *out_iso_face_vtx_idx     = iso_face_vtx_idx;
  *out_iso_face_vtx         = iso_face_vtx;
  *out_iso_n_face_parent    = iso_n_face_parent;
  *out_iso_face_parent_idx  = iso_face_parent_idx;
  *out_iso_face_parent      = iso_face_parent;
  *out_iso_face_parent_gnum = iso_face_parent_gnum;
 
  PDM_free(iso_face_def);
  PDM_free(elt_n_child);
}


static void
_debug_ngon
(
 int     n_face,
 int     n_edge,
 int    *face_edge_idx,
 int    *face_edge,
 int    *edge_vtx,
 double *vtx_coord,
 double *vtx_field,
 double  isovalue
 )
{
  int *face_vtx = NULL;
  PDM_compute_face_vtx_from_face_and_edge(n_face,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          &face_vtx);
  int n_vtx = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_vtx = PDM_MAX(n_vtx, face_vtx[i]);
  }

  double *_vtx_field = NULL; 
  PDM_malloc(_vtx_field, n_vtx, double);
  for (int i = 0; i < n_vtx; i++) {
    _vtx_field[i] = vtx_field[i] - isovalue;
  }

  const char *name1 = "debug_ngon_faces.vtk";
  PDM_vtk_write_polydata_field(name1,
                               n_vtx,
                               vtx_coord,
                               NULL,
                               n_face,
                               face_edge_idx,
                               face_vtx,
                               NULL,
                               NULL,
                               NULL,
                               "field",
                               _vtx_field);
  PDM_free(face_vtx);

  const char   *name2 = "debug_ngon_edges.vtk";
  const char   *field_name [] = {"field"};
  const double *field_value[] = {_vtx_field};
  PDM_vtk_write_std_elements_ho_with_vtx_field(name2,
                                               1,
                                               n_vtx,
                                               vtx_coord,
                                               NULL,
                                               PDM_MESH_NODAL_BAR2,
                                               n_edge,
                                               edge_vtx,
                                               NULL,
                                               0,
                                               NULL,
                                               NULL,
                                               1,
                                               field_name,
                                               field_value);
  PDM_free(_vtx_field);
}


static void
_debug_ngon_cell
(
 int     cell_face_n,
 int    *cell_face,
 int    *face_edge_idx,
 int    *face_edge,
 int    *edge_vtx,
 double *vtx_coord,
 double *vtx_field,
 double  isovalue
 )
{
  int *_face_vtx_idx = NULL;
  PDM_malloc(_face_vtx_idx, cell_face_n+1, int);
  _face_vtx_idx[0] = 0;
  for (int i = 0; i < cell_face_n; i++) {
    int i_face = PDM_ABS(cell_face[i]) - 1;
    _face_vtx_idx[i+1] = _face_vtx_idx[i] + face_edge_idx[i_face+1] - face_edge_idx[i_face];
  }

  int *_face_edge = NULL;
  PDM_malloc(_face_edge, _face_vtx_idx[cell_face_n], int);
  for (int i = 0; i < cell_face_n; i++) {
    int i_face = PDM_ABS(cell_face[i]) - 1;
    int k = _face_vtx_idx[i];
    for (int j = face_edge_idx[i_face]; j < face_edge_idx[i_face+1]; j++) {
      _face_edge[k++] = face_edge[j] * PDM_SIGN(cell_face[i]);
    }
  }
  // PDM_log_trace_connectivity_int(_face_vtx_idx, _face_edge, cell_face_n, "_face_edge : ");

  int *_face_vtx = NULL;
  PDM_compute_face_vtx_from_face_and_edge(cell_face_n,
                                          _face_vtx_idx,
                                          _face_edge,
                                          edge_vtx,
                                          &_face_vtx);
  PDM_free(_face_edge);

  // PDM_log_trace_connectivity_int(_face_vtx_idx, _face_vtx, cell_face_n, "_face_vtx : ");
  int *order_in_unique = NULL;
  PDM_malloc(order_in_unique, _face_vtx_idx[cell_face_n], int);
  int _n_vtx = PDM_inplace_unique_int_with_order_in_unique(_face_vtx,
                                                           order_in_unique,
                                                           0,
                                                           _face_vtx_idx[cell_face_n]-1);

  double *_vtx_coord = NULL;
  double *_vtx_field = NULL;
  PDM_malloc(_vtx_coord, _n_vtx * 3, double);
  PDM_malloc(_vtx_field, _n_vtx    , double);
  for (int i = 0; i < _n_vtx; i++) {
    int i_vtx = _face_vtx[i] - 1;
    memcpy(&_vtx_coord[3*i], &vtx_coord[3*i_vtx], sizeof(double) * 3);
    _vtx_field[i] = vtx_field[i_vtx] - isovalue;
  }


  // log_trace("_n_vtx = %d\n", _n_vtx);
  // PDM_log_trace_array_int(_face_vtx,       _n_vtx,                     "_face_vtx (unique) : ");
  // PDM_log_trace_array_int(order_in_unique, _face_vtx_idx[cell_face_n], "order_in_unique    : ");

  for (int i = 0; i < _face_vtx_idx[cell_face_n]; i++) {
    _face_vtx[i] = order_in_unique[i] + 1;
  }

  const char *name = "debug_ngon_cell.vtk";
  PDM_vtk_write_polydata_field(name,
                               _n_vtx,
                               _vtx_coord,
                               NULL,
                               cell_face_n,
                               _face_vtx_idx,
                               _face_vtx,
                               NULL,
                               NULL,
                               NULL,
                               "field",
                               _vtx_field);
  PDM_free(_face_vtx_idx);
  PDM_free(_face_vtx);
  PDM_free(order_in_unique);
  PDM_free(_vtx_coord);
  PDM_free(_vtx_field);
}


static void
_trace_isopolygon_in_cell
(
 int   cell_face_n,
 int  *cell_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 double *vtx_coord, // only for debug
 double *vtx_field, // only for debug
 double  isovalue,  // only for debug
 int  *is_active_face,
 int  *edge_to_iso_vtx,
 int  *i_edge_in_cell,
 int  *cell_edge_face,
 int  *cell_edge,
 int  *is_used_edge,
 int  *face_tag,
 int  *iso_n_face,
 int **iso_face_vtx_idx,
 int **iso_face_vtx,
 int  *tmp_iso_n_face,
 int  *iso_n_edge,
 int **iso_edge_vtx,
 int **iso_edge_parent_idx,
 int **iso_edge_parent,
 int  *tmp_iso_n_edge
 )
{
  int dbg = 0;

  // Collect active edges incident to current cell
  int cell_edge_n = 0;

  for (int i_face = 0; i_face < cell_face_n; i_face++) {
    int face_id = PDM_ABS(cell_face[i_face]) - 1;

    if (!is_active_face[face_id]) continue;

    for (int i_edge = face_edge_idx[face_id]; i_edge < face_edge_idx[face_id+1]; i_edge++) {
      int edge_id = PDM_ABS(face_edge[i_edge]) - 1;
      if (edge_to_iso_vtx[edge_id] == 0) continue;

      if (i_edge_in_cell[edge_id] < 0) {
        i_edge_in_cell[edge_id] = cell_edge_n;
        cell_edge[cell_edge_n]  = edge_id;
        cell_edge_n++;
      }

      if (PDM_SIGN(cell_face[i_face]) * PDM_SIGN(face_edge[i_edge]) * edge_to_iso_vtx[edge_id] > 0) {
        cell_edge_face[2*i_edge_in_cell[edge_id]  ] = i_face; // face for which current edge's isovtx is key (== IPISE in isoap)
        cell_edge_face[2*i_edge_in_cell[edge_id]+1] = i_edge - face_edge_idx[face_id]; // position of current edge in this face (== IVISE in isoap)
      }

    } // End of loop on edges incident to current face
  } // End of loop on faces incident to current cell

  if (dbg) {
    PDM_log_trace_array_int(cell_edge, cell_edge_n, "  cell_edge : ");
  }
  if (cell_edge_n < 3) {//== 0) {
    // Current cell is not traversed by the isosurface
    return;
  }



  // Trace isopolygons in current cell
  int n_used_edge = 0;

  int start_edge = 0;
  int safeguard1 = 0;
  while (n_used_edge < cell_edge_n) {
    safeguard1++;
    if (safeguard1 > cell_edge_n) {
      PDM_error(__FILE__, __LINE__, 0, "Failed to use all cell edges\n");
    }
    // Start a new isopolygon
    if (*iso_n_face >= *tmp_iso_n_face) {
      *tmp_iso_n_face = PDM_MAX(*iso_n_face+1, 2*(*tmp_iso_n_face));
      PDM_realloc(*iso_face_vtx_idx, *iso_face_vtx_idx, *tmp_iso_n_face + 1, int);
    }

    for (int i = 0; i < cell_edge_n; i++) {
      if (!is_used_edge[i]) {
        start_edge = i;
        break;
      }
    }
    n_used_edge++;
    assert(!is_used_edge[start_edge]);
    is_used_edge[start_edge] = 1;

    int *ifv = *iso_face_vtx + (*iso_face_vtx_idx)[*iso_n_face];
    int s_ifv = 0;

    ifv[s_ifv++] = PDM_ABS(edge_to_iso_vtx[cell_edge[start_edge]]);

    int current_cell_edge = start_edge;
    int current_edge      = cell_edge[current_cell_edge];
    int current_dest      = -1;

    if (dbg) {
      log_trace("start new isopoly from edge %d (#%d in cell, isovtx = %d)\n",
                current_edge, current_cell_edge, ifv[0]);
    }

    // Complete isopolygon
    int is_complete = 0;
    int safeguard2  = 0;
    while (!is_complete) {
      safeguard2++;
      if (safeguard2 > cell_edge_n) {
        PDM_error(__FILE__, __LINE__, 0, "Failed to complete isopolygon\n");
      }
      // Get face for which current edge's isovtx is "key"
      int i_face = cell_edge_face[2*current_cell_edge];
      int face_id = PDM_ABS(cell_face[i_face]) - 1;
      int *fe = face_edge + face_edge_idx[face_id];
      int face_edge_n = face_edge_idx[face_id+1] - face_edge_idx[face_id];
      if (dbg) {
        log_trace("  face_id : %d (sign = %d)\n", face_id, PDM_SIGN(cell_face[i_face]));
      }


      // Get position of current edge in this face
      int i_edge = cell_edge_face[2*current_cell_edge+1];

      if (dbg) {
        log_trace("    i_edge : %d (fe = %d)\n", i_edge, fe[i_edge]);
      }

      // Get next *active* edge in this face
      // (Be careful: adjacent edges in a face are not necessarily contiguous in the face_edge connectivity table!)

      // destination vertex of current edge in this face
      if (PDM_SIGN(cell_face[i_face]) * PDM_SIGN(fe[i_edge]) > 0) {
        current_dest = edge_vtx[2*current_edge+1];
      }
      else {
        current_dest = edge_vtx[2*current_edge  ];
      }


      if (dbg) {
        log_trace("    current_dest : %d\n", current_dest);
      }

      int found_next = 0;
      // "external" loop to find the next *active* edge in current face
      for (int k = 0; k < face_edge_n; k++) {
        if (dbg) {
          log_trace("     k = %d\n", k);
        }

        // "internal" loop to find the (unique) edge in current face that starts at vertex current_dest
        for (int j = 0; j < face_edge_n; j++) {
          int edge_id = PDM_ABS(fe[j]) - 1;//(i_edge+j)%face_edge_n]) - 1;
          if (edge_id == current_edge) continue;
          if (dbg) {
            log_trace("      j = %d, edge_id = %d, isovtx : %d\n", j, edge_id, edge_to_iso_vtx[edge_id]);
          }
          int j_orig, j_dest;
          if (PDM_SIGN(cell_face[i_face]) * PDM_SIGN(fe[j]) > 0) {
            j_orig = edge_vtx[2*edge_id  ];
            j_dest = edge_vtx[2*edge_id+1];
          }
          else {
            j_orig = edge_vtx[2*edge_id+1];
            j_dest = edge_vtx[2*edge_id  ];
          }

          if (dbg) {
            log_trace("      j_orig/dest : %d %d\n", j_orig, j_dest);
          }

          if (j_orig == current_dest) {
            // we found next edge
            current_dest = j_dest;
            current_edge = edge_id;
            if (dbg) {
              log_trace("      next edge found : %d\n", current_edge);
            }

            if (edge_to_iso_vtx[edge_id] != 0) {
              // this is an active edge
              current_cell_edge = i_edge_in_cell[edge_id];
              found_next        = 1;
              if (dbg) {
                log_trace("      this is an active edge :)\n");
              }


              int iso_edge_orig = ifv[s_ifv-1];
              int iso_edge_dest = -1;
              if (current_cell_edge == start_edge) {
                is_complete = 1;
                iso_edge_dest = ifv[0];
              }
              else {
                // assert(is_used_edge[current_cell_edge] == 0);
                if (is_used_edge[current_cell_edge]) {
                  int i_wrong_edge = cell_edge[current_cell_edge];
                  if (0) {
                    _debug_ngon_cell(cell_face_n,
                                     cell_face,
                                     face_edge_idx,
                                     face_edge,
                                     edge_vtx,
                                     vtx_coord,
                                     vtx_field,
                                     isovalue);
                  }
                  PDM_error(__FILE__, __LINE__, 0, "Edge %d : (%d %d) already used\n",
                            i_wrong_edge,
                            edge_vtx[2*i_wrong_edge  ],
                            edge_vtx[2*i_wrong_edge+1]);
                }
                is_used_edge[current_cell_edge] = 1;
                n_used_edge++;

                if (dbg) {
                  log_trace("      n_used_edge = %d / %d\n", n_used_edge, cell_edge_n);
                }

                iso_edge_dest = PDM_ABS(edge_to_iso_vtx[current_edge]);
                if (iso_edge_dest != ifv[s_ifv-1]) {
                  ifv[s_ifv++] = iso_edge_dest;
                }
              }

              /* Add iso_edge if necessary */
              if (face_tag != NULL) {
                if (face_tag[face_id] > 0) {
                  if (iso_edge_dest > 0) {
                    // add 1 edge
                    if (*iso_n_edge >= *tmp_iso_n_edge) {
                      *tmp_iso_n_edge = PDM_MAX(*iso_n_edge+1, 2*(*tmp_iso_n_edge));
                      PDM_realloc(*iso_edge_parent_idx, *iso_edge_parent_idx, *tmp_iso_n_edge + 1, int);
                      PDM_realloc(*iso_edge_parent    , *iso_edge_parent,     *tmp_iso_n_edge    , int); // Warning : fails if more than one parent
                      PDM_realloc(*iso_edge_vtx       , *iso_edge_vtx,        *tmp_iso_n_edge * 2, int);
                    }
                    (*iso_edge_parent_idx)[*iso_n_edge+1] = (*iso_edge_parent_idx)[*iso_n_edge] + 1;
                    (*iso_edge_parent)[(*iso_edge_parent_idx)[*iso_n_edge]] = face_id + 1;

                    (*iso_edge_vtx)[2*(*iso_n_edge)  ] = iso_edge_orig;
                    (*iso_edge_vtx)[2*(*iso_n_edge)+1] = iso_edge_dest;
                    (*iso_n_edge)++;
                  }
                }
              }

            }
            break;
          }
        } // end of search for next active edge

        if (found_next) break;
      }

      assert(found_next);
    } // end of current isopolygon

    // make sure last isovtx is not identical to first isovtx (can happen in degenerate cases)
    if (ifv[s_ifv-1] == ifv[0]) {
      s_ifv--;
    }

    // last bnd edge if


    if (s_ifv > 2) {
      // Non-degenerate isopolygon
      if (dbg) {
        log_trace("complete isopoly : ");
        PDM_log_trace_array_int(ifv, s_ifv, "");
        log_trace("\n");
      }
      (*iso_face_vtx_idx)[*iso_n_face+1] = (*iso_face_vtx_idx)[*iso_n_face] + s_ifv;
      (*iso_n_face)++;
    }
    else {
      // TODO: cancel iso_edges??
    }

  } // end of current cell


  // Reset
  for (int i_edge = 0; i_edge < cell_edge_n; i_edge++) {
    i_edge_in_cell[cell_edge[i_edge]] = -1;
    is_used_edge[i_edge] = 0;
  }
}


static void
_isosurface_ngon_single_part
(
  int          mesh_dimension,
  int          n_isovalues,
  double      *isovalues,
  double       tol,
  int          n_cell,
  int          n_face,
  int          n_edge,
  int          n_vtx,
  int         *cell_face_idx,
  int         *cell_face,
  int         *face_edge_idx,
  int         *face_edge,
  int         *edge_vtx,
  double      *vtx_coord,
  double      *vtx_field,
  int         *face_tag,
  PDM_g_num_t *cell_ln_to_gn,
  int         *out_iso_n_vtx,
  double     **out_iso_vtx_coord,
  int        **out_iso_vtx_parent_idx,
  int        **out_iso_vtx_parent,
  double     **out_iso_vtx_parent_weight,
  int        **out_iso_vtx_parent_edge,
  int        **out_isovalue_vtx_idx,
  int         *out_iso_n_face,
  int        **out_iso_face_vtx_idx,
  int        **out_iso_face_vtx,
  int        **out_iso_face_parent_idx,
  int        **out_iso_face_parent,
  int        **out_isovalue_face_idx,
  int         *out_iso_n_edge,
  int        **out_iso_edge_vtx,
  int        **out_iso_edge_parent_idx,
  int        **out_iso_edge_parent,
  int        **out_isovalue_edge_idx
)
{
  int dbg = 0;
  int is_3d = (mesh_dimension == 3);

  /* Count isosurface vertices */
  int iso_n_vtx = 0;
  for (int i_edge = 0; i_edge < n_edge; i_edge++) {
    int i_vtx0 = edge_vtx[2*i_edge  ] - 1;
    int i_vtx1 = edge_vtx[2*i_edge+1] - 1;

    iso_n_vtx += _isosurface_cross_any_level_ngon(vtx_field[i_vtx0],
                                                  vtx_field[i_vtx1],
                                                  n_isovalues,
                                                  isovalues,
                                                  tol);
  }

  double *iso_vtx_coord         = NULL;
  int    *iso_vtx_parent_idx    = NULL;
  int    *iso_vtx_parent        = NULL;
  int    *iso_vtx_parent_edge   = NULL;
  double *iso_vtx_parent_weight = NULL;
  PDM_malloc(iso_vtx_coord        , iso_n_vtx * 3, double);
  PDM_malloc(iso_vtx_parent_idx   , iso_n_vtx + 1, int   );
  PDM_malloc(iso_vtx_parent       , iso_n_vtx * 2, int   );
  PDM_malloc(iso_vtx_parent_edge  , iso_n_vtx    , int   );
  PDM_malloc(iso_vtx_parent_weight, iso_n_vtx * 2, double);
  iso_vtx_parent_idx[0] = 0;

  if (iso_n_vtx == 0) {
    /* Empty isosurface => early return */
    *out_iso_n_vtx             = 0;
    *out_iso_vtx_coord         = iso_vtx_coord;
    *out_iso_vtx_parent_idx    = iso_vtx_parent_idx;
    *out_iso_vtx_parent        = iso_vtx_parent;
    *out_iso_vtx_parent_edge   = iso_vtx_parent_edge;
    *out_iso_vtx_parent_weight = iso_vtx_parent_weight;

    if (is_3d) {
      *out_iso_n_face          = 0;
      *out_iso_face_parent_idx = PDM_array_zeros_int(1);
      PDM_malloc(*out_iso_face_parent, 0, int);
      *out_iso_face_vtx_idx    = PDM_array_zeros_int(1);
      PDM_malloc(*out_iso_face_vtx, 0, int);
      *out_isovalue_face_idx   = PDM_array_zeros_int(n_isovalues + 1);
    }

    if (face_tag != NULL) {
      *out_iso_n_edge          = 0;
      *out_iso_edge_parent_idx = PDM_array_zeros_int(1);
      PDM_malloc(*out_iso_edge_parent, 0, int);
      PDM_malloc(*out_iso_edge_vtx   , 0, int);
      *out_isovalue_edge_idx   = PDM_array_zeros_int(n_isovalues + 1);
    }

    return;
  }


  iso_n_vtx = 0;

  int max_cell_face_n = 0;
  for (int i_cell = 0; i_cell < n_cell; i_cell++) {
    max_cell_face_n = PDM_MAX(max_cell_face_n, cell_face_idx[i_cell+1] - cell_face_idx[i_cell]);
  }

  int max_face_edge_n = 0;
  for (int i_face = 0; i_face < n_face; i_face++) {
    int face_edge_n = face_edge_idx[i_face+1] - face_edge_idx[i_face];
    max_face_edge_n = PDM_MAX(max_face_edge_n, face_edge_n);
  }


  /* Loop over isovalues */
  int *edge_to_iso_vtx = NULL;
  int *is_active_face  = NULL;
  PDM_malloc(edge_to_iso_vtx, n_edge, int);
  PDM_malloc(is_active_face , n_face, int);

  int *i_edge_in_cell   = PDM_array_const_int(n_edge, -1);
  int s_cell_edge       = max_cell_face_n * max_face_edge_n;
  int *cell_edge        = NULL;
  int *cell_edge_face   = NULL;
  PDM_malloc(cell_edge     , s_cell_edge  , int);
  PDM_malloc(cell_edge_face, s_cell_edge*2, int);
  int *is_used_edge     = PDM_array_zeros_int(s_cell_edge);
  // TODO (optim) : malloc a single, larger buffer for work arrays

  int s_iso_face           = n_cell;
  int s_iso_face_vtx       = n_cell * s_cell_edge;
  int *iso_face_vtx_idx    = NULL;
  int *iso_face_vtx        = NULL;
  int *iso_face_parent_idx = NULL;
  int *iso_face_parent     = NULL;
  if (is_3d) {
    PDM_malloc(iso_face_vtx_idx   , s_iso_face + 1, int);
    PDM_malloc(iso_face_vtx       , s_iso_face_vtx, int);
    PDM_malloc(iso_face_parent_idx, s_iso_face + 1, int);
    PDM_malloc(iso_face_parent    , s_iso_face    , int);
    iso_face_parent_idx[0] = 0;
    iso_face_vtx_idx   [0] = 0;
}

  int *isovalue_face_idx = NULL;
  if (is_3d) {
    PDM_malloc(isovalue_face_idx, n_isovalues + 1, int);
    isovalue_face_idx[0] = 0;
  }

  int iso_n_face = 0;


  int s_iso_edge           = 0;
  int *iso_edge_vtx        = NULL;
  int *iso_edge_parent_idx = NULL;
  int *iso_edge_parent     = NULL;

  int *isovalue_edge_idx = NULL;

  int iso_n_edge = 0;
  if (face_tag != NULL) {
    PDM_malloc(isovalue_edge_idx, n_isovalues + 1, int);
    isovalue_edge_idx[0] = 0;
    for (int i_face = 0; i_face < n_face; i_face++) {
      s_iso_edge += n_isovalues * (face_tag[i_face] > 0);
    }
    PDM_malloc(iso_edge_vtx       , s_iso_edge * 2, int);
    PDM_malloc(iso_edge_parent_idx, s_iso_edge + 1, int);
    PDM_malloc(iso_edge_parent    , s_iso_edge * 2, int);

    iso_edge_parent_idx[0] = 0;
  }

  int *vtx_to_iso_vtx = NULL;
  PDM_malloc(vtx_to_iso_vtx, n_vtx, int);

  int *isovalue_vtx_idx = NULL;
  PDM_malloc(isovalue_vtx_idx, n_isovalues + 1, int);
  isovalue_vtx_idx[0] = 0;

  for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {

    /* Identify active edges for current isovalue */
    // maybe optimize reset?
    PDM_array_reset_int(edge_to_iso_vtx, n_edge, 0);
    PDM_array_reset_int(vtx_to_iso_vtx,  n_vtx,  0);

    int isovalue_n_active_edge = 0;

    for (int i_edge = 0; i_edge < n_edge; i_edge++) {
      int i_vtx0 = edge_vtx[2*i_edge  ] - 1;
      int i_vtx1 = edge_vtx[2*i_edge+1] - 1;

      double val0 = vtx_field[i_vtx0] - isovalues[i_isovalue];
      double val1 = vtx_field[i_vtx1] - isovalues[i_isovalue];

      if (_isosurface_cross_0_level_ngon(val0, val1, tol)) {

        isovalue_n_active_edge++;

        double t = val0 / (val0 - val1);
        for (int i = 0; i < 3; i++) {
          iso_vtx_coord[3*iso_n_vtx+i] = (1-t)*vtx_coord[3*i_vtx0+i] + t*vtx_coord[3*i_vtx1+i];
        }

        // Reduce t == 0/1 to a single parent
        int iso_vtx_id = 0;
        if (t < ISOSURFACE_EPS_T) { // if (PDM_ABS(v0) < ISOSURFACE_EPS) ??
          if (vtx_to_iso_vtx[i_vtx0] == 0) {
            iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]] = i_vtx0 + 1;
            iso_vtx_parent_edge  [                   iso_n_vtx ] = i_edge;
            iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]] = 1.;
            iso_vtx_parent_idx[iso_n_vtx+1] = iso_vtx_parent_idx[iso_n_vtx] + 1;
            iso_vtx_id = ++iso_n_vtx;

            vtx_to_iso_vtx[i_vtx0] = iso_vtx_id;
          }
          else {
            iso_vtx_id = vtx_to_iso_vtx[i_vtx0];
          }
        }
        else if (t > 1 - ISOSURFACE_EPS_T) { // if (PDM_ABS(v1) < ISOSURFACE_EPS) ??
          if (vtx_to_iso_vtx[i_vtx1] == 0) {
            iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]] = i_vtx1 + 1;
            iso_vtx_parent_edge  [                   iso_n_vtx ] = i_edge;
            iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]] = 1.;
            iso_vtx_parent_idx[iso_n_vtx+1] = iso_vtx_parent_idx[iso_n_vtx] + 1;
            iso_vtx_id = ++iso_n_vtx;

            vtx_to_iso_vtx[i_vtx1] = iso_vtx_id;
          }
          else {
            iso_vtx_id = vtx_to_iso_vtx[i_vtx1];
          }
        }
        else {
          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]  ] = i_vtx0 + 1;
          iso_vtx_parent_edge  [                   iso_n_vtx   ] = i_edge;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]  ] = 1. - t;

          iso_vtx_parent       [iso_vtx_parent_idx[iso_n_vtx]+1] = i_vtx1 + 1;
          iso_vtx_parent_weight[iso_vtx_parent_idx[iso_n_vtx]+1] = t;
          iso_vtx_parent_idx[iso_n_vtx+1] = iso_vtx_parent_idx[iso_n_vtx] + 2;

          iso_vtx_id = ++iso_n_vtx;
        }

        if (dbg) {
          log_trace("iso vtx %d on edge %d : t = %f, xyz = %f %f %f, parent(s) :",
                    iso_vtx_id, i_edge, t,
                    iso_vtx_coord[3*(iso_vtx_id-1)  ],
                    iso_vtx_coord[3*(iso_vtx_id-1)+1],
                    iso_vtx_coord[3*(iso_vtx_id-1)+2]);
          for (int idx = iso_vtx_parent_idx[iso_vtx_id-1]; idx < iso_vtx_parent_idx[iso_vtx_id]; idx++) {
            log_trace(" %d", iso_vtx_parent[idx]);
          }
          log_trace(")\n");
        }

        int sign = val0 > val1 ? 1 : -1;
        edge_to_iso_vtx[i_edge] = sign*iso_vtx_id;
      }
    } // End of loop on edges

    if (dbg && 0) {
      const char *field_name [] = {"edge_to_iso_vtx"};
      const int  *field_value[] = {edge_to_iso_vtx};

      char filename[999];
      sprintf(filename, "active_edges_%d.vtk", i_isovalue);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 vtx_coord,
                                 NULL,
                                 PDM_MESH_NODAL_BAR2,
                                 n_edge,
                                 edge_vtx,
                                 NULL,
                                 1,
                                 field_name,
                                 field_value);
    }

    // Skip isovalue if zero active edge
    if (isovalue_n_active_edge == 0) {
      isovalue_vtx_idx [i_isovalue+1] = iso_n_vtx;
      isovalue_face_idx[i_isovalue+1] = iso_n_face;
      if (face_tag != NULL) {
        isovalue_edge_idx[i_isovalue+1] = iso_n_edge;
      }
      continue;
    }

    /* Identify active faces for current isovalue (with at least one active edge) */
    PDM_array_reset_int(is_active_face, n_face, 0);
    for (int i_face = 0; i_face < n_face; i_face++) {
      for (int i_edge = face_edge_idx[i_face]; i_edge < face_edge_idx[i_face+1]; i_edge++) {
        int edge_id = PDM_ABS(face_edge[i_edge]) - 1;
        if (edge_to_iso_vtx[edge_id] != 0) {
          is_active_face[i_face] = 1;
          break;
        }
      }
    }

    /* Polygonize isovalue */
    if (dbg && 0) {
      _debug_ngon(n_face,
                  n_edge,
                  face_edge_idx,
                  face_edge,
                  edge_vtx,
                  vtx_coord,
                  vtx_field,
                  isovalues[i_isovalue]);
    }
    for (int i_cell = 0; i_cell < n_cell; i_cell++) {

      int cell_face_n = cell_face_idx[i_cell+1] - cell_face_idx[i_cell];

      int tmp_s_iso_face = s_iso_face;
      int tmp_iso_n_face = iso_n_face;
      if (dbg) {
        log_trace("cell %d ("PDM_FMT_G_NUM"):\n", i_cell, cell_ln_to_gn[i_cell]);
      }
      _trace_isopolygon_in_cell(cell_face_n,
                                cell_face + cell_face_idx[i_cell],
                                face_edge_idx,
                                face_edge,
                                edge_vtx,
                                vtx_coord,
                                vtx_field,
                                isovalues[i_isovalue],
                                is_active_face,
                                edge_to_iso_vtx,
                                i_edge_in_cell,
                                cell_edge_face,
                                cell_edge,
                                is_used_edge,
                                face_tag,
                                &iso_n_face,
                                &iso_face_vtx_idx,
                                &iso_face_vtx,
                                &s_iso_face,
                                &iso_n_edge,
                                &iso_edge_vtx,
                                &iso_edge_parent_idx,
                                &iso_edge_parent,
                                &s_iso_edge);

      if (iso_n_face > tmp_s_iso_face) {
        // TODO: handle multiple parents??
        PDM_realloc(iso_face_parent_idx, iso_face_parent_idx, s_iso_face + 1, int);
        PDM_realloc(iso_face_parent    , iso_face_parent    , s_iso_face    , int);
      }
      for (int i = tmp_iso_n_face; i < iso_n_face; i++) {
        iso_face_parent_idx[i+1] = iso_face_parent_idx[i] + 1;
        iso_face_parent[iso_face_parent_idx[i]] = i_cell + 1;
      }

    } // End of loop on cells

    isovalue_vtx_idx [i_isovalue+1] = iso_n_vtx;
    if (is_3d) {
      isovalue_face_idx[i_isovalue+1] = iso_n_face;
    }
    if (face_tag != NULL) {
      isovalue_edge_idx[i_isovalue+1] = iso_n_edge;
    }

  } // End of loop on isovalues


  /* Free memory */
  PDM_free(edge_to_iso_vtx);
  PDM_free(vtx_to_iso_vtx );
  PDM_free(is_active_face );
  PDM_free(i_edge_in_cell );
  PDM_free(cell_edge      );
  PDM_free(cell_edge_face );
  PDM_free(is_used_edge   );

  /* Outputs */
  *out_iso_n_vtx = iso_n_vtx;
  PDM_realloc(iso_vtx_coord        , *out_iso_vtx_coord        ,                           iso_n_vtx * 3, double);
  PDM_realloc(iso_vtx_parent_idx   , *out_iso_vtx_parent_idx   ,                           iso_n_vtx + 1, int   );
  PDM_realloc(iso_vtx_parent       , *out_iso_vtx_parent       , (*out_iso_vtx_parent_idx)[iso_n_vtx]   , int   );
  PDM_realloc(iso_vtx_parent_edge  , *out_iso_vtx_parent_edge  ,                           iso_n_vtx    , int   );
  PDM_realloc(iso_vtx_parent_weight, *out_iso_vtx_parent_weight, (*out_iso_vtx_parent_idx)[iso_n_vtx]   , double);
  *out_isovalue_vtx_idx = isovalue_vtx_idx;

  if (is_3d) {
    *out_iso_n_face = iso_n_face;
    PDM_realloc(iso_face_parent_idx, *out_iso_face_parent_idx,   iso_n_face + 1                      , int);
    PDM_realloc(iso_face_parent    , *out_iso_face_parent    , (*out_iso_face_parent_idx)[iso_n_face], int);
    PDM_realloc(iso_face_vtx_idx   , *out_iso_face_vtx_idx   ,   iso_n_face + 1                      , int);
    PDM_realloc(iso_face_vtx       , *out_iso_face_vtx       , (*out_iso_face_vtx_idx)   [iso_n_face], int);
    *out_isovalue_face_idx = isovalue_face_idx;
  }

  if (face_tag != NULL) { // TODO dim 2
    *out_iso_n_edge          = iso_n_edge;
    *out_iso_edge_parent_idx = NULL;
    *out_iso_edge_parent     = NULL;
    *out_iso_edge_vtx        = NULL;
    PDM_realloc(iso_edge_parent_idx, *out_iso_edge_parent_idx,   iso_n_edge + 1                      , int);
    PDM_realloc(iso_edge_parent    , *out_iso_edge_parent    , (*out_iso_edge_parent_idx)[iso_n_edge], int);
    PDM_realloc(iso_edge_vtx       , *out_iso_edge_vtx       ,   iso_n_edge * 2                      , int);
    *out_isovalue_edge_idx = isovalue_edge_idx;
  }
}


static void
_generate_gnum_from_parents
(
 PDM_MPI_Comm  comm,
 int           n_part,
 int          *n_entity,
 PDM_g_num_t **parent_gnum,
 int           nuplet,
 PDM_g_num_t **gnum
)
{
  PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,  // unused,
                                             n_part,
                                             PDM_FALSE,
                                             1., // unused,
                                             comm,
                                             PDM_OWNERSHIP_USER);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gen_gnum, i_part, n_entity[i_part], parent_gnum[i_part]);
  }

  PDM_gnum_set_parents_nuplet(gen_gnum, nuplet);

  PDM_gnum_compute(gen_gnum);

  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum[i_part] = PDM_gnum_get(gen_gnum, i_part);
    // PDM_free(parent_gnum[i_part]);
  }
  // PDM_free(parent_gnum);
  PDM_gnum_free(gen_gnum);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_isosurface_marching_algo
(
  PDM_isosurface_t        *isos,
  int                      id_iso
)
{
  int debug      = 1;
  int debug_visu = 1;
  double t_start, t_end;

  if (debug==1) {
    log_trace("id_iso = %d\n", id_iso);
  }

  // > Comm info
  PDM_MPI_Comm comm = isos->comm;
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // > Isosurface information
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];
  int   n_isovalues = _iso->n_isovalues;
  double *isovalues = _iso->  isovalues;

  double **vtx_field = isos->extract_field;

  // > Mesh information
  int n_part = isos->n_part;

  // > Allocate isosurface vertices
  int          *iso_n_vtx             = PDM_array_zeros_int(n_part);
  double      **iso_vtx_coord         = NULL;
  int         **iso_vtx_parent_idx    = NULL;
  int         **iso_vtx_parent        = NULL;
  double      **iso_vtx_parent_weight = NULL;
  PDM_g_num_t **iso_vtx_parent_gnum   = NULL;
  int         **isovalue_vtx_idx      = NULL;
  PDM_malloc(iso_vtx_coord        , n_part, double      *);
  PDM_malloc(iso_vtx_parent_idx   , n_part, int         *);
  PDM_malloc(iso_vtx_parent       , n_part, int         *);
  PDM_malloc(iso_vtx_parent_weight, n_part, double      *);
  PDM_malloc(iso_vtx_parent_gnum  , n_part, PDM_g_num_t *);
  PDM_malloc(isovalue_vtx_idx     , n_part, int         *);

  // > Allocate isosurface edges
  int          *iso_n_edge           = PDM_array_zeros_int(n_part);
  int         **iso_edge_vtx         = NULL;
  int         **iso_edge_bnd_tag_idx = NULL;
  int         **iso_edge_bnd_tag     = NULL;
  PDM_malloc(iso_edge_vtx        , n_part, int *);
  PDM_malloc(iso_edge_bnd_tag_idx, n_part, int *);
  PDM_malloc(iso_edge_bnd_tag    , n_part, int *);
  int           iso_n_edge_group     = 0;
  int         **iso_edge_group_idx   = NULL;
  int         **iso_edge_group_lnum  = NULL;
  PDM_g_num_t **iso_edge_parent_gnum = NULL;
  int         **iso_edge_parent_idx  = NULL;
  int         **iso_edge_parent      = NULL;
  int         **isovalue_edge_idx    = NULL;
  PDM_malloc(iso_edge_group_idx  , n_part, int         *);
  PDM_malloc(iso_edge_group_lnum , n_part, int         *);
  PDM_malloc(iso_edge_parent_gnum, n_part, PDM_g_num_t *);
  PDM_malloc(iso_edge_parent_idx , n_part, int         *);
  PDM_malloc(iso_edge_parent     , n_part, int         *);
  PDM_malloc(isovalue_edge_idx   , n_part, int         *);

  // > Allocate face edges
  // TODO: transformer en part_mesh_nodal (tetra sera que TRI et QUAD)
  //       quid des autres elements que les TETRA ? POLY_2D ? (MR gitlab)
  int          *iso_n_face           = NULL;
  int         **iso_face_vtx_idx     = NULL;
  int         **iso_face_vtx         = NULL;
  PDM_g_num_t **iso_face_parent_gnum = NULL;
  int         **iso_face_parent_idx  = NULL;
  int         **iso_face_parent      = NULL;
  int         **isovalue_face_idx    = NULL;
  if (isos->entry_mesh_dim==3) {
    iso_n_face = PDM_array_zeros_int(n_part);
    PDM_malloc(iso_face_vtx_idx    , n_part, int         *);
    PDM_malloc(iso_face_vtx        , n_part, int         *);
    PDM_malloc(iso_face_parent_gnum, n_part, PDM_g_num_t *);
    PDM_malloc(iso_face_parent_idx , n_part, int         *);
    PDM_malloc(iso_face_parent     , n_part, int         *);
    PDM_malloc(isovalue_face_idx   , n_part, int         *);
  }


  for (int i_part = 0; i_part < n_part; i_part++) {

    // > Get vtx info
    int          n_vtx     = isos->extract_n_vtx    [i_part];
    double      *vtx_coord = isos->extract_vtx_coord[i_part];
    PDM_g_num_t *vtx_gnum  = isos->extract_vtx_gnum [i_part];
    

    /*
     * Build active edges from entry mesh (i.e. crossing or on at least one isosurface)
     */
    int  n_edge           = 0;
    int *edge_vtx         = NULL;
    int *edge_bnd_tag_idx = NULL;
    int *edge_bnd_tag     = NULL;
    int *edge_parent_idx  = NULL;
    int *edge_parent      = NULL;
    int  n_crossings      = 0;

    // > Active edges for each sections
    int   n_section = 2;
    int **elt_edge  = NULL;
    PDM_malloc(elt_edge, n_section, int *);
    for (int i_section = 0; i_section < n_section; i_section++) {
      elt_edge[i_section] = NULL;
    }

    PDM_Mesh_nodal_elt_t sections_elt_t   [2] = {PDM_MESH_NODAL_TRIA3,
                                                 PDM_MESH_NODAL_TETRA4};
    int                  sections_n_elt   [2] = {isos->extract_n_tri     [i_part],
                                                 isos->extract_n_tetra   [i_part]};
    int                 *sections_elt_vtx [2] = {isos->extract_tri_vtx   [i_part],
                                                 isos->extract_tetra_vtx [i_part]};
    int                 *sections_elt_tag [2] = {isos->extract_tri_tag   [i_part],
                                                 NULL};

    if (debug_visu) {
      char name[999];
      const char   *section_name[] = {"tria", "tetra"};
      const char   *field_name  [] = {"field"};
      const double *field_value [] = {vtx_field[i_part]};
      for (int i_section = 0; i_section < 2; i_section++) {
        sprintf(name, "extract_%s_%d.vtk", section_name[i_section], i_rank*n_part + i_part);
        PDM_vtk_write_std_elements_ho_with_vtx_field(name,
                                                     1,
                                                     n_vtx,
                                                     vtx_coord,
                                                     vtx_gnum,
                                                     sections_elt_t  [i_section],
                                                     sections_n_elt  [i_section],
                                                     sections_elt_vtx[i_section],
                                                     NULL,
                                                     0,
                                                     NULL,
                                                     NULL,
                                                     1,
                                                     field_name,
                                                     field_value);
      }
    }

    t_start = PDM_MPI_Wtime();
    _build_active_edges(n_section,
                        sections_elt_t,
                        sections_n_elt,
                        sections_elt_vtx,
                        sections_elt_tag,
                        n_isovalues,
                        isovalues,
                        isos->ISOSURFACE_EPS,
                        vtx_field[i_part],
                       &n_edge,
                       &edge_vtx,
                       &edge_bnd_tag_idx,
                       &edge_bnd_tag,
                       &edge_parent_idx,
                       &edge_parent,
                        elt_edge,
                       &n_crossings);
    // TODO: merge iso_edge_parent with parallel

    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
      // should we do this here??
      for (int i = 0; i < iso_edge_parent_idx[i_part][iso_n_edge[i_part]]; i++) {
        iso_edge_parent[i_part][i] = isos->extract_tri_lnum[i_part][iso_edge_parent[i_part][i]-1] + 1;
      }
    }

    t_end = PDM_MPI_Wtime();

    if (debug==1) {
      log_trace("\n");
      log_trace("Build active edges : %.3fs  (%d edges)\n", t_end - t_start, n_edge);
      int n_bnd_tag_tot = edge_bnd_tag_idx[n_edge];
      PDM_log_trace_array_int(edge_bnd_tag_idx, n_edge+1     , "edge_bnd_tag_idx ::");
      PDM_log_trace_array_int(edge_bnd_tag    , n_bnd_tag_tot, "edge_bnd_tag     ::");
      int n_parent_tot = edge_parent_idx[n_edge];
      PDM_log_trace_array_int(edge_parent_idx, n_edge+1    , "edge_parent_idx ::");
      PDM_log_trace_array_int(edge_parent    , n_parent_tot, "edge_parent     ::");
    }

    if (debug_visu==1) {
      char outname[999];
      sprintf(outname, "marching_sections_active_edges_%d_%d.vtk", i_rank, i_part);
      PDM_vtk_write_std_elements(outname,
                                 n_vtx,
                                 vtx_coord,
                                 NULL,
                                 PDM_MESH_NODAL_BAR2,
                                 n_edge,
                                 edge_vtx,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }



    /*
     * Build vertices on actives edges
     */
    iso_vtx_coord        [i_part] = NULL;
    iso_vtx_parent_idx   [i_part] = NULL;
    iso_vtx_parent       [i_part] = NULL;
    iso_vtx_parent_weight[i_part] = NULL;
    iso_vtx_parent_gnum  [i_part] = NULL;
    int *vtx_to_iso_vtx           = NULL;
    int *iso_vtx_to_edge          = NULL;

    // > Build iso-vertices on active edges
    iso_n_vtx[i_part] = 0;
    t_start = PDM_MPI_Wtime();

    _build_iso_vtx(n_isovalues, 
                   isovalues,
                   isos->ISOSURFACE_EPS,
                   n_vtx,
                   vtx_gnum,
                   vtx_coord,
                   vtx_field[i_part],
                   n_edge,
                   edge_vtx,
                   n_crossings,
                  &vtx_to_iso_vtx,
                  &iso_n_vtx            [i_part],
                  &iso_vtx_coord        [i_part],
                  &iso_vtx_parent_idx   [i_part],
                  &iso_vtx_parent       [i_part],
                  &iso_vtx_parent_weight[i_part],
                  &iso_vtx_parent_gnum  [i_part],
                  &iso_vtx_to_edge,
                  &isovalue_vtx_idx     [i_part]);

    if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
      // should we do this here??
      for (int i = 0; i < iso_vtx_parent_idx[i_part][iso_n_vtx[i_part]]; i++) {
        iso_vtx_parent[i_part][i] = isos->extract_vtx_lnum[i_part][iso_vtx_parent[i_part][i]-1] + 1;
      }
    }

    t_end = PDM_MPI_Wtime();

    if (debug==1) {
      log_trace("\n");
      log_trace("Build active vertices : %.3fs  (%d vertices)\n", t_end - t_start, iso_n_vtx[i_part]);
      
      int n_parent_tot = iso_vtx_parent_idx[i_part][iso_n_vtx[i_part]];
      PDM_log_trace_array_int (iso_vtx_parent_idx [i_part],   iso_n_vtx[i_part] , "iso_vtx_parent_idx ::");
      PDM_log_trace_array_int (iso_vtx_parent     [i_part],   n_parent_tot      , "iso_vtx_parent     ::");
      PDM_log_trace_array_long(iso_vtx_parent_gnum[i_part], 3*iso_n_vtx[i_part] , "iso_vtx_parent_gnum::");
      PDM_log_trace_array_int (iso_vtx_to_edge            ,   iso_n_vtx[i_part] , "iso_vtx_to_edge    ::");
      PDM_log_trace_array_int (vtx_to_iso_vtx             ,   n_vtx             , "vtx_to_iso_vtx     ::");
      PDM_log_trace_array_int (isovalue_vtx_idx   [i_part],   n_isovalues+1     , "isovalue_vtx_idx   ::");
    }



    /*
     * Build active faces from entry mesh (i.e. on at least one isosurface)
     */
    int  n_face           = 0;
    int *face_vtx_idx     = NULL;
    int *face_vtx         = NULL;
    int *face_parent_idx  = NULL;
    int *face_parent      = NULL;

    // > Active faces for each sections
    int *elt_face = NULL;

    if (isos->entry_mesh_dim==3) {
      t_start = PDM_MPI_Wtime();
      _build_active_faces(isos->extract_n_tetra  [i_part],
                          isos->extract_tetra_vtx[i_part],
                          n_isovalues,
                          isovalues,
                          isos->ISOSURFACE_EPS,
                          vtx_field[i_part],
                          vtx_to_iso_vtx,
                          &n_face,
                          &face_vtx_idx,
                          &face_vtx,
                          &face_parent_idx,
                          &face_parent,
                          &elt_face);

      if (isos->extract_kind == PDM_EXTRACT_PART_KIND_LOCAL) {
        // should we do this here??
        for (int i = 0; i < iso_face_parent_idx[i_part][iso_n_face[i_part]]; i++) {
          iso_face_parent[i_part][i] = isos->extract_tetra_lnum[i_part][iso_face_parent[i_part][i]-1] + 1;
        }
      }


      t_end = PDM_MPI_Wtime();
      if (debug==1) {
        log_trace("\n");
        log_trace("Build active faces : %.3fs  (%d faces)\n", t_end - t_start, n_face);
        int face_vtx_size = face_vtx_idx[n_face];
        PDM_log_trace_array_int(face_vtx_idx, n_face+1     , "face_vtx_idx ::");
        PDM_log_trace_array_int(face_vtx    , face_vtx_size, "face_vtx     ::");
        int face_parent_size = face_parent_idx[n_face];
        PDM_log_trace_array_int(face_parent_idx, n_face+1        , "face_parent_idx ::");
        PDM_log_trace_array_int(face_parent    , face_parent_size, "face_parent     ::");
      }
    }


    // > Allocate edge
    iso_edge_vtx         [i_part] = NULL;
    iso_edge_bnd_tag_idx [i_part] = NULL;
    iso_edge_bnd_tag     [i_part] = NULL;
    iso_edge_parent_idx  [i_part] = NULL;
    iso_edge_parent      [i_part] = NULL;
    iso_edge_parent_gnum [i_part] = NULL;
    PDM_malloc(isovalue_edge_idx[i_part], n_isovalues + 1, int);
    isovalue_edge_idx    [i_part][0] = 0;

    // > Allocate face
    if (isos->entry_mesh_dim==3) {
      iso_face_vtx_idx    [i_part] = NULL;
      iso_face_vtx        [i_part] = NULL;
      iso_face_parent_idx [i_part] = NULL;
      iso_face_parent     [i_part] = NULL;
      iso_face_parent_gnum[i_part] = NULL;
      PDM_malloc(isovalue_face_idx[i_part], n_isovalues + 1, int);
      isovalue_face_idx   [i_part][0] = 0;
    }

    /*
     * Build isosurface mesh
     */
    int iso_n_edge_bnd_tag = 0;
    int iso_n_edge_parent  = 0;
    int iso_n_face_parent  = 0;
    for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
      
      int *edge_to_iso_vtx = PDM_array_zeros_int(n_edge);
      int i_beg_iso_vtx = isovalue_vtx_idx[i_part][i_isovalue  ];
      int i_end_iso_vtx = isovalue_vtx_idx[i_part][i_isovalue+1];

      for (int i_iso_vtx=i_beg_iso_vtx; i_iso_vtx<i_end_iso_vtx; ++i_iso_vtx) {
        int cross_edge = iso_vtx_to_edge[i_iso_vtx];
        if (cross_edge!=0) {
          edge_to_iso_vtx[iso_vtx_to_edge[i_iso_vtx]-1] = i_iso_vtx+1;
        }
      }
      int *_vtx_to_iso_vtx = NULL;
      PDM_malloc(_vtx_to_iso_vtx, n_vtx, int);
      memcpy(_vtx_to_iso_vtx, vtx_to_iso_vtx, n_vtx * sizeof(int));
      for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
        if (_vtx_to_iso_vtx[i_vtx] <= i_beg_iso_vtx || i_end_iso_vtx < _vtx_to_iso_vtx[i_vtx]) {
          _vtx_to_iso_vtx[i_vtx] = 0;
        }
      }
      if (debug) {
        log_trace("\n");
        log_trace("Building i_isovalue %d\n", i_isovalue);
        PDM_log_trace_array_int (edge_to_iso_vtx, n_edge, "edge_to_iso_vtx ::");
        PDM_log_trace_array_int ( vtx_to_iso_vtx, n_vtx , " vtx_to_iso_vtx ::");
        PDM_log_trace_array_int (_vtx_to_iso_vtx, n_vtx , "_vtx_to_iso_vtx ::");
      }

      t_start = PDM_MPI_Wtime();
      _contouring_triangles(isos->extract_n_tri   [i_part],
                            isos->extract_tri_vtx [i_part],
                            isos->extract_tri_gnum[i_part],
                            isos->extract_tri_tag [i_part],
                            n_edge,
                            edge_bnd_tag_idx,
                            edge_bnd_tag,
                            edge_parent_idx,
                            edge_parent,
                            elt_edge[0],
                            vtx_gnum,
                            vtx_field[i_part],
                            _vtx_to_iso_vtx,
                            edge_to_iso_vtx,
                            &iso_n_edge          [i_part],
                            &iso_edge_vtx        [i_part],
                            &iso_edge_parent_gnum[i_part],
                            &iso_n_edge_bnd_tag,
                            &iso_edge_bnd_tag_idx[i_part],
                            &iso_edge_bnd_tag    [i_part],
                            &iso_n_edge_parent,
                            &iso_edge_parent_idx [i_part],
                            &iso_edge_parent     [i_part]);
      if (isos->entry_mesh_dim==3) {
        _contouring_tetrahedra(isos->extract_n_tetra   [i_part],
                               isos->extract_tetra_vtx [i_part],
                               isos->extract_tetra_gnum[i_part],
                               n_face,
                               face_vtx_idx,
                               face_vtx,
                               face_parent_idx,
                               face_parent,
                               elt_edge[1],
                               elt_face,
                               vtx_gnum,
                               vtx_field[i_part],
                               _vtx_to_iso_vtx,
                               edge_to_iso_vtx,
                               &iso_n_face          [i_part],
                               &iso_face_vtx_idx    [i_part],
                               &iso_face_vtx        [i_part],
                               &iso_face_parent_gnum[i_part],
                               &iso_n_face_parent,
                               &iso_face_parent_idx [i_part],
                               &iso_face_parent     [i_part]);
      }
      
      PDM_free(edge_to_iso_vtx);
      PDM_free(_vtx_to_iso_vtx);

      isovalue_edge_idx[i_part][i_isovalue+1] = iso_n_edge[i_part];
      if (isos->entry_mesh_dim==3) {
        isovalue_face_idx[i_part][i_isovalue+1] = iso_n_face[i_part];
      }
      
      t_end = PDM_MPI_Wtime();
      if (isos->entry_mesh_dim==2) {
        printf("   isovalue %2d : %.3fs  (%6d vertices, %6d edges)\n",
               i_isovalue,
               t_end - t_start,
               iso_n_vtx[i_part],
               isovalue_edge_idx[i_part][i_isovalue+1] - isovalue_edge_idx[i_part][i_isovalue]);
        fflush(stdout);
      }
      else if (isos->entry_mesh_dim==3) {
        printf("   isovalue %2d : %.3fs  (%6d vertices, %6d edges, %6d faces)\n",
               i_isovalue,
               t_end - t_start,
               iso_n_vtx[i_part],
               isovalue_edge_idx[i_part][i_isovalue+1] - isovalue_edge_idx[i_part][i_isovalue],
               isovalue_face_idx[i_part][i_isovalue+1] - isovalue_face_idx[i_part][i_isovalue]);
        fflush(stdout);
      }

    } // End of loop on isovalues



    /*
     * Convert bnd tag to group
     * Warning edges can be in multiple group
     */

    if (debug) {
      PDM_log_trace_array_int(iso_edge_bnd_tag_idx[i_part], iso_n_edge[i_part], "iso_edge_bnd_tag_idx ::");
      PDM_log_trace_array_int(iso_edge_bnd_tag    [i_part], iso_edge_bnd_tag_idx[i_part][iso_n_edge[i_part]], "iso_edge_bnd_tag ::");
    }
    iso_n_edge_group = isos->extract_tri_n_group[i_part];

    // > Delete duplicate in edge tags
    int *iso_edge_bnd_tag_unique_n = PDM_array_zeros_int(iso_n_edge[i_part]);
    int i_write_tag = 0;
    for (int i_edge=0; i_edge<iso_n_edge[i_part]; ++i_edge) {
      int i_beg_bnd_tag = iso_edge_bnd_tag_idx[i_part][i_edge  ];
      int i_end_bnd_tag = iso_edge_bnd_tag_idx[i_part][i_edge+1];
      int n_unique_tag = PDM_inplace_unique(iso_edge_bnd_tag[i_part], i_beg_bnd_tag, i_end_bnd_tag-1);
      for (int i_tag=0; i_tag<n_unique_tag; ++i_tag) {
        iso_edge_bnd_tag[i_part][i_write_tag++] = iso_edge_bnd_tag[i_part][i_beg_bnd_tag+i_tag];
        iso_edge_bnd_tag_unique_n[i_edge]++;
      }
    }
    PDM_free(iso_edge_bnd_tag_idx[i_part]);
    iso_edge_bnd_tag_idx[i_part] = PDM_array_new_idx_from_sizes_int(iso_edge_bnd_tag_unique_n, iso_n_edge[i_part]);
    PDM_realloc(iso_edge_bnd_tag[i_part], iso_edge_bnd_tag[i_part], iso_edge_bnd_tag_idx[i_part][iso_n_edge[i_part]], int);
    if (debug) {
      PDM_log_trace_array_int(iso_edge_bnd_tag_unique_n   , iso_n_edge[i_part], "iso_edge_bnd_tag_unique_n ::");
      PDM_log_trace_array_int(iso_edge_bnd_tag_idx[i_part], iso_n_edge[i_part], "iso_edge_bnd_tag_idx ::");
      PDM_log_trace_array_int(iso_edge_bnd_tag    [i_part], iso_edge_bnd_tag_idx[i_part][iso_n_edge[i_part]], "iso_edge_bnd_tag ::");
    }
    PDM_free(iso_edge_bnd_tag_unique_n);


    // > Count number of entities in each group
    int *iso_edge_group_n = PDM_array_zeros_int(iso_n_edge_group);
    for (int i_tag=0; i_tag<iso_edge_bnd_tag_idx[i_part][iso_n_edge[i_part]]; ++i_tag) {
      int i_group = iso_edge_bnd_tag[i_part][i_tag]-1;
      iso_edge_group_n[i_group]++;
    }

    // > Fill groups with entities
    iso_edge_group_idx [i_part] = PDM_array_new_idx_from_sizes_int(iso_edge_group_n, iso_n_edge_group);
    iso_edge_group_lnum[i_part] = PDM_array_zeros_int(iso_edge_group_idx[i_part][iso_n_edge_group]);
    PDM_array_reset_int(iso_edge_group_n, iso_n_edge_group, 0);
    for (int i_edge=0; i_edge<iso_n_edge[i_part]; ++i_edge) {
      int i_beg_bnd_tag = iso_edge_bnd_tag_idx[i_part][i_edge  ];
      int i_end_bnd_tag = iso_edge_bnd_tag_idx[i_part][i_edge+1];

      for (int i_tag=i_beg_bnd_tag; i_tag<i_end_bnd_tag; ++i_tag) {
        int i_group = iso_edge_bnd_tag[i_part][i_tag]-1;
        int i_write = iso_edge_group_idx[i_part][i_group]+iso_edge_group_n[i_group];

        iso_edge_group_lnum[i_part][i_write] = i_edge+1;
        iso_edge_group_n[i_group]++;
      }
    }
    PDM_free(iso_edge_group_n);
    PDM_free(iso_edge_bnd_tag_idx[i_part]);
    PDM_free(iso_edge_bnd_tag    [i_part]);

    if (debug) {
      PDM_log_trace_array_int(iso_edge_group_idx [i_part], iso_n_edge_group+1, "iso_edge_group_idx ::");
      PDM_log_trace_array_int(iso_edge_group_lnum[i_part], iso_edge_group_idx[i_part][iso_n_edge_group], "iso_edge_group_lnum ::");
      for (int i_group=0; i_group<iso_n_edge_group; ++i_group) {
        log_trace("i_group = %d\n", i_group+1);
        int i_beg      = iso_edge_group_idx[i_part][i_group];
        int n_in_group = iso_edge_group_idx[i_part][i_group+1]-iso_edge_group_idx[i_part][i_group];
        PDM_log_trace_array_int(&iso_edge_group_lnum[i_part][i_beg], n_in_group, "iso_edge_group_lnum ::");
      }
    }




    /* 
     * Visu
     */
    if (debug_visu) {

      log_trace("\n");
      log_trace("iso_n_edge = %d\n", iso_n_edge[i_part]);
      PDM_log_trace_array_int(iso_edge_vtx[i_part], iso_n_edge[i_part]*2, "iso_edge_vtx     ::");
      char outname[999];

      log_trace("\n");
      log_trace("iso_n_vertices = %d \n", iso_n_vtx[i_part]);

      int *iso_edge_isovalue = NULL;
      PDM_malloc(iso_edge_isovalue, iso_n_edge[i_part], int);
      for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
        for (int i = isovalue_edge_idx[i_part][i_isovalue]; i < isovalue_edge_idx[i_part][i_isovalue+1]; i++) {
          iso_edge_isovalue[i] = i_isovalue;
        }
      }

      const char *field_name [] = {"i_isovalue"};
      const int  *field_value[] = {iso_edge_isovalue};

      sprintf(outname, "marching_sections_iso_edges_%d.vtk", i_part);
      PDM_vtk_write_std_elements(outname,
                                 iso_n_vtx[i_part],
                                 iso_vtx_coord[i_part],
                                 NULL,
                                 PDM_MESH_NODAL_BAR2,
                                 iso_n_edge[i_part],
                                 iso_edge_vtx[i_part],
                                 NULL,
                                 1,
                                 field_name,
                                 field_value);
      PDM_free(iso_edge_isovalue);
    } // End visu


    // > Free temp arrays
    for (int i_section = 0; i_section < n_section; i_section++) {
      PDM_free(elt_edge[i_section]);
    }
    if (isos->entry_mesh_dim==3) {
      PDM_free(elt_face);
    }

    PDM_free(elt_edge);
    PDM_free(edge_vtx);
    PDM_free(edge_bnd_tag_idx);
    PDM_free(edge_bnd_tag);
    PDM_free(edge_parent_idx);
    PDM_free(edge_parent);
    PDM_free(iso_edge_bnd_tag_idx);
    PDM_free(iso_edge_bnd_tag);
    PDM_free(face_vtx_idx);
    PDM_free(face_vtx);
    PDM_free(face_parent_idx);
    PDM_free(face_parent);
    PDM_free(vtx_to_iso_vtx);
    PDM_free(iso_vtx_to_edge);
  } // End of loop on partitions


  /*
   * Build vertices/edges/faces gnum
   */
  PDM_g_num_t **iso_vtx_gnum  = NULL;
  PDM_malloc(iso_vtx_gnum , n_part, PDM_g_num_t *);
  _generate_gnum_from_parents(isos->comm,
                              isos->n_part,
                              iso_n_vtx,
                              iso_vtx_parent_gnum,
                              3,
                              iso_vtx_gnum);

  PDM_g_num_t **iso_edge_gnum = NULL;
  PDM_malloc(iso_edge_gnum, n_part, PDM_g_num_t *);
  _generate_gnum_from_parents(isos->comm,
                              isos->n_part,
                              iso_n_edge,
                              iso_edge_parent_gnum,
                              2,
                              iso_edge_gnum);

  PDM_g_num_t **iso_face_gnum = NULL;
  if (isos->entry_mesh_dim==3) {
    PDM_malloc(iso_face_gnum, n_part, PDM_g_num_t *);
    _generate_gnum_from_parents(isos->comm,
                                isos->n_part,
                                iso_n_face,
                                iso_face_parent_gnum,
                                3,
                                iso_face_gnum);
  }

  /**
   * Now we have generated iso gnums, we can use these array to store parent_gnum
   */
  if (isos->entry_is_part == 0) {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      // > Vertices
      int i_write = 0;
      for (int i_vtx = 0; i_vtx < iso_n_vtx[i_part]; i_vtx++) {
        iso_vtx_parent_gnum  [i_part][i_write++] = iso_vtx_parent_gnum[i_part][3*i_vtx  ];
        if (iso_vtx_parent_gnum[i_part][3*i_vtx+1]>0) {
          iso_vtx_parent_gnum[i_part][i_write++] = iso_vtx_parent_gnum[i_part][3*i_vtx+1];
        }
      }
      PDM_realloc(iso_vtx_parent_gnum[i_part], iso_vtx_parent_gnum[i_part], i_write, PDM_g_num_t);

      // > Edges
      i_write = 0;
      for (int i_edge = 0; i_edge < iso_n_edge[i_part]; i_edge++) {
        int i_beg_parent = iso_edge_parent_idx[i_part][i_edge  ];
        int i_end_parent = iso_edge_parent_idx[i_part][i_edge+1];
        for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
          int tri_lnum = iso_edge_parent[i_part][i_parent];
          iso_edge_parent_gnum[i_part][i_write++] = isos->extract_tri_gnum[i_part][tri_lnum-1];
        }
      }
      PDM_realloc(iso_edge_parent_gnum[i_part], iso_edge_parent_gnum[i_part], i_write, PDM_g_num_t);
      

      // > Faces
      if (isos->entry_mesh_dim==3) {
        i_write = 0;
        for (int i_face = 0; i_face < iso_n_face[i_part]; i_face++) {
          int i_beg_parent = iso_face_parent_idx[i_part][i_face  ];
          int i_end_parent = iso_face_parent_idx[i_part][i_face+1];
          for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
            int tetra_lnum = iso_face_parent[i_part][i_parent];
            iso_face_parent_gnum[i_part][i_write++] = isos->extract_tetra_gnum[i_part][tetra_lnum-1];
          }
        }
        PDM_realloc(iso_face_parent_gnum[i_part], iso_face_parent_gnum[i_part], i_write, PDM_g_num_t);
      }
    }
  }
  else {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(iso_vtx_parent_gnum [i_part]);
      PDM_free(iso_edge_parent_gnum[i_part]);
      if (isos->entry_mesh_dim==3) {
        PDM_free(iso_face_parent_gnum[i_part]);
      }
    }
    PDM_free(iso_vtx_parent_gnum );
    PDM_free(iso_edge_parent_gnum);
    PDM_free(iso_face_parent_gnum);
  }



  /*
   * Build edge groups gnum
   */
  PDM_g_num_t **iso_edge_group_gnum = NULL;
  PDM_malloc(iso_edge_group_gnum, n_part, PDM_g_num_t *);
  for (int i_part=0; i_part<n_part; i_part++) {
    PDM_malloc(iso_edge_group_gnum[i_part], iso_edge_group_idx[i_part][iso_n_edge_group], PDM_g_num_t);
  }
  PDM_gen_gnum_t *gen_gnum_edge_group = PDM_gnum_create(3, n_part, PDM_FALSE, 1., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_parents_nuplet(gen_gnum_edge_group, 1);

  PDM_g_num_t   **_iso_edge_group_gnum = NULL;
  PDM_malloc(_iso_edge_group_gnum, n_part, PDM_g_num_t *);
  for (int i_group=0; i_group<iso_n_edge_group; ++i_group) {
    for (int i_part=0; i_part<n_part; i_part++) {
      int i_beg_group = iso_edge_group_idx[i_part][i_group  ];
      int i_end_group = iso_edge_group_idx[i_part][i_group+1];
      int n_edge_in_group = i_end_group - i_beg_group;
      PDM_malloc(_iso_edge_group_gnum[i_part], n_edge_in_group, PDM_g_num_t);
      int i_write = 0;
      for (int i_edge=i_beg_group; i_edge<i_end_group; ++i_edge) {
        int lnum = iso_edge_group_lnum[i_part][i_edge];
        _iso_edge_group_gnum[i_part][i_write++] = iso_edge_gnum[i_part][lnum-1];
      }
      PDM_gnum_set_from_parents(gen_gnum_edge_group, i_part, n_edge_in_group, _iso_edge_group_gnum[i_part]);
    }
    PDM_gnum_compute(gen_gnum_edge_group);

    for (int i_part=0; i_part<n_part; i_part++) {
      PDM_free(_iso_edge_group_gnum[i_part]);
      int i_beg_group = iso_edge_group_idx[i_part][i_group  ];
      int i_end_group = iso_edge_group_idx[i_part][i_group+1];
      _iso_edge_group_gnum[i_part] = PDM_gnum_get(gen_gnum_edge_group, i_part);
      int i_read = 0;
      for (int i_edge=i_beg_group; i_edge<i_end_group; ++i_edge) {
        iso_edge_group_gnum[i_part][i_edge] = _iso_edge_group_gnum[i_part][i_read++];
      }
      PDM_free(_iso_edge_group_gnum[i_part]);
    }
  }
  PDM_gnum_free(gen_gnum_edge_group);
  PDM_free(_iso_edge_group_gnum);

  /*
   * Store isosurface in part_mesh_nodal
   */
  // isos->n_part = n_part; // isos->iso_n_part

  // Vertices
  _iso->iso_n_entity         [PDM_MESH_ENTITY_VTX] = iso_n_vtx;
  _iso->iso_entity_gnum      [PDM_MESH_ENTITY_VTX] = iso_vtx_gnum;
  _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_VTX] = iso_vtx_parent_idx;
  _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_VTX] = isovalue_vtx_idx;
  if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
    _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_VTX] = iso_vtx_parent_gnum;
  }
  _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_VTX] = iso_vtx_parent;
  _iso->iso_vtx_coord         = iso_vtx_coord;
  _iso->iso_vtx_parent_weight = iso_vtx_parent_weight;

  // Edges
  _iso->iso_n_entity         [PDM_MESH_ENTITY_EDGE] = iso_n_edge;
  _iso->iso_entity_gnum      [PDM_MESH_ENTITY_EDGE] = iso_edge_gnum;
  _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_EDGE] = iso_edge_parent_idx;
  _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_EDGE] = isovalue_edge_idx;
  if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
    _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_EDGE] = iso_edge_parent_gnum;
  }
  _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_EDGE] = iso_edge_parent;
  _iso->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = iso_edge_vtx;
  _iso->iso_n_edge_group    = iso_n_edge_group;
  _iso->iso_edge_group_idx  = iso_edge_group_idx;
  _iso->iso_edge_group_lnum = iso_edge_group_lnum;
  _iso->iso_edge_group_gnum = iso_edge_group_gnum;

  // Faces
  if (isos->entry_mesh_dim==3) {
    _iso->iso_n_entity         [PDM_MESH_ENTITY_FACE] = iso_n_face;
    _iso->iso_entity_gnum      [PDM_MESH_ENTITY_FACE] = iso_face_gnum;
    _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_FACE] = iso_face_parent_idx;
    _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_FACE] = isovalue_face_idx;
    if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
      _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_FACE] = iso_face_parent_gnum;
    }
    _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_FACE] = iso_face_parent;
    _iso->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX] = iso_face_vtx_idx;
    _iso->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX] = iso_face_vtx;
  }

  // Ownerships
  PDM_malloc(_iso->iso_owner_vtx_coord        , isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_vtx_parent_weight, isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_edge_bnd         , isos->n_part, PDM_ownership_t);
  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    PDM_malloc(_iso->iso_owner_gnum               [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_parent_lnum        [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_parent_idx         [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_isovalue_entity_idx[i_entity], isos->n_part, PDM_ownership_t);
  }
  PDM_malloc(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_part, PDM_ownership_t);

  for (int i_part=0; i_part<isos->n_part; ++i_part) {
    _iso->iso_owner_vtx_coord        [i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_vtx_parent_weight[i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_edge_bnd         [i_part] = PDM_OWNERSHIP_KEEP;
    for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
      _iso->iso_owner_gnum               [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_parent_lnum        [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_parent_idx         [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_isovalue_entity_idx[i_entity][i_part] = PDM_OWNERSHIP_KEEP;
    }

    _iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part] = PDM_OWNERSHIP_KEEP;
  }

}



void
PDM_isosurface_ngon_algo
(
  PDM_isosurface_t        *isos,
  int                      id_iso
)
{
  /**
   * TODO (à discuter avec mon rayon de soleil):
   *   - isoler et factoriser tronc commun 2d/3d
   *   - gérer 2d:
   *     - trianguler faces (besoin de face_vtx, garder les edges ou juste les actives?)
   *     - factoriser avec contouring_triangles
   */
  // int debug      = 1;
  int debug_visu = 0;
  double t_start, t_end;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  PDM_MPI_Comm_size(isos->comm, &n_rank);

  // > Isosurface information
  _isosurface_t *_iso = &isos->isosurfaces[id_iso];
  int   n_isovalues = _iso->n_isovalues;
  double *isovalues = _iso->  isovalues;

  double **vtx_field = isos->extract_field;

  // > Allocate isosurface vertices
  int          *iso_n_vtx             = NULL;
  double      **iso_vtx_coord         = NULL;
  PDM_g_num_t **iso_vtx_gnum          = NULL;
  int         **iso_vtx_parent_idx    = NULL;
  int         **iso_vtx_parent        = NULL;
  double      **iso_vtx_parent_weight = NULL;
  PDM_g_num_t **iso_vtx_parent_gnum   = NULL;
  int         **isovalue_vtx_idx      = NULL;
  PDM_malloc(iso_n_vtx            , isos->n_part, int          );
  PDM_malloc(iso_vtx_coord        , isos->n_part, double      *);
  PDM_malloc(iso_vtx_gnum         , isos->n_part, PDM_g_num_t *);
  PDM_malloc(iso_vtx_parent_idx   , isos->n_part, int         *);
  PDM_malloc(iso_vtx_parent       , isos->n_part, int         *);
  PDM_malloc(iso_vtx_parent_weight, isos->n_part, double      *);
  PDM_malloc(iso_vtx_parent_gnum  , isos->n_part, PDM_g_num_t *);
  PDM_malloc(isovalue_vtx_idx     , isos->n_part, int         *);


  // > Allocate isosurface edges
  int          *iso_n_edge           = NULL;
  int         **iso_edge_vtx         = NULL;
  PDM_g_num_t **iso_edge_gnum        = NULL;
  int         **iso_edge_parent_idx  = NULL;
  int         **iso_edge_parent      = NULL;
  PDM_g_num_t **iso_edge_parent_gnum = NULL;
  int         **isovalue_edge_idx    = NULL;
  int         **iso_edge_group_idx   = NULL;
  int         **iso_edge_group_lnum  = NULL;
  PDM_g_num_t **iso_edge_group_gnum   = NULL;
  PDM_malloc(iso_n_edge          , isos->n_part, int          );
  PDM_malloc(iso_edge_vtx        , isos->n_part, int         *);
  PDM_malloc(iso_edge_gnum       , isos->n_part, PDM_g_num_t *);
  PDM_malloc(iso_edge_parent_idx , isos->n_part, int         *);
  PDM_malloc(iso_edge_parent     , isos->n_part, int         *);
  PDM_malloc(iso_edge_parent_gnum, isos->n_part, PDM_g_num_t *);
  PDM_malloc(isovalue_edge_idx   , isos->n_part, int         *);
  PDM_malloc(iso_edge_group_idx  , isos->n_part, int         *);
  PDM_malloc(iso_edge_group_lnum , isos->n_part, int         *);
  PDM_malloc(iso_edge_group_gnum , isos->n_part, PDM_g_num_t *);


  // > Allocate face edges
  int          *iso_n_face           = NULL;
  int         **iso_face_vtx_idx     = NULL;
  int         **iso_face_vtx         = NULL;
  PDM_g_num_t **iso_face_gnum        = NULL;
  int         **iso_face_parent_idx  = NULL;
  int         **iso_face_parent      = NULL;
  PDM_g_num_t **iso_face_parent_gnum = NULL;
  int         **isovalue_face_idx    = NULL;
  if (isos->entry_mesh_dim == 3) {
    PDM_malloc(iso_n_face          , isos->n_part, int          );
    PDM_malloc(iso_face_vtx_idx    , isos->n_part, int         *);
    PDM_malloc(iso_face_vtx        , isos->n_part, int         *);
    PDM_malloc(iso_face_gnum       , isos->n_part, PDM_g_num_t *);
    PDM_malloc(iso_face_parent_idx , isos->n_part, int         *);
    PDM_malloc(iso_face_parent     , isos->n_part, int         *);
    PDM_malloc(iso_face_parent_gnum, isos->n_part, PDM_g_num_t *);
    PDM_malloc(isovalue_face_idx   , isos->n_part, int         *);
  }

  PDM_part_mesh_t *pmesh = isos->extract_pmesh;
  assert(pmesh != NULL);

  int *iso_edge_group_n = NULL;
  int  n_surface        = 0;
  if (isos->entry_mesh_dim==3) {
    n_surface = isos->n_group_face;
    PDM_malloc(iso_edge_group_n, n_surface, int);
  }

  t_start = PDM_MPI_Wtime();

  /* Generate edges if necessary */
  int          *pn_face        = NULL;
  int         **pface_vtx_idx  = NULL;
  int         **pface_vtx      = NULL;
  int          *pn_vtx         = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  int          *pn_edge        = NULL;
  int         **pedge_vtx      = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;

  if (!isos->we_have_edges) {
    PDM_malloc(pn_face       , isos->n_part, int          );
    PDM_malloc(pn_vtx        , isos->n_part, int          );
    PDM_malloc(pface_vtx_idx , isos->n_part, int         *);
    PDM_malloc(pface_vtx     , isos->n_part, int         *);
    PDM_malloc(pface_ln_to_gn, isos->n_part, PDM_g_num_t *);
    PDM_malloc(pvtx_ln_to_gn , isos->n_part, PDM_g_num_t *);

    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      pn_face[i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE);
      pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pmesh,
                                                   i_part,
                                                   PDM_MESH_ENTITY_VTX);
      PDM_part_mesh_connectivity_get(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     &pface_vtx    [i_part],
                                     &pface_vtx_idx[i_part],
                                     PDM_OWNERSHIP_BAD_VALUE);

      PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &pface_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_BAD_VALUE);

      PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &pvtx_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_BAD_VALUE);
    }

    PDM_compute_face_edge_from_face_vtx(isos->comm,
                                        isos->n_part,
                                        pn_face,
                                        pn_vtx,
                                        pface_vtx_idx,
                                        pface_vtx,
                                        pface_ln_to_gn,
                                        pvtx_ln_to_gn,
                                        &pface_edge_idx,
                                        &pface_edge,
                                        &pn_edge,
                                        &pedge_vtx,
                                        &pedge_ln_to_gn);
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(pedge_ln_to_gn[i_part]);
    }
    PDM_free(pn_face       );
    PDM_free(pn_vtx        );
    PDM_free(pface_vtx_idx );
    PDM_free(pface_vtx     );
    PDM_free(pface_ln_to_gn);
    PDM_free(pvtx_ln_to_gn );
    PDM_free(pedge_ln_to_gn);
  }

  /* Build isosurface */
  for (int i_part = 0; i_part < isos->n_part; i_part++) {

    int n_cell = 0;
    int n_face = 0;
    int n_edge = 0;
    int n_vtx  = 0;

    int    *cell_face_idx = NULL;
    int    *cell_face     = NULL;
    int    *face_edge_idx = NULL;
    int    *face_edge     = NULL;
    int    *edge_vtx      = NULL;
    double *vtx_coord     = NULL;

    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    n_cell = PDM_part_mesh_n_entity_get(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL);
    n_face = PDM_part_mesh_n_entity_get(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE);
    n_vtx  = PDM_part_mesh_n_entity_get(pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_connectivity_get(pmesh,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   &cell_face,
                                   &cell_face_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    if (isos->we_have_edges) {
      n_edge = PDM_part_mesh_n_entity_get(pmesh,
                                          i_part,
                                          PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_connectivity_get(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                     &face_edge,
                                     &face_edge_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);

      int *edge_vtx_idx = NULL;
      PDM_part_mesh_connectivity_get(pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                     &edge_vtx,
                                     &edge_vtx_idx,
                                     PDM_OWNERSHIP_BAD_VALUE);
    }
    else {
      n_edge        = pn_edge       [i_part];
      face_edge_idx = pface_edge_idx[i_part];
      face_edge     = pface_edge    [i_part];
      edge_vtx      = pedge_vtx     [i_part];
    }

    PDM_part_mesh_vtx_coord_get(pmesh,
                                i_part,
                                &vtx_coord,
                                PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    // group_face -> face_tag
    int *face_tag = NULL;
    if (n_surface > 0) {
      face_tag = PDM_array_zeros_int(n_face);
      for (int i_surface = 0; i_surface < n_surface; i_surface++) {

        int          surface_face_n    = 0;
        int         *surface_face_lnum = NULL;
        PDM_g_num_t *surface_face_gnum = NULL;
        PDM_part_mesh_bound_get(pmesh,
                                i_part,
                                i_surface,
                                PDM_BOUND_TYPE_FACE,
                                &surface_face_n,
                                &surface_face_lnum,
                                &surface_face_gnum,
                                PDM_OWNERSHIP_BAD_VALUE);

        for (int i = 0; i < surface_face_n; i++) {
          face_tag[surface_face_lnum[i]-1] = i_surface + 1;
        }
      }
    }


    /* Build isosurface */
    iso_n_vtx            [i_part] = 0;
    iso_vtx_coord        [i_part] = NULL;
    iso_vtx_parent_idx   [i_part] = NULL;
    iso_vtx_parent       [i_part] = NULL;
    iso_vtx_parent_weight[i_part] = NULL;
    isovalue_vtx_idx     [i_part] = NULL;
    iso_n_edge           [i_part] = 0;
    iso_edge_vtx         [i_part] = NULL;
    iso_edge_parent_idx  [i_part] = NULL;
    iso_edge_parent      [i_part] = NULL;
    isovalue_edge_idx    [i_part] = NULL;
    if (isos->entry_mesh_dim == 3) {
      iso_n_face         [i_part] = 0;
      iso_face_vtx_idx   [i_part] = NULL;
      iso_face_vtx       [i_part] = NULL;
      iso_face_parent    [i_part] = NULL;
      iso_face_parent_idx[i_part] = NULL;
      isovalue_face_idx  [i_part] = NULL;
    }
    int *iso_vtx_parent_edge = NULL;

    if (debug_visu) {
      int *face_vtx = NULL;
      if (isos->we_have_edges) {
        PDM_compute_face_vtx_from_face_and_edge(n_face,
                                                face_edge_idx,
                                                face_edge,
                                                edge_vtx,
                                                &face_vtx);
      }
      else {
        PDM_part_mesh_connectivity_get(pmesh,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       &face_vtx,
                                       &face_edge_idx,
                                       PDM_OWNERSHIP_BAD_VALUE);
      }

      char name[999];
      sprintf(name, "extract_faces_%d.vtk", id_iso);
      PDM_vtk_write_polydata_field(name,
                                   n_vtx,
                                   vtx_coord,
                                   NULL,
                                   n_face,
                                   face_edge_idx,
                                   face_vtx,
                                   NULL,
                                   NULL,
                                   NULL,
                                   "field",
                                   // isos->extract_field[id_iso][i_part]);
                                   isos->extract_field[i_part]);
      if (isos->we_have_edges) {
        PDM_free(face_vtx);
      }
    }

    _isosurface_ngon_single_part(isos->entry_mesh_dim,
                                 n_isovalues,
                                 isovalues,
                                 isos->ISOSURFACE_EPS,
                                 n_cell,
                                 n_face,
                                 n_edge,
                                 n_vtx,
                                 cell_face_idx,
                                 cell_face,
                                 face_edge_idx,
                                 face_edge,
                                 edge_vtx,
                                 vtx_coord,
                                 vtx_field[i_part],
                                 face_tag,
                                 cell_ln_to_gn,
                                 &iso_n_vtx            [i_part],
                                 &iso_vtx_coord        [i_part],
                                 &iso_vtx_parent_idx   [i_part],
                                 &iso_vtx_parent       [i_part],
                                 &iso_vtx_parent_weight[i_part],
                                 &iso_vtx_parent_edge,
                                 &isovalue_vtx_idx     [i_part],
                                 &iso_n_face           [i_part],
                                 &iso_face_vtx_idx     [i_part],
                                 &iso_face_vtx         [i_part],
                                 &iso_face_parent_idx  [i_part],
                                 &iso_face_parent      [i_part],
                                 &isovalue_face_idx    [i_part],
                                 &iso_n_edge           [i_part],
                                 &iso_edge_vtx         [i_part],
                                 &iso_edge_parent_idx  [i_part],
                                 &iso_edge_parent      [i_part],
                                 &isovalue_edge_idx    [i_part]);

    // Setup parent gnums
    int *edge_n_child = NULL;
    PDM_calloc(edge_n_child, n_edge, int);
    PDM_malloc(iso_vtx_parent_gnum[i_part], iso_n_vtx[i_part] * 3, PDM_g_num_t);
    for (int i = 0; i < iso_n_vtx[i_part]; i++) {
      iso_vtx_parent_gnum  [i_part][3*i  ] =  vtx_ln_to_gn[iso_vtx_parent[i_part][iso_vtx_parent_idx[i_part][i]] - 1];
      if (iso_vtx_parent_idx[i_part][i+1] - iso_vtx_parent_idx[i_part][i] > 1) {
        iso_vtx_parent_gnum[i_part][3*i+1] =  vtx_ln_to_gn[iso_vtx_parent[i_part][iso_vtx_parent_idx[i_part][i]+1] - 1];
        iso_vtx_parent_gnum[i_part][3*i+2] =--edge_n_child[iso_vtx_parent_edge[i]];
      }
      else {
        iso_vtx_parent_gnum[i_part][3*i+1] = 0;
        iso_vtx_parent_gnum[i_part][3*i+2] = 0;
      }
    }
    free(edge_n_child);
    free(iso_vtx_parent_edge);


    // TODO: handle case n_parent > 1
    int *face_n_child = NULL;
    PDM_calloc(face_n_child, n_face, int);
    PDM_malloc(iso_edge_parent_gnum[i_part], 2*iso_n_edge[i_part], PDM_g_num_t);
    for (int i = 0; i < iso_n_edge[i_part]; i++) {
      int face_lnum = iso_edge_parent[i_part][iso_edge_parent_idx[i_part][i]];
      iso_edge_parent_gnum[i_part][2*i  ] =  face_ln_to_gn[face_lnum-1];
      iso_edge_parent_gnum[i_part][2*i+1] =--face_n_child [face_lnum-1];
    }
    free(face_n_child);


    // TODO: handle case n_parent > 1
    if (isos->entry_mesh_dim == 3) {
      int *cell_n_child = NULL;
      PDM_calloc(cell_n_child, n_cell, int);
      PDM_malloc(iso_face_parent_gnum[i_part], 2*iso_n_face[i_part], PDM_g_num_t);
      for (int i = 0; i < iso_n_face[i_part]; i++) {
        int cell_lnum = iso_face_parent[i_part][iso_face_parent_idx[i_part][i]];
        iso_face_parent_gnum[i_part][2*i  ] =  cell_ln_to_gn[cell_lnum-1];
        iso_face_parent_gnum[i_part][2*i+1] =--cell_n_child [cell_lnum-1];
      }
      free(cell_n_child);
    }


    // Groups
    if (n_surface > 0) {
      PDM_array_reset_int(iso_edge_group_n, n_surface, 0);
      for (int i_edge = 0; i_edge < iso_n_edge[i_part]; i_edge++) {
        int i_face = iso_edge_parent[i_part][iso_edge_parent_idx[i_part][i_edge]] - 1;
        assert(face_tag[i_face] > 0);
        iso_edge_group_n[face_tag[i_face]-1]++;
      }

      iso_edge_group_idx [i_part] = PDM_array_new_idx_from_sizes_int(iso_edge_group_n, n_surface);
      PDM_malloc(iso_edge_group_lnum[i_part], iso_edge_group_idx[i_part][n_surface], int);
      PDM_array_reset_int(iso_edge_group_n, n_surface, 0);
      for (int i_edge = 0; i_edge < iso_n_edge[i_part]; i_edge++) {
        int i_face    = iso_edge_parent[i_part][iso_edge_parent_idx[i_part][i_edge]] - 1;
        int i_surface = face_tag[i_face]-1;
        // TODO: handle case with face on multiple surfaces?
        iso_edge_group_lnum[i_part][iso_edge_group_idx[i_part][i_surface] + iso_edge_group_n[i_surface]] = i_edge+1;
        iso_edge_group_n[i_surface]++;
      }
      PDM_free(face_tag);
    }

    if (!isos->we_have_edges) {
      PDM_free(pface_edge_idx[i_part]);
      PDM_free(pface_edge    [i_part]);
      PDM_free(pedge_vtx     [i_part]);
    }
  }
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge    );
  PDM_free(pn_edge       );
  PDM_free(pedge_vtx     );
  PDM_free(iso_edge_group_n);

  t_end = PDM_MPI_Wtime();
  log_trace("Contouring ngon : %.3fs \n", t_end - t_start);

  /* Generate global IDs */
  PDM_MPI_Barrier(isos->comm);
  t_start = PDM_MPI_Wtime();


  // TODO: use "fast gnums"?
  _generate_gnum_from_parents(isos->comm,
                              isos->n_part,
                              iso_n_vtx,
                              iso_vtx_parent_gnum,
                              3,
                              iso_vtx_gnum);
  _generate_gnum_from_parents(isos->comm,
                              isos->n_part,
                              iso_n_edge,
                              iso_edge_parent_gnum,
                              2,
                              iso_edge_gnum);
  if (isos->entry_mesh_dim == 3) {
    _generate_gnum_from_parents(isos->comm,
                                isos->n_part,
                                iso_n_face,
                                iso_face_parent_gnum,
                                2,
                                iso_face_gnum);
  }
  


  /*
   * For vertices, parent gnum can be oversized for gnum computation, so we need to remove gaps
   */
  if (isos->entry_is_part == 0) {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      int i_write = 0;
      for (int i = 0; i < iso_n_vtx[i_part]; i++) {
        iso_vtx_parent_gnum[i_part][i_write++] = iso_vtx_parent_gnum[i_part][3*i];
        if (iso_vtx_parent_idx[i_part][i+1] - iso_vtx_parent_idx[i_part][i] > 1) {
          iso_vtx_parent_gnum[i_part][i_write++] = iso_vtx_parent_gnum[i_part][3*i+1];
        }
      }
      PDM_realloc(iso_vtx_parent_gnum[i_part], iso_vtx_parent_gnum[i_part], i_write, PDM_g_num_t);

      i_write = 0;
      for (int i = 0; i < iso_n_edge[i_part]; i++) {
        iso_edge_parent_gnum[i_part][i_write++] = iso_edge_parent_gnum[i_part][2*i];
      }
      PDM_realloc(iso_edge_parent_gnum[i_part], iso_edge_parent_gnum[i_part], i_write, PDM_g_num_t);

      if (isos->entry_mesh_dim == 3) {
        i_write = 0;
        for (int i = 0; i < iso_n_face[i_part]; i++) {
          iso_face_parent_gnum[i_part][i_write++] = iso_face_parent_gnum[i_part][2*i];
        }
        PDM_realloc(iso_face_parent_gnum[i_part], iso_face_parent_gnum[i_part], i_write, PDM_g_num_t);
      }
    }
  }
  else {
    // > We can free parent gnum now that isosurface gnum are computed
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(iso_vtx_parent_gnum [i_part]);
      PDM_free(iso_edge_parent_gnum[i_part]);
      if (isos->entry_mesh_dim == 3) {
        PDM_free(iso_face_parent_gnum[i_part]);
      }
    }
    PDM_free(iso_vtx_parent_gnum);
    PDM_free(iso_edge_parent_gnum);
    PDM_free(iso_face_parent_gnum);
  }


  if (n_surface > 0) {
    PDM_g_num_t **iso_edge_group_parent_gnum = NULL;
    PDM_malloc(iso_edge_group_parent_gnum, isos->n_part, PDM_g_num_t *);
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_malloc(iso_edge_group_parent_gnum[i_part], iso_edge_group_idx[i_part][n_surface], PDM_g_num_t);
      PDM_malloc(iso_edge_group_gnum       [i_part], iso_edge_group_idx[i_part][n_surface], PDM_g_num_t);
      for (int i = 0; i < iso_edge_group_idx[i_part][n_surface]; i++) {
        iso_edge_group_parent_gnum[i_part][i] = iso_edge_gnum[i_part][iso_edge_group_lnum[i_part][i] - 1];
      }
    }

    for (int i_surface = 0; i_surface < n_surface; i_surface++) {
      PDM_gen_gnum_t *gen_gnum = PDM_gnum_create(3,  // unused,
                                                 isos->n_part,
                                                 PDM_FALSE,
                                                 1., // unused,
                                                 isos->comm,
                                                 PDM_OWNERSHIP_KEEP);

      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        PDM_gnum_set_from_parents(gen_gnum,
                                  i_part,
                                  iso_edge_group_idx[i_part][i_surface+1] - iso_edge_group_idx[i_part][i_surface],
                                  &iso_edge_group_parent_gnum[i_part][iso_edge_group_idx[i_part][i_surface]]);
      }

      PDM_gnum_set_parents_nuplet(gen_gnum, 1);
      PDM_gnum_compute(gen_gnum);


      for (int i_part = 0; i_part < isos->n_part; i_part++) {
        PDM_g_num_t *gnum = PDM_gnum_get(gen_gnum, i_part);
        memcpy(iso_edge_group_gnum[i_part] + iso_edge_group_idx[i_part][i_surface],
               gnum,
               sizeof(PDM_g_num_t) * (iso_edge_group_idx[i_part][i_surface+1] - iso_edge_group_idx[i_part][i_surface]));
      }
      PDM_gnum_free(gen_gnum);
    }
    for (int i_part = 0; i_part < isos->n_part; i_part++) {
      PDM_free(iso_edge_group_parent_gnum[i_part]);
    }
    PDM_free(iso_edge_group_parent_gnum);
  }

  t_end = PDM_MPI_Wtime();
  log_trace("Generate global IDs : %.3fs \n", t_end - t_start);


  /*
   * Store isosurface in struct
   */
  // Vertices
  _iso->iso_n_entity         [PDM_MESH_ENTITY_VTX] = iso_n_vtx;
  _iso->iso_entity_gnum      [PDM_MESH_ENTITY_VTX] = iso_vtx_gnum;
  _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_VTX] = iso_vtx_parent_idx;
  _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_VTX] = isovalue_vtx_idx;

  if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
    _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_VTX] = iso_vtx_parent_gnum;
  }
  _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_VTX] = iso_vtx_parent;
  _iso->iso_vtx_coord         = iso_vtx_coord;
  _iso->iso_vtx_parent_weight = iso_vtx_parent_weight;

  // Edges
  _iso->iso_n_entity         [PDM_MESH_ENTITY_EDGE] = iso_n_edge;
  _iso->iso_entity_gnum      [PDM_MESH_ENTITY_EDGE] = iso_edge_gnum;
  _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_EDGE] = iso_edge_parent_idx;
  _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_EDGE] = isovalue_edge_idx;
  if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
    _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_EDGE] = iso_edge_parent_gnum;
  }
  _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_EDGE] = iso_edge_parent;
  _iso->iso_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = iso_edge_vtx;
  _iso->iso_n_edge_group    = n_surface;
  _iso->iso_edge_group_idx  = iso_edge_group_idx;
  _iso->iso_edge_group_lnum = iso_edge_group_lnum;
  _iso->iso_edge_group_gnum = iso_edge_group_gnum;

  // Faces
  if (isos->entry_mesh_dim==3) {
    _iso->iso_n_entity         [PDM_MESH_ENTITY_FACE] = iso_n_face;
    _iso->iso_entity_gnum      [PDM_MESH_ENTITY_FACE] = iso_face_gnum;
    _iso->iso_entity_parent_idx[PDM_MESH_ENTITY_FACE] = iso_face_parent_idx;
    _iso->isovalue_entity_idx  [PDM_MESH_ENTITY_FACE] = isovalue_face_idx;
    if (isos->extract_kind==PDM_EXTRACT_PART_KIND_REEQUILIBRATE) {
      _iso->iso_entity_parent_gnum[PDM_MESH_ENTITY_FACE] = iso_face_parent_gnum;
    }
    _iso->iso_entity_parent_lnum[PDM_MESH_ENTITY_FACE] = iso_face_parent;
    _iso->iso_connec_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX] = iso_face_vtx_idx;
    _iso->iso_connec    [PDM_CONNECTIVITY_TYPE_FACE_VTX] = iso_face_vtx;
  }

  // Ownerships
  PDM_malloc(_iso->iso_owner_vtx_coord        , isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_vtx_parent_weight, isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_edge_bnd         , isos->n_part, PDM_ownership_t);
  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    PDM_malloc(_iso->iso_owner_gnum               [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_parent_lnum        [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_parent_idx         [i_entity], isos->n_part, PDM_ownership_t);
    PDM_malloc(_iso->iso_owner_isovalue_entity_idx[i_entity], isos->n_part, PDM_ownership_t);
  }
  PDM_malloc(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX], isos->n_part, PDM_ownership_t);
  PDM_malloc(_iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX], isos->n_part, PDM_ownership_t);

  for (int i_part=0; i_part<isos->n_part; ++i_part) {
    _iso->iso_owner_vtx_coord        [i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_vtx_parent_weight[i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_edge_bnd         [i_part] = PDM_OWNERSHIP_KEEP;
    for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
      _iso->iso_owner_gnum               [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_parent_lnum        [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_parent_idx         [i_entity][i_part] = PDM_OWNERSHIP_KEEP;
      _iso->iso_owner_isovalue_entity_idx[i_entity][i_part] = PDM_OWNERSHIP_KEEP;
    }

    _iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part] = PDM_OWNERSHIP_KEEP;
    _iso->iso_owner_connec[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part] = PDM_OWNERSHIP_KEEP;
  }

}

#ifdef  __cplusplus
}
#endif

