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
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_to_part.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"
#include "pdm_vtk.h"

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


static const int bar_pairs[] = {
  0, 1
};

static const int tria_pairs[] = {
  1, 2,
  2, 0,
  0, 1
};

static const int quad_pairs[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0
};

static const int tetra_pairs[] = {
  0, 1,
  0, 2,
  0, 3,
  1, 2,
  1, 3,
  2, 3
};

static const int pyramid_pairs[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  0, 4,
  1, 4,
  2, 4,
  3, 4
};

static const int prism_pairs[] = {
  1, 2,
  2, 0,
  0, 1,
  4, 5,
  5, 3,
  3, 4,
  0, 3,
  1, 4,
  2, 5
};

static const int hexa_pairs[] = {
  0, 1,
  1, 2,
  2, 3,
  3, 0,
  4, 5,
  5, 6,
  6, 7,
  7, 4,
  0, 4,
  1, 5,
  2, 6,
  3, 7
};

static const int *elt_pairs[] = {
  NULL,          // PDM_MESH_NODAL_POINT
  bar_pairs,     // PDM_MESH_NODAL_BAR2
  tria_pairs,    // PDM_MESH_NODAL_TRIA3
  quad_pairs,    // PDM_MESH_NODAL_QUAD4
  NULL,          // PDM_MESH_NODAL_POLY_2D
  tetra_pairs,   // PDM_MESH_NODAL_TETRA4
  pyramid_pairs, // PDM_MESH_NODAL_PYRAMID5
  prism_pairs,   // PDM_MESH_NODAL_PRISM6
  hexa_pairs     // PDM_MESH_NODAL_HEXA8
};


/* From split operator  >>> */
static int _pattern_permutation_2d[8] = {
  0, 0, 1, 0, 2, 2, 1, 0
};


static const double ISOSURFACE_EPS = 1e-6;
static inline int
_cross_0_level
(
  const double v0,
  const double v1
)
{
  return (PDM_ABS(v0) > ISOSURFACE_EPS) && (PDM_ABS(v1) > ISOSURFACE_EPS) && (v0*v1 < 0);
}


static inline int
_cross_any_level
(
  const double v0,
  const double v1,
  const int    n_isovalues,
  const double isovalues[]
)
{
  int n_crossings = 0;
  for (int i = 0; i < n_isovalues; i++) {
    n_crossings += _cross_0_level(v0 - isovalues[i], v1 - isovalues[i]);
  }

  return n_crossings;
}



static void
_count_active_vertices
(
  const int               n_isovalues,
  const double           *isovalues,
  const int               n_vtx,
  const double           *vtx_field,
  int                    *n_vtx_on_vtx
)
{
  *n_vtx_on_vtx = 0;
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
      if (PDM_ABS(vtx_field[i_vtx] - isovalues[i_isovalue]) <= ISOSURFACE_EPS) {
        ++*n_vtx_on_vtx;
      }
    }
  }
}


static void
_build_active_edges
(
  PDM_part_mesh_nodal_t  *pmn,
  const int               i_part,
  const int               n_isovalues,
  const double           *isovalues,
  double                 *vtx_field,
  int                    *n_edge,
  int                   **edge_vtx,
  int                   **elt_edge,
  int                    *n_crossings
)
{
  int debug = 1;

  *n_crossings = 0;
  int n_section = PDM_part_mesh_nodal_n_section_get(pmn);

  // Prepare edge hash table
  const int max_key = 1024; // ?


  if (debug==1) {
    int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    PDM_log_trace_array_double(vtx_field, n_vtx, "vtx_field ::");
  }



  int key_count[max_key];
  int key_idx  [max_key+1];

  PDM_array_reset_int(key_count, max_key, 0);

  /* First loop to count */
  for (int i_section = 0; i_section < n_section; i_section++) {

    // Get section type
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(pmn, i_section);

    // Get table of vertex pairs for current element type
    const int *pairs = elt_pairs[t_elt];
    if (pairs == NULL) {
      // Skip irrelevant sections
      continue;
    }

    // Get section information : number of elts and elt->vtx connectivity
    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, i_section, i_part);

    int         *connec;
    PDM_g_num_t *ln_to_gn;
    int         *parent_num;
    PDM_g_num_t *parent_entity_g_num;

    PDM_part_mesh_nodal_section_std_get(pmn,
                                        i_section,
                                        i_part,
                                        &connec,
                                        &ln_to_gn,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_BAD_VALUE);

    // Increment key occurrences
    int n_pair    = PDM_n_nedge_elt_per_elmt(t_elt);
    int elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int *_connec = connec + elt_n_vtx*i_elt;
      for (int i_pair = 0; i_pair < n_pair; i_pair++) {
        int i_vtx0 = _connec[pairs[2*i_pair  ]];
        int i_vtx1 = _connec[pairs[2*i_pair+1]];

        if (_cross_any_level(vtx_field[i_vtx0-1],
                             vtx_field[i_vtx1-1],
                             n_isovalues,
                             isovalues)) {
          // This is an active edge
          int key = (i_vtx0 + i_vtx1 - 2) % max_key;
          key_count[key]++;
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

  int *key_edge  = malloc(sizeof(int) * key_idx[max_key]);
  int *_edge_vtx = malloc(sizeof(int) * key_idx[max_key] * 2);
  *n_edge = 0;

  /* Second loop to fill */
  for (int i_section = 0; i_section < n_section; i_section++) {

    // Get section type
    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(pmn, i_section);

    // Get table of vertex pairs for current element type
    const int *pairs = elt_pairs[t_elt];
    if (pairs == NULL) {
      // Skip irrelevant sections
      continue;
    }

    // Get section information : number of elts and elt->vtx connectivity
    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, i_section, i_part);

    int         *connec;
    PDM_g_num_t *ln_to_gn;
    int         *parent_num;
    PDM_g_num_t *parent_entity_g_num;

    PDM_part_mesh_nodal_section_std_get(pmn,
                                        i_section,
                                        i_part,
                                        &connec,
                                        &ln_to_gn,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_BAD_VALUE);

    int n_pair    = PDM_n_nedge_elt_per_elmt(t_elt);
    int elt_n_vtx = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, 1);

    elt_edge[i_section] = malloc(sizeof(int) * n_pair*n_elt);

    for (int i_elt = 0; i_elt < n_elt; i_elt++) {
      int *_connec = connec + elt_n_vtx*i_elt;
      for (int i_pair = 0; i_pair < n_pair; i_pair++) {
        int i_vtx0 = _connec[pairs[2*i_pair  ]];
        int i_vtx1 = _connec[pairs[2*i_pair+1]];

        int edge_id = 0;
        int _n_crossings = _cross_any_level(vtx_field[i_vtx0-1],
                                            vtx_field[i_vtx1-1],
                                            n_isovalues,
                                            isovalues);
        if (_n_crossings > 0) {
          int key = (i_vtx0 + i_vtx1 - 2) % max_key;

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

          if (edge_id == 0) {
            // Create new edge
            *n_crossings += _n_crossings;
            _edge_vtx[2*(*n_edge)  ] = i_vtx0;
            _edge_vtx[2*(*n_edge)+1] = i_vtx1;
            key_edge[key_idx[key] + key_count[key]++] = (*n_edge);
            edge_id = ++(*n_edge);
          }
        } // End if active edge

        elt_edge[i_section][n_pair*i_elt+i_pair] = edge_id;
      } // End of loop on pairs
    } // End of loop on elements
  } // End of loop on sections
  free(key_edge);

  *edge_vtx = realloc(_edge_vtx, sizeof(int) * (*n_edge) * 2);
}



static void
_build_iso_vtx_on_active_edge
(
  const double        isovalue,
  const int           n_vtx,
  const PDM_g_num_t  *vtx_gnum,
  const double       *vtx_coord,
  const double       *vtx_field,
  const int           n_edge,
  const int          *edge_vtx,
        int         **out_vtx_to_iso_vtx,
        int         **out_edge_to_iso_vtx,
        int          *out_iso_n_vtx,
        double      **out_iso_vtx_coord,
        int         **out_iso_vtx_parent_idx,
        int         **out_iso_vtx_parent,
        double      **out_iso_vtx_parent_weight,
        PDM_g_num_t **out_iso_vtx_parent_gnum
)
{
  /*
   * TODO: 
   *   - put iso vtx from mesh vtx next to iso vtx
   *     from mesh edge close in memory (in build active edge directly ?)
   */
  int         *vtx_to_iso_vtx        = *out_vtx_to_iso_vtx;
  int         *edge_to_iso_vtx       = *out_edge_to_iso_vtx;
  int          _iso_n_vtx            = *out_iso_n_vtx;
  double      *iso_vtx_coord         = *out_iso_vtx_coord;
  int         *iso_vtx_parent_idx    = *out_iso_vtx_parent_idx;
  int         *iso_vtx_parent        = *out_iso_vtx_parent;
  double      *iso_vtx_parent_weight = *out_iso_vtx_parent_weight;
  PDM_g_num_t *iso_vtx_parent_gnum   = *out_iso_vtx_parent_gnum;


  /*
   * Keep link between mesh vertices and iso vertices
   */  
  PDM_array_reset_int(vtx_to_iso_vtx, n_vtx, 0);
  int tmp_iso_n_vtx = _iso_n_vtx;
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    vtx_to_iso_vtx[i_vtx] = 0;
    if (PDM_ABS(vtx_field[i_vtx] - isovalue) <= ISOSURFACE_EPS) {
      vtx_to_iso_vtx[i_vtx] = ++tmp_iso_n_vtx;
    }
  }


  /*
   * Fill iso vertices infos
   */
  iso_vtx_parent_idx[0] = 0;
 
  // > Iso vertices on mesh vertices
  for (int i_vtx = 0; i_vtx < n_vtx; i_vtx++) {
    if (vtx_to_iso_vtx[i_vtx] > 0) {
      memcpy(iso_vtx_coord + 3*_iso_n_vtx,
             vtx_coord     + 3*i_vtx,
             sizeof(double) * 3);
      iso_vtx_parent_idx   [                   _iso_n_vtx+1] = iso_vtx_parent_idx[_iso_n_vtx] + 1;
      iso_vtx_parent       [iso_vtx_parent_idx[_iso_n_vtx ]] = i_vtx + 1;
      iso_vtx_parent_weight[iso_vtx_parent_idx[_iso_n_vtx ]] = 1.;
      iso_vtx_parent_gnum  [                 2*_iso_n_vtx  ] = vtx_gnum[i_vtx];
      iso_vtx_parent_gnum  [                 2*_iso_n_vtx+1] = vtx_gnum[i_vtx];
      _iso_n_vtx++;
    }
  }
  
  // > Iso vertices on mesh edges
  PDM_array_reset_int(edge_to_iso_vtx, n_edge, 0);

  for (int i_edge = 0; i_edge < n_edge; i_edge++) {
    int i_vtx0 = edge_vtx[2*i_edge  ] - 1;
    int i_vtx1 = edge_vtx[2*i_edge+1] - 1;

    double val0 = vtx_field[i_vtx0] - isovalue;
    double val1 = vtx_field[i_vtx1] - isovalue;

    if (_cross_0_level(val0, val1)) {
      double t = val0 / (val0 - val1);
      for (int i = 0; i < 3; i++) {
        iso_vtx_coord[3*_iso_n_vtx+i] = (1-t)*vtx_coord[3*i_vtx0+i] + t*vtx_coord[3*i_vtx1+i];
      }
      iso_vtx_parent_idx[_iso_n_vtx+1] = iso_vtx_parent_idx[_iso_n_vtx] + 2;

      iso_vtx_parent       [iso_vtx_parent_idx[_iso_n_vtx]  ] = i_vtx0 + 1;
      iso_vtx_parent_weight[iso_vtx_parent_idx[_iso_n_vtx]  ] = 1. - t;

      iso_vtx_parent       [iso_vtx_parent_idx[_iso_n_vtx]+1] = i_vtx1 + 1;
      iso_vtx_parent_weight[iso_vtx_parent_idx[_iso_n_vtx]+1] = t;
      
      iso_vtx_parent_gnum  [2*_iso_n_vtx  ] = vtx_gnum[i_vtx0];
      iso_vtx_parent_gnum  [2*_iso_n_vtx+1] = vtx_gnum[i_vtx1];
      
      edge_to_iso_vtx[i_edge] = ++_iso_n_vtx;
    }

  } // End of loop on edges


  /*
   * Output
   */
  *out_vtx_to_iso_vtx        = vtx_to_iso_vtx;
  *out_edge_to_iso_vtx       = edge_to_iso_vtx;
  *out_iso_n_vtx             = _iso_n_vtx;
  *out_iso_vtx_coord         = iso_vtx_coord;
  *out_iso_vtx_parent_idx    = iso_vtx_parent_idx;
  *out_iso_vtx_parent        = iso_vtx_parent;
  *out_iso_vtx_parent_weight = iso_vtx_parent_weight;
  *out_iso_vtx_parent_gnum   = iso_vtx_parent_gnum;
}



static void
_contouring_triangles
(
  int           n_elt,
  int          *elt_vtx,
  PDM_g_num_t  *elt_gnum,
  PDM_g_num_t   max_elt_gnum,
  int          *parent_num,
  int          *elt_edge,
  PDM_g_num_t  *vtx_gnum,
  double       *vtx_field,
  int          *vtx_to_iso_vtx,
  int          *edge_to_iso_vtx,
  int          *out_iso_n_edge,
  int         **out_iso_edge_vtx,
  int         **out_iso_edge_parent,
  PDM_g_num_t **out_iso_edge_parent_gnum
)
{
  /*
   * TODO:
   *   - is parent_gnum = (edge_gnum, edge_gnum) or (edge_gnum, max_edge_num + vtx_gnum) or (max_edge_num + vtx_gnum, max_edge_num + vtx_gnum)
   *      --> how is that ???
   */

  /* First loop to count */
  int iso_n_edge = *out_iso_n_edge;
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {

    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 3; i_edge++) {
      char bit = 0;
      if (elt_edge[3*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[3*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges

    if (pattern > 0 && pattern < 7) {
      iso_n_edge++;
    }

  } // End of loop on elements

  /* Allocate */
  int         *iso_edge_vtx         = realloc(*out_iso_edge_vtx,         sizeof(int        ) * iso_n_edge * 2);
  int         *iso_edge_parent      = realloc(*out_iso_edge_parent,      sizeof(int        ) * iso_n_edge    );
  PDM_g_num_t *iso_edge_parent_gnum = realloc(*out_iso_edge_parent_gnum, sizeof(PDM_g_num_t) * iso_n_edge * 2);

  iso_n_edge = *out_iso_n_edge;


  /* Second loop to fill */
  for (int i_elt = 0; i_elt < n_elt; i_elt++) {

    unsigned char pattern = 0;
    for (int i_edge = 0; i_edge < 3; i_edge++) {
      char bit = 0;
      if (elt_edge[3*i_elt+i_edge] != 0) {
        bit = (edge_to_iso_vtx[PDM_ABS(elt_edge[3*i_elt+i_edge]) - 1] != 0);
      }
      pattern |= (bit << i_edge);
    } // End of loop on edges

    if (pattern == 0 || pattern == 7) {
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
      case 1: case 2: case 4: {
        int vtx_on_vtx0 = vtx_to_iso_vtx[perm_elt_vtx[0]];
        assert(vtx_on_vtx0);
        int edge_id0 = PDM_ABS(perm_elt_edge[0]) - 1;
        int vtx_on_edge0 = edge_to_iso_vtx[edge_id0];
        assert(vtx_on_edge0);

        if (vtx_field[perm_elt_vtx[1]] > vtx_field[perm_elt_vtx[2]]) {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_vtx0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_edge0;

          iso_edge_parent_gnum[2*iso_n_edge  ] = max_elt_gnum + vtx_gnum[perm_elt_vtx[2]];
          iso_edge_parent_gnum[2*iso_n_edge+1] = elt_gnum[edge_id0];
        }
        else {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_edge0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_vtx0;

          iso_edge_parent_gnum[2*iso_n_edge  ] = elt_gnum[edge_id0];
          iso_edge_parent_gnum[2*iso_n_edge+1] = max_elt_gnum + vtx_gnum[perm_elt_vtx[2]];
        }
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

          iso_edge_parent_gnum[2*iso_n_edge  ] = elt_gnum[edge_id1];
          iso_edge_parent_gnum[2*iso_n_edge+1] = elt_gnum[edge_id0];
        }
        else {
          iso_edge_vtx[2*iso_n_edge  ] = vtx_on_edge0;
          iso_edge_vtx[2*iso_n_edge+1] = vtx_on_edge1;

          iso_edge_parent_gnum[2*iso_n_edge  ] = elt_gnum[edge_id0];
          iso_edge_parent_gnum[2*iso_n_edge+1] = elt_gnum[edge_id1];
        }
        break;
      }

      default: {
        break;
      }
    }

    iso_edge_parent[iso_n_edge++] = i_elt + 1;

  } // End of loop on elements


  if (parent_num) {
    for (int i = *out_iso_n_edge; i < iso_n_edge; i++) {
      iso_edge_parent[i] = parent_num[iso_edge_parent[i] - 1];
    }
  }


  *out_iso_n_edge           = iso_n_edge;
  *out_iso_edge_vtx         = iso_edge_vtx;
  *out_iso_edge_parent      = iso_edge_parent;
  *out_iso_edge_parent_gnum = iso_edge_parent_gnum;
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
  int debug = 1;
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

  // > Part mesh nodal
  PDM_part_mesh_nodal_t *pmn = isos->pmesh_nodal;

  // > Isosurface information
  int   n_isovalues = isos->n_isovalues[id_iso];
  double *isovalues = isos->  isovalues[id_iso];

  double **vtx_field = isos->field[id_iso];

  // > Mesh information
  int mesh_dim  = PDM_part_mesh_nodal_mesh_dimension_get(pmn);
  int n_part    = PDM_part_mesh_nodal_n_part_get        (pmn);
  int n_section = PDM_part_mesh_nodal_n_section_get     (pmn);


  // > Allocate isosurface vertices
  int          *iso_n_vtx             = PDM_array_zeros_int(n_part);
  double      **iso_vtx_coord         = malloc(sizeof(double      *) * n_part);
  int         **iso_vtx_parent_idx    = malloc(sizeof(int         *) * n_part);
  int         **iso_vtx_parent        = malloc(sizeof(int         *) * n_part);
  double      **iso_vtx_parent_weight = malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **iso_vtx_parent_gnum   = malloc(sizeof(PDM_g_num_t *) * n_part);

  // > Allocate isosurface edges
  int          *iso_n_edge           = PDM_array_zeros_int(n_part);
  int         **iso_edge_vtx         = malloc(sizeof(int         *) * n_part);
  int         **iso_edge_parent      = malloc(sizeof(int         *) * n_part);
  PDM_g_num_t **iso_edge_parent_gnum = malloc(sizeof(PDM_g_num_t *) * n_part);
  int         **isovalue_edge_idx    = malloc(sizeof(int         *) * n_part);

  // > Allocate face edges
  // TODO: transformer en part_mesh_nodal (tetra sera que TRI et QUAD)
  //       quid des autres elements que les TETRA ? POLY_2D ? (MR gitlab)
  int          *iso_n_face           = PDM_array_zeros_int(n_part);
  int         **iso_face_vtx_idx     = malloc(sizeof(int         *) * n_part);
  int         **iso_face_vtx         = malloc(sizeof(int         *) * n_part);
  int         **iso_face_parent      = malloc(sizeof(int         *) * n_part);
  PDM_g_num_t **iso_face_parent_gnum = malloc(sizeof(PDM_g_num_t *) * n_part);
  int         **isovalue_face_idx    = malloc(sizeof(int         *) * n_part);



  for (int i_part = 0; i_part < n_part; i_part++) {

    // > Get vtx info
    int          n_vtx     = PDM_part_mesh_nodal_n_vtx_get    (pmn, i_part);
    PDM_g_num_t *vtx_gnum  = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);
    double      *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
    
    // > Active edges for each sections
    int **elt_edge = malloc(sizeof(int *) * n_section);
    for (int i_section = 0; i_section < n_section; i_section++) {
      elt_edge[i_section] = NULL;
    }


    /*
     * Build active edges from entry mesh (i.e. crossing at least one isosurface)
     */
    t_start = PDM_MPI_Wtime();
    
    int  n_edge       = 0;
    int *edge_vtx     = NULL;
    int  n_crossings  = 0;
    
    _build_active_edges(pmn,
                        i_part,
                        n_isovalues,
                        isovalues,
                        vtx_field[i_part],
                        &n_edge,
                        &edge_vtx,
                        elt_edge,
                        &n_crossings);

    t_end = PDM_MPI_Wtime();
    printf("   build active edges : %.3fs  (%d active edges)\n", 
           t_end - t_start,
           n_edge);
    fflush(stdout);


    if (debug==1) {
      log_trace("n_active_edge = %d\n", n_edge);
      log_trace("n_crossings   = %d\n", n_crossings);

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
     * Count active vertices from entry mesh (i.e. intersecting at least one isosurface)
     */
    t_start = PDM_MPI_Wtime();

    int n_vtx_on_vtx = 0;
    _count_active_vertices( n_isovalues,
                            isovalues,
                            n_vtx,
                            vtx_field[i_part],
                           &n_vtx_on_vtx);

    t_end = PDM_MPI_Wtime();
    printf("   count active vertices : %.3fs  (%d active vertices)\n", 
           t_end - t_start,
           n_vtx_on_vtx);
    fflush(stdout);

    if (debug==1) {
      log_trace("n_vtx_on_vtx  = %d\n", n_vtx_on_vtx);
    }


    // > Allocate vertex
    iso_n_vtx[i_part] = n_vtx_on_vtx + n_crossings;
    iso_vtx_coord        [i_part] = malloc(sizeof(double     ) *  iso_n_vtx[i_part] * 3);
    iso_vtx_parent_idx   [i_part] = malloc(sizeof(int        ) * (iso_n_vtx[i_part] + 1));
    iso_vtx_parent       [i_part] = malloc(sizeof(int        ) * (n_vtx_on_vtx + 2*n_crossings));
    iso_vtx_parent_weight[i_part] = malloc(sizeof(double     ) * (n_vtx_on_vtx + 2*n_crossings));
    iso_vtx_parent_gnum  [i_part] = malloc(sizeof(PDM_g_num_t) *  iso_n_vtx[i_part] * 2);
    int *vtx_to_iso_vtx           = malloc(sizeof(int        ) * n_vtx);
    int *edge_to_iso_vtx          = malloc(sizeof(int        ) * n_edge); // TODO: why set at 0 in proto


    // > Allocate edge
    iso_edge_vtx        [i_part] = NULL;
    iso_edge_parent     [i_part] = NULL;
    iso_edge_parent_gnum[i_part] = NULL;
    isovalue_edge_idx   [i_part] = malloc(sizeof(int) * (n_isovalues + 1));
    isovalue_edge_idx   [i_part][0] = 0;


    // > Allocate face
    iso_face_vtx_idx    [i_part] = NULL;
    iso_face_vtx        [i_part] = NULL;
    iso_face_parent     [i_part] = NULL;
    iso_face_parent_gnum[i_part] = NULL;
    isovalue_face_idx   [i_part] = malloc(sizeof(int) * (n_isovalues + 1));
    isovalue_face_idx   [i_part][0] = 0;



    /*
     * Build isosurface mesh
     */
    iso_n_vtx[i_part] = 0;
    for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
      
      if (debug) {
        log_trace("Building i_isovalue %d\n", i_isovalue);
      }
      t_start = PDM_MPI_Wtime();


      /* Build iso-vertices on active edges */
      _build_iso_vtx_on_active_edge(isovalues[i_isovalue],
                                    n_vtx,
                                    vtx_gnum,
                                    vtx_coord,
                                    vtx_field[i_part],
                                    n_edge,
                                    edge_vtx,
                                   &vtx_to_iso_vtx,
                                   &edge_to_iso_vtx,
                                   &(iso_n_vtx[i_part]),
                                   &iso_vtx_coord[i_part],
                                   &iso_vtx_parent_idx[i_part],
                                   &iso_vtx_parent[i_part],
                                   &iso_vtx_parent_weight[i_part],
                                   &iso_vtx_parent_gnum[i_part]
                                   );
      if (debug) {
        log_trace("\n");
        log_trace("Vertices \n");
        log_trace("iso_n_vtx[i_part] = %d \n", iso_n_vtx[i_part]);
        // PDM_log_trace_array_double(iso_vtx_coord[i_part], iso_n_vtx[i_part]*3, "iso_vtx_coord ::");
      }

      for (int i_section = 0; i_section < n_section; i_section++) {

        // Get section type
        PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_section_elt_type_get(pmn, i_section);

        // Get section information : number of elts and elt->vtx connectivity
        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, i_section, i_part);

        int         *connec;
        PDM_g_num_t *gnum;
        int         *parent_num;
        PDM_g_num_t *parent_entity_g_num;
        PDM_part_mesh_nodal_section_std_get(pmn,
                                            i_section,
                                            i_part,
                                            &connec,
                                            &gnum,
                                            &parent_num,
                                            &parent_entity_g_num,
                                            PDM_OWNERSHIP_BAD_VALUE);



        switch (t_elt) {

          case PDM_MESH_NODAL_TRIA3: {
            PDM_g_num_t l_max_elt_gnum = PDM_array_max_gnum(gnum, n_elt);
            PDM_g_num_t   max_elt_gnum = 0;
            PDM_MPI_Allreduce (&l_max_elt_gnum, &max_elt_gnum, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, isos->comm);
            
            _contouring_triangles(n_elt,
                                  connec,
                                  gnum,
                                  max_elt_gnum,
                                  parent_num,
                                  elt_edge[i_section],
                                  vtx_gnum,
                                  vtx_field[i_part],
                                  vtx_to_iso_vtx,
                                  edge_to_iso_vtx,
                                  &iso_n_edge          [i_part],
                                  &iso_edge_vtx        [i_part],
                                  &iso_edge_parent     [i_part],
                                  &iso_edge_parent_gnum[i_part]);
            break;
          }

          default: {
            printf("Warning: Contouring not yet implemented for elt type %d\n", t_elt);
            continue;
          }
        }

      } // End of loop on sections
      
      isovalue_edge_idx[i_part][i_isovalue+1] = iso_n_edge[i_part];
      
      t_end = PDM_MPI_Wtime();
      printf("   isovalue %2d : %.3fs  (%6d vertices, %6d edges, %6d faces)\n",
             i_isovalue,
             t_end - t_start,
             iso_n_vtx[i_part],
             isovalue_edge_idx[i_part][i_isovalue+1] - isovalue_edge_idx[i_part][i_isovalue],
             isovalue_face_idx[i_part][i_isovalue+1] - isovalue_face_idx[i_part][i_isovalue]);
      fflush(stdout);

    } // End of loop on isovalues

    /* 
     * Visu
     */
    if (debug) {
      char outname[999];

      log_trace("n_vertices = %d \n", iso_n_vtx[i_part]);

      int *iso_edge_isovalue = malloc(sizeof(int) * iso_n_edge[i_part]);
      for (int i_isovalue = 0; i_isovalue < n_isovalues; i_isovalue++) {
        for (int i = isovalue_edge_idx[i_part][i_isovalue]; i < isovalue_edge_idx[i_part][i_isovalue+1]; i++) {
          iso_edge_isovalue[i] = i_isovalue;
        }
      }

      const char *field_name [] = {"parent", "i_isovalue"};
      const int  *field_value[] = {iso_edge_parent[i_part], iso_edge_isovalue};

      sprintf(outname, "marching_sections_iso_edges_%d.vtk", i_part);
      PDM_vtk_write_std_elements(outname,
                                 iso_n_vtx[i_part],
                                 iso_vtx_coord[i_part],
                                 NULL,
                                 PDM_MESH_NODAL_BAR2,
                                 iso_n_edge[i_part],
                                 iso_edge_vtx[i_part],
                                 NULL,
                                 2,
                                 field_name,
                                 field_value);
      free(iso_edge_isovalue);
    } // End visu


    // > Free temp arrays
    for (int i_section = 0; i_section < n_section; i_section++) {
      free(elt_edge[i_section]);
    }
    free(elt_edge);
    free(edge_vtx);
    free(vtx_to_iso_vtx);
    free(edge_to_iso_vtx);

  } // End of loop on partitions


  /*
   * Build vertices gnum
   */
  PDM_gen_gnum_t *gen_gnum_vtx  = PDM_gnum_create(mesh_dim-1, n_part, PDM_FALSE, 1., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gen_gnum_t *gen_gnum_edge = PDM_gnum_create(mesh_dim-1, n_part, PDM_FALSE, 1., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_parents_nuplet(gen_gnum_vtx , 2);
  PDM_gnum_set_parents_nuplet(gen_gnum_edge, 2);

  for (int i_part=0; i_part<n_part; i_part++) {
    PDM_gnum_set_from_parents(gen_gnum_vtx , i_part, iso_n_vtx [i_part], iso_vtx_parent_gnum [i_part]);
    PDM_gnum_set_from_parents(gen_gnum_edge, i_part, iso_n_edge[i_part], iso_edge_parent_gnum[i_part]);
    PDM_log_trace_array_long(iso_vtx_parent_gnum [i_part], iso_n_vtx [i_part], "iso_vtx_parent_gnum  :: ");
    PDM_log_trace_array_long(iso_edge_parent_gnum[i_part], iso_n_edge[i_part], "iso_edge_parent_gnum :: ");
  }
  PDM_gnum_compute(gen_gnum_vtx);
  PDM_gnum_compute(gen_gnum_edge);

  PDM_g_num_t **iso_vtx_gnum  = malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **iso_edge_gnum = malloc(sizeof(PDM_g_num_t *) * n_part);
  for (int i_part=0; i_part<n_part; i_part++) {
    iso_vtx_gnum [i_part] = PDM_gnum_get(gen_gnum_vtx , i_part);
    iso_edge_gnum[i_part] = PDM_gnum_get(gen_gnum_edge, i_part);
  }

  for (int i_part=0; i_part<n_part; i_part++) {
    free(iso_vtx_parent_gnum[i_part]);
    free(iso_edge_parent_gnum[i_part]);
    free(iso_face_parent_gnum[i_part]);
  }
  free(iso_vtx_parent_gnum);
  free(iso_edge_parent_gnum);
  free(iso_face_parent_gnum);
  PDM_gnum_free(gen_gnum_vtx );
  PDM_gnum_free(gen_gnum_edge);

  /*
   * Store isosurface in part_mesh_nodal
   *
   * TODO:
   */
  isos->iso_n_part = n_part;

  isos->iso_n_vtx            [id_iso] = iso_n_vtx;
  isos->iso_vtx_coord        [id_iso] = iso_vtx_coord;
  isos->iso_vtx_gnum         [id_iso] = iso_vtx_gnum;
  isos->iso_vtx_lparent_idx  [id_iso] = iso_vtx_parent_idx;
  isos->iso_vtx_lparent      [id_iso] = iso_vtx_parent;
  isos->iso_vtx_parent_weight[id_iso] = iso_vtx_parent_weight;
  
  isos->iso_n_edge       [id_iso] = iso_n_edge;
  isos->iso_edge_vtx     [id_iso] = iso_edge_vtx;
  isos->iso_edge_gnum    [id_iso] = iso_edge_gnum;
  isos->iso_edge_lparent  [id_iso] = iso_edge_parent;
  isos->isovalue_edge_idx[id_iso] = isovalue_edge_idx;
  
  isos->iso_n_face      [id_iso] = iso_n_face;
  isos->iso_face_vtx_idx[id_iso] = iso_face_vtx_idx;
  isos->iso_face_vtx    [id_iso] = iso_face_vtx;
  // isos->iso_face_gnum   [id_iso] = iso_face_gnum;
  isos->iso_face_lparent [id_iso] = iso_face_parent;
  isos->isovalue_face_idx[id_iso] = isovalue_face_idx;

  isos->iso_owner_vtx_coord         [id_iso] = malloc(sizeof(PDM_ownership_t  ) * n_part);
  isos->iso_owner_vtx_parent_weight [id_iso] = malloc(sizeof(PDM_ownership_t  ) * n_part);
  isos->iso_owner_gnum              [id_iso] = malloc(sizeof(PDM_ownership_t *) * n_part);
  isos->iso_owner_connec            [id_iso] = malloc(sizeof(PDM_ownership_t *) * n_part);
  isos->iso_owner_lparent           [id_iso] = malloc(sizeof(PDM_ownership_t *) * n_part);
  for (int i_part=0; i_part<n_part; ++i_part) {
    isos->iso_owner_vtx_coord         [id_iso][i_part] = PDM_OWNERSHIP_KEEP;
    isos->iso_owner_vtx_parent_weight [id_iso][i_part] = PDM_OWNERSHIP_KEEP;

    isos->iso_owner_gnum              [id_iso][i_part] = malloc(sizeof(PDM_ownership_t *) * PDM_MESH_ENTITY_MAX );
    isos->iso_owner_connec            [id_iso][i_part] = malloc(sizeof(PDM_ownership_t *) * PDM_MESH_ENTITY_MAX );
    isos->iso_owner_lparent           [id_iso][i_part] = malloc(sizeof(PDM_ownership_t *) * PDM_MESH_ENTITY_MAX );
    for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
      isos->iso_owner_gnum   [id_iso][i_part][i_entity] = PDM_OWNERSHIP_KEEP;
      isos->iso_owner_connec [id_iso][i_part][i_entity] = PDM_OWNERSHIP_KEEP;
      isos->iso_owner_lparent[id_iso][i_part][i_entity] = PDM_OWNERSHIP_KEEP;
    }
  }

}

#ifdef  __cplusplus
}
#endif

