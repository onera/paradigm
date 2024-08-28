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

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"
#include "pdm_unique.h"
#include "pdm_binary_search.h"
#include "pdm_order.h"
#include "pdm_error.h"
#include "pdm_part_extension.h"
#include "pdm_part_to_part.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_extension_algorithm.h"
#include "pdm_domain_utils.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_vtk.h"

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
static void
_hexa_ngon_to_nodal
(
 int   n_cell,
 int  *cell_face_idx,
 int  *cell_face,
 int  *face_vtx_idx,
 int  *face_vtx,
 int **cell_vtx
 )
{
  int debug = 1;

  PDM_malloc(*cell_vtx, n_cell * 8, int);

  for (int i_cell = 0; i_cell < n_cell; i_cell++) {

    if (debug==1) {
      log_trace("\nCell %d\n", i_cell+1);
      for (int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; idx_face++) {
        int iface = PDM_ABS(cell_face[idx_face]) - 1;
        log_trace("  face %d: %d %d %d %d\n",
                  cell_face[idx_face],
                  face_vtx[4*iface], face_vtx[4*iface+1], face_vtx[4*iface+2], face_vtx[4*iface+3]);
      }
    }

    int *cv = (*cell_vtx) + 8*i_cell;
    int *cf = cell_face + cell_face_idx[i_cell];
    assert(cell_face_idx[i_cell+1] - cell_face_idx[i_cell] == 6);

    int i_face = PDM_ABS(cf[0]) - 1;
    assert(face_vtx_idx[i_face+1] - face_vtx_idx[i_face] == 4);
    if (cf[0] < 0) {
      for (int i = 0; i < 4; i++) {
        cv[  i] = face_vtx[face_vtx_idx[i_face] + i];
      }
    }
    else {
      for (int i = 0; i < 4; i++) {
        cv[3-i] = face_vtx[face_vtx_idx[i_face] + i];
      }
    }

    int found[2] = {0, 0};
    for (int i = 1; i < 6; i++) {
      i_face = PDM_ABS(cf[i]) - 1;
      assert(face_vtx_idx[i_face+1] - face_vtx_idx[i_face] == 4);
      int *fv = face_vtx + face_vtx_idx[i_face];

      int skip = 0;
      for (int j = 0; j < 4; j++) {
        int i_vtx0, i_vtx1;
        if (cf[i] < 0) {
          i_vtx0 = fv[(j+1)%4];
          i_vtx1 = fv[j      ];
        }
        else {
          i_vtx0 = fv[j      ];
          i_vtx1 = fv[(j+1)%4];
        }

        if (!found[0]) {
          if (i_vtx0 == cv[0] && i_vtx1 == cv[1]) {
            if (cf[i] < 0) {
              cv[4] = fv[(j+2)%4];
              cv[5] = fv[(j+3)%4];
            }
            else {
              cv[4] = fv[(j+3)%4];
              cv[5] = fv[(j+2)%4];
            }
            found[0] = 1;
            skip = 1;
            break;
          }
        }
      }
      if (skip) continue;

      if (!found[1]) {
        for (int j = 0; j < 4; j++) {
          int i_vtx0, i_vtx1;
          if (cf[i] < 0) {
            i_vtx0 = fv[j      ];
            i_vtx1 = fv[(j+1)%4];
          }
          else {
            i_vtx0 = fv[(j+1)%4];
            i_vtx1 = fv[j      ];
          }

          if (i_vtx0 == cv[3] && i_vtx1 == cv[2]) {
            if (cf[i] < 0) {
              cv[6] = fv[(j+2)%4];
              cv[7] = fv[(j+3)%4];
            }
            else {
              cv[6] = fv[(j+3)%4];
              cv[7] = fv[(j+2)%4];
            }
            found[1] = 1;
            break;
          }
        }
      }

    }

    assert(found[0] && found[1]);

  }
}


static
void
_update_propagating_graph_for_depth
(
  int    depth,
  int    n_part,
  int   *pn_entity,
  int  **pn_entity_by_depth,
  int  **pentity_to_entity_idx,
  int  **pentity_to_entity_triplet,
  int  **pentity_to_entity_interface,
  int ***new_pentity_to_entity_idx,
  int ***new_pentity_to_entity_triplet,
  int ***new_pentity_to_entity_interface
)
{
  int **_new_pentity_to_entity_idx       = *new_pentity_to_entity_idx;
  int **_new_pentity_to_entity_triplet   = *new_pentity_to_entity_triplet;
  int **_new_pentity_to_entity_interface = *new_pentity_to_entity_interface;

  for(int i_part=0; i_part<n_part; ++i_part) {

    // > Remove previous graph
    PDM_free(_new_pentity_to_entity_idx      [i_part]);
    PDM_free(_new_pentity_to_entity_triplet  [i_part]);
    PDM_free(_new_pentity_to_entity_interface[i_part]);
    
    int beg_entity = pn_entity[i_part] - pn_entity_by_depth[depth][i_part];
    int shift = 0;
    for(int k = 0; k < depth; ++k) {
      shift += pn_entity_by_depth[k][i_part];
    }

    int *pentity_to_entity_n = NULL;
    PDM_malloc(pentity_to_entity_n               , pn_entity[i_part]  , int);
    PDM_malloc(_new_pentity_to_entity_idx[i_part], pn_entity[i_part]+1, int);
    _new_pentity_to_entity_idx[i_part][0] = 0;
    
    for(int i=0; i<pn_entity[i_part]; ++i) {
      pentity_to_entity_n               [i  ] = 0;
      _new_pentity_to_entity_idx[i_part][i+1] = _new_pentity_to_entity_idx[i_part][i];
    }

    for(int i=0; i<pn_entity_by_depth[depth][i_part]; ++i) {
      int idx_write = beg_entity+i;
      int n_elmt = pentity_to_entity_idx[i_part][shift+i+1] - pentity_to_entity_idx[i_part][shift+i];
      _new_pentity_to_entity_idx[i_part][idx_write+1] = _new_pentity_to_entity_idx[i_part][idx_write] + n_elmt;
    }

    int n_new_size = _new_pentity_to_entity_idx[i_part][pn_entity[i_part]]/3;
    PDM_malloc(_new_pentity_to_entity_triplet  [i_part], 3 * n_new_size, int);
    PDM_malloc(_new_pentity_to_entity_interface[i_part],     n_new_size, int);

    for(int i = 0; i < pn_entity_by_depth[depth][i_part]; ++i) {

      for(int j = pentity_to_entity_idx[i_part][shift+i]/3; j < pentity_to_entity_idx[i_part][shift+i+1]/3; ++j) {

        int idx_write = _new_pentity_to_entity_idx[i_part][beg_entity+i]/3+
                             pentity_to_entity_n          [beg_entity+i];

        _new_pentity_to_entity_triplet  [i_part][3*idx_write  ] = pentity_to_entity_triplet  [i_part][3*j  ];
        _new_pentity_to_entity_triplet  [i_part][3*idx_write+1] = pentity_to_entity_triplet  [i_part][3*j+1];
        _new_pentity_to_entity_triplet  [i_part][3*idx_write+2] = pentity_to_entity_triplet  [i_part][3*j+2];
        _new_pentity_to_entity_interface[i_part][  idx_write  ] = pentity_to_entity_interface[i_part][  j  ];

        pentity_to_entity_n[beg_entity+i] ++;
      }
    }

    PDM_free(pentity_to_entity_n);
  } // End loop on partitions
}



// static
// void
// _concat_int_idx_array_current_with_extended
// (
//   int            ln_part_tot,
//   int           *pn_entity,
//   int           *pn_entity_extended,
//   int          **pentity_idx_array,
//   int          **pentity_extended_idx_array
// )
// {
//   for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
//     int pn_concat_entity  = pn_entity [i_part] +1 + pn_entity_extended [i_part] +1;
//     int prev_last_idx = pentity_idx_array[i_part][pn_entity[i_part]];
//     PDM_realloc(pentity_idx_array[i_part], pentity_idx_array[i_part], pn_concat_entity+1, int);
//     for(int i_entity = 0; i_entity < pn_entity_extended[i_part]; ++i_entity) {
//       pentity_idx_array[i_part][pn_entity[i_part]+1+i_entity] = pentity_extended_idx_array[i_part][i_entity+1]+prev_last_idx;
//     }
//   }
// }

// static
// void
// _concat_int_array_current_with_extended_from_idx
// (
//   int            ln_part_tot,
//   int           *pn_entity,
//   int           *pn_entity_extended,
//   int          **pentity_array_idx,
//   int          **pentity_array,
//   int          **pentity_extended_array_idx,
//   int          **pentity_extended_array
// )
// {
//   for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
//     int          array_size = pentity_array_idx         [i_part][pn_entity         [i_part]];
//     int extended_array_size = pentity_extended_array_idx[i_part][pn_entity_extended[i_part]];
//     PDM_realloc(pentity_array[i_part], pentity_array[i_part], array_size+extended_array_size, int);
//     for(int i_entity = 0; i_entity < extended_array_size; ++i_entity) {
//       pentity_array[i_part][array_size+i_entity] = pentity_extended_array[i_part][i_entity];
//     }
//   }
// }

static
void
_concat_int_array_current_with_extended
(
  int            ln_part_tot,
  int           *pn_entity,
  int           *pn_entity_extended,
  int          **pentity_array,
  int          **pentity_extended_array
)
{
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int pn_concat_entity  = pn_entity [i_part] + pn_entity_extended [i_part];
    PDM_realloc(pentity_array[i_part], pentity_array[i_part], pn_concat_entity, int);
    for(int i_entity = 0; i_entity < pn_entity_extended[i_part]; ++i_entity) {
      pentity_array[i_part][pn_entity[i_part]+i_entity] = pentity_extended_array[i_part][i_entity];
    }
  }
}

// static
// void
// _concat_gnum_array_current_with_extended
// (
//   int            ln_part_tot,
//   int           *pn_entity,
//   int           *pn_entity_extended,
//   PDM_g_num_t  **pentity_array,
//   PDM_g_num_t  **pentity_extended_array
// )
// {
//   for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
//     int pn_concat_entity  = pn_entity [i_part] + pn_entity_extended [i_part];
//     PDM_realloc(pentity_array[i_part], pentity_array[i_part], pn_concat_entity, int);
//     for(int i_entity = 0; i_entity < pn_entity_extended[i_part]; ++i_entity) {
//       pentity_array[i_part][pn_entity[i_part]+i_entity] = pentity_extended_array[i_part][i_entity];
//     }
//   }
// }


static
void
_concat_ln_to_gn_current_with_extended
(
  int            ln_part_tot,
  int           *pn_entity,
  int           *pn_entity_extended,
  PDM_g_num_t  **pentity_ln_to_gn,
  PDM_g_num_t  **pentity_extended_ln_to_gn,
  PDM_g_num_t   *shift_by_domain_entity
)
{
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int pn_concat_entity  = pn_entity [i_part] + pn_entity_extended [i_part];
    PDM_realloc(pentity_ln_to_gn[i_part], pentity_ln_to_gn[i_part], pn_concat_entity, PDM_g_num_t);
    for(int i_entity = 0; i_entity < pn_entity_extended[i_part]; ++i_entity) {
      pentity_ln_to_gn[i_part][pn_entity[i_part]+i_entity] = pentity_extended_ln_to_gn[i_part][i_entity];
      *shift_by_domain_entity = PDM_MAX(*shift_by_domain_entity, pentity_extended_ln_to_gn[i_part][i_entity]);
    }
  }
}


static
void
_concat_coords_current_with_extended
(
  int            ln_part_tot,
  int           *pn_entity,
  int           *pn_entity_extended,
  double       **pentity_coords,
  double       **pentity_extended_coords
)
{
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int pn_concat_entity  = pn_entity [i_part] + pn_entity_extended [i_part];
    PDM_realloc(pentity_coords[i_part], pentity_coords[i_part], 3 * pn_concat_entity, double);
    for(int i_entity = 0; i_entity < 3 * pn_entity_extended[i_part]; ++i_entity) {
      pentity_coords[i_part][3*pn_entity[i_part]+i_entity] = pentity_extended_coords[i_part][i_entity];
    }
  }
}


static
void
_concat_connectivity_with_extended
(
  int            ln_part_tot,
  int           *pn_entity1,
  int           *pn_entity1_extended,
  int          **pentity1_entity2_idx,
  int          **pentity1_entity2,
  int          **pextended_entity1_entity2_idx,
  int          **pextended_entity1_entity2
)
{
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int pn_concat_entity  = pn_entity1 [i_part] + pn_entity1_extended [i_part];
    int pn_concat_entity1_entity2_idx = pentity1_entity2_idx [i_part][pn_entity1[i_part]] + pextended_entity1_entity2_idx[i_part][pn_entity1_extended[i_part]];

    PDM_realloc(pentity1_entity2_idx[i_part], pentity1_entity2_idx[i_part], pn_concat_entity+1           , int);
    PDM_realloc(pentity1_entity2    [i_part], pentity1_entity2    [i_part], pn_concat_entity1_entity2_idx, int);

    for(int i_face = 0; i_face < pn_entity1_extended[i_part]; ++i_face) {
      int ln_vtx = pextended_entity1_entity2_idx[i_part][i_face+1] - pextended_entity1_entity2_idx[i_part][i_face];
      pentity1_entity2_idx [i_part][pn_entity1[i_part]+i_face+1] = pentity1_entity2_idx [i_part][pn_entity1[i_part]+i_face] + ln_vtx;
    }

    /* Concatenate graphe in other array */
    for(int i = 0; i < pextended_entity1_entity2_idx[i_part][pn_entity1_extended[i_part]]; ++i) {
      pentity1_entity2[i_part][pentity1_entity2_idx [i_part][pn_entity1[i_part]]+i] = pextended_entity1_entity2[i_part][i];
    }

  }
}

static
void
_concat_full_with_extended
(
  int             ln_part_tot,
  int            *pfull_n_entity_extended,
  int            *pn_entity_extended,
  int            *pn_entity_extended_old,
  PDM_g_num_t   **pentity_extended_ln_to_gn,
  int           **pentity_extended_to_pentity_idx,
  int           **pentity_extended_to_pentity_triplet,
  int           **pentity_extended_to_pentity_interface,
  PDM_g_num_t   **pfull_entity_extended_ln_to_gn,
  int           **pfull_entity_extended_to_pentity_idx,
  int           **pfull_entity_extended_to_pentity_triplet,
  int           **pfull_entity_extended_to_pentity_interface

)
{
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int size_entity_entity = (pfull_entity_extended_to_pentity_idx[i_part][pn_entity_extended_old[i_part]] + pentity_extended_to_pentity_idx[i_part][pn_entity_extended[i_part]])/3;
    PDM_realloc(pfull_entity_extended_ln_to_gn            [i_part], pfull_entity_extended_ln_to_gn            [i_part], pfull_n_entity_extended[i_part]  , PDM_g_num_t);
    PDM_realloc(pfull_entity_extended_to_pentity_idx      [i_part], pfull_entity_extended_to_pentity_idx      [i_part], pfull_n_entity_extended[i_part]+1, int        );
    PDM_realloc(pfull_entity_extended_to_pentity_triplet  [i_part], pfull_entity_extended_to_pentity_triplet  [i_part], 3 * size_entity_entity           , int        );
    PDM_realloc(pfull_entity_extended_to_pentity_interface[i_part], pfull_entity_extended_to_pentity_interface[i_part],     size_entity_entity           , int        );

    for(int i_entity = 0; i_entity < pn_entity_extended[i_part]; ++i_entity) {
      int ln_entity = pentity_extended_to_pentity_idx[i_part][i_entity+1] - pentity_extended_to_pentity_idx[i_part][i_entity];
      pfull_entity_extended_to_pentity_idx[i_part][pn_entity_extended_old[i_part]+i_entity+1] = pfull_entity_extended_to_pentity_idx[i_part][pn_entity_extended_old[i_part]+i_entity] + ln_entity;

      /* ln_to_gn */
      pfull_entity_extended_ln_to_gn      [i_part][pn_entity_extended_old[i_part]+i_entity] = pentity_extended_ln_to_gn[i_part][i_entity];
    }

    for(int i = 0; i < pentity_extended_to_pentity_idx[i_part][pn_entity_extended[i_part]]/3; ++i) {
      int idx_write = pfull_entity_extended_to_pentity_idx[i_part][pn_entity_extended_old[i_part]]/3 + i;
      pfull_entity_extended_to_pentity_triplet  [i_part][3*idx_write  ] = pentity_extended_to_pentity_triplet  [i_part][3*i  ];
      pfull_entity_extended_to_pentity_triplet  [i_part][3*idx_write+1] = pentity_extended_to_pentity_triplet  [i_part][3*i+1];
      pfull_entity_extended_to_pentity_triplet  [i_part][3*idx_write+2] = pentity_extended_to_pentity_triplet  [i_part][3*i+2];
      pfull_entity_extended_to_pentity_interface[i_part][  idx_write  ] = pentity_extended_to_pentity_interface[i_part][  i  ];
    }

  }
}



static
void
_compute_offset
(
  PDM_part_extension_t *part_ext
)
{
  int **pn_cell       = NULL;
  int **pn_face       = NULL;
  int **pn_edge       = NULL;
  int **pn_vtx        = NULL;
  int **pn_edge_group = NULL;
  int **pn_face_group = NULL;
  PDM_malloc(pn_cell      , part_ext->n_domain, int *);
  PDM_malloc(pn_face      , part_ext->n_domain, int *);
  PDM_malloc(pn_edge      , part_ext->n_domain, int *);
  PDM_malloc(pn_vtx       , part_ext->n_domain, int *);
  PDM_malloc(pn_edge_group, part_ext->n_domain, int *);
  PDM_malloc(pn_face_group, part_ext->n_domain, int *);

  PDM_g_num_t ***cell_ln_to_gn       = NULL;
  PDM_g_num_t ***face_ln_to_gn       = NULL;
  PDM_g_num_t ***edge_ln_to_gn       = NULL;
  PDM_g_num_t ***vtx_ln_to_gn        = NULL;
  PDM_g_num_t ***face_group_ln_to_gn = NULL;
  PDM_g_num_t ***edge_group_ln_to_gn = NULL;
  PDM_malloc(cell_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(face_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(edge_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(vtx_ln_to_gn       , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(edge_group_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(face_group_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_malloc(pn_cell      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_face      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_edge      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_vtx       [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_edge_group[i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_face_group[i_domain], part_ext->n_part[i_domain], int);

    PDM_malloc(cell_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(face_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(edge_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(vtx_ln_to_gn       [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(face_group_ln_to_gn[i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(edge_group_ln_to_gn[i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      pn_cell            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_cell;
      pn_face            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pn_edge            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_vtx             [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge_group      [i_domain][i_part] = 0;
      pn_face_group      [i_domain][i_part] = 0;
      if (part_ext->parts[i_domain][i_part].n_edge_group > 0) {
        pn_edge_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_bound_idx[part_ext->parts[i_domain][i_part].n_edge_group];
      }
      if (part_ext->parts[i_domain][i_part].n_face_group > 0) {
        pn_face_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_idx[part_ext->parts[i_domain][i_part].n_face_group];
      }

      cell_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
      face_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      edge_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      vtx_ln_to_gn       [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      edge_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_bound_ln_to_gn;
      face_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_ln_to_gn;

    }
  }

  assert(part_ext->shift_by_domain_cell == NULL);
  part_ext->shift_by_domain_cell = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                         part_ext->n_part,
                                                                         pn_cell,
                                                                         cell_ln_to_gn,
                                                                         part_ext->comm);

  assert(part_ext->shift_by_domain_face == NULL);
  part_ext->shift_by_domain_face = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                         part_ext->n_part,
                                                                         pn_face,
                                                                         face_ln_to_gn,
                                                                         part_ext->comm);

  assert(part_ext->shift_by_domain_edge == NULL);
  part_ext->shift_by_domain_edge = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                         part_ext->n_part,
                                                                         pn_edge,
                                                                         edge_ln_to_gn,
                                                                         part_ext->comm);

  assert(part_ext->shift_by_domain_vtx == NULL);
  part_ext->shift_by_domain_vtx = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                        part_ext->n_part,
                                                                        pn_vtx,
                                                                        vtx_ln_to_gn,
                                                                        part_ext->comm);

  // Not used but keep for memory / debug
  // assert(part_ext->shift_by_domain_edge_group == NULL);
  // part_ext->shift_by_domain_edge_group = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
  //                                                                              part_ext->n_part,
  //                                                                              pn_edge_group,
  //                                                                              edge_group_ln_to_gn,
  //                                                                              part_ext->comm);

  // assert(part_ext->shift_by_domain_face_group == NULL);
  // part_ext->shift_by_domain_face_group = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
  //                                                                              part_ext->n_part,
  //                                                                              pn_face_group,
  //                                                                              face_group_ln_to_gn,
  //                                                                              part_ext->comm);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(pn_cell      [i_domain]);
    PDM_free(pn_face      [i_domain]);
    PDM_free(pn_edge      [i_domain]);
    PDM_free(pn_vtx       [i_domain]);
    PDM_free(pn_edge_group[i_domain]);
    PDM_free(pn_face_group[i_domain]);

    PDM_free(cell_ln_to_gn      [i_domain]);
    PDM_free(face_ln_to_gn      [i_domain]);
    PDM_free(edge_ln_to_gn      [i_domain]);
    PDM_free(vtx_ln_to_gn       [i_domain]);
    PDM_free(edge_group_ln_to_gn[i_domain]);
    PDM_free(face_group_ln_to_gn[i_domain]);
  }

  PDM_free(pn_cell      );
  PDM_free(pn_face      );
  PDM_free(pn_edge      );
  PDM_free(pn_vtx       );
  PDM_free(pn_edge_group);
  PDM_free(pn_face_group);

  PDM_free(cell_ln_to_gn      );
  PDM_free(face_ln_to_gn      );
  PDM_free(edge_ln_to_gn      );
  PDM_free(vtx_ln_to_gn       );
  PDM_free(edge_group_ln_to_gn);
  PDM_free(face_group_ln_to_gn);
}



static
void
_offset_parts_by_domain
(
  PDM_part_extension_t *part_ext,
  int                   sens
)
{
  int **pn_cell       = NULL;
  int **pn_face       = NULL;
  int **pn_edge       = NULL;
  int **pn_vtx        = NULL;
  int **pn_edge_group = NULL;
  int **pn_face_group = NULL;
  PDM_malloc(pn_cell      , part_ext->n_domain, int *);
  PDM_malloc(pn_face      , part_ext->n_domain, int *);
  PDM_malloc(pn_edge      , part_ext->n_domain, int *);
  PDM_malloc(pn_vtx       , part_ext->n_domain, int *);
  PDM_malloc(pn_edge_group, part_ext->n_domain, int *);
  PDM_malloc(pn_face_group, part_ext->n_domain, int *);

  PDM_g_num_t ***cell_ln_to_gn       = NULL;
  PDM_g_num_t ***face_ln_to_gn       = NULL;
  PDM_g_num_t ***edge_ln_to_gn       = NULL;
  PDM_g_num_t ***vtx_ln_to_gn        = NULL;
  PDM_g_num_t ***edge_group_ln_to_gn = NULL;
  PDM_g_num_t ***face_group_ln_to_gn = NULL;
  PDM_malloc(cell_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(face_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(edge_ln_to_gn      , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(vtx_ln_to_gn       , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(edge_group_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(face_group_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_malloc(pn_cell      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_face      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_edge      [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_vtx       [i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_edge_group[i_domain], part_ext->n_part[i_domain], int);
    PDM_malloc(pn_face_group[i_domain], part_ext->n_part[i_domain], int);

    PDM_malloc(cell_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(face_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(edge_ln_to_gn      [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(vtx_ln_to_gn       [i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(edge_group_ln_to_gn[i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);
    PDM_malloc(face_group_ln_to_gn[i_domain], part_ext->n_part[i_domain], PDM_g_num_t *);

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      pn_cell            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_cell;
      pn_face            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pn_edge            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_vtx             [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge_group      [i_domain][i_part] = 0;
      pn_face_group      [i_domain][i_part] = 0;
      if (part_ext->parts[i_domain][i_part].n_edge_group > 0) {
        pn_edge_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_bound_idx[part_ext->parts[i_domain][i_part].n_edge_group];
      }
      if (part_ext->parts[i_domain][i_part].n_face_group > 0) {
        pn_face_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_idx[part_ext->parts[i_domain][i_part].n_face_group];
      }

      cell_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
      face_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      edge_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      vtx_ln_to_gn       [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      edge_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_bound_ln_to_gn;
      face_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_ln_to_gn;

    }
  }

  PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                part_ext->n_part,
                                pn_cell,
                                cell_ln_to_gn,
                                part_ext->shift_by_domain_cell,
                                sens);


  PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                part_ext->n_part,
                                pn_face,
                                face_ln_to_gn,
                                part_ext->shift_by_domain_face,
                                sens);

  PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                part_ext->n_part,
                                pn_edge,
                                edge_ln_to_gn,
                                part_ext->shift_by_domain_edge,
                                sens);

  PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                part_ext->n_part,
                                pn_vtx,
                                vtx_ln_to_gn,
                                part_ext->shift_by_domain_vtx,
                                sens);

  // > Group gnum is only exchanged so shift seems to be unnecessary
  // PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
  //                               part_ext->n_part,
  //                               pn_edge_group,
  //                               edge_group_ln_to_gn,
  //                               part_ext->shift_by_domain_edge_group,
  //                               sens);

  // PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
  //                               part_ext->n_part,
  //                               pn_face_group,
  //                               face_group_ln_to_gn,
  //                               part_ext->shift_by_domain_face_group,
  //                               sens);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(pn_cell      [i_domain]);
    PDM_free(pn_face      [i_domain]);
    PDM_free(pn_edge      [i_domain]);
    PDM_free(pn_vtx       [i_domain]);
    PDM_free(pn_edge_group[i_domain]);
    PDM_free(pn_face_group[i_domain]);

    PDM_free(cell_ln_to_gn      [i_domain]);
    PDM_free(face_ln_to_gn      [i_domain]);
    PDM_free(edge_ln_to_gn      [i_domain]);
    PDM_free(vtx_ln_to_gn       [i_domain]);
    PDM_free(edge_group_ln_to_gn[i_domain]);
    PDM_free(face_group_ln_to_gn[i_domain]);
  }

  PDM_free(pn_cell      );
  PDM_free(pn_face      );
  PDM_free(pn_edge      );
  PDM_free(pn_vtx       );
  PDM_free(pn_edge_group);
  PDM_free(pn_face_group);

  PDM_free(cell_ln_to_gn      );
  PDM_free(face_ln_to_gn      );
  PDM_free(edge_ln_to_gn      );
  PDM_free(vtx_ln_to_gn       );
  PDM_free(edge_group_ln_to_gn);
  PDM_free(face_group_ln_to_gn);
}

static
void
_compute_other_part_domain_interface
(
 PDM_part_extension_t *part_ext
)
{
  int have_edge = 0;
  int have_face = 0;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      if(part_ext->parts[i_domain][i_part].n_edge > 0 &&
         part_ext->parts[i_domain][i_part].edge_ln_to_gn != NULL) {
        have_edge = 1;
      }
      if(part_ext->parts[i_domain][i_part].n_face > 0 &&
         part_ext->parts[i_domain][i_part].face_ln_to_gn != NULL) {
        have_face = 1;
      }
    }
  }
  part_ext->have_edge = have_edge;
  part_ext->have_face = have_face;

  if(part_ext->pdi == NULL) {
    return;
  }

  int is_describe_vtx  = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_VTX );
  int is_describe_edge = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_EDGE);
  int is_describe_face = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_FACE);

  int is_describe_vtx_l  = is_describe_vtx;
  int is_describe_edge_l = is_describe_edge;
  int is_describe_face_l = is_describe_face;
  PDM_MPI_Allreduce(&is_describe_vtx_l , &is_describe_vtx , 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_edge_l, &is_describe_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_face_l, &is_describe_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);


  // En gros noeud centré avec toutes les connectivités
  if(is_describe_vtx == 1 &&
     (is_describe_edge == 0 || is_describe_face == 0) &&
     (have_edge       == 1 || (have_face == 1 && part_ext->dim==3))) {

    // Rebuild domaine_interface in distributed frame
    // PDM_domain_interface* dintrf = PDM_part_domain_interface_to_domain_interface()


    int          **pn_vtx         = NULL;
    int          **pn_edge        = NULL;
    int          **pn_face        = NULL;
    PDM_g_num_t ***pvtx_ln_to_gn  = NULL;
    PDM_g_num_t ***pedge_ln_to_gn = NULL;
    PDM_g_num_t ***pface_ln_to_gn = NULL;
    int         ***pedge_vtx_idx  = NULL;
    int         ***pedge_vtx      = NULL;
    int         ***pface_edge_idx = NULL;
    int         ***pface_edge     = NULL;
    int         ***pface_vtx_idx  = NULL;
    int         ***pface_vtx      = NULL;
    PDM_malloc(pn_vtx        , part_ext->n_domain, int          *);
    PDM_malloc(pn_edge       , part_ext->n_domain, int          *);
    PDM_malloc(pn_face       , part_ext->n_domain, int          *);
    PDM_malloc(pvtx_ln_to_gn , part_ext->n_domain, PDM_g_num_t **);
    PDM_malloc(pedge_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);
    PDM_malloc(pface_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);
    PDM_malloc(pedge_vtx_idx , part_ext->n_domain, int         **);
    PDM_malloc(pedge_vtx     , part_ext->n_domain, int         **);
    PDM_malloc(pface_edge_idx, part_ext->n_domain, int         **);
    PDM_malloc(pface_edge    , part_ext->n_domain, int         **);
    PDM_malloc(pface_vtx_idx , part_ext->n_domain, int         **);
    PDM_malloc(pface_vtx     , part_ext->n_domain, int         **);

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_malloc(pn_vtx          [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pn_edge         [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pn_face         [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pvtx_ln_to_gn   [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pedge_ln_to_gn  [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pface_ln_to_gn  [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pedge_vtx_idx   [i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pedge_vtx       [i_domain], part_ext->n_domain, int         *);
      if(have_edge == 1) {
        PDM_malloc(pface_edge_idx[i_domain], part_ext->n_domain, int         *);
        PDM_malloc(pface_edge    [i_domain], part_ext->n_domain, int         *);
      }
      else {
        PDM_malloc(pface_vtx_idx [i_domain], part_ext->n_domain, int         *);
        PDM_malloc(pface_vtx     [i_domain], part_ext->n_domain, int         *);
      }

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_vtx          [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        pn_edge         [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
        pn_face         [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        pvtx_ln_to_gn   [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        pedge_ln_to_gn  [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
        pface_ln_to_gn  [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
        PDM_malloc(pedge_vtx[i_domain][i_part], 2 * pn_edge[i_domain][i_part], int);
        if(have_edge == 1) {
          pface_edge_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge_idx;
          pface_edge    [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge;
        }
        else {
          pface_vtx_idx [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_vtx_idx;
          pface_vtx     [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_vtx;
          // pedge_vtx     [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_vtx;
        }


        int *_edge_vtx = part_ext->parts[i_domain][i_part].edge_vtx;
        for(int i_edge = 0; i_edge < pn_edge[i_domain][i_part]; ++i_edge) {
          pedge_vtx     [i_domain][i_part][2*i_edge  ] =  _edge_vtx[2*i_edge  ];
          pedge_vtx     [i_domain][i_part][2*i_edge+1] = -_edge_vtx[2*i_edge+1];
          // printf("i_edge = %i (%i)- i_vtx1 = %i | i_vtx2 = %i \n", i_edge, (int) pedge_ln_to_gn[i_domain][i_part][i_edge], _edge_vtx[2*i_edge  ], -_edge_vtx[2*i_edge+1]);
        }

        int _nedge = pn_edge[i_domain][i_part];
        PDM_malloc(pedge_vtx_idx[i_domain][i_part], _nedge+1, int);

        pedge_vtx_idx [i_domain][i_part][0] = 0;
        for(int i = 0; i < _nedge; ++i) {
          pedge_vtx_idx [i_domain][i_part][i+1] = pedge_vtx_idx [i_domain][i_part][i] + 2;
        }
      }
    }

    if(is_describe_edge == 0) {

      // Translate
      PDM_part_domain_interface_add(part_ext->pdi,
                                    PDM_BOUND_TYPE_VTX,
                                    PDM_BOUND_TYPE_EDGE,
                                    part_ext->n_part,
                                    pn_vtx,
                                    pvtx_ln_to_gn,
                                    pn_edge,
                                    pedge_ln_to_gn,
                                    pedge_vtx_idx,
                                    pedge_vtx,
                                    1); // Connectivity_is_signed
    }


    if(part_ext->dim==3 && have_face == 1 &&  is_describe_face == 0) {
      // Translate
      if(have_edge == 1) {
        PDM_part_domain_interface_add(part_ext->pdi,
                                      PDM_BOUND_TYPE_EDGE,
                                      PDM_BOUND_TYPE_FACE,
                                      part_ext->n_part,
                                      pn_edge,
                                      pedge_ln_to_gn,
                                      pn_face,
                                      pface_ln_to_gn,
                                      pface_edge_idx,
                                      pface_edge,
                                      1);// Connectivity_is_signed
      } else {
        PDM_part_domain_interface_add(part_ext->pdi,
                                      PDM_BOUND_TYPE_VTX,
                                      PDM_BOUND_TYPE_FACE,
                                      part_ext->n_part,
                                      pn_vtx,
                                      pvtx_ln_to_gn,
                                      pn_face,
                                      pface_ln_to_gn,
                                      pface_vtx_idx,
                                      pface_vtx,
                                      0);// Connectivity_is_signed
      }
    }


    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        PDM_free(pedge_vtx_idx [i_domain][i_part]);
        PDM_free(pedge_vtx     [i_domain][i_part]);
      }
      PDM_free(pn_vtx        [i_domain]);
      PDM_free(pn_edge       [i_domain]);
      PDM_free(pn_face       [i_domain]);
      PDM_free(pvtx_ln_to_gn [i_domain]);
      PDM_free(pedge_ln_to_gn[i_domain]);
      PDM_free(pface_ln_to_gn[i_domain]);
      PDM_free(pedge_vtx_idx [i_domain]);
      PDM_free(pedge_vtx     [i_domain]);
      if(have_edge == 1) {
        PDM_free(pface_edge_idx[i_domain]);
        PDM_free(pface_edge    [i_domain]);
      }
      else {
        PDM_free(pface_vtx_idx[i_domain]);
        PDM_free(pface_vtx    [i_domain]);
      }
    }
    PDM_free(pn_vtx        );
    PDM_free(pn_edge       );
    PDM_free(pn_face       );
    PDM_free(pvtx_ln_to_gn );
    PDM_free(pedge_ln_to_gn);
    PDM_free(pface_ln_to_gn);
    PDM_free(pedge_vtx_idx );
    PDM_free(pedge_vtx     );
    PDM_free(pface_edge_idx);
    PDM_free(pface_edge    );
    PDM_free(pface_vtx_idx );
    PDM_free(pface_vtx     );

  } else if (is_describe_face == 1) {

    // assert(is_describe_vtx == 0);

    // Faire la méthode de Julien face2vtx

  }
}

// UNUSED BUT TO KEEP
// static
// void
// _build_part_extension_graph_to_new
// (
//  PDM_part_extension_t *part_ext,
//  int                  *pn_init_entity,
//  int                  *pn_entity,
//  PDM_g_num_t         **pentity_gnum,
//  int                 **pentity_ancstr_strd,
//  PDM_g_num_t         **pentity_ancstr,
//  int                 **pentity_path_itrf_strd,
//  int                 **pentity_path_itrf,
//  int                ***out_pentity_to_entity_idx,
//  int                ***out_pentity_to_entity_trplt,
//  int                ***out_pentity_to_entity_path_itrf_idx,
//  int                ***out_pentity_to_entity_path_itrf
// )
// {
//   int debug = 0;
//   if (debug==1) {
//     log_trace("\n\n _build_part_extension_graph_to_new\n\n");
//   }


//   /**
//    * Build part_to_part between new entities all ancestors
//    */
//   int **pentity_ancstr_idx = NULL;
//   PDM_malloc(pentity_ancstr_idx, part_ext->ln_part_tot, int *);
//   for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
//     pentity_ancstr_idx[i_part] = PDM_array_new_idx_from_sizes_int(pentity_ancstr_strd[i_part], pn_entity[i_part]);
//   }

//   PDM_part_to_part_t *ptp = PDM_part_to_part_create(
//     (const PDM_g_num_t **) pentity_gnum,
//                            pn_entity,
//                            part_ext->ln_part_tot,
//     (const PDM_g_num_t **) pentity_gnum,
//                            pn_init_entity,
//                            part_ext->ln_part_tot,
//     (const int         **) pentity_ancstr_idx,
//     (const PDM_g_num_t **) pentity_ancstr,
//                            part_ext->comm);

//   int  *n_ref_lnum2 = NULL;
//   int **ref_lnum2   = NULL;
//   PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

//   int         **gnum1_come_from_idx = NULL;
//   PDM_g_num_t **gnum1_come_from     = NULL;
//   PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);


//   /**
//    * Build triplet and send it with path itrf
//    */
//   int i_rank;
//   PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
//   int **pentity_trplt = NULL;
//   PDM_malloc(pentity_trplt, part_ext->ln_part_tot, int *);
//   for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
//     PDM_malloc(pentity_trplt[i_part], 3*pn_entity[i_part], int);
//     for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
//       pentity_trplt[i_part][3*i_entity  ] = i_rank;
//       pentity_trplt[i_part][3*i_entity+1] = i_part;
//       pentity_trplt[i_part][3*i_entity+2] = i_entity;
//     }
//   }

//   int request_exch_ancstr = 0;
//   int **rcvd_trplt      = NULL;
//   PDM_part_to_part_iexch( ptp,
//                           PDM_MPI_COMM_KIND_P2P,
//                           PDM_STRIDE_CST_INTERLACED,
//                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
//                           3,
//                           sizeof(int),
//         (const int  **)   NULL,
//         (const void **)   pentity_trplt,
//                           NULL,
//               (void ***) &rcvd_trplt,
//                          &request_exch_ancstr);
//   PDM_part_to_part_iexch_wait(ptp, request_exch_ancstr);

//   int **rcvd_path_itrf_strd = NULL;
//   int **rcvd_path_itrf      = NULL;
//   PDM_part_to_part_iexch( ptp,
//                           PDM_MPI_COMM_KIND_P2P,
//                           PDM_STRIDE_VAR_INTERLACED,
//                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
//                           1,
//                           1*sizeof(int),
//         (const int  **)   pentity_path_itrf_strd,
//         (const void **)   pentity_path_itrf,
//                          &rcvd_path_itrf_strd,
//               (void ***) &rcvd_path_itrf,
//                          &request_exch_ancstr);
//   PDM_part_to_part_iexch_wait(ptp, request_exch_ancstr);

//   if (debug==1) {
//     for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
//       log_trace("IPART = %d \n", i_part);
//       int i_read_path = 0;
//       for (int i_ref=0; i_ref<n_ref_lnum2[i_part]; ++i_ref) {
//         int ref_lnum = ref_lnum2[i_part][i_ref];
//         log_trace("gnum "PDM_FMT_G_NUM" \n", pentity_gnum[i_part][ref_lnum-1]);
//         int i_beg_cf = gnum1_come_from_idx[i_part][ref_lnum-1];
//         int i_end_cf = gnum1_come_from_idx[i_part][ref_lnum  ];
//         log_trace("i_beg_cf = %d; i_end_cf = %d \n", i_beg_cf, i_end_cf);
//         for (int i_cf=i_beg_cf; i_cf<i_end_cf; ++i_cf) {
//           log_trace("i_cf = %d; \n", i_cf);
//           log_trace("\t connected to %d %d %d (%d/%d) with path ",
//                             rcvd_trplt[i_part][3*i_cf  ],
//                             rcvd_trplt[i_part][3*i_cf+1],
//                             rcvd_trplt[i_part][3*i_cf+2],
//                             i_cf    -i_beg_cf,
//                             i_end_cf-i_beg_cf);
//           for (int i_path=0; i_path<rcvd_path_itrf_strd[i_part][i_cf]; ++i_path) {
//             log_trace("%d ", rcvd_path_itrf[i_part][i_read_path++]);
//           }
//           log_trace("\n");
//         }
//       }
//     }
//   }



//   /**
//    * Set result
//    */
//   int **rcvd_trplt_idx     = NULL;
//   int **rcvd_path_itrf_idx = NULL;
//   PDM_malloc(rcvd_trplt_idx    , part_ext->ln_part_tot, int *);
//   PDM_malloc(rcvd_path_itrf_idx, part_ext->ln_part_tot, int *);
//   for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

//     // > Get strd for all init entities
//     int *_rcvd_trplt_strd     = PDM_array_zeros_int(pn_init_entity[i_part]);
//     for (int i_ref=0; i_ref<n_ref_lnum2[i_part]; ++i_ref) {
//       int ref_lnum    = ref_lnum2[i_part][i_ref];
//       int n_come_from = gnum1_come_from_idx[i_part][i_ref+1]-gnum1_come_from_idx[i_part][i_ref];
//       _rcvd_trplt_strd[ref_lnum-1] = n_come_from;
//     }
//     int n_rcvd = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]];
//     rcvd_trplt_idx    [i_part] = PDM_array_new_idx_from_sizes_int(_rcvd_trplt_strd    , pn_init_entity[i_part]);
//     rcvd_path_itrf_idx[i_part] = PDM_array_new_idx_from_sizes_int(rcvd_path_itrf_strd[i_part], n_rcvd);
    
//     PDM_free(_rcvd_trplt_strd);
//     PDM_free(rcvd_path_itrf_strd[i_part]);
//     PDM_free(pentity_trplt[i_part]);
//     PDM_free(pentity_ancstr_idx[i_part]);
//   }
//   PDM_free(rcvd_path_itrf_strd);
//   PDM_free(pentity_trplt);
//   PDM_free(pentity_ancstr_idx);
//   PDM_part_to_part_free(ptp);

//   *out_pentity_to_entity_idx           = rcvd_trplt_idx;
//   *out_pentity_to_entity_trplt         = rcvd_trplt;
//   *out_pentity_to_entity_path_itrf_idx = rcvd_path_itrf_idx;
//   *out_pentity_to_entity_path_itrf     = rcvd_path_itrf;
// }


static
void
_build_part_extension_graph_to_old
(
  PDM_part_extension_t  *part_ext,
  int                   *init_n_entity,
  PDM_g_num_t          **init_entity_gnum,
  int                   *init_n_entity_group,
  int                  **init_entity_group_tag,
  PDM_g_num_t          **init_entity_group_gnum,
  int                   *new_n_entity,
  PDM_g_num_t          **new_entity_gnum,
  int                  **new_entity_ancstr_strd,
  PDM_g_num_t          **new_entity_ancstr,
  PDM_g_num_t           *shift_entity,
  // int                  **new_entity_path_itrf_strd,
  // int                  **new_entity_path_itrf,
  PDM_g_num_t         ***out_pentity_to_entity_ancstr,
  int                 ***out_pentity_to_entity_trplt,
  int                 ***out_pentity_group_idx,
  int                 ***out_pentity_group,
  PDM_g_num_t         ***out_pentity_group_gnum,
  PDM_mesh_entities_t    mesh_entity
)
{
  int debug = 0;
  if (debug==1) {
    log_trace("\n\n _build_part_extension_graph_to_old\n\n");
  }

  int have_group = init_entity_group_tag!=NULL;

  /**
   * Build part_to_part between new entities all ancestors
   * assuming ancestor = gnum when missing
   */
  int         **l_entity_ancstr_strd = NULL;
  int         **l_entity_ancstr_idx  = NULL;
  PDM_g_num_t **l_entity_ancstr      = NULL;
  PDM_malloc(l_entity_ancstr_strd, part_ext->ln_part_tot, int         *);
  PDM_malloc(l_entity_ancstr_idx , part_ext->ln_part_tot, int         *);
  PDM_malloc(l_entity_ancstr     , part_ext->ln_part_tot, PDM_g_num_t *);
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_malloc(l_entity_ancstr[i_part], new_n_entity[i_part], PDM_g_num_t);
    l_entity_ancstr_strd[i_part] = PDM_array_copy_int(new_entity_ancstr_strd[i_part], new_n_entity[i_part]);

    int i_read_ancstr = 0;
    for (int i_entity=0; i_entity<new_n_entity[i_part]; ++i_entity) {
      if (l_entity_ancstr_strd[i_part][i_entity]==0) {
        l_entity_ancstr     [i_part][i_entity] = new_entity_gnum[i_part][i_entity];
        l_entity_ancstr_strd[i_part][i_entity] = 1;
      }
      else {
        l_entity_ancstr[i_part][i_entity] = new_entity_ancstr[i_part][i_read_ancstr++];
      }
    }
    l_entity_ancstr_idx[i_part] = PDM_array_new_idx_from_sizes_int(l_entity_ancstr_strd[i_part], new_n_entity[i_part]);
  
    if (debug==1) {
      log_trace("IPART = %d \n", i_part);
      log_trace("\t init_n_entity = %d \n", init_n_entity[i_part]);
      PDM_log_trace_array_long(init_entity_gnum[i_part], init_n_entity[i_part], "\t init_entity_gnum ::");
      
      log_trace("\t new_n_entity = %d \n", new_n_entity[i_part]);
      PDM_log_trace_array_long(new_entity_gnum     [i_part], new_n_entity[i_part]  , "\t new_entity_gnum      ::");
      PDM_log_trace_array_int (l_entity_ancstr_strd[i_part], new_n_entity[i_part]  , "\t l_entity_ancstr_strd ::");
      PDM_log_trace_array_int (l_entity_ancstr_idx [i_part], new_n_entity[i_part]+1, "\t l_entity_ancstr_idx  ::");
      PDM_log_trace_array_long(l_entity_ancstr     [i_part], new_n_entity[i_part]  , "\t l_entity_ancstr      ::");

    }

    PDM_free(l_entity_ancstr_strd[i_part]);
  }
  PDM_free(l_entity_ancstr_strd);

  PDM_part_to_part_t *ptp = PDM_part_to_part_create(
    (const PDM_g_num_t **) new_entity_gnum, 
                           new_n_entity,
                           part_ext->ln_part_tot,
    (const PDM_g_num_t **) init_entity_gnum, 
                           init_n_entity,
                           part_ext->ln_part_tot,
    (const int         **) l_entity_ancstr_idx,
    (const PDM_g_num_t **) l_entity_ancstr,
                           part_ext->comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);


  /**
   * Build triplet and send it with path itrf
   */
  int i_rank;
  PDM_MPI_Comm_rank(part_ext->comm, &i_rank);

  int **pentity_trplt = NULL;
  PDM_malloc(pentity_trplt, part_ext->ln_part_tot, int *);
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_malloc(pentity_trplt[i_part], 3*init_n_entity[i_part], int);
    for (int i_ref=0; i_ref<n_ref_lnum2[i_part]; ++i_ref) {
      int ref_lnum = ref_lnum2[i_part][i_ref];
      pentity_trplt[i_part][3*(ref_lnum-1)  ] = i_rank;
      pentity_trplt[i_part][3*(ref_lnum-1)+1] = i_part;
      pentity_trplt[i_part][3*(ref_lnum-1)+2] = ref_lnum-1;
    }
  }

  int request_exch = 0;
  int **rcvd_trplt      = NULL;
  PDM_part_to_part_reverse_iexch( ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_CST_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                          3,
                          sizeof(int),
        (const int  **)   NULL,
        (const void **)   pentity_trplt,
                          NULL,
              (void ***) &rcvd_trplt,
                         &request_exch);
  PDM_part_to_part_reverse_iexch_wait(ptp, request_exch);

  int         **rcvd_group_tag  = NULL;
  PDM_g_num_t **rcvd_group_gnum = NULL;
  if (have_group) {
    PDM_part_to_part_reverse_iexch( ptp,
                            PDM_MPI_COMM_KIND_P2P,
                            PDM_STRIDE_CST_INTERLACED,
                            PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                            1,
                            sizeof(int),
          (const int  **)   NULL,
          (const void **)   init_entity_group_tag,
                            NULL,
                (void ***) &rcvd_group_tag,
                           &request_exch);
    PDM_part_to_part_reverse_iexch_wait(ptp, request_exch);

    PDM_part_to_part_reverse_iexch( ptp,
                            PDM_MPI_COMM_KIND_P2P,
                            PDM_STRIDE_CST_INTERLACED,
                            PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                            1,
                            sizeof(PDM_g_num_t),
          (const int  **)   NULL,
          (const void **)   init_entity_group_gnum,
                            NULL,
                (void ***) &rcvd_group_gnum,
                           &request_exch);
    PDM_part_to_part_reverse_iexch_wait(ptp, request_exch); 
  }


  if (debug==1) {
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      log_trace("IPART = %d \n", i_part);
      // int i_read_path = 0;
      for (int i_entity=0; i_entity<new_n_entity[i_part]; ++i_entity) {
        log_trace("gnum "PDM_FMT_G_NUM" \n", new_entity_gnum[i_part][i_entity]);
        log_trace("\t from gnum "PDM_FMT_G_NUM" \n", l_entity_ancstr[i_part][i_entity]);
        // log_trace("   via path ");
        // for (int i_itrf=0; i_itrf<new_entity_path_itrf_strd[i_part][i_entity]; ++i_itrf) {
        //   log_trace("%d ", new_entity_path_itrf[i_part][i_read_path]);
        //   i_read_path++;
        // }
        log_trace("\n");
        log_trace("\t received triplet %d %d %d\n", 
                  rcvd_trplt[i_part][3*i_entity  ], 
                  rcvd_trplt[i_part][3*i_entity+1], 
                  rcvd_trplt[i_part][3*i_entity+2]);
        if (have_group) {
          log_trace("\t get tag %d (with gnum "PDM_FMT_G_NUM"\n", 
              rcvd_group_tag [i_part][i_entity],
              rcvd_group_gnum[i_part][i_entity]);
        }
      }
    }
  }

  
  /**
   * Transform tag into group
   */
  int         **entity_group_idx  = NULL;
  int         **entity_group      = NULL;
  PDM_g_num_t **entity_group_gnum = NULL;
  if (have_group) {
    PDM_malloc(entity_group_idx , part_ext->ln_part_tot, int         *);
    PDM_malloc(entity_group     , part_ext->ln_part_tot, int         *);
    PDM_malloc(entity_group_gnum, part_ext->ln_part_tot, PDM_g_num_t *);

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      int *entity_group_n = NULL;
      PDM_calloc(entity_group_n, init_n_entity_group[i_part], int);

      // > Count number of entities in group
      for (int i_entity=0; i_entity<new_n_entity[i_part]; ++i_entity) {
        int group_tag = rcvd_group_tag[i_part][i_entity];
        if (group_tag>0) {
          entity_group_n[group_tag-1]++;
        }
      }

      // > Set idx and alloc
      entity_group_idx[i_part] = PDM_array_new_idx_from_sizes_int(entity_group_n, init_n_entity_group[i_part]);
      PDM_array_reset_int(entity_group_n, init_n_entity_group[i_part], 0);

      int n_entity_in_group = entity_group_idx[i_part][init_n_entity_group[i_part]];
      PDM_malloc(entity_group     [i_part], n_entity_in_group, int        );
      PDM_malloc(entity_group_gnum[i_part], n_entity_in_group, PDM_g_num_t);

      // > Fill
      for (int i_entity=0; i_entity<new_n_entity[i_part]; ++i_entity) {
        int group_tag = rcvd_group_tag[i_part][i_entity];
        if (group_tag>0) {
          int i_write = entity_group_idx[i_part][group_tag-1] + entity_group_n[group_tag-1];
          entity_group     [i_part][i_write] = init_n_entity[i_part]+i_entity+1;
          entity_group_gnum[i_part][i_write] = rcvd_group_gnum[i_part][i_entity];
          entity_group_n[group_tag-1]++;
        }
      }

      if (debug == 1) {
        log_trace("n_group = %d\n", init_n_entity_group[i_part]);
        int size_group = entity_group_idx[i_part][init_n_entity_group[i_part]];
        PDM_log_trace_connectivity_int(entity_group_idx [i_part],
                                       entity_group     [i_part],
                                       init_n_entity_group[i_part], "entity_group");
        PDM_log_trace_array_long(entity_group_gnum[i_part], size_group                 , "entity_group_gnum ::");
      }

      free(entity_group_n);
      free(rcvd_group_tag [i_part]);
      free(rcvd_group_gnum[i_part]);
    }
    free(rcvd_group_tag);
    free(rcvd_group_gnum);
  }

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(l_entity_ancstr_idx[i_part]);
    PDM_free(pentity_trplt      [i_part]);
  }
  PDM_free(l_entity_ancstr_idx);
  PDM_free(pentity_trplt);
  PDM_part_to_part_free(ptp);

  /*
   * We reapply shift after part_to_part to avoid merge same gnum in part_to_part (but from different domain)
   */
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    for (int i_entity=0; i_entity<new_n_entity[i_part]; ++i_entity) {
      int shift = shift_entity[part_ext->lpart_to_dom[i_part]];
      l_entity_ancstr[i_part][i_entity] -= shift;
    }
  }

  *out_pentity_to_entity_trplt  = rcvd_trplt;
  *out_pentity_to_entity_ancstr = l_entity_ancstr;
  // > Set ownership
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
      part_ext->ownership_border_ln_to_gn_ancstr[mesh_entity][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
    }
  }

  if (have_group) {
    *out_pentity_group_idx       = entity_group_idx;
    *out_pentity_group           = entity_group;
    *out_pentity_group_gnum      = entity_group_gnum;
    // > Set ownership
    for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
      for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
        part_ext->ownership_border_group[mesh_entity][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
      }
    }
  }
}



static
void
_get_block_data_base_on_part_ext
(
 PDM_part_extension_t  *part_ext,
 int                    dn_entity,
 PDM_g_num_t           *dentity_gnum,
 int                   *dentity_ancstr_strd,
 PDM_g_num_t           *dentity_ancstr,
 int                   *dentity_path_itrf_strd,
 int                   *dentity_path_itrf,
 int                   *pn_entity,
 PDM_g_num_t          **pentity_gnum,
 int                 ***out_pentity_ancstr_strd,
 // int                 ***out_pentity_ancstr_idx,
 PDM_g_num_t         ***out_pentity_ancstr,
 // int                 ***out_pentity_path_itrf_strd,
 int                 ***out_pentity_path_itrf_idx,
 int                 ***out_pentity_path_itrf,
 PDM_mesh_entities_t    mesh_entity
)
{
  int debug = 0;

  /**
   * From data_base get ancestor and interface path on new entities on partition
   */

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(dentity_gnum,
                                                                        dn_entity, 
                                                 (const PDM_g_num_t **) pentity_gnum,
                                                 (const int          *) pn_entity,
                                                                        part_ext->ln_part_tot,
                                                                        part_ext->comm);
  
  // > Var stride here cause some entities can be missing from init db
  int         **pentity_ancstr_strd = NULL;
  // int         **pentity_ancstr_idx  = NULL;
  PDM_g_num_t **pentity_ancstr      = NULL;
  PDM_block_to_part_exch(btp,
                         1 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         dentity_ancstr_strd,
                         dentity_ancstr,
                        &pentity_ancstr_strd,
         (void ***)     &pentity_ancstr);

  int **pentity_path_itrf_strd = NULL;
  int **pentity_path_itrf_idx  = NULL;
  int **pentity_path_itrf      = NULL;
  PDM_block_to_part_exch(btp,
                         1 * sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         dentity_path_itrf_strd,
                         dentity_path_itrf,
                        &pentity_path_itrf_strd,
         (void ***)     &pentity_path_itrf);
  PDM_block_to_part_free(btp);


  if (debug==1) {
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      int len_tot_path_itrf = 0;
      int n_ancstr          = 0;
      for(int i_entity = 0; i_entity < pn_entity[i_part]; ++i_entity) {
        len_tot_path_itrf += pentity_path_itrf_strd[i_part][i_entity];
        n_ancstr          += pentity_ancstr_strd   [i_part][i_entity];
      }
      log_trace("i_part = %d\n", i_part);
      log_trace("\t pn_entity = %d\n", pn_entity[i_part]);
      PDM_log_trace_array_long(pentity_gnum          [i_part], pn_entity[i_part], "\t pentity_gnum          [i_part]");
      PDM_log_trace_array_int (pentity_ancstr_strd   [i_part], pn_entity[i_part], "\t pentity_ancstr_strd   [i_part]");
      PDM_log_trace_array_long(pentity_ancstr        [i_part], n_ancstr         , "\t pentity_ancstr        [i_part]");
      PDM_log_trace_array_int (pentity_path_itrf_strd[i_part], pn_entity[i_part], "\t pentity_path_itrf_strd[i_part]");
      PDM_log_trace_array_int (pentity_path_itrf     [i_part], len_tot_path_itrf, "\t pentity_path_itrf     [i_part]");
    }
  }


  // PDM_malloc(pentity_ancstr_idx, part_ext->ln_part_tot, int *);
  PDM_malloc(pentity_path_itrf_idx, part_ext->ln_part_tot, int *);
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    // pentity_ancstr_idx   [i_part] = PDM_array_new_idx_from_sizes_int(pentity_ancstr_strd   [i_part], pn_entity[i_part]);
    pentity_path_itrf_idx[i_part] = PDM_array_new_idx_from_sizes_int(pentity_path_itrf_strd[i_part], pn_entity[i_part]);
    // PDM_free(pentity_ancstr_strd   [i_part]);
    PDM_free(pentity_path_itrf_strd[i_part]);
  }
  // PDM_free(pentity_ancstr_strd);
  PDM_free(pentity_path_itrf_strd);



  *out_pentity_ancstr_strd    = pentity_ancstr_strd;
  // *out_pentity_ancstr_idx     = pentity_ancstr_idx;
  *out_pentity_ancstr         = pentity_ancstr;
  // *out_pentity_path_itrf_strd = pentity_path_itrf_strd;
  *out_pentity_path_itrf_idx  = pentity_path_itrf_idx;
  *out_pentity_path_itrf      = pentity_path_itrf;
  // > Set ownership
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
      part_ext->ownership_border_graph    [mesh_entity][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
      part_ext->ownership_border_path_itrf[mesh_entity][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
    }
  }

}


static
void
_store_extended_entity
(
  PDM_part_extension_t      *part_ext,
  int                       *pn_init_entity1,
  int                       *pn_entity1,
  PDM_g_num_t              **pentity1_gnum,
  int                      **pentity1_pentity2_idx,
  int                      **pentity1_pentity2,
  int                      **out_pn_entity1,
  PDM_g_num_t             ***out_pentity1_gnum,
  int                     ***out_pentity1_pentity2_idx,
  int                     ***out_pentity1_pentity2,
  PDM_mesh_entities_t        mesh_entity,
  PDM_connectivity_type_t    connectivity_type
)
{
  int debug = 0;
  if (debug==1) {
    log_trace("\n\n_store_extended_entity ::\n");
  }

  int has_connectivity = pentity1_pentity2_idx!=NULL;


  int          *_out_pn_entity1            = NULL;
  PDM_g_num_t **_out_pentity1_gnum         = NULL;
  int         **_out_pentity1_pentity2_idx = NULL;
  int         **_out_pentity1_pentity2     = NULL;
  PDM_malloc(_out_pn_entity1   , part_ext->ln_part_tot, int          );
  PDM_malloc(_out_pentity1_gnum, part_ext->ln_part_tot, PDM_g_num_t *);
  if (has_connectivity) {
    PDM_malloc(_out_pentity1_pentity2_idx, part_ext->ln_part_tot, int *);
    PDM_malloc(_out_pentity1_pentity2    , part_ext->ln_part_tot, int *);
  }

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    _out_pn_entity1[i_part] = pn_entity1[i_part]-pn_init_entity1[i_part];
    PDM_malloc(_out_pentity1_gnum[i_part], _out_pn_entity1[i_part], PDM_g_num_t);
    int connectivity_size = 0;
    if (has_connectivity) {
      connectivity_size = pentity1_pentity2_idx[i_part][pn_entity1     [i_part]]
                        - pentity1_pentity2_idx[i_part][pn_init_entity1[i_part]];
      PDM_malloc(_out_pentity1_pentity2_idx[i_part],_out_pn_entity1[i_part]+1, int);
      PDM_malloc(_out_pentity1_pentity2    [i_part], connectivity_size       , int);
      _out_pentity1_pentity2_idx[i_part][0] = 0;
    }

    int i_write      = 0;
    int i_write_conn = 0;
    for (int i_entity1=pn_init_entity1[i_part];
             i_entity1<pn_entity1     [i_part]; ++i_entity1) {
      if (has_connectivity) {
        int i_beg_entity2 = pentity1_pentity2_idx[i_part][i_entity1  ];
        int i_end_entity2 = pentity1_pentity2_idx[i_part][i_entity1+1];
        int n_entity2 = i_end_entity2-i_beg_entity2;
        for (int i_entity2=i_beg_entity2; i_entity2<i_end_entity2; ++i_entity2) {
          _out_pentity1_pentity2[i_part][i_write_conn] = pentity1_pentity2[i_part][i_entity2];
          i_write_conn++;
        }
        _out_pentity1_pentity2_idx[i_part][i_write+1] = _out_pentity1_pentity2_idx[i_part][i_write]+n_entity2;
      }
      _out_pentity1_gnum        [i_part][i_write  ] = pentity1_gnum[i_part][i_entity1];
      i_write++;
    }
    
    if (debug==1) {
      log_trace("\ni_part = %d\n", i_part);
      log_trace("\t _out_pn_entity1 = %d\n", _out_pn_entity1[i_part]);
      PDM_log_trace_array_long(_out_pentity1_gnum        [i_part], _out_pn_entity1[i_part], "\t out_pentity1_gnum         ::");
      if (has_connectivity) {
        PDM_log_trace_array_int (_out_pentity1_pentity2_idx[i_part], _out_pn_entity1[i_part], "\t out_pentity1_pentity2_idx ::");
        PDM_log_trace_array_int (_out_pentity1_pentity2    [i_part], connectivity_size      , "\t out_pentity1_pentity2     ::");
      }
    }
  }


  *out_pn_entity1    = _out_pn_entity1;
  *out_pentity1_gnum = _out_pentity1_gnum;
  // > Set ownership
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
      part_ext->ownership_border_ln_to_gn[mesh_entity][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
    }
  }

  if (has_connectivity) {
    *out_pentity1_pentity2_idx = _out_pentity1_pentity2_idx;
    *out_pentity1_pentity2     = _out_pentity1_pentity2;
    // > Set ownership
    for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
      for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
        part_ext->ownership_border_connectivity[connectivity_type][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
      }
    }
  }


}



static
void
_store_extended_vtx
(
 PDM_part_extension_t *part_ext,
 int                  *pn_init_entity1,
 int                  *pn_entity1,
 PDM_g_num_t         **pentity1_gnum,
 double              **pentity1_coord,
 int                 **out_pn_entity1,
 PDM_g_num_t        ***out_pentity1_gnum,
 double             ***out_pentity1_coord

)
{
  int debug = 0;
  if (debug==1) {
    log_trace("\n\n_store_extended_vtx ::\n");
  }


  int          *_out_pn_entity1     = NULL;
  PDM_g_num_t **_out_pentity1_gnum  = NULL;
  double      **_out_pentity1_coord = NULL;
  PDM_malloc(_out_pn_entity1    , part_ext->ln_part_tot, int          );
  PDM_malloc(_out_pentity1_gnum , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(_out_pentity1_coord, part_ext->ln_part_tot, double      *);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    _out_pn_entity1    [i_part] = pn_entity1[i_part]-pn_init_entity1[i_part];
    PDM_malloc(_out_pentity1_gnum [i_part],   _out_pn_entity1[i_part], PDM_g_num_t);
    PDM_malloc(_out_pentity1_coord[i_part], 3*_out_pn_entity1[i_part], double     );

    int i_write = 0;
    for (int i_entity1=pn_init_entity1[i_part];
             i_entity1<pn_entity1     [i_part]; ++i_entity1) {
      _out_pentity1_gnum [i_part][  i_write  ] = pentity1_gnum [i_part][  i_entity1  ];
      _out_pentity1_coord[i_part][3*i_write  ] = pentity1_coord[i_part][3*i_entity1  ];
      _out_pentity1_coord[i_part][3*i_write+1] = pentity1_coord[i_part][3*i_entity1+1];
      _out_pentity1_coord[i_part][3*i_write+2] = pentity1_coord[i_part][3*i_entity1+2];
      i_write++;
    }
    
    if (debug==1) {
      log_trace("\ni_part = %d\n", i_part);
      log_trace("\t _out_pn_entity1 = %d\n", _out_pn_entity1[i_part]);
      PDM_log_trace_array_long(_out_pentity1_gnum        [i_part], _out_pn_entity1[i_part], "\t out_pentity1_gnum         ::");
    }
  }

  *out_pn_entity1     = _out_pn_entity1;
  *out_pentity1_gnum  = _out_pentity1_gnum;
  *out_pentity1_coord = _out_pentity1_coord;
  // > Set ownership
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for (int i_part = 0; i_part < part_ext->n_part[i_dom]; i_part++) {
      part_ext->ownership_border_ln_to_gn [PDM_MESH_ENTITY_VTX][i_dom][i_part] = PDM_OWNERSHIP_KEEP;
      part_ext->ownership_border_vtx_coord                     [i_dom][i_part] = PDM_OWNERSHIP_KEEP;
    }
  }
}


static
void
_setup_domain_interface_in_block_frame
(
 PDM_part_extension_t *part_ext
)
{
  int debug = 0;

  part_ext->dom_itrf = NULL;
  if(part_ext->pdi == NULL) {
    return;
  }

  int n_interface = 0;
  if(part_ext->pdi != NULL) {
    n_interface = PDM_part_domain_interface_n_interface_get(part_ext->pdi);
  }

  int is_describe_vtx  = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_VTX );
  int is_describe_edge = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_EDGE);
  int is_describe_face = PDM_part_domain_interface_exist_get(part_ext->pdi, PDM_BOUND_TYPE_FACE);

  int is_describe_vtx_l  = is_describe_vtx;
  int is_describe_edge_l = is_describe_edge;
  int is_describe_face_l = is_describe_face;
  PDM_MPI_Allreduce(&is_describe_vtx_l , &is_describe_vtx , 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_edge_l, &is_describe_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
  PDM_MPI_Allreduce(&is_describe_face_l, &is_describe_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

  /*
   * Do shortcut
   */

  int          **pn_vtx         = NULL;
  int          **pn_edge        = NULL;
  int          **pn_face        = NULL;
  PDM_g_num_t ***pvtx_ln_to_gn  = NULL;
  PDM_g_num_t ***pedge_ln_to_gn = NULL;
  PDM_g_num_t ***pface_ln_to_gn = NULL;
  PDM_malloc(pn_vtx        , part_ext->n_domain, int          *);
  PDM_malloc(pn_edge       , part_ext->n_domain, int          *);
  PDM_malloc(pn_face       , part_ext->n_domain, int          *);
  PDM_malloc(pvtx_ln_to_gn , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(pedge_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(pface_ln_to_gn, part_ext->n_domain, PDM_g_num_t **);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_malloc(pn_vtx        [i_domain], part_ext->n_domain, int          );
    PDM_malloc(pn_edge       [i_domain], part_ext->n_domain, int          );
    PDM_malloc(pn_face       [i_domain], part_ext->n_domain, int          );
    PDM_malloc(pvtx_ln_to_gn [i_domain], part_ext->n_domain, PDM_g_num_t *);
    PDM_malloc(pedge_ln_to_gn[i_domain], part_ext->n_domain, PDM_g_num_t *);
    PDM_malloc(pface_ln_to_gn[i_domain], part_ext->n_domain, PDM_g_num_t *);

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      pn_vtx        [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_face       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pvtx_ln_to_gn [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      pedge_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      pface_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
    }
  }


  if(is_describe_face) {
    int **is_face_on_itrf = NULL;
    PDM_part_domain_interface_to_domain_interface(part_ext->pdi,
                                                  PDM_BOUND_TYPE_FACE,
                                                  part_ext->n_part,
                                                  pn_face,
                                                  pface_ln_to_gn,
                                                  &part_ext->dom_itrf,
                                                  &is_face_on_itrf);
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(is_face_on_itrf[i_part]);
    }
    PDM_free(is_face_on_itrf);
  }

  if(is_describe_edge) {
    int **is_edge_on_itrf = NULL;
    PDM_part_domain_interface_to_domain_interface(part_ext->pdi,
                                                  PDM_BOUND_TYPE_EDGE,
                                                  part_ext->n_part,
                                                  pn_edge,
                                                  pedge_ln_to_gn,
                                                  &part_ext->dom_itrf,
                                                  &is_edge_on_itrf);
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(is_edge_on_itrf[i_part]);
    }
    PDM_free(is_edge_on_itrf);
  }

  if(is_describe_vtx) {
    int **is_vtx_on_itrf = NULL;
    PDM_part_domain_interface_to_domain_interface(part_ext->pdi,
                                                  PDM_BOUND_TYPE_VTX,
                                                  part_ext->n_part,
                                                  pn_vtx,
                                                  pvtx_ln_to_gn,
                                                  &part_ext->dom_itrf,
                                                  &is_vtx_on_itrf);
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(is_vtx_on_itrf[i_part]);
    }
    PDM_free(is_vtx_on_itrf);
  }

  /*
   * Free all
   */
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(pn_vtx        [i_domain]);
    PDM_free(pn_edge       [i_domain]);
    PDM_free(pn_face       [i_domain]);
    PDM_free(pvtx_ln_to_gn [i_domain]);
    PDM_free(pedge_ln_to_gn[i_domain]);
    PDM_free(pface_ln_to_gn[i_domain]);
  }
  PDM_free(pn_vtx        );
  PDM_free(pn_edge       );
  PDM_free(pn_face       );
  PDM_free(pvtx_ln_to_gn );
  PDM_free(pedge_ln_to_gn);
  PDM_free(pface_ln_to_gn);

  /*
   * At this stage we have for all interface the information for all entities
   *   During the part_extension problem we need this information to merge / amend all extend entity
   *   We make a distributed vision that more convenient to do all this kind of operation
   */
  part_ext->n_interface = n_interface;

  // for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
  //   part_ext->ptb_itrf[i] = NULL; // (PDM_part_to_block_t **) PDM_malloc( n_interface * sizeof(PDM_part_to_block_t **));
  //   part_ext->opp_gnum[i] = NULL; // (PDM_g_num_t         **) PDM_malloc( n_interface * sizeof(PDM_g_num_t         **));
  //   part_ext->opp_sens[i] = NULL; // (PDM_g_num_t         **) PDM_malloc( n_interface * sizeof(PDM_g_num_t         **));

  //   // for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
  //   //   part_ext->ptb_itrf[i][i_itrf] = NULL;
  //   //   part_ext->opp_gnum[i][i_itrf] = NULL;
  //   // }
  // }

  // PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
  //                                     PDM_BOUND_TYPE_VTX,
  //                                     part_ext->shift_by_domain_vtx,
  //                                     &part_ext->ptb_itrf[PDM_BOUND_TYPE_VTX],
  //                                     &part_ext->opp_gnum[PDM_BOUND_TYPE_VTX],
  //                                     &part_ext->opp_sens[PDM_BOUND_TYPE_VTX]);


  // if(is_describe_edge) {
  //   PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
  //                                       PDM_BOUND_TYPE_EDGE,
  //                                       part_ext->shift_by_domain_edge,
  //                                      &part_ext->ptb_itrf[PDM_BOUND_TYPE_EDGE],
  //                                      &part_ext->opp_gnum[PDM_BOUND_TYPE_EDGE],
  //                                      &part_ext->opp_sens[PDM_BOUND_TYPE_EDGE]);
  // }


  // if(is_describe_face) {
  //   PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
  //                                       PDM_BOUND_TYPE_FACE,
  //                                       part_ext->shift_by_domain_face,
  //                                      &part_ext->ptb_itrf[PDM_BOUND_TYPE_FACE],
  //                                      &part_ext->opp_gnum[PDM_BOUND_TYPE_FACE],
  //                                      &part_ext->opp_sens[PDM_BOUND_TYPE_FACE]);
  // }


  assert(part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] == 0   );
  assert(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] == NULL);
  assert(part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] == NULL);
  assert(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] == NULL);
  assert(part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX] == NULL);

  if(part_ext->dom_itrf != NULL) {
    PDM_domain_interface_make_flat_view2(part_ext->dom_itrf,
                                         PDM_BOUND_TYPE_VTX,
                                         part_ext->shift_by_domain_vtx,
                                         &part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX],
                                         &part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX],
                                         &part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX],
                                         &part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX],
                                         &part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX]);
    if(is_describe_edge) {
      PDM_domain_interface_make_flat_view2(part_ext->dom_itrf,
                                           PDM_BOUND_TYPE_EDGE,
                                           part_ext->shift_by_domain_edge,
                                           &part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE],
                                           &part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE]);
    }
    if(is_describe_face) {
      PDM_domain_interface_make_flat_view2(part_ext->dom_itrf,
                                           PDM_BOUND_TYPE_FACE,
                                           part_ext->shift_by_domain_face,
                                           &part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_FACE],
                                           &part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_FACE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_FACE],
                                           &part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_FACE]);
    }
    // PDM_free(part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX]);
  } else {
    part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = 0;
    part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX] = NULL;

    part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE] = 0;
    part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE] = NULL;

    part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_FACE] = 0;
    part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_FACE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_FACE] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_FACE] = NULL;
  }


  if(debug == 1) {
    int n_data_vtx = 0;
    for(int i = 0; i < part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX]; ++i) {
      n_data_vtx += part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX][i];
    }
    log_trace("\n");
    log_trace("Vertices\n");
    PDM_log_trace_array_long(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX], "dvtx_itrf_blk_gnum            ::");
    PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX], "dvtx_itrf_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX], 2 * n_data_vtx                                  , "dvtx_itrf_gnum_and_itrf_data  ::");
    PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX], n_data_vtx                                      , "dvtx_itrf_gnum_and_itrf_sens  ::");
    
    if(is_describe_edge) {
      int n_data_edge = 0;
      for(int i = 0; i < part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_EDGE]; ++i) {
        n_data_edge += part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE][i];
      }
      log_trace("\n");
      log_trace("Edges\n");
      PDM_log_trace_array_long(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_EDGE], "dedge_itrf_blk_gnum            ::");
      PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_EDGE], "dedge_itrf_gnum_and_itrf_strid ::");
      PDM_log_trace_array_long(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE], 2 * n_data_edge                                  , "dedge_itrf_gnum_and_itrf_data  ::");
      PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE], n_data_edge                                      , "dedge_itrf_gnum_and_itrf_sens  ::");
    }

    if(is_describe_face) {
      int n_data_face = 0;
      for(int i = 0; i < part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_FACE]; ++i) {
        n_data_face += part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE][i];
      }
      log_trace("\n");
      log_trace("Faces\n");
      PDM_log_trace_array_long(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_FACE], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_FACE], "dface_itrf_blk_gnum            ::");
      PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_FACE], "dface_itrf_gnum_and_itrf_strid ::");
      PDM_log_trace_array_long(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_FACE], 2 * n_data_face                                  , "dface_itrf_gnum_and_itrf_data  ::");
      PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_FACE], n_data_face                                      , "dface_itrf_gnum_and_itrf_sens  ::");
    }
  }


  // PDM_free(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX]);
  // PDM_free(part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX]);
  // PDM_free(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX]);
  // PDM_free(part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX]);
}

static
void
_build_bound_graph
(
 PDM_part_extension_t   *part_ext
)
{
  /*
   * Dans tous les cas on cherche a obtenir le graphe entre le propagateur et l'entité principale :
   *    - PDM_EXTEND_FROM_VTX
   *    - PDM_EXTEND_FROM_EDGE
   *    - PDM_EXTEND_FROM_FACE
   */
  int debug = 0;

  int          **pn_entity1            = NULL;
  PDM_g_num_t ***pentity1_ln_to_gn     = NULL;
  int         ***pentity1_bnd_proc_idx = NULL;
  int         ***pentity1_bnd_part_idx = NULL;
  int         ***pentity1_bnd          = NULL;
  PDM_malloc(pn_entity1           , part_ext->n_domain, int          *);
  PDM_malloc(pentity1_ln_to_gn    , part_ext->n_domain, PDM_g_num_t **);
  PDM_malloc(pentity1_bnd_proc_idx, part_ext->n_domain, int         **);
  PDM_malloc(pentity1_bnd_part_idx, part_ext->n_domain, int         **);
  PDM_malloc(pentity1_bnd         , part_ext->n_domain, int         **);
  PDM_bound_type_t bound_type = PDM_BOUND_TYPE_MAX;
  if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
    bound_type = PDM_BOUND_TYPE_VTX;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_malloc(pn_entity1           [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pentity1_ln_to_gn    [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pentity1_bnd_proc_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd_part_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd         [i_domain], part_ext->n_domain, int         *);
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1           [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        pentity1_ln_to_gn    [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        pentity1_bnd_proc_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx;
        pentity1_bnd_part_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx;
        pentity1_bnd         [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_part_bound;
      }
    }
  } else if(part_ext->extend_type == PDM_EXTEND_FROM_EDGE) {
    bound_type = PDM_BOUND_TYPE_EDGE;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_malloc(pn_entity1           [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pentity1_ln_to_gn    [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pentity1_bnd_proc_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd_part_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd         [i_domain], part_ext->n_domain, int         *);
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1           [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
        pentity1_ln_to_gn    [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
        pentity1_bnd_proc_idx[i_domain][i_part] = NULL;
        pentity1_bnd_part_idx[i_domain][i_part] = NULL;
        pentity1_bnd         [i_domain][i_part] = NULL;
      }
    }
  } else if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
    bound_type = PDM_BOUND_TYPE_FACE;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_malloc(pn_entity1           [i_domain], part_ext->n_domain, int          );
      PDM_malloc(pentity1_ln_to_gn    [i_domain], part_ext->n_domain, PDM_g_num_t *);
      PDM_malloc(pentity1_bnd_proc_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd_part_idx[i_domain], part_ext->n_domain, int         *);
      PDM_malloc(pentity1_bnd         [i_domain], part_ext->n_domain, int         *);
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1           [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        pentity1_ln_to_gn    [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
        pentity1_bnd_proc_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_part_bound_proc_idx;
        pentity1_bnd_part_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_part_bound_part_idx;
        pentity1_bnd         [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_part_bound;
      }
    }
  }

  /*
   * 1st step :
   *   - Connectivity between partition have two kind :
   *      + by partitionning interface (same domain)
   *      + by domaine interface (other domain or periodic or multidomain)
   *   - This step give rebuild a connectivity graphe with both contribution
   *      + The first is deduce by global numbering
   *      + The second is deduce by all domain_interface give by the user
   */
  int **pentity_bound_to_pentity_bound_idx       = NULL;
  int **pentity_bound_to_pentity_bound_triplet   = NULL;
  int **pentity_bound_to_pentity_bound_interface = NULL;
  PDM_part_extension_build_entity1_graph(part_ext->pdi,
                                         bound_type,
                                         part_ext->n_domain,
                                         part_ext->n_part,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         pentity1_bnd_proc_idx,
                                         pentity1_bnd_part_idx,
                                         pentity1_bnd,
                                         NULL,
                                         part_ext->user_defined_bnd_graph,
                                         &pentity_bound_to_pentity_bound_idx,
                                         &pentity_bound_to_pentity_bound_triplet,
                                         &pentity_bound_to_pentity_bound_interface,
                                         part_ext->comm);

  if(debug == 1) {
    int l_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        PDM_log_trace_array_int(pentity_bound_to_pentity_bound_idx      [l_part], pn_entity1[i_domain][i_part], "pentity_bound_to_pentity_bound_idx ::");
        PDM_log_trace_array_int(pentity_bound_to_pentity_bound_triplet  [l_part], pentity_bound_to_pentity_bound_idx[l_part][pn_entity1[i_domain][i_part]], "pentity_bound_to_pentity_bound_triplet ::");
        PDM_log_trace_array_int(pentity_bound_to_pentity_bound_interface[l_part], pentity_bound_to_pentity_bound_idx[l_part][pn_entity1[i_domain][i_part]]/3, "pentity_bound_to_pentity_bound_interface ::");
        l_part++;
      }
    }
  }

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(pn_entity1           [i_domain]);
    PDM_free(pentity1_ln_to_gn    [i_domain]);
    PDM_free(pentity1_bnd_proc_idx[i_domain]);
    PDM_free(pentity1_bnd_part_idx[i_domain]);
    PDM_free(pentity1_bnd         [i_domain]);
  }
  PDM_free(pn_entity1);
  PDM_free(pentity1_ln_to_gn);
  PDM_free(pentity1_bnd_proc_idx);
  PDM_free(pentity1_bnd_part_idx);
  PDM_free(pentity1_bnd);

  part_ext->pinit_entity_bound_to_pentity_bound_idx       = pentity_bound_to_pentity_bound_idx;
  part_ext->pinit_entity_bound_to_pentity_bound_triplet   = pentity_bound_to_pentity_bound_triplet;
  part_ext->pinit_entity_bound_to_pentity_bound_interface = pentity_bound_to_pentity_bound_interface;

}


static
void
_build_rotation_matrix
(
  double  *rotation_direction,
  double   rotation_angle,
  int      sgn_itrf,
  double **rotation_matrix
)
{
  double angle = -sgn_itrf*rotation_angle;

  if (PDM_ABS(rotation_direction[0])>1e-15 && 
      PDM_ABS(rotation_direction[1])<1e-15 &&
      PDM_ABS(rotation_direction[2])<1e-15) {
    rotation_matrix[0][0] = 1.;
    rotation_matrix[0][1] = 0.;
    rotation_matrix[0][2] = 0.;

    rotation_matrix[1][0] = 0.;
    rotation_matrix[1][1] = cos(angle);
    rotation_matrix[1][2] =-sin(angle);

    rotation_matrix[2][0] = 0.;
    rotation_matrix[2][1] = sin(angle);
    rotation_matrix[2][2] = cos(angle);
  }
  else if (PDM_ABS(rotation_direction[0])<1e-15 && 
           PDM_ABS(rotation_direction[1])>1e-15 &&
           PDM_ABS(rotation_direction[2])<1e-15) {
    rotation_matrix[0][0] = cos(angle);
    rotation_matrix[0][1] = 0.;
    rotation_matrix[0][2] = sin(angle);

    rotation_matrix[1][0] = 0.;
    rotation_matrix[1][1] = 1.;
    rotation_matrix[1][2] = 0.;

    rotation_matrix[2][0] =-sin(angle);
    rotation_matrix[2][1] = 0.;
    rotation_matrix[2][2] = cos(angle);
  }
  else if (PDM_ABS(rotation_direction[0])<1e-15 && 
           PDM_ABS(rotation_direction[1])<1e-15 &&
           PDM_ABS(rotation_direction[2])>1e-15) {
    rotation_matrix[0][0] = cos(angle);
    rotation_matrix[0][1] =-sin(angle);
    rotation_matrix[0][2] = 0.;

    rotation_matrix[1][0] = sin(angle);
    rotation_matrix[1][1] = cos(angle);
    rotation_matrix[1][2] = 0.;

    rotation_matrix[2][0] = 0.;
    rotation_matrix[2][1] = 0.;
    rotation_matrix[2][2] = 1.;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Don't know how to build rotation matrix around axis (%.3f,%.3f,%.3f)\n",
                                  rotation_direction[0],
                                  rotation_direction[1],
                                  rotation_direction[2]);
  } 
}


static
void
_build_homogeneous_matrix
(
  double **rotation_matrix,
  double  *rotation_center,
  double **homogeneous_matrix
)
{
  double apply_center_matrix[3] = {0, 0., 0.};

  for (int i=0; i<3; ++i) {
    apply_center_matrix[i] = rotation_matrix[0][i]*rotation_center[0]
                           + rotation_matrix[1][i]*rotation_center[1]
                           + rotation_matrix[2][i]*rotation_center[2];

    homogeneous_matrix[i][0] = rotation_matrix[i][0];
    homogeneous_matrix[i][1] = rotation_matrix[i][1];
    homogeneous_matrix[i][2] = rotation_matrix[i][2];
    homogeneous_matrix[i][3] = rotation_center[i] - apply_center_matrix[i];
  }

  homogeneous_matrix[3][0] = 0.;
  homogeneous_matrix[3][1] = 0.;
  homogeneous_matrix[3][2] = 0.;
  homogeneous_matrix[3][3] = 1.;
}


static
void
_apply_homogeneous_matrix
(
  double **homogeneous_matrix,
  double  *vector
)
{

  double _vector[4] = {0., 0., 0., 0.};
  double _result[4] = {0., 0., 0., 0.};

  _vector[0] = vector[0];
  _vector[1] = vector[1];
  _vector[2] = vector[2];
  _vector[3] = 1.;

  for (int i=0; i<4; ++i) {
    _result[i] = homogeneous_matrix[0][i]*_vector[0]
               + homogeneous_matrix[1][i]*_vector[1]
               + homogeneous_matrix[2][i]*_vector[2]
               + homogeneous_matrix[3][i]*_vector[3];
  }

  for (int i=0; i<3; ++i) {
    vector[i] = _result[i];
  }
}


static
void
_exchange_coord_and_apply_transform
(
PDM_part_extension_t  *part_ext,
int                   *pn_vtx_extended,
PDM_g_num_t          **pvtx_extended_ln_to_gn,
int                   *pn_vtx,
double               **pvtx_coords,
int                  **pvtx_extended_to_pvtx_idx,
int                  **pvtx_extended_to_pvtx_triplet,
int                  **pvtx_extended_to_pvtx_interface,
double              ***pvtx_extended_coords_out
)
{
  int debug = 0;

  // > Unique graph
  int **_pvtx_extended_to_pvtx_idx   = NULL;
  int **_pvtx_extended_to_pvtx_trplt = NULL;
  int **_pvtx_extended_to_pvtx_itrf  = NULL;
  PDM_malloc(_pvtx_extended_to_pvtx_idx  , part_ext->ln_part_tot, int *);
  PDM_malloc(_pvtx_extended_to_pvtx_trplt, part_ext->ln_part_tot, int *);
  PDM_malloc(_pvtx_extended_to_pvtx_itrf , part_ext->ln_part_tot, int *);
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_malloc(_pvtx_extended_to_pvtx_idx  [i_part],   pn_vtx_extended[i_part]+1, int);
    PDM_malloc(_pvtx_extended_to_pvtx_trplt[i_part], 3*pn_vtx_extended[i_part]  , int);
    PDM_malloc(_pvtx_extended_to_pvtx_itrf [i_part],   pn_vtx_extended[i_part]  , int);
    _pvtx_extended_to_pvtx_idx  [i_part][0] = 0;
    for (int i_vtx=0; i_vtx<pn_vtx_extended[i_part]; ++i_vtx) {
      int i_read = pvtx_extended_to_pvtx_idx[i_part][i_vtx]/3;
      _pvtx_extended_to_pvtx_trplt[i_part][3*i_vtx  ] =  pvtx_extended_to_pvtx_triplet  [i_part][3*i_read  ];
      _pvtx_extended_to_pvtx_trplt[i_part][3*i_vtx+1] =  pvtx_extended_to_pvtx_triplet  [i_part][3*i_read+1];
      _pvtx_extended_to_pvtx_trplt[i_part][3*i_vtx+2] =  pvtx_extended_to_pvtx_triplet  [i_part][3*i_read+2];
      _pvtx_extended_to_pvtx_itrf [i_part][  i_vtx  ] =  pvtx_extended_to_pvtx_interface[i_part][  i_read  ];
      _pvtx_extended_to_pvtx_idx  [i_part][  i_vtx+1] = _pvtx_extended_to_pvtx_idx      [i_part][  i_vtx   ]+3;
    }
  }

  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_extended_ln_to_gn,
                                                                          (const int          *) pn_vtx_extended,
                                                                          part_ext->ln_part_tot,
                                                                          (const int          *) pn_vtx,
                                                                          part_ext->ln_part_tot,
                                                                          (const int         **) _pvtx_extended_to_pvtx_idx,
                                                                          (const int         **) NULL,
                                                                          (const int         **) _pvtx_extended_to_pvtx_trplt,
                                                                          part_ext->comm);
  /*
   *
   */
  int exch_request = -1;
  double      **pextract_vtx_coords           = NULL;
  PDM_part_to_part_reverse_iexch(ptp_vtx,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 3 * sizeof(double),
                                 NULL,
                (const void **)  pvtx_coords,
                                 NULL,
                    (void ***)   &pextract_vtx_coords,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp_vtx, exch_request);

  PDM_part_to_part_free(ptp_vtx);

  /*
   * Apply transformation if any
   */
  int n_interface = 0;
  if(part_ext->pdi != NULL) {
    n_interface = PDM_part_domain_interface_n_interface_get(part_ext->pdi);
  }
  double  **translation_vector = NULL;
  double ***rotation_matrix    = NULL;
  double ***homogeneous_matrix = NULL;
  double  **rotation_direction = NULL;
  double  **rotation_center    = NULL;
  double   *rotation_angle     = NULL;
  PDM_malloc(translation_vector, n_interface, double * );
  PDM_malloc(rotation_matrix   , n_interface, double **);
  PDM_malloc(homogeneous_matrix, n_interface, double **);
  PDM_malloc(rotation_direction, n_interface, double * );
  PDM_malloc(rotation_center   , n_interface, double * );
  PDM_malloc(rotation_angle    , n_interface, double   );
  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    translation_vector[i_interf] = NULL;
    PDM_part_domain_interface_translation_get(part_ext->pdi, i_interf, &translation_vector[i_interf]);

    rotation_matrix   [i_interf] = NULL;
    homogeneous_matrix[i_interf] = NULL;
    PDM_part_domain_interface_rotation_get(part_ext->pdi,
                                           i_interf,
                                           &rotation_direction[i_interf],
                                           &rotation_center   [i_interf],
                                           &rotation_angle    [i_interf]);
    if(rotation_center[i_interf] != NULL) {
      
      // > Allocate tmp and transfo matrix
      PDM_malloc(rotation_matrix   [i_interf], 3, double *);
      PDM_malloc(homogeneous_matrix[i_interf], 4, double *);
      for(int k = 0; k < 3; ++k) {
        PDM_malloc(rotation_matrix[i_interf][k], 3, double);
      }
      for(int k = 0; k < 4; ++k) {
        PDM_malloc(homogeneous_matrix[i_interf][k], 4, double);
      }
    }
  }

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    for(int i_vtx = 0; i_vtx < pn_vtx_extended[i_part]; ++i_vtx) {
      int i_interface   = PDM_ABS (_pvtx_extended_to_pvtx_itrf[i_part][i_vtx]);
      int sgn_interface = PDM_SIGN(_pvtx_extended_to_pvtx_itrf[i_part][i_vtx]);
      if(i_interface != 0 && translation_vector[PDM_ABS(i_interface)-1] != NULL) {
        for(int k = 0; k < 3; ++k) {
          pextract_vtx_coords[i_part][3*i_vtx+k] += sgn_interface * translation_vector[PDM_ABS(i_interface)-1][k];
        }
      }
      if(i_interface != 0 && rotation_direction[PDM_ABS(i_interface)-1] != NULL) {

        // > Build matrix
        _build_rotation_matrix(rotation_direction[PDM_ABS(i_interface-1)],
                               rotation_angle    [PDM_ABS(i_interface-1)],
                               sgn_interface,
                               rotation_matrix   [PDM_ABS(i_interface-1)]);

        _build_homogeneous_matrix(rotation_matrix   [PDM_ABS(i_interface-1)],
                                  rotation_center   [PDM_ABS(i_interface-1)],
                                  homogeneous_matrix[PDM_ABS(i_interface-1)]);

        // > Apply matrix
        _apply_homogeneous_matrix(homogeneous_matrix[PDM_ABS(i_interface-1)],
                                 &pextract_vtx_coords[i_part][3*i_vtx]);
      }
    }

    PDM_free(_pvtx_extended_to_pvtx_idx  [i_part]);
    PDM_free(_pvtx_extended_to_pvtx_trplt[i_part]);
    PDM_free(_pvtx_extended_to_pvtx_itrf [i_part]);
  }
  PDM_free(_pvtx_extended_to_pvtx_idx);
  PDM_free(_pvtx_extended_to_pvtx_trplt);
  PDM_free(_pvtx_extended_to_pvtx_itrf);


  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    if(translation_vector[i_interf] != NULL) {
      PDM_free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        PDM_free(rotation_matrix[i_interf][k]);
      }
      for(int k = 0; k < 4; ++k) {
        PDM_free(homogeneous_matrix[i_interf][k]);
      }
      PDM_free(rotation_matrix[i_interf]);
      PDM_free(homogeneous_matrix[i_interf]);
    }
  }
  PDM_free(translation_vector);
  PDM_free(homogeneous_matrix);
  PDM_free(rotation_matrix);
  PDM_free(rotation_direction);
  PDM_free(rotation_center);
  PDM_free(rotation_angle);

  /*
   * Petit vtk en légende
   */
  if(debug == 1) {
    int i_rank;
    PDM_MPI_Comm_rank(part_ext->comm, &i_rank);

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      char filename[999];
      sprintf(filename, "extended_vtx_coords_%i_%i.vtk", i_part, i_rank);

      PDM_vtk_write_point_cloud(filename,
                                pn_vtx_extended       [i_part],
                                pextract_vtx_coords   [i_part],
                                pvtx_extended_ln_to_gn[i_part],
                                NULL);
    }
  }

  *pvtx_extended_coords_out = pextract_vtx_coords;
}

static
void
_part_extension_3d
(
 PDM_part_extension_t *part_ext
)
{

  int visu  = 0;
  int debug = 0;

  /*
   * 2 possibilities :
   *   - With face_vtx
   *   - With face_edge + edge_vtx
   */

  /* Size */
  int *pn_vtx       = NULL;
  int *pn_edge      = NULL;
  int *pn_face      = NULL;
  int *pn_cell      = NULL;
  int *pn_init_vtx  = NULL;
  int *pn_init_edge = NULL;
  int *pn_init_face = NULL;
  int *pn_init_cell = NULL;
  PDM_malloc(pn_vtx      , part_ext->ln_part_tot, int );
  PDM_malloc(pn_edge     , part_ext->ln_part_tot, int );
  PDM_malloc(pn_face     , part_ext->ln_part_tot, int );
  PDM_malloc(pn_cell     , part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_vtx , part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_edge, part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_face, part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_cell, part_ext->ln_part_tot, int );

  /* ln_to_gn */
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_g_num_t **pcell_ln_to_gn = NULL;
  PDM_malloc(pvtx_ln_to_gn , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pedge_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pface_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pcell_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);

  /* Connectivity */
  int **pcell_face_idx = NULL;
  int **pcell_face     = NULL;
  int **pcell_vtx_idx  = NULL;
  int **pcell_vtx      = NULL;
  int **pface_edge_idx = NULL;
  int **pface_edge     = NULL;
  int **pface_vtx_idx  = NULL;
  int **pface_vtx      = NULL;
  int **pedge_vtx_idx  = NULL;
  int **pedge_vtx      = NULL;
  PDM_malloc(pcell_face_idx, part_ext->ln_part_tot, int *);
  PDM_malloc(pcell_face    , part_ext->ln_part_tot, int *);
  PDM_malloc(pcell_vtx_idx , part_ext->ln_part_tot, int *);
  PDM_malloc(pcell_vtx     , part_ext->ln_part_tot, int *);
  PDM_malloc(pface_edge_idx, part_ext->ln_part_tot, int *);
  PDM_malloc(pface_edge    , part_ext->ln_part_tot, int *);
  PDM_malloc(pface_vtx_idx , part_ext->ln_part_tot, int *);
  PDM_malloc(pface_vtx     , part_ext->ln_part_tot, int *);
  PDM_malloc(pedge_vtx_idx , part_ext->ln_part_tot, int *);
  PDM_malloc(pedge_vtx     , part_ext->ln_part_tot, int *);

  /* Groups */
  int          *pn_face_group    = NULL;
  int         **pface_group_tag  = NULL;
  PDM_g_num_t **pface_group_gnum = NULL;
  PDM_malloc(pn_face_group   , part_ext->ln_part_tot, int          );
  PDM_malloc(pface_group_tag , part_ext->ln_part_tot, int         *);
  PDM_malloc(pface_group_gnum, part_ext->ln_part_tot, PDM_g_num_t *);

  /* Coordinates */
  double **pvtx_coords = NULL;
  PDM_malloc(pvtx_coords, part_ext->ln_part_tot, double *);

  int **pcell_alrdy_sent = NULL;
  PDM_malloc(pcell_alrdy_sent, part_ext->ln_part_tot, int *);

  int lpart = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      pn_vtx          [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_init_vtx     [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      if(part_ext->have_edge == 1) {
        pn_edge       [lpart] = part_ext->parts[i_domain][i_part].n_edge;
        pn_init_edge  [lpart] = part_ext->parts[i_domain][i_part].n_edge;
      }
      pn_face         [lpart] = part_ext->parts[i_domain][i_part].n_face;
      pn_init_face    [lpart] = part_ext->parts[i_domain][i_part].n_face;
      pn_cell         [lpart] = part_ext->parts[i_domain][i_part].n_cell;
      pn_init_cell    [lpart] = part_ext->parts[i_domain][i_part].n_cell;


      /* Copy to realloc after all step */
      PDM_malloc(pvtx_ln_to_gn   [lpart], pn_vtx [lpart], PDM_g_num_t);
      if(part_ext->have_edge == 1) {
        PDM_malloc(pedge_ln_to_gn[lpart], pn_edge[lpart], PDM_g_num_t);
      }
      PDM_malloc(pface_ln_to_gn  [lpart], pn_face[lpart], PDM_g_num_t);
      PDM_malloc(pcell_ln_to_gn  [lpart], pn_cell[lpart], PDM_g_num_t);
      PDM_malloc(pcell_alrdy_sent[lpart], pn_face[lpart], int        );


      for(int i_cell = 0; i_cell < pn_cell[lpart]; ++i_cell) {
        pcell_ln_to_gn  [lpart][i_cell] = part_ext->parts[i_domain][i_part].cell_ln_to_gn[i_cell];
        pcell_alrdy_sent[lpart][i_cell] = 0;
      }
      for(int i_face = 0; i_face < pn_face[lpart]; ++i_face) {
        pface_ln_to_gn[lpart][i_face] = part_ext->parts[i_domain][i_part].face_ln_to_gn[i_face];
      }
      if(part_ext->have_edge == 1) {
        for(int i_edge = 0; i_edge < pn_edge[lpart]; ++i_edge) {
          pedge_ln_to_gn[lpart][i_edge] = part_ext->parts[i_domain][i_part].edge_ln_to_gn[i_edge];
        }
      }
      for(int i_vtx = 0; i_vtx < pn_vtx[lpart]; ++i_vtx) {
        pvtx_ln_to_gn[lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn[i_vtx];
      }

      /* cell_face connectivity */
      PDM_malloc(pcell_face_idx[lpart],                                                 pn_cell[lpart]+1, int);
      PDM_malloc(pcell_face    [lpart], part_ext->parts[i_domain][i_part].cell_face_idx[pn_cell[lpart]] , int);

      for(int i_cell = 0; i_cell < pn_cell[lpart]+1; ++i_cell) {
        pcell_face_idx [lpart][i_cell] = part_ext->parts[i_domain][i_part].cell_face_idx[i_cell];
      }

      for(int idx = 0; idx < pcell_face_idx[i_part][pn_cell[lpart]]; ++idx) {
        pcell_face [lpart][idx] = part_ext->parts[i_domain][i_part].cell_face[idx];
      }

      /* face_vtx connectivity */
      if(part_ext->have_edge == 0) {
        PDM_malloc(pface_vtx_idx[lpart],                                                pn_face[lpart]+1, int);
        PDM_malloc(pface_vtx    [lpart], part_ext->parts[i_domain][i_part].face_vtx_idx[pn_face[lpart]] , int);
        for(int i_face = 0; i_face < pn_face[lpart]+1; ++i_face) {
          pface_vtx_idx [lpart][i_face] = part_ext->parts[i_domain][i_part].face_vtx_idx[i_face];
        }

        for(int idx = 0; idx < pface_vtx_idx[i_part][pn_face[lpart]]; ++idx) {
          pface_vtx [lpart][idx] = part_ext->parts[i_domain][i_part].face_vtx[idx];
        }
      }

      if(part_ext->have_edge == 1) {
        /* face_edge connectivity */
        PDM_malloc(pface_edge_idx[lpart],                                                 pn_face[lpart]+1, int);
        PDM_malloc(pface_edge    [lpart], part_ext->parts[i_domain][i_part].face_edge_idx[pn_face[lpart]] , int);
        for(int i_face = 0; i_face < pn_face[lpart]+1; ++i_face) {
          pface_edge_idx [lpart][i_face] = part_ext->parts[i_domain][i_part].face_edge_idx[i_face];
        }

        for(int idx = 0; idx < pface_edge_idx[i_part][pn_face[lpart]]; ++idx) {
          pface_edge [lpart][idx] = part_ext->parts[i_domain][i_part].face_edge[idx];
        }
        
        /* edge_vtx connectivity */
        PDM_malloc(pedge_vtx_idx[lpart],     pn_edge[lpart]+1, int);
        PDM_malloc(pedge_vtx    [lpart], 2 * pn_edge[lpart]  , int);
        for(int i_edge = 0; i_edge < pn_edge[lpart]+1; ++i_edge) {
          pedge_vtx_idx [lpart][i_edge] = 2 * i_edge;
        }

        for(int idx = 0; idx < pedge_vtx_idx[i_part][pn_edge[lpart]]; ++idx) {
          pedge_vtx [lpart][idx] = part_ext->parts[i_domain][i_part].edge_vtx[idx];
        }
      }

      PDM_calloc(pface_group_tag [lpart], pn_face[lpart], int        );
      PDM_calloc(pface_group_gnum[lpart], pn_face[lpart], PDM_g_num_t);

      /* Coordinnates */
      PDM_malloc(pvtx_coords[lpart], 3 * pn_vtx [lpart], double);
      for(int i_vtx = 0; i_vtx < 3 * pn_vtx[lpart]; ++i_vtx) {
        pvtx_coords   [lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx[i_vtx];
      }


      // > Transform group into tag

      int          _n_group           = part_ext->parts[i_domain][i_part].n_face_group;
      int         *_group_entity_idx  = part_ext->parts[i_domain][i_part].face_bound_idx;
      int         *_group_entity      = part_ext->parts[i_domain][i_part].face_bound;
      PDM_g_num_t *_group_entity_gnum = part_ext->parts[i_domain][i_part].face_bound_ln_to_gn;
      pn_face_group[lpart] = _n_group;
      for (int i_group=0; i_group<_n_group; ++i_group) {
        for (int i_entity=_group_entity_idx[i_group  ];
                 i_entity<_group_entity_idx[i_group+1];
                 ++i_entity) {
          int lnum = _group_entity     [i_entity];
          int gnum = _group_entity_gnum[i_entity];
          pface_group_tag [lpart][lnum-1] = i_group+1;
          pface_group_gnum[lpart][lnum-1] = gnum;
        }
      }

      lpart++;
    }
  }


  /*
   * On va etendre la partition avec le graphe de base tant que de nouveau elements apparaissent
   *   -> Il faut ajuster la taille du graphe en fonction des nouvelles entités (juste une rallonge)
   *   -> On doit également alimenter un tableau pour le lien avec les entités de la recursion d'après
   */
  int           *pn_vtx_extended_old                  = NULL;
  int           *pfull_n_vtx_extended                 = NULL;
  PDM_g_num_t  **pfull_vtx_extended_ln_to_gn          = NULL;
  int          **pfull_vtx_extended_to_pvtx_idx       = NULL;
  int          **pfull_vtx_extended_to_pvtx_triplet   = NULL;
  int          **pfull_vtx_extended_to_pvtx_interface = NULL;
  PDM_malloc(pn_vtx_extended_old                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_vtx_extended                , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_vtx_extended_ln_to_gn         , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_interface, part_ext->ln_part_tot, int         *);

  int           *pn_edge_extended_old                   = NULL;
  int           *pfull_n_edge_extended                  = NULL;
  PDM_g_num_t  **pfull_edge_extended_ln_to_gn           = NULL;
  int          **pfull_edge_extended_to_pedge_idx       = NULL;
  int          **pfull_edge_extended_to_pedge_triplet   = NULL;
  int          **pfull_edge_extended_to_pedge_interface = NULL;
  PDM_malloc(pn_edge_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_edge_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_edge_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_edge_extended_to_pedge_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_interface, part_ext->ln_part_tot, int         *);

  int           *pn_face_extended_old                   = NULL;
  int           *pfull_n_face_extended                  = NULL;
  PDM_g_num_t  **pfull_face_extended_ln_to_gn           = NULL;
  int          **pfull_face_extended_to_pface_idx       = NULL;
  int          **pfull_face_extended_to_pface_triplet   = NULL;
  int          **pfull_face_extended_to_pface_interface = NULL;
  PDM_malloc(pn_face_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_face_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_face_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_face_extended_to_pface_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_face_extended_to_pface_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_face_extended_to_pface_interface, part_ext->ln_part_tot, int         *);

  int           *pn_cell_extended_old                   = NULL;
  int           *pfull_n_cell_extended                  = NULL;
  PDM_g_num_t  **pfull_cell_extended_ln_to_gn           = NULL;
  int          **pfull_cell_extended_to_pcell_idx       = NULL;
  int          **pfull_cell_extended_to_pcell_triplet   = NULL;
  int          **pfull_cell_extended_to_pcell_interface = NULL;
  PDM_malloc(pn_cell_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_cell_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_cell_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_cell_extended_to_pcell_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_cell_extended_to_pcell_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_cell_extended_to_pcell_interface, part_ext->ln_part_tot, int         *);


  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

    pfull_n_vtx_extended[i_part] = 0;
    PDM_malloc(pfull_vtx_extended_to_pvtx_idx      [i_part], pfull_n_vtx_extended[i_part]+1, int);
    PDM_malloc(pfull_vtx_extended_ln_to_gn         [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_vtx_extended_to_pvtx_interface[i_part], 0, int        );
    pfull_vtx_extended_to_pvtx_idx[i_part][0] = 0;

    if (part_ext->have_edge==1) {
      pfull_n_edge_extended[i_part] = 0;
      PDM_malloc(pfull_edge_extended_to_pedge_idx      [i_part], pfull_n_edge_extended[i_part]+1, int);
      PDM_malloc(pfull_edge_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
      PDM_malloc(pfull_edge_extended_to_pedge_triplet  [i_part], 0, int        );
      PDM_malloc(pfull_edge_extended_to_pedge_interface[i_part], 0, int        );
      pfull_edge_extended_to_pedge_idx[i_part][0] = 0;
    }

    pfull_n_face_extended[i_part] = 0;
    PDM_malloc(pfull_face_extended_to_pface_idx      [i_part], pfull_n_face_extended[i_part]+1, int);
    PDM_malloc(pfull_face_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_face_extended_to_pface_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_face_extended_to_pface_interface[i_part], 0, int        );
    pfull_face_extended_to_pface_idx[i_part][0] = 0;

    pfull_n_cell_extended[i_part] = 0;
    PDM_malloc(pfull_cell_extended_to_pcell_idx      [i_part], pfull_n_cell_extended[i_part]+1, int);
    PDM_malloc(pfull_cell_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_cell_extended_to_pcell_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_cell_extended_to_pcell_interface[i_part], 0, int        );
    pfull_cell_extended_to_pcell_idx[i_part][0] = 0;
  }

  PDM_g_num_t shift_by_domain_vtx  = part_ext->shift_by_domain_vtx [part_ext->n_domain];
  PDM_g_num_t shift_by_domain_edge = part_ext->shift_by_domain_edge[part_ext->n_domain];
  PDM_g_num_t shift_by_domain_face = part_ext->shift_by_domain_face[part_ext->n_domain];
  PDM_g_num_t shift_by_domain_cell = part_ext->shift_by_domain_cell[part_ext->n_domain];



  /**
   * Cell link between interface
   */
  int          prev_dcell_itrf_n_blk               = 0;
  PDM_g_num_t *prev_dcell_itrf_blk_gnum            = NULL;
  int         *prev_dcell_itrf_blk_ancstr_strd     = NULL;
  PDM_g_num_t *prev_dcell_itrf_blk_ancstr          = NULL;
  int         *prev_dcell_itrf_blk_path_itrf_strd  = NULL;
  int         *prev_dcell_itrf_blk_path_itrf       = NULL;
  int         *prev_dcell_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *prev_dcell_itrf_gnum_and_itrf_data  = NULL;


  /**
   * Face link between interface
   */
  int          prev_dface_itrf_n_blk               = part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_FACE];
  PDM_g_num_t *prev_dface_itrf_blk_gnum            = part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_FACE];
  int         *prev_dface_itrf_blk_ancstr_strd     = PDM_array_zeros_int(prev_dface_itrf_n_blk);
  PDM_g_num_t *prev_dface_itrf_blk_ancstr          = PDM_array_zeros_gnum(0);
  int         *prev_dface_itrf_blk_path_itrf_strd  = PDM_array_zeros_int(prev_dface_itrf_n_blk);
  int         *prev_dface_itrf_blk_path_itrf       = PDM_array_zeros_int(0);
  int         *prev_dface_itrf_gnum_and_itrf_strid = part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE];
  PDM_g_num_t *prev_dface_itrf_gnum_and_itrf_data  = part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_FACE];
  int         *prev_dface_itrf_gnum_and_itrf_sens  = part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_FACE];


  /**
   * Edge link between interface
   */
  int          prev_dedge_itrf_n_blk               = 0;
  PDM_g_num_t *prev_dedge_itrf_blk_gnum            = NULL;
  int         *prev_dedge_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *prev_dedge_itrf_gnum_and_itrf_data  = NULL;
  int         *prev_dedge_itrf_gnum_and_itrf_sens  = NULL;
  int         *prev_dedge_itrf_blk_ancstr_strd     = NULL;
  PDM_g_num_t *prev_dedge_itrf_blk_ancstr          = NULL;
  int         *prev_dedge_itrf_blk_path_itrf_strd  = NULL;
  int         *prev_dedge_itrf_blk_path_itrf       = NULL;
  if(part_ext->have_edge == 1) {
    prev_dedge_itrf_n_blk               = part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_blk_gnum            = part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_strid = part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_data  = part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_sens  = part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_blk_ancstr_strd     = PDM_array_zeros_int(prev_dedge_itrf_n_blk);
    prev_dedge_itrf_blk_ancstr          = PDM_array_zeros_gnum(0);
    prev_dedge_itrf_blk_path_itrf_strd  = PDM_array_zeros_int(prev_dedge_itrf_n_blk);
    prev_dedge_itrf_blk_path_itrf       = PDM_array_zeros_int(0);
  }


  /**
   * Vtx link between interface
   */
  int          prev_dvtx_itrf_n_blk               = part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX];
  PDM_g_num_t *prev_dvtx_itrf_blk_gnum            = part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX];
  int         *prev_dvtx_itrf_gnum_and_itrf_strid = part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX];
  PDM_g_num_t *prev_dvtx_itrf_gnum_and_itrf_data  = part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX];
  int         *prev_dvtx_itrf_blk_ancstr_strd     = PDM_array_zeros_int(prev_dvtx_itrf_n_blk);
  PDM_g_num_t *prev_dvtx_itrf_blk_ancstr          = PDM_array_zeros_gnum(0);
  int         *prev_dvtx_itrf_blk_path_itrf_strd  = PDM_array_zeros_int(prev_dvtx_itrf_n_blk);
  int         *prev_dvtx_itrf_blk_path_itrf       = PDM_array_zeros_int(0);

  int **pcurr_entity_bound_to_pentity_bound_idx       = part_ext->pinit_entity_bound_to_pentity_bound_idx;
  int **pcurr_entity_bound_to_pentity_bound_triplet   = part_ext->pinit_entity_bound_to_pentity_bound_triplet;
  int **pcurr_entity_bound_to_pentity_bound_interface = part_ext->pinit_entity_bound_to_pentity_bound_interface;


  int **pn_vtx_extended_by_depth  = NULL;
  int **pn_edge_extended_by_depth = NULL;
  int **pn_face_extended_by_depth = NULL;
  int **pn_cell_extended_by_depth = NULL;
  PDM_malloc(pn_vtx_extended_by_depth , part_ext->depth, int *);
  PDM_malloc(pn_edge_extended_by_depth, part_ext->depth, int *);
  PDM_malloc(pn_face_extended_by_depth, part_ext->depth, int *);
  PDM_malloc(pn_cell_extended_by_depth, part_ext->depth, int *);
  for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {
    PDM_malloc(pn_vtx_extended_by_depth   [i_depth], part_ext->ln_part_tot, int);
    if(part_ext->have_edge == 1) {
      PDM_malloc(pn_edge_extended_by_depth[i_depth], part_ext->ln_part_tot, int);
    }
    PDM_malloc(pn_face_extended_by_depth  [i_depth], part_ext->ln_part_tot, int);
    PDM_malloc(pn_cell_extended_by_depth  [i_depth], part_ext->ln_part_tot, int);

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      pn_vtx_extended_by_depth [i_depth][i_part] = 0;
      if(part_ext->have_edge == 1) {
        pn_edge_extended_by_depth[i_depth][i_part] = 0;
      }
      pn_face_extended_by_depth[i_depth][i_part] = 0;
      pn_cell_extended_by_depth[i_depth][i_part] = 0;
    }
  }


  if (part_ext->extend_type==PDM_EXTEND_FROM_EDGE && part_ext->have_edge==0) {
    PDM_error(__FILE__, __LINE__, 0, "part_extension with extend_type %d asked but edge seems to miss\n", part_ext->extend_type);
  }



  int i_depth   = 0;
  int step      = 0;
  while(i_depth < part_ext->depth) {


    int i_rank;
    PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
    if (i_rank==0) printf("Computing DEPTH %d (step = %i) \n", i_depth, step);

    if(debug == 1) {
      log_trace("\n\n\n >> DEPTH %d step = %i\n", i_depth, step);
    }
    double t_start = PDM_MPI_Wtime();



    /* Use descending connectivity to deduce connectivity and extend cells */
    int          *pn_cell_extended                  = NULL;
    PDM_g_num_t **pcell_extended_ln_to_gn           = NULL;
    int         **pcell_extended_alrdy_sent         = NULL;
    int         **pcell_extended_to_pcell_idx       = NULL;
    int         **pcell_extended_to_pcell_triplet   = NULL;
    int         **pcell_extended_to_pcell_interface = NULL;


    if (part_ext->extend_type==PDM_EXTEND_FROM_VTX) {
      if(part_ext->have_edge == 1) {
        for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
          PDM_compute_face_vtx_from_face_and_edge(pn_face       [i_part],
                                                  pface_edge_idx[i_part],
                                                  pface_edge    [i_part],
                                                  pedge_vtx     [i_part],
                                                 &pface_vtx     [i_part]);
          pface_vtx_idx[i_part] = pface_edge_idx[i_part];

          PDM_combine_connectivity(pn_cell       [i_part],
                                   pcell_face_idx[i_part],
                                   pcell_face    [i_part],
                                   pface_vtx_idx [i_part],
                                   pface_vtx     [i_part],
                                  &pcell_vtx_idx [i_part],
                                  &pcell_vtx     [i_part]);

          PDM_free(pface_vtx[i_part]);
        }
      }
      else {
        for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
          PDM_combine_connectivity(pn_cell       [i_part],
                                   pcell_face_idx[i_part],
                                   pcell_face    [i_part],
                                   pface_vtx_idx [i_part],
                                   pface_vtx     [i_part],
                                  &pcell_vtx_idx [i_part],
                                  &pcell_vtx     [i_part]);
        }
      }
    }



    /*
     * Hook the most leading entity connexion
     */
    int          next_dcell_itrf_n_blk               = 0;
    PDM_g_num_t *next_dcell_itrf_blk_gnum            = NULL;
    int         *next_dcell_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dcell_itrf_blk_ancstr          = NULL;
    int         *next_dcell_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dcell_itrf_blk_path_itrf       = NULL;
    int         *next_dcell_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dcell_itrf_gnum_and_itrf_data  = NULL;

    if (debug==1) {
      log_trace("\n\n");
      log_trace("========================================= \n");
      log_trace("PDM_part_extension_entity1_to_entity2 beg \n");
    }
    if(part_ext->extend_type==PDM_EXTEND_FROM_VTX) {
      PDM_part_extension_entity1_to_entity2(shift_by_domain_cell,
                                            part_ext->ln_part_tot,
                                            pn_vtx,
                                            pvtx_ln_to_gn,
                                            pcurr_entity_bound_to_pentity_bound_idx,
                                            pcurr_entity_bound_to_pentity_bound_triplet,
                                            pcurr_entity_bound_to_pentity_bound_interface,
                                            pn_cell,
                                            pcell_ln_to_gn,
                                            pcell_alrdy_sent,
                                            pcell_vtx_idx,
                                            pcell_vtx,
                                            prev_dcell_itrf_n_blk,
                                            prev_dcell_itrf_blk_gnum,
                                            prev_dcell_itrf_blk_ancstr_strd,
                                            prev_dcell_itrf_blk_ancstr,
                                            prev_dcell_itrf_blk_path_itrf_strd,
                                            prev_dcell_itrf_blk_path_itrf,
                                            prev_dcell_itrf_gnum_and_itrf_strid,
                                            prev_dcell_itrf_gnum_and_itrf_data,
                                            &pn_cell_extended,
                                            &pcell_extended_ln_to_gn,
                                            &pcell_extended_alrdy_sent,
                                            &pcell_extended_to_pcell_idx,
                                            &pcell_extended_to_pcell_triplet,
                                            &pcell_extended_to_pcell_interface,
                                            &next_dcell_itrf_n_blk,
                                            &next_dcell_itrf_blk_gnum,
                                            &next_dcell_itrf_blk_ancstr_strd,
                                            &next_dcell_itrf_blk_ancstr,
                                            &next_dcell_itrf_blk_path_itrf_strd,
                                            &next_dcell_itrf_blk_path_itrf,
                                            &next_dcell_itrf_gnum_and_itrf_strid,
                                            &next_dcell_itrf_gnum_and_itrf_data,
                                            part_ext->comm);
      
      for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
        PDM_free(pcell_vtx_idx[i_part]);
        PDM_free(pcell_vtx    [i_part]);
      }
    }
    else if (part_ext->extend_type==PDM_EXTEND_FROM_FACE) {
      PDM_part_extension_entity1_to_entity2(shift_by_domain_cell,
                                            part_ext->ln_part_tot,
                                            pn_face,
                                            pface_ln_to_gn,
                                            pcurr_entity_bound_to_pentity_bound_idx,
                                            pcurr_entity_bound_to_pentity_bound_triplet,
                                            pcurr_entity_bound_to_pentity_bound_interface,
                                            pn_cell,
                                            pcell_ln_to_gn,
                                            pcell_alrdy_sent,
                                            pcell_face_idx,
                                            pcell_face,
                                            prev_dcell_itrf_n_blk,
                                            prev_dcell_itrf_blk_gnum,
                                            prev_dcell_itrf_blk_ancstr_strd,
                                            prev_dcell_itrf_blk_ancstr,
                                            prev_dcell_itrf_blk_path_itrf_strd,
                                            prev_dcell_itrf_blk_path_itrf,
                                            prev_dcell_itrf_gnum_and_itrf_strid,
                                            prev_dcell_itrf_gnum_and_itrf_data,
                                            &pn_cell_extended,
                                            &pcell_extended_ln_to_gn,
                                            &pcell_extended_alrdy_sent,
                                            &pcell_extended_to_pcell_idx,
                                            &pcell_extended_to_pcell_triplet,
                                            &pcell_extended_to_pcell_interface,
                                            &next_dcell_itrf_n_blk,
                                            &next_dcell_itrf_blk_gnum,
                                            &next_dcell_itrf_blk_ancstr_strd,
                                            &next_dcell_itrf_blk_ancstr,
                                            &next_dcell_itrf_blk_path_itrf_strd,
                                            &next_dcell_itrf_blk_path_itrf,
                                            &next_dcell_itrf_gnum_and_itrf_strid,
                                            &next_dcell_itrf_gnum_and_itrf_data,
                                            part_ext->comm);
    }

    if (debug==1) {
      log_trace("\n\n");
      log_trace("PDM_part_extension_entity1_to_entity2 end \n");
      log_trace("========================================= \n");
    }

    PDM_free(prev_dcell_itrf_blk_gnum           );
    PDM_free(prev_dcell_itrf_blk_ancstr_strd    );
    PDM_free(prev_dcell_itrf_blk_ancstr         );
    PDM_free(prev_dcell_itrf_blk_path_itrf_strd );
    PDM_free(prev_dcell_itrf_blk_path_itrf      );
    PDM_free(prev_dcell_itrf_gnum_and_itrf_strid);
    PDM_free(prev_dcell_itrf_gnum_and_itrf_data );
    prev_dcell_itrf_n_blk               = next_dcell_itrf_n_blk;
    prev_dcell_itrf_blk_gnum            = next_dcell_itrf_blk_gnum;
    prev_dcell_itrf_blk_ancstr_strd     = next_dcell_itrf_blk_ancstr_strd;
    prev_dcell_itrf_blk_ancstr          = next_dcell_itrf_blk_ancstr;
    prev_dcell_itrf_blk_path_itrf_strd  = next_dcell_itrf_blk_path_itrf_strd;
    prev_dcell_itrf_blk_path_itrf       = next_dcell_itrf_blk_path_itrf;
    prev_dcell_itrf_gnum_and_itrf_strid = next_dcell_itrf_gnum_and_itrf_strid;
    prev_dcell_itrf_gnum_and_itrf_data  = next_dcell_itrf_gnum_and_itrf_data;


    /*
     * Update with descending connectivity :
     */

    // > Vertex to vertex link
    int          *pn_vtx_extended                   = NULL;
    PDM_g_num_t **pvtx_extended_ln_to_gn            = NULL;
    int         **pvtx_extended_to_pvtx_idx         = NULL;
    int         **pvtx_extended_to_pvtx_triplet     = NULL;
    int         **pvtx_extended_to_pvtx_interface   = NULL;
    
    // > Edge to edge link
    int          *pn_edge_extended                  = NULL;
    PDM_g_num_t **pedge_extended_ln_to_gn           = NULL;
    int         **pedge_extended_to_pedge_idx       = NULL;
    int         **pedge_extended_to_pedge_triplet   = NULL;
    int         **pedge_extended_to_pedge_interface = NULL;
    int         **pedge_extended_to_pedge_sens      = NULL;
    
    // > Face to face link
    int          *pn_face_extended                  = NULL;
    PDM_g_num_t **pface_extended_ln_to_gn           = NULL;
    int         **pface_extended_to_pface_idx       = NULL;
    int         **pface_extended_to_pface_triplet   = NULL;
    int         **pface_extended_to_pface_interface = NULL;
    int         **pface_extended_to_pface_sens      = NULL;
    
    // > Extended connectivities
    int         **pextended_cell_face_idx           = NULL;
    int         **pextended_cell_face               = NULL;
    int         **pextended_face_edge_idx           = NULL;
    int         **pextended_face_edge               = NULL;
    int         **pextended_edge_vtx_idx            = NULL;
    int         **pextended_edge_vtx                = NULL;
    int         **pextended_face_vtx_idx            = NULL;
    int         **pextended_face_vtx                = NULL;

    // > Vertices next data_base 
    int          next_dvtx_itrf_n_blk               = 0;
    PDM_g_num_t *next_dvtx_itrf_blk_gnum            = NULL;
    int         *next_dvtx_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dvtx_itrf_blk_ancstr          = NULL;
    int         *next_dvtx_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dvtx_itrf_blk_path_itrf       = NULL;
    int         *next_dvtx_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dvtx_itrf_gnum_and_itrf_data  = NULL;

    // > Edges next data_base 
    int          next_dedge_itrf_n_blk               = 0;
    PDM_g_num_t *next_dedge_itrf_blk_gnum            = NULL;
    int         *next_dedge_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dedge_itrf_blk_ancstr          = NULL;
    int         *next_dedge_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dedge_itrf_blk_path_itrf       = NULL;
    int         *next_dedge_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dedge_itrf_gnum_and_itrf_data  = NULL;
    int         *next_dedge_itrf_gnum_and_itrf_sens  = NULL;

    // > Faces next data_base 
    int          next_dface_itrf_n_blk               = 0;
    PDM_g_num_t *next_dface_itrf_blk_gnum            = NULL;
    int         *next_dface_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dface_itrf_blk_ancstr          = NULL;
    int         *next_dface_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dface_itrf_blk_path_itrf       = NULL;
    int         *next_dface_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dface_itrf_gnum_and_itrf_data  = NULL;
    int         *next_dface_itrf_gnum_and_itrf_sens  = NULL;


    if (debug==1) {
      log_trace("\n\n\n\n");
      log_trace("=================================================================================\n");
      log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Cell->Face)\n");
    }
    PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                     part_ext->n_interface,
                                                                     shift_by_domain_face,
                                                                     prev_dface_itrf_n_blk,
                                                                     prev_dface_itrf_blk_gnum,
                                                                     prev_dface_itrf_blk_ancstr_strd,
                                                                     prev_dface_itrf_blk_ancstr,
                                                                     prev_dface_itrf_blk_path_itrf_strd,
                                                                     prev_dface_itrf_blk_path_itrf,
                                                                     prev_dface_itrf_gnum_and_itrf_strid,
                                                                     prev_dface_itrf_gnum_and_itrf_data,
                                                                     prev_dface_itrf_gnum_and_itrf_sens,
                                                                     pn_cell,
                                                                     pcell_ln_to_gn,
                                                                     pn_face,
                                                                     pface_ln_to_gn,
                                                                     pcell_face_idx,
                                                                     pcell_face,
                                                                     pn_cell_extended,
                                                                     pcell_extended_ln_to_gn,
                                                                     pcell_extended_to_pcell_idx,
                                                                     pcell_extended_to_pcell_triplet,
                                                                     pcell_extended_to_pcell_interface,
                                                                     NULL,
                                                                     1,
                                                                     &pn_face_extended,
                                                                     &pface_extended_ln_to_gn,
                                                                     &pextended_cell_face_idx,
                                                                     &pextended_cell_face,
                                                                     &pface_extended_to_pface_idx,
                                                                     &pface_extended_to_pface_triplet,
                                                                     &pface_extended_to_pface_interface,
                                                                     &pface_extended_to_pface_sens,
                                                                     &next_dface_itrf_n_blk,
                                                                     &next_dface_itrf_blk_gnum,
                                                                     &next_dface_itrf_blk_ancstr_strd,
                                                                     &next_dface_itrf_blk_ancstr,
                                                                     &next_dface_itrf_blk_path_itrf_strd,
                                                                     &next_dface_itrf_blk_path_itrf,
                                                                     &next_dface_itrf_gnum_and_itrf_strid,
                                                                     &next_dface_itrf_gnum_and_itrf_data,
                                                                     &next_dface_itrf_gnum_and_itrf_sens,
                                                                     part_ext->comm);
    if (debug==1) {
      log_trace("\n\n");
      log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Cell->Face)\n");
      log_trace("=================================================================================\n");
    }

    PDM_free(prev_dface_itrf_blk_gnum           );
    PDM_free(prev_dface_itrf_blk_ancstr_strd    );
    PDM_free(prev_dface_itrf_blk_ancstr         );
    PDM_free(prev_dface_itrf_blk_path_itrf_strd );
    PDM_free(prev_dface_itrf_blk_path_itrf      );
    PDM_free(prev_dface_itrf_gnum_and_itrf_strid);
    PDM_free(prev_dface_itrf_gnum_and_itrf_data );
    PDM_free(prev_dface_itrf_gnum_and_itrf_sens );
    prev_dface_itrf_n_blk               = next_dface_itrf_n_blk;
    prev_dface_itrf_blk_gnum            = next_dface_itrf_blk_gnum;
    prev_dface_itrf_blk_ancstr_strd     = next_dface_itrf_blk_ancstr_strd;
    prev_dface_itrf_blk_ancstr          = next_dface_itrf_blk_ancstr;
    prev_dface_itrf_blk_path_itrf_strd  = next_dface_itrf_blk_path_itrf_strd;
    prev_dface_itrf_blk_path_itrf       = next_dface_itrf_blk_path_itrf;
    prev_dface_itrf_gnum_and_itrf_strid = next_dface_itrf_gnum_and_itrf_strid;
    prev_dface_itrf_gnum_and_itrf_data  = next_dface_itrf_gnum_and_itrf_data;
    prev_dface_itrf_gnum_and_itrf_sens  = next_dface_itrf_gnum_and_itrf_sens;

    part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_FACE] = prev_dface_itrf_n_blk;
    part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_FACE] = prev_dface_itrf_blk_gnum;
    part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_FACE] = prev_dface_itrf_gnum_and_itrf_strid;
    part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_FACE] = prev_dface_itrf_gnum_and_itrf_data;
    part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_FACE] = prev_dface_itrf_gnum_and_itrf_sens;


    if(part_ext->have_edge == 1) {
      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("=================================================================================\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Face->Edge)\n");
      }
      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_edge,
                                                                       prev_dedge_itrf_n_blk,
                                                                       prev_dedge_itrf_blk_gnum,
                                                                       prev_dedge_itrf_blk_ancstr_strd,
                                                                       prev_dedge_itrf_blk_ancstr,
                                                                       prev_dedge_itrf_blk_path_itrf_strd,
                                                                       prev_dedge_itrf_blk_path_itrf,
                                                                       prev_dedge_itrf_gnum_and_itrf_strid,
                                                                       prev_dedge_itrf_gnum_and_itrf_data,
                                                                       prev_dedge_itrf_gnum_and_itrf_sens,
                                                                       pn_face,
                                                                       pface_ln_to_gn,
                                                                       pn_edge,
                                                                       pedge_ln_to_gn,
                                                                       pface_edge_idx,
                                                                       pface_edge,
                                                                       pn_face_extended,
                                                                       pface_extended_ln_to_gn,
                                                                       pface_extended_to_pface_idx,
                                                                       pface_extended_to_pface_triplet,
                                                                       pface_extended_to_pface_interface,
                                                                       pface_extended_to_pface_sens,
                                                                       1,
                                                                       &pn_edge_extended,
                                                                       &pedge_extended_ln_to_gn,
                                                                       &pextended_face_edge_idx,
                                                                       &pextended_face_edge,
                                                                       &pedge_extended_to_pedge_idx,
                                                                       &pedge_extended_to_pedge_triplet,
                                                                       &pedge_extended_to_pedge_interface,
                                                                       &pedge_extended_to_pedge_sens,
                                                                       &next_dedge_itrf_n_blk,
                                                                       &next_dedge_itrf_blk_gnum,
                                                                       &next_dedge_itrf_blk_ancstr_strd,
                                                                       &next_dedge_itrf_blk_ancstr,
                                                                       &next_dedge_itrf_blk_path_itrf_strd,
                                                                       &next_dedge_itrf_blk_path_itrf,
                                                                       &next_dedge_itrf_gnum_and_itrf_strid,
                                                                       &next_dedge_itrf_gnum_and_itrf_data,
                                                                       &next_dedge_itrf_gnum_and_itrf_sens,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Face->Edge)\n");
        log_trace("=================================================================================\n");
      }

      PDM_free(prev_dedge_itrf_blk_gnum);
      PDM_free(prev_dedge_itrf_blk_ancstr_strd);
      PDM_free(prev_dedge_itrf_blk_ancstr);
      PDM_free(prev_dedge_itrf_blk_path_itrf_strd);
      PDM_free(prev_dedge_itrf_blk_path_itrf);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_data);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_sens);
      prev_dedge_itrf_n_blk               = next_dedge_itrf_n_blk;
      prev_dedge_itrf_blk_gnum            = next_dedge_itrf_blk_gnum;
      prev_dedge_itrf_blk_ancstr_strd     = next_dedge_itrf_blk_ancstr_strd;
      prev_dedge_itrf_blk_ancstr          = next_dedge_itrf_blk_ancstr;
      prev_dedge_itrf_blk_path_itrf_strd  = next_dedge_itrf_blk_path_itrf_strd;
      prev_dedge_itrf_blk_path_itrf       = next_dedge_itrf_blk_path_itrf;
      prev_dedge_itrf_gnum_and_itrf_strid = next_dedge_itrf_gnum_and_itrf_strid;
      prev_dedge_itrf_gnum_and_itrf_data  = next_dedge_itrf_gnum_and_itrf_data;
      prev_dedge_itrf_gnum_and_itrf_sens  = next_dedge_itrf_gnum_and_itrf_sens;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_blk_gnum;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_data;
      part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_sens;


      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("================================================================================\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Edge->Vtx)\n");
      }

      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_vtx, // Attention il va evoluer lui
                                                                       prev_dvtx_itrf_n_blk,
                                                                       prev_dvtx_itrf_blk_gnum,
                                                                       prev_dvtx_itrf_blk_ancstr_strd,
                                                                       prev_dvtx_itrf_blk_ancstr,
                                                                       prev_dvtx_itrf_blk_path_itrf_strd,
                                                                       prev_dvtx_itrf_blk_path_itrf,
                                                                       prev_dvtx_itrf_gnum_and_itrf_strid,
                                                                       prev_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       pn_edge,
                                                                       pedge_ln_to_gn,
                                                                       pn_vtx,
                                                                       pvtx_ln_to_gn,
                                                                       pedge_vtx_idx,
                                                                       pedge_vtx,
                                                                       pn_edge_extended,
                                                                       pedge_extended_ln_to_gn,
                                                                       pedge_extended_to_pedge_idx,
                                                                       pedge_extended_to_pedge_triplet,
                                                                       pedge_extended_to_pedge_interface,
                                                                       pedge_extended_to_pedge_sens,
                                                                       1,
                                                                      &pn_vtx_extended,
                                                                      &pvtx_extended_ln_to_gn,
                                                                      &pextended_edge_vtx_idx,
                                                                      &pextended_edge_vtx,
                                                                      &pvtx_extended_to_pvtx_idx,
                                                                      &pvtx_extended_to_pvtx_triplet,
                                                                      &pvtx_extended_to_pvtx_interface,
                                                                       NULL,
                                                                      &next_dvtx_itrf_n_blk,
                                                                      &next_dvtx_itrf_blk_gnum,
                                                                      &next_dvtx_itrf_blk_ancstr_strd,
                                                                      &next_dvtx_itrf_blk_ancstr,
                                                                      &next_dvtx_itrf_blk_path_itrf_strd,
                                                                      &next_dvtx_itrf_blk_path_itrf,
                                                                      &next_dvtx_itrf_gnum_and_itrf_strid,
                                                                      &next_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Edge->Vtx)\n");
        log_trace("================================================================================\n");
      }

      PDM_free(prev_dvtx_itrf_blk_gnum);
      PDM_free(prev_dvtx_itrf_blk_ancstr_strd);
      PDM_free(prev_dvtx_itrf_blk_ancstr);
      PDM_free(prev_dvtx_itrf_blk_path_itrf_strd);
      PDM_free(prev_dvtx_itrf_blk_path_itrf);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_data);
      prev_dvtx_itrf_n_blk               = next_dvtx_itrf_n_blk;
      prev_dvtx_itrf_blk_gnum            = next_dvtx_itrf_blk_gnum;
      prev_dvtx_itrf_blk_ancstr_strd     = next_dvtx_itrf_blk_ancstr_strd;
      prev_dvtx_itrf_blk_ancstr          = next_dvtx_itrf_blk_ancstr;
      prev_dvtx_itrf_blk_path_itrf_strd  = next_dvtx_itrf_blk_path_itrf_strd;
      prev_dvtx_itrf_blk_path_itrf       = next_dvtx_itrf_blk_path_itrf;
      prev_dvtx_itrf_gnum_and_itrf_strid = next_dvtx_itrf_gnum_and_itrf_strid;
      prev_dvtx_itrf_gnum_and_itrf_data  = next_dvtx_itrf_gnum_and_itrf_data;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_gnum;
      // part_ext->dentity_itrf_blk_ancstr         [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_ancstr;
      // part_ext->dentity_itrf_blk_path_itrf_strd [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf_strd;
      // part_ext->dentity_itrf_blk_path_itrf      [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_data;

    }
    else {

      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("================================================================================ \n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Face->Vtx) \n");
      }

      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_vtx, // Attention il va evoluer lui
                                                                       prev_dvtx_itrf_n_blk,
                                                                       prev_dvtx_itrf_blk_gnum,
                                                                       prev_dvtx_itrf_blk_ancstr_strd,
                                                                       prev_dvtx_itrf_blk_ancstr,
                                                                       prev_dvtx_itrf_blk_path_itrf_strd,
                                                                       prev_dvtx_itrf_blk_path_itrf,
                                                                       prev_dvtx_itrf_gnum_and_itrf_strid,
                                                                       prev_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       pn_face,
                                                                       pface_ln_to_gn,
                                                                       pn_vtx,
                                                                       pvtx_ln_to_gn,
                                                                       pface_vtx_idx,
                                                                       pface_vtx,
                                                                       pn_face_extended,
                                                                       pface_extended_ln_to_gn,
                                                                       pface_extended_to_pface_idx,
                                                                       pface_extended_to_pface_triplet,
                                                                       pface_extended_to_pface_interface,
                                                                       pface_extended_to_pface_sens,
                                                                       1,
                                                                      &pn_vtx_extended,
                                                                      &pvtx_extended_ln_to_gn,
                                                                      &pextended_face_vtx_idx,
                                                                      &pextended_face_vtx,
                                                                      &pvtx_extended_to_pvtx_idx,
                                                                      &pvtx_extended_to_pvtx_triplet,
                                                                      &pvtx_extended_to_pvtx_interface,
                                                                       NULL,
                                                                      &next_dvtx_itrf_n_blk,
                                                                      &next_dvtx_itrf_blk_gnum,
                                                                      &next_dvtx_itrf_blk_ancstr_strd,
                                                                      &next_dvtx_itrf_blk_ancstr,
                                                                      &next_dvtx_itrf_blk_path_itrf_strd,
                                                                      &next_dvtx_itrf_blk_path_itrf,
                                                                      &next_dvtx_itrf_gnum_and_itrf_strid,
                                                                      &next_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Face->Vtx) \n");
        log_trace("================================================================================ \n");
      }

      PDM_free(prev_dvtx_itrf_blk_gnum           );
      PDM_free(prev_dvtx_itrf_blk_ancstr_strd    );
      PDM_free(prev_dvtx_itrf_blk_ancstr         );
      PDM_free(prev_dvtx_itrf_blk_path_itrf_strd );
      PDM_free(prev_dvtx_itrf_blk_path_itrf      );
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_data );
      prev_dvtx_itrf_n_blk               = next_dvtx_itrf_n_blk;
      prev_dvtx_itrf_blk_gnum            = next_dvtx_itrf_blk_gnum;
      prev_dvtx_itrf_blk_ancstr_strd     = next_dvtx_itrf_blk_ancstr_strd;
      prev_dvtx_itrf_blk_ancstr          = next_dvtx_itrf_blk_ancstr;
      prev_dvtx_itrf_blk_path_itrf_strd  = next_dvtx_itrf_blk_path_itrf_strd;
      prev_dvtx_itrf_blk_path_itrf       = next_dvtx_itrf_blk_path_itrf;
      prev_dvtx_itrf_gnum_and_itrf_strid = next_dvtx_itrf_gnum_and_itrf_strid;
      prev_dvtx_itrf_gnum_and_itrf_data  = next_dvtx_itrf_gnum_and_itrf_data;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_gnum;
      // part_ext->dentity_itrf_blk_ancstr         [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_ancstr;
      // part_ext->dentity_itrf_blk_path_itrf_strd [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf_strd;
      // part_ext->dentity_itrf_blk_path_itrf      [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_data;


    }




    /*
     * Hook coordinates
     */
    double **pvtx_extended_coords = NULL;
    _exchange_coord_and_apply_transform(part_ext,
                                        pn_vtx_extended,
                                        pvtx_extended_ln_to_gn,
                                        pn_vtx,
                                        pvtx_coords,
                                        pvtx_extended_to_pvtx_idx,
                                        pvtx_extended_to_pvtx_triplet,
                                        pvtx_extended_to_pvtx_interface,
                                        &pvtx_extended_coords);

    /*
     * Concatenate all information to continue recursion
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      if(debug == 1) {
        log_trace("  -> i_part = %d \n", i_part);
      }

      /* Update size */
      pn_vtx_extended_old              [i_part]  = pfull_n_vtx_extended[i_part];
      pfull_n_vtx_extended             [i_part] +=      pn_vtx_extended[i_part];
      pn_vtx_extended_by_depth[i_depth][i_part] +=      pn_vtx_extended[i_part];
      
      if(part_ext->have_edge == 1) {
        pn_edge_extended_old              [i_part]  = pfull_n_edge_extended[i_part];
        pfull_n_edge_extended             [i_part] +=      pn_edge_extended[i_part];
        pn_edge_extended_by_depth[i_depth][i_part] +=      pn_edge_extended[i_part];
      }

      pn_face_extended_old              [i_part]  = pfull_n_face_extended[i_part];
      pfull_n_face_extended             [i_part] +=      pn_face_extended[i_part];
      pn_face_extended_by_depth[i_depth][i_part] +=      pn_face_extended[i_part];

      pn_cell_extended_old              [i_part]  = pfull_n_cell_extended[i_part];
      pfull_n_cell_extended             [i_part] +=      pn_cell_extended[i_part];
      pn_cell_extended_by_depth[i_depth][i_part] +=      pn_cell_extended[i_part];

    }





    /* Cells */
    _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                           pn_cell,
                                           pn_cell_extended,
                                           pcell_ln_to_gn,
                                           pcell_extended_ln_to_gn,
                                           &shift_by_domain_cell);

    _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                       pn_cell,
                                       pn_cell_extended,
                                       pcell_face_idx,
                                       pcell_face,
                                       pextended_cell_face_idx,
                                       pextended_cell_face);

    _concat_int_array_current_with_extended(part_ext->ln_part_tot,
                                            pn_cell,
                                            pn_cell_extended,
                                            pcell_alrdy_sent,
                                            pcell_extended_alrdy_sent);

    /* Faces */
    _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                           pn_face,
                                           pn_face_extended,
                                           pface_ln_to_gn,
                                           pface_extended_ln_to_gn,
                                           &shift_by_domain_face);

    if(part_ext->have_edge == 1) {
      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_face,
                                         pn_face_extended,
                                         pface_edge_idx,
                                         pface_edge,
                                         pextended_face_edge_idx,
                                         pextended_face_edge);
    }
    else {
      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_face,
                                         pn_face_extended,
                                         pface_vtx_idx,
                                         pface_vtx,
                                         pextended_face_vtx_idx,
                                         pextended_face_vtx);
    }

    /* Edges */
    if(part_ext->have_edge == 1) {

      _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                             pn_edge,
                                             pn_edge_extended,
                                             pedge_ln_to_gn,
                                             pedge_extended_ln_to_gn,
                                             &shift_by_domain_edge);

      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_edge,
                                         pn_edge_extended,
                                         pedge_vtx_idx,
                                         pedge_vtx,
                                         pextended_edge_vtx_idx,
                                         pextended_edge_vtx);
    }


    /* Vertices */
    _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                           pn_vtx,
                                           pn_vtx_extended,
                                           pvtx_ln_to_gn,
                                           pvtx_extended_ln_to_gn,
                                           &shift_by_domain_vtx);

    _concat_coords_current_with_extended(part_ext->ln_part_tot,
                                         pn_vtx,
                                         pn_vtx_extended,
                                         pvtx_coords,
                                         pvtx_extended_coords);

    _concat_full_with_extended(part_ext->ln_part_tot,
                               pfull_n_vtx_extended,
                               pn_vtx_extended,
                               pn_vtx_extended_old,
                               pvtx_extended_ln_to_gn,
                               pvtx_extended_to_pvtx_idx,
                               pvtx_extended_to_pvtx_triplet,
                               pvtx_extended_to_pvtx_interface,
                               pfull_vtx_extended_ln_to_gn,
                               pfull_vtx_extended_to_pvtx_idx,
                               pfull_vtx_extended_to_pvtx_triplet,
                               pfull_vtx_extended_to_pvtx_interface);

    if (part_ext->have_edge==1) {
      _concat_full_with_extended(part_ext->ln_part_tot,
                                 pfull_n_edge_extended,
                                 pn_edge_extended,
                                 pn_edge_extended_old,
                                 pedge_extended_ln_to_gn,
                                 pedge_extended_to_pedge_idx,
                                 pedge_extended_to_pedge_triplet,
                                 pedge_extended_to_pedge_interface,
                                 pfull_edge_extended_ln_to_gn,
                                 pfull_edge_extended_to_pedge_idx,
                                 pfull_edge_extended_to_pedge_triplet,
                                 pfull_edge_extended_to_pedge_interface);
    }

    _concat_full_with_extended(part_ext->ln_part_tot,
                               pfull_n_face_extended,
                               pn_face_extended,
                               pn_face_extended_old,
                               pface_extended_ln_to_gn,
                               pface_extended_to_pface_idx,
                               pface_extended_to_pface_triplet,
                               pface_extended_to_pface_interface,
                               pfull_face_extended_ln_to_gn,
                               pfull_face_extended_to_pface_idx,
                               pfull_face_extended_to_pface_triplet,
                               pfull_face_extended_to_pface_interface);

    _concat_full_with_extended(part_ext->ln_part_tot,
                               pfull_n_cell_extended,
                               pn_cell_extended,
                               pn_cell_extended_old,
                               pcell_extended_ln_to_gn,
                               pcell_extended_to_pcell_idx,
                               pcell_extended_to_pcell_triplet,
                               pcell_extended_to_pcell_interface,
                               pfull_cell_extended_ln_to_gn,
                               pfull_cell_extended_to_pcell_idx,
                               pfull_cell_extended_to_pcell_triplet,
                               pfull_cell_extended_to_pcell_interface);




    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      /*   */
      int pn_concat_vtx  = pn_vtx [i_part] + pn_vtx_extended [i_part];
      int pn_concat_edge = pn_edge[i_part];
      if (part_ext->have_edge==1) {
        pn_concat_edge += pn_edge_extended[i_part];
      }
      int pn_concat_face = pn_face[i_part] + pn_face_extended[i_part];
      int pn_concat_cell = pn_cell[i_part] + pn_cell_extended[i_part];

      int size_vtx_vtx   = (pfull_vtx_extended_to_pvtx_idx  [i_part][pn_vtx_extended_old [i_part]] + pvtx_extended_to_pvtx_idx  [i_part][pn_vtx_extended [i_part]])/3;
      int size_edge_edge = 0;
      if (part_ext->have_edge==1) {
        size_edge_edge = (pfull_edge_extended_to_pedge_idx[i_part][pn_edge_extended_old[i_part]] + pedge_extended_to_pedge_idx[i_part][pn_edge_extended[i_part]])/3;
      }
      int size_face_face = (pfull_face_extended_to_pface_idx[i_part][pn_face_extended_old[i_part]] + pface_extended_to_pface_idx[i_part][pn_face_extended[i_part]])/3;
      int size_cell_cell = (pfull_cell_extended_to_pcell_idx[i_part][pn_cell_extended_old[i_part]] + pcell_extended_to_pcell_idx[i_part][pn_cell_extended[i_part]])/3;

      if(debug == 1) {
        log_trace("\n");
        log_trace("i_part = %d\n", i_part);
        log_trace("Vertices connection :: \n");
        PDM_log_trace_array_long(pfull_vtx_extended_ln_to_gn           [i_part], pfull_n_vtx_extended[i_part]  , "pfull_vtx_extended_ln_to_gn         ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_idx        [i_part], pfull_n_vtx_extended[i_part]+1, "pfull_vtx_extended_to_pvtx_idx      ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_triplet    [i_part], 3 * size_vtx_vtx              , "pfull_vtx_extended_to_pvtx_triplet  ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_interface  [i_part],     size_vtx_vtx              , "pfull_vtx_extended_to_pvtx_interface::");
        if (part_ext->have_edge==1) {
          log_trace("Edges connection :: \n");
          PDM_log_trace_array_long(pfull_edge_extended_ln_to_gn            [i_part], pfull_n_edge_extended[i_part]  , "pfull_edge_extended_ln_to_gn         ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_idx        [i_part], pfull_n_edge_extended[i_part]+1, "pfull_edge_extended_to_pedge_idx      ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_triplet    [i_part], 3 * size_edge_edge             , "pfull_edge_extended_to_pedge_triplet  ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_interface  [i_part],     size_edge_edge             , "pfull_edge_extended_to_pedge_interface::");
        }
        log_trace("Faces connection :: \n");
        PDM_log_trace_array_long(pfull_face_extended_ln_to_gn          [i_part], pfull_n_face_extended[i_part]  , "pfull_face_extended_ln_to_gn           ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_idx      [i_part], pfull_n_face_extended[i_part]+1, "pfull_face_extended_to_pface_idx       ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_triplet  [i_part], 3 * size_face_face             , "pfull_face_extended_to_pface_triplet   ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_interface[i_part],     size_face_face             , "pfull_face_extended_to_pface_interface ::");
        log_trace("Cells connection :: \n");
        PDM_log_trace_array_long(pfull_cell_extended_ln_to_gn          [i_part], pfull_n_cell_extended[i_part]  , "pfull_cell_extended_ln_to_gn           ::");
        PDM_log_trace_array_int (pfull_cell_extended_to_pcell_idx      [i_part], pfull_n_cell_extended[i_part]+1, "pfull_cell_extended_to_pcell_idx       ::");
        PDM_log_trace_array_int (pfull_cell_extended_to_pcell_triplet  [i_part], 3 * size_cell_cell             , "pfull_cell_extended_to_pcell_triplet   ::");
        PDM_log_trace_array_int (pfull_cell_extended_to_pcell_interface[i_part],     size_cell_cell             , "pfull_cell_extended_to_pcell_interface ::");
      }


      // printf("part_ext->extend_type = %i \n", part_ext->extend_type);

      /*
       * Update graphe - Only extend idx since no connextion is create AT this stage
       */
      int pn_concat_entity_extended = 0;
      int pn_entity_extended        = 0;
      int pn_entity                 = 0;
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
        pn_concat_entity_extended = pn_concat_vtx;
        pn_entity_extended        = pn_vtx_extended[i_part];
        pn_entity                 = pn_vtx         [i_part];
      } else if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        pn_concat_entity_extended = pn_concat_face;
        pn_entity_extended        = pn_face_extended[i_part];
        pn_entity                 = pn_face         [i_part];
      } else {
        PDM_error(__FILE__, __LINE__, 0, "part_extension with extend_type %d for 3d mesh\n", part_ext->extend_type);
      }
      // log_trace("pn_concat_entity_extended = %i \n", pn_concat_entity_extended);

      /*
       * Only extend the index array since connexion is freeze for one step
       */
      PDM_realloc(pcurr_entity_bound_to_pentity_bound_idx[i_part], pcurr_entity_bound_to_pentity_bound_idx[i_part], pn_concat_entity_extended+1, int);
      for(int i = 0; i < pn_entity_extended; ++i) {
        pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity];
      }




      if(debug == 1) {
        log_trace("\n");
        log_trace("i_part = %d\n", i_part);
        log_trace("pn_concat_vtx = %d\n", pn_concat_vtx);
        PDM_log_trace_array_long(pvtx_ln_to_gn[i_part], pn_concat_vtx, "pvtx_ln_to_gn ::");

        log_trace("\n");
        log_trace("CELL_FACE\n");
        int cell_face_size = pcell_face_idx[i_part][pn_concat_cell];
        log_trace("pn_concat_cell = %d\n", pn_concat_cell);
        PDM_log_trace_array_int (pcell_face_idx[i_part], pn_concat_cell+1, "pcell_face_idx ::");
        PDM_log_trace_array_int (pcell_face    [i_part], cell_face_size  , "pcell_face     ::");
        PDM_log_trace_array_long(pcell_ln_to_gn[i_part], pn_concat_cell  , "pcell_ln_to_gn ::");

        if (part_ext->have_edge==1) {
          log_trace("\n");
          log_trace("FACE_EDGE\n");
          int face_edge_size = pface_edge_idx[i_part][pn_concat_face];
          log_trace("pn_concat_face = %d\n", pn_concat_face);
          PDM_log_trace_array_int (pface_edge_idx[i_part], pn_concat_face+1, "pface_edge_idx ::");
          PDM_log_trace_array_int (pface_edge    [i_part], face_edge_size  , "pface_edge     ::");
          PDM_log_trace_array_long(pface_ln_to_gn[i_part], pn_concat_face  , "pface_ln_to_gn ::");
          
          log_trace("\n");
          log_trace("EDGE_VTX\n");
          int edge_vtx_size = pedge_vtx_idx[i_part][pn_concat_edge];
          log_trace("pn_concat_edge = %d\n", pn_concat_edge);
          PDM_log_trace_array_int (pedge_vtx_idx [i_part], pn_concat_edge+1, "pedge_vtx_idx  ::");
          PDM_log_trace_array_int (pedge_vtx     [i_part], edge_vtx_size   , "pedge_vtx      ::");
          PDM_log_trace_array_long(pedge_ln_to_gn[i_part], pn_concat_edge  , "pedge_ln_to_gn ::");
        } 
        else {
          log_trace("\n");
          log_trace("FACE_VTX\n");
          int face_vtx_size = pface_vtx_idx[i_part][pn_concat_face];
          log_trace("pn_concat_face = %d\n", pn_concat_face);
          PDM_log_trace_array_int (pface_vtx_idx [i_part], pn_concat_face+1, "pface_vtx_idx  ::");
          PDM_log_trace_array_int (pface_vtx     [i_part], face_vtx_size   , "pface_vtx      ::");
          PDM_log_trace_array_long(pface_ln_to_gn[i_part], pn_concat_face  , "pface_ln_to_gn ::");
        }
      }




      if(visu == 1) {
        char filename[999];

        int *_pface_vtx_idx = NULL;
        int *_pface_vtx     = NULL;
        if (part_ext->have_edge==1) {
          sprintf(filename, "out_edge_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
          PDM_vtk_write_std_elements(filename,
                                     pn_concat_vtx,
                                     pvtx_coords   [i_part],
                                     pvtx_ln_to_gn [i_part],
                                     PDM_MESH_NODAL_BAR2,
                                     pn_concat_edge,
                                     pedge_vtx     [i_part],
                                     pedge_ln_to_gn[i_part],
                                     0,
                                     NULL,
                                     NULL);
          PDM_compute_face_vtx_from_face_and_edge(pn_concat_face,
                                                  pface_edge_idx[i_part],
                                                  pface_edge    [i_part],
                                                  pedge_vtx     [i_part],
                                                &_pface_vtx);
          _pface_vtx_idx = pface_edge_idx[i_part];

        }
        else {
          _pface_vtx_idx = pface_vtx_idx[i_part];
          _pface_vtx     = pface_vtx    [i_part];
        }


        sprintf(filename, "out_face_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               pn_concat_vtx,
                               pvtx_coords   [i_part],
                               pvtx_ln_to_gn [i_part],
                               pn_concat_face,
                              _pface_vtx_idx,
                              _pface_vtx,
                               pface_ln_to_gn[i_part],
                               NULL);


        int *cell_vtx = NULL;
        _hexa_ngon_to_nodal(pn_concat_cell,
                            pcell_face_idx[i_part],
                            pcell_face    [i_part],
                            _pface_vtx_idx,
                            _pface_vtx,
                           &cell_vtx);
        sprintf(filename, "out_cell_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
        PDM_vtk_write_std_elements(filename,
                                   pn_concat_vtx,
                                   pvtx_coords   [i_part],
                                   pvtx_ln_to_gn [i_part],
                                   PDM_MESH_NODAL_HEXA8,
                                   pn_concat_cell,
                                   cell_vtx,
                                   pcell_ln_to_gn[i_part],
                                   0,
                                   NULL,
                                   NULL);

        if (part_ext->have_edge==1) {
          PDM_free(_pface_vtx);
        }
        PDM_free(cell_vtx);
      } // End visu



    } /* End loop part */



    /*
     * Free coords
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(pvtx_extended_coords[i_part]);
    }
    PDM_free(pvtx_extended_coords);

    /*
     * Update shift_by_domain_***
     */
    PDM_g_num_t _shift_by_domain_cell = shift_by_domain_cell;
    PDM_MPI_Allreduce(&_shift_by_domain_cell, &shift_by_domain_cell, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    PDM_g_num_t _shift_by_domain_face = shift_by_domain_face;
    PDM_MPI_Allreduce(&_shift_by_domain_face, &shift_by_domain_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    if (part_ext->have_edge==1) {
      PDM_g_num_t _shift_by_domain_edge = shift_by_domain_edge;
      PDM_MPI_Allreduce(&_shift_by_domain_edge, &shift_by_domain_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
    }

    PDM_g_num_t _shift_by_domain_vtx = shift_by_domain_vtx;
    PDM_MPI_Allreduce(&_shift_by_domain_vtx, &shift_by_domain_vtx, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    PDM_g_num_t _pn_cell_extended_tot = 0;
    PDM_g_num_t  pn_cell_extended_tot = 0;
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      _pn_cell_extended_tot += pn_cell_extended[i_part];
    }

    PDM_MPI_Allreduce(&_pn_cell_extended_tot, &pn_cell_extended_tot, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, part_ext->comm);

    if(debug == 1) {
      log_trace("pn_cell_extended_tot = %i (local = %i ) \n", pn_cell_extended_tot, _pn_cell_extended_tot);
    }


    /*
     * Free
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(pvtx_extended_ln_to_gn         [i_part]);
      PDM_free(pvtx_extended_to_pvtx_idx      [i_part]);
      PDM_free(pvtx_extended_to_pvtx_triplet  [i_part]);
      PDM_free(pvtx_extended_to_pvtx_interface[i_part]);

      PDM_free(pextended_cell_face_idx[i_part]);
      PDM_free(pextended_cell_face    [i_part]);

      if (part_ext->have_edge==1) {
        PDM_free(pedge_extended_ln_to_gn          [i_part]);
        PDM_free(pedge_extended_to_pedge_idx      [i_part]);
        PDM_free(pedge_extended_to_pedge_triplet  [i_part]);
        PDM_free(pedge_extended_to_pedge_interface[i_part]);
        if (pedge_extended_to_pedge_sens!=NULL) {
          PDM_free(pedge_extended_to_pedge_sens   [i_part]);
        }

        PDM_free(pextended_face_edge_idx[i_part]);
        PDM_free(pextended_face_edge    [i_part]);
        PDM_free(pextended_edge_vtx_idx [i_part]);
        PDM_free(pextended_edge_vtx     [i_part]);
      }
      else {
        PDM_free(pextended_face_vtx_idx[i_part]);
        PDM_free(pextended_face_vtx    [i_part]);
      }
      
      PDM_free(pface_extended_ln_to_gn          [i_part]);
      PDM_free(pface_extended_to_pface_idx      [i_part]);
      PDM_free(pface_extended_to_pface_triplet  [i_part]);
      PDM_free(pface_extended_to_pface_interface[i_part]);
      if (pface_extended_to_pface_sens!=NULL) {
        PDM_free(pface_extended_to_pface_sens   [i_part]);
      }

      PDM_free(pcell_extended_ln_to_gn          [i_part]);
      PDM_free(pcell_extended_to_pcell_idx      [i_part]);
      PDM_free(pcell_extended_to_pcell_triplet  [i_part]);
      PDM_free(pcell_extended_to_pcell_interface[i_part]);
      PDM_free(pcell_extended_alrdy_sent        [i_part]);
    }

    PDM_free(pvtx_extended_ln_to_gn);
    PDM_free(pvtx_extended_to_pvtx_idx);
    PDM_free(pvtx_extended_to_pvtx_triplet);
    PDM_free(pvtx_extended_to_pvtx_interface);

    PDM_free(pedge_extended_ln_to_gn);
    PDM_free(pedge_extended_to_pedge_idx);
    PDM_free(pedge_extended_to_pedge_triplet);
    PDM_free(pedge_extended_to_pedge_interface);
    PDM_free(pedge_extended_to_pedge_sens);

    PDM_free(pface_extended_ln_to_gn);
    PDM_free(pface_extended_to_pface_idx);
    PDM_free(pface_extended_to_pface_triplet);
    PDM_free(pface_extended_to_pface_interface);
    PDM_free(pface_extended_to_pface_sens);

    PDM_free(pcell_extended_ln_to_gn);
    PDM_free(pcell_extended_to_pcell_idx);
    PDM_free(pcell_extended_to_pcell_triplet);
    PDM_free(pcell_extended_to_pcell_interface);
    PDM_free(pcell_extended_alrdy_sent);

    PDM_free(pextended_cell_face_idx);
    PDM_free(pextended_cell_face);
    PDM_free(pextended_face_vtx_idx);
    PDM_free(pextended_face_vtx);
    PDM_free(pextended_face_edge_idx);
    PDM_free(pextended_face_edge);
    PDM_free(pextended_edge_vtx_idx);
    PDM_free(pextended_edge_vtx);


    /*
     * A chaque étape :
     *   - On garde le même graphe entre les entitiés, mais on agrandit le tableau idx (pour être cohérent avec le part_to_part )
     *   - A l'issu d'une étape, il faut swap le graph avec celui de la nouvelle depth
     *   - Pour avoir l'historique complet on peut agglomerer tout les graphe de chaque depth to have the full one
     */
    if(pn_cell_extended_tot == 0) {
      // Change graph
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {

        _update_propagating_graph_for_depth(i_depth,
                                            part_ext->ln_part_tot,
                                            pn_vtx,
                                            pn_vtx_extended_by_depth,
                                            pfull_vtx_extended_to_pvtx_idx,
                                            pfull_vtx_extended_to_pvtx_triplet,
                                            pfull_vtx_extended_to_pvtx_interface,
                                           &pcurr_entity_bound_to_pentity_bound_idx,
                                           &pcurr_entity_bound_to_pentity_bound_triplet,
                                           &pcurr_entity_bound_to_pentity_bound_interface);
      }
      else if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        _update_propagating_graph_for_depth(i_depth,
                                            part_ext->ln_part_tot,
                                            pn_face,
                                            pn_face_extended_by_depth,
                                            pfull_face_extended_to_pface_idx,
                                            pfull_face_extended_to_pface_triplet,
                                            pfull_face_extended_to_pface_interface,
                                           &pcurr_entity_bound_to_pentity_bound_idx,
                                           &pcurr_entity_bound_to_pentity_bound_triplet,
                                           &pcurr_entity_bound_to_pentity_bound_interface);
      }

      i_depth++;
    }
    step++;


    /*
     * Update for next step
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      
      for (int i_cell=0; i_cell<pn_cell[i_part]; ++i_cell) {
        pcell_alrdy_sent[i_part][i_cell] = 0;
      }

      /* Update size */
      pn_vtx [i_part] += pn_vtx_extended [i_part];
      if (part_ext->have_edge==1) {
        pn_edge[i_part] += pn_edge_extended[i_part];
      }
      pn_face[i_part] += pn_face_extended[i_part];
      pn_cell[i_part] += pn_cell_extended[i_part];
    }
    PDM_free(pn_vtx_extended);
    if (part_ext->have_edge==1) {
      PDM_free(pn_edge_extended);
    }
    PDM_free(pn_face_extended);
    PDM_free(pn_cell_extended);

    // if(step > 6) {
    //   abort();
    // }

    double t_end = PDM_MPI_Wtime();

    if (i_rank==0) printf(" (%.3fs)\n", t_end - t_start);
  }


  /**
   * Copy extended entities on part_extension structure
   */
  // > Cells
  _store_extended_entity(part_ext, 
                         pn_init_cell,
                         pn_cell,
                         pcell_ln_to_gn,
                         pcell_face_idx,
                         pcell_face,
                        &part_ext->n_cell_border,
                        &part_ext->border_cell_ln_to_gn,
                        &part_ext->border_cell_face_idx,
                        &part_ext->border_cell_face,
                         PDM_MESH_ENTITY_CELL, 
                         PDM_CONNECTIVITY_TYPE_CELL_FACE);

  // > Faces
  if(part_ext->have_edge == 1) {
    _store_extended_entity(part_ext, 
                           pn_init_face,
                           pn_face,
                           pface_ln_to_gn,
                           pface_edge_idx,
                           pface_edge,
                          &part_ext->n_face_border,
                          &part_ext->border_face_ln_to_gn,
                          &part_ext->border_face_edge_idx,
                          &part_ext->border_face_edge,
                           PDM_MESH_ENTITY_FACE, 
                           PDM_CONNECTIVITY_TYPE_FACE_EDGE);

    _store_extended_entity(part_ext, 
                           pn_init_edge,
                           pn_edge,
                           pedge_ln_to_gn,
                           pedge_vtx_idx,
                           pedge_vtx,
                          &part_ext->n_edge_border,
                          &part_ext->border_edge_ln_to_gn,
                          &part_ext->border_edge_vtx_idx,
                          &part_ext->border_edge_vtx,
                           PDM_MESH_ENTITY_EDGE, 
                           PDM_CONNECTIVITY_TYPE_EDGE_VTX);
  }
  else {
    _store_extended_entity(part_ext, 
                           pn_init_face,
                           pn_face,
                           pface_ln_to_gn,
                           pface_vtx_idx,
                           pface_vtx,
                          &part_ext->n_face_border,
                          &part_ext->border_face_ln_to_gn,
                          &part_ext->border_face_vtx_idx,
                          &part_ext->border_face_vtx,
                           PDM_MESH_ENTITY_FACE, 
                           PDM_CONNECTIVITY_TYPE_FACE_VTX);
  }

  // _store_extended_entity(part_ext, 
  //                        pn_init_vtx,
  //                        pn_vtx,
  //                        pvtx_ln_to_gn,
  //                        NULL,
  //                        NULL,
  //                       &part_ext->n_vtx_border,
  //                       &part_ext->border_vtx_ln_to_gn,
  //                        NULL,
  //                        NULL);
  _store_extended_vtx(part_ext, 
                         pn_init_vtx,
                         pn_vtx,
                         pvtx_ln_to_gn,
                         pvtx_coords,
                        &part_ext->n_vtx_border,
                        &part_ext->border_vtx_ln_to_gn,
                        &part_ext->border_vtx);


  /**
   * Need to get distributed data_base info on partition
   */

  // > Cells
  int         **pcell_ancstr_strd    = NULL;
  PDM_g_num_t **pcell_ancstr         = NULL;
  // int         **pcell_path_itrf_strd = NULL;
  // int         **pcell_path_itrf      = NULL;
  _get_block_data_base_on_part_ext( part_ext,
                                    prev_dcell_itrf_n_blk,
                                    prev_dcell_itrf_blk_gnum,
                                    prev_dcell_itrf_blk_ancstr_strd,
                                    prev_dcell_itrf_blk_ancstr,
                                    prev_dcell_itrf_blk_path_itrf_strd,
                                    prev_dcell_itrf_blk_path_itrf,
                                    part_ext->n_cell_border,
                                    part_ext->border_cell_ln_to_gn,
                                   &pcell_ancstr_strd,
                                   &pcell_ancstr,
                                   // &pcell_path_itrf_strd,
                                   &part_ext->cell_cell_path_itrf_idx,
                                   &part_ext->cell_cell_path_itrf,
                                    PDM_MESH_ENTITY_CELL);


  // > Faces
  int         **pface_ancstr_strd    = NULL;
  PDM_g_num_t **pface_ancstr         = NULL;
  // int         **pface_path_itrf_strd = NULL;
  // int         **pface_path_itrf      = NULL;
  _get_block_data_base_on_part_ext( part_ext,
                                    prev_dface_itrf_n_blk,
                                    prev_dface_itrf_blk_gnum,
                                    prev_dface_itrf_blk_ancstr_strd,
                                    prev_dface_itrf_blk_ancstr,
                                    prev_dface_itrf_blk_path_itrf_strd,
                                    prev_dface_itrf_blk_path_itrf,
                                    part_ext->n_face_border,
                                    part_ext->border_face_ln_to_gn,
                                   &pface_ancstr_strd,
                                   &pface_ancstr,
                                   // &pface_path_itrf_strd,
                                   &part_ext->face_face_path_itrf_idx,
                                   &part_ext->face_face_path_itrf,
                                    PDM_MESH_ENTITY_FACE);


  // > Edges
  int         **pedge_ancstr_strd    = NULL;
  PDM_g_num_t **pedge_ancstr         = NULL;
  // int         **pedge_path_itrf_strd = NULL;
  // int         **pedge_path_itrf      = NULL;
  if(part_ext->have_edge == 1) {
    _get_block_data_base_on_part_ext( part_ext,
                                      prev_dedge_itrf_n_blk,
                                      prev_dedge_itrf_blk_gnum,
                                      prev_dedge_itrf_blk_ancstr_strd,
                                      prev_dedge_itrf_blk_ancstr,
                                      prev_dedge_itrf_blk_path_itrf_strd,
                                      prev_dedge_itrf_blk_path_itrf,
                                      part_ext->n_edge_border,
                                      part_ext->border_edge_ln_to_gn,
                                     &pedge_ancstr_strd,
                                     &pedge_ancstr,
                                     // &pedge_path_itrf_strd,
                                     &part_ext->edge_edge_path_itrf_idx,
                                     &part_ext->edge_edge_path_itrf,
                                      PDM_MESH_ENTITY_EDGE);
  }

  // > Vertex
  int         **pvtx_ancstr_strd    = NULL;
  PDM_g_num_t **pvtx_ancstr         = NULL;
  // int         **pvtx_path_itrf_strd = NULL;
  // int         **pvtx_path_itrf      = NULL;
  _get_block_data_base_on_part_ext( part_ext,
                                    prev_dvtx_itrf_n_blk,
                                    prev_dvtx_itrf_blk_gnum,
                                    prev_dvtx_itrf_blk_ancstr_strd,
                                    prev_dvtx_itrf_blk_ancstr,
                                    prev_dvtx_itrf_blk_path_itrf_strd,
                                    prev_dvtx_itrf_blk_path_itrf,
                                    part_ext->n_vtx_border,
                                    part_ext->border_vtx_ln_to_gn,
                                   &pvtx_ancstr_strd,
                                   &pvtx_ancstr,
                                   // &pvtx_path_itrf_strd,
                                   &part_ext->vtx_vtx_path_itrf_idx,
                                   &part_ext->vtx_vtx_path_itrf,
                                    PDM_MESH_ENTITY_VTX);


  /**
   * Use ancestor info to build graph between original partition and extended partition
   */
  // > Cells
  _build_part_extension_graph_to_old(part_ext,
                                     pn_init_cell,
                                     pcell_ln_to_gn,
                                     NULL,
                                     NULL,
                                     NULL,
                                     part_ext->n_cell_border,
                                     part_ext->border_cell_ln_to_gn,
                                     pcell_ancstr_strd,
                                     pcell_ancstr,
                                     // pcell_path_itrf_strd,
                                     // pcell_path_itrf,
                                     part_ext->shift_by_domain_cell,
                                    &part_ext->border_cell_ln_to_gn_ancstr,
                                    &part_ext->cell_cell_extended2,
                                     NULL,
                                     NULL,
                                     NULL,
                                     PDM_MESH_ENTITY_CELL);

  // > Faces
  _build_part_extension_graph_to_old(part_ext,
                                     pn_init_face,
                                     pface_ln_to_gn,
                                     pn_face_group,
                                     pface_group_tag,
                                     pface_group_gnum,
                                     part_ext->n_face_border,
                                     part_ext->border_face_ln_to_gn,
                                     pface_ancstr_strd,
                                     pface_ancstr,
                                     // pface_path_itrf_strd,
                                     // pface_path_itrf,
                                     part_ext->shift_by_domain_face,
                                    &part_ext->border_face_ln_to_gn_ancstr,
                                    &part_ext->face_face_extended,
                                    &part_ext->border_face_group_idx,
                                    &part_ext->border_face_group,
                                    &part_ext->border_face_group_ln_to_gn,
                                     PDM_MESH_ENTITY_FACE);

  // > Edges
  if(part_ext->have_edge == 1) {
    _build_part_extension_graph_to_old(part_ext,
                                       pn_init_edge,
                                       pedge_ln_to_gn,
                                       NULL,
                                       NULL,
                                       NULL,
                                       part_ext->n_edge_border,
                                       part_ext->border_edge_ln_to_gn,
                                       pedge_ancstr_strd,
                                       pedge_ancstr,
                                       // pedge_path_itrf_strd,
                                       // pedge_path_itrf,
                                       part_ext->shift_by_domain_edge,
                                      &part_ext->border_edge_ln_to_gn_ancstr,
                                      &part_ext->edge_edge_extended,
                                       NULL,
                                       NULL,
                                       NULL,
                                       PDM_MESH_ENTITY_EDGE);
  }

  // > Vertex
  _build_part_extension_graph_to_old(part_ext,
                                     pn_init_vtx,
                                     pvtx_ln_to_gn,
                                     NULL,
                                     NULL,
                                     NULL,
                                     part_ext->n_vtx_border,
                                     part_ext->border_vtx_ln_to_gn,
                                     pvtx_ancstr_strd,
                                     pvtx_ancstr,
                                     // pvtx_path_itrf_strd,
                                     // pvtx_path_itrf,
                                     part_ext->shift_by_domain_vtx,
                                    &part_ext->border_vtx_ln_to_gn_ancstr,
                                    &part_ext->vtx_vtx_extended,
                                     NULL,
                                     NULL,
                                     NULL,
                                     PDM_MESH_ENTITY_VTX);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pcell_ancstr_strd   [i_part]);
    PDM_free(pcell_ancstr        [i_part]);
    // PDM_free(pcell_path_itrf_strd[i_part]);
    // PDM_free(pcell_path_itrf     [i_part]);

    PDM_free(pface_ancstr_strd   [i_part]);
    PDM_free(pface_ancstr        [i_part]);
    // PDM_free(pface_path_itrf_strd[i_part]);
    // PDM_free(pface_path_itrf     [i_part]);

    if(part_ext->have_edge == 1) {
      PDM_free(pedge_ancstr_strd   [i_part]);
      PDM_free(pedge_ancstr        [i_part]);
      // PDM_free(pedge_path_itrf_strd[i_part]);
      // PDM_free(pedge_path_itrf     [i_part]);
    }

    PDM_free(pvtx_ancstr_strd    [i_part]);
    PDM_free(pvtx_ancstr         [i_part]);
    // PDM_free(pvtx_path_itrf_strd [i_part]);
    // PDM_free(pvtx_path_itrf      [i_part]);
  }
  PDM_free(pcell_ancstr_strd);
  PDM_free(pcell_ancstr);
  // PDM_free(pcell_path_itrf_strd);
  // PDM_free(pcell_path_itrf);

  PDM_free(pface_ancstr_strd);
  PDM_free(pface_ancstr);
  // PDM_free(pface_path_itrf_strd);
  // PDM_free(pface_path_itrf);
  
  if(part_ext->have_edge == 1) {
    PDM_free(pedge_ancstr_strd);
    PDM_free(pedge_ancstr);
    // PDM_free(pedge_path_itrf_strd);
    // PDM_free(pedge_path_itrf);
  }
  
  PDM_free(pvtx_ancstr_strd);
  PDM_free(pvtx_ancstr);
  // PDM_free(pvtx_path_itrf_strd);
  // PDM_free(pvtx_path_itrf);



  PDM_free(prev_dcell_itrf_blk_gnum           );
  PDM_free(prev_dcell_itrf_blk_ancstr_strd    );
  PDM_free(prev_dcell_itrf_blk_ancstr         );
  PDM_free(prev_dcell_itrf_blk_path_itrf_strd );
  PDM_free(prev_dcell_itrf_blk_path_itrf      );
  PDM_free(prev_dcell_itrf_gnum_and_itrf_strid);
  PDM_free(prev_dcell_itrf_gnum_and_itrf_data );

  PDM_free(prev_dface_itrf_blk_ancstr_strd    );
  PDM_free(prev_dface_itrf_blk_ancstr         );
  PDM_free(prev_dface_itrf_blk_path_itrf_strd );
  PDM_free(prev_dface_itrf_blk_path_itrf      );

  if (part_ext->have_edge==1) {
    PDM_free(prev_dedge_itrf_blk_ancstr_strd    );
    PDM_free(prev_dedge_itrf_blk_ancstr         );
    PDM_free(prev_dedge_itrf_blk_path_itrf_strd );
    PDM_free(prev_dedge_itrf_blk_path_itrf      );
  }
  
  PDM_free(prev_dvtx_itrf_blk_ancstr_strd);
  PDM_free(prev_dvtx_itrf_blk_ancstr);
  PDM_free(prev_dvtx_itrf_blk_path_itrf_strd);
  PDM_free(prev_dvtx_itrf_blk_path_itrf);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pvtx_ln_to_gn      [i_part]);
    if (part_ext->have_edge==1) {
      PDM_free(pedge_ln_to_gn     [i_part]);
    }
    PDM_free(pface_ln_to_gn     [i_part]);
    PDM_free(pcell_ln_to_gn     [i_part]);
    PDM_free(pcell_alrdy_sent   [i_part]);
    PDM_free(pcell_face_idx     [i_part]);
    PDM_free(pcell_face         [i_part]);
    if (part_ext->have_edge==1) {
      PDM_free(pface_edge_idx   [i_part]);
      PDM_free(pface_edge       [i_part]);
      PDM_free(pedge_vtx_idx    [i_part]);
      PDM_free(pedge_vtx        [i_part]);
    }
    else {
      PDM_free(pface_vtx_idx    [i_part]);
      PDM_free(pface_vtx        [i_part]);
    }
    PDM_free(pface_group_tag    [i_part]);
    PDM_free(pface_group_gnum   [i_part]);
    PDM_free(pvtx_coords        [i_part]);
  }

  PDM_free(pn_vtx);
  PDM_free(pn_init_vtx);
  PDM_free(pvtx_ln_to_gn);
  PDM_free(pn_edge);
  PDM_free(pn_init_edge);
  PDM_free(pedge_ln_to_gn);
  PDM_free(pn_face);
  PDM_free(pn_init_face);
  PDM_free(pface_ln_to_gn);
  PDM_free(pn_cell);
  PDM_free(pn_init_cell);
  PDM_free(pcell_ln_to_gn);
  PDM_free(pcell_alrdy_sent);
  PDM_free(pcell_face_idx);
  PDM_free(pcell_face);
  PDM_free(pcell_vtx_idx);
  PDM_free(pcell_vtx);
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge);
  PDM_free(pface_vtx_idx);
  PDM_free(pface_vtx);
  PDM_free(pedge_vtx_idx);
  PDM_free(pedge_vtx);
  PDM_free(pn_face_group);
  PDM_free(pface_group_tag);
  PDM_free(pface_group_gnum);
  PDM_free(pvtx_coords);


  for(int i = 0; i < part_ext->depth; ++i) {
    PDM_free(pn_vtx_extended_by_depth [i]);
    if (part_ext->have_edge==1) {
      PDM_free(pn_edge_extended_by_depth[i]);
    }
    PDM_free(pn_face_extended_by_depth[i]);
    PDM_free(pn_cell_extended_by_depth[i]);
  }
  PDM_free(pn_vtx_extended_by_depth );
  PDM_free(pn_edge_extended_by_depth);
  PDM_free(pn_face_extended_by_depth);
  PDM_free(pn_cell_extended_by_depth);

  /*
   * To keep
   */
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pfull_vtx_extended_ln_to_gn         [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_idx      [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_triplet  [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_interface[i_part]);

    if (part_ext->have_edge==1) {
      PDM_free(pfull_edge_extended_ln_to_gn          [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_idx      [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_triplet  [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_interface[i_part]);
    }

    PDM_free(pfull_face_extended_ln_to_gn          [i_part]);
    PDM_free(pfull_face_extended_to_pface_idx      [i_part]);
    PDM_free(pfull_face_extended_to_pface_triplet  [i_part]);
    PDM_free(pfull_face_extended_to_pface_interface[i_part]);

    PDM_free(pfull_cell_extended_ln_to_gn          [i_part]);
    PDM_free(pfull_cell_extended_to_pcell_idx      [i_part]);
    PDM_free(pfull_cell_extended_to_pcell_triplet  [i_part]);
    PDM_free(pfull_cell_extended_to_pcell_interface[i_part]);

  }
  PDM_free(pn_vtx_extended_old                 );
  PDM_free(pfull_n_vtx_extended                );
  PDM_free(pfull_vtx_extended_ln_to_gn         );
  PDM_free(pfull_vtx_extended_to_pvtx_idx      );
  PDM_free(pfull_vtx_extended_to_pvtx_triplet  );
  PDM_free(pfull_vtx_extended_to_pvtx_interface);

  PDM_free(pn_edge_extended_old                  );
  PDM_free(pfull_n_edge_extended                 );
  PDM_free(pfull_edge_extended_ln_to_gn          );
  PDM_free(pfull_edge_extended_to_pedge_idx      );
  PDM_free(pfull_edge_extended_to_pedge_triplet  );
  PDM_free(pfull_edge_extended_to_pedge_interface);

  PDM_free(pn_face_extended_old);
  PDM_free(pfull_n_face_extended                 );
  PDM_free(pfull_face_extended_ln_to_gn          );
  PDM_free(pfull_face_extended_to_pface_idx      );
  PDM_free(pfull_face_extended_to_pface_triplet  );
  PDM_free(pfull_face_extended_to_pface_interface);

  PDM_free(pn_cell_extended_old);
  PDM_free(pfull_n_cell_extended                 );
  PDM_free(pfull_cell_extended_ln_to_gn          );
  PDM_free(pfull_cell_extended_to_pcell_idx      );
  PDM_free(pfull_cell_extended_to_pcell_triplet  );
  PDM_free(pfull_cell_extended_to_pcell_interface);

}




static
void
_part_extension_2d
(
 PDM_part_extension_t *part_ext
)
{
  int visu  = 0;
  int debug = 0;

  /*
   * 2 possibilities :
   *   - With face_vtx
   *   - With face_edge + edge_vtx
   */

  /* Size */
  int *pn_vtx       = NULL;
  int *pn_edge      = NULL;
  int *pn_face      = NULL;
  int *pn_init_vtx  = NULL;
  int *pn_init_edge = NULL;
  int *pn_init_face = NULL;
  PDM_malloc(pn_vtx      , part_ext->ln_part_tot, int );
  PDM_malloc(pn_edge     , part_ext->ln_part_tot, int );
  PDM_malloc(pn_face     , part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_vtx , part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_edge, part_ext->ln_part_tot, int );
  PDM_malloc(pn_init_face, part_ext->ln_part_tot, int );

  /* ln_to_gn */
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_malloc(pvtx_ln_to_gn , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pedge_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pface_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);

  /* Connectivity */
  int **pface_edge_idx = NULL;
  int **pface_edge     = NULL;
  int **pface_vtx_idx  = NULL;
  int **pface_vtx      = NULL;
  int **pedge_vtx_idx  = NULL;
  int **pedge_vtx      = NULL;
  PDM_malloc(pface_edge_idx, part_ext->ln_part_tot, int *);
  PDM_malloc(pface_edge    , part_ext->ln_part_tot, int *);
  PDM_malloc(pface_vtx_idx , part_ext->ln_part_tot, int *);
  PDM_malloc(pface_vtx     , part_ext->ln_part_tot, int *);
  PDM_malloc(pedge_vtx_idx , part_ext->ln_part_tot, int *);
  PDM_malloc(pedge_vtx     , part_ext->ln_part_tot, int *);

  /* Groups */
  int          *pn_edge_group    = NULL;
  int         **pedge_group_tag  = NULL;
  PDM_g_num_t **pedge_group_gnum = NULL;
  PDM_malloc(pn_edge_group   , part_ext->ln_part_tot, int          );
  PDM_malloc(pedge_group_tag , part_ext->ln_part_tot, int         *);
  PDM_malloc(pedge_group_gnum, part_ext->ln_part_tot, PDM_g_num_t *);


  /* Coordinates */
  double **pvtx_coords = NULL;
  PDM_malloc(pvtx_coords, part_ext->ln_part_tot, double *);

  int **pface_alrdy_sent = NULL;
  PDM_malloc(pface_alrdy_sent, part_ext->ln_part_tot, int *);

  int lpart = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      // Init at null to avoid wrong free
      pface_edge_idx[lpart] = NULL;
      pface_edge    [lpart] = NULL;
      pface_vtx_idx [lpart] = NULL;
      pface_vtx     [lpart] = NULL;
      pedge_vtx_idx [lpart] = NULL;
      pedge_vtx     [lpart] = NULL;

      pn_vtx        [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_init_vtx   [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      if(part_ext->have_edge == 1) {
        pn_edge     [lpart] = part_ext->parts[i_domain][i_part].n_edge;
        pn_init_edge[lpart] = part_ext->parts[i_domain][i_part].n_edge;
      }
      pn_face       [lpart] = part_ext->parts[i_domain][i_part].n_face;
      pn_init_face  [lpart] = part_ext->parts[i_domain][i_part].n_face;

      /* Copy to realloc after all step */
      PDM_malloc(pvtx_ln_to_gn [lpart], pn_vtx [lpart], PDM_g_num_t);
      PDM_malloc(pface_ln_to_gn[lpart], pn_face[lpart], PDM_g_num_t);
      if(part_ext->have_edge == 1) {
        PDM_malloc(pface_edge_idx[lpart],                                                 pn_face[lpart]+1, int);
        PDM_malloc(pface_edge    [lpart], part_ext->parts[i_domain][i_part].face_edge_idx[pn_face[lpart]] , int);
        PDM_malloc(pedge_ln_to_gn[lpart],   pn_edge[lpart]   , PDM_g_num_t);
        PDM_malloc(pedge_vtx_idx [lpart],   pn_edge[lpart]+1 , int        );
        PDM_malloc(pedge_vtx     [lpart], 2*pn_edge[lpart]   , int        );
        pedge_vtx_idx    [lpart][0] = 0;

        PDM_calloc(pedge_group_tag [lpart],  pn_edge[lpart]   , int        );
        PDM_calloc(pedge_group_gnum[lpart],  pn_edge[lpart]   , PDM_g_num_t);
      }
      else {
        PDM_malloc(pface_vtx_idx[lpart],                                                pn_face[lpart]+1, int);
        PDM_malloc(pface_vtx    [lpart], part_ext->parts[i_domain][i_part].face_vtx_idx[pn_face[lpart]] , int);
      }
      PDM_malloc(pvtx_coords     [lpart], 3 * pn_vtx [lpart], double);
      PDM_malloc(pface_alrdy_sent[lpart],     pn_face[lpart], int   );


      for(int i_face = 0; i_face < pn_face[lpart]; ++i_face) {
        pface_ln_to_gn  [lpart][i_face] = part_ext->parts[i_domain][i_part].face_ln_to_gn[i_face];
        pface_alrdy_sent[lpart][i_face] = 0;
      }
      if(part_ext->have_edge == 1) {
        for(int i_edge = 0; i_edge < pn_edge[lpart]; ++i_edge) {
          pedge_ln_to_gn[lpart][i_edge] = part_ext->parts[i_domain][i_part].edge_ln_to_gn[i_edge];
        }
      }
      for(int i_vtx = 0; i_vtx < pn_vtx[lpart]; ++i_vtx) {
        pvtx_ln_to_gn[lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn[i_vtx];
      }

      if(part_ext->have_edge == 1) {
        for(int i_face = 0; i_face < pn_face[lpart]+1; ++i_face) {
          pface_edge_idx [lpart][i_face] = part_ext->parts[i_domain][i_part].face_edge_idx[i_face];
        }
        for(int idx = 0; idx < pface_edge_idx[i_part][pn_face[lpart]]; ++idx) {
          pface_edge [lpart][idx] = part_ext->parts[i_domain][i_part].face_edge[idx];
        }
        for(int i_edge = 0; i_edge < pn_edge[lpart]; ++i_edge) {
          pedge_vtx_idx[lpart][  i_edge+1] = (i_edge+1)*2;
          pedge_vtx    [lpart][2*i_edge  ] = part_ext->parts[i_domain][i_part].edge_vtx[2*i_edge  ];
          pedge_vtx    [lpart][2*i_edge+1] = part_ext->parts[i_domain][i_part].edge_vtx[2*i_edge+1];
        }
      }
      else {
        for(int i_face = 0; i_face < pn_face[lpart]+1; ++i_face) {
          pface_vtx_idx [lpart][i_face] = part_ext->parts[i_domain][i_part].face_vtx_idx[i_face];
        }
        for(int idx = 0; idx < pface_vtx_idx[i_part][pn_face[lpart]]; ++idx) {
          pface_vtx [lpart][idx] = part_ext->parts[i_domain][i_part].face_vtx[idx];
        }
      }

      for(int i_vtx = 0; i_vtx < 3 * pn_vtx[lpart]; ++i_vtx) {
        pvtx_coords   [lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx[i_vtx];
      }

      // > Transform group into tag
      if (part_ext->have_edge) {
        

        int          _n_group           = part_ext->parts[i_domain][i_part].n_edge_group;
        int         *_group_entity_idx  = part_ext->parts[i_domain][i_part].edge_bound_idx;
        int         *_group_entity      = part_ext->parts[i_domain][i_part].edge_bound;
        PDM_g_num_t *_group_entity_gnum = part_ext->parts[i_domain][i_part].edge_bound_ln_to_gn;
        pn_edge_group[lpart] = _n_group;
        for (int i_group=0; i_group<_n_group; ++i_group) {
          for (int i_entity=_group_entity_idx[i_group  ];
                   i_entity<_group_entity_idx[i_group+1];
                   ++i_entity) {
            int lnum = _group_entity     [i_entity];
            int gnum = _group_entity_gnum[i_entity];
            pedge_group_tag [lpart][lnum-1] = i_group+1;
            pedge_group_gnum[lpart][lnum-1] = gnum;
          }
        }
      }

      lpart++;
    }
  }


  /*
   * On va etendre la partition avec le graphe de base tant que de nouveau elements apparaissent
   *   -> Il faut ajuster la taille du graphe en fonction des nouvelles entités (juste une rallonge)
   *   -> On doit également alimenter un tableau pour le lien avec les entités de la recursion d'après
   */
  int           *pn_vtx_extended_old                  = NULL;
  int           *pfull_n_vtx_extended                 = NULL;
  PDM_g_num_t  **pfull_vtx_extended_ln_to_gn          = NULL;
  int          **pfull_vtx_extended_to_pvtx_idx       = NULL;
  int          **pfull_vtx_extended_to_pvtx_triplet   = NULL;
  int          **pfull_vtx_extended_to_pvtx_interface = NULL;
  PDM_malloc(pn_vtx_extended_old                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_vtx_extended                , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_vtx_extended_ln_to_gn         , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_interface, part_ext->ln_part_tot, int         *);

  int           *pn_edge_extended_old                   = NULL;
  int           *pfull_n_edge_extended                  = NULL;
  PDM_g_num_t  **pfull_edge_extended_ln_to_gn           = NULL;
  int          **pfull_edge_extended_to_pedge_idx       = NULL;
  int          **pfull_edge_extended_to_pedge_triplet   = NULL;
  int          **pfull_edge_extended_to_pedge_interface = NULL;
  PDM_malloc(pn_edge_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_edge_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_edge_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_edge_extended_to_pedge_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_interface, part_ext->ln_part_tot, int         *);

  int           *pn_face_extended_old                   = NULL;
  int           *pfull_n_face_extended                  = NULL;
  PDM_g_num_t  **pfull_face_extended_ln_to_gn           = NULL;
  int          **pfull_face_extended_to_pface_idx       = NULL;
  int          **pfull_face_extended_to_pface_triplet   = NULL;
  int          **pfull_face_extended_to_pface_interface = NULL;
  PDM_malloc(pn_face_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_face_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_face_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_face_extended_to_pface_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_face_extended_to_pface_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_face_extended_to_pface_interface, part_ext->ln_part_tot, int         *);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

    pfull_n_vtx_extended[i_part] = 0;
    PDM_malloc(pfull_vtx_extended_to_pvtx_idx      [i_part], pfull_n_vtx_extended[i_part]+1, int);
    PDM_malloc(pfull_vtx_extended_ln_to_gn         [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_vtx_extended_to_pvtx_interface[i_part], 0, int        );
    pfull_vtx_extended_to_pvtx_idx[i_part][0] = 0;

    if (part_ext->have_edge==1) {
      pfull_n_edge_extended[i_part] = 0;
      PDM_malloc(pfull_edge_extended_to_pedge_idx      [i_part], pfull_n_edge_extended[i_part]+1, int);
      PDM_malloc(pfull_edge_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
      PDM_malloc(pfull_edge_extended_to_pedge_triplet  [i_part], 0, int        );
      PDM_malloc(pfull_edge_extended_to_pedge_interface[i_part], 0, int        );
      pfull_edge_extended_to_pedge_idx[i_part][0] = 0;
    }

    pfull_n_face_extended[i_part] = 0;
    PDM_malloc(pfull_face_extended_to_pface_idx      [i_part], pfull_n_face_extended[i_part]+1, int);
    PDM_malloc(pfull_face_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_face_extended_to_pface_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_face_extended_to_pface_interface[i_part], 0, int        );
    pfull_face_extended_to_pface_idx[i_part][0] = 0;
  }

  PDM_g_num_t shift_by_domain_vtx  = part_ext->shift_by_domain_vtx [part_ext->n_domain];
  PDM_g_num_t shift_by_domain_edge = part_ext->shift_by_domain_edge[part_ext->n_domain];
  PDM_g_num_t shift_by_domain_face = part_ext->shift_by_domain_face[part_ext->n_domain];


  /**
   * Face link between interface
   */
  int          prev_dface_itrf_n_blk               = 0;
  PDM_g_num_t *prev_dface_itrf_blk_gnum            = NULL;
  int         *prev_dface_itrf_blk_ancstr_strd     = NULL;
  PDM_g_num_t *prev_dface_itrf_blk_ancstr          = NULL;
  int         *prev_dface_itrf_blk_path_itrf_strd  = NULL;
  int         *prev_dface_itrf_blk_path_itrf       = NULL;
  int         *prev_dface_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *prev_dface_itrf_gnum_and_itrf_data  = NULL;


  /**
   * Edge link between interface
   */
  int          prev_dedge_itrf_n_blk               = 0;
  PDM_g_num_t *prev_dedge_itrf_blk_gnum            = NULL;
  int         *prev_dedge_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *prev_dedge_itrf_gnum_and_itrf_data  = NULL;
  int         *prev_dedge_itrf_gnum_and_itrf_sens  = NULL;
  int         *prev_dedge_itrf_blk_ancstr_strd     = NULL;
  PDM_g_num_t *prev_dedge_itrf_blk_ancstr          = NULL;
  int         *prev_dedge_itrf_blk_path_itrf_strd  = NULL;
  int         *prev_dedge_itrf_blk_path_itrf       = NULL;
  if(part_ext->have_edge == 1) {
    prev_dedge_itrf_n_blk               = part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_blk_gnum            = part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_strid = part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_data  = part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_gnum_and_itrf_sens  = part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE];
    prev_dedge_itrf_blk_ancstr_strd     = PDM_array_zeros_int(prev_dedge_itrf_n_blk);
    prev_dedge_itrf_blk_ancstr          = PDM_array_zeros_gnum(0);
    prev_dedge_itrf_blk_path_itrf_strd  = PDM_array_zeros_int(prev_dedge_itrf_n_blk);
    prev_dedge_itrf_blk_path_itrf       = PDM_array_zeros_int(0);
  }


  /**
   * Vtx link between interface
   */
  int          prev_dvtx_itrf_n_blk               = part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX];
  PDM_g_num_t *prev_dvtx_itrf_blk_gnum            = part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX];
  int         *prev_dvtx_itrf_gnum_and_itrf_strid = part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX];
  PDM_g_num_t *prev_dvtx_itrf_gnum_and_itrf_data  = part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX];
  int         *prev_dvtx_itrf_blk_ancstr_strd     = PDM_array_zeros_int(prev_dvtx_itrf_n_blk);
  PDM_g_num_t *prev_dvtx_itrf_blk_ancstr          = PDM_array_zeros_gnum(0);
  int         *prev_dvtx_itrf_blk_path_itrf_strd  = PDM_array_zeros_int(prev_dvtx_itrf_n_blk);
  int         *prev_dvtx_itrf_blk_path_itrf       = PDM_array_zeros_int(0);

  int **pcurr_entity_bound_to_pentity_bound_idx       = part_ext->pinit_entity_bound_to_pentity_bound_idx;
  int **pcurr_entity_bound_to_pentity_bound_triplet   = part_ext->pinit_entity_bound_to_pentity_bound_triplet;
  int **pcurr_entity_bound_to_pentity_bound_interface = part_ext->pinit_entity_bound_to_pentity_bound_interface;


  int **pn_vtx_extended_by_depth   = NULL;
  int **pn_edge_extended_by_depth  = NULL;
  int **pn_face_extended_by_depth  = NULL;
  PDM_malloc(pn_vtx_extended_by_depth , part_ext->depth, int *);
  PDM_malloc(pn_edge_extended_by_depth, part_ext->depth, int *);
  PDM_malloc(pn_face_extended_by_depth, part_ext->depth, int *);
  for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {
    PDM_malloc(pn_vtx_extended_by_depth   [i_depth], part_ext->ln_part_tot, int);
    if(part_ext->have_edge == 1) {
      PDM_malloc(pn_edge_extended_by_depth[i_depth], part_ext->ln_part_tot, int);
    }
    PDM_malloc(pn_face_extended_by_depth  [i_depth], part_ext->ln_part_tot, int);

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      pn_vtx_extended_by_depth [i_depth][i_part] = 0;
      if(part_ext->have_edge == 1) {
        pn_edge_extended_by_depth[i_depth][i_part] = 0;
      }
      pn_face_extended_by_depth[i_depth][i_part] = 0;
    }
  }


  if (part_ext->extend_type==PDM_EXTEND_FROM_EDGE && part_ext->have_edge==0) {
    PDM_error(__FILE__, __LINE__, 0, "part_extension with extend_type %d asked but edge seems to miss\n", part_ext->extend_type);
  }


  int i_depth   = 0;
  int step      = 0;
  while(i_depth < part_ext->depth) {



    int i_rank;
    PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
    if (i_rank==0) printf("Computing DEPTH %d (step = %i) \n", i_depth, step);

    if (debug==1) {
      log_trace("\n\n\n >> DEPTH %d step = %i\n", i_depth, step);
    }
    double t_start = PDM_MPI_Wtime();

    /* Use descending connectivity to deduce connectivity and extend_face */
    int          *pn_face_extended                  = NULL;
    PDM_g_num_t **pface_extended_ln_to_gn           = NULL;
    int         **pface_extended_alrdy_sent         = NULL;
    int         **pface_extended_to_pface_idx       = NULL;
    int         **pface_extended_to_pface_triplet   = NULL;
    int         **pface_extended_to_pface_interface = NULL;



    if(part_ext->have_edge == 1 && part_ext->extend_type==PDM_EXTEND_FROM_VTX) {
      /**
       * From face_edge and edge_vtx, get face_vtx to compute extended faces.
       * Then, use extended faces to rebuild face_edge_extended and edge_vtx_extended
       */
      for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
        PDM_compute_face_vtx_from_face_and_edge(pn_face[i_part],
                                                pface_edge_idx[i_part],
                                                pface_edge[i_part],
                                                pedge_vtx[i_part],
                                               &pface_vtx[i_part]);
        pface_vtx_idx[i_part] = pface_edge_idx[i_part];
      }
    }

    /*
     * Hook the most leading entity connexion
     */
    int          next_dface_itrf_n_blk               = 0;
    PDM_g_num_t *next_dface_itrf_blk_gnum            = NULL;
    int         *next_dface_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dface_itrf_blk_ancstr          = NULL;
    int         *next_dface_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dface_itrf_blk_path_itrf       = NULL;
    int         *next_dface_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dface_itrf_gnum_and_itrf_data  = NULL;
    
    if (debug==1) {
      log_trace("\n\n");
      log_trace("========================================= \n");
      log_trace("PDM_part_extension_entity1_to_entity2 beg \n");
    }

    if(part_ext->extend_type==PDM_EXTEND_FROM_VTX) {
      PDM_part_extension_entity1_to_entity2(shift_by_domain_face, // Attention il va evoluer lui
                                            part_ext->ln_part_tot,
                                            pn_vtx,
                                            pvtx_ln_to_gn,
                                            pcurr_entity_bound_to_pentity_bound_idx,
                                            pcurr_entity_bound_to_pentity_bound_triplet,
                                            pcurr_entity_bound_to_pentity_bound_interface,
                                            pn_face,
                                            pface_ln_to_gn,
                                            pface_alrdy_sent,
                                            pface_vtx_idx,
                                            pface_vtx,
                                            prev_dface_itrf_n_blk,
                                            prev_dface_itrf_blk_gnum,
                                            prev_dface_itrf_blk_ancstr_strd,
                                            prev_dface_itrf_blk_ancstr,
                                            prev_dface_itrf_blk_path_itrf_strd,
                                            prev_dface_itrf_blk_path_itrf,
                                            prev_dface_itrf_gnum_and_itrf_strid,
                                            prev_dface_itrf_gnum_and_itrf_data,
                                            &pn_face_extended,
                                            &pface_extended_ln_to_gn,
                                            &pface_extended_alrdy_sent,
                                            &pface_extended_to_pface_idx,
                                            &pface_extended_to_pface_triplet,
                                            &pface_extended_to_pface_interface,
                                            &next_dface_itrf_n_blk,
                                            &next_dface_itrf_blk_gnum,
                                            &next_dface_itrf_blk_ancstr_strd,
                                            &next_dface_itrf_blk_ancstr,
                                            &next_dface_itrf_blk_path_itrf_strd,
                                            &next_dface_itrf_blk_path_itrf,
                                            &next_dface_itrf_gnum_and_itrf_strid,
                                            &next_dface_itrf_gnum_and_itrf_data,
                                            part_ext->comm);
    }
    else if (part_ext->extend_type==PDM_EXTEND_FROM_EDGE) {
      PDM_part_extension_entity1_to_entity2(shift_by_domain_face, // Attention il va evoluer lui
                                            part_ext->ln_part_tot,
                                            pn_edge,
                                            pedge_ln_to_gn,
                                            pcurr_entity_bound_to_pentity_bound_idx,
                                            pcurr_entity_bound_to_pentity_bound_triplet,
                                            pcurr_entity_bound_to_pentity_bound_interface,
                                            pn_face,
                                            pface_ln_to_gn,
                                            pface_alrdy_sent,
                                            pface_edge_idx,
                                            pface_edge,
                                            prev_dface_itrf_n_blk,
                                            prev_dface_itrf_blk_gnum,
                                            prev_dface_itrf_blk_ancstr_strd,
                                            prev_dface_itrf_blk_ancstr,
                                            prev_dface_itrf_blk_path_itrf_strd,
                                            prev_dface_itrf_blk_path_itrf,
                                            prev_dface_itrf_gnum_and_itrf_strid,
                                            prev_dface_itrf_gnum_and_itrf_data,
                                            &pn_face_extended,
                                            &pface_extended_ln_to_gn,
                                            &pface_extended_alrdy_sent,
                                            &pface_extended_to_pface_idx,
                                            &pface_extended_to_pface_triplet,
                                            &pface_extended_to_pface_interface,
                                            &next_dface_itrf_n_blk,
                                            &next_dface_itrf_blk_gnum,
                                            &next_dface_itrf_blk_ancstr_strd,
                                            &next_dface_itrf_blk_ancstr,
                                            &next_dface_itrf_blk_path_itrf_strd,
                                            &next_dface_itrf_blk_path_itrf,
                                            &next_dface_itrf_gnum_and_itrf_strid,
                                            &next_dface_itrf_gnum_and_itrf_data,
                                            part_ext->comm);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "part_extension with extend_type %d is invalid in 2d\n", part_ext->extend_type);
    }


    if (debug==1) {
      log_trace("\n\n");
      log_trace("PDM_part_extension_entity1_to_entity2 end \n");
      log_trace("========================================= \n");
    }
    // if(step == 1) {
    //   exit(1);
    // }

    PDM_free(prev_dface_itrf_blk_gnum           );
    PDM_free(prev_dface_itrf_blk_ancstr_strd    );
    PDM_free(prev_dface_itrf_blk_ancstr         );
    PDM_free(prev_dface_itrf_blk_path_itrf_strd );
    PDM_free(prev_dface_itrf_blk_path_itrf      );
    PDM_free(prev_dface_itrf_gnum_and_itrf_strid);
    PDM_free(prev_dface_itrf_gnum_and_itrf_data );
    prev_dface_itrf_n_blk               = next_dface_itrf_n_blk;
    prev_dface_itrf_blk_gnum            = next_dface_itrf_blk_gnum;
    prev_dface_itrf_blk_ancstr_strd     = next_dface_itrf_blk_ancstr_strd;
    prev_dface_itrf_blk_ancstr          = next_dface_itrf_blk_ancstr;
    prev_dface_itrf_blk_path_itrf_strd  = next_dface_itrf_blk_path_itrf_strd;
    prev_dface_itrf_blk_path_itrf       = next_dface_itrf_blk_path_itrf;
    prev_dface_itrf_gnum_and_itrf_strid = next_dface_itrf_gnum_and_itrf_strid;
    prev_dface_itrf_gnum_and_itrf_data  = next_dface_itrf_gnum_and_itrf_data;

    if(part_ext->have_edge == 1 && pface_vtx != NULL ) {
      for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
        PDM_free(pface_vtx[i_part]);
      }
    }

    /*
     * Update with descending connectivity :
     *   - Mandatory because we need to iterate the connectivity face_vtx (but with the new faces)
     */

    // > Vertex to vertex link
    int          *pn_vtx_extended                   = NULL;
    PDM_g_num_t **pvtx_extended_ln_to_gn            = NULL;
    int         **pvtx_extended_to_pvtx_idx         = NULL;
    int         **pvtx_extended_to_pvtx_triplet     = NULL;
    int         **pvtx_extended_to_pvtx_interface   = NULL;
    
    // > Edge to edge link
    int          *pn_edge_extended                  = NULL;
    PDM_g_num_t **pedge_extended_ln_to_gn           = NULL;
    int         **pedge_extended_to_pedge_idx       = NULL;
    int         **pedge_extended_to_pedge_triplet   = NULL;
    int         **pedge_extended_to_pedge_interface = NULL;
    int         **pedge_extended_to_pedge_sens      = NULL;
    
    // > Extended connectivities
    int         **pextended_face_edge_idx           = NULL;
    int         **pextended_face_edge               = NULL;
    int         **pextended_edge_vtx_idx            = NULL;
    int         **pextended_edge_vtx                = NULL;
    int         **pextended_face_vtx_idx            = NULL;
    int         **pextended_face_vtx                = NULL;

    // > Vertices next data_base 
    int          next_dvtx_itrf_n_blk               = 0;
    PDM_g_num_t *next_dvtx_itrf_blk_gnum            = NULL;
    int         *next_dvtx_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dvtx_itrf_blk_ancstr          = NULL;
    int         *next_dvtx_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dvtx_itrf_blk_path_itrf       = NULL;
    int         *next_dvtx_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dvtx_itrf_gnum_and_itrf_data  = NULL;

    // > Edges next data_base 
    int          next_dedge_itrf_n_blk               = 0;
    PDM_g_num_t *next_dedge_itrf_blk_gnum            = NULL;
    int         *next_dedge_itrf_blk_ancstr_strd     = NULL;
    PDM_g_num_t *next_dedge_itrf_blk_ancstr          = NULL;
    int         *next_dedge_itrf_blk_path_itrf_strd  = NULL;
    int         *next_dedge_itrf_blk_path_itrf       = NULL;
    int         *next_dedge_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dedge_itrf_gnum_and_itrf_data  = NULL;
    int         *next_dedge_itrf_gnum_and_itrf_sens  = NULL;


    if(part_ext->have_edge == 1) {
      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("=================================================================================\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Face->Edge)\n");
      }
      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_edge,
                                                                       prev_dedge_itrf_n_blk,
                                                                       prev_dedge_itrf_blk_gnum,
                                                                       prev_dedge_itrf_blk_ancstr_strd,
                                                                       prev_dedge_itrf_blk_ancstr,
                                                                       prev_dedge_itrf_blk_path_itrf_strd,
                                                                       prev_dedge_itrf_blk_path_itrf,
                                                                       prev_dedge_itrf_gnum_and_itrf_strid,
                                                                       prev_dedge_itrf_gnum_and_itrf_data,
                                                                       prev_dedge_itrf_gnum_and_itrf_sens,
                                                                       pn_face,
                                                                       pface_ln_to_gn,
                                                                       pn_edge,
                                                                       pedge_ln_to_gn,
                                                                       pface_edge_idx,
                                                                       pface_edge,
                                                                       pn_face_extended,
                                                                       pface_extended_ln_to_gn,
                                                                       pface_extended_to_pface_idx,
                                                                       pface_extended_to_pface_triplet,
                                                                       pface_extended_to_pface_interface,
                                                                       NULL,
                                                                       1,
                                                                       &pn_edge_extended,
                                                                       &pedge_extended_ln_to_gn,
                                                                       &pextended_face_edge_idx,
                                                                       &pextended_face_edge,
                                                                       &pedge_extended_to_pedge_idx,
                                                                       &pedge_extended_to_pedge_triplet,
                                                                       &pedge_extended_to_pedge_interface,
                                                                       &pedge_extended_to_pedge_sens,
                                                                       &next_dedge_itrf_n_blk,
                                                                       &next_dedge_itrf_blk_gnum,
                                                                       &next_dedge_itrf_blk_ancstr_strd,
                                                                       &next_dedge_itrf_blk_ancstr,
                                                                       &next_dedge_itrf_blk_path_itrf_strd,
                                                                       &next_dedge_itrf_blk_path_itrf,
                                                                       &next_dedge_itrf_gnum_and_itrf_strid,
                                                                       &next_dedge_itrf_gnum_and_itrf_data,
                                                                       &next_dedge_itrf_gnum_and_itrf_sens,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Face->Edge)\n");
        log_trace("=================================================================================\n");
      }

      PDM_free(prev_dedge_itrf_blk_gnum);
      PDM_free(prev_dedge_itrf_blk_ancstr_strd);
      PDM_free(prev_dedge_itrf_blk_ancstr);
      PDM_free(prev_dedge_itrf_blk_path_itrf_strd);
      PDM_free(prev_dedge_itrf_blk_path_itrf);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_data);
      PDM_free(prev_dedge_itrf_gnum_and_itrf_sens);
      prev_dedge_itrf_n_blk               = next_dedge_itrf_n_blk;
      prev_dedge_itrf_blk_gnum            = next_dedge_itrf_blk_gnum;
      prev_dedge_itrf_blk_ancstr_strd     = next_dedge_itrf_blk_ancstr_strd;
      prev_dedge_itrf_blk_ancstr          = next_dedge_itrf_blk_ancstr;
      prev_dedge_itrf_blk_path_itrf_strd  = next_dedge_itrf_blk_path_itrf_strd;
      prev_dedge_itrf_blk_path_itrf       = next_dedge_itrf_blk_path_itrf;
      prev_dedge_itrf_gnum_and_itrf_strid = next_dedge_itrf_gnum_and_itrf_strid;
      prev_dedge_itrf_gnum_and_itrf_data  = next_dedge_itrf_gnum_and_itrf_data;
      prev_dedge_itrf_gnum_and_itrf_sens  = next_dedge_itrf_gnum_and_itrf_sens;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_blk_gnum;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_data;
      part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_EDGE] = prev_dedge_itrf_gnum_and_itrf_sens;

      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("================================================================================\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Edge->Vtx)\n");
      }

      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_vtx,
                                                                       prev_dvtx_itrf_n_blk,
                                                                       prev_dvtx_itrf_blk_gnum,
                                                                       prev_dvtx_itrf_blk_ancstr_strd,
                                                                       prev_dvtx_itrf_blk_ancstr,
                                                                       prev_dvtx_itrf_blk_path_itrf_strd,
                                                                       prev_dvtx_itrf_blk_path_itrf,
                                                                       prev_dvtx_itrf_gnum_and_itrf_strid,
                                                                       prev_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       pn_edge,
                                                                       pedge_ln_to_gn,
                                                                       pn_vtx,
                                                                       pvtx_ln_to_gn,
                                                                       pedge_vtx_idx,
                                                                       pedge_vtx,
                                                                       pn_edge_extended,
                                                                       pedge_extended_ln_to_gn,
                                                                       pedge_extended_to_pedge_idx,
                                                                       pedge_extended_to_pedge_triplet,
                                                                       pedge_extended_to_pedge_interface,
                                                                       pedge_extended_to_pedge_sens,
                                                                       1,
                                                                      &pn_vtx_extended,
                                                                      &pvtx_extended_ln_to_gn,
                                                                      &pextended_edge_vtx_idx,
                                                                      &pextended_edge_vtx,
                                                                      &pvtx_extended_to_pvtx_idx,
                                                                      &pvtx_extended_to_pvtx_triplet,
                                                                      &pvtx_extended_to_pvtx_interface,
                                                                       NULL,
                                                                      &next_dvtx_itrf_n_blk,
                                                                      &next_dvtx_itrf_blk_gnum,
                                                                      &next_dvtx_itrf_blk_ancstr_strd,
                                                                      &next_dvtx_itrf_blk_ancstr,
                                                                      &next_dvtx_itrf_blk_path_itrf_strd,
                                                                      &next_dvtx_itrf_blk_path_itrf,
                                                                      &next_dvtx_itrf_gnum_and_itrf_strid,
                                                                      &next_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Edge->Vtx)\n");
        log_trace("================================================================================\n");
      }
      
      PDM_free(prev_dvtx_itrf_blk_gnum);
      PDM_free(prev_dvtx_itrf_blk_ancstr_strd);
      PDM_free(prev_dvtx_itrf_blk_ancstr);
      PDM_free(prev_dvtx_itrf_blk_path_itrf_strd);
      PDM_free(prev_dvtx_itrf_blk_path_itrf);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_data);
      prev_dvtx_itrf_n_blk               = next_dvtx_itrf_n_blk;
      prev_dvtx_itrf_blk_gnum            = next_dvtx_itrf_blk_gnum;
      prev_dvtx_itrf_blk_ancstr_strd     = next_dvtx_itrf_blk_ancstr_strd;
      prev_dvtx_itrf_blk_ancstr          = next_dvtx_itrf_blk_ancstr;
      prev_dvtx_itrf_blk_path_itrf_strd  = next_dvtx_itrf_blk_path_itrf_strd;
      prev_dvtx_itrf_blk_path_itrf       = next_dvtx_itrf_blk_path_itrf;
      prev_dvtx_itrf_gnum_and_itrf_strid = next_dvtx_itrf_gnum_and_itrf_strid;
      prev_dvtx_itrf_gnum_and_itrf_data  = next_dvtx_itrf_gnum_and_itrf_data;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_gnum;
      // part_ext->dentity_itrf_blk_ancstr         [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_ancstr;
      // part_ext->dentity_itrf_blk_path_itrf_strd [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf_strd;
      // part_ext->dentity_itrf_blk_path_itrf      [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_data;

    }
    else {

      if (debug==1) {
        log_trace("\n\n\n\n");
        log_trace("================================================================================ \n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 beg (Face->Vtx) \n");
      }
      PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2(part_ext->ln_part_tot,
                                                                       part_ext->n_interface,
                                                                       shift_by_domain_vtx,
                                                                       prev_dvtx_itrf_n_blk,
                                                                       prev_dvtx_itrf_blk_gnum,
                                                                       prev_dvtx_itrf_blk_ancstr_strd,
                                                                       prev_dvtx_itrf_blk_ancstr,
                                                                       prev_dvtx_itrf_blk_path_itrf_strd,
                                                                       prev_dvtx_itrf_blk_path_itrf,
                                                                       prev_dvtx_itrf_gnum_and_itrf_strid,
                                                                       prev_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       pn_face,
                                                                       pface_ln_to_gn,
                                                                       pn_vtx,
                                                                       pvtx_ln_to_gn,
                                                                       pface_vtx_idx,
                                                                       pface_vtx,
                                                                       pn_face_extended,
                                                                       pface_extended_ln_to_gn,
                                                                       pface_extended_to_pface_idx,
                                                                       pface_extended_to_pface_triplet,
                                                                       pface_extended_to_pface_interface,
                                                                       NULL,
                                                                       1,
                                                                      &pn_vtx_extended,
                                                                      &pvtx_extended_ln_to_gn,
                                                                      &pextended_face_vtx_idx,
                                                                      &pextended_face_vtx,
                                                                      &pvtx_extended_to_pvtx_idx,
                                                                      &pvtx_extended_to_pvtx_triplet,
                                                                      &pvtx_extended_to_pvtx_interface,
                                                                       NULL,
                                                                      &next_dvtx_itrf_n_blk,
                                                                      &next_dvtx_itrf_blk_gnum,
                                                                      &next_dvtx_itrf_blk_ancstr_strd,
                                                                      &next_dvtx_itrf_blk_ancstr,
                                                                      &next_dvtx_itrf_blk_path_itrf_strd,
                                                                      &next_dvtx_itrf_blk_path_itrf,
                                                                      &next_dvtx_itrf_gnum_and_itrf_strid,
                                                                      &next_dvtx_itrf_gnum_and_itrf_data,
                                                                       NULL,
                                                                       part_ext->comm);
      if (debug==1) {
        log_trace("\n\n");
        log_trace("PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2 end (Face->Vtx) \n");
        log_trace("================================================================================ \n");
      }

      PDM_free(prev_dvtx_itrf_blk_gnum           );
      PDM_free(prev_dvtx_itrf_blk_ancstr_strd    );
      PDM_free(prev_dvtx_itrf_blk_ancstr         );
      PDM_free(prev_dvtx_itrf_blk_path_itrf_strd );
      PDM_free(prev_dvtx_itrf_blk_path_itrf      );
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_strid);
      PDM_free(prev_dvtx_itrf_gnum_and_itrf_data );
      prev_dvtx_itrf_n_blk               = next_dvtx_itrf_n_blk;
      prev_dvtx_itrf_blk_gnum            = next_dvtx_itrf_blk_gnum;
      prev_dvtx_itrf_blk_ancstr_strd     = next_dvtx_itrf_blk_ancstr_strd;
      prev_dvtx_itrf_blk_ancstr          = next_dvtx_itrf_blk_ancstr;
      prev_dvtx_itrf_blk_path_itrf_strd  = next_dvtx_itrf_blk_path_itrf_strd;
      prev_dvtx_itrf_blk_path_itrf       = next_dvtx_itrf_blk_path_itrf;
      prev_dvtx_itrf_gnum_and_itrf_strid = next_dvtx_itrf_gnum_and_itrf_strid;
      prev_dvtx_itrf_gnum_and_itrf_data  = next_dvtx_itrf_gnum_and_itrf_data;

      part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_n_blk;
      part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_gnum;
      // part_ext->dentity_itrf_blk_ancstr         [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_ancstr;
      // part_ext->dentity_itrf_blk_path_itrf_strd [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf_strd;
      // part_ext->dentity_itrf_blk_path_itrf      [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_blk_path_itrf;
      part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_strid;
      part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = prev_dvtx_itrf_gnum_and_itrf_data;

    }

    /*
     * Hook coordinates
     */
    double **pvtx_extended_coords = NULL;
    _exchange_coord_and_apply_transform(part_ext,
                                        pn_vtx_extended,
                                        pvtx_extended_ln_to_gn,
                                        pn_vtx,
                                        pvtx_coords,
                                        pvtx_extended_to_pvtx_idx,
                                        pvtx_extended_to_pvtx_triplet,
                                        pvtx_extended_to_pvtx_interface,
                                        &pvtx_extended_coords);

    /*
     * Concatenate all information to continue recursion
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      if (debug==1) {
        log_trace("  -> i_part = %d \n", i_part);
      }

      /* Update size */
      pn_vtx_extended_old              [i_part]  = pfull_n_vtx_extended[i_part];
      pfull_n_vtx_extended             [i_part] +=      pn_vtx_extended[i_part];
      pn_vtx_extended_by_depth[i_depth][i_part] +=      pn_vtx_extended[i_part];
      
      if(part_ext->have_edge == 1) {
        pn_edge_extended_old              [i_part]  = pfull_n_edge_extended[i_part];
        pfull_n_edge_extended             [i_part] +=      pn_edge_extended[i_part];
        pn_edge_extended_by_depth[i_depth][i_part] +=      pn_edge_extended[i_part];
      }

      pn_face_extended_old              [i_part]  = pfull_n_face_extended[i_part];
      pfull_n_face_extended             [i_part] +=      pn_face_extended[i_part];
      pn_face_extended_by_depth[i_depth][i_part] +=      pn_face_extended[i_part];

    }


    /* Faces */
    _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                           pn_face,
                                           pn_face_extended,
                                           pface_ln_to_gn,
                                           pface_extended_ln_to_gn,
                                           &shift_by_domain_face);

    if(part_ext->have_edge == 1) {
      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_face,
                                         pn_face_extended,
                                         pface_edge_idx,
                                         pface_edge,
                                         pextended_face_edge_idx,
                                         pextended_face_edge);
    }
    else {
      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_face,
                                         pn_face_extended,
                                         pface_vtx_idx,
                                         pface_vtx,
                                         pextended_face_vtx_idx,
                                         pextended_face_vtx);
    }

    _concat_int_array_current_with_extended(part_ext->ln_part_tot,
                                            pn_face,
                                            pn_face_extended,
                                            pface_alrdy_sent,
                                            pface_extended_alrdy_sent);

    /* Edges */
    if(part_ext->have_edge == 1) {

      _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                             pn_edge,
                                             pn_edge_extended,
                                             pedge_ln_to_gn,
                                             pedge_extended_ln_to_gn,
                                             &shift_by_domain_edge);

      _concat_connectivity_with_extended(part_ext->ln_part_tot,
                                         pn_edge,
                                         pn_edge_extended,
                                         pedge_vtx_idx,
                                         pedge_vtx,
                                         pextended_edge_vtx_idx,
                                         pextended_edge_vtx);
    }


    /* Vertices */
    _concat_ln_to_gn_current_with_extended(part_ext->ln_part_tot,
                                           pn_vtx,
                                           pn_vtx_extended,
                                           pvtx_ln_to_gn,
                                           pvtx_extended_ln_to_gn,
                                           &shift_by_domain_vtx);

    _concat_coords_current_with_extended(part_ext->ln_part_tot,
                                         pn_vtx,
                                         pn_vtx_extended,
                                         pvtx_coords,
                                         pvtx_extended_coords);

    _concat_full_with_extended(part_ext->ln_part_tot,
                               pfull_n_vtx_extended,
                               pn_vtx_extended,
                               pn_vtx_extended_old,
                               pvtx_extended_ln_to_gn,
                               pvtx_extended_to_pvtx_idx,
                               pvtx_extended_to_pvtx_triplet,
                               pvtx_extended_to_pvtx_interface,
                               pfull_vtx_extended_ln_to_gn,
                               pfull_vtx_extended_to_pvtx_idx,
                               pfull_vtx_extended_to_pvtx_triplet,
                               pfull_vtx_extended_to_pvtx_interface);

    if (part_ext->have_edge==1) {
      _concat_full_with_extended(part_ext->ln_part_tot,
                                 pfull_n_edge_extended,
                                 pn_edge_extended,
                                 pn_edge_extended_old,
                                 pedge_extended_ln_to_gn,
                                 pedge_extended_to_pedge_idx,
                                 pedge_extended_to_pedge_triplet,
                                 pedge_extended_to_pedge_interface,
                                 pfull_edge_extended_ln_to_gn,
                                 pfull_edge_extended_to_pedge_idx,
                                 pfull_edge_extended_to_pedge_triplet,
                                 pfull_edge_extended_to_pedge_interface);
    }

    _concat_full_with_extended(part_ext->ln_part_tot,
                               pfull_n_face_extended,
                               pn_face_extended,
                               pn_face_extended_old,
                               pface_extended_ln_to_gn,
                               pface_extended_to_pface_idx,
                               pface_extended_to_pface_triplet,
                               pface_extended_to_pface_interface,
                               pfull_face_extended_ln_to_gn,
                               pfull_face_extended_to_pface_idx,
                               pfull_face_extended_to_pface_triplet,
                               pfull_face_extended_to_pface_interface);


    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      /*   */
      int pn_concat_vtx  = pn_vtx [i_part] + pn_vtx_extended [i_part];
      int pn_concat_edge = pn_edge[i_part];
      if (part_ext->have_edge==1) {
        pn_concat_edge += pn_edge_extended[i_part];
      }
      int pn_concat_face = pn_face[i_part] + pn_face_extended[i_part];

      int size_vtx_vtx   = (pfull_vtx_extended_to_pvtx_idx  [i_part][pn_vtx_extended_old [i_part]] + pvtx_extended_to_pvtx_idx  [i_part][pn_vtx_extended [i_part]])/3;
      int size_edge_edge = 0;
      if (part_ext->have_edge==1) {
        size_edge_edge = (pfull_edge_extended_to_pedge_idx[i_part][pn_edge_extended_old[i_part]] + pedge_extended_to_pedge_idx[i_part][pn_edge_extended[i_part]])/3;
      }
      int size_face_face = (pfull_face_extended_to_pface_idx[i_part][pn_face_extended_old[i_part]] + pface_extended_to_pface_idx[i_part][pn_face_extended[i_part]])/3;

      if(debug == 1) {
        log_trace("\n");
        log_trace("i_part = %d\n", i_part);
        log_trace("Vertices connection :: \n");
        PDM_log_trace_array_long(pfull_vtx_extended_ln_to_gn           [i_part], pfull_n_vtx_extended[i_part]  , "pfull_vtx_extended_ln_to_gn         ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_idx        [i_part], pfull_n_vtx_extended[i_part]+1, "pfull_vtx_extended_to_pvtx_idx      ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_triplet    [i_part], 3 * size_vtx_vtx              , "pfull_vtx_extended_to_pvtx_triplet  ::");
        PDM_log_trace_array_int (pfull_vtx_extended_to_pvtx_interface  [i_part],     size_vtx_vtx              , "pfull_vtx_extended_to_pvtx_interface::");
        if (part_ext->have_edge==1) {
          log_trace("Edges connection :: \n");
          PDM_log_trace_array_long(pfull_edge_extended_ln_to_gn            [i_part], pfull_n_edge_extended[i_part]  , "pfull_edge_extended_ln_to_gn         ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_idx        [i_part], pfull_n_edge_extended[i_part]+1, "pfull_edge_extended_to_pedge_idx      ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_triplet    [i_part], 3 * size_edge_edge             , "pfull_edge_extended_to_pedge_triplet  ::");
          PDM_log_trace_array_int (pfull_edge_extended_to_pedge_interface  [i_part],     size_edge_edge             , "pfull_edge_extended_to_pedge_interface::");
        }
        log_trace("Faces connection :: \n");
        PDM_log_trace_array_long(pfull_face_extended_ln_to_gn          [i_part], pfull_n_face_extended[i_part]  , "pfull_face_extended_ln_to_gn           ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_idx      [i_part], pfull_n_face_extended[i_part]+1, "pfull_face_extended_to_pface_idx       ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_triplet  [i_part], 3 * size_face_face             , "pfull_face_extended_to_pface_triplet   ::");
        PDM_log_trace_array_int (pfull_face_extended_to_pface_interface[i_part],     size_face_face             , "pfull_face_extended_to_pface_interface ::");
      }

      /*
       * Update graphe - Only extend idx since no connextion is create AT this stage
       */
      int pn_concat_entity_extended = 0;
      int pn_entity_extended        = 0;
      int pn_entity                 = 0;
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
        pn_concat_entity_extended = pn_concat_vtx;
        pn_entity_extended        = pn_vtx_extended[i_part];
        pn_entity                 = pn_vtx         [i_part];
      } else if(part_ext->extend_type == PDM_EXTEND_FROM_EDGE) {
        pn_concat_entity_extended = pn_concat_edge;
        pn_entity_extended        = pn_edge_extended[i_part];
        pn_entity                 = pn_edge         [i_part];
      } else {
        PDM_error(__FILE__, __LINE__, 0, "part_extension with extend_type %d for 3d mesh\n", part_ext->extend_type);
      }
      // log_trace("pn_concat_entity_extended = %i \n", pn_concat_entity_extended);

      /*
       * Only extend the index array since connexion is freeze for one step
       */
      PDM_realloc(pcurr_entity_bound_to_pentity_bound_idx[i_part], pcurr_entity_bound_to_pentity_bound_idx[i_part], pn_concat_entity_extended+1, int);
      for(int i = 0; i < pn_entity_extended; ++i) {
        pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity];
      }



      if(visu == 1) {

        if (part_ext->have_edge==1) {
          char filename[999];
          sprintf(filename, "out_edge_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
          PDM_vtk_write_std_elements(filename,
                                     pn_concat_vtx,
                                     pvtx_coords   [i_part],
                                     pvtx_ln_to_gn [i_part],
                                     PDM_MESH_NODAL_BAR2,
                                     pn_concat_edge,
                                     pedge_vtx     [i_part],
                                     pedge_ln_to_gn[i_part],
                                     0,
                                     NULL,
                                     NULL);
        }


        int *_pface_vtx_idx = NULL;
        int *_pface_vtx     = NULL;
        if (part_ext->have_edge==1) {
          if(debug == 1) {
            log_trace("\n");
            log_trace("i_part = %d\n", i_part);
            log_trace("pn_concat_vtx = %d\n", pn_concat_vtx);
            PDM_log_trace_array_long(pvtx_ln_to_gn[i_part], pn_concat_vtx, "pvtx_ln_to_gn ::");

            log_trace("\n");
            log_trace("FACE_EDGE\n");
            int face_edge_size = pface_edge_idx[i_part][pn_concat_face];
            log_trace("pn_concat_face = %d\n", pn_concat_face);
            PDM_log_trace_array_int (pface_edge_idx[i_part], pn_concat_face+1, "pface_edge_idx ::");
            PDM_log_trace_array_int (pface_edge    [i_part], face_edge_size  , "pface_edge     ::");
            PDM_log_trace_array_long(pface_ln_to_gn[i_part], pn_concat_face  , "pface_ln_to_gn ::");

            log_trace("\n");
            log_trace("EDGE_VTX\n");
            int edge_vtx_size = pedge_vtx_idx[i_part][pn_concat_edge];
            log_trace("pn_concat_edge = %d\n", pn_concat_edge);
            PDM_log_trace_array_int (pedge_vtx_idx [i_part], pn_concat_edge+1, "pedge_vtx_idx  ::");
            PDM_log_trace_array_int (pedge_vtx     [i_part], edge_vtx_size   , "pedge_vtx      ::");
            PDM_log_trace_array_long(pedge_ln_to_gn[i_part], pn_concat_edge  , "pedge_ln_to_gn ::");
          }

          PDM_compute_face_vtx_from_face_and_edge(pn_concat_face,
                                                  pface_edge_idx[i_part],
                                                  pface_edge    [i_part],
                                                  pedge_vtx     [i_part],
                                                &_pface_vtx);
          _pface_vtx_idx = pface_edge_idx[i_part];

        }
        else {
          _pface_vtx_idx = pface_vtx_idx[i_part];
          _pface_vtx     = pface_vtx    [i_part];
        }


        if(debug == 1) {
          log_trace("\n");
          log_trace("i_part = %d\n", i_part);
          log_trace("pn_concat_vtx = %d\n", pn_concat_vtx);
          PDM_log_trace_array_long(pvtx_ln_to_gn[i_part], pn_concat_vtx, "pvtx_ln_to_gn ::");

          log_trace("\n");
          log_trace("FACE_VTX\n");
          int face_vtx_size = _pface_vtx_idx[pn_concat_face];
          log_trace("pn_concat_face = %d\n", pn_concat_face);
          PDM_log_trace_array_int (_pface_vtx_idx         , pn_concat_face+1, "pface_vtx_idx  ::");
          PDM_log_trace_array_int (_pface_vtx             , face_vtx_size   , "pface_vtx      ::");
          PDM_log_trace_array_long( pface_ln_to_gn[i_part], pn_concat_face  , "pface_ln_to_gn ::");
        }

        PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
        char filename[999];
        sprintf(filename, "out_face_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               pn_concat_vtx,
                               pvtx_coords   [i_part],
                               pvtx_ln_to_gn [i_part],
                               pn_concat_face,
                              _pface_vtx_idx,
                              _pface_vtx,
                               pface_ln_to_gn[i_part],
                               NULL);

        if (part_ext->have_edge==1) {
          PDM_free(_pface_vtx);
        }
      } // End visu

    } /* End loop part */

    /*
     * Free coords
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(pvtx_extended_coords[i_part]);
    }
    PDM_free(pvtx_extended_coords);

    /*
     * Update shift_by_domain_face
     */
    PDM_g_num_t _shift_by_domain_face = shift_by_domain_face;
    PDM_MPI_Allreduce(&_shift_by_domain_face, &shift_by_domain_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    if (part_ext->have_edge==1) {
      PDM_g_num_t _shift_by_domain_edge = shift_by_domain_edge;
      PDM_MPI_Allreduce(&_shift_by_domain_edge, &shift_by_domain_edge, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);
    }

    PDM_g_num_t _shift_by_domain_vtx = shift_by_domain_vtx;
    PDM_MPI_Allreduce(&_shift_by_domain_vtx, &shift_by_domain_vtx, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    PDM_g_num_t _pn_face_extended_tot = 0;
    PDM_g_num_t  pn_face_extended_tot = 0;
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      _pn_face_extended_tot += pn_face_extended[i_part];
    }

    PDM_MPI_Allreduce(&_pn_face_extended_tot, &pn_face_extended_tot, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, part_ext->comm);

    // log_trace("pn_face_extended_tot = %i (local = %i ) \n", pn_face_extended_tot, _pn_face_extended_tot);

    /*
     * Free
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(pvtx_extended_ln_to_gn         [i_part]);
      PDM_free(pvtx_extended_to_pvtx_idx      [i_part]);
      PDM_free(pvtx_extended_to_pvtx_triplet  [i_part]);
      PDM_free(pvtx_extended_to_pvtx_interface[i_part]);

      if (part_ext->have_edge==1) {
        PDM_free(pedge_extended_ln_to_gn          [i_part]);
        PDM_free(pedge_extended_to_pedge_idx      [i_part]);
        PDM_free(pedge_extended_to_pedge_triplet  [i_part]);
        PDM_free(pedge_extended_to_pedge_interface[i_part]);
        if (pedge_extended_to_pedge_sens!=NULL) {
          PDM_free(pedge_extended_to_pedge_sens   [i_part]);
        }

        PDM_free(pextended_face_edge_idx[i_part]);
        PDM_free(pextended_face_edge    [i_part]);
        PDM_free(pextended_edge_vtx_idx [i_part]);
        PDM_free(pextended_edge_vtx     [i_part]);
      }
      else {
        PDM_free(pextended_face_vtx_idx[i_part]);
        PDM_free(pextended_face_vtx    [i_part]);
      }
      
      PDM_free(pface_extended_ln_to_gn          [i_part]);
      PDM_free(pface_extended_to_pface_idx      [i_part]);
      PDM_free(pface_extended_to_pface_triplet  [i_part]);
      PDM_free(pface_extended_to_pface_interface[i_part]);
      PDM_free(pface_extended_alrdy_sent        [i_part]);
    }

    PDM_free(pvtx_extended_ln_to_gn);
    PDM_free(pvtx_extended_to_pvtx_idx);
    PDM_free(pvtx_extended_to_pvtx_triplet);
    PDM_free(pvtx_extended_to_pvtx_interface);

    PDM_free(pedge_extended_ln_to_gn);
    PDM_free(pedge_extended_to_pedge_idx);
    PDM_free(pedge_extended_to_pedge_triplet);
    PDM_free(pedge_extended_to_pedge_interface);
    PDM_free(pedge_extended_to_pedge_sens);

    PDM_free(pextended_face_vtx_idx);
    PDM_free(pextended_face_vtx);
    PDM_free(pextended_face_edge_idx);
    PDM_free(pextended_face_edge);
    PDM_free(pextended_edge_vtx_idx);
    PDM_free(pextended_edge_vtx);
    
    PDM_free(pface_extended_ln_to_gn);
    PDM_free(pface_extended_to_pface_idx);
    PDM_free(pface_extended_to_pface_triplet);
    PDM_free(pface_extended_to_pface_interface);
    PDM_free(pface_extended_alrdy_sent);


    /*
     * A chaque étape :
     *   - On garde le même graphe entre les entitiés, mais on agrandit le tableau idx (pour être cohérent avec le part_to_part )
     *   - A l'issu d'une étape, il faut swap le graph avec celui de la nouvelle depth
     *   - Pour avoir l'historique complet on peut agglomerer tout les graphe de chaque depth to have the full one
     */
    if(pn_face_extended_tot == 0) {
      // Change graph
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {

        _update_propagating_graph_for_depth(i_depth,
                                            part_ext->ln_part_tot,
                                            pn_vtx,
                                            pn_vtx_extended_by_depth,
                                            pfull_vtx_extended_to_pvtx_idx,
                                            pfull_vtx_extended_to_pvtx_triplet,
                                            pfull_vtx_extended_to_pvtx_interface,
                                           &pcurr_entity_bound_to_pentity_bound_idx,
                                           &pcurr_entity_bound_to_pentity_bound_triplet,
                                           &pcurr_entity_bound_to_pentity_bound_interface);
      }
      else if(part_ext->extend_type == PDM_EXTEND_FROM_EDGE) {
        _update_propagating_graph_for_depth(i_depth,
                                            part_ext->ln_part_tot,
                                            pn_edge,
                                            pn_edge_extended_by_depth,
                                            pfull_edge_extended_to_pedge_idx,
                                            pfull_edge_extended_to_pedge_triplet,
                                            pfull_edge_extended_to_pedge_interface,
                                           &pcurr_entity_bound_to_pentity_bound_idx,
                                           &pcurr_entity_bound_to_pentity_bound_triplet,
                                           &pcurr_entity_bound_to_pentity_bound_interface);
      }

      i_depth++;
    }
    step++;


    /*
     * Update for next step
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      for (int i_face=0; i_face<pn_face[i_part]; ++i_face) {
        pface_alrdy_sent[i_part][i_face] = 0;
      }

      /* Update size */
      pn_vtx [i_part] += pn_vtx_extended [i_part];
      if (part_ext->have_edge==1) {
        pn_edge[i_part] += pn_edge_extended[i_part];
      }
      pn_face[i_part] += pn_face_extended[i_part];
    }
    PDM_free(pn_vtx_extended);
    if (part_ext->have_edge==1) {
      PDM_free(pn_edge_extended);
    }
    PDM_free(pn_face_extended);

    // if(step > 1) {
    //   abort();
    // }

    double t_end = PDM_MPI_Wtime();

    if (i_rank==0) printf(" (%.3fs)\n", t_end - t_start);
  }


  /**
   * Copy extended entities on part_extension structure
   */
  // > Faces
  if(part_ext->have_edge == 1) {
    _store_extended_entity(part_ext, 
                           pn_init_face,
                           pn_face,
                           pface_ln_to_gn,
                           pface_edge_idx,
                           pface_edge,
                          &part_ext->n_face_border,
                          &part_ext->border_face_ln_to_gn,
                          &part_ext->border_face_edge_idx,
                          &part_ext->border_face_edge,
                           PDM_MESH_ENTITY_FACE, 
                           PDM_CONNECTIVITY_TYPE_FACE_EDGE);

    _store_extended_entity(part_ext, 
                           pn_init_edge,
                           pn_edge,
                           pedge_ln_to_gn,
                           pedge_vtx_idx,
                           pedge_vtx,
                          &part_ext->n_edge_border,
                          &part_ext->border_edge_ln_to_gn,
                          &part_ext->border_edge_vtx_idx,
                          &part_ext->border_edge_vtx,
                           PDM_MESH_ENTITY_EDGE, 
                           PDM_CONNECTIVITY_TYPE_EDGE_VTX);
  }
  else {
    _store_extended_entity(part_ext, 
                           pn_init_face,
                           pn_face,
                           pface_ln_to_gn,
                           pface_vtx_idx,
                           pface_vtx,
                          &part_ext->n_face_border,
                          &part_ext->border_face_ln_to_gn,
                          &part_ext->border_face_vtx_idx,
                          &part_ext->border_face_vtx,
                           PDM_MESH_ENTITY_FACE, 
                           PDM_CONNECTIVITY_TYPE_FACE_VTX);
  }

  // _store_extended_entity(part_ext, 
  //                        pn_init_vtx,
  //                        pn_vtx,
  //                        pvtx_ln_to_gn,
  //                        NULL,
  //                        NULL,
  //                       &part_ext->n_vtx_border,
  //                       &part_ext->border_vtx_ln_to_gn,
  //                        NULL,
  //                        NULL);
  _store_extended_vtx(part_ext, 
                         pn_init_vtx,
                         pn_vtx,
                         pvtx_ln_to_gn,
                         pvtx_coords,
                        &part_ext->n_vtx_border,
                        &part_ext->border_vtx_ln_to_gn,
                        &part_ext->border_vtx);


  /**
   * Need to get distributed data_base info on partition
   */

  // > Faces
  int         **pface_ancstr_strd    = NULL;
  PDM_g_num_t **pface_ancstr         = NULL;
  // int         **pface_path_itrf_strd = NULL;
  // int         **pface_path_itrf      = NULL;
  _get_block_data_base_on_part_ext( part_ext,
                                    prev_dface_itrf_n_blk,
                                    prev_dface_itrf_blk_gnum,
                                    prev_dface_itrf_blk_ancstr_strd,
                                    prev_dface_itrf_blk_ancstr,
                                    prev_dface_itrf_blk_path_itrf_strd,
                                    prev_dface_itrf_blk_path_itrf,
                                    part_ext->n_face_border,
                                    part_ext->border_face_ln_to_gn,
                                   &pface_ancstr_strd,
                                   &pface_ancstr,
                                   // &pface_path_itrf_strd,
                                   &part_ext->face_face_path_itrf_idx,
                                   &part_ext->face_face_path_itrf,
                                    PDM_MESH_ENTITY_FACE);


  // > Edges
  int         **pedge_ancstr_strd    = NULL;
  PDM_g_num_t **pedge_ancstr         = NULL;
  // int         **pedge_path_itrf_strd = NULL;
  // int         **pedge_path_itrf      = NULL;
  if(part_ext->have_edge == 1) {
    _get_block_data_base_on_part_ext( part_ext,
                                      prev_dedge_itrf_n_blk,
                                      prev_dedge_itrf_blk_gnum,
                                      prev_dedge_itrf_blk_ancstr_strd,
                                      prev_dedge_itrf_blk_ancstr,
                                      prev_dedge_itrf_blk_path_itrf_strd,
                                      prev_dedge_itrf_blk_path_itrf,
                                      part_ext->n_edge_border,
                                      part_ext->border_edge_ln_to_gn,
                                     &pedge_ancstr_strd,
                                     &pedge_ancstr,
                                     // &pedge_path_itrf_strd,
                                     &part_ext->edge_edge_path_itrf_idx,
                                     &part_ext->edge_edge_path_itrf,
                                      PDM_MESH_ENTITY_EDGE);
  }

  // > Vertex
  int         **pvtx_ancstr_strd    = NULL;
  PDM_g_num_t **pvtx_ancstr         = NULL;
  // int         **pvtx_path_itrf_strd = NULL;
  // int         **pvtx_path_itrf      = NULL;
  _get_block_data_base_on_part_ext( part_ext,
                                    prev_dvtx_itrf_n_blk,
                                    prev_dvtx_itrf_blk_gnum,
                                    prev_dvtx_itrf_blk_ancstr_strd,
                                    prev_dvtx_itrf_blk_ancstr,
                                    prev_dvtx_itrf_blk_path_itrf_strd,
                                    prev_dvtx_itrf_blk_path_itrf,
                                    part_ext->n_vtx_border,
                                    part_ext->border_vtx_ln_to_gn,
                                   &pvtx_ancstr_strd,
                                   &pvtx_ancstr,
                                   // &pvtx_path_itrf_strd,
                                   &part_ext->vtx_vtx_path_itrf_idx,
                                   &part_ext->vtx_vtx_path_itrf,
                                    PDM_MESH_ENTITY_VTX);


  /**
   * Use ancestor info to build graph between original partition and extended partition
   */
  // > Faces
  _build_part_extension_graph_to_old(part_ext,
                                     pn_init_face,
                                     pface_ln_to_gn,
                                     NULL,
                                     NULL,
                                     NULL,
                                     part_ext->n_face_border,
                                     part_ext->border_face_ln_to_gn,
                                     pface_ancstr_strd,
                                     pface_ancstr,
                                     // pface_path_itrf_strd,
                                     // pface_path_itrf,
                                     part_ext->shift_by_domain_face,
                                    &part_ext->border_face_ln_to_gn_ancstr,
                                    &part_ext->face_face_extended,
                                     NULL,
                                     NULL,
                                     NULL,
                                     PDM_MESH_ENTITY_FACE);

  // > Edges
  if(part_ext->have_edge == 1) {
    _build_part_extension_graph_to_old(part_ext,
                                       pn_init_edge,
                                       pedge_ln_to_gn,
                                       pn_edge_group,
                                       pedge_group_tag,
                                       pedge_group_gnum,
                                       part_ext->n_edge_border,
                                       part_ext->border_edge_ln_to_gn,
                                       pedge_ancstr_strd,
                                       pedge_ancstr,
                                       // pedge_path_itrf_strd,
                                       // pedge_path_itrf,
                                       part_ext->shift_by_domain_edge,
                                      &part_ext->border_edge_ln_to_gn_ancstr,
                                      &part_ext->edge_edge_extended,
                                      &part_ext->border_edge_group_idx,
                                      &part_ext->border_edge_group,
                                      &part_ext->border_edge_group_ln_to_gn,
                                       PDM_MESH_ENTITY_EDGE);
  }

  // > Vertex
  _build_part_extension_graph_to_old(part_ext,
                                     pn_init_vtx,
                                     pvtx_ln_to_gn,
                                     NULL,
                                     NULL,
                                     NULL,
                                     part_ext->n_vtx_border,
                                     part_ext->border_vtx_ln_to_gn,
                                     pvtx_ancstr_strd,
                                     pvtx_ancstr,
                                     // pvtx_path_itrf_strd,
                                     // pvtx_path_itrf,
                                     part_ext->shift_by_domain_vtx,
                                    &part_ext->border_vtx_ln_to_gn_ancstr,
                                    &part_ext->vtx_vtx_extended,
                                     NULL,
                                     NULL,
                                     NULL,
                                     PDM_MESH_ENTITY_VTX);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pface_ancstr_strd   [i_part]);
    PDM_free(pface_ancstr        [i_part]);
    // PDM_free(pface_path_itrf_strd[i_part]);
    // PDM_free(pface_path_itrf     [i_part]);

    if(part_ext->have_edge == 1) {
      PDM_free(pedge_ancstr_strd   [i_part]);
      PDM_free(pedge_ancstr        [i_part]);
      // PDM_free(pedge_path_itrf_strd[i_part]);
      // PDM_free(pedge_path_itrf     [i_part]);
    }

    PDM_free(pvtx_ancstr_strd    [i_part]);
    PDM_free(pvtx_ancstr         [i_part]);
    // PDM_free(pvtx_path_itrf_strd [i_part]);
    // PDM_free(pvtx_path_itrf      [i_part]);
  }
  PDM_free(pface_ancstr_strd);
  PDM_free(pface_ancstr);
  // PDM_free(pface_path_itrf_strd);
  // PDM_free(pface_path_itrf);
  
  if(part_ext->have_edge == 1) {
    PDM_free(pedge_ancstr_strd);
    PDM_free(pedge_ancstr);
    // PDM_free(pedge_path_itrf_strd);
    // PDM_free(pedge_path_itrf);
  }
  
  PDM_free(pvtx_ancstr_strd);
  PDM_free(pvtx_ancstr);
  // PDM_free(pvtx_path_itrf_strd);
  // PDM_free(pvtx_path_itrf);



  PDM_free(prev_dface_itrf_blk_gnum           );
  PDM_free(prev_dface_itrf_blk_ancstr_strd    );
  PDM_free(prev_dface_itrf_blk_ancstr         );
  PDM_free(prev_dface_itrf_blk_path_itrf_strd );
  PDM_free(prev_dface_itrf_blk_path_itrf      );
  PDM_free(prev_dface_itrf_gnum_and_itrf_strid);
  PDM_free(prev_dface_itrf_gnum_and_itrf_data );

  if (part_ext->have_edge==1) {
    PDM_free(prev_dedge_itrf_blk_ancstr_strd    );
    PDM_free(prev_dedge_itrf_blk_ancstr         );
    PDM_free(prev_dedge_itrf_blk_path_itrf_strd );
    PDM_free(prev_dedge_itrf_blk_path_itrf      );
  }
  
  PDM_free(prev_dvtx_itrf_blk_ancstr_strd);
  PDM_free(prev_dvtx_itrf_blk_ancstr);
  PDM_free(prev_dvtx_itrf_blk_path_itrf_strd);
  PDM_free(prev_dvtx_itrf_blk_path_itrf);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pvtx_ln_to_gn      [i_part]);
    PDM_free(pface_ln_to_gn     [i_part]);
    PDM_free(pface_alrdy_sent   [i_part]);
    if (part_ext->have_edge==1) {
      PDM_free(pedge_ln_to_gn   [i_part]);
      PDM_free(pface_edge_idx   [i_part]);
      PDM_free(pface_edge       [i_part]);
      PDM_free(pedge_vtx_idx    [i_part]);
      PDM_free(pedge_vtx        [i_part]);
      PDM_free(pedge_group_tag  [i_part]);
      PDM_free(pedge_group_gnum [i_part]);
    }
    else {
      PDM_free(pface_vtx_idx    [i_part]);
      PDM_free(pface_vtx        [i_part]);
    }
    PDM_free(pvtx_coords        [i_part]);
  }

  PDM_free(pn_vtx);
  PDM_free(pn_init_vtx);
  PDM_free(pvtx_ln_to_gn);
  PDM_free(pn_edge);
  PDM_free(pn_init_edge);
  PDM_free(pedge_ln_to_gn);
  PDM_free(pn_face);
  PDM_free(pn_init_face);
  PDM_free(pface_ln_to_gn);
  PDM_free(pface_alrdy_sent);
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge);
  PDM_free(pface_vtx_idx);
  PDM_free(pface_vtx);
  PDM_free(pedge_vtx_idx);
  PDM_free(pedge_vtx);
  PDM_free(pn_edge_group);
  PDM_free(pedge_group_tag);
  PDM_free(pedge_group_gnum);
  PDM_free(pvtx_coords);


  for(int i = 0; i < part_ext->depth; ++i) {
    PDM_free(pn_vtx_extended_by_depth [i]);
    if (part_ext->have_edge==1) {
      PDM_free(pn_edge_extended_by_depth[i]);
    }
    PDM_free(pn_face_extended_by_depth[i]);
  }
  PDM_free(pn_vtx_extended_by_depth );
  PDM_free(pn_edge_extended_by_depth);
  PDM_free(pn_face_extended_by_depth);

  /*
   * To keep
   */
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pfull_vtx_extended_ln_to_gn         [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_idx      [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_triplet  [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_interface[i_part]);

    if (part_ext->have_edge==1) {
      PDM_free(pfull_edge_extended_ln_to_gn          [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_idx      [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_triplet  [i_part]);
      PDM_free(pfull_edge_extended_to_pedge_interface[i_part]);
    }

    PDM_free(pfull_face_extended_ln_to_gn          [i_part]);
    PDM_free(pfull_face_extended_to_pface_idx      [i_part]);
    PDM_free(pfull_face_extended_to_pface_triplet  [i_part]);
    PDM_free(pfull_face_extended_to_pface_interface[i_part]);

  }
  PDM_free(pn_vtx_extended_old                 );
  PDM_free(pfull_n_vtx_extended                );
  PDM_free(pfull_vtx_extended_ln_to_gn         );
  PDM_free(pfull_vtx_extended_to_pvtx_idx      );
  PDM_free(pfull_vtx_extended_to_pvtx_triplet  );
  PDM_free(pfull_vtx_extended_to_pvtx_interface);

  PDM_free(pn_edge_extended_old                  );
  PDM_free(pfull_n_edge_extended                 );
  PDM_free(pfull_edge_extended_ln_to_gn          );
  PDM_free(pfull_edge_extended_to_pedge_idx      );
  PDM_free(pfull_edge_extended_to_pedge_triplet  );
  PDM_free(pfull_edge_extended_to_pedge_interface);

  PDM_free(pn_face_extended_old);
  PDM_free(pfull_n_face_extended                 );
  PDM_free(pfull_face_extended_ln_to_gn          );
  PDM_free(pfull_face_extended_to_pface_idx      );
  PDM_free(pfull_face_extended_to_pface_triplet  );
  PDM_free(pfull_face_extended_to_pface_interface);

}

static
void
_part_extension_1d
(
 PDM_part_extension_t *part_ext
)
{
  /* Size */
  int *pn_vtx  = NULL;
  int *pn_edge = NULL;
  PDM_malloc(pn_vtx , part_ext->ln_part_tot, int );
  PDM_malloc(pn_edge, part_ext->ln_part_tot, int );

  /* ln_to_gn */
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_malloc(pvtx_ln_to_gn , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pedge_ln_to_gn, part_ext->ln_part_tot, PDM_g_num_t *);

  /* Connectivity */
  int **pedge_vtx_idx  = NULL;
  int **pedge_vtx      = NULL;
  PDM_malloc(pedge_vtx_idx, part_ext->ln_part_tot, int *);
  PDM_malloc(pedge_vtx    , part_ext->ln_part_tot, int *);

  /* Coordinates */
  double **pvtx_coords = NULL;
  PDM_malloc(pvtx_coords, part_ext->ln_part_tot, double *);

  int lpart = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      pn_vtx        [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge       [lpart] = part_ext->parts[i_domain][i_part].n_edge;

      /* Copy to realloc after all step */
      PDM_malloc(pvtx_ln_to_gn [lpart], pn_vtx [lpart], PDM_g_num_t);
      PDM_malloc(pedge_ln_to_gn[lpart], pn_edge[lpart], PDM_g_num_t);

      /* ln_to_gn */
      for(int i_edge = 0; i_edge < pn_edge[lpart]; ++i_edge) {
        pedge_ln_to_gn[lpart][i_edge] = part_ext->parts[i_domain][i_part].edge_ln_to_gn[i_edge];
      }

      for(int i_vtx = 0; i_vtx < pn_vtx[lpart]; ++i_vtx) {
        pvtx_ln_to_gn[lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn[i_vtx];
      }

      PDM_malloc(pedge_vtx_idx[lpart],   pn_edge[lpart]+1 , int);
      PDM_malloc(pedge_vtx    [lpart], 2*pn_edge[lpart]   , int);
      for(int i_edge = 0; i_edge < pn_edge[lpart]+1; ++i_edge) {
        pedge_vtx_idx [lpart][i_edge] = 2 * i_edge;
      }

      for(int idx = 0; idx < pedge_vtx_idx[i_part][pn_edge[lpart]]; ++idx) {
        pedge_vtx [lpart][idx] = part_ext->parts[i_domain][i_part].edge_vtx[idx];
      }

      PDM_malloc(pvtx_coords[lpart], 3 * pn_vtx [lpart], double);
      for(int i_vtx = 0; i_vtx < 3 * pn_vtx[lpart]; ++i_vtx) {
        pvtx_coords   [lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx[i_vtx];
      }

      lpart++;
    }
  }

  int           *pn_vtx_extended_old                  = NULL;
  int           *pfull_n_vtx_extended                 = NULL;
  PDM_g_num_t  **pfull_vtx_extended_ln_to_gn          = NULL;
  int          **pfull_vtx_extended_to_pvtx_idx       = NULL;
  int          **pfull_vtx_extended_to_pvtx_triplet   = NULL;
  int          **pfull_vtx_extended_to_pvtx_interface = NULL;
  PDM_malloc(pn_vtx_extended_old                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_vtx_extended                , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_vtx_extended_ln_to_gn         , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_vtx_extended_to_pvtx_interface, part_ext->ln_part_tot, int         *);

  int           *pn_edge_extended_old                   = NULL;
  int           *pfull_n_edge_extended                  = NULL;
  PDM_g_num_t  **pfull_edge_extended_ln_to_gn           = NULL;
  int          **pfull_edge_extended_to_pedge_idx       = NULL;
  int          **pfull_edge_extended_to_pedge_triplet   = NULL;
  int          **pfull_edge_extended_to_pedge_interface = NULL;
  PDM_malloc(pn_edge_extended_old                  , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_n_edge_extended                 , part_ext->ln_part_tot, int          );
  PDM_malloc(pfull_edge_extended_ln_to_gn          , part_ext->ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pfull_edge_extended_to_pedge_idx      , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_triplet  , part_ext->ln_part_tot, int         *);
  PDM_malloc(pfull_edge_extended_to_pedge_interface, part_ext->ln_part_tot, int         *);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_malloc(pfull_vtx_extended_to_pvtx_idx      [i_part], pfull_n_vtx_extended[i_part]+1, int);
    PDM_malloc(pfull_vtx_extended_ln_to_gn         [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_vtx_extended_to_pvtx_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_vtx_extended_to_pvtx_interface[i_part], 0, int        );
    pfull_n_vtx_extended          [i_part] = 0;
    pfull_vtx_extended_to_pvtx_idx[i_part][0] = 0;

    PDM_malloc(pfull_edge_extended_to_pedge_idx      [i_part], pfull_n_edge_extended[i_part]+1, int);
    PDM_malloc(pfull_edge_extended_ln_to_gn          [i_part], 0, PDM_g_num_t);
    PDM_malloc(pfull_edge_extended_to_pedge_triplet  [i_part], 0, int        );
    PDM_malloc(pfull_edge_extended_to_pedge_interface[i_part], 0, int        );
    pfull_n_edge_extended           [i_part] = 0;
    pfull_edge_extended_to_pedge_idx[i_part][0] = 0;
  }



  /*
   * Free all
   */
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pfull_vtx_extended_to_pvtx_idx        [i_part]);
    PDM_free(pfull_edge_extended_to_pedge_idx      [i_part]);
    PDM_free(pfull_vtx_extended_ln_to_gn           [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_triplet    [i_part]);
    PDM_free(pfull_vtx_extended_to_pvtx_interface  [i_part]);
    PDM_free(pfull_edge_extended_ln_to_gn          [i_part]);
    PDM_free(pfull_edge_extended_to_pedge_triplet  [i_part]);
    PDM_free(pfull_edge_extended_to_pedge_interface[i_part]);

  }
  PDM_free(pn_vtx_extended_old                 );
  PDM_free(pfull_n_vtx_extended                );
  PDM_free(pfull_vtx_extended_ln_to_gn         );
  PDM_free(pfull_vtx_extended_to_pvtx_idx      );
  PDM_free(pfull_vtx_extended_to_pvtx_triplet  );
  PDM_free(pfull_vtx_extended_to_pvtx_interface);
  PDM_free(pn_edge_extended_old                 );
  PDM_free(pfull_n_edge_extended                );
  PDM_free(pfull_edge_extended_ln_to_gn         );
  PDM_free(pfull_edge_extended_to_pedge_idx      );
  PDM_free(pfull_edge_extended_to_pedge_triplet  );
  PDM_free(pfull_edge_extended_to_pedge_interface);

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    PDM_free(pvtx_ln_to_gn [i_part]);
    PDM_free(pedge_ln_to_gn[i_part]);
    PDM_free(pedge_vtx_idx [i_part]);
    PDM_free(pedge_vtx     [i_part]);
    PDM_free(pvtx_coords   [i_part]);
  }

  PDM_free(pn_vtx        );
  PDM_free(pn_edge       );
  PDM_free(pvtx_ln_to_gn );
  PDM_free(pedge_ln_to_gn);
  PDM_free(pedge_vtx_idx );
  PDM_free(pedge_vtx     );
  PDM_free(pvtx_coords   );
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_part_extension_compute2
(
        PDM_part_extension_t *part_ext,
  const int                   dim
)
{
  // TODO : mv dim in create but break API
  part_ext->dim          = dim;
  part_ext->compute_kind = 1;

  part_ext->ln_part_tot = 0;
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < part_ext->n_part[i_dom]; ++i_part) {
      part_ext->ln_part_tot += 1;
    }
  }

  part_ext->lpart_to_dom = NULL;
  PDM_malloc(part_ext->lpart_to_dom, part_ext->ln_part_tot, int);
  int l_part = 0;
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < part_ext->n_part[i_dom]; ++i_part) {
      part_ext->lpart_to_dom[l_part] = i_dom;
      l_part += 1;
    }
  }

  _compute_offset(part_ext);

  /*
   * Warm up all part_domain_interface
   *   - Identify all domain interface we have and make missing one
   *   - Go to distributed frame
   *   - Create the extend_type part_domain_interface to setup the first step of the algorithme
   */
  _compute_other_part_domain_interface(part_ext);

  /*
   * Let's do an block frame domain_interface
   */
  _setup_domain_interface_in_block_frame(part_ext);

  /* Manage shift */
  _offset_parts_by_domain(part_ext, 1);

  /*
   * Create the initial graphe with the extend_kind
   */
  _build_bound_graph(part_ext);



  /* Manage dim */
  if(dim == 3) {
    _part_extension_3d(part_ext);
  } else if(dim == 2) {
    _part_extension_2d(part_ext);
  } else if(dim == 1) {
    _part_extension_1d(part_ext);
  } else  {
    PDM_error(__FILE__, __LINE__, 0, "Wrong dim size in PDM_part_extension_compute2 : %d ( Should be >=1 )\n", (int) dim);
  }


  /* Manage loop for depth AND multiple transformation */

  /* Manage unshift */
  _offset_parts_by_domain(part_ext, -1);
  // _offset_results_by_domain(part_ext);


  /* Free data - Should be usefull to redo all part_interface after all the process */
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    if(part_ext->dentity_itrf_blk_gnum[i] != NULL) {
      PDM_free(part_ext->dentity_itrf_blk_gnum           [i]);
      PDM_free(part_ext->dentity_itrf_gnum_and_itrf_strid[i]);
      PDM_free(part_ext->dentity_itrf_gnum_and_itrf_data [i]);
      // PDM_free(part_ext->dentity_ancstr_strd[i]);
      // PDM_free(part_ext->dentity_ancstr     [i]);
      // PDM_free(part_ext->dentity_path_itrf_strd[i]);
      // PDM_free(part_ext->dentity_path_itrf     [i]);
      PDM_free(part_ext->dentity_itrf_gnum_and_itrf_sens [i]);
    }
  }

  if(part_ext->dom_itrf != NULL) {
    PDM_domain_interface_free(part_ext->dom_itrf);
  }


}

/**
 *
 * \brief Free a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_part_extension_free
(
 PDM_part_extension_t *part_ext
)
{
  if (part_ext == NULL) {
    return;
  }

  if(part_ext->n_tot_part_by_domain != NULL) {
    PDM_free(part_ext->n_tot_part_by_domain);
    part_ext->n_tot_part_by_domain = NULL;
  }

  if(part_ext->neighbor_idx != NULL) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        PDM_free(part_ext->neighbor_idx       [i_part+shift_part]);
        PDM_free(part_ext->neighbor_desc      [i_part+shift_part]);
        PDM_free(part_ext->neighbor_interface [i_part+shift_part]);

        PDM_free(part_ext->dist_neighbor_cell_n           [i_part+shift_part]);
        PDM_free(part_ext->dist_neighbor_cell_idx         [i_part+shift_part]);
        PDM_free(part_ext->dist_neighbor_cell_desc        [i_part+shift_part]);
        PDM_free(part_ext->dist_neighbor_cell_interface   [i_part+shift_part]);
        PDM_free(part_ext->unique_order_dist_neighbor_cell[i_part+shift_part]);

        PDM_free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        PDM_free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        PDM_free(part_ext->entity_cell_opp    [i_part+shift_part]);

        PDM_free(part_ext->border_cell_list    [i_part+shift_part]);

        // if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        PDM_free(part_ext->entity_cell_idx[i_part+shift_part]);
        PDM_free(part_ext->entity_cell_n  [i_part+shift_part]);
        PDM_free(part_ext->entity_cell    [i_part+shift_part]);
        // }

        PDM_free(part_ext->opp_interface_and_gnum_face[i_part+shift_part]);
        PDM_free(part_ext->opp_interface_and_gnum_edge[i_part+shift_part]);
        PDM_free(part_ext->opp_interface_and_gnum_vtx[i_part+shift_part]);

        PDM_free(part_ext->cur_interface_face[i_part+shift_part]);
        PDM_free(part_ext->cur_interface_edge[i_part+shift_part]);
        PDM_free(part_ext->cur_interface_vtx[i_part+shift_part]);

        PDM_free(part_ext->cur_sens_face[i_part+shift_part]);
        PDM_free(part_ext->cur_sens_edge[i_part+shift_part]);
        PDM_free(part_ext->cur_sens_vtx [i_part+shift_part]);

        if (part_ext->face_face_extended_idx!=NULL) {
          PDM_free(part_ext->face_face_extended_idx[i_part+shift_part]);
          PDM_free(part_ext->face_face_extended    [i_part+shift_part]);
        }
        if (part_ext->edge_edge_extended_idx!=NULL) {
          PDM_free(part_ext->edge_edge_extended_idx[i_part+shift_part]);
          PDM_free(part_ext->edge_edge_extended    [i_part+shift_part]);
        }
        if (part_ext->vtx_vtx_extended_idx!=NULL) {
          PDM_free(part_ext->vtx_vtx_extended_idx  [i_part+shift_part]);
          PDM_free(part_ext->vtx_vtx_extended      [i_part+shift_part]);
        }

        if(part_ext->owner_vtx_part_bound == 1) {
          PDM_free(part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx);
          PDM_free(part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx);
          PDM_free(part_ext->parts[i_domain][i_part].vtx_part_bound);
        }

        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {

          if(part_ext->border_cell_face_idx != NULL) {
            PDM_free(part_ext->border_cell_face_idx[i_part+shift_part]);
            PDM_free(part_ext->border_cell_face    [i_part+shift_part]);
          }

          if(part_ext->border_face_edge_idx != NULL) {
            PDM_free(part_ext->border_face_edge_idx [i_part+shift_part]);
            PDM_free(part_ext->border_face_edge     [i_part+shift_part]);
          }

          if(part_ext->border_edge_vtx_idx != NULL) {
            PDM_free(part_ext->border_edge_vtx_idx[i_part+shift_part]);
            PDM_free(part_ext->border_edge_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_vtx_idx != NULL) {
            PDM_free(part_ext->border_face_vtx_idx[i_part+shift_part]);
            PDM_free(part_ext->border_face_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_group_idx != NULL) {
            PDM_free(part_ext->border_face_group_idx[i_part+shift_part]);
            PDM_free(part_ext->border_face_group    [i_part+shift_part]);
          }

          if(part_ext->border_vtx != NULL) {
            PDM_free(part_ext->border_vtx[i_part+shift_part]);
          }

          if(part_ext->border_cell_ln_to_gn != NULL) {
            PDM_free(part_ext->border_cell_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_face_ln_to_gn != NULL) {
            PDM_free(part_ext->border_face_ln_to_gn[i_part+shift_part]);
            PDM_free(part_ext->face_face_interface  [i_part+shift_part]);
          }

          if(part_ext->border_edge_ln_to_gn != NULL) {
            PDM_free(part_ext->border_edge_ln_to_gn[i_part+shift_part]);
            PDM_free(part_ext->edge_edge_interface[i_part+shift_part]);
          }

          if(part_ext->border_vtx_ln_to_gn != NULL) {
            PDM_free(part_ext->border_vtx_ln_to_gn[i_part+shift_part]);
            PDM_free(part_ext->vtx_vtx_interface  [i_part+shift_part]);
          }

          if(part_ext->border_face_group_ln_to_gn != NULL) {
            PDM_free(part_ext->border_face_group_ln_to_gn[i_part+shift_part]);
          }

        }

        PDM_free(part_ext->cell_cell_idx[i_part+shift_part]);
        PDM_free(part_ext->cell_cell    [i_part+shift_part]);

      }

      shift_part += part_ext->n_part[i_domain];
    }
    PDM_free(part_ext->n_entity_bound);
  }
  PDM_free(part_ext->neighbor_idx);
  PDM_free(part_ext->neighbor_desc);
  PDM_free(part_ext->neighbor_interface);

  if( (part_ext->cell_cell_extended != NULL ||
      part_ext->face_face_extended != NULL ||
      part_ext->edge_edge_extended != NULL ||
      part_ext->  vtx_vtx_extended != NULL) && part_ext->compute_kind != 0 ) {
    int shift_part = 0;

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

        if (part_ext->ownership_border_graph[PDM_MESH_ENTITY_CELL][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          // PDM_free(part_ext->cell_cell_extended_idx [i_part+shift_part]);
          PDM_free(part_ext->cell_cell_extended2    [i_part+shift_part]);
        }
        if (part_ext->ownership_border_path_itrf[PDM_MESH_ENTITY_CELL][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->cell_cell_path_itrf_idx[i_part+shift_part]);
          PDM_free(part_ext->cell_cell_path_itrf    [i_part+shift_part]);
        }

        if (part_ext->ownership_border_graph[PDM_MESH_ENTITY_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          // PDM_free(part_ext->face_face_extended_idx [i_part+shift_part]);
          PDM_free(part_ext->face_face_extended     [i_part+shift_part]);
        }
        if (part_ext->ownership_border_path_itrf[PDM_MESH_ENTITY_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->face_face_path_itrf_idx[i_part+shift_part]);
          PDM_free(part_ext->face_face_path_itrf    [i_part+shift_part]);
        }

        if (part_ext->ownership_border_graph[PDM_MESH_ENTITY_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          // PDM_free(part_ext->edge_edge_extended_idx [i_part+shift_part]);
          PDM_free(part_ext->edge_edge_extended     [i_part+shift_part]);
        }
        if (part_ext->ownership_border_path_itrf[PDM_MESH_ENTITY_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->edge_edge_path_itrf_idx[i_part+shift_part]);
          PDM_free(part_ext->edge_edge_path_itrf    [i_part+shift_part]);
        }

        if (part_ext->ownership_border_graph[PDM_MESH_ENTITY_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          // PDM_free(part_ext->vtx_vtx_extended_idx [i_part+shift_part]);
          PDM_free(part_ext->vtx_vtx_extended     [i_part+shift_part]);
        }
        if (part_ext->ownership_border_path_itrf[PDM_MESH_ENTITY_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->vtx_vtx_path_itrf_idx[i_part+shift_part]);
          PDM_free(part_ext->vtx_vtx_path_itrf    [i_part+shift_part]);
        }
      }
      shift_part += part_ext->n_part[i_domain];
    }
    // PDM_free(part_ext->cell_cell_extended_idx);
    PDM_free(part_ext->cell_cell_extended2);
    PDM_free(part_ext->cell_cell_path_itrf_idx);
    PDM_free(part_ext->cell_cell_path_itrf);

    // PDM_free(part_ext->face_face_extended_idx);
    PDM_free(part_ext->face_face_extended);
    PDM_free(part_ext->face_face_path_itrf_idx);
    PDM_free(part_ext->face_face_path_itrf);

    // PDM_free(part_ext->edge_edge_extended_idx);
    PDM_free(part_ext->edge_edge_extended);
    PDM_free(part_ext->edge_edge_path_itrf_idx);
    PDM_free(part_ext->edge_edge_path_itrf);

    // PDM_free(part_ext->vtx_vtx_extended_idx);
    PDM_free(part_ext->vtx_vtx_extended);
    PDM_free(part_ext->vtx_vtx_path_itrf_idx);
    PDM_free(part_ext->vtx_vtx_path_itrf);
  }


  if (part_ext->n_vtx_border!=NULL) {    
    int shift_part = 0;

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        
        // > Cell
        if (part_ext->ownership_border_ln_to_gn[PDM_MESH_ENTITY_CELL][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_cell_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_cell_face_idx[i_part+shift_part]);
          PDM_free(part_ext->border_cell_face    [i_part+shift_part]);
        }
        if (part_ext->ownership_border_ln_to_gn_ancstr[PDM_MESH_ENTITY_CELL][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_cell_ln_to_gn_ancstr[i_part+shift_part]);
        }

        // > Face
        if (part_ext->ownership_border_ln_to_gn[PDM_MESH_ENTITY_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_face_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_face_edge_idx[i_part+shift_part]);
          PDM_free(part_ext->border_face_edge    [i_part+shift_part]);
        }
        if (part_ext->ownership_border_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_face_vtx_idx[i_part+shift_part]);
          PDM_free(part_ext->border_face_vtx    [i_part+shift_part]);
        }
        if (part_ext->ownership_border_group[PDM_MESH_ENTITY_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_face_group_idx     [i_part+shift_part]);
          PDM_free(part_ext->border_face_group         [i_part+shift_part]);
          PDM_free(part_ext->border_face_group_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_ln_to_gn_ancstr[PDM_MESH_ENTITY_FACE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_face_ln_to_gn_ancstr[i_part+shift_part]);
        }

        // > Edge
        if (part_ext->ownership_border_ln_to_gn[PDM_MESH_ENTITY_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_edge_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_edge_vtx_idx[i_part+shift_part]);
          PDM_free(part_ext->border_edge_vtx    [i_part+shift_part]);
        }
        if (part_ext->ownership_border_group[PDM_MESH_ENTITY_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_edge_group_idx     [i_part+shift_part]);
          PDM_free(part_ext->border_edge_group         [i_part+shift_part]);
          PDM_free(part_ext->border_edge_group_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_ln_to_gn_ancstr[PDM_MESH_ENTITY_EDGE][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_edge_ln_to_gn_ancstr[i_part+shift_part]);
        }

        // > Vtx
        if (part_ext->ownership_border_ln_to_gn[PDM_MESH_ENTITY_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_vtx_ln_to_gn[i_part+shift_part]);
        }
        if (part_ext->ownership_border_vtx_coord[i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_vtx[i_part+shift_part]);
        }
        if (part_ext->ownership_border_ln_to_gn_ancstr[PDM_MESH_ENTITY_VTX][i_domain][i_part]==PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->border_vtx_ln_to_gn_ancstr[i_part+shift_part]);
        }

      }
      shift_part += part_ext->n_part[i_domain];
    }

    // > Cell
    PDM_free(part_ext->n_cell_border);
    PDM_free(part_ext->border_cell_ln_to_gn);
    PDM_free(part_ext->border_cell_ln_to_gn_ancstr);
    PDM_free(part_ext->border_cell_face_idx);
    PDM_free(part_ext->border_cell_face    );

    // > Face
    PDM_free(part_ext->n_face_border);
    PDM_free(part_ext->border_face_ln_to_gn);
    PDM_free(part_ext->border_face_ln_to_gn_ancstr);
    PDM_free(part_ext->border_face_edge_idx);
    PDM_free(part_ext->border_face_edge);
    PDM_free(part_ext->border_face_vtx_idx);
    PDM_free(part_ext->border_face_vtx);
    PDM_free(part_ext->border_face_group_idx);
    PDM_free(part_ext->border_face_group);
    PDM_free(part_ext->border_face_group_ln_to_gn);

    // > Edge
    PDM_free(part_ext->n_edge_border);
    PDM_free(part_ext->border_edge_ln_to_gn);
    PDM_free(part_ext->border_edge_ln_to_gn_ancstr);
    PDM_free(part_ext->border_edge_vtx_idx);
    PDM_free(part_ext->border_edge_vtx);
    PDM_free(part_ext->border_edge_group_idx);
    PDM_free(part_ext->border_edge_group);
    PDM_free(part_ext->border_edge_group_ln_to_gn);

    // > Vtx
    PDM_free(part_ext->n_vtx_border);
    PDM_free(part_ext->border_vtx_ln_to_gn);
    PDM_free(part_ext->border_vtx_ln_to_gn_ancstr);
    PDM_free(part_ext->border_vtx);

    PDM_free(part_ext->lpart_to_dom);


  }



  if(part_ext->compute_kind == 0) {
    PDM_free(part_ext->n_cell);
    PDM_free(part_ext->n_cell_border);
    PDM_free(part_ext->border_cell_list);

    PDM_free(part_ext->cell_cell_idx);
    PDM_free(part_ext->cell_cell);

    if(part_ext->face_face_extended_idx != NULL) {
      PDM_free(part_ext->face_face_extended);
      PDM_free(part_ext->face_face_extended_idx);
      PDM_free(part_ext->face_face_interface);
    }

    if(part_ext->opp_interface_and_gnum_face != NULL) {
      PDM_free(part_ext->opp_interface_and_gnum_face);
      PDM_free(part_ext->cur_interface_face);
      PDM_free(part_ext->cur_sens_face);
      PDM_free(part_ext->n_cur_interface_face);
    }
    if(part_ext->opp_interface_and_gnum_edge != NULL) {
      PDM_free(part_ext->opp_interface_and_gnum_edge);
      PDM_free(part_ext->cur_interface_edge);
      PDM_free(part_ext->cur_sens_edge);
      PDM_free(part_ext->n_cur_interface_edge);
    }
    if(part_ext->opp_interface_and_gnum_vtx != NULL) {
      PDM_free(part_ext->opp_interface_and_gnum_vtx);
      PDM_free(part_ext->cur_interface_vtx);
      PDM_free(part_ext->cur_sens_vtx);
      PDM_free(part_ext->n_cur_interface_vtx);
    }


    if(part_ext->edge_edge_extended_idx != NULL) {
      PDM_free(part_ext->edge_edge_extended_idx);
      PDM_free(part_ext->edge_edge_extended);
      PDM_free(part_ext->edge_edge_interface);
      PDM_free(part_ext->edge_edge_path_itrf_idx);
      PDM_free(part_ext->edge_edge_path_itrf);
    }

    if(part_ext->vtx_vtx_extended_idx != NULL) {
      PDM_free(part_ext->vtx_vtx_extended_idx);
      PDM_free(part_ext->vtx_vtx_extended);
      PDM_free(part_ext->vtx_vtx_interface);
      PDM_free(part_ext->vtx_vtx_path_itrf_idx);
      PDM_free(part_ext->vtx_vtx_path_itrf);
    }

    /* Peu import l'ownership on free car on rend à l'utilisateur l'interface i_domain / i_part */
    if(part_ext->border_cell_face_idx != NULL) {
      PDM_free(part_ext->border_cell_face_idx);
      PDM_free(part_ext->border_cell_face);
    }

    if(part_ext->border_face_edge_idx != NULL) {
      PDM_free(part_ext->border_face_edge_idx);
      PDM_free(part_ext->border_face_edge);
    }

    if(part_ext->border_edge_vtx_idx != NULL) {
      PDM_free(part_ext->border_edge_vtx_idx);
      PDM_free(part_ext->border_edge_vtx);
    }

    if(part_ext->border_face_vtx_idx != NULL) {
      PDM_free(part_ext->border_face_vtx_idx);
      PDM_free(part_ext->border_face_vtx);
    }

    if(part_ext->border_face_group_idx != NULL) {
      PDM_free(part_ext->border_face_group_idx);
      PDM_free(part_ext->border_face_group);
    }

    if(part_ext->border_vtx != NULL) {
      PDM_free(part_ext->border_vtx);
    }

    if(part_ext->border_cell_ln_to_gn != NULL) {
      PDM_free(part_ext->border_cell_ln_to_gn);
    }

    if(part_ext->border_face_ln_to_gn != NULL) {
      PDM_free(part_ext->border_face_ln_to_gn);
    }

    if(part_ext->border_edge_ln_to_gn != NULL) {
      PDM_free(part_ext->border_edge_ln_to_gn);
    }

    if(part_ext->border_vtx_ln_to_gn != NULL) {
      PDM_free(part_ext->border_vtx_ln_to_gn);
    }

    if(part_ext->border_face_group_ln_to_gn != NULL) {
      PDM_free(part_ext->border_face_group_ln_to_gn);
    }


    for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {

      int shift_part = 0;
      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
          PDM_free(part_ext->cell_cell_extended_idx         [i_depth][i_part+shift_part]);
          PDM_free(part_ext->cell_cell_extended_n           [i_depth][i_part+shift_part]);
          PDM_free(part_ext->cell_cell_extended             [i_depth][i_part+shift_part]);
          PDM_free(part_ext->unique_order_cell_cell_extended[i_depth][i_part+shift_part]);
          PDM_free(part_ext->cell_cell_interface            [i_depth][i_part+shift_part]);
        }
        shift_part += part_ext->n_part[i_domain];
      }
      PDM_free(part_ext->cell_cell_extended_idx           [i_depth]);
      PDM_free(part_ext->cell_cell_extended_n             [i_depth]);
      PDM_free(part_ext->cell_cell_extended               [i_depth]);
      PDM_free(part_ext->unique_order_cell_cell_extended  [i_depth]);
      PDM_free(part_ext->cell_cell_interface              [i_depth]);
      PDM_free(part_ext->n_unique_order_cell_cell_extended[i_depth]);
    }

    PDM_free(part_ext->cell_cell_extended_idx);
    PDM_free(part_ext->cell_cell_extended_n);
    PDM_free(part_ext->cell_cell_extended);
    PDM_free(part_ext->cell_cell_path_itrf_idx);
    PDM_free(part_ext->cell_cell_path_itrf);
    // PDM_free(part_ext->cell_cell_extended);
    PDM_free(part_ext->unique_order_cell_cell_extended);
    PDM_free(part_ext->cell_cell_interface);
    PDM_free(part_ext->n_unique_order_cell_cell_extended);

    /* Only shortcut of user data */
    PDM_free(part_ext->entity_cell_idx    );
    PDM_free(part_ext->entity_cell        );
    PDM_free(part_ext->entity_cell_n);

    PDM_free(part_ext->dist_neighbor_cell_n   );
    PDM_free(part_ext->dist_neighbor_cell_idx );
    PDM_free(part_ext->dist_neighbor_cell_desc);
    PDM_free(part_ext->dist_neighbor_cell_interface);
    PDM_free(part_ext->unique_order_dist_neighbor_cell);
    PDM_free(part_ext->n_unique_order_dist_neighbor_cell);

    /* Allocated by distant neightbor */
    PDM_free(part_ext->entity_cell_opp_idx);
    PDM_free(part_ext->entity_cell_opp_n  );
    PDM_free(part_ext->entity_cell_opp    );

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        PDM_free(part_ext->cell_cell_extended_pruned_idx[i_part+shift_part]);
        PDM_free(part_ext->cell_cell_extended_pruned    [i_part+shift_part]);
        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
          PDM_free(part_ext->cell_cell_interface_pruned   [i_part+shift_part]);
        }
      }
      shift_part += part_ext->n_part[i_domain];
    }
    PDM_free(part_ext->cell_cell_extended_pruned_idx);
    PDM_free(part_ext->cell_cell_extended_pruned    );
    PDM_free(part_ext->cell_cell_interface_pruned   );
  }


  if(part_ext->pinit_entity_bound_to_pentity_bound_idx != NULL) {
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_idx      [i_part]);
      PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_triplet  [i_part]);
      PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_interface[i_part]);
    }
    PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_idx      );
    PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_triplet  );
    PDM_free(part_ext->pinit_entity_bound_to_pentity_bound_interface);
  }


  PDM_free(part_ext->n_part_idx);
  PDM_free(part_ext->n_part_g_idx);
  PDM_free(part_ext->n_part);
  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(part_ext->parts[i_domain]);
  }
  PDM_free(part_ext->parts);

  PDM_free(part_ext->shift_by_domain_cell);
  PDM_free(part_ext->shift_by_domain_face);
  PDM_free(part_ext->shift_by_domain_edge);
  PDM_free(part_ext->shift_by_domain_vtx);

  PDM_free(part_ext->shift_by_domain_edge_group);
  PDM_free(part_ext->shift_by_domain_face_group);

  if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
    PDM_free(part_ext->composed_interface_idx);
    PDM_free(part_ext->composed_interface);
    PDM_free(part_ext->composed_ln_to_gn_sorted);
  }


  /**
   * Free ownerships
   */
  for (int i_entity=0; i_entity<PDM_MESH_ENTITY_MAX; ++i_entity) {
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_free(part_ext->ownership_border_ln_to_gn       [i_entity][i_domain]);
      PDM_free(part_ext->ownership_border_ln_to_gn_ancstr[i_entity][i_domain]);
      PDM_free(part_ext->ownership_border_group          [i_entity][i_domain]);
      PDM_free(part_ext->ownership_border_graph          [i_entity][i_domain]);
      PDM_free(part_ext->ownership_border_path_itrf      [i_entity][i_domain]);
    }
    PDM_free(part_ext->ownership_border_ln_to_gn       [i_entity]);
    PDM_free(part_ext->ownership_border_ln_to_gn_ancstr[i_entity]);
    PDM_free(part_ext->ownership_border_group          [i_entity]);
    PDM_free(part_ext->ownership_border_graph          [i_entity]);
    PDM_free(part_ext->ownership_border_path_itrf      [i_entity]);
  }
  PDM_free(part_ext->ownership_border_ln_to_gn);
  PDM_free(part_ext->ownership_border_ln_to_gn_ancstr);
  PDM_free(part_ext->ownership_border_group);
  PDM_free(part_ext->ownership_border_graph);
  PDM_free(part_ext->ownership_border_path_itrf);

  for (int i_conn=0; i_conn<PDM_CONNECTIVITY_TYPE_MAX; ++i_conn) {
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      PDM_free(part_ext->ownership_border_connectivity[i_conn][i_domain]);
    }
    PDM_free(part_ext->ownership_border_connectivity[i_conn]);
  }
  PDM_free(part_ext->ownership_border_connectivity);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    PDM_free(part_ext->ownership_border_vtx_coord[i_domain]);
  }
  PDM_free(part_ext->ownership_border_vtx_coord);


  PDM_free(part_ext);
}


/**
 *
 * \brief Get extended connectivity
 *
 * \param [in]  part_ext            \p PDM_part_extension_t structure instance
 * \param [in]  i_domain            Domain identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  connectivity_type   Connectivity type
 * \param [out] connect_idx         Connectivity index
 * \param [out] connect             Connectivity
 *
 * \return Number of leading entities
 *
 */
int
PDM_part_extension_connectivity_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect_idx,
 int                     **connect
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(connectivity_type)
  {
    case PDM_CONNECTIVITY_TYPE_CELL_FACE:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *connect_idx = part_ext->border_cell_face_idx         [shift_part+i_part];
      *connect     = part_ext->border_cell_face             [shift_part+i_part];

      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "cell_face");

    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_edge_idx  [shift_part+i_part];
      *connect     = part_ext->border_face_edge      [shift_part+i_part];
      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "face_edge");
    }
    break;

    case PDM_CONNECTIVITY_TYPE_FACE_VTX:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *connect_idx = part_ext->border_face_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_face_vtx       [shift_part+i_part];
    }
    break;

    case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *connect_idx = part_ext->border_edge_vtx_idx   [shift_part+i_part];
      *connect     = part_ext->border_edge_vtx       [shift_part+i_part];
      // PDM_log_trace_connectivity_int(*connect_idx, *connect, n_entity, "edge_vtx");
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  part_ext->ownership_border_connectivity[connectivity_type][i_domain][i_part] = PDM_OWNERSHIP_USER;

  return n_entity;
}

int
PDM_part_extension_connectivity_get2
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_connectivity_type_t   connectivity_type,
  int                     **connect_idx,
  int                     **connect,
  PDM_ownership_t           ownership
)
{
  if(part_ext->compute_kind == 0) {
    return PDM_part_extension_connectivity_get(part_ext,
                                               i_domain,
                                               i_part,
                                               connectivity_type,
                                               connect_idx,
                                               connect);
  }

  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  if (part_ext->ownership_border_connectivity[connectivity_type][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(connectivity_type)
    {
      case PDM_CONNECTIVITY_TYPE_CELL_FACE:
      {
        n_entity     = part_ext->n_cell_border       [shift_part+i_part];
        *connect_idx = part_ext->border_cell_face_idx[shift_part+i_part];
        *connect     = part_ext->border_cell_face    [shift_part+i_part];
      }
      break;

      case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
      {
        n_entity     = part_ext->n_face_border       [shift_part+i_part];
        *connect_idx = part_ext->border_face_edge_idx[shift_part+i_part];
        *connect     = part_ext->border_face_edge    [shift_part+i_part];
      }
      break;

      case PDM_CONNECTIVITY_TYPE_FACE_VTX:
      {
        n_entity     = part_ext->n_face_border      [shift_part+i_part];
        *connect_idx = part_ext->border_face_vtx_idx[shift_part+i_part];
        *connect     = part_ext->border_face_vtx    [shift_part+i_part];
      }
      break;

      case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
      {
        n_entity     = part_ext->n_edge_border      [shift_part+i_part];
        *connect_idx = part_ext->border_edge_vtx_idx[shift_part+i_part];
        *connect     = part_ext->border_edge_vtx    [shift_part+i_part];
      }
      break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
      break;

    }
    part_ext->ownership_border_connectivity[connectivity_type][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_connectivity_get2 connectivity_type %d seems not to be computed.\n", connectivity_type);
  }

  return n_entity;
}


/**
 *
 * \brief Get global ids
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] ln_to_gn     Global ids (size = \ref n_elt)
 *
 * \return  n_elt  Number of elements
 *
 */
int
PDM_part_extension_ln_to_gn_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  PDM_g_num_t             **ln_to_gn
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      int n_cell   = part_ext->parts[i_domain][i_part].n_cell;
      n_entity     = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *ln_to_gn    = part_ext->border_cell_ln_to_gn         [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_cell :: ");
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face   = part_ext->parts[i_domain][i_part].n_face;
      n_entity     = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *ln_to_gn    = part_ext->border_face_ln_to_gn  [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_face :: ");
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge   = part_ext->parts[i_domain][i_part].n_edge;
      n_entity     = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *ln_to_gn    = part_ext->border_edge_ln_to_gn  [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      int n_vtx   = part_ext->parts[i_domain][i_part].n_vtx;
      n_entity     = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      *ln_to_gn    = part_ext->border_vtx_ln_to_gn  [shift_part+i_part];
      // PDM_log_trace_array_int(*ln_to_gn, n_entity, "extent_vtx :: ");
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;
  }

  part_ext->ownership_border_ln_to_gn[mesh_entity][i_domain][i_part] = PDM_OWNERSHIP_USER;

  return n_entity;
}
int
PDM_part_extension_ln_to_gn_get2
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  PDM_g_num_t             **ln_to_gn,
  PDM_ownership_t           ownership
)
{

  if(part_ext->compute_kind == 0) {
    return PDM_part_extension_ln_to_gn_get(part_ext,
                                            i_domain,
                                            i_part,
                                            mesh_entity,
                                            ln_to_gn);
  }

  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  if (part_ext->ownership_border_ln_to_gn[mesh_entity][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(mesh_entity)
    {
      case PDM_MESH_ENTITY_CELL:
      {
        n_entity     = part_ext->n_cell_border[shift_part+i_part];
        *ln_to_gn    = part_ext->border_cell_ln_to_gn[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_FACE:
      {
        n_entity     = part_ext->n_face_border[shift_part+i_part];
        *ln_to_gn    = part_ext->border_face_ln_to_gn[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_EDGE:
      {
        n_entity     = part_ext->n_edge_border[shift_part+i_part];
        *ln_to_gn    = part_ext->border_edge_ln_to_gn[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_VTX:
      {
        n_entity     = part_ext->n_vtx_border[shift_part+i_part];
        *ln_to_gn    = part_ext->border_vtx_ln_to_gn[shift_part+i_part];
      }
      break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown mesh_entity \n");
      break;
    }
    part_ext->ownership_border_ln_to_gn[mesh_entity][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_ln_to_gn_get2 mesh_entity %d seems not to be computed.\n", mesh_entity);
  }

  return n_entity;
}


int
PDM_part_extension_ancestor_ln_to_gn_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  PDM_g_num_t             **ancestor_ln_to_gn,
  PDM_ownership_t           ownership
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  if (part_ext->ownership_border_ln_to_gn_ancstr[mesh_entity][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(mesh_entity)
    {
      case PDM_MESH_ENTITY_CELL:
      {
        n_entity           = part_ext->n_cell_border              [shift_part+i_part];
        *ancestor_ln_to_gn = part_ext->border_cell_ln_to_gn_ancstr[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_FACE:
      {
        n_entity           = part_ext->n_face_border              [shift_part+i_part];
        *ancestor_ln_to_gn = part_ext->border_face_ln_to_gn_ancstr[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_EDGE:
      {
        n_entity           = part_ext->n_edge_border              [shift_part+i_part];
        *ancestor_ln_to_gn = part_ext->border_edge_ln_to_gn_ancstr[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_VTX:
      {
        n_entity           = part_ext->n_vtx_border              [shift_part+i_part];
        *ancestor_ln_to_gn = part_ext->border_vtx_ln_to_gn_ancstr[shift_part+i_part];
      }
      break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown mesh_entity \n");
      break;

    }
    part_ext->ownership_border_ln_to_gn_ancstr[mesh_entity][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_ancestor_ln_to_gn_get mesh_entity %d seems not to be computed.\n", mesh_entity);
  }

  return n_entity;
}


/**
 *
 * \brief Get global ids
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [out] ln_to_gn     Global ids (size = \ref n_elt)
 *
 * \return  n_elt  Number of elements
 *
 */

int
PDM_part_extension_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                     **interface_no
)
{
  int n_entity = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      int n_cell    = part_ext->parts[i_domain][i_part].n_cell;
      n_entity      = part_ext->cell_cell_extended_pruned_idx[shift_part+i_part][n_cell];
      *interface_no = part_ext->cell_cell_interface_pruned   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face    = part_ext->parts[i_domain][i_part].n_face;
      n_entity      = part_ext->face_face_extended_idx[shift_part+i_part][n_face];
      *interface_no = part_ext->face_face_interface   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge    = part_ext->parts[i_domain][i_part].n_edge;
      n_entity      = part_ext->edge_edge_extended_idx[shift_part+i_part][n_edge];
      *interface_no = part_ext->edge_edge_interface   [shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      int n_vtx   =  part_ext->parts[i_domain][i_part].n_vtx;
      n_entity      = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
      *interface_no = part_ext->vtx_vtx_interface   [shift_part+i_part];
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

  }

  return n_entity;
}


void
PDM_part_extension_graph_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  // int                     **pentity_to_entity_idx,
  int                     **pentity_to_entity,
  PDM_ownership_t           ownership
)
{
  int shift_part = part_ext->n_part_idx[i_domain];
  if (part_ext->ownership_border_graph[mesh_entity][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(mesh_entity)
    {
      case PDM_MESH_ENTITY_CELL:
      {
        // *pentity_to_entity_idx = part_ext->cell_cell_extended_idx2[shift_part+i_part];
        *pentity_to_entity     = part_ext->cell_cell_extended2    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_FACE:
      {
        // *pentity_to_entity_idx = part_ext->face_face_extended_idx[shift_part+i_part];
        *pentity_to_entity     = part_ext->face_face_extended    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_EDGE:
      {
        // *pentity_to_entity_idx = part_ext->edge_edge_extended_idx[shift_part+i_part];
        *pentity_to_entity     = part_ext->edge_edge_extended    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_VTX:
      {
        // *pentity_to_entity_idx = part_ext->vtx_vtx_extended_idx[shift_part+i_part];
        *pentity_to_entity     = part_ext->vtx_vtx_extended    [shift_part+i_part];
        break;
      }

      default :
      {
        PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
        break;
      }
    }
    part_ext->ownership_border_graph[mesh_entity][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_graph_get mesh_entity %d seems not to be computed.\n", mesh_entity);
  }
}


int
PDM_part_extension_path_interface_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                     **path_itrf_idx,
  int                     **path_itrf,
  PDM_ownership_t           ownership
)
{
  int shift_part = part_ext->n_part_idx[i_domain];
  int n_entity = 0;
  if (part_ext->ownership_border_path_itrf[mesh_entity][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(mesh_entity)
    {
      case PDM_MESH_ENTITY_CELL:
      {
        n_entity       = part_ext->n_cell_border[shift_part+i_part];
        *path_itrf_idx = part_ext->cell_cell_path_itrf_idx[shift_part+i_part];
        *path_itrf     = part_ext->cell_cell_path_itrf    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_FACE:
      {
        n_entity       = part_ext->n_face_border[shift_part+i_part];
        *path_itrf_idx = part_ext->face_face_path_itrf_idx[shift_part+i_part];
        *path_itrf     = part_ext->face_face_path_itrf    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_EDGE:
      {
        n_entity       = part_ext->n_edge_border[shift_part+i_part];
        *path_itrf_idx = part_ext->edge_edge_path_itrf_idx[shift_part+i_part];
        *path_itrf     = part_ext->edge_edge_path_itrf    [shift_part+i_part];
        break;
      }

      case PDM_MESH_ENTITY_VTX:
      {
        n_entity       = part_ext->n_vtx_border[shift_part+i_part];
        *path_itrf_idx = part_ext->vtx_vtx_path_itrf_idx[shift_part+i_part];
        *path_itrf     = part_ext->vtx_vtx_path_itrf    [shift_part+i_part];
        break;
      }

      default :
      {
        PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
        break;
      }
    }
    part_ext->ownership_border_path_itrf[mesh_entity][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_path_interface_get mesh_entity %d seems not to be computed.\n", mesh_entity);
  }
  return n_entity;
}


int
PDM_part_extension_group_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                     **group_entity_idx,
  int                     **group_entity,
  PDM_g_num_t             **group_entity_ln_to_gn
)
{
  int n_group = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  switch(mesh_entity)
  {
    case PDM_MESH_ENTITY_CELL:
    {
      PDM_error(__FILE__, __LINE__, 0, "PDM_MESH_ENTITY_CELL has no groups\n");
    }
    break;

    case PDM_MESH_ENTITY_FACE:
    {
      int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
      n_group                = n_face_group;
      *group_entity_idx      = part_ext->border_face_group_idx     [shift_part+i_part];
      *group_entity          = part_ext->border_face_group         [shift_part+i_part];
      *group_entity_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_EDGE:
    {
      int n_edge_group = part_ext->parts[i_domain][i_part].n_edge_group;
      n_group                = n_edge_group;
      *group_entity_idx      = part_ext->border_edge_group_idx     [shift_part+i_part];
      *group_entity          = part_ext->border_edge_group         [shift_part+i_part];
      *group_entity_ln_to_gn = part_ext->border_edge_group_ln_to_gn[shift_part+i_part];
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      PDM_error(__FILE__, __LINE__, 0, "PDM_MESH_ENTITY_VTX has no groups\n");
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown mesh_entity \n");
    break;

  }
  
  part_ext->ownership_border_group[mesh_entity][i_domain][i_part] = PDM_OWNERSHIP_USER;

  return n_group;
}

int
PDM_part_extension_group_get2
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                     **group_entity_idx,
  int                     **group_entity,
  PDM_g_num_t             **group_entity_ln_to_gn,
  PDM_ownership_t           ownership
)
{

  if(part_ext->compute_kind == 0) {
    return PDM_part_extension_group_get(part_ext,
                                            i_domain,
                                            i_part,
                                            mesh_entity,
                                            group_entity_idx,
                                            group_entity,
                                            group_entity_ln_to_gn);
  }

  int n_group = -1;
  int shift_part = part_ext->n_part_idx[i_domain];
  if (part_ext->ownership_border_group[mesh_entity][i_domain][i_part]!=PDM_OWNERSHIP_BAD_VALUE) {
    switch(mesh_entity)
    {
      case PDM_MESH_ENTITY_CELL:
      {
        PDM_error(__FILE__, __LINE__, 0, "PDM_MESH_ENTITY_CELL has no groups\n");
      }
      break;

      case PDM_MESH_ENTITY_FACE:
      {
        int n_face_group = part_ext->parts[i_domain][i_part].n_face_group;
        n_group                = n_face_group;
        *group_entity_idx      = part_ext->border_face_group_idx     [shift_part+i_part];
        *group_entity          = part_ext->border_face_group         [shift_part+i_part];
        *group_entity_ln_to_gn = part_ext->border_face_group_ln_to_gn[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_EDGE:
      {
        int n_edge_group = part_ext->parts[i_domain][i_part].n_edge_group;
        n_group                = n_edge_group;
        *group_entity_idx      = part_ext->border_edge_group_idx     [shift_part+i_part];
        *group_entity          = part_ext->border_edge_group         [shift_part+i_part];
        *group_entity_ln_to_gn = part_ext->border_edge_group_ln_to_gn[shift_part+i_part];
      }
      break;

      case PDM_MESH_ENTITY_VTX:
      {
        PDM_error(__FILE__, __LINE__, 0, "PDM_MESH_ENTITY_VTX has no groups\n");
      }
      break;

    default :
      PDM_error(__FILE__, __LINE__, 0, "Unknown mesh_entity \n");
      break;

    }
    part_ext->ownership_border_group[mesh_entity][i_domain][i_part] = ownership;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_extension_group_get2 mesh_entity %d seems not to be computed.\n", mesh_entity);
  }

  return n_group;
}


/**
 *
 * \brief Get vertex coordinates
 *
 * \param [in]  part_ext     Pointer to \ref PDM_part_extension_t object
 * \param [in]  i_domain     Id of current domain
 * \param [in]  i_part       Id of current partition
 * \param [out] vtx_coord    Vertex coordinates (size = \ref n_vtx * 3)
 *
 * \return  n_vtx  Number of vertices
 *
 */
int
PDM_part_extension_vtx_coord_get
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  double                  **vtx_coord
)
{
  int shift_part     = part_ext->n_part_idx[i_domain];
  int n_vtx          = part_ext->parts[i_domain][i_part].n_vtx;
  int n_vtx_extended = part_ext->vtx_vtx_extended_idx[shift_part+i_part][n_vtx];
  *vtx_coord = part_ext->border_vtx[shift_part+i_part];

  part_ext->ownership_border_vtx_coord[i_domain][i_part] = PDM_OWNERSHIP_USER;

  return n_vtx_extended;
}
int
PDM_part_extension_vtx_coord_get2
(
  PDM_part_extension_t     *part_ext,
  int                       i_domain,
  int                       i_part,
  double                  **vtx_coord,
  PDM_ownership_t           ownership
)
{
  if(part_ext-> compute_kind == 0) {
    return PDM_part_extension_vtx_coord_get(part_ext,
                                            i_domain,
                                            i_part,
                                            vtx_coord);
  }

  int shift_part     = part_ext->n_part_idx[i_domain];
  int n_vtx_extended = part_ext->n_vtx_border[shift_part+i_part];
  *vtx_coord = part_ext->border_vtx[shift_part+i_part];
  
  part_ext->ownership_border_vtx_coord[i_domain][i_part] = ownership;

  return n_vtx_extended;
}


int
PDM_part_extension_composed_interface_get
(
 PDM_part_extension_t     *part_ext,
 int                     **composed_interface_idx,
 int                     **composed_interface,
 PDM_g_num_t             **composed_ln_to_gn_sorted
)
{
  *composed_interface_idx   = part_ext->composed_interface_idx;
  *composed_interface       = part_ext->composed_interface;
  *composed_ln_to_gn_sorted = part_ext->composed_ln_to_gn_sorted;

  return part_ext->n_composed_interface;
}


/**
 *
 * \brief Create part_to_part from interior and extended elements
 *
 * \param [out]  ptp                             Part to part structure
 * \param [in]   n_part                          Number of partitions
 * \param [in]   n_int_cell                      Number of interior elements
 * \param [in]   int_cell_ln_to_gn               gnum of interior elements
 * \param [in]   n_ext_cell                      Number of extended elements
 * \param [in]   ext_cell_ln_to_gn               gnum of extended elements
 * \param [out]  n_selected_cell_to_send         Number of elements selected for send
 * \param [out]  selected_cell_to_send           Local numbering of elements selected for send
 *
 */

/* TODO : to generalize for SONICS for vertices
 */

void
PDM_part_to_part_create_from_extension
(
       PDM_part_to_part_t **ptp,
 const int                  n_part,
       int                 *n_int_cell,
 const PDM_g_num_t        **int_cell_ln_to_gn,
       int                 *n_ext_cell,
 const PDM_g_num_t        **ext_cell_ln_to_gn,
       int                **n_selected_cell_to_send,
       int               ***selected_cell_to_send,
 const PDM_MPI_Comm         comm
)
{

  PDM_g_num_t _max_g_num = -1;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i = 0; i < n_ext_cell[i_part]; ++i) {
      _max_g_num = PDM_MAX(_max_g_num, ext_cell_ln_to_gn[i_part][i]);
    }
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      _max_g_num = PDM_MAX(_max_g_num, int_cell_ln_to_gn[i_part][i]);
    }
  }

  PDM_g_num_t max_g_num = 0;
  PDM_MPI_Allreduce(&_max_g_num, &max_g_num, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_g_num_t* distrib_cell = PDM_compute_uniform_entity_distribution(comm, max_g_num);

  PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_MERGE,
                                                                   1.,
                                             (PDM_g_num_t **)      ext_cell_ln_to_gn,
                                                                   distrib_cell,
                                                                   n_ext_cell,
                                                                   n_part,
                                                                   comm);

  int  block_n_elt = PDM_part_to_block_n_elt_block_get(ptb);
  int *block_n = NULL;
  PDM_malloc(block_n, block_n_elt, int);
  for(int i = 0; i < block_n_elt; ++i) {
    block_n[i] = 1;
  }

  PDM_g_num_t* distrib_adapt = PDM_part_to_block_adapt_partial_block_to_block(ptb, &block_n, max_g_num);
  PDM_free(distrib_adapt);

  PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_cell,
                              (const PDM_g_num_t **)  int_cell_ln_to_gn,
                                                      n_int_cell,
                                                      n_part,
                                                      comm);

  int   stride_one = 1;
  int **is_ext_cell = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         block_n,
                         NULL,
           (void ***)    &is_ext_cell);

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);
  PDM_free(block_n);
  PDM_free(distrib_cell);

  int          *_n_selected_cell_to_send        = NULL;
  int         **_selected_cell_to_send_idx      = NULL;
  int         **_selected_cell_to_send          = NULL;
  PDM_g_num_t **_selected_cell_to_send_ln_to_gn = NULL;
  PDM_malloc(_n_selected_cell_to_send       , n_part, int          );
  PDM_malloc(_selected_cell_to_send_idx     , n_part, int         *);
  PDM_malloc(_selected_cell_to_send         , n_part, int         *);
  PDM_malloc(_selected_cell_to_send_ln_to_gn, n_part, PDM_g_num_t *);

  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Compute buffer size
    int n_ext_to_send = 0;
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      n_ext_to_send += is_ext_cell[i_part][i];
    }

    _n_selected_cell_to_send       [i_part] = n_ext_to_send;
    PDM_malloc(_selected_cell_to_send_idx     [i_part], n_ext_to_send+1, int        );
    PDM_malloc(_selected_cell_to_send         [i_part], n_ext_to_send  , int        );
    PDM_malloc(_selected_cell_to_send_ln_to_gn[i_part], n_ext_to_send  , PDM_g_num_t);

    int idx_write = 0;
    _selected_cell_to_send_idx[i_part][0] = 0;
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      if(is_ext_cell[i_part][i] == 1) {
        _selected_cell_to_send_idx     [i_part][idx_write+1] = _selected_cell_to_send_idx[i_part][idx_write] + is_ext_cell[i_part][i];
        _selected_cell_to_send         [i_part][idx_write]   = i+1;
        _selected_cell_to_send_ln_to_gn[i_part][idx_write]   = int_cell_ln_to_gn[i_part][i];
        idx_write++;
      }
    }

  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(is_ext_cell[i_part]);
  }
  PDM_free(is_ext_cell);

  PDM_part_to_part_t  *_ptp = PDM_part_to_part_create((const PDM_g_num_t **) _selected_cell_to_send_ln_to_gn,
                                                                             _n_selected_cell_to_send,
                                                                             n_part,
                                                      (const PDM_g_num_t **) ext_cell_ln_to_gn,
                                                                             n_ext_cell,
                                                                             n_part,
                                                              (const int **) _selected_cell_to_send_idx,
                                                      (const PDM_g_num_t **) _selected_cell_to_send_ln_to_gn,
                                                                             comm);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(_selected_cell_to_send_idx     [i_part]);
    PDM_free(_selected_cell_to_send_ln_to_gn[i_part]);
  }
  PDM_free(_selected_cell_to_send_idx);
  PDM_free(_selected_cell_to_send_ln_to_gn);

  *n_selected_cell_to_send = _n_selected_cell_to_send;
  *selected_cell_to_send   = _selected_cell_to_send;
  *ptp = _ptp;
}




/**
 *
 * \brief Set connectivity
 *
 * \param [in]  part_ext           \p PDM_part_extension_t structure instance
 * \param [in]  i_domain           Domain identifier
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  connect            Connectivity
 *
 */

void
PDM_part_extension_connectivity_set
(
 PDM_part_extension_t    *part_ext,
 int                      i_domain,
 int                      i_part,
 PDM_connectivity_type_t  connectivity_type,
 int                     *connect_idx,
 int                     *connect
 )
{
  part_ext->has_connectivity[connectivity_type] = PDM_TRUE;

  switch (connectivity_type) {
    case PDM_CONNECTIVITY_TYPE_EDGE_VTX: {
      part_ext->parts[i_domain][i_part].edge_vtx = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_VTX: {
      part_ext->parts[i_domain][i_part].face_vtx_idx = connect_idx;
      part_ext->parts[i_domain][i_part].face_vtx     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_EDGE: {
      part_ext->parts[i_domain][i_part].face_edge_idx = connect_idx;
      part_ext->parts[i_domain][i_part].face_edge     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_CELL_FACE: {
      part_ext->parts[i_domain][i_part].cell_face_idx = connect_idx;
      part_ext->parts[i_domain][i_part].cell_face     = connect;
      break;
    }

    case PDM_CONNECTIVITY_TYPE_FACE_CELL: {
      part_ext->parts[i_domain][i_part].face_cell = connect;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Connectivity type %d not yet supported\n",
                connectivity_type);
      break;
    }

  }
}


/**
 *
 * \brief Set global ids
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [in]  n_entity     Local number of entities
 * \param [in]  ln_to_gn     Global ids (size = \p n_entity)
 *
 */

void
PDM_part_extension_ln_to_gn_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       mesh_entity,
 int                       n_entity,
 PDM_g_num_t              *ln_to_gn
)
{
  switch (mesh_entity) {
    case PDM_MESH_ENTITY_VTX: {
      part_ext->parts[i_domain][i_part].n_vtx         = n_entity;
      part_ext->parts[i_domain][i_part].vtx_ln_to_gn  = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_EDGE: {
      part_ext->parts[i_domain][i_part].n_edge        = n_entity;
      part_ext->parts[i_domain][i_part].edge_ln_to_gn = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].n_face        = n_entity;
      part_ext->parts[i_domain][i_part].face_ln_to_gn = ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_CELL: {
      part_ext->parts[i_domain][i_part].n_cell        = n_entity;
      part_ext->parts[i_domain][i_part].cell_ln_to_gn = ln_to_gn;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Invalid entity type %d\n", mesh_entity);
      break;
    }

  }
}


/**
 *
 * \brief Set vertex coordinates
 *
 * \param [in]  part_ext     \p PDM_part_extension_t structure instance
 * \param [in]  i_domain     Domain identifier
 * \param [in]  i_part       Partition identifier
 * \param [in]  vtx_coord    Vertex coordinates (size = 3 * *n_vtx*)
 *
 */

void
PDM_part_extension_vtx_coord_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 double                   *vtx_coord
)
{
  part_ext->parts[i_domain][i_part].vtx = vtx_coord;
}


/**
 *
 * \brief Set the connection graph between partitions for the requested entity type
 *
 * \param [in]  multipart             \p PDM_part_extension_t structure instance
 * \param [in]  i_domain              Domain identifier
 * \param [in]  i_part                Partition identifier
 * \param [in]  entity_type           Type of mesh entity
 * \param [in]  part_bound_proc_idx   Partitioning boundary entities index from process (size = *n_rank* + 1)
 * \param [in]  part_bound_part_idx   Partitioning boundary entities index from partition (size = *n_total_part* + 1)
 * \param [in]  part_bound            Partitioning boundary entities (size = 4 * *n_entity_part_bound* = \p part_bound_proc_idx[*n_rank])
 */

void
PDM_part_extension_part_bound_graph_set
(
 PDM_part_extension_t *part_ext,
 int                   i_domain,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 int                  *part_bound_proc_idx,
 int                  *part_bound_part_idx,
 int                  *part_bound
)
{
  switch (entity_type) {

    case PDM_MESH_ENTITY_VTX: {
      part_ext->parts[i_domain][i_part].vtx_part_bound_proc_idx  = part_bound_proc_idx;
      part_ext->parts[i_domain][i_part].vtx_part_bound_part_idx  = part_bound_part_idx;
      part_ext->parts[i_domain][i_part].vtx_part_bound           = part_bound;
      if (part_bound_proc_idx!=NULL) {
        part_ext->user_defined_bnd_graph = 1;
      }
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = part_bound_proc_idx;
      part_ext->parts[i_domain][i_part].face_part_bound_part_idx = part_bound_part_idx;
      part_ext->parts[i_domain][i_part].face_part_bound          = part_bound;
      if (part_bound_proc_idx!=NULL) {
        part_ext->user_defined_bnd_graph = 1;
      }
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Entity type %d not yet supported\n",
                entity_type);
      break;
    }

  }
}


/**
 *
 * \brief Set group description
 *
 * \param [in]  part_ext               \p PDM_part_extension_t structure instance
 * \param [in]  i_domain               Domain identifier
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 * \param [in]  group_entity_idx       Index for group->entity connectivity (size = \p n_group)
 * \param [in]  group_entity           Group->entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  group_entity_ln_to_gn  Group->entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 *
 */

void
PDM_part_extension_group_set
(
 PDM_part_extension_t     *part_ext,
 int                       i_domain,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       n_group,
 int                      *group_entity_idx,
 int                      *group_entity,
 PDM_g_num_t              *group_entity_ln_to_gn
)
{
  switch (entity_type) {

    case PDM_MESH_ENTITY_EDGE: {
      part_ext->parts[i_domain][i_part].n_edge_group        = n_group;
      part_ext->parts[i_domain][i_part].edge_bound_idx      = group_entity_idx;
      part_ext->parts[i_domain][i_part].edge_bound          = group_entity;
      part_ext->parts[i_domain][i_part].edge_bound_ln_to_gn = group_entity_ln_to_gn;
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].n_face_group        = n_group;
      part_ext->parts[i_domain][i_part].face_bound_idx      = group_entity_idx;
      part_ext->parts[i_domain][i_part].face_bound          = group_entity;
      part_ext->parts[i_domain][i_part].face_bound_ln_to_gn = group_entity_ln_to_gn;
      break;
    }

    default: {
      PDM_error(__FILE__, __LINE__, 0, "Entity type %d not yet supported\n",
                entity_type);
      break;
    }
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

