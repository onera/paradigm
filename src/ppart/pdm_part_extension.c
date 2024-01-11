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

static
void
_compute_offset
(
  PDM_part_extension_t *part_ext
)
{
  int **pn_cell       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_edge       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_vtx        = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face_group = malloc(part_ext->n_domain * sizeof(int *));

  PDM_g_num_t ***cell_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***edge_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***vtx_ln_to_gn        = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_group_ln_to_gn = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    pn_cell      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_edge      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_vtx       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face_group[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));

    cell_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    edge_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    vtx_ln_to_gn       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_group_ln_to_gn[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      pn_cell            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_cell;
      pn_face            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pn_edge            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_vtx             [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_face_group      [i_domain][i_part] = 0;
      if (part_ext->parts[i_domain][i_part].n_face_group > 0) {
        pn_face_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_idx[part_ext->parts[i_domain][i_part].n_face_group];
      }

      cell_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
      face_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      edge_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      vtx_ln_to_gn       [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      face_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_group_ln_to_gn;

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

  assert(part_ext->shift_by_domain_face_group == NULL);
  part_ext->shift_by_domain_face_group = PDM_compute_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                                                               part_ext->n_part,
                                                                               pn_face_group,
                                                                               face_group_ln_to_gn,
                                                                               part_ext->comm);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(pn_cell      [i_domain]);
    free(pn_face      [i_domain]);
    free(pn_edge      [i_domain]);
    free(pn_vtx       [i_domain]);
    free(pn_face_group[i_domain]);

    free(cell_ln_to_gn      [i_domain]);
    free(face_ln_to_gn      [i_domain]);
    free(edge_ln_to_gn      [i_domain]);
    free(vtx_ln_to_gn       [i_domain]);
    free(face_group_ln_to_gn[i_domain]);
  }

  free(pn_cell      );
  free(pn_face      );
  free(pn_edge      );
  free(pn_vtx       );
  free(pn_face_group);

  free(cell_ln_to_gn      );
  free(face_ln_to_gn      );
  free(edge_ln_to_gn      );
  free(vtx_ln_to_gn       );
  free(face_group_ln_to_gn);
}



static
void
_offset_parts_by_domain
(
  PDM_part_extension_t *part_ext,
  int                   sens
)
{
  int **pn_cell       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_edge       = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_vtx        = malloc(part_ext->n_domain * sizeof(int *));
  int **pn_face_group = malloc(part_ext->n_domain * sizeof(int *));

  PDM_g_num_t ***cell_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***edge_ln_to_gn       = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***vtx_ln_to_gn        = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***face_group_ln_to_gn = malloc(part_ext->n_domain * sizeof(PDM_g_num_t **));

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    pn_cell      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_edge      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_vtx       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));
    pn_face_group[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(int *));

    cell_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    edge_ln_to_gn      [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    vtx_ln_to_gn       [i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));
    face_group_ln_to_gn[i_domain] = malloc(part_ext->n_part[i_domain] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {

      pn_cell            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_cell;
      pn_face            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pn_edge            [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_vtx             [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_face_group      [i_domain][i_part] = 0;
      if (part_ext->parts[i_domain][i_part].n_face_group > 0) {
        pn_face_group      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_bound_idx[part_ext->parts[i_domain][i_part].n_face_group];
      }

      cell_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].cell_ln_to_gn;
      face_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      edge_ln_to_gn      [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      vtx_ln_to_gn       [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      face_group_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_group_ln_to_gn;

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

  PDM_offset_ln_to_gn_by_domain(part_ext->n_domain,
                                part_ext->n_part,
                                pn_face_group,
                                face_group_ln_to_gn,
                                part_ext->shift_by_domain_face_group,
                                sens);

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(pn_cell      [i_domain]);
    free(pn_face      [i_domain]);
    free(pn_edge      [i_domain]);
    free(pn_vtx       [i_domain]);
    free(pn_face_group[i_domain]);

    free(cell_ln_to_gn      [i_domain]);
    free(face_ln_to_gn      [i_domain]);
    free(edge_ln_to_gn      [i_domain]);
    free(vtx_ln_to_gn       [i_domain]);
    free(face_group_ln_to_gn[i_domain]);
  }

  free(pn_cell      );
  free(pn_face      );
  free(pn_edge      );
  free(pn_vtx       );
  free(pn_face_group);

  free(cell_ln_to_gn      );
  free(face_ln_to_gn      );
  free(edge_ln_to_gn      );
  free(vtx_ln_to_gn       );
  free(face_group_ln_to_gn);
}

static
void
_compute_other_part_domain_interface
(
 PDM_part_extension_t *part_ext
)
{
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

  // En gros noeud centré avec toutes les connectivités
  if(is_describe_vtx == 1 &&
     (is_describe_edge == 0 || is_describe_face == 0) &&
     (have_edge       == 1 && have_face == 1)) {

    // Rebuild domaine_interface in distributed frame
    // PDM_domain_interface* dintrf = PDM_part_domain_interface_to_domain_interface()


    int          **pn_vtx         = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    int          **pn_edge        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    int          **pn_face        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
    PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
    int         ***pedge_vtx_idx  = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pedge_vtx      = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pface_edge_idx = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));
    int         ***pface_edge     = (int         ***) malloc( part_ext->n_domain * sizeof(int         **));

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      pn_vtx        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pn_edge       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pn_face       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pvtx_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pedge_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pface_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      pedge_vtx_idx [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pedge_vtx     [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pface_edge_idx[i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
      pface_edge    [i_domain] = (int         **) malloc( part_ext->n_domain * sizeof(int         *));

      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_vtx        [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        pn_edge       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
        pn_face       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        pvtx_ln_to_gn [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
        pedge_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
        pface_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
        pface_edge_idx[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge_idx;
        pface_edge    [i_domain][i_part] = part_ext->parts[i_domain][i_part].face_edge;
        // pedge_vtx     [i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_vtx;
        pedge_vtx     [i_domain][i_part] = (int *) malloc((2 * pn_edge[i_domain][i_part]) * sizeof(int));


        int *_edge_vtx = part_ext->parts[i_domain][i_part].edge_vtx;
        for(int i_edge = 0; i_edge < pn_edge[i_domain][i_part]; ++i_edge) {
          pedge_vtx     [i_domain][i_part][2*i_edge  ] =  _edge_vtx[2*i_edge  ];
          pedge_vtx     [i_domain][i_part][2*i_edge+1] = -_edge_vtx[2*i_edge+1];
          // printf("i_edge = %i (%i)- i_vtx1 = %i | i_vtx2 = %i \n", i_edge, (int) pedge_ln_to_gn[i_domain][i_part][i_edge], _edge_vtx[2*i_edge  ], -_edge_vtx[2*i_edge+1]);
        }

        int _nedge = pn_edge[i_domain][i_part];
        pedge_vtx_idx [i_domain][i_part] = malloc((_nedge+1) * sizeof(int));

        pedge_vtx_idx [i_domain][i_part][0] = 0;
        for(int i = 0; i < _nedge; ++i) {
          pedge_vtx_idx [i_domain][i_part][i+1] = pedge_vtx_idx [i_domain][i_part][i] + 2;
        }
      }
    }

    if(is_describe_edge == 0) {

      // Translate
      printf("Translate vtx to edge ... \n");
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
      printf("Translate vtx to edge END \n");
    }


    if(have_face == 1 &&  is_describe_face == 0) {

      // Translate
      printf("Translate edge to face ... \n");
      // Translate
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
    }


    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(pedge_vtx_idx [i_domain][i_part]);
        free(pedge_vtx     [i_domain][i_part]);
      }
      free(pn_vtx        [i_domain]);
      free(pn_edge       [i_domain]);
      free(pn_face       [i_domain]);
      free(pvtx_ln_to_gn [i_domain]);
      free(pedge_ln_to_gn[i_domain]);
      free(pface_ln_to_gn[i_domain]);
      free(pedge_vtx_idx [i_domain]);
      free(pedge_vtx     [i_domain]);
      free(pface_edge_idx[i_domain]);
      free(pface_edge    [i_domain]);
    }
    free(pn_vtx        );
    free(pn_edge       );
    free(pn_face       );
    free(pvtx_ln_to_gn );
    free(pedge_ln_to_gn);
    free(pface_ln_to_gn);
    free(pedge_vtx_idx );
    free(pedge_vtx     );
    free(pface_edge_idx );
    free(pface_edge     );

  } else if (is_describe_face == 1) {

    // assert(is_describe_vtx == 0);

    // Faire la méthode de Julien face2vtx

  }
}

static
void
_setup_domain_interface_in_block_frame
(
 PDM_part_extension_t *part_ext
)
{
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

  int          **pn_vtx         = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
  int          **pn_edge        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
  int          **pn_face        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    pn_vtx        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
    pn_edge       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
    pn_face       [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
    pvtx_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
    pface_ln_to_gn[i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      pn_vtx        [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
      pn_face       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
      pvtx_ln_to_gn [i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      pedge_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      pface_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
    }
  }


  part_ext->dom_itrf = NULL;
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
      free(is_face_on_itrf[i_part]);
    }
    free(is_face_on_itrf);
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
      free(is_edge_on_itrf[i_part]);
    }
    free(is_edge_on_itrf);
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
      free(is_vtx_on_itrf[i_part]);
    }
    free(is_vtx_on_itrf);
  }

  /*
   * Free all
   */
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(pn_vtx        [i_domain]);
    free(pn_edge       [i_domain]);
    free(pn_face       [i_domain]);
    free(pvtx_ln_to_gn [i_domain]);
    free(pedge_ln_to_gn[i_domain]);
    free(pface_ln_to_gn[i_domain]);
  }
  free(pn_vtx        );
  free(pn_edge       );
  free(pn_face       );
  free(pvtx_ln_to_gn );
  free(pedge_ln_to_gn);
  free(pface_ln_to_gn);

  /*
   * At this stage we have for all interface the information for all entities
   *   During the part_extension problem we need this information to merge / amend all extend entity
   *   We make a distributed vision that more convenient to do all this kind of operation
   */
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    part_ext->ptb_itrf[i] = NULL; // (PDM_part_to_block_t **) malloc( n_interface * sizeof(PDM_part_to_block_t **));
    part_ext->opp_gnum[i] = NULL; // (PDM_g_num_t         **) malloc( n_interface * sizeof(PDM_g_num_t         **));

    // for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
    //   part_ext->ptb_itrf[i][i_itrf] = NULL;
    //   part_ext->opp_gnum[i][i_itrf] = NULL;
    // }
  }

  part_ext->n_interface = n_interface;

  PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
                                      PDM_BOUND_TYPE_VTX,
                                      part_ext->shift_by_domain_vtx,
                                      &part_ext->ptb_itrf[PDM_BOUND_TYPE_VTX],
                                      &part_ext->opp_gnum[PDM_BOUND_TYPE_VTX]);


  if(is_describe_edge) {
    PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
                                        PDM_BOUND_TYPE_EDGE,
                                        part_ext->shift_by_domain_edge,
                                       &part_ext->ptb_itrf[PDM_BOUND_TYPE_EDGE],
                                       &part_ext->opp_gnum[PDM_BOUND_TYPE_EDGE]);
  }


  if(is_describe_face) {
    PDM_domain_interface_make_flat_view(part_ext->dom_itrf,
                                        PDM_BOUND_TYPE_FACE,
                                        part_ext->shift_by_domain_face,
                                       &part_ext->ptb_itrf[PDM_BOUND_TYPE_FACE],
                                       &part_ext->opp_gnum[PDM_BOUND_TYPE_FACE]);
  }


}


static
void
_part_extension_3d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);
}

static
void
_part_extension_2d
(
 PDM_part_extension_t *part_ext
)
{

  /*
   * Dans tous les cas on cherche a obtenir le graphe entre le propagateur et l'entité principale :
   *    - PDM_EXTEND_FROM_VTX
   *    - PDM_EXTEND_FROM_EDGE
   */
  int          **pn_entity1        = (int         ** ) malloc( part_ext->n_domain * sizeof(int          *));
  PDM_g_num_t ***pentity1_ln_to_gn = (PDM_g_num_t ***) malloc( part_ext->n_domain * sizeof(PDM_g_num_t **));
  PDM_bound_type_t bound_type = PDM_BOUND_TYPE_MAX;
  if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
    bound_type = PDM_BOUND_TYPE_VTX;

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      pn_entity1        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pentity1_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_vtx;
        pentity1_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      }
    }

  } else if(part_ext->extend_type == PDM_EXTEND_FROM_EDGE) {
    bound_type = PDM_BOUND_TYPE_EDGE;

    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      pn_entity1        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pentity1_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_edge;
        pentity1_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
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
  int **pentity1_extented_to_pentity1_idx       = NULL;
  int **pentity1_extented_to_pentity1_triplet   = NULL;
  int **pentity1_extented_to_pentity1_interface = NULL;
  PDM_part_extension_build_entity1_graph(part_ext->pdi,
                                         bound_type,
                                         part_ext->n_domain,
                                         part_ext->n_part,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         NULL,
                                         &pentity1_extented_to_pentity1_idx,
                                         &pentity1_extented_to_pentity1_triplet,
                                         &pentity1_extented_to_pentity1_interface,
                                         part_ext->comm);


  /*
   * 2 possibilities :
   *   - With face_vtx
   *   - With face_edge + edge_vtx
   */
  if(part_ext->have_edge == 1) {
    printf("Part extension with edges not implemented yet !");
    abort();
  } else {

  }


  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    free(pentity1_extented_to_pentity1_idx      [i_part]);
    free(pentity1_extented_to_pentity1_triplet  [i_part]);
    free(pentity1_extented_to_pentity1_interface[i_part]);
  }
  free(pentity1_extented_to_pentity1_idx      );
  free(pentity1_extented_to_pentity1_triplet  );
  free(pentity1_extented_to_pentity1_interface);


  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(pn_entity1       [i_domain]);
    free(pentity1_ln_to_gn[i_domain]);
  }
  free(pn_entity1);
  free(pentity1_ln_to_gn);

}

static
void
_part_extension_1d
(
 PDM_part_extension_t *part_ext
)
{
  PDM_UNUSED(part_ext);

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
  part_ext->compute_kind  = 1;

  part_ext->ln_part_tot = 0;
  for(int i_dom = 0; i_dom < part_ext->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < part_ext->n_part[i_dom]; ++i_part) {
      part_ext->ln_part_tot += 1;
    }
  }

  _compute_offset(part_ext);

  /*
   * Warm up all part_domain_interface
   *   - Identify all domain interface we have and make missing one
   *   - Go to distributed frame
   *   - Create the extend_type part_domain_interface to setup the first step of the algorithme
   */
  printf("_compute_other_part_domain_interface \n");
  _compute_other_part_domain_interface(part_ext);
  printf("_compute_other_part_domain_interface end\n");

  /*
   * Let's do an block frame domain_interface
   */
  printf("_setup_domain_interface_in_block_frame \n");
  _setup_domain_interface_in_block_frame(part_ext);
  printf("_setup_domain_interface_in_block_frame end\n");


  /* Manage shift */
  printf("_offset_parts_by_domain \n");
  _offset_parts_by_domain(part_ext, 1);
  printf("_offset_parts_by_domain end \n");

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
    if(part_ext->ptb_itrf[i] != NULL) {
      for(int i_itrf = 0; i_itrf < part_ext->n_interface; ++i_itrf) {
        if(part_ext->ptb_itrf[i][i_itrf] != NULL) {
          PDM_part_to_block_free(part_ext->ptb_itrf[i][i_itrf]);
          free(part_ext->opp_gnum[i][i_itrf]);
        }
      }
    }
    free(part_ext->ptb_itrf[i]);
    free(part_ext->opp_gnum[i]);
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
    free(part_ext->n_tot_part_by_domain);
    part_ext->n_tot_part_by_domain = NULL;
  }

  if(part_ext->neighbor_idx != NULL) {
    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->neighbor_idx       [i_part+shift_part]);
        free(part_ext->neighbor_desc      [i_part+shift_part]);
        free(part_ext->neighbor_interface [i_part+shift_part]);

        free(part_ext->dist_neighbor_cell_n           [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_idx         [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_desc        [i_part+shift_part]);
        free(part_ext->dist_neighbor_cell_interface   [i_part+shift_part]);
        free(part_ext->unique_order_dist_neighbor_cell[i_part+shift_part]);

        free(part_ext->entity_cell_opp_idx[i_part+shift_part]);
        free(part_ext->entity_cell_opp_n  [i_part+shift_part]);
        free(part_ext->entity_cell_opp    [i_part+shift_part]);

        free(part_ext->border_cell_list    [i_part+shift_part]);

        // if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
        free(part_ext->entity_cell_idx[i_part+shift_part]);
        free(part_ext->entity_cell_n  [i_part+shift_part]);
        free(part_ext->entity_cell    [i_part+shift_part]);
        // }


        if(part_ext->opp_interface_and_gnum_face[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_face[i_part+shift_part]);
        }

        if(part_ext->opp_interface_and_gnum_edge[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_edge[i_part+shift_part]);
        }

        if(part_ext->opp_interface_and_gnum_vtx[i_part+shift_part] != NULL) {
          free(part_ext->opp_interface_and_gnum_vtx[i_part+shift_part]);
        }

        if(part_ext->cur_interface_face[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_face[i_part+shift_part]);
        }

        if(part_ext->cur_interface_edge[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_edge[i_part+shift_part]);
        }

        if(part_ext->cur_interface_vtx[i_part+shift_part] != NULL) {
          free(part_ext->cur_interface_vtx[i_part+shift_part]);
        }

        if(part_ext->cur_sens_face[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_face[i_part+shift_part]);
        }

        if(part_ext->cur_sens_edge[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_edge[i_part+shift_part]);
        }

        if(part_ext->cur_sens_vtx[i_part+shift_part] != NULL) {
          free(part_ext->cur_sens_vtx[i_part+shift_part]);
        }


        if(part_ext->face_face_extended_idx != NULL) {
          free(part_ext->face_face_extended_idx[i_part+shift_part]);
        }

        if(part_ext->face_face_extended != NULL) {
          free(part_ext->face_face_extended[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended_idx != NULL) {
          free(part_ext->edge_edge_extended_idx[i_part+shift_part]);
        }

        if(part_ext->edge_edge_extended != NULL) {
          free(part_ext->edge_edge_extended[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended_idx != NULL) {
          free(part_ext->vtx_vtx_extended_idx[i_part+shift_part]);
        }

        if(part_ext->vtx_vtx_extended != NULL) {
          free(part_ext->vtx_vtx_extended[i_part+shift_part]);
        }

        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {

          if(part_ext->border_cell_face_idx != NULL) {
            free(part_ext->border_cell_face_idx[i_part+shift_part]);
            free(part_ext->border_cell_face    [i_part+shift_part]);
          }

          if(part_ext->border_face_edge_idx != NULL) {
            free(part_ext->border_face_edge_idx [i_part+shift_part]);
            free(part_ext->border_face_edge     [i_part+shift_part]);
          }

          if(part_ext->border_edge_vtx_idx != NULL) {
            free(part_ext->border_edge_vtx_idx[i_part+shift_part]);
            free(part_ext->border_edge_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_vtx_idx != NULL) {
            free(part_ext->border_face_vtx_idx[i_part+shift_part]);
            free(part_ext->border_face_vtx    [i_part+shift_part]);
          }

          if(part_ext->border_face_group_idx != NULL) {
            free(part_ext->border_face_group_idx[i_part+shift_part]);
            free(part_ext->border_face_group    [i_part+shift_part]);
          }

          if(part_ext->border_vtx != NULL) {
            free(part_ext->border_vtx[i_part+shift_part]);
          }

          if(part_ext->border_cell_ln_to_gn != NULL) {
            free(part_ext->border_cell_ln_to_gn[i_part+shift_part]);
          }

          if(part_ext->border_face_ln_to_gn != NULL) {
            free(part_ext->border_face_ln_to_gn[i_part+shift_part]);
            free(part_ext->face_face_interface  [i_part+shift_part]);
          }

          if(part_ext->border_edge_ln_to_gn != NULL) {
            free(part_ext->border_edge_ln_to_gn[i_part+shift_part]);
            free(part_ext->edge_edge_interface[i_part+shift_part]);
          }

          if(part_ext->border_vtx_ln_to_gn != NULL) {
            free(part_ext->border_vtx_ln_to_gn[i_part+shift_part]);
            free(part_ext->vtx_vtx_interface  [i_part+shift_part]);
          }

          if(part_ext->border_face_group_ln_to_gn != NULL) {
            free(part_ext->border_face_group_ln_to_gn[i_part+shift_part]);
          }

        }

        free(part_ext->cell_cell_idx[i_part+shift_part]);
        free(part_ext->cell_cell    [i_part+shift_part]);

      }

      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->n_entity_bound);
  }
  free(part_ext->neighbor_idx);
  free(part_ext->neighbor_desc);
  free(part_ext->neighbor_interface);

  if(part_ext->compute_kind == 0) {
    free(part_ext->n_cell);
    free(part_ext->n_cell_border);
    free(part_ext->border_cell_list);

    free(part_ext->cell_cell_idx);
    free(part_ext->cell_cell);

    if(part_ext->face_face_extended_idx != NULL) {
      free(part_ext->face_face_extended_idx);
      free(part_ext->face_face_extended);
      free(part_ext->face_face_interface);
    }

    if(part_ext->opp_interface_and_gnum_face != NULL) {
      free(part_ext->opp_interface_and_gnum_face);
      free(part_ext->cur_interface_face);
      free(part_ext->cur_sens_face);
      free(part_ext->n_cur_interface_face);
    }
    if(part_ext->opp_interface_and_gnum_edge != NULL) {
      free(part_ext->opp_interface_and_gnum_edge);
      free(part_ext->cur_interface_edge);
      free(part_ext->cur_sens_edge);
      free(part_ext->n_cur_interface_edge);
    }
    if(part_ext->opp_interface_and_gnum_vtx != NULL) {
      free(part_ext->opp_interface_and_gnum_vtx);
      free(part_ext->cur_interface_vtx);
      free(part_ext->cur_sens_vtx);
      free(part_ext->n_cur_interface_vtx);
    }


    if(part_ext->edge_edge_extended_idx != NULL) {
      free(part_ext->edge_edge_extended_idx);
      free(part_ext->edge_edge_extended);
      free(part_ext->edge_edge_interface);
    }

    if(part_ext->vtx_vtx_extended_idx != NULL) {
      free(part_ext->vtx_vtx_extended_idx);
      free(part_ext->vtx_vtx_extended);
      free(part_ext->vtx_vtx_interface);
    }

    /* Peu import l'ownership on free car on rend à l'utilisateur l'interface i_domain / i_part */
    if(part_ext->border_cell_face_idx != NULL) {
      free(part_ext->border_cell_face_idx);
      free(part_ext->border_cell_face);
    }

    if(part_ext->border_face_edge_idx != NULL) {
      free(part_ext->border_face_edge_idx);
      free(part_ext->border_face_edge);
    }

    if(part_ext->border_edge_vtx_idx != NULL) {
      free(part_ext->border_edge_vtx_idx);
      free(part_ext->border_edge_vtx);
    }

    if(part_ext->border_face_vtx_idx != NULL) {
      free(part_ext->border_face_vtx_idx);
      free(part_ext->border_face_vtx);
    }

    if(part_ext->border_face_group_idx != NULL) {
      free(part_ext->border_face_group_idx);
      free(part_ext->border_face_group);
    }

    if(part_ext->border_vtx != NULL) {
      free(part_ext->border_vtx);
    }

    if(part_ext->border_cell_ln_to_gn != NULL) {
      free(part_ext->border_cell_ln_to_gn);
    }

    if(part_ext->border_face_ln_to_gn != NULL) {
      free(part_ext->border_face_ln_to_gn);
    }

    if(part_ext->border_edge_ln_to_gn != NULL) {
      free(part_ext->border_edge_ln_to_gn);
    }

    if(part_ext->border_vtx_ln_to_gn != NULL) {
      free(part_ext->border_vtx_ln_to_gn);
    }

    if(part_ext->border_face_group_ln_to_gn != NULL) {
      free(part_ext->border_face_group_ln_to_gn);
    }


    for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {

      int shift_part = 0;
      for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
        for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
          free(part_ext->cell_cell_extended_idx         [i_depth][i_part+shift_part]);
          free(part_ext->cell_cell_extended_n           [i_depth][i_part+shift_part]);
          free(part_ext->cell_cell_extended             [i_depth][i_part+shift_part]);
          free(part_ext->unique_order_cell_cell_extended[i_depth][i_part+shift_part]);
          free(part_ext->cell_cell_interface            [i_depth][i_part+shift_part]);
        }
        shift_part += part_ext->n_part[i_domain];
      }
      free(part_ext->cell_cell_extended_idx           [i_depth]);
      free(part_ext->cell_cell_extended_n             [i_depth]);
      free(part_ext->cell_cell_extended               [i_depth]);
      free(part_ext->unique_order_cell_cell_extended  [i_depth]);
      free(part_ext->cell_cell_interface              [i_depth]);
      free(part_ext->n_unique_order_cell_cell_extended[i_depth]);
    }

    free(part_ext->cell_cell_extended_idx);
    free(part_ext->cell_cell_extended_n);
    free(part_ext->cell_cell_extended);
    free(part_ext->unique_order_cell_cell_extended);
    free(part_ext->cell_cell_interface);
    free(part_ext->n_unique_order_cell_cell_extended);

    /* Only shortcut of user data */
    free(part_ext->entity_cell_idx    );
    free(part_ext->entity_cell        );
    free(part_ext->entity_cell_n);

    free(part_ext->dist_neighbor_cell_n   );
    free(part_ext->dist_neighbor_cell_idx );
    free(part_ext->dist_neighbor_cell_desc);
    free(part_ext->dist_neighbor_cell_interface);
    free(part_ext->unique_order_dist_neighbor_cell);
    free(part_ext->n_unique_order_dist_neighbor_cell);

    /* Allocated by distant neightbor */
    free(part_ext->entity_cell_opp_idx);
    free(part_ext->entity_cell_opp_n  );
    free(part_ext->entity_cell_opp    );

    int shift_part = 0;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        free(part_ext->cell_cell_extended_pruned_idx[i_part+shift_part]);
        free(part_ext->cell_cell_extended_pruned    [i_part+shift_part]);
        if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
          free(part_ext->cell_cell_interface_pruned   [i_part+shift_part]);
        }
      }
      shift_part += part_ext->n_part[i_domain];
    }
    free(part_ext->cell_cell_extended_pruned_idx);
    free(part_ext->cell_cell_extended_pruned    );
    free(part_ext->cell_cell_interface_pruned   );
  }

  free(part_ext->n_part_idx);
  free(part_ext->n_part_g_idx);
  free(part_ext->n_part);
  part_ext->n_part = NULL;

  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    free(part_ext->parts[i_domain]);
  }
  free(part_ext->parts);

  if(part_ext->shift_by_domain_cell != NULL){
    free(part_ext->shift_by_domain_cell);
  }
  if(part_ext->shift_by_domain_face != NULL){
    free(part_ext->shift_by_domain_face);
  }
  if(part_ext->shift_by_domain_edge != NULL){
    free(part_ext->shift_by_domain_edge);
  }
  if(part_ext->shift_by_domain_vtx != NULL){
    free(part_ext->shift_by_domain_vtx);
  }

  if(part_ext->shift_by_domain_face_group != NULL){
    free(part_ext->shift_by_domain_face_group);
  }

  if(part_ext->owner == PDM_OWNERSHIP_KEEP) {
    if(part_ext->composed_interface_idx != NULL) {
      free(part_ext->composed_interface_idx);
    }

    if(part_ext->composed_interface != NULL) {
      free(part_ext->composed_interface);
    }

    if(part_ext->composed_ln_to_gn_sorted != NULL) {
      free(part_ext->composed_ln_to_gn_sorted);
    }
  }

  free(part_ext);
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
      abort();
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
      abort();
    }
    break;

    case PDM_MESH_ENTITY_VTX:
    {
      abort();
    }
    break;

  default :
    PDM_error(__FILE__, __LINE__, 0, "Unknown connectivity_type \n");
    break;

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
  int *block_n = (int *) malloc( block_n_elt * sizeof(int));
  for(int i = 0; i < block_n_elt; ++i) {
    block_n[i] = 1;
  }

  PDM_g_num_t* distrib_adapt = PDM_part_to_block_adapt_partial_block_to_block(ptb, &block_n, max_g_num);
  free(distrib_adapt);

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
  free(block_n);
  free(distrib_cell);

  int          *_n_selected_cell_to_send        = (int         *)  malloc(n_part * sizeof(int         ));
  int         **_selected_cell_to_send_idx      = (int         **) malloc(n_part * sizeof(int         *));
  int         **_selected_cell_to_send          = (int         **) malloc(n_part * sizeof(int         *));
  PDM_g_num_t **_selected_cell_to_send_ln_to_gn = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Compute buffer size
    int n_ext_to_send = 0;
    for(int i = 0; i < n_int_cell[i_part]; ++i) {
      n_ext_to_send += is_ext_cell[i_part][i];
    }

    _n_selected_cell_to_send       [i_part] = n_ext_to_send;
    _selected_cell_to_send_idx     [i_part] = malloc((n_ext_to_send+1) * sizeof(int        ));
    _selected_cell_to_send         [i_part] = malloc( n_ext_to_send    * sizeof(int        ));
    _selected_cell_to_send_ln_to_gn[i_part] = malloc( n_ext_to_send    * sizeof(PDM_g_num_t));

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
    free(is_ext_cell[i_part]);
  }
  free(is_ext_cell);

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
    free(_selected_cell_to_send_idx     [i_part]);
    free(_selected_cell_to_send_ln_to_gn[i_part]);
  }
  free(_selected_cell_to_send_idx);
  free(_selected_cell_to_send_ln_to_gn);

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
      break;
    }

    case PDM_MESH_ENTITY_FACE: {
      part_ext->parts[i_domain][i_part].face_part_bound_proc_idx = part_bound_proc_idx;
      part_ext->parts[i_domain][i_part].face_part_bound_part_idx = part_bound_part_idx;
      part_ext->parts[i_domain][i_part].face_part_bound          = part_bound;
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

