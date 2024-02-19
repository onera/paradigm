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
  part_ext->n_interface = n_interface;

  // for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
  //   part_ext->ptb_itrf[i] = NULL; // (PDM_part_to_block_t **) malloc( n_interface * sizeof(PDM_part_to_block_t **));
  //   part_ext->opp_gnum[i] = NULL; // (PDM_g_num_t         **) malloc( n_interface * sizeof(PDM_g_num_t         **));
  //   part_ext->opp_sens[i] = NULL; // (PDM_g_num_t         **) malloc( n_interface * sizeof(PDM_g_num_t         **));

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
  } else {
    part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX] = 0;
    part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX] = NULL;
    part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX] = NULL;
  }

  int n_data = 0;
  for(int i = 0; i < part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX]; ++i) {
    n_data += part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX][i];
  }

  if(1 == 1) {
    PDM_log_trace_array_long(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX], "dvtx_itrf_blk_gnum            ::");
    PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX], part_ext->dentity_itrf_n_blk[PDM_BOUND_TYPE_VTX], "dvtx_itrf_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX], 2 * n_data                                      , "dvtx_itrf_gnum_and_itrf_data  ::");
    PDM_log_trace_array_int (part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX], n_data                                          , "dvtx_itrf_gnum_and_itrf_sens  ::");
  }


  // free(part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX]);
  // free(part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX]);
  // free(part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX]);
  // free(part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX]);
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
  } else if(part_ext->extend_type == PDM_EXTEND_FROM_FACE) {
    bound_type = PDM_BOUND_TYPE_FACE;
    for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
      pn_entity1        [i_domain] = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
      pentity1_ln_to_gn [i_domain] = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
      for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
        pn_entity1       [i_domain][i_part] = part_ext->parts[i_domain][i_part].n_face;
        pentity1_ln_to_gn[i_domain][i_part] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
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
                                         NULL,
                                         &pentity_bound_to_pentity_bound_idx,
                                         &pentity_bound_to_pentity_bound_triplet,
                                         &pentity_bound_to_pentity_bound_interface,
                                         part_ext->comm);

  if(1 == 1) {
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
    free(pn_entity1       [i_domain]);
    free(pentity1_ln_to_gn[i_domain]);
  }
  free(pn_entity1);
  free(pentity1_ln_to_gn);

  part_ext->pinit_entity_bound_to_pentity_bound_idx       = pentity_bound_to_pentity_bound_idx;
  part_ext->pinit_entity_bound_to_pentity_bound_triplet   = pentity_bound_to_pentity_bound_triplet;
  part_ext->pinit_entity_bound_to_pentity_bound_interface = pentity_bound_to_pentity_bound_interface;

}

static
void
_exchange_coord_and_apply_transform
(
PDM_part_extension_t  *part_ext,
int                   *pn_vtx_extented,
PDM_g_num_t          **pvtx_extented_ln_to_gn,
int                   *pn_vtx,
double               **pvtx_coords,
int                  **pvtx_extented_to_pvtx_idx,
int                  **pvtx_extented_to_pvtx_triplet,
int                  **pvtx_extented_to_pvtx_interface,
double              ***pvtx_extented_coords_out
)
{
  PDM_part_to_part_t* ptp_vtx = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pvtx_extented_ln_to_gn,
                                                                          (const int          *) pn_vtx_extented,
                                                                          part_ext->ln_part_tot,
                                                                          (const int          *) pn_vtx,
                                                                          part_ext->ln_part_tot,
                                                                          (const int         **) pvtx_extented_to_pvtx_idx,
                                                                          (const int         **) NULL,
                                                                          (const int         **) pvtx_extented_to_pvtx_triplet,
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
  double  **translation_vector = malloc(n_interface * sizeof(double *  ));
  double ***rotation_matrix    = malloc(n_interface * sizeof(double ** ));
  double  **rotation_direction = malloc(n_interface * sizeof(double *  ));
  double  **rotation_center    = malloc(n_interface * sizeof(double *  ));
  double   *rotation_angle     = malloc(n_interface * sizeof(double    ));
  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    translation_vector[i_interf] = NULL;
    PDM_part_domain_interface_translation_get(part_ext->pdi, i_interf, &translation_vector[i_interf]);

    rotation_matrix[i_interf] = NULL;
    PDM_part_domain_interface_rotation_get   (part_ext->pdi,
                                              i_interf,
                                              &rotation_direction[i_interf],
                                              &rotation_center   [i_interf],
                                              &rotation_angle    [i_interf]);

    if(rotation_center    [i_interf] != NULL) {
      rotation_matrix[i_interf] = malloc(3 * sizeof(double *));
      for(int k = 0; k < 3; ++k) {
        rotation_matrix[i_interf][k] = malloc(3 * sizeof(double));
      }
    }
  }

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
      int i_interface   = PDM_ABS (pvtx_extented_to_pvtx_interface[i_part][i_vtx]);
      int sgn_interface = PDM_SIGN(pvtx_extented_to_pvtx_interface[i_part][i_vtx]);
      if(i_interface != 0 && translation_vector[PDM_ABS(i_interface)-1] != NULL) {
        for(int k = 0; k < 3; ++k) {
          pextract_vtx_coords[i_part][3*i_vtx+k] += sgn_interface * translation_vector[PDM_ABS(i_interface)-1][k];
        }
      }
    }
  }


  for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
    if(translation_vector[i_interf] != NULL) {
      free(translation_vector[i_interf]);
    }
    if(rotation_center    [i_interf] != NULL) {
      for(int k = 0; k < 3; ++k) {
        free(rotation_matrix[i_interf][k]);
      }
      free(rotation_matrix[i_interf]);
    }
  }
  free(translation_vector);
  free(rotation_matrix);
  free(rotation_direction);
  free(rotation_center);
  free(rotation_angle);

  /*
   * Petit vtk en légende
   */
  if(1 == 1) {
    int i_rank;
    PDM_MPI_Comm_rank(part_ext->comm, &i_rank);

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      char filename[999];
      sprintf(filename, "extented_vtx_coords_%i_%i.vtk", i_part, i_rank);

      PDM_vtk_write_point_cloud(filename,
                                pn_vtx_extented       [i_part],
                                pextract_vtx_coords   [i_part],
                                pvtx_extented_ln_to_gn[i_part],
                                NULL);
    }
  }

  *pvtx_extented_coords_out = pextract_vtx_coords;
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
   * 2 possibilities :
   *   - With face_vtx
   *   - With face_edge + edge_vtx
   */
  if(part_ext->have_edge == 1) {
    printf("Part extension with edges not implemented yet !");
    abort();
  } else {

  }

  int           *pn_vtx         = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
  int           *pn_edge        = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
  int           *pn_face        = (int         * ) malloc( part_ext->n_domain * sizeof(int          ));
  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t **) malloc( part_ext->n_domain * sizeof(PDM_g_num_t *));
  int          **pface_vtx_idx  = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
  int          **pface_vtx      = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
  int          **pedge_vtx_idx  = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
  int          **pedge_vtx      = (int         **) malloc( part_ext->n_domain * sizeof(int         *));
  double       **pvtx_coords    = (double      **) malloc( part_ext->n_domain * sizeof(double      *));

  int lpart = 0;
  for(int i_domain = 0; i_domain < part_ext->n_domain; ++i_domain) {
    for(int i_part = 0; i_part < part_ext->n_part[i_domain]; ++i_part) {
      pn_vtx        [lpart] = part_ext->parts[i_domain][i_part].n_vtx;
      pn_edge       [lpart] = part_ext->parts[i_domain][i_part].n_edge;
      pn_face       [lpart] = part_ext->parts[i_domain][i_part].n_face;
      pvtx_ln_to_gn [lpart] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn;
      pedge_ln_to_gn[lpart] = part_ext->parts[i_domain][i_part].edge_ln_to_gn;
      pface_ln_to_gn[lpart] = part_ext->parts[i_domain][i_part].face_ln_to_gn;
      pface_vtx_idx [lpart] = part_ext->parts[i_domain][i_part].face_vtx_idx;
      pface_vtx     [lpart] = part_ext->parts[i_domain][i_part].face_vtx;
      pedge_vtx_idx [lpart] = NULL;
      pedge_vtx     [lpart] = part_ext->parts[i_domain][i_part].edge_vtx;

      /* Copy to realloc after all step */
      pvtx_ln_to_gn [lpart] = malloc((pn_vtx [lpart]  ) * sizeof(PDM_g_num_t));
      pedge_ln_to_gn[lpart] = malloc((pn_edge[lpart]  ) * sizeof(PDM_g_num_t));
      pface_ln_to_gn[lpart] = malloc((pn_face[lpart]  ) * sizeof(PDM_g_num_t));
      pface_vtx_idx [lpart] = malloc((pn_face[lpart]+1) * sizeof(int        ));
      pface_vtx     [lpart] = malloc(part_ext->parts[i_domain][i_part].face_vtx_idx[pn_face[lpart]]   * sizeof(int));
      pvtx_coords   [lpart] = malloc(3 * pn_vtx [lpart] * sizeof(double));


      for(int i_face = 0; i_face < pn_face[lpart]; ++i_face) {
        pface_ln_to_gn[lpart][i_face] = part_ext->parts[i_domain][i_part].face_ln_to_gn[i_face];
      }
      for(int i_edge = 0; i_edge < pn_edge[lpart]; ++i_edge) {
        pedge_ln_to_gn[lpart][i_edge] = part_ext->parts[i_domain][i_part].edge_ln_to_gn[i_edge];
      }
      for(int i_vtx = 0; i_vtx < pn_vtx[lpart]; ++i_vtx) {
        pvtx_ln_to_gn[lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx_ln_to_gn[i_vtx];
      }

      for(int i_face = 0; i_face < pn_face[lpart]+1; ++i_face) {
        pface_vtx_idx [lpart][i_face] = part_ext->parts[i_domain][i_part].face_vtx_idx[i_face];
      }

      for(int idx = 0; idx < pface_vtx_idx[i_part][pn_face[lpart]]; ++idx) {
        pface_vtx [lpart][idx] = part_ext->parts[i_domain][i_part].face_vtx[idx];
      }

      for(int i_vtx = 0; i_vtx < 3 * pn_vtx[lpart]; ++i_vtx) {
        pvtx_coords   [lpart][i_vtx] = part_ext->parts[i_domain][i_part].vtx[i_vtx];
      }

      lpart++;
    }
  }


  /*
   * On va etendre la partition avec le graphe de base tant que de nouveau elements apparaissent
   *   -> Il faut ajuster la taille du graphe en fonction des nouvelles entités (juste une rallonge)
   *   -> On doit également alimenter un tableau pour le lien avec les entités de la recursion d'après
   */

  PDM_g_num_t *prev_dentity2_elt_gnum           = NULL;
  PDM_g_num_t *prev_dentity2_orig_gnum_and_itrf = NULL;
  PDM_g_num_t *prev_distrib_extented_entity2    = NULL;


  int           *pfull_n_vtx_extented                 = malloc(part_ext->ln_part_tot * sizeof(int         *));
  PDM_g_num_t  **pfull_vtx_extented_ln_to_gn          = malloc(part_ext->ln_part_tot * sizeof(PDM_g_num_t *));
  int          **pfull_vtx_extented_to_pvtx_idx       = malloc(part_ext->ln_part_tot * sizeof(int         *));
  int          **pfull_vtx_extented_to_pvtx_triplet   = malloc(part_ext->ln_part_tot * sizeof(int         *));
  int          **pfull_vtx_extented_to_pvtx_interface = malloc(part_ext->ln_part_tot * sizeof(int         *));

  int           *pfull_n_face_extented                  = malloc(part_ext->ln_part_tot * sizeof(int         *));
  PDM_g_num_t  **pfull_face_extented_ln_to_gn           = malloc(part_ext->ln_part_tot * sizeof(PDM_g_num_t *));
  int          **pfull_face_extented_to_pface_idx       = malloc(part_ext->ln_part_tot * sizeof(int         *));
  int          **pfull_face_extented_to_pface_triplet   = malloc(part_ext->ln_part_tot * sizeof(int         *));
  int          **pfull_face_extented_to_pface_interface = malloc(part_ext->ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    pfull_n_vtx_extented [i_part] = 0;
    pfull_n_face_extented[i_part] = 0;

    pfull_vtx_extented_to_pvtx_idx  [i_part] = malloc((pfull_n_vtx_extented [i_part]+1) * sizeof(int));
    pfull_face_extented_to_pface_idx[i_part] = malloc((pfull_n_face_extented[i_part]+1) * sizeof(int));
    pfull_vtx_extented_to_pvtx_idx  [i_part][0] = 0;
    pfull_face_extented_to_pface_idx[i_part][0] = 0;

    pfull_vtx_extented_ln_to_gn           [i_part] = malloc(0 * sizeof(PDM_g_num_t ));
    pfull_vtx_extented_to_pvtx_triplet    [i_part] = malloc(0 * sizeof(int         ));
    pfull_vtx_extented_to_pvtx_interface  [i_part] = malloc(0 * sizeof(int         ));

    pfull_face_extented_ln_to_gn          [i_part] = malloc(0 * sizeof(PDM_g_num_t ));
    pfull_face_extented_to_pface_triplet  [i_part] = malloc(0 * sizeof(int         ));
    pfull_face_extented_to_pface_interface[i_part] = malloc(0 * sizeof(int         ));
  }

  PDM_g_num_t shift_by_domain_face = part_ext->shift_by_domain_face[part_ext->n_domain];
  PDM_g_num_t shift_by_domain_vtx  = part_ext->shift_by_domain_vtx [part_ext->n_domain];


  int          prev_dface_itrf_n_blk               = 0;
  PDM_g_num_t *prev_dface_itrf_blk_gnum            = NULL;
  int         *prev_dface_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *prev_dface_itrf_gnum_and_itrf_data  = NULL;

  int **pcurr_entity_bound_to_pentity_bound_idx       = part_ext->pinit_entity_bound_to_pentity_bound_idx;
  int **pcurr_entity_bound_to_pentity_bound_triplet   = part_ext->pinit_entity_bound_to_pentity_bound_triplet;
  int **pcurr_entity_bound_to_pentity_bound_interface = part_ext->pinit_entity_bound_to_pentity_bound_interface;

  // int ***pcurr_entity_bound_to_pentity_bound_idx       = malloc( (part_ext->depth+1) * sizeof(int         **));
  // int ***pcurr_entity_bound_to_pentity_bound_triplet   = malloc( (part_ext->depth+1) * sizeof(int         **));
  // int ***pcurr_entity_bound_to_pentity_bound_interface = malloc( (part_ext->depth+1) * sizeof(int         **));

  // /*
  //  * Init
  //  */
  // pcurr_entity_bound_to_pentity_bound_idx      [0] = part_ext->pinit_entity_bound_to_pentity_bound_idx;
  // pcurr_entity_bound_to_pentity_bound_triplet  [0] = part_ext->pinit_entity_bound_to_pentity_bound_triplet;
  // pcurr_entity_bound_to_pentity_bound_interface[0] = part_ext->pinit_entity_bound_to_pentity_bound_interface;

  // for(int i_depth = 1; i_depth < part_ext->depth+1; ++i_depth) {
  //   pcurr_entity_bound_to_pentity_bound_idx      [i_depth] = NULL;
  //   pcurr_entity_bound_to_pentity_bound_triplet  [i_depth] = NULL;
  //   pcurr_entity_bound_to_pentity_bound_interface[i_depth] = NULL;
  // }


  int       **pn_vtx_extented_by_depth   = malloc(part_ext->depth * sizeof(int         *));
  int       **pn_face_extented_by_depth  = malloc(part_ext->depth * sizeof(int         *));
  for(int i_depth = 0; i_depth < part_ext->depth; ++i_depth) {
    pn_vtx_extented_by_depth [i_depth] = malloc(part_ext->ln_part_tot * sizeof(int));
    pn_face_extented_by_depth[i_depth] = malloc(part_ext->ln_part_tot * sizeof(int));

    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      pn_vtx_extented_by_depth [i_depth][i_part] = 0;
      pn_face_extented_by_depth[i_depth][i_part] = 0;
    }

  }


  int i_depth   = 0;
  // int converged = 0;
  int step      = 0;
  while(i_depth < part_ext->depth) {
  // // while(converged == 0) {
  // while(step < 2) {

    /* Use descending connectivity to deduce connectivity and extend_face */
    int          *pn_face_extented                  = NULL;
    PDM_g_num_t **pface_extented_ln_to_gn           = NULL;
    int         **pface_extented_to_pface_idx       = NULL;
    int         **pface_extented_to_pface_triplet   = NULL;
    int         **pface_extented_to_pface_interface = NULL;

    log_trace(" --------------------- step = %i / i_depth = %i \n", step, i_depth);

    /*
     * TODO :
     *   - Keep an block array containaing blk_gnum -> (orig_gnum, interface)
     *   - Use this information to remove alreay faces in the mesh (at the second step of the algorithm)
     */
    int          next_dface_itrf_n_blk               = 0;
    PDM_g_num_t *next_dface_itrf_blk_gnum            = NULL;
    int         *next_dface_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dface_itrf_gnum_and_itrf_data  = NULL;
    log_trace(" PDM_part_extension_entity1_to_entity2 beg \n");
    PDM_part_extension_entity1_to_entity2(shift_by_domain_face, // Attention il va evoluer lui
                                          part_ext->ln_part_tot,
                                          pn_vtx,
                                          pvtx_ln_to_gn,
                                          // pcurr_entity_bound_to_pentity_bound_idx      [i_depth],       // pvtx_to_pvtx_idx,
                                          // pcurr_entity_bound_to_pentity_bound_triplet  [i_depth],   // pvtx_to_pvtx_triplet,
                                          // pcurr_entity_bound_to_pentity_bound_interface[i_depth], // pvtx_to_pvtx_interface,
                                          pcurr_entity_bound_to_pentity_bound_idx      ,       // pvtx_to_pvtx_idx,
                                          pcurr_entity_bound_to_pentity_bound_triplet  ,   // pvtx_to_pvtx_triplet,
                                          pcurr_entity_bound_to_pentity_bound_interface, // pvtx_to_pvtx_interface,
                                          pn_face,
                                          pface_ln_to_gn,
                                          pface_vtx_idx,
                                          pface_vtx,
                                          prev_dface_itrf_n_blk,
                                          prev_dface_itrf_blk_gnum,
                                          prev_dface_itrf_gnum_and_itrf_strid,
                                          prev_dface_itrf_gnum_and_itrf_data,
                                          &pn_face_extented,
                                          &pface_extented_ln_to_gn,
                                          &pface_extented_to_pface_idx,
                                          &pface_extented_to_pface_triplet,
                                          &pface_extented_to_pface_interface,
                                          &next_dface_itrf_n_blk,
                                          &next_dface_itrf_blk_gnum,
                                          &next_dface_itrf_gnum_and_itrf_strid,
                                          &next_dface_itrf_gnum_and_itrf_data,
                                          part_ext->comm);
    log_trace(" PDM_part_extension_entity1_to_entity2 end \n");
    // if(step == 1) {
    //   exit(1);
    // }

    free(prev_dface_itrf_blk_gnum           );
    free(prev_dface_itrf_gnum_and_itrf_strid);
    free(prev_dface_itrf_gnum_and_itrf_data );
    prev_dface_itrf_n_blk               = next_dface_itrf_n_blk;
    prev_dface_itrf_blk_gnum            = next_dface_itrf_blk_gnum;
    prev_dface_itrf_gnum_and_itrf_strid = next_dface_itrf_gnum_and_itrf_strid;
    prev_dface_itrf_gnum_and_itrf_data  = next_dface_itrf_gnum_and_itrf_data;

    /*
     * Update with descending connectivity :
     *   - Mandatory because we need to iterate the connectivity face_vtx (but with the new faces)
     */
    int                           *pn_vtx_extented                 = NULL;
    PDM_g_num_t                  **pvtx_extented_ln_to_gn          = NULL;
    int                          **pextented_face_vtx_idx          = NULL;
    int                          **pextented_face_vtx              = NULL;
    int                          **pvtx_extented_to_pvtx_idx       = NULL;
    int                          **pvtx_extented_to_pvtx_triplet   = NULL;
    int                          **pvtx_extented_to_pvtx_interface = NULL;


    int          next_dvtx_itrf_n_blk               = 0;
    PDM_g_num_t *next_dvtx_itrf_blk_gnum            = NULL;
    int         *next_dvtx_itrf_gnum_and_itrf_strid = NULL;
    PDM_g_num_t *next_dvtx_itrf_gnum_and_itrf_data  = NULL;
    int         *next_dvtx_itrf_gnum_and_itrf_sens  = NULL;
    log_trace(" PDM_part_extension_pentity1_entity2_to_extented_pentity1_entity2 beg \n");
    PDM_part_extension_pentity1_entity2_to_extented_pentity1_entity2(part_ext->ln_part_tot,
                                                                     part_ext->n_interface,
                                                                     shift_by_domain_vtx, // Attention il va evoluer lui
                                                                     // part_ext->ptb_itrf[PDM_BOUND_TYPE_VTX],
                                                                     // part_ext->opp_gnum[PDM_BOUND_TYPE_VTX],
                                                                     // part_ext->opp_sens[PDM_BOUND_TYPE_VTX],
                                                                     // New
                                                                     part_ext->dentity_itrf_n_blk              [PDM_BOUND_TYPE_VTX],
                                                                     part_ext->dentity_itrf_blk_gnum           [PDM_BOUND_TYPE_VTX],
                                                                     part_ext->dentity_itrf_gnum_and_itrf_strid[PDM_BOUND_TYPE_VTX],
                                                                     part_ext->dentity_itrf_gnum_and_itrf_data [PDM_BOUND_TYPE_VTX],
                                                                     part_ext->dentity_itrf_gnum_and_itrf_sens [PDM_BOUND_TYPE_VTX],
                                                                     pn_face,
                                                                     pface_ln_to_gn,
                                                                     pn_vtx,
                                                                     pvtx_ln_to_gn,
                                                                     pface_vtx_idx,
                                                                     pface_vtx,
                                                                     pn_face_extented,
                                                                     pface_extented_ln_to_gn,
                                                                     pface_extented_to_pface_idx,
                                                                     pface_extented_to_pface_triplet,
                                                                     pface_extented_to_pface_interface,
                                                                     &pn_vtx_extented,
                                                                     &pvtx_extented_ln_to_gn,
                                                                     &pextented_face_vtx_idx,
                                                                     &pextented_face_vtx,
                                                                     &pvtx_extented_to_pvtx_idx,
                                                                     &pvtx_extented_to_pvtx_triplet,
                                                                     &pvtx_extented_to_pvtx_interface,
                                                                     // NEW
                                                                     &next_dvtx_itrf_n_blk,
                                                                     &next_dvtx_itrf_blk_gnum,
                                                                     &next_dvtx_itrf_gnum_and_itrf_strid,
                                                                     &next_dvtx_itrf_gnum_and_itrf_data,
                                                                     &next_dvtx_itrf_gnum_and_itrf_sens,
                                                                     part_ext->comm);
    log_trace(" PDM_part_extension_pentity1_entity2_to_extented_pentity1_entity2 end \n");

    if(1 == 1) {
      for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
        PDM_log_trace_array_long(pvtx_extented_ln_to_gn[i_part], pn_vtx_extented[i_part], "pvtx_extented_ln_to_gn ::");
        PDM_log_trace_connectivity_int(pextented_face_vtx_idx[i_part], pextented_face_vtx[i_part], pn_face_extented[i_part], "pextented_face_vtx ::");
      }
    }

    /*
     * Hook coordinates
     */
    double **pvtx_extented_coords = NULL;
    _exchange_coord_and_apply_transform(part_ext,
                                        pn_vtx_extented,
                                        pvtx_extented_ln_to_gn,
                                        pn_vtx,
                                        pvtx_coords,
                                        pvtx_extented_to_pvtx_idx,
                                        pvtx_extented_to_pvtx_triplet,
                                        pvtx_extented_to_pvtx_interface,
                                        &pvtx_extented_coords);

    /*
     * Concatenate all information to continue recursion
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      /* Update size */
      int pn_vtx_extented_old  = pfull_n_vtx_extented [i_part];
      int pn_face_extented_old = pfull_n_face_extented[i_part];
      pfull_n_vtx_extented [i_part] += pn_vtx_extented [i_part];
      pfull_n_face_extented[i_part] += pn_face_extented[i_part];

      pn_vtx_extented_by_depth [i_depth][i_part] += pn_vtx_extented [i_part];
      pn_face_extented_by_depth[i_depth][i_part] += pn_face_extented[i_part];

      /*   */
      int pn_concat_vtx  = pn_vtx [i_part] + pn_vtx_extented [i_part];
      int pn_concat_edge = pn_edge[i_part]; // + pn_edge_extented[i_part];
      int pn_concat_face = pn_face[i_part] + pn_face_extented[i_part];

      int pn_concat_face_vtx_idx = pface_vtx_idx [i_part][pn_face[i_part]] + pextented_face_vtx_idx[i_part][pn_face_extented[i_part]];

      /* Realloc */
      pvtx_ln_to_gn [i_part] = realloc(pvtx_ln_to_gn [i_part],     pn_concat_vtx          * sizeof(PDM_g_num_t));
      pvtx_coords   [i_part] = realloc(pvtx_coords   [i_part], 3 * pn_concat_vtx          * sizeof(double     ));
      pedge_ln_to_gn[i_part] = realloc(pedge_ln_to_gn[i_part],     pn_concat_edge         * sizeof(PDM_g_num_t));
      pface_ln_to_gn[i_part] = realloc(pface_ln_to_gn[i_part],     pn_concat_face         * sizeof(PDM_g_num_t));
      pface_vtx_idx [i_part] = realloc(pface_vtx_idx [i_part],     (pn_concat_face+1)     * sizeof(int        ));
      pface_vtx     [i_part] = realloc(pface_vtx     [i_part],     pn_concat_face_vtx_idx * sizeof(int        ));

      /* Concatenate */
      for(int i_face = 0; i_face < pn_face_extented[i_part]; ++i_face) {
        pface_ln_to_gn[i_part][pn_face[i_part]+i_face] = pface_extented_ln_to_gn[i_part][i_face];
        shift_by_domain_face = PDM_MAX(shift_by_domain_face, pface_extented_ln_to_gn[i_part][i_face]);
      }

      // for(int i_edge = ; i_edge < pn_edge_extented[i_part]; ++i_edge) {
      //   pedge_ln_to_gn[i_part][pn_edge[i_part]+i_edge] = pedge_extented_ln_to_gn[i_part][i_edge];
      // }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
        pvtx_ln_to_gn[i_part][pn_vtx[i_part]+i_vtx] = pvtx_extented_ln_to_gn[i_part][i_vtx];
        shift_by_domain_vtx = PDM_MAX(shift_by_domain_vtx, pvtx_extented_ln_to_gn[i_part][i_vtx]);
      }

      for(int i_vtx = 0; i_vtx < 3 * pn_vtx_extented[i_part]; ++i_vtx) {
        pvtx_coords[i_part][3*pn_vtx[i_part]+i_vtx] = pvtx_extented_coords[i_part][i_vtx];
      }

      /* Concatenate graphe in other array */
      for(int i_face = 0; i_face < pn_face_extented[i_part]; ++i_face) {
        int ln_vtx = pextented_face_vtx_idx[i_part][i_face+1] - pextented_face_vtx_idx[i_part][i_face];
        pface_vtx_idx [i_part][pn_face[i_part]+i_face+1] = pface_vtx_idx [i_part][pn_face[i_part]+i_face] + ln_vtx;
      }

      /* Concatenate graphe in other array */
      for(int i = 0; i < pextented_face_vtx_idx[i_part][pn_face_extented[i_part]]; ++i) {
        pface_vtx[i_part][pface_vtx_idx [i_part][pn_face[i_part]]+i] = pextented_face_vtx[i_part][i];
      }

      int size_vtx_vtx = (pfull_vtx_extented_to_pvtx_idx[i_part][pn_vtx_extented_old] + pvtx_extented_to_pvtx_idx[i_part][pn_vtx_extented[i_part]])/3;
      pfull_vtx_extented_ln_to_gn         [i_part] = realloc(pfull_vtx_extented_ln_to_gn         [i_part],  pfull_n_vtx_extented[i_part]    * sizeof(PDM_g_num_t));
      pfull_vtx_extented_to_pvtx_idx      [i_part] = realloc(pfull_vtx_extented_to_pvtx_idx      [i_part], (pfull_n_vtx_extented[i_part]+1) * sizeof(int        ));
      pfull_vtx_extented_to_pvtx_triplet  [i_part] = realloc(pfull_vtx_extented_to_pvtx_triplet  [i_part], 3 * size_vtx_vtx                 * sizeof(int        ));
      pfull_vtx_extented_to_pvtx_interface[i_part] = realloc(pfull_vtx_extented_to_pvtx_interface[i_part],     size_vtx_vtx                 * sizeof(int        ));

      int size_face_face = (pfull_face_extented_to_pface_idx[i_part][pn_face_extented_old] + pface_extented_to_pface_idx[i_part][pn_face_extented[i_part]])/3;
      pfull_face_extented_ln_to_gn          [i_part] = realloc(pfull_face_extented_ln_to_gn          [i_part],  pfull_n_face_extented[i_part]    * sizeof(PDM_g_num_t));
      pfull_face_extented_to_pface_idx      [i_part] = realloc(pfull_face_extented_to_pface_idx      [i_part], (pfull_n_face_extented[i_part]+1) * sizeof(int        ));
      pfull_face_extented_to_pface_triplet  [i_part] = realloc(pfull_face_extented_to_pface_triplet  [i_part], 3 * size_face_face                * sizeof(int        ));
      pfull_face_extented_to_pface_interface[i_part] = realloc(pfull_face_extented_to_pface_interface[i_part],     size_face_face                * sizeof(int        ));

      for(int i_face = 0; i_face < pn_face_extented[i_part]; ++i_face) {
        int ln_face = pface_extented_to_pface_idx[i_part][i_face+1] - pface_extented_to_pface_idx[i_part][i_face];
        pfull_face_extented_to_pface_idx[i_part][pn_face_extented_old+i_face+1] = pfull_face_extented_to_pface_idx[i_part][pn_face_extented_old+i_face] + ln_face;
      }

      for(int i = 0; i < pface_extented_to_pface_idx[i_part][pn_face_extented[i_part]]/3; ++i) {
        int idx_write = pfull_face_extented_to_pface_idx[i_part][pn_face_extented_old]/3 + i;
        pfull_face_extented_to_pface_triplet  [i_part][3*idx_write  ] = pface_extented_to_pface_triplet  [i_part][3*i  ];
        pfull_face_extented_to_pface_triplet  [i_part][3*idx_write+1] = pface_extented_to_pface_triplet  [i_part][3*i+1];
        pfull_face_extented_to_pface_triplet  [i_part][3*idx_write+2] = pface_extented_to_pface_triplet  [i_part][3*i+2];
        pfull_face_extented_to_pface_interface[i_part][  idx_write  ] = pface_extented_to_pface_interface[i_part][  i  ];
      }

      for(int i_vtx = 0; i_vtx < pn_vtx_extented[i_part]; ++i_vtx) {
        int ln_vtx = pvtx_extented_to_pvtx_idx[i_part][i_vtx+1] - pvtx_extented_to_pvtx_idx[i_part][i_vtx];
        pfull_vtx_extented_to_pvtx_idx[i_part][pn_vtx_extented_old+i_vtx+1] = pfull_vtx_extented_to_pvtx_idx[i_part][pn_vtx_extented_old+i_vtx] + ln_vtx;
      }

      for(int i = 0; i < pvtx_extented_to_pvtx_idx[i_part][pn_vtx_extented[i_part]]/3; ++i) {
        int idx_write = pfull_vtx_extented_to_pvtx_idx[i_part][pn_vtx_extented_old]/3 + i;
        pfull_vtx_extented_to_pvtx_triplet  [i_part][3*idx_write  ] = pvtx_extented_to_pvtx_triplet  [i_part][3*i  ];
        pfull_vtx_extented_to_pvtx_triplet  [i_part][3*idx_write+1] = pvtx_extented_to_pvtx_triplet  [i_part][3*i+1];
        pfull_vtx_extented_to_pvtx_triplet  [i_part][3*idx_write+2] = pvtx_extented_to_pvtx_triplet  [i_part][3*i+2];
        pfull_vtx_extented_to_pvtx_interface[i_part][  idx_write  ] = pvtx_extented_to_pvtx_interface[i_part][  i  ];
      }

      if(1 == 1) {
        PDM_log_trace_array_int(pfull_vtx_extented_to_pvtx_idx[i_part], pfull_n_vtx_extented[i_part]+1  , "pfull_vtx_extented_to_pvtx_idx     ::");
        PDM_log_trace_array_int(pfull_face_extented_to_pface_triplet  [i_part], 3 * size_face_face, "pfull_face_extented_to_pface_triplet   ::");
        PDM_log_trace_array_int(pfull_face_extented_to_pface_interface[i_part],     size_face_face, "pfull_face_extented_to_pface_interface ::");
        PDM_log_trace_array_int(pfull_vtx_extented_to_pvtx_triplet    [i_part], 3 * size_vtx_vtx  , "pfull_vtx_extented_to_pvtx_triplet     ::");
        PDM_log_trace_array_int(pfull_vtx_extented_to_pvtx_interface  [i_part],     size_vtx_vtx  , "pfull_vtx_extented_to_pvtx_interface   ::");
      }

      /*
       * Update graphe - Only extend idx since no connextion is create AT this stage
       */
      int pn_concat_entity_extented = 0;
      int pn_entity_extented        = 0;
      int pn_entity                 = 0;
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
        pn_concat_entity_extented = pn_concat_vtx;
        pn_entity_extented        = pn_vtx_extented[i_part];
        pn_entity                 = pn_vtx         [i_part];
      } else if(part_ext->extend_type == PDM_EXTEND_FROM_EDGE) {
        pn_concat_entity_extented = pn_concat_edge;
        abort();
        // pn_entity_extented        = pn_edge_extented[i_part];
        // pn_entity                 = pn_edge         [i_part];
      }
      log_trace("pn_concat_entity_extented = %i \n", pn_concat_entity_extented);
      pcurr_entity_bound_to_pentity_bound_idx[i_part] = realloc(pcurr_entity_bound_to_pentity_bound_idx[i_part], (pn_concat_entity_extented+1) * sizeof(int));
      for(int i = 0; i < pn_entity_extented; ++i) {
        pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity];
      }

      // pcurr_entity_bound_to_pentity_bound_idx[i_depth][i_part] = realloc(pcurr_entity_bound_to_pentity_bound_idx[i_depth][i_part], (pn_concat_entity_extented+1) * sizeof(int));
      // for(int i = 0; i < pn_entity_extented; ++i) {
      //   pcurr_entity_bound_to_pentity_bound_idx[i_depth][i_part][pn_entity+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_depth][i_part][pn_entity];
      // }

      if(1 == 1) {
        int i_rank;
        PDM_MPI_Comm_rank(part_ext->comm, &i_rank);
        char filename[999];
        sprintf(filename, "out_face_vtx_step=%i_%i_%i.vtk", step, i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               pn_concat_vtx,
                               pvtx_coords   [i_part],
                               pvtx_ln_to_gn [i_part],
                               pn_concat_face,
                               pface_vtx_idx [i_part],
                               pface_vtx     [i_part],
                               pface_ln_to_gn[i_part],
                               NULL);

      }
    } /* End loop part */

    /*
     * Free coords
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      free(pvtx_extented_coords[i_part]);
    }
    free(pvtx_extented_coords);

    /*
     * Update shift_by_domain_face
     */
    PDM_g_num_t _shift_by_domain_face = shift_by_domain_face;
    PDM_MPI_Allreduce(&_shift_by_domain_face, &shift_by_domain_face, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    PDM_g_num_t _shift_by_domain_vtx = shift_by_domain_vtx;
    PDM_MPI_Allreduce(&_shift_by_domain_vtx, &shift_by_domain_vtx, 1, PDM_MPI_INT, PDM_MPI_MAX, part_ext->comm);

    PDM_g_num_t _pn_face_extented_tot = 0;
    PDM_g_num_t  pn_face_extented_tot = 0;
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      _pn_face_extented_tot += pn_face_extented[i_part];
    }

    PDM_MPI_Allreduce(&_pn_face_extented_tot, &pn_face_extented_tot, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, part_ext->comm);

    log_trace("pn_face_extented_tot = %i \n", pn_face_extented_tot);

    /*
     * Free
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      free(pvtx_extented_ln_to_gn           [i_part]);
      free(pextented_face_vtx_idx           [i_part]);
      free(pextented_face_vtx               [i_part]);
      free(pvtx_extented_to_pvtx_idx        [i_part]);
      free(pvtx_extented_to_pvtx_triplet    [i_part]);
      free(pvtx_extented_to_pvtx_interface  [i_part]);
      free(pface_extented_ln_to_gn          [i_part]);
      free(pface_extented_to_pface_idx      [i_part]);
      free(pface_extented_to_pface_triplet  [i_part]);
      free(pface_extented_to_pface_interface[i_part]);
    }
    free(pvtx_extented_ln_to_gn         );
    free(pextented_face_vtx_idx         );
    free(pextented_face_vtx             );
    free(pvtx_extented_to_pvtx_idx      );
    free(pvtx_extented_to_pvtx_triplet  );
    free(pvtx_extented_to_pvtx_interface);
    free(pface_extented_ln_to_gn          );
    free(pface_extented_to_pface_idx      );
    free(pface_extented_to_pface_triplet  );
    free(pface_extented_to_pface_interface);


    /*
     * A chaque étape :
     *   - On garde le même graphe entre les entitiés, mais on agrandit le tableau idx (pour être cohérent avec le part_to_part )
     *   - A l'issu d'une étape, il faut swap le graph avec celui de la nouvelle depth
     *   - Pour avoir l'historique complet on peut agglomerer tout les graphe de chaque depth to have the full one
     */


    if(pn_face_extented_tot == 0) {
      // Change graph
      if(part_ext->extend_type == PDM_EXTEND_FROM_VTX) {
        // pcurr_entity_bound_to_pentity_bound_idx       = pfull_vtx_extented_to_pvtx_idx;
        pcurr_entity_bound_to_pentity_bound_triplet   = pfull_vtx_extented_to_pvtx_triplet;
        pcurr_entity_bound_to_pentity_bound_interface = pfull_vtx_extented_to_pvtx_interface;
        //

        // Ne marche pas si depth == 3 il faut uniquement le detph-1
        if(i_depth == 2) {
          abort();
        }

        for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
          int pn_vtx_tot = pn_vtx [i_part] + pn_vtx_extented [i_part];

          // Deja fait plus haut
          // pcurr_entity_bound_to_pentity_bound_idx[i_part] = realloc(pcurr_entity_bound_to_pentity_bound_idx[i_part], (pn_vtx_tot+1) * sizeof(int));
          // for(int i = 0; i < pn_entity_extented; ++i) {
          //   pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_entity];
          // }
          // int *pcurr_entity_bound_to_pentity_bound_n = PDM_array_zeros_int(pn_vtx_tot);
          pcurr_entity_bound_to_pentity_bound_idx[i_part][0] = 0;

          for(int i = 0; i < pn_vtx [i_part]; ++i) {
            pcurr_entity_bound_to_pentity_bound_idx[i_part][i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][i];
          }

          PDM_log_trace_array_int(pfull_vtx_extented_to_pvtx_idx[i_part], pn_vtx_extented_by_depth[i_depth][i_part]+1, "pfull_vtx_extented_to_pvtx_idx ::");

          int pn_vtx_old_depth = pn_vtx [i_part] - pn_vtx_extented_by_depth[i_depth][i_part];
          for(int i = 0; i < pn_vtx_extented_by_depth[i_depth][i_part]; ++i) {
            int ln_curr_centity = pfull_vtx_extented_to_pvtx_idx[i_part][i+1] - pfull_vtx_extented_to_pvtx_idx[i_part][i];
            pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_vtx_old_depth+i+1] = pcurr_entity_bound_to_pentity_bound_idx[i_part][pn_vtx_old_depth+i] + ln_curr_centity;
          }

          // free(pcurr_entity_bound_to_pentity_bound_n);
        }

      } else {
        abort();
      }

      i_depth++;
    }
    step++;


    /*
     * Update for next step
     */
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {

      /* Update size */
      pn_vtx [i_part] += pn_vtx_extented [i_part];
      pn_edge[i_part] += 0; // + pn_edge_extented[i_part];
      pn_face[i_part] += pn_face_extented[i_part];
    }
    free(pn_vtx_extented);
    free(pn_face_extented);

    if(step > 2) {
      abort();
    }

  }

  free(prev_dface_itrf_blk_gnum           );
  free(prev_dface_itrf_gnum_and_itrf_strid);
  free(prev_dface_itrf_gnum_and_itrf_data );

  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    free(pvtx_ln_to_gn [i_part]);
    free(pedge_ln_to_gn[i_part]);
    free(pface_ln_to_gn[i_part]);
    free(pface_vtx_idx [i_part]);
    free(pface_vtx     [i_part]);
    free(pvtx_coords   [i_part]);
  }

  free(pn_vtx        );
  free(pn_edge       );
  free(pn_face       );
  free(pvtx_ln_to_gn );
  free(pedge_ln_to_gn);
  free(pface_ln_to_gn);
  free(pface_vtx_idx );
  free(pface_vtx     );
  free(pedge_vtx_idx );
  free(pedge_vtx     );
  free(pvtx_coords   );

  /*
   * To keep
   */
  for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
    free(pfull_vtx_extented_ln_to_gn         [i_part]);
    free(pfull_vtx_extented_to_pvtx_idx      [i_part]);
    free(pfull_vtx_extented_to_pvtx_triplet  [i_part]);
    free(pfull_vtx_extented_to_pvtx_interface[i_part]);

    free(pfull_face_extented_ln_to_gn          [i_part]);
    free(pfull_face_extented_to_pface_idx      [i_part]);
    free(pfull_face_extented_to_pface_triplet  [i_part]);
    free(pfull_face_extented_to_pface_interface[i_part]);

  }
  free(pfull_n_vtx_extented                );
  free(pfull_vtx_extented_ln_to_gn         );
  free(pfull_vtx_extented_to_pvtx_idx      );
  free(pfull_vtx_extented_to_pvtx_triplet  );
  free(pfull_vtx_extented_to_pvtx_interface);

  free(pfull_n_face_extented                 );
  free(pfull_face_extented_ln_to_gn          );
  free(pfull_face_extented_to_pface_idx      );
  free(pfull_face_extented_to_pface_triplet  );
  free(pfull_face_extented_to_pface_interface);



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
    // if(part_ext->ptb_itrf[i] != NULL) {
    //   for(int i_itrf = 0; i_itrf < part_ext->n_interface; ++i_itrf) {
    //     if(part_ext->ptb_itrf[i][i_itrf] != NULL) {
    //       PDM_part_to_block_free(part_ext->ptb_itrf[i][i_itrf]);
    //       free(part_ext->opp_gnum[i][i_itrf]);
    //       free(part_ext->opp_sens[i][i_itrf]);
    //     }
    //   }
    // }
    // free(part_ext->ptb_itrf[i]);
    // free(part_ext->opp_sens[i]);
    // free(part_ext->opp_gnum[i]);

    if(part_ext->dentity_itrf_blk_gnum[i] != NULL) {
      free(part_ext->dentity_itrf_blk_gnum           [i]);
      free(part_ext->dentity_itrf_gnum_and_itrf_strid[i]);
      free(part_ext->dentity_itrf_gnum_and_itrf_data [i]);
      free(part_ext->dentity_itrf_gnum_and_itrf_sens [i]);
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


  if(part_ext->pinit_entity_bound_to_pentity_bound_idx != NULL) {
    for(int i_part = 0; i_part < part_ext->ln_part_tot; ++i_part) {
      free(part_ext->pinit_entity_bound_to_pentity_bound_idx      [i_part]);
      free(part_ext->pinit_entity_bound_to_pentity_bound_triplet  [i_part]);
      free(part_ext->pinit_entity_bound_to_pentity_bound_interface[i_part]);
    }
    free(part_ext->pinit_entity_bound_to_pentity_bound_idx      );
    free(part_ext->pinit_entity_bound_to_pentity_bound_triplet  );
    free(part_ext->pinit_entity_bound_to_pentity_bound_interface);
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

