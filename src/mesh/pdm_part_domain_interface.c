/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_distant_neighbor.h"
#include "pdm_order.h"

#include "pdm_part_domain_interface.h"
#include "pdm_part_domain_interface_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/

static inline
int
_is_same_quadruplet
(
int iproc1, int ipart1, int ielt1, int iinterf1,
int iproc2, int ipart2, int ielt2, int iinterf2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        if(iinterf1 == iinterf2){
          return 1;
        }
      }
    }
  }
  return 0;
}

static
void
_unique_quadruplet
(
  int   n_entity,
  int  *neighbor_entity_idx,
  int  *neighbor_entity,
  int **unique_neighbor_entity_idx,
  int **unique_neighbor_entity_n,
  int **unique_neighbor_entity
)
{

  int* _unique_neighbor_entity_idx = malloc( (n_entity + 1) * sizeof(int));
  int* _unique_neighbor_entity_n   = malloc( (n_entity    ) * sizeof(int));
  int* _unique_neighbor_entity     = malloc( 4 * neighbor_entity_idx[n_entity] * sizeof(int));
  int* order                       = malloc(     neighbor_entity_idx[n_entity] * sizeof(int)); // Suralloc

  _unique_neighbor_entity_idx[0] = 0;
  for(int i_entity = 0; i_entity < n_entity; ++i_entity) {

    int beg       = neighbor_entity_idx[i_entity];
    int n_connect = neighbor_entity_idx[i_entity+1] - beg;

    PDM_order_lnum_s(&neighbor_entity[4*beg], 4, order, n_connect);

    _unique_neighbor_entity_n  [i_entity  ] = 0;
    _unique_neighbor_entity_idx[i_entity+1] = _unique_neighbor_entity_idx[i_entity];

    int last_proc  = -1;
    int last_part  = -1;
    int last_elmt  = -1;
    int last_inte  = -40;
    for(int i = 0; i < n_connect; ++i) {
      int old_order   = order[i];
      int curr_proc   = neighbor_entity[4*(beg+old_order)  ];
      int curr_part   = neighbor_entity[4*(beg+old_order)+1];
      int curr_entity = neighbor_entity[4*(beg+old_order)+2];
      int curr_inte   = neighbor_entity[4*(beg+old_order)+3];
      int is_same  = _is_same_quadruplet(last_proc, last_part, last_elmt  , last_inte,
                                         curr_proc, curr_part, curr_entity, curr_inte);

      if(is_same == 0){ // N'est pas le meme
        // idx_unique++;
        last_proc = curr_proc;
        last_part = curr_part;
        last_elmt = curr_entity;
        last_inte = curr_inte;

        int beg_write = 4 * _unique_neighbor_entity_idx[i_entity+1];
        // printf("beg_write = %i | curr_proc = %i | curr_part = %i | curr_entity = %i \n", beg_write, curr_proc, curr_part, curr_entity);
        _unique_neighbor_entity[beg_write  ] = curr_proc;
        _unique_neighbor_entity[beg_write+1] = curr_part;
        _unique_neighbor_entity[beg_write+2] = curr_entity;
        _unique_neighbor_entity[beg_write+3] = curr_inte;

        /* Increment the new counter */
        _unique_neighbor_entity_idx[i_entity+1]++;
        _unique_neighbor_entity_n  [i_entity  ]++;
      }
    }
  }

  _unique_neighbor_entity = realloc(_unique_neighbor_entity, 4 * neighbor_entity_idx[n_entity] * sizeof(int));

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  free(order);
}

static
int
_concatenate_neighbor
(
  int             n_part,
  int            *n_entity,
  int           **neighbor_idx,
  int           **init_neighbor_idx,
  int           **init_neighbor_desc,
  int           **next_neighbor_n,
  int           **next_neighbor_desc,
  int          ***concat_neighbor_idx,
  int          ***concat_neighbor_desc
)
{

  int **_concat_neighbor_idx  = malloc(n_part * sizeof(int *));
  int **_concat_neighbor_desc = malloc(n_part * sizeof(int *));

  int is_same = 1;
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int   n_elmt           = n_entity    [i_part];
    int  *_neighbor_idx    = neighbor_idx[i_part];

    int  *_init_neighbor_idx  = init_neighbor_idx [i_part];
    int  *_init_neighbor_desc = init_neighbor_desc[i_part];

    int  *_next_neighbor_n    = next_neighbor_n   [i_part];
    int  *_next_neighbor_desc = next_neighbor_desc[i_part];


    int n_tot_next = 0;
    for(int i = 0; i < n_elmt; ++i) {
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        n_tot_next += _next_neighbor_n[idx];
      }
    }

    int n_concat_idx_tot = _init_neighbor_idx[n_elmt] + n_tot_next;
    _concat_neighbor_idx [i_part] = PDM_array_zeros_int(n_elmt+1);
    _concat_neighbor_desc[i_part] = malloc(4 * n_concat_idx_tot * sizeof(int));

    int idx_read  = 0;
    _concat_neighbor_idx[i_part][0] = 0;
    for(int i = 0; i < n_elmt; ++i) {
      _concat_neighbor_idx[i_part][i+1] = _concat_neighbor_idx[i_part][i];

      /* Current part */
      for(int idx = _init_neighbor_idx[i]; idx < _init_neighbor_idx[i+1]; ++idx) {
        int idx_write = _concat_neighbor_idx[i_part][i+1]++;

        _concat_neighbor_desc[i_part][4*idx_write  ] = _init_neighbor_desc[4*idx  ];
        _concat_neighbor_desc[i_part][4*idx_write+1] = _init_neighbor_desc[4*idx+1];
        _concat_neighbor_desc[i_part][4*idx_write+2] = _init_neighbor_desc[4*idx+2];
        _concat_neighbor_desc[i_part][4*idx_write+3] = _init_neighbor_desc[4*idx+3];

      }

      /* Opp part */
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        for(int k = 0; k < _next_neighbor_n[idx]; ++k) {
          int idx_write = _concat_neighbor_idx[i_part][i+1]++;

          _concat_neighbor_desc[i_part][4*idx_write  ] = _next_neighbor_desc[4*idx_read  ];
          _concat_neighbor_desc[i_part][4*idx_write+1] = _next_neighbor_desc[4*idx_read+1];
          _concat_neighbor_desc[i_part][4*idx_write+2] = _next_neighbor_desc[4*idx_read+2];
          _concat_neighbor_desc[i_part][4*idx_write+3] = _next_neighbor_desc[4*idx_read+3];
          idx_read++;
        }
      }
    }

    /* Unique */
    int *_unique_concat_neighbor_idx  = NULL;
    int *_unique_concat_neighbor_n    = NULL;
    int *_unique_concat_neighbor_desc = NULL;
    _unique_quadruplet(n_elmt,
                       _concat_neighbor_idx [i_part],
                       _concat_neighbor_desc[i_part],
                       &_unique_concat_neighbor_idx,
                       &_unique_concat_neighbor_n,
                       &_unique_concat_neighbor_desc);


    for(int i = 0; i < n_elmt; ++i) {
      int n_old_connect = _init_neighbor_idx[i+1] - _init_neighbor_idx[i];
      if(n_old_connect != _unique_concat_neighbor_n[i]) {
        is_same = 0;
        break;
      }
    }

    // printf("is_same = %i\n", is_same);

    free(_concat_neighbor_idx [i_part]);
    free(_concat_neighbor_desc[i_part]);
    free(_unique_concat_neighbor_n);
    _concat_neighbor_idx [i_part] = _unique_concat_neighbor_idx;
    _concat_neighbor_desc[i_part] = _unique_concat_neighbor_desc;

    /* Debug */
    // PDM_log_trace_graph_nuplet_int(_concat_neighbor_idx[i_part], _concat_neighbor_desc[i_part], 4, n_elmt, "_concat_neighbor_desc :");


  }

  *concat_neighbor_idx  = _concat_neighbor_idx;
  *concat_neighbor_desc = _concat_neighbor_desc;

  return is_same;
}


static
void
_exchange_and_sort_neighbor
(
  PDM_MPI_Comm    comm,
  int             n_part,
  int            *n_entity,
  int           **neighbor_idx,
  int           **neighbor_desc,
  int           **neighbor_interface,
  int           **init_neighbor_idx,
  int           **init_neighbor_desc,
  int          ***all_neighbor_idx,
  int          ***all_neighbor_desc
)
{

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(comm,
                                                           n_part,
                                                           n_entity,
                                                           neighbor_idx,
                                                           neighbor_desc);

  int **prev_neighbor_idx  = init_neighbor_idx;
  int **prev_neighbor_desc = init_neighbor_desc;

  int is_same    = 0;
  int i_step     = 0;
  int first_step = 1;
  int **concat_neighbor_opp_idx = NULL;
  int **concat_neighbor_opp     = NULL;
  while(is_same != 1) {

    /*
     * Conpute stride
     */
    int **prev_neighbor_n = malloc(n_part * sizeof(int *));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      prev_neighbor_n[i_part] = malloc(n_entity[i_part] * sizeof(int));
      for(int i = 0; i < n_entity[i_part]; ++i) {
        prev_neighbor_n[i_part][i] = prev_neighbor_idx[i_part][i+1] - prev_neighbor_idx[i_part][i];
      }
    }

    int **next_neighbor_opp_n = NULL;
    int **next_neighbor_opp   = NULL;
    PDM_distant_neighbor_exch(dn,
                              4 * sizeof(int),
                              PDM_STRIDE_VAR_INTERLACED,
                              -1,
                              prev_neighbor_n,
                    (void **) prev_neighbor_desc,
                             &next_neighbor_opp_n,
                   (void ***)&next_neighbor_opp);


    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(prev_neighbor_n[i_part]);
    }
    free(prev_neighbor_n);

    is_same = _concatenate_neighbor(n_part,
                                    n_entity,
                                    neighbor_idx,
                                    prev_neighbor_idx,
                                    prev_neighbor_desc,
                                    next_neighbor_opp_n,
                                    next_neighbor_opp,
                                    &concat_neighbor_opp_idx,
                                    &concat_neighbor_opp);

    int is_same_l = is_same;
    PDM_MPI_Allreduce (&is_same_l, &is_same, 1, PDM_MPI_INT, PDM_MPI_MIN, comm);

    /* Init next step */
    if(!first_step) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        free(prev_neighbor_idx  [i_part]);
        free(prev_neighbor_desc [i_part]);
      }
      free(prev_neighbor_idx);
      free(prev_neighbor_desc);
    }
    prev_neighbor_idx  = concat_neighbor_opp_idx;
    prev_neighbor_desc = concat_neighbor_opp;


    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(next_neighbor_opp_n[i_part]);
      free(next_neighbor_opp  [i_part]);
    }
    free(next_neighbor_opp_n);
    free(next_neighbor_opp);

    first_step = 0;
    i_step++;
    if(i_step > 50) {
      abort();
    }
  }

  PDM_distant_neighbor_free(dn);

  /*
   * Filter by removing self and already direct neighbor
   */
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int **filter_neighbor_idx  = malloc(n_part * sizeof(int *));
  int **filter_neighbor_desc = malloc(n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *_concat_neighbor_opp_idx = concat_neighbor_opp_idx[i_part];
    int *_concat_neighbor_opp     = concat_neighbor_opp    [i_part];

    int *_neighbor_idx   = neighbor_idx      [i_part];
    int *_neighbor_desc  = neighbor_desc     [i_part];
    int *_neighbor_intrf = neighbor_interface[i_part];

    filter_neighbor_idx [i_part] = malloc(     (n_entity[i_part] + 1)                       * sizeof(int));
    filter_neighbor_desc[i_part] = malloc( 4 * (_concat_neighbor_opp_idx[n_entity[i_part]]) * sizeof(int));

    int *_filter_neighbor_idx = filter_neighbor_idx [i_part];
    int *_filter_neighbor     = filter_neighbor_desc[i_part];

    _filter_neighbor_idx[0] = 0;
    for(int i = 0; i < n_entity[i_part]; ++i) {

      _filter_neighbor_idx[i+1] = _filter_neighbor_idx[i];

      for(int idx = _concat_neighbor_opp_idx[i]; idx < _concat_neighbor_opp_idx[i+1]; ++idx) {

        int curr_proc = _concat_neighbor_opp[4*idx  ];
        int curr_part = _concat_neighbor_opp[4*idx+1];
        int curr_elmt = _concat_neighbor_opp[4*idx+2];
        int curr_inte = _concat_neighbor_opp[4*idx+3];

        /* Rm if current elemt is inside */
        int is_define_in_direct_neight = 0;

        if(curr_proc == i_rank && curr_part == i_part && curr_elmt == i) {
          continue;
        }


        /* Brut force */
        for(int idx2 = _neighbor_idx[i]; idx2 < _neighbor_idx[i+1]; ++idx2) {
          int opp_proc = _neighbor_desc[3*idx2  ];
          int opp_part = _neighbor_desc[3*idx2+1];
          int opp_elmt = _neighbor_desc[3*idx2+2];
          int opp_inte = curr_inte; // Normal because i can come from another interface

          is_define_in_direct_neight = _is_same_quadruplet(curr_proc, curr_part, curr_elmt, curr_inte,
                                                           opp_proc , opp_part , opp_elmt , opp_inte);

          if(is_define_in_direct_neight) {
            break;
          }
        }

        if(is_define_in_direct_neight) {
          continue;
        }

        int idx_write = _filter_neighbor_idx[i+1]++;
        _filter_neighbor[4*idx_write  ] = _concat_neighbor_opp[4*idx  ];
        _filter_neighbor[4*idx_write+1] = _concat_neighbor_opp[4*idx+1];
        _filter_neighbor[4*idx_write+2] = _concat_neighbor_opp[4*idx+2];
        _filter_neighbor[4*idx_write+3] = _concat_neighbor_opp[4*idx+3];

      }

      /* On re-rajoute les neighbor */
      for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
        int idx_write = _filter_neighbor_idx[i+1]++;
        _filter_neighbor[4*idx_write  ] = _neighbor_desc[3*idx  ];
        _filter_neighbor[4*idx_write+1] = _neighbor_desc[3*idx+1];
        _filter_neighbor[4*idx_write+2] = _neighbor_desc[3*idx+2];
        _filter_neighbor[4*idx_write+3] = _neighbor_intrf[idx];
      }



    }

    /*
     * Realloc
     */
    filter_neighbor_desc[i_part] = realloc(filter_neighbor_desc[i_part], 4 * (_filter_neighbor_idx[n_entity[i_part]]) * sizeof(int));
    free(_concat_neighbor_opp_idx);
    free(_concat_neighbor_opp);

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(_filter_neighbor_idx, filter_neighbor_desc[i_part], 4, n_entity[i_part], "filter_neighbor_desc OOOO :");
    }


  }

  free(concat_neighbor_opp_idx);
  free(concat_neighbor_opp    );

  *all_neighbor_idx  = filter_neighbor_idx;
  *all_neighbor_desc = filter_neighbor_desc;

}


/*============================================================================
 * Public function definitions
 *============================================================================*/
PDM_part_domain_interface_t*
PDM_part_domain_interface_create
(
const int                    n_interface,
const int                    n_domain,
const int                   *n_part,
PDM_domain_interface_mult_t  multidomain_interface,
PDM_ownership_t              ownership,
PDM_MPI_Comm                 comm
)
{

  PDM_part_domain_interface_t *dom_intrf = (PDM_part_domain_interface_t *) malloc (sizeof(PDM_part_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  dom_intrf->n_part            = malloc(n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->n_part[i_domain] = n_part[i_domain];
  }
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_pn_vtx        = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_vtx_ln_to_gn  = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx_idx   = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));

  dom_intrf->interface_pn_edge        = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_edge_ln_to_gn  = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_edge       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_edge       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_edge_idx   = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_edge       = (int         ****) malloc(n_domain * sizeof(int         ***));

  dom_intrf->interface_pn_face       = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_face_ln_to_gn = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face_idx  = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_face      = (int         ****) malloc(n_domain * sizeof(int         ***));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->interface_pn_vtx       [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_vtx_ln_to_gn [i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx_idx  [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    dom_intrf->interface_pn_edge       [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_edge_ln_to_gn [i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_edge      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_edge      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_edge_idx  [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_edge      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    dom_intrf->interface_pn_face      [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_face_ln_to_gn[i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face_idx [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      dom_intrf->interface_pn_vtx       [i_domain][i_part] = (int          *) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
      dom_intrf->interface_sgn_vtx      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_vtx      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_dom_vtx      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));

      dom_intrf->interface_pn_edge       [i_domain][i_part] = (int          *) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_edge_ln_to_gn [i_domain][i_part] = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
      dom_intrf->interface_sgn_edge      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_edge      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_edge_idx  [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_dom_edge      [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));

      dom_intrf->interface_pn_face      [i_domain][i_part] = (int          *) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = (PDM_g_num_t **) malloc(n_interface * sizeof(PDM_g_num_t *));
      dom_intrf->interface_sgn_face     [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_face     [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_ids_face_idx [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));
      dom_intrf->interface_dom_face     [i_domain][i_part] = (int         **) malloc(n_interface * sizeof(int         *));

      for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
        dom_intrf->interface_pn_vtx       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interf] = NULL;

        dom_intrf->interface_pn_edge       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_edge_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_edge      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_edge      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_edge_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_edge      [i_domain][i_part][i_interf] = NULL;

        dom_intrf->interface_pn_face       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_face_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_face      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_face      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_face_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_face      [i_domain][i_part][i_interf] = NULL;
      }
    }
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    dom_intrf->is_result[i] = 0;
  }

  dom_intrf->translation_vect   = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_direction = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_center    = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_angle     = (double  *) malloc(n_interface * sizeof(double  ));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    dom_intrf->translation_vect  [i_interface] = NULL;
    dom_intrf->rotation_direction[i_interface] = NULL;
    dom_intrf->rotation_center   [i_interface] = NULL;
    dom_intrf->rotation_angle    [i_interface] = 0.;
  }

  dom_intrf->interface_describe_by_face = 0;
  dom_intrf->interface_describe_by_edge = 0;
  dom_intrf->interface_describe_by_vtx  = 0;

  return dom_intrf;
}

void
PDM_part_domain_interface_set
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind,
 int                           i_domain,
 int                           i_part,
 int                           i_interface,
 int                           interface_pn,
 PDM_g_num_t                  *interface_ln_to_gn,
 int                          *interface_sgn,
 int                          *interface_ids,
 int                          *interface_ids_idx,
 int                          *interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_pn_vtx      [i_domain][i_part][i_interface] = interface_pn;
    dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part][i_interface] = interface_ln_to_gn;
    dom_intrf->interface_sgn_vtx     [i_domain][i_part][i_interface] = interface_sgn;
    dom_intrf->interface_ids_vtx     [i_domain][i_part][i_interface] = interface_ids;
    dom_intrf->interface_ids_vtx_idx [i_domain][i_part][i_interface] = interface_ids_idx;
    dom_intrf->interface_dom_vtx     [i_domain][i_part][i_interface] = interface_dom;
    dom_intrf->interface_describe_by_vtx  = 1;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    dom_intrf->interface_pn_edge      [i_domain][i_part][i_interface] = interface_pn;
    dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface] = interface_ln_to_gn;
    dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface] = interface_sgn;
    dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface] = interface_ids;
    dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface] = interface_ids_idx;
    dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface] = interface_dom;
    dom_intrf->interface_describe_by_edge  = 1;
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_pn_face      [i_domain][i_part][i_interface] = interface_pn;
    dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface] = interface_ln_to_gn;
    dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface] = interface_sgn;
    dom_intrf->interface_ids_face     [i_domain][i_part][i_interface] = interface_ids;
    dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface] = interface_ids_idx;
    dom_intrf->interface_dom_face     [i_domain][i_part][i_interface] = interface_dom;
    dom_intrf->interface_describe_by_face  = 1;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

void
PDM_part_domain_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind,
 int                            i_domain,
 int                            i_part,
 int                            i_interface,
 int                           *interface_pn,
 PDM_g_num_t                  **interface_ln_to_gn,
 int                          **interface_sgn,
 int                          **interface_ids,
 int                          **interface_ids_idx,
 int                          **interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    *interface_pn       = dom_intrf->interface_pn_vtx      [i_domain][i_part][i_interface];
    *interface_ln_to_gn = dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part][i_interface];
    *interface_sgn      = dom_intrf->interface_sgn_vtx     [i_domain][i_part][i_interface];
    *interface_ids      = dom_intrf->interface_ids_vtx     [i_domain][i_part][i_interface];
    *interface_ids_idx  = dom_intrf->interface_ids_vtx_idx [i_domain][i_part][i_interface];
    *interface_dom      = dom_intrf->interface_dom_vtx     [i_domain][i_part][i_interface];
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    *interface_pn       = dom_intrf->interface_pn_edge      [i_domain][i_part][i_interface];
    *interface_ln_to_gn = dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface];
    *interface_sgn      = dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface];
    *interface_ids      = dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface];
    *interface_ids_idx  = dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface];
    *interface_dom      = dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface];
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    *interface_pn       = dom_intrf->interface_pn_face      [i_domain][i_part][i_interface];
    *interface_ln_to_gn = dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface];
    *interface_sgn      = dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface];
    *interface_ids      = dom_intrf->interface_ids_face     [i_domain][i_part][i_interface];
    *interface_ids_idx  = dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface];
    *interface_dom      = dom_intrf->interface_dom_face     [i_domain][i_part][i_interface];
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

int
PDM_part_domain_interface_n_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf
)
{
  return dom_intrf->n_interface;
}


int
PDM_part_domain_interface_exist_get
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind
)
{
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    return dom_intrf->interface_describe_by_vtx;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    return dom_intrf->interface_describe_by_edge;
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    return dom_intrf->interface_describe_by_face;
  }
  return 0;
}

void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
)
{

  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {

    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      if(dom_intrf->ownership == PDM_OWNERSHIP_KEEP) {
        if(dom_intrf->interface_pn_vtx       [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_vtx       [i_domain][i_part]);
          dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part]);
          dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_vtx      [i_domain][i_part]);
          dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx      [i_domain][i_part]);
          dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part]);
          dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_vtx      [i_domain][i_part]);
          dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
        };

        if(dom_intrf->interface_pn_edge      [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_edge      [i_domain][i_part]);
          dom_intrf->interface_pn_edge      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_edge_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_edge     [i_domain][i_part]);
          dom_intrf->interface_sgn_edge     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_edge     [i_domain][i_part]);
          dom_intrf->interface_ids_edge     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_edge_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_edge_idx [i_domain][i_part]);
          dom_intrf->interface_ids_edge_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_edge     [i_domain][i_part]);
          dom_intrf->interface_dom_edge     [i_domain][i_part] = NULL;
        };

        if(dom_intrf->interface_pn_face      [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_face      [i_domain][i_part]);
          dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_face_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_face     [i_domain][i_part]);
          dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face     [i_domain][i_part]);
          dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face_idx [i_domain][i_part]);
          dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_face     [i_domain][i_part]);
          dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
        };



      }
    }

    free(dom_intrf->interface_pn_vtx       [i_domain]);
    free(dom_intrf->interface_vtx_ln_to_gn [i_domain]);
    free(dom_intrf->interface_sgn_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx_idx  [i_domain]);
    free(dom_intrf->interface_dom_vtx      [i_domain]);

    free(dom_intrf->interface_pn_edge      [i_domain]);
    free(dom_intrf->interface_edge_ln_to_gn[i_domain]);
    free(dom_intrf->interface_sgn_edge     [i_domain]);
    free(dom_intrf->interface_ids_edge     [i_domain]);
    free(dom_intrf->interface_ids_edge_idx [i_domain]);
    free(dom_intrf->interface_dom_edge     [i_domain]);

    free(dom_intrf->interface_pn_face      [i_domain]);
    free(dom_intrf->interface_face_ln_to_gn[i_domain]);
    free(dom_intrf->interface_sgn_face     [i_domain]);
    free(dom_intrf->interface_ids_face     [i_domain]);
    free(dom_intrf->interface_ids_face_idx [i_domain]);
    free(dom_intrf->interface_dom_face     [i_domain]);

  }
  free(dom_intrf->n_part);

  free(dom_intrf->interface_pn_vtx       );
  free(dom_intrf->interface_vtx_ln_to_gn );
  free(dom_intrf->interface_sgn_vtx      );
  free(dom_intrf->interface_ids_vtx      );
  free(dom_intrf->interface_ids_vtx_idx  );
  free(dom_intrf->interface_dom_vtx      );

  free(dom_intrf->interface_pn_edge       );
  free(dom_intrf->interface_edge_ln_to_gn );
  free(dom_intrf->interface_sgn_edge      );
  free(dom_intrf->interface_ids_edge      );
  free(dom_intrf->interface_ids_edge_idx  );
  free(dom_intrf->interface_dom_edge      );

  free(dom_intrf->interface_pn_face      );
  free(dom_intrf->interface_face_ln_to_gn);
  free(dom_intrf->interface_sgn_face     );
  free(dom_intrf->interface_ids_face     );
  free(dom_intrf->interface_ids_face_idx );
  free(dom_intrf->interface_dom_face     );

  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    if(dom_intrf->translation_vect[i_interface]   != NULL) {
      free(dom_intrf->translation_vect[i_interface]);
      dom_intrf->translation_vect[i_interface] = NULL;
    }
    if(dom_intrf->rotation_direction[i_interface]   != NULL) {
      free(dom_intrf->rotation_direction[i_interface]);
      dom_intrf->rotation_direction[i_interface] = NULL;
    }
    if(dom_intrf->rotation_center[i_interface]   != NULL) {
      free(dom_intrf->rotation_center[i_interface]);
      dom_intrf->rotation_center[i_interface] = NULL;
    }
  }

  free(dom_intrf->translation_vect  );
  free(dom_intrf->rotation_direction);
  free(dom_intrf->rotation_center   );
  free(dom_intrf->rotation_angle    );

  free(dom_intrf);
}




void
PDM_part_domain_interface_translation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
  const double                       *vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->translation_vect[i_interface] == NULL);

  dom_intrf->translation_vect[i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->translation_vect[i_interface][i] = vect[i];
  }

}

void
PDM_part_domain_interface_rotation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
  const double                       *direction,
  const double                       *center,
  const double                        angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->rotation_direction[i_interface] == NULL);
  assert(dom_intrf->rotation_center   [i_interface] == NULL);

  dom_intrf->rotation_direction[i_interface] = (double *) malloc( 3 * sizeof(double));
  dom_intrf->rotation_center   [i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->rotation_direction[i_interface][i] = direction[i];
    dom_intrf->rotation_center   [i_interface][i] = center   [i];
  }
  dom_intrf->rotation_angle[i_interface] = angle;
}



void
PDM_part_domain_interface_translation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
        double                      **vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->translation_vect[i_interface] != NULL){

    *vect = (double *) malloc( 3 * sizeof(double));
    double* _vect = *vect;

    for(int i = 0; i < 3; ++i) {
      _vect[i] = dom_intrf->translation_vect[i_interface][i];
    }
  } else {
    *vect = NULL;
  }
}

void
PDM_part_domain_interface_rotation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->rotation_direction[i_interface] != NULL) {
    assert(dom_intrf->rotation_center   [i_interface] != NULL);

    *direction = (double *) malloc( 3 * sizeof(double));
    *center    = (double *) malloc( 3 * sizeof(double));
    double *_direction = *direction;
    double *_center    = *center   ;

    for(int i = 0; i < 3; ++i) {
      _direction[i] = dom_intrf->rotation_direction[i_interface][i];
      _center   [i] = dom_intrf->rotation_center   [i_interface][i];
    }
    *angle = dom_intrf->rotation_angle[i_interface];
  } else {
    *direction = NULL;
    *center    = NULL;
    *angle     = 0;
  }
}


void
PDM_part_domain_interface_as_graph
(
  PDM_part_domain_interface_t    *dom_intrf,
  PDM_bound_type_t                interface_kind,
  int                           **n_entity,
  PDM_g_num_t                  ***entity_ln_to_gn,
  int                          ***neighbor_entity_idx,
  int                          ***neighbor_entity_desc,
  int                            *n_g_interface,
  int                           **all_composed_id_idx,
  int                           **all_composed_id,
  PDM_g_num_t                   **all_composed_ln_to_gn_sorted
)
{
  PDM_UNUSED(entity_ln_to_gn);

  int i_rank;
  PDM_MPI_Comm_rank(dom_intrf->comm, &i_rank);

  // assert(dom_intrf->n_domain == 1); // TODO --> shift of gnum AND part_id
  if(dom_intrf->n_domain > 1) {
    printf("WARNING : PDM_part_domain_interface_as_graph not general is n_domain > 1");
  }

  int* n_tot_part_by_domain = (int *) malloc( dom_intrf->n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain) {

    n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = dom_intrf->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, dom_intrf->comm);
  }

  /* Si multidomain on fait un shift et tt roule */
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain) {
    n_part_loc_all_domain += dom_intrf->n_part[i_domain];
  }

  int n_interface = PDM_part_domain_interface_n_interface_get(dom_intrf);

  int **neighbor_n         = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_idx       = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_desc      = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_interface = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int  *n_entity_bound     = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  int **neighbor_opp_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_opp_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_opp_desc = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );

  /*
   * Loop over all interfaces to create distant neighbor structure
   */
  int shift_part   = 0;
  int shift_part_g = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {
    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      n_entity_bound[i_part+shift_part] = n_entity[i_domain][i_part];
      neighbor_idx  [i_part+shift_part] = (int *) malloc( (n_entity_bound[i_part+shift_part]+1) * sizeof(int) );
      neighbor_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_n   = neighbor_n    [i_part+shift_part];
      int* _neighbor_idx = neighbor_idx  [i_part+shift_part];

      neighbor_opp_idx  [i_part+shift_part] = (int *) malloc( (n_entity_bound[i_part+shift_part]+1) * sizeof(int) );
      neighbor_opp_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_opp_n   = neighbor_opp_n    [i_part+shift_part];
      int* _neighbor_opp_idx = neighbor_opp_idx  [i_part+shift_part];

      int           *interface_pn       = malloc(n_interface * sizeof(int          ));
      PDM_g_num_t  **interface_ln_to_gn = malloc(n_interface * sizeof(PDM_g_num_t *));
      int          **interface_sgn      = malloc(n_interface * sizeof(int         *));
      int          **interface_ids      = malloc(n_interface * sizeof(int         *));
      int          **interface_ids_idx  = malloc(n_interface * sizeof(int         *));
      int          **interface_dom      = malloc(n_interface * sizeof(int         *));
      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
        PDM_part_domain_interface_get(dom_intrf,
                                      interface_kind,
                                      i_domain,
                                      i_part,
                                      i_interface,
                                      &interface_pn[i_interface],
                                      &interface_ln_to_gn[i_interface],
                                      &interface_sgn[i_interface],
                                      &interface_ids[i_interface],
                                      &interface_ids_idx[i_interface],
                                      &interface_dom[i_interface]);
      }

      /*
       * First step : Count interface to add in distant neighbor due to connectivity betwenn domain
       */
      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

        // log_trace("-------------------------------- i_interface = %i  -------------------------------- \n", i_interface);
        // PDM_log_trace_array_int(interface_sgn[i_interface], interface_pn[i_interface], "interface_sgn :: ");

        for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

          // Search the first in list that is in current part/proc
          // int i_proc_cur   = -1;
          // int i_part_cur   = -1;
          int i_entity_cur = -1;
          int found        = 0;
          int idx_current  = -1;
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
              // i_proc_cur   = i_proc_opp;
              // i_part_cur   = i_part_opp;
              i_entity_cur = i_entity_opp;
              idx_current  = j;
              found = 1;
            }
          }

          if(!found) {
            continue;
          }

          // Il manque une notion de direction sinon on sait pas dans quelle sens va le raccord

          assert(found == 1);

          // log_trace("i_proc_cur = %i | i_part_cur = %i | i_entity_cur = %i | sgn = %i \n", i_proc_cur, i_part_cur, i_entity_cur, interface_sgn[i_interface][idx_entity]);

          // Only add the opposite part of the graph
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            // int i_proc_opp   = interface_ids[i_interface][3*j  ];
            // int i_part_opp   = interface_ids[i_interface][3*j+1];
            // int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(idx_current != j) {
              // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
              _neighbor_n[i_entity_cur] += 1;
              _neighbor_opp_n[i_entity_cur] += 1;
            // } else {
            //   _neighbor_opp_n[i_entity_cur] += 1;
            }
          }
        }
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      _neighbor_opp_idx[0] = 0;
      for(int i_entity = 0; i_entity < n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_idx    [i_entity+1] = _neighbor_idx    [i_entity] + _neighbor_n    [i_entity];
        _neighbor_opp_idx[i_entity+1] = _neighbor_opp_idx[i_entity] + _neighbor_opp_n[i_entity];
        _neighbor_n    [i_entity] = 0;
        _neighbor_opp_n[i_entity] = 0;
      }

      neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx    [n_entity_bound[i_part+shift_part]] * sizeof(int) );
      neighbor_opp_desc [i_part+shift_part] = (int *) malloc( 4 * _neighbor_opp_idx[n_entity_bound[i_part+shift_part]] * sizeof(int) );
      neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[n_entity_bound[i_part+shift_part]] * sizeof(int) );
      int* _neighbor_desc      = neighbor_desc     [i_part+shift_part];
      int* _neighbor_opp_desc  = neighbor_opp_desc [i_part+shift_part];
      int* _neighbor_interface = neighbor_interface[i_part+shift_part];

      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
        for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

          // Search the first in list that is in current part/proc
          // int i_proc_cur   = -1;
          // int i_part_cur   = -1;
          int i_entity_cur = -1;
          int found        = 0;
          int idx_current  = -1;
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
              // i_proc_cur   = i_proc_opp;
              // i_part_cur   = i_part_opp;
              i_entity_cur = i_entity_opp;
              idx_current  = j;
              found = 1;
            }
          }

          if(!found) {
            continue;
          }

          assert(found == 1);

          // Only add the opposite part of the graph
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(idx_current != j) {
              // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
              int idx_write = _neighbor_idx[i_entity_cur] + _neighbor_n[i_entity_cur]++;
              _neighbor_desc[3*idx_write  ] = i_proc_opp;              // i_proc_opp;
              _neighbor_desc[3*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
              _neighbor_desc[3*idx_write+2] = i_entity_opp;            // i_entity_opp
              _neighbor_interface[idx_write] = (i_interface+1) * interface_sgn[i_interface][idx_entity];

              idx_write = _neighbor_opp_idx[i_entity_cur] + _neighbor_opp_n[i_entity_cur]++;
              _neighbor_opp_desc[4*idx_write  ] = i_proc_opp;              // i_proc_opp;
              _neighbor_opp_desc[4*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
              _neighbor_opp_desc[4*idx_write+2] = i_entity_opp;            // i_entity_opp
              _neighbor_opp_desc[4*idx_write+3] = (i_interface+1) * interface_sgn[i_interface][idx_entity];


            // } else {
            //   int idx_write = _neighbor_opp_idx[i_entity_cur] + _neighbor_opp_n[i_entity_cur]++;
            //   _neighbor_opp_desc[4*idx_write  ] = i_proc_opp;              // i_proc_opp;
            //   _neighbor_opp_desc[4*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
            //   _neighbor_opp_desc[4*idx_write+2] = i_entity_opp;            // i_entity_opp
            //   _neighbor_opp_desc[4*idx_write+3] = - (i_interface+1) * interface_sgn[i_interface][idx_entity];
            }
          }
        }
      }

      if(0 == 1) {
        PDM_log_trace_graph_nuplet_int(_neighbor_idx, _neighbor_desc, 3, n_entity_bound[i_part+shift_part], "_neighbor graph :");
        PDM_log_trace_graph_nuplet_int(_neighbor_opp_idx, _neighbor_opp_desc, 4, n_entity_bound[i_part+shift_part], "_neighbor_opp_desc :");
      }

      free(interface_pn      );
      free(interface_ln_to_gn);
      free(interface_sgn     );
      free(interface_ids     );
      free(interface_ids_idx );
      free(interface_dom     );

    }
    shift_part   += dom_intrf->n_part[i_domain];
    shift_part_g += n_tot_part_by_domain[i_domain];
  }


  /*
   * Au moment ou on propage il faut appliqué le signe de la transformation du raccord maybe ?
   */
  _exchange_and_sort_neighbor(dom_intrf->comm,
                              n_part_loc_all_domain,
                              n_entity_bound,
                              neighbor_idx,
                              neighbor_desc,
                              neighbor_interface,
                              neighbor_opp_idx,
                              neighbor_opp_desc,
                              neighbor_entity_idx,
                              neighbor_entity_desc);


  /*
   * At this stage we have interface composition AND
   * All combination is a real new connection
   * We need to fused all combinaison in a sigle interface and keep the original link
   */

  /*
   * Traduce sign of all interface
   */
  int *sgn_interf_to_interf = malloc(2 * n_interface * sizeof(int));

  int j = 0;
  for(int i = n_interface; i > 0 ; --i) {
    sgn_interf_to_interf[j++] = -i;
  }
  for(int i = 0; i < n_interface ; ++i) {
    sgn_interf_to_interf[n_interface+i] = (i+1);
  }
  // PDM_log_trace_array_int(sgn_interf_to_interf, 2 * n_interface, "sgn_interf_to_interf ::");

  int i_composed_interface = 0;
  int max_composed  = 0;
  int max_lcomposed = 0;
  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    int n_elmt = n_entity_bound[i_part];
    int *_neighbor_entity_idx  = (*neighbor_entity_idx)[i_part];
    max_composed += _neighbor_entity_idx[n_elmt];

    for(int i = 0; i < n_elmt; ++i) {
      int dn = _neighbor_entity_idx[i+1] - _neighbor_entity_idx[i];
      max_lcomposed = PDM_MAX(max_lcomposed, dn);
    }
  }

  int         *composed_id_idx         = malloc(    (max_composed+1) * sizeof(int        ));
  int         *composed_id             = malloc(     max_composed    * sizeof(int        ));
  int         *composed_id_tmp         = malloc(     max_lcomposed   * sizeof(int        ));
  PDM_g_num_t *composed_key            = malloc(     max_composed    * sizeof(PDM_g_num_t));
  int         *composed_key_update_idx = malloc( 2 * max_composed    * sizeof(PDM_g_num_t));

  composed_id_idx[0] = 0;
  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

    int n_elmt = n_entity_bound[i_part];
    int *_neighbor_entity_idx  = (*neighbor_entity_idx )[i_part];
    int *_neighbor_entity_desc = (*neighbor_entity_desc)[i_part];

    int *_filter_neighbor_entity_idx = malloc((n_elmt+1) * sizeof(int));

    _filter_neighbor_entity_idx[0] = 0;
    for(int i = 0; i < n_elmt; ++i) {

      _filter_neighbor_entity_idx[i+1] = _filter_neighbor_entity_idx[i];
      /* Exclude void conectivity rapidly */
      if(_neighbor_entity_idx[i] == _neighbor_entity_idx[i+1]) {
        continue;
      }

      int idx_fisrt = _neighbor_entity_idx[i];
      int first_proc = _neighbor_entity_desc[4*idx_fisrt  ];
      int first_part = _neighbor_entity_desc[4*idx_fisrt+1];
      int first_elmt = _neighbor_entity_desc[4*idx_fisrt+2];
      int first_inte = _neighbor_entity_desc[4*idx_fisrt+3];

      int pos_first = PDM_binary_search_int(first_inte, sgn_interf_to_interf, 2 * n_interface);
      assert(pos_first != -1);

      int composed_interface = n_interface + (pos_first+1); // n_interface is here to not conflict with previous one

      int i_tmp = 0;
      composed_id_tmp[i_tmp++] = first_inte;

      for(int idx = _neighbor_entity_idx[i]+1; idx < _neighbor_entity_idx[i+1]; ++idx){

        int next_proc = _neighbor_entity_desc[4*idx  ];
        int next_part = _neighbor_entity_desc[4*idx+1];
        int next_elmt = _neighbor_entity_desc[4*idx+2];
        int next_inte = _neighbor_entity_desc[4*idx+3];

        if(first_proc == next_proc &&
           first_part == next_part &&
           first_elmt == next_elmt &&
           first_inte != next_inte) {

          // Donc même element mais interface différente --> On compose

          int pos_next = PDM_binary_search_int(next_inte, sgn_interf_to_interf, 2 * n_interface);
          assert(pos_next != -1);

          composed_interface += (pos_next+1) * ( 2 * n_interface);

          // printf("Hit !!!! (%i %i %i %i) -> %i / (%i %i %i %i) -> %i \n",
          //        first_proc, first_part, first_elmt, first_inte, pos_first,
          //        next_proc , next_part , next_elmt , next_inte , pos_next);
          // printf(" --> Composed  :  %i \n", composed_interface);

          composed_id_tmp[i_tmp++] = next_inte;

        } else {

          /* Create the new array - inplace */
          composed_id_idx[i_composed_interface+1] = composed_id_idx[i_composed_interface];
          int idx_write = _filter_neighbor_entity_idx[i+1]++;

          _neighbor_entity_desc[4*idx_write  ] = first_proc;
          _neighbor_entity_desc[4*idx_write+1] = first_part;
          _neighbor_entity_desc[4*idx_write+2] = first_elmt;

          if(composed_interface == n_interface + (pos_first+1)) {
            _neighbor_entity_desc[4*idx_write+3] = first_inte;
          } else {
            _neighbor_entity_desc[4*idx_write+3] = composed_interface;

            for(int k = 0; k < i_tmp; ++k ) {
              int idx_write_composed = composed_id_idx[i_composed_interface+1]++;
              composed_id[idx_write_composed] = composed_id_tmp[k];
            }

            composed_key_update_idx[2*i_composed_interface  ] = i_part;
            composed_key_update_idx[2*i_composed_interface+1] = idx_write;
            composed_key           [i_composed_interface++] = composed_interface;
          }

          first_proc = next_proc;
          first_part = next_part;
          first_elmt = next_elmt;
          first_inte = next_inte;
          pos_first = PDM_binary_search_int(first_inte, sgn_interf_to_interf, 2 * n_interface);
          composed_interface = n_interface + (pos_first+1);

          /* Reset */
          i_tmp = 0;
          composed_id_tmp[i_tmp++] = first_inte;

        }
      }

      /*
       * Management of last
       */
      if(composed_interface == n_interface + (pos_first+1)) {
        int idx_write = _filter_neighbor_entity_idx[i+1]++;

        _neighbor_entity_desc[4*idx_write  ] = first_proc;
        _neighbor_entity_desc[4*idx_write+1] = first_part;
        _neighbor_entity_desc[4*idx_write+2] = first_elmt;
        _neighbor_entity_desc[4*idx_write+3] = first_inte;
      }

    }

    free(_neighbor_entity_idx);
    (*neighbor_entity_idx)[i_part] = _filter_neighbor_entity_idx;

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(_filter_neighbor_entity_idx, _neighbor_entity_desc, 4, n_elmt, "_neighbor_entity_desc :");
    }
  }

  /*
   * Realloc
   */
  free(composed_id_tmp);
  composed_id_idx = realloc(composed_id_idx, (i_composed_interface+1)              * sizeof(int        ));
  composed_id     = realloc(composed_id    , composed_id_idx[i_composed_interface] * sizeof(int        ));
  composed_key    = realloc(composed_key   , (i_composed_interface+1)              * sizeof(PDM_g_num_t));

  /*
   * Generate table to give from interface number the composed one
   *   - We need the g_id of new interface and the associate composition
   */
  // PDM_gen_gnum_t* gen_gnum_interf = PDM_gnum_create(3, 1, PDM_FALSE, 1.e-6, dom_intrf->comm, PDM_OWNERSHIP_USER);

  // PDM_gnum_set_from_parents(gen_gnum_interf, 0, i_composed_interface, composed_key);
  // PDM_gnum_compute(gen_gnum_interf);

  // PDM_g_num_t* composed_id_gnum = PDM_gnum_get(gen_gnum_interf, 0);

  // PDM_gnum_free(gen_gnum_interf);
  PDM_g_num_t* composed_id_gnum = composed_key;

  if(0 == 1) {
    PDM_log_trace_array_long(composed_key    , i_composed_interface, "composed_id :: ");
    PDM_log_trace_array_long(composed_id_gnum, i_composed_interface, "composed_id_gnum :: ");
    PDM_log_trace_connectivity_int(composed_id_idx, composed_id, i_composed_interface, "composed_id :: ");
  }


  // No update because we need FOR all entity the same interface global numbering so we keep the rules
  // for(int i = 0; i < i_composed_interface; ++i) {
  //   int i_part     = composed_key_update_idx[2*i  ];
  //   int idx_update = composed_key_update_idx[2*i+1];
  //   (*neighbor_entity_desc[i_part])[4*idx_update+3] = n_interface + composed_id_gnum[i];
  // }

  /*
   * Panic verbose
   */
  if(0 == 1) {
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

      int n_elmt = n_entity_bound[i_part];
      int *_neighbor_entity_idx  = (*neighbor_entity_idx )[i_part];
      int *_neighbor_entity_desc = (*neighbor_entity_desc)[i_part];
      PDM_log_trace_graph_nuplet_int(_neighbor_entity_idx, _neighbor_entity_desc, 4, n_elmt, "_neighbor_entity_desc :");
    }
  }
  free(composed_key_update_idx);

  int n_rank = -1;
  PDM_MPI_Comm_size(dom_intrf->comm, &n_rank);

  int *composed_id_n = malloc( i_composed_interface * sizeof(int        ));
  PDM_g_num_t max_loc = 0;
  for(int i = 0; i < i_composed_interface; ++i) {
    max_loc = PDM_MAX(max_loc, composed_id_gnum[i]);
    composed_id_n[i] = composed_id_idx[i+1] - composed_id_idx[i];
  }
  free(composed_id_idx);
  PDM_g_num_t max_glob = -1;
  PDM_MPI_Allreduce(&max_loc, &max_glob, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, dom_intrf->comm);

  PDM_g_num_t* distrib_interf = malloc( (n_rank + 1) * sizeof(PDM_g_num_t));
  distrib_interf[0] = 0;
  for(int i = 1; i < n_rank+1; ++i) {
    distrib_interf[i] = max_glob;
  }

  /*
   * Each proc can create a global id but all procs is suceptible to know the decomposition of all interface
   *  --> Unified
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                   1.,
                                                                   &composed_id_gnum,
                                                                   distrib_interf,
                                                                   &i_composed_interface,
                                                                   1,
                                                                   dom_intrf->comm);

  // free(composed_id_gnum);
  int *_composed_id_n = NULL;
  int *_composed_id   = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         &composed_id_n,
             (void **)   &composed_id,
                         &_composed_id_n,
             (void **)   &_composed_id);

  free(composed_id);
  free(composed_id_n);

  int _n_g_interface = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *ptb_composed_ln_to_gn_sorted = PDM_part_to_block_block_gnum_get(ptb);

  if(i_rank > 0) {
    assert(_n_g_interface == 0);
    free(_composed_id);
  }


  free(distrib_interf);

  free(composed_key);
  free(sgn_interf_to_interf);


  /*
   *  Broadcast to all
   */
  PDM_MPI_Bcast(&_n_g_interface  , 1                                , PDM_MPI_INT, 0, dom_intrf->comm);

  PDM_g_num_t *_composed_ln_to_gn_sorted     = malloc(_n_g_interface * sizeof(PDM_g_num_t));

  if(i_rank == 0) {
    for(int i = 0; i < _n_g_interface; ++i) {
      _composed_ln_to_gn_sorted[i] =  ptb_composed_ln_to_gn_sorted[i];
    }
  }

  int *_composed_id_idx = malloc( (_n_g_interface+1) * sizeof(int));

  if(i_rank == 0) {
    _composed_id_idx[0] = 0;
    for(int i = 0; i < _n_g_interface; ++i) {
      _composed_id_idx[i+1] = _composed_id_idx[i] + _composed_id_n[i];
    }
  }
  free(_composed_id_n);

  PDM_part_to_block_free(ptb);

  PDM_MPI_Bcast(_composed_ln_to_gn_sorted    , _n_g_interface, PDM__PDM_MPI_G_NUM, 0, dom_intrf->comm);

  PDM_MPI_Bcast(_composed_id_idx, _n_g_interface+1                , PDM_MPI_INT, 0, dom_intrf->comm);


  if(i_rank != 0) {
    _composed_id = malloc(_composed_id_idx[_n_g_interface] * sizeof(int));
  }


  PDM_MPI_Bcast(_composed_id    , _composed_id_idx[_n_g_interface], PDM_MPI_INT, 0, dom_intrf->comm);


  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    free(neighbor_n        [i_part]);
    free(neighbor_idx      [i_part]);
    free(neighbor_desc     [i_part]);
    free(neighbor_interface[i_part]);
    free(neighbor_opp_n    [i_part]);
    free(neighbor_opp_idx  [i_part]);
    free(neighbor_opp_desc [i_part]);
  }
  free(neighbor_n        );
  free(neighbor_idx      );
  free(neighbor_desc     );
  free(neighbor_interface);
  free(neighbor_opp_n    );
  free(neighbor_opp_idx  );
  free(neighbor_opp_desc );

  free(n_entity_bound);
  free(n_tot_part_by_domain);


  *n_g_interface                = _n_g_interface;
  *all_composed_id_idx          = _composed_id_idx;
  *all_composed_id              = _composed_id;
  *all_composed_ln_to_gn_sorted = _composed_ln_to_gn_sorted;

}

void
PDM_part_domain_interface_translate
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind1,
 PDM_bound_type_t               interface_kind2,
 int                           *n_part,
 int                          **pn_entity1,
 int                          **pn_entity2,
 PDM_g_num_t                 ***entity2_ln_to_gn,
 int                         ***entity2_entity1_idx,
 int                         ***entity2_entity1
)
{
  PDM_UNUSED(dom_intrf);
  PDM_UNUSED(interface_kind1);
  PDM_UNUSED(interface_kind2);
  PDM_UNUSED(n_part);
  PDM_UNUSED(pn_entity1);
  PDM_UNUSED(pn_entity2);
  PDM_UNUSED(entity2_ln_to_gn);
  PDM_UNUSED(entity2_entity1_idx);
  PDM_UNUSED(entity2_entity1);

  int         **pdi_neighbor_idx         = NULL;
  int         **pdi_neighbor             = NULL;
  int           n_composed_interface     = 0;
  int          *composed_interface_idx   = NULL;
  int          *composed_interface       = NULL;
  PDM_g_num_t  *composed_ln_to_gn_sorted = NULL;

  PDM_part_domain_interface_as_graph(dom_intrf,
                                     interface_kind1,
                                     pn_entity1,
                                     NULL,
                                     &pdi_neighbor_idx,
                                     &pdi_neighbor,
                                     &n_composed_interface,
                                     &composed_interface_idx,
                                     &composed_interface,
                                     &composed_ln_to_gn_sorted);
  free(composed_interface_idx);
  free(composed_interface);
  free(composed_ln_to_gn_sorted);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    n_part_loc_all_domain += n_part[i_dom];
  }

  int          *n_entity1          = (int          *) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int         **neighbor_interface = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int         **neighbor_idx       = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );
  int         **neighbor_desc      = (int         **) malloc( n_part_loc_all_domain * sizeof(int         *) );

  int shift_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      neighbor_idx      [i_part] = pdi_neighbor_idx[i_part];
      neighbor_desc     [i_part] = malloc( 3 * (neighbor_idx [i_part][pn_entity1[i_dom][i_part]]) * sizeof(int));
      neighbor_interface[i_part] = malloc(     (neighbor_idx [i_part][pn_entity1[i_dom][i_part]]) * sizeof(int));

      /* Copy */
      for(int i = 0; i < neighbor_idx [i_part][pn_entity1[i_dom][i_part]]; ++i) {
        neighbor_desc     [i_part][3*i  ] = pdi_neighbor[i_part][4*i  ];
        neighbor_desc     [i_part][3*i+1] = pdi_neighbor[i_part][4*i+1];
        neighbor_desc     [i_part][3*i+2] = pdi_neighbor[i_part][4*i+2];
        neighbor_interface[i_part][  i  ] = pdi_neighbor[i_part][4*i+3];
      }
      PDM_log_trace_graph_nuplet_int(neighbor_idx[i_part], neighbor_desc[i_part], 3, pn_entity1[i_dom][i_part], "neighbor_desc (debug) :");
      free(pdi_neighbor[i_part]);

      n_entity1[i_part+shift_part] = pn_entity1[i_dom][i_part];

    }
    shift_part += n_part[i_dom];
  }
  free(pdi_neighbor_idx);
  free(pdi_neighbor);

  /*
   * Prepare exchange by transform with local interface information - We change frame
   */
  int** entity1_is_dom_intrf = malloc(n_part_loc_all_domain * sizeof(int *));
  int** entity2_is_dom_intrf = malloc(n_part_loc_all_domain * sizeof(int *));

  shift_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      int spart = shift_part + i_part;
      entity1_is_dom_intrf[spart] = malloc(n_entity1        [spart] * sizeof(int));
      entity2_is_dom_intrf[spart] = malloc(pn_entity2[i_dom][i_part] * sizeof(int));

      for(int i = 0; i < n_entity1[spart]; ++i) {
        entity1_is_dom_intrf[spart][i] = 0;
      }
      for(int i = 0; i < pn_entity2[i_dom][i_part]; ++i) {
        entity2_is_dom_intrf[spart][i] = 0;
      }

      int *_entity1_is_dom_intrf = entity1_is_dom_intrf[spart];
      int *_entity2_is_dom_intrf = entity2_is_dom_intrf[spart];

      int *_entity2_entity1_idx = entity2_entity1_idx[i_dom][i_part];
      int *_entity2_entity1     = entity2_entity1    [i_dom][i_part];

      /* Loop over graph to flag all entity1 that have an interface */
      for(int i_entity1 = 0; i_entity1 < pn_entity1[i_dom][i_part]; ++i_entity1) {
        for(int idx_neight = neighbor_idx[i_part][i_entity1]; idx_neight < neighbor_idx[i_part][i_entity1+1]; ++idx_neight) {
          _entity1_is_dom_intrf[i_entity1] += 1;
        }
      }

      PDM_log_trace_array_int(_entity1_is_dom_intrf, n_entity1[spart], "_entity1_is_dom_intrf ::");

    }
    shift_part += n_part[i_dom];
  }


  /*
   *  Create distant neighbor for exchange between entity1
   */
  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(dom_intrf->comm,
                                                           n_part_loc_all_domain,
                                                           n_entity1,
                                                           neighbor_idx,
                                                           neighbor_desc);




  PDM_distant_neighbor_free(dn);



  /*
   *  Compute the transpose connectivity to post-treat exch
   */


  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    free(neighbor_interface[i_part]);
    free(neighbor_idx      [i_part]);
    free(neighbor_desc     [i_part]);
  }
  free(neighbor_interface);
  free(neighbor_idx      );
  free(neighbor_desc     );
  free(n_entity1);



}














