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
#include "pdm_part_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_distant_neighbor.h"
#include "pdm_order.h"

#include "pdm_domain_interface.h"
#include "pdm_domain_interface_priv.h"
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

  int *_unique_neighbor_entity_idx = NULL;
  int *_unique_neighbor_entity_n   = NULL;
  int *_unique_neighbor_entity     = NULL;
  int *order                       = NULL;
  PDM_malloc(_unique_neighbor_entity_idx, n_entity + 1                     , int);
  PDM_malloc(_unique_neighbor_entity_n  , n_entity                         , int);
  PDM_malloc(_unique_neighbor_entity    , 4 * neighbor_entity_idx[n_entity], int);
  PDM_malloc(order                      ,     neighbor_entity_idx[n_entity], int); // Suralloc

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

  PDM_realloc(_unique_neighbor_entity ,_unique_neighbor_entity , 4 * neighbor_entity_idx[n_entity] ,int);

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  PDM_free(order);
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

  int **_concat_neighbor_idx  = NULL;
  int **_concat_neighbor_desc = NULL;
  PDM_malloc(_concat_neighbor_idx , n_part, int *);
  PDM_malloc(_concat_neighbor_desc, n_part, int *);

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
    PDM_malloc(_concat_neighbor_desc[i_part], 4 * n_concat_idx_tot, int);

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

    PDM_free(_concat_neighbor_idx [i_part]);
    PDM_free(_concat_neighbor_desc[i_part]);
    PDM_free(_unique_concat_neighbor_n);
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
    int **prev_neighbor_n;
    PDM_malloc(prev_neighbor_n, n_part, int *);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_malloc(prev_neighbor_n[i_part], n_entity[i_part], int);
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
      PDM_free(prev_neighbor_n[i_part]);
    }
    PDM_free(prev_neighbor_n);

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
        PDM_free(prev_neighbor_idx  [i_part]);
        PDM_free(prev_neighbor_desc [i_part]);
      }
      PDM_free(prev_neighbor_idx);
      PDM_free(prev_neighbor_desc);
    }
    prev_neighbor_idx  = concat_neighbor_opp_idx;
    prev_neighbor_desc = concat_neighbor_opp;


    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_free(next_neighbor_opp_n[i_part]);
      PDM_free(next_neighbor_opp  [i_part]);
    }
    PDM_free(next_neighbor_opp_n);
    PDM_free(next_neighbor_opp);

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

  int **filter_neighbor_idx  = NULL;
  int **filter_neighbor_desc = NULL;
  PDM_malloc(filter_neighbor_idx , n_part, int *);
  PDM_malloc(filter_neighbor_desc, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *_concat_neighbor_opp_idx = concat_neighbor_opp_idx[i_part];
    int *_concat_neighbor_opp     = concat_neighbor_opp    [i_part];

    int *_neighbor_idx   = neighbor_idx      [i_part];
    int *_neighbor_desc  = neighbor_desc     [i_part];
    int *_neighbor_intrf = neighbor_interface[i_part];

    PDM_malloc(filter_neighbor_idx [i_part],     (n_entity[i_part] + 1)                      , int);
    PDM_malloc(filter_neighbor_desc[i_part], 4 * (_concat_neighbor_opp_idx[n_entity[i_part]]), int);

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
    PDM_realloc(filter_neighbor_desc[i_part] ,filter_neighbor_desc[i_part] , 4 * (_filter_neighbor_idx[n_entity[i_part]]) ,int);
    PDM_free(_concat_neighbor_opp_idx);
    PDM_free(_concat_neighbor_opp);

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(_filter_neighbor_idx, filter_neighbor_desc[i_part], 4, n_entity[i_part], "filter_neighbor_desc OOOO :");
    }


  }

  PDM_free(concat_neighbor_opp_idx);
  PDM_free(concat_neighbor_opp    );

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

  PDM_part_domain_interface_t *dom_intrf;
  PDM_malloc(dom_intrf, 1, PDM_part_domain_interface_t);
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  PDM_malloc(dom_intrf->n_part, n_domain, int);
  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->n_part[i_domain] = n_part[i_domain];
  }
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  PDM_malloc(dom_intrf->interface_pn_vtx      , n_domain, int          **);
  PDM_malloc(dom_intrf->interface_vtx_ln_to_gn, n_domain, PDM_g_num_t ***);
  PDM_malloc(dom_intrf->interface_sgn_vtx     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_sens_vtx    , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_vtx     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_vtx_idx , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_dom_vtx     , n_domain, int         ***);

  PDM_malloc(dom_intrf->interface_pn_edge      , n_domain, int          **);
  PDM_malloc(dom_intrf->interface_edge_ln_to_gn, n_domain, PDM_g_num_t ***);
  PDM_malloc(dom_intrf->interface_sgn_edge     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_sens_edge    , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_edge     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_edge_idx , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_dom_edge     , n_domain, int         ***);

  PDM_malloc(dom_intrf->interface_pn_face      , n_domain, int          **);
  PDM_malloc(dom_intrf->interface_face_ln_to_gn, n_domain, PDM_g_num_t ***);
  PDM_malloc(dom_intrf->interface_sgn_face     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_sens_face    , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_face     , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_ids_face_idx , n_domain, int         ***);
  PDM_malloc(dom_intrf->interface_dom_face     , n_domain, int         ***);

  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    PDM_malloc(dom_intrf->interface_pn_vtx      [i_domain], n_part[i_domain], int          *);
    PDM_malloc(dom_intrf->interface_vtx_ln_to_gn[i_domain], n_part[i_domain], PDM_g_num_t **);
    PDM_malloc(dom_intrf->interface_sgn_vtx     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_sens_vtx    [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_vtx     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_vtx_idx [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_dom_vtx     [i_domain], n_part[i_domain], int         **);

    PDM_malloc(dom_intrf->interface_pn_edge      [i_domain], n_part[i_domain], int          *);
    PDM_malloc(dom_intrf->interface_edge_ln_to_gn[i_domain], n_part[i_domain], PDM_g_num_t **);
    PDM_malloc(dom_intrf->interface_sgn_edge     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_sens_edge    [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_edge     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_edge_idx [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_dom_edge     [i_domain], n_part[i_domain], int         **);

    PDM_malloc(dom_intrf->interface_pn_face      [i_domain], n_part[i_domain], int          *);
    PDM_malloc(dom_intrf->interface_face_ln_to_gn[i_domain], n_part[i_domain], PDM_g_num_t **);
    PDM_malloc(dom_intrf->interface_sgn_face     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_sens_face    [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_face     [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_ids_face_idx [i_domain], n_part[i_domain], int         **);
    PDM_malloc(dom_intrf->interface_dom_face     [i_domain], n_part[i_domain], int         **);

    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      PDM_malloc(dom_intrf->interface_pn_vtx      [i_domain][i_part], n_interface, int          );
      PDM_malloc(dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part], n_interface, PDM_g_num_t *);
      PDM_malloc(dom_intrf->interface_sgn_vtx     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_sens_vtx    [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_vtx     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_vtx_idx [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_dom_vtx     [i_domain][i_part], n_interface, int         *);

      PDM_malloc(dom_intrf->interface_pn_edge      [i_domain][i_part], n_interface, int          );
      PDM_malloc(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part], n_interface, PDM_g_num_t *);
      PDM_malloc(dom_intrf->interface_sgn_edge     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_sens_edge    [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_edge     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_edge_idx [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_dom_edge     [i_domain][i_part], n_interface, int         *);

      PDM_malloc(dom_intrf->interface_pn_face      [i_domain][i_part], n_interface, int          );
      PDM_malloc(dom_intrf->interface_face_ln_to_gn[i_domain][i_part], n_interface, PDM_g_num_t *);
      PDM_malloc(dom_intrf->interface_sgn_face     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_sens_face    [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_face     [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_ids_face_idx [i_domain][i_part], n_interface, int         *);
      PDM_malloc(dom_intrf->interface_dom_face     [i_domain][i_part], n_interface, int         *);

      for(int i_interf = 0; i_interf < n_interface; ++i_interf) {
        dom_intrf->interface_pn_vtx       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sens_vtx     [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interf] = NULL;

        dom_intrf->interface_pn_edge       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_edge_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_edge      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sens_edge     [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_edge      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_edge_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_edge      [i_domain][i_part][i_interf] = NULL;

        dom_intrf->interface_pn_face       [i_domain][i_part][i_interf] = 0;
        dom_intrf->interface_face_ln_to_gn [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sgn_face      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_sens_face     [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_face      [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_ids_face_idx  [i_domain][i_part][i_interf] = NULL;
        dom_intrf->interface_dom_face      [i_domain][i_part][i_interf] = NULL;
      }
    }
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    dom_intrf->is_result[i] = 0;
  }

  PDM_malloc(dom_intrf->translation_vect  , n_interface, double *);
  PDM_malloc(dom_intrf->rotation_direction, n_interface, double *);
  PDM_malloc(dom_intrf->rotation_center   , n_interface, double *);
  PDM_malloc(dom_intrf->rotation_angle    , n_interface, double  );

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
 int                          *interface_sens,
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
    dom_intrf->interface_sens_vtx    [i_domain][i_part][i_interface] = interface_sens;
    dom_intrf->interface_ids_vtx     [i_domain][i_part][i_interface] = interface_ids;
    dom_intrf->interface_ids_vtx_idx [i_domain][i_part][i_interface] = interface_ids_idx;
    dom_intrf->interface_dom_vtx     [i_domain][i_part][i_interface] = interface_dom;
    dom_intrf->interface_describe_by_vtx  = 1;
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    dom_intrf->interface_pn_edge      [i_domain][i_part][i_interface] = interface_pn;
    dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface] = interface_ln_to_gn;
    dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface] = interface_sgn;
    dom_intrf->interface_sens_edge    [i_domain][i_part][i_interface] = interface_sens;
    dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface] = interface_ids;
    dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface] = interface_ids_idx;
    dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface] = interface_dom;
    dom_intrf->interface_describe_by_edge  = 1;
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_pn_face      [i_domain][i_part][i_interface] = interface_pn;
    dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface] = interface_ln_to_gn;
    dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface] = interface_sgn;
    dom_intrf->interface_sens_face    [i_domain][i_part][i_interface] = interface_sens;
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
 int                          **interface_sens,
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
    *interface_sens     = dom_intrf->interface_sens_vtx    [i_domain][i_part][i_interface];
    *interface_ids      = dom_intrf->interface_ids_vtx     [i_domain][i_part][i_interface];
    *interface_ids_idx  = dom_intrf->interface_ids_vtx_idx [i_domain][i_part][i_interface];
    *interface_dom      = dom_intrf->interface_dom_vtx     [i_domain][i_part][i_interface];
  } else if (interface_kind == PDM_BOUND_TYPE_EDGE) {
    *interface_pn       = dom_intrf->interface_pn_edge      [i_domain][i_part][i_interface];
    *interface_ln_to_gn = dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface];
    *interface_sgn      = dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface];
    *interface_sens     = dom_intrf->interface_sens_edge    [i_domain][i_part][i_interface];
    *interface_ids      = dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface];
    *interface_ids_idx  = dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface];
    *interface_dom      = dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface];
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    *interface_pn       = dom_intrf->interface_pn_face      [i_domain][i_part][i_interface];
    *interface_ln_to_gn = dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface];
    *interface_sgn      = dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface];
    *interface_sens     = dom_intrf->interface_sens_face    [i_domain][i_part][i_interface];
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
          PDM_free(dom_intrf->interface_pn_vtx       [i_domain][i_part]);
          dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part]);
          dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sgn_vtx      [i_domain][i_part]);
          dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sens_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sens_vtx      [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sens_vtx      [i_domain][i_part]);
          dom_intrf->interface_sens_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_vtx      [i_domain][i_part]);
          dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part]);
          dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_dom_vtx      [i_domain][i_part]);
          dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
        };

        if(dom_intrf->interface_pn_edge      [i_domain][i_part] != NULL) {
          PDM_free(dom_intrf->interface_pn_edge      [i_domain][i_part]);
          dom_intrf->interface_pn_edge      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_edge_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_edge_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sgn_edge     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sgn_edge     [i_domain][i_part]);
          dom_intrf->interface_sgn_edge     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sens_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sens_edge     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sens_edge     [i_domain][i_part]);
          dom_intrf->interface_sens_edge     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_edge     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_edge     [i_domain][i_part]);
          dom_intrf->interface_ids_edge     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_edge_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_edge_idx [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_edge_idx [i_domain][i_part]);
          dom_intrf->interface_ids_edge_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_edge     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_dom_edge     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_dom_edge     [i_domain][i_part]);
          dom_intrf->interface_dom_edge     [i_domain][i_part] = NULL;
        };

        if(dom_intrf->interface_pn_face      [i_domain][i_part] != NULL) {
          PDM_free(dom_intrf->interface_pn_face      [i_domain][i_part]);
          dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_face_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sgn_face     [i_domain][i_part]);
          dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sens_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_sens_face     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_sens_face     [i_domain][i_part]);
          dom_intrf->interface_sens_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_face     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_face     [i_domain][i_part]);
          dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_ids_face_idx [i_domain][i_part]);
          dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            PDM_free(dom_intrf->interface_dom_face     [i_domain][i_part][i_interface]);
          }
          PDM_free(dom_intrf->interface_dom_face     [i_domain][i_part]);
          dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
        };



      }
    }

    PDM_free(dom_intrf->interface_pn_vtx       [i_domain]);
    PDM_free(dom_intrf->interface_vtx_ln_to_gn [i_domain]);
    PDM_free(dom_intrf->interface_sgn_vtx      [i_domain]);
    PDM_free(dom_intrf->interface_sens_vtx     [i_domain]);
    PDM_free(dom_intrf->interface_ids_vtx      [i_domain]);
    PDM_free(dom_intrf->interface_ids_vtx_idx  [i_domain]);
    PDM_free(dom_intrf->interface_dom_vtx      [i_domain]);

    PDM_free(dom_intrf->interface_pn_edge      [i_domain]);
    PDM_free(dom_intrf->interface_edge_ln_to_gn[i_domain]);
    PDM_free(dom_intrf->interface_sgn_edge     [i_domain]);
    PDM_free(dom_intrf->interface_sens_edge    [i_domain]);
    PDM_free(dom_intrf->interface_ids_edge     [i_domain]);
    PDM_free(dom_intrf->interface_ids_edge_idx [i_domain]);
    PDM_free(dom_intrf->interface_dom_edge     [i_domain]);

    PDM_free(dom_intrf->interface_pn_face      [i_domain]);
    PDM_free(dom_intrf->interface_face_ln_to_gn[i_domain]);
    PDM_free(dom_intrf->interface_sgn_face     [i_domain]);
    PDM_free(dom_intrf->interface_sens_face    [i_domain]);
    PDM_free(dom_intrf->interface_ids_face     [i_domain]);
    PDM_free(dom_intrf->interface_ids_face_idx [i_domain]);
    PDM_free(dom_intrf->interface_dom_face     [i_domain]);

  }
  PDM_free(dom_intrf->n_part);

  PDM_free(dom_intrf->interface_pn_vtx       );
  PDM_free(dom_intrf->interface_vtx_ln_to_gn );
  PDM_free(dom_intrf->interface_sgn_vtx      );
  PDM_free(dom_intrf->interface_sens_vtx     );
  PDM_free(dom_intrf->interface_ids_vtx      );
  PDM_free(dom_intrf->interface_ids_vtx_idx  );
  PDM_free(dom_intrf->interface_dom_vtx      );

  PDM_free(dom_intrf->interface_pn_edge       );
  PDM_free(dom_intrf->interface_edge_ln_to_gn );
  PDM_free(dom_intrf->interface_sgn_edge      );
  PDM_free(dom_intrf->interface_sens_edge     );
  PDM_free(dom_intrf->interface_ids_edge      );
  PDM_free(dom_intrf->interface_ids_edge_idx  );
  PDM_free(dom_intrf->interface_dom_edge      );

  PDM_free(dom_intrf->interface_pn_face      );
  PDM_free(dom_intrf->interface_face_ln_to_gn);
  PDM_free(dom_intrf->interface_sgn_face     );
  PDM_free(dom_intrf->interface_sens_face    );
  PDM_free(dom_intrf->interface_ids_face     );
  PDM_free(dom_intrf->interface_ids_face_idx );
  PDM_free(dom_intrf->interface_dom_face     );

  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    if(dom_intrf->translation_vect[i_interface]   != NULL) {
      PDM_free(dom_intrf->translation_vect[i_interface]);
      dom_intrf->translation_vect[i_interface] = NULL;
    }
    if(dom_intrf->rotation_direction[i_interface]   != NULL) {
      PDM_free(dom_intrf->rotation_direction[i_interface]);
      dom_intrf->rotation_direction[i_interface] = NULL;
    }
    if(dom_intrf->rotation_center[i_interface]   != NULL) {
      PDM_free(dom_intrf->rotation_center[i_interface]);
      dom_intrf->rotation_center[i_interface] = NULL;
    }
  }

  PDM_free(dom_intrf->translation_vect  );
  PDM_free(dom_intrf->rotation_direction);
  PDM_free(dom_intrf->rotation_center   );
  PDM_free(dom_intrf->rotation_angle    );

  PDM_free(dom_intrf);
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

  PDM_malloc(dom_intrf->translation_vect[i_interface], 3, double);

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

  PDM_malloc(dom_intrf->rotation_direction[i_interface], 3, double);
  PDM_malloc(dom_intrf->rotation_center   [i_interface], 3, double);

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

    PDM_malloc(*vect, 3 , double);
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

    PDM_malloc(*direction, 3, double);
    PDM_malloc(*center   , 3, double);
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

  int *n_tot_part_by_domain;
  PDM_malloc(n_tot_part_by_domain, dom_intrf->n_domain, int);
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

  int **neighbor_n         = NULL;
  int **neighbor_idx       = NULL;
  int **neighbor_desc      = NULL;
  int **neighbor_interface = NULL;
  int  *n_entity_bound     = NULL;
  PDM_malloc(neighbor_n        , n_part_loc_all_domain, int *);
  PDM_malloc(neighbor_idx      , n_part_loc_all_domain, int *);
  PDM_malloc(neighbor_desc     , n_part_loc_all_domain, int *);
  PDM_malloc(neighbor_interface, n_part_loc_all_domain, int *);
  PDM_malloc(n_entity_bound    , n_part_loc_all_domain, int  );

  int **neighbor_opp_n    = NULL;
  int **neighbor_opp_idx  = NULL;
  int **neighbor_opp_desc = NULL;
  PDM_malloc(neighbor_opp_n   , n_part_loc_all_domain, int *);
  PDM_malloc(neighbor_opp_idx , n_part_loc_all_domain, int *);
  PDM_malloc(neighbor_opp_desc, n_part_loc_all_domain, int *);

  /*
   * Loop over all interfaces to create distant neighbor structure
   */
  int shift_part   = 0;
  int shift_part_g = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {
    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      n_entity_bound[i_part+shift_part] = n_entity[i_domain][i_part];
      PDM_malloc(neighbor_idx  [i_part+shift_part], n_entity_bound[i_part+shift_part] + 1, int);
      neighbor_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_n   = neighbor_n    [i_part+shift_part];
      int* _neighbor_idx = neighbor_idx  [i_part+shift_part];

      PDM_malloc(neighbor_opp_idx  [i_part+shift_part], n_entity_bound[i_part+shift_part] + 1, int);
      neighbor_opp_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_opp_n   = neighbor_opp_n    [i_part+shift_part];
      int* _neighbor_opp_idx = neighbor_opp_idx  [i_part+shift_part];

      int          *interface_pn       = NULL;
      PDM_g_num_t **interface_ln_to_gn = NULL;
      int         **interface_sgn      = NULL;
      int         **interface_sens     = NULL;
      int         **interface_ids      = NULL;
      int         **interface_ids_idx  = NULL;
      int         **interface_dom      = NULL;
      PDM_malloc(interface_pn      , n_interface, int          );
      PDM_malloc(interface_ln_to_gn, n_interface, PDM_g_num_t *);
      PDM_malloc(interface_sgn     , n_interface, int         *);
      PDM_malloc(interface_sens    , n_interface, int         *);
      PDM_malloc(interface_ids     , n_interface, int         *);
      PDM_malloc(interface_ids_idx , n_interface, int         *);
      PDM_malloc(interface_dom     , n_interface, int         *);
      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
        PDM_part_domain_interface_get(dom_intrf,
                                      interface_kind,
                                      i_domain,
                                      i_part,
                                      i_interface,
                                      &interface_pn      [i_interface],
                                      &interface_ln_to_gn[i_interface],
                                      &interface_sgn     [i_interface],
                                      &interface_sens    [i_interface],
                                      &interface_ids     [i_interface],
                                      &interface_ids_idx [i_interface],
                                      &interface_dom     [i_interface]);
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

      PDM_malloc(neighbor_desc     [i_part+shift_part], 3 * _neighbor_idx    [n_entity_bound[i_part+shift_part]], int);
      PDM_malloc(neighbor_opp_desc [i_part+shift_part], 4 * _neighbor_opp_idx[n_entity_bound[i_part+shift_part]], int);
      PDM_malloc(neighbor_interface[i_part+shift_part],     _neighbor_idx[n_entity_bound[i_part+shift_part]]    , int);
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

      PDM_free(interface_pn      );
      PDM_free(interface_ln_to_gn);
      PDM_free(interface_sgn     );
      PDM_free(interface_sens    );
      PDM_free(interface_ids     );
      PDM_free(interface_ids_idx );
      PDM_free(interface_dom     );

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
  int *sgn_interf_to_interf = NULL;
  PDM_malloc(sgn_interf_to_interf, 2 * n_interface, int);

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

  int         *composed_id_idx         = NULL;
  int         *composed_id             = NULL;
  int         *composed_id_tmp         = NULL;
  PDM_g_num_t *composed_key            = NULL;
  int         *composed_key_update_idx = NULL;
  PDM_malloc(composed_id_idx        ,    (max_composed+1), int        );
  PDM_malloc(composed_id            ,     max_composed   , int        );
  PDM_malloc(composed_id_tmp        ,     max_lcomposed  , int        );
  PDM_malloc(composed_key           ,     max_composed   , PDM_g_num_t);
  PDM_malloc(composed_key_update_idx, 2 * max_composed   , int        );

  composed_id_idx[0] = 0;
  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {

    int n_elmt = n_entity_bound[i_part];
    int *_neighbor_entity_idx  = (*neighbor_entity_idx )[i_part];
    int *_neighbor_entity_desc = (*neighbor_entity_desc)[i_part];

    int *_filter_neighbor_entity_idx;
    PDM_malloc(_filter_neighbor_entity_idx, n_elmt + 1, int);

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

    PDM_free(_neighbor_entity_idx);
    (*neighbor_entity_idx)[i_part] = _filter_neighbor_entity_idx;

    if(0 == 1) {
      PDM_log_trace_graph_nuplet_int(_filter_neighbor_entity_idx, _neighbor_entity_desc, 4, n_elmt, "_neighbor_entity_desc :");
    }
  }

  /*
   * Realloc
   */
  PDM_free(composed_id_tmp);
  PDM_realloc(composed_id_idx ,composed_id_idx , (i_composed_interface+1)              ,int        );
  PDM_realloc(composed_id     ,composed_id     , composed_id_idx[i_composed_interface] ,int        );
  PDM_realloc(composed_key    ,composed_key    , (i_composed_interface+1)              ,PDM_g_num_t);

  /*
   * Generate table to give from interface number the composed one
   *   - We need the g_id of new interface and the associate composition
   */
  PDM_g_num_t* composed_id_gnum = composed_key;

  if(0 == 1) {
    PDM_log_trace_array_long(composed_key    , i_composed_interface, "composed_id :: ");
    PDM_log_trace_array_long(composed_id_gnum, i_composed_interface, "composed_id_gnum :: ");
    PDM_log_trace_connectivity_int(composed_id_idx, composed_id, i_composed_interface, "composed_id :: ");
  }

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
  PDM_free(composed_key_update_idx);

  int n_rank = -1;
  PDM_MPI_Comm_size(dom_intrf->comm, &n_rank);

  int *composed_id_n;
  PDM_malloc(composed_id_n, i_composed_interface, int        );
  PDM_g_num_t max_loc = 0;
  for(int i = 0; i < i_composed_interface; ++i) {
    max_loc = PDM_MAX(max_loc, composed_id_gnum[i]);
    composed_id_n[i] = composed_id_idx[i+1] - composed_id_idx[i];
  }
  PDM_free(composed_id_idx);
  PDM_g_num_t max_glob = -1;
  PDM_MPI_Allreduce(&max_loc, &max_glob, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, dom_intrf->comm);

  PDM_g_num_t *distrib_interf;
  PDM_malloc(distrib_interf, n_rank + 1, PDM_g_num_t);
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

  //PDM_free(composed_id_gnum);
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

  PDM_free(composed_id);
  PDM_free(composed_id_n);

  int _n_g_interface = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *ptb_composed_ln_to_gn_sorted = PDM_part_to_block_block_gnum_get(ptb);

  if(i_rank > 0) {
    assert(_n_g_interface == 0);
    PDM_free(_composed_id);
  }


  PDM_free(distrib_interf);

  PDM_free(composed_key);
  PDM_free(sgn_interf_to_interf);


  /*
   *  Broadcast to all
   */
  PDM_MPI_Bcast(&_n_g_interface  , 1, PDM_MPI_INT, 0, dom_intrf->comm);

  PDM_g_num_t *_composed_ln_to_gn_sorted;
  PDM_malloc(_composed_ln_to_gn_sorted, _n_g_interface, PDM_g_num_t);

  if(i_rank == 0) {
    for(int i = 0; i < _n_g_interface; ++i) {
      _composed_ln_to_gn_sorted[i] =  ptb_composed_ln_to_gn_sorted[i];
    }
  }

  int *_composed_id_idx;
  PDM_malloc(_composed_id_idx, _n_g_interface + 1, int);

  if(i_rank == 0) {
    _composed_id_idx[0] = 0;
    for(int i = 0; i < _n_g_interface; ++i) {
      _composed_id_idx[i+1] = _composed_id_idx[i] + _composed_id_n[i];
    }
  }
  PDM_free(_composed_id_n);

  PDM_part_to_block_free(ptb);

  PDM_MPI_Bcast(_composed_ln_to_gn_sorted    , _n_g_interface, PDM__PDM_MPI_G_NUM, 0, dom_intrf->comm);

  PDM_MPI_Bcast(_composed_id_idx, _n_g_interface+1                , PDM_MPI_INT, 0, dom_intrf->comm);


  if(i_rank != 0) {
    PDM_malloc(_composed_id, _composed_id_idx[_n_g_interface], int);
  }


  PDM_MPI_Bcast(_composed_id    , _composed_id_idx[_n_g_interface], PDM_MPI_INT, 0, dom_intrf->comm);


  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    PDM_free(neighbor_n        [i_part]);
    PDM_free(neighbor_idx      [i_part]);
    PDM_free(neighbor_desc     [i_part]);
    PDM_free(neighbor_interface[i_part]);
    PDM_free(neighbor_opp_n    [i_part]);
    PDM_free(neighbor_opp_idx  [i_part]);
    PDM_free(neighbor_opp_desc [i_part]);
  }
  PDM_free(neighbor_n        );
  PDM_free(neighbor_idx      );
  PDM_free(neighbor_desc     );
  PDM_free(neighbor_interface);
  PDM_free(neighbor_opp_n    );
  PDM_free(neighbor_opp_idx  );
  PDM_free(neighbor_opp_desc );

  PDM_free(n_entity_bound);
  PDM_free(n_tot_part_by_domain);


  *n_g_interface                = _n_g_interface;
  *all_composed_id_idx          = _composed_id_idx;
  *all_composed_id              = _composed_id;
  *all_composed_ln_to_gn_sorted = _composed_ln_to_gn_sorted;

}

void
PDM_part_domain_interface_view_by_part
(
  PDM_part_domain_interface_t   *pdi,
  PDM_bound_type_t               interface_kind,
  int                           *pn_entity,
  PDM_g_num_t                  **pentity_ln_to_gn,
  int                          **pn_entity_num_out,
  int                         ***pentity_num_out,
  int                         ***pentity_opp_location_idx_out,
  int                         ***pentity_opp_location_out,
  int                         ***pentity_opp_interface_idx_out,
  int                         ***pentity_opp_interface_out,
  int                         ***pentity_opp_sens_out,
  PDM_g_num_t                 ***pentity_opp_gnum_out
)
{

  int i_rank;
  PDM_MPI_Comm_rank(pdi->comm, &i_rank);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < pdi->n_domain; ++i_dom) {
    n_part_loc_all_domain += pdi->n_part[i_dom];
  }

  int  *pn_entity_num              = NULL;
  int **pentity_num                = NULL;
  int **pentity_opp_location_idx   = NULL;
  int **pentity_opp_location       = NULL;
  int **pentity_opp_interface_idx  = NULL;
  int **pentity_opp_interface      = NULL;
  int **pentity_opp_sens           = NULL;
  PDM_malloc(pn_entity_num            , n_part_loc_all_domain, int  );
  PDM_malloc(pentity_num              , n_part_loc_all_domain, int *);
  PDM_malloc(pentity_opp_location_idx , n_part_loc_all_domain, int *);
  PDM_malloc(pentity_opp_location     , n_part_loc_all_domain, int *);
  PDM_malloc(pentity_opp_interface_idx, n_part_loc_all_domain, int *);
  PDM_malloc(pentity_opp_interface    , n_part_loc_all_domain, int *);
  PDM_malloc(pentity_opp_sens         , n_part_loc_all_domain, int *);


  int s_part = 0;
  for(int i_dom = 0; i_dom < pdi->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < pdi->n_part[i_dom]; ++i_part) {

      pentity_num  [s_part+i_part] = PDM_array_zeros_int(pn_entity[s_part+i_part]);
      int *pentity_cur_n   = PDM_array_zeros_int(pn_entity[s_part+i_part]);
      int *_pentity_num = pentity_num[s_part+i_part];

      /*  */
      for(int i_interface = 0; i_interface < pdi->n_interface; ++i_interface) {

        int           ln_interface        = 0;
        PDM_g_num_t  *pinterface_ln_to_gn = NULL;
        int          *pinterface_sgn      = NULL;
        int          *pinterface_sens     = NULL;
        int          *pinterface_ids      = NULL;
        int          *pinterface_ids_idx  = NULL;
        int          *pinterface_dom      = NULL;

        PDM_part_domain_interface_get(pdi,
                                      interface_kind,
                                      i_dom,
                                      i_part,
                                      i_interface,
                                      &ln_interface,
                                      &pinterface_ln_to_gn,
                                      &pinterface_sgn,
                                      &pinterface_sens,
                                      &pinterface_ids,
                                      &pinterface_ids_idx,
                                      &pinterface_dom);

        if(0 == 1) {
          PDM_log_trace_array_int (pinterface_sgn     ,   ln_interface, "pinterface_sgn      ::");
          PDM_log_trace_array_int (pinterface_sens    ,   ln_interface, "pinterface_sens     ::");
          PDM_log_trace_array_int (pinterface_dom     , 2*ln_interface, "pinterface_dom      ::");
          PDM_log_trace_array_long(pinterface_ln_to_gn,   ln_interface, "pinterface_ln_to_gn ::");
          PDM_log_trace_array_int (pinterface_ids     , 3 * pinterface_ids_idx[ln_interface], "pinterface_ids ::");
        }

        /*
         * First step : Count interface to add in distant neighbor due to connectivity betwenn domain
         */
        for(int idx_entity = 0; idx_entity < ln_interface; ++idx_entity) {

          // Search the first in list that is in current part/proc
          int i_entity_cur = -1;
          int found        = 0;
          for(int j = pinterface_ids_idx[idx_entity]; j < pinterface_ids_idx[idx_entity+1]; ++j) {
            int i_proc_opp   = pinterface_ids[3*j  ];
            int i_part_opp   = pinterface_ids[3*j+1];
            int i_entity_opp = pinterface_ids[3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == s_part + i_part && found == 0) {
              i_entity_cur = i_entity_opp;
              found = 1;
              int n_opp = pinterface_ids_idx[idx_entity+1] - pinterface_ids_idx[idx_entity];
              pentity_cur_n[i_entity_cur] += (n_opp - 1);
            }
          }

          assert(found == 1);
        }
      } /* End i_interface */

      if(1 == 0) {
        PDM_log_trace_array_int (_pentity_num , pn_entity[s_part+i_part], "_pentity_num  ::");
        PDM_log_trace_array_int (pentity_cur_n, pn_entity[s_part+i_part], "pentity_cur_n ::");
      }

      /* Creation de l'indirection inverse */
      int *pentity_cur_idx = PDM_array_const_int(pn_entity[s_part+i_part], -1);

      pn_entity_num[s_part+i_part] = 0;
      for(int i = 0; i < pn_entity[s_part+i_part]; ++i ) {
        if(pentity_cur_n[i] > 0) {
          _pentity_num[pn_entity_num[s_part+i_part]] = i;
          pentity_cur_idx[i] = pn_entity_num[s_part+i_part];
          pn_entity_num[s_part+i_part]++;
        }
      }

      PDM_malloc(pentity_opp_location_idx[s_part+i_part], pn_entity_num[s_part+i_part]+1, int);
      int *_pentity_opp_location_idx = pentity_opp_location_idx[s_part+i_part];

      _pentity_opp_location_idx[0] = 0;
      for(int idx_entity = 0; idx_entity < pn_entity_num[s_part+i_part]; ++idx_entity) {
        int i_entity = _pentity_num[idx_entity];
        _pentity_opp_location_idx[idx_entity+1] = _pentity_opp_location_idx[idx_entity] + pentity_cur_n[i_entity];
        pentity_cur_n[i_entity] = 0;
      }

      int n_location = _pentity_opp_location_idx[pn_entity_num[s_part+i_part]];
      PDM_malloc(pentity_opp_location [s_part+i_part], 3 * n_location, int);
      PDM_malloc(pentity_opp_interface[s_part+i_part],     n_location, int);
      PDM_malloc(pentity_opp_sens     [s_part+i_part],     n_location, int);
      int *_pentity_opp_location  = pentity_opp_location [s_part+i_part];
      int *_pentity_opp_interface = pentity_opp_interface[s_part+i_part];
      int *_pentity_opp_sens      = pentity_opp_sens     [s_part+i_part];

      /* Fill */
      for(int i_interface = 0; i_interface < pdi->n_interface; ++i_interface) {

        int           ln_interface        = 0;
        PDM_g_num_t  *pinterface_ln_to_gn = NULL;
        int          *pinterface_sgn      = NULL;
        int          *pinterface_sens     = NULL;
        int          *pinterface_ids      = NULL;
        int          *pinterface_ids_idx  = NULL;
        int          *pinterface_dom      = NULL;

        PDM_part_domain_interface_get(pdi,
                                      interface_kind,
                                      i_dom,
                                      i_part,
                                      i_interface,
                                      &ln_interface,
                                      &pinterface_ln_to_gn,
                                      &pinterface_sgn,
                                      &pinterface_sens,
                                      &pinterface_ids,
                                      &pinterface_ids_idx,
                                      &pinterface_dom);

        for(int idx_entity = 0; idx_entity < ln_interface; ++idx_entity) {

          // Search the first in list that is in current part/proc
          int i_entity_cur = -1;
          int found        = 0;
          int idx_current  = -1;
          for(int j = pinterface_ids_idx[idx_entity]; j < pinterface_ids_idx[idx_entity+1]; ++j) {
            int i_proc_opp   = pinterface_ids[3*j  ];
            int i_part_opp   = pinterface_ids[3*j+1];
            int i_entity_opp = pinterface_ids[3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == s_part+i_part && found == 0) {
              i_entity_cur = i_entity_opp;
              found = 1;
              idx_current  = j;
            }
          }

          assert(found == 1);

          for(int j = pinterface_ids_idx[idx_entity]; j < pinterface_ids_idx[idx_entity+1]; ++j) {
            int i_proc_opp   = pinterface_ids[3*j  ];
            int i_part_opp   = pinterface_ids[3*j+1];
            int i_entity_opp = pinterface_ids[3*j+2];

            if(idx_current != j) {
              int idx_write = _pentity_opp_location_idx[pentity_cur_idx[i_entity_cur]] + pentity_cur_n[i_entity_cur]++;

              _pentity_opp_location [3*idx_write  ] = i_proc_opp;
              _pentity_opp_location [3*idx_write+1] = i_part_opp;
              _pentity_opp_location [3*idx_write+2] = i_entity_opp;
              _pentity_opp_interface[  idx_write  ] = (i_interface+1)*pinterface_sgn [idx_entity];
              _pentity_opp_sens     [  idx_write  ] = pinterface_sens[idx_entity];
            }
          }
        }
      } /* End i_interface */

      PDM_free(pentity_cur_n);
      PDM_free(pentity_cur_idx);

    }
    s_part += pdi->n_part[i_dom];
  }

  /*
   *  Le _pentity_opp_interface_idx apparaitra quand on fera de la combinaison
   */

  /*
   * Exchange to have all opposit gnum
   */

  int         **part1_to_part2_idx         = NULL;
  int         **part1_to_part2_triplet_idx = NULL;
  int         **part1_to_part2_triplet     = NULL;
  int         **part1_to_part2_interface   = NULL;
  PDM_g_num_t **part1_to_part2_gnum        = NULL;
  PDM_malloc(part1_to_part2_idx,         n_part_loc_all_domain, int         *);
  // PDM_malloc(part1_to_part2_triplet_idx, n_part_loc_all_domain, int         *);
  PDM_malloc(part1_to_part2_triplet,     n_part_loc_all_domain, int         *);
  PDM_malloc(part1_to_part2_interface,   n_part_loc_all_domain, int         *);
  PDM_malloc(part1_to_part2_gnum,        n_part_loc_all_domain, PDM_g_num_t *);

  int li_part = 0;
  for(int i_dom = 0; i_dom < pdi->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < pdi->n_part[i_dom]; ++i_part) {

      PDM_malloc(part1_to_part2_idx[li_part],(pn_entity[li_part] + 1) ,int);
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pn_entity[li_part]);

      for(int idx_entity = 0; idx_entity < pn_entity_num[li_part]; ++idx_entity) {
        int i_entity = pentity_num[li_part][idx_entity];
        int n_opp = pentity_opp_location_idx[li_part][idx_entity+1] - pentity_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }

      /* Setup idx */
      for(int i_entity = 0; i_entity < pn_entity[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pn_entity[li_part]];
      PDM_malloc(part1_to_part2_triplet  [li_part],n_connect_tot  , int        );
      PDM_malloc(part1_to_part2_interface[li_part],n_connect_tot/3, int        );
      PDM_malloc(part1_to_part2_gnum     [li_part],n_connect_tot/3, PDM_g_num_t);

      for(int i = 0; i < n_connect_tot; ++i) {
        part1_to_part2_triplet  [li_part][i] = -10000;
      }
      // PDM_log_trace_array_int(part1_to_part2_idx      [li_part], pn_entity[li_part], "part1_to_part2_idx       ::");

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity_num[li_part]; ++idx_entity) {
        int i_entity = pentity_num[li_part][idx_entity];
        for(int idx_opp = pentity_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pentity_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pentity_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pentity_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pentity_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pentity_opp_interface[li_part][  idx_opp  ];

          part1_to_part2_gnum     [li_part][idx_write  ] = pentity_ln_to_gn     [li_part][  i_entity ];

        }
      }

      // PDM_log_trace_array_int(part1_to_part2_triplet  [li_part], n_connect_tot     , "part1_to_part2_triplet   ::");
      // PDM_log_trace_array_int(part1_to_part2_interface[li_part], n_connect_tot/3   , "part1_to_part2_interface ::");

      PDM_free(part1_to_part2_n);

      li_part += 1;
    }
  }

  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity_ln_to_gn,
                                                                      (const int          *) pn_entity,
                                                                      n_part_loc_all_domain,
                                                                      (const int          *) pn_entity,
                                                                      n_part_loc_all_domain,
                                                                      (const int         **) part1_to_part2_idx,
                                                                      (const int         **) part1_to_part2_triplet_idx,
                                                                      (const int         **) part1_to_part2_triplet,
                                                                      pdi->comm);

  int exch_request = -1;
  PDM_g_num_t **pextract_gnum_opp = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                                 NULL,
                (const void **)  pentity_ln_to_gn,
                                 NULL,
                    (void ***)   &pextract_gnum_opp,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  if(0 == 1) { // Usefull to know how many data is transfer
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");
    }

    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "extract_lnum2 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_ref_lnum2  [i_part], "gnum1_come_from ::");
    }

    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      int n_send = part1_to_part2_idx[i_part][pn_entity[i_part]]/3;
      log_trace("n_send = %i \n", n_send);
      PDM_log_trace_array_long(pextract_gnum_opp[i_part], n_send, "pextract_gnum_opp      : ");
    }
  }


  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    PDM_free(part1_to_part2_gnum     [i_part]);
    PDM_free(part1_to_part2_idx      [i_part]);
    PDM_free(part1_to_part2_triplet  [i_part]);
    PDM_free(part1_to_part2_interface[i_part]);
  }
  PDM_free(part1_to_part2_gnum     );
  PDM_free(part1_to_part2_idx      );
  PDM_free(part1_to_part2_triplet  );
  PDM_free(part1_to_part2_interface);

  PDM_part_to_part_free(ptp);

  *pn_entity_num_out             = pn_entity_num;
  *pentity_num_out               = pentity_num;
  *pentity_opp_location_idx_out  = pentity_opp_location_idx;
  *pentity_opp_location_out      = pentity_opp_location;
  *pentity_opp_interface_idx_out = pentity_opp_interface_idx;
  *pentity_opp_interface_out     = pentity_opp_interface;
  *pentity_opp_sens_out          = pentity_opp_sens;
  *pentity_opp_gnum_out          = pextract_gnum_opp;
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
  PDM_free(composed_interface_idx);
  PDM_free(composed_interface);
  PDM_free(composed_ln_to_gn_sorted);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    n_part_loc_all_domain += n_part[i_dom];
  }

  int  *n_entity1          = NULL;
  int **neighbor_interface = NULL;
  int **neighbor_idx       = NULL;
  int **neighbor_desc      = NULL;
  PDM_malloc(n_entity1         , n_part_loc_all_domain, int          );
  PDM_malloc(neighbor_interface, n_part_loc_all_domain, int         *);
  PDM_malloc(neighbor_idx      , n_part_loc_all_domain, int         *);
  PDM_malloc(neighbor_desc     , n_part_loc_all_domain, int         *);

  int shift_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      neighbor_idx      [i_part] = pdi_neighbor_idx[i_part];
      PDM_malloc(neighbor_desc     [i_part], 3 * (neighbor_idx [i_part][pn_entity1[i_dom][i_part]]), int);
      PDM_malloc(neighbor_interface[i_part],     (neighbor_idx [i_part][pn_entity1[i_dom][i_part]]), int);

      /* Copy */
      for(int i = 0; i < neighbor_idx [i_part][pn_entity1[i_dom][i_part]]; ++i) {
        neighbor_desc     [i_part][3*i  ] = pdi_neighbor[i_part][4*i  ];
        neighbor_desc     [i_part][3*i+1] = pdi_neighbor[i_part][4*i+1];
        neighbor_desc     [i_part][3*i+2] = pdi_neighbor[i_part][4*i+2];
        neighbor_interface[i_part][  i  ] = pdi_neighbor[i_part][4*i+3];
      }
      PDM_log_trace_graph_nuplet_int(neighbor_idx[i_part], neighbor_desc[i_part], 3, pn_entity1[i_dom][i_part], "neighbor_desc (debug) :");
      PDM_free(pdi_neighbor[i_part]);

      n_entity1[i_part+shift_part] = pn_entity1[i_dom][i_part];

    }
    shift_part += n_part[i_dom];
  }
  PDM_free(pdi_neighbor_idx);
  PDM_free(pdi_neighbor);

  /*
   * Prepare exchange by transform with local interface information - We change frame
   */
  int* *entity1_is_dom_intrf = NULL;
  PDM_malloc(entity1_is_dom_intrf, n_part_loc_all_domain, int *);

  shift_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      int spart = shift_part + i_part;
      PDM_malloc(entity1_is_dom_intrf[spart], n_entity1[spart], int);
      // PDM_malloc(entity2_is_dom_intrf[spart],pn_entity2[i_dom][i_part] ,int);

      for(int i = 0; i < n_entity1[spart]; ++i) {
        entity1_is_dom_intrf[spart][i] = 0;
      }
      int *_entity1_is_dom_intrf = entity1_is_dom_intrf[spart];

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
    PDM_free(neighbor_interface[i_part]);
    PDM_free(neighbor_idx      [i_part]);
    PDM_free(neighbor_desc     [i_part]);
  }
  PDM_free(neighbor_interface);
  PDM_free(neighbor_idx      );
  PDM_free(neighbor_desc     );
  PDM_free(n_entity1);



}

void
PDM_part_domain_interface_to_domain_interface
(
  PDM_part_domain_interface_t    *dom_intrf,
  PDM_bound_type_t                interface_kind1,
  int                            *n_part,
  int                           **pn_entity1,
  PDM_g_num_t                  ***entity1_ln_to_gn,
  PDM_domain_interface_t        **ditrf_out,
  int                          ***is_entity1_on_itrf_out
)
{
  int debug = 0;

  int i_rank;
  PDM_MPI_Comm_rank(dom_intrf->comm, &i_rank);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    n_part_loc_all_domain += n_part[i_dom];
  }

  int n_interface = PDM_part_domain_interface_n_interface_get(dom_intrf);
  int          **pn_interface               = NULL;
  PDM_g_num_t ***interface_entity1_ln_to_gn = NULL;
  PDM_g_num_t ***interface_ln_to_gn         = NULL;
  int         ***interface_sgn              = NULL;
  int         ***interface_sens             = NULL;
  int         ***interface_dom              = NULL;
  PDM_malloc(pn_interface              , n_interface, int          *);
  PDM_malloc(interface_entity1_ln_to_gn, n_interface, PDM_g_num_t **);
  PDM_malloc(interface_ln_to_gn        , n_interface, PDM_g_num_t **);
  PDM_malloc(interface_sgn             , n_interface, int         **);
  PDM_malloc(interface_sens            , n_interface, int         **);
  PDM_malloc(interface_dom             , n_interface, int         **);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    PDM_malloc(pn_interface              [i_interface], n_part_loc_all_domain, int          );
    PDM_malloc(interface_entity1_ln_to_gn[i_interface], n_part_loc_all_domain, PDM_g_num_t *);
    PDM_malloc(interface_ln_to_gn        [i_interface], n_part_loc_all_domain, PDM_g_num_t *);
    PDM_malloc(interface_sgn             [i_interface], n_part_loc_all_domain, int         *);
    PDM_malloc(interface_sens            [i_interface], n_part_loc_all_domain, int         *);
    PDM_malloc(interface_dom             [i_interface], n_part_loc_all_domain, int         *);
  }

  /*
   * Re-Create a domain interface
   */
  PDM_domain_interface_t* ditrf = *ditrf_out;
  if(ditrf == NULL) {
    ditrf = PDM_domain_interface_create(n_interface,
                                        dom_intrf->n_domain,
                                        PDM_DOMAIN_INTERFACE_MULT_YES,
                                        PDM_OWNERSHIP_KEEP,
                                        dom_intrf->comm);
  }


  /*
   * Tag all entity concerns by interface
   */
  int **is_entity1_on_itrf = NULL;
  PDM_malloc(is_entity1_on_itrf, n_part_loc_all_domain, int *);

  /* Fix output */
  *is_entity1_on_itrf_out = is_entity1_on_itrf;
  *ditrf_out              = ditrf;

  int s_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      is_entity1_on_itrf[s_part+i_part] = PDM_array_zeros_int(pn_entity1[i_dom][i_part]);

      /*  */
      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

        int           ln_interface        = 0;
        PDM_g_num_t  *pinterface_ln_to_gn = NULL;
        int          *pinterface_sgn      = NULL;
        int          *pinterface_sens     = NULL;
        int          *pinterface_ids      = NULL;
        int          *pinterface_ids_idx  = NULL;
        int          *pinterface_dom      = NULL;

        PDM_part_domain_interface_get(dom_intrf,
                                      interface_kind1,
                                      i_dom,
                                      i_part,
                                      i_interface,
                                      &ln_interface,
                                      &pinterface_ln_to_gn,
                                      &pinterface_sgn,
                                      &pinterface_sens,
                                      &pinterface_ids,
                                      &pinterface_ids_idx,
                                      &pinterface_dom);

        if(debug == 1) {
          PDM_log_trace_array_long(pinterface_ln_to_gn,   ln_interface, "pinterface_ln_to_gn ::");
          PDM_log_trace_array_int (pinterface_sgn     ,   ln_interface, "pinterface_sgn      ::");
          PDM_log_trace_array_int (pinterface_sens    ,   ln_interface, "pinterface_sens     ::");
          PDM_log_trace_array_int (pinterface_dom     , 2*ln_interface, "pinterface_dom      ::");
          PDM_log_trace_array_long(pinterface_ln_to_gn,   ln_interface, "pinterface_ln_to_gn ::");
          PDM_log_trace_graph_nuplet_int(pinterface_ids_idx, pinterface_ids, 3, ln_interface, "pinterface_ids ::");
        }

        PDM_malloc(interface_entity1_ln_to_gn[i_interface][s_part+i_part], ln_interface, PDM_g_num_t);
        PDM_malloc(interface_dom             [i_interface][s_part+i_part], ln_interface, int        );
        PDM_g_num_t* _interface_entity1_ln_to_gn = interface_entity1_ln_to_gn[i_interface][s_part+i_part];
        int        * _interface_dom              = interface_dom             [i_interface][s_part+i_part];

        int idx_write = 0;
        for(int idx_entity = 0; idx_entity < ln_interface; ++idx_entity) {

          int found        = 0;
          for(int j = pinterface_ids_idx[idx_entity]; j < pinterface_ids_idx[idx_entity+1]; ++j) {
            int i_proc_opp   = pinterface_ids[3*j  ];
            int i_part_opp   = pinterface_ids[3*j+1];
            int i_entity_opp = pinterface_ids[3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == s_part+i_part && found == 0) {
              found = 1;
              is_entity1_on_itrf         [s_part+i_part][i_entity_opp] = 1;
              _interface_entity1_ln_to_gn[idx_write] = entity1_ln_to_gn[i_dom][i_part][i_entity_opp];
              _interface_dom             [idx_write] = i_dom;
              idx_write++;
            }
          }
          assert(found == 1);
        }

        assert(idx_write == ln_interface);

        pn_interface      [i_interface][s_part+i_part] = ln_interface;
        interface_ln_to_gn[i_interface][s_part+i_part] = pinterface_ln_to_gn;
        interface_sgn     [i_interface][s_part+i_part] = pinterface_sgn;
        interface_sens    [i_interface][s_part+i_part] = pinterface_sens;

      }

      // PDM_log_trace_array_int(is_entity1_on_itrf[s_part+i_part], pn_entity1[i_dom][i_part], "is_entity1_on_itrf ::");
    }
    s_part += n_part[i_dom];
  }


  int          *dinterface_dn  = NULL;
  int         **dinterface_dom = NULL;
  PDM_g_num_t **dinterface_ids = NULL;
  PDM_malloc(dinterface_dn , n_interface, int           );
  PDM_malloc(dinterface_dom, n_interface, int         * );
  PDM_malloc(dinterface_ids, n_interface, PDM_g_num_t * );


  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        interface_ln_to_gn[i_interface],
                                                        NULL,
                                                        pn_interface[i_interface],
                                                        n_part_loc_all_domain,
                                                        dom_intrf->comm);

    int* *stride_one;
    PDM_malloc(stride_one, n_part_loc_all_domain, int *);
    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      stride_one[i_part] = PDM_array_const_int(pn_interface[i_interface][i_part], 1);
    }


    PDM_g_num_t *dentity1_gnum = NULL;
    int         *dentity1_sgn  = NULL;
    int         *dentity1_sens = NULL;
    int         *dentity1_dom  = NULL;
    int         *dblk_strid    = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                 (void **) interface_entity1_ln_to_gn[i_interface],
                           &dblk_strid,
                 (void **) &dentity1_gnum);
    PDM_free(dblk_strid);

    PDM_part_to_block_exch(ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                 (void **) interface_sgn[i_interface],
                           &dblk_strid,
                 (void **) &dentity1_sgn);
    PDM_free(dblk_strid);

    PDM_part_to_block_exch(ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                 (void **) interface_sens[i_interface],
                           &dblk_strid,
                 (void **) &dentity1_sens);
    PDM_free(dblk_strid);

    PDM_part_to_block_exch(ptb,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           stride_one,
                 (void **) interface_dom[i_interface],
                           &dblk_strid,
                 (void **) &dentity1_dom);


    for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
      PDM_free(interface_entity1_ln_to_gn[i_interface][i_part]);
      //PDM_free(interface_sgn             [i_interface][i_part]);
      PDM_free(interface_dom             [i_interface][i_part]);
      PDM_free(stride_one[i_part]);
    }

    PDM_free(pn_interface              [i_interface]);
    PDM_free(interface_entity1_ln_to_gn[i_interface]);
    PDM_free(interface_ln_to_gn        [i_interface]);
    PDM_free(interface_sgn             [i_interface]);
    PDM_free(interface_sens            [i_interface]);
    PDM_free(interface_dom             [i_interface]);

    int          n_gnum     = PDM_part_to_block_n_elt_block_get  (ptb);
    // PDM_g_num_t* block_gnum = PDM_part_to_block_block_gnum_get   (ptb);

    int n_data = 0;
    for(int i = 0; i < n_gnum; ++i) {
      n_data += dblk_strid[i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int (dblk_strid   , n_gnum, "dblk_strid    ::");
      PDM_log_trace_array_long(dentity1_gnum, n_data, "dentity1_gnum ::");
      PDM_log_trace_array_int (dentity1_sgn , n_data, "dentity1_sgn  ::");
      PDM_log_trace_array_int (dentity1_sens, n_data, "dentity1_sens ::");
      PDM_log_trace_array_int (dentity1_dom , n_data, "dentity1_dom  ::");
    }

    dinterface_dn [i_interface] = n_gnum;
    PDM_malloc(dinterface_ids[i_interface], 2 * n_gnum, PDM_g_num_t);
    PDM_malloc(dinterface_dom[i_interface], 2 * n_gnum, int        );

    /*
     * We can have multiple occurence for the same gnum (triple point for exemple)
     *
     */
    int max_blk_strid = 0;
    for(int i = 0; i < n_gnum; ++i) {
      max_blk_strid = PDM_MAX(max_blk_strid, dblk_strid[i]);
    }
    int         *order             = NULL;
    int         *tmp_dentity1_sgn  = NULL;
    int         *tmp_dentity1_sens = NULL;
    int         *tmp_dentity1_dom  = NULL;
    PDM_g_num_t *tmp_dentity1_gnum = NULL;
    PDM_malloc(order            , max_blk_strid, int        );
    PDM_malloc(tmp_dentity1_sgn , max_blk_strid, int        );
    PDM_malloc(tmp_dentity1_sens, max_blk_strid, int        );
    PDM_malloc(tmp_dentity1_dom , max_blk_strid, int        );
    PDM_malloc(tmp_dentity1_gnum, max_blk_strid, PDM_g_num_t);

    int *dblk_strid_unique;
    PDM_malloc(dblk_strid_unique, n_gnum, int);

    int idx_read  = 0;
    int idx_write = 0;
    for(int i = 0; i < n_gnum; ++i) {
      int beg = idx_read;
      int end = beg + dblk_strid[i];
      int n_ldata = end - beg;
      for(int k = 0; k < n_ldata; ++k) {
        order[k] = k;
        tmp_dentity1_gnum[k] = dentity1_gnum[idx_read+k];
        tmp_dentity1_sgn [k] = dentity1_sgn [idx_read+k];
        tmp_dentity1_sens[k] = dentity1_sens[idx_read+k];
        tmp_dentity1_dom [k] = dentity1_dom [idx_read+k];
      }

      PDM_sort_long(tmp_dentity1_gnum, order, n_ldata);

      PDM_g_num_t first = tmp_dentity1_gnum[0];
      dblk_strid_unique[i] = 1;

      dentity1_gnum[idx_write] = tmp_dentity1_gnum[0];
      dentity1_sgn [idx_write] = tmp_dentity1_sgn [order[0]];
      dentity1_sens[idx_write] = tmp_dentity1_sens[order[0]];
      dentity1_dom [idx_write] = tmp_dentity1_dom [order[0]];
      idx_write++;
      for(int k = 1; k < n_ldata; ++k) {
        if(first != tmp_dentity1_gnum[k]) {
          dentity1_gnum[idx_write] = tmp_dentity1_gnum[k];
          dentity1_sgn [idx_write] = tmp_dentity1_sgn [order[k]];
          dentity1_sens[idx_write] = tmp_dentity1_sens[order[k]];
          dentity1_dom [idx_write] = tmp_dentity1_dom [order[k]];
          idx_write++;

          dblk_strid_unique[i]++;
          first = tmp_dentity1_gnum[k];
        }
      }

      idx_read += dblk_strid[i];
    }

    PDM_free(order);
    PDM_free(tmp_dentity1_sgn );
    PDM_free(tmp_dentity1_sens);
    PDM_free(tmp_dentity1_dom );
    PDM_free(tmp_dentity1_gnum);
    PDM_free(dblk_strid);
    dblk_strid = dblk_strid_unique;

    n_data = 0;
    for(int i = 0; i < n_gnum; ++i) {
      n_data += dblk_strid[i];
    }

    if(0 == 1) {
      PDM_log_trace_array_int (dblk_strid   , n_gnum, "dblk_strid    (unique)::");
      PDM_log_trace_array_long(dentity1_gnum, n_data, "dentity1_gnum (unique)::");
      PDM_log_trace_array_int (dentity1_sgn , n_data, "dentity1_sgn  (unique)::");
      PDM_log_trace_array_int (dentity1_sens, n_data, "dentity1_sens (unique)::");
      PDM_log_trace_array_int (dentity1_dom , n_data, "dentity1_dom  (unique)::");
    }


    idx_read = 0;
    for(int i = 0; i < n_gnum; ++i) {

      assert(dblk_strid[i] % 2 == 0);
      for(int k = 0; k < dblk_strid[i]/2; ++k) {

        int sgn1 = dentity1_sgn[idx_read  ];
        int sgn2 = dentity1_sgn[idx_read+1];

        int sens1 = dentity1_sens[idx_read  ];
        int sens2 = dentity1_sens[idx_read+1];

        // assert(sens1 == sens2);

        PDM_g_num_t gnum1 = dentity1_gnum[idx_read  ];
        PDM_g_num_t gnum2 = dentity1_gnum[idx_read+1];

        int dom1 = dentity1_dom[idx_read  ];
        int dom2 = dentity1_dom[idx_read+1];

        assert(sgn1 * sgn2 == -1); // Forcement de sgn opposé

        if(sgn1 == -1) {
          // Swap
          dentity1_gnum[idx_read  ] = gnum2;
          dentity1_gnum[idx_read+1] = gnum1;

          dentity1_dom[idx_read  ] = dom2;
          dentity1_dom[idx_read+1] = dom1;

          dentity1_sens[idx_read  ] = sens2;
          dentity1_sens[idx_read+1] = sens1;

          gnum1 = dentity1_gnum[idx_read  ];
          gnum2 = dentity1_gnum[idx_read+1];

          dom1 = dentity1_dom[idx_read  ];
          dom2 = dentity1_dom[idx_read+1];

          sens1 = dentity1_sens[idx_read  ];
          sens2 = dentity1_sens[idx_read+1];

        }

        /**
         * For pair entities, multiple case can happen
         *   - entities        have same sens (possible with edge for example)
         *   - entities do not have same sens (normally true for faces)
         *      in this case multiple possibility can happen 
         *        - sens1 =  1 ; sens2 = -1 (typically information coming from pdm)
         *        - sens1 = -1 ; sens2 = -1 (if user do not have information of which entity
         *                                   is with the good sens, but knows that pair have
         *                                   opposite sens)
         * In both cases, we set first entity with + sens
         * and if they have opposite sens, we set the sens 2 to -1.
         */
        dinterface_ids[i_interface][2*i  ] =         gnum1;
        dinterface_ids[i_interface][2*i+1] = sens2 * gnum2;
        // dinterface_ids[i_interface][2*i  ] = sens1 * gnum1;
        // dinterface_ids[i_interface][2*i+1] = sens2 * gnum2;
        // dinterface_ids[i_interface][2*i+1] = gnum2;

        dinterface_dom[i_interface][2*i  ] = dom1;
        dinterface_dom[i_interface][2*i+1] = dom2;


        idx_read += 2;

      }
    }

    PDM_free(dentity1_gnum);
    PDM_free(dentity1_sgn);
    PDM_free(dentity1_sens);
    PDM_free(dentity1_dom);
    PDM_free(dblk_strid);
    PDM_free(stride_one);

    // PDM_log_trace_array_long(dinterface_ids[i_interface], 2 * dinterface_dn[i_interface], "PDM_part_domain_interface_to_domain_interface : dinterface_ids ::");

    PDM_part_to_block_free(ptb);

  }

  ditrf->is_result[interface_kind1] = 1;
  PDM_domain_interface_set(ditrf,
                           interface_kind1,
                           dinterface_dn,
                           dinterface_ids,
                           dinterface_dom);

  PDM_free(pn_interface);
  PDM_free(interface_ln_to_gn);
  PDM_free(interface_entity1_ln_to_gn);
  PDM_free(interface_sgn             );
  PDM_free(interface_sens            );
  PDM_free(interface_dom             );

  // log_trace("PDM_part_domain_interface_to_domain_interface end\n");

}


void
PDM_part_domain_interface_add
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind1,
 PDM_bound_type_t               interface_kind2,
 int                           *n_part,
 int                          **pn_entity1,
 PDM_g_num_t                 ***entity1_ln_to_gn,
 int                          **pn_entity2,
 PDM_g_num_t                 ***entity2_ln_to_gn,
 int                         ***entity2_entity1_idx,
 int                         ***entity2_entity1,
 int                            connectivity_is_signed
)
{
  int debug = 0;

  // log_trace("PDM_part_domain_interface_to_domain_interface %i to %i \n", interface_kind1, interface_kind2);
  int i_rank;
  PDM_MPI_Comm_rank(dom_intrf->comm, &i_rank);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    n_part_loc_all_domain += n_part[i_dom];
  }

  int n_interface = PDM_part_domain_interface_n_interface_get(dom_intrf);

  /*
   * Re-Create a domain interface
   */
  PDM_domain_interface_t  *ditrf              = NULL;
  int                    **is_entity1_on_itrf = NULL;
  PDM_part_domain_interface_to_domain_interface(dom_intrf,
                                                interface_kind1,
                                                n_part,
                                                pn_entity1,
                                                entity1_ln_to_gn,
                                                &ditrf,
                                                &is_entity1_on_itrf);

  /*
   * Extract from part and prepare the way of block
   */
  int          *dn_entity1                  = NULL;
  int          *dn_entity2                  = NULL;
  PDM_g_num_t **dfilter_entity2_entity1     = NULL;
  int         **dfilter_entity2_entity1_idx = NULL;
  PDM_malloc(dn_entity1                 , dom_intrf->n_domain, int          );
  PDM_malloc(dn_entity2                 , dom_intrf->n_domain, int          );
  PDM_malloc(dfilter_entity2_entity1    , dom_intrf->n_domain, PDM_g_num_t *);
  PDM_malloc(dfilter_entity2_entity1_idx, dom_intrf->n_domain, int         *);

  int s_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {

    int          *n_filter_entity2           = NULL;
    PDM_g_num_t **filter_entity2_entity1     = NULL;
    PDM_g_num_t **filter_entity2_ln_to_gn    = NULL;
    int         **filter_entity2_entity1_idx = NULL;
    int         **filter_entity2_entity1_n   = NULL;
    PDM_malloc(n_filter_entity2          , n_part[i_dom], int          );
    PDM_malloc(filter_entity2_entity1    , n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(filter_entity2_ln_to_gn   , n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(filter_entity2_entity1_idx, n_part[i_dom], int         *);
    PDM_malloc(filter_entity2_entity1_n  , n_part[i_dom], int         *);
    PDM_g_num_t max_gnum1 = 0;
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      int n_entity1 = pn_entity1[i_dom][i_part];
      int n_entity2 = pn_entity2[i_dom][i_part];
      int *_pentity2_entity1_idx = entity2_entity1_idx [i_dom][i_part];
      int *_pentity2_entity1     = entity2_entity1     [i_dom][i_part];
      int *_is_entity1_on_itrf   = is_entity1_on_itrf  [s_part+i_part];

      PDM_malloc(filter_entity2_entity1    [i_part],_pentity2_entity1_idx[n_entity2], PDM_g_num_t);
      PDM_malloc(filter_entity2_entity1_idx[i_part], (n_entity2 + 1)                , int        );
      PDM_malloc(filter_entity2_entity1_n  [i_part], (n_entity2    )                , int        );

      PDM_g_num_t *_filter_entity2_entity1     = filter_entity2_entity1    [i_part];
      int         *_filter_entity2_entity1_idx = filter_entity2_entity1_idx[i_part];
      int         *_filter_entity2_entity1_n   = filter_entity2_entity1_n  [i_part];

      n_filter_entity2[i_part] = 0;



      for(int i = 0; i < n_entity1; ++i) {
        max_gnum1 = PDM_MAX(max_gnum1, entity1_ln_to_gn[i_dom][i_part][i]);
      }

      _filter_entity2_entity1_idx[0] = 0;
      for(int i = 0; i < n_entity2; ++i) {

        int have_interf = 1;
        for(int idx_entity1 = _pentity2_entity1_idx[i]; idx_entity1 < _pentity2_entity1_idx[i+1]; ++idx_entity1) {
          int i_entity1 = PDM_ABS(_pentity2_entity1[idx_entity1])-1;
          if(_is_entity1_on_itrf[i_entity1] == 0) {
            have_interf = 0;
            break;
          }
        }

        if(have_interf) {
          _filter_entity2_entity1_idx[n_filter_entity2[i_part]+1] = _filter_entity2_entity1_idx[n_filter_entity2[i_part]];
          _filter_entity2_entity1_n[n_filter_entity2[i_part]] = 0;
          for(int idx_entity1 = _pentity2_entity1_idx[i]; idx_entity1 < _pentity2_entity1_idx[i+1]; ++idx_entity1) {
            int sgn       = PDM_SIGN(_pentity2_entity1[idx_entity1]);
            int i_entity1 = PDM_ABS (_pentity2_entity1[idx_entity1])-1;
            int idx_write = _filter_entity2_entity1_idx[n_filter_entity2[i_part]+1]++;
            _filter_entity2_entity1_n[n_filter_entity2[i_part]]++;
            _filter_entity2_entity1[idx_write] = sgn * entity1_ln_to_gn[i_dom][i_part][i_entity1];
          }
          n_filter_entity2[i_part]++;
        } else {
          _filter_entity2_entity1_idx[n_filter_entity2[i_part]+1] = _filter_entity2_entity1_idx[n_filter_entity2[i_part]];
          _filter_entity2_entity1_n[n_filter_entity2[i_part]] = 0;
          n_filter_entity2[i_part]++;
        }
      }

      /*
       * Realloc
       */
      if(debug == 1) {
        PDM_log_trace_connectivity_long(_filter_entity2_entity1_idx, _filter_entity2_entity1, n_filter_entity2[i_part], "_filter_entity2_entity1 ::");
      }

      PDM_realloc(filter_entity2_entity1    [i_part] ,filter_entity2_entity1    [i_part] , _filter_entity2_entity1_idx[n_filter_entity2[i_part]] ,PDM_g_num_t);
      PDM_realloc(filter_entity2_entity1_idx[i_part] ,filter_entity2_entity1_idx[i_part] ,  (n_filter_entity2[i_part] + 1)                       ,int        );
      filter_entity2_ln_to_gn   [i_part] = entity2_ln_to_gn[i_dom][i_part];

      assert(n_filter_entity2[i_part] == n_entity2);
    }


    /*
     * Return to block vision
     */
    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        filter_entity2_ln_to_gn,
                                                        NULL,
                                                        n_filter_entity2,
                                                        n_part[i_dom],
                                                        dom_intrf->comm);

    int         *dfilter_entity2_entity1_n = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           filter_entity2_entity1_n,
                 (void **) filter_entity2_entity1,
                           &dfilter_entity2_entity1_n,
                 (void **) &dfilter_entity2_entity1[i_dom]);


    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      PDM_free(filter_entity2_entity1    [i_part]);
      PDM_free(filter_entity2_entity1_idx[i_part]);
      PDM_free(filter_entity2_entity1_n  [i_part]);

    }
    PDM_free(filter_entity2_entity1    );
    PDM_free(filter_entity2_entity1_idx);
    PDM_free(filter_entity2_entity1_n  );
    PDM_free(filter_entity2_ln_to_gn   );
    PDM_free(n_filter_entity2);

    dn_entity2[i_dom]                  = PDM_part_to_block_n_elt_block_get  (ptb);
    dfilter_entity2_entity1_idx[i_dom] = PDM_array_new_idx_from_sizes_int(dfilter_entity2_entity1_n, dn_entity2[i_dom] );

    if(debug == 1) {
      PDM_log_trace_connectivity_long(dfilter_entity2_entity1_idx[i_dom], dfilter_entity2_entity1[i_dom], dn_entity2[i_dom] , "dfilter_entity2_entity1 ::");
    }

    PDM_free(dfilter_entity2_entity1_n);

    PDM_part_to_block_free(ptb);

    // Compute dn_entity1
    PDM_g_num_t n_g_entity1 = -1;
    PDM_MPI_Allreduce(&max_gnum1, &n_g_entity1, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ditrf->comm);

    dn_entity1[i_dom] = PDM_compute_uniform_dn_entity(ditrf->comm, n_g_entity1);

    s_part += n_part[i_dom];
  }

  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    PDM_free(is_entity1_on_itrf[i_part]);
  }
  PDM_free(is_entity1_on_itrf);

  /*
   * Translate in distributed
   */
  if(0 == 1) {
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      PDM_log_trace_array_long(ditrf->interface_ids_vtx[i_interface], 2 *ditrf->interface_dn_vtx[i_interface], "interface_ids_vtx ::");
      PDM_log_trace_array_int (ditrf->interface_dom_vtx[i_interface], 2 *ditrf->interface_dn_vtx[i_interface], "interface_dom_vtx ::");
    }
  }

  /*
   * Management of cases
   */
  int          *kind1_interface_dn   = NULL;
  // int         **kind1_dinterface_dom = dinterface_dom;
  int         **kind1_dinterface_dom = NULL;
  PDM_g_num_t **kind1_dinterface_ids = NULL;
  if(interface_kind1 == PDM_BOUND_TYPE_VTX) {
    kind1_interface_dn   = ditrf->interface_dn_vtx;
    kind1_dinterface_dom = ditrf->interface_dom_vtx;
    kind1_dinterface_ids = ditrf->interface_ids_vtx;
  } else if(interface_kind1 == PDM_BOUND_TYPE_EDGE) {
    kind1_interface_dn   = ditrf->interface_dn_edge;
    kind1_dinterface_dom = ditrf->interface_dom_edge;
    kind1_dinterface_ids = ditrf->interface_ids_edge;
  }

  int          **kind2_interface_dn   = NULL;
  int         ***kind2_dinterface_dom = NULL;
  PDM_g_num_t ***kind2_dinterface_ids = NULL;
  if(interface_kind2 == PDM_BOUND_TYPE_EDGE) {
    kind2_interface_dn   = &ditrf->interface_dn_edge;
    kind2_dinterface_dom = &ditrf->interface_dom_edge;
    kind2_dinterface_ids = &ditrf->interface_ids_edge;
  } else if(interface_kind2 == PDM_BOUND_TYPE_FACE) {
    kind2_interface_dn   = &ditrf->interface_dn_face;
    kind2_dinterface_dom = &ditrf->interface_dom_face;
    kind2_dinterface_ids = &ditrf->interface_ids_face;
  }

  PDM_domain_interface_translate_entity1_entity2(ditrf->n_domain,
                                                 ditrf->n_interface,
                                                 dn_entity1,
                                                 dn_entity2,
                                                 kind1_interface_dn,
                                                 kind1_dinterface_dom,
                                                 kind1_dinterface_ids,
                                                 dfilter_entity2_entity1_idx,
                                                 dfilter_entity2_entity1,
                                                 connectivity_is_signed,
                                                 ditrf->comm,
                                                 kind2_interface_dn,
                                                 kind2_dinterface_ids,
                                                 kind2_dinterface_dom);

  if(debug == 1) {
    for(int i = 0; i < ditrf->n_interface; ++i) {
      PDM_log_trace_array_long((*kind2_dinterface_ids)[i], 2 *  (*kind2_interface_dn)[i], "kind2_dinterface_ids ::");
    }
  }

  ditrf->is_result[interface_kind2] = 1;

  PDM_ddomain_interface_to_pdomain_interface(ditrf->comm,
                                             ditrf->n_interface,
                                             ditrf->n_domain,
                                             ditrf->multidomain_intrf,
                                             interface_kind2,
                                             *(kind2_interface_dn),
                                             *(kind2_dinterface_ids),
                                             *(kind2_dinterface_dom),
                                             n_part,
                                             pn_entity2,
                                             entity2_ln_to_gn,
                                             dom_intrf);
  PDM_free(dn_entity1);
  PDM_free(dn_entity2);

  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    PDM_free(dfilter_entity2_entity1_idx[i_dom]);
    PDM_free(dfilter_entity2_entity1    [i_dom]);
  }
  PDM_free(dfilter_entity2_entity1_idx);
  PDM_free(dfilter_entity2_entity1);


  PDM_domain_interface_free(ditrf);

  // log_trace("PDM_part_domain_interface_to_domain_interface %i to %i \n", interface_kind1, interface_kind2);

}


void
PDM_part_domain_interface_face2vtx
(
 PDM_part_domain_interface_t   *dom_intrf,
 int                           *n_part,
 int                          **pn_face,
 PDM_g_num_t                 ***pface_ln_to_gn,
 int                          **pn_vtx,
 PDM_g_num_t                 ***pvtx_ln_to_gn,
 int                         ***pface_vtx_idx,
 int                         ***pface_vtx
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dom_intrf->comm, &i_rank);

  int n_part_loc_all_domain = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    n_part_loc_all_domain += n_part[i_dom];
  }

  int n_interface = PDM_part_domain_interface_n_interface_get(dom_intrf);
  int          **pn_interface            = NULL;
  PDM_g_num_t ***interface_face_ln_to_gn = NULL;
  PDM_g_num_t ***interface_ln_to_gn      = NULL;
  int         ***interface_sgn           = NULL;
  int         ***interface_sens          = NULL;
  int         ***interface_dom           = NULL;
  PDM_malloc(pn_interface           , n_interface, int          *);
  PDM_malloc(interface_face_ln_to_gn, n_interface, PDM_g_num_t **);
  PDM_malloc(interface_ln_to_gn     , n_interface, PDM_g_num_t **);
  PDM_malloc(interface_sgn          , n_interface, int         **);
  PDM_malloc(interface_sens         , n_interface, int         **);
  PDM_malloc(interface_dom          , n_interface, int         **);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    PDM_malloc(pn_interface           [i_interface], n_part_loc_all_domain, int          );
    PDM_malloc(interface_face_ln_to_gn[i_interface], n_part_loc_all_domain, PDM_g_num_t *);
    PDM_malloc(interface_ln_to_gn     [i_interface], n_part_loc_all_domain, PDM_g_num_t *);
    PDM_malloc(interface_sgn          [i_interface], n_part_loc_all_domain, int         *);
    PDM_malloc(interface_sens         [i_interface], n_part_loc_all_domain, int         *);
    PDM_malloc(interface_dom          [i_interface], n_part_loc_all_domain, int         *);
  }

  /*
   * Re-Create a domain interface
   */
  PDM_domain_interface_t  *ditrf           = NULL;
  int                    **is_face_on_itrf = NULL;
  PDM_part_domain_interface_to_domain_interface(dom_intrf,
                                                PDM_BOUND_TYPE_FACE,
                                                n_part,
                                                pn_face,
                                                pface_ln_to_gn,
                                                &ditrf,
                                                &is_face_on_itrf);

  /*
   * Extract from part and prepare the way of block
   */
  int          *dn_vtx               = NULL;
  int          *dn_face              = NULL;
  PDM_g_num_t **dfilter_face_vtx     = NULL;
  int         **dfilter_face_vtx_idx = NULL;
  PDM_malloc(dn_vtx              , dom_intrf->n_domain, int          );
  PDM_malloc(dn_face             , dom_intrf->n_domain, int          );
  PDM_malloc(dfilter_face_vtx    , dom_intrf->n_domain, PDM_g_num_t *);
  PDM_malloc(dfilter_face_vtx_idx, dom_intrf->n_domain, int         *);

  int s_part = 0;
  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {

    int          *n_filter_face        = NULL;
    PDM_g_num_t **filter_face_vtx      = NULL;
    PDM_g_num_t **filter_face_ln_to_gn = NULL;
    int         **filter_face_vtx_idx  = NULL;
    int         **filter_face_vtx_n    = NULL;
    PDM_malloc(n_filter_face       , n_part[i_dom], int          );
    PDM_malloc(filter_face_vtx     , n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(filter_face_ln_to_gn, n_part[i_dom], PDM_g_num_t *);
    PDM_malloc(filter_face_vtx_idx , n_part[i_dom], int         *);
    PDM_malloc(filter_face_vtx_n   , n_part[i_dom], int         *);
    PDM_g_num_t max_gnum_vtx = 0;
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      int n_face = pn_face[i_dom][i_part];
      int n_vtx  = pn_vtx [i_dom][i_part];
      int *_pface_vtx_idx   = pface_vtx_idx [i_dom][i_part];
      int *_pface_vtx       = pface_vtx     [i_dom][i_part];
      int *_is_face_on_itrf = is_face_on_itrf  [s_part+i_part];

      PDM_malloc(filter_face_vtx    [i_part], _pface_vtx_idx[n_face], PDM_g_num_t);
      PDM_malloc(filter_face_vtx_idx[i_part],  n_face + 1           , int        );
      PDM_malloc(filter_face_vtx_n  [i_part],  n_face               , int        );

      PDM_g_num_t *_filter_face_vtx     = filter_face_vtx    [i_part];
      int         *_filter_face_vtx_idx = filter_face_vtx_idx[i_part];
      int         *_filter_face_vtx_n   = filter_face_vtx_n  [i_part];

      n_filter_face[i_part] = 0;

      for(int i = 0; i < n_vtx; ++i) {
        max_gnum_vtx = PDM_MAX(max_gnum_vtx, pvtx_ln_to_gn[i_dom][i_part][i]);
      }

      /* Count */
      _filter_face_vtx_idx[0] = 0;
      for(int i = 0; i < n_face; ++i) {

        if(_is_face_on_itrf[i] == 0) {
          _filter_face_vtx_idx[n_filter_face[i_part]+1] = _filter_face_vtx_idx[n_filter_face[i_part]];
          _filter_face_vtx_n[n_filter_face[i_part]] = 0;
          n_filter_face[i_part]++;
        } else {
          _filter_face_vtx_idx[n_filter_face[i_part]+1] = _filter_face_vtx_idx[n_filter_face[i_part]];
          _filter_face_vtx_n[n_filter_face[i_part]] = 0;
          for(int idx_vtx = _pface_vtx_idx[i]; idx_vtx < _pface_vtx_idx[i+1]; ++idx_vtx) {
            int sgn   = PDM_SIGN(_pface_vtx[idx_vtx]);
            int i_vtx = PDM_ABS (_pface_vtx[idx_vtx])-1;
            int idx_write = _filter_face_vtx_idx[n_filter_face[i_part]+1]++;
            _filter_face_vtx_n[n_filter_face[i_part]]++;
            _filter_face_vtx[idx_write] = sgn * pvtx_ln_to_gn[i_dom][i_part][i_vtx];
          }
          n_filter_face[i_part]++;
        }

      }

      /*
       * Realloc
       */
      if(0 == 1) {
        PDM_log_trace_connectivity_long(_filter_face_vtx_idx, _filter_face_vtx, n_filter_face[i_part], "_filter_face_vtx ::");
      }

      PDM_realloc(filter_face_vtx     [i_part] ,filter_face_vtx     [i_part] , _filter_face_vtx_idx[n_filter_face[i_part]] ,PDM_g_num_t);
      PDM_realloc(filter_face_vtx_idx [i_part] ,filter_face_vtx_idx [i_part] ,  (n_filter_face[i_part] + 1)                       ,int        );
      filter_face_ln_to_gn[i_part] = pface_ln_to_gn[i_dom][i_part];

      assert(n_filter_face[i_part] == n_face);

    }

    /*
     * Return to block vision
     */
    PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        filter_face_ln_to_gn,
                                                        NULL,
                                                        n_filter_face,
                                                        n_part[i_dom],
                                                        dom_intrf->comm);
    int         *dfilter_face_vtx_n = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           filter_face_vtx_n,
                 (void **) filter_face_vtx,
                           &dfilter_face_vtx_n,
                 (void **) &dfilter_face_vtx[i_dom]);


    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      PDM_free(filter_face_vtx    [i_part]);
      PDM_free(filter_face_vtx_idx[i_part]);
      PDM_free(filter_face_vtx_n  [i_part]);

    }
    PDM_free(filter_face_vtx    );
    PDM_free(filter_face_vtx_idx);
    PDM_free(filter_face_vtx_n  );
    PDM_free(filter_face_ln_to_gn);
    PDM_free(n_filter_face);

    dn_face[i_dom]                  = PDM_part_to_block_n_elt_block_get  (ptb);
    dfilter_face_vtx_idx[i_dom] = PDM_array_new_idx_from_sizes_int(dfilter_face_vtx_n, dn_face[i_dom] );

    if(0 == 1) {
      PDM_log_trace_connectivity_long(dfilter_face_vtx_idx[i_dom], dfilter_face_vtx[i_dom], dn_face[i_dom] , "dfilter_face_vtx ::");
    }

    PDM_free(dfilter_face_vtx_n);

    PDM_part_to_block_free(ptb);

    // Compute dn_vtx
    PDM_g_num_t n_g_vtx = -1;
    PDM_MPI_Allreduce(&max_gnum_vtx, &n_g_vtx, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, ditrf->comm);

    dn_vtx[i_dom] = PDM_compute_uniform_dn_entity(ditrf->comm, n_g_vtx);

    s_part += n_part[i_dom];


  }

  for(int i_part = 0; i_part < n_part_loc_all_domain; ++i_part) {
    PDM_free(is_face_on_itrf[i_part]);
  }
  PDM_free(is_face_on_itrf);

  PDM_domain_interface_translate_face2vtx(ditrf,
                                          dn_vtx,
                                          dn_face,
                                          dfilter_face_vtx_idx,
                                          dfilter_face_vtx);

  /*
   * Translate
   */
  PDM_ddomain_interface_to_pdomain_interface(ditrf->comm,
                                             ditrf->n_interface,
                                             ditrf->n_domain,
                                             ditrf->multidomain_intrf,
                                             PDM_BOUND_TYPE_VTX,
                                             ditrf->interface_dn_vtx,
                                             ditrf->interface_ids_vtx,
                                             ditrf->interface_dom_vtx,
                                             n_part,
                                             pn_vtx,
                                             pvtx_ln_to_gn,
                                             dom_intrf);
  ditrf->is_result[PDM_BOUND_TYPE_VTX] = 1;
  PDM_free(dn_vtx);
  PDM_free(dn_face);

  for(int i_dom = 0; i_dom < dom_intrf->n_domain; ++i_dom) {
    PDM_free(dfilter_face_vtx_idx[i_dom]);
    PDM_free(dfilter_face_vtx    [i_dom]);
  }
  PDM_free(dfilter_face_vtx_idx);
  PDM_free(dfilter_face_vtx);


  PDM_domain_interface_free(ditrf);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    PDM_free(pn_interface           [i_interface]);
    PDM_free(interface_face_ln_to_gn[i_interface]);
    PDM_free(interface_ln_to_gn     [i_interface]);
    PDM_free(interface_sgn          [i_interface]);
    PDM_free(interface_sens         [i_interface]);
    PDM_free(interface_dom          [i_interface]);
  }
  PDM_free(pn_interface);
  PDM_free(interface_ln_to_gn);
  PDM_free(interface_face_ln_to_gn);
  PDM_free(interface_sgn             );
  PDM_free(interface_sens            );
  PDM_free(interface_dom             );

}






