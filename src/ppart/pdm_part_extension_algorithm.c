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
#include "pdm_part_extension_algorithm.h"
#include "pdm_part_extension_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_domain_interface.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
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
_concatenate_part_graph
(
 int           pn_entity1,
 int          *part1_to_part2_idx,
 int          *part1_to_part2_interface,
 int          *bentity1_entity2_n,
 PDM_g_num_t  *bentity1_entity2_gnum,
 int          *bentity1_entity2_triplet,
 int          *bentity1_entity2_interface_n,
 int          *bentity1_entity2_interface_tot_n,
 int          *bentity1_entity2_interface,
 int          *pnext_bentity1_entity2_n,
 PDM_g_num_t  *pnext_bentity1_entity2_gnum,
 int          *pnext_bentity1_entity2_triplet,
 int          *pnext_bentity1_entity2_interface_n,
 int          *pnext_bentity1_entity2_interface_tot_n,
 int          *pnext_bentity1_entity2_interface,
 int         **concat_bentity1_entity2_n_out,
 PDM_g_num_t **concat_bentity1_entity2_gnum_out,
 int         **concat_bentity1_entity2_triplet_out,
 int         **concat_bentity1_entity2_interface_n_out,
 int         **concat_bentity1_entity2_interface_tot_n_out,
 int         **concat_bentity1_entity2_interface_out
)
{
  PDM_UNUSED(part1_to_part2_interface);

  int         *concat_bentity1_entity2_n               = malloc(pn_entity1 * sizeof(int));;
  int         *concat_bentity1_entity2_interface_tot_n = malloc(pn_entity1 * sizeof(int));;

  int n_concat_entity1_entity2      = 0;
  int n_concat_entity1_entity2_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {
    concat_bentity1_entity2_n              [i] = 0;
    concat_bentity1_entity2_interface_tot_n[i] = 0;

    /* From previous in layout of part2 */
    concat_bentity1_entity2_n              [i] += bentity1_entity2_n[i];
    concat_bentity1_entity2_interface_tot_n[i] += bentity1_entity2_interface_tot_n[i];

    /* From recv in layout part1_to_part2 */
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {
      concat_bentity1_entity2_n              [i] += pnext_bentity1_entity2_n              [idx];
      concat_bentity1_entity2_interface_tot_n[i] += pnext_bentity1_entity2_interface_tot_n[idx];
    }

    // Update global count
    n_concat_entity1_entity2      += concat_bentity1_entity2_n              [i];
    n_concat_entity1_entity2_itrf += concat_bentity1_entity2_interface_tot_n[i];
  }

  /* Allocate */
  PDM_g_num_t *concat_bentity1_entity2_gnum        = malloc(    n_concat_entity1_entity2      * sizeof(PDM_g_num_t));
  int         *concat_bentity1_entity2_triplet     = malloc(3 * n_concat_entity1_entity2      * sizeof(int        ));
  int         *concat_bentity1_entity2_interface_n = malloc(    n_concat_entity1_entity2      * sizeof(int        ));
  int         *concat_bentity1_entity2_interface   = malloc(    n_concat_entity1_entity2_itrf * sizeof(int        ));

  n_concat_entity1_entity2      = 0;
  n_concat_entity1_entity2_itrf = 0;
  int idx_readb      = 0;
  int idx_readb_itrf = 0;
  int idx_readp      = 0;
  int idx_readp_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {

    /* From previous in layout of part2 */
    for(int k = 0; k < bentity1_entity2_n[i]; ++k) {
      concat_bentity1_entity2_gnum       [  n_concat_entity1_entity2  ] = bentity1_entity2_gnum       [  idx_readb  ];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2  ] = bentity1_entity2_triplet    [3*idx_readb  ];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+1] = bentity1_entity2_triplet    [3*idx_readb+1];
      concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+2] = bentity1_entity2_triplet    [3*idx_readb+2];
      concat_bentity1_entity2_interface_n[  n_concat_entity1_entity2  ] = bentity1_entity2_interface_n[  idx_readb  ];

      for(int p = 0; p <  bentity1_entity2_interface_n[  idx_readb  ]; ++p) {
        concat_bentity1_entity2_interface[n_concat_entity1_entity2_itrf++] = bentity1_entity2_interface[idx_readb_itrf++];
      }
      idx_readb++;
      n_concat_entity1_entity2++;
    }

    /* From recv in layout part1_to_part2 */
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {

      /* From previous in layout of part2 */
      for(int k = 0; k < pnext_bentity1_entity2_n[idx]; ++k) {
        concat_bentity1_entity2_gnum       [  n_concat_entity1_entity2  ] = pnext_bentity1_entity2_gnum       [  idx_readp  ];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2  ] = pnext_bentity1_entity2_triplet    [3*idx_readp  ];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+1] = pnext_bentity1_entity2_triplet    [3*idx_readp+1];
        concat_bentity1_entity2_triplet    [3*n_concat_entity1_entity2+2] = pnext_bentity1_entity2_triplet    [3*idx_readp+2];
        concat_bentity1_entity2_interface_n[  n_concat_entity1_entity2  ] = pnext_bentity1_entity2_interface_n[  idx_readp  ];

        for(int p = 0; p <  pnext_bentity1_entity2_interface_n[  idx_readp  ]; ++p) {
          concat_bentity1_entity2_interface[n_concat_entity1_entity2_itrf++] = pnext_bentity1_entity2_interface[idx_readp_itrf++];
        }
        idx_readp++;
        n_concat_entity1_entity2++;
      }
    }
  }


  *concat_bentity1_entity2_n_out               = concat_bentity1_entity2_n;
  *concat_bentity1_entity2_gnum_out            = concat_bentity1_entity2_gnum;
  *concat_bentity1_entity2_triplet_out         = concat_bentity1_entity2_triplet;
  *concat_bentity1_entity2_interface_n_out     = concat_bentity1_entity2_interface_n;
  *concat_bentity1_entity2_interface_tot_n_out = concat_bentity1_entity2_interface_tot_n;
  *concat_bentity1_entity2_interface_out       = concat_bentity1_entity2_interface;

}

static
void
_dump_graph_info
(
 int           pn_entity1,
 int          *bentity1_entity2_n,
 PDM_g_num_t  *bentity1_entity2_gnum,
 int          *bentity1_entity2_triplet,
 int          *bentity1_entity2_interface_n,
 int          *bentity1_entity2_interface_tot_n,
 int          *bentity1_entity2_interface
)
{

  PDM_UNUSED(bentity1_entity2_interface);

  log_trace(" -------------------------- graph info -------------------------- \n");
  int idx_readb      = 0;
  int idx_readb_itrf = 0;
  for(int i = 0; i < pn_entity1; ++i) {
    log_trace("i_entity = %i - bentity1_entity2_n = %i - bentity1_entity2_interface_tot_n = %i\n", i, bentity1_entity2_n[i], bentity1_entity2_interface_tot_n[i]);
    /* From previous in layout of part2 */
    for(int k = 0; k < bentity1_entity2_n[i]; ++k) {
      // log_trace("\t bentity1_entity2_interface_n = %i \n", bentity1_entity2_interface_n[idx_readb]);
      log_trace("\t (gnum=%i, [%i/%i/%i] ", bentity1_entity2_gnum[idx_readb],
                                              bentity1_entity2_triplet[3*idx_readb],
                                              bentity1_entity2_triplet[3*idx_readb+1],
                                              bentity1_entity2_triplet[3*idx_readb+2]);
      if(bentity1_entity2_interface_n[idx_readb] > 0) {
        log_trace(" - Interface : ");
      }
      for(int p = 0; p < bentity1_entity2_interface_n[idx_readb]; ++p) {
        log_trace("%i ", bentity1_entity2_interface[idx_readb_itrf++]);
      }
      log_trace("\n");

      idx_readb++;
    }
  }


}



static
void
exchange_and_concat_part_graph
(
 PDM_part_to_part_t  *ptp,
 int                  n_part_tot,
 int                 *pn_entity1,
 int                **part1_to_part2_idx,
 int                **part1_to_part2_interface,
 int                **bentity1_entity2_n,
 PDM_g_num_t        **bentity1_entity2_gnum,
 int                **bentity1_entity2_triplet,
 int                **bentity1_entity2_interface_n,
 int                **bentity1_entity2_interface_tot_n,
 int                **bentity1_entity2_interface,
 int               ***pnext_bentity1_entity2_n_out,
 PDM_g_num_t       ***pnext_bentity1_entity2_gnum_out,
 int               ***pnext_bentity1_entity2_triplet_out,
 int               ***pnext_bentity1_entity2_interface_n_out,
 int               ***pnext_bentity1_entity2_interface_tot_n_out,
 int               ***pnext_bentity1_entity2_interface_out
)
{
  int exch_request = 0;
  int         **pnext_bentity1_entity2_n    = NULL;
  PDM_g_num_t **pnext_bentity1_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_gnum,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pnext_bentity1_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 3 * sizeof(int),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_triplet,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  /*
   * Exchange stride or interface
   */
  int **pnext_bentity1_entity2_interface_n = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(int),
                (const int **)   bentity1_entity2_n,
                (const void **)  bentity1_entity2_interface_n,
                                 &pnext_bentity1_entity2_n,
                    (void ***)   &pnext_bentity1_entity2_interface_n,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  int **pnext_bentity1_entity2_interface_tot_n = NULL;
  int **pnext_bentity1_entity2_interface       = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(int),
                (const int **)   bentity1_entity2_interface_tot_n,
                (const void **)  bentity1_entity2_interface,
                                 &pnext_bentity1_entity2_interface_tot_n,
                    (void ***)   &pnext_bentity1_entity2_interface,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  /*
   * Tout a été échanger avec l'interface courante, donc on rajoute à l'arrivé la provenace
   */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    //

    int n_interface_tot = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int idx = part1_to_part2_idx[i_part][i]/3; idx < part1_to_part2_idx[i_part][i+1]/3; ++idx) {
        n_interface_tot += pnext_bentity1_entity2_interface_tot_n[i_part][idx];
        if(part1_to_part2_interface[i_part][idx] != 0) {
          n_interface_tot += pnext_bentity1_entity2_n[i_part][idx];
        }
      }
    }

    int *pupdate_bentity1_entity2_interface = malloc(n_interface_tot * sizeof(int));

    n_interface_tot = 0;
    int idx_read      = 0;
    int idx_read_itrf = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int idx = part1_to_part2_idx[i_part][i]/3; idx < part1_to_part2_idx[i_part][i+1]/3; ++idx) {


        for(int k = 0; k < pnext_bentity1_entity2_n[i_part][idx]; ++k) {


          for(int p = 0; p < pnext_bentity1_entity2_interface_n[i_part][idx_read]; ++p) {
            pupdate_bentity1_entity2_interface[n_interface_tot] = pnext_bentity1_entity2_interface[i_part][idx_read_itrf++];
            n_interface_tot++;
          }

          if(part1_to_part2_interface[i_part][idx] != 0) {
            pnext_bentity1_entity2_interface_tot_n[i_part][idx] += 1;
            pnext_bentity1_entity2_interface_n    [i_part][idx_read] += 1;
            pupdate_bentity1_entity2_interface[n_interface_tot++] = -part1_to_part2_interface[i_part][idx];
          }
          idx_read++;
        }

      }
    }

    free(pnext_bentity1_entity2_interface[i_part]);
    pnext_bentity1_entity2_interface[i_part] = pupdate_bentity1_entity2_interface;


  }



  *pnext_bentity1_entity2_n_out               = pnext_bentity1_entity2_n;
  *pnext_bentity1_entity2_gnum_out            = pnext_bentity1_entity2_gnum;
  *pnext_bentity1_entity2_triplet_out         = pnext_bentity1_entity2_triplet;
  *pnext_bentity1_entity2_interface_n_out     = pnext_bentity1_entity2_interface_n;
  *pnext_bentity1_entity2_interface_tot_n_out = pnext_bentity1_entity2_interface_tot_n;
  *pnext_bentity1_entity2_interface_out       = pnext_bentity1_entity2_interface;
}



static
void
_update_in_part2_layout
(
 int   pn_entity1,
 int  *part1_to_part2_idx,
 int  *pnext_bentity1_entity2_n,
 int  *pnext_bentity1_entity2_interface_tot_n,
 int **update_bentity1_entity2_n_out,
 int **update_bentity1_entity2_interface_tot_n_out
)
{
  int *update_bentity1_entity2_n               = malloc(pn_entity1 * sizeof(int));
  int *update_bentity1_entity2_interface_tot_n = malloc(pn_entity1 * sizeof(int));
  for(int i = 0; i < pn_entity1; ++i) {
    update_bentity1_entity2_n              [i] = 0;
    update_bentity1_entity2_interface_tot_n[i] = 0;
    for(int idx = part1_to_part2_idx[i]/3; idx < part1_to_part2_idx[i+1]/3; ++idx) {
      update_bentity1_entity2_n              [i] += pnext_bentity1_entity2_n              [idx];
      update_bentity1_entity2_interface_tot_n[i] += pnext_bentity1_entity2_interface_tot_n[idx];
    }
  }
  *update_bentity1_entity2_n_out               = update_bentity1_entity2_n;
  *update_bentity1_entity2_interface_tot_n_out = update_bentity1_entity2_interface_tot_n;
}



static
void
_recurse_and_filter
(
 PDM_part_to_part_t  *ptp,
 int                  n_part_tot,
 int                 *pn_entity1,
 int                **part1_to_part2_idx,
 int                **part1_to_part2_interface,
 int                **pentity1_entity2_idx,
 int                **pentity1_entity2,
 PDM_g_num_t        **pentity2_ln_to_gn,
 int                **pextract_entity2_n,
 PDM_g_num_t        **pextract_entity2_gnum,
 int                **pextract_entity2_triplet,
 PDM_MPI_Comm         comm
)
{

  PDM_UNUSED(pextract_entity2_n);
  PDM_UNUSED(pextract_entity2_gnum);
  PDM_UNUSED(pextract_entity2_triplet);


  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   * We wan to build the complete graphe betwenn entity1 and entity2
   * But at each step we need to fix it with transformation
   * We need to build :
   *     - bentity1_entity2_n
   *     - bentity1_entity2_gnum
   */
  int         **bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));
  int         **bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));

  /* Init */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    bentity1_entity2_n              [i_part] = malloc(pn_entity1[i_part] * sizeof(int));
    bentity1_entity2_interface_tot_n[i_part] = malloc(pn_entity1[i_part] * sizeof(int));

    int         *_pentity1_entity2_idx             = pentity1_entity2_idx            [i_part];
    int         *_pentity1_entity2                 = pentity1_entity2                [i_part];
    PDM_g_num_t *_pentity2_ln_to_gn                = pentity2_ln_to_gn               [i_part];
    // int         *_part1_to_part2_idx               = part1_to_part2_idx              [i_part];
    int         *_bentity1_entity2_n               = bentity1_entity2_n              [i_part];
    int         *_bentity1_entity2_interface_tot_n = bentity1_entity2_interface_tot_n[i_part];

    int n_entity1_entity2      = 0;
    int n_entity1_entity2_itrf = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      _bentity1_entity2_n              [i] = 0;
      _bentity1_entity2_interface_tot_n[i] = 0;
      _bentity1_entity2_n[i] += _pentity1_entity2_idx[i+1] - _pentity1_entity2_idx[i];
      n_entity1_entity2 += _bentity1_entity2_n[i];
    }

    bentity1_entity2_gnum       [i_part] = malloc(    n_entity1_entity2      * sizeof(PDM_g_num_t));
    bentity1_entity2_triplet    [i_part] = malloc(3 * n_entity1_entity2      * sizeof(int        ));
    bentity1_entity2_interface_n[i_part] = malloc(    n_entity1_entity2      * sizeof(int        ));
    bentity1_entity2_interface  [i_part] = malloc(    n_entity1_entity2_itrf * sizeof(int        ));
    PDM_g_num_t *_bentity1_entity2_gnum        = bentity1_entity2_gnum       [i_part];
    int         *_bentity1_entity2_triplet     = bentity1_entity2_triplet    [i_part];
    int         *_bentity1_entity2_interface_n = bentity1_entity2_interface_n[i_part];
    int         *_bentity1_entity2_interface   = bentity1_entity2_interface  [i_part];

    n_entity1_entity2 = 0;
    for(int i = 0; i < pn_entity1[i_part]; ++i) {
      for(int j = _pentity1_entity2_idx[i]; j < _pentity1_entity2_idx[i+1]; ++j) {
        int i_entity2 = PDM_ABS(_pentity1_entity2[j])-1;
        _bentity1_entity2_gnum       [n_entity1_entity2] = _pentity2_ln_to_gn[i_entity2];
        _bentity1_entity2_interface_n[n_entity1_entity2] = 0;

        _bentity1_entity2_interface_tot_n[i] += _bentity1_entity2_interface_n[n_entity1_entity2];

        _bentity1_entity2_triplet    [3*n_entity1_entity2  ] = i_rank;
        _bentity1_entity2_triplet    [3*n_entity1_entity2+1] = i_part;
        _bentity1_entity2_triplet    [3*n_entity1_entity2+2] = i_entity2;

        n_entity1_entity2++;
      }
    }

    PDM_log_trace_array_int (_bentity1_entity2_n              , pn_entity1[i_part]    , "_bentity1_entity2_n           ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface_tot_n, pn_entity1[i_part]    , "_bentity1_entity2_n           ::");
    PDM_log_trace_array_long(_bentity1_entity2_gnum           , n_entity1_entity2     , "_bentity1_entity2_gnum        ::");
    PDM_log_trace_array_int (_bentity1_entity2_triplet        , 3 * n_entity1_entity2 , "_bentity1_entity2_gnum        ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface_n    , n_entity1_entity2     , "_bentity1_entity2_interface_n ::");
    PDM_log_trace_array_int (_bentity1_entity2_interface      , n_entity1_entity2_itrf, "_bentity1_entity2_interface   ::");
  }

  /*
   * Previous concat
   *
   */
  int         **prev_concat_bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **prev_concat_bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **prev_concat_bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **prev_concat_bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    prev_concat_bentity1_entity2_n              [i_part] = bentity1_entity2_n              [i_part];
    prev_concat_bentity1_entity2_gnum           [i_part] = bentity1_entity2_gnum           [i_part];
    prev_concat_bentity1_entity2_triplet        [i_part] = bentity1_entity2_triplet        [i_part];
    prev_concat_bentity1_entity2_interface_n    [i_part] = bentity1_entity2_interface_n    [i_part];
    prev_concat_bentity1_entity2_interface_tot_n[i_part] = bentity1_entity2_interface_tot_n[i_part];
    prev_concat_bentity1_entity2_interface      [i_part] = bentity1_entity2_interface      [i_part];
  }

  int         **concat_bentity1_entity2_n               = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **concat_bentity1_entity2_gnum            = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **concat_bentity1_entity2_triplet         = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface_n     = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface_tot_n = malloc(n_part_tot * sizeof(int         *));
  int         **concat_bentity1_entity2_interface       = malloc(n_part_tot * sizeof(int         *));

  for(int i_step = 0; i_step < 2; ++i_step)  {

    log_trace("------- STEP %i ------- \n", i_step);
    /*
     * Il faudrait refaire l'échange qu'on fait en gnum_come_from avec ce layout pour l'appeler en part2
     * On doit avoir le même resultats
     */
    int         **pnext_bentity1_entity2_n               = NULL;
    PDM_g_num_t **pnext_bentity1_entity2_gnum            = NULL;
    int         **pnext_bentity1_entity2_triplet         = NULL;
    int         **pnext_bentity1_entity2_interface_n     = NULL;
    int         **pnext_bentity1_entity2_interface_tot_n = NULL;
    int         **pnext_bentity1_entity2_interface       = NULL;

    for(int i_part = 0; i_part < n_part_tot; ++i_part) {
      PDM_log_trace_array_int (bentity1_entity2_n              [i_part], pn_entity1[i_part]    , "bentity1_entity2_n               ::");
      PDM_log_trace_array_int (bentity1_entity2_interface_tot_n[i_part], pn_entity1[i_part]    , "bentity1_entity2_interface_tot_n ::");
    }
    exchange_and_concat_part_graph(ptp,
                                   n_part_tot,
                                   pn_entity1,
                                   part1_to_part2_idx,
                                   part1_to_part2_interface,
                                   bentity1_entity2_n,
                                   bentity1_entity2_gnum,
                                   bentity1_entity2_triplet,
                                   bentity1_entity2_interface_n,
                                   bentity1_entity2_interface_tot_n,
                                   bentity1_entity2_interface,
                                   &pnext_bentity1_entity2_n,
                                   &pnext_bentity1_entity2_gnum,
                                   &pnext_bentity1_entity2_triplet,
                                   &pnext_bentity1_entity2_interface_n,
                                   &pnext_bentity1_entity2_interface_tot_n,
                                   &pnext_bentity1_entity2_interface);

    for(int i_part = 0; i_part < n_part_tot; ++i_part) {

      // Il faut rajouter l'interface number
      _concatenate_part_graph(pn_entity1[i_part],
                              part1_to_part2_idx                          [i_part],
                              part1_to_part2_interface                    [i_part],
                              prev_concat_bentity1_entity2_n              [i_part],
                              prev_concat_bentity1_entity2_gnum           [i_part],
                              prev_concat_bentity1_entity2_triplet        [i_part],
                              prev_concat_bentity1_entity2_interface_n    [i_part],
                              prev_concat_bentity1_entity2_interface_tot_n[i_part],
                              prev_concat_bentity1_entity2_interface      [i_part],
                              pnext_bentity1_entity2_n                    [i_part],
                              pnext_bentity1_entity2_gnum                 [i_part],
                              pnext_bentity1_entity2_triplet              [i_part],
                              pnext_bentity1_entity2_interface_n          [i_part],
                              pnext_bentity1_entity2_interface_tot_n      [i_part],
                              pnext_bentity1_entity2_interface            [i_part],
                              &concat_bentity1_entity2_n                  [i_part],
                              &concat_bentity1_entity2_gnum               [i_part],
                              &concat_bentity1_entity2_triplet            [i_part],
                              &concat_bentity1_entity2_interface_n        [i_part],
                              &concat_bentity1_entity2_interface_tot_n    [i_part],
                              &concat_bentity1_entity2_interface          [i_part]);

      int *update_bentity1_entity2_n               = NULL;
      int *update_bentity1_entity2_interface_tot_n = NULL;
      _update_in_part2_layout(pn_entity1                            [i_part],
                              part1_to_part2_idx                    [i_part],
                              pnext_bentity1_entity2_n              [i_part],
                              pnext_bentity1_entity2_interface_tot_n[i_part],
                              &update_bentity1_entity2_n,
                              &update_bentity1_entity2_interface_tot_n);

      free(pnext_bentity1_entity2_n              [i_part]);
      free(pnext_bentity1_entity2_interface_tot_n[i_part]);
      pnext_bentity1_entity2_n              [i_part] = update_bentity1_entity2_n;
      pnext_bentity1_entity2_interface_tot_n[i_part] = update_bentity1_entity2_interface_tot_n;

      /*
       * Preapre for next step
       */
      free(prev_concat_bentity1_entity2_n              [i_part]);
      free(prev_concat_bentity1_entity2_gnum           [i_part]);
      free(prev_concat_bentity1_entity2_triplet        [i_part]);
      free(prev_concat_bentity1_entity2_interface_n    [i_part]);
      free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
      free(prev_concat_bentity1_entity2_interface      [i_part]);

      prev_concat_bentity1_entity2_n              [i_part] = concat_bentity1_entity2_n              [i_part];
      prev_concat_bentity1_entity2_gnum           [i_part] = concat_bentity1_entity2_gnum           [i_part];
      prev_concat_bentity1_entity2_triplet        [i_part] = concat_bentity1_entity2_triplet        [i_part];
      prev_concat_bentity1_entity2_interface_n    [i_part] = concat_bentity1_entity2_interface_n    [i_part];
      prev_concat_bentity1_entity2_interface_tot_n[i_part] = concat_bentity1_entity2_interface_tot_n[i_part];
      prev_concat_bentity1_entity2_interface      [i_part] = concat_bentity1_entity2_interface      [i_part];

      concat_bentity1_entity2_n              [i_part] = NULL;
      concat_bentity1_entity2_gnum           [i_part] = NULL;
      concat_bentity1_entity2_triplet        [i_part] = NULL;
      concat_bentity1_entity2_interface_n    [i_part] = NULL;
      concat_bentity1_entity2_interface_tot_n[i_part] = NULL;
      concat_bentity1_entity2_interface      [i_part] = NULL;

      if(i_step > 0) {
        free(bentity1_entity2_n              [i_part]);
        free(bentity1_entity2_gnum           [i_part]);
        free(bentity1_entity2_triplet        [i_part]);
        free(bentity1_entity2_interface_n    [i_part]);
        free(bentity1_entity2_interface_tot_n[i_part]);
        free(bentity1_entity2_interface      [i_part]);
      }

      bentity1_entity2_n              [i_part] = pnext_bentity1_entity2_n              [i_part];
      bentity1_entity2_gnum           [i_part] = pnext_bentity1_entity2_gnum           [i_part];
      bentity1_entity2_triplet        [i_part] = pnext_bentity1_entity2_triplet        [i_part];
      bentity1_entity2_interface_n    [i_part] = pnext_bentity1_entity2_interface_n    [i_part];
      bentity1_entity2_interface_tot_n[i_part] = pnext_bentity1_entity2_interface_tot_n[i_part];
      bentity1_entity2_interface      [i_part] = pnext_bentity1_entity2_interface      [i_part];

      pnext_bentity1_entity2_n              [i_part] = NULL;
      pnext_bentity1_entity2_gnum           [i_part] = NULL;
      pnext_bentity1_entity2_triplet        [i_part] = NULL;
      pnext_bentity1_entity2_interface_n    [i_part] = NULL;
      pnext_bentity1_entity2_interface_tot_n[i_part] = NULL;
      pnext_bentity1_entity2_interface      [i_part] = NULL;


      _dump_graph_info(pn_entity1[i_part],
                       prev_concat_bentity1_entity2_n              [i_part],
                       prev_concat_bentity1_entity2_gnum           [i_part],
                       prev_concat_bentity1_entity2_triplet        [i_part],
                       prev_concat_bentity1_entity2_interface_n    [i_part],
                       prev_concat_bentity1_entity2_interface_tot_n[i_part],
                       prev_concat_bentity1_entity2_interface      [i_part]);

    }

    /*
     * Swap ptr
     */

    free(pnext_bentity1_entity2_n              );
    free(pnext_bentity1_entity2_gnum           );
    free(pnext_bentity1_entity2_triplet        );
    free(pnext_bentity1_entity2_interface_n    );
    free(pnext_bentity1_entity2_interface_tot_n);
    free(pnext_bentity1_entity2_interface      );

  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(prev_concat_bentity1_entity2_n              [i_part]);
    free(prev_concat_bentity1_entity2_gnum           [i_part]);
    free(prev_concat_bentity1_entity2_triplet        [i_part]);
    free(prev_concat_bentity1_entity2_interface_n    [i_part]);
    free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
    free(prev_concat_bentity1_entity2_interface      [i_part]);
  }
  free(prev_concat_bentity1_entity2_n              );
  free(prev_concat_bentity1_entity2_gnum           );
  free(prev_concat_bentity1_entity2_triplet        );
  free(prev_concat_bentity1_entity2_interface_n    );
  free(prev_concat_bentity1_entity2_interface_tot_n);
  free(prev_concat_bentity1_entity2_interface      );


  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   free(concat_bentity1_entity2_n              [i_part]);
  //   free(concat_bentity1_entity2_gnum           [i_part]);
  //   free(concat_bentity1_entity2_triplet        [i_part]);
  //   free(concat_bentity1_entity2_interface_n    [i_part]);
  //   free(concat_bentity1_entity2_interface_tot_n[i_part]);
  //   free(concat_bentity1_entity2_interface      [i_part]);
  // }
  free(concat_bentity1_entity2_n              );
  free(concat_bentity1_entity2_gnum           );
  free(concat_bentity1_entity2_triplet        );
  free(concat_bentity1_entity2_interface_n    );
  free(concat_bentity1_entity2_interface_tot_n);
  free(concat_bentity1_entity2_interface      );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(bentity1_entity2_gnum           [i_part]);
    free(bentity1_entity2_n              [i_part]);
    free(bentity1_entity2_interface_n    [i_part]);
    free(bentity1_entity2_interface      [i_part]);
    free(bentity1_entity2_triplet        [i_part]);
    free(bentity1_entity2_interface_tot_n[i_part]);
  }
  free(bentity1_entity2_gnum           );
  free(bentity1_entity2_n              );
  free(bentity1_entity2_interface_n    );
  free(bentity1_entity2_interface      );
  free(bentity1_entity2_triplet        );
  free(bentity1_entity2_interface_tot_n);

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_part_extension_build_entity1_graph
(
  PDM_part_domain_interface_t   *pdi,
  PDM_bound_type_t               entity1_bound,
  int                            n_domain,
  int                           *n_part,
  int                          **pn_entity1_in,
  PDM_g_num_t                 ***pentity1_ln_to_gn_in,
  int                         ***pentity1_hint_in,
  int                         ***pentity1_extented_to_pentity1_idx_out,
  int                         ***pentity1_extented_to_pentity1_triplet_out,
  int                         ***pentity1_extented_to_pentity1_interface_out,
  PDM_MPI_Comm                   comm
)
{

  int n_part_tot = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_part_tot += n_part[i_domain];
  }

  int          *pn_entity1            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_hint         = NULL;

  if(pentity1_hint_in != NULL) {
    pentity1_hint = malloc(n_part_tot * sizeof(int         *));
  }

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot] = pentity1_ln_to_gn_in   [i_dom][i_part];

      if(pentity1_hint_in != NULL && pentity1_hint_in[i_dom] != NULL) {
        pentity1_hint[ln_part_tot] = pentity1_hint_in[i_dom][i_part];
      }

      ln_part_tot += 1;
    }
  }

  int          *pn_entity1_num             = NULL;
  int         **pentity1_num               = NULL;
  int         **pentity1_opp_location_idx  = NULL;
  int         **pentity1_opp_location      = NULL;
  int         **pentity1_opp_interface_idx = NULL;
  int         **pentity1_opp_interface     = NULL;
  int         **pentity1_opp_sens          = NULL;
  PDM_g_num_t **pentity1_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity1_bound,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         &pn_entity1_num,
                                         &pentity1_num,
                                         &pentity1_opp_location_idx,
                                         &pentity1_opp_location,
                                         &pentity1_opp_interface_idx,
                                         &pentity1_opp_interface,
                                         &pentity1_opp_sens,
                                         &pentity1_opp_gnum);

  /*
   * Prepare creation of part_to_part
   */
  int         **part1_to_part2_idx              = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_triplet          = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_interface        = malloc(n_part_tot * sizeof(int         *));

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(comm, n_part_tot);

  int **ppart_entity1_proc_idx = NULL;
  int **ppart_entity1_part_idx = NULL;
  int **ppart_entity1          = NULL; // (i_entity1, i_proc, i_part, i_entity1_opp)
  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      NULL,
                                      n_part_tot,
                                      pn_entity1,
               (const PDM_g_num_t **) pentity1_ln_to_gn,
               (const int         **) pentity1_hint,
                                      &ppart_entity1_proc_idx,
                                      &ppart_entity1_part_idx,
                                      &ppart_entity1,
                                      NULL);

  free(part_distribution);

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);
  int n_g_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_g_part_tot += n_part_g[i_dom];
  }

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = (int *) malloc( (n_rank) * sizeof(int));
    PDM_MPI_Allgather(&li_part,
                      1,
                      PDM_MPI_INT,
                      n_part_shift,
                      1,
                      PDM_MPI_INT,
                      comm);

    if(0 == 1) {
      PDM_log_trace_array_int(n_part_shift, n_rank+1, "n_part_shift ::");
    }

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      part1_to_part2_idx[li_part] = malloc((pn_entity1[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pn_entity1[li_part]);

      int n_entity_bound = ppart_entity1_part_idx[i_part][n_g_part_tot];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_entity1[i_part][4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        int n_opp = pentity1_opp_location_idx[li_part][idx_entity+1] - pentity1_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }


      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pn_entity1[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      // printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_entity1[i_part][4*idx_entity]-1;
        int i_proc_opp   = ppart_entity1[i_part][4*idx_entity+1];
        int i_part_opp   = ppart_entity1[i_part][4*idx_entity+2]-1;
        int i_entity_opp = ppart_entity1[i_part][4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp + n_part_shift[i_proc_opp];
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        for(int idx_opp = pentity1_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pentity1_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pentity1_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pentity1_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pentity1_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pentity1_opp_interface[li_part][  idx_opp  ];
        }
      }

      free(part1_to_part2_n);

      li_part += 1;
    }

    free(n_part_shift);
  }

  free(n_part_g);

  *pentity1_extented_to_pentity1_idx_out       = part1_to_part2_idx;
  *pentity1_extented_to_pentity1_triplet_out   = part1_to_part2_triplet;
  *pentity1_extented_to_pentity1_interface_out = part1_to_part2_interface;

  free(pn_entity1          );
  free(pentity1_ln_to_gn   );

  if(pentity1_hint_in != NULL) {
    free(pentity1_hint);
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity1_num             [i_part]);
    free(pentity1_opp_location_idx[i_part]);
    free(pentity1_opp_location    [i_part]);
    free(pentity1_opp_interface   [i_part]);
    free(pentity1_opp_sens        [i_part]);
    free(pentity1_opp_gnum        [i_part]);
    free(ppart_entity1            [i_part]);
    free(ppart_entity1_proc_idx   [i_part]);
    free(ppart_entity1_part_idx   [i_part]);
  }
  free(pn_entity1_num            );
  free(pentity1_num              );
  free(pentity1_opp_location_idx );
  free(pentity1_opp_location     );
  free(pentity1_opp_interface_idx);
  free(pentity1_opp_interface    );
  free(pentity1_opp_sens         );
  free(pentity1_opp_gnum         );
  free(ppart_entity1             );
  free(ppart_entity1_proc_idx    );
  free(ppart_entity1_part_idx    );


}


/*
 * Translate and post-treated link by interface with on entity to another
 *  Exemple of use :
 *     - entity1 to face
 *     - entity1 to edge
 *     - face to cell
 */
void
PDM_part_extension_entity1_to_entity2
(
  PDM_g_num_t                   shift_by_domain_entity2,
  int                           n_part,
  int                          *pn_entity1,
  PDM_g_num_t                 **pentity1_ln_to_gn,
  int                         **pentity1_to_pentity1_idx,
  int                         **pentity1_to_pentity1_triplet,
  int                         **pentity1_to_pentity1_interface,
  int                          *pn_entity2,
  PDM_g_num_t                 **pentity2_ln_to_gn,
  int                         **pentity2_alrdy_sent,
  PDM_g_num_t                 **pentity2_ancstr,
  int                         **pentity2_path_itrf_idx,
  int                         **pentity2_path_itrf,
  int                         **pentity2_entity1_idx,
  int                         **pentity2_entity1,
  int                           prev_dentity2_itrf_n_blk,
  PDM_g_num_t                  *prev_dentity2_itrf_blk_gnum,
  int                          *prev_dentity2_itrf_gnum_and_itrf_strid,
  PDM_g_num_t                  *prev_dentity2_itrf_gnum_and_itrf_data,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_alrdy_sent_out,
  PDM_g_num_t                ***pentity2_extented_ancstr_out,
  int                        ***pentity2_extented_path_itrf_idx_out,
  int                        ***pentity2_extented_path_itrf_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  int                          *next_dentity2_itrf_n_blk_out,
  PDM_g_num_t                 **next_dentity2_itrf_blk_gnum_out,
  int                         **next_dentity2_itrf_gnum_and_itrf_strid_out,
  PDM_g_num_t                 **next_dentity2_itrf_gnum_and_itrf_data_out,
  PDM_MPI_Comm                  comm
)
{
  int debug = 1;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int         **pentity1_entity2_idx  = NULL;
  int         **pentity1_entity2      = NULL;
  PDM_part_connectivity_transpose(n_part,
                                  pn_entity2,
                                  pn_entity1,
                                  pentity2_entity1_idx,
                                  pentity2_entity1,
                                  &pentity1_entity2_idx,
                                  &pentity1_entity2);

  if(1 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      PDM_log_trace_array_long(pentity1_ln_to_gn[i_part], pn_entity1[i_part], "pentity1_ln_to_gn ::");

      log_trace("\n");
      PDM_log_trace_array_long(pentity2_ln_to_gn     [i_part],                                pn_entity2[i_part]  , "pentity2_ln_to_gn      ::");
      PDM_log_trace_array_int (pentity2_alrdy_sent   [i_part],                                pn_entity2[i_part]  , "pentity2_alrdy_sent    ::");
      PDM_log_trace_array_long(pentity2_ancstr       [i_part],                                pn_entity2[i_part]  , "pentity2_ancstr        ::");
      PDM_log_trace_array_int (pentity2_path_itrf_idx[i_part],                                pn_entity2[i_part]+1, "pentity2_path_itrf_idx ::");
      PDM_log_trace_array_int (pentity2_path_itrf    [i_part], pentity2_path_itrf_idx[i_part][pn_entity2[i_part]] , "pentity2_path_itrf     ::");

      log_trace("\n");
      PDM_log_trace_array_int(pentity1_to_pentity1_idx      [i_part],                                  pn_entity1[i_part]+1 , "pentity1_to_pentity1_idx       ::");
      PDM_log_trace_array_int(pentity1_to_pentity1_triplet  [i_part], pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]  , "pentity1_to_pentity1_triplet   ::");
      PDM_log_trace_array_int(pentity1_to_pentity1_interface[i_part], pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3, "pentity1_to_pentity1_interface ::");
    }
  }

  /*
   * Create part_to_part to exchange all data in opposit part
   * Exchange informations:
   *   - face gnum
   *   - face triplet
   *   - face number of ancestors
   *   - interface path which created faces
   *   - face ancestor (original gnum)
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------\n");
    log_trace("PTP to exchange entities2\n");
  }
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_ln_to_gn,
                                                                      (const int          *) pn_entity1,
                                                                      n_part,
                                                                      (const int          *) pn_entity1,
                                                                      n_part,
                                                                      (const int         **) pentity1_to_pentity1_idx,
                                                                      (const int         **) NULL,
                                                                      (const int         **) pentity1_to_pentity1_triplet,
                                                                      comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  if (1 ==0) {
    log_trace("\n");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int (ref_lnum2          [i_part], n_ref_lnum2[i_part]                             , "ref_lnum2           ::");
      PDM_log_trace_array_int (gnum1_come_from_idx[i_part], n_ref_lnum2[i_part]+1                           , "gnum1_come_from_idx ::");
      PDM_log_trace_array_long(gnum1_come_from    [i_part], gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]], "gnum1_come_from     ::");
    }
  }

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_triplet_n      = malloc(n_part * sizeof(int         *));
  int         **gnum1_com_from_triplet_send   = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_gnum_send      = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **gnum1_com_from_gnum_ancstr    = malloc(n_part * sizeof(PDM_g_num_t *));
  int         **gnum1_com_from_triplet_n_itrf = malloc(n_part * sizeof(int         *));
  int         **gnum1_com_from_gnum_n_itrf    = malloc(n_part * sizeof(int         *));
  int         **gnum1_com_from_gnum_itrf      = malloc(n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx     = pentity1_entity2_idx    [i_part];

    /* Count while excluding faces:
     *   - that have already been sent at previous step
     */
    int interf_entity_1 = 0;
    int n_send_part2 = 0;
    int n_send_ancestors = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        interf_entity_1 = pentity1_to_pentity1_interface[i_part][k];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          int i_entity2 = PDM_ABS(pentity1_entity2[i_part][idx_entity2])-1;
          if (pentity2_alrdy_sent[i_part][i_entity2]==0) {
            n_send_part2 += 1;
            n_send_ancestors += pentity2_path_itrf_idx[i_part][i_entity2+1]-pentity2_path_itrf_idx[i_part][i_entity2];
          }
        }
      }
    }

    /* Allocate */
    gnum1_com_from_triplet_n     [i_part] = malloc(    n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_triplet_send  [i_part] = malloc(3 * n_send_part2      * sizeof(int        ));
    gnum1_com_from_gnum_send     [i_part] = malloc(    n_send_part2      * sizeof(PDM_g_num_t));
    gnum1_com_from_gnum_ancstr   [i_part] = malloc(    n_send_part2      * sizeof(PDM_g_num_t));
    gnum1_com_from_triplet_n_itrf[i_part] = malloc(    n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_gnum_n_itrf   [i_part] = malloc(    n_send_part2      * sizeof(int        ));
    gnum1_com_from_gnum_itrf     [i_part] = malloc(    n_send_ancestors  * sizeof(int        ));

    int         *_gnum1_com_from_triplet_n      = gnum1_com_from_triplet_n     [i_part];
    int         *_gnum1_com_from_triplet_send   = gnum1_com_from_triplet_send  [i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_send      = gnum1_com_from_gnum_send     [i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_ancstr    = gnum1_com_from_gnum_ancstr   [i_part];
    int         *_gnum1_com_from_triplet_n_itrf = gnum1_com_from_triplet_n_itrf[i_part];
    int         *_gnum1_com_from_gnum_n_itrf    = gnum1_com_from_gnum_n_itrf   [i_part];
    int         *_gnum1_com_from_gnum_itrf      = gnum1_com_from_gnum_itrf     [i_part];

    int         *pentity2_alrdy_sent_tmp      = malloc(    pn_entity2[i_part]* sizeof(int));
    for (int i_entity2=0; i_entity2<pn_entity2[i_part]; ++i_entity2) {
      pentity2_alrdy_sent_tmp[i_entity2] = pentity2_alrdy_sent[i_part][i_entity2];
    }


    int l_n_send_part2     = 0;
    int l_n_send_ancestors = 0;
    n_send_ancestors = 0;
    n_send_part2     = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) { // Loop on referenced entities of part2
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) { // Loop on entities of part1 referenced by entity2
        interf_entity_1 = pentity1_to_pentity1_interface[i_part][k];
        
        // log_trace("    via interf_entity1 = %d, with entity = %d\n", interf_entity_1, gnum1_come_from[i_part][k]);
        l_n_send_part2     = 0;
        l_n_send_ancestors = 0;
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          // boucle sur les faces incidente aux vertex part2
          int i_entity2 = PDM_ABS(pentity1_entity2[i_part][idx_entity2])-1;
          // log_trace("        get i_entity2 = %d : interf_entity2 = %d\n", i_entity2, pentity2_path_itrf[i_part][i_entity2]);
          // log_trace("        get i_entity2 = %d\n", i_entity2);
          // log_trace("             alrdy_sent = %d\n", pentity2_alrdy_sent[i_part][i_entity2]);

          // assert (pentity2_alrdy_sent[i_part][i_entity2]==0);
          if (pentity2_alrdy_sent[i_part][i_entity2]==0) {
            // log_trace("            SAVING i_entity2 = (%d %d %d)\n", i_rank, i_part, i_entity2);
            _gnum1_com_from_gnum_send    [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
            _gnum1_com_from_gnum_ancstr  [  n_send_part2  ] = pentity2_ancstr  [i_part][i_entity2];
            _gnum1_com_from_triplet_send [3*n_send_part2  ] = i_rank;
            _gnum1_com_from_triplet_send [3*n_send_part2+1] = i_part;
            _gnum1_com_from_triplet_send [3*n_send_part2+2] = i_entity2;

            _gnum1_com_from_gnum_n_itrf  [  n_send_part2  ] = pentity2_path_itrf_idx[i_part][i_entity2+1]-pentity2_path_itrf_idx[i_part][i_entity2];
            for (int i_ancestor=pentity2_path_itrf_idx[i_part][i_entity2]; i_ancestor<pentity2_path_itrf_idx[i_part][i_entity2+1]; ++i_ancestor) {
              _gnum1_com_from_gnum_itrf[n_send_ancestors] = pentity2_path_itrf[i_part][i_ancestor];
              // log_trace("                i_ancstr = %d --> itrf = %d ; ancstr = %d\n", i_ancestor, pentity2_path_itrf[i_part][i_ancestor], pentity2_ancstr[i_part][i_ancestor]);
              l_n_send_ancestors++;
              n_send_ancestors++;
            }
            pentity2_alrdy_sent_tmp[i_entity2] = 1 ;
            l_n_send_part2++;
            n_send_part2++;
            // }
          }
        }
        _gnum1_com_from_triplet_n     [k] = l_n_send_part2;
        _gnum1_com_from_triplet_n_itrf[k] = l_n_send_ancestors;

      }
    }

    for (int i_entity2=0; i_entity2<pn_entity2[i_part]; ++i_entity2) {
      pentity2_alrdy_sent[i_part][i_entity2] = pentity2_alrdy_sent_tmp[i_entity2];
    }
    free(pentity2_alrdy_sent_tmp);

    if(1 == 1) {
      // Ce tableau contient pour chaque vertex d'un interface, ses faces associées.
      log_trace("\n");
      PDM_log_trace_array_int (_gnum1_com_from_triplet_n    , n_ref_lnum2[i_part], "_gnum1_com_from_triplet_n    ::");
      PDM_log_trace_array_long(_gnum1_com_from_gnum_send    , n_send_part2       , "_gnum1_com_from_gnum_send    ::");
      PDM_log_trace_array_long(_gnum1_com_from_gnum_ancstr  , n_send_part2       , "_gnum1_com_from_gnum_ancstr  ::");
      PDM_log_trace_array_int (_gnum1_com_from_triplet_send , n_send_part2*3     , "_gnum1_com_from_triplet_send ::");
      log_trace("\n");
      PDM_log_trace_array_int (_gnum1_com_from_triplet_n_itrf , n_ref_lnum2[i_part], "_gnum1_com_from_triplet_n_itrf  ::");
      PDM_log_trace_array_int (_gnum1_com_from_gnum_n_itrf    , n_send_part2       , "_gnum1_com_from_gnum_n_itrf     ::");
      PDM_log_trace_array_int (_gnum1_com_from_gnum_itrf      , n_send_ancestors   , "_gnum1_com_from_gnum_itrf       ::");
    }


  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity1_entity2_idx[i_part]);
    free(pentity1_entity2    [i_part]);
  }
  free(pentity1_entity2_idx);
  free(pentity1_entity2    );


  int exch_request = -1;
  int         **pextract_entity2_n      = NULL;
  int         **pextract_entity2_itrf_n = NULL;

  PDM_g_num_t **pextract_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  PDM_g_num_t **pextract_entity2_ancstr = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_ancstr,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_ancstr,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int  **)  gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_triplet_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity2_n_itrf   = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(int),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_n_itrf  ,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_n_itrf,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity2_path_itrf = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(int),
                (const int **)   gnum1_com_from_triplet_n_itrf,
                (const void **)  gnum1_com_from_gnum_itrf,
                                 &pextract_entity2_itrf_n,
                    (void ***)   &pextract_entity2_path_itrf,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(gnum1_com_from_triplet_n     [i_part]);
    free(gnum1_com_from_triplet_send  [i_part]);
    free(gnum1_com_from_gnum_send     [i_part]);
    free(gnum1_com_from_gnum_ancstr   [i_part]);
    free(gnum1_com_from_triplet_n_itrf[i_part]);
    free(gnum1_com_from_gnum_n_itrf   [i_part]);
    free(gnum1_com_from_gnum_itrf     [i_part]);
  }
  free(gnum1_com_from_triplet_n     );
  free(gnum1_com_from_triplet_n_itrf);
  free(gnum1_com_from_triplet_send  );
  free(gnum1_com_from_gnum_send     );
  free(gnum1_com_from_gnum_n_itrf   );
  free(gnum1_com_from_gnum_itrf     );
  free(gnum1_com_from_gnum_ancstr   );

  PDM_part_to_part_free(ptp);



  
  // > Compute n_entity2 for every entity1 and n_interface for each entity2 of every entity1
  int **pextract_entity2_idx           = malloc(n_part * sizeof(int *));
  int **pextract_entity2_path_itrf_idx = malloc(n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    
    pextract_entity2_idx          [i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity2_n      [i_part], n_part1_to_part2);
    pextract_entity2_path_itrf_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity2_n_itrf [i_part], pextract_entity2_idx[i_part][n_part1_to_part2]);
    
    if ((debug==1)&&(1 == 1)) {
      PDM_log_trace_array_int(pextract_entity2_idx[i_part], n_part1_to_part2+1, "pextract_entity2_idx        ::");
      PDM_log_trace_array_int(pextract_entity2_path_itrf_idx[i_part], pextract_entity2_idx[i_part][n_part1_to_part2]+1, "pextract_entity2_path_itrf_idx ::");
    }
  }


  if ((debug==1)&&(1 == 1)) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_part1_to_part2                 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
      int n_part1_to_part2_recv_tot        =                                        pextract_entity2_idx[i_part][n_part1_to_part2];
      int n_part1_to_part2_recv_tot_ancstr = pextract_entity2_path_itrf_idx[i_part][pextract_entity2_idx[i_part][n_part1_to_part2]];

      log_trace("\n");
      PDM_log_trace_array_long(pextract_entity2_gnum     [i_part], n_part1_to_part2_recv_tot       , "_pextract_entity2_gnum     ::");
      PDM_log_trace_array_long(pextract_entity2_ancstr   [i_part], n_part1_to_part2_recv_tot       , "_pextract_entity2_ancstr   ::");
      PDM_log_trace_array_int (pextract_entity2_n_itrf   [i_part], n_part1_to_part2_recv_tot       , "_pextract_entity2_n_itrf   ::");
      PDM_log_trace_array_int (pextract_entity2_path_itrf[i_part], n_part1_to_part2_recv_tot_ancstr, "_pextract_entity2_path_itrf::");
    }


    log_trace("\n");
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_gnum           , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_gnum            ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_gnum_and_itrf_strid, prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_gnum_and_itrf_strid ::");

    int n_data = 0;
    for(int i = 0; i < prev_dentity2_itrf_n_blk; ++i) {
      n_data += prev_dentity2_itrf_gnum_and_itrf_strid[i];
    }
    PDM_log_trace_array_long(prev_dentity2_itrf_gnum_and_itrf_data , 2 * n_data              , "prev_dentity2_itrf_gnum_and_itrf_data  ::");
  }







  /***************************************************************
   * We have to check that composition (gnum, interface(s?)) doesn't already exist
   * So we ask block interface data_base about entity candidate.
   *
   *
   * Introducing pextract_entity2_kind:
   *    0/ Invalid
   *       a/ Passing through an interface that created it         :-1
   *    1/ Internal
   *       a/ Already know in current partition                    : 0
   *       b/ New in current partition                             : 2
   *    2/ From interface
   *       a/ Know by relation table and know in current partition : 1
   *       b/ Know by relation table but not local                 : 3
   *       c/ New                                                  : 4
   *
   * We try to store entities in following order : [2/3/4]
   *
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("---------------------------\n");
    log_trace("BTP to find valid entities2\n");
  }


  /*
   * Need to sort pentity2_gnum to search if a candidate entity
   * doesn't already exist in partition.
   *
   */
  PDM_g_num_t **pentity2_gnum_sorted = malloc( n_part * sizeof(PDM_g_num_t *));
  
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pentity2_gnum_sorted[i_part] = malloc( pn_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {
      pentity2_gnum_sorted[i_part][i_entity2] = pentity2_ln_to_gn[i_part][i_entity2];
    }
    PDM_sort_long(pentity2_gnum_sorted[i_part], NULL, pn_entity2[i_part]);
  }



  /*
   * Get entity2 informations from distributed data_base
   *
   */
  int *pextract_entity1_entity2_n_elmt = malloc(n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    pextract_entity1_entity2_n_elmt[i_part] = pextract_entity2_idx[i_part][n_part1_to_part2];
  }
  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(prev_dentity2_itrf_blk_gnum,
                                                                        prev_dentity2_itrf_n_blk,
                                                (const PDM_g_num_t **)  pextract_entity2_gnum,
                                                                        pextract_entity1_entity2_n_elmt,
                                                                        n_part,
                                                                        comm);

  int         **prev_pentity2_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t **prev_pentity2_itrf_gnum_and_itrf_data  = NULL;

  PDM_block_to_part_exch(btp,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         prev_dentity2_itrf_gnum_and_itrf_strid,
                         prev_dentity2_itrf_gnum_and_itrf_data,
                         &prev_pentity2_itrf_gnum_and_itrf_strid,
         (void ***)      &prev_pentity2_itrf_gnum_and_itrf_data);
  PDM_block_to_part_free(btp);


  if ((debug==1)&&(1 == 1)) {
    log_trace("\n");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_data = 0;
      for(int i = 0; i < pextract_entity1_entity2_n_elmt[i_part]; ++i) {
        n_data += prev_pentity2_itrf_gnum_and_itrf_strid[i_part][i];
      }
      PDM_log_trace_array_int(prev_pentity2_itrf_gnum_and_itrf_strid[i_part], pextract_entity1_entity2_n_elmt[i_part], "prev_pentity2_itrf_gnum_and_itrf_strid ::");
      PDM_log_trace_array_long(prev_pentity2_itrf_gnum_and_itrf_data[i_part], 2 * n_data                             , "prev_pentity2_itrf_gnum_and_itrf_data  ::");
    }
  }

  int          *pn_new_entity2        = malloc(n_part * sizeof(int          ));
  int          *pn_new_entity2_kind4  = malloc(n_part * sizeof(int          ));
  int          *pn_new_entity2_n_itrf = malloc(n_part * sizeof(int          ));
  int         **pextract_entity2_kind = malloc(n_part * sizeof(int         *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_new_entity2       [i_part] = 0;
    pn_new_entity2_kind4 [i_part] = 0;
    pn_new_entity2_n_itrf[i_part] = 0;

    int n_part1_to_part2 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = pextract_entity2_idx[i_part][n_part1_to_part2];

    // > Partitioned data
    PDM_g_num_t *_pentity2_gnum_sorted                  = pentity2_gnum_sorted                 [i_part];
    PDM_g_num_t *_pextract_entity2_gnum                 = pextract_entity2_gnum                [i_part];
    int         *_pextract_entity2_n_itrf               = pextract_entity2_n_itrf              [i_part];

    // > Partitioned data from data_base
    PDM_g_num_t *_prev_pentity2_itrf_gnum_and_itrf_data = prev_pentity2_itrf_gnum_and_itrf_data[i_part];

    pextract_entity2_kind[i_part] = malloc(n_part1_to_part2_recv_tot * sizeof(int));

    int idx_read_data = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {

      if(pentity1_to_pentity1_interface[i_part][i] != 0) {

        for(int j = pextract_entity2_idx[i_part][i]; j < pextract_entity2_idx[i_part][i+1]; ++j) {
          int cur_itrf     = PDM_ABS (pentity1_to_pentity1_interface[i_part][i]);
          int sgn_cur_itrf = PDM_SIGN(pentity1_to_pentity1_interface[i_part][i]);

          int keep = 1;
          PDM_g_num_t gnum_opp = 0;
          for(int k = 0; k < prev_pentity2_itrf_gnum_and_itrf_strid[i_part][j]; ++k) {

            int         opp_itrf     = PDM_ABS (_prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data+1]);
            int         opp_sgn_itrf = PDM_SIGN(_prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data+1]);

            if((cur_itrf==opp_itrf) && (sgn_cur_itrf==- opp_sgn_itrf)) {
              keep = 0;
              gnum_opp = _prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data];
            }

            idx_read_data++;
          }

          if(keep == 1) { // Entity has no relation in data_base, so its new
            pextract_entity2_kind[i_part][j] = 4;
            pn_new_entity2       [i_part]   += 1;
            pn_new_entity2_kind4 [i_part]   += 1;
            pn_new_entity2_n_itrf[i_part]   += _pextract_entity2_n_itrf[j]+1;
          }
          else { // Entity has relation in data_base, so part can alrdy have it, or not
            assert(gnum_opp != 0);
            int pos_int = PDM_binary_search_long(gnum_opp, _pentity2_gnum_sorted, pn_entity2[i_part]);
            if( pos_int != -1) { // Found in partition, no need to create it
              pextract_entity2_kind[i_part][j] = 1;
            } else { // Not found in partition, need to create it
              pextract_entity2_kind[i_part][j] = 3;
              pn_new_entity2       [i_part]   += 1;
              pn_new_entity2_n_itrf[i_part]   += _pextract_entity2_n_itrf[j]+1;
            }
          }
        }
      }
      else{
        for(int j = pextract_entity2_idx[i_part][i]; j < pextract_entity2_idx[i_part][i+1]; ++j) { // For joins search if entity alrdy on partition
          PDM_g_num_t gnum_cur = _pextract_entity2_gnum[j];
          int pos_int = PDM_binary_search_long(gnum_cur, _pentity2_gnum_sorted, pn_entity2[i_part]);
          if( pos_int != -1) {
            pextract_entity2_kind[i_part][j] = 0;
          } else {
            pextract_entity2_kind[i_part][j] = 2;
            pn_new_entity2       [i_part]   += 1;
            pn_new_entity2_n_itrf[i_part]   += _pextract_entity2_n_itrf[j]; // Not +1 because don't care of interface=0
          }

          for(int k = 0; k < prev_pentity2_itrf_gnum_and_itrf_strid[i_part][j]; ++k) {
            idx_read_data++;
          }
        }

      }
    }

    log_trace("pn_new_entity2        = %d\n",pn_new_entity2       [i_part]);
    log_trace("pn_new_entity2_kind4  = %d\n",pn_new_entity2_kind4 [i_part]);
    log_trace("pn_new_entity2_n_itrf = %d\n",pn_new_entity2_n_itrf[i_part]);
    PDM_log_trace_array_int(pextract_entity2_kind[i_part], pextract_entity1_entity2_n_elmt[i_part], "pextract_entity2_kind ::");
    printf("ATTTENTION L'ORDRE DANS LE TABLEAU TRIE DOIT ETRE RETROUVE DANS LA POSITION !!!!!!\n");

    free(prev_pentity2_itrf_gnum_and_itrf_strid[i_part]);
    free(prev_pentity2_itrf_gnum_and_itrf_data [i_part]);
  }

  free(prev_pentity2_itrf_gnum_and_itrf_strid);
  free(prev_pentity2_itrf_gnum_and_itrf_data );







  /*
   * At this stage :
   *   - We have the connectivity exchange
   *   - We receive face that can be already define. 
   *       --> They can be merged using ancestor face and interfaces it went through,
   *           but they must be sorted, to be sure that computed gnum is the same for faces in conflict:
   *           gnum = nuplet(ancestor, sort(interf1, interf2, ...))
   *
   * Post-treatment
   *   - Remove duplicate (same gnum or same equivalent transformation )
   *   - Keep link between entity2 (extented and current part)
   *   - Create new global numbering for all new entities
   */

  if (debug==1) {
    log_trace("\n\n");
    log_trace("------------------------------------------\n");
    log_trace("compute gnum for new entities2 from interf\n");
  }

  /*
   * Compute gnum of extracted entities2 using nuplet(ancestor, sort(interf1, interf2, ...))
   */

  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     n_part,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);


  /*
   * Count the max number of ancestors for entities generating a new gnum
   */
  int l_n_max_ancstr = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2          = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = pextract_entity2_idx[i_part][n_part1_to_part2];

    for (int i_entity=0; i_entity<pextract_entity2_idx[i_part][n_part1_to_part2]; ++i_entity) {
      if (pextract_entity2_kind[i_part][i_entity] == 4) {
        l_n_max_ancstr = PDM_MAX(pextract_entity2_n_itrf[i_part][i_entity], l_n_max_ancstr);
      }
    }
  }
  log_trace("\n-->l_n_max_ancstr = %d\n", l_n_max_ancstr);
  
  int   n_max_ancstr_t[1];
  int l_n_max_ancstr_t[1] = {l_n_max_ancstr};
  PDM_MPI_Allreduce(l_n_max_ancstr_t, n_max_ancstr_t, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  int   n_max_ancstr = n_max_ancstr_t[0];
  log_trace("\n-->n_max_ancstr = %d\n", n_max_ancstr);
  // TODO: faire un from nuplet avec stride variable
  
  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2+n_max_ancstr);


  // > Count number of entities2 got by an interface and number of valid face (faces not passing though same interface)
  PDM_g_num_t **pentity2_only_by_interface_nuplet = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pentity2_only_by_interface_parent = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int         *_pextract_entity2_idx           = pextract_entity2_idx          [i_part];
    PDM_g_num_t *_pextract_entity2_gnum          = pextract_entity2_gnum         [i_part];
    PDM_g_num_t *_pextract_entity2_ancstr        = pextract_entity2_ancstr       [i_part];
    int         *_pextract_entity2_n_itrf        = pextract_entity2_n_itrf       [i_part];
    int         *_pextract_entity2_path_itrf_idx = pextract_entity2_path_itrf_idx[i_part];
    int         *_pextract_entity2_path_itrf     = pextract_entity2_path_itrf    [i_part];
    int         *_pextract_entity2_kind          = pextract_entity2_kind         [i_part];

    int n_part1_to_part2                 =  pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot        = _pextract_entity2_idx[n_part1_to_part2];
    int n_part1_to_part2_recv_tot_ancstr = _pextract_entity2_path_itrf_idx[_pextract_entity2_idx[n_part1_to_part2]];



    // > Build nuplet for entities2 got by an interface
    int *tmp_itrf_ancstr = malloc((n_max_ancstr+1) * sizeof(int));
    // printf("NOT SURE pentity2_only_by_interface_parent IS USEFUL (info of ancestor in nuplet) or for update graph ?\n");
    pentity2_only_by_interface_nuplet[i_part] = malloc((2+n_max_ancstr)* pn_new_entity2_kind4[i_part] * sizeof(PDM_g_num_t));
    pentity2_only_by_interface_parent[i_part] = malloc( 2              * pn_new_entity2_kind4[i_part] * sizeof(PDM_g_num_t));
    PDM_g_num_t *_pentity2_only_by_interface_nuplet = pentity2_only_by_interface_nuplet[i_part];
    PDM_g_num_t *_pentity2_only_by_interface_parent = pentity2_only_by_interface_parent[i_part];

    int idx_write = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
        if (_pextract_entity2_kind[j]==4) {
          _pentity2_only_by_interface_parent[              2 *idx_write  ] = _pextract_entity2_gnum  [j];
          _pentity2_only_by_interface_parent[              2 *idx_write+1] =  pentity1_to_pentity1_interface[i_part][i];

          // > Fill with previous interface, those who has no previous or less than other, fill with current interface
          _pentity2_only_by_interface_nuplet[(n_max_ancstr+2)*idx_write  ] = _pextract_entity2_ancstr[j];
          int l_i_ancstr = 0;
          for (int i_ancstr=_pextract_entity2_path_itrf_idx[j]; i_ancstr<_pextract_entity2_path_itrf_idx[j+1]; ++i_ancstr) {
            tmp_itrf_ancstr[l_i_ancstr++] = _pextract_entity2_path_itrf[i_ancstr];
          }
          for (int i_ancstr=_pextract_entity2_n_itrf[j]; i_ancstr<n_max_ancstr+1; ++i_ancstr){
            tmp_itrf_ancstr[l_i_ancstr++] = pentity1_to_pentity1_interface[i_part][i];
          }

          // > Sort interface so that entity with same ancstr will give same gnum
          PDM_sort_int(tmp_itrf_ancstr, NULL, n_max_ancstr+1);
          
          for (int i_ancstr=0; i_ancstr<n_max_ancstr+1; ++i_ancstr){
            _pentity2_only_by_interface_nuplet[(n_max_ancstr+2)*idx_write+1+i_ancstr] = tmp_itrf_ancstr[i_ancstr];
          }
          idx_write++;
        }
      }
    }
    assert(idx_write == pn_new_entity2_kind4[i_part]);


    if(1 == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(pentity2_only_by_interface_nuplet[i_part], (n_max_ancstr+2)* pn_new_entity2_kind4 [i_part], "pentity2_only_by_interface_nuplet ::");
      PDM_log_trace_array_long(pentity2_only_by_interface_parent[i_part],               2 * pn_new_entity2_kind4 [i_part], "pentity2_only_by_interface_parent ::");
      // PDM_log_trace_array_long(pentity2_interface         [i_part],     pn_entity2_only_by_interface[i_part], "pentity2_interface          ::");
    }


    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              pn_new_entity2_kind4 [i_part],
                              pentity2_only_by_interface_nuplet[i_part]);

  }

  PDM_gnum_compute(gen_gnum_entity2);


  /* Update new entities2 with computed gnum */
  PDM_g_num_t **extented_entity2_ln_to_gn  = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pnew_entity2_gnum          = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pnew_entity2_ancstr        = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pnew_entity2_parent_t      = malloc(n_part * sizeof(PDM_g_num_t *));
  int         **pnew_entity2_n_itrf        = malloc(n_part * sizeof(int         *));
  int         **pnew_entity2_kind          = malloc(n_part * sizeof(int         *));
  int         **pnew_entity2_path_itrf_idx = malloc(n_part * sizeof(int         *));
  int         **pnew_entity2_path_itrf     = malloc(n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int         *_pextract_entity2_idx           = pextract_entity2_idx          [i_part];
    PDM_g_num_t *_pextract_entity2_gnum          = pextract_entity2_gnum         [i_part];
    PDM_g_num_t *_pextract_entity2_ancstr        = pextract_entity2_ancstr       [i_part];
    int         *_pextract_entity2_n_itrf        = pextract_entity2_n_itrf       [i_part];
    int         *_pextract_entity2_path_itrf_idx = pextract_entity2_path_itrf_idx[i_part];
    int         *_pextract_entity2_path_itrf     = pextract_entity2_path_itrf    [i_part];
    int         *_pextract_entity2_kind          = pextract_entity2_kind         [i_part];

    int n_part1_to_part2          =  pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];
    int n_new_entity2             =  pn_new_entity2       [i_part];
    int n_new_entity2_n_itrf      =  pn_new_entity2_n_itrf[i_part];
    // log_trace(" n_valid_entity2         = %d\n",  n_valid_entity2        );
    // log_trace(" n_valid_entity2_n_itrf  = %d\n",  n_valid_entity2_n_itrf );
    
    // log_trace("Looking for size = %d\n", _pextract_entity2_path_itrf_idx[n_part1_to_part2_recv_tot]);
    // int n_itrf_tot = _pextract_entity2_path_itrf_idx[n_part1_to_part2_recv_tot]+n_part1_to_part2_recv_tot;
    // log_trace("n_itrf_tot = %d\n", n_itrf_tot);

    extented_entity2_ln_to_gn [i_part] = PDM_gnum_get(gen_gnum_entity2, i_part);
    pnew_entity2_gnum         [i_part] = malloc( n_new_entity2       * sizeof(PDM_g_num_t *));
    pnew_entity2_ancstr       [i_part] = malloc( n_new_entity2       * sizeof(PDM_g_num_t *));
    pnew_entity2_parent_t     [i_part] = malloc( n_new_entity2*3     * sizeof(PDM_g_num_t *));
    pnew_entity2_n_itrf       [i_part] = malloc( n_new_entity2       * sizeof(int         *));
    pnew_entity2_kind         [i_part] = malloc( n_new_entity2       * sizeof(int         *));
    pnew_entity2_path_itrf_idx[i_part] = malloc((n_new_entity2+1)    * sizeof(int         *));
    pnew_entity2_path_itrf    [i_part] = malloc( n_new_entity2_n_itrf* sizeof(int         *));
    PDM_g_num_t *_pnew_entity2_gnum          = pnew_entity2_gnum    [i_part];
    PDM_g_num_t *_pnew_entity2_ancstr        = pnew_entity2_ancstr  [i_part];
    PDM_g_num_t *_pnew_entity2_parent_t      = pnew_entity2_parent_t[i_part];
    int         *_pnew_entity2_n_itrf        = pnew_entity2_n_itrf  [i_part];
    int         *_pnew_entity2_kind          = pnew_entity2_kind    [i_part];
    int         *_pnew_entity2_path_itrf     = pnew_entity2_path_itrf    [i_part];
    int         *_pnew_entity2_path_itrf_idx = pnew_entity2_path_itrf_idx[i_part];

    int idx_read = 0;
    int idx_write = 0;
    int interf_entity_1 = 0;
    _pnew_entity2_path_itrf_idx[0] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      interf_entity_1 = pentity1_to_pentity1_interface[i_part][i];
      if(interf_entity_1 != 0) {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {

          if (_pextract_entity2_kind[j]==4) {
            _pnew_entity2_path_itrf_idx[  idx_write+1] = _pnew_entity2_path_itrf_idx[idx_write]+_pextract_entity2_n_itrf  [j]+1;
            _pnew_entity2_gnum         [  idx_write  ] =  extented_entity2_ln_to_gn[i_part][idx_read++] + shift_by_domain_entity2;
            _pnew_entity2_ancstr       [  idx_write  ] = _pextract_entity2_ancstr  [j];
            _pnew_entity2_kind         [  idx_write  ] = _pextract_entity2_kind    [j];
            _pnew_entity2_parent_t     [3*idx_write  ] =  pextract_entity2_triplet[i_part][3*j  ];
            _pnew_entity2_parent_t     [3*idx_write+1] =  pextract_entity2_triplet[i_part][3*j+1];
            _pnew_entity2_parent_t     [3*idx_write+2] =  pextract_entity2_triplet[i_part][3*j+2];
            _pnew_entity2_n_itrf       [  idx_write  ] = _pextract_entity2_n_itrf  [j]+1;
            int l_i_ancstr = _pnew_entity2_path_itrf_idx[idx_write];
            for (int i_ancstr=_pextract_entity2_path_itrf_idx[j]; i_ancstr<_pextract_entity2_path_itrf_idx[j+1]; ++i_ancstr) {
              _pnew_entity2_path_itrf[l_i_ancstr++] = _pextract_entity2_path_itrf[i_ancstr];
            }
            _pnew_entity2_path_itrf[l_i_ancstr] = interf_entity_1;
            idx_write++;
          } else if (_pextract_entity2_kind[j]==3) {
            _pnew_entity2_path_itrf_idx[  idx_write+1] = _pnew_entity2_path_itrf_idx[idx_write]+_pextract_entity2_n_itrf  [j]+1;
            _pnew_entity2_gnum         [  idx_write  ] = _pextract_entity2_gnum    [j];
            _pnew_entity2_ancstr       [  idx_write  ] = _pextract_entity2_ancstr  [j];
            _pnew_entity2_kind         [  idx_write  ] = _pextract_entity2_kind    [j];
            _pnew_entity2_parent_t     [3*idx_write  ] =  pextract_entity2_triplet[i_part][3*j  ];
            _pnew_entity2_parent_t     [3*idx_write+1] =  pextract_entity2_triplet[i_part][3*j+1];
            _pnew_entity2_parent_t     [3*idx_write+2] =  pextract_entity2_triplet[i_part][3*j+2];
            _pnew_entity2_n_itrf       [  idx_write  ] = _pextract_entity2_n_itrf  [j]+1;
            int l_i_ancstr = _pnew_entity2_path_itrf_idx[idx_write];
            for (int i_ancstr=_pextract_entity2_path_itrf_idx[j]; i_ancstr<_pextract_entity2_path_itrf_idx[j+1]; ++i_ancstr) {
              _pnew_entity2_path_itrf[l_i_ancstr++] = _pextract_entity2_path_itrf[i_ancstr];
            }
            _pnew_entity2_path_itrf[l_i_ancstr] = interf_entity_1;
            idx_write++;
          }
        }
      }
      else {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
          PDM_g_num_t gnum_cur = _pextract_entity2_ancstr[j];

          if (_pextract_entity2_kind[j]==2) {
            _pnew_entity2_path_itrf_idx[  idx_write+1] = _pnew_entity2_path_itrf_idx[idx_write]+_pextract_entity2_n_itrf  [j];
            _pnew_entity2_gnum         [  idx_write  ] = _pextract_entity2_gnum    [j];
            _pnew_entity2_ancstr       [  idx_write  ] = _pextract_entity2_ancstr  [j];
            _pnew_entity2_kind         [  idx_write  ] = _pextract_entity2_kind    [j];
            _pnew_entity2_parent_t     [3*idx_write  ] =  pextract_entity2_triplet[i_part][3*j  ];
            _pnew_entity2_parent_t     [3*idx_write+1] =  pextract_entity2_triplet[i_part][3*j+1];
            _pnew_entity2_parent_t     [3*idx_write+2] =  pextract_entity2_triplet[i_part][3*j+2];
            _pnew_entity2_n_itrf       [  idx_write  ] = _pextract_entity2_n_itrf  [j];
            int l_i_ancstr = _pnew_entity2_path_itrf_idx[idx_write];
            for (int i_ancstr=_pextract_entity2_path_itrf_idx[j]; i_ancstr<_pextract_entity2_path_itrf_idx[j+1]; ++i_ancstr) {
              _pnew_entity2_path_itrf[l_i_ancstr++] = _pextract_entity2_path_itrf[i_ancstr];
            }
            idx_write++;
          }
        }
      }
    }

    if(1 == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(_pnew_entity2_gnum    , n_new_entity2        , "_pnew_entity2_gnum    (Update) :: ");
      PDM_log_trace_array_long(_pnew_entity2_ancstr  , n_new_entity2        , "_pnew_entity2_ancstr  (Update) :: ");
      PDM_log_trace_array_long(_pnew_entity2_parent_t, n_new_entity2*3      , "_pnew_entity2_ancstr_t(Update) :: ");
      PDM_log_trace_array_int (_pnew_entity2_n_itrf  , n_new_entity2        , "_pnew_entity2_n_itrf  (Update) :: ");
      PDM_log_trace_array_int (_pnew_entity2_kind    , n_new_entity2        , "_pnew_entity2_kind    (Update) :: ");
      PDM_log_trace_array_int (_pnew_entity2_path_itrf_idx, n_new_entity2+1      , "_pnew_entity2_path_itrf_idx(Update) :: ");
      PDM_log_trace_array_int (_pnew_entity2_path_itrf    , n_new_entity2_n_itrf , "_pnew_entity2_path_itrf    (Update) :: ");
    }
  }


  /*
   * Create part_to_block to keep link between the new extent entity2 link with origin (CAUTION TO SHIFT) :
   *   - Reminder : The relation is bijective
   *
   *
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------------------\n");
    log_trace("complete data_base with new entities2\n");
  }

  /*
   *  Creation d'une table à partir de pentity2_only_by_interface_parent et extented_entity2_ln_to_gn (Attention au shift)
   */
  PDM_g_num_t **gnum_itrf_link                     = malloc((n_part+1) * sizeof(PDM_g_num_t *));
  PDM_g_num_t **gnum_and_itrf_data                 = malloc((n_part+1) * sizeof(PDM_g_num_t *));
  int         **gnum_and_itrf_stri                 = malloc((n_part+1) * sizeof(int         *));
  double      **gnum_and_itrf_weight               = malloc((n_part+1) * sizeof(double      *));
  int          *pn_entity2_only_by_interface_twice = malloc((n_part+1) * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    PDM_g_num_t* _pentity2_only_by_interface_parent = pentity2_only_by_interface_parent[i_part];
    PDM_g_num_t* _pnew_entity2_gnum  = pnew_entity2_gnum [i_part];

    PDM_g_num_t *_pextract_entity2_gnum    = pextract_entity2_gnum   [i_part];
    int         *_pextract_entity2_idx     = pextract_entity2_idx    [i_part];

    int _pn_new_entity2_kind4 = pn_new_entity2_kind4[i_part];
    pn_entity2_only_by_interface_twice[i_part] = 2 * pn_new_entity2_kind4 [i_part];

    // Symétrise information
    gnum_itrf_link      [i_part] = malloc(2 * _pn_new_entity2_kind4 * sizeof(PDM_g_num_t));
    gnum_and_itrf_data  [i_part] = malloc(4 * _pn_new_entity2_kind4 * sizeof(PDM_g_num_t));
    gnum_and_itrf_stri  [i_part] = malloc(2 * _pn_new_entity2_kind4 * sizeof(int        ));
    gnum_and_itrf_weight[i_part] = malloc(2 * _pn_new_entity2_kind4 * sizeof(double     ));

    for(int i = 0; i < _pn_new_entity2_kind4; ++i) {

      PDM_g_num_t orig_gnum = _pentity2_only_by_interface_parent[2*i  ];
      PDM_g_num_t i_itrf    = _pentity2_only_by_interface_parent[2*i+1];

      PDM_g_num_t new_gnum = extented_entity2_ln_to_gn[i_part][i]+ shift_by_domain_entity2;;

      int j = _pn_new_entity2_kind4+i;

      /* Remplissage */
      gnum_itrf_link[i_part][i] = orig_gnum;
      gnum_itrf_link[i_part][j] = new_gnum;

      gnum_and_itrf_data[i_part][2*i  ] = new_gnum;
      gnum_and_itrf_data[i_part][2*i+1] = -_pentity2_only_by_interface_parent[2*i+1];

      gnum_and_itrf_data[i_part][2*j  ] =  orig_gnum;
      gnum_and_itrf_data[i_part][2*j+1] = _pentity2_only_by_interface_parent[2*i+1];

      gnum_and_itrf_stri[i_part][i] = 1;
      gnum_and_itrf_stri[i_part][j] = 1;

      gnum_and_itrf_weight[i_part][i] = 1.;
      gnum_and_itrf_weight[i_part][j] = 1.;

    }
  }

  /*
   * Merge with previous
   */
  gnum_itrf_link                    [n_part] = prev_dentity2_itrf_blk_gnum;
  gnum_and_itrf_data                [n_part] = prev_dentity2_itrf_gnum_and_itrf_data;
  gnum_and_itrf_stri                [n_part] = prev_dentity2_itrf_gnum_and_itrf_strid;
  gnum_and_itrf_weight              [n_part] = malloc(prev_dentity2_itrf_n_blk * sizeof(double));
  pn_entity2_only_by_interface_twice[n_part] = prev_dentity2_itrf_n_blk;
  for(int i = 0; i < prev_dentity2_itrf_n_blk; ++i) {
    gnum_and_itrf_weight[n_part][i] = 1.;
  }

  PDM_part_to_block_t* ptb_itrf = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           gnum_itrf_link,
                                                           gnum_and_itrf_weight,
                                                           pn_entity2_only_by_interface_twice,
                                                           n_part+1,
                                                           comm);

  free(gnum_and_itrf_weight[n_part]);

  int         *dentity2_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *dentity2_gnum_and_itrf_data  = NULL;
  PDM_part_to_block_exch(ptb_itrf,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         gnum_and_itrf_stri,
          (void **)      gnum_and_itrf_data,
                         &dentity2_gnum_and_itrf_strid,
          (void **)      &dentity2_gnum_and_itrf_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(gnum_itrf_link      [i_part]);
    free(gnum_and_itrf_data  [i_part]);
    free(gnum_and_itrf_stri  [i_part]);
    free(gnum_and_itrf_weight[i_part]);
  }
  free(gnum_itrf_link      );
  free(gnum_and_itrf_data  );
  free(gnum_and_itrf_stri  );
  free(gnum_and_itrf_weight);
  free(pn_entity2_only_by_interface_twice);

  int n_blk_elt_itrf = PDM_part_to_block_n_elt_block_get(ptb_itrf);
  PDM_g_num_t *next_dentity2_elt_gnum_ptp_itrf = PDM_part_to_block_block_gnum_get   (ptb_itrf);

  PDM_g_num_t *dentity2_itrf_blk_gnum = (PDM_g_num_t *) malloc(n_blk_elt_itrf * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_blk_elt_itrf; ++i) {
    dentity2_itrf_blk_gnum[i] = next_dentity2_elt_gnum_ptp_itrf[i];
  }

  if(1 == 0) {
    log_trace("\n");
    log_trace("\n");
    int n_data = 0;
    for(int i = 0; i < n_blk_elt_itrf; ++i) {
      n_data += dentity2_gnum_and_itrf_strid[i];
    }
    PDM_log_trace_array_long(next_dentity2_elt_gnum_ptp_itrf,  n_blk_elt_itrf, "next_dentity2_elt_gnum       ::");
    PDM_log_trace_array_int (dentity2_gnum_and_itrf_strid,     n_blk_elt_itrf, "dentity2_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(dentity2_gnum_and_itrf_data , 2 * n_data        , "dentity2_gnum_and_itrf_data  ::");
  }

  int max_strid = 0;
  for(int i = 0; i < n_blk_elt_itrf; ++i) {
    max_strid = PDM_MAX(max_strid, dentity2_gnum_and_itrf_strid[i]);
  }
  int *data_order = malloc(max_strid * sizeof(int));

  /* Unique incomming data */
  int idx_read      = 0;
  int idx_blk_write = 0;
  for(int i = 0; i < n_blk_elt_itrf; ++i) {

    int n_data = dentity2_gnum_and_itrf_strid[i];

    int n_unique = PDM_order_inplace_unique_long(n_data, 2, &dentity2_gnum_and_itrf_data[idx_read], data_order);

    /* Copy */
    for(int i_unique = 0; i_unique < n_unique; ++i_unique) {
      dentity2_gnum_and_itrf_data[2*idx_blk_write  ] = dentity2_gnum_and_itrf_data[(idx_read+2*i_unique)  ];
      dentity2_gnum_and_itrf_data[2*idx_blk_write+1] = dentity2_gnum_and_itrf_data[(idx_read+2*i_unique)+1];
      idx_blk_write++;
    }

    dentity2_gnum_and_itrf_strid[i] = n_unique;
    idx_read += 2 * n_data;
  }
  free(data_order);

  dentity2_gnum_and_itrf_data = realloc(dentity2_gnum_and_itrf_data, 2 * idx_blk_write * sizeof(PDM_g_num_t));

  if(1 == 1) {
    log_trace("\n");
    PDM_log_trace_array_long(next_dentity2_elt_gnum_ptp_itrf,  n_blk_elt_itrf, "next_dentity2_elt_gnum       ::");
    PDM_log_trace_array_int (dentity2_gnum_and_itrf_strid,     n_blk_elt_itrf, "dentity2_gnum_and_itrf_strid (Unique) ::");
    PDM_log_trace_array_long(dentity2_gnum_and_itrf_data , 2 * idx_blk_write , "dentity2_gnum_and_itrf_data  (Unique) ::");
  }

  PDM_part_to_block_free(ptb_itrf);

  /*
   * Keep pointer for output
   */
  *next_dentity2_itrf_n_blk_out               = n_blk_elt_itrf;
  *next_dentity2_itrf_blk_gnum_out            = dentity2_itrf_blk_gnum;
  *next_dentity2_itrf_gnum_and_itrf_strid_out = dentity2_gnum_and_itrf_strid;
  *next_dentity2_itrf_gnum_and_itrf_data_out  = dentity2_gnum_and_itrf_data;

  PDM_gnum_free(gen_gnum_entity2);







  /*
   * Remove entity2 duplicates
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------------------\n");
    log_trace("unique on new entities2\n");
  }

  int          *pn_entity2_extented                     = malloc(n_part * sizeof(int          ));
  int         **pentity2_extented_to_pentity2_idx       = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_to_pentity2_triplet   = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **pentity2_extented_ln_to_gn              = malloc(n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pentity2_extented_ancstr                = malloc(n_part * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_alrdy_sent            = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_n_itrf                = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_path_itrf_idx         = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_path_itrf             = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_to_pentity2_interface = malloc(n_part * sizeof(int         *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int         *_pextract_entity2_idx        = pextract_entity2_idx  [i_part];
    PDM_g_num_t *_pextract_entity2_gnum       = pextract_entity2_gnum [i_part];
    PDM_g_num_t *_pnew_entity2_gnum           = pnew_entity2_gnum     [i_part];
    PDM_g_num_t *_pnew_entity2_ancstr         = pnew_entity2_ancstr   [i_part];
    int         *_pnew_entity2_n_itrf         = pnew_entity2_n_itrf   [i_part];
    int         *_pnew_entity2_kind           = pnew_entity2_kind     [i_part];
    int         *_pnew_entity2_path_itrf_idx  = pnew_entity2_path_itrf_idx [i_part];
    int         *_pnew_entity2_path_itrf      = pnew_entity2_path_itrf     [i_part];

    int n_part1_to_part2          = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];
    int n_new_entity2             =  pn_new_entity2       [i_part];
    int n_new_entity2_n_itrf      =  pn_new_entity2_n_itrf[i_part];

    // All gnum has be unified / shift w/r of interface, only sort along gnum
    int *order = malloc(n_new_entity2 * sizeof(int));
    for(int i_entity2 = 0; i_entity2 < n_new_entity2; ++i_entity2) {
      order[i_entity2] = i_entity2;
    }

    int n_unique = PDM_inplace_unique_long_and_order(_pnew_entity2_gnum,
                                                     order,
                                                     0,
                                                     n_new_entity2-1);

    pentity2_extented_ln_to_gn  [i_part] = malloc(n_unique * sizeof(PDM_g_num_t));
    pentity2_extented_ancstr    [i_part] = malloc(n_unique * sizeof(PDM_g_num_t));
    pentity2_extented_alrdy_sent[i_part] = malloc(n_unique * sizeof(int        ));
    int n_unique2 = 0;
    int size_itrf_array = 0; 
    for(int i = 0; i < n_unique; ++i) {
      // log_trace("i = %d \n", i);

      int pos = PDM_binary_search_long(_pnew_entity2_gnum[i], pentity2_gnum_sorted[i_part], pn_entity2[i_part]);
      if(pos == -1) { // Stockage que si n'existe pas dans la partition
        // log_trace("gnum_old = %d ; gnum_new = %d ; ancstr = %d ; n_itrf= %d ; \n",  _pextract_entity2_gnum[i],  _pnew_entity2_gnum[i],  _pnew_entity2_ancstr[order[i]], _pnew_entity2_n_itrf  [order[i]]);

        pentity2_extented_ln_to_gn  [i_part][n_unique2] = _pnew_entity2_gnum  [i];
        pentity2_extented_alrdy_sent[i_part][n_unique2] = 0;
        pentity2_extented_ancstr    [i_part][n_unique2] = _pnew_entity2_ancstr[order[i]];
        if (_pnew_entity2_kind[order[i]] == 2) {
          pentity2_extented_alrdy_sent[i_part][n_unique2] = 1;
        }
        size_itrf_array += _pnew_entity2_n_itrf  [order[i]];
        n_unique2++;
      }
    }
    pn_entity2_extented[i_part] = n_unique2;

    // Extract all data
    // log_trace("size_itrf_array = %d \n", size_itrf_array);
    pentity2_extented_to_pentity2_idx      [i_part] = malloc( (     n_unique2 + 1) * sizeof(int        ));
    pentity2_extented_to_pentity2_triplet  [i_part] = malloc( ( 3 * n_unique2    ) * sizeof(int        ));
    // pentity2_extented_n_itrf               [i_part] = malloc( (     n_unique2    ) * sizeof(PDM_g_num_t));
    pentity2_extented_path_itrf_idx        [i_part] = malloc( (     n_unique2 + 1) * sizeof(int));
    pentity2_extented_path_itrf            [i_part] = malloc( ( size_itrf_array  ) * sizeof(int));
    pentity2_extented_to_pentity2_interface[i_part] = malloc( (     n_unique2    ) * sizeof(int        ));
    /* Count */
    pentity2_extented_to_pentity2_idx    [i_part][0] = 0;
    int idx_write = 0;
    pentity2_extented_path_itrf_idx[i_part][idx_write] = 0;
    for(int i_entity2 = 0; i_entity2 < n_unique; ++i_entity2) {
      int pos = PDM_binary_search_long(_pnew_entity2_gnum[i_entity2], pentity2_gnum_sorted[i_part], pn_entity2[i_part]);

      if(pos == -1) {
        int old_pos = order[i_entity2];
        pentity2_extented_to_pentity2_idx[i_part][idx_write+1] = pentity2_extented_to_pentity2_idx[i_part][idx_write] + 3;

        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write  ] = pnew_entity2_parent_t[i_part][3*old_pos  ];
        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write+1] = pnew_entity2_parent_t[i_part][3*old_pos+1];
        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write+2] = pnew_entity2_parent_t[i_part][3*old_pos+2];

        // Save the link (gnum, interface)
        if ((_pnew_entity2_kind[old_pos] == 3)||(_pnew_entity2_kind[old_pos] == 4)) {
          int i_last_interf = pnew_entity2_path_itrf_idx[i_part][old_pos+1]-1;
          pentity2_extented_to_pentity2_interface[i_part][idx_write] =  pnew_entity2_path_itrf    [i_part][i_last_interf];
        }
        else if (_pnew_entity2_kind[old_pos] == 2) {
          pentity2_extented_to_pentity2_interface[i_part][idx_write] = 0;
        }
        else {
          PDM_error(__FILE__, __LINE__, 0, "Kind should be 2 or 4 (here = %d).\n", _pnew_entity2_kind[old_pos]);
        }

        pentity2_extented_path_itrf_idx[i_part][idx_write+1] = pentity2_extented_path_itrf_idx[i_part][idx_write] + pnew_entity2_n_itrf  [i_part][  old_pos  ];
        int l_i_ancstr = 0;
        int i_beg = pentity2_extented_path_itrf_idx[i_part][idx_write];
        // log_trace("    i_beg = %d \n", i_beg);

        for (int i_ancstr=_pnew_entity2_path_itrf_idx[old_pos]; i_ancstr<_pnew_entity2_path_itrf_idx[old_pos+1]; ++i_ancstr) {
          // log_trace("    i_ancstr = %d --> _pnew_entity2_path_itrf[i_ancstr] = %d\n", i_ancstr, _pnew_entity2_path_itrf[i_ancstr]);
          pentity2_extented_path_itrf[i_part][i_beg+l_i_ancstr] = _pnew_entity2_path_itrf[i_ancstr];
          l_i_ancstr ++;
        }

        // entities which came from interf=0 has already be sent (cause if dont they wouldn't be here)
        idx_write++;
      }
    }
    assert(idx_write == n_unique2);
    n_unique = n_unique2;

    int n_triplet = pentity2_extented_to_pentity2_idx[i_part][n_unique];


    // // > Build interface path of unique new entity2
    // pentity2_extented_path_itrf_idx[i_part] = PDM_array_new_idx_from_sizes_int(pentity2_extented_n_itrf[i_part], n_unique);
    
    if(1 == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(pentity2_extented_ln_to_gn             [i_part], n_unique        , "pentity2_extented_ln_to_gn      (UNIQUE)  ::");
      PDM_log_trace_array_long(pentity2_extented_ancstr               [i_part], n_unique        , "pentity2_extented_ancstr        (UNIQUE)  ::");
      PDM_log_trace_array_int (pentity2_extented_alrdy_sent           [i_part], n_unique        , "pentity2_extented_alrdy_sent    (UNIQUE)  ::");
      // PDM_log_trace_array_int (pentity2_extented_n_itrf               [i_part], n_unique        , "pentity2_extented_n_itrf        (UNIQUE)  ::");
      PDM_log_trace_array_int (pentity2_extented_path_itrf_idx        [i_part], n_unique+1      , "pentity2_extented_path_itrf_idx (UNIQUE)  ::");
      PDM_log_trace_array_int (pentity2_extented_path_itrf            [i_part], size_itrf_array , "pentity2_extented_path_itrf     (UNIQUE)  ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_idx      [i_part], n_unique+1      , "pentity2_extented_to_pentity2_idx      ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_interface[i_part], n_unique        , "pentity2_extented_to_pentity2_interface::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_triplet  [i_part], n_triplet       , "pentity2_extented_to_pentity2_triplet  ::");
    }

    // /*
    //  * Tag faces that have been sent
    //  */
    // for (int i_face=0; i_face<n_unique; ++i_face) {
    //   int i_entity2 = pentity2_extented_to_pentity2_triplet[i_part][3*i_face+2];
    //   pentity2_alrdy_sent[i_part][i_entity2] = 1;
    // }

    free(order);

  }



  /*
   * To keep after all
   */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    // free(pentity2_interface         [i_part]);
    free(pentity2_only_by_interface_nuplet[i_part]);
    free(pentity2_only_by_interface_parent[i_part]);
  }
  free(pn_new_entity2       );
  free(pn_new_entity2_kind4 );
  free(pn_new_entity2_n_itrf);
  // free(pentity2_interface          );
  free(pentity2_only_by_interface_nuplet);
  free(pentity2_only_by_interface_parent);


  // for(int i_part = 0; i_part < n_part; ++i_part) {
  //   free(pentity2_extented_n_itrf[i_part]);
  // }
  // free(pentity2_extented_n_itrf);



  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pentity2_extented_ancstr_out                = pentity2_extented_ancstr;
  *pentity2_extented_alrdy_sent_out            = pentity2_extented_alrdy_sent;
  *pentity2_extented_path_itrf_idx_out         = pentity2_extented_path_itrf_idx;
  *pentity2_extented_path_itrf_out             = pentity2_extented_path_itrf; //TODO: use this array instead of pentity2_extented_path_itrf ????
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_pentity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_to_pentity2_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_to_pentity2_interface;

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pextract_entity2_idx[i_part]);
  }
  free(pextract_entity2_idx);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pextract_entity2_n         [i_part]);
    free(pextract_entity2_gnum      [i_part]);
    free(pnew_entity2_path_itrf  [i_part]);
    free(pextract_entity2_triplet   [i_part]);
  }
  free(pextract_entity2_n         );
  free(pextract_entity2_gnum      );
  free(pnew_entity2_path_itrf);
  free(pextract_entity2_triplet   );

}

/*
 * Translate and post-treated link by interface with on entity to another
 *  Exemple of use :
 *     - entity1 to face
 *     - entity1 to edge
 *     - face to cell
 */
void
PDM_part_extension_interface_by_entity1_to_interface_by_entity2
(
  PDM_part_domain_interface_t  *pdi,
  PDM_bound_type_t              entity1_bound,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_entity2,
  int                          *n_part,
  int                         **pn_entity1_in,
  PDM_g_num_t                ***pentity1_ln_to_gn_in,
  int                        ***pentity1_hint_in,
  int                         **pn_entity2_in,
  PDM_g_num_t                ***pentity2_ln_to_gn_in,
  int                        ***pentity2_entity1_idx_in,
  int                        ***pentity2_entity1_in,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                  comm
)
{
  int n_part_tot = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    n_part_tot += n_part[i_domain];
  }

  int          *pn_entity1            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_hint         = NULL;
  int          *pn_entity2            = malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ln_to_gn     = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_entity1_idx  = malloc(n_part_tot * sizeof(int         *));
  int         **pentity2_entity1      = malloc(n_part_tot * sizeof(int         *));

  if(pentity1_hint_in != NULL) {
    pentity1_hint = malloc(n_part_tot * sizeof(int         *));
  }

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot] = pentity1_ln_to_gn_in   [i_dom][i_part];
      pn_entity2          [ln_part_tot] = pn_entity2_in          [i_dom][i_part];
      pentity2_ln_to_gn   [ln_part_tot] = pentity2_ln_to_gn_in   [i_dom][i_part];
      pentity2_entity1_idx[ln_part_tot] = pentity2_entity1_idx_in[i_dom][i_part];
      pentity2_entity1    [ln_part_tot] = pentity2_entity1_in    [i_dom][i_part];

      if(pentity1_hint_in != NULL && pentity1_hint_in[i_dom] != NULL) {
        pentity1_hint[ln_part_tot] = pentity1_hint_in[i_dom][i_part];
      }

      ln_part_tot += 1;
    }
  }


  int         **pentity1_entity2_idx  = NULL;
  int         **pentity1_entity2      = NULL;
  PDM_part_connectivity_transpose(ln_part_tot,
                                  pn_entity2,
                                  pn_entity1,
                                  pentity2_entity1_idx,
                                  pentity2_entity1,
                                  &pentity1_entity2_idx,
                                  &pentity1_entity2);


  int          *pn_entity1_num             = NULL;
  int         **pentity1_num               = NULL;
  int         **pentity1_opp_location_idx  = NULL;
  int         **pentity1_opp_location      = NULL;
  int         **pentity1_opp_interface_idx = NULL;
  int         **pentity1_opp_interface     = NULL;
  int         **pentity1_opp_sens          = NULL;
  PDM_g_num_t **pentity1_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity1_bound,
                                         pn_entity1,
                                         pentity1_ln_to_gn,
                                         &pn_entity1_num,
                                         &pentity1_num,
                                         &pentity1_opp_location_idx,
                                         &pentity1_opp_location,
                                         &pentity1_opp_interface_idx,
                                         &pentity1_opp_interface,
                                         &pentity1_opp_sens,
                                         &pentity1_opp_gnum);

  /*
   * Prepare creation of part_to_part
   */
  int         **part1_to_part2_idx              = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_triplet_idx      = NULL; //malloc(n_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet          = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_interface        = malloc(n_part_tot * sizeof(int         *));

  int         **part1_to_part2_entity2_n           = malloc(n_part_tot * sizeof(int         *));
  int         **part1_to_part2_entity2_triplet     = malloc(n_part_tot * sizeof(int         *));

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(comm, n_part_tot);

  int **ppart_entity1_proc_idx = NULL;
  int **ppart_entity1_part_idx = NULL;
  int **ppart_entity1          = NULL; // (i_entity1, i_proc, i_part, i_entity1_opp)
  PDM_part_generate_entity_graph_comm(comm,
                                      part_distribution,
                                      NULL,
                                      n_part_tot,
                                      pn_entity1,
               (const PDM_g_num_t **) pentity1_ln_to_gn,
               (const int         **) pentity1_hint,
                                      &ppart_entity1_proc_idx,
                                      &ppart_entity1_part_idx,
                                      &ppart_entity1,
                                      NULL);

  free(part_distribution);

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);
  int n_g_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_g_part_tot += n_part_g[i_dom];
  }

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = (int *) malloc( (n_rank) * sizeof(int));
    PDM_MPI_Allgather(&li_part,
                      1,
                      PDM_MPI_INT,
                      n_part_shift,
                      1,
                      PDM_MPI_INT,
                      comm);

    if(0 == 1) {
      PDM_log_trace_array_int(n_part_shift, n_rank+1, "n_part_shift ::");
    }

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      part1_to_part2_idx[li_part] = malloc((pn_entity1[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pn_entity1[li_part]);

      int n_entity_bound = ppart_entity1_part_idx[i_part][n_g_part_tot];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_entity1[i_part][4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        int n_opp = pentity1_opp_location_idx[li_part][idx_entity+1] - pentity1_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }


      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pn_entity1[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      // printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_entity1[i_part][4*idx_entity]-1;
        int i_proc_opp   = ppart_entity1[i_part][4*idx_entity+1];
        int i_part_opp   = ppart_entity1[i_part][4*idx_entity+2]-1;
        int i_entity_opp = ppart_entity1[i_part][4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp + n_part_shift[i_proc_opp];
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
        int i_entity = pentity1_num[li_part][idx_entity];
        for(int idx_opp = pentity1_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pentity1_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pentity1_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pentity1_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pentity1_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pentity1_opp_interface[li_part][  idx_opp  ];
        }
      }


      /*
       * Creation des buffers d'envoi des connectivités
       */
      int n_send = part1_to_part2_idx[li_part][pn_entity1[li_part]]/3;
      part1_to_part2_entity2_n    [li_part] = malloc(n_send * sizeof(int));
      int n_send_entity2     = 0;

      int *_pentity1_entity2_idx     = pentity1_entity2_idx    [li_part];
      int *_part1_to_part2_entity2_n = part1_to_part2_entity2_n[li_part];

      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {
          _part1_to_part2_entity2_n    [idx] = 0;
          for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2 ) {
            // int i_entity2 = PDM_ABS(pentity1_entity2[li_part][idx_entity2])-1;
            _part1_to_part2_entity2_n    [idx] += 1;
          }
          n_send_entity2 += _part1_to_part2_entity2_n    [idx];
        }
      }


      part1_to_part2_entity2_triplet    [li_part] = malloc(3 * n_send_entity2     * sizeof(int        ));
      int         *_part1_to_part2_entity2_triplet   = part1_to_part2_entity2_triplet  [li_part];

      /* Fill loop */
      n_send_entity2     = 0;
      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          /* Copy */
          for(int idx_entity2 = _pentity1_entity2_idx[i_entity]; idx_entity2 < _pentity1_entity2_idx[i_entity+1]; ++idx_entity2 ) {
            int i_entity2 = PDM_ABS(pentity1_entity2[li_part][idx_entity2])-1;

            _part1_to_part2_entity2_triplet[3*n_send_entity2  ] = i_rank;
            _part1_to_part2_entity2_triplet[3*n_send_entity2+1] = li_part;
            _part1_to_part2_entity2_triplet[3*n_send_entity2+2] = i_entity2;
            n_send_entity2++;

          }
        }
      }

      // PDM_log_trace_array_int(part1_to_part2_triplet[li_part] , n_connect_tot, "part1_to_part2_triplet ::");

      free(part1_to_part2_n);

      li_part += 1;
    }

    free(n_part_shift);
  }

  free(n_part_g);

  // TODO : Faire un part_to_block pour hook TOUTES les faces qui pointe sur le même sommet par exemple


  /*
   * Create part_to_part to exchange all data in opposit part
   */
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_ln_to_gn,
                                                                      (const int          *) pn_entity1,
                                                                      n_part_tot,
                                                                      (const int          *) pn_entity1,
                                                                      n_part_tot,
                                                                      (const int         **) part1_to_part2_idx,
                                                                      (const int         **) part1_to_part2_triplet_idx,
                                                                      (const int         **) part1_to_part2_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_triplet_n    = malloc(n_part_tot * sizeof(int         *));
  int         **gnum1_com_from_triplet_send = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_gnum_send    = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx     = pentity1_entity2_idx    [i_part];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
      }
    }

    /* Allocate */
    gnum1_com_from_triplet_n   [i_part] = malloc(    n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_triplet_send[i_part] = malloc(3 * n_send_part2      * sizeof(int        ));
    gnum1_com_from_gnum_send   [i_part] = malloc(    n_send_part2      * sizeof(PDM_g_num_t));

    int         *_gnum1_com_from_triplet_n    = gnum1_com_from_triplet_n   [i_part];
    int         *_gnum1_com_from_triplet_send = gnum1_com_from_triplet_send[i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_send    = gnum1_com_from_gnum_send   [i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {

          int i_entity2 = PDM_ABS(pentity1_entity2[i_part][idx_entity2])-1;
          _gnum1_com_from_gnum_send   [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
          _gnum1_com_from_triplet_send[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_triplet_send[3*n_send_part2+1] = i_part;
          _gnum1_com_from_triplet_send[3*n_send_part2+2] = i_entity2;
          n_send_part2++;
        }
      }
    }

    PDM_log_trace_array_long(_gnum1_com_from_gnum_send, n_send_part2, "_gnum1_com_from_gnum_send ::");


  }

  int exch_request = -1;
  int         **pextract_entity2_n    = NULL;
  PDM_g_num_t **pextract_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_gnum_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int  **)  gnum1_com_from_triplet_n,
                (const void **)  gnum1_com_from_triplet_send,
                                 &pextract_entity2_n,
                    (void ***)   &pextract_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(gnum1_com_from_triplet_n   [i_part]);
    free(gnum1_com_from_triplet_send[i_part]);
    free(gnum1_com_from_gnum_send   [i_part]);
  }


  if(1 == 1) { // Usefull to know how many data is transfer
    for(int i_part = 0; i_part < n_part_tot; ++i_part) {

      int n_triplet = part1_to_part2_idx[i_part][pn_entity1[i_part]];
      PDM_log_trace_array_int(part1_to_part2_idx       [i_part], pn_entity1[i_part]+1, "part1_to_part2_idx ::");
      PDM_log_trace_array_int(part1_to_part2_triplet   [i_part], n_triplet           , "part1_to_part2_triplet ::");
      PDM_log_trace_array_int(part1_to_part2_interface [i_part], n_triplet/3         , "part1_to_part2_interface ::");
      PDM_log_trace_graph_nuplet_int(part1_to_part2_idx[i_part],
                                     part1_to_part2_triplet[i_part], 1,
                                     pn_entity1[i_part], "part1_to_part2_triplet ::");


      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");

      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "extract_lnum2 :");
      PDM_log_trace_connectivity_long(gnum1_come_from_idx[i_part],
                                      gnum1_come_from    [i_part],
                                      n_ref_lnum2  [i_part], "gnum1_come_from ::");
    }
  }


  int **pextract_entity2_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    pextract_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity2_n[i_part], n_part1_to_part2);
  }


  _recurse_and_filter(ptp,
                      n_part_tot,
                      pn_entity1,
                      part1_to_part2_idx,
                      part1_to_part2_interface,
                      pentity1_entity2_idx,
                      pentity1_entity2,
                      pentity2_ln_to_gn,
                      pextract_entity2_n,
                      pextract_entity2_gnum,
                      pextract_entity2_triplet,
                      comm);

  free(gnum1_com_from_triplet_n   );
  free(gnum1_com_from_triplet_send);
  free(gnum1_com_from_gnum_send   );


  /*
   * Post-treatment
   *   - Remove duplicate (same gnum or same equivalent transformation )
   *   - Keep link between entity2 (extented and current part)
   *   - Create new global numbering for all new entities
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     n_part_tot,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2);

  int          *pn_entity2_only_by_interface        = malloc(ln_part_tot * sizeof(int          ));
  int         **pentity2_interface                  = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_ln_to_gn_only_by_interface = malloc(ln_part_tot * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_entity2_gnum    = pextract_entity2_gnum   [i_part];
    int         *_pextract_entity2_idx     = pextract_entity2_idx    [i_part];

    int n_part1_to_part2          = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];
    pn_entity2_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        pn_entity2_only_by_interface[i_part] += pextract_entity2_n[i_part][i];
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_entity2_gnum, n_part1_to_part2_recv_tot, "_pextract_entity2_gnum ::");
    }

    pentity2_ln_to_gn_only_by_interface[i_part] = malloc(2 * pn_entity2_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    pentity2_interface                 [i_part] = malloc(    n_part1_to_part2_recv_tot            * sizeof(int        ));
    PDM_g_num_t *_pentity2_ln_to_gn_only_by_interface = pentity2_ln_to_gn_only_by_interface[i_part];

    pn_entity2_only_by_interface[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
          int idx_write = pn_entity2_only_by_interface[i_part]++;
          _pentity2_ln_to_gn_only_by_interface[2*idx_write  ] = _pextract_entity2_gnum[j];
          _pentity2_ln_to_gn_only_by_interface[2*idx_write+1] = part1_to_part2_interface[i_part][i];
        }
      }

      for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
        pentity2_interface[i_part][j] = part1_to_part2_interface[i_part][i];
      }
    }

    if(1 == 0) {
      PDM_log_trace_array_long(pentity2_ln_to_gn_only_by_interface[i_part], 2 * pn_entity2_only_by_interface[i_part], "pentity2_ln_to_gn_only_by_interface ::");
    }


    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              pn_entity2_only_by_interface[i_part],
                              pentity2_ln_to_gn_only_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_entity2);

  /* Update */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    PDM_g_num_t* extented_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    PDM_g_num_t *_pextract_entity2_gnum    = pextract_entity2_gnum   [i_part];
    int         *_pextract_entity2_idx     = pextract_entity2_idx    [i_part];

    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;

    int idx_read = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(part1_to_part2_interface[i_part][i] != 0) {
        for(int j = _pextract_entity2_idx[i]; j < _pextract_entity2_idx[i+1]; ++j) {
          _pextract_entity2_gnum[j] = extented_entity2_ln_to_gn[idx_read++] + shift_by_domain_entity2[n_domain];
        }
      }
    }

    if(1 == 0) {
      PDM_log_trace_array_long(_pextract_entity2_gnum, _pextract_entity2_idx[n_part1_to_part2], "_pextract_entity2_gnum (Update) : ");
    }

  }


  PDM_gnum_free(gen_gnum_entity2);

  /*
   * Local post-treatment
   *   - Remove duplicate
   *   - Remove same entity but wih different path
   */
  int          *pn_entity2_extented                     = malloc(ln_part_tot * sizeof(int          ));
  int         **pentity2_extented_to_pentity2_idx       = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_to_pentity2_triplet   = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_extented_ln_to_gn              = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **extented_entity2_orig_gnum              = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_to_pentity2_interface = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_entity2_gnum = pextract_entity2_gnum[i_part];
    int         *_pextract_entity2_idx  = pextract_entity2_idx [i_part];

    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];

    // All gnum has be unified / shift w/r of interface, only sort along gnum
    int *order = malloc(n_part1_to_part2_recv_tot * sizeof(int));
    for(int i_entity2 = 0; i_entity2 < n_part1_to_part2_recv_tot; ++i_entity2) {
      order[i_entity2] = i_entity2;
    }

    // Maybe faire inplace sort and revert with order after ?
    PDM_g_num_t* sorted_pentity2_ln_to_gn = malloc(pn_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i = 0; i < pn_entity2[i_part]; ++i) {
      sorted_pentity2_ln_to_gn[i] = pentity2_ln_to_gn[i_part][i];
    }
    PDM_sort_long(sorted_pentity2_ln_to_gn, NULL, pn_entity2[i_part]);

    int n_unique = PDM_inplace_unique_long_and_order(_pextract_entity2_gnum,
                                                     order,
                                                     0,
                                                     n_part1_to_part2_recv_tot-1);

    pentity2_extented_ln_to_gn[i_part] = malloc(n_unique * sizeof(PDM_g_num_t));
    int n_unique2 = 0;
    for(int i = 0; i < n_unique; ++i) {
      int pos = PDM_binary_search_long(_pextract_entity2_gnum[i], sorted_pentity2_ln_to_gn, pn_entity2[i_part]);
      if(pos == -1) {
        pentity2_extented_ln_to_gn[i_part][n_unique2] = _pextract_entity2_gnum[i];
        n_unique2++;
      }
    }
    pn_entity2_extented[i_part] = n_unique2;

    // Extract all data
    pentity2_extented_to_pentity2_idx      [i_part] = malloc( (     n_unique2 + 1) * sizeof(int        ));
    pentity2_extented_to_pentity2_triplet  [i_part] = malloc( ( 3 * n_unique2    ) * sizeof(int        ));
    extented_entity2_orig_gnum             [i_part] = malloc( (     n_unique2    ) * sizeof(PDM_g_num_t));
    pentity2_extented_to_pentity2_interface[i_part] = malloc( (     n_unique2    ) * sizeof(int        ));

    /* Count */
    pentity2_extented_to_pentity2_idx    [i_part][0] = 0;
    int idx_write = 0;
    for(int i_entity2 = 0; i_entity2 < n_unique; ++i_entity2) {
      int pos = PDM_binary_search_long(_pextract_entity2_gnum[i_entity2], sorted_pentity2_ln_to_gn, pn_entity2[i_part]);

      if(pos == -1) {
        int old_pos = order[i_entity2];
        pentity2_extented_to_pentity2_idx[i_part][idx_write+1] = pentity2_extented_to_pentity2_idx[i_part][idx_write] + 3;

        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write  ] = pextract_entity2_triplet[i_part][3*old_pos  ];
        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write+1] = pextract_entity2_triplet[i_part][3*old_pos+1];
        pentity2_extented_to_pentity2_triplet[i_part][3*idx_write+2] = pextract_entity2_triplet[i_part][3*old_pos+2];

        // Save the link (gnum, interface)
        pentity2_extented_to_pentity2_interface[i_part][idx_write] = pentity2_interface    [i_part][  old_pos  ];
        extented_entity2_orig_gnum             [i_part][idx_write] = _pextract_entity2_gnum        [  i_entity2]; // Because it's sorted already
        idx_write++;
      }
    }
    assert(idx_write == n_unique2);
    n_unique = n_unique2;

    int n_triplet = pentity2_extented_to_pentity2_idx[i_part][n_unique];

    if(0 == 1) {
      PDM_log_trace_array_long(_pextract_entity2_gnum                         , n_unique   , "_pextract_entity2_gnum (UNIQUE)        ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_idx      [i_part], n_unique+1 , "pentity2_extented_to_pentity2_idx      ::");
      PDM_log_trace_array_long(extented_entity2_orig_gnum             [i_part], n_triplet/3, "extented_entity2_orig_gnum             ::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_interface[i_part], n_triplet/3, "pentity2_extented_to_pentity2_interface::");
      PDM_log_trace_array_int (pentity2_extented_to_pentity2_triplet  [i_part], n_triplet  , "pentity2_extented_to_pentity2_triplet  ::");
    }

    free(order);
    free(sorted_pentity2_ln_to_gn);

  }

  /*
   * A partir de toute ces infos on peut rajouter tout dans le PDM_part_domain_interface_set
   *   -> A généraliser pour des cellules
   *   -> Il manque les échanges pour les domaine
   *   -> Faire également un PDM_part_domain_interface_set_replace pour updater avec les nouvelles infos
   *   -> C'est cette info qu'on utilisera pour le rang 2 pour unifier les faces deja ramener
   */


  /*
   * Puis il faut également sortir toutes les infos pour construire le ptp pour faire la cascade
   *  part1 = part_extented
   *  part2 = part_current
   */
  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_pentity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_to_pentity2_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_to_pentity2_interface;



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(extented_entity2_orig_gnum[i_part]);
  }
  free(extented_entity2_orig_gnum);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_ln_to_gn_only_by_interface[i_part]);
  }

  free(pn_entity2_only_by_interface);
  free(pentity2_ln_to_gn_only_by_interface);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity2_idx[i_part]);
    free(pentity2_interface  [i_part]);
  }
  free(pextract_entity2_idx);
  free(pentity2_interface);



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity2_n        [i_part]);
    free(pextract_entity2_gnum     [i_part]);
    free(pextract_entity2_triplet  [i_part]);
    free(pentity1_entity2_idx      [i_part]);
    free(pentity1_entity2          [i_part]);
  }
  free(pextract_entity2_n        );
  free(pextract_entity2_gnum     );
  free(pextract_entity2_triplet  );
  free(pentity1_entity2_idx      );
  free(pentity1_entity2          );


  PDM_part_to_part_free(ptp);


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(ppart_entity1_proc_idx[i_part]);
    free(ppart_entity1_part_idx[i_part]);
    free(ppart_entity1         [i_part]);
  }
  free(ppart_entity1_proc_idx);
  free(ppart_entity1_part_idx);
  free(ppart_entity1         );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(part1_to_part2_idx            [i_part] );
    free(part1_to_part2_triplet        [i_part]);
    free(part1_to_part2_interface      [i_part]);
    free(part1_to_part2_entity2_n      [i_part]);
    free(part1_to_part2_entity2_triplet[i_part]);

  }

  free(part1_to_part2_idx            );
  free(part1_to_part2_triplet        );
  free(part1_to_part2_interface      );
  free(part1_to_part2_entity2_n      );
  free(part1_to_part2_entity2_triplet);

  free(pn_entity1          );
  free(pentity1_ln_to_gn   );
  free(pn_entity2          );
  free(pentity2_ln_to_gn   );
  free(pentity2_entity1_idx);
  free(pentity2_entity1    );

  if(pentity1_hint_in != NULL) {
    free(pentity1_hint);
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity1_num             [i_part]);
    free(pentity1_opp_location_idx[i_part]);
    free(pentity1_opp_location    [i_part]);
    free(pentity1_opp_interface   [i_part]);
    free(pentity1_opp_sens        [i_part]);
    free(pentity1_opp_gnum        [i_part]);
  }
  free(pn_entity1_num            );
  free(pentity1_num              );
  free(pentity1_opp_location_idx );
  free(pentity1_opp_location     );
  free(pentity1_opp_interface_idx);
  free(pentity1_opp_interface    );
  free(pentity1_opp_sens         );
  free(pentity1_opp_gnum         );

}


void
PDM_part_extension_pentity1_entity2_to_extented_pentity1_entity2
(
 int                    n_part,
 int                    n_interface,
 PDM_g_num_t            shift_by_domain_entity2,
 int                    prev_dentity2_itrf_n_blk,
 PDM_g_num_t           *prev_dentity2_itrf_blk_gnum,
 int                   *prev_dentity2_itrf_gnum_and_itrf_strid,
 PDM_g_num_t           *prev_dentity2_itrf_gnum_and_itrf_data,
 int                   *prev_dentity2_itrf_gnum_and_itrf_sens,
 int                   *pn_entity1,
 PDM_g_num_t          **pentity1_ln_to_gn,
 int                   *pn_entity2,
 PDM_g_num_t          **pentity2_ln_to_gn,
 int                  **pentity1_entity2_idx,
 int                  **pentity1_entity2,
 int                   *pn_entity1_extented,
 PDM_g_num_t          **pentity1_extented_ln_to_gn,
 int                  **pentity1_extented_to_pentity1_idx,
 int                  **pentity1_extented_to_pentity1_triplet,
 int                  **pentity1_extented_to_pentity1_interface,
 int                  **pn_entity2_extented_out,
 PDM_g_num_t         ***pentity2_extented_ln_to_gn_out,
 int                 ***pextented_entity1_entity2_idx_out,
 int                 ***pextented_entity1_entity2_out,
 int                 ***pentity2_extented_to_pentity2_idx_out,
 int                 ***pentity2_extented_to_pentity2_triplet_out,
 int                 ***pentity2_extented_to_pentity2_interface_out,
 int                   *next_dentity2_itrf_n_blk_out,
 PDM_g_num_t          **next_dentity2_itrf_blk_gnum_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_strid_out,
 PDM_g_num_t          **next_dentity2_itrf_gnum_and_itrf_data_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_sens_out,
 PDM_MPI_Comm           comm
)
{
  PDM_UNUSED(n_interface);
  PDM_UNUSED(pentity1_ln_to_gn);

  int debug = 1;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Creation du ptp between extented part and current partition for entity1
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------\n");
    log_trace("PTP to exchange entities2\n");
  }

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pentity1_extented_ln_to_gn[i_part], pn_entity1_extented[i_part], "pentity1_extended_ln_to_gn ::");
      PDM_log_trace_array_int(pentity1_extented_to_pentity1_idx[i_part], pn_entity1_extented[i_part]+1, "pentity1_extended_to_pentity1_idx ::");
      PDM_log_trace_array_int(pentity1_extented_to_pentity1_triplet[i_part], pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]], "pentity1_extended_to_pentity1_triplet ::");
      PDM_log_trace_array_int(pentity1_extented_to_pentity1_interface[i_part], pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3, "pentity1_extended_to_pentity1_interface ::");

    }
  }

  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_extented_ln_to_gn,
                                                                      (const int          *) pn_entity1_extented,
                                                                      n_part,
                                                                      (const int          *) pn_entity1,
                                                                      n_part,
                                                                      (const int         **) pentity1_extented_to_pentity1_idx,
                                                                      (const int         **) NULL,
                                                                      (const int         **) pentity1_extented_to_pentity1_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Remplissage des buffers et envoi
   */
  int         **gnum1_com_from_entity1_entity2_n       = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_entity1_entity2         = malloc(n_part * sizeof(PDM_g_num_t *));
  int         **gnum1_com_from_entity1_entity2_triplet = malloc(n_part * sizeof(int         *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx = pentity1_entity2_idx[i_part];
    int *_pentity1_entity2     = pentity1_entity2    [i_part];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
      }
    }

    /* Allocate */
    gnum1_com_from_entity1_entity2_n      [i_part] = malloc(     n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_entity1_entity2        [i_part] = malloc(     n_send_part2      * sizeof(PDM_g_num_t));
    gnum1_com_from_entity1_entity2_triplet[i_part] = malloc( 3 * n_send_part2      * sizeof(int        ));

    int         *_gnum1_com_from_triplet_n               = gnum1_com_from_entity1_entity2_n      [i_part];
    PDM_g_num_t *_gnum1_com_from_entity1_entity2         = gnum1_com_from_entity1_entity2        [i_part];
    int         *_gnum1_com_from_entity1_entity2_triplet = gnum1_com_from_entity1_entity2_triplet[i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
          _gnum1_com_from_entity1_entity2        [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+1] = i_part;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+2] = i_entity2;
          n_send_part2++;
        }
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(_gnum1_com_from_entity1_entity2        , n_send_part2  , "_gnum1_com_from_gnum_send               ::");
      PDM_log_trace_array_int (_gnum1_com_from_entity1_entity2_triplet, n_send_part2*3, "_gnum1_com_from_entity1_entity2_triplet ::");
    }

  }

  /*
   * Exchange:
   *     - vertices gnum    of face from initial domain referenced
   *     - vertices triplet of face from initial domain referenced
   */
  int exch_request = -1;
  int         **pextract_entity1_entity2_n    = NULL;
  PDM_g_num_t **pextract_entity1_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity1_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2_triplet,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(gnum1_com_from_entity1_entity2_n      [i_part]);
    free(gnum1_com_from_entity1_entity2        [i_part]);
    free(gnum1_com_from_entity1_entity2_triplet[i_part]);
  }
  free(gnum1_com_from_entity1_entity2_n      );
  free(gnum1_com_from_entity1_entity2        );
  free(gnum1_com_from_entity1_entity2_triplet);
  PDM_part_to_part_free(ptp);


  int **pextract_entity1_entity2_idx = malloc(n_part * sizeof(int *));
  int *pextract_entity1_entity2_n_elmt = malloc(n_part * sizeof(int *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    pextract_entity1_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity1_entity2_n[i_part], n_part1_to_part2);

    pextract_entity1_entity2_n_elmt[i_part] = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      pextract_entity1_entity2_n_elmt[i_part] += pextract_entity1_entity2_n[i_part][i];
    }
  }

  /*
   * Prepare unification for all incoming entity2 by interface
   *
   */

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(prev_dentity2_itrf_blk_gnum,
                                                                        prev_dentity2_itrf_n_blk,
                                                (const PDM_g_num_t **)  pextract_entity1_entity2_gnum,
                                                                        pextract_entity1_entity2_n_elmt,
                                                                        n_part,
                                                                        comm);

  int         **prev_pentity2_itrf_gnum_and_itrf_strid = NULL;
  PDM_g_num_t **prev_pentity2_itrf_gnum_and_itrf_data  = NULL;

  PDM_block_to_part_exch(btp,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         prev_dentity2_itrf_gnum_and_itrf_strid,
                         prev_dentity2_itrf_gnum_and_itrf_data,
                         &prev_pentity2_itrf_gnum_and_itrf_strid,
         (void ***)      &prev_pentity2_itrf_gnum_and_itrf_data);
  PDM_block_to_part_free(btp);



  /*
   * Need to sort pentity2_gnum to search if a candidate entity
   * doesn't already exist in partition.
   *
   */
  int         **entity2_order            = malloc( n_part * sizeof(int         *));
  PDM_g_num_t **pentity2_ln_to_gn_sorted = malloc( n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {
    entity2_order           [i_part] = malloc( pn_entity2[i_part] * sizeof(int        ));
    pentity2_ln_to_gn_sorted[i_part] = malloc( pn_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {
      entity2_order           [i_part][i_entity2] = i_entity2;
      pentity2_ln_to_gn_sorted[i_part][i_entity2] = pentity2_ln_to_gn[i_part][i_entity2];
    }
    PDM_sort_long(pentity2_ln_to_gn_sorted[i_part], entity2_order[i_part], pn_entity2[i_part]);
  }



  /*
   * Post-treatment - same as PDM_part_extension_entity1_to_entity2
   */

  /*
   * pextract_entity2_kind :
   *    1/ Internal
   *       a/ Already know in current partition                    : 0
   *       b/ New in current partition                             : 2
   *    2/ From interface
   *       a/ Know by relation table and know in current partition : 1
   *       b/ Know by relation table but not local                 : 3
   *       c/ New                                                  : 4
   * On essaye après d'organiser les nouvelles entités dans l'ordre suivant :
   *    [2/3/4]
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("--------------------\n");
    log_trace("find valid entities2\n");
  }

  int n_extented_kind = 5;
  int         **pextract_entity2_kind           = malloc(n_part * sizeof(int *));
  int         **pextract_entity2_translate      = malloc(n_part * sizeof(int *));
  int         **pextract_entity2_interface      = malloc(n_part * sizeof(int *));
  int         **recv_buffer_to_sort_kind_order  = malloc(n_part * sizeof(int *));
  int         **pentity2_extented_by_kind_idx   = malloc(n_part * sizeof(int *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

   int n_data = 0;
    for(int i = 0; i < pextract_entity1_entity2_n_elmt[i_part]; ++i) {
      n_data += prev_pentity2_itrf_gnum_and_itrf_strid[i_part][i];
    }

    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    if(1 == 1) {
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum          [i_part], pextract_entity1_entity2_n_elmt[i_part], "pextract_entity1_entity2_gnum          ::");
      PDM_log_trace_array_int (prev_pentity2_itrf_gnum_and_itrf_strid [i_part], pextract_entity1_entity2_n_elmt[i_part], "prev_pentity2_itrf_gnum_and_itrf_strid ::");
      PDM_log_trace_array_int (pentity1_extented_to_pentity1_interface[i_part],     n_part1_to_part2                   , "pentity1_extented_to_pentity1_interface::");
      PDM_log_trace_array_long(prev_pentity2_itrf_gnum_and_itrf_data  [i_part], 2 * n_data                             , "prev_pentity2_itrf_gnum_and_itrf_data  ::");
    }

    PDM_g_num_t *_prev_pentity2_itrf_gnum_and_itrf_data = prev_pentity2_itrf_gnum_and_itrf_data[i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum         = pextract_entity1_entity2_gnum        [i_part];
    pextract_entity2_kind     [i_part] = malloc(pextract_entity1_entity2_n_elmt[i_part] * sizeof(int        ));
    pextract_entity2_translate[i_part] = malloc(pextract_entity1_entity2_n_elmt[i_part] * sizeof(int));
    pextract_entity2_interface[i_part] = malloc(pextract_entity1_entity2_n_elmt[i_part] * sizeof(int));

    int         *_pextract_entity2_translate = pextract_entity2_translate[i_part];
    int         *_pextract_entity2_interface = pextract_entity2_interface[i_part];
    PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted  [i_part];


    int idx_read      = 0;
    int idx_read_data = 0;
    for(int i = 0; i < n_part1_to_part2; ++i) {
      if(pentity1_extented_to_pentity1_interface[i_part][i] != 0) {
        for(int j = 0; j < pextract_entity1_entity2_n[i_part][i]; ++j) {

          int cur_itrf     = PDM_ABS (pentity1_extented_to_pentity1_interface[i_part][i]);
          int sgn_cur_itrf = PDM_SIGN(pentity1_extented_to_pentity1_interface[i_part][i]);

          // log_trace(" ----------------- gnum = ("PDM_FMT_G_NUM",%i) \n", _pextract_entity1_entity2_gnum[idx_read], pentity1_extented_to_pentity1_interface[i_part][i]);

          int         keep = 1;
          PDM_g_num_t gnum_opp = 0;
          PDM_g_num_t gnum_cur = 0;
          for(int k = 0; k < prev_pentity2_itrf_gnum_and_itrf_strid[i_part][idx_read]; ++k) {
            // log_trace("\t ("PDM_FMT_G_NUM"/"PDM_FMT_G_NUM") \n", _prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data], _prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data+1]);

            int         opp_itrf     = PDM_ABS (_prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data+1]);
            int         opp_sgn_itrf = PDM_SIGN(_prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data+1]);

            if(cur_itrf == opp_itrf && sgn_cur_itrf == opp_sgn_itrf && gnum_opp == 0) { // Sign a unifier avec l'autre
              keep = 0;

              // Si on connait une correspondance
              //  -> Soit on a la correspondance chez nous (on cherche dans le ln_to_gn )
              //  -> Soit c'est nouveau chez nous, donc on doit recreer un numero local
              // _pextract_entity2_translate[idx_read_data]
              gnum_cur     = _pextract_entity1_entity2_gnum        [idx_read];
              gnum_opp     = _prev_pentity2_itrf_gnum_and_itrf_data[2*idx_read_data];
            }

            idx_read_data++;
          }

          if(keep == 1) {
            // log_trace("Not found = ("PDM_FMT_G_NUM",%i) \n", _pextract_entity1_entity2_gnum[idx_read], pentity1_extented_to_pentity1_interface[i_part][i]);

            _pextract_entity2_translate  [idx_read] = -1;
            _pextract_entity2_interface  [idx_read] = pentity1_extented_to_pentity1_interface[i_part][i];
            pextract_entity2_kind[i_part][idx_read] = 4;

          } else {
            // Search in local
            // log_trace("Local found = ("PDM_FMT_G_NUM",%i) / gnum_opp = "PDM_FMT_G_NUM" \n", _pextract_entity1_entity2_gnum[idx_read], pentity1_extented_to_pentity1_interface[i_part][i], gnum_opp);

            assert(gnum_opp != 0);
            assert(gnum_cur != 0);

            // CAUTION SENS NEED TO BE MANAGE !!!!!!!
            int pos_int = PDM_binary_search_long(gnum_opp, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);
            // log_trace("pos_int = %i \n", pos_int);
            if( pos_int != -1) {
              _pextract_entity2_translate  [idx_read] = entity2_order[i_part][pos_int];
              _pextract_entity2_interface  [idx_read] = pentity1_extented_to_pentity1_interface[i_part][i];
              pextract_entity2_kind[i_part][idx_read] = 1;
            } else {
              _pextract_entity1_entity2_gnum[idx_read] = gnum_opp;
              _pextract_entity2_translate  [idx_read] = -1;
              _pextract_entity2_interface  [idx_read] = pentity1_extented_to_pentity1_interface[i_part][i];
              pextract_entity2_kind[i_part][idx_read] = 3;
            }

          }

          idx_read++;
        }
      } else {
        for(int j = 0; j < pextract_entity1_entity2_n[i_part][i]; ++j) {

          PDM_g_num_t gnum_cur = _pextract_entity1_entity2_gnum        [idx_read];

          int pos_int = PDM_binary_search_long(gnum_cur, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);
          if( pos_int != -1) {
            _pextract_entity2_translate  [idx_read] = entity2_order[i_part][pos_int];
            _pextract_entity2_interface  [idx_read] = pentity1_extented_to_pentity1_interface[i_part][i];
            pextract_entity2_kind[i_part][idx_read] = 0;
          } else {
            _pextract_entity2_translate  [idx_read] = -1;
            _pextract_entity2_interface  [idx_read] = pentity1_extented_to_pentity1_interface[i_part][i];
            pextract_entity2_kind[i_part][idx_read] = 2;
          }

          for(int k = 0; k < prev_pentity2_itrf_gnum_and_itrf_strid[i_part][idx_read]; ++k) {
            idx_read_data++;
          }

          idx_read++;
        }
      }
    }

    // PDM_log_trace_array_int(pextract_entity2_kind[i_part], pextract_entity1_entity2_n_elmt[i_part], "pextract_entity2_kind ::");

    free(prev_pentity2_itrf_gnum_and_itrf_strid[i_part]);
    free(prev_pentity2_itrf_gnum_and_itrf_data [i_part]);

    /* Sort by kind and keep order */
    recv_buffer_to_sort_kind_order[i_part] = malloc(pextract_entity1_entity2_n_elmt[i_part] * sizeof(int));
    for(int i = 0; i < pextract_entity1_entity2_n_elmt[i_part]; ++i) {
      recv_buffer_to_sort_kind_order[i_part][i] = i;
    }

    /* Sort, keep order and setup idx */
    PDM_sort_int(pextract_entity2_kind[i_part], recv_buffer_to_sort_kind_order[i_part], pextract_entity1_entity2_n_elmt[i_part]);

    pentity2_extented_by_kind_idx[i_part] = malloc((n_extented_kind+1) * sizeof(int));
    int *pentity2_extented_by_kind_n = PDM_array_zeros_int(n_extented_kind);

    for(int i = 0; i < pextract_entity1_entity2_n_elmt[i_part]; ++i) {
      pentity2_extented_by_kind_n[pextract_entity2_kind[i_part][i]]++;
    }

    pentity2_extented_by_kind_idx[i_part][0] = 0;
    for(int i = 0; i < n_extented_kind; ++i) {
      pentity2_extented_by_kind_idx[i_part][i+1] = pentity2_extented_by_kind_idx[i_part][i] + pentity2_extented_by_kind_n[i];
    }

    PDM_log_trace_array_int(pentity2_extented_by_kind_idx[i_part], n_extented_kind, "pentity2_extented_by_kind_idx ::");

    free(pentity2_extented_by_kind_n);
  }

  free(prev_pentity2_itrf_gnum_and_itrf_strid);
  free(prev_pentity2_itrf_gnum_and_itrf_data );
  free(pextract_entity1_entity2_n_elmt);

  /*
   * Create new global numbering for all new entities
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     n_part,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2);

  /* Create des nouveaux gnum */
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_interface  = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int *_pentity2_extented_by_kind_idx = pentity2_extented_by_kind_idx[i_part];

    /* Prepare copy */
    int n_gnum = _pentity2_extented_by_kind_idx[n_extented_kind] - _pentity2_extented_by_kind_idx[n_extented_kind-1];
    pentity2_extented_ln_to_gn_by_interface[i_part] = malloc(2 * n_gnum * sizeof(PDM_g_num_t));

    PDM_g_num_t *_pextract_entity1_entity2_gnum           = pextract_entity1_entity2_gnum          [i_part];
    int         *_recv_buffer_to_sort_kind_order          = recv_buffer_to_sort_kind_order         [i_part];
    PDM_g_num_t *_pentity2_extented_ln_to_gn_by_interface = pentity2_extented_ln_to_gn_by_interface[i_part];
    int         *_pextract_entity2_translate              = pextract_entity2_translate             [i_part];
    int         *_pextract_entity2_interface              = pextract_entity2_interface             [i_part];

    int idx_write = 0;
    for(int i = _pentity2_extented_by_kind_idx[n_extented_kind-1]; i < _pentity2_extented_by_kind_idx[n_extented_kind]; ++i) {
      int idx_read = _recv_buffer_to_sort_kind_order[i];
      _pentity2_extented_ln_to_gn_by_interface[2*idx_write  ] = _pextract_entity1_entity2_gnum[idx_read];
      _pentity2_extented_ln_to_gn_by_interface[2*idx_write+1] = _pextract_entity2_interface   [idx_read];
      idx_write++;
    }

    if (1==1) {
      log_trace("\n");
      log_trace("ngum = %d\n", n_gnum);
      log_trace("pextract_entity1_entity2_n_elmt[i_part] = %d\n", pextract_entity1_entity2_n_elmt[i_part]);
      PDM_log_trace_array_long(_pextract_entity1_entity2_gnum, pextract_entity1_entity2_n_elmt[i_part], "pextract_entity1_entity2_gnum :: ");
      PDM_log_trace_array_long(_pextract_entity2_translate, pextract_entity1_entity2_n_elmt[i_part] , "_pextract_entity2_translate :: ");
      PDM_log_trace_array_long(_pentity2_extented_ln_to_gn_by_interface, 2 * n_gnum, "_pentity2_extented_ln_to_gn_by_interface :: ");
    }

    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              n_gnum,
                              _pentity2_extented_ln_to_gn_by_interface);

  }

  PDM_gnum_compute(gen_gnum_entity2);


  /* Post-treatment */
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int         *_pentity2_extented_by_kind_idx = pentity2_extented_by_kind_idx          [i_part];
    int n_gnum = _pentity2_extented_by_kind_idx[n_extented_kind] - _pentity2_extented_by_kind_idx[n_extented_kind-1];

    PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    PDM_g_num_t *_pextract_entity1_entity2_gnum           = pextract_entity1_entity2_gnum          [i_part];
    int         *_recv_buffer_to_sort_kind_order          = recv_buffer_to_sort_kind_order         [i_part];

    // PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, n_gnum, "extented_from_itrf_entity2_ln_to_gn :: ");

    /* Update array */
    int idx_write = 0;
    for(int i = _pentity2_extented_by_kind_idx[n_extented_kind-1]; i < _pentity2_extented_by_kind_idx[n_extented_kind]; ++i) {
      int idx_read =_recv_buffer_to_sort_kind_order[i];
      _pextract_entity1_entity2_gnum[idx_read] = extented_from_itrf_entity2_ln_to_gn[idx_write++] + shift_by_domain_entity2;
    }
  }

  /*
   * At this stage all gnum are generated and properly shifted
   * You need to update all incoming conectivity
   * For each kind (if new), we unique all entries in order to setup a new local order
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("--------------------------------\n");
    log_trace("update connectivity of entities1\n");
  }
  int          *pn_entity2_extented              = malloc(n_part * sizeof(int          ));
  PDM_g_num_t **pentity2_extented_ln_to_gn       = malloc(n_part * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_triplet        = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_interface      = malloc(n_part * sizeof(int         *));
  int         **pentity2_extented_to_entity2_idx = malloc(n_part * sizeof(int         *));
  int         **pextract_entity1_entity2         = malloc(n_part * sizeof(int         *));


  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_entity2_extented[i_part] = 0;



    int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
    int         *_pentity2_extented_by_kind_idx    = pentity2_extented_by_kind_idx   [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
    int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];
    int         *_recv_buffer_to_sort_kind_order   = recv_buffer_to_sort_kind_order  [i_part];
    int         *_pextract_entity2_translate       = pextract_entity2_translate      [i_part];

    int beg_entry = _pentity2_extented_by_kind_idx[2];
    int end_entry = _pentity2_extented_by_kind_idx[5];

    int n_entry = end_entry - beg_entry;

    int idx_write = 0;
    PDM_g_num_t *tmp_ln_to_gn = malloc(n_entry * sizeof(PDM_g_num_t));
    for(int i = beg_entry; i < end_entry; ++i) {
      int idx_read = _recv_buffer_to_sort_kind_order[i];
      tmp_ln_to_gn[idx_write++] = _pextract_entity1_entity2_gnum[idx_read];
    }
    int *unique_order_entity2 = malloc(n_entry * sizeof(int));

    int n_unique = PDM_inplace_unique_long2(tmp_ln_to_gn,
                                            unique_order_entity2,
                                            0,
                                            n_entry-1);


    for(int i = 0; i < n_entry; ++i) {
      unique_order_entity2[i] += pn_entity2_extented[i_part];
    }
    pn_entity2_extented[i_part] += n_unique;
    free(tmp_ln_to_gn);


    /* Allocate array */
    pentity2_extented_to_entity2_idx[i_part] = malloc(   (pn_entity2_extented[i_part]+1) * sizeof(int        ));
    pentity2_extented_ln_to_gn      [i_part] = malloc(    pn_entity2_extented[i_part]    * sizeof(PDM_g_num_t));
    pentity2_extented_triplet       [i_part] = malloc(3 * pn_entity2_extented[i_part]    * sizeof(int        ));
    pentity2_extented_interface     [i_part] = malloc(    pn_entity2_extented[i_part]    * sizeof(int        ));

    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;

    pextract_entity1_entity2[i_part] = malloc( _pextract_entity1_entity2_idx[n_part1_to_part2] * sizeof(int));
    int *_pextract_entity1_entity2  = pextract_entity1_entity2 [i_part];

    int         *_pentity2_extented_to_entity2_idx = pentity2_extented_to_entity2_idx[i_part];
    PDM_g_num_t *_pentity2_extented_ln_to_gn       = pentity2_extented_ln_to_gn      [i_part];
    int         *_pentity2_extented_triplet        = pentity2_extented_triplet       [i_part];
    int         *_pentity2_extented_interface      = pentity2_extented_interface     [i_part];
    int         *_pextract_entity2_interface       = pextract_entity2_interface      [i_part];

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented[i_part]+1; ++i_entity2) {
      _pentity2_extented_to_entity2_idx[i_entity2] = 3 * i_entity2;
    }

    /* Setup all new connexion and update connectivity */
    for(int i_kind = 0; i_kind < n_extented_kind; ++i_kind) {

      if(i_kind == 0 || i_kind == 1) {

        for(int i = _pentity2_extented_by_kind_idx[i_kind]; i < _pentity2_extented_by_kind_idx[i_kind+1]; ++i) {
          int idx_read = _recv_buffer_to_sort_kind_order[i];
          _pextract_entity1_entity2[idx_read] = _pextract_entity2_translate[idx_read]+1; // TODO : SIGN
        }

        // juste update
        continue;
      }

      int         *_unique_order_entity2          =  unique_order_entity2;

      for(int i = _pentity2_extented_by_kind_idx[i_kind]; i < _pentity2_extented_by_kind_idx[i_kind+1]; ++i) {
        int idx_read = _recv_buffer_to_sort_kind_order[i];
        int i_unique = _unique_order_entity2[i-_pentity2_extented_by_kind_idx[2]];

        log_trace("[%i] - i_unique = %i / %i \n", i, i_unique, pn_entity2_extented[i_part]);

        /* Can be erase multiple time */
        _pentity2_extented_ln_to_gn [  i_unique  ] = _pextract_entity1_entity2_gnum   [  idx_read  ];
        _pentity2_extented_triplet  [3*i_unique  ] = _pextract_entity1_entity2_triplet[3*idx_read  ];
        _pentity2_extented_triplet  [3*i_unique+1] = _pextract_entity1_entity2_triplet[3*idx_read+1];
        _pentity2_extented_triplet  [3*i_unique+2] = _pextract_entity1_entity2_triplet[3*idx_read+2];
        _pentity2_extented_interface[  i_unique  ] = _pextract_entity2_interface      [  idx_read  ];

        _pextract_entity1_entity2[idx_read] = (pn_entity2[i_part] + i_unique + 1); // TODO : SIGN

      }

    }

    if(0 == 1) {
      PDM_log_trace_array_int (pentity2_extented_to_entity2_idx[i_part],     pn_entity2_extented[i_part]+1, "pentity2_extented_to_entity2_idx ::" );
      PDM_log_trace_array_int (pentity2_extented_triplet       [i_part], 3 * pn_entity2_extented[i_part]  , "pentity2_extented_triplet    ::" );
      PDM_log_trace_array_int (pentity2_extented_interface     [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_interface  ::" );
      PDM_log_trace_array_long(pentity2_extented_ln_to_gn      [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_ln_to_gn   ::" );
    }


    free(unique_order_entity2);

  }



  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pextract_entity2_translate[i_part]);
    free(pextract_entity2_interface[i_part]);
  }
  free(pextract_entity2_translate);
  free(pextract_entity2_interface);


  /* Update new block of interface */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------\n");
    log_trace("update entity1 interfaces\n");
  }

  PDM_g_num_t **gnum_itrf_link                     = malloc((n_part+1) * sizeof(PDM_g_num_t *));
  PDM_g_num_t **gnum_and_itrf_data                 = malloc((n_part+1) * sizeof(PDM_g_num_t *));
  int         **gnum_and_itrf_stri                 = malloc((n_part+1) * sizeof(int         *));
  double      **gnum_and_itrf_weight               = malloc((n_part+1) * sizeof(double      *));
  int          *pn_entity2_only_by_interface_twice = malloc((n_part+1) * sizeof(int          ));


  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Only new kind by interface need to update table
    int         *_pentity2_extented_by_kind_idx = pentity2_extented_by_kind_idx          [i_part];
    int n_gnum = _pentity2_extented_by_kind_idx[n_extented_kind] - _pentity2_extented_by_kind_idx[n_extented_kind-1];

    PDM_g_num_t* _pentity2_extented_ln_to_gn_by_interface = pentity2_extented_ln_to_gn_by_interface[i_part];
    PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    pn_entity2_only_by_interface_twice[i_part] = 2 * n_gnum;

    gnum_itrf_link      [i_part] = malloc(2 * n_gnum * sizeof(PDM_g_num_t));
    gnum_and_itrf_data  [i_part] = malloc(4 * n_gnum * sizeof(PDM_g_num_t));
    gnum_and_itrf_stri  [i_part] = malloc(2 * n_gnum * sizeof(int        ));
    gnum_and_itrf_weight[i_part] = malloc(2 * n_gnum * sizeof(double     ));

    for(int i = 0; i < n_gnum; ++i) {

      PDM_g_num_t orig_gnum = _pentity2_extented_ln_to_gn_by_interface[2*i  ];
      PDM_g_num_t i_itrf    = _pentity2_extented_ln_to_gn_by_interface[2*i+1];

      PDM_g_num_t new_gnum = extented_from_itrf_entity2_ln_to_gn[i] + shift_by_domain_entity2;

      int j = n_gnum+i;

      /* Remplissage */
      gnum_itrf_link[i_part][i] = orig_gnum;
      gnum_itrf_link[i_part][j] = new_gnum;

      gnum_and_itrf_data[i_part][2*i  ] = new_gnum;
      gnum_and_itrf_data[i_part][2*i+1] = i_itrf;

      gnum_and_itrf_data[i_part][2*j  ] = orig_gnum;
      gnum_and_itrf_data[i_part][2*j+1] = -i_itrf;

      gnum_and_itrf_stri[i_part][i] = 1;
      gnum_and_itrf_stri[i_part][j] = 1;

      gnum_and_itrf_weight[i_part][i] = 1.;
      gnum_and_itrf_weight[i_part][j] = 1.;

    }
  }

  /*
   * Merge with previous
   */
  gnum_itrf_link                    [n_part] = prev_dentity2_itrf_blk_gnum;
  gnum_and_itrf_data                [n_part] = prev_dentity2_itrf_gnum_and_itrf_data;
  gnum_and_itrf_stri                [n_part] = prev_dentity2_itrf_gnum_and_itrf_strid;
  gnum_and_itrf_weight              [n_part] = malloc(prev_dentity2_itrf_n_blk * sizeof(double));
  pn_entity2_only_by_interface_twice[n_part] = prev_dentity2_itrf_n_blk;
  for(int i = 0; i < prev_dentity2_itrf_n_blk; ++i) {
    gnum_and_itrf_weight[n_part][i] = 1.;
  }


  PDM_part_to_block_t* ptb_itrf = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           gnum_itrf_link,
                                                           gnum_and_itrf_weight,
                                                           pn_entity2_only_by_interface_twice,
                                                           n_part+1,
                                                           comm);
  free(gnum_and_itrf_weight[n_part]);

  int         *dentity2_gnum_and_itrf_strid = NULL;
  PDM_g_num_t *dentity2_gnum_and_itrf_data  = NULL;
  PDM_part_to_block_exch(ptb_itrf,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         gnum_and_itrf_stri,
          (void **)      gnum_and_itrf_data,
                         &dentity2_gnum_and_itrf_strid,
          (void **)      &dentity2_gnum_and_itrf_data);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(gnum_itrf_link      [i_part]);
    free(gnum_and_itrf_data  [i_part]);
    free(gnum_and_itrf_stri  [i_part]);
    free(gnum_and_itrf_weight[i_part]);
  }
  free(gnum_itrf_link      );
  free(gnum_and_itrf_data  );
  free(gnum_and_itrf_stri  );
  free(gnum_and_itrf_weight);
  free(pn_entity2_only_by_interface_twice);

  int n_blk_elt_itrf = PDM_part_to_block_n_elt_block_get(ptb_itrf);
  PDM_g_num_t *next_dentity2_elt_gnum_ptp_itrf = PDM_part_to_block_block_gnum_get   (ptb_itrf);

  PDM_g_num_t *dentity2_itrf_blk_gnum = (PDM_g_num_t *) malloc(n_blk_elt_itrf * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_blk_elt_itrf; ++i) {
    dentity2_itrf_blk_gnum[i] = next_dentity2_elt_gnum_ptp_itrf[i];
  }

  if(0 == 1) {
    PDM_log_trace_array_long(next_dentity2_elt_gnum_ptp_itrf       ,     n_blk_elt_itrf, "next_dentity2_elt_gnum      ::");
    int n_data = 0;
    for(int i = 0; i < n_blk_elt_itrf; ++i) {
      n_data += dentity2_gnum_and_itrf_strid[i];
    }
    PDM_log_trace_array_int (dentity2_gnum_and_itrf_strid,     n_blk_elt_itrf, "dentity2_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(dentity2_gnum_and_itrf_data , 2 * n_data        , "dentity2_gnum_and_itrf_data  ::");
  }

  int max_strid = 0;
  for(int i = 0; i < n_blk_elt_itrf; ++i) {
    max_strid = PDM_MAX(max_strid, dentity2_gnum_and_itrf_strid[i]);
  }
  int *data_order = malloc(max_strid * sizeof(int));

  /* Unique incomming data */
  int idx_read      = 0;
  int idx_blk_write = 0;
  for(int i = 0; i < n_blk_elt_itrf; ++i) {

    int n_data = dentity2_gnum_and_itrf_strid[i];

    int n_unique = PDM_order_inplace_unique_long(n_data, 2, &dentity2_gnum_and_itrf_data[idx_read], data_order);

    /* Copy */
    for(int i_unique = 0; i_unique < n_unique; ++i_unique) {
      dentity2_gnum_and_itrf_data[2*idx_blk_write  ] = dentity2_gnum_and_itrf_data[(idx_read+2*i_unique)  ];
      dentity2_gnum_and_itrf_data[2*idx_blk_write+1] = dentity2_gnum_and_itrf_data[(idx_read+2*i_unique)+1];
      idx_blk_write++;
    }

    dentity2_gnum_and_itrf_strid[i] = n_unique;
    idx_read += 2 * n_data;
  }
  free(data_order);

  dentity2_gnum_and_itrf_data = realloc(dentity2_gnum_and_itrf_data, 2 * idx_blk_write * sizeof(PDM_g_num_t));

  if(0 == 1) {
    PDM_log_trace_array_int (dentity2_gnum_and_itrf_strid,     n_blk_elt_itrf, "dentity2_gnum_and_itrf_strid (Unique) ::");
    PDM_log_trace_array_long(dentity2_gnum_and_itrf_data , 2 * idx_blk_write , "dentity2_gnum_and_itrf_data  (Unique) ::");
  }

  PDM_part_to_block_free(ptb_itrf);

  /*
   * Keep pointer for output
   */
  *next_dentity2_itrf_n_blk_out               = n_blk_elt_itrf;
  *next_dentity2_itrf_blk_gnum_out            = dentity2_itrf_blk_gnum;
  *next_dentity2_itrf_gnum_and_itrf_strid_out = dentity2_gnum_and_itrf_strid;
  *next_dentity2_itrf_gnum_and_itrf_data_out  = dentity2_gnum_and_itrf_data;


  PDM_gnum_free(gen_gnum_entity2);

  /* Free all working array */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_extented_by_kind_idx          [i_part]);
    free(pextract_entity2_kind                  [i_part]);
    free(recv_buffer_to_sort_kind_order         [i_part]);
    free(pentity2_extented_ln_to_gn_by_interface[i_part]);
  }

  free(pentity2_extented_by_kind_idx);
  free(pextract_entity2_kind);
  free(recv_buffer_to_sort_kind_order);
  free(pentity2_extented_ln_to_gn_by_interface);

  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pextented_entity1_entity2_idx_out           = pextract_entity1_entity2_idx;
  *pextented_entity1_entity2_out               = pextract_entity1_entity2;
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_entity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_interface;

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(entity2_order           [i_part]);
    free(pentity2_ln_to_gn_sorted[i_part]);
  }
  free(entity2_order           );
  free(pentity2_ln_to_gn_sorted);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pextract_entity1_entity2_n      [i_part]);
    free(pextract_entity1_entity2_gnum   [i_part]);
    free(pextract_entity1_entity2_triplet[i_part]);
    // free(pextract_entity1_entity2_idx [i_part]);
  }
  free(pextract_entity1_entity2_n      );
  free(pextract_entity1_entity2_gnum   );
  free(pextract_entity1_entity2_triplet);


  return;
  // abort();

  // /*
  //  * We have two kind of extented :
  //  *   - From partition
  //  *   - From interfaces (request new gnum génération ...)
  //  */
  // int          *pn_entity2_extented_by_interface         = malloc(n_part * sizeof(int          ));
  // int          *pn_entity2_extented_by_partition         = malloc(n_part * sizeof(int          ));
  // // PDM_g_num_t **pentity2_extented_ln_to_gn_by_interface  = malloc(n_part * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **pentity2_extented_ln_to_gn_by_partition  = malloc(n_part * sizeof(PDM_g_num_t *));
  // int         **pentity2_extented_triplet_by_interface   = malloc(n_part * sizeof(int         *));
  // int         **pentity2_extented_triplet_by_partition   = malloc(n_part * sizeof(int         *));
  // int         **pentity2_extented_interface_by_interface = malloc(n_part * sizeof(int         *));
  // int         **pentity2_extented_interface_by_partition = malloc(n_part * sizeof(int         *));

  // for(int i_part = 0; i_part < n_part; ++i_part) {

  //   int         *_pextract_entity1_interface = pentity1_extented_to_pentity1_interface[i_part];
  //   int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
  //   PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
  //   int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

  //   PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];

  //   pn_entity2_extented_by_interface[i_part] = 0;
  //   pn_entity2_extented_by_partition[i_part] = 0;

  //   /* Count */
  //   int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
  //   for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
  //     if(_pextract_entity1_interface[i_ref] != 0) {

  //       int i_itrf = _pextract_entity1_interface[i_ref];

  //       /* Blinder d'assert */
  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         /* Old manner */
  //         // PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         // PDM_g_num_t gnum_to_find[2] = {entity2_g_num, i_itrf};

  //         // /* First look on interface array */
  //         // int pos = -1; // PDM_order_binary_search_long(gnum_to_find, pentity2_ext_opp_gnum_and_itrf[i_part], 2, pn_entity2_ext_opp_gnum_and_itrf[i_part]);
  //         // abort();
  //         // // log_trace("gnum_to_find = %i/%i --> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);

  //         // if(pos == -1) {
  //         //   pn_entity2_extented_by_interface[i_part]++;
  //         // }

  //         if(pextract_entity2_kind[i_part][idx_entity1] == 1) {
  //           pn_entity2_extented_by_interface[i_part]++;
  //         }

  //       }
  //     } else { // Second case : Interior entity1

  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

  //         if(pos_int == -1) { // Not in current partition
  //           pn_entity2_extented_by_partition[i_part]++;
  //         }
  //       }
  //     }
  //   }

  //   /* Fill */
  //   pentity2_extented_ln_to_gn_by_interface [i_part] = malloc(2 * pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
  //   pentity2_extented_ln_to_gn_by_partition [i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
  //   pentity2_extented_triplet_by_interface  [i_part] = malloc(3 * pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
  //   pentity2_extented_triplet_by_partition  [i_part] = malloc(3 * pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
  //   pentity2_extented_interface_by_interface[i_part] = malloc(    pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
  //   pentity2_extented_interface_by_partition[i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(int        ));

  //   pn_entity2_extented_by_interface[i_part] = 0;
  //   pn_entity2_extented_by_partition[i_part] = 0;
  //   for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
  //     if(_pextract_entity1_interface[i_ref] != 0) {

  //       /* Blinder d'assert */
  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};
  //         // Old manner
  //         // int pos = -1; // PDM_order_binary_search_long(gnum_to_find, pentity2_ext_opp_gnum_and_itrf[i_part], 2, pn_entity2_ext_opp_gnum_and_itrf[i_part]);
  //         // abort();
  //         // if(pos == -1) {
  //         if(pextract_entity2_kind[i_part][idx_entity1] == 1) {
  //           int idx_write = pn_entity2_extented_by_interface[i_part]++;
  //           pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write  ] = entity2_g_num;
  //           pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write+1] = _pextract_entity1_interface[i_ref];

  //           pentity2_extented_triplet_by_interface  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
  //           pentity2_extented_triplet_by_interface  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
  //           pentity2_extented_triplet_by_interface  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
  //           pentity2_extented_interface_by_interface[i_part][  idx_write  ] = _pextract_entity1_interface[i_ref];
  //         }
  //       }
  //     } else { // Second case : Interior entity1

  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

  //         if(pos_int == -1) { // Not in current partition
  //           int idx_write = pn_entity2_extented_by_partition[i_part]++;
  //           pentity2_extented_ln_to_gn_by_partition [i_part][idx_write] = entity2_g_num;
  //           pentity2_extented_triplet_by_partition  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
  //           pentity2_extented_triplet_by_partition  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
  //           pentity2_extented_triplet_by_partition  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
  //           pentity2_extented_interface_by_partition[i_part][  idx_write  ] = 0; // Because interior
  //         }
  //       }
  //     }
  //   }

  //   // if(1 == 1) {
  //   //   PDM_log_trace_array_long(pentity2_extented_ln_to_gn_by_interface[i_part],
  //   //                            2 * pn_entity2_extented_by_interface[i_part],
  //   //                            "pentity2_extented_ln_to_gn_by_interface ::");
  //   // }

  //   PDM_gnum_set_from_parents(gen_gnum_entity2,
  //                             i_part,
  //                             pn_entity2_extented_by_interface[i_part],
  //                             pentity2_extented_ln_to_gn_by_interface[i_part]);

  // }

  // PDM_gnum_compute(gen_gnum_entity2);

  // /*
  //  * Unify entity1_entity2
  //  */
  // // int **pextract_entity1_entity2 = malloc(n_part * sizeof(int *));
  // for(int i_part = 0; i_part < n_part; ++i_part) {

  //   int         *_pextract_entity1_interface  = pentity1_extented_to_pentity1_interface[i_part];

  //   int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
  //   PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
  //   int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

  //   PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];
  //   int         *_entity2_order              = entity2_order           [i_part];

  //   int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;


  //   pextract_entity1_entity2[i_part] = malloc( _pextract_entity1_entity2_idx[n_part1_to_part2] * sizeof(int));
  //   int *_pextract_entity1_entity2  = pextract_entity1_entity2 [i_part];

  //   PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

  //   if(1 == 0) {
  //     PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn :: ");
  //   }

  //   // Erase and realloc :
  //   pentity2_extented_ln_to_gn_by_interface[i_part] = realloc(pentity2_extented_ln_to_gn_by_interface[i_part], pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
  //     pentity2_extented_ln_to_gn_by_interface[i_part][i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2] + shift_by_domain_entity2;
  //   }

  //   PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn ::");

  //   /* Sort unique the new gnum to unify */
  //   PDM_g_num_t *extented_from_itrf_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
  //     extented_from_itrf_entity2_ln_to_gn_sorted[i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2];
  //   }
  //   // PDM_sort_long(extented_from_itrf_entity2_ln_to_gn_sorted, extented_from_itrf_entity2_order, pn_entity2_extented_by_interface[i_part]);
  //   pn_entity2_extented_by_interface[i_part] = PDM_inplace_unique_long(extented_from_itrf_entity2_ln_to_gn_sorted,
  //                                                                      NULL,
  //                                                                      0,
  //                                                                      pn_entity2_extented_by_interface[i_part]-1);
  //   extented_from_itrf_entity2_ln_to_gn_sorted = realloc(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
  //   int         *extented_from_itrf_entity2_order = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
  //     extented_from_itrf_entity2_order[i_entity2] = -1;
  //   }

  //   printf("pn_entity2_extented_by_interface[i_part] = %i \n", pn_entity2_extented_by_interface[i_part]);


  //   /* Sort unique the new gnum to unify */
  //   PDM_g_num_t *extented_from_part_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
  //     extented_from_part_entity2_ln_to_gn_sorted[i_entity2] = pentity2_extented_ln_to_gn_by_partition[i_part][i_entity2];
  //     // extented_from_part_entity2_order          [i_entity2] = i_entity2;
  //   }
  //   // PDM_sort_long(extented_from_part_entity2_ln_to_gn_sorted, extented_from_part_entity2_order, pn_entity2_extented_by_partition[i_part]);

  //   pn_entity2_extented_by_partition[i_part] = PDM_inplace_unique_long(extented_from_part_entity2_ln_to_gn_sorted,
  //                                                                      NULL,
  //                                                                      0,
  //                                                                      pn_entity2_extented_by_partition[i_part]-1);
  //   int *extented_from_part_entity2_order = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
  //   extented_from_part_entity2_ln_to_gn_sorted = realloc(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
  //     extented_from_part_entity2_order[i_entity2] = -1;
  //   }

  //   if(1 == 1) {
  //     printf("pn_entity2_extented_by_partition[i_part] = %i \n", pn_entity2_extented_by_partition[i_part]);
  //     PDM_log_trace_array_long(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part], "extented_from_part_entity2_ln_to_gn_sorted ::");
  //     PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn_sorted ::");
  //   }

  //   /* */
  //   int pn_entity2_extented_by_interface2 = 0; // To read in extented_from_itrf_entity2_ln_to_gn
  //   for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {

  //     log_trace("_pextract_entity1_interface[i_ref] = %i \n", _pextract_entity1_interface[i_ref] );

  //     /*
  //      * First case :
  //      *   - entity1 is move by interface
  //      */
  //     if(_pextract_entity1_interface[i_ref] != 0) {
  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};

  //         int pos = -1; // PDM_order_binary_search_long(gnum_to_find, pentity2_ext_opp_gnum_and_itrf[i_part], 2, pn_entity2_ext_opp_gnum_and_itrf[i_part]);
  //         abort();
  //         // log_trace("\t Search %i/%i in _pentity2_opp_gnum_and_itrf -> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);

  //         if(pos == -1) {
  //           /*
  //            * Subcase :
  //            *   - entity2 is not in table of interface : it's a new entity2
  //            */
  //           PDM_g_num_t entity2_extented_g_num = extented_from_itrf_entity2_ln_to_gn[pn_entity2_extented_by_interface2++];
  //           int pos_new = PDM_binary_search_long(entity2_extented_g_num,
  //                                                extented_from_itrf_entity2_ln_to_gn_sorted,
  //                                                pn_entity2_extented_by_interface[i_part]);

  //           assert(pos_new != -1);
  //           int shift = pn_entity2[i_part] + pn_entity2_extented_by_partition[i_part];

  //           if(extented_from_itrf_entity2_order[pos_new] == -1) {
  //             extented_from_itrf_entity2_order[pos_new] = 1;

  //             pentity2_extented_triplet_by_interface  [i_part][3*pos_new  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
  //             pentity2_extented_triplet_by_interface  [i_part][3*pos_new+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
  //             pentity2_extented_triplet_by_interface  [i_part][3*pos_new+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
  //             pentity2_extented_interface_by_interface[i_part][  pos_new  ] = _pextract_entity1_interface[i_ref];

  //           }

  //           // _pextract_entity1_entity2[idx_entity1] = ( shift + extented_from_itrf_entity2_order[pos_new] + 1); // ATTENTION SIGN
  //           _pextract_entity1_entity2[idx_entity1] = ( shift + pos_new + 1); // ATTENTION SIGN

  //           log_trace("\t Translate cas 1 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

  //         } else {
  //           /*
  //            * Subcase :
  //            *   - entity2 is in table of interface
  //            */
  //           abort();
  //           // PDM_g_num_t cur_gnum = pentity2_ext_cur_gnum[i_part][pos];
  //           // int pos2 = PDM_binary_search_long(cur_gnum, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);
  //           // // int pos2      = entity2_opp_position[i_part][pos ];
  //           // _pextract_entity1_entity2[idx_entity1] = ( pos2 + 1);

  //           // log_trace("\t Translate cas 2 : %i  ---> (idx_entity1 = %i) --> %i (pos=%i/pos2=%i) \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1], pos, pos2);


  //         }
  //       }
  //     } else { // Second case : Interior entity1

  //       for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

  //         PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
  //         int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

  //         if(pos_int != -1) {                   // entity2 is already in current partition
  //           assert(pos_int != -1);
  //           int orig_pos = _entity2_order[pos_int];
  //           _pextract_entity1_entity2[idx_entity1] = ( orig_pos + 1);

  //           log_trace("\t Translate cas 3 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

  //         } else {                              // Not in current partition
  //           int pos_ext = PDM_binary_search_long(entity2_g_num,
  //                                                extented_from_part_entity2_ln_to_gn_sorted,
  //                                                pn_entity2_extented_by_partition[i_part]);
  //           assert(pos_ext != -1);

  //           if(extented_from_part_entity2_order[pos_ext] == -1) {
  //             extented_from_part_entity2_order[pos_ext] = 1;
  //             pentity2_extented_ln_to_gn_by_partition [i_part][  pos_ext  ] = entity2_g_num;

  //             pentity2_extented_triplet_by_partition  [i_part][3*pos_ext  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
  //             pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
  //             pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
  //             pentity2_extented_interface_by_partition[i_part][  pos_ext  ] = 0;

  //           }

  //           // _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + extented_from_part_entity2_order[pos_ext] + 1);
  //           _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + pos_ext + 1);

  //           log_trace("\t Translate cas 4 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

  //         }
  //       }
  //     }
  //   }

  //   // assert(pn_entity2_extented_by_interface2 == pn_entity2_extented_by_interface[i_part]);

  //   if(1 == 1) {
  //     PDM_log_trace_connectivity_int(_pextract_entity1_entity2_idx,
  //                                    _pextract_entity1_entity2,
  //                                    n_part1_to_part2,
  //                                    "pextract_entity1_entity2 ::");
  //   }

  //   free(extented_from_itrf_entity2_order          );
  //   free(extented_from_itrf_entity2_ln_to_gn_sorted);
  //   free(extented_from_part_entity2_order          );
  //   free(extented_from_part_entity2_ln_to_gn_sorted);

  // }

  // PDM_gnum_free(gen_gnum_entity2);


  // /*
  //  * Reconstruction
  //  */
  // // int          *pn_entity2_extented              = malloc(n_part * sizeof(int          ));
  // // PDM_g_num_t **pentity2_extented_ln_to_gn       = malloc(n_part * sizeof(PDM_g_num_t *));
  // // int         **pentity2_extented_triplet        = malloc(n_part * sizeof(int         *));
  // // int         **pentity2_extented_interface      = malloc(n_part * sizeof(int         *));
  // // int         **pentity2_extented_to_entity2_idx = malloc(n_part * sizeof(int         *));

  // for(int i_part = 0; i_part < n_part; ++i_part) {

  //   pn_entity2_extented[i_part] = pn_entity2_extented_by_partition[i_part] + pn_entity2_extented_by_interface[i_part];
  //   pentity2_extented_ln_to_gn      [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(PDM_g_num_t));
  //   pentity2_extented_triplet       [i_part] = malloc(3 * pn_entity2_extented[i_part]      * sizeof(int        ));
  //   pentity2_extented_interface     [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(int        ));
  //   pentity2_extented_to_entity2_idx[i_part] = malloc( (  pn_entity2_extented[i_part] + 1) * sizeof(int        ));

  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
  //     pentity2_extented_ln_to_gn [i_part][  i_entity2  ] = pentity2_extented_ln_to_gn_by_partition [i_part][  i_entity2  ];
  //     pentity2_extented_triplet  [i_part][3*i_entity2  ] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2  ];
  //     pentity2_extented_triplet  [i_part][3*i_entity2+1] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+1];
  //     pentity2_extented_triplet  [i_part][3*i_entity2+2] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+2];
  //     pentity2_extented_interface[i_part][  i_entity2  ] = pentity2_extented_interface_by_partition[i_part][  i_entity2  ];
  //   }

  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
  //     int idx_write = pn_entity2_extented_by_partition[i_part]+i_entity2;
  //     pentity2_extented_ln_to_gn [i_part][  idx_write  ] = pentity2_extented_ln_to_gn_by_interface [i_part][  i_entity2  ];
  //     pentity2_extented_triplet  [i_part][3*idx_write  ] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2  ];
  //     pentity2_extented_triplet  [i_part][3*idx_write+1] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+1];
  //     pentity2_extented_triplet  [i_part][3*idx_write+2] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+2];
  //     pentity2_extented_interface[i_part][  idx_write  ] = pentity2_extented_interface_by_interface[i_part][  i_entity2  ];
  //   }

  //   for(int i_entity2 = 0; i_entity2 < pn_entity2_extented[i_part]+1; ++i_entity2) {
  //     pentity2_extented_to_entity2_idx[i_part][i_entity2] = 3 * i_entity2;
  //   }

  //   if(1 == 1) {
  //     PDM_log_trace_array_int (pentity2_extented_to_entity2_idx[i_part],     pn_entity2_extented[i_part]+1, "pentity2_extented_to_entity2_idx ::" );
  //     PDM_log_trace_array_int (pentity2_extented_triplet       [i_part], 3 * pn_entity2_extented[i_part]  , "pentity2_extented_triplet    ::" );
  //     PDM_log_trace_array_int (pentity2_extented_interface     [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_interface  ::" );
  //     PDM_log_trace_array_long(pentity2_extented_ln_to_gn      [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_ln_to_gn   ::" );
  //   }

  // }

  // *pn_entity2_extented_out                     = pn_entity2_extented;
  // *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  // *pextented_entity1_entity2_idx_out           = pextract_entity1_entity2_idx;
  // *pextented_entity1_entity2_out               = pextract_entity1_entity2;
  // *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_entity2_idx;
  // *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_triplet;
  // *pentity2_extented_to_pentity2_interface_out = pentity2_extented_interface;

  // for(int i_part = 0; i_part < n_part; ++i_part) {
  //   free(entity2_order           [i_part]);
  //   free(pentity2_ln_to_gn_sorted[i_part]);
  // }
  // free(entity2_order           );
  // free(pentity2_ln_to_gn_sorted);

  // for(int i_part = 0; i_part < n_part; ++i_part) {
  //   free(pentity2_extented_ln_to_gn_by_interface [i_part]);
  //   free(pentity2_extented_ln_to_gn_by_partition [i_part]);
  //   free(pentity2_extented_triplet_by_interface  [i_part]);
  //   free(pentity2_extented_triplet_by_partition  [i_part]);
  //   free(pentity2_extented_interface_by_interface[i_part]);
  //   free(pentity2_extented_interface_by_partition[i_part]);
  // }
  // free(pentity2_extented_ln_to_gn_by_interface );
  // free(pentity2_extented_ln_to_gn_by_partition );
  // free(pentity2_extented_triplet_by_interface  );
  // free(pentity2_extented_triplet_by_partition  );
  // free(pentity2_extented_interface_by_interface);
  // free(pentity2_extented_interface_by_partition);

  // free(pn_entity2_extented_by_interface );
  // free(pn_entity2_extented_by_partition );

  // for(int i_part = 0; i_part < n_part; ++i_part) {
  //   free(pextract_entity1_entity2_n      [i_part]);
  //   free(pextract_entity1_entity2_gnum   [i_part]);
  //   free(pextract_entity1_entity2_triplet[i_part]);
  //   // free(pextract_entity1_entity2_idx [i_part]);
  // }
  // free(pextract_entity1_entity2_n      );
  // free(pextract_entity1_entity2_gnum   );
  // free(pextract_entity1_entity2_triplet);

}



void
PDM_part_extension_pconnectivity_to_extented_pconnectivity
(
  PDM_part_domain_interface_t    *pdi,
  PDM_bound_type_t                entity2_bound,
  int                             n_domain,
  PDM_g_num_t                    *shift_by_domain_entity2,
  int                            *n_part,
  int                           **pn_entity1_in,
  PDM_g_num_t                  ***pentity1_ln_to_gn_in,
  int                           **pn_entity2_in,
  PDM_g_num_t                  ***pentity2_ln_to_gn_in,
  int                          ***pentity1_entity2_idx_in,
  int                          ***pentity1_entity2_in,
  int                            *pn_entity1_extented,
  PDM_g_num_t                   **pentity1_extented_ln_to_gn,
  int                           **pentity1_extented_to_pentity1_idx,
  int                           **pentity1_extented_to_pentity1_triplet,
  int                           **pentity1_extented_to_pentity1_interface,
  int                           **pn_entity2_extented_out,
  PDM_g_num_t                  ***pentity2_extented_ln_to_gn_out,
  int                          ***pextented_entity1_entity2_idx_out,
  int                          ***pextented_entity1_entity2_out,
  int                          ***pentity2_extented_to_pentity2_idx_out,
  int                          ***pentity2_extented_to_pentity2_triplet_out,
  int                          ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                    comm
)
{
  log_trace("---  PDM_part_extension_pconnectivity_to_extented_pconnectivity --- \n");

  // Maybe supoose on entry already concat domain in part and juste give shift ?
  PDM_UNUSED(pdi);
  PDM_UNUSED(entity2_bound);
  PDM_UNUSED(n_domain);
  PDM_UNUSED(shift_by_domain_entity2);
  PDM_UNUSED(n_part);
  PDM_UNUSED(pn_entity1_in);
  PDM_UNUSED(pentity1_ln_to_gn_in);
  PDM_UNUSED(pn_entity2_in);
  PDM_UNUSED(pentity2_ln_to_gn_in);
  PDM_UNUSED(pentity1_entity2_idx_in);
  PDM_UNUSED(pentity1_entity2_in);
  PDM_UNUSED(pn_entity1_extented);
  PDM_UNUSED(pentity1_extented_ln_to_gn);
  PDM_UNUSED(pentity1_extented_to_pentity1_idx);
  PDM_UNUSED(pentity1_extented_to_pentity1_triplet);
  PDM_UNUSED(pentity1_extented_to_pentity1_interface);
  PDM_UNUSED(pn_entity2_extented_out);
  PDM_UNUSED(pentity2_extented_ln_to_gn_out);
  PDM_UNUSED(pextented_entity1_entity2_idx_out);
  PDM_UNUSED(pextented_entity1_entity2_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_idx_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_triplet_out);
  PDM_UNUSED(pentity2_extented_to_pentity2_interface_out);
  PDM_UNUSED(comm);

  PDM_domain_interface_t  *ditrf              = NULL;
  int                    **is_entity2_on_itrf = NULL;
  PDM_part_domain_interface_to_domain_interface(pdi,
                                                entity2_bound,
                                                n_part,
                                                pn_entity2_in,
                                                pentity2_ln_to_gn_in,
                                                &ditrf,
                                                &is_entity2_on_itrf);

  /*
   * Create flat view to have in flat frame all information for couple iterface / gnum
   */
  PDM_part_to_block_t **ptb_interface_entity2      = NULL;
  PDM_g_num_t         **interface_entity2_opp_gnum = NULL;
  int                 **interface_entity2_opp_sens = NULL;
  PDM_domain_interface_make_flat_view(ditrf,
                                      entity2_bound,
                                      shift_by_domain_entity2,
                                      &ptb_interface_entity2,
                                      &interface_entity2_opp_gnum,
                                      &interface_entity2_opp_sens);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  int n_interface = PDM_part_domain_interface_n_interface_get(pdi);

  for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
    free(interface_entity2_opp_sens[i_itrf]);
  }
  free(interface_entity2_opp_sens);


  int n_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_part_tot += n_part[i_dom];
  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    free(is_entity2_on_itrf[i_part]);
  }
  free(is_entity2_on_itrf);

  PDM_domain_interface_free(ditrf);

  int          *pn_entity1           = (int          * ) malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity1_ln_to_gn    = (PDM_g_num_t ** ) malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity1_entity2_idx = (int         ** ) malloc(n_part_tot * sizeof(int         *));
  int         **pentity1_entity2     = (int         ** ) malloc(n_part_tot * sizeof(int         *));

  int          *pn_entity2        = (int          * ) malloc(n_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ln_to_gn = (PDM_g_num_t ** ) malloc(n_part_tot * sizeof(PDM_g_num_t *));

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      pn_entity1          [ln_part_tot+i_part] = pn_entity1_in          [i_dom][i_part];
      pentity1_ln_to_gn   [ln_part_tot+i_part] = pentity1_ln_to_gn_in   [i_dom][i_part];
      pentity1_entity2_idx[ln_part_tot+i_part] = pentity1_entity2_idx_in[i_dom][i_part];
      pentity1_entity2    [ln_part_tot+i_part] = pentity1_entity2_in    [i_dom][i_part];
      pn_entity2          [ln_part_tot+i_part] = pn_entity2_in          [i_dom][i_part];
      pentity2_ln_to_gn   [ln_part_tot+i_part] = pentity2_ln_to_gn_in   [i_dom][i_part];

    }
    ln_part_tot += n_part[i_dom];
  }

  /*
   * Preparation of array to manage interface at post-treatment
   */
  int          *pn_entity2_num             = NULL;
  int         **pentity2_num               = NULL;
  int         **pentity2_opp_location_idx  = NULL;
  int         **pentity2_opp_location      = NULL;
  int         **pentity2_opp_interface_idx = NULL;
  int         **pentity2_opp_interface     = NULL;
  int         **pentity2_opp_sens          = NULL;
  PDM_g_num_t **pentity2_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         entity2_bound,
                                         pn_entity2,
                                         pentity2_ln_to_gn,
                                         &pn_entity2_num,
                                         &pentity2_num,
                                         &pentity2_opp_location_idx,
                                         &pentity2_opp_location,
                                         &pentity2_opp_interface_idx,
                                         &pentity2_opp_interface,
                                         &pentity2_opp_sens,
                                         &pentity2_opp_gnum);

  /*
   *  Il faut sortir le interface ln_to_gn + le shift puis garder un multiblock to part avec les doublé + les sens
   *
   */
  // Supprmier ce bout et utiliser le flat view APRès l'écahnge

  int          *pn_entity2_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **entity2_opp_position         = malloc(ln_part_tot * sizeof(int         *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    // int *_pentity2_opp_location = pentity2_opp_location_idx[i_part];
    int *_pentity2_opp_location_idx = pentity2_opp_location_idx[i_part];
    int n_connect_tot = _pentity2_opp_location_idx[pn_entity2_num[i_part]];

    // On garde le lien courant (indirect sort of pentity2_num) + également keep le sens_opp
    pentity2_opp_gnum_and_itrf[i_part] = malloc(2 * n_connect_tot * sizeof(PDM_g_num_t));
    PDM_g_num_t* _pentity2_opp_gnum_and_itrf = pentity2_opp_gnum_and_itrf[i_part];

    // log_trace("pn_entity2_num[i_part] = %i \n", pn_entity2_num[i_part]);
    // PDM_log_trace_array_int(_pentity2_opp_location_idx, pn_entity2_num[i_part]+1, "_pentity2_opp_location_idx ::");

    int *tmp_opp_position = malloc(n_connect_tot * sizeof(int));
    for(int idx_entity = 0; idx_entity < pn_entity2_num[i_part]; ++idx_entity) {
      int i_entity = pentity2_num[i_part][idx_entity];
      for(int idx_opp = _pentity2_opp_location_idx[idx_entity]; idx_opp < _pentity2_opp_location_idx[idx_entity+1]; ++idx_opp) {

        _pentity2_opp_gnum_and_itrf[2*idx_opp  ] = pentity2_opp_gnum     [i_part][idx_opp];
        _pentity2_opp_gnum_and_itrf[2*idx_opp+1] = pentity2_opp_interface[i_part][idx_opp];

        tmp_opp_position[idx_opp] = i_entity;
      }
    }

    if(0 == 1) {
      PDM_log_trace_array_long(pentity2_opp_gnum_and_itrf[i_part], 2 * n_connect_tot, "pentity2_opp_gnum_and_itrf ::");
      PDM_log_trace_array_int(tmp_opp_position, n_connect_tot, "tmp_opp_position ::");
    }

    entity2_opp_position        [i_part] = malloc(n_connect_tot * sizeof(int));

    int* order = malloc(n_connect_tot * sizeof(int));
    pn_entity2_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_connect_tot,
                                                                                   2,
                                                                                   pentity2_opp_gnum_and_itrf[i_part],
                                                                                   order);


    if(0 == 1) {
      PDM_log_trace_array_int(order, n_connect_tot, "order = ");
      PDM_log_trace_array_long(pentity2_opp_gnum_and_itrf[i_part],
                               2 * pn_entity2_opp_gnum_and_itrf[i_part], "pentity2_opp_gnum_and_itrf ::");
      log_trace("pn_entity2_opp_gnum_and_itrf[i_part] = %i \n", pn_entity2_opp_gnum_and_itrf[i_part]);
      log_trace("n_connect_tot = %i \n", n_connect_tot);
    }

    for(int i = 0; i < pn_entity2_opp_gnum_and_itrf[i_part]; ++i) {
      entity2_opp_position[i_part][i] = tmp_opp_position[order[i]];
    }


    free(order);
    free(tmp_opp_position);

    if(0 == 1) {
      PDM_log_trace_array_int(entity2_opp_position[i_part],
                              pn_entity2_opp_gnum_and_itrf[i_part], "entity2_opp_position = ");
      PDM_log_trace_array_int(pentity2_num[i_part], pn_entity2_num[i_part], "pentity2_num ::");
    }
  }


  if(0 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pentity2_num[i_part], pn_entity2_num[i_part], "pentity2_num ::");
      PDM_log_trace_graph_nuplet_int(pentity2_opp_location_idx[i_part],
                                     pentity2_opp_location    [i_part],
                                     3,
                                     pn_entity2_num[i_part], "pentity2_opp_location ::");
      PDM_log_trace_array_int (pentity2_opp_interface[i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_interface ::");
      PDM_log_trace_array_int (pentity2_opp_sens     [i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_sens ::");
      PDM_log_trace_array_long(pentity2_opp_gnum     [i_part], pentity2_opp_location_idx[i_part][pn_entity2_num[i_part]], "pentity2_opp_gnum ::");
    }
  }

  /*
   * Creation du part_to_part entre entity2
   */
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_extented_ln_to_gn,
                                                                      (const int          *) pn_entity1_extented,
                                                                      ln_part_tot,
                                                                      (const int          *) pn_entity1,
                                                                      ln_part_tot,
                                                                      (const int         **) pentity1_extented_to_pentity1_idx,
                                                                      (const int         **) NULL,
                                                                      (const int         **) pentity1_extented_to_pentity1_triplet,
                                                                      comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  /*
   * Prepare buffer
   */
  int         **gnum1_com_from_entity1_entity2_n       = malloc(n_part_tot * sizeof(int         *));
  PDM_g_num_t **gnum1_com_from_entity1_entity2         = malloc(n_part_tot * sizeof(PDM_g_num_t *));
  int         **gnum1_com_from_entity1_entity2_triplet = malloc(n_part_tot * sizeof(int         *));
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx = pentity1_entity2_idx[i_part];
    int *_pentity1_entity2     = pentity1_entity2    [i_part];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        n_send_part2 += _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
      }
    }

    /* Allocate */
    gnum1_com_from_entity1_entity2_n      [i_part] = malloc(     n_gnum1_come_from * sizeof(int        ));
    gnum1_com_from_entity1_entity2        [i_part] = malloc(     n_send_part2      * sizeof(PDM_g_num_t));
    gnum1_com_from_entity1_entity2_triplet[i_part] = malloc( 3 * n_send_part2      * sizeof(int        ));

    int         *_gnum1_com_from_triplet_n               = gnum1_com_from_entity1_entity2_n      [i_part];
    PDM_g_num_t *_gnum1_com_from_entity1_entity2         = gnum1_com_from_entity1_entity2        [i_part];
    int         *_gnum1_com_from_entity1_entity2_triplet = gnum1_com_from_entity1_entity2_triplet[i_part];

    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        _gnum1_com_from_triplet_n[k] = _pentity1_entity2_idx[i_entity1+1] - _pentity1_entity2_idx[i_entity1];
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          int i_entity2 = PDM_ABS(_pentity1_entity2[idx_entity2])-1;
          _gnum1_com_from_entity1_entity2        [  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2  ] = i_rank;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+1] = i_part;
          _gnum1_com_from_entity1_entity2_triplet[3*n_send_part2+2] = i_entity2;
          n_send_part2++;
        }
      }
    }
  }


  int exch_request = -1;
  int         **pextract_entity1_entity2_n    = NULL;
  PDM_g_num_t **pextract_entity1_entity2_gnum = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_gnum,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  int **pextract_entity1_entity2_triplet = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 3 * sizeof(int),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2_triplet,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_triplet,
                                 &exch_request);
  PDM_part_to_part_reverse_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(gnum1_com_from_entity1_entity2_n      [i_part]);
    free(gnum1_com_from_entity1_entity2        [i_part]);
    free(gnum1_com_from_entity1_entity2_triplet[i_part]);
  }
  free(gnum1_com_from_entity1_entity2_n      );
  free(gnum1_com_from_entity1_entity2        );
  free(gnum1_com_from_entity1_entity2_triplet);

  PDM_part_to_part_free(ptp);

  int **pextract_entity1_entity2_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    pextract_entity1_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity1_entity2_n[i_part], n_part1_to_part2);
  }


  /*
   * Get for all receive connectivity information about sens and opp gnum
   * CAUTION we need to manage SENS and SIGN
   * We concatenate for all partition
   */
  PDM_g_num_t ***query_itrf_gnum    = malloc(n_interface * sizeof(PDM_g_num_t **));
  int          **query_itrf_n       = malloc(n_interface * sizeof(int          *));
  int         ***recv_itrf_opp_n    = malloc(n_interface * sizeof(int         **));
  PDM_g_num_t ***recv_itrf_opp_gnum = malloc(n_interface * sizeof(PDM_g_num_t **));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    query_itrf_gnum   [i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
    recv_itrf_opp_n   [i_interface] = malloc(ln_part_tot * sizeof(int         *));
    recv_itrf_opp_gnum[i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
    query_itrf_n      [i_interface] = malloc(ln_part_tot * sizeof(PDM_g_num_t  ));

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      query_itrf_n      [i_interface][i_part] = 0;
    }
  }


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int         *_pextract_entity1_interface    = pentity1_extented_to_pentity1_interface[i_part];
    int         *_pextract_entity1_entity2_idx  = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum = pextract_entity1_entity2_gnum   [i_part];

    /* Count */
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {
        int i_itrf = PDM_ABS(_pextract_entity1_interface[i_ref])-1;
        query_itrf_n[i_itrf][i_part] += _pextract_entity1_entity2_idx[i_ref+1]-_pextract_entity1_entity2_idx[i_ref];
      }
    }

    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      query_itrf_gnum   [i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(PDM_g_num_t));
      recv_itrf_opp_n   [i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(int        ));
      recv_itrf_opp_gnum[i_interface][i_part] = malloc(query_itrf_n[i_interface][i_part] * sizeof(PDM_g_num_t));
      query_itrf_n[i_interface][i_part] = 0;
    }

    /* Fill */
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {
        int i_itrf = PDM_ABS(_pextract_entity1_interface[i_ref])-1;
        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {
          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int idx_write = query_itrf_n[i_itrf][i_part]++;
          query_itrf_gnum[i_itrf][i_part][idx_write] = entity2_g_num;
        }
      }
    }
  }


  if(1 == 1) {
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
        PDM_log_trace_array_long(query_itrf_gnum[i_interface][i_part], query_itrf_n[i_interface][i_part], "query_itrf_gnum ::");
      }
    }
  }

  PDM_g_num_t ***recv_interface_entity2_opp_gnum   = malloc(n_interface * sizeof(PDM_g_num_t **));
  int          **n_recv_interface_entity2_opp_gnum = malloc(n_interface * sizeof(int          *));
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    n_recv_interface_entity2_opp_gnum[i_interface] = malloc(ln_part_tot * sizeof(int));

    int n_gnum            = PDM_part_to_block_n_elt_block_get(ptb_interface_entity2[i_interface]);
    PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get (ptb_interface_entity2[i_interface]);

    PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(blk_gnum,
                                                                          n_gnum,
                                                  (const PDM_g_num_t **)  query_itrf_gnum[i_interface],
                                                                          query_itrf_n   [i_interface],
                                                                          ln_part_tot,
                                                                          comm);
    int         **recv_stride = NULL;
    int *send_stride = PDM_array_const_int(n_gnum, 1);
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           send_stride,
                           interface_entity2_opp_gnum[i_interface],
                           &recv_stride,
           (void ***)      &recv_interface_entity2_opp_gnum[i_interface]);
    free(send_stride);

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(recv_stride[i_part], query_itrf_n[i_interface][i_part], "recv_stride ::");

      int n_recv_tot = 0;
      for(int i = 0; i < query_itrf_n[i_interface][i_part]; ++i) {
        n_recv_tot += recv_stride[i_part][i];
      }
      n_recv_interface_entity2_opp_gnum[i_interface][i_part] = n_recv_tot;
      free(recv_stride[i_part]);
      PDM_log_trace_array_long(recv_interface_entity2_opp_gnum[i_interface][i_part], n_recv_tot, "recv_interface_entity2_opp_gnum ::");
    }
    free(recv_stride);

    PDM_block_to_part_free(btp);
  }

  for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
    free(interface_entity2_opp_gnum[i_itrf]);
    PDM_part_to_block_free(ptb_interface_entity2[i_itrf]);
  }
  free(interface_entity2_opp_gnum);
  free(ptb_interface_entity2);

  int          *pn_entity2_ext_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_ext_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_recv_tot = 0;
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      n_recv_tot += n_recv_interface_entity2_opp_gnum[i_interface][i_part];
    }
    pentity2_ext_opp_gnum_and_itrf[i_part] = malloc(2 * n_recv_tot * sizeof(PDM_g_num_t));

    n_recv_tot = 0;
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      for(int i = 0; i < n_recv_interface_entity2_opp_gnum[i_interface][i_part]; ++i) {
        PDM_g_num_t gnum = PDM_ABS (recv_interface_entity2_opp_gnum[i_interface][i_part][i]);
        int         sgn  = PDM_SIGN(recv_interface_entity2_opp_gnum[i_interface][i_part][i]);
        pentity2_ext_opp_gnum_and_itrf[i_part][2*n_recv_tot  ] = gnum;
        pentity2_ext_opp_gnum_and_itrf[i_part][2*n_recv_tot+1] = -sgn * (i_interface+1); // Avec le sgn ?
        n_recv_tot++;
      }
    }

    int *order = malloc(n_recv_tot * sizeof(int));
    pn_entity2_ext_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_recv_tot,
                                                                                   2,
                                                                                   pentity2_ext_opp_gnum_and_itrf[i_part],
                                                                                   order);
    pentity2_ext_opp_gnum_and_itrf[i_part] = realloc(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part] * sizeof(PDM_g_num_t));


    free(order);


    PDM_log_trace_array_long(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part], "pentity2_opp_gnum_and_itrf ::");

  }


  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      free(query_itrf_gnum   [i_interface][i_part]);
      free(recv_itrf_opp_n   [i_interface][i_part]);
      free(recv_itrf_opp_gnum[i_interface][i_part]);
      free(recv_interface_entity2_opp_gnum[i_interface][i_part]);
    }
    free(query_itrf_gnum   [i_interface]);
    free(recv_itrf_opp_n   [i_interface]);
    free(recv_itrf_opp_gnum[i_interface]);
    free(query_itrf_n      [i_interface]);
    free(recv_interface_entity2_opp_gnum  [i_interface]);
    free(n_recv_interface_entity2_opp_gnum[i_interface]);
  }
  free(recv_interface_entity2_opp_gnum  );
  free(n_recv_interface_entity2_opp_gnum);

  free(query_itrf_gnum   );
  free(query_itrf_n      );
  free(recv_itrf_opp_n   );
  free(recv_itrf_opp_gnum);




  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     ln_part_tot,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2);


  /*
   * Prepare for unify
   */
  int         **entity2_order            = malloc( ln_part_tot * sizeof(int         *));
  PDM_g_num_t **pentity2_ln_to_gn_sorted = malloc( ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    /*
     * Sort current part entity2_ln_to_gn
     */
    entity2_order           [i_part] = malloc( pn_entity2[i_part] * sizeof(int        ));
    pentity2_ln_to_gn_sorted[i_part] = malloc( pn_entity2[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2[i_part]; ++i_entity2) {
      pentity2_ln_to_gn_sorted[i_part][i_entity2] = pentity2_ln_to_gn[i_part][i_entity2];
      entity2_order           [i_part][i_entity2] = i_entity2;
    }
    PDM_sort_long(pentity2_ln_to_gn_sorted[i_part], entity2_order[i_part], pn_entity2[i_part]);
  }


  /*
   * We have two kind of extented :
   *   - From partition
   *   - From interfaces (request new gnum génération ...)
   */
  int          *pn_entity2_extented_by_interface         = malloc(ln_part_tot * sizeof(int          ));
  int          *pn_entity2_extented_by_partition         = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_interface  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_partition  = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_triplet_by_interface   = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_triplet_by_partition   = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface_by_interface = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface_by_partition = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* _pentity2_opp_gnum_and_itrf = pentity2_opp_gnum_and_itrf[i_part];

    int         *_pextract_entity1_interface = pentity1_extented_to_pentity1_interface[i_part];
    // PDM_g_num_t *_pextract_entity1_gnum      = pextract_entity1_gnum     [i_part];
    // int         *_pextract_entity1_triplet   = pextract_entity1_triplet  [i_part];

    int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
    int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

    PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];

    pn_entity2_extented_by_interface[i_part] = 0;
    pn_entity2_extented_by_partition[i_part] = 0;

    /* Count */
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, pn_entity2_opp_gnum_and_itrf[i_part]);
          // log_trace("gnum_to_find = %i/%i --> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);
          if(pos == -1) {
            pn_entity2_extented_by_interface[i_part]++;
          }
          // printf("Search (%i/%i) -> pos = %i \n", entity2_g_num, _pextract_entity1_interface[i_ref], pos);
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int == -1) { // Not in current partition
            pn_entity2_extented_by_partition[i_part]++;
          }
        }
      }
    }

    /* Fill */
    pentity2_extented_ln_to_gn_by_interface [i_part] = malloc(2 * pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    pentity2_extented_ln_to_gn_by_partition [i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    pentity2_extented_triplet_by_interface  [i_part] = malloc(3 * pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    pentity2_extented_triplet_by_partition  [i_part] = malloc(3 * pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
    pentity2_extented_interface_by_interface[i_part] = malloc(    pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    pentity2_extented_interface_by_partition[i_part] = malloc(    pn_entity2_extented_by_partition[i_part] * sizeof(int        ));

    pn_entity2_extented_by_interface[i_part] = 0;
    pn_entity2_extented_by_partition[i_part] = 0;
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {
      if(_pextract_entity1_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, pn_entity2_opp_gnum_and_itrf[i_part]);
          if(pos == -1) {
            int idx_write = pn_entity2_extented_by_interface[i_part]++;
            pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write  ] = entity2_g_num;
            pentity2_extented_ln_to_gn_by_interface[i_part][2*idx_write+1] = _pextract_entity1_interface[i_ref];

            pentity2_extented_triplet_by_interface  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
            pentity2_extented_triplet_by_interface  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
            pentity2_extented_triplet_by_interface  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
            pentity2_extented_interface_by_interface[i_part][  idx_write  ] = _pextract_entity1_interface[i_ref];
          }
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int == -1) { // Not in current partition
            int idx_write = pn_entity2_extented_by_partition[i_part]++;
            pentity2_extented_ln_to_gn_by_partition [i_part][idx_write] = entity2_g_num;
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
            pentity2_extented_triplet_by_partition  [i_part][3*idx_write+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
            pentity2_extented_interface_by_partition[i_part][  idx_write  ] = 0; // Because interior
          }
        }
      }
    }

    PDM_log_trace_array_long(pentity2_extented_ln_to_gn_by_interface[i_part],
                             2 * pn_entity2_extented_by_interface[i_part],
                             "pentity2_extented_ln_to_gn_by_interface ::");

    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              pn_entity2_extented_by_interface[i_part],
                              pentity2_extented_ln_to_gn_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_entity2);

  /*
   * Unify entity1_entity2
   */
  int **pextract_entity1_entity2 = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pentity2_opp_gnum_and_itrf  = pentity2_opp_gnum_and_itrf             [i_part];
    int         *_pextract_entity1_interface  = pentity1_extented_to_pentity1_interface[i_part];
    // PDM_g_num_t *_pextract_entity1_gnum      = pextract_entity1_gnum     [i_part];

    int         *_pextract_entity1_entity2_idx     = pextract_entity1_entity2_idx    [i_part];
    PDM_g_num_t *_pextract_entity1_entity2_gnum    = pextract_entity1_entity2_gnum   [i_part];
    int         *_pextract_entity1_entity2_triplet = pextract_entity1_entity2_triplet[i_part];

    int           n_entity2_opp_position     = pn_entity2_opp_gnum_and_itrf[i_part];

    PDM_g_num_t *_pentity2_ln_to_gn_sorted   = pentity2_ln_to_gn_sorted[i_part];
    int         *_entity2_order              = entity2_order           [i_part];

    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;


    pextract_entity1_entity2[i_part] = malloc( _pextract_entity1_entity2_idx[n_part1_to_part2] * sizeof(int));
    int *_pextract_entity1_entity2  = pextract_entity1_entity2 [i_part];

    PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    // Erase and realloc :
    pentity2_extented_ln_to_gn_by_interface[i_part] = realloc(pentity2_extented_ln_to_gn_by_interface[i_part], pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      pentity2_extented_ln_to_gn_by_interface[i_part][i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2] + shift_by_domain_entity2[n_domain];
    }

    // PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn ::");

    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_itrf_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_ln_to_gn_sorted[i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2];
    }
    // PDM_sort_long(extented_from_itrf_entity2_ln_to_gn_sorted, extented_from_itrf_entity2_order, pn_entity2_extented_by_interface[i_part]);
    pn_entity2_extented_by_interface[i_part] = PDM_inplace_unique_long(extented_from_itrf_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_interface[i_part]-1);
    extented_from_itrf_entity2_ln_to_gn_sorted = realloc(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part] * sizeof(PDM_g_num_t));
    int         *extented_from_itrf_entity2_order = malloc( pn_entity2_extented_by_interface[i_part] * sizeof(int        ));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_order[i_entity2] = -1;
    }

    printf("pn_entity2_extented_by_interface[i_part] = %i \n", pn_entity2_extented_by_interface[i_part]);


    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_part_entity2_ln_to_gn_sorted = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      extented_from_part_entity2_ln_to_gn_sorted[i_entity2] = pentity2_extented_ln_to_gn_by_partition[i_part][i_entity2];
      // extented_from_part_entity2_order          [i_entity2] = i_entity2;
    }
    // PDM_sort_long(extented_from_part_entity2_ln_to_gn_sorted, extented_from_part_entity2_order, pn_entity2_extented_by_partition[i_part]);

    pn_entity2_extented_by_partition[i_part] = PDM_inplace_unique_long(extented_from_part_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_partition[i_part]-1);
    int *extented_from_part_entity2_order = malloc( pn_entity2_extented_by_partition[i_part] * sizeof(int        ));
    extented_from_part_entity2_ln_to_gn_sorted = realloc(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part] * sizeof(PDM_g_num_t));
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      extented_from_part_entity2_order[i_entity2] = -1;
    }

    printf("pn_entity2_extented_by_partition[i_part] = %i \n", pn_entity2_extented_by_partition[i_part]);


    PDM_log_trace_array_long(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part], "extented_from_part_entity2_ln_to_gn_sorted ::");
    PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn_sorted ::");



    /* */
    int pn_entity2_extented_by_interface2 = 0; // To read in extented_from_itrf_entity2_ln_to_gn
    for(int i_ref = 0; i_ref < n_part1_to_part2; ++i_ref) {

      log_trace("_pextract_entity1_interface[i_ref] = %i \n", _pextract_entity1_interface[i_ref] );

      /*
       * First case :
       *   - entity1 is move by interface
       */
      if(_pextract_entity1_interface[i_ref] != 0) {
        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          PDM_g_num_t gnum_to_find[2] = {entity2_g_num, _pextract_entity1_interface[i_ref]};

          int pos = PDM_order_binary_search_long(gnum_to_find, _pentity2_opp_gnum_and_itrf, 2, n_entity2_opp_position);
          // log_trace("\t Search %i/%i in _pentity2_opp_gnum_and_itrf -> pos = %i \n", gnum_to_find[0], gnum_to_find[1], pos);

          if(pos == -1) {
            /*
             * Subcase :
             *   - entity2 is not in table of interface : it's a new entity2
             */
            PDM_g_num_t entity2_extented_g_num = extented_from_itrf_entity2_ln_to_gn[pn_entity2_extented_by_interface2++];
            int pos_new = PDM_binary_search_long(entity2_extented_g_num,
                                                 extented_from_itrf_entity2_ln_to_gn_sorted,
                                                 pn_entity2_extented_by_interface[i_part]);

            assert(pos_new != -1);
            int shift = pn_entity2[i_part] + pn_entity2_extented_by_partition[i_part];

            if(extented_from_itrf_entity2_order[pos_new] == -1) {
              extented_from_itrf_entity2_order[pos_new] = 1;

              pentity2_extented_ln_to_gn_by_interface [i_part][  pos_new  ] = entity2_extented_g_num + shift_by_domain_entity2[n_domain];

              pentity2_extented_triplet_by_interface  [i_part][3*pos_new  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
              pentity2_extented_triplet_by_interface  [i_part][3*pos_new+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
              pentity2_extented_triplet_by_interface  [i_part][3*pos_new+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
              pentity2_extented_interface_by_interface[i_part][  pos_new  ] = _pextract_entity1_interface[i_ref];

            }

            // _pextract_entity1_entity2[idx_entity1] = ( shift + extented_from_itrf_entity2_order[pos_new] + 1); // ATTENTION SIGN
            _pextract_entity1_entity2[idx_entity1] = ( shift + pos_new + 1); // ATTENTION SIGN

            log_trace("\t Translate cas 1 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          } else {
            /*
             * Subcase :
             *   - entity2 is in table of interface
             */
            int pos2      = entity2_opp_position[i_part][pos ];
            _pextract_entity1_entity2[idx_entity1] = ( pos2 + 1);

            // log_trace("\t Translate cas 2 : %i  ---> (idx_entity1 = %i) --> %i (pos=%i/pos2=%i) \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1], pos, pos2);


          }
        }
      } else { // Second case : Interior entity1

        for(int idx_entity1 = _pextract_entity1_entity2_idx[i_ref]; idx_entity1 < _pextract_entity1_entity2_idx[i_ref+1]; ++idx_entity1) {

          PDM_g_num_t entity2_g_num = _pextract_entity1_entity2_gnum[idx_entity1];
          int pos_int = PDM_binary_search_long(entity2_g_num, _pentity2_ln_to_gn_sorted, pn_entity2[i_part]);

          if(pos_int != -1) {                   // entity2 is already in current partition
            assert(pos_int != -1);
            int orig_pos = _entity2_order[pos_int];
            _pextract_entity1_entity2[idx_entity1] = ( orig_pos + 1);

            log_trace("\t Translate cas 3 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          } else {                              // Not in current partition
            int pos_ext = PDM_binary_search_long(entity2_g_num,
                                                 extented_from_part_entity2_ln_to_gn_sorted,
                                                 pn_entity2_extented_by_partition[i_part]);
            assert(pos_ext != -1);

            if(extented_from_part_entity2_order[pos_ext] == -1) {
              extented_from_part_entity2_order[pos_ext] = 1;
              pentity2_extented_ln_to_gn_by_partition [i_part][  pos_ext  ] = entity2_g_num;

              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext  ] = _pextract_entity1_entity2_triplet[3*idx_entity1  ];
              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+1] = _pextract_entity1_entity2_triplet[3*idx_entity1+1];
              pentity2_extented_triplet_by_partition  [i_part][3*pos_ext+2] = _pextract_entity1_entity2_triplet[3*idx_entity1+2];
              pentity2_extented_interface_by_partition[i_part][  pos_ext  ] = 0;

            }

            // _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + extented_from_part_entity2_order[pos_ext] + 1);
            _pextract_entity1_entity2[idx_entity1] = ( pn_entity2[i_part] + pos_ext + 1);

            log_trace("\t Translate cas 4 : %i  ---> (idx_entity1 = %i) --> %i \n", entity2_g_num, idx_entity1, _pextract_entity1_entity2[idx_entity1]);

          }
        }
      }
    }

    // assert(pn_entity2_extented_by_interface2 == pn_entity2_extented_by_interface[i_part]);

    if(1 == 1) {
      PDM_log_trace_connectivity_int(_pextract_entity1_entity2_idx,
                                     _pextract_entity1_entity2,
                                     n_part1_to_part2,
                                     "pextract_entity1_entity2 ::");
    }

    free(extented_from_itrf_entity2_order          );
    free(extented_from_itrf_entity2_ln_to_gn_sorted);
    free(extented_from_part_entity2_order          );
    free(extented_from_part_entity2_ln_to_gn_sorted);

  }


  PDM_gnum_free(gen_gnum_entity2);

  /*
   * Reconstruction
   */
  int          *pn_entity2_extented              = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pentity2_extented_ln_to_gn       = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  int         **pentity2_extented_triplet        = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_interface      = malloc(ln_part_tot * sizeof(int         *));
  int         **pentity2_extented_to_entity2_idx = malloc(ln_part_tot * sizeof(int         *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    pn_entity2_extented[i_part] = pn_entity2_extented_by_partition[i_part] + pn_entity2_extented_by_interface[i_part];
    pentity2_extented_ln_to_gn      [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(PDM_g_num_t));
    pentity2_extented_triplet       [i_part] = malloc(3 * pn_entity2_extented[i_part]      * sizeof(int        ));
    pentity2_extented_interface     [i_part] = malloc(    pn_entity2_extented[i_part]      * sizeof(int        ));
    pentity2_extented_to_entity2_idx[i_part] = malloc( (  pn_entity2_extented[i_part] + 1) * sizeof(int        ));

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      pentity2_extented_ln_to_gn [i_part][  i_entity2  ] = pentity2_extented_ln_to_gn_by_partition [i_part][  i_entity2  ];
      pentity2_extented_triplet  [i_part][3*i_entity2  ] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2  ];
      pentity2_extented_triplet  [i_part][3*i_entity2+1] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+1];
      pentity2_extented_triplet  [i_part][3*i_entity2+2] = pentity2_extented_triplet_by_partition  [i_part][3*i_entity2+2];
      pentity2_extented_interface[i_part][  i_entity2  ] = pentity2_extented_interface_by_partition[i_part][  i_entity2  ];
    }

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      int idx_write = pn_entity2_extented_by_partition[i_part]+i_entity2;
      pentity2_extented_ln_to_gn [i_part][  idx_write  ] = pentity2_extented_ln_to_gn_by_interface [i_part][  i_entity2  ];
      pentity2_extented_triplet  [i_part][3*idx_write  ] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2  ];
      pentity2_extented_triplet  [i_part][3*idx_write+1] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+1];
      pentity2_extented_triplet  [i_part][3*idx_write+2] = pentity2_extented_triplet_by_interface  [i_part][3*i_entity2+2];
      pentity2_extented_interface[i_part][  idx_write  ] = pentity2_extented_interface_by_interface[i_part][  i_entity2  ];
    }

    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented[i_part]+1; ++i_entity2) {
      pentity2_extented_to_entity2_idx[i_part][i_entity2] = 3 * i_entity2;
    }

    if(1 == 1) {
      PDM_log_trace_array_int (pentity2_extented_to_entity2_idx[i_part],     pn_entity2_extented[i_part]+1, "pentity2_extented_to_entity2_idx ::" );
      PDM_log_trace_array_int (pentity2_extented_triplet       [i_part], 3 * pn_entity2_extented[i_part]  , "pentity2_extented_triplet    ::" );
      PDM_log_trace_array_int (pentity2_extented_interface     [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_interface  ::" );
      PDM_log_trace_array_long(pentity2_extented_ln_to_gn      [i_part],     pn_entity2_extented[i_part]  , "pentity2_extented_ln_to_gn   ::" );
    }

  }

  *pn_entity2_extented_out                     = pn_entity2_extented;
  *pentity2_extented_ln_to_gn_out              = pentity2_extented_ln_to_gn;
  *pextented_entity1_entity2_idx_out           = pextract_entity1_entity2_idx;
  *pextented_entity1_entity2_out               = pextract_entity1_entity2;
  *pentity2_extented_to_pentity2_idx_out       = pentity2_extented_to_entity2_idx;
  *pentity2_extented_to_pentity2_triplet_out   = pentity2_extented_triplet;
  *pentity2_extented_to_pentity2_interface_out = pentity2_extented_interface;



  // for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
  //   free(pextract_entity1_entity2[i_part]);
  // }
  // free(pextract_entity1_entity2);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_entity1_entity2_n      [i_part]);
    free(pextract_entity1_entity2_gnum   [i_part]);
    free(pextract_entity1_entity2_triplet[i_part]);
    // free(pextract_entity1_entity2_idx [i_part]);
  }
  free(pextract_entity1_entity2_n      );
  free(pextract_entity1_entity2_gnum   );
  free(pextract_entity1_entity2_triplet);
  // free(pextract_entity1_entity2_idx );


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_extented_ln_to_gn_by_interface [i_part]);
    free(pentity2_extented_ln_to_gn_by_partition [i_part]);
    free(pentity2_extented_triplet_by_interface  [i_part]);
    free(pentity2_extented_triplet_by_partition  [i_part]);
    free(pentity2_extented_interface_by_interface[i_part]);
    free(pentity2_extented_interface_by_partition[i_part]);
  }
  free(pentity2_extented_ln_to_gn_by_interface );
  free(pentity2_extented_ln_to_gn_by_partition );
  free(pentity2_extented_triplet_by_interface  );
  free(pentity2_extented_triplet_by_partition  );
  free(pentity2_extented_interface_by_interface);
  free(pentity2_extented_interface_by_partition);

  free(pn_entity2_extented_by_interface );
  free(pn_entity2_extented_by_partition );

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(entity2_order           [i_part]);
    free(pentity2_ln_to_gn_sorted[i_part]);
  }
  free(entity2_order           );
  free(pentity2_ln_to_gn_sorted);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_opp_gnum_and_itrf[i_part]);
    free(pentity2_ext_opp_gnum_and_itrf[i_part]);
    free(entity2_opp_position      [i_part]);
  }
  free(entity2_opp_position        );
  free(pentity2_opp_gnum_and_itrf  );
  free(pentity2_ext_opp_gnum_and_itrf  );
  free(pn_entity2_opp_gnum_and_itrf);
  free(pn_entity2_ext_opp_gnum_and_itrf);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pentity2_num             [i_part]);
    free(pentity2_opp_location_idx[i_part]);
    free(pentity2_opp_location    [i_part]);
    free(pentity2_opp_interface   [i_part]);
    free(pentity2_opp_sens        [i_part]);
    free(pentity2_opp_gnum        [i_part]);
  }
  free(pn_entity2_num            );
  free(pentity2_num              );
  free(pentity2_opp_location_idx );
  free(pentity2_opp_location     );
  free(pentity2_opp_interface_idx);
  free(pentity2_opp_interface    );
  free(pentity2_opp_sens         );
  free(pentity2_opp_gnum         );

  free(pn_entity1          );
  free(pentity1_ln_to_gn   );
  free(pentity1_entity2_idx);
  free(pentity1_entity2    );
  free(pn_entity2          );
  free(pentity2_ln_to_gn   );

}



#ifdef __cplusplus
}
#endif /* __cplusplus */

