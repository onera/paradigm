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
_tell_me_more_data_base
(
  int          n_entity,
  PDM_g_num_t *pentity_gnum,
  int         *pentity_ancstr_strd,
  PDM_g_num_t *pentity_ancstr,
  int         *pentity_path_itrf_strd,
  int         *pentity_path_itrf,
  int         *pentity_data_strd,
  PDM_g_num_t *pentity_data,
  PDM_g_num_t  wanted_gnum
)
{
  int i_read_ancstr    = 0;
  int i_read_path_itrf = 0;
  int i_read_data      = 0;
  for (int i_entity=0; i_entity<n_entity; ++i_entity) {
    if (pentity_gnum[i_entity]==wanted_gnum) {
      log_trace("\ndata_base:: entity gnum = "PDM_FMT_G_NUM" :\n", pentity_gnum[i_entity]);

      // > print ancestor
      int l_i_read_ancstr    = i_read_ancstr   ;
      log_trace("\t ancestor = ");
      for (int i_ancstr=0; i_ancstr<pentity_ancstr_strd[i_entity]; ++i_ancstr) {
        log_trace("%d ", pentity_ancstr[l_i_read_ancstr]);
        l_i_read_ancstr++;
      }
      log_trace("\n");

      // > print path itrf
      int l_i_read_path_itrf = i_read_path_itrf;
      log_trace("\t path itrf = ");
      for (int i_path_itrf=0; i_path_itrf<pentity_path_itrf_strd[i_entity]; ++i_path_itrf) {
        log_trace("%d ", pentity_path_itrf[l_i_read_path_itrf]);
        l_i_read_path_itrf++;
      }
      log_trace("\n");

      // > print data
      int l_i_read_data = i_read_data;
      log_trace("\t is linked with :\n");
      for (int i_data=0; i_data<pentity_data_strd[i_entity]; ++i_data) {
        log_trace("gnum = "PDM_FMT_G_NUM" via itrf = "PDM_FMT_G_NUM" \n",
                        pentity_data[2*l_i_read_data  ],
                        pentity_data[2*l_i_read_data+1]);
        l_i_read_data++;
      }
      log_trace("\n");
    }
    else {
      i_read_data      += pentity_data_strd     [i_entity];
      i_read_ancstr    += pentity_ancstr_strd   [i_entity];
      i_read_path_itrf += pentity_path_itrf_strd[i_entity];
    }
  }
}



static
void
_tell_me_more_valid_entities
(
  int          pn_entity1,
  PDM_g_num_t *pentity1_gnum,
  int         *pentity1_to_pentity1_idx,
  int         *pentity1_entity2_idx,
  PDM_g_num_t *pentity1_entity2_gnum,
  int         *pentity1_entity2_sign,
  int         *pentity1_entity2_kind,
  int         *pentity1_entity2_lnum,
  int         *pentity1_entity2_itrf,
  int         *pentity1_entity2_sens,
  PDM_g_num_t  wanted_gnum
)
{
  for (int i_entity1=0; i_entity1<pn_entity1; ++i_entity1) {
    if (pentity1_gnum[i_entity1]==wanted_gnum) {
      log_trace("\nvalid_entities:: entity1 gnum = "PDM_FMT_G_NUM" connectivity :\n", pentity1_gnum[i_entity1]);
      int i_beg1 = pentity1_to_pentity1_idx[i_entity1  ]/3;
      int i_end1 = pentity1_to_pentity1_idx[i_entity1+1]/3;
      int n_read1 = i_end1-i_beg1;
      for (int i_read1=i_beg1; i_read1<i_end1; ++i_read1) {
        log_trace("\t candidate %d/%d \n", i_read1-i_beg1+1, n_read1);
        int i_beg = pentity1_entity2_idx[i_read1  ];
        int i_end = pentity1_entity2_idx[i_read1+1];
        for (int i_entity2=i_beg; i_entity2<i_end; ++i_entity2) {
          log_trace("\t\t entity2 gnum = "PDM_FMT_G_NUM" ", pentity1_entity2_gnum[i_entity2]);
          if (pentity1_entity2_sign!=NULL) {
            log_trace(" with sign = %d", pentity1_entity2_sign[i_entity2]);
          }
          log_trace(" became: \n");
          log_trace("\t\t\t   kind = %d\n", pentity1_entity2_kind[i_entity2]);
          log_trace("\t\t\t   lnum = %d\n", pentity1_entity2_lnum[i_entity2]);
          log_trace("\t\t\t   itrf = %d\n", pentity1_entity2_itrf[i_entity2]);
          if (pentity1_entity2_sens!=NULL) {
            log_trace("\t\t\t   sens = %d\n", pentity1_entity2_sens[i_entity2]);
          }
        }
      }
    }
  }
}

static
void
_tell_me_more_received
(
  int          pn_entity1,
  PDM_g_num_t *pentity1_gnum,
  int         *pentity1_to_pentity1_idx,
  int         *pentity1_entity2_idx,
  PDM_g_num_t *pentity1_entity2_gnum,
  int         *pentity1_entity2_trplt,
  int         *pentity1_entity2_sign,
  PDM_g_num_t  wanted_gnum
)
{
  for (int i_entity1=0; i_entity1<pn_entity1; ++i_entity1) {
    if (pentity1_gnum[i_entity1]==wanted_gnum) {
      log_trace("\nreceived:: entity1 gnum = "PDM_FMT_G_NUM" connectivity :\n", pentity1_gnum[i_entity1]);
      int i_beg1 = pentity1_to_pentity1_idx[i_entity1  ]/3;
      int i_end1 = pentity1_to_pentity1_idx[i_entity1+1]/3;
      int n_read1 = i_end1-i_beg1;
      for (int i_read1=i_beg1; i_read1<i_end1; ++i_read1) {
        log_trace("\t candidate %d/%d \n", i_read1-i_beg1+1, n_read1);
        int i_beg2 = pentity1_entity2_idx[i_read1  ];
        int i_end2 = pentity1_entity2_idx[i_read1+1];
        for (int i_entity2=i_beg2; i_entity2<i_end2; ++i_entity2) {
          log_trace("\t\t entity2 gnum = "PDM_FMT_G_NUM" ", pentity1_entity2_gnum[i_entity2]);
          log_trace("(%d, %d, %d) ", pentity1_entity2_trplt[3*i_entity2  ]
                                   , pentity1_entity2_trplt[3*i_entity2+1]
                                   , pentity1_entity2_trplt[3*i_entity2+2]);
          if (pentity1_entity2_sign!=NULL) {
            log_trace(" with sign = %d", pentity1_entity2_sign[i_entity2]);
          }
          log_trace("\n");
        }
      }
    }
  }
}

static
void
_tell_me_more_connectivity
(
  int          pn_entity1,
  PDM_g_num_t *pentity1_gnum,
  int         *pentity1_entity2_idx,
  int         *pentity1_entity2,
  PDM_g_num_t *pentity2_gnum,
  PDM_g_num_t  wanted_gnum
)
{
  for (int i_entity1=0; i_entity1<pn_entity1; ++i_entity1) {
    if (pentity1_gnum[i_entity1]==wanted_gnum) {
      log_trace("\nconnectivity:: entity1 gnum = "PDM_FMT_G_NUM" connectivity :\n", pentity1_gnum[i_entity1]);
      int i_beg = pentity1_entity2_idx[i_entity1  ];
      int i_end = pentity1_entity2_idx[i_entity1+1];
      log_trace("\t entity2 lnum = ");
      for (int i_entity2=i_beg; i_entity2<i_end; ++i_entity2) {
        log_trace("%d ", pentity1_entity2[i_entity2]);
      }
      log_trace("\n");
      log_trace("\t entity2 gnum = ");
      for (int i_entity2=i_beg; i_entity2<i_end; ++i_entity2) {
        log_trace("%d ", PDM_SIGN(pentity1_entity2[i_entity2])*pentity2_gnum[PDM_ABS(pentity1_entity2[i_entity2])-1]);
      }
      log_trace("\n");
    }
  }
}

static
void
_tell_me_more_link
(
  int          pn_entity,
  PDM_g_num_t *pentity_gnum,
  int         *pentity_to_pentity_idx,
  int         *pentity_to_pentity_trplt,
  int         *pentity_to_pentity_itrf,
  int         *pentity_to_pentity_sens,
  PDM_g_num_t  wanted_gnum
)
{
  for (int i_entity=0; i_entity<pn_entity; ++i_entity) {
    if (pentity_gnum[i_entity]==wanted_gnum) {
      log_trace("\nlink:: gnum = "PDM_FMT_G_NUM" connected to \n", pentity_gnum[i_entity]);
      int i_beg = pentity_to_pentity_idx[i_entity  ]/3;
      int i_end = pentity_to_pentity_idx[i_entity+1]/3;
      for (int i_trplt=i_beg; i_trplt<i_end; ++i_trplt) {
        log_trace("\t triplet = (%d,%d,%d)\n", pentity_to_pentity_trplt[3*i_trplt  ],
                                               pentity_to_pentity_trplt[3*i_trplt+1],
                                               pentity_to_pentity_trplt[3*i_trplt+2]);
        log_trace("\t itrf    = %d\n", pentity_to_pentity_itrf[i_trplt]);
        if (pentity_to_pentity_sens!=NULL) {
        log_trace("\t sens    = %d\n", pentity_to_pentity_sens[i_trplt]);
        }
      }
    }
  }
}


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

  int *concat_bentity1_entity2_n               = NULL;
  int *concat_bentity1_entity2_interface_tot_n = NULL;
  PDM_malloc(concat_bentity1_entity2_n              , pn_entity1, int);
  PDM_malloc(concat_bentity1_entity2_interface_tot_n, pn_entity1, int);

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
  PDM_g_num_t *concat_bentity1_entity2_gnum        = NULL;
  int         *concat_bentity1_entity2_triplet     = NULL;
  int         *concat_bentity1_entity2_interface_n = NULL;
  int         *concat_bentity1_entity2_interface   = NULL;
  PDM_malloc(concat_bentity1_entity2_gnum       ,     n_concat_entity1_entity2     , PDM_g_num_t);
  PDM_malloc(concat_bentity1_entity2_triplet    , 3 * n_concat_entity1_entity2     , int        );
  PDM_malloc(concat_bentity1_entity2_interface_n,     n_concat_entity1_entity2     , int        );
  PDM_malloc(concat_bentity1_entity2_interface  ,     n_concat_entity1_entity2_itrf, int        );

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

    int *pupdate_bentity1_entity2_interface = NULL;
    PDM_malloc(pupdate_bentity1_entity2_interface, n_interface_tot, int);

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

    PDM_free(pnext_bentity1_entity2_interface[i_part]);
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
  int *update_bentity1_entity2_n               = NULL;
  int *update_bentity1_entity2_interface_tot_n = NULL;
  PDM_malloc(update_bentity1_entity2_n              , pn_entity1, int);
  PDM_malloc(update_bentity1_entity2_interface_tot_n, pn_entity1, int);
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

  int* _unique_neighbor_entity_idx = NULL;
  int* _unique_neighbor_entity_n   = NULL;
  int* _unique_neighbor_entity     = NULL;
  int* order                       = NULL;
  PDM_malloc(_unique_neighbor_entity_idx, n_entity + 1, int);
  PDM_malloc(_unique_neighbor_entity_n  , n_entity    , int);
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

  PDM_realloc(_unique_neighbor_entity, _unique_neighbor_entity, 4 * neighbor_entity_idx[n_entity], int);

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  PDM_free(order);
}


static inline
int
_is_same_quintuplet
(
int iproc1, int ipart1, int ielt1, int iinterf1, int isens1,
int iproc2, int ipart2, int ielt2, int iinterf2, int isens2
)
{
  if(iproc1 == iproc2){
    if(ipart1 == ipart2){
      if(ielt1 == ielt2){
        if(iinterf1 == iinterf2){
          if(isens1 == isens2){
            return 1;
          }
        }
      }
    }
  }
  return 0;
}


static
void
_unique_quintuplet
(
  int   n_entity,
  int  *neighbor_entity_idx,
  int  *neighbor_entity,
  int **unique_neighbor_entity_idx,
  int **unique_neighbor_entity_n,
  int **unique_neighbor_entity
)
{

  int* _unique_neighbor_entity_idx = NULL;
  int* _unique_neighbor_entity_n   = NULL;
  int* _unique_neighbor_entity     = NULL;
  int* order                       = NULL;
  PDM_malloc(_unique_neighbor_entity_idx, n_entity + 1, int);
  PDM_malloc(_unique_neighbor_entity_n  , n_entity    , int);
  PDM_malloc(_unique_neighbor_entity    , 5 * neighbor_entity_idx[n_entity], int);
  PDM_malloc(order                      ,     neighbor_entity_idx[n_entity], int); // Suralloc

  _unique_neighbor_entity_idx[0] = 0;
  for(int i_entity = 0; i_entity < n_entity; ++i_entity) {

    int beg       = neighbor_entity_idx[i_entity];
    int n_connect = neighbor_entity_idx[i_entity+1] - beg;

    PDM_order_lnum_s(&neighbor_entity[5*beg], 5, order, n_connect);

    _unique_neighbor_entity_n  [i_entity  ] = 0;
    _unique_neighbor_entity_idx[i_entity+1] = _unique_neighbor_entity_idx[i_entity];

    int last_proc  = -1;
    int last_part  = -1;
    int last_elmt  = -1;
    int last_inte  = -40;
    int last_sens  = -40;
    for(int i = 0; i < n_connect; ++i) {
      int old_order   = order[i];
      int curr_proc   = neighbor_entity[5*(beg+old_order)  ];
      int curr_part   = neighbor_entity[5*(beg+old_order)+1];
      int curr_entity = neighbor_entity[5*(beg+old_order)+2];
      int curr_inte   = neighbor_entity[5*(beg+old_order)+3];
      int curr_sens   = neighbor_entity[5*(beg+old_order)+4];
      int is_same  = _is_same_quintuplet(last_proc, last_part, last_elmt  , last_inte, last_sens,
                                         curr_proc, curr_part, curr_entity, curr_inte, curr_sens);

      if(is_same == 0){ // N'est pas le meme
        // idx_unique++;
        last_proc = curr_proc;
        last_part = curr_part;
        last_elmt = curr_entity;
        last_inte = curr_inte;
        last_sens = curr_sens;

        int beg_write = 5 * _unique_neighbor_entity_idx[i_entity+1];
        // printf("beg_write = %i | curr_proc = %i | curr_part = %i | curr_entity = %i \n", beg_write, curr_proc, curr_part, curr_entity);
        _unique_neighbor_entity[beg_write  ] = curr_proc;
        _unique_neighbor_entity[beg_write+1] = curr_part;
        _unique_neighbor_entity[beg_write+2] = curr_entity;
        _unique_neighbor_entity[beg_write+3] = curr_inte;
        _unique_neighbor_entity[beg_write+4] = curr_sens;

        /* Increment the new counter */
        _unique_neighbor_entity_idx[i_entity+1]++;
        _unique_neighbor_entity_n  [i_entity  ]++;
      }
    }
  }

  PDM_realloc(_unique_neighbor_entity, _unique_neighbor_entity, 5 * neighbor_entity_idx[n_entity], int);

  *unique_neighbor_entity_idx = _unique_neighbor_entity_idx;
  *unique_neighbor_entity_n   = _unique_neighbor_entity_n;
  *unique_neighbor_entity     = _unique_neighbor_entity;
  PDM_free(order);
}




static
void
_find_twin_interface
(
  int            n_part,
  int            d_db_n_entity, // Distributed data_base
  PDM_g_num_t   *d_db_entity_gnum,
  int           *d_db_entity_strd,
  PDM_g_num_t   *d_db_entity_data,
  int           *d_db_entity_sens,
  int           *pn_entity,
  PDM_g_num_t  **pentity_gnum,
  int          **pentity_itrf,
  int          **pentity_twin_itrf_idx,
  int          **pentity_twin_itrf_itrf,
  int          **pentity_twin_itrf_sens,
  PDM_g_num_t  **pentity_twin_itrf_gnum,
  int          **pentity_twin_itrf_lnum,
  int          **p_db_n_entity_out,
  PDM_g_num_t ***p_db_entity_gnum_out,
  PDM_g_num_t ***p_db_entity_data_out,
  int         ***p_db_entity_sens_out,
  PDM_MPI_Comm   comm
)
{
  int debug      = 0;
  int debug_loop = 0;

  if (debug==1) {
    log_trace("\n");
    log_trace("------------------------\n");
    log_trace("BEG _find_twin_interface\n");
  }


  int has_sens = 0;
  if (d_db_entity_sens!=NULL) {
    has_sens = 1;
  }


  int *pn_twin_itrf = NULL;
  PDM_malloc(pn_twin_itrf, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_twin_itrf[i_part] = pentity_twin_itrf_idx[i_part][pn_entity[i_part]];
  }

  /*
   * Get informations from data_base for each element
   *
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(d_db_entity_gnum,
                                                   d_db_n_entity,
                           (const PDM_g_num_t **)  pentity_twin_itrf_gnum,
                                                   pn_twin_itrf,
                                                   n_part,
                                                   comm);

  int         **p_db_entity_strd = NULL;
  PDM_g_num_t **p_db_entity_data = NULL;
  PDM_block_to_part_exch(btp,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         d_db_entity_strd,
                         d_db_entity_data,
                        &p_db_entity_strd,
         (void ***)     &p_db_entity_data);

  int **p_db_entity_strd2= NULL;
  int **p_db_entity_sens = NULL;
  if (has_sens==1) {
    PDM_block_to_part_exch(btp,
                           1 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           d_db_entity_strd,
                           d_db_entity_sens,
                          &p_db_entity_strd2,
           (void ***)     &p_db_entity_sens);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_free(p_db_entity_strd2[i_part]);
    }
    PDM_free(p_db_entity_strd2);
  }
  PDM_block_to_part_free(btp);


  int          *_p_db_n_entity_out    = *p_db_n_entity_out;
  PDM_g_num_t **_p_db_entity_gnum_out = *p_db_entity_gnum_out;
  PDM_g_num_t **_p_db_entity_data_out = *p_db_entity_data_out;
  int         **_p_db_entity_sens_out = NULL;
  if (has_sens==1) {
    _p_db_entity_sens_out = *p_db_entity_sens_out;
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {

    if (debug==1) {
      PDM_log_trace_array_long(_p_db_entity_gnum_out[i_part],   _p_db_n_entity_out[i_part], "_p_db_entity_gnum_out ::");
      PDM_log_trace_array_long(_p_db_entity_data_out[i_part], 2*_p_db_n_entity_out[i_part], "_p_db_entity_data_out ::");
      if (has_sens==1) {
        PDM_log_trace_array_int (_p_db_entity_sens_out[i_part], _p_db_n_entity_out[i_part], "_p_db_entity_data_out ::");
      }
    }

    int *p_db_entity2_idx = PDM_array_new_idx_from_sizes_int(p_db_entity_strd[i_part], pn_twin_itrf[i_part]);

    // > Realloc entry args
    if (debug==1) {log_trace("\n_p_db_n_entity_out[i_part] = %d\n", _p_db_n_entity_out[i_part]);}
    if (debug==1) {log_trace("\n pn_twin_itrf     [i_part] = %d\n",  pn_twin_itrf     [i_part]);}
    PDM_realloc(_p_db_entity_gnum_out[i_part], _p_db_entity_gnum_out[i_part],   (_p_db_n_entity_out[i_part]+pn_twin_itrf[i_part]), PDM_g_num_t);
    PDM_realloc(_p_db_entity_data_out[i_part], _p_db_entity_data_out[i_part], 2*(_p_db_n_entity_out[i_part]+pn_twin_itrf[i_part]), PDM_g_num_t);
    if (has_sens==1) {
      PDM_realloc(_p_db_entity_sens_out[i_part], _p_db_entity_sens_out[i_part],   (_p_db_n_entity_out[i_part]+pn_twin_itrf[i_part]), int      );
    }

    int n_twin_itrf = 0;
    int idx_write = _p_db_n_entity_out[i_part];
    for(int i = 0; i < _p_db_n_entity_out[i_part]; ++i) {
      int i_pos = pentity_twin_itrf_lnum[i_part][i];
      if (debug_loop==1) {log_trace("\ni = %d --> i_pos = %d\n", i, i_pos);}

      int cur_gnum = pentity_gnum[i_part][i_pos];
      int cur_itrf = pentity_itrf[i_part][i_pos];
      if (debug_loop==1) {log_trace("\t -> gnum = %d ; interf = %d \n", cur_gnum, cur_itrf);}

      for(int j = pentity_twin_itrf_idx[i_part][i_pos]; j < pentity_twin_itrf_idx[i_part][i_pos+1]; ++j) {
        for(int k = p_db_entity2_idx[j]; k < p_db_entity2_idx[j+1]; ++k) {
          PDM_g_num_t opp_gnum = p_db_entity_data[i_part][2*k  ];
          int         opp_itrf = p_db_entity_data[i_part][2*k+1];
          if (debug_loop==1) {log_trace("\t\t candidate %d/%d: (%d,%d) ", k-p_db_entity2_idx[j]+1,
                                                        p_db_entity2_idx[j+1]-p_db_entity2_idx[j],
                                                        opp_gnum, opp_itrf);}
          if (cur_itrf==-opp_itrf) {
            if (debug_loop==1) {log_trace("is valid (idx_write=%d)", idx_write);}
            _p_db_entity_gnum_out[i_part][  idx_write  ] =  cur_gnum;
            _p_db_entity_data_out[i_part][2*idx_write  ] =  opp_gnum;
            _p_db_entity_data_out[i_part][2*idx_write+1] =  pentity_twin_itrf_itrf[i_part][j];
            if (has_sens==1) {
              _p_db_entity_sens_out[i_part][ idx_write ] =  pentity_twin_itrf_sens[i_part][j];
            }

            idx_write++;
            n_twin_itrf++;
          }
          if (debug_loop==1) {log_trace("\n");}
        }
      }
    }
    _p_db_n_entity_out[i_part] += n_twin_itrf; // Some interface can be wrong
    PDM_realloc(_p_db_entity_gnum_out[i_part], _p_db_entity_gnum_out[i_part],   _p_db_n_entity_out[i_part], PDM_g_num_t);
    PDM_realloc(_p_db_entity_data_out[i_part], _p_db_entity_data_out[i_part], 2*_p_db_n_entity_out[i_part], PDM_g_num_t);
    if (has_sens==1) {
      PDM_realloc(_p_db_entity_sens_out[i_part], _p_db_entity_sens_out[i_part],   _p_db_n_entity_out[i_part], int        );
    }

    PDM_free(p_db_entity2_idx);

    if (debug==1) {
      PDM_log_trace_array_long(_p_db_entity_gnum_out[i_part],   _p_db_n_entity_out[i_part], "_p_db_entity_gnum_out ::");
      PDM_log_trace_array_long(_p_db_entity_data_out[i_part], 2*_p_db_n_entity_out[i_part], "_p_db_entity_data_out ::");
      if (has_sens==1) {
        PDM_log_trace_array_int (_p_db_entity_sens_out[i_part],   _p_db_n_entity_out[i_part], "_p_db_entity_sens_out ::");
      }
    }
    
    PDM_free(p_db_entity_data[i_part]);
    if (has_sens==1) {
      PDM_free(p_db_entity_sens[i_part]);
    }
    PDM_free(p_db_entity_strd[i_part]);
  }
  
  PDM_free(p_db_entity_data);
  if (has_sens==1) {
    PDM_free(p_db_entity_sens);
  }
  PDM_free(p_db_entity_strd);
  PDM_free(pn_twin_itrf);

  if (debug==1) {
    log_trace("END _find_twin_interface\n");
    log_trace("------------------------\n");
    log_trace("\n");
  }

}

static
void
_find_valid_entities
(
  int            n_part,
  int           *pn_entity1,
  int          **pentity1_itrf,
  int          **pentity1_entity2_idx,
  PDM_g_num_t  **pentity1_entity2_gnum,
  int          **pentity1_entity2_triplet,
  int           *pentity1_entity2_ntot,
  int           *pn_entity2,
  PDM_g_num_t  **pentity2_gnum,
  int            d_db_n_entity2, // Distributed data_base
  PDM_g_num_t   *d_db_entity2_gnum,
  int           *d_db_entity2_ancstr_strd,
  PDM_g_num_t   *d_db_entity2_ancstr,
  int           *d_db_entity2_path_itrf_strd,
  int           *d_db_entity2_path_itrf,
  int           *d_db_entity2_strd,
  PDM_g_num_t   *d_db_entity2_data,
  int           *d_db_entity2_sens,
  int         ***pentity2_kind_out,
  int         ***pentity2_lnum_out,
  int         ***pentity2_itrf_out,
  int         ***pentity2_sens_out,
  PDM_g_num_t ***pentity2_ancstr_out,
  int         ***pentity2_path_itrf_strd_out,
  int         ***pentity2_path_itrf_out,
  int         ***pentity2_kind_idx_out,
  int         ***pentity2_kind_ordr_out,
  int         ***pentity2_twin_itrf_idx_out,
  PDM_g_num_t ***pentity2_twin_itrf_gnum_out,
  int         ***pentity2_twin_itrf_itrf_out,
  int         ***pentity2_twin_itrf_sens_out,
  PDM_MPI_Comm   comm
)
{
  /*
   * TODO:
   *   - add pentity1_entity2_gnum sorted in args
   */

  int debug      = 0;
  int debug_loop = 0;
  if (debug==1) {
    log_trace("\n");
    log_trace("------------------------\n");
    log_trace("BEG _find_valid_entities\n");
    PDM_log_trace_array_int(pentity1_entity2_ntot, n_part, "pn_entity2 ::");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity1_entity2_idx[i_part], pn_entity1[i_part]+1, "pentity1_entity2_idx ::");
    }
    
    log_trace("\n");
    int n_ancstr = 0;
    int n_itrf   = 0;
    for (int i_entity=0; i_entity<d_db_n_entity2; ++i_entity) {
      n_ancstr+=d_db_entity2_ancstr_strd   [i_entity];
      n_itrf  +=d_db_entity2_path_itrf_strd[i_entity];
    }
    PDM_log_trace_array_long(d_db_entity2_gnum          , d_db_n_entity2, "d_db_entity2_gnum           ::");
    PDM_log_trace_array_int (d_db_entity2_ancstr_strd   , d_db_n_entity2, "d_db_entity2_ancstr_strd    ::");
    PDM_log_trace_array_long(d_db_entity2_ancstr        , n_ancstr      , "d_db_entity2_ancstr         ::");
    PDM_log_trace_array_int (d_db_entity2_path_itrf_strd, d_db_n_entity2, "d_db_entity2_path_itrf_strd ::");
    PDM_log_trace_array_int (d_db_entity2_path_itrf     , n_itrf        , "d_db_entity2_path_itrf      ::");
  }
  


  /*
   * Get informations from data_base for each element
   *
   */
  int has_sens = 0;
  if (d_db_entity2_sens!=NULL) {
    has_sens = 1;
  }


  if (debug==1) {
    log_trace("\n");
    log_trace("has_sens = %d\n", has_sens);
    log_trace("\n");
  }

  PDM_block_to_part_t* btp = PDM_block_to_part_create_from_sparse_block(d_db_entity2_gnum,
                                                                        d_db_n_entity2,
                                                (const PDM_g_num_t **)  pentity1_entity2_gnum,
                                                                        pentity1_entity2_ntot,
                                                                        n_part,
                                                                        comm);

  int         **p_db_entity2_strd = NULL;
  PDM_g_num_t **p_db_entity2_data = NULL;
  PDM_block_to_part_exch(btp,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         d_db_entity2_strd,
                         d_db_entity2_data,
                        &p_db_entity2_strd,
         (void ***)     &p_db_entity2_data);

  int **p_db_entity2_strd2 = NULL;
  int **p_db_entity2_sens  = NULL;
  if (has_sens==1) {
    PDM_block_to_part_exch(btp,
                           1 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           d_db_entity2_strd,
                           d_db_entity2_sens,
                          &p_db_entity2_strd2,
           (void ***)     &p_db_entity2_sens);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_free(p_db_entity2_strd2[i_part]);
    }
    PDM_free(p_db_entity2_strd2);
  }

  // > Var stride here cause some entities can be missing from init db
  int         **p_db_entity2_ancstr_strd = NULL;
  PDM_g_num_t **p_db_entity2_ancstr = NULL;
  PDM_block_to_part_exch(btp,
                         1 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         d_db_entity2_ancstr_strd,
                         d_db_entity2_ancstr,
                        &p_db_entity2_ancstr_strd,
         (void ***)     &p_db_entity2_ancstr);

  int **p_db_entity2_path_itrf_strd = NULL;
  int **p_db_entity2_path_itrf      = NULL;
  PDM_block_to_part_exch(btp,
                         1 * sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         d_db_entity2_path_itrf_strd,
                         d_db_entity2_path_itrf,
                        &p_db_entity2_path_itrf_strd,
         (void ***)     &p_db_entity2_path_itrf);
  PDM_block_to_part_free(btp);


  /*
   * Need to sort pentity_gnum to search if a candidate entity
   * doesn't already exist in partition.
   *
   */
  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      int n_data   = 0;
      int n_ancstr = 0;
      int n_itrf   = 0;
      for (int i_entity = 0; i_entity < pentity1_entity2_ntot[i_part]; ++i_entity) {
        n_data   += p_db_entity2_strd[i_part][i_entity];
        n_ancstr += p_db_entity2_ancstr_strd[i_part][i_entity];
        n_itrf   += p_db_entity2_path_itrf_strd[i_part][i_entity];
      }
      PDM_log_trace_array_long(pentity2_gnum              [i_part], pn_entity2           [i_part], "pentity2_gnum               ::");
      PDM_log_trace_array_int (p_db_entity2_strd          [i_part], pentity1_entity2_ntot[i_part], "p_db_entity2_strd           ::");
      PDM_log_trace_array_long(p_db_entity2_data          [i_part], 2*n_data                     , "p_db_entity2_data           ::");
      if (has_sens==1) {
        PDM_log_trace_array_int (p_db_entity2_sens          [i_part],   n_data                     , "p_db_entity2_sens           ::");
      }
      PDM_log_trace_array_int (p_db_entity2_ancstr_strd   [i_part], pentity1_entity2_ntot[i_part], "p_db_entity2_ancstr_strd    ::");
      PDM_log_trace_array_long(p_db_entity2_ancstr        [i_part], n_ancstr                     , "p_db_entity2_ancstr         ::");
      PDM_log_trace_array_int (p_db_entity2_path_itrf_strd[i_part], pentity1_entity2_ntot[i_part], "p_db_entity2_path_itrf_strd ::");
      PDM_log_trace_array_int (p_db_entity2_path_itrf     [i_part], n_itrf                       , "p_db_entity2_path_itrf      ::");
    }
  }

  int         **pentity2_ordr        = NULL;
  PDM_g_num_t **pentity2_gnum_sorted = NULL;
  PDM_malloc(pentity2_ordr       , n_part, int         *);
  PDM_malloc(pentity2_gnum_sorted, n_part, PDM_g_num_t *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_malloc(pentity2_ordr       [i_part], pn_entity2[i_part], int        );
    PDM_malloc(pentity2_gnum_sorted[i_part], pn_entity2[i_part], PDM_g_num_t);
    for(int i_entity = 0; i_entity < pn_entity2[i_part]; ++i_entity) {
      pentity2_ordr       [i_part][i_entity] = i_entity;
      pentity2_gnum_sorted[i_part][i_entity] = pentity2_gnum[i_part][i_entity];
    }
    PDM_sort_long(pentity2_gnum_sorted[i_part], pentity2_ordr[i_part], pn_entity2[i_part]);
  }
  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      int n_data = 0;
      for (int i_entity = 0; i_entity < pentity1_entity2_ntot[i_part]; ++i_entity) {
        n_data += p_db_entity2_ancstr_strd[i_part][i_entity];
      }
      PDM_log_trace_array_long(pentity1_entity2_gnum   [i_part], pentity1_entity2_ntot  [i_part], "pentity1_entity2_gnum ::");
      PDM_log_trace_array_int (pentity1_entity2_triplet[i_part], 3*pentity1_entity2_ntot[i_part], "pentity1_entity2_trplt::");
      PDM_log_trace_array_int (pentity2_ordr           [i_part], pn_entity2             [i_part], "pentity2_ordr         ::");
      PDM_log_trace_array_long(pentity2_gnum_sorted    [i_part], pn_entity2             [i_part], "pentity2_gnum_sorted  ::");
    }
  }



  /*
   * Post-treatment - find if new entity requested (parent, itrf) not already in data_base
   *
   * Possible kind :
   *    1/ Internal
   *       a/ Already know in current partition                    : 0
   *       b/ New in current partition                             : 2
   *    2/ From interface
   *       a/ Know by relation table and know in current partition : 1
   *       b/ Know by relation table but not local                 : 3
   *       c/ Unknown by relation table but from other part        : 4
   *       d/ New                                                  : 5
   * On essaye après d'organiser les nouvelles entités dans l'ordre suivant :
   *    [2/3/4]
   */

  int n_kind = 6;

  // > Allocate
  int **pentity2_kind = NULL;
  int **pentity2_lnum = NULL;
  int **pentity2_itrf = NULL;
  int **pentity2_sens = NULL;
  PDM_malloc(pentity2_kind, n_part, int *);
  PDM_malloc(pentity2_lnum, n_part, int *);
  PDM_malloc(pentity2_itrf, n_part, int *);
  if (has_sens==1) {
    PDM_malloc(pentity2_sens, n_part, int *);
  }

  int         **pentity2_twin_itrf_idx  = NULL;
  PDM_g_num_t **pentity2_twin_itrf_gnum = NULL;
  int         **pentity2_twin_itrf_itrf = NULL;
  int         **pentity2_twin_itrf_sens = NULL;
  PDM_malloc(pentity2_twin_itrf_idx , n_part, int         *);
  PDM_malloc(pentity2_twin_itrf_gnum, n_part, PDM_g_num_t *);
  PDM_malloc(pentity2_twin_itrf_itrf, n_part, int         *);
  if (has_sens==1) {
    PDM_malloc(pentity2_twin_itrf_sens, n_part, int       *);
  }

  PDM_g_num_t **pentity2_ancstr         = NULL;
  int         **pentity2_path_itrf_strd = NULL;
  int         **pentity2_path_itrf      = NULL;
  PDM_malloc(pentity2_ancstr        , n_part, PDM_g_num_t *);
  PDM_malloc(pentity2_path_itrf_strd, n_part, int         *);
  PDM_malloc(pentity2_path_itrf     , n_part, int         *);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    
    // > Allocate
    //TODO: fonciton qui alloue plusieurs tableaux en même temps
    pentity2_kind[i_part] = PDM_array_const_int(pentity1_entity2_ntot[i_part], -1);
    pentity2_lnum[i_part] = PDM_array_const_int(pentity1_entity2_ntot[i_part],  0);
    pentity2_itrf[i_part] = PDM_array_const_int(pentity1_entity2_ntot[i_part], -1);
    if (has_sens==1) {
      pentity2_sens[i_part] = PDM_array_const_int(pentity1_entity2_ntot[i_part],  1);
    }

    pentity2_ancstr        [i_part] = PDM_array_const_gnum(pentity1_entity2_ntot[i_part],  0);
    pentity2_path_itrf_strd[i_part] = PDM_array_const_int (pentity1_entity2_ntot[i_part],  0);
  
    PDM_malloc(pentity2_twin_itrf_idx[i_part], pentity1_entity2_ntot[i_part]+1, int);
    pentity2_twin_itrf_idx[i_part][0] = 0;

    int len_path_itrf_tot = 0;
    int idx_read_data  = 0; // read index in received entity2 data_base infos
    int lidx_read_data = 0;
    int i_read_ancstr = 0; // read index in received entity2 data_base ancstr infos

    for (int i = 0; i < pn_entity1[i_part]; ++i) {
      if (debug_loop==1) {log_trace("\n\ni_entity1 = %d -> interf = %d \n", i, pentity1_itrf[i_part][i]);}

      if(pentity1_itrf[i_part][i] != 0) {

        for(int j = pentity1_entity2_idx[i_part][i]; j < pentity1_entity2_idx[i_part][i+1]; ++j) {
          PDM_g_num_t opp_gnum     = 0;
          int         opp_itrf     = 0;
          int         opp_itrf_sgn = 0;
          int         opp_sens     = 1;
          
          // Entity2 informations
          PDM_g_num_t cur_gnum         = pentity1_entity2_gnum[i_part][j];
          int         cur_trplt_proc   = pentity1_entity2_triplet[i_part][3*j  ];
          int         cur_trplt_part   = pentity1_entity2_triplet[i_part][3*j+1];
          int         cur_trplt_entity = pentity1_entity2_triplet[i_part][3*j+2];
          int         cur_itrf     = PDM_ABS (pentity1_itrf[i_part][i]);
          int         cur_itrf_sgn = PDM_SIGN(pentity1_itrf[i_part][i]);
          if (debug_loop==1) {log_trace("\n\t j_entity2 =%d -> gnum = ("PDM_FMT_G_NUM",%i) trplt = (%d, %d, %d)\n", j, cur_gnum, cur_itrf_sgn*cur_itrf, cur_trplt_proc, cur_trplt_part, cur_trplt_entity);}

          // Candidate for entity2 (from data_base)
          int known_in_db = 0;
          lidx_read_data = idx_read_data;
          for(int k = 0; k < p_db_entity2_strd[i_part][j]; ++k) {
          
            // Candidate informations
            opp_gnum     =          p_db_entity2_data[i_part][2*lidx_read_data  ];
            opp_itrf     = PDM_ABS (p_db_entity2_data[i_part][2*lidx_read_data+1]);
            opp_itrf_sgn = PDM_SIGN(p_db_entity2_data[i_part][2*lidx_read_data+1]);
            if (has_sens==1) {
              opp_sens   =          p_db_entity2_sens[i_part][  lidx_read_data  ];
            }
            if (debug_loop==1) {log_trace("\t\t candidate %d/%d: ("PDM_FMT_G_NUM"/%d) \n", k+1, p_db_entity2_strd[i_part][j], opp_gnum, opp_itrf_sgn*opp_itrf);}
            
            // Check if candidate have same interface that entity2, if so it (entity2, itrf) alrdy know
            // NOTE: je commente && gnum_opp == 0 (normalement on peut pas avoir 2 candidats valides différents?)
            if(cur_itrf == opp_itrf && cur_itrf_sgn == -opp_itrf_sgn) {// && gnum_opp == 0) { // Sign a unifier avec l'autre
              known_in_db = 1;
              break;
            }
  
            lidx_read_data++;
          }


          if (known_in_db==0) {
            if (debug_loop==1) {log_trace("\t\t ===> Not found = ("PDM_FMT_G_NUM",%i), j = %d \n", cur_gnum, cur_itrf_sgn*cur_itrf, j);}
            
            if (cur_trplt_part==i_part) { // From the same partition (domain)
              pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
              pentity2_kind[i_part][j] = 5;
              
              if (p_db_entity2_ancstr_strd[i_part][j]==0) { // not in db now, wo it means that it is an ancstr
                pentity2_ancstr[i_part][j] = cur_gnum;
              }
              else {
                pentity2_ancstr[i_part][j] = p_db_entity2_ancstr[i_part][i_read_ancstr];
              }
              pentity2_path_itrf_strd[i_part][j] = p_db_entity2_path_itrf_strd[i_part][j];
              
              len_path_itrf_tot += p_db_entity2_path_itrf_strd[i_part][j];

              // > Find other interfaces of the entity
              pentity2_twin_itrf_idx[i_part][j+1] = pentity2_twin_itrf_idx[i_part][j];
              lidx_read_data = idx_read_data;
              for(int k = 0; k < p_db_entity2_strd[i_part][j]; ++k) {
                opp_gnum     =          p_db_entity2_data[i_part][2*lidx_read_data  ];
                opp_itrf     = PDM_ABS (p_db_entity2_data[i_part][2*lidx_read_data+1]);
                opp_itrf_sgn = PDM_SIGN(p_db_entity2_data[i_part][2*lidx_read_data+1]);
                if (debug_loop==1) {log_trace("\t\t\t lidx_read_data: %d \n", lidx_read_data);}
                if (debug_loop==1) {log_trace("\t\t\t twin itrf candidate: (%d,%d) \n", opp_gnum, opp_itrf_sgn*opp_itrf);}
                if (cur_itrf!=opp_itrf) {
                  if (debug_loop==1) {log_trace("\t\t\t twin itrf: (%d,%d) \n", opp_gnum, opp_itrf_sgn*opp_itrf);}

                  pentity2_twin_itrf_idx[i_part][j+1]++;
                }
                lidx_read_data++;
              }
            }

            else { // From a different partition (domain)
              if (debug_loop==1) {log_trace("\t\t ===> From other part = ("PDM_FMT_G_NUM",part=%d)\n", cur_gnum, cur_trplt_part);}

              // PDM_log_trace_array_long(pentity2_gnum_sorted  [i_part], pn_entity2[i_part], "pentity2_gnum_sorted   ::");
              int pos_int = PDM_binary_search_long(cur_gnum, pentity2_gnum_sorted[i_part], pn_entity2[i_part]);
              if( pos_int != -1) {
                if (debug_loop==1) {log_trace("\t\t ===> And known in partition = (i_pos=%d) sens = %d\n", pos_int, opp_sens);}
                assert((has_sens==0)||(opp_sens==-1||opp_sens==1));
                pentity2_lnum[i_part][j] = opp_sens*(pentity2_ordr[i_part][pos_int]+1);
                pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
                pentity2_kind[i_part][j] = 1; //Bluff
                if (has_sens==1) {
                  pentity2_sens[i_part][j] = opp_sens; //Bluff2
                }
              }
              else {

                //TODO: vérifier si dans les candidats y'a pas un sommet dans la partition -> pb comment savoir que c'est le bon ?

                // Is there any candidate alrdy in partition. If yes, hope its to good one
                int pos_int_c = -1;
                opp_sens  =  0;
                lidx_read_data = idx_read_data;
                for(int k = 0; k < p_db_entity2_strd[i_part][j]; ++k) {
                
                  // Candidate informations
                  if (has_sens==1) {
                    opp_sens = p_db_entity2_sens[i_part][lidx_read_data];
                  }
                  else {
                    opp_sens = 1;
                  }
                  opp_gnum = p_db_entity2_data  [i_part][2*lidx_read_data  ];
                  if (debug_loop==1) {log_trace("\t\t\t candidate %d/%d: ("PDM_FMT_G_NUM") \n", k+1, p_db_entity2_strd[i_part][j], opp_gnum);}
                  pos_int_c = PDM_binary_search_long(opp_gnum, pentity2_gnum_sorted[i_part], pn_entity2[i_part]);
                  if (pos_int_c!=-1) {
                    break;
                  }
                  lidx_read_data++;
                }
                lidx_read_data = idx_read_data+p_db_entity2_strd[i_part][j];

                if( pos_int_c != -1) {
                  assert((has_sens==0)||(opp_sens==-1||opp_sens==1));
                  if (debug_loop==1) {log_trace("\t\t ===> Twin known in partition = (i_pos=%d) sens = %d\n", pos_int_c, opp_sens);}
                  pentity2_lnum[i_part][j] = opp_sens*(pentity2_ordr[i_part][pos_int_c]+1);
                  pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
                  pentity2_kind[i_part][j] = 1; //Bluff
                  // pentity2_sens[i_part][j] = opp_sens; //Bluff2 ?
                }

                else { 
                  /**
                   * Entity from other partition (domain) is unknown in current partition:
                   *   -> get the entity with the original gnum
                   */

                  if (debug_loop==1) {log_trace("\t\t ===> Unknown in partition = (i_pos=%d)\n", pos_int);}
                  pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
                  pentity2_kind[i_part][j] = 4;

                  if (p_db_entity2_ancstr_strd[i_part][j]==0) { // not in db now, wo it means that it is an ancstr
                    pentity2_ancstr[i_part][j] = cur_gnum;
                  }
                  else {
                    pentity2_ancstr[i_part][j] = p_db_entity2_ancstr[i_part][i_read_ancstr];
                  }

                  pentity2_path_itrf_strd[i_part][j] = p_db_entity2_path_itrf_strd[i_part][j];
                  len_path_itrf_tot += p_db_entity2_path_itrf_strd[i_part][j];
                  // TODO: do we need to get twin interface ? -> probably not
                }

              }
              pentity2_twin_itrf_idx[i_part][j+1] = pentity2_twin_itrf_idx[i_part][j];
            }
          }

          else { // Known in data_base

            // Search if entity2 is known in partition
            int pos_int = PDM_binary_search_long(opp_gnum, pentity2_gnum_sorted[i_part], pn_entity2[i_part]);
            if( pos_int != -1) {
              // opp_sens = -1;
              // if (has_sens==1) {
              //   opp_sens = p_db_entity2_sens[i_part][lidx_read_data];
              // }
              // else {
              //   opp_sens = 1;
              // }
              if (debug_loop==1) {log_trace("\t\t ===> Local found = ("PDM_FMT_G_NUM",%i) sens = %d -> kind = 1\n", opp_gnum, opp_itrf_sgn*opp_itrf, opp_sens);}

              assert((has_sens==0)||(opp_sens==-1||opp_sens==1));
              pentity2_lnum[i_part][j] = opp_sens*(pentity2_ordr[i_part][pos_int]+1);
              pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
              pentity2_kind[i_part][j] = 1;
                // pentity2_sens[i_part][j] = opp_sens;
            }
            else {
              if (debug_loop==1) {log_trace("\t\t ===> Not found = ("PDM_FMT_G_NUM",%i) -> kind = 3 (sens = %d)\n", opp_gnum, opp_itrf_sgn*opp_itrf, opp_sens);}
              pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
              pentity2_kind[i_part][j] = 3;
              
              // entity2 take the opposite gnum from db
              pentity1_entity2_gnum[i_part][j] = opp_gnum;
              if (has_sens==1) {
                pentity2_sens[i_part][j] = opp_sens;
              }
            }

            // > Other interfaces of the entity not interesting
            pentity2_twin_itrf_idx[i_part][j+1] = pentity2_twin_itrf_idx[i_part][j];
          }

          // Increment reading
          idx_read_data += p_db_entity2_strd[i_part][j];
          i_read_ancstr+=p_db_entity2_ancstr_strd[i_part][j];
          // log_trace("i_read_ancstr = %d\n", i_read_ancstr)
        }

      }
      else {

        for(int j = pentity1_entity2_idx[i_part][i]; j < pentity1_entity2_idx[i_part][i+1]; ++j) {
         
          // Entity2 informations
          PDM_g_num_t cur_gnum = pentity1_entity2_gnum[i_part][j];
          if (debug_loop==1) {log_trace("\n\t j_entity2 =%d -> gnum = ("PDM_FMT_G_NUM",%i) \n", j, cur_gnum, pentity1_itrf[i_part][i]);}
          
          // Search if entity2 is known in partition
          int pos_int = PDM_binary_search_long(cur_gnum, pentity2_gnum_sorted[i_part], pn_entity2[i_part]);
        
          if( pos_int != -1) {
            // log_trace("cur_gnum = %d -> pos_int = %d\n", cur_gnum, pos_int);
            pentity2_lnum[i_part][j] = pentity2_ordr[i_part][pos_int]+1;
            pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
            pentity2_kind[i_part][j] = 0;
            if (debug_loop==1) {log_trace("\t\t ===> Local found = ("PDM_FMT_G_NUM",%d) -> kind = 0\n", cur_gnum, pentity1_itrf[i_part][i]);}
          } else {
            if (debug_loop==1) {log_trace("\t\t ===> Not found -> kind = 2\n");}
            pentity2_itrf[i_part][j] = pentity1_itrf[i_part][i];
            pentity2_kind[i_part][j] = 2;
          }
          
          // > Other interfaces of the entity not interesting
          pentity2_twin_itrf_idx[i_part][j+1] = pentity2_twin_itrf_idx[i_part][j];
          
          // Increment reading
          idx_read_data += p_db_entity2_strd[i_part][j];
          i_read_ancstr+=p_db_entity2_ancstr_strd[i_part][j];
          // log_trace("i_read_ancstr = %d\n", i_read_ancstr)
        }
      }
    }

    if (debug==1) {
      PDM_log_trace_array_long(pentity1_entity2_gnum  [i_part], pentity1_entity2_ntot[i_part], "pentity1_entity2_gnum   ::");
      PDM_log_trace_array_int (pentity2_lnum          [i_part], pentity1_entity2_ntot[i_part], "pentity2_lnum           ::");
      PDM_log_trace_array_int (pentity2_itrf          [i_part], pentity1_entity2_ntot[i_part], "pentity2_itrf           ::");
      PDM_log_trace_array_int (pentity2_kind          [i_part], pentity1_entity2_ntot[i_part], "pentity2_kind           ::");
      if (has_sens==1) {
        PDM_log_trace_array_int (pentity2_sens          [i_part], pentity1_entity2_ntot[i_part], "pentity2_sens           ::");
      }
      PDM_log_trace_array_long(pentity2_ancstr        [i_part], pentity1_entity2_ntot[i_part], "pentity2_ancstr         ::");
      PDM_log_trace_array_int (pentity2_path_itrf_strd[i_part], pentity1_entity2_ntot[i_part], "pentity2_path_itrf_strd ::");
    }


    // > Allocate
    int n_twin_interface = pentity2_twin_itrf_idx[i_part][pentity1_entity2_ntot[i_part]];
    PDM_malloc(pentity2_twin_itrf_gnum[i_part], n_twin_interface, PDM_g_num_t);
    PDM_malloc(pentity2_twin_itrf_itrf[i_part], n_twin_interface, int        );
    if (has_sens==1) {
      PDM_malloc(pentity2_twin_itrf_sens[i_part], n_twin_interface, int        );
    }

    PDM_malloc(pentity2_path_itrf[i_part], len_path_itrf_tot, int        );
    
    idx_read_data  = 0; // read index in received entity2 data_base infos
    int idx_read_itrf  = 0; // read index in received entity2 data_base infos
    int idx_write_data = 0;
    int idx_write_itrf = 0;

    for (int i = 0; i < pn_entity1[i_part]; ++i) {
      if (debug_loop==1) {log_trace("\n\ni_entity1 = %d -> interf = %d \n", i, pentity1_itrf[i_part][i]);}

      for(int j = pentity1_entity2_idx[i_part][i]; j < pentity1_entity2_idx[i_part][i+1]; ++j) {

        if (pentity2_kind[i_part][j]==5) {
          PDM_g_num_t cur_gnum = pentity1_entity2_gnum[i_part][j];
          int         cur_itrf = PDM_ABS(pentity1_itrf[i_part][i]);
          if (debug_loop==1) {log_trace("\n\t j_entity2 =%d -> gnum = ("PDM_FMT_G_NUM",%i) \n", j, cur_gnum, cur_itrf);}

          for(int k = 0; k < p_db_entity2_strd[i_part][j]; ++k) {
          
            // Candidate informations
            PDM_g_num_t opp_gnum     =          p_db_entity2_data[i_part][2*idx_read_data  ];
            int         opp_itrf     = PDM_ABS (p_db_entity2_data[i_part][2*idx_read_data+1]);
            int         opp_itrf_sgn = PDM_SIGN(p_db_entity2_data[i_part][2*idx_read_data+1]);
            int         opp_sens     = 0;
            if (has_sens==1) {
              opp_sens = p_db_entity2_sens[i_part][  idx_read_data  ];
            }

            if (debug_loop==1) {log_trace("\t\t candidate %d%d: ("PDM_FMT_G_NUM"/"PDM_FMT_G_NUM") \n", k+1, p_db_entity2_strd[i_part][j], opp_gnum, opp_itrf_sgn*opp_itrf);}
            if (cur_itrf!=opp_itrf) {
              if (debug_loop==1) {log_trace("\t\t with = (%d,%d)\n", opp_gnum, opp_itrf);}
              if (debug_loop==1) {log_trace("\t\t idx_write_data = %d\n", idx_write_data);}
              pentity2_twin_itrf_gnum[i_part][idx_write_data] = opp_gnum;
              pentity2_twin_itrf_itrf[i_part][idx_write_data] = opp_itrf_sgn*opp_itrf;
              if (has_sens==1) {
                pentity2_twin_itrf_sens[i_part][idx_write_data] = opp_sens;
              }
              idx_write_data++;
            }
            idx_read_data++;
          }
          for(int k = 0; k < p_db_entity2_path_itrf_strd[i_part][j]; ++k) {
            if (debug_loop==1) {log_trace("idx_write_itrf = %d -> interf = %d \n", idx_write_itrf, p_db_entity2_path_itrf[i_part][idx_read_itrf]);}
            pentity2_path_itrf[i_part][idx_write_itrf++] = p_db_entity2_path_itrf[i_part][idx_read_itrf];
            idx_read_itrf++;
          }

        }
        else {
          idx_read_data += p_db_entity2_strd          [i_part][j];
          idx_read_itrf += p_db_entity2_path_itrf_strd[i_part][j];
        }
      }
    }

    if (debug==1) {
      PDM_log_trace_array_int (pentity2_twin_itrf_idx [i_part], pentity1_entity2_ntot[i_part]+1, "pentity2_twin_itrf_idx  ::");
      PDM_log_trace_array_long(pentity2_twin_itrf_gnum[i_part], n_twin_interface               , "pentity2_twin_itrf_gnum ::");
      PDM_log_trace_array_int (pentity2_twin_itrf_itrf[i_part], n_twin_interface               , "pentity2_twin_itrf_itrf ::");
      if (has_sens==1) {
        PDM_log_trace_array_int (pentity2_twin_itrf_sens[i_part], n_twin_interface               , "pentity2_twin_itrf_sens ::");
      }
      PDM_log_trace_array_int (pentity2_path_itrf     [i_part], len_path_itrf_tot              , "pentity2_path_itrf      ::");
    }
  }


  // > Free
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity2_ordr       [i_part]);
    PDM_free(pentity2_gnum_sorted[i_part]);

    PDM_free(p_db_entity2_strd[i_part]);
    PDM_free(p_db_entity2_data[i_part]);
    if (has_sens==1) {
      PDM_free(p_db_entity2_sens[i_part]);
    }

    PDM_free(p_db_entity2_ancstr_strd[i_part]);
    PDM_free(p_db_entity2_ancstr[i_part]);

    PDM_free(p_db_entity2_path_itrf_strd[i_part]);
    PDM_free(p_db_entity2_path_itrf[i_part]);
  }
  PDM_free(pentity2_ordr       );
  PDM_free(pentity2_gnum_sorted);
 
  PDM_free(p_db_entity2_strd);
  PDM_free(p_db_entity2_data);
  if (has_sens==1) {
    PDM_free(p_db_entity2_sens);
  }

  PDM_free(p_db_entity2_ancstr_strd);
  PDM_free(p_db_entity2_ancstr);

  PDM_free(p_db_entity2_path_itrf_strd);
  PDM_free(p_db_entity2_path_itrf);



  /*
   * Sort by kind and keep order
   */
  int **pentity2_kind_idx  = NULL;
  int **pentity2_kind_ordr = NULL;
  PDM_malloc(pentity2_kind_idx , n_part, int *);
  PDM_malloc(pentity2_kind_ordr, n_part, int *);

  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Initialize
    int *pentity2_kind_cp = NULL;
    PDM_malloc(pentity2_kind_cp          , pentity1_entity2_ntot[i_part], int);
    PDM_malloc(pentity2_kind_ordr[i_part], pentity1_entity2_ntot[i_part], int);
    for(int i = 0; i < pentity1_entity2_ntot[i_part]; ++i) {
      pentity2_kind_ordr[i_part][i] = i;
      pentity2_kind_cp[i] = pentity2_kind[i_part][i];
    }

    // Sort
    PDM_sort_int(pentity2_kind_cp, pentity2_kind_ordr[i_part], pentity1_entity2_ntot[i_part]);

    // Count and generate idx
    int *pentity2_kind_n      = PDM_array_zeros_int(n_kind);
    pentity2_kind_idx[i_part] = PDM_array_zeros_int(n_kind+1);
    for(int i = 0; i < pentity1_entity2_ntot[i_part]; ++i) {
      pentity2_kind_n[pentity2_kind[i_part][i]]++;
    }
    for(int i = 0; i < n_kind; ++i) {
      pentity2_kind_idx[i_part][i+1] = pentity2_kind_idx[i_part][i] + pentity2_kind_n[i];
    }


    if (debug==1) {
      PDM_log_trace_array_int( pentity2_kind_n          , n_kind            , "pentity2_kind_n    ::");
      PDM_log_trace_array_int(pentity2_kind_idx [i_part], n_kind+1          , "pentity2_kind_idx  ::");
      PDM_log_trace_array_int(pentity2_kind_ordr[i_part], pentity1_entity2_ntot[i_part], "pentity2_kind_ordr ::");
    }

    PDM_free(pentity2_kind_n);
    PDM_free(pentity2_kind_cp);
  }



  // > Return
  *pentity2_kind_out      = pentity2_kind;
  *pentity2_lnum_out      = pentity2_lnum;
  *pentity2_itrf_out      = pentity2_itrf;
  if (has_sens==1){
    *pentity2_sens_out      = pentity2_sens;
  }
  *pentity2_ancstr_out    = pentity2_ancstr;

  *pentity2_path_itrf_strd_out  = pentity2_path_itrf_strd;
  *pentity2_path_itrf_out       = pentity2_path_itrf;
  
  *pentity2_kind_idx_out  = pentity2_kind_idx;
  *pentity2_kind_ordr_out = pentity2_kind_ordr;

  *pentity2_twin_itrf_idx_out  = pentity2_twin_itrf_idx;
  *pentity2_twin_itrf_gnum_out = pentity2_twin_itrf_gnum;
  *pentity2_twin_itrf_itrf_out = pentity2_twin_itrf_itrf;
  if (has_sens==1) {
    *pentity2_twin_itrf_sens_out = pentity2_twin_itrf_sens;
  }
  else {
    PDM_free(pentity2_twin_itrf_sens);
  }

  if (debug==1) {
    log_trace("END _find_valid_entities\n");
    log_trace("------------------------\n");
    log_trace("\n");
  }

}



static
void
_enrich_block_interface
(
  int           n_part,
  int          *pn_entity_cur,
  PDM_g_num_t **cur_pentity_itrf_gnum,
  PDM_g_num_t **cur_pentity_itrf_data,
  int         **cur_pentity_itrf_sens,
  int           pn_entity_prev,
  PDM_g_num_t  *prev_pentity_itrf_gnum,
  int          *prev_pentity_itrf_strd,
  PDM_g_num_t  *prev_pentity_itrf_data,
  int          *prev_pentity_itrf_sens,
  int          *dn_entity_next,
  PDM_g_num_t **next_dentity_itrf_gnum,
  int         **next_dentity_itrf_strd,
  PDM_g_num_t **next_dentity_itrf_data,
  int         **next_dentity_itrf_sens,
  PDM_MPI_Comm  comm
)
{
  int debug      = 0;
  // int debug_loop = 0;

  if (debug==1) {
    log_trace("\n");
    log_trace("---------------------------\n");
    log_trace("BEG _enrich_block_interface\n");
  }

  int has_sens = 0;
  if (prev_pentity_itrf_sens!=NULL) {
    has_sens = 1;
  }

  /*
   * Concatenate array as one "domain"
   */
  int          *pn_entity                = NULL;
  PDM_g_num_t **concat_pentity_itrf_gnum = NULL;
  int         **concat_pentity_itrf_strd = NULL;
  PDM_g_num_t **concat_pentity_itrf_data = NULL;
  int         **concat_pentity_itrf_sens = NULL;
  double      **concat_pentity_itrf_wght = NULL;
  PDM_malloc(pn_entity               , n_part+1, int          );
  PDM_malloc(concat_pentity_itrf_gnum, n_part+1, PDM_g_num_t *);
  PDM_malloc(concat_pentity_itrf_strd, n_part+1, int         *);
  PDM_malloc(concat_pentity_itrf_data, n_part+1, PDM_g_num_t *);
  PDM_malloc(concat_pentity_itrf_sens, n_part+1, int         *);
  PDM_malloc(concat_pentity_itrf_wght, n_part+1, double      *);
    
  for(int i_part = 0; i_part < n_part; ++i_part) {
    if(debug == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(cur_pentity_itrf_gnum[i_part],   pn_entity_cur[i_part], "cur_pentity_itrf_gnum ::");
      PDM_log_trace_array_long(cur_pentity_itrf_data[i_part], 2*pn_entity_cur[i_part], "cur_pentity_itrf_data ::");
      if (has_sens==1) {
        PDM_log_trace_array_int(cur_pentity_itrf_sens[i_part],   pn_entity_cur[i_part], "cur_pentity_itrf_sens ::");
      }
    }

    pn_entity[i_part] = 2*pn_entity_cur[i_part];

    PDM_malloc(concat_pentity_itrf_gnum  [i_part],     pn_entity[i_part], PDM_g_num_t);
    PDM_malloc(concat_pentity_itrf_strd  [i_part],     pn_entity[i_part], int        );
    PDM_malloc(concat_pentity_itrf_data  [i_part], 2 * pn_entity[i_part], PDM_g_num_t);
    if (has_sens==1) {
      PDM_malloc(concat_pentity_itrf_sens[i_part],     pn_entity[i_part], int        );
    }
    PDM_malloc(concat_pentity_itrf_wght  [i_part],     pn_entity[i_part], double     );

    for(int i = 0; i < pn_entity_cur[i_part]; ++i) {
      PDM_g_num_t orig_gnum = cur_pentity_itrf_data[i_part][2*i  ];
      PDM_g_num_t i_itrf    = cur_pentity_itrf_data[i_part][2*i+1];
      PDM_g_num_t new_gnum  = cur_pentity_itrf_gnum[i_part][  i  ];
      int j = pn_entity_cur[i_part]+i;

      /* Remplissage */
      concat_pentity_itrf_gnum[i_part][i] = orig_gnum;
      concat_pentity_itrf_gnum[i_part][j] = new_gnum;

      concat_pentity_itrf_data[i_part][2*i  ] = new_gnum;
      concat_pentity_itrf_data[i_part][2*i+1] =-i_itrf;

      concat_pentity_itrf_data[i_part][2*j  ] = orig_gnum;
      concat_pentity_itrf_data[i_part][2*j+1] = i_itrf;

      concat_pentity_itrf_strd[i_part][i] = 1;
      concat_pentity_itrf_strd[i_part][j] = 1;

      if (has_sens==1) {
        concat_pentity_itrf_sens[i_part][i] = cur_pentity_itrf_sens[i_part][i];
        concat_pentity_itrf_sens[i_part][j] = cur_pentity_itrf_sens[i_part][i];
      }

      concat_pentity_itrf_wght[i_part][i] = 1.;
      concat_pentity_itrf_wght[i_part][j] = 1.;

    }

    if(debug == 1) {
      log_trace("\n");
      PDM_log_trace_array_long  (concat_pentity_itrf_gnum[i_part],   pn_entity[i_part], "concat_pentity_itrf_gnum ::");
      PDM_log_trace_array_int   (concat_pentity_itrf_strd[i_part],   pn_entity[i_part], "concat_pentity_itrf_strd ::");
      PDM_log_trace_array_long  (concat_pentity_itrf_data[i_part], 2*pn_entity[i_part], "concat_pentity_itrf_data ::");
      if (has_sens==1) {
        PDM_log_trace_array_int   (concat_pentity_itrf_sens[i_part],   pn_entity[i_part], "concat_pentity_itrf_sens ::");
      }
    }
  }

  concat_pentity_itrf_gnum[n_part] = prev_pentity_itrf_gnum;
  concat_pentity_itrf_strd[n_part] = prev_pentity_itrf_strd;
  concat_pentity_itrf_data[n_part] = prev_pentity_itrf_data;
  if (has_sens==1) {
    concat_pentity_itrf_sens[n_part] = prev_pentity_itrf_sens;
  }
  concat_pentity_itrf_wght[n_part] = PDM_array_const_double(pn_entity_prev, 1.);

  pn_entity[n_part] = pn_entity_prev;

  if(debug == 1) {
    int n_data = 0;
    for (int i=0; i<pn_entity_prev;++i) {
      n_data += concat_pentity_itrf_strd[n_part][i];
    }
    log_trace("\n");
    PDM_log_trace_array_long  (concat_pentity_itrf_gnum[n_part], pn_entity_prev, "concat_pentity_itrf_gnum ::");
    PDM_log_trace_array_int   (concat_pentity_itrf_strd[n_part], pn_entity_prev, "concat_pentity_itrf_strd ::");
    PDM_log_trace_array_long  (concat_pentity_itrf_data[n_part], 2*n_data      , "concat_pentity_itrf_data ::");
    if (has_sens==1) {
      PDM_log_trace_array_int   (concat_pentity_itrf_sens[n_part],   n_data      , "concat_pentity_itrf_sens ::");
    }
    // PDM_log_trace_array_double(concat_pentity_itrf_wght[n_part], pn_entity_prev, "concat_pentity_itrf_wght ::");
    PDM_log_trace_array_int   (pn_entity, n_part+1, "pn_entity ::");
  }

  /*
   * Merge informations with PTB
   */
  PDM_part_to_block_t* ptb_itrf = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           concat_pentity_itrf_gnum,
                                                           concat_pentity_itrf_wght,
                                                           pn_entity,
                                                           n_part+1,
                                                           comm);

  // > Exchange data
  int         *dentity_itrf_strd = NULL;
  PDM_g_num_t *dentity_itrf_data = NULL;
  PDM_part_to_block_exch(ptb_itrf,
                         2 * sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         concat_pentity_itrf_strd,
          (void **)      concat_pentity_itrf_data,
                        &dentity_itrf_strd,
          (void **)     &dentity_itrf_data);

  int *dentity_itrf_strd2 = NULL;
  int *dentity_itrf_sens  = NULL;
  if (has_sens==1) {
    PDM_part_to_block_exch(ptb_itrf,
                           1 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           concat_pentity_itrf_strd,
            (void **)      concat_pentity_itrf_sens,
                          &dentity_itrf_strd2,
            (void **)     &dentity_itrf_sens);
    PDM_free(dentity_itrf_strd2);
  }

  PDM_free(concat_pentity_itrf_wght[n_part]);

  // > Free local arrays
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(concat_pentity_itrf_gnum[i_part]);
    PDM_free(concat_pentity_itrf_strd[i_part]);
    PDM_free(concat_pentity_itrf_data[i_part]);
    if (has_sens==1) {
      PDM_free(concat_pentity_itrf_sens[i_part]);
    }
    PDM_free(concat_pentity_itrf_wght[i_part]);
  }
  PDM_free(concat_pentity_itrf_gnum);
  PDM_free(concat_pentity_itrf_strd);
  PDM_free(concat_pentity_itrf_data);
  PDM_free(concat_pentity_itrf_sens);
  PDM_free(concat_pentity_itrf_wght);
  PDM_free(pn_entity);


  // TODO comprendre pourquoi il répartit de manière cheloue
  /*
   * Post-treatment of PTB
   */
  int dn_entity_itrf = PDM_part_to_block_n_elt_block_get(ptb_itrf);
  PDM_g_num_t *next_dentity_itrf_gnum_ptp = PDM_part_to_block_block_gnum_get(ptb_itrf);
  

  PDM_g_num_t *dentity_itrf_gnum = NULL;
  PDM_malloc(dentity_itrf_gnum, dn_entity_itrf, PDM_g_num_t);
  for(int i = 0; i < dn_entity_itrf; ++i) {
    dentity_itrf_gnum[i] = next_dentity_itrf_gnum_ptp[i];
  }

  // > Compute max stride
  int max_strid = 0;
  for(int i = 0; i < dn_entity_itrf; ++i) {
    max_strid = PDM_MAX(max_strid, dentity_itrf_strd[i]);
  }
  int *data_order = NULL;
  PDM_malloc(data_order, max_strid, int);

  int *dentity_itrf_sens_cp = NULL; // cause of data_order
  if (has_sens==1) { 
    int n_sens = 0;
    for (int i_entity=0; i_entity<dn_entity_itrf; ++i_entity) {
      n_sens += dentity_itrf_strd[i_entity];
    }
    dentity_itrf_sens_cp = PDM_array_copy_int(dentity_itrf_sens, n_sens);
  }

  // > Unique incomming data
  int idx_read  = 0;
  int idx_write = 0;
  for(int i = 0; i < dn_entity_itrf; ++i) {

    int n_data   = dentity_itrf_strd[i];
    int n_unique = PDM_order_inplace_unique_and_order_long(n_data, 2, &dentity_itrf_data[2*idx_read], data_order);

    /* Copy */
    for(int i_unique = 0; i_unique < n_unique; ++i_unique) {
      dentity_itrf_data[2*idx_write  ] = dentity_itrf_data[2*(idx_read+i_unique)  ];
      dentity_itrf_data[2*idx_write+1] = dentity_itrf_data[2*(idx_read+i_unique)+1];
      if (has_sens==1) {
        dentity_itrf_sens[  idx_write  ] = dentity_itrf_sens_cp[ idx_read+  data_order [i_unique]   ];
      }
      idx_write++;
    }

    dentity_itrf_strd[i] = n_unique;
    idx_read += n_data;
  }
  PDM_free(data_order);
  PDM_free(dentity_itrf_sens_cp);
  

  PDM_realloc(dentity_itrf_data, dentity_itrf_data, 2 * idx_write, PDM_g_num_t);
  if (has_sens==1) {
    PDM_realloc(dentity_itrf_sens, dentity_itrf_sens,   idx_write, int        );
  }

  if(debug == 1) {
    log_trace("\n");
    log_trace("dn_entity_itrf = %d\n", dn_entity_itrf);
    PDM_log_trace_array_long(dentity_itrf_gnum,     dn_entity_itrf, "dentity_itrf_gnum (Unique) ::");
    PDM_log_trace_array_int (dentity_itrf_strd,     dn_entity_itrf, "dentity_itrf_strd (Unique) ::");
    PDM_log_trace_array_long(dentity_itrf_data, 2 * idx_write     , "dentity_itrf_data (Unique) ::");
    if (has_sens==1) {
      PDM_log_trace_array_int (dentity_itrf_sens,     idx_write     , "dentity_itrf_sens (Unique) ::");
    }
  }

  PDM_part_to_block_free(ptb_itrf);
  


  /*
   * Keep pointer for output
   */
  *dn_entity_next         = dn_entity_itrf;
  *next_dentity_itrf_gnum = dentity_itrf_gnum;
  *next_dentity_itrf_strd = dentity_itrf_strd;
  *next_dentity_itrf_data = dentity_itrf_data;
  if (has_sens==1) {
    *next_dentity_itrf_sens = dentity_itrf_sens;
  }
  if (debug==1) {
    log_trace("END _enrich_block_interface\n");
    log_trace("---------------------------\n");
    log_trace("\n");
  }
}


static
void
_enrich_block_ancestors
(
  int           n_part,
  int          *pn_entity_cur,
  PDM_g_num_t **cur_pentity_itrf_gnum,
  PDM_g_num_t **cur_pentity_itrf_ancstr,
  int         **cur_pentity_itrf_path_itrf_strd,
  int         **cur_pentity_itrf_path_itrf,
  int           dn_entity_prev,
  PDM_g_num_t  *prev_dentity_itrf_gnum,
  int          *prev_dentity_itrf_ancstr_strd,
  PDM_g_num_t  *prev_dentity_itrf_ancstr,
  int          *prev_dentity_itrf_path_itrf_strd,
  int          *prev_dentity_itrf_path_itrf,
  int           dn_entity_next,
  PDM_g_num_t  *next_dentity_itrf_gnum,
  int         **next_dentity_itrf_ancstr_strd,
  PDM_g_num_t **next_dentity_itrf_ancstr,
  int         **next_dentity_itrf_path_itrf_strd,
  int         **next_dentity_itrf_path_itrf,
  PDM_MPI_Comm  comm
)
{

  int debug = 0;
  if (debug==1) {
    log_trace("\n");
    log_trace("---------------------------\n");
    log_trace("BEG _enrich_block_ancestors\n");
  }

  /**
   * Concatenate arrays for part_to_part
   */
  int          *concat_pn_entity                   = NULL;
  PDM_g_num_t **concat_pentity_itrf_gnum           = NULL;
  int         **concat_pentity_itrf_ancstr_strd    = NULL;
  PDM_g_num_t **concat_pentity_itrf_ancstr         = NULL;
  int         **concat_pentity_itrf_len_path_strd  = NULL;
  int         **concat_pentity_itrf_path_itrf_strd = NULL;
  int         **concat_pentity_itrf_path_itrf      = NULL;
  int         **concat_pentity_to_cur_dentity_idx  = NULL;
  PDM_malloc(concat_pn_entity                  , n_part+1, int          );
  PDM_malloc(concat_pentity_itrf_gnum          , n_part+1, PDM_g_num_t *);
  PDM_malloc(concat_pentity_itrf_ancstr_strd   , n_part+1, int         *);
  PDM_malloc(concat_pentity_itrf_ancstr        , n_part+1, PDM_g_num_t *);
  PDM_malloc(concat_pentity_itrf_len_path_strd , n_part+1, int         *);
  PDM_malloc(concat_pentity_itrf_path_itrf_strd, n_part+1, int         *);
  PDM_malloc(concat_pentity_itrf_path_itrf     , n_part+1, int         *);
  PDM_malloc(concat_pentity_to_cur_dentity_idx , n_part+1, int         *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    concat_pn_entity                  [i_part] = pn_entity_cur[i_part];
    concat_pentity_itrf_gnum          [i_part] = cur_pentity_itrf_gnum[i_part];
    concat_pentity_itrf_ancstr_strd   [i_part] = PDM_array_const_int(pn_entity_cur[i_part], 1);
    concat_pentity_itrf_ancstr        [i_part] = cur_pentity_itrf_ancstr[i_part];
    concat_pentity_itrf_len_path_strd [i_part] = cur_pentity_itrf_path_itrf_strd[i_part];
    concat_pentity_itrf_path_itrf_strd[i_part] = cur_pentity_itrf_path_itrf_strd[i_part];
    concat_pentity_itrf_path_itrf     [i_part] = cur_pentity_itrf_path_itrf[i_part];
    concat_pentity_to_cur_dentity_idx [i_part] = PDM_array_new_idx_from_const_stride_int(1, pn_entity_cur[i_part]);
    
    if(debug == 1) {
      int n_ancstr = 0;
      int n_itrf   = 0;
      for (int i=0; i<concat_pn_entity[i_part];++i) {
        n_ancstr += concat_pentity_itrf_ancstr_strd   [i_part][i];
        n_itrf   += concat_pentity_itrf_path_itrf_strd[i_part][i];
      }

      PDM_log_trace_array_long(concat_pentity_itrf_gnum          [i_part], concat_pn_entity[i_part], "concat_pentity_itrf_gnum          [i_part] ::");
      PDM_log_trace_array_int (concat_pentity_itrf_ancstr_strd   [i_part], concat_pn_entity[i_part], "concat_pentity_itrf_ancstr_strd   [i_part] ::");
      PDM_log_trace_array_long(concat_pentity_itrf_ancstr        [i_part], n_ancstr                , "concat_pentity_itrf_ancstr        [i_part] ::");
      PDM_log_trace_array_int (concat_pentity_itrf_path_itrf_strd[i_part], concat_pn_entity[i_part], "concat_pentity_itrf_path_itrf_strd[i_part] ::");
      PDM_log_trace_array_int (concat_pentity_itrf_path_itrf     [i_part], n_itrf                  , "concat_pentity_itrf_path_itrf     [i_part] ::");
    }
  }
  concat_pn_entity                  [n_part] = dn_entity_prev;
  concat_pentity_itrf_gnum          [n_part] = prev_dentity_itrf_gnum;
  concat_pentity_itrf_ancstr_strd   [n_part] = prev_dentity_itrf_ancstr_strd;
  concat_pentity_itrf_ancstr        [n_part] = prev_dentity_itrf_ancstr;
  concat_pentity_itrf_len_path_strd [n_part] = prev_dentity_itrf_path_itrf_strd;
  concat_pentity_itrf_path_itrf_strd[n_part] = prev_dentity_itrf_path_itrf_strd;
  concat_pentity_itrf_path_itrf     [n_part] = prev_dentity_itrf_path_itrf;
  concat_pentity_to_cur_dentity_idx [n_part] = PDM_array_new_idx_from_const_stride_int(1, dn_entity_prev);
  if(debug == 1) {
    int n_ancstr = 0;
    int n_itrf   = 0;
    for (int i=0; i<concat_pn_entity[n_part];++i) {
      n_ancstr += concat_pentity_itrf_ancstr_strd   [n_part][i];
      n_itrf   += concat_pentity_itrf_path_itrf_strd[n_part][i];
    }

    PDM_log_trace_array_long(concat_pentity_itrf_gnum          [n_part], concat_pn_entity[n_part], "concat_pentity_itrf_gnum          [n_part] ::");
    PDM_log_trace_array_int (concat_pentity_itrf_ancstr_strd   [n_part], concat_pn_entity[n_part], "concat_pentity_itrf_ancstr_strd   [i_part] ::");
    PDM_log_trace_array_long(concat_pentity_itrf_ancstr        [n_part], n_ancstr                , "concat_pentity_itrf_ancstr        [n_part] ::");
    PDM_log_trace_array_int (concat_pentity_itrf_path_itrf_strd[n_part], concat_pn_entity[n_part], "concat_pentity_itrf_path_itrf_strd[n_part] ::");
    PDM_log_trace_array_int (concat_pentity_itrf_path_itrf     [n_part], n_itrf                  , "concat_pentity_itrf_path_itrf     [n_part] ::");
  }
  


  /**
   * part_to_part to get ancestors and path interface
   */
  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **)concat_pentity_itrf_gnum,
                                                    concat_pn_entity,
                                                    n_part+1,
                                                   (const PDM_g_num_t **) &next_dentity_itrf_gnum,
                                                   &dn_entity_next,
                                                    1,
                                                   (const int **) concat_pentity_to_cur_dentity_idx,
                                                   (const PDM_g_num_t **) concat_pentity_itrf_gnum,
                                                   comm);

  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2 = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int           request_exch_ancstr = 0;
  int         **rcvd_ancstr_strd = NULL;
  PDM_g_num_t **rcvd_ancstr = NULL;
  PDM_part_to_part_iexch( ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_VAR_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2 ,
                          1,
                          1*sizeof(PDM_g_num_t),
        (const int  **)   concat_pentity_itrf_ancstr_strd,
        (const void **)   concat_pentity_itrf_ancstr,
                         &rcvd_ancstr_strd,
              (void ***) &rcvd_ancstr,
                         &request_exch_ancstr);
  PDM_part_to_part_iexch_wait(ptp, request_exch_ancstr);

  int request_exch_path_itrf = 0;
  int **rcvd_path_itrf_strd = NULL;
  int **rcvd_path_itrf = NULL;
  PDM_part_to_part_iexch( ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_VAR_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2 ,
                          1,
                          1*sizeof(int),
        (const int  **)   concat_pentity_itrf_path_itrf_strd,
        (const void **)   concat_pentity_itrf_path_itrf,
                         &rcvd_path_itrf_strd,
              (void ***) &rcvd_path_itrf,
                         &request_exch_path_itrf);
  PDM_part_to_part_iexch_wait(ptp, request_exch_path_itrf);

  

  if (debug==1) {
    log_trace("dn_entity_next = %d\n", dn_entity_next);
    log_trace("n_ref_lnum2[0] = %d\n", n_ref_lnum2[0]);
    PDM_log_trace_array_int (ref_lnum2[0], n_ref_lnum2[0], "ref_lnum2 ::");
    int n_ancstr = 0;
    int n_itrf   = 0;
    for (int i=0; i<n_ref_lnum2[0]; ++i) {
      n_ancstr += rcvd_ancstr_strd[0][i];
      n_itrf   += rcvd_path_itrf_strd[0][i];
    }
    PDM_log_trace_array_int (rcvd_ancstr_strd   [0], n_ref_lnum2[0], "rcvd_ancstr_strd    ::");
    PDM_log_trace_array_long(next_dentity_itrf_gnum, dn_entity_next, "next_dentity_itrf_gnum ::");
    PDM_log_trace_array_int (rcvd_ancstr_strd   [0], n_ref_lnum2[0], "rcvd_ancstr_strd    ::");
    PDM_log_trace_array_long(rcvd_ancstr        [0], n_ancstr      , "rcvd_ancstr         ::");
    PDM_log_trace_array_int (rcvd_path_itrf_strd[0], n_ref_lnum2[0], "rcvd_path_itrf_strd ::");
    PDM_log_trace_array_int (rcvd_path_itrf     [0], n_itrf        , "rcvd_path_itrf      ::");
  }


  /**
   * Post treat part_to_part to get ancestors and path interface on db referencial
   */
  int         *_next_dentity_itrf_ancstr_strd    = PDM_array_zeros_int(dn_entity_next);
  PDM_g_num_t *_next_dentity_itrf_ancstr         = rcvd_ancstr[0];
  int         *_next_dentity_itrf_path_itrf_strd = PDM_array_zeros_int(dn_entity_next);
  int         *_next_dentity_itrf_path_itrf      = rcvd_path_itrf[0];

  for (int i_ref_entity=0; i_ref_entity<n_ref_lnum2[0]; ++i_ref_entity) {
    int i_entity = ref_lnum2[0][i_ref_entity]-1;
    _next_dentity_itrf_ancstr_strd   [i_entity] = rcvd_ancstr_strd   [0][i_ref_entity];
    _next_dentity_itrf_path_itrf_strd[i_entity] = rcvd_path_itrf_strd[0][i_ref_entity];
  }

  *next_dentity_itrf_ancstr_strd    = _next_dentity_itrf_ancstr_strd;
  *next_dentity_itrf_ancstr         = _next_dentity_itrf_ancstr;
  *next_dentity_itrf_path_itrf_strd = _next_dentity_itrf_path_itrf_strd;
  *next_dentity_itrf_path_itrf      = _next_dentity_itrf_path_itrf;


  if (debug==1) {
    log_trace("\n");
    log_trace("END _enrich_block_ancestors\n");
    log_trace("---------------------------\n");
  }




  /**
   * Free
   */
  PDM_free(rcvd_ancstr_strd[0]);
  PDM_free(rcvd_path_itrf_strd[0]);
  PDM_free(rcvd_ancstr_strd);
  PDM_free(rcvd_ancstr);
  PDM_free(rcvd_path_itrf_strd);
  PDM_free(rcvd_path_itrf);
  PDM_part_to_part_free(ptp);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(concat_pentity_itrf_ancstr_strd  [i_part]);
    PDM_free(concat_pentity_to_cur_dentity_idx[i_part]);
  }
  PDM_free(concat_pentity_to_cur_dentity_idx[n_part]);

  PDM_free(concat_pn_entity);
  PDM_free(concat_pentity_itrf_gnum);
  PDM_free(concat_pentity_itrf_ancstr_strd);
  PDM_free(concat_pentity_itrf_ancstr);
  PDM_free(concat_pentity_itrf_len_path_strd);
  PDM_free(concat_pentity_itrf_path_itrf_strd);
  PDM_free(concat_pentity_itrf_path_itrf);
  PDM_free(concat_pentity_to_cur_dentity_idx);
}


static
void
_compute_gnum_from_ancestor_and_itrfs
(
  int            n_part,
  PDM_g_num_t    shift,
  int           *pn_entity1,
  int          **pentity1_entity2_idx,
  int          **pentity1_entity2_kind,
  PDM_g_num_t  **pentity1_entity2_gnum,
  int          **pentity1_entity2_itrf,
  int           *pentity1_entity2_ntot,
  PDM_g_num_t  **pentity1_entity2_ancstr,
  int          **pentity1_entity2_path_itrf_strd,
  int          **pentity1_entity2_path_itrf,
  int          **pentity1_entity2_kind_idx,
  int          **pentity1_entity2_kind_order,
  PDM_g_num_t ***pentity1_entity2_gnum_out,
  PDM_g_num_t ***pnew_entity2_gnum_out,
  PDM_g_num_t ***pnew_entity2_ancstr_out,
  int         ***pnew_entity2_path_itrf_strd_out,
  int         ***pnew_entity2_path_itrf_out,
  PDM_g_num_t ***pnew_entity2_parent_nuplet_out,
  PDM_MPI_Comm   comm
)
{
  /*
   * TODO:
   *   - compute gnum from nuplet
   *      (cause not sure that actual version is robust:
   *       can new entities have different number of previous interface ?)
   *   - add shift on new entities
   */

  int debug      = 0;
  int debug_loop = 0;

  if (debug==1) {
    log_trace("\n");
    log_trace("-----------------------------------------\n");
    log_trace("BEG _compute_gnum_from_ancestor_and_itrfs\n");
    log_trace("\n");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_part1_to_part2 = pentity1_entity2_idx[i_part][pn_entity1[i_part]];
      int n_data = 0;
      for(int i_entity = 0; i_entity < n_part1_to_part2; ++i_entity) {
        n_data +=pentity1_entity2_path_itrf_strd[i_part][i_entity];
      }
      PDM_log_trace_array_long(pentity1_entity2_gnum[i_part], n_part1_to_part2, "pentity1_entity2_gnum ::");
      PDM_log_trace_array_int (pentity1_entity2_itrf[i_part], n_part1_to_part2, "pentity1_entity2_itrf ::");
      PDM_log_trace_array_int (pentity1_entity2_kind[i_part], n_part1_to_part2, "pentity1_entity2_kind ::");
      PDM_log_trace_array_long(pentity1_entity2_ancstr[i_part], n_part1_to_part2, "pentity1_entity2_ancstr ::");
      PDM_log_trace_array_int (pentity1_entity2_path_itrf_strd[i_part], n_part1_to_part2, "pentity1_entity2_path_itrf_strd ::");
      PDM_log_trace_array_int (pentity1_entity2_path_itrf     [i_part], n_data, "pentity1_entity2_path_itrf ::");
    }
  }


  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3,
                                                     n_part,
                                                     PDM_TRUE,
                                                     1.e-6,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);


  /*
   * Count the max number of ancestors for entities generating a new gnum
   */
  int l_max_len_path = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    for (int i_entity=0; i_entity<pentity1_entity2_idx[i_part][pn_entity1[i_part]]; ++i_entity) {
      if (pentity1_entity2_kind[i_part][i_entity] == 5) {
        int len_path = pentity1_entity2_path_itrf_strd[i_part][i_entity];
        l_max_len_path = PDM_MAX(len_path, l_max_len_path);
      }
    }
  }

  // > Reduce
  int   max_len_path_t[1];
  int l_max_len_path_t[1] = {l_max_len_path};
  if (debug==1) {log_trace("\n-->l_max_len_path = %d\n", l_max_len_path);}
  PDM_MPI_Allreduce(l_max_len_path_t, max_len_path_t, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  int   max_len_path = max_len_path_t[0];
  if (debug==1) {log_trace("\n-->max_len_path = %d\n", max_len_path);}
  
  PDM_gnum_set_parents_nuplet(gen_gnum_entity2, 2+max_len_path);


  /*
   * Build nuplets for gnum computation
   */
  PDM_g_num_t **pentity2_ancstr_nuplet = NULL;
  PDM_g_num_t **pentity2_parent_nuplet = NULL;
  PDM_malloc(pentity2_ancstr_nuplet, n_part, PDM_g_num_t *);
  PDM_malloc(pentity2_parent_nuplet, n_part, PDM_g_num_t *);
  
  int **pentity1_entity2_path_itrf_idx = NULL;
  PDM_malloc(pentity1_entity2_path_itrf_idx, n_part, int *);
 
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pentity1_entity2_path_itrf_idx[i_part] = PDM_array_new_idx_from_sizes_int(pentity1_entity2_path_itrf_strd[i_part], pentity1_entity2_ntot[i_part]);

    PDM_g_num_t *tmp_ancstr_nplt = NULL;
    PDM_malloc(tmp_ancstr_nplt, max_len_path+1, PDM_g_num_t);
    int  n_entity2_kind4 = pentity1_entity2_kind_idx[i_part][6]-pentity1_entity2_kind_idx[i_part][5];
    
    PDM_malloc(pentity2_ancstr_nuplet[i_part], (2+max_len_path)* n_entity2_kind4, PDM_g_num_t);
    PDM_malloc(pentity2_parent_nuplet[i_part],  2              * n_entity2_kind4, PDM_g_num_t);
    
    int idx_write = 0;
    for(int i = pentity1_entity2_kind_idx[i_part][5]; i < pentity1_entity2_kind_idx[i_part][6]; ++i) {
      int j = pentity1_entity2_kind_order[i_part][i];
      if (debug_loop==1) {
        log_trace("j = %d :: gnum = "PDM_FMT_G_NUM" ; itrf = %d\n", j, pentity1_entity2_gnum[i_part][j], pentity1_entity2_itrf[i_part][j]);
        log_trace("idx = %d\n", pentity1_entity2_path_itrf_idx[i_part][j]);
      }
      pentity2_parent_nuplet[i_part][2*idx_write  ] = pentity1_entity2_gnum[i_part][j];
      pentity2_parent_nuplet[i_part][2*idx_write+1] = pentity1_entity2_itrf[i_part][j]; // TODO: why this "-" ??? not symmetric with vtx why ??

      // > Fill with previous interface, those who has no previous or less than other, fill with current interface
      int l_i_ancstr = 0;
      if (debug_loop==1) log_trace("beg l_i_ancstr = %d \n", l_i_ancstr);
      for (int i_ancstr=pentity1_entity2_path_itrf_idx[i_part][j]; i_ancstr<pentity1_entity2_path_itrf_idx[i_part][j+1]; ++i_ancstr) {
        tmp_ancstr_nplt[l_i_ancstr++] = pentity1_entity2_path_itrf[i_part][i_ancstr];
      }
      if (debug_loop==1) log_trace("    l_i_ancstr = %d \n", l_i_ancstr);
      int len_path = pentity1_entity2_path_itrf_idx[i_part][j+1]-pentity1_entity2_path_itrf_idx[i_part][j];
      tmp_ancstr_nplt[l_i_ancstr++] = pentity1_entity2_itrf[i_part][j];
      if (debug_loop==1) log_trace("    l_i_ancstr = %d \n", l_i_ancstr);
      for (int i_ancstr=len_path+1; i_ancstr<max_len_path+1; ++i_ancstr){
        tmp_ancstr_nplt[l_i_ancstr++] = 0; // works ?? --> if not really need a compute gnum from idx
      }
      if (debug_loop==1) log_trace("end l_i_ancstr = %d\n", l_i_ancstr);

      // > Sort interface so that entity with same ancstr will give same gnum
      PDM_sort_long(tmp_ancstr_nplt, NULL, max_len_path+1);
      if (debug_loop==1) {
        PDM_log_trace_array_long(tmp_ancstr_nplt, max_len_path+1, "----> nuplet ::");
      }
      
      // log_trace("Writing case =  ");
      // log_trace(" %d ", (max_len_path+2)*idx_write);
      if (debug_loop==1) log_trace("\t nuplet =  ");
      if (debug_loop==1) log_trace(" %d ", pentity1_entity2_ancstr[i_part][j]);
      pentity2_ancstr_nuplet[i_part][(max_len_path+2)*idx_write  ] = pentity1_entity2_ancstr[i_part][j];
      for (int i_ancstr=0; i_ancstr<max_len_path+1; ++i_ancstr){
        // log_trace(" %d ", (max_len_path+2)*idx_write+1+i_ancstr);
        if (debug_loop==1) log_trace(" %d ", tmp_ancstr_nplt[i_ancstr]);
        pentity2_ancstr_nuplet[i_part][(max_len_path+2)*idx_write+1+i_ancstr] = tmp_ancstr_nplt[i_ancstr];
      }
      if (debug_loop==1) log_trace("\n");
      idx_write++;
    }
    assert(idx_write == n_entity2_kind4);
    PDM_free(tmp_ancstr_nplt);


    if(debug == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(pentity2_ancstr_nuplet[i_part], (max_len_path+2)* n_entity2_kind4, "pentity2_ancstr_nuplet ::");
      PDM_log_trace_array_long(pentity2_parent_nuplet[i_part],               2 * n_entity2_kind4, "pentity2_parent_nuplet ::");
    }


    PDM_gnum_set_from_parents(gen_gnum_entity2,
                              i_part,
                              n_entity2_kind4,
                              pentity2_ancstr_nuplet[i_part]);
  }

  PDM_gnum_compute(gen_gnum_entity2);


  /*
   * Replace new kind4 entities gnum
   */
  PDM_g_num_t **pnew_entity2_gnum           = NULL;
  PDM_g_num_t **pnew_entity2_ancstr         = NULL;
  int         **pnew_entity2_path_itrf_strd = NULL;
  int         **pnew_entity2_path_itrf      = NULL;
  PDM_malloc(pnew_entity2_gnum          , n_part, PDM_g_num_t *);
  PDM_malloc(pnew_entity2_ancstr        , n_part, PDM_g_num_t *);
  PDM_malloc(pnew_entity2_path_itrf_strd, n_part, int         *);
  PDM_malloc(pnew_entity2_path_itrf     , n_part, int         *);

  for(int i_part = 0; i_part < n_part; ++i_part) {


    int  n_entity2_kind4 = pentity1_entity2_kind_idx[i_part][6]-pentity1_entity2_kind_idx[i_part][5];

    int len_path_itrf_tot = pentity1_entity2_path_itrf_idx[i_part][pentity1_entity2_ntot[i_part]] + n_entity2_kind4;
    if (debug_loop==1) log_trace("len_path_itrf_tot = %d \n", pentity1_entity2_path_itrf_idx[i_part][n_entity2_kind4]);
    if (debug_loop==1) log_trace("len_path_itrf_tot = %d \n", n_entity2_kind4);
    if (debug_loop==1) log_trace("len_path_itrf_tot = %d \n", len_path_itrf_tot);
    PDM_malloc(pnew_entity2_gnum          [i_part], n_entity2_kind4  , PDM_g_num_t);
    PDM_malloc(pnew_entity2_ancstr        [i_part], n_entity2_kind4  , PDM_g_num_t);
    PDM_malloc(pnew_entity2_path_itrf_strd[i_part], n_entity2_kind4  , int        );
    PDM_malloc(pnew_entity2_path_itrf     [i_part], len_path_itrf_tot, int        );

    PDM_g_num_t* kind4_entity2_new_gnum = PDM_gnum_get(gen_gnum_entity2, i_part);

    /* Update array */
    int idx_write = 0;
    int idx_write_data = 0;
    int idx_write_dbg = 0;
    for(int i = pentity1_entity2_kind_idx[i_part][5]; i < pentity1_entity2_kind_idx[i_part][6]; ++i) {
      int idx_read = pentity1_entity2_kind_order[i_part][i];
      if (debug_loop==1) {
        log_trace("new gnum "PDM_FMT_G_NUM" (path = ", kind4_entity2_new_gnum[idx_write] + shift);
        for (int i_path=0; i_path<2+max_len_path; ++i_path) {
          log_trace("%d ", pentity2_ancstr_nuplet[i_part][idx_write_dbg++]);
        }
        log_trace(")\n");
      }
      pentity1_entity2_gnum[i_part][idx_read]  = kind4_entity2_new_gnum[idx_write] + shift;
      pnew_entity2_gnum    [i_part][idx_write] = pentity1_entity2_gnum  [i_part][idx_read];
      
      pnew_entity2_ancstr  [i_part][idx_write] = pentity1_entity2_ancstr[i_part][idx_read];
      pnew_entity2_path_itrf_strd[i_part][idx_write] = pentity1_entity2_path_itrf_strd[i_part][idx_read]+1;
      int i_beg_read_data = pentity1_entity2_path_itrf_idx[i_part][idx_read  ];
      int i_end_read_data = pentity1_entity2_path_itrf_idx[i_part][idx_read+1];
      for (int i_read_data=i_beg_read_data; i_read_data<i_end_read_data; ++i_read_data) {
        pnew_entity2_path_itrf[i_part][idx_write_data] = pentity1_entity2_path_itrf[i_part][i_read_data];
        idx_write_data++;
      }
      pnew_entity2_path_itrf[i_part][idx_write_data++] = pentity1_entity2_itrf[i_part][idx_read];
      idx_write++;
    }

    if (debug==1) {
      int n_itrf = 0;
      for (int i=0; i<n_entity2_kind4;++i) {
        n_itrf += pnew_entity2_path_itrf_strd[i_part][i];
      }
      log_trace("\n");
      PDM_log_trace_array_long(pentity1_entity2_gnum [i_part],   pentity1_entity2_ntot[i_part], "pentity1_entity2_gnum (new) ::");
      PDM_log_trace_array_long(pnew_entity2_gnum     [i_part],   n_entity2_kind4              , "pnew_entity2_gnum           ::");
      PDM_log_trace_array_long(pnew_entity2_ancstr   [i_part],   n_entity2_kind4              , "pnew_entity2_ancstr         ::");
      PDM_log_trace_array_int (pnew_entity2_path_itrf_strd[i_part],   n_entity2_kind4              , "pnew_entity2_path_itrf_strd::");
      PDM_log_trace_array_int (pnew_entity2_path_itrf[i_part],   n_itrf              , "pnew_entity2_path_itrf::");
      PDM_log_trace_array_long(pentity2_parent_nuplet[i_part], 2*n_entity2_kind4              , "pentity2_parent_nuplet      ::");
    }
  }

  *pentity1_entity2_gnum_out       = pentity1_entity2_gnum;
  *pnew_entity2_gnum_out           = pnew_entity2_gnum;
  *pnew_entity2_ancstr_out         = pnew_entity2_ancstr;
  *pnew_entity2_path_itrf_strd_out = pnew_entity2_path_itrf_strd;
  *pnew_entity2_path_itrf_out      = pnew_entity2_path_itrf;
  *pnew_entity2_parent_nuplet_out  = pentity2_parent_nuplet;

  if (debug==1) {
    log_trace("END _compute_gnum_from_ancestor_and_itrfs\n");
    log_trace("-----------------------------------------\n");
    log_trace("\n");
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity2_ancstr_nuplet[i_part]);
    PDM_free(pentity1_entity2_path_itrf_idx[i_part]);
  }
  PDM_free(pentity2_ancstr_nuplet);
  PDM_free(pentity1_entity2_path_itrf_idx);
  PDM_gnum_free(gen_gnum_entity2);
}


static void 
_unique_entities_and_update_connectivity
(
  int            n_part,
  int           *pn_entity,
  int          **pentity_sign,
  int          **pentity_lnum,
  int          **pentity_itrf,
  int          **pentity_sens,
  int          **pentity_triplet,
  PDM_g_num_t  **pentity_gnum,
  int          **pentity_kind,
  int          **pentity_kind_idx,
  int          **pentity_kind_order,
  int          **out_pn_entity,
  int         ***out_pconnec_to_entity,
  PDM_g_num_t ***out_pentity_gnum,
  int         ***out_pentity_alrdy_sent,
  int         ***out_pentity_to_entity_idx,
  int         ***out_pentity_to_entity_trplt,
  int         ***out_pentity_to_entity_itrf,
  int         ***out_pentity_to_entity_sens,
  int            keep_all_parent
)
{
  int debug = 0;

  int has_sens = 0;
  if (pentity_sens!=NULL) {
    has_sens = 1;
  }

  /**
   * Prepare output
   */
  int          *_out_pn_entity               = NULL;
  PDM_g_num_t **_out_pentity_gnum            = NULL;
  int         **_out_pentity_alrdy_sent      = NULL;
  int         **_out_pentity_to_entity_idx   = NULL;
  int         **_out_pentity_to_entity_trplt = NULL;
  int         **_out_pentity_to_entity_itrf  = NULL;
  int         **_out_pentity_to_entity_sens  = NULL;
  int         **_out_pconnec_to_entity       = NULL;
  PDM_malloc(_out_pn_entity              , n_part, int          );
  PDM_malloc(_out_pentity_gnum           , n_part, PDM_g_num_t *);
  PDM_malloc(_out_pentity_alrdy_sent     , n_part, int         *);
  PDM_malloc(_out_pentity_to_entity_idx  , n_part, int         *);
  PDM_malloc(_out_pentity_to_entity_trplt, n_part, int         *);
  PDM_malloc(_out_pentity_to_entity_itrf , n_part, int         *);
  if (has_sens==1){
    PDM_malloc(_out_pentity_to_entity_sens, n_part, int *);
  }
  if (out_pconnec_to_entity!=NULL) {
    _out_pconnec_to_entity = *out_pconnec_to_entity;
  }
  
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_entity = pentity_kind_idx[i_part][6];
    
    if (debug==1) {
      log_trace("\n");
      PDM_log_trace_array_long(pentity_gnum    [i_part],   n_entity, "pentity_gnum   ::");
      PDM_log_trace_array_int (pentity_triplet [i_part], 3*n_entity, "pentity_trplt  ::");
    }

    /**
     * Order gnum for searching
     */
    int beg_entry = pentity_kind_idx[i_part][2];
    int end_entry = pentity_kind_idx[i_part][6];
    int n_entry   = end_entry - beg_entry;

    int i_write = 0;
    int         *order    = NULL;
    PDM_g_num_t *tmp_gnum = NULL;
    PDM_malloc(order   , n_entry, int        );
    PDM_malloc(tmp_gnum, n_entry, PDM_g_num_t);

    for(int i = beg_entry; i < end_entry; ++i) {
      int idx_read = pentity_kind_order[i_part][i];
      tmp_gnum[i_write++] = pentity_gnum[i_part][idx_read];
    }

    int n_unique = PDM_inplace_unique_long2(tmp_gnum,
                                            order,
                                            0,
                                            n_entry-1);
    _out_pn_entity[i_part] = n_unique;

    if (debug==1) {
      log_trace("\n");
      PDM_log_trace_array_int (order   , n_entry , "order               ::");
      PDM_log_trace_array_long(tmp_gnum, n_unique, "tmp_gnum (n_unique) ::");
    }


    /**
     * Count :
     *   - how many triplet have created entity
     *   - new len_path for each new entity in partition
     */
    int  pn_entity_to_entity_tot = 0;
    int *pentity_to_entity_n = PDM_array_zeros_int(_out_pn_entity[i_part]);
    for(int i = beg_entry; i < end_entry; ++i) {
      int idx_read = pentity_kind_order[i_part][i];
      int pos = PDM_binary_search_long(pentity_gnum[i_part][idx_read], tmp_gnum, n_unique);
      
      if (keep_all_parent==1 || (keep_all_parent==0 && pentity_to_entity_n[pos]==0)) {        
        pentity_to_entity_n[pos]++;
        pn_entity_to_entity_tot++;
      }

    }
    _out_pentity_to_entity_idx[i_part] = PDM_array_new_idx_from_sizes_int(pentity_to_entity_n, _out_pn_entity[i_part]);
        
    if (debug==1) {
      log_trace("\n");
      log_trace("pn_entity_to_entity_tot = %d\n",pn_entity_to_entity_tot);
      PDM_log_trace_array_int( pentity_to_entity_n              , _out_pn_entity[i_part]  , "p_entity_to_entity_n   ::");
      PDM_log_trace_array_int(_out_pentity_to_entity_idx[i_part], _out_pn_entity[i_part]+1, "p_entity_to_entity_idx ::");
      log_trace("\n");
    }


    /**
     * Allocate array 
     */
    PDM_malloc(_out_pentity_gnum            [i_part],   _out_pn_entity[i_part] , PDM_g_num_t);
    PDM_malloc(_out_pentity_alrdy_sent      [i_part],   _out_pn_entity[i_part] , int        );
    PDM_malloc(_out_pentity_to_entity_trplt [i_part], 3*pn_entity_to_entity_tot, int        );
    PDM_malloc(_out_pentity_to_entity_itrf  [i_part],   pn_entity_to_entity_tot, int        );
    if (has_sens==1){
      PDM_malloc(_out_pentity_to_entity_sens[i_part],   pn_entity_to_entity_tot, int        );
    }
    int *pentity_to_entity_i_write = PDM_array_zeros_int(_out_pn_entity[i_part]);


    /**
     * Update connectivity for entities of kind 0 or 1
     */
    if (out_pconnec_to_entity!=NULL) {
      for(int i_entity = pentity_kind_idx[i_part][0]; i_entity < pentity_kind_idx[i_part][2]; ++i_entity) {
        int idx_read = pentity_kind_order[i_part][i_entity];
        int sign     = pentity_sign      [i_part][idx_read];
        int sens = 1;
        if (has_sens==1) {
          sens = pentity_sens[i_part][idx_read];
        }
        _out_pconnec_to_entity[i_part][idx_read] = sign*sens*pentity_lnum[i_part][idx_read];
      }
    }

    /**
     * Unique and update connectivity for entities of kind 2 to 5
     */
    for(int i = 0; i < n_entry; ++i) {
      int i_unique = order[i];
      
      int i_beg = pentity_kind_idx  [i_part][2];
      int i_pos = pentity_kind_order[i_part][i_beg+i];
 
      _out_pentity_gnum[i_part][i_unique] = tmp_gnum[i_unique];
      if (pentity_kind[i_part][i_pos]==2){
        _out_pentity_alrdy_sent[i_part][i_unique] = 1;
      } else {
        _out_pentity_alrdy_sent[i_part][i_unique] = 0;
      }

      // > Entity2 to entity2 triplet
      int sens = 1;
      if (has_sens==1) {
        sens = pentity_sens[i_part][i_pos];
      }
      i_write = _out_pentity_to_entity_idx[i_part][i_unique]+pentity_to_entity_i_write[i_unique];
      _out_pentity_to_entity_itrf [i_part][  i_write  ] = pentity_itrf   [i_part][   i_pos   ];
      _out_pentity_to_entity_trplt[i_part][3*i_write  ] = pentity_triplet[i_part][3*(i_pos)  ];
      _out_pentity_to_entity_trplt[i_part][3*i_write+1] = pentity_triplet[i_part][3*(i_pos)+1];
      _out_pentity_to_entity_trplt[i_part][3*i_write+2] = pentity_triplet[i_part][3*(i_pos)+2];
      if (has_sens==1){
        _out_pentity_to_entity_sens [i_part][  i_write  ] = sens;
      }
      if (keep_all_parent==1) {
        pentity_to_entity_i_write[i_unique]++;
      }
      
      if (out_pconnec_to_entity!=NULL) {
        int sign = pentity_sign[i_part][i_pos];
        _out_pconnec_to_entity[i_part][i_pos] = sign * sens * (pn_entity[i_part] + i_unique +1);
      }
    }

    PDM_free(pentity_to_entity_i_write);
    PDM_free(pentity_to_entity_n);
    PDM_free(tmp_gnum);
    PDM_free(order);


    if (debug==1) {
      log_trace("\n");
      log_trace("pn_entity[i_part] = %d\n", _out_pn_entity[i_part]);
      PDM_log_trace_array_long(_out_pentity_gnum           [i_part],  _out_pn_entity[i_part]    , "pentity_gnum            ::");
      log_trace("\n");
      PDM_log_trace_array_int (_out_pentity_alrdy_sent     [i_part],  _out_pn_entity[i_part]    , "pentity_alrdy_sent      ::");
      log_trace("\n");
      PDM_log_trace_array_int (_out_pentity_to_entity_idx  [i_part], (_out_pn_entity[i_part]+1) , "pentity_to_entity_idx   ::");
      PDM_log_trace_array_int (_out_pentity_to_entity_trplt[i_part], 3*pn_entity_to_entity_tot  , "pentity_to_entity_trplt ::");
      PDM_log_trace_array_int (_out_pentity_to_entity_itrf [i_part],   pn_entity_to_entity_tot  , "pentity_to_entity_itrf  ::");
      if (out_pconnec_to_entity!=NULL) {
        PDM_log_trace_array_int(_out_pconnec_to_entity[i_part], n_entity , "pconnec_to_entity  ::");
      }
    }


    /**
     * Unique quadruplet for unique parent if keep_all_parent active
     */
    if (keep_all_parent==1) {
      if (has_sens==1) {
        int *tmp_array = NULL;
        PDM_malloc(tmp_array, 5*pn_entity_to_entity_tot, int);
        for(int i = 0; i < pn_entity_to_entity_tot; ++i) {
          tmp_array[5*i  ] = _out_pentity_to_entity_trplt[i_part][3*i  ];
          tmp_array[5*i+1] = _out_pentity_to_entity_trplt[i_part][3*i+1];
          tmp_array[5*i+2] = _out_pentity_to_entity_trplt[i_part][3*i+2];
          tmp_array[5*i+3] = _out_pentity_to_entity_itrf [i_part][  i  ];
          tmp_array[5*i+4] = _out_pentity_to_entity_sens [i_part][  i  ];
        }

        int  *unique_entity_to_entity_idx = NULL;
        int  *unique_entity_to_entity_n   = NULL;
        int  *unique_entity_to_entity     = NULL;
        _unique_quintuplet(_out_pn_entity[i_part],
                           _out_pentity_to_entity_idx[i_part],
                           tmp_array,
                          &unique_entity_to_entity_idx,
                          &unique_entity_to_entity_n,
                          &unique_entity_to_entity);


        int n_unique_trplt_tot = unique_entity_to_entity_idx[_out_pn_entity[i_part]];
        if (debug==1) {
          log_trace("\n");
          log_trace("n_unique_trplt_tot = %d\n", n_unique_trplt_tot);
          PDM_log_trace_array_int (unique_entity_to_entity_idx, _out_pn_entity[i_part]+1, "unique_entity_to_entity_idx ::");
        }
        PDM_realloc(_out_pentity_to_entity_trplt[i_part], _out_pentity_to_entity_trplt[i_part], 3*n_unique_trplt_tot, int);
        PDM_realloc(_out_pentity_to_entity_itrf [i_part], _out_pentity_to_entity_itrf [i_part],   n_unique_trplt_tot, int);
        for(int i = 0; i < n_unique_trplt_tot; ++i) {
          _out_pentity_to_entity_trplt[i_part][3*i  ] = unique_entity_to_entity[5*i  ];
          _out_pentity_to_entity_trplt[i_part][3*i+1] = unique_entity_to_entity[5*i+1];
          _out_pentity_to_entity_trplt[i_part][3*i+2] = unique_entity_to_entity[5*i+2];
          _out_pentity_to_entity_itrf [i_part][  i  ] = unique_entity_to_entity[5*i+3];
          _out_pentity_to_entity_sens [i_part][  i  ] = unique_entity_to_entity[5*i+4];
        }

        PDM_free(_out_pentity_to_entity_idx[i_part]);
        _out_pentity_to_entity_idx[i_part] = unique_entity_to_entity_idx;
        PDM_free(unique_entity_to_entity_n);
        PDM_free(unique_entity_to_entity);
        PDM_free(tmp_array);
      
        if (debug==1) {
          PDM_log_trace_array_int(_out_pentity_to_entity_idx  [i_part], (_out_pn_entity[i_part]+1), "pentity_to_entity_idx   (unique)::");
          PDM_log_trace_array_int(_out_pentity_to_entity_trplt[i_part], 3*n_unique_trplt_tot      , "pentity_to_entity_trplt (unique)::");
          PDM_log_trace_array_int(_out_pentity_to_entity_itrf [i_part],   n_unique_trplt_tot      , "pentity_to_entity_itrf  (unique)::");
          PDM_log_trace_array_int(_out_pentity_to_entity_sens [i_part],   n_unique_trplt_tot      , "pentity_to_entity_sens  (unique)::");
        }

      }
      else {
        int *tmp_array = NULL;
        PDM_malloc(tmp_array, 4*pn_entity_to_entity_tot, int);
        for(int i = 0; i < pn_entity_to_entity_tot; ++i) {
          tmp_array[4*i  ] = _out_pentity_to_entity_trplt[i_part][3*i  ];
          tmp_array[4*i+1] = _out_pentity_to_entity_trplt[i_part][3*i+1];
          tmp_array[4*i+2] = _out_pentity_to_entity_trplt[i_part][3*i+2];
          tmp_array[4*i+3] = _out_pentity_to_entity_itrf [i_part][  i  ];
        }

        int  *unique_entity_to_entity_idx = NULL;
        int  *unique_entity_to_entity_n   = NULL;
        int  *unique_entity_to_entity     = NULL;
        _unique_quadruplet(_out_pn_entity[i_part],
                           _out_pentity_to_entity_idx[i_part],
                           tmp_array,
                          &unique_entity_to_entity_idx,
                          &unique_entity_to_entity_n,
                          &unique_entity_to_entity);


        int n_unique_trplt_tot = unique_entity_to_entity_idx[_out_pn_entity[i_part]];
        if (debug==1) {
          log_trace("\n");
          log_trace("n_unique_trplt_tot = %d\n", n_unique_trplt_tot);
          PDM_log_trace_array_int (unique_entity_to_entity_idx, _out_pn_entity[i_part]+1, "unique_entity_to_entity_idx ::");
        }
        PDM_realloc(_out_pentity_to_entity_trplt[i_part], _out_pentity_to_entity_trplt[i_part], 3*n_unique_trplt_tot, int);
        PDM_realloc(_out_pentity_to_entity_itrf [i_part], _out_pentity_to_entity_itrf [i_part],   n_unique_trplt_tot, int);
        for(int i = 0; i < n_unique_trplt_tot; ++i) {
          _out_pentity_to_entity_trplt[i_part][3*i  ] = unique_entity_to_entity[4*i  ];
          _out_pentity_to_entity_trplt[i_part][3*i+1] = unique_entity_to_entity[4*i+1];
          _out_pentity_to_entity_trplt[i_part][3*i+2] = unique_entity_to_entity[4*i+2];
          _out_pentity_to_entity_itrf [i_part][  i  ] = unique_entity_to_entity[4*i+3];
        }

        PDM_free(_out_pentity_to_entity_idx[i_part]);
        _out_pentity_to_entity_idx[i_part] = unique_entity_to_entity_idx;
        PDM_free(unique_entity_to_entity_n);
        PDM_free(unique_entity_to_entity);
        PDM_free(tmp_array);
      
        if (debug==1) {
          PDM_log_trace_array_int(_out_pentity_to_entity_idx  [i_part], (_out_pn_entity[i_part]+1), "pentity_to_entity_idx   (unique)::");
          PDM_log_trace_array_int(_out_pentity_to_entity_trplt[i_part], 3*n_unique_trplt_tot      , "pentity_to_entity_trplt (unique)::");
          PDM_log_trace_array_int(_out_pentity_to_entity_itrf [i_part],   n_unique_trplt_tot      , "pentity_to_entity_itrf  (unique)::");
        }

      }
    }

    for(int i = 0; i < n_unique+1; ++i) {
      _out_pentity_to_entity_idx[i_part][i]*=3;
    }
  }


  /**
   * Output
   */
  *out_pn_entity               = _out_pn_entity;
  *out_pentity_gnum            = _out_pentity_gnum;
  *out_pentity_alrdy_sent      = _out_pentity_alrdy_sent;
  *out_pentity_to_entity_idx   = _out_pentity_to_entity_idx;
  *out_pentity_to_entity_trplt = _out_pentity_to_entity_trplt;
  *out_pentity_to_entity_itrf  = _out_pentity_to_entity_itrf;
  if (has_sens==1) {
    *out_pentity_to_entity_sens  = _out_pentity_to_entity_sens;
  }
  if (out_pconnec_to_entity!=NULL) {
    *out_pconnec_to_entity       = _out_pconnec_to_entity;
  }
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
  int         **bentity1_entity2_n               = NULL;
  PDM_g_num_t **bentity1_entity2_gnum            = NULL;
  int         **bentity1_entity2_triplet         = NULL;
  int         **bentity1_entity2_interface_tot_n = NULL;
  int         **bentity1_entity2_interface       = NULL;
  int         **bentity1_entity2_interface_n     = NULL;
  PDM_malloc(bentity1_entity2_n              , n_part_tot, int         *);
  PDM_malloc(bentity1_entity2_gnum           , n_part_tot, PDM_g_num_t *);
  PDM_malloc(bentity1_entity2_triplet        , n_part_tot, int         *);
  PDM_malloc(bentity1_entity2_interface_tot_n, n_part_tot, int         *);
  PDM_malloc(bentity1_entity2_interface      , n_part_tot, int         *);
  PDM_malloc(bentity1_entity2_interface_n    , n_part_tot, int         *);

  /* Init */
  for(int i_part = 0; i_part < n_part_tot; ++i_part) {

    PDM_malloc(bentity1_entity2_n              [i_part], pn_entity1[i_part], int);
    PDM_malloc(bentity1_entity2_interface_tot_n[i_part], pn_entity1[i_part], int);

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

    PDM_malloc(bentity1_entity2_gnum       [i_part],     n_entity1_entity2     , PDM_g_num_t);
    PDM_malloc(bentity1_entity2_triplet    [i_part], 3 * n_entity1_entity2     , int        );
    PDM_malloc(bentity1_entity2_interface_n[i_part],     n_entity1_entity2     , int        );
    PDM_malloc(bentity1_entity2_interface  [i_part],     n_entity1_entity2_itrf, int        );
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
  int         **prev_concat_bentity1_entity2_n               = NULL;
  PDM_g_num_t **prev_concat_bentity1_entity2_gnum            = NULL;
  int         **prev_concat_bentity1_entity2_triplet         = NULL;
  int         **prev_concat_bentity1_entity2_interface_n     = NULL;
  int         **prev_concat_bentity1_entity2_interface_tot_n = NULL;
  int         **prev_concat_bentity1_entity2_interface       = NULL;
  PDM_malloc(prev_concat_bentity1_entity2_n              , n_part_tot, int         *);
  PDM_malloc(prev_concat_bentity1_entity2_gnum           , n_part_tot, PDM_g_num_t *);
  PDM_malloc(prev_concat_bentity1_entity2_triplet        , n_part_tot, int         *);
  PDM_malloc(prev_concat_bentity1_entity2_interface_n    , n_part_tot, int         *);
  PDM_malloc(prev_concat_bentity1_entity2_interface_tot_n, n_part_tot, int         *);
  PDM_malloc(prev_concat_bentity1_entity2_interface      , n_part_tot, int         *);

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    prev_concat_bentity1_entity2_n              [i_part] = bentity1_entity2_n              [i_part];
    prev_concat_bentity1_entity2_gnum           [i_part] = bentity1_entity2_gnum           [i_part];
    prev_concat_bentity1_entity2_triplet        [i_part] = bentity1_entity2_triplet        [i_part];
    prev_concat_bentity1_entity2_interface_n    [i_part] = bentity1_entity2_interface_n    [i_part];
    prev_concat_bentity1_entity2_interface_tot_n[i_part] = bentity1_entity2_interface_tot_n[i_part];
    prev_concat_bentity1_entity2_interface      [i_part] = bentity1_entity2_interface      [i_part];
  }

  int         **concat_bentity1_entity2_n               = NULL;
  PDM_g_num_t **concat_bentity1_entity2_gnum            = NULL;
  int         **concat_bentity1_entity2_triplet         = NULL;
  int         **concat_bentity1_entity2_interface_n     = NULL;
  int         **concat_bentity1_entity2_interface_tot_n = NULL;
  int         **concat_bentity1_entity2_interface       = NULL;
  PDM_malloc(concat_bentity1_entity2_n              , n_part_tot, int         *);
  PDM_malloc(concat_bentity1_entity2_gnum           , n_part_tot, PDM_g_num_t *);
  PDM_malloc(concat_bentity1_entity2_triplet        , n_part_tot, int         *);
  PDM_malloc(concat_bentity1_entity2_interface_n    , n_part_tot, int         *);
  PDM_malloc(concat_bentity1_entity2_interface_tot_n, n_part_tot, int         *);
  PDM_malloc(concat_bentity1_entity2_interface      , n_part_tot, int         *);

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

      PDM_free(pnext_bentity1_entity2_n              [i_part]);
      PDM_free(pnext_bentity1_entity2_interface_tot_n[i_part]);
      pnext_bentity1_entity2_n              [i_part] = update_bentity1_entity2_n;
      pnext_bentity1_entity2_interface_tot_n[i_part] = update_bentity1_entity2_interface_tot_n;

      /*
       * Preapre for next step
       */
      PDM_free(prev_concat_bentity1_entity2_n              [i_part]);
      PDM_free(prev_concat_bentity1_entity2_gnum           [i_part]);
      PDM_free(prev_concat_bentity1_entity2_triplet        [i_part]);
      PDM_free(prev_concat_bentity1_entity2_interface_n    [i_part]);
      PDM_free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
      PDM_free(prev_concat_bentity1_entity2_interface      [i_part]);

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
        PDM_free(bentity1_entity2_n              [i_part]);
        PDM_free(bentity1_entity2_gnum           [i_part]);
        PDM_free(bentity1_entity2_triplet        [i_part]);
        PDM_free(bentity1_entity2_interface_n    [i_part]);
        PDM_free(bentity1_entity2_interface_tot_n[i_part]);
        PDM_free(bentity1_entity2_interface      [i_part]);
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

    PDM_free(pnext_bentity1_entity2_n              );
    PDM_free(pnext_bentity1_entity2_gnum           );
    PDM_free(pnext_bentity1_entity2_triplet        );
    PDM_free(pnext_bentity1_entity2_interface_n    );
    PDM_free(pnext_bentity1_entity2_interface_tot_n);
    PDM_free(pnext_bentity1_entity2_interface      );

  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    PDM_free(prev_concat_bentity1_entity2_n              [i_part]);
    PDM_free(prev_concat_bentity1_entity2_gnum           [i_part]);
    PDM_free(prev_concat_bentity1_entity2_triplet        [i_part]);
    PDM_free(prev_concat_bentity1_entity2_interface_n    [i_part]);
    PDM_free(prev_concat_bentity1_entity2_interface_tot_n[i_part]);
    PDM_free(prev_concat_bentity1_entity2_interface      [i_part]);
  }
  PDM_free(prev_concat_bentity1_entity2_n              );
  PDM_free(prev_concat_bentity1_entity2_gnum           );
  PDM_free(prev_concat_bentity1_entity2_triplet        );
  PDM_free(prev_concat_bentity1_entity2_interface_n    );
  PDM_free(prev_concat_bentity1_entity2_interface_tot_n);
  PDM_free(prev_concat_bentity1_entity2_interface      );


  // for(int i_part = 0; i_part < n_part_tot; ++i_part) {
  //   PDM_free(concat_bentity1_entity2_n              [i_part]);
  //   PDM_free(concat_bentity1_entity2_gnum           [i_part]);
  //   PDM_free(concat_bentity1_entity2_triplet        [i_part]);
  //   PDM_free(concat_bentity1_entity2_interface_n    [i_part]);
  //   PDM_free(concat_bentity1_entity2_interface_tot_n[i_part]);
  //   PDM_free(concat_bentity1_entity2_interface      [i_part]);
  // }
  PDM_free(concat_bentity1_entity2_n              );
  PDM_free(concat_bentity1_entity2_gnum           );
  PDM_free(concat_bentity1_entity2_triplet        );
  PDM_free(concat_bentity1_entity2_interface_n    );
  PDM_free(concat_bentity1_entity2_interface_tot_n);
  PDM_free(concat_bentity1_entity2_interface      );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    PDM_free(bentity1_entity2_gnum           [i_part]);
    PDM_free(bentity1_entity2_n              [i_part]);
    PDM_free(bentity1_entity2_interface_n    [i_part]);
    PDM_free(bentity1_entity2_interface      [i_part]);
    PDM_free(bentity1_entity2_triplet        [i_part]);
    PDM_free(bentity1_entity2_interface_tot_n[i_part]);
  }
  PDM_free(bentity1_entity2_gnum           );
  PDM_free(bentity1_entity2_n              );
  PDM_free(bentity1_entity2_interface_n    );
  PDM_free(bentity1_entity2_interface      );
  PDM_free(bentity1_entity2_triplet        );
  PDM_free(bentity1_entity2_interface_tot_n);

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
  int                         ***pentity1_bnd_proc_idx_in,
  int                         ***pentity1_bnd_part_idx_in,
  int                         ***pentity1_bnd_in,
  int                         ***pentity1_hint_in,
  int                            user_defined_bnd_graph,
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

  int          *pn_entity1            = NULL;
  PDM_g_num_t **pentity1_ln_to_gn     = NULL;
  int         **pentity1_hint         = NULL;
  PDM_malloc(pn_entity1       , n_part_tot, int          );
  PDM_malloc(pentity1_ln_to_gn, n_part_tot, PDM_g_num_t *);
  if(pentity1_hint_in != NULL) {
    PDM_malloc(pentity1_hint  , n_part_tot, int         *);
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
  int **part1_to_part2_idx       = NULL;
  int **part1_to_part2_triplet   = NULL;
  int **part1_to_part2_interface = NULL;
  PDM_malloc(part1_to_part2_idx      , n_part_tot, int *);
  PDM_malloc(part1_to_part2_triplet  , n_part_tot, int *);
  PDM_malloc(part1_to_part2_interface, n_part_tot, int *);

  int dom_itrf_defined = 1;
  if (pn_entity1_num==NULL) { // if pdi==NULL
    dom_itrf_defined = 0;
  }

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(comm, n_part_tot);

  int **ppart_entity1_proc_idx = *pentity1_bnd_proc_idx_in;
  int **ppart_entity1_part_idx = *pentity1_bnd_part_idx_in;
  int **ppart_entity1          = *pentity1_bnd_in;
  if (user_defined_bnd_graph==0) {
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

    PDM_free(part_distribution);
  }

  int* n_part_g = NULL;
  PDM_malloc(n_part_g, n_domain, int);
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);
  int n_g_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_g_part_tot += n_part_g[i_dom];
  }

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = NULL;
    PDM_malloc(n_part_shift, n_rank, int);
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
      PDM_malloc(part1_to_part2_idx[li_part], pn_entity1[li_part] + 1, int);
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pn_entity1[li_part]);

      int n_entity_bound = ppart_entity1_part_idx[li_part][n_g_part_tot];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_entity1[li_part][4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From domain interface */
      if (dom_itrf_defined==1) {
        for(int idx_entity = 0; idx_entity < pn_entity1_num[li_part]; ++idx_entity) {
          int i_entity = pentity1_num[li_part][idx_entity];
          int n_opp = pentity1_opp_location_idx[li_part][idx_entity+1] - pentity1_opp_location_idx[li_part][idx_entity];
          part1_to_part2_n[i_entity] += n_opp;
        }
      }


      for(int i_entity = 0; i_entity < pn_entity1[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pn_entity1[li_part]];
      PDM_malloc(part1_to_part2_triplet  [li_part], n_connect_tot  , int);
      PDM_malloc(part1_to_part2_interface[li_part], n_connect_tot/3, int);

      // printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_entity1[li_part][4*idx_entity]-1;
        int i_proc_opp   = ppart_entity1[li_part][4*idx_entity+1];
        int i_part_opp   = ppart_entity1[li_part][4*idx_entity+2]-1;
        int i_entity_opp = ppart_entity1[li_part][4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp + n_part_shift[i_proc_opp];
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From domain interface */
      if (dom_itrf_defined==1) {
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
      }

      PDM_free(part1_to_part2_n);

      li_part += 1;
    }

    PDM_free(n_part_shift);
  }

  PDM_free(n_part_g);

  *pentity1_extented_to_pentity1_idx_out       = part1_to_part2_idx;
  *pentity1_extented_to_pentity1_triplet_out   = part1_to_part2_triplet;
  *pentity1_extented_to_pentity1_interface_out = part1_to_part2_interface;

  PDM_free(pn_entity1          );
  PDM_free(pentity1_ln_to_gn   );

  if(pentity1_hint_in != NULL) {
    PDM_free(pentity1_hint);
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    if (dom_itrf_defined==1) {
      PDM_free(pentity1_num             [i_part]);
      PDM_free(pentity1_opp_location_idx[i_part]);
      PDM_free(pentity1_opp_location    [i_part]);
      PDM_free(pentity1_opp_interface   [i_part]);
      PDM_free(pentity1_opp_sens        [i_part]);
      PDM_free(pentity1_opp_gnum        [i_part]);
    }
    else if (dom_itrf_defined==0) {
      PDM_free(ppart_entity1            [i_part]);
      PDM_free(ppart_entity1_proc_idx   [i_part]);
      PDM_free(ppart_entity1_part_idx   [i_part]);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Invalid dom_itrf_defined value (%d) \n", dom_itrf_defined);
    }
  }
  if (dom_itrf_defined==1) {
    PDM_free(pn_entity1_num            );
    PDM_free(pentity1_num              );
    PDM_free(pentity1_opp_location_idx );
    PDM_free(pentity1_opp_location     );
    PDM_free(pentity1_opp_interface_idx);
    PDM_free(pentity1_opp_interface    );
    PDM_free(pentity1_opp_sens         );
    PDM_free(pentity1_opp_gnum         );
  }
  else if (dom_itrf_defined==0) {
    PDM_free(ppart_entity1             );
    PDM_free(ppart_entity1_proc_idx    );
    PDM_free(ppart_entity1_part_idx    );
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Invalid dom_itrf_defined value (%d) \n", dom_itrf_defined);
  }


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
  int                         **pentity2_entity1_idx,
  int                         **pentity2_entity1,
  int                           prev_dentity2_itrf_n_blk,
  PDM_g_num_t                  *prev_dentity2_itrf_blk_gnum,
  int                          *prev_dentity2_itrf_blk_ancstr_strd,
  PDM_g_num_t                  *prev_dentity2_itrf_blk_ancstr,
  int                          *prev_dentity2_itrf_blk_path_itrf_strd,
  int                          *prev_dentity2_itrf_blk_path_itrf,
  int                          *prev_dentity2_itrf_gnum_and_itrf_strid,
  PDM_g_num_t                  *prev_dentity2_itrf_gnum_and_itrf_data,
  int                         **pn_entity2_extended_out,
  PDM_g_num_t                ***pentity2_extended_ln_to_gn_out,
  int                        ***pentity2_extended_alrdy_sent_out,
  int                        ***pentity2_extended_to_pentity2_idx_out,
  int                        ***pentity2_extended_to_pentity2_triplet_out,
  int                        ***pentity2_extended_to_pentity2_interface_out,
  int                          *next_dentity2_itrf_n_blk_out,
  PDM_g_num_t                 **next_dentity2_itrf_blk_gnum_out,
  int                         **next_dentity2_itrf_blk_ancstr_strd_out,
  PDM_g_num_t                 **next_dentity2_itrf_blk_ancstr_out,
  int                         **next_dentity2_itrf_blk_path_itrf_strd_out,
  int                         **next_dentity2_itrf_blk_path_itrf_out,
  int                         **next_dentity2_itrf_gnum_and_itrf_strid_out,
  PDM_g_num_t                 **next_dentity2_itrf_gnum_and_itrf_data_out,
  PDM_MPI_Comm                  comm
)
{
  int debug = 0;

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

  if(debug == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      PDM_log_trace_array_int (pentity2_entity1[i_part], pentity2_entity1_idx[i_part][pn_entity2[i_part]], "pentity2_entity1 ::");

      log_trace("\n");
      PDM_log_trace_array_long(pentity1_ln_to_gn[i_part], pn_entity1[i_part], "pentity1_ln_to_gn ::");

      log_trace("\n");
      PDM_log_trace_array_int(pentity1_to_pentity1_idx      [i_part],                                  pn_entity1[i_part]+1 , "pentity1_to_pentity1_idx       ::");
      PDM_log_trace_array_int(pentity1_to_pentity1_triplet  [i_part], pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]  , "pentity1_to_pentity1_triplet   ::");
      PDM_log_trace_array_int(pentity1_to_pentity1_interface[i_part], pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3, "pentity1_to_pentity1_interface ::");

      log_trace("\n");
      PDM_log_trace_array_long(pentity2_ln_to_gn     [i_part],                                pn_entity2[i_part]  , "pentity2_ln_to_gn      ::");
      PDM_log_trace_array_int (pentity2_alrdy_sent   [i_part],                                pn_entity2[i_part]  , "pentity2_alrdy_sent    ::");
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

  if (debug == 1) {
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
  int         **gnum1_com_from_triplet_n      = NULL;
  int         **gnum1_com_from_triplet_send   = NULL;
  PDM_g_num_t **gnum1_com_from_gnum_send      = NULL;
  PDM_malloc(gnum1_com_from_triplet_n   , n_part, int         *);
  PDM_malloc(gnum1_com_from_triplet_send, n_part, int         *);
  PDM_malloc(gnum1_com_from_gnum_send   , n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][ n_ref_lnum2[i_part]];
    int *_pentity1_entity2_idx     = pentity1_entity2_idx    [i_part];

    /* Count while excluding faces:
     *   - that have already been sent at previous step
     */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        for(int idx_entity2 = _pentity1_entity2_idx[i_entity1]; idx_entity2 < _pentity1_entity2_idx[i_entity1+1]; ++idx_entity2) {
          int i_entity2 = PDM_ABS(pentity1_entity2[i_part][idx_entity2])-1;
          if (pentity2_alrdy_sent[i_part][i_entity2]==0) {
            n_send_part2 += 1;
          }
        }
      }
    }

    /* Allocate */
    PDM_malloc(gnum1_com_from_triplet_n   [i_part],     n_gnum1_come_from, int        );
    PDM_malloc(gnum1_com_from_triplet_send[i_part], 3 * n_send_part2     , int        );
    PDM_malloc(gnum1_com_from_gnum_send   [i_part],     n_send_part2     , PDM_g_num_t);

    int         *_gnum1_com_from_triplet_n    = gnum1_com_from_triplet_n     [i_part];
    int         *_gnum1_com_from_triplet_send = gnum1_com_from_triplet_send  [i_part];
    PDM_g_num_t *_gnum1_com_from_gnum_send    = gnum1_com_from_gnum_send     [i_part];

    int *pentity2_alrdy_sent_tmp      = NULL;
    PDM_malloc(pentity2_alrdy_sent_tmp, pn_entity2[i_part], int);
    for (int i_entity2=0; i_entity2<pn_entity2[i_part]; ++i_entity2) {
      pentity2_alrdy_sent_tmp[i_entity2] = pentity2_alrdy_sent[i_part][i_entity2];
    }


    int l_n_send_part2     = 0;
    n_send_part2     = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) { // Loop on referenced entities of part2
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) { // Loop on entities of part1 referenced by entity2
        
        // log_trace("    via interf_entity1 = %d, with entity = %d\n", interf_entity_1, gnum1_come_from[i_part][k]);
        l_n_send_part2     = 0;
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
            _gnum1_com_from_triplet_send [3*n_send_part2  ] = i_rank;
            _gnum1_com_from_triplet_send [3*n_send_part2+1] = i_part;
            _gnum1_com_from_triplet_send [3*n_send_part2+2] = i_entity2;

            pentity2_alrdy_sent_tmp[i_entity2] = 1 ;
            l_n_send_part2++;
            n_send_part2++;
            // }
          }
        }
        _gnum1_com_from_triplet_n[k] = l_n_send_part2;

      }
    }

    for (int i_entity2=0; i_entity2<pn_entity2[i_part]; ++i_entity2) {
      pentity2_alrdy_sent[i_part][i_entity2] = pentity2_alrdy_sent_tmp[i_entity2];
    }
    PDM_free(pentity2_alrdy_sent_tmp);

    if(debug == 1) {
      // Ce tableau contient pour chaque vertex d'un interface, ses faces associées.
      log_trace("\n");
      PDM_log_trace_array_int (_gnum1_com_from_triplet_n    , n_ref_lnum2[i_part], "_gnum1_com_from_triplet_n    ::");
      PDM_log_trace_array_long(_gnum1_com_from_gnum_send    , n_send_part2       , "_gnum1_com_from_gnum_send    ::");
      PDM_log_trace_array_int (_gnum1_com_from_triplet_send , n_send_part2*3     , "_gnum1_com_from_triplet_send ::");
    }


  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity1_entity2_idx[i_part]);
    PDM_free(pentity1_entity2    [i_part]);
  }
  PDM_free(pentity1_entity2_idx);
  PDM_free(pentity1_entity2    );


  int exch_request = -1;
  int         **pextract_entity2_n = NULL;

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

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(gnum1_com_from_triplet_n     [i_part]);
    PDM_free(gnum1_com_from_triplet_send  [i_part]);
    PDM_free(gnum1_com_from_gnum_send     [i_part]);
  }
  PDM_free(gnum1_com_from_triplet_n     );
  PDM_free(gnum1_com_from_triplet_send  );
  PDM_free(gnum1_com_from_gnum_send     );

  PDM_part_to_part_free(ptp);



  
  // > Compute n_entity2 for every entity1 and n_interface for each entity2 of every entity1
  int  *pn_entity1_entity2             = NULL;
  int **pextract_entity2_idx           = NULL;
  int  *pextract_entity1_entity2_ntot  = NULL;
  PDM_malloc(pn_entity1_entity2           , n_part, int  );
  PDM_malloc(pextract_entity2_idx         , n_part, int *);
  PDM_malloc(pextract_entity1_entity2_ntot, n_part, int  );;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    
    pn_entity1_entity2            [i_part] = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
    pextract_entity2_idx          [i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity2_n      [i_part],                              pn_entity1_entity2[i_part]);
    pextract_entity1_entity2_ntot [i_part] = pextract_entity2_idx[i_part][pn_entity1_entity2[i_part]];
    
    if (debug==1) {
      PDM_log_trace_array_int(pextract_entity2_idx          [i_part],                              n_part1_to_part2 +1, "pextract_entity2_idx           ::");
    }
  }


  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_part1_to_part2                 = pentity1_to_pentity1_idx[i_part][pn_entity1[i_part]]/3;
      int n_part1_to_part2_recv_tot        =                                        pextract_entity2_idx[i_part][n_part1_to_part2];

      log_trace("\n");
      PDM_log_trace_array_long(pextract_entity2_gnum     [i_part], n_part1_to_part2_recv_tot       , "_pextract_entity2_gnum     ::");
    }


    log_trace("\n");
    log_trace("------------------\n");
    log_trace("previous data_base\n");
    int n_itrf = 0;
    int n_data = 0;
    int n_ancstr = 0;
    for(int i = 0; i < prev_dentity2_itrf_n_blk; ++i) {
      n_itrf += prev_dentity2_itrf_blk_path_itrf_strd [i];
      n_data += prev_dentity2_itrf_gnum_and_itrf_strid[i];
      n_ancstr += prev_dentity2_itrf_blk_ancstr_strd[i];
    }
    log_trace("\n");
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_gnum           , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_gnum            ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_ancstr_strd    , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_ancstr_strd     ::");
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_ancstr         , n_ancstr                , "prev_dentity2_itrf_blk_ancstr          ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_path_itrf_strd , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_path_itrf_strd  ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_path_itrf      , n_itrf                  , "prev_dentity2_itrf_blk_path_itrf       ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_gnum_and_itrf_strid, prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(prev_dentity2_itrf_gnum_and_itrf_data , 2 * n_data              , "prev_dentity2_itrf_gnum_and_itrf_data  ::");

  }







  /*
   * We have to check that composition (gnum, interface(s?)) doesn't already exist
   * So we ask block interface data_base about entity candidate.
   *
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("---------------------------\n");
    log_trace("BTP to find valid entities2\n");

    _tell_me_more_data_base(prev_dentity2_itrf_n_blk, 
                            prev_dentity2_itrf_blk_gnum, 
                            prev_dentity2_itrf_blk_ancstr_strd, 
                            prev_dentity2_itrf_blk_ancstr, 
                            prev_dentity2_itrf_blk_path_itrf_strd, 
                            prev_dentity2_itrf_blk_path_itrf, 
                            prev_dentity2_itrf_gnum_and_itrf_strid, 
                            prev_dentity2_itrf_gnum_and_itrf_data, 
                            278);
  }

  int         **pextract_entity2_kind = NULL;
  int         **pextract_entity2_lnum = NULL;
  int         **pextract_entity2_itrf = NULL;
  PDM_g_num_t **pextract_entity2_ancstr         = NULL;
  int         **pextract_entity2_path_itrf_strd = NULL;
  int         **pextract_entity2_path_itrf      = NULL;
  int         **pextract_entity2_kind_idx  = NULL;
  int         **pextract_entity2_kind_ordr = NULL;
  int         **pextract_entity2_twin_itrf_idx  = NULL;
  PDM_g_num_t **pextract_entity2_twin_itrf_gnum = NULL;
  int         **pextract_entity2_twin_itrf_itrf = NULL;


  int n_kind = 6;
  _find_valid_entities(n_part,
                       pn_entity1_entity2,
                       pentity1_to_pentity1_interface,
                       pextract_entity2_idx,
                       pextract_entity2_gnum,
                       pextract_entity2_triplet,
                       pextract_entity1_entity2_ntot,
                       pn_entity2,
                       pentity2_ln_to_gn,
                       prev_dentity2_itrf_n_blk, // Distributed data_base
                       prev_dentity2_itrf_blk_gnum,
                       prev_dentity2_itrf_blk_ancstr_strd,
                       prev_dentity2_itrf_blk_ancstr,
                       prev_dentity2_itrf_blk_path_itrf_strd,
                       prev_dentity2_itrf_blk_path_itrf,
                       prev_dentity2_itrf_gnum_and_itrf_strid,
                       prev_dentity2_itrf_gnum_and_itrf_data,
                       NULL,
                      &pextract_entity2_kind,
                      &pextract_entity2_lnum,
                      &pextract_entity2_itrf,
                       NULL,
                      &pextract_entity2_ancstr,
                      &pextract_entity2_path_itrf_strd,
                      &pextract_entity2_path_itrf,
                      &pextract_entity2_kind_idx,
                      &pextract_entity2_kind_ordr,
                      &pextract_entity2_twin_itrf_idx,
                      &pextract_entity2_twin_itrf_gnum,
                      &pextract_entity2_twin_itrf_itrf,
                       NULL,
                       comm);


  if (debug==1) {
    // PDM_log_trace_array_int(pn_new_entity2_n_itrf, n_part  , "pn_new_entity2_n_itrf ::");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pextract_entity2_gnum[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_gnum ::");
      PDM_log_trace_array_int (pextract_entity2_kind[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_kind ::");
      PDM_log_trace_array_int (pextract_entity2_lnum[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_lnum ::");
      PDM_log_trace_array_int (pextract_entity2_itrf[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_itrf ::");
      
      PDM_log_trace_array_long(pextract_entity2_ancstr[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_ancstr ::");
      PDM_log_trace_array_int (pextract_entity2_path_itrf_strd[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_path_itrf_strd ::");
      
      PDM_log_trace_array_int (pextract_entity2_kind_idx [i_part], n_kind+1, "pextract_entity2_kind_idx  ::");
      PDM_log_trace_array_int (pextract_entity2_kind_ordr[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_kind_ordr ::");
    }
  }



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
   * Compute gnum of extracted entities2 of kind4 using nuplet(ancestor, sort(interf1, interf2, ...))
   */
  PDM_g_num_t **pnew_entity2_gnum           = NULL;
  PDM_g_num_t **pnew_entity2_ancstr         = NULL;
  int         **pnew_entity2_path_itrf_strd = NULL;
  int         **pnew_entity2_path_itrf      = NULL;
  PDM_g_num_t **pnew_entity2_parent_nuplet  = NULL;
  _compute_gnum_from_ancestor_and_itrfs(n_part,
                                        shift_by_domain_entity2,
                                        pn_entity1_entity2,
                                        pextract_entity2_idx,
                                        pextract_entity2_kind,
                                        pextract_entity2_gnum,
                                        pextract_entity2_itrf,
                                        pextract_entity1_entity2_ntot,
                                        pextract_entity2_ancstr,
                                        pextract_entity2_path_itrf_strd,
                                        pextract_entity2_path_itrf,
                                        pextract_entity2_kind_idx,
                                        pextract_entity2_kind_ordr,
                                       &pextract_entity2_gnum,
                                       &pnew_entity2_gnum,
                                       &pnew_entity2_ancstr,
                                       &pnew_entity2_path_itrf_strd,
                                       &pnew_entity2_path_itrf,
                                       &pnew_entity2_parent_nuplet,
                                        comm);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    if (debug==1) {
      int n_new_num = pextract_entity2_kind_idx[i_part][6]-pextract_entity2_kind_idx[i_part][5];
      PDM_log_trace_array_long(pextract_entity2_gnum[i_part], pextract_entity1_entity2_ntot[i_part], "pextract_entity2_gnum             ::");
      log_trace("\n");
      PDM_log_trace_array_long(pnew_entity2_gnum  [i_part], n_new_num, "pextract_entity2_gnum   ::");
      PDM_log_trace_array_long(pnew_entity2_ancstr[i_part], n_new_num, "pextract_entity2_ancstr ::");
    }
    PDM_free(pextract_entity2_ancstr[i_part]);
    PDM_free(pextract_entity2_path_itrf_strd[i_part]);
    PDM_free(pextract_entity2_path_itrf[i_part]);
  }
  PDM_free(pextract_entity2_ancstr);
  PDM_free(pextract_entity2_path_itrf_strd);
  PDM_free(pextract_entity2_path_itrf);
  PDM_free(pn_entity1_entity2);


  /*
   * Keep entities of kind 2,3,4
   */

  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------------------------\n");
    log_trace("unique entities while updating connectivity\n");
  }

  int keep_all_parent = 1;

  int          *pn_entity2_extended                = NULL;
  PDM_g_num_t **pextended_entity2_gnum             = NULL;
  // PDM_g_num_t **pextended_entity2_ancstr           = NULL;
  int         **pextended_entity2_alrdy_sent       = NULL;
  // int         **pextended_entity2_path_itrf_idx    = NULL;
  // int         **pextended_entity2_path_itrf        = NULL;
  int         **pextended_entity2_to_entity2_idx   = NULL;
  int         **pextended_entity2_to_entity2_trplt = NULL;
  int         **pextended_entity2_to_entity2_itrf  = NULL;
  int         **pextended_entity2_to_entity2_sens  = NULL;

  _unique_entities_and_update_connectivity( n_part,
                                            pn_entity2,
                                            NULL,
                                            pextract_entity2_lnum,
                                            pextract_entity2_itrf,
                                            NULL,
                                            pextract_entity2_triplet,
                                            pextract_entity2_gnum,
                                            pextract_entity2_kind,
                                            pextract_entity2_kind_idx,
                                            pextract_entity2_kind_ordr,
                                           &pn_entity2_extended,
                                            NULL,
                                           &pextended_entity2_gnum,
                                           &pextended_entity2_alrdy_sent,
                                           &pextended_entity2_to_entity2_idx,
                                           &pextended_entity2_to_entity2_trplt,
                                           &pextended_entity2_to_entity2_itrf,
                                           &pextended_entity2_to_entity2_sens,
                                            keep_all_parent);
  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_triplet_data = pextended_entity2_to_entity2_idx[i_part][pn_entity2_extended[i_part]];
      PDM_log_trace_array_long(pextended_entity2_gnum            [i_part], pn_entity2_extended[i_part]  , "pextended_entity2_gnum             ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_idx  [i_part], pn_entity2_extended[i_part]+1, "pextended_entity2_to_entity2_idx   ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_trplt[i_part], n_triplet_data               , "pextended_entity2_to_entity2_trplt ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_itrf [i_part], n_triplet_data/3             , "pextended_entity2_to_entity2_itrf  ::");

      _tell_me_more_link(pn_entity2_extended               [i_part],
                         pextended_entity2_gnum            [i_part],
                         pextended_entity2_to_entity2_idx  [i_part],
                         pextended_entity2_to_entity2_trplt[i_part],
                         pextended_entity2_to_entity2_itrf [i_part],
                         NULL,
           (PDM_g_num_t) 282);
    }
  }


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pextract_entity2_triplet   [i_part]);
  }
  PDM_free(pextract_entity2_triplet   );

  /**
   * FOR fix ansctr issue:
   * put ancstr in DB
   * get ansctr and path in find_valid_entities
   * for kind4 build new path and new ancstr at this stage (use it for new gnum)
   * enrich block interf with this info
   * should be good
   */




  /*
   * Update new block of interface
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------------------\n");
    log_trace("complete data_base with new entities2\n");
  }

  // > Prepare array for _enrich_block_interface
  int  *p_db_n_entity2   = NULL;
  PDM_malloc(p_db_n_entity2, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    p_db_n_entity2  [i_part] = pextract_entity2_kind_idx[i_part][6]-pextract_entity2_kind_idx[i_part][5];
  }
  
  int          dn_entity2_itrf    = 0;
  PDM_g_num_t *dentity2_itrf_gnum = NULL;
  int         *dentity2_itrf_ancstr_strd    = NULL;
  PDM_g_num_t *dentity2_itrf_ancstr         = NULL;
  int         *dentity2_itrf_path_itrf_strd = NULL;
  int         *dentity2_itrf_path_itrf      = NULL;
  int         *dentity2_itrf_strd = NULL;
  PDM_g_num_t *dentity2_itrf_data = NULL;
  _enrich_block_interface(n_part,
                          p_db_n_entity2,
                          pnew_entity2_gnum,
                          pnew_entity2_parent_nuplet,
                          NULL,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_gnum_and_itrf_strid,
                          prev_dentity2_itrf_gnum_and_itrf_data,
                          NULL,
                         &dn_entity2_itrf,
                         &dentity2_itrf_gnum,
                         &dentity2_itrf_strd,
                         &dentity2_itrf_data,
                          NULL,
                          comm);

  if(debug == 1) {
    int n_data = 0;
    for (int i=0; i<dn_entity2_itrf;++i) {
      n_data += dentity2_itrf_strd[i];
    }
    log_trace("\n");
    log_trace("dn_entity2_itrf = %d\n", dn_entity2_itrf);
    PDM_log_trace_array_long(dentity2_itrf_gnum, dn_entity2_itrf, "dentity_itrf_gnum (Unique) ::");
    PDM_log_trace_array_int (dentity2_itrf_strd, dn_entity2_itrf, "dentity_itrf_strd (Unique) ::");
    PDM_log_trace_array_long(dentity2_itrf_data, 2*n_data       , "dentity_itrf_data (Unique) ::");
  }


  /*
   * Transfer interface from prev to next
   */
  int          *p_db_n_entity2_with_twin    = PDM_array_copy_int(p_db_n_entity2, n_part);
  int         **pentity2_twin_itrf_lnum     = NULL;
  PDM_g_num_t **pnew_entity2_gnum_with_twin = NULL;
  PDM_malloc(pentity2_twin_itrf_lnum    , n_part, int         *);
  PDM_malloc(pnew_entity2_gnum_with_twin, n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int i_write = 0;
    PDM_malloc(pentity2_twin_itrf_lnum[i_part], p_db_n_entity2[i_part], int);;
    for(int i = pextract_entity2_kind_idx[i_part][5]; i < pextract_entity2_kind_idx[i_part][6]; ++i) {
      int i_pos = pextract_entity2_kind_ordr[i_part][i];
      pentity2_twin_itrf_lnum[i_part][i_write] = i_pos;
      i_write++;
    }
    pnew_entity2_gnum_with_twin[i_part] = PDM_array_copy_gnum(pnew_entity2_gnum[i_part], p_db_n_entity2[i_part]);
  }


  _find_twin_interface(n_part,
                       dn_entity2_itrf, // Distributed data_base
                       dentity2_itrf_gnum,
                       dentity2_itrf_strd,
                       dentity2_itrf_data,
                       NULL,
                       pextract_entity1_entity2_ntot,
                       pextract_entity2_gnum,
                       pextract_entity2_itrf,
                       pextract_entity2_twin_itrf_idx,
                       pextract_entity2_twin_itrf_itrf,
                       NULL,
                       pextract_entity2_twin_itrf_gnum,
                       pentity2_twin_itrf_lnum,
                      &p_db_n_entity2_with_twin,
                      &pnew_entity2_gnum_with_twin,
                      &pnew_entity2_parent_nuplet,
                       NULL,
                       comm);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity2_twin_itrf_lnum[i_part]);
  }
  PDM_free(pentity2_twin_itrf_lnum);

  int          dn_entity2_itrf_with_twin    = 0;
  PDM_g_num_t *dentity2_itrf_gnum_with_twin = NULL;
  int         *dentity2_itrf_strd_with_twin = NULL;
  PDM_g_num_t *dentity2_itrf_data_with_twin = NULL;

  _enrich_block_interface(n_part,
                          p_db_n_entity2_with_twin,
                          pnew_entity2_gnum_with_twin,
                          pnew_entity2_parent_nuplet,
                          NULL,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_gnum_and_itrf_strid,
                          prev_dentity2_itrf_gnum_and_itrf_data,
                          NULL,
                         &dn_entity2_itrf_with_twin,
                         &dentity2_itrf_gnum_with_twin,
                         &dentity2_itrf_strd_with_twin,
                         &dentity2_itrf_data_with_twin,
                          NULL,
                          comm);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnew_entity2_gnum_with_twin[i_part]);
  }
  PDM_free(p_db_n_entity2_with_twin);
  PDM_free(pnew_entity2_gnum_with_twin);

  _enrich_block_ancestors(n_part,
                          p_db_n_entity2,
                          pnew_entity2_gnum,
                          pnew_entity2_ancstr,
                          pnew_entity2_path_itrf_strd,
                          pnew_entity2_path_itrf,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_blk_ancstr_strd,
                          prev_dentity2_itrf_blk_ancstr,
                          prev_dentity2_itrf_blk_path_itrf_strd,
                          prev_dentity2_itrf_blk_path_itrf,
                          dn_entity2_itrf_with_twin,
                          dentity2_itrf_gnum_with_twin,
                         &dentity2_itrf_ancstr_strd,
                         &dentity2_itrf_ancstr,
                         &dentity2_itrf_path_itrf_strd,
                         &dentity2_itrf_path_itrf,
                          comm);

    
  if(debug == 1) {
    int n_data = 0;
    int n_itrf = 0;
    int n_ancstr = 0;
    for (int i=0; i<dn_entity2_itrf_with_twin;++i) {
      n_data += dentity2_itrf_strd_with_twin[i];
      n_itrf += dentity2_itrf_path_itrf_strd[i];
      n_ancstr += dentity2_itrf_ancstr_strd[i];
    }
    log_trace("\n");
    log_trace("dn_entity2_itrf_with_twin = %d\n", dn_entity2_itrf_with_twin);
    PDM_log_trace_array_long(dentity2_itrf_gnum_with_twin          ,     dn_entity2_itrf_with_twin, "dentity_itrf_gnum_with_twin (Unique) ::");
    PDM_log_trace_array_int (dentity2_itrf_ancstr_strd             ,     dn_entity2_itrf_with_twin, "dentity2_itrf_ancstr_strd    (with twin) ::");
    PDM_log_trace_array_long(dentity2_itrf_ancstr                  ,     n_ancstr                 , "dentity2_itrf_ancstr         (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_path_itrf_strd          ,     dn_entity2_itrf_with_twin, "dentity2_itrf_path_itrf_strd (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_path_itrf               ,     n_itrf                   , "dentity2_itrf_path_itrf      (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_strd_with_twin          ,     dn_entity2_itrf_with_twin, "dentity_itrf_strd_with_twin (Unique) ::");
    PDM_log_trace_array_long(dentity2_itrf_data_with_twin          ,     2*n_data                 , "dentity_itrf_data_with_twin (Unique) ::");
  }

  /*
   * Keep pointer for output
   */
  *next_dentity2_itrf_n_blk_out               = dn_entity2_itrf_with_twin;
  *next_dentity2_itrf_blk_gnum_out            = dentity2_itrf_gnum_with_twin;
  *next_dentity2_itrf_blk_ancstr_strd_out     = dentity2_itrf_ancstr_strd;
  *next_dentity2_itrf_blk_ancstr_out          = dentity2_itrf_ancstr;
  *next_dentity2_itrf_blk_path_itrf_strd_out  = dentity2_itrf_path_itrf_strd;
  *next_dentity2_itrf_blk_path_itrf_out       = dentity2_itrf_path_itrf;
  *next_dentity2_itrf_gnum_and_itrf_strid_out = dentity2_itrf_strd_with_twin;
  *next_dentity2_itrf_gnum_and_itrf_data_out  = dentity2_itrf_data_with_twin;

  // > Free tmp data_base
  PDM_free(dentity2_itrf_gnum);
  PDM_free(dentity2_itrf_strd);
  PDM_free(dentity2_itrf_data);
  






  for(int i_part = 0; i_part < n_part; ++i_part) {

    PDM_free(pnew_entity2_gnum[i_part]);
    PDM_free(pnew_entity2_ancstr[i_part]);
    PDM_free(pnew_entity2_path_itrf_strd[i_part]);
    PDM_free(pnew_entity2_path_itrf[i_part]);
    PDM_free(pnew_entity2_parent_nuplet[i_part]);

    PDM_free(pextract_entity2_idx       [i_part]);
    PDM_free(pextract_entity2_n         [i_part]);
    PDM_free(pextract_entity2_gnum      [i_part]);

    PDM_free(pextract_entity2_kind      [i_part]);
    PDM_free(pextract_entity2_lnum      [i_part]);
    PDM_free(pextract_entity2_itrf      [i_part]);

    PDM_free(pextract_entity2_twin_itrf_idx[i_part]);
    PDM_free(pextract_entity2_twin_itrf_gnum[i_part]);
    PDM_free(pextract_entity2_twin_itrf_itrf[i_part]);
    
    PDM_free(pextract_entity2_kind_idx [i_part]);
    PDM_free(pextract_entity2_kind_ordr[i_part]);
  }
  PDM_free(p_db_n_entity2);
  PDM_free(pnew_entity2_gnum);
  PDM_free(pnew_entity2_ancstr);
  PDM_free(pnew_entity2_path_itrf_strd);
  PDM_free(pnew_entity2_path_itrf);
  PDM_free(pnew_entity2_parent_nuplet);

  PDM_free(pextract_entity2_idx);
  PDM_free(pextract_entity2_n         );
  PDM_free(pextract_entity2_gnum      );
  PDM_free(pextract_entity2_kind      );
  PDM_free(pextract_entity2_lnum      );
  PDM_free(pextract_entity2_itrf      );

  PDM_free(pextract_entity2_twin_itrf_idx);
  PDM_free(pextract_entity2_twin_itrf_gnum);
  PDM_free(pextract_entity2_twin_itrf_itrf);

  PDM_free(pextract_entity2_kind_idx);
  PDM_free(pextract_entity2_kind_ordr);

  PDM_free(pextract_entity1_entity2_ntot);


  // > Return
  *pn_entity2_extended_out                     = pn_entity2_extended;
  *pentity2_extended_ln_to_gn_out              = pextended_entity2_gnum;
  *pentity2_extended_alrdy_sent_out            = pextended_entity2_alrdy_sent;
  *pentity2_extended_to_pentity2_idx_out       = pextended_entity2_to_entity2_idx;
  *pentity2_extended_to_pentity2_triplet_out   = pextended_entity2_to_entity2_trplt;
  *pentity2_extended_to_pentity2_interface_out = pextended_entity2_to_entity2_itrf;
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

  int          *pn_entity1            = NULL;
  PDM_g_num_t **pentity1_ln_to_gn     = NULL;
  int         **pentity1_hint         = NULL;
  int          *pn_entity2            = NULL;
  PDM_g_num_t **pentity2_ln_to_gn     = NULL;
  int         **pentity2_entity1_idx  = NULL;
  int         **pentity2_entity1      = NULL;
  PDM_malloc(pn_entity1          , n_part_tot, int          );
  PDM_malloc(pentity1_ln_to_gn   , n_part_tot, PDM_g_num_t *);
  PDM_malloc(pn_entity2          , n_part_tot, int          );
  PDM_malloc(pentity2_ln_to_gn   , n_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity2_entity1_idx, n_part_tot, int         *);
  PDM_malloc(pentity2_entity1    , n_part_tot, int         *);
  if(pentity1_hint_in != NULL) {
    PDM_malloc(pentity1_hint     , n_part_tot, int         *);
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
  int **part1_to_part2_idx             = NULL;
  int **part1_to_part2_triplet_idx     = NULL;
  int **part1_to_part2_triplet         = NULL;
  int **part1_to_part2_interface       = NULL;
  int **part1_to_part2_entity2_n       = NULL;
  int **part1_to_part2_entity2_triplet = NULL;
  PDM_malloc(part1_to_part2_idx            , n_part_tot, int *);
  PDM_malloc(part1_to_part2_triplet        , n_part_tot, int *);
  PDM_malloc(part1_to_part2_interface      , n_part_tot, int *);
  PDM_malloc(part1_to_part2_entity2_n      , n_part_tot, int *);
  PDM_malloc(part1_to_part2_entity2_triplet, n_part_tot, int *);


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

  PDM_free(part_distribution);

  int* n_part_g = NULL;
  PDM_malloc(n_part_g, n_domain, int);
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);
  int n_g_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_g_part_tot += n_part_g[i_dom];
  }

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    int *n_part_shift = NULL;
    PDM_malloc(n_part_shift, n_rank, int);
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
      PDM_malloc(part1_to_part2_idx[li_part], pn_entity1[li_part] + 1, int);
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
      PDM_malloc(part1_to_part2_triplet  [li_part], n_connect_tot  , int);
      PDM_malloc(part1_to_part2_interface[li_part], n_connect_tot/3, int);

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
      PDM_malloc(part1_to_part2_entity2_n[li_part], n_send, int);
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


      PDM_malloc(part1_to_part2_entity2_triplet[li_part], 3 * n_send_entity2, int);
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

      PDM_free(part1_to_part2_n);

      li_part += 1;
    }

    PDM_free(n_part_shift);
  }

  PDM_free(n_part_g);

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
  int         **gnum1_com_from_triplet_n    = NULL;
  int         **gnum1_com_from_triplet_send = NULL;
  PDM_g_num_t **gnum1_com_from_gnum_send    = NULL;
  PDM_malloc(gnum1_com_from_triplet_n   , n_part_tot, int         *);
  PDM_malloc(gnum1_com_from_triplet_send, n_part_tot, int         *);
  PDM_malloc(gnum1_com_from_gnum_send   , n_part_tot, PDM_g_num_t *);
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
    PDM_malloc(gnum1_com_from_triplet_n   [i_part],     n_gnum1_come_from, int        );
    PDM_malloc(gnum1_com_from_triplet_send[i_part], 3 * n_send_part2     , int        );
    PDM_malloc(gnum1_com_from_gnum_send   [i_part],     n_send_part2     , PDM_g_num_t);

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
    PDM_free(gnum1_com_from_triplet_n   [i_part]);
    PDM_free(gnum1_com_from_triplet_send[i_part]);
    PDM_free(gnum1_com_from_gnum_send   [i_part]);
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


  int **pextract_entity2_idx = NULL;
  PDM_malloc(pextract_entity2_idx, ln_part_tot, int *);
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

  PDM_free(gnum1_com_from_triplet_n   );
  PDM_free(gnum1_com_from_triplet_send);
  PDM_free(gnum1_com_from_gnum_send   );


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

  int          *pn_entity2_only_by_interface        = NULL;
  int         **pentity2_interface                  = NULL;
  PDM_g_num_t **pentity2_ln_to_gn_only_by_interface = NULL;
  PDM_malloc(pn_entity2_only_by_interface       , ln_part_tot, int          );
  PDM_malloc(pentity2_interface                 , ln_part_tot, int         *);
  PDM_malloc(pentity2_ln_to_gn_only_by_interface, ln_part_tot, PDM_g_num_t *);

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

    PDM_malloc(pentity2_ln_to_gn_only_by_interface[i_part], 2 * pn_entity2_only_by_interface[i_part], PDM_g_num_t);
    PDM_malloc(pentity2_interface                 [i_part],     n_part1_to_part2_recv_tot           , int        );
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
  int          *pn_entity2_extented                     = NULL;
  int         **pentity2_extented_to_pentity2_idx       = NULL;
  int         **pentity2_extented_to_pentity2_triplet   = NULL;
  PDM_g_num_t **pentity2_extented_ln_to_gn              = NULL;
  PDM_g_num_t **extented_entity2_orig_gnum              = NULL;
  int         **pentity2_extented_to_pentity2_interface = NULL;
  PDM_malloc(pn_entity2_extented                    , ln_part_tot, int          );
  PDM_malloc(pentity2_extented_to_pentity2_idx      , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_to_pentity2_triplet  , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_ln_to_gn             , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(extented_entity2_orig_gnum             , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity2_extented_to_pentity2_interface, ln_part_tot, int         *);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t *_pextract_entity2_gnum = pextract_entity2_gnum[i_part];
    int         *_pextract_entity2_idx  = pextract_entity2_idx [i_part];

    int n_part1_to_part2 = part1_to_part2_idx[i_part][pn_entity1[i_part]]/3;
    int n_part1_to_part2_recv_tot = _pextract_entity2_idx[n_part1_to_part2];

    // All gnum has be unified / shift w/r of interface, only sort along gnum
    int *order = NULL;
    PDM_malloc(order, n_part1_to_part2_recv_tot, int);
    for(int i_entity2 = 0; i_entity2 < n_part1_to_part2_recv_tot; ++i_entity2) {
      order[i_entity2] = i_entity2;
    }

    // Maybe faire inplace sort and revert with order after ?
    PDM_g_num_t* sorted_pentity2_ln_to_gn = NULL;
    PDM_malloc(sorted_pentity2_ln_to_gn, pn_entity2[i_part], PDM_g_num_t);
    for(int i = 0; i < pn_entity2[i_part]; ++i) {
      sorted_pentity2_ln_to_gn[i] = pentity2_ln_to_gn[i_part][i];
    }
    PDM_sort_long(sorted_pentity2_ln_to_gn, NULL, pn_entity2[i_part]);

    int n_unique = PDM_inplace_unique_long_and_order(_pextract_entity2_gnum,
                                                     order,
                                                     0,
                                                     n_part1_to_part2_recv_tot-1);

    PDM_malloc(pentity2_extented_ln_to_gn[i_part], n_unique, PDM_g_num_t);
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
    PDM_malloc(pentity2_extented_to_pentity2_idx      [i_part],     n_unique2 + 1, int        );
    PDM_malloc(pentity2_extented_to_pentity2_triplet  [i_part], 3 * n_unique2    , int        );
    PDM_malloc(extented_entity2_orig_gnum             [i_part],     n_unique2    , PDM_g_num_t);
    PDM_malloc(pentity2_extented_to_pentity2_interface[i_part],     n_unique2    , int        );

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

    PDM_free(order);
    PDM_free(sorted_pentity2_ln_to_gn);

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
    PDM_free(extented_entity2_orig_gnum[i_part]);
  }
  PDM_free(extented_entity2_orig_gnum);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pentity2_ln_to_gn_only_by_interface[i_part]);
  }

  PDM_free(pn_entity2_only_by_interface);
  PDM_free(pentity2_ln_to_gn_only_by_interface);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_entity2_idx[i_part]);
    PDM_free(pentity2_interface  [i_part]);
  }
  PDM_free(pextract_entity2_idx);
  PDM_free(pentity2_interface);



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_entity2_n        [i_part]);
    PDM_free(pextract_entity2_gnum     [i_part]);
    PDM_free(pextract_entity2_triplet  [i_part]);
    PDM_free(pentity1_entity2_idx      [i_part]);
    PDM_free(pentity1_entity2          [i_part]);
  }
  PDM_free(pextract_entity2_n        );
  PDM_free(pextract_entity2_gnum     );
  PDM_free(pextract_entity2_triplet  );
  PDM_free(pentity1_entity2_idx      );
  PDM_free(pentity1_entity2          );


  PDM_part_to_part_free(ptp);


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    PDM_free(ppart_entity1_proc_idx[i_part]);
    PDM_free(ppart_entity1_part_idx[i_part]);
    PDM_free(ppart_entity1         [i_part]);
  }
  PDM_free(ppart_entity1_proc_idx);
  PDM_free(ppart_entity1_part_idx);
  PDM_free(ppart_entity1         );


  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    PDM_free(part1_to_part2_idx            [i_part] );
    PDM_free(part1_to_part2_triplet        [i_part]);
    PDM_free(part1_to_part2_interface      [i_part]);
    PDM_free(part1_to_part2_entity2_n      [i_part]);
    PDM_free(part1_to_part2_entity2_triplet[i_part]);
  }

  PDM_free(part1_to_part2_idx            );
  PDM_free(part1_to_part2_triplet        );
  PDM_free(part1_to_part2_interface      );
  PDM_free(part1_to_part2_entity2_n      );
  PDM_free(part1_to_part2_entity2_triplet);

  PDM_free(pn_entity1          );
  PDM_free(pentity1_ln_to_gn   );
  PDM_free(pn_entity2          );
  PDM_free(pentity2_ln_to_gn   );
  PDM_free(pentity2_entity1_idx);
  PDM_free(pentity2_entity1    );

  if(pentity1_hint_in != NULL) {
    PDM_free(pentity1_hint);
  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pentity1_num             [i_part]);
    PDM_free(pentity1_opp_location_idx[i_part]);
    PDM_free(pentity1_opp_location    [i_part]);
    PDM_free(pentity1_opp_interface   [i_part]);
    PDM_free(pentity1_opp_sens        [i_part]);
    PDM_free(pentity1_opp_gnum        [i_part]);
  }
  PDM_free(pn_entity1_num            );
  PDM_free(pentity1_num              );
  PDM_free(pentity1_opp_location_idx );
  PDM_free(pentity1_opp_location     );
  PDM_free(pentity1_opp_interface_idx);
  PDM_free(pentity1_opp_interface    );
  PDM_free(pentity1_opp_sens         );
  PDM_free(pentity1_opp_gnum         );

}


void
PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2
(
 int                    n_part,
 int                    n_interface,
 PDM_g_num_t            shift_by_domain_entity2,
 int                    prev_dentity2_itrf_n_blk,
 PDM_g_num_t           *prev_dentity2_itrf_blk_gnum,
 int                   *prev_dentity2_itrf_blk_ancstr_strd,
 PDM_g_num_t           *prev_dentity2_itrf_blk_ancstr,
 int                   *prev_dentity2_itrf_blk_path_itrf_strd,
 int                   *prev_dentity2_itrf_blk_path_itrf,
 int                   *prev_dentity2_itrf_gnum_and_itrf_strid,
 PDM_g_num_t           *prev_dentity2_itrf_gnum_and_itrf_data,
 int                   *prev_dentity2_itrf_gnum_and_itrf_sens,
 int                   *pn_entity1,
 PDM_g_num_t          **pentity1_ln_to_gn,
 int                   *pn_entity2,
 PDM_g_num_t          **pentity2_ln_to_gn,
 int                  **pentity1_entity2_idx,
 int                  **pentity1_entity2,
 int                   *pn_entity1_extended,
 PDM_g_num_t          **pentity1_extended_ln_to_gn,
 int                  **pentity1_extended_to_pentity1_idx,
 int                  **pentity1_extended_to_pentity1_triplet,
 int                  **pentity1_extended_to_pentity1_interface,
 int                  **pentity1_extended_to_pentity1_sens,
 int                    keep_all_parent,
 int                  **pn_entity2_extended_out,
 PDM_g_num_t         ***pentity2_extended_ln_to_gn_out,
 int                 ***pextended_entity1_entity2_idx_out,
 int                 ***pextended_entity1_entity2_out,
 int                 ***pentity2_extended_to_pentity2_idx_out,
 int                 ***pentity2_extended_to_pentity2_triplet_out,
 int                 ***pentity2_extended_to_pentity2_interface_out,
 int                 ***pentity2_extended_to_pentity2_sens_out,
 int                   *next_dentity2_itrf_n_blk_out,
 PDM_g_num_t          **next_dentity2_itrf_blk_gnum_out,
 int                  **next_dentity2_itrf_blk_ancstr_strd_out,
 PDM_g_num_t          **next_dentity2_itrf_blk_ancstr_out,
 int                  **next_dentity2_itrf_blk_path_itrf_strd_out,
 int                  **next_dentity2_itrf_blk_path_itrf_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_strid_out,
 PDM_g_num_t          **next_dentity2_itrf_gnum_and_itrf_data_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_sens_out,
 PDM_MPI_Comm           comm
)
{
  PDM_UNUSED(n_interface);
  PDM_UNUSED(pentity1_ln_to_gn);

  int debug = 0;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int entity1_has_sens = 0;
  if (pentity1_extended_to_pentity1_sens!=NULL) {
    entity1_has_sens = 1;
  }
  int entity2_has_sens = 0;
  if (prev_dentity2_itrf_gnum_and_itrf_sens!=NULL) {
    entity2_has_sens = 1;
  }


  /*
   * Creation du ptp between extended part and current partition for entity1
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------\n");
    log_trace("PTP to exchange entities2\n");
    log_trace("\n");
    log_trace("entity1_has_sens = %d\n", entity1_has_sens);
    log_trace("entity2_has_sens = %d\n", entity2_has_sens);
  }

  if(debug == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      // > Faces
      _tell_me_more_connectivity(pn_entity1[i_part],
                                 pentity1_ln_to_gn[i_part],
                                 pentity1_entity2_idx[i_part],
                                 pentity1_entity2[i_part],
                                 pentity2_ln_to_gn[i_part],
                   (PDM_g_num_t) 45);
      
      // > Edges
      _tell_me_more_connectivity(pn_entity1[i_part],
                                 pentity1_ln_to_gn[i_part],
                                 pentity1_entity2_idx[i_part],
                                 pentity1_entity2[i_part],
                                 pentity2_ln_to_gn[i_part],
                   (PDM_g_num_t) 97);

      log_trace("\n");
      PDM_log_trace_array_long(pentity1_extended_ln_to_gn[i_part], pn_entity1_extended[i_part], "pentity1_extended_ln_to_gn ::");
      
      log_trace("\n");
      PDM_log_trace_array_int(pentity1_extended_to_pentity1_idx      [i_part],                                           pn_entity1_extended[i_part]+1 , "pentity1_extended_to_pentity1_idx       ::");
      PDM_log_trace_array_int(pentity1_extended_to_pentity1_triplet  [i_part], pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]  , "pentity1_extended_to_pentity1_triplet   ::");
      PDM_log_trace_array_int(pentity1_extended_to_pentity1_interface[i_part], pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3, "pentity1_extended_to_pentity1_interface ::");
      if (entity1_has_sens==1) {
        PDM_log_trace_array_int(pentity1_extended_to_pentity1_sens     [i_part], pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3, "pentity1_extended_to_pentity1_sens      ::");
      }

      log_trace("\n");
      PDM_log_trace_array_long(pentity2_ln_to_gn[i_part], pn_entity2[i_part]  , "pentity2_ln_to_gn      ::");
      
      int pentity1_entity2_size = pentity1_entity2_idx[i_part][pn_entity1[i_part]];
      log_trace("\n");
      log_trace("pn_entity1 = %d\n", pn_entity1[i_part]);
      PDM_log_trace_array_int(pentity1_entity2_idx[i_part], pn_entity1[i_part]+1 , "pentity1_entity2_idx ::");
      PDM_log_trace_array_int(pentity1_entity2    [i_part], pentity1_entity2_size, "pentity1_entity2     ::");
    }
  }

  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pentity1_extended_ln_to_gn,
                                                                      (const int          *) pn_entity1_extended,
                                                                                             n_part,
                                                                      (const int          *) pn_entity1,
                                                                                             n_part,
                                                                      (const int         **) pentity1_extended_to_pentity1_idx,
                                                                      (const int         **) NULL,
                                                                      (const int         **) pentity1_extended_to_pentity1_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  int         **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &gnum1_come_from_idx, &gnum1_come_from);

  int **gnum1_come_from_sens = NULL;
  if (entity1_has_sens==1) {

    int exch_request = -1;
    PDM_part_to_part_iexch(ptp,
                           PDM_MPI_COMM_KIND_P2P,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                           1,
                           sizeof(int),
          (const int **)   NULL,
          (const void **)  pentity1_extended_to_pentity1_sens,
                           NULL,
              (void ***)   &gnum1_come_from_sens,
                           &exch_request);
    PDM_part_to_part_iexch_wait(ptp, exch_request);

  }


  /*
   * Remplissage des buffers et envoi
   */
  int         **gnum1_com_from_entity1_entity2_n       = NULL;
  PDM_g_num_t **gnum1_com_from_entity1_entity2         = NULL;
  int         **gnum1_com_from_entity1_entity2_sign    = NULL;
  int         **gnum1_com_from_entity1_entity2_triplet = NULL;
  PDM_malloc(gnum1_com_from_entity1_entity2_n      , n_part, int         *);
  PDM_malloc(gnum1_com_from_entity1_entity2        , n_part, PDM_g_num_t *);
  PDM_malloc(gnum1_com_from_entity1_entity2_sign   , n_part, int         *);
  PDM_malloc(gnum1_com_from_entity1_entity2_triplet, n_part, int         *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_gnum1_come_from = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]];

    /* Count */
    int n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        for(int idx_entity2 = pentity1_entity2_idx[i_part][i_entity1]; idx_entity2 < pentity1_entity2_idx[i_part][i_entity1+1]; ++idx_entity2) {
          n_send_part2 += 1;
        }
      }
    }

    /* Allocate */
    PDM_malloc(gnum1_com_from_entity1_entity2_n      [i_part],     n_gnum1_come_from, int        );
    PDM_malloc(gnum1_com_from_entity1_entity2        [i_part],     n_send_part2     , PDM_g_num_t);
    PDM_malloc(gnum1_com_from_entity1_entity2_sign   [i_part],     n_send_part2     , int        );
    PDM_malloc(gnum1_com_from_entity1_entity2_triplet[i_part], 3 * n_send_part2     , int        );


    n_send_part2 = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      int i_entity1 = ref_lnum2[i_part][i_ref]-1;
      for(int k = gnum1_come_from_idx[i_part][i_ref]; k < gnum1_come_from_idx[i_part][i_ref+1]; ++k) {
        int sens = 1;
        if (entity1_has_sens==1) {
          sens = gnum1_come_from_sens[i_part][k];
        }
        gnum1_com_from_entity1_entity2_n[i_part][k] = pentity1_entity2_idx[i_part][i_entity1+1] - pentity1_entity2_idx[i_part][i_entity1];
        if (entity2_has_sens==0) { // > typically face_vtx: to inverse orientation we need to change vtx order and not vtx sign
          if (sens==1) {
            for(int idx_entity2 = pentity1_entity2_idx[i_part][i_entity1]; idx_entity2 < pentity1_entity2_idx[i_part][i_entity1+1]; ++idx_entity2) {
              int i_entity2      = PDM_ABS (pentity1_entity2[i_part][idx_entity2])-1;
              int i_entity2_sign = PDM_SIGN(pentity1_entity2[i_part][idx_entity2]);
              gnum1_com_from_entity1_entity2        [i_part][  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
              gnum1_com_from_entity1_entity2_sign   [i_part][  n_send_part2  ] = i_entity2_sign;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2  ] = i_rank;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+1] = i_part;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+2] = i_entity2;
              n_send_part2++;
            }
          }
          else {
            for(int idx_entity2 = pentity1_entity2_idx[i_part][i_entity1+1]-1; idx_entity2 >= pentity1_entity2_idx[i_part][i_entity1]; --idx_entity2) {
              int i_entity2      = PDM_ABS (pentity1_entity2[i_part][idx_entity2])-1;
              int i_entity2_sign = PDM_SIGN(pentity1_entity2[i_part][idx_entity2]);
              gnum1_com_from_entity1_entity2        [i_part][  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
              gnum1_com_from_entity1_entity2_sign   [i_part][  n_send_part2  ] = i_entity2_sign;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2  ] = i_rank;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+1] = i_part;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+2] = i_entity2;
              n_send_part2++;
            }
          }
        }
        else { // > typically face_edge: to inverse orientation we need to change edge sign and not edge order
          if (sens==1) {
            for(int idx_entity2 = pentity1_entity2_idx[i_part][i_entity1]; idx_entity2 < pentity1_entity2_idx[i_part][i_entity1+1]; ++idx_entity2) {
              int i_entity2      = PDM_ABS (pentity1_entity2[i_part][idx_entity2])-1;
              int i_entity2_sign = PDM_SIGN(pentity1_entity2[i_part][idx_entity2]);
              gnum1_com_from_entity1_entity2        [i_part][  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
              gnum1_com_from_entity1_entity2_sign   [i_part][  n_send_part2  ] = i_entity2_sign;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2  ] = i_rank;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+1] = i_part;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+2] = i_entity2;
              n_send_part2++;
            }
          }
          else {
            for(int idx_entity2 = pentity1_entity2_idx[i_part][i_entity1]; idx_entity2 < pentity1_entity2_idx[i_part][i_entity1+1]; ++idx_entity2) {
              int i_entity2      = PDM_ABS (pentity1_entity2[i_part][idx_entity2])-1;
              int i_entity2_sign = PDM_SIGN(pentity1_entity2[i_part][idx_entity2]);
              gnum1_com_from_entity1_entity2        [i_part][  n_send_part2  ] = pentity2_ln_to_gn[i_part][i_entity2];
              gnum1_com_from_entity1_entity2_sign   [i_part][  n_send_part2  ] = -i_entity2_sign;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2  ] = i_rank;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+1] = i_part;
              gnum1_com_from_entity1_entity2_triplet[i_part][3*n_send_part2+2] = i_entity2;
              n_send_part2++;
            }
          }
        }


      }
    }

    if(debug == 1) {
      log_trace("\n");
      PDM_log_trace_array_long(gnum1_com_from_entity1_entity2        [i_part], n_send_part2  , "gnum1_com_from_gnum_send               ::");
      PDM_log_trace_array_int (gnum1_com_from_entity1_entity2_sign   [i_part], n_send_part2  , "gnum1_com_from_gnum_send_sign          ::");
      PDM_log_trace_array_int (gnum1_com_from_entity1_entity2_triplet[i_part], n_send_part2*3, "gnum1_com_from_entity1_entity2_triplet ::");
    }

  }

  /*
   * Exchange:
   *     - vertices gnum    of face from initial domain referenced
   *     - vertices triplet of face from initial domain referenced
   */
  int exch_request = -1;
  int         **pextract_entity1_entity2_n = NULL;

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

  int **pextract_entity1_entity2_sign = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(int),
                (const int **)   gnum1_com_from_entity1_entity2_n,
                (const void **)  gnum1_com_from_entity1_entity2_sign,
                                 &pextract_entity1_entity2_n,
                    (void ***)   &pextract_entity1_entity2_sign,
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
    if (entity1_has_sens==1) {
      PDM_free(gnum1_come_from_sens                [i_part]);
    }
    PDM_free(gnum1_com_from_entity1_entity2_n      [i_part]);
    PDM_free(gnum1_com_from_entity1_entity2        [i_part]);
    PDM_free(gnum1_com_from_entity1_entity2_sign   [i_part]);
    PDM_free(gnum1_com_from_entity1_entity2_triplet[i_part]);
  }
  if (entity1_has_sens==1) {
    PDM_free(gnum1_come_from_sens);
  }
  PDM_free(gnum1_com_from_entity1_entity2_n);
  PDM_free(gnum1_com_from_entity1_entity2);
  PDM_free(gnum1_com_from_entity1_entity2_sign);
  PDM_free(gnum1_com_from_entity1_entity2_triplet);
  PDM_part_to_part_free(ptp);



  /*
   * Prepare unification for all incoming entity2 by interface
   *
   */


  if(debug == 1) {
    log_trace("\n");
    log_trace("------------------\n");
    log_trace("previous data_base\n");
    log_trace("prev_dentity2_itrf_n_blk = %d\n", prev_dentity2_itrf_n_blk);
    int n_itrf   = 0;
    int n_data   = 0;
    int n_ancstr = 0;
    for (int i=0; i<prev_dentity2_itrf_n_blk; ++i) {
      n_itrf   += prev_dentity2_itrf_blk_path_itrf_strd [i];
      n_data   += prev_dentity2_itrf_gnum_and_itrf_strid[i];
      n_ancstr += prev_dentity2_itrf_blk_ancstr_strd    [i];
    }
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_gnum           , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_gnum            ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_path_itrf_strd , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_path_itrf_strd  ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_ancstr_strd    , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_ancstr_strd     ::");
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_ancstr         , n_ancstr                , "prev_dentity2_itrf_blk_ancstr          ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_path_itrf_strd , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_blk_path_itrf_strd  ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_blk_path_itrf      , n_itrf                  , "prev_dentity2_itrf_blk_path_itrf       ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_gnum_and_itrf_strid, prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_gnum_and_itrf_strid ::");
    PDM_log_trace_array_long(prev_dentity2_itrf_gnum_and_itrf_data , 2*n_data                , "prev_dentity2_itrf_gnum_and_itrf_data  ::");
    if (entity2_has_sens==1) {
      PDM_log_trace_array_int (prev_dentity2_itrf_gnum_and_itrf_sens ,   n_data                , "prev_dentity2_itrf_gnum_and_itrf_sens  ::");
    }
  }

  /**
   * 
   * TODO: bypass if empty
   */


  /*
   * We have to check that composition (gnum, interface(s?)) doesn't already exist
   * So we ask block interface data_base about entity candidate.
   *
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("---------------------------\n");
    log_trace("BTP to find valid entities2\n");
  }


  int  *pn_entity1_entity2            = NULL;
  int **pextract_entity1_entity2_idx  = NULL;
  int  *pextract_entity1_entity2_ntot = NULL;
  PDM_malloc(pn_entity1_entity2           , n_part, int  );
  PDM_malloc(pextract_entity1_entity2_idx , n_part, int *);
  PDM_malloc(pextract_entity1_entity2_ntot, n_part, int  );;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_entity1_entity2[i_part] = pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3;

    pextract_entity1_entity2_idx [i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity1_entity2_n[i_part], pn_entity1_entity2[i_part]);
    pextract_entity1_entity2_ntot[i_part] = pextract_entity1_entity2_idx[i_part][pn_entity1_entity2[i_part]];
    if (debug==1) {
      PDM_log_trace_array_int(pextract_entity1_entity2_idx[i_part], pn_entity1_entity2[i_part] +1, "pextract_entity1_entity2_idx ::");
    }
  }



  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      log_trace("i_part = %d\n", i_part);
      // > Faces
      _tell_me_more_received(pn_entity1_extended          [i_part],
                             pentity1_extended_ln_to_gn   [i_part],
                             pentity1_extended_to_pentity1_idx   [i_part],
                             pextract_entity1_entity2_idx [i_part],
                             pextract_entity1_entity2_gnum[i_part],
                             pextract_entity1_entity2_triplet[i_part],
                             pextract_entity1_entity2_sign[i_part],
                             45);
      // > Edges
      _tell_me_more_received(pn_entity1_extended          [i_part],
                             pentity1_extended_ln_to_gn   [i_part],
                             pentity1_extended_to_pentity1_idx   [i_part],
                             pextract_entity1_entity2_idx [i_part],
                             pextract_entity1_entity2_gnum[i_part],
                             pextract_entity1_entity2_triplet[i_part],
                             pextract_entity1_entity2_sign[i_part],
                             99);
    }
  }


  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_part1_to_part2          = pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3;
      int n_part1_to_part2_recv_tot = pextract_entity1_entity2_idx[i_part][n_part1_to_part2];

      log_trace("\n");
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum[i_part], n_part1_to_part2_recv_tot, "_pextract_entity1_entity2_gnum ::");
      PDM_log_trace_array_int (pextract_entity1_entity2_sign[i_part], n_part1_to_part2_recv_tot, "_pextract_entity1_entity2_sign ::");
    }
  }

  int         **pextract_entity2_kind = NULL;
  int         **pextract_entity2_lnum = NULL;
  int         **pextract_entity2_itrf = NULL;
  int         **pextract_entity2_sens = NULL;
  PDM_g_num_t **pextract_entity2_ancstr         = NULL;
  int         **pextract_entity2_path_itrf_strd = NULL;
  int         **pextract_entity2_path_itrf      = NULL;
  int         **pextract_entity2_kind_idx  = NULL;
  int         **pextract_entity2_kind_ordr = NULL;
  int         **pextract_entity2_twin_itrf_idx  = NULL;
  PDM_g_num_t **pextract_entity2_twin_itrf_gnum = NULL;
  int         **pextract_entity2_twin_itrf_itrf = NULL;
  int         **pextract_entity2_twin_itrf_sens = NULL;


  int n_kind = 6;
  _find_valid_entities(n_part,
                       pn_entity1_entity2,
                       pentity1_extended_to_pentity1_interface,
                       pextract_entity1_entity2_idx,
                       pextract_entity1_entity2_gnum,
                       pextract_entity1_entity2_triplet,
                       pextract_entity1_entity2_ntot,
                       pn_entity2,
                       pentity2_ln_to_gn,
                       prev_dentity2_itrf_n_blk, // Distributed data_base
                       prev_dentity2_itrf_blk_gnum,
                       prev_dentity2_itrf_blk_ancstr_strd,
                       prev_dentity2_itrf_blk_ancstr,
                       prev_dentity2_itrf_blk_path_itrf_strd,
                       prev_dentity2_itrf_blk_path_itrf,
                       prev_dentity2_itrf_gnum_and_itrf_strid,
                       prev_dentity2_itrf_gnum_and_itrf_data,
                       prev_dentity2_itrf_gnum_and_itrf_sens,
                      &pextract_entity2_kind,
                      &pextract_entity2_lnum,
                      &pextract_entity2_itrf,
                      &pextract_entity2_sens,
                      &pextract_entity2_ancstr,
                      &pextract_entity2_path_itrf_strd,
                      &pextract_entity2_path_itrf,
                      &pextract_entity2_kind_idx,
                      &pextract_entity2_kind_ordr,
                      &pextract_entity2_twin_itrf_idx,
                      &pextract_entity2_twin_itrf_gnum,
                      &pextract_entity2_twin_itrf_itrf,
                      &pextract_entity2_twin_itrf_sens,
                       comm);
  if (debug==1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("\n");
      log_trace("i_part = %d\n", i_part);

      // > Faces
      if (entity2_has_sens==1) {
        _tell_me_more_valid_entities(pn_entity1_extended          [i_part],
                                     pentity1_extended_ln_to_gn   [i_part],
                                     pentity1_extended_to_pentity1_idx   [i_part],
                                     pextract_entity1_entity2_idx [i_part],
                                     pextract_entity1_entity2_gnum[i_part],
                                     pextract_entity1_entity2_sign[i_part],
                                     pextract_entity2_kind[i_part],
                                     pextract_entity2_lnum[i_part],
                                     pextract_entity2_itrf[i_part],
                                     pextract_entity2_sens[i_part],
                                     45);
      } 
      else {
        _tell_me_more_valid_entities(pn_entity1_extended          [i_part],
                                     pentity1_extended_ln_to_gn   [i_part],
                                     pentity1_extended_to_pentity1_idx   [i_part],
                                     pextract_entity1_entity2_idx [i_part],
                                     pextract_entity1_entity2_gnum[i_part],
                                     NULL,
                                     pextract_entity2_kind[i_part],
                                     pextract_entity2_lnum[i_part],
                                     pextract_entity2_itrf[i_part],
                                     NULL,
                                     45);
      }


      int n_part1_to_part2          = pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3;
      int n_part1_to_part2_recv_tot = pextract_entity1_entity2_idx[i_part][n_part1_to_part2];

      log_trace("\n");
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum[i_part], n_part1_to_part2_recv_tot, "_pextract_entity1_entity2_gnum ::");
      PDM_log_trace_array_int (pextract_entity1_entity2_sign[i_part], n_part1_to_part2_recv_tot, "_pextract_entity1_entity2_sign ::");
    }
  }


  /*
   * TODO: use this function to compute gnum
   *       need to build path_itrf_idx, path_itrf, ancestor arrays
   *
   */

  PDM_g_num_t **p_db_entity_gnum           = NULL;
  PDM_g_num_t **p_db_entity_ancstr         = NULL;
  int         **p_db_entity_path_itrf_strd = NULL;
  int         **p_db_entity_path_itrf      = NULL;
  PDM_g_num_t **p_db_entity_data           = NULL;
  

  for(int i_part = 0; i_part < n_part; ++i_part) {
    if(debug == 1) {
      int n_part1_to_part2 = pextract_entity1_entity2_idx[i_part][pn_entity1_entity2[i_part]];
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum[i_part], n_part1_to_part2, "pextract_entity2_gnum ::" );
      PDM_log_trace_array_int (pextract_entity2_itrf        [i_part], n_part1_to_part2, "pextract_entity2_itrf ::" );
    }
  }
  

   _compute_gnum_from_ancestor_and_itrfs(n_part,
                                         shift_by_domain_entity2,
                                         pn_entity1_entity2,
                                         pextract_entity1_entity2_idx,
                                         pextract_entity2_kind,
                                         pextract_entity1_entity2_gnum,
                                         pextract_entity2_itrf,
                                         pextract_entity1_entity2_ntot,
                                         pextract_entity2_ancstr,
                                         pextract_entity2_path_itrf_strd,
                                         pextract_entity2_path_itrf,
                                         pextract_entity2_kind_idx,
                                         pextract_entity2_kind_ordr,
                                        &pextract_entity1_entity2_gnum,
                                        &p_db_entity_gnum,
                                        &p_db_entity_ancstr,
                                        &p_db_entity_path_itrf_strd,
                                        &p_db_entity_path_itrf,
                                        &p_db_entity_data,
                                         comm);

  if(debug == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_part1_to_part2 = pextract_entity1_entity2_idx[i_part][pn_entity1_entity2[i_part]];
      PDM_log_trace_array_int (pextract_entity1_entity2_idx [i_part], pn_entity1_entity2[i_part], "pextract_entity1_entity2_idx ::" );
      PDM_log_trace_array_int (pextract_entity2_kind        [i_part], n_part1_to_part2, "pextract_entity2_kind ::" );
      PDM_log_trace_array_int (pextract_entity2_lnum        [i_part], n_part1_to_part2, "pextract_entity2_lnum ::" );
      PDM_log_trace_array_int (pextract_entity2_itrf        [i_part], n_part1_to_part2, "pextract_entity2_itrf ::" );
      if (entity2_has_sens==1) {
        PDM_log_trace_array_int (pextract_entity2_sens        [i_part], n_part1_to_part2, "pextract_entity2_sens ::" );
      }
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum[i_part], n_part1_to_part2, "pextract_entity2_gnum ::" );
      PDM_log_trace_array_long(pextract_entity2_ancstr      [i_part], n_part1_to_part2, "pextract_entity2_ancstr ::" );
    }
  }

  int  *p_db_n_entity2   = NULL;
  int **p_db_entity_sens = NULL;
  PDM_malloc(p_db_n_entity2  , n_part, int  );
  PDM_malloc(p_db_entity_sens, n_part, int *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    p_db_n_entity2  [i_part] = pextract_entity2_kind_idx[i_part][n_kind] - pextract_entity2_kind_idx[i_part][n_kind-1];
    p_db_entity_sens[i_part] = PDM_array_const_int(p_db_n_entity2[i_part], 1);

    PDM_free(pextract_entity2_ancstr[i_part]);
    PDM_free(pextract_entity2_path_itrf_strd[i_part]);
    PDM_free(pextract_entity2_path_itrf[i_part]);
  }
  PDM_free(pextract_entity2_ancstr);
  PDM_free(pextract_entity2_path_itrf_strd);
  PDM_free(pextract_entity2_path_itrf);
  PDM_free(pn_entity1_entity2);







  /**
   * At this stage all gnum are generated and properly shifted
   * You need to update all incoming conectivity
   * For each kind (if new), we unique all entries in order to setup a new local order
   * 
   * --> Keep entities of kind 2,3,4 while updating connectivity
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("--------------------------------\n");
    log_trace("update connectivity of entities1\n");
  }

  // int keep_all_parent = 1;

  int          *pn_entity2_extended                = NULL;
  PDM_g_num_t **pextended_entity2_gnum             = NULL;
  int         **pextended_entity2_alrdy_sent       = NULL;
  int         **pextended_entity2_to_entity2_idx   = NULL;
  int         **pextended_entity2_to_entity2_trplt = NULL;
  int         **pextended_entity2_to_entity2_itrf  = NULL;
  int         **pextended_entity2_to_entity2_sens  = NULL;

  // > Allocate connectivity array
  int **pextended_entity1_entity2 = NULL;
  PDM_malloc(pextended_entity1_entity2, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_part1_to_part2 = pentity1_extended_to_pentity1_idx[i_part][pn_entity1_extended[i_part]]/3;
    pextended_entity1_entity2[i_part] = PDM_array_const_int(pextract_entity1_entity2_idx[i_part][n_part1_to_part2], 0);

    if(debug == 1) {
      int n_entity = pextract_entity2_kind_idx[i_part][6];
      PDM_log_trace_array_long(pextract_entity1_entity2_gnum[i_part], n_entity, "pextract_entity1_entity2_gnum   ::" );
    }
  }

  _unique_entities_and_update_connectivity( n_part,
                                            pn_entity2,
                                            pextract_entity1_entity2_sign,
                                            pextract_entity2_lnum,
                                            pextract_entity2_itrf,
                                            pextract_entity2_sens,
                                            pextract_entity1_entity2_triplet,
                                            pextract_entity1_entity2_gnum,
                                            pextract_entity2_kind,
                                            pextract_entity2_kind_idx,
                                            pextract_entity2_kind_ordr,
                                           &pn_entity2_extended,
                                           &pextended_entity1_entity2,
                                           &pextended_entity2_gnum,
                                           &pextended_entity2_alrdy_sent,
                                           &pextended_entity2_to_entity2_idx,
                                           &pextended_entity2_to_entity2_trplt,
                                           &pextended_entity2_to_entity2_itrf,
                                           &pextended_entity2_to_entity2_sens,
                                            keep_all_parent);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pextended_entity2_alrdy_sent[i_part]);

    if(debug == 1) {
      log_trace("\n");
      log_trace("i_part = %d\n", i_part);
      int n_part1_to_part2 = pextract_entity1_entity2_idx[i_part][pn_entity1_extended[i_part]];
      int n_triplet_data = pextended_entity2_to_entity2_idx[i_part][pn_entity2_extended[i_part]];
      PDM_log_trace_array_long(pextended_entity2_gnum            [i_part], pn_entity2_extended[i_part]  , "pextended_entity2_gnum             ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_idx  [i_part], pn_entity2_extended[i_part]+1, "pextended_entity2_to_entity2_idx   ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_trplt[i_part], n_triplet_data               , "pextended_entity2_to_entity2_trplt ::");
      PDM_log_trace_array_int (pextended_entity2_to_entity2_itrf [i_part], n_triplet_data/3             , "pextended_entity2_to_entity2_itrf  ::");
      if (entity2_has_sens==1) {
        PDM_log_trace_array_int (pextended_entity2_to_entity2_sens [i_part], n_triplet_data/3           , "pextended_entity2_to_entity2_sens  ::");
      }
      PDM_log_trace_array_int (pextract_entity1_entity2_idx      [i_part], pn_entity1_extended[i_part]  , "pextended_entity1_entity2_idx      ::");
      PDM_log_trace_array_int (pextended_entity1_entity2         [i_part], n_part1_to_part2             , "pextended_entity1_entity2          ::");

      if (entity2_has_sens==1) {
        _tell_me_more_link(pn_entity2_extended               [i_part],
                           pextended_entity2_gnum            [i_part],
                           pextended_entity2_to_entity2_idx  [i_part],
                           pextended_entity2_to_entity2_trplt[i_part],
                           pextended_entity2_to_entity2_itrf [i_part],
                           pextended_entity2_to_entity2_sens [i_part],
             (PDM_g_num_t) 99);
      }
      else {
        _tell_me_more_link(pn_entity2_extended               [i_part],
                           pextended_entity2_gnum            [i_part],
                           pextended_entity2_to_entity2_idx  [i_part],
                           pextended_entity2_to_entity2_trplt[i_part],
                           pextended_entity2_to_entity2_itrf [i_part],
                           NULL,
             (PDM_g_num_t) 99);
      }
    }
  }
  PDM_free(pextended_entity2_alrdy_sent);
  



  /*
   * Update new block of interface
   */
  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------------\n");
    log_trace("update entity1 interfaces\n");
  }

  if(debug == 1) {
    int n_data = 0;
    for (int i=0; i<prev_dentity2_itrf_n_blk;++i) {
      n_data += prev_dentity2_itrf_gnum_and_itrf_strid[i];
    }
    log_trace("\n");
    PDM_log_trace_array_long(prev_dentity2_itrf_blk_gnum           , prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_gnum (Unique) ::");
    PDM_log_trace_array_int (prev_dentity2_itrf_gnum_and_itrf_strid, prev_dentity2_itrf_n_blk, "prev_dentity2_itrf_strf (Unique) ::");
    PDM_log_trace_array_long(prev_dentity2_itrf_gnum_and_itrf_data , 2*n_data                , "prev_dentity2_itrf_data (Unique) ::");
  }


  
  int          dn_entity2_itrf    = 0;
  PDM_g_num_t *dentity2_itrf_gnum = NULL;
  int         *dentity2_itrf_ancstr_strd    = NULL;
  PDM_g_num_t *dentity2_itrf_ancstr         = NULL;
  int         *dentity2_itrf_path_itrf_strd = NULL;
  int         *dentity2_itrf_path_itrf      = NULL;
  int         *dentity2_itrf_strd = NULL;
  PDM_g_num_t *dentity2_itrf_data = NULL;
  int         *dentity2_itrf_sens = NULL;

  if(debug == 1) {
    log_trace("\n");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("p_db_n_entity2[i_part] = %d\n", p_db_n_entity2[i_part]);
      PDM_log_trace_array_long(p_db_entity_gnum  [i_part],     p_db_n_entity2[i_part], "p_db_entity_gnum   (Unique) ::");
      PDM_log_trace_array_long(p_db_entity_ancstr[i_part],     p_db_n_entity2[i_part], "p_db_entity_ancstr (Unique) ::");
      PDM_log_trace_array_long(p_db_entity_data  [i_part], 2 * p_db_n_entity2[i_part], "p_db_entity_data   (Unique) ::");
    }
  }
  _enrich_block_interface(n_part,
                          p_db_n_entity2,
                          p_db_entity_gnum,
                          p_db_entity_data,
                          p_db_entity_sens,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_gnum_and_itrf_strid,
                          prev_dentity2_itrf_gnum_and_itrf_data,
                          prev_dentity2_itrf_gnum_and_itrf_sens,
                         &dn_entity2_itrf,
                         &dentity2_itrf_gnum,
                         &dentity2_itrf_strd,
                         &dentity2_itrf_data,
                         &dentity2_itrf_sens,
                          comm);

  if(debug == 1) {
    int n_data = 0;
    for (int i=0; i<dn_entity2_itrf;++i) {
      n_data += dentity2_itrf_strd[i];
    }
    log_trace("\n");
    log_trace("dn_entity2_itrf = %d\n", dn_entity2_itrf);
    PDM_log_trace_array_long(dentity2_itrf_gnum,     dn_entity2_itrf, "dentity2_itrf_gnum (Unique) ::");
    PDM_log_trace_array_int (dentity2_itrf_strd,     dn_entity2_itrf, "dentity2_itrf_strd (Unique) ::");
    PDM_log_trace_array_long(dentity2_itrf_data, 2 * n_data         , "dentity2_itrf_data (Unique) ::");
    if (entity2_has_sens==1) {
      PDM_log_trace_array_int (dentity2_itrf_sens,     n_data         , "dentity2_itrf_sens (Unique) ::");
    }
  }




  /*
   * Transfer interface from prev to next
   * TODO: gérer la stride variable ??
   * TODO: bypass if no twin.
   */

  if (debug==1) {
    log_trace("\n\n");
    log_trace("--------------------\n");
    log_trace("find twin interfaces\n");
  }
  int          *p_db_n_entity2_with_twin   = PDM_array_copy_int(p_db_n_entity2, n_part);
  int         **pentity2_twin_itrf_lnum    = NULL;
  PDM_g_num_t **p_db_entity_gnum_with_twin = NULL;
  PDM_malloc(pentity2_twin_itrf_lnum   , n_part, int         *);
  PDM_malloc(p_db_entity_gnum_with_twin, n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int i_write = 0;
    PDM_malloc(pentity2_twin_itrf_lnum[i_part], p_db_n_entity2[i_part], int);;
    for(int i = pextract_entity2_kind_idx[i_part][5]; i < pextract_entity2_kind_idx[i_part][6]; ++i) {
      int i_pos = pextract_entity2_kind_ordr[i_part][i];
      pentity2_twin_itrf_lnum[i_part][i_write] = i_pos;
      i_write++;
    }
    p_db_entity_gnum_with_twin[i_part] = PDM_array_copy_gnum(p_db_entity_gnum[i_part], p_db_n_entity2[i_part]);
  }


  _find_twin_interface(n_part,
                       dn_entity2_itrf, // Distributed data_base
                       dentity2_itrf_gnum,
                       dentity2_itrf_strd,
                       dentity2_itrf_data,
                       dentity2_itrf_sens,
                       pextract_entity1_entity2_ntot,
                       pextract_entity1_entity2_gnum,
                       pextract_entity2_itrf,
                       pextract_entity2_twin_itrf_idx,
                       pextract_entity2_twin_itrf_itrf,
                       pextract_entity2_twin_itrf_sens,
                       pextract_entity2_twin_itrf_gnum,
                       pentity2_twin_itrf_lnum,
                      &p_db_n_entity2_with_twin,
                      &p_db_entity_gnum_with_twin,
                      &p_db_entity_data,
                      &p_db_entity_sens,
                       comm);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pentity2_twin_itrf_lnum[i_part]);
  }
  PDM_free(pentity2_twin_itrf_lnum);




  if (debug==1) {
    log_trace("\n\n");
    log_trace("-------------------\n");
    log_trace("enrich db with twin\n");
  }
  int          dn_entity2_itrf_with_twin    = 0;
  PDM_g_num_t *dentity2_itrf_gnum_with_twin = NULL;
  int         *dentity2_itrf_strd_with_twin = NULL;
  PDM_g_num_t *dentity2_itrf_data_with_twin = NULL;
  int         *dentity2_itrf_sens_with_twin = NULL;

  if(debug == 1) {
    log_trace("\n");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      log_trace("p_db_n_entity2_with_twin[i_part] = %d\n", p_db_n_entity2_with_twin[i_part]);
      PDM_log_trace_array_long(p_db_entity_gnum_with_twin[i_part],     p_db_n_entity2_with_twin[i_part], "p_db_entity_gnum (with_twin) ::");
      PDM_log_trace_array_long(p_db_entity_data          [i_part], 2 * p_db_n_entity2_with_twin[i_part], "p_db_entity_data (with_twin) ::");
      if (entity2_has_sens==1) {
        PDM_log_trace_array_int (p_db_entity_sens          [i_part],     p_db_n_entity2_with_twin[i_part], "p_db_entity_sens (with_twin) ::");
      }
    }
  }
  _enrich_block_interface(n_part,
                          p_db_n_entity2_with_twin,
                          p_db_entity_gnum_with_twin,
                          p_db_entity_data,
                          p_db_entity_sens,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_gnum_and_itrf_strid,
                          prev_dentity2_itrf_gnum_and_itrf_data,
                          prev_dentity2_itrf_gnum_and_itrf_sens,
                         &dn_entity2_itrf_with_twin,
                         &dentity2_itrf_gnum_with_twin,
                         &dentity2_itrf_strd_with_twin,
                         &dentity2_itrf_data_with_twin,
                         &dentity2_itrf_sens_with_twin,
                          comm);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(p_db_entity_gnum_with_twin[i_part]);
  }
  PDM_free(p_db_n_entity2_with_twin);
  PDM_free(p_db_entity_gnum_with_twin);
  
  _enrich_block_ancestors(n_part,
                          p_db_n_entity2,
                          p_db_entity_gnum,
                          p_db_entity_ancstr,
                          p_db_entity_path_itrf_strd,
                          p_db_entity_path_itrf,
                          prev_dentity2_itrf_n_blk,
                          prev_dentity2_itrf_blk_gnum,
                          prev_dentity2_itrf_blk_ancstr_strd,
                          prev_dentity2_itrf_blk_ancstr,
                          prev_dentity2_itrf_blk_path_itrf_strd,
                          prev_dentity2_itrf_blk_path_itrf,
                          dn_entity2_itrf_with_twin,
                          dentity2_itrf_gnum_with_twin,
                         &dentity2_itrf_ancstr_strd,
                         &dentity2_itrf_ancstr,
                         &dentity2_itrf_path_itrf_strd,
                         &dentity2_itrf_path_itrf,
                          comm);
  if(debug == 1) {
    int n_data   = 0;
    int n_itrf   = 0;
    int n_ancstr = 0;
    for (int i=0; i<dn_entity2_itrf_with_twin;++i) {
      n_data   += dentity2_itrf_strd_with_twin[i];
      n_itrf   += dentity2_itrf_path_itrf_strd[i];
      n_ancstr += dentity2_itrf_ancstr_strd[i];
    }
    log_trace("\n");
    log_trace("dn_entity2_itrf_with_twin = %d\n", dn_entity2_itrf_with_twin);
    PDM_log_trace_array_long(dentity2_itrf_gnum_with_twin,     dn_entity2_itrf_with_twin, "dentity2_itrf_gnum           (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_ancstr_strd   ,     dn_entity2_itrf_with_twin, "dentity2_itrf_ancstr_strd    (with twin) ::");
    PDM_log_trace_array_long(dentity2_itrf_ancstr        ,     n_ancstr                 , "dentity2_itrf_ancstr         (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_path_itrf_strd,     dn_entity2_itrf_with_twin, "dentity2_itrf_path_itrf_strd (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_path_itrf     ,     n_itrf                   , "dentity2_itrf_path_itrf      (with twin) ::");
    PDM_log_trace_array_int (dentity2_itrf_strd_with_twin,     dn_entity2_itrf_with_twin, "dentity2_itrf_strd           (with twin) ::");
    PDM_log_trace_array_long(dentity2_itrf_data_with_twin, 2 * n_data                   , "dentity2_itrf_data           (with twin) ::");
    if (entity2_has_sens==1) {
      PDM_log_trace_array_int (dentity2_itrf_sens_with_twin,     n_data                   , "dentity2_itrf_sens           (with twin) ::");
    }
  }

  /*
   * Keep pointer for output
   */
  *next_dentity2_itrf_n_blk_out               = dn_entity2_itrf_with_twin;
  *next_dentity2_itrf_blk_gnum_out            = dentity2_itrf_gnum_with_twin;
  *next_dentity2_itrf_blk_ancstr_strd_out     = dentity2_itrf_ancstr_strd;
  *next_dentity2_itrf_blk_ancstr_out          = dentity2_itrf_ancstr;
  *next_dentity2_itrf_blk_path_itrf_strd_out  = dentity2_itrf_path_itrf_strd;
  *next_dentity2_itrf_blk_path_itrf_out       = dentity2_itrf_path_itrf;
  *next_dentity2_itrf_gnum_and_itrf_strid_out = dentity2_itrf_strd_with_twin;
  *next_dentity2_itrf_gnum_and_itrf_data_out  = dentity2_itrf_data_with_twin;
  if (entity2_has_sens==1) {
    *next_dentity2_itrf_gnum_and_itrf_sens_out  = dentity2_itrf_sens_with_twin;
  }
  // > Free previous data_base
  PDM_free(dentity2_itrf_gnum);
  PDM_free(dentity2_itrf_strd);
  PDM_free(dentity2_itrf_data);
  PDM_free(dentity2_itrf_sens);

  // > Free tmp arrays used to enrich data_base

  /* Free all working array */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(p_db_entity_gnum[i_part]);
    PDM_free(p_db_entity_ancstr[i_part]);
    PDM_free(p_db_entity_path_itrf_strd[i_part]);
    PDM_free(p_db_entity_path_itrf[i_part]);
    PDM_free(p_db_entity_data[i_part]);
    PDM_free(p_db_entity_sens[i_part]);

    PDM_free(pextract_entity1_entity2_n      [i_part]);
    PDM_free(pextract_entity1_entity2_gnum   [i_part]);
    PDM_free(pextract_entity1_entity2_sign   [i_part]);
    PDM_free(pextract_entity1_entity2_triplet[i_part]);

    PDM_free(pextract_entity2_kind[i_part]);
    PDM_free(pextract_entity2_lnum[i_part]);
    PDM_free(pextract_entity2_itrf[i_part]);
    if (entity2_has_sens==1) {
      PDM_free(pextract_entity2_sens[i_part]);
    }

    PDM_free(pextract_entity2_twin_itrf_idx[i_part]);
    PDM_free(pextract_entity2_twin_itrf_gnum[i_part]);
    PDM_free(pextract_entity2_twin_itrf_itrf[i_part]);
    if (entity2_has_sens==1) {
      PDM_free(pextract_entity2_twin_itrf_sens[i_part]);
    }
    
    PDM_free(pextract_entity2_kind_idx [i_part]);
    PDM_free(pextract_entity2_kind_ordr[i_part]);
  }
  PDM_free(p_db_n_entity2);
  PDM_free(p_db_entity_ancstr);
  PDM_free(p_db_entity_path_itrf_strd);
  PDM_free(p_db_entity_path_itrf);
  PDM_free(p_db_entity_data);
  PDM_free(p_db_entity_gnum);
  PDM_free(p_db_entity_sens);

  PDM_free(pextract_entity1_entity2_n);
  PDM_free(pextract_entity1_entity2_gnum);
  PDM_free(pextract_entity1_entity2_sign);
  PDM_free(pextract_entity1_entity2_triplet);

  PDM_free(pextract_entity2_kind);
  PDM_free(pextract_entity2_lnum);
  PDM_free(pextract_entity2_itrf);
  if (entity2_has_sens==1) {
    PDM_free(pextract_entity2_sens);
  }

  PDM_free(pextract_entity2_twin_itrf_idx);
  PDM_free(pextract_entity2_twin_itrf_gnum);
  PDM_free(pextract_entity2_twin_itrf_itrf);
  PDM_free(pextract_entity2_twin_itrf_sens);
  
  PDM_free(pextract_entity2_kind_idx);
  PDM_free(pextract_entity2_kind_ordr);

  PDM_free(pextract_entity1_entity2_ntot);
  

  for (int i_part=0; i_part<n_part; ++i_part) {
    int i_write = 0;
    for (int i_entity=0; i_entity<pn_entity1_extended[i_part]; ++i_entity) {
      int i_first_entity = pentity1_extended_to_pentity1_idx[i_part][i_entity]/3;
      int i_beg = pextract_entity1_entity2_idx[i_part][i_first_entity  ];
      int i_end = pextract_entity1_entity2_idx[i_part][i_first_entity+1];
      int n_entity2 = i_end-i_beg;
      for (int i_entity2=i_beg; i_entity2<i_end; ++i_entity2) {
              pextended_entity1_entity2[i_part][i_write++] = pextended_entity1_entity2[i_part][i_entity2];
      }
      pextract_entity1_entity2_idx[i_part][i_entity+1] = pextract_entity1_entity2_idx[i_part][i_entity]+n_entity2;
    }
    int entity1_entity2_size = pextract_entity1_entity2_idx[i_part][pn_entity1_extended[i_part]];
    PDM_realloc(pextended_entity1_entity2   [i_part], pextended_entity1_entity2   [i_part],     entity1_entity2_size       , int);
    PDM_realloc(pextract_entity1_entity2_idx[i_part], pextract_entity1_entity2_idx[i_part], (pn_entity1_extended[i_part]+1), int);

  }


  // > Return
  *pn_entity2_extended_out                     = pn_entity2_extended;
  *pentity2_extended_ln_to_gn_out              = pextended_entity2_gnum;
  *pextended_entity1_entity2_idx_out           = pextract_entity1_entity2_idx;
  *pextended_entity1_entity2_out               = pextended_entity1_entity2;
  *pentity2_extended_to_pentity2_idx_out       = pextended_entity2_to_entity2_idx;
  *pentity2_extended_to_pentity2_triplet_out   = pextended_entity2_to_entity2_trplt;
  if (entity2_has_sens==1) {
    *pentity2_extended_to_pentity2_sens_out      = pextended_entity2_to_entity2_sens;
  }
  *pentity2_extended_to_pentity2_interface_out = pextended_entity2_to_entity2_itrf;

  return;
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
    PDM_free(interface_entity2_opp_sens[i_itrf]);
  }
  PDM_free(interface_entity2_opp_sens);


  int n_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    n_part_tot += n_part[i_dom];
  }

  for(int i_part = 0; i_part < n_part_tot; ++i_part) {
    PDM_free(is_entity2_on_itrf[i_part]);
  }
  PDM_free(is_entity2_on_itrf);

  PDM_domain_interface_free(ditrf);

  int          *pn_entity1           = NULL;
  PDM_g_num_t **pentity1_ln_to_gn    = NULL;
  int         **pentity1_entity2_idx = NULL;
  int         **pentity1_entity2     = NULL;
  PDM_malloc(pn_entity1          , n_part_tot, int          );
  PDM_malloc(pentity1_ln_to_gn   , n_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity1_entity2_idx, n_part_tot, int         *);
  PDM_malloc(pentity1_entity2    , n_part_tot, int         *);

  int          *pn_entity2        = NULL;
  PDM_g_num_t **pentity2_ln_to_gn = NULL;
  PDM_malloc(pn_entity2       , n_part_tot, int          );
  PDM_malloc(pentity2_ln_to_gn, n_part_tot, PDM_g_num_t *);

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

  int          *pn_entity2_opp_gnum_and_itrf = NULL;
  PDM_g_num_t **pentity2_opp_gnum_and_itrf   = NULL;
  int         **entity2_opp_position         = NULL;
  PDM_malloc(pn_entity2_opp_gnum_and_itrf, ln_part_tot, int          );
  PDM_malloc(pentity2_opp_gnum_and_itrf  , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(entity2_opp_position        , ln_part_tot, int         *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    // int *_pentity2_opp_location = pentity2_opp_location_idx[i_part];
    int *_pentity2_opp_location_idx = pentity2_opp_location_idx[i_part];
    int n_connect_tot = _pentity2_opp_location_idx[pn_entity2_num[i_part]];

    // On garde le lien courant (indirect sort of pentity2_num) + également keep le sens_opp
    PDM_malloc(pentity2_opp_gnum_and_itrf[i_part], 2 * n_connect_tot, PDM_g_num_t);
    PDM_g_num_t* _pentity2_opp_gnum_and_itrf = pentity2_opp_gnum_and_itrf[i_part];

    // log_trace("pn_entity2_num[i_part] = %i \n", pn_entity2_num[i_part]);
    // PDM_log_trace_array_int(_pentity2_opp_location_idx, pn_entity2_num[i_part]+1, "_pentity2_opp_location_idx ::");

    int *tmp_opp_position = NULL;
    PDM_malloc(tmp_opp_position, n_connect_tot, int);
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

    PDM_malloc(entity2_opp_position[i_part], n_connect_tot, int);

    int* order = NULL;
    PDM_malloc(order, n_connect_tot, int);
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


    PDM_free(order);
    PDM_free(tmp_opp_position);

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
  int         **gnum1_com_from_entity1_entity2_n       = NULL;
  PDM_g_num_t **gnum1_com_from_entity1_entity2         = NULL;
  int         **gnum1_com_from_entity1_entity2_triplet = NULL;
  PDM_malloc(gnum1_com_from_entity1_entity2_n      , n_part_tot, int         *);
  PDM_malloc(gnum1_com_from_entity1_entity2        , n_part_tot, PDM_g_num_t *);
  PDM_malloc(gnum1_com_from_entity1_entity2_triplet, n_part_tot, int         *);
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
    PDM_malloc(gnum1_com_from_entity1_entity2_n      [i_part],     n_gnum1_come_from, int        );
    PDM_malloc(gnum1_com_from_entity1_entity2        [i_part],     n_send_part2     , PDM_g_num_t);
    PDM_malloc(gnum1_com_from_entity1_entity2_triplet[i_part], 3 * n_send_part2     , int        );

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
    PDM_free(gnum1_com_from_entity1_entity2_n      [i_part]);
    PDM_free(gnum1_com_from_entity1_entity2        [i_part]);
    PDM_free(gnum1_com_from_entity1_entity2_triplet[i_part]);
  }
  PDM_free(gnum1_com_from_entity1_entity2_n      );
  PDM_free(gnum1_com_from_entity1_entity2        );
  PDM_free(gnum1_com_from_entity1_entity2_triplet);

  PDM_part_to_part_free(ptp);

  int **pextract_entity1_entity2_idx = NULL;
  PDM_malloc(pextract_entity1_entity2_idx, ln_part_tot, int *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    int n_part1_to_part2 = pentity1_extented_to_pentity1_idx[i_part][pn_entity1_extented[i_part]]/3;
    pextract_entity1_entity2_idx[i_part] = PDM_array_new_idx_from_sizes_int(pextract_entity1_entity2_n[i_part], n_part1_to_part2);
  }


  /*
   * Get for all receive connectivity information about sens and opp gnum
   * CAUTION we need to manage SENS and SIGN
   * We concatenate for all partition
   */
  PDM_g_num_t ***query_itrf_gnum    = NULL;
  int          **query_itrf_n       = NULL;
  int         ***recv_itrf_opp_n    = NULL;
  PDM_g_num_t ***recv_itrf_opp_gnum = NULL;
  PDM_malloc(query_itrf_gnum   , n_interface, PDM_g_num_t **);
  PDM_malloc(query_itrf_n      , n_interface, int          *);
  PDM_malloc(recv_itrf_opp_n   , n_interface, int         **);
  PDM_malloc(recv_itrf_opp_gnum, n_interface, PDM_g_num_t **);

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    PDM_malloc(query_itrf_gnum   [i_interface], ln_part_tot, PDM_g_num_t *);
    PDM_malloc(query_itrf_n      [i_interface], ln_part_tot, int          );
    PDM_malloc(recv_itrf_opp_n   [i_interface], ln_part_tot, int         *);
    PDM_malloc(recv_itrf_opp_gnum[i_interface], ln_part_tot, PDM_g_num_t *);

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      query_itrf_n[i_interface][i_part] = 0;
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
      PDM_malloc(query_itrf_gnum   [i_interface][i_part], query_itrf_n[i_interface][i_part], PDM_g_num_t);
      PDM_malloc(recv_itrf_opp_n   [i_interface][i_part], query_itrf_n[i_interface][i_part], int        );
      PDM_malloc(recv_itrf_opp_gnum[i_interface][i_part], query_itrf_n[i_interface][i_part], PDM_g_num_t);
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

  PDM_g_num_t ***recv_interface_entity2_opp_gnum   = NULL;
  int          **n_recv_interface_entity2_opp_gnum = NULL;
  PDM_malloc(recv_interface_entity2_opp_gnum  , n_interface, PDM_g_num_t **);
  PDM_malloc(n_recv_interface_entity2_opp_gnum, n_interface, int          *);
  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

    PDM_malloc(n_recv_interface_entity2_opp_gnum[i_interface], ln_part_tot, int);

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
    PDM_free(send_stride);

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(recv_stride[i_part], query_itrf_n[i_interface][i_part], "recv_stride ::");

      int n_recv_tot = 0;
      for(int i = 0; i < query_itrf_n[i_interface][i_part]; ++i) {
        n_recv_tot += recv_stride[i_part][i];
      }
      n_recv_interface_entity2_opp_gnum[i_interface][i_part] = n_recv_tot;
      PDM_free(recv_stride[i_part]);
      PDM_log_trace_array_long(recv_interface_entity2_opp_gnum[i_interface][i_part], n_recv_tot, "recv_interface_entity2_opp_gnum ::");
    }
    PDM_free(recv_stride);

    PDM_block_to_part_free(btp);
  }

  for(int i_itrf = 0; i_itrf < n_interface; ++i_itrf) {
    PDM_free(interface_entity2_opp_gnum[i_itrf]);
    PDM_part_to_block_free(ptb_interface_entity2[i_itrf]);
  }
  PDM_free(interface_entity2_opp_gnum);
  PDM_free(ptb_interface_entity2);

  int          *pn_entity2_ext_opp_gnum_and_itrf = NULL;
  PDM_g_num_t **pentity2_ext_opp_gnum_and_itrf   = NULL;
  PDM_malloc(pn_entity2_ext_opp_gnum_and_itrf, ln_part_tot, int          );
  PDM_malloc(pentity2_ext_opp_gnum_and_itrf  , ln_part_tot, PDM_g_num_t *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int n_recv_tot = 0;
    for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
      n_recv_tot += n_recv_interface_entity2_opp_gnum[i_interface][i_part];
    }
    PDM_malloc(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * n_recv_tot, PDM_g_num_t);

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

    int *order = NULL;
    PDM_malloc(order, n_recv_tot, int);
    pn_entity2_ext_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_recv_tot,
                                                                                   2,
                                                                                   pentity2_ext_opp_gnum_and_itrf[i_part],
                                                                                   order);
    PDM_realloc(pentity2_ext_opp_gnum_and_itrf[i_part], pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part], PDM_g_num_t);


    PDM_free(order);


    PDM_log_trace_array_long(pentity2_ext_opp_gnum_and_itrf[i_part], 2 * pn_entity2_ext_opp_gnum_and_itrf[i_part], "pentity2_opp_gnum_and_itrf ::");

  }


  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_free(query_itrf_gnum   [i_interface][i_part]);
      PDM_free(recv_itrf_opp_n   [i_interface][i_part]);
      PDM_free(recv_itrf_opp_gnum[i_interface][i_part]);
      PDM_free(recv_interface_entity2_opp_gnum[i_interface][i_part]);
    }
    PDM_free(query_itrf_gnum   [i_interface]);
    PDM_free(recv_itrf_opp_n   [i_interface]);
    PDM_free(recv_itrf_opp_gnum[i_interface]);
    PDM_free(query_itrf_n      [i_interface]);
    PDM_free(recv_interface_entity2_opp_gnum  [i_interface]);
    PDM_free(n_recv_interface_entity2_opp_gnum[i_interface]);
  }
  PDM_free(recv_interface_entity2_opp_gnum  );
  PDM_free(n_recv_interface_entity2_opp_gnum);

  PDM_free(query_itrf_gnum   );
  PDM_free(query_itrf_n      );
  PDM_free(recv_itrf_opp_n   );
  PDM_free(recv_itrf_opp_gnum);




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
  int         **entity2_order            = NULL;
  PDM_g_num_t **pentity2_ln_to_gn_sorted = NULL;
  PDM_malloc(entity2_order           , ln_part_tot, int         *);
  PDM_malloc(pentity2_ln_to_gn_sorted, ln_part_tot, PDM_g_num_t *);
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    /*
     * Sort current part entity2_ln_to_gn
     */
    PDM_malloc(entity2_order           [i_part], pn_entity2[i_part], int        );
    PDM_malloc(pentity2_ln_to_gn_sorted[i_part], pn_entity2[i_part], PDM_g_num_t);
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
  int          *pn_entity2_extented_by_interface         = NULL;
  int          *pn_entity2_extented_by_partition         = NULL;
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_interface  = NULL;
  PDM_g_num_t **pentity2_extented_ln_to_gn_by_partition  = NULL;
  int         **pentity2_extented_triplet_by_interface   = NULL;
  int         **pentity2_extented_triplet_by_partition   = NULL;
  int         **pentity2_extented_interface_by_interface = NULL;
  int         **pentity2_extented_interface_by_partition = NULL;
  PDM_malloc(pn_entity2_extented_by_interface        , ln_part_tot, int          );
  PDM_malloc(pn_entity2_extented_by_partition        , ln_part_tot, int          );
  PDM_malloc(pentity2_extented_ln_to_gn_by_interface , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity2_extented_ln_to_gn_by_partition , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity2_extented_triplet_by_interface  , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_triplet_by_partition  , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_interface_by_interface, ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_interface_by_partition, ln_part_tot, int         *);

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
    PDM_malloc(pentity2_extented_ln_to_gn_by_interface [i_part], 2 * pn_entity2_extented_by_interface[i_part], PDM_g_num_t);
    PDM_malloc(pentity2_extented_ln_to_gn_by_partition [i_part],     pn_entity2_extented_by_partition[i_part], PDM_g_num_t);
    PDM_malloc(pentity2_extented_triplet_by_interface  [i_part], 3 * pn_entity2_extented_by_interface[i_part], int        );
    PDM_malloc(pentity2_extented_triplet_by_partition  [i_part], 3 * pn_entity2_extented_by_partition[i_part], int        );
    PDM_malloc(pentity2_extented_interface_by_interface[i_part],     pn_entity2_extented_by_interface[i_part], int        );
    PDM_malloc(pentity2_extented_interface_by_partition[i_part],     pn_entity2_extented_by_partition[i_part], int        );

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
  int **pextract_entity1_entity2 = NULL;
  PDM_malloc(pextract_entity1_entity2, ln_part_tot, int *);
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


    PDM_malloc(pextract_entity1_entity2[i_part], _pextract_entity1_entity2_idx[n_part1_to_part2], int);
    int *_pextract_entity1_entity2  = pextract_entity1_entity2 [i_part];

    PDM_g_num_t* extented_from_itrf_entity2_ln_to_gn = PDM_gnum_get(gen_gnum_entity2, i_part);

    // Erase and realloc :
    PDM_realloc(pentity2_extented_ln_to_gn_by_interface[i_part], pentity2_extented_ln_to_gn_by_interface[i_part], pn_entity2_extented_by_interface[i_part], PDM_g_num_t);
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      pentity2_extented_ln_to_gn_by_interface[i_part][i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2] + shift_by_domain_entity2[n_domain];
    }

    // PDM_log_trace_array_long(extented_from_itrf_entity2_ln_to_gn, pn_entity2_extented_by_interface[i_part], "extented_from_itrf_entity2_ln_to_gn ::");

    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_itrf_entity2_ln_to_gn_sorted = NULL;
    PDM_malloc(extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part], PDM_g_num_t);
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_ln_to_gn_sorted[i_entity2] = extented_from_itrf_entity2_ln_to_gn[i_entity2];
    }
    // PDM_sort_long(extented_from_itrf_entity2_ln_to_gn_sorted, extented_from_itrf_entity2_order, pn_entity2_extented_by_interface[i_part]);
    pn_entity2_extented_by_interface[i_part] = PDM_inplace_unique_long(extented_from_itrf_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_interface[i_part]-1);
    PDM_realloc(extented_from_itrf_entity2_ln_to_gn_sorted, extented_from_itrf_entity2_ln_to_gn_sorted, pn_entity2_extented_by_interface[i_part], PDM_g_num_t);
    int *extented_from_itrf_entity2_order = NULL;
    PDM_malloc(extented_from_itrf_entity2_order, pn_entity2_extented_by_interface[i_part], int);
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_interface[i_part]; ++i_entity2) {
      extented_from_itrf_entity2_order[i_entity2] = -1;
    }

    printf("pn_entity2_extented_by_interface[i_part] = %i \n", pn_entity2_extented_by_interface[i_part]);


    /* Sort unique the new gnum to unify */
    PDM_g_num_t *extented_from_part_entity2_ln_to_gn_sorted = NULL;
    PDM_malloc(extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part], PDM_g_num_t);
    for(int i_entity2 = 0; i_entity2 < pn_entity2_extented_by_partition[i_part]; ++i_entity2) {
      extented_from_part_entity2_ln_to_gn_sorted[i_entity2] = pentity2_extented_ln_to_gn_by_partition[i_part][i_entity2];
      // extented_from_part_entity2_order          [i_entity2] = i_entity2;
    }
    // PDM_sort_long(extented_from_part_entity2_ln_to_gn_sorted, extented_from_part_entity2_order, pn_entity2_extented_by_partition[i_part]);

    pn_entity2_extented_by_partition[i_part] = PDM_inplace_unique_long(extented_from_part_entity2_ln_to_gn_sorted,
                                                                       NULL,
                                                                       0,
                                                                       pn_entity2_extented_by_partition[i_part]-1);
    int *extented_from_part_entity2_order = NULL;
    PDM_malloc(extented_from_part_entity2_order, pn_entity2_extented_by_partition[i_part], int);
    PDM_realloc(extented_from_part_entity2_ln_to_gn_sorted, extented_from_part_entity2_ln_to_gn_sorted, pn_entity2_extented_by_partition[i_part], PDM_g_num_t);
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

    PDM_free(extented_from_itrf_entity2_order          );
    PDM_free(extented_from_itrf_entity2_ln_to_gn_sorted);
    PDM_free(extented_from_part_entity2_order          );
    PDM_free(extented_from_part_entity2_ln_to_gn_sorted);

  }


  PDM_gnum_free(gen_gnum_entity2);

  /*
   * Reconstruction
   */
  int          *pn_entity2_extented              = NULL;
  PDM_g_num_t **pentity2_extented_ln_to_gn       = NULL;
  int         **pentity2_extented_triplet        = NULL;
  int         **pentity2_extented_interface      = NULL;
  int         **pentity2_extented_to_entity2_idx = NULL;
  PDM_malloc(pn_entity2_extented             , ln_part_tot, int          );
  PDM_malloc(pentity2_extented_ln_to_gn      , ln_part_tot, PDM_g_num_t *);
  PDM_malloc(pentity2_extented_triplet       , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_interface     , ln_part_tot, int         *);
  PDM_malloc(pentity2_extented_to_entity2_idx, ln_part_tot, int         *);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    pn_entity2_extented[i_part] = pn_entity2_extented_by_partition[i_part] + pn_entity2_extented_by_interface[i_part];
    PDM_malloc(pentity2_extented_ln_to_gn      [i_part],     pn_entity2_extented[i_part]    , PDM_g_num_t);
    PDM_malloc(pentity2_extented_triplet       [i_part], 3 * pn_entity2_extented[i_part]    , int        );
    PDM_malloc(pentity2_extented_interface     [i_part],     pn_entity2_extented[i_part]    , int        );
    PDM_malloc(pentity2_extented_to_entity2_idx[i_part],     pn_entity2_extented[i_part] + 1, int        );

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
  //   PDM_free(pextract_entity1_entity2[i_part]);
  // }
  // PDM_free(pextract_entity1_entity2);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pextract_entity1_entity2_n      [i_part]);
    PDM_free(pextract_entity1_entity2_gnum   [i_part]);
    PDM_free(pextract_entity1_entity2_triplet[i_part]);
    // PDM_free(pextract_entity1_entity2_idx [i_part]);
  }
  PDM_free(pextract_entity1_entity2_n      );
  PDM_free(pextract_entity1_entity2_gnum   );
  PDM_free(pextract_entity1_entity2_triplet);
  // PDM_free(pextract_entity1_entity2_idx );


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pentity2_extented_ln_to_gn_by_interface [i_part]);
    PDM_free(pentity2_extented_ln_to_gn_by_partition [i_part]);
    PDM_free(pentity2_extented_triplet_by_interface  [i_part]);
    PDM_free(pentity2_extented_triplet_by_partition  [i_part]);
    PDM_free(pentity2_extented_interface_by_interface[i_part]);
    PDM_free(pentity2_extented_interface_by_partition[i_part]);
  }
  PDM_free(pentity2_extented_ln_to_gn_by_interface );
  PDM_free(pentity2_extented_ln_to_gn_by_partition );
  PDM_free(pentity2_extented_triplet_by_interface  );
  PDM_free(pentity2_extented_triplet_by_partition  );
  PDM_free(pentity2_extented_interface_by_interface);
  PDM_free(pentity2_extented_interface_by_partition);

  PDM_free(pn_entity2_extented_by_interface );
  PDM_free(pn_entity2_extented_by_partition );

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(entity2_order           [i_part]);
    PDM_free(pentity2_ln_to_gn_sorted[i_part]);
  }
  PDM_free(entity2_order           );
  PDM_free(pentity2_ln_to_gn_sorted);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pentity2_opp_gnum_and_itrf[i_part]);
    PDM_free(pentity2_ext_opp_gnum_and_itrf[i_part]);
    PDM_free(entity2_opp_position      [i_part]);
  }
  PDM_free(entity2_opp_position        );
  PDM_free(pentity2_opp_gnum_and_itrf  );
  PDM_free(pentity2_ext_opp_gnum_and_itrf  );
  PDM_free(pn_entity2_opp_gnum_and_itrf);
  PDM_free(pn_entity2_ext_opp_gnum_and_itrf);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_free(pentity2_num             [i_part]);
    PDM_free(pentity2_opp_location_idx[i_part]);
    PDM_free(pentity2_opp_location    [i_part]);
    PDM_free(pentity2_opp_interface   [i_part]);
    PDM_free(pentity2_opp_sens        [i_part]);
    PDM_free(pentity2_opp_gnum        [i_part]);
  }
  PDM_free(pn_entity2_num            );
  PDM_free(pentity2_num              );
  PDM_free(pentity2_opp_location_idx );
  PDM_free(pentity2_opp_location     );
  PDM_free(pentity2_opp_interface_idx);
  PDM_free(pentity2_opp_interface    );
  PDM_free(pentity2_opp_sens         );
  PDM_free(pentity2_opp_gnum         );

  PDM_free(pn_entity1          );
  PDM_free(pentity1_ln_to_gn   );
  PDM_free(pentity1_entity2_idx);
  PDM_free(pentity1_entity2    );
  PDM_free(pn_entity2          );
  PDM_free(pentity2_ln_to_gn   );

}



#ifdef __cplusplus
}
#endif /* __cplusplus */

