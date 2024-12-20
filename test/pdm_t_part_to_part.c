#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_part_to_part.h"
#include "pdm_array.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  assert(n_rank == 2);

  int n_part1 = 3;
  int n_part2 = 3;

  int *n_elt1;
  PDM_malloc(n_elt1,n_part1,int);
  PDM_g_num_t **gnum_elt1;
  PDM_malloc(gnum_elt1,n_part1,PDM_g_num_t *);
  int **part1_to_part2_idx;
  PDM_malloc(part1_to_part2_idx,n_part1,int *);
  PDM_g_num_t **part1_to_part2;
  PDM_malloc(part1_to_part2,n_part1,PDM_g_num_t *);

  int *n_elt2 = n_elt1;
  PDM_g_num_t **gnum_elt2 = gnum_elt1;


  if (i_rank == 0) {
    n_elt1[0] = 1;
    n_elt1[1] = 4;
    n_elt1[2] = 2;
  } else {
    n_elt1[0] = 2;
    n_elt1[1] = 6;
    n_elt1[2] = 2;
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_malloc(gnum_elt1[i],n_elt1[i],PDM_g_num_t);
    PDM_malloc(part1_to_part2_idx[i],(n_elt1[i] + 1),int);
  }


  if (i_rank == 0) {
    gnum_elt1[0][0] = 9;

    gnum_elt1[1][0] = 15;
    gnum_elt1[1][1] = 13;
    gnum_elt1[1][2] = 5;
    gnum_elt1[1][3] = 2;

    gnum_elt1[2][0] = 3;
    gnum_elt1[2][1] = 8;

    part1_to_part2_idx[0][0] = 0;
    part1_to_part2_idx[0][1] = 0;

    part1_to_part2_idx[1][0] = 0;
    part1_to_part2_idx[1][1] = 2;
    part1_to_part2_idx[1][2] = 4;
    part1_to_part2_idx[1][3] = 5;
    part1_to_part2_idx[1][4] = 6;

    part1_to_part2_idx[2][0] = 0;
    part1_to_part2_idx[2][1] = 1;
    part1_to_part2_idx[2][2] = 1;
  }

  else {
    gnum_elt1[0][0] = 12;
    gnum_elt1[0][1] = 16;

    gnum_elt1[1][0] = 11;
    gnum_elt1[1][1] = 10;
    gnum_elt1[1][2] = 17;
    gnum_elt1[1][3] = 14;
    gnum_elt1[1][4] = 1;
    gnum_elt1[1][5] = 4;

    gnum_elt1[2][0] = 7;
    gnum_elt1[2][1] = 6;

    part1_to_part2_idx[0][0] = 0;
    part1_to_part2_idx[0][1] = 0;
    part1_to_part2_idx[0][2] = 0;

    part1_to_part2_idx[1][0] = 0;
    part1_to_part2_idx[1][1] = 0;
    part1_to_part2_idx[1][2] = 0;
    part1_to_part2_idx[1][3] = 2;
    part1_to_part2_idx[1][4] = 4;
    part1_to_part2_idx[1][5] = 6;
    part1_to_part2_idx[1][6] = 8;

    part1_to_part2_idx[2][0] = 0;
    part1_to_part2_idx[2][1] = 1;
    part1_to_part2_idx[2][2] = 2;
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_malloc(part1_to_part2[i],part1_to_part2_idx[i][n_elt1[i]],PDM_g_num_t);
  }

  if (i_rank == 0) {
    // part1_to_part2[0][0] = ;

    part1_to_part2[1][0] = 9;
    part1_to_part2[1][1] = 11;
    part1_to_part2[1][2] = 9;
    part1_to_part2[1][3] = 10;
    part1_to_part2[1][4] = 11;
    part1_to_part2[1][5] = 10;

    part1_to_part2[2][0] = 9;
  }

  else {
    // part1_to_part2[0][0] = ;

    part1_to_part2[1][0] = 12;
    part1_to_part2[1][1] = 11;
    part1_to_part2[1][2] = 12;
    part1_to_part2[1][3] = 10;
    part1_to_part2[1][4] = 16;
    part1_to_part2[1][5] = 11;
    part1_to_part2[1][6] = 16;
    part1_to_part2[1][7] = 10;

    part1_to_part2[2][0] = 12;
    part1_to_part2[2][1] = 16;
  }

  /*
   *  Create Part-to-part object
   */
  PDM_part_to_part_t *ptp = PDM_part_to_part_create ((const PDM_g_num_t **) gnum_elt1,
                                                     n_elt1,
                                                     n_part1,
                                                     (const PDM_g_num_t **)gnum_elt2,
                                                     n_elt2,
                                                     n_part2,
                                                     (const int **) part1_to_part2_idx,
                                                     (const PDM_g_num_t **) part1_to_part2,
                                                     comm);
  /*
   *  Check Part2 ref/unref/gnum1_come_from
   */
  int  *n_ref_num2 = NULL;
  int **ref_num2   = NULL;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_ref_num2,
                                  &ref_num2);

  int  *n_unref_num2 = NULL;
  int **unref_num2   = NULL;
  PDM_part_to_part_ref_lnum2_get (ptp,
                                  &n_unref_num2,
                                  &unref_num2);

  int **gnum1_come_from_idx = NULL;
  PDM_g_num_t **gnum1_come_from = NULL;
  PDM_part_to_part_gnum1_come_from_get (ptp,
                                        &gnum1_come_from_idx,
                                        &gnum1_come_from);
  // for (int i = 0; i < n_part2; i++) {

  //   log_trace("\npart2 %d\n", i);
  //   PDM_log_trace_array_int(ref_num2[i], n_ref_num2[i], "referenced (l_num) : ");
  //   log_trace("unreferenced (g_num) : ");
  //   for (int j = 0; j < n_unref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" ", gnum_elt2[i][unref_num2[i][j]-1]);
  //   }
  //   log_trace("\n");
  //   log_trace("referenced (g_num) -> gnum1_come_from :\n");
  //   for (int j = 0; j < n_ref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" -> ", gnum_elt2[i][ref_num2[i][j]-1]);
  //     for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
  //       log_trace(PDM_FMT_G_NUM" ", gnum1_come_from[i][k]);
  //     }
  //     log_trace("\n");
  //   }


  // }


  int **part1_stride;
  PDM_malloc(part1_stride,n_part1,int         *);
  PDM_g_num_t **part1_data;
  PDM_malloc(part1_data,n_part1,PDM_g_num_t *);

  for (int i = 0; i < n_part1; i++) {

    // log_trace("\npart1 %d\n", i);

    PDM_malloc(part1_stride[i],n_elt1[i],int);

    int s_part1_data = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      part1_stride[i][j] = (int) gnum_elt1[i][j];
      // part1_stride[i][j] = 1;
      s_part1_data += part1_stride[i][j];
    }
    // PDM_log_trace_array_int(part1_stride[i], n_elt1[i], "part1_stride : ");

    // log_trace("g_num -> data:\n");
    PDM_malloc(part1_data[i],s_part1_data,PDM_g_num_t);
    int idx = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      // int idx0 = idx;
      for (int k = 1; k <= part1_stride[i][j]; k++) {
        part1_data[i][idx++] = k;
      }
      // log_trace(PDM_FMT_G_NUM, gnum_elt1[i][j]);
      // PDM_log_trace_array_long(part1_data[i] + idx0,
      //                          part1_stride[i][j],
      //                          " -> ");
    }
    // for (int j = 0; j < n_elt1[i]; j++) {
    //   part1_data[i][j] = gnum_elt1[i][j];
    // }
  }


  int         **part2_stride = NULL;
  PDM_g_num_t **part2_data   = NULL;
  int request;

  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_VAR_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                          0,
                          sizeof(PDM_g_num_t),
           (const int **) part1_stride,
          (const void **) part1_data,
                          &part2_stride,
               (void ***) &part2_data,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);

  // log_trace("\n\n---- Check iexch ----\n");
  // for (int i = 0; i < n_part2; i++) {

  //   log_trace("\npart2 %d\n", i);
  //   PDM_log_trace_array_int(part2_stride[i], gnum1_come_from_idx[i][n_ref_num2[i]], "stride : ");
  //   log_trace("referenced (g_num) -> data :\n");
  //   int idx = 0;
  //   for (int j = 0; j < n_ref_num2[i]; j++) {
  //     log_trace(PDM_FMT_G_NUM" :\n", gnum_elt2[i][ref_num2[i][j]-1]);
  //     for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
  //       log_trace("    "PDM_FMT_G_NUM" -> ", gnum1_come_from[i][k]);
  //       for (int l = 0; l < part2_stride[i][k]; l++) {
  //         log_trace(PDM_FMT_G_NUM" ", part2_data[i][idx++]);
  //       }
  //       log_trace("\n");
  //     }
  //     log_trace("\n");
  //   }


  // }




  /*
   *  Exchange an interleaved, constant-stride field
   */
  // log_trace("\n\n---- Exchange an interleaved, constant-stride field ----\n");
  PDM_g_num_t **part1_field;
  PDM_malloc(part1_field,n_part1,PDM_g_num_t *);
  for (int i = 0; i < n_part1; i++) {
    // int n = part1_to_part2_idx[i][n_elt1[i]];
    // PDM_malloc(part1_field[i],n * 2,PDM_g_num_t);

    // for (int j = 0; j < n_elt1[i]; j++) {
    //   for (int k = part1_to_part2_idx[i][j]; k < part1_to_part2_idx[i][j+1]; k++) {
    //     part1_field[i][k  ] = gnum_elt1[i][j];
    //     part1_field[i][k+n] = part1_to_part2[i][k];
    //   }
    // }
    int n = n_elt1[i];
    PDM_malloc(part1_field[i],n * 2,PDM_g_num_t);

    for (int j = 0; j < n_elt1[i]; j++) {
      part1_field[i][j  ] = gnum_elt1[i][j];
      part1_field[i][j+n] = gnum_elt1[i][j]+1;
    }

    // log_trace("\npart1 %d\n", i);
    // PDM_log_trace_array_long(part1_field[i],     n, "  part1_field (1st comp) : ");
    // PDM_log_trace_array_long(part1_field[i] + n, n, "  part1_field (2nd comp) : ");
  }



  PDM_g_num_t **part2_field = NULL;
  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_CST_INTERLEAVED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,//PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                          2,
                          sizeof(PDM_g_num_t),
                          NULL,
          (const void **) part1_field,
                          NULL,
               (void ***) &part2_field,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);


  // for (int i = 0; i < n_part2; i++) {
  //   log_trace("\npart2 %d\n", i);
  //   int n = gnum1_come_from_idx[i][n_ref_num2[i]];
  //   PDM_log_trace_array_long(part2_field[i],     n, "  part2_field (1st comp) : ");
  //   PDM_log_trace_array_long(part2_field[i] + n, n, "  part2_field (2nd comp) : ");
  // }



  /* Reverse */
  for (int i = 0; i < n_part1; i++) {
    PDM_free(part1_field[i]);
  }
  PDM_free(part1_field);

  PDM_part_to_part_reverse_iexch (ptp,
                                  PDM_MPI_COMM_KIND_P2P,
                                  PDM_STRIDE_CST_INTERLEAVED,
                                  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                  2,
                                  sizeof(PDM_g_num_t),
                                  NULL,
                  (const void **) part2_field,
                                  NULL,
                       (void ***) &part1_field,
                                  &request);

  PDM_part_to_part_reverse_iexch_wait (ptp,
                                       request);

  // log_trace("Reverse\n");
  // for (int i = 0; i < n_part1; i++) {
  //   int n = part1_to_part2_idx[i][n_elt1[i]];
  //   log_trace("\npart1 %d\n", i);
  //   PDM_log_trace_array_long(part1_field[i],     n, "  part1_field (1st comp) : ");
  //   PDM_log_trace_array_long(part1_field[i] + n, n, "  part1_field (2nd comp) : ");
  // }


  /*
   *  Exchange an interlaced, constant-stride field
   */

  for (int i = 0; i < n_part1; i++) {
    PDM_free(part1_field[i]);
  }
  PDM_free(part1_field);

  PDM_malloc(part1_field,n_part1,PDM_g_num_t *);
  for (int i = 0; i < n_part1; i++) {
    int n = n_elt1[i];
    PDM_malloc(part1_field[i],n * 2,PDM_g_num_t);
    for (int j = 0; j < n_elt1[i]; j++) {
      part1_field[i][2*j  ] = gnum_elt1[i][j];
      part1_field[i][2*j+1] = gnum_elt1[i][j]+1;
    }
  }

  PDM_g_num_t **part2_field_p2p = NULL;
  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_P2P,
                          PDM_STRIDE_CST_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                          2,
                          sizeof(PDM_g_num_t),
                          NULL,
          (const void **) part1_field,
                          NULL,
               (void ***) &part2_field_p2p,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);

  /* Reverse */
  PDM_g_num_t **part1_field_p2p = NULL;
  PDM_part_to_part_reverse_iexch (ptp,
                                  PDM_MPI_COMM_KIND_P2P,
                                  PDM_STRIDE_CST_INTERLACED,
                                  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                  2,
                                  sizeof(PDM_g_num_t),
                                  NULL,
                  (const void **) part2_field_p2p,
                                  NULL,
                       (void ***) &part1_field_p2p,
                                  &request);

  PDM_part_to_part_reverse_iexch_wait (ptp,
                                       request);

  PDM_g_num_t **part2_field_coll = NULL;
  PDM_part_to_part_iexch (ptp,
                          PDM_MPI_COMM_KIND_COLLECTIVE,
                          PDM_STRIDE_CST_INTERLACED,
                          PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                          2,
                          sizeof(PDM_g_num_t),
                          NULL,
          (const void **) part1_field,
                          NULL,
               (void ***) &part2_field_coll,
                          &request);

  PDM_part_to_part_iexch_wait (ptp,
                               request);

  /* Reverse */
  PDM_g_num_t **part1_field_coll = NULL;
  PDM_part_to_part_reverse_iexch (ptp,
                                  PDM_MPI_COMM_KIND_COLLECTIVE,
                                  PDM_STRIDE_CST_INTERLACED,
                                  PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                  2,
                                  sizeof(PDM_g_num_t),
                                  NULL,
                  (const void **) part2_field_coll,
                                  NULL,
                       (void ***) &part1_field_coll,
                                  &request);

  PDM_part_to_part_reverse_iexch_wait (ptp,
                                       request);

  for (int i = 0; i < n_part2; i++) {
    for (int j = 0; j < n_ref_num2[i]; j++) {
      for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
        assert(part2_field_p2p[i][2*k  ] == part2_field_coll[i][2*k  ]);
        assert(part2_field_p2p[i][2*k+1] == part2_field_coll[i][2*k+1]);
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    for (int j = 0; j < n_elt1[i]; j++) {
      for (int k = part1_to_part2_idx[i][j]; k < part1_to_part2_idx[i][j+1]; k++) {
        assert(part1_field_p2p[i][2*k  ] == part1_field_coll[i][2*k  ]);
        assert(part1_field_p2p[i][2*k+1] == part1_field_coll[i][2*k+1]);
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_free(part1_field_p2p[i]);
  }
  PDM_free(part1_field_p2p);

  for (int i = 0; i < n_part2; i++) {
    PDM_free(part2_field_p2p[i]);
  }
  PDM_free(part2_field_p2p);

  for (int i = 0; i < n_part1; i++) {
    PDM_free(part1_field_coll[i]);
  }
  PDM_free(part1_field_coll);

  for (int i = 0; i < n_part2; i++) {
    PDM_free(part2_field_coll[i]);
  }
  PDM_free(part2_field_coll);

  // log_trace("==== P1 -> P2 ====\n");
  /* 2 consecutive iexch in stride var with same stride */
  for (int ipart = 0; ipart < n_part1; ipart++) {
    int s_part1_data = 0;
    PDM_realloc(part1_stride[ipart] ,part1_stride[ipart] , n_elt1[ipart],int);
    for (int i = 0; i < n_elt1[ipart]; i++) {
      part1_stride[ipart][i] = (int) (gnum_elt1[ipart][i] % 2) + 1;
      s_part1_data += part1_stride[ipart][i];
    }

    PDM_realloc(part1_data[ipart] ,part1_data[ipart] , s_part1_data,PDM_g_num_t);
    int idx = 0;
    for (int i = 0; i < n_elt1[ipart]; i++) {
      for (int j = 0; j < part1_stride[ipart][i]; j++) {
        part1_data[ipart][idx++] = gnum_elt1[ipart][i];
      }
    }
  }

  for (int i = 0; i < n_part2; i++) {
    PDM_free(part2_stride[i]);
    PDM_free(part2_data  [i]);
  }
  PDM_free(part2_stride);
  PDM_free(part2_data);


  part2_stride = NULL;
  // log_trace("1\n");
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
         (const int  **) part1_stride,
         (const void **) part1_data,
                         &part2_stride,
              (void ***) &part2_data,
                         &request);
  PDM_part_to_part_iexch_wait (ptp, request);

  // log_trace("2\n");
  PDM_g_num_t **part2_data2 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
         (const int  **) part1_stride,
         (const void **) part1_data,
                         &part2_stride,
              (void ***) &part2_data2,
                         &request);
  PDM_part_to_part_iexch_wait (ptp, request);

  for (int i = 0; i < n_part2; i++) {
    int idx = 0;
    for (int j = 0; j < n_ref_num2[i]; j++) {
      // int id2 = ref_num2[i][j] - 1;
      for (int k = gnum1_come_from_idx[i][j]; k < gnum1_come_from_idx[i][j+1]; k++) {
        // log_trace("gnum2 "PDM_FMT_G_NUM", gnum1 "PDM_FMT_G_NUM", expected stride = %d, got %d\n",
        //           gnum_elt2[i][id2], gnum1_come_from[i][k], (int) (gnum1_come_from[i][k]%2) + 1,
        //           part2_stride[i][k]);
        for (int l = 0; l < part2_stride[i][k]; l++) {
          // log_trace("  "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
          //           part2_data2[i][idx], part2_data[i][idx]);
          assert(part2_data2[i][idx] == part2_data[i][idx]);
          idx++;
        }
      }
    }

    PDM_free(part2_data2[i]);
  }
  PDM_free(part2_data2);



  // log_trace("==== P1 <- P2 ====\n");
  /* 2 consecutive reverse iexch in stride var with same stride */
  for (int ipart = 0; ipart < n_part2; ipart++) {
    int s_part2_data = 0;
    PDM_realloc(part2_stride[ipart] ,part2_stride[ipart] , n_elt2[ipart],int);
    for (int i = 0; i < n_elt2[ipart]; i++) {
      part2_stride[ipart][i] = (int) (gnum_elt2[ipart][i] % 2) + 1;
      s_part2_data += part2_stride[ipart][i];
    }

    PDM_realloc(part2_data[ipart] ,part2_data[ipart] , s_part2_data,PDM_g_num_t);
    int idx = 0;
    for (int i = 0; i < n_elt2[ipart]; i++) {
      for (int j = 0; j < part2_stride[ipart][i]; j++) {
        part2_data[ipart][idx++] = gnum_elt2[ipart][i];
      }
    }
  }

  for (int i = 0; i < n_part1; i++) {
    PDM_free(part1_stride[i]);
    PDM_free(part1_data  [i]);
  }
  PDM_free(part1_stride);
  PDM_free(part1_data);

  part1_stride = NULL;
  // log_trace("1\n");
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                 (const int  **) part2_stride,
                 (const void **) part2_data,
                                 &part1_stride,
                      (void ***) &part1_data,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait (ptp, request);
  // for (int i = 0; i < n_part1; i++) {
  //   log_trace("part1 %d\n", i);
  //   PDM_log_trace_array_int(part1_stride[i], part1_to_part2_idx[i][n_elt1[i]], "part1_stride : ");
  // }
  // log_trace("2\n");

  // for (int i = 0; i < n_part1; i++) {
  //  PDM_free(part1_stride[i]);
  // }
  //PDM_free(part1_stride);
  // part1_stride = NULL;

  PDM_g_num_t **part1_data2 = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_VAR_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                 (const int  **) part2_stride,
                 (const void **) part2_data,
                                 &part1_stride,
                      (void ***) &part1_data2,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait (ptp, request);
  // for (int i = 0; i < n_part1; i++) {
  //   log_trace("part1 %d\n", i);
  //   PDM_log_trace_array_int(part1_stride[i], part1_to_part2_idx[i][n_elt1[i]], "part1_stride : ");
  // }

  for (int i = 0; i < n_part1; i++) {
    int idx = 0;
    for (int j = 0; j < n_elt1[i]; j++) {
      for (int k = part1_to_part2_idx[i][j]; k < part1_to_part2_idx[i][j+1]; k++) {
        // log_trace("gnum1 "PDM_FMT_G_NUM", gnum2 "PDM_FMT_G_NUM", expected stride = %d, got %d\n",
        //           gnum_elt1[i][j], part1_to_part2[i][k], (int) (part1_to_part2[i][k]%2) + 1,
        //           part1_stride[i][k]);
        for (int l = 0; l < part1_stride[i][k]; l++) {
          // log_trace("  "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n",
          //           part1_data2[i][idx], part1_data[i][idx]);
          assert(part1_data2[i][idx] == part1_data[i][idx]);
          idx++;
        }
      }
    }

    PDM_free(part1_data2[i]);
  }
  PDM_free(part1_data2);

  /*
   *  Free memory
   */
  PDM_part_to_part_free (ptp);

  for (int i = 0; i < n_part1; i++) {
    PDM_free(gnum_elt1[i]);
    PDM_free(part1_to_part2_idx[i]);
    PDM_free(part1_to_part2[i]);

    PDM_free(part1_stride[i]);
    PDM_free(part1_data[i]);
    PDM_free(part1_field[i]);
  }

  for (int i = 0; i < n_part2; i++) {
    PDM_free(part2_stride[i]);
    PDM_free(part2_data[i]);
    PDM_free(part2_field[i]);
  }

  PDM_free(n_elt1);
  PDM_free(gnum_elt1);
  PDM_free(part1_to_part2_idx);
  PDM_free(part1_to_part2);

  PDM_free(part1_stride);
  PDM_free(part1_data);
  PDM_free(part1_field);
  PDM_free(part2_stride);
  PDM_free(part2_data);
  PDM_free(part2_field);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
