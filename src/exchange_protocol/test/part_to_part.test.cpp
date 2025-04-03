#include <stddef.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_part_to_part.h"




MPI_TEST_CASE("[pdm_part_to_part] - 2p - part_to_part", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  /*
   *                   p0                              p1
   *
   *                       2        3        3        1
   * part1 =               x--------x        x--------x
   *
   *
   * part2 =      x--------x--------x        x--------x--------x
   *              10       20       30       30       40       50
   */
  std::vector<std::vector<PDM_g_num_t>> vgnum_elt1 = {{2, 3}, {3, 1}};
  int n_part1 = 1;

  std::vector<std::vector<PDM_g_num_t>> vgnum_elt2 = {{10, 20, 30}, {30, 40, 50}};
  int n_part2 = 1;

  std::vector<std::vector<int>>         vpart1_to_part2_idx = {{0 , 1, 3  }, {0, 2, 4}};
  std::vector<std::vector<PDM_g_num_t>> vpart1_to_part2     = {{20, 30, 30}, {30, 30, 40, 30}};


  int n_elt1 = vgnum_elt1[i_rank].size();
  int n_elt2 = vgnum_elt2[i_rank].size();

  PDM_g_num_t *gnum_elt1 = vgnum_elt1[i_rank].data();
  PDM_g_num_t *gnum_elt2 = vgnum_elt2[i_rank].data();

  int         *part1_to_part2_idx = vpart1_to_part2_idx[i_rank].data();
  PDM_g_num_t *part1_to_part2     = vpart1_to_part2    [i_rank].data();
  PDM_part_to_part_t* ptp = PDM_part_to_part_create((const PDM_g_num_t **) &gnum_elt1,
                                                    (const int          *) &n_elt1,
                                                                           n_part1,
                                                    (const PDM_g_num_t **) &gnum_elt2,
                                                    (const int          *) &n_elt2,
                                                                           n_part2,
                                                    (const int         **) &part1_to_part2_idx,
                                                    (const PDM_g_num_t **) &part1_to_part2,
                                                                           pdm_comm);

  int  *pn_ref_lnum2 = NULL;
  int **pref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp,
                                 &pn_ref_lnum2,
                                 &pref_lnum2);

  int  n_ref_lnum2 = pn_ref_lnum2[0];
  int *ref_lnum2   = pref_lnum2  [0];

  int  *pn_unref_lnum2 = NULL;
  int **punref_lnum2   = NULL;
  PDM_part_to_part_unref_lnum2_get(ptp,
                                   &pn_unref_lnum2,
                                   &punref_lnum2);
  int  n_unref_lnum2 = pn_unref_lnum2[0];
  int *unref_lnum2   = punref_lnum2  [0];

  int         **pgnum1_come_from_idx = NULL;
  PDM_g_num_t **pgnum1_come_from     = NULL;
  PDM_part_to_part_gnum1_come_from_get(ptp, &pgnum1_come_from_idx, &pgnum1_come_from);
  int         *gnum1_come_from_idx = pgnum1_come_from_idx[0];
  PDM_g_num_t *gnum1_come_from     = pgnum1_come_from    [0];

  if(0 == 1) {
    PDM_log_trace_array_int (ref_lnum2  , n_ref_lnum2  , "ref_lnum2   ::");
    PDM_log_trace_array_int (unref_lnum2, n_unref_lnum2, "unref_lnum2 ::");
    PDM_log_trace_array_int (gnum1_come_from_idx, n_ref_lnum2+1, "gnum1_come_from_idx ::");
    PDM_log_trace_array_long(gnum1_come_from, gnum1_come_from_idx[n_ref_lnum2], "gnum1_come_from ::");
  }

  /*
   * Check
   */
  int         p0_expected_ref_lnum2          [2] = {2, 3};
  int         p0_expected_unref_lnum2        [1] = {1};
  int         p0_expected_gnum1_come_from_idx[3] = {0, 1, 3};
  PDM_g_num_t p0_expected_gnum1_come_from    [3] = {2, 1, 3};

  MPI_CHECK_EQ_C_ARRAY(0, ref_lnum2          , p0_expected_ref_lnum2          , 2);
  MPI_CHECK_EQ_C_ARRAY(0, unref_lnum2        , p0_expected_unref_lnum2        , 1);
  MPI_CHECK_EQ_C_ARRAY(0, gnum1_come_from_idx, p0_expected_gnum1_come_from_idx, 3);
  MPI_CHECK_EQ_C_ARRAY(0, gnum1_come_from    , p0_expected_gnum1_come_from    , 3);

  int         p1_expected_ref_lnum2          [2] = {1, 2};
  int         p1_expected_unref_lnum2        [1] = {3};
  int         p1_expected_gnum1_come_from_idx[3] = {0, 2, 3};
  PDM_g_num_t p1_expected_gnum1_come_from    [3] = {1, 3, 1};

  MPI_CHECK_EQ_C_ARRAY(1, ref_lnum2          , p1_expected_ref_lnum2          , 2);
  MPI_CHECK_EQ_C_ARRAY(1, unref_lnum2        , p1_expected_unref_lnum2        , 1);
  MPI_CHECK_EQ_C_ARRAY(1, gnum1_come_from_idx, p1_expected_gnum1_come_from_idx, 3);
  MPI_CHECK_EQ_C_ARRAY(1, gnum1_come_from    , p1_expected_gnum1_come_from    , 3);

  /*
   * ************************************************************************************************
   * 1/ Exchange gnum1 to part1->part2 with ORDER_PART1
   */
  PDM_g_num_t **tmp_recv_part2_to_part1_gnum_elt2 = NULL;
  int request = -1;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
                         NULL,
       (const void **)   &gnum_elt1,
                         NULL,
           (void ***)    &tmp_recv_part2_to_part1_gnum_elt2,
                         &request);
  PDM_part_to_part_iexch_wait(ptp, request);

  PDM_g_num_t *recv_part2_to_part1_gnum_elt2 = tmp_recv_part2_to_part1_gnum_elt2[0];
  PDM_free(tmp_recv_part2_to_part1_gnum_elt2);

  if(0 == 1) {
    PDM_log_trace_array_long(recv_part2_to_part1_gnum_elt2, gnum1_come_from_idx[n_ref_lnum2], "recv_part2_to_part1_gnum_elt2 ::");
  }

  /*
   *  En tout rigeur le 3 pourrait venir du rang 0 et du rang 1 mais en interne on prends le premier venu --> See pdm_part_migrate
   */
  PDM_g_num_t p0_expected_recv_part2_to_part1_gnum_elt2[3] = {2, 1, 3};
  PDM_g_num_t p1_expected_recv_part2_to_part1_gnum_elt2[3] = {1, 3, 1};

  MPI_CHECK_EQ_C_ARRAY(0, recv_part2_to_part1_gnum_elt2, p0_expected_recv_part2_to_part1_gnum_elt2, 3);
  MPI_CHECK_EQ_C_ARRAY(1, recv_part2_to_part1_gnum_elt2, p1_expected_recv_part2_to_part1_gnum_elt2, 3);

  /*
   * ************************************************************************************************
   * 2/ Exchange gnum1 to part2->part1 with ORDER_PART2
   */
  PDM_g_num_t **tmp_recv_part1_to_part2_gnum_elt2 = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(PDM_g_num_t),
                                 NULL,
               (const void **)   &gnum_elt2,
                                 NULL,
                   (void ***)    &tmp_recv_part1_to_part2_gnum_elt2,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait(ptp, request);

  PDM_g_num_t *recv_part1_to_part2_gnum_elt2 = tmp_recv_part1_to_part2_gnum_elt2[0];
  PDM_free(tmp_recv_part1_to_part2_gnum_elt2);

  if(0 == 1) {
    PDM_log_trace_array_long(recv_part1_to_part2_gnum_elt2, part1_to_part2_idx[n_elt1], "recv_part1_to_part2_gnum_elt2 ::");
  }

  PDM_g_num_t p0_expected_recv_part1_to_part2_gnum_elt2[3] = {20, 30, 30};
  PDM_g_num_t p1_expected_recv_part1_to_part2_gnum_elt2[4] = {30, 30, 40, 30};

  MPI_CHECK_EQ_C_ARRAY(0, recv_part1_to_part2_gnum_elt2, p0_expected_recv_part1_to_part2_gnum_elt2, 3);
  MPI_CHECK_EQ_C_ARRAY(1, recv_part1_to_part2_gnum_elt2, p1_expected_recv_part1_to_part2_gnum_elt2, 3);

  /*
   * ************************************************************************************************
   * 3/ Exchange gnum1 to part2->part1 with ORDER_PART1_TO_PART2 : (Smart move : On recupère le resulat de l'échange d'avant :p )
   */
  PDM_g_num_t **tmp_check_recv_part2_to_part1_gnum_elt2 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(PDM_g_num_t),
                         NULL,
       (const void **)   &recv_part1_to_part2_gnum_elt2,
                         NULL,
           (void ***)    &tmp_check_recv_part2_to_part1_gnum_elt2,
                         &request);
  PDM_part_to_part_iexch_wait(ptp, request);

  PDM_g_num_t *check_recv_part2_to_part1_gnum_elt2 = tmp_check_recv_part2_to_part1_gnum_elt2[0];
  PDM_free(tmp_check_recv_part2_to_part1_gnum_elt2);

  if(0 == 1) {
    PDM_log_trace_array_long(check_recv_part2_to_part1_gnum_elt2, gnum1_come_from_idx[n_ref_lnum2], "check_recv_part2_to_part1_gnum_elt2 ::");
  }

  /*
   *  En tout rigeur le 3 pourrait venir du rang 0 et du rang 1 mais en interne on prends le premier venu --> See pdm_part_migrate
   */
  PDM_g_num_t p0_expected_check_recv_part2_to_part1_gnum_elt2[3] = {20, 30, 30};
  PDM_g_num_t p1_expected_check_recv_part2_to_part1_gnum_elt2[3] = {30, 30, 40};

  MPI_CHECK_EQ_C_ARRAY(0, check_recv_part2_to_part1_gnum_elt2, p0_expected_check_recv_part2_to_part1_gnum_elt2, 3);
  MPI_CHECK_EQ_C_ARRAY(1, check_recv_part2_to_part1_gnum_elt2, p1_expected_check_recv_part2_to_part1_gnum_elt2, 3);

  /*
   * ************************************************************************************************
   * 4/ Exchange gnum1 to part2->part1 with ORDER_PART2 (Smart move : On recupère le resulat de l'échange d'avant :p )
   */
  PDM_g_num_t **tmp_check_recv_part1_to_part2_gnum_elt2 = NULL;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM,
                                 1,
                                 sizeof(PDM_g_num_t),
                                 NULL,
               (const void **)   &recv_part2_to_part1_gnum_elt2,
                                 NULL,
                   (void ***)    &tmp_check_recv_part1_to_part2_gnum_elt2,
                                 &request);
  PDM_part_to_part_reverse_iexch_wait(ptp, request);

  PDM_g_num_t *check_recv_part1_to_part2_gnum_elt2 = tmp_check_recv_part1_to_part2_gnum_elt2[0];
  PDM_free(tmp_check_recv_part1_to_part2_gnum_elt2);

  if(0 == 1) {
    PDM_log_trace_array_long(check_recv_part1_to_part2_gnum_elt2, part1_to_part2_idx[n_elt1], "check_recv_part1_to_part2_gnum_elt2 ::");
  }

  PDM_g_num_t p0_expected_check_recv_part1_to_part2_gnum_elt2[3] = {2, 3, 3};
  PDM_g_num_t p1_expected_check_recv_part1_to_part2_gnum_elt2[4] = {3, 3, 1, 1};

  MPI_CHECK_EQ_C_ARRAY(0, check_recv_part1_to_part2_gnum_elt2, p0_expected_check_recv_part1_to_part2_gnum_elt2, 3);
  MPI_CHECK_EQ_C_ARRAY(1, check_recv_part1_to_part2_gnum_elt2, p1_expected_check_recv_part1_to_part2_gnum_elt2, 3);

  /*
   * Free
   */
  PDM_free(recv_part1_to_part2_gnum_elt2);
  PDM_free(recv_part2_to_part1_gnum_elt2);
  PDM_free(check_recv_part2_to_part1_gnum_elt2);
  PDM_free(check_recv_part1_to_part2_gnum_elt2);

  PDM_part_to_part_free(ptp);
}

