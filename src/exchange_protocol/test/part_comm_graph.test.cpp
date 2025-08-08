#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_algorithm.h"
#include "pdm_sort.h"
#include "pdm_vtk.h"
#include <functional>


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};
  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();

  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  // PDM_log_trace_array_int(lowner_bound, 3, "lowner_bound ::");

  static int lowner_bound_expected_p0[3] = {1, 1, 1};
  static int lowner_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, 3);

  // ---------------------------------------------------------------------------
  // Exchange stride cst
  std::vector<std::vector<int>> vsend_cst_data = {{-3, -2, -1}, {10, 20, 30}};
  int *send_cst_data = vsend_cst_data[i_rank].data();

  int **tmp_recv_cst_data = NULL;
  PDM_part_comm_graph_exch(pcg,
                           sizeof(int),
                           PDM_STRIDE_CST_INTERLACED,
                           1,
                           NULL,
            (void **)      &send_cst_data,
                           NULL,
           (void ***)      &tmp_recv_cst_data);
  int *recv_cst_data = tmp_recv_cst_data[0];
  free(tmp_recv_cst_data);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_cst_data, n_entity_bound, "recv_cst_data ::");
  }

  static int recv_cst_data_expected_p0[3] = {10, 20, 30};
  static int recv_cst_data_expected_p1[3] = {-3, -2, -1};

  MPI_CHECK_EQ_C_ARRAY(0, recv_cst_data, recv_cst_data_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, recv_cst_data, recv_cst_data_expected_p1, n_entity_bound);

  free(recv_cst_data);

  // ---------------------------------------------------------------------------
  // Exch variable
  std::vector<std::vector<int>> vsend_strid = {{0, 2   , 1}, {2      , 0, 1 }};
  std::vector<std::vector<int>> vsend_data  = {{  -2, 2,-1}, {10, -10,    30}};
  int *send_strid = vsend_strid[i_rank].data();
  int *send_data  = vsend_data[i_rank].data();

  int **tmp_recv_data = NULL;
  int **tmp_recv_stri = NULL;
  PDM_part_comm_graph_exch(pcg,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           -1,
                           &send_strid,
            (void **)      &send_data,
                           &tmp_recv_stri,
           (void ***)      &tmp_recv_data);
  int *recv_stri = tmp_recv_stri[0];
  int *recv_data = tmp_recv_data[0];
  free(tmp_recv_data);
  free(tmp_recv_stri);

  int n_recv_tot = 0;
  for(int i = 0; i < n_entity_bound; ++i) {
    n_recv_tot += recv_stri[i];
  }

  if(0 == 1) {
    log_trace("n_recv_tot = %i \n", n_recv_tot);
    PDM_log_trace_array_int(send_strid, n_entity_bound, "send_strid ::");
    PDM_log_trace_array_int(recv_data, n_recv_tot, "recv_data ::");
  }

  static int recv_stri_expected_p0[3] = {2, 0, 1};
  static int recv_stri_expected_p1[3] = {0, 2, 1};

  MPI_CHECK_EQ_C_ARRAY(0, recv_stri, recv_stri_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, recv_stri, recv_stri_expected_p1, n_entity_bound);

  static int recv_data_expected_p0[3] = {10, -10, 30};
  static int recv_data_expected_p1[3] = {-2, 2, -1};

  MPI_CHECK_EQ_C_ARRAY(0, recv_data, recv_data_expected_p0, n_recv_tot);
  MPI_CHECK_EQ_C_ARRAY(1, recv_data, recv_data_expected_p1, n_recv_tot);

  free(recv_stri);
  free(recv_data);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2p - Persistent exchange", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};
  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();

  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  // PDM_log_trace_array_int(lowner_bound, 3, "lowner_bound ::");

  static int lowner_bound_expected_p0[3] = {1, 1, 1};
  static int lowner_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, 3);

  // ---------------------------------------------------------------------------
  // Exchange stride cst
  std::vector<std::vector<int>> vsend_cst_data = {{-3, -2, -1}, {10, 20, 30}};
  int *send_cst_data = vsend_cst_data[i_rank].data();

  int **tmp_recv_cst_data = NULL;
  int req0 = PDM_part_comm_graph_exch_init(pcg,
                                            PDM_MPI_COMM_KIND_P2P,
                                            sizeof(int),
                                            PDM_STRIDE_CST_INTERLACED,
                                            1,
                                            NULL,
                             (void **)      &send_cst_data,
                                            NULL,
                            (void ***)      &tmp_recv_cst_data);
  PDM_part_comm_graph_exch_start(pcg, req0);
  PDM_part_comm_graph_exch_wait(pcg, req0);

  int *recv_cst_data = tmp_recv_cst_data[0];

  if(0 == 1) {
    PDM_log_trace_array_int(recv_cst_data, n_entity_bound, "recv_cst_data ::");
  }

  static int recv_cst_data_expected_p0[3] = {10, 20, 30};
  static int recv_cst_data_expected_p1[3] = {-3, -2, -1};

  MPI_CHECK_EQ_C_ARRAY(0, recv_cst_data, recv_cst_data_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, recv_cst_data, recv_cst_data_expected_p1, n_entity_bound);

  for(int i = 0; i < static_cast<int>(vsend_cst_data[i_rank].size()); ++i) {
    vsend_cst_data[i_rank][i] += 10;
  }

  PDM_part_comm_graph_exch_start(pcg, req0);
  PDM_part_comm_graph_exch_wait(pcg, req0);

  if(0 == 1) {
    PDM_log_trace_array_int(recv_cst_data, n_entity_bound, "recv_cst_data ::");
  }

  static int recv_cst_data_expected2_p0[3] = {20, 30, 40};
  static int recv_cst_data_expected2_p1[3] = { 7,  8,  9};

  MPI_CHECK_EQ_C_ARRAY(0, recv_cst_data, recv_cst_data_expected2_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, recv_cst_data, recv_cst_data_expected2_p1, n_entity_bound);

  free(recv_cst_data);
  free(tmp_recv_cst_data);

  PDM_part_comm_graph_exch_free(pcg, req0);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2p - iexch", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};
  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();

  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  // PDM_log_trace_array_int(lowner_bound, 3, "lowner_bound ::");

  static int lowner_bound_expected_p0[3] = {1, 1, 1};
  static int lowner_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, 3);

  // ---------------------------------------------------------------------------
  // Exchange stride cst
  std::vector<std::vector<int>> vsend_cst_data = {{-3, -2, -1}, {10, 20, 30}};
  int *send_cst_data = vsend_cst_data[i_rank].data();

  static int recv_cst_data_expected_p0[3] = {10, 20, 30};
  static int recv_cst_data_expected_p1[3] = {-3, -2, -1};

  std::vector<PDM_mpi_comm_kind_t> lexch_type = {PDM_MPI_COMM_KIND_P2P,
                                                 PDM_MPI_COMM_KIND_COLLECTIVE,
                                                 PDM_MPI_COMM_KIND_WIN_RMA};
  int n_type_exch = lexch_type.size();
  int n_try = 4;

  for(int i_try = 0; i_try < n_try; ++i_try) {
    for(int i_type_exch = 0; i_type_exch < n_type_exch; ++i_type_exch) {

#ifndef HAVE_MPI_COLLECTIVE_INIT_FUNC
      if(lexch_type[i_type_exch] == PDM_MPI_COMM_KIND_COLLECTIVE) {
        continue;
      }
#endif

      int **tmp_recv_cst_data = NULL;
      int req0 = PDM_part_comm_graph_iexch(pcg,
                                           lexch_type[i_type_exch],
                                           sizeof(int),
                                           PDM_STRIDE_CST_INTERLACED,
                                           1,
                                           NULL,
                            (void **)      &send_cst_data,
                                           NULL,
                           (void ***)      &tmp_recv_cst_data);
      PDM_part_comm_graph_exch_wait(pcg, req0);

      int *recv_cst_data = tmp_recv_cst_data[0];

      if(0 == 1) {
        PDM_log_trace_array_int(recv_cst_data, n_entity_bound, "recv_cst_data ::");
      }

      MPI_CHECK_EQ_C_ARRAY(0, recv_cst_data, recv_cst_data_expected_p0, n_entity_bound);
      MPI_CHECK_EQ_C_ARRAY(1, recv_cst_data, recv_cst_data_expected_p1, n_entity_bound);

      free(recv_cst_data);
      free(tmp_recv_cst_data);
    }
  }

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2p - allreduce ", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};
  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();

  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  // ------------------ MAX / INT ------------------
  std::vector<std::vector<int>> vpdata = {{1, 1, 1, 1, 1, 1, 1, 1, 1},
                                          {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}};
  int *pdata = vpdata[i_rank].data();

  PDM_part_comm_graph_all_reduce(pcg,
                                 PDM_MPI_INT,
                                 1,
                                 PDM_MPI_MAX,
            ( unsigned char **)  &pdata);

  if(0 == 1) {
    PDM_log_trace_array_int(pdata, vn_elt[i_rank], "pdata ::");
  }

  int expexted_max_int_p0[9]  = {1, 1, 2, 1, 1, 2, 1, 1, 2};
  int expexted_max_int_p1[12] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};

  MPI_CHECK_EQ_C_ARRAY(0, pdata, expexted_max_int_p0,  9);
  MPI_CHECK_EQ_C_ARRAY(1, pdata, expexted_max_int_p1, 12);

  // ------------------ MIN / INT ------------------
  for(int i = 0; i < vn_elt[i_rank]; ++i) {
    vpdata[i_rank][i] = i_rank+1;
  }

  PDM_part_comm_graph_all_reduce(pcg,
                                 PDM_MPI_INT,
                                 1,
                                 PDM_MPI_MIN,
            ( unsigned char **)  &pdata);

  if(0 == 1) {
    PDM_log_trace_array_int(pdata, vn_elt[i_rank], "pdata ::");
  }

  int expexted_min_int_p0[9]  = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  int expexted_min_int_p1[12] = {1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2};

  MPI_CHECK_EQ_C_ARRAY(0, pdata, expexted_min_int_p0,  9);
  MPI_CHECK_EQ_C_ARRAY(1, pdata, expexted_min_int_p1, 12);

  // ------------------ SUM / INT ------------------
  for(int i = 0; i < vn_elt[i_rank]; ++i) {
    vpdata[i_rank][i] = i_rank+1;
  }

  PDM_part_comm_graph_all_reduce(pcg,
                                 PDM_MPI_INT,
                                 1,
                                 PDM_MPI_SUM,
            ( unsigned char **)  &pdata);

  if(0 == 1) {
    PDM_log_trace_array_int(pdata, vn_elt[i_rank], "pdata ::");
  }

  int expexted_sum_int_p0[9]  = {1, 1, 3, 1, 1, 3, 1, 1, 3};
  int expexted_sum_int_p1[12] = {3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2};

  MPI_CHECK_EQ_C_ARRAY(0, pdata, expexted_sum_int_p0,  9);
  MPI_CHECK_EQ_C_ARRAY(1, pdata, expexted_sum_int_p1, 12);

  // ------------------ SUM / INT STRIDED ----------
  std::vector<std::vector<int>> vpdata_strided = {{1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2},
                                                  {2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3}};
  int *pdata_strided = vpdata_strided[i_rank].data();

  for(int i = 0; i < 2*vn_elt[i_rank]; ++i) {
    vpdata_strided[i_rank][i] = i_rank+1;
  }

  PDM_part_comm_graph_all_reduce(pcg,
                                 PDM_MPI_INT,
                                 2,
                                 PDM_MPI_SUM,
             ( unsigned char **) &pdata_strided);

  if(0 == 1) {
    PDM_log_trace_array_int(pdata_strided, 2*vn_elt[i_rank], "pdata_strided ::");
  }

  int expexted_sum_int_strided_p0[18] = {1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 3, 3,};
  int expexted_sum_int_strided_p1[24] = {3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2};

  MPI_CHECK_EQ_C_ARRAY(0, pdata_strided, expexted_sum_int_strided_p0, 18);
  MPI_CHECK_EQ_C_ARRAY(1, pdata_strided, expexted_sum_int_strided_p1, 24);

  PDM_part_comm_graph_free(pcg);

}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 3p", 3) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);


  /*
   *  RANK2
   *    |++++|++++|        |++++|++++|++++| 18
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3        3    4    5    6
   *
   *  RANK0                        RANK1
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */


  /* Part */
  // std::vector<int> vn_elt = {9, 12, 18};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  std::vector<int> vn_entity_bound = {6, 7, 7};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9,
                                                  7, 2, 1, 1,
                                                  8, 2, 1, 2,
                                                  9, 2, 1, 3},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9,
                                                  9, 2, 1, 3,
                                                  10, 2, 1, 4,
                                                  11, 2, 1, 5,
                                                  12, 2, 1, 6},
                                                 {1, 0, 1, 7,
                                                  2, 0, 1, 8,
                                                  3, 0, 1, 9,
                                                  4, 1, 1, 10,
                                                  5, 1, 1, 11,
                                                  3, 1, 1, 9,
                                                  6, 1, 1, 12}};


  int  n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound   = ventity_bound  [i_rank].data();

  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  // PDM_log_trace_array_int(lowner_bound, n_entity_bound, "lowner_bound ::");

  static int lowner_bound_expected_p0[6] = {1, 1, 1, 1, 1, 1};
  static int lowner_bound_expected_p1[7] = {0, 0, 0, 0, 1, 1, 1 };
  static int lowner_bound_expected_p2[7] = {0, 0, 0, 0, 0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(2, lowner_bound, lowner_bound_expected_p2, n_entity_bound);


  PDM_part_comm_graph_free(pcg);

}

MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2p - order ", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  // int n_elt1 = vn_elt[i_rank];
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};

  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();


  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  // PDM_log_trace_array_int(lowner_bound, 3, "lowner_bound ::");

  static int lowner_bound_expected_p0[3] = {1, 1, 1};
  static int lowner_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, 3);


  // std::vector<std::vector<int>> vold_to_new = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
  //                                              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
  std::vector<std::vector<int>> vold_to_new = {{8, 7, 6, 5, 4, 3, 2, 1, 0},
                                               {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};

  /* sur p0  :
    entity[0]= 1 -> entity[vold_to_new[0]] = 9
    entity[1]= 2 -> entity[vold_to_new[1]] = 8
                  ...
    entity[8]=1->entity[vold_to_new[8]] = 1

    meme chose sur p1

  */
  int *old_to_new = vold_to_new[i_rank].data();

  // PDM_log_trace_array_int(entity_bound, 4 * n_entity_bound, "entity_bound (Avant) ::");

  PDM_part_comm_graph_reorder(pcg,
                              &old_to_new);

  // PDM_log_trace_array_int(entity_bound, 4 * n_entity_bound, "entity_bound ::");

  static int entity_bound_reorder_p0[12] = {7 , 1, 1, 12,   4, 1, 1, 8,   1, 1, 1, 4};
  static int entity_bound_reorder_p1[12] = {12, 0, 1,  7,   8, 0, 1, 4,   4, 0, 1, 1};

  MPI_CHECK_EQ_C_ARRAY(0, entity_bound, entity_bound_reorder_p0, 12);
  MPI_CHECK_EQ_C_ARRAY(1, entity_bound, entity_bound_reorder_p1, 12);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p", 2) {

  // Corresponds to a QUAD of n_vtx_seg = 3
  /*

       p0                 p1
           6                 6
     3 |+++++++| 6     3 |+++++++| 6
       |       |         |       |
   3   |       |  7  2   |       |   7
       |   4   |         |   4   |
     2 |+++++++| 5     2 |+++++++| 5
       |       |         |       |
   1   |       |  5  1   |       |   5
       |       |         |       |
     1 |+++++++| 4     1 |+++++++| 4
           2                 3

  */

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<int> vn_entity1      = {6, 6};
  std::vector<int> vn_entity2      = {7, 7};
  std::vector<std::vector<int>> ventity_bound = {{4, 1, 1, 1,
                                                  5, 1, 1, 2,
                                                  6, 1, 1, 3 },
                                                 {1, 0, 1, 4,
                                                  2, 0, 1, 5,
                                                  3, 0, 1, 6 }};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 2, 4, 6, 8, 10, 12, 14},
                                                        {0, 2, 4, 6, 8, 10, 12, 14}};

  // entity2s are the faces of each mesh (edges) they are struct made of entity1s

                                                  /*  1       2       3       4       5       6       7    */
  std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   5, 4,   3, 6,   6, 5},
                                                    {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 2; // nombre de faces de bords attendu

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[8] = {5, 1, 1, 1,   7, 1, 1, 2};
  static int entity_bound_reorder_p1[8] = {1, 0, 1, 5,   2, 0, 1, 7};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 8);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 8);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}



MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - invert sense", 2) {

  // Correspond to a QUAD of n_vtx_seg = 3

/*
       p0                 p1
           6                 6
     3 |+++++++| 6     3 |+++++++| 6
       |       |         |       |
   3   |       |  7  2   |       |   7
       |   4   |         |   4   |
     2 |+++++++| 5     2 |+++++++| 5
       |       |         |       |
   1   |       |  5  1   |       |   5
       |       |         |       |
     1 |+++++++| 4     1 |+++++++| 4
           2                 3

*/
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<int> vn_entity1      = {6, 6};
  std::vector<int> vn_entity2      = {7, 7};
  std::vector<std::vector<int>> ventity_bound = {{4, 1, 1, 1,
                                                  5, 1, 1, 2,
                                                  6, 1, 1, 3 },
                                                 {1, 0, 1, 4,
                                                  2, 0, 1, 5,
                                                  3, 0, 1, 6 }};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 2, 4, 6, 8, 10, 12, 14},
                                                        {0, 2, 4, 6, 8, 10, 12, 14}};

  // std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   5, 4,   3, 6,   6, 5},
  //                                                   {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

                                                                                 /*inverted        inverted   */
  //                                                                                |----|          |----|
  std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   4, 5,   3, 6,   5, 6},
                                                    {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 2;

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[8] = {5, 1, 1, -1,   7, 1, 1, -2};
  static int entity_bound_reorder_p1[8] = {1, 0, 1, -5,   2, 0, 1, -7};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 8);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 8);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}




MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - 3D ", 2) {

  // Corresponds to a HEXA of n_vtx_seg = 3
  // here we only represent the boundary between the 2 parts
  //

  //              --------- +18              9+--------
  //                       /|                /|
  //                      / |               / |
  //                  15 /  |              /  |
  //              ----- +   |            6+---|---
  //                   /|   |            /|   |
  //                  / |17 |           / | 4 |
  //              12 /- |-- +17        /  |  8+-------
  //             -- +   |  /|        3+-- |--/|--
  //                |   | / |         |   | / |
  //                |15 |/  |         | 2 |/  |
  //             -- |-- +14 |         |  5+-- |----
  //                |  /|   |         |  /|   |
  //                | / |16 |         | / | 3 |
  //             11 |/ -| --+16       |/  |  7+-------
  //            --- +   |  /         2+-- |- /----
  //                |   | /           | 1 | /
  //                |14 |/            |   |/
  //   x          --|-- +             |  4+------
  //   ^  y         |  /13            |  /
  //   | +          | /               | /
  //   |/           |/                |/
  //   +--->z   --- +                1+------
  //              10

  //                        z=0.5 plane
  // all normals of boundary faces are z-positive


  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // Keep for debug
  std::vector<std::vector<double>> vvtx_coords = {{0.0, 0.0, 0.0,    /*1*/
                                                   0.5, 0.0, 0.0,    /*2*/
                                                   1.0, 0.0, 0.0,    /*3*/
                                                   0.0, 0.5, 0.0,    /*4*/
                                                   0.5, 0.5, 0.0,    /*5*/
                                                   1.0, 0.5, 0.0,    /*6*/
                                                   0.0, 1.0, 0.0,    /*7*/
                                                   0.5, 1.0, 0.0,    /*8*/
                                                   1.0, 1.0, 0.0,    /*9*/
                                                   0.0, 0.0, 0.5,    /*10*/
                                                   0.5, 0.0, 0.5,    /*11*/
                                                   1.0, 0.0, 0.5,    /*12*/
                                                   0.0, 0.5, 0.5,    /*13*/
                                                   0.5, 0.5, 0.5,    /*14*/
                                                   1.0, 0.5, 0.5,    /*15*/
                                                   0.0, 1.0, 0.5,    /*16*/
                                                   0.5, 1.0, 0.5,    /*17*/
                                                   1.0, 1.0, 0.5 },  /*18*/

                                                  {0.0, 0.0, 0.5,    /*1*/
                                                   0.5, 0.0, 0.5,    /*2*/
                                                   1.0, 0.0, 0.5,    /*3*/
                                                   0.0, 0.5, 0.5,    /*4*/
                                                   0.5, 0.5, 0.5,    /*5*/
                                                   1.0, 0.5, 0.5,    /*6*/
                                                   0.0, 1.0, 0.5,    /*7*/
                                                   0.5, 1.0, 0.5,    /*8*/
                                                   1.0, 1.0, 0.5,    /*9*/
                                                   0.0, 0.0, 1.0,    /*10*/
                                                   0.5, 0.0, 1.0,    /*11*/
                                                   1.0, 0.0, 1.0,    /*12*/
                                                   0.0, 0.5, 1.0,    /*13*/
                                                   0.5, 0.5, 1.0,    /*14*/
                                                   1.0, 0.5, 1.0,    /*15*/
                                                   0.0, 1.0, 1.0,    /*16*/
                                                   0.5, 1.0, 1.0,    /*17*/
                                                   1.0, 1.0, 1.0}};  /*18*/

  std::vector<int> vn_entity_bound = {9 ,  9};
  std::vector<int> vn_entity1      = {18, 18};
  std::vector<int> vn_entity2      = {20, 20};
  std::vector<std::vector<int>> ventity_bound = {{10, 1, 1, 1,
                                                  11, 1, 1, 2,
                                                  12, 1, 1, 3,
                                                  13, 1, 1, 4,
                                                  14, 1, 1, 5,
                                                  15, 1, 1, 6,
                                                  16, 1, 1, 7,
                                                  17, 1, 1, 8,
                                                  18, 1, 1, 9 },
                                                 {1, 0, 1, 10,
                                                  2, 0, 1, 11,
                                                  3, 0, 1, 12,
                                                  4, 0, 1, 13,
                                                  5, 0, 1, 14,
                                                  6, 0, 1, 15,
                                                  7, 0, 1, 16,
                                                  8, 0, 1, 17,
                                                  9, 0, 1, 18}};


  /* ici les entity2 sont des quads -> donc composés de 4 noeuds à chaque fois*/
  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80},
                                                        {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80}};

  std::vector<std::vector<int>> ventity2_entity1 = {{2, 1, 4, 5,      /*1 */
                                                     2, 5, 6, 3,      /*2 */
                                                     1, 2, 11, 10,    /*3 */
                                                     7, 8, 5, 4,      /*4 */
                                                     1, 10, 13, 4,    /*5 */
                                                     11, 2, 3, 12,    /*6 */
                                                     8, 9, 6, 5,      /*7 */
                                                     11, 2, 5, 14,    /*8 */
                                                     14, 5, 4, 13,    /*9 */
                                                     12, 3, 6, 15,    /*10 */
                                                     15, 6, 5, 14,    /*11 */
                                                     4, 13, 16, 7,    /*12 */
                                                     14, 5, 8, 17,    /*13 */
                                                     14, 13, 10, 11,  /*14  bound*/
                                                     17, 8, 7, 16,    /*15 */
                                                     15, 6, 9, 18,    /*16 */
                                                     15, 14, 11, 12,  /*17  bound*/
                                                     18, 9, 8, 17,    /*18 */
                                                     17, 16, 13, 14,  /*19  bound*/
                                                     18, 17, 14, 15}, /*20  bound*/

                                                    {5, 4, 1, 2,      /*1   bound*/
                                                     6, 5, 2, 3,      /*2   bound*/
                                                     8, 7, 4, 5,      /*3   bound*/
                                                     11, 10, 1, 2,    /*4*/
                                                     9, 8, 5, 6,      /*5   bound*/
                                                     4, 1, 10, 13,    /*6*/
                                                     12, 11, 2, 3,    /*7*/
                                                     11, 2, 5, 14,    /*8*/
                                                     14, 5, 4, 13,    /*9*/
                                                     15, 12, 3, 6,    /*10*/
                                                     15, 6, 5, 14,    /*11*/
                                                     16, 7, 4, 13,    /*12*/
                                                     14, 5, 8, 17,    /*13*/
                                                     14, 13, 10, 11,  /*14*/
                                                     17, 8, 7, 16,    /*15*/
                                                     18, 15, 6, 9,    /*16*/
                                                     15, 14, 11, 12,  /*17*/
                                                     18, 9, 8, 17,    /*18*/
                                                     17, 16, 13, 14,  /*19*/
                                                     18, 17, 14, 15}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 4; // nombre de faces de bords attendu

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[16] = {14, 1, 1,  1,   17, 1, 1,  2,   19, 1, 1,  3,   20, 1, 1, 5};
  static int entity_bound_reorder_p1[16] = { 1, 0, 1, 14,    2, 0, 1, 17,    3, 0, 1, 19,    5, 0, 1,20};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 16);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");

      int *face_tag = (int *) malloc(pn_entity2 * sizeof(int));

      for(int i = 0; i < pn_entity2; ++i) {
        face_tag[i] = -1;
      }

      for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {
        int i_face = pentity2_graph[i_part][4*idx]-1;
        int t_rank = pentity2_graph[i_part][4*idx+1];
        face_tag[i_face] = t_rank;
      }

      const char* field_name[] = {"face_tag", 0 };
      const int*  field     [] = {face_tag};

      char filename[999];
      sprintf(filename, "out_face_graph_i_part=%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 pn_entity1,
                                 vvtx_coords[i_rank].data(),
                                 NULL,
                                 PDM_MESH_NODAL_QUAD4,
                                 pn_entity2,
                                 entity2_entity1,
                                 NULL,
                                 1,
                                 field_name,
                                 (const int **)   &field);


      free(face_tag);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}



MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - 3D - invert face_vtx", 2) {

  // Correspond to a HEXA of n_vtx_seg = 3

  // here we only represent the boundary between the 2 parts
  //

  //              --------- +18              9+--------
  //                       /|                /|
  //                      / |               / |
  //                  15 /  |              /  |
  //              ----- +   |            6+---|---
  //                   /|   |            /|   |
  //                  / |17 |           / | 4 |
  //              12 /- |-- +17        /  |  8+-------
  //             -- +   |  /|        3+-- |--/|--
  //                |   | / |         |   | / |
  //                |15 |/  |         | 2 |/  |
  //             -- |-- +14 |         |  5+-- |----
  //                |  /|   |         |  /|   |
  //                | / |16 |         | / | 3 |
  //             11 |/ -| --+16       |/  |  7+-------
  //            --- +   |  /         2+-- |- /----
  //                |   | /           | 1 | /
  //                |14 |/            |   |/
  //   x          --|-- +             |  4+------
  //   ^  y         |  /13            |  /
  //   | +          | /               | /
  //   |/           |/                |/
  //   +--->z   --- +                1+------
  //              10

  //                        z=0.5 plane
  // p1 : all normals of boundary faces are z-negative
  // p2 :  "    "     "     "      "     "  z-positive


  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // Keep for debug
  std::vector<std::vector<double>> vvtx_coords = {{0.0, 0.0, 0.0,
                                                   0.5, 0.0, 0.0,
                                                   1.0, 0.0, 0.0,
                                                   0.0, 0.5, 0.0,
                                                   0.5, 0.5, 0.0,
                                                   1.0, 0.5, 0.0,
                                                   0.0, 1.0, 0.0,
                                                   0.5, 1.0, 0.0,
                                                   1.0, 1.0, 0.0,
                                                   0.0, 0.0, 0.5,
                                                   0.5, 0.0, 0.5,
                                                   1.0, 0.0, 0.5,
                                                   0.0, 0.5, 0.5,
                                                   0.5, 0.5, 0.5,
                                                   1.0, 0.5, 0.5,
                                                   0.0, 1.0, 0.5,
                                                   0.5, 1.0, 0.5,
                                                   1.0, 1.0, 0.5 },
                                                  {0.0, 0.0, 0.5,
                                                   0.5, 0.0, 0.5,
                                                   1.0, 0.0, 0.5,
                                                   0.0, 0.5, 0.5,
                                                   0.5, 0.5, 0.5,
                                                   1.0, 0.5, 0.5,
                                                   0.0, 1.0, 0.5,
                                                   0.5, 1.0, 0.5,
                                                   1.0, 1.0, 0.5,
                                                   0.0, 0.0, 1.0,
                                                   0.5, 0.0, 1.0,
                                                   1.0, 0.0, 1.0,
                                                   0.0, 0.5, 1.0,
                                                   0.5, 0.5, 1.0,
                                                   1.0, 0.5, 1.0,
                                                   0.0, 1.0, 1.0,
                                                   0.5, 1.0, 1.0,
                                                   1.0, 1.0, 1.0}};

  std::vector<int> vn_entity_bound = {9 ,  9};
  std::vector<int> vn_entity1      = {18, 18};
  std::vector<int> vn_entity2      = {20, 20};
  std::vector<std::vector<int>> ventity_bound = {{10, 1, 1, 1,
                                                  11, 1, 1, 2,
                                                  12, 1, 1, 3,
                                                  13, 1, 1, 4,
                                                  14, 1, 1, 5,
                                                  15, 1, 1, 6,
                                                  16, 1, 1, 7,
                                                  17, 1, 1, 8,
                                                  18, 1, 1, 9 },
                                                 {1, 0, 1, 10,
                                                  2, 0, 1, 11,
                                                  3, 0, 1, 12,
                                                  4, 0, 1, 13,
                                                  5, 0, 1, 14,
                                                  6, 0, 1, 15,
                                                  7, 0, 1, 16,
                                                  8, 0, 1, 17,
                                                  9, 0, 1, 18}};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80},
                                                        {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80}};

  std::vector<std::vector<int>> ventity2_entity1 = {{2, 1, 4, 5,      // 1
                                                     2, 5, 6, 3,      // 2
                                                     1, 2, 11, 10,    // 3
                                                     7, 8, 5, 4,      // 4
                                                     1, 10, 13, 4,    // 5
                                                     11, 2, 3, 12,    // 6
                                                     8, 9, 6, 5,      // 7
                                                     11, 2, 5, 14,    // 8
                                                     14, 5, 4, 13,    // 9
                                                     12, 3, 6, 15,    // 10
                                                     15, 6, 5, 14,    // 11
                                                     4, 13, 16, 7,    // 12
                                                     14, 5, 8, 17,    // 13
                                                     11, 10, 13, 14,  // 14 original : 14, 13, 10, 11,
                                                     17, 8, 7, 16,    // 15
                                                     15, 6, 9, 18,    // 16
                                                     12, 11, 14, 15,  // 17 original : 15, 14, 11, 12
                                                     18, 9, 8, 17,    // 18
                                                     14, 13, 16, 17,  // 19 origianl : 17, 16, 13, 14
                                                     15, 14, 17, 18}, // 20 original : 18, 17, 14, 15
                                                    {5, 4, 1, 2,
                                                     6, 5, 2, 3,
                                                     8, 7, 4, 5,
                                                     11, 10, 1, 2,
                                                     9, 8, 5, 6,
                                                     4, 1, 10, 13,
                                                     12, 11, 2, 3,
                                                     11, 2, 5, 14,
                                                     14, 5, 4, 13,
                                                     15, 12, 3, 6,
                                                     15, 6, 5, 14,
                                                     16, 7, 4, 13,
                                                     14, 5, 8, 17,
                                                     14, 13, 10, 11,
                                                     17, 8, 7, 16,
                                                     18, 15, 6, 9,
                                                     15, 14, 11, 12,
                                                     18, 9, 8, 17,
                                                     17, 16, 13, 14,
                                                     18, 17, 14, 15}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 4;

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[16] = {14, 1, 1,  -1, 17, 1, 1,  -2, 19, 1, 1, -3, 20, 1, 1, -5};
  static int entity_bound_reorder_p1[16] = { 1, 0, 1, -14,  2, 0, 1, -17,  3, 0, 1,-19,  5, 0, 1,-20};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 16);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");

      int *face_tag = (int *) malloc(pn_entity2 * sizeof(int));

      for(int i = 0; i < pn_entity2; ++i) {
        face_tag[i] = -1;
      }

      for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {
        int i_face = pentity2_graph[i_part][4*idx]-1;
        int t_rank = pentity2_graph[i_part][4*idx+1];
        face_tag[i_face] = t_rank;
      }

      const char* field_name[] = {"face_tag", 0 };
      const int*  field     [] = {face_tag};

      char filename[999];
      sprintf(filename, "out_face_graph_i_part=%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 pn_entity1,
                                 vvtx_coords[i_rank].data(),
                                 NULL,
                                 PDM_MESH_NODAL_QUAD4,
                                 pn_entity2,
                                 entity2_entity1,
                                 NULL,
                                 1,
                                 field_name,
                                 (const int **)   &field);


      free(face_tag);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 1 perio - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *              9 +---+---+---+ 12
   *                |           |
   *              5 +   rank 1  + 8
   *                |           |
   *                +---+---+---+
   *                1   2   3   4
   * interface -1                   interface +1
   *                9  10  11  12
   *                +---+---+---+
   *                |           |
   *              5 +   rank 0  + 8
   *                |           |
   *              1 +---+---+---+ 4
   *
   * --- rank 0 ---
   *  1 -> (0, 1,  4) through interface -1
   *
   *  4 -> (0, 1,  1) through interface  1
   *
   *  5 -> (0, 1,  8) through interface -1
   *
   *  8 -> (0, 1,  5) through interface  1
   *
   *  9 -> (1, 1,  1) through interface 0
   *    -> (0, 1, 12) through interface -1
   *    -> (1, 1,  4) through interface -1
   *
   * 10 -> (1, 1,  2) through interface 0
   *
   * 11 -> (1, 1,  3) through interface 0
   *
   * 12 -> (1, 1,  4) through interface 0
   * 12 -> (0, 1,  9) through interface 1
   * 12 -> (1, 1,  1) through interface 1
   *
   * owners : 1, 5, 9, 10, 11, 12
   *
   *
   * --- rank 1 ---
   *  1 -> (0, 1,  9) through interface 0
   *    -> (1, 1,  4) through interface -1
   *    -> (0, 1, 12) through interface -1
   *
   *  2 -> (0, 1, 10) through interface 0
   *
   *  3 -> (0, 1, 11) through interface 0
   *
   *  4 -> (0, 1, 12) through interface 0
   *    -> (1, 1,  1) through interface 1
   *    -> (0, 1,  9) through interface 1
   *
   *  5 -> (1, 1,  8) through interface -1
   *
   *  8 -> (1, 1,  5) through interface  1
   *
   *  9 -> (1, 1, 12) through interface -1
   *
   * 12 -> (1, 1,  9) through interface  1
   *
   * owners : 5, 9
   *
   */

  /* Part */
  int n_part = 1;

  /* Comm graph */
  std::vector<int> vn_entity_bound = {12, 12};
  std::vector<std::vector<int>> ventity_bound = {{1,  0, 1,  4,
                                                  4,  0, 1,  1,
                                                  5,  0, 1,  8,
                                                  8,  0, 1,  5,
                                                  9,  1, 1,  1,
                                                  9,  0, 1, 12,
                                                  9,  1, 1,  4,
                                                  10, 1, 1,  2,
                                                  11, 1, 1,  3,
                                                  12, 1, 1,  4,
                                                  12, 0, 1,  9,
                                                  12, 1, 1,  1},
                                                 {1,  0, 1,  9,
                                                  1,  1, 1,  4,
                                                  1,  0, 1, 12,
                                                  2,  0, 1, 10,
                                                  3,  0, 1, 11,
                                                  4,  0, 1, 12,
                                                  4,  1, 1,  1,
                                                  4,  0, 1,  9,
                                                  5,  1, 1,  8,
                                                  8,  1, 1,  5,
                                                  9,  1, 1, 12,
                                                  12, 1, 1,  9}};
  std::vector<std::vector<int>> ventity_interface = {{-1,  1, -1, 1, 0, -1, -1, 0,  0, 0,  1, 1},
                                                     { 0, -1, -1, 0, 0,  0,  1, 1, -1, 1, -1, 1}};

  int n_entity_bound = vn_entity_bound  [i_rank];
  int *entity_bound  = ventity_bound    [i_rank].data();
  int *entity_nuplet = ventity_interface[i_rank].data();

  PDM_part_comm_graph_t *pcg = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                                      &n_entity_bound,
                                                                      &entity_bound,
                                                                      PDM_OWNERSHIP_USER,
                                                                      1,
                                                                      &entity_nuplet,
                                                                      PDM_OWNERSHIP_USER,
                                                                      PDM_TRUE,
                                                                      pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);
  // PDM_log_trace_array_int(lowner_bound, n_entity_bound, "lowner_bound ::");

  static int lowner_bound_expected_p0[12] = {1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0};
  static int lowner_bound_expected_p1[12] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, n_entity_bound);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 2 perio - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *                 interface +2
   *
   *                6  4----5----6
   *                |\  \        |
   *                | \  \  rank |
   *                |  \  \   1  |
   *  interface -1  |   \  \     |  interface +1
   *                4    5  2    3
   *                |     \  \   |
   *                |      \  \  |
   *                | rank  \  \ |
   *                |  0     \  \|
   *                1----2----3  1
   *
   *                 interface -2
   *
   * --- rank 0 ---
   * 1 -> (1, 1, 1) through interface -1
   *   -> (0, 1, 3) through interface -1
   *   -> (0, 1, 6) through interface -2
   *   -> (1, 1, 4) through interface -2
   *   -> (1, 1, 6) through interface  3
   *
   * 2 -> (1, 1, 5) through interface -2
   *
   * 3 -> (1, 1, 1) through interface  0
   *   -> (0, 1, 1) through interface +1
   *   -> (1, 1, 4) through interface -4
   *   -> (0, 1, 6) through interface -4
   *   -> (1, 1, 6) through interface -2
   *
   * 4 -> (1, 1, 3) through interface -1
   *
   * 5 -> (1, 1, 2) through interface  0
   *
   * 6 -> (1, 1, 4) through interface  0
   *   -> (1, 1, 6) through interface -1
   *   -> (0, 1, 1) through interface +2
   *   -> (0, 1, 3) through interface +4
   *   -> (1, 1, 1) through interface +4
   *
   * owners : 1, 2, 4, 5
   *
   * --- rank 1 ---
   * 1 -> (0, 1, 1) through interface +1
   *   -> (0, 1, 3) through interface  0
   *   -> (0, 1, 6) through interface -4
   *   -> (1, 1, 4) through interface -4
   *   -> (1, 1, 6) through interface -2
   *
   * 2 -> (0, 1, 5) through interface  0
   *
   * 3 -> (0, 1, 4) through interface +1
   *
   * 4 -> (0, 1, 6) through interface  0
   *   -> (1, 1, 6) through interface -1
   *   -> (0, 1, 1) through interface +2
   *   -> (1, 1, 1) through interface +4
   *   -> (0, 1, 3) through interface +4
   *
   * 5 -> (0, 1, 2) through interface +2
   *
   * 6 -> (0, 1, 1) through interface -3
   *   -> (1, 1, 1) through interface +2
   *   -> (0, 1, 3) through interface +2
   *   -> (1, 1, 4) through interface +1
   *   -> (0, 1, 6) through interface +1
   *
   * owners : ∅
   *
   */

  /* Part */
  int n_part = 1;

  /* Comm graph */
  std::vector<int> vn_entity_bound = {18, 18};
  std::vector<std::vector<int>> ventity_bound = {{1, 1, 1, 1,
                                                  1, 0, 1, 3,
                                                  1, 0, 1, 6,
                                                  1, 1, 1, 4,
                                                  1, 1, 1, 6,
                                                  2, 1, 1, 5,
                                                  3, 1, 1, 1,
                                                  3, 0, 1, 1,
                                                  3, 1, 1, 4,
                                                  3, 0, 1, 6,
                                                  3, 1, 1, 6,
                                                  4, 1, 1, 3,
                                                  5, 1, 1, 2,
                                                  6, 1, 1, 4,
                                                  6, 1, 1, 6,
                                                  6, 0, 1, 1,
                                                  6, 0, 1, 3,
                                                  6, 1, 1, 1},
                                                 {1, 0, 1, 1,
                                                  1, 0, 1, 3,
                                                  1, 0, 1, 6,
                                                  1, 1, 1, 4,
                                                  1, 1, 1, 6,
                                                  2, 0, 1, 5,
                                                  3, 0, 1, 4,
                                                  4, 0, 1, 6,
                                                  4, 1, 1, 6,
                                                  4, 0, 1, 1,
                                                  4, 1, 1, 1,
                                                  4, 0, 1, 3,
                                                  5, 0, 1, 2,
                                                  6, 0, 1, 1,
                                                  6, 1, 1, 1,
                                                  6, 0, 1, 3,
                                                  6, 1, 1, 4,
                                                  6, 0, 1, 6}};
  std::vector<std::vector<int>> ventity_interface = {{-1, -1, -2, -2,  3, -2, 0, 1, -4, -4, -2, -1, 0,  0, -1, 2, 4, 4},
                                                     { 1,  0, -4, -4, -2,  0, 1, 0, -1,  2,  4,  4, 2, -3,  2, 2, 1, 1}};

  int n_entity_bound = vn_entity_bound  [i_rank];
  int *entity_bound  = ventity_bound    [i_rank].data();
  int *entity_nuplet = ventity_interface[i_rank].data();

  PDM_part_comm_graph_t *pcg = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                                      &n_entity_bound,
                                                                      &entity_bound,
                                                                      PDM_OWNERSHIP_USER,
                                                                      1,
                                                                      &entity_nuplet,
                                                                      PDM_OWNERSHIP_USER,
                                                                      PDM_TRUE,
                                                                      pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);
  // PDM_log_trace_array_int(lowner_bound, n_entity_bound, "lowner_bound ::");

  static int lowner_bound_expected_p0[18] = {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
  static int lowner_bound_expected_p1[18] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, n_entity_bound);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - 1 part - 1 perio axi - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *                             2 ⎽----⎽ 1
   *                           ⟋ rank ⟋ /
   *                         ⟋   1  ⟋  /
   *                       5 ⎽----⎽ 4  /
   *                  5 ⎽----⎽ 6  /  ⟋ 3
   *                ⟋ rank ⟋ /  / ⟋
   *              ⟋   0  ⟋  /\ /⟋
   *            3 ⎽----⎽ 4  /  6
   *              \    /  ⟋ 2
   * interface -1  \  / ⟋  interface +1
   *                \/⟋
   *                1
   *
   *
   * --- rank 0 ---
   * 1 -> (0, 1, 1) through interface -1
   *   -> (0, 1, 1) through interface +1
   *
   * 2 -> (1, 1, 6) through interface  0
   *   -> (0, 1, 2) through interface -1
   *   -> (0, 1, 2) through interface +1
   *
   * 3 -> (0, 1, 4) through interface -1
   *
   * 4 -> (0, 1, 3) through interface +1
   *
   * 5 -> (0, 1, 6) through interface -1
   *   -> (1, 1, 5) through interface  0
   *   -> (1, 1, 4) through interface -1
   *
   * 6 -> (0, 1, 5) through interface +1
   *   -> (1, 1, 4) through interface  0
   *   -> (1, 1, 5) through interface +1
   *
   * owners : 1, 2, 3, 5
   *
   *
   * --- rank 1 ---
   * 1 -> (1, 1, 2) through interface +1
   *
   * 2 -> (1, 1, 1) through interface -1
   *
   * 3 -> (1, 1, 3) through interface -1
   *   -> (1, 1, 3) through interface +1
   *
   * 4 -> (1, 1, 5) through interface +1
   *   -> (0, 1, 6) through interface  0
   *   -> (0, 1, 5) through interface +1
   *
   * 5 -> (1, 1, 4) through interface -1
   *   -> (0, 1, 5) through interface  0
   *   -> (0, 1, 6) through interface -1
   *
   * 6 -> (0, 1, 2) through interface  0
   *   -> (1, 1, 6) through interface -1
   *   -> (1, 1, 6) through interface +1
   *
   * owners : 1, 3
   *
   */

  /* Part */
  int n_part = 1;

  /* Comm graph */
  std::vector<int> vn_entity_bound = {13, 13};
  std::vector<std::vector<int>> ventity_bound = {{1, 0, 1, 1,
                                                  1, 0, 1, 1,
                                                  2, 1, 1, 6,
                                                  2, 0, 1, 2,
                                                  2, 0, 1, 2,
                                                  3, 0, 1, 4,
                                                  4, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  5, 1, 1, 5,
                                                  5, 1, 1, 4,
                                                  6, 0, 1, 5,
                                                  6, 1, 1, 4,
                                                  6, 1, 1, 5},
                                                 {1, 1, 1, 2,
                                                  2, 1, 1, 1,
                                                  3, 1, 1, 3,
                                                  3, 1, 1, 3,
                                                  4, 1, 1, 5,
                                                  4, 0, 1, 6,
                                                  4, 0, 1, 5,
                                                  5, 1, 1, 4,
                                                  5, 0, 1, 5,
                                                  5, 0, 1, 6,
                                                  6, 0, 1, 2,
                                                  6, 1, 1, 6,
                                                  6, 1, 1, 6}};

  std::vector<std::vector<int>> ventity_interface = {{-1,  1,  0, -1,  1, -1,  1, -1,  0, -1,  1,  0,  1},
                                                     { 1, -1, -1,  1,  1,  0,  1, -1,  0, -1,  0, -1,  1}};

  int n_entity_bound = vn_entity_bound  [i_rank];
  int *entity_bound  = ventity_bound    [i_rank].data();
  int *entity_nuplet = ventity_interface[i_rank].data();

  PDM_part_comm_graph_t *pcg = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                                      &n_entity_bound,
                                                                      &entity_bound,
                                                                      PDM_OWNERSHIP_USER,
                                                                      1,
                                                                      &entity_nuplet,
                                                                      PDM_OWNERSHIP_USER,
                                                                      PDM_TRUE,
                                                                      pdm_comm);

  const int* lowner_bound = PDM_part_comm_graph_owner_get(pcg, 0);

  if(1 == 0) {
    PDM_log_trace_array_int(lowner_bound, n_entity_bound, "lowner_bound ::");
  }

  static int lowner_bound_expected_p0[13] = {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0};
  static int lowner_bound_expected_p1[13] = {1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lowner_bound, lowner_bound_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, lowner_bound, lowner_bound_expected_p1, n_entity_bound);

  PDM_part_comm_graph_free(pcg);
}


MPI_TEST_CASE("[PDM_part_comm_graph] - gather strided data", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *    |++++|++++| 9    9 |++++|++++|++++| 12
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++| 6    5 |++++|++++|++++| 8
   *    |    |    |        |    |    |    |
   *    |    |    |        |    |    |    |
   *    |++++|++++|        |++++|++++|++++|
   *   1     2    3       1     2    3    4
   */

  /* Part */
  std::vector<int> vn_elt = {9, 12};
  int n_part = 1;

  /* Graphe comm */
  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<std::vector<int>> ventity_bound = {{3, 1, 1, 1,
                                                  6, 1, 1, 5,
                                                  9, 1, 1, 9},
                                                 {1, 0, 1, 3,
                                                  5, 0, 1, 6,
                                                  9, 0, 1, 9}};

  int n_entity_bound = vn_entity_bound[i_rank];
  int *entity_bound  = ventity_bound  [i_rank].data();


  PDM_part_comm_graph_t* pcg = PDM_part_comm_graph_create(n_part,
                                                          &n_entity_bound,
                                                          &entity_bound,
                                                          PDM_OWNERSHIP_USER,
                                                          pdm_comm);


  int n_vtx = vn_elt[i_rank];
  std::vector<std::vector<int   >> vtx_data_n   = {{1, 0, 1, 0, 0, 0, 1, 0, 1},
                                                   {1, 0, 0, 1, 2, 0, 0, 0, 1, 0, 0, 1}};
  std::vector<std::vector<double>> vtx_data     = {{0.,0.,0.,
                                                    0.,0.,2.,
                                                    2.,0.,0.,
                                                    2.,0.,2.,},
                                                   {0.,1.,2.,
                                                    0.,1.,6.,
                                                    2.,1.,1., -2.,1.,-1.,
                                                    2.,1.,2.,
                                                    6.,1.,2.}};

  int    **gather_vtx_data_n = NULL;
  double **gather_vtx_data   = NULL;
  PDM_part_comm_graph_gather_strided_data(pcg,
                                          3*sizeof(double),
                                          PDM_STRIDE_CST_INTERLACED,
                                          &n_vtx,
                               (int   **) &vtx_data_n[i_rank],
                               (void  **) &vtx_data  [i_rank],
                                          &gather_vtx_data_n,
                               (void ***) &gather_vtx_data);

  std::vector<std::vector<int   >> expctd_vtx_data_n   = {{1, 0, 2, 0, 0, 2, 1, 0, 2},
                                                          {2, 0, 0, 1, 2, 0, 0, 0, 2, 0, 0, 1}};
  std::vector<std::vector<double>> expctd_vtx_data     = {{0.,0.,0.,
                                                           0.,0.,2.,  0.,1., 2.,
                                                           2.,1.,1., -2.,1.,-1.,
                                                           2.,0.,0.,
                                                           2.,0.,2., 2.,1.,2.},
                                                          {0.,1.,2., 0.,0.,2.,
                                                           0.,1.,6.,
                                                           2.,1.,1., -2.,1.,-1.,
                                                           2.,1.,2.,  2.,0., 2.,
                                                           6.,1.,2.}};

  int i_read = 0;
  for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
    CHECK(gather_vtx_data_n[0][i_vtx] == expctd_vtx_data_n[i_rank][i_vtx]);
    for (int i_data=0; i_data<gather_vtx_data_n[0][i_vtx]; ++i_data) {
      CHECK(gather_vtx_data[0][3*i_read  ] == doctest::Approx(expctd_vtx_data[i_rank][3*i_read  ]).epsilon(0.01));
      CHECK(gather_vtx_data[0][3*i_read+1] == doctest::Approx(expctd_vtx_data[i_rank][3*i_read+1]).epsilon(0.01));
      CHECK(gather_vtx_data[0][3*i_read+2] == doctest::Approx(expctd_vtx_data[i_rank][3*i_read+2]).epsilon(0.01));
      i_read++;
    }
  }

  PDM_free(gather_vtx_data_n[0]);
  PDM_free(gather_vtx_data  [0]);
  PDM_free(gather_vtx_data_n);
  PDM_free(gather_vtx_data);

  PDM_part_comm_graph_free(pcg);
}
