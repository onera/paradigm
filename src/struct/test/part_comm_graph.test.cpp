#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_comm_graph.h"
#include "pdm_logging.h"
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

  PDM_part_comm_graph_t* ptpgc = PDM_part_comm_graph_create(n_part,
                                                                            &n_entity_bound,
                                                                            &entity_bound,
                                                                            pdm_comm);

  const int* lower_bound = PDM_part_comm_graph_owner_get(ptpgc, 0);

  // PDM_log_trace_array_int(lower_bound, 3, "lower_bound ::");

  static int lower_bound_expected_p0[3] = {1, 1, 1};
  static int lower_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lower_bound, lower_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lower_bound, lower_bound_expected_p1, 3);

  // ---------------------------------------------------------------------------
  // Exchange stride cst
  std::vector<std::vector<int>> vsend_cst_data = {{-3, -2, -1}, {10, 20, 30}};
  int *send_cst_data = vsend_cst_data[i_rank].data();

  int **tmp_recv_cst_data = NULL;
  PDM_part_comm_graph_exch(ptpgc,
                                   sizeof(int),
                                   PDM_STRIDE_CST_INTERLACED,
                                   1,
                                   NULL,
                    (void **)      &send_cst_data,
                                   NULL,
                   (void ***)      &tmp_recv_cst_data);
  int *recv_cst_data = tmp_recv_cst_data[0];
  free(tmp_recv_cst_data);

  // std::vector<int> part1_flags(n_elt1);
  // int **tmp_recv_cst_data = NULL;
  // PDM_part_comm_graph_exch(ptpgc,
  //                                  sizeof(int),
  //                                  PDM_STRIDE_CST_INTERLACED,
  //                                  PDM_PART_TO_PART_DATA_DEF_ORDER_PART,
  //                                  1,
  //                   (void **)      &part1_flags,
  //                  (void ***)      &part1_flags);
  // int *recv_cst_data = tmp_recv_cst_data[0];
  // free(tmp_recv_cst_data);

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
  PDM_part_comm_graph_exch(ptpgc,
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

  // ---------------------------------------------------------------------------
  // Exchange cst raw




  PDM_part_comm_graph_free(ptpgc);
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

  PDM_part_comm_graph_t* ptpgc = PDM_part_comm_graph_create(n_part,
                                                                            &n_entity_bound,
                                                                            &entity_bound,
                                                                            pdm_comm);

  const int* lower_bound = PDM_part_comm_graph_owner_get(ptpgc, 0);

  // PDM_log_trace_array_int(lower_bound, n_entity_bound, "lower_bound ::");

  static int lower_bound_expected_p0[6] = {1, 1, 1, 1, 1, 1};
  static int lower_bound_expected_p1[7] = {0, 0, 0, 0, 1, 1, 1 };
  static int lower_bound_expected_p2[7] = {0, 0, 0, 0, 0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lower_bound, lower_bound_expected_p0, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(1, lower_bound, lower_bound_expected_p1, n_entity_bound);
  MPI_CHECK_EQ_C_ARRAY(2, lower_bound, lower_bound_expected_p2, n_entity_bound);


  PDM_part_comm_graph_free(ptpgc);

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
  int *entity_bound  = ventity_bound         [i_rank].data();


  PDM_part_comm_graph_t* ptpgc = PDM_part_comm_graph_create(n_part,
                                                                            &n_entity_bound,
                                                                            &entity_bound,
                                                                            pdm_comm);

  const int* lower_bound = PDM_part_comm_graph_owner_get(ptpgc, 0);

  // PDM_log_trace_array_int(lower_bound, 3, "lower_bound ::");

  static int lower_bound_expected_p0[3] = {1, 1, 1};
  static int lower_bound_expected_p1[3] = {0, 0, 0};

  MPI_CHECK_EQ_C_ARRAY(0, lower_bound, lower_bound_expected_p0, 3);
  MPI_CHECK_EQ_C_ARRAY(1, lower_bound, lower_bound_expected_p1, 3);


  // std::vector<std::vector<int>> vold_to_new = {{0, 1, 2, 3, 4, 5, 6, 7, 8},
  //                                              {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
  std::vector<std::vector<int>> vold_to_new = {{8, 7, 6, 5, 4, 3, 2, 1, 0},
                                               {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
  int *old_to_new = vold_to_new[i_rank].data();

  // PDM_log_trace_array_int(entity_bound, 4 * n_entity_bound, "entity_bound (Avant) ::");

  PDM_part_comm_graph_reorder(ptpgc,
                                      &entity_bound,
                                      &old_to_new);

  // PDM_log_trace_array_int(entity_bound, 4 * n_entity_bound, "entity_bound ::");

  static int entity_bound_reorder_p0[12] = {7 , 1, 1, 12, 4, 1, 1, 8, 1, 1, 1, 4};
  static int entity_bound_reorder_p1[12] = {12, 0, 1,  7, 8, 0, 1, 4, 4, 0, 1, 1};

  MPI_CHECK_EQ_C_ARRAY(0, entity_bound, entity_bound_reorder_p0, 12);
  MPI_CHECK_EQ_C_ARRAY(1, entity_bound, entity_bound_reorder_p1, 12);

  PDM_part_comm_graph_free(ptpgc);
}
