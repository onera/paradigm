#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include <vector>

MPI_TEST_CASE("[pdm_gnum] - 1p - from_parent",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<PDM_g_num_t> parent_gnum = {8, 1, 12, 8, 4};
  int n_elmt = parent_gnum.size();

  PDM_gnum_set_from_parents(gen_gnum,
                            0,
                            n_elmt,
                            parent_gnum.data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected[5] = {3, 1, 4, 3, 2};

  CHECK_EQ_C_ARRAY(ln_to_gn, ln_to_gn_expected, n_elmt);

  PDM_gnum_free(gen_gnum);
}

MPI_TEST_CASE("[pdm_gnum] - 2p - from_parent", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<std::vector<PDM_g_num_t>> parent_gnum = {{8, 1, 12}, {8, 4}};
  int n_elmt = parent_gnum[i_rank].size();

  PDM_gnum_set_from_parents(gen_gnum,
                            0,
                            n_elmt,
                            parent_gnum[i_rank].data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[3] = {3, 1, 4};
  PDM_g_num_t ln_to_gn_expected_p1[2] = {3, 2};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, n_elmt);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, n_elmt);

  PDM_gnum_free(gen_gnum);

}



MPI_TEST_CASE("[pdm_gnum] - 1p - from_parent nuplet = 2 ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<PDM_g_num_t> parent_gnum = {8, 0, 1, 1, 12, 0, 8, -4, 4, 0};
  int n_elmt = parent_gnum.size()/2;

  PDM_gnum_set_parents_nuplet(gen_gnum, 2);

  PDM_gnum_set_from_parents(gen_gnum,
                            0,
                            n_elmt,
                            parent_gnum.data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected[5] = {3, 1, 5, 4, 2};

  CHECK_EQ_C_ARRAY(ln_to_gn, ln_to_gn_expected, n_elmt);

  PDM_gnum_free(gen_gnum);
}

MPI_TEST_CASE("[pdm_gnum] - 2p - from_parent nuplet = 2 ", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<std::vector<PDM_g_num_t>> parent_gnum = {{8, 0, 1, 1, 12, 0}, {8, -4, 4, 0}};
  int n_elmt = parent_gnum[i_rank].size()/2;

  PDM_gnum_set_parents_nuplet(gen_gnum, 2);

  PDM_gnum_set_from_parents(gen_gnum,
                            0,
                            n_elmt,
                            parent_gnum[i_rank].data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[3] = {3, 1, 5};
  PDM_g_num_t ln_to_gn_expected_p1[2] = {4, 2};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, n_elmt);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, n_elmt);


  PDM_gnum_free(gen_gnum);

}


MPI_TEST_CASE("[pdm_gnum] - 1p - from_coords",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<double> coords      = {0. ,  0. ,  0.,
                                     0.5,  0.5,  0.,
                                     -1.,  0. ,  0.,
                                      1.,  1. ,  1.,
                                     0.8, -1. , 0.5,
                                      1.,  1. ,  1.,}; // Merge with 4
  std::vector<double> char_length = {1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6, 1.e-6};
  int n_elmt = coords.size()/3;

  PDM_gnum_set_from_coords(gen_gnum,
                           0,
                           n_elmt,
                           coords.data(),
                           char_length.data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected[6] = {3, 4, 1, 5, 2, 5};

  CHECK_EQ_C_ARRAY(ln_to_gn, ln_to_gn_expected, n_elmt);

  PDM_gnum_free(gen_gnum);
}


MPI_TEST_CASE("[pdm_gnum] - 2p - from_coords",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_FALSE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<std::vector<double>> coords      = {{0. ,  0. ,  0.,
                                                  0.5,  0.5,  0.,
                                                  -1.,  0. ,  0.,
                                                   1.,  1. ,  1.},
                                                  {0.8, -1. , 0.5,
                                                   1.,  1. ,  1.,}}; // Merge with 4
  std::vector<std::vector<double>> char_length = {{1.e-6, 1.e-6, 1.e-6, 1.e-6}, {1.e-6, 1.e-6}};
  int n_elmt = coords[i_rank].size()/3;

  PDM_gnum_set_from_coords(gen_gnum,
                           0,
                           n_elmt,
                           coords     [i_rank].data(),
                           char_length[i_rank].data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[4] = {3, 4, 1, 6};
  PDM_g_num_t ln_to_gn_expected_p1[2] = {2, 5};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, n_elmt);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, n_elmt);

  PDM_gnum_free(gen_gnum);
}

MPI_TEST_CASE("[pdm_gnum] - 2p - from_coords - merge ",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);

  std::vector<std::vector<double>> coords      = {{0. ,  0. ,  0.,
                                                  0.5,  0.5,  0.,
                                                  -1.,  0. ,  0.,
                                                   1.,  1. ,  1.},
                                                  {0.8, -1. , 0.5,
                                                   1.,  1. ,  1.,}}; // Merge with 4
  std::vector<std::vector<double>> char_length = {{1.e-6, 1.e-6, 1.e-6, 1.e-6}, {1.e-6, 1.e-6}};
  int n_elmt = coords[i_rank].size()/3;

  PDM_gnum_set_from_coords(gen_gnum,
                           0,
                           n_elmt,
                           coords     [i_rank].data(),
                           char_length[i_rank].data());

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[4] = {3, 4, 1, 5};
  PDM_g_num_t ln_to_gn_expected_p1[2] = {2, 5};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, n_elmt);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, n_elmt);

  PDM_gnum_free(gen_gnum);
}



MPI_TEST_CASE("[pdm_gnum] - 2p - from_part_comm_graph 2p", 2) {
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

  PDM_part_comm_graph_t* pgc = PDM_part_comm_graph_create(n_part,
                                                            &n_entity_bound,
                                                            &entity_bound,
                                                            pdm_comm);

  std::vector<int> pn_elmt = {9, 12};
  const int n_elmt = pn_elmt[i_rank];


  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);


  PDM_gnum_set_from_part_comm_graph(gen_gnum,
                                    &n_elmt,
                                    pgc);

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);


  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[9]  = {1, 2, 7, 3, 4, 8, 5, 6, 9 };
  PDM_g_num_t ln_to_gn_expected_p1[12] = {7, 10, 11, 12, 8, 13, 14, 15, 9, 16, 17, 18};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, pn_elmt[i_rank]);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, pn_elmt[i_rank]);

  PDM_gnum_free(gen_gnum);

  PDM_part_comm_graph_free(pgc);
}



MPI_TEST_CASE("[pdm_gnum] - 2p - from_entity_graph 2p", 2) {
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
  // int n_part = 1;

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

  std::vector<int> pn_elmt = {9, 12};
  const int n_elmt = pn_elmt[i_rank];


  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1.e-6,
                                             pdm_comm,
                                             PDM_OWNERSHIP_KEEP);


  PDM_gnum_set_from_entity_graph(gen_gnum,
                                 0,
                                 n_elmt,
                                 n_entity_bound,
                                 entity_bound);

  PDM_gnum_compute(gen_gnum);

  PDM_g_num_t* ln_to_gn = PDM_gnum_get(gen_gnum, 0);


  if(0 == 1) {
    PDM_log_trace_array_long(ln_to_gn, n_elmt, "ln_to_gn");
  }

  PDM_g_num_t ln_to_gn_expected_p0[9]  = {1, 2, 7, 3, 4, 8, 5, 6, 9 };
  PDM_g_num_t ln_to_gn_expected_p1[12] = {7, 10, 11, 12, 8, 13, 14, 15, 9, 16, 17, 18};

  MPI_CHECK_EQ_C_ARRAY(0, ln_to_gn, ln_to_gn_expected_p0, pn_elmt[i_rank]);
  MPI_CHECK_EQ_C_ARRAY(1, ln_to_gn, ln_to_gn_expected_p1, pn_elmt[i_rank]);

  PDM_gnum_free(gen_gnum);

}
