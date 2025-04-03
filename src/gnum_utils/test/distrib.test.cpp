#include <memory>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
#include "pdm_mem_tool.h"


MPI_TEST_CASE("[PDM_distrib_compute] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<int> pdn_elmt = {4, 8};

  int dn_elmt = pdn_elmt[i_rank];
  std::vector<PDM_g_num_t> distrib(n_rank+1);
  PDM_distrib_compute(dn_elmt, distrib.data(), -1, pdm_comm);

  PDM_g_num_t expected_distrib[3] = {0, 4, 12};

  CHECK_EQ_C_ARRAY(distrib.data(), expected_distrib, n_rank+1);

}

MPI_TEST_CASE("[PDM_compute_entity_distribution] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  std::vector<int> pdn_elmt = {4, 8};
  int dn_elmt = pdn_elmt[i_rank];

  PDM_g_num_t* distrib = PDM_compute_entity_distribution(pdm_comm, dn_elmt);

  PDM_free(distrib);
}

MPI_TEST_CASE("[PDM_compute_uniform_entity_distribution] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  PDM_g_num_t n_g_tot = 100;

  PDM_g_num_t* distrib = PDM_compute_uniform_entity_distribution(pdm_comm, n_g_tot);

  PDM_g_num_t expected_distrib[3] = {0, 50, 100};

  CHECK_EQ_C_ARRAY(distrib, expected_distrib, n_rank+1);

  PDM_free(distrib);
}



MPI_TEST_CASE("[PDM_compute_uniform_entity_distribution_from_partition] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  int n_part_in = 1;
  std::vector<std::vector<PDM_g_num_t>> pln_to_gn = {{40, 12, 1, 21}, {40, 1, 2, 3, 5}};

  int          n_elmt   = pln_to_gn[i_rank].size();
  PDM_g_num_t *ln_to_gn = pln_to_gn[i_rank].data();

  PDM_g_num_t *distrib = PDM_compute_uniform_entity_distribution_from_partition(pdm_comm,
                                                                                n_part_in,
                                                                                &n_elmt,
                                                        (const PDM_g_num_t **)  &ln_to_gn);

  // PDM_log_trace_array_long(distrib, n_rank+1, "distrib ::");

  PDM_g_num_t expected_distrib[3] = {0, 20, 40};

  CHECK_EQ_C_ARRAY(distrib, expected_distrib, n_rank+1);

  PDM_free(distrib);
}


MPI_TEST_CASE("[PDM_distrib_weight] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  int n_part_in = 1;
  std::vector<std::vector<PDM_g_num_t>> pln_to_gn = {{40 , 12, 1,  21}, {40 , 1 , 2   , 3 , 5}};
  std::vector<std::vector<double>>      pweight   = {{10., 1., 1., 1.}, {10., 1., 100., 1., 1.}};

  int          n_elmt   = pln_to_gn[i_rank].size();
  PDM_g_num_t *ln_to_gn = pln_to_gn[i_rank].data();
  double      *weight   = pweight  [i_rank].data();

  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol          = 0.10;
  PDM_g_num_t* distrib = NULL;
  PDM_distrib_weight(    sampling_factor,
                         n_rank,
                         n_part_in,
                         &n_elmt,
 (const PDM_g_num_t **)  &ln_to_gn,
 (const double      **)  &weight,
                         n_iter_max,
                         tol,
                         pdm_comm,
                         &distrib);

  // PDM_log_trace_array_long(distrib, n_rank+1, "distrib ::");

  PDM_g_num_t expected_distrib[3] = {0, 6, 40};

  CHECK_EQ_C_ARRAY(distrib, expected_distrib, n_rank+1);

  PDM_free(distrib);
}
