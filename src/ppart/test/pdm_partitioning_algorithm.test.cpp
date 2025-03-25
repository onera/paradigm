#include <stddef.h>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_logging.h"
#include "pdm_partitioning_algorithm.h"
#include <vector>



MPI_TEST_CASE("[PDM_part_generate_entity_graph_comm] - 1p - n_part=3", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_rank = 0;
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<PDM_g_num_t> part_distribution = {0, 3};
  int n_part = 3;
  int n_part_tot = part_distribution[n_rank];

  std::vector<std::vector<PDM_g_num_t>> vpentity_ln_to_gn = {{1, 2, 3},
                                                             {3, 4, 5},
                                                             {3, 5, 7, 8}};
  std::vector<int> vpn_entity(n_part);

  PDM_g_num_t **pentity_ln_to_gn = NULL;
  PDM_malloc(pentity_ln_to_gn, n_part, PDM_g_num_t *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    pentity_ln_to_gn[i_part] = vpentity_ln_to_gn[i_part].data();
    vpn_entity      [i_part] = vpentity_ln_to_gn[i_part].size();
  }

  int **pproc_bound_idx  = NULL;
  int **ppart_bound_idx  = NULL;
  int **pentity_bound    = NULL;
  int **pentity_priority = NULL;
  PDM_part_generate_entity_graph_comm(pdm_comm,
                                      part_distribution.data(),
                                      NULL,
                                      n_part,
                                      vpn_entity.data(),
              (const PDM_g_num_t **)  pentity_ln_to_gn,
                                      NULL,
                                      &pproc_bound_idx,
                                      &ppart_bound_idx,
                                      &pentity_bound,
                                      &pentity_priority);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      PDM_log_trace_array_int(pproc_bound_idx [i_part], n_rank+1    , "pproc_bound_idx ::");
      PDM_log_trace_array_int(ppart_bound_idx [i_part], n_part_tot+1, "ppart_bound_idx ::");
      PDM_log_trace_array_int(pentity_bound   [i_part], 3 * ppart_bound_idx[i_part][n_part_tot], "pentity_bound ::");
      PDM_log_trace_array_int(pentity_priority[i_part], ppart_bound_idx[i_part][n_part_tot], "pentity_priority ::");
    }
  }

  // Check
  int n0_pproc_bound_idx [2] = {0, 2};
  int n0_ppart_bound_idx [4] = {0, 0, 1, 2};
  int n0_pentity_bound   [6] = {3, 0, 2, 1, 3, 0};
  int n0_pentity_priority[2] = {-1, -1};

  int n1_pproc_bound_idx [2] = {0, 3};
  int n1_ppart_bound_idx [4] = {0, 1, 1, 3};
  int n1_pentity_bound   [9] = {1, 0, 1, 3, 1, 0, 3, 1, 3};
  int n1_pentity_priority[3] = {2, -1, 1};

  int n2_pproc_bound_idx [2] = {0, 3};
  int n2_ppart_bound_idx [4] = {0, 1, 3, 3};
  int n2_pentity_bound   [9] = {1, 0, 1, 3, 1, 0, 2, 1, 2};
  int n2_pentity_priority[3] = {2, 1, -1};

  MPI_CHECK_EQ_C_ARRAY(0, pproc_bound_idx [0], n0_pproc_bound_idx , 2);
  MPI_CHECK_EQ_C_ARRAY(0, ppart_bound_idx [0], n0_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_bound   [0], n0_pentity_bound   , 6);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_priority[0], n0_pentity_priority, 2);

  MPI_CHECK_EQ_C_ARRAY(0, pproc_bound_idx [1], n1_pproc_bound_idx , 2);
  MPI_CHECK_EQ_C_ARRAY(0, ppart_bound_idx [1], n1_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_bound   [1], n1_pentity_bound   , 9);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_priority[1], n1_pentity_priority, 3);

  MPI_CHECK_EQ_C_ARRAY(0, pproc_bound_idx [2], n2_pproc_bound_idx , 2);
  MPI_CHECK_EQ_C_ARRAY(0, ppart_bound_idx [2], n2_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_bound   [2], n2_pentity_bound   , 9);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_priority[2], n2_pentity_priority, 3);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pproc_bound_idx [i_part]);
    PDM_free(ppart_bound_idx [i_part]);
    PDM_free(pentity_bound   [i_part]);
    PDM_free(pentity_priority[i_part]);
  }
  PDM_free(pproc_bound_idx );
  PDM_free(ppart_bound_idx );
  PDM_free(pentity_bound   );
  PDM_free(pentity_priority);

  PDM_free(pentity_ln_to_gn);
}


MPI_TEST_CASE("[PDM_part_generate_entity_graph_comm] - 2p - n_part_tot=3", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_rank = 0;
  int i_rank = 0;
  PDM_MPI_Comm_size(pdm_comm, &n_rank);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<PDM_g_num_t> part_distribution = {0, 2, 3};
  int n_part = part_distribution[i_rank+1]-part_distribution[i_rank];
  int n_part_tot = part_distribution[n_rank];

  std::vector<std::vector<std::vector<PDM_g_num_t>>> vpentity_ln_to_gn = {{{1, 2, 3},
                                                                          {3, 4, 5}},
                                                                         {{3, 5, 7, 8}}};
  std::vector<int> vpn_entity(n_part);

  PDM_g_num_t **pentity_ln_to_gn = NULL;
  PDM_malloc(pentity_ln_to_gn, n_part, PDM_g_num_t *);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    pentity_ln_to_gn[i_part] = vpentity_ln_to_gn[i_rank][i_part].data();
    vpn_entity      [i_part] = vpentity_ln_to_gn[i_rank][i_part].size();
  }

  int **pproc_bound_idx  = NULL;
  int **ppart_bound_idx  = NULL;
  int **pentity_bound    = NULL;
  int **pentity_priority = NULL;
  PDM_part_generate_entity_graph_comm(pdm_comm,
                                      part_distribution.data(),
                                      NULL,
                                      n_part,
                                      vpn_entity.data(),
              (const PDM_g_num_t **)  pentity_ln_to_gn,
                                      NULL,
                                      &pproc_bound_idx,
                                      &ppart_bound_idx,
                                      &pentity_bound,
                                      &pentity_priority);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      PDM_log_trace_array_int(pproc_bound_idx [i_part], n_rank+1    , "pproc_bound_idx ::");
      PDM_log_trace_array_int(ppart_bound_idx [i_part], n_part_tot+1, "ppart_bound_idx ::");
      PDM_log_trace_array_int(pentity_bound   [i_part], 3 * ppart_bound_idx[i_part][n_part_tot], "pentity_bound ::");
      PDM_log_trace_array_int(pentity_priority[i_part], ppart_bound_idx[i_part][n_part_tot], "pentity_priority ::");
    }
  }

  // Check
  int p0_n0_pproc_bound_idx [3] = {0, 1, 2};
  int p0_n0_ppart_bound_idx [4] = {0, 0, 1, 2};
  int p0_n0_pentity_bound   [6] = {3, 0, 2, 1, 3, 1};
  int p0_n0_pentity_priority[2] = {-1, -1};

  int p0_n1_pproc_bound_idx [3] = {0, 1, 3};
  int p0_n1_ppart_bound_idx [4] = {0, 1, 1, 3};
  int p0_n1_pentity_bound   [9] = {1, 0, 1, 3, 1, 1, 1, 1, 3};
  int p0_n1_pentity_priority[3] = {0, -1, 1};

  int p1_n2_pproc_bound_idx [3] = {0, 3, 3};
  int p1_n2_ppart_bound_idx [4] = {0, 1, 3, 3};
  int p1_n2_pentity_bound   [9] = {1, 0, 1, 3, 1, 0, 2, 1, 2};
  int p1_n2_pentity_priority[3] = {0, 1, -1};

  MPI_CHECK_EQ_C_ARRAY(0, pproc_bound_idx [0], p0_n0_pproc_bound_idx , 3);
  MPI_CHECK_EQ_C_ARRAY(0, ppart_bound_idx [0], p0_n0_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_bound   [0], p0_n0_pentity_bound   , 6);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_priority[0], p0_n0_pentity_priority, 2);

  MPI_CHECK_EQ_C_ARRAY(0, pproc_bound_idx [1], p0_n1_pproc_bound_idx , 3);
  MPI_CHECK_EQ_C_ARRAY(0, ppart_bound_idx [1], p0_n1_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_bound   [1], p0_n1_pentity_bound   , 9);
  MPI_CHECK_EQ_C_ARRAY(0, pentity_priority[1], p0_n1_pentity_priority, 3);

  MPI_CHECK_EQ_C_ARRAY(1, pproc_bound_idx [0], p1_n2_pproc_bound_idx , 3);
  MPI_CHECK_EQ_C_ARRAY(1, ppart_bound_idx [0], p1_n2_ppart_bound_idx , 4);
  MPI_CHECK_EQ_C_ARRAY(1, pentity_bound   [0], p1_n2_pentity_bound   , 9);
  MPI_CHECK_EQ_C_ARRAY(1, pentity_priority[0], p1_n2_pentity_priority, 3);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pproc_bound_idx [i_part]);
    PDM_free(ppart_bound_idx [i_part]);
    PDM_free(pentity_bound   [i_part]);
    PDM_free(pentity_priority[i_part]);
  }
  PDM_free(pproc_bound_idx );
  PDM_free(ppart_bound_idx );
  PDM_free(pentity_bound   );
  PDM_free(pentity_priority);

  PDM_free(pentity_ln_to_gn);
}
