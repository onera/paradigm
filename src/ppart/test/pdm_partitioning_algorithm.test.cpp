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




MPI_TEST_CASE("[PDM_compute_face_edge_from_face_vtx] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_rank = 0;
  int i_rank = 0;
  PDM_MPI_Comm_size(pdm_comm, &n_rank);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  int n_part = 1;

  // Test cases generated on 2 procs with test : test/pdm_t_compute_part_edges -n 3 (with n_vtx_y=n_vtx_z=2)
  std::vector<std::vector<int>> vpface_vtx_idx = {{0, 4, 8, 12, 16, 20, 24},
                                                  {0, 4, 8, 12, 16, 20, 24}};

  std::vector<std::vector<int>> vpface_vtx = {{2, 1, 3, 4, 1, 2, 6, 5, 1, 5, 7, 3, 6, 2, 4, 8, 8, 4, 3, 7, 8, 7, 5, 6},
                                              {1, 3, 4, 2, 5, 1, 2, 6, 5, 1, 3, 7, 6, 2, 4, 8, 8, 4, 3, 7, 8, 7, 5, 6}};

  std::vector<std::vector<PDM_g_num_t>> vpface_ln_to_gn = {{1, 3, 4, 6, 7, 10},
                                                           {2, 5, 6, 8, 9, 11}};

  std::vector<std::vector<PDM_g_num_t>> vpvtx_ln_to_gn = {{1, 2, 4, 5, 7, 8, 10, 11},
                                                          {2, 3, 5, 6, 8, 9, 11, 12}};

  int pn_face = vpface_ln_to_gn[i_rank].size();
  int pn_vtx  = vpvtx_ln_to_gn [i_rank].size();

  int         *pface_vtx_idx  = vpface_vtx_idx [i_rank].data();
  int         *pface_vtx      = vpface_vtx     [i_rank].data();
  PDM_g_num_t *pface_ln_to_gn = vpface_ln_to_gn[i_rank].data();
  PDM_g_num_t *pvtx_ln_to_gn  = vpvtx_ln_to_gn [i_rank].data();

  int           *pn_edge        = NULL;
  int          **pface_edge_idx = NULL;
  int          **pface_edge     = NULL;
  int          **pedge_vtx      = NULL;
  PDM_g_num_t  **pedge_ln_to_gn = NULL;
  PDM_compute_face_edge_from_face_vtx(pdm_comm,
                                      n_part,
                                      &pn_face,
                                      &pn_vtx,
                                      &pface_vtx_idx,
                                      &pface_vtx,
                                      &pface_ln_to_gn,
                                      &pvtx_ln_to_gn,
                                      &pface_edge_idx,
                                      &pface_edge,
                                      &pn_edge,
                                      &pedge_vtx,
                                      &pedge_ln_to_gn);

  int         p0_expected_pface_edge_idx[ 7] = {0, 4, 8, 12, 16, 20, 24};
  int         p0_expected_pface_edge    [24] = {1, 2, 3, 5, -1, 4, 6, 8, -4, -2, 7, 10, -6, -3, 9, 11, -9, -7, -5, 12, -12, -11, -10, -8};
  int         p0_expected_pedge_vtx     [24] = {2, 1, 1, 3, 4, 2, 5, 1, 3, 4, 2, 6, 7, 3, 6, 5, 4, 8, 5, 7, 8, 6, 7, 8};
  PDM_g_num_t p0_expected_pedge_ln_to_gn[12] = {1, 2, 4, 5, 6, 8, 11, 12, 13, 14, 17, 18};

  int         p1_expected_pface_edge_idx[ 7] = {0, 4, 8, 12, 16, 20, 24};
  int         p1_expected_pface_edge    [24] = {-2, 1, 3, 5, -4, -1, 6, 8, -4, -2, 7, 10, -6, -3, 9, 11, -9, -5, 7, 12, -12, -11, -8, 10};
  int         p1_expected_pedge_vtx     [24] = {2, 1, 3, 1, 4, 2, 1, 5, 3, 4, 2, 6, 3, 7, 6, 5, 4, 8, 7, 5, 8, 6, 7, 8};
  PDM_g_num_t p1_expected_pedge_ln_to_gn[12] = {3, 4, 7, 8, 9, 10, 13, 15, 16, 17, 19, 20};

  CHECK(pn_edge[0] == 12);

  MPI_CHECK_EQ_C_ARRAY(0, pface_edge_idx[0], p0_expected_pface_edge_idx,  7);
  MPI_CHECK_EQ_C_ARRAY(0, pface_edge    [0], p0_expected_pface_edge    , 24);
  MPI_CHECK_EQ_C_ARRAY(0, pedge_vtx     [0], p0_expected_pedge_vtx     , 24);
  MPI_CHECK_EQ_C_ARRAY(0, pedge_ln_to_gn[0], p0_expected_pedge_ln_to_gn, 12);

  MPI_CHECK_EQ_C_ARRAY(1, pface_edge_idx[0], p1_expected_pface_edge_idx,  7);
  MPI_CHECK_EQ_C_ARRAY(1, pface_edge    [0], p1_expected_pface_edge    , 24);
  MPI_CHECK_EQ_C_ARRAY(1, pedge_vtx     [0], p1_expected_pedge_vtx     , 24);
  MPI_CHECK_EQ_C_ARRAY(1, pedge_ln_to_gn[0], p1_expected_pedge_ln_to_gn, 12);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pface_edge_idx[i_part]);
    PDM_free(pface_edge    [i_part]);
    PDM_free(pedge_vtx     [i_part]);
    PDM_free(pedge_ln_to_gn[i_part]);
  }
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge    );
  PDM_free(pedge_vtx     );
  PDM_free(pedge_ln_to_gn);
  PDM_free(pn_edge);
}
