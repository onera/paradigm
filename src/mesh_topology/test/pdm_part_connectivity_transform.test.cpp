#include <stdlib.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm_doctest.h"
#include "pdm_mem_tool.h"
#include "pdm_part_connectivity_transform.h"

/*
 *  Use case
 *     ./paradigm/test/pdm_t_partitioning_dcube -n 3 -n_part 1 -pt-scotch (Générateur dcube_gen sur 8 cellules)
 */

MPI_TEST_CASE("[pdm_part_connectivity_transform] - 1p - pdm_connectivity_transform ",1) {

  const int n_cell = 8;
  const int n_face = 36;

  int face_cell_idx[37] = {0 , 1, 2, 3, 4, 6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26,
                           28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 44, 45, 46, 47, 48};
  int face_cell    [48] = {1, 6, 7, 4, 1, 5, 6, 2, 7, 3, 4, 8, 5, 2, 3, 8, 1, 7, 5, 3, 1, 6, 7,
                           4, 5, 2, 3, 8, 6, 4, 2, 8, 1, 5, 6, 2, 1, 7, 5, 3, 6, 4, 2, 8, 7, 3, 4, 8};

  int cell_face_idx[n_cell+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
  int cell_face    [48]       = { 1,  5,  13, 17, 25, 29,
                                 -6, 10, -19, 23, 28, 32,
                                 -7, 11,  16, 20,-30, 34,
                                  4,  8, -18, 22,-31, 35,
                                 -5,  9,  15, 19, 26, 30,
                                  2,  6, -17, 21, 27, 31,
                                  3,  7, 14, 18, -29, 33,
                                 -8, 12,-20, 24, -32, 36};
  int* cell_cell_idx;
  int* cell_cell;
  PDM_combine_connectivity(n_cell,
                           cell_face_idx,
                           cell_face,
                           face_cell_idx,
                           face_cell,
                           &cell_cell_idx,
                           &cell_cell);

  // PDM_log_trace_array_int(cell_cell_idx, n_cell+1         , "cell_cell_idx::");
  // PDM_log_trace_array_int(cell_cell, cell_cell_idx[n_cell], "cell_cell::");

  int cell_cell_idx_expected[n_cell+1] = {0, 4, 8, 12, 16, 20, 24, 28, 32};
  int cell_cell_expected[32]           = {1, 5, 6, 7, 2, 5, 6, 8, 3, 5, 7, 8, 4, 6, 7, 8, 1, 2, 3, 5, 1, 2, 4, 6, 1, 3, 4, 7, 2, 3, 4, 8};

  CHECK_EQ_C_ARRAY(cell_cell_idx, cell_cell_idx_expected, n_cell+1);
  CHECK_EQ_C_ARRAY(cell_cell    , cell_cell_expected    , cell_cell_idx[n_cell]);

  free(cell_cell_idx);
  free(cell_cell);


  int* cell_face_from_transpose_idx;
  int* cell_face_from_transpose;
  PDM_connectivity_transpose(n_face,
                             n_cell,
                             face_cell_idx,
                             face_cell,
                             &cell_face_from_transpose_idx,
                             &cell_face_from_transpose);

  int cell_face_from_transpose_idx_expected[n_cell+1] = {0, 6, 12, 18, 24, 30, 36, 42, 48};
  int cell_face_from_transpose_expected    [48]       = {1, 5, 13, 17, 25, 29, 6, 10, 19, 23, 28, 32, 7, 11, 16, 20, 30, 34, 4, 8, 18, 22, 31, 35, 5, 9, 15, 19, 26, 30, 2, 6, 17, 21, 27, 31, 3, 7, 14, 18, 29, 33, 8, 12, 20, 24, 32, 36};

  CHECK_EQ_C_ARRAY(cell_face_from_transpose_idx, cell_face_from_transpose_idx_expected, n_cell+1);
  CHECK_EQ_C_ARRAY(cell_face_from_transpose    , cell_face_from_transpose_expected    , cell_face_from_transpose_idx_expected[n_cell]);

  PDM_free(cell_face_from_transpose_idx);
  PDM_free(cell_face_from_transpose);

  int* face_cell_from_transpose_idx;
  int* face_cell_from_transpose;
  PDM_connectivity_transpose(n_cell,
                             n_face,
                             cell_face_idx,
                             cell_face,
                             &face_cell_from_transpose_idx,
                             &face_cell_from_transpose);

  int face_cell_from_transpose_idx_expected[n_face+1] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 44, 45, 46, 47, 48};
  int face_cell_from_transpose_expected    [48]       = {1, 6, 7, 4, 1, -5, -2, 6, -3, 7, 4, -8, 5, 2, 3, 8, 1, 7, 5, 3, 1, -6, -4, 7, -2, 5, 3, -8, 6, 4, 2, 8, 1, 5, 6, 2, 1, -7, -3, 5, -4, 6, 2, -8, 7, 3, 4, 8};

  CHECK_EQ_C_ARRAY(face_cell_from_transpose_idx, face_cell_from_transpose_idx_expected, n_face+1);
  CHECK_EQ_C_ARRAY(face_cell_from_transpose    , face_cell_from_transpose_expected    , face_cell_from_transpose_idx_expected[n_face]);

  PDM_free(face_cell_from_transpose_idx);
  PDM_free(face_cell_from_transpose);

}



MPI_TEST_CASE("[PDM_compute_dface_vtx_from_edges] - 2p",2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<std::vector<int>> dface_edge_idx = {{0, 4, 8, 12, 16, 20, 24},
                                                  {0, 4, 8, 12, 16, 20}};

  std::vector<std::vector<PDM_g_num_t>> dface_edge = {{1, 2, 4, 6, -4, 3, 7, 9, -1, 5, 8, 12, -5, -2, 11, 14, -8, -3, 10, 15, -8, -4, 13, 17},
                                                      {-13, -11, -6, 18, -10, -7, 16, 19, -16, -9, 13, 20, -18, -17, -14, -12, -20, -19, -15, 17}};

  std::vector<std::vector<PDM_g_num_t>> dedge_vtx = {{2, 1, 1, 4, 3, 2, 5, 2, 7, 1, 4, 5, 6, 3, 2, 8, 5, 6, 3, 9},
                                                     {10, 4, 8, 7, 5, 11, 7, 10, 9, 8, 6, 12, 11, 8, 10, 11, 12, 9, 11, 12}};

  int dn_face = dface_edge_idx[i_rank].size()-1;
  int dn_edge = dedge_vtx     [i_rank].size()/2;

  PDM_g_num_t *dface_vtx = NULL;
  PDM_compute_dface_vtx_from_edges(pdm_comm,
                                   dn_face,
                                   dn_edge,
                                   dface_edge_idx[i_rank].data(),
                                   dface_edge    [i_rank].data(),
                                   dedge_vtx     [i_rank].data(),
                                   &dface_vtx);

  PDM_g_num_t p0_expected_dface_vtx[24] = {2, 1, 4, 5, 2, 5, 6, 3, 1, 2, 8, 7, 1, 7, 10, 4, 8, 2, 3, 9, 8, 2, 5, 11};
  PDM_g_num_t p1_expected_dface_vtx[20] = {11, 5, 4, 10, 9, 3, 6, 12, 12, 6, 5, 11, 11, 10, 7, 8, 12, 11, 8, 9};

  MPI_CHECK_EQ_C_ARRAY(0, dface_vtx, p0_expected_dface_vtx, 24);
  MPI_CHECK_EQ_C_ARRAY(1, dface_vtx, p1_expected_dface_vtx, 20);


  PDM_free(dface_vtx);
}




MPI_TEST_CASE("[PDM_graph_compress] ", 1) {

  int n_entity = 5;
  std::vector<int> graph_idx = {0, 5, 10, 15, 20, 26};
  std::vector<int> graph     = {4, 5, 2, 1, 1,
                                1, 2, 3, 5, 3,
                                4, 3, 2, 4, 5,
                                5, 4, 1, 3, 5,
                                1, 2, 3, 4, 5, 5};

  PDM_graph_compress(n_entity, graph_idx.data(), graph.data());

  // PDM_log_trace_array_int(graph_idx.data(), n_entity+1, "graph_idx ::");
  // PDM_log_trace_array_int(graph.data(), graph_idx[n_entity], "graph ::");

  int expected_graph_idx[6 ] = {0, 3, 6, 9, 12, 16};
  int expected_graph    [16] = {1, 3, 4,
                                0, 2, 4,
                                1, 3, 4,
                                0, 2, 4,
                                0, 1, 2, 3};

  MPI_CHECK_EQ_C_ARRAY(0, graph_idx, expected_graph_idx, 6);
  MPI_CHECK_EQ_C_ARRAY(0, graph    , expected_graph    , 16);


}
