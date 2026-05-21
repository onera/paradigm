#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include <vector>
#include "pdm.h"
#include "pdm_dmesh.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_mem_tool.h"


MPI_TEST_CASE("[PDM_dmesh] PDM_dmesh_find_topological_ridges", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  /**
   *                12
   *              / |
   *             /  |
   *           11   |
   *            | 3 |
   *       6 ---|-- 7 ------ 8 ------ 9 ------ 10
   *      /     |  /        /        /        /
   *     /   1  | /   2    /   4    /   5    /
   *    /       |/        /        /        /
   *   1 ------ 2 ------ 3 ------ 4 ------ 5
   */

  // Faces
  PDM_g_num_t distrib_face[3] = {0, 3, 5};

  std::vector<std::vector<int>> vdface_vtx_idx = {
    {0, 4, 8, 12},
    {0, 4, 8}
  };

  std::vector<std::vector<PDM_g_num_t>> vdface_vtx = {
    {
      1, 2, 7,  6,
      2, 3, 8,  7,
      2, 7, 12, 11
    },
    {
      3, 4, 9,  8,
      4, 5, 10, 9
    }
  };

  int         *dface_vtx_idx = vdface_vtx_idx[i_rank].data();
  PDM_g_num_t *dface_vtx     = vdface_vtx    [i_rank].data();

  // Groups
  int n_group_face = 3;

  std::vector<std::vector<int>> vdgroup_face_idx = {
    {0, 1, 3, 4},
    {0, 0, 1, 1}
  };

  std::vector<std::vector<PDM_g_num_t>> vdgroup_face = {
    {1, 2, 4, 3},
    {5}
  };

  int         *dgroup_face_idx = vdgroup_face_idx[i_rank].data();
  PDM_g_num_t *dgroup_face     = vdgroup_face    [i_rank].data();


  // Find ridges
  PDM_g_num_t *distrib_ridge         = NULL;
  PDM_g_num_t *dridge_vtx            = NULL;
  int          n_group_ridge         = 0;
  int         *dgroup_edge_idx       = NULL;
  PDM_g_num_t *dgroup_edge           = NULL;
  int         *dridge_face_group_idx = NULL;
  int         *dridge_face_group     = NULL;
  PDM_dmesh_find_topological_ridges(pdm_comm,
                                    distrib_face,
                                    dface_vtx_idx,
                                    dface_vtx,
                                    n_group_face,
                                    dgroup_face_idx,
                                    dgroup_face,
                                    &distrib_ridge,
                                    &dridge_vtx,
                                    &n_group_ridge,
                                    &dgroup_edge_idx,
                                    &dgroup_edge,
                                    &dridge_face_group_idx,
                                    &dridge_face_group);

  // Check result
  int dn_ridge = (int) (distrib_ridge[i_rank+1] - distrib_ridge[i_rank]);

  PDM_g_num_t exp_distrib_ridge[3] = {0, 6, 14};
  std::vector<std::vector<PDM_g_num_t>> vexp_dridge_vtx = {
    {1, 2, 2, 3, 6, 1, 3, 4, 2, 7, 4, 5},
    {7, 6, 11, 2, 8, 7, 5, 10, 9, 8, 7, 12, 10, 9, 12, 11}
  };

  int exp_n_group_ridge = 4;

  std::vector<std::vector<int>> vexp_dgroup_edge_idx = {
    {0, 2, 5, 5, 6},
    {0, 1, 5, 8, 8}
  };
  std::vector<std::vector<PDM_g_num_t>> vexp_dgroup_edge = {
    {1, 3, 2, 4, 6, 5},
    {7, 9, 10, 11, 13, 8, 12, 14}
  };

  std::vector<std::vector<int>> vexp_dridge_face_group_idx = {
    {0, 1, 2, 3, 4, 7, 8},
    {0, 1, 2, 3, 4, 5, 6, 7, 8}
  };
  std::vector<std::vector<int>> vexp_dridge_face_group = {
    {1, 2, 1, 2, 1, 2, 3, 2},
    {1, 3, 2, 2, 2, 3, 2, 3}
  };


  MPI_CHECK_EQ_C_ARRAY(i_rank, distrib_ridge, exp_distrib_ridge, n_rank+1);
  CHECK(n_group_ridge == exp_n_group_ridge);
  MPI_CHECK_EQ_C_ARRAY(i_rank, dgroup_edge_idx, vexp_dgroup_edge_idx[i_rank].data(), n_group_ridge+1);
  MPI_CHECK_EQ_C_ARRAY(i_rank, dgroup_edge,     vexp_dgroup_edge    [i_rank].data(), dgroup_edge_idx[n_group_ridge]);
  MPI_CHECK_EQ_C_ARRAY(i_rank, dridge_face_group_idx, vexp_dridge_face_group_idx[i_rank].data(), dn_ridge+1);
  MPI_CHECK_EQ_C_ARRAY(i_rank, dridge_face_group,     vexp_dridge_face_group    [i_rank].data(), dridge_face_group_idx[dn_ridge]);

  // Free memory
  PDM_free(distrib_ridge);
  PDM_free(dridge_vtx);
  PDM_free(dgroup_edge_idx);
  PDM_free(dgroup_edge);
  PDM_free(dridge_face_group_idx);
  PDM_free(dridge_face_group);
}