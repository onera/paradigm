#include <stddef.h>
#include <vector>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"



MPI_TEST_CASE("[pdm_multipart] - 2p - 1 domain, with dpart_id", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  const int n_domain = 1;
  const int n_part   = 1;
  const PDM_split_dual_t split_method     = PDM_SPLIT_DUAL_WITH_IMPLICIT; // not relevant here
  const PDM_part_size_t  part_size_method = PDM_PART_SIZE_HOMOGENEOUS;    // not relevant here
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                split_method,
                                                part_size_method,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);


  // Set block-distributed mesh
  int dn_face = 2;
  int dn_edge = 6;
  int dn_vtx[] = {5, 4};

  int dface_edge_idx[] = {0, 4, 8};
  std::vector<PDM_g_num_t> dface_edge[] = {
    {1, 8, -3, -7, 2, 9, -4, -8},
    {3, 11, -5, -10, 4, 12, -6, -11}
  };

  int dedge_vtx_idx[] = {0, 2, 4, 6, 8, 10, 12};
  std::vector<PDM_g_num_t> dedge_vtx[] = {
    {1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9},
    {1, 4, 2, 5, 3, 6, 4, 7, 5, 8, 6, 9}
  };

  std::vector<double> dvtx_coord[] = {
    {
      0, 0, 0,
      1, 0, 0,
      2, 0, 0,
      0, 1, 0,
      1, 1, 0
    },
    {
      2, 1, 0,
      0, 2, 0,
      1, 2, 0,
      2, 2, 0
    }
  };


  int n_edge_group = 0;
  int dedge_group_idx[] = {0};
  PDM_multipart_block_set(mpart,
                          0,
                          dn_face,
                          dn_edge,
                          dn_vtx[test_rank],
                          n_edge_group,
                          dface_edge_idx,
                          dface_edge[test_rank].data(),
                          NULL,
                          dedge_vtx_idx,
                          dedge_vtx[test_rank].data(),
                          dvtx_coord[test_rank].data(),
                          dedge_group_idx,
                          NULL);

  // Set desired partitioning
  PDM_g_num_t distrib_face[] = {0, 2, 4};

  PDM_g_num_t gn_part = 2;
  int *dpart_id = NULL;
  PDM_malloc(dpart_id, dn_face, int);
  for (int i_face = 0; i_face < dn_face; i_face++) {
    dpart_id[i_face] = (i_face + distrib_face[test_rank] + 1) % gn_part;
  }

  PDM_multipart_dpart_id_set(mpart,
                             0,
                             dpart_id,
                             PDM_OWNERSHIP_USER);

  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *face_ln_to_gn = NULL;
  int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_CELL, // /!\ FACE = CELL here
                                               &face_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);

  // Check face global IDs
  static PDM_g_num_t expected_face_ln_to_gn0[] = {2, 4};
  static PDM_g_num_t expected_face_ln_to_gn1[] = {1, 3};

  MPI_CHECK_EQ_C_ARRAY(0, face_ln_to_gn, expected_face_ln_to_gn0, n_face);
  MPI_CHECK_EQ_C_ARRAY(1, face_ln_to_gn, expected_face_ln_to_gn1, n_face);

  // Free memory
  PDM_free(dpart_id);
  PDM_multipart_free(mpart);
}