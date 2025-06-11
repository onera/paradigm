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



MPI_TEST_CASE("[pdm_multipart] - 2p - 2D, 1 domain, with dpart_id", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 2D, 1 domain, with dpart_id\n");

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
    {
      1, 8, -3, 7,
      2, 9, -4, -8
    },
    {
      3, 11, 5, 10,
      4, 12, 6, -11
    }
  };

  int dedge_vtx_idx[] = {0, 2, 4, 6, 8, 10, 12};
  std::vector<PDM_g_num_t> dedge_vtx[] = {
    {
      1, 2,
      2, 3,
      4, 5,
      5, 6,
      8, 7,
      9, 8
    },
    {
      4, 1,
      2, 5,
      3, 6,
      7, 4,
      5, 8,
      6, 9
    }
  };

  std::vector<PDM_g_num_t> dedge_face[] = {
    {
      1, 0,
      2, 0,
      3, 1,
      4, 2,
      3, 0,
      4, 0
    },
    {
      1, 0,
      1, 2,
      2, 0,
      3, 0,
      3, 4,
      4, 0
    }
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

  SUBCASE("From dface_edge") {
    if (test_rank) printf("  From dface_edge\n");
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
  }

  SUBCASE("From dedge_face") {
    if (test_rank) printf("  From dedge_face\n");
    PDM_multipart_block_set(mpart,
                            0,
                            dn_face,
                            dn_edge,
                            dn_vtx[test_rank],
                            n_edge_group,
                            NULL,
                            NULL,
                            dedge_face[test_rank].data(),
                            dedge_vtx_idx,
                            dedge_vtx[test_rank].data(),
                            dvtx_coord[test_rank].data(),
                            dedge_group_idx,
                            NULL);
  }

  // Set desired partitioning
  std::vector<int> dpart_id[] = {
    {0, 1},
    {0, 1}
  };

  PDM_multipart_dpart_id_set(mpart,
                             0,
                             dpart_id[test_rank].data(),
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

  int expected_n_face = 2;
  CHECK(n_face == expected_n_face);

  // Check face global IDs
  static PDM_g_num_t expected_face_ln_to_gn0[] = {1, 3};
  static PDM_g_num_t expected_face_ln_to_gn1[] = {2, 4};

  MPI_CHECK_EQ_C_ARRAY(0, face_ln_to_gn, expected_face_ln_to_gn0, n_face);
  MPI_CHECK_EQ_C_ARRAY(1, face_ln_to_gn, expected_face_ln_to_gn1, n_face);

  // Free memory
  PDM_multipart_free(mpart);
}


MPI_TEST_CASE("[pdm_multipart] - 2p - 1D, 1 domain, with dpart_id", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 1D, 1 domain, with dpart_id\n");

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
  int dn_edge[] = {3, 2};
  int dn_vtx [] = {3, 3};

  std::vector<PDM_g_num_t> dedge_vtx[] = {
    {
      1, 2,
      2, 3,
      3, 4
    },
    {
      4, 5,
      5, 6
    }
  };

  std::vector<double> dvtx_coord[] = {
    {
      0, 0, 0,
      1, 0, 0,
      2, 0, 0
    },
    {
      3, 0, 0,
      4, 0, 0,
      5, 0, 0
    }
  };

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_USER,
                                        0,
                                        0,
                                        dn_edge[test_rank],
                                        dn_vtx [test_rank],
                                        pdm_comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord[test_rank].data(),
                          PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                            dedge_vtx[test_rank].data(),
                            NULL,
                            PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(mpart,
                          0,
                          dmesh);


  // Set desired partitioning
  std::vector<int> dpart_id[] = {
    {1, 0, 1},
    {0, 1}
  };

  PDM_multipart_dpart_id_set(mpart,
                             0,
                             dpart_id[test_rank].data(),
                             PDM_OWNERSHIP_USER);

  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *edge_ln_to_gn = NULL;
  int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_EDGE,
                                               &edge_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);

  int expected_n_edge[] = {2, 3};
  CHECK(n_edge == expected_n_edge[test_rank]);

  // Check edge global IDs
  static PDM_g_num_t expected_edge_ln_to_gn0[] = {2, 4};
  static PDM_g_num_t expected_edge_ln_to_gn1[] = {1, 3, 5};

  MPI_CHECK_EQ_C_ARRAY(0, edge_ln_to_gn, expected_edge_ln_to_gn0, n_edge);
  MPI_CHECK_EQ_C_ARRAY(1, edge_ln_to_gn, expected_edge_ln_to_gn1, n_edge);

  // Free memory
  PDM_dmesh_free(dmesh);
  PDM_multipart_free(mpart);
}



MPI_TEST_CASE("[pdm_multipart] - 2p - 2D, 1 domain, one empty block", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 2D, 1 domain, one empty block\n");

  const int n_domain = 1;
  const int n_part   = 1;
  const PDM_part_size_t part_size_method = PDM_PART_SIZE_HOMOGENEOUS;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;


#ifdef PDM_HAVE_PARMETIS
  SUBCASE("PDM_SPLIT_DUAL_WITH_PARMETIS") {
    split_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_PARMETIS\n");
  }
#endif
#ifdef PDM_HAVE_PTSCOTCH
  SUBCASE("PDM_SPLIT_DUAL_WITH_PTSCOTCH") {
    split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_PTSCOTCH\n");
  }
#endif
  SUBCASE("PDM_SPLIT_DUAL_WITH_HILBERT") {
    split_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_HILBERT\n");
  }
  SUBCASE("PDM_SPLIT_DUAL_WITH_IMPLICIT") {
    split_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_IMPLICIT\n");
  }

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                split_method,
                                                part_size_method,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);


  // Set block-distributed mesh
  int dn_face[] = {4,  0};
  int dn_edge[] = {12, 0};
  int dn_vtx [] = {9,  0};

  std::vector<int> dface_edge_idx[] = {
    {0, 4, 8, 12, 16},
    {0}
  };
  std::vector<PDM_g_num_t> dface_edge[] = {
    {
      1, 8, -3, 7,
      2, 9, -4, -8,
      3, 11, 5, 10,
      4, 12, 6, -11
    },
    {}
  };

  std::vector<PDM_g_num_t> dedge_vtx[] = {
    {
      1, 2,
      2, 3,
      4, 5,
      5, 6,
      8, 7,
      9, 8,
      4, 1,
      2, 5,
      3, 6,
      7, 4,
      5, 8,
      6, 9
    },
    {}
  };

  std::vector<double> dvtx_coord[] = {
    {
      0, 0, 0,
      1, 0, 0,
      2, 0, 0,
      0, 1, 0,
      1, 1, 0,
      2, 1, 0,
      0, 2, 0,
      1, 2, 0,
      2, 2, 0
    },
    {}
  };

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_USER,
                                        0,
                                        dn_face[test_rank],
                                        dn_edge[test_rank],
                                        dn_vtx [test_rank],
                                        pdm_comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord[test_rank].data(),
                          PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             dface_edge    [test_rank].data(),
                             dface_edge_idx[test_rank].data(),
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             dedge_vtx[test_rank].data(),
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(mpart,
                          0,
                          dmesh);


  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *face_ln_to_gn = NULL;
  int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_FACE,
                                               &face_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);
  PDM_UNUSED(n_face);

  // Free memory
  PDM_dmesh_free(dmesh);
  PDM_multipart_free(mpart);
}



MPI_TEST_CASE("[pdm_multipart] - 2p - 1D, 1 domain, one empty block", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 1D, 1 domain, one empty block\n");

  const int n_domain = 1;
  const int n_part   = 1;
  const PDM_part_size_t part_size_method = PDM_PART_SIZE_HOMOGENEOUS;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;


#ifdef PDM_HAVE_PARMETIS
  SUBCASE("PDM_SPLIT_DUAL_WITH_PARMETIS") {
    split_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_PARMETIS\n");
  }
#endif
#ifdef PDM_HAVE_PTSCOTCH
  SUBCASE("PDM_SPLIT_DUAL_WITH_PTSCOTCH") {
    split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_PTSCOTCH\n");
  }
#endif
  SUBCASE("PDM_SPLIT_DUAL_WITH_HILBERT") {
    split_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_HILBERT\n");
  }
  SUBCASE("PDM_SPLIT_DUAL_WITH_IMPLICIT") {
    split_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;
    if (test_rank == 0) printf("  PDM_SPLIT_DUAL_WITH_IMPLICIT\n");
  }

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                split_method,
                                                part_size_method,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);


  // Set block-distributed mesh
  int dn_edge[] = {5, 0};
  int dn_vtx [] = {6, 0};

  std::vector<PDM_g_num_t> dedge_vtx[] = {
    {
      1, 2,
      2, 3,
      3, 4,
      4, 5,
      5, 6
    },
    {}
  };

  std::vector<double> dvtx_coord[] = {
    {
      0, 0, 0,
      1, 0, 0,
      2, 0, 0,
      3, 0, 0,
      4, 0, 0,
      5, 0, 0
    },
    {}
  };

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_USER,
                                        0,
                                        0,
                                        dn_edge[test_rank],
                                        dn_vtx [test_rank],
                                        pdm_comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord[test_rank].data(),
                          PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                             dedge_vtx[test_rank].data(),
                             NULL,
                             PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(mpart,
                          0,
                          dmesh);


  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *edge_ln_to_gn = NULL;
  int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_EDGE,
                                               &edge_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);
  PDM_UNUSED(n_edge);

  // Free memory
  PDM_dmesh_free(dmesh);
  PDM_multipart_free(mpart);
}