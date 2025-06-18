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


MPI_TEST_CASE("[pdm_multipart] - 2p - 3D, 1 domain, with dpart_id", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 3D, 1 domain, with dpart_id\n");

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
  int dn_cell = 4;
  int dn_face = 18;
  int dn_vtx[] = {14, 13};

  int dcell_face_idx[] = {0, 6, 12, 18, 24};
  std::vector<PDM_g_num_t> dcell_face[] = {
    {
      1, 5, 13, 17, 25, 29,
      -5, 9, 15, 19, 26, 30,
      2, 6, -17, 21, 27, 31,
      -6, 10, -19, 23, 28, 32
    },
    {
      3, 7, 14, 18, -29, 33,
      -7, 11, 16, 20, -30, 34,
      4, 8, -18, 22, -31, 35,
      -8, 12, -20, 24, -32, 36
    }
  };

  int dface_vtx_idx[] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72};
  std::vector<PDM_g_num_t> dface_vtx[] = {
    {
      10, 13, 4, 1,
      13, 16, 7, 4,
      19, 22, 13, 10,
      22, 25, 16, 13,
      2, 5, 14, 11,
      5, 8, 17, 14,
      11, 14, 23, 20,
      14, 17, 26, 23,
      3, 6, 15, 12,
      6, 9, 18, 15,
      12, 15, 24, 21,
      15, 18, 27, 24,
      2, 11, 10, 1,
      11, 20, 19, 10,
      3, 12, 11, 2,
      12, 21, 20, 11,
      4, 13, 14, 5,
      13, 22, 23, 14
    },
    {
      5, 14, 15, 6,
      14, 23, 24, 15,
      7, 16, 17, 8,
      16, 25, 26, 17,
      8, 17, 18, 9,
      17, 26, 27, 18,
      4, 5, 2, 1,
      5, 6, 3, 2,
      7, 8, 5, 4,
      8, 9, 6, 5,
      10, 11, 14, 13,
      11, 12, 15, 14,
      13, 14, 17, 16,
      14, 15, 18, 17,
      19, 20, 23, 22,
      20, 21, 24, 23,
      22, 23, 26, 25,
      23, 24, 27, 26
    }
  };

  std::vector<PDM_g_num_t> dface_cell[] = {
    {
      1, 0,
      3, 0,
      5, 0,
      7, 0,
      1, 2,
      3, 4,
      5, 6,
      7, 8,
      2, 0,
      4, 0,
      6, 0,
      8, 0,
      1, 0,
      5, 0,
      2, 0,
      6, 0,
      1, 3,
      5, 7
    },
    {
      2, 4,
      6, 8,
      3, 0,
      7, 0,
      4, 0,
      8, 0,
      1, 0,
      2, 0,
      3, 0,
      4, 0,
      1, 5,
      2, 6,
      3, 7,
      4, 8,
      5, 0,
      6, 0,
      7, 0,
      8, 0
    }
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
      2, 2, 0,
      0, 0, 1,
      1, 0, 1,
      2, 0, 1,
      0, 1, 1,
      1, 1, 1
    },
    {
      2, 1, 1,
      0, 2, 1,
      1, 2, 1,
      2, 2, 1,
      0, 0, 2,
      1, 0, 2,
      2, 0, 2,
      0, 1, 2,
      1, 1, 2,
      2, 1, 2,
      0, 2, 2,
      1, 2, 2,
      2, 2, 2
    }
  };


  int n_face_group = 0;
  int dface_group_idx[] = {0};

  SUBCASE("From dcell_face") {
    if (test_rank) printf("  From dcell_face\n");
    PDM_multipart_block_set(mpart,
                            0,
                            dn_cell,
                            dn_face,
                            dn_vtx[test_rank],
                            n_face_group,
                            dcell_face_idx,
                            dcell_face[test_rank].data(),
                            NULL,
                            dface_vtx_idx,
                            dface_vtx[test_rank].data(),
                            dvtx_coord[test_rank].data(),
                            dface_group_idx,
                            NULL);
  }

  SUBCASE("From dface_cell") {
    if (test_rank) printf("  From dface_cell\n");
    PDM_multipart_block_set(mpart,
                            0,
                            dn_cell,
                            dn_face,
                            dn_vtx[test_rank],
                            n_face_group,
                            NULL,
                            NULL,
                            dface_cell[test_rank].data(),
                            dface_vtx_idx,
                            dface_vtx[test_rank].data(),
                            dvtx_coord[test_rank].data(),
                            dface_group_idx,
                            NULL);
  }

  // Set desired partitioning
  std::vector<int> dpart_id[] = {
    {0, 1, 0, 1},
    {0, 1, 0, 1}
  };

  PDM_multipart_dpart_id_set(mpart,
                             0,
                             dpart_id[test_rank].data(),
                             PDM_OWNERSHIP_USER);

  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *cell_ln_to_gn = NULL;
  int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_CELL,
                                               &cell_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);

  int expected_n_cell = 4;
  CHECK(n_cell == expected_n_cell);

  // Check cell global IDs
  static PDM_g_num_t expected_cell_ln_to_gn0[] = {1, 3, 5, 7};
  static PDM_g_num_t expected_cell_ln_to_gn1[] = {2, 4, 6, 8};

  MPI_CHECK_EQ_C_ARRAY(0, cell_ln_to_gn, expected_cell_ln_to_gn0, n_cell);
  MPI_CHECK_EQ_C_ARRAY(1, cell_ln_to_gn, expected_cell_ln_to_gn1, n_cell);

  // Free memory
  PDM_multipart_free(mpart);
}



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



MPI_TEST_CASE("[pdm_multipart] - 2p - 3D, 1 domain, one empty block", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  if (test_rank == 0) printf("RUN [pdm_multipart] - 2p - 3D, 1 domain, one empty block\n");

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
  int dn_cell[] = {8,  0};
  int dn_face[] = {36, 0};
  int dn_vtx [] = {27, 0};

  std::vector<int> dcell_face_idx[] = {
    {0, 6, 12, 18, 24, 30, 36, 42, 48},
    {0}
  };

  std::vector<PDM_g_num_t> dcell_face[] = {
    {
      1, 5, 13, 17, 25, 29,
      -5, 9, 15, 19, 26, 30,
      2, 6, -17, 21, 27, 31,
      -6, 10, -19, 23, 28, 32,
      3, 7, 14, 18, -29, 33,
      -7, 11, 16, 20, -30, 34,
      4, 8, -18, 22, -31, 35,
      -8, 12, -20, 24, -32, 36
    },
    {}
  };

  std::vector<int> dface_vtx_idx[] = {
    {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120, 124, 128, 132, 136, 140, 144},
    {0}
  };

  std::vector<PDM_g_num_t> dface_vtx[] = {
    {
      10, 13, 4, 1,
      13, 16, 7, 4,
      19, 22, 13, 10,
      22, 25, 16, 13,
      2, 5, 14, 11,
      5, 8, 17, 14,
      11, 14, 23, 20,
      14, 17, 26, 23,
      3, 6, 15, 12,
      6, 9, 18, 15,
      12, 15, 24, 21,
      15, 18, 27, 24,
      2, 11, 10, 1,
      11, 20, 19, 10,
      3, 12, 11, 2,
      12, 21, 20, 11,
      4, 13, 14, 5,
      13, 22, 23, 14,
      5, 14, 15, 6,
      14, 23, 24, 15,
      7, 16, 17, 8,
      16, 25, 26, 17,
      8, 17, 18, 9,
      17, 26, 27, 18,
      4, 5, 2, 1,
      5, 6, 3, 2,
      7, 8, 5, 4,
      8, 9, 6, 5,
      10, 11, 14, 13,
      11, 12, 15, 14,
      13, 14, 17, 16,
      14, 15, 18, 17,
      19, 20, 23, 22,
      20, 21, 24, 23,
      22, 23, 26, 25,
      23, 24, 27, 26
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
      2, 2, 0,
      0, 0, 1,
      1, 0, 1,
      2, 0, 1,
      0, 1, 1,
      1, 1, 1,
      2, 1, 1,
      0, 2, 1,
      1, 2, 1,
      2, 2, 1,
      0, 0, 2,
      1, 0, 2,
      2, 0, 2,
      0, 1, 2,
      1, 1, 2,
      2, 1, 2,
      0, 2, 2,
      1, 2, 2,
      2, 2, 2
    },
    {}
  };

  PDM_dmesh_t *dmesh = PDM_dmesh_create(PDM_OWNERSHIP_USER,
                                        dn_cell[test_rank],
                                        dn_face[test_rank],
                                        0,
                                        dn_vtx [test_rank],
                                        pdm_comm);

  PDM_dmesh_vtx_coord_set(dmesh,
                          dvtx_coord[test_rank].data(),
                          PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_CELL_FACE,
                             dcell_face    [test_rank].data(),
                             dcell_face_idx[test_rank].data(),
                             PDM_OWNERSHIP_USER);

  PDM_dmesh_connectivity_set(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             dface_vtx    [test_rank].data(),
                             dface_vtx_idx[test_rank].data(),
                             PDM_OWNERSHIP_USER);

  PDM_multipart_dmesh_set(mpart,
                          0,
                          dmesh);
              
  // Compute partitioning
  PDM_multipart_compute(mpart);

  // Get partitioned mesh
  PDM_g_num_t *cell_ln_to_gn = NULL;
  int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                               0,
                                               0,
                                               PDM_MESH_ENTITY_CELL,
                                               &cell_ln_to_gn,
                                               PDM_OWNERSHIP_KEEP);
  PDM_UNUSED(n_cell);

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