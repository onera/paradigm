#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_logging.h"

// double coord_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};

MPI_TEST_CASE("decomposes hexa ",1) {
  const PDM_g_num_t n_vtx            = 12;
  const PDM_g_num_t n_cell           = 2;
  const int         n_hexa_section_1 = 2;
  PDM_g_num_t connec_hexa_1[16] = {1, 2, 5, 4, 7, 8, 11, 10, // First
                                   2, 3, 6, 5, 8, 9, 12, 11};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int dmesh_nodal_id = PDM_DMesh_nodal_create(pdm_comm, n_vtx, n_cell);

  int hexa_section_1 = PDM_DMesh_nodal_section_add(dmesh_nodal_id, PDM_MESH_NODAL_HEXA8);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal_id,
                                  hexa_section_1,
                                  n_hexa_section_1,
                                  connec_hexa_1);
  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_decompose_faces_get_size(dmesh_nodal_id, &n_face_elt_tot, &n_sum_vtx_face_tot);

  printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 12 );
  CHECK( n_sum_vtx_face_tot == 48 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_dmesh_nodal_decompose_faces(dmesh_nodal_id,
                                  dcell_face_vtx_idx.data(),
                                  dcell_face_vtx.data(),
                                  delmt_face_cell.data(),
                                  NULL);
  // dcell_face_vtx_idx[0] = 0;
  // PDM_dmesh_nodal_decompose_faces(dmesh_nodal_id,
  //                                 dcell_face_vtx_idx,
  //                                 dcell_face_vtx,
  //                                 delmt_face_cell,
  //                                 NULL);

  // PDM_log_trace_array_long(delmt_face_cell, n_face_elt_tot, "delmt_face_cell:: ");
  // PDM_log_trace_array_int(dcell_face_vtx_idx, n_face_elt_tot+1, "dcell_face_vtx_idx:: ");
  // PDM_log_trace_array_long(dcell_face_vtx, n_sum_vtx_face_tot, "dcell_face_vtx:: ");

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {4, 5, 2, 1, 11, 10, 7, 8,
                                                          7, 10, 4, 1, 10, 11, 5, 4,
                                                          5, 11, 8, 2, 2, 8, 7, 1,
                                                          5, 6, 3, 2, 12, 11, 8, 9,
                                                          8, 11, 5, 2, 11, 12, 6, 5,
                                                          6, 12, 9, 3, 3, 9, 8, 2};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);

  PDM_DMesh_nodal_cell_face_compute(dmesh_nodal_id);

  // free(delmt_face_cell);
  // free(dcell_face_vtx_idx);
  // free(dcell_face_vtx);
  PDM_g_num_t* dface_cell;
  PDM_g_num_t* dface_vtx;
  int*         dface_vtx_idx;
  int dn_face = PDM_DMesh_nodal_face_cell_get(dmesh_nodal_id, &dface_cell);
  PDM_DMesh_nodal_face_vtx_get(dmesh_nodal_id, &dface_vtx_idx, &dface_vtx);

  CHECK( dn_face == 11);

  PDM_g_num_t dface_cell_expected[22]    = {1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0 };
  int         dface_vtx_idx_expected[12] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44};
  PDM_g_num_t dface_vtx_expected[44]     = {4, 5, 2, 1, 5, 6, 3, 2, 2, 8, 7, 1, 3, 9, 8, 2, 7, 10, 4, 1, 8, 11,
                                            5, 2, 6, 12, 9, 3, 10, 11, 5, 4, 11, 12, 6, 5, 11, 10, 7, 8, 12, 11, 8, 9};

  MPI_CHECK_EQ_C_ARRAY(0, dface_cell   , dface_cell_expected   , 2*dn_face             );
  MPI_CHECK_EQ_C_ARRAY(0, dface_vtx_idx, dface_vtx_idx_expected, dn_face+1             );
  MPI_CHECK_EQ_C_ARRAY(0, dface_vtx    , dface_vtx_expected    , dface_vtx_idx[dn_face]);

  // printf(" dface_vtx_idx[dn_face]::%i\n",  dface_vtx_idx[dn_face]);
  // PDM_log_trace_array_long(dface_cell, 2*dn_face, "dface_cell:: ");
  // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1, "dface_vtx_idx:: ");
  // PDM_log_trace_array_long(dface_vtx, dface_vtx_idx[dn_face], "dface_vtx:: ");

  PDM_DMesh_nodal_free(dmesh_nodal_id);
}


// double coord_tetra_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_tetra_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_tetra_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};

MPI_TEST_CASE("decomposes tetra ",1) {
  const PDM_g_num_t n_vtx             = 12;
  const PDM_g_num_t n_cell            = 10;
  const int         n_tetra_section_1 = 10;
  PDM_g_num_t connec_tetra_1[40] = {1, 2, 4, 7,
                                    2, 5, 4, 11,
                                    4, 7, 11, 10,
                                    2, 7, 8, 11,
                                    2, 4, 7, 11,
                                    2, 6, 5, 11,
                                    2, 9, 3, 6,
                                    11, 12, 9, 6,
                                    9, 11, 2, 8,
                                    9, 2, 11, 6};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int dmesh_nodal_id = PDM_DMesh_nodal_create(pdm_comm, n_vtx, n_cell);

  int tetra_section_1 = PDM_DMesh_nodal_section_add(dmesh_nodal_id, PDM_MESH_NODAL_TETRA4);

  PDM_DMesh_nodal_section_std_set(dmesh_nodal_id,
                                  tetra_section_1,
                                  n_tetra_section_1,
                                  connec_tetra_1);
  int n_face_elt_tot     = -1;
  int n_sum_vtx_face_tot = -1;
  PDM_dmesh_nodal_decompose_faces_get_size(dmesh_nodal_id, &n_face_elt_tot, &n_sum_vtx_face_tot);

  printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  CHECK( n_face_elt_tot     == 40  );
  CHECK( n_sum_vtx_face_tot == 120 );

  // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);

  dcell_face_vtx_idx[0] = 0;
  PDM_dmesh_nodal_decompose_faces(dmesh_nodal_id,
                                  dcell_face_vtx_idx.data(),
                                  dcell_face_vtx.data(),
                                  delmt_face_cell.data(),
                                  NULL);
  // dcell_face_vtx_idx[0] = 0;
  // PDM_dmesh_nodal_decompose_faces(dmesh_nodal_id,
  //                                 dcell_face_vtx_idx,
  //                                 dcell_face_vtx,
  //                                 delmt_face_cell,
  //                                 NULL);

  PDM_log_trace_array_long(delmt_face_cell.data(), n_face_elt_tot, "delmt_face_cell:: ");
  PDM_log_trace_array_int(dcell_face_vtx_idx.data(), n_face_elt_tot+1, "dcell_face_vtx_idx:: ");
  PDM_log_trace_array_long(dcell_face_vtx.data(), n_sum_vtx_face_tot, "dcell_face_vtx:: ");

  std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10};
  std::vector<int>         dcell_face_vtx_idx_expected = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120};
  std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {1, 2, 4, 1, 7, 2, 1, 4, 7, 2, 7, 4, 2, 5, 4, 2, 11, 5, 2, 4, 11, 5, 11, 4, 4, 7, 11, 4, 10, 7, 4, 11, 10, 7, 10, 11, 2, 7, 8, 2, 11, 7, 2, 8, 11, 7, 11, 8, 2, 4, 7, 2, 11, 4, 2, 7, 11, 4, 11, 7, 2, 6, 5, 2, 11, 6, 2, 5, 11, 6, 11, 5, 2, 9, 3, 2, 6, 9, 2, 3, 6, 9, 6, 3, 11, 12, 9, 11, 6, 12, 11, 9, 6, 12, 6, 9, 9, 11, 2, 9, 8, 11, 9, 2, 8, 11, 8, 2, 9, 2, 11, 9, 6, 2, 9, 11, 6, 2, 6, 11};

  CHECK( delmt_face_cell    == delmt_face_cell_expected);
  CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  CHECK( dcell_face_vtx     == dcell_face_vtx_expected);

  PDM_DMesh_nodal_cell_face_compute(dmesh_nodal_id);

  // free(delmt_face_cell);
  // free(dcell_face_vtx_idx);
  // free(dcell_face_vtx);
  PDM_g_num_t* dface_cell;
  PDM_g_num_t* dface_vtx;
  int*         dface_vtx_idx;
  int dn_face = PDM_DMesh_nodal_face_cell_get(dmesh_nodal_id, &dface_cell);
  PDM_DMesh_nodal_face_vtx_get(dmesh_nodal_id, &dface_vtx_idx, &dface_vtx);

  CHECK( dn_face == 30);

  PDM_g_num_t dface_cell_expected[60]   = {1, 0, 1, 0, 2, 0, 7, 0, 1, 0, 5, 1, 6, 0, 7, 0, 10, 7, 5, 2, 4, 0, 6, 2, 7, 0, 9, 0, 6, 10, 5, 4, 2, 0, 4, 9, 3, 0, 9, 10, 3, 5, 6, 0, 3, 0, 8, 10, 4, 0, 8, 0, 3, 0, 9, 0, 8, 0, 8, 0};
  int         dface_vtx_idx_expected[61] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90};
  PDM_g_num_t dface_vtx_expected[90]     = {1, 2, 4, 1, 7, 2, 2, 5, 4, 2, 3, 6, 1, 4, 7, 2, 4, 7, 2, 6, 5, 2, 9, 3, 9, 6, 2, 2, 11, 4, 2, 7, 8, 2, 5, 11, 9, 6, 3, 9, 2, 8, 2, 11, 6, 2, 7, 11, 5, 11, 4, 2, 8, 11, 4, 10, 7, 9, 11, 2, 4, 7, 11, 6, 11, 5, 4, 11, 10, 11, 9, 6, 7, 11, 8, 12, 6, 9, 7, 10, 11, 9, 8, 11, 11, 6, 12, 11, 12, 9};

  MPI_CHECK_EQ_C_ARRAY(0, dface_cell   , dface_cell_expected   , 2*dn_face             );
  MPI_CHECK_EQ_C_ARRAY(0, dface_vtx_idx, dface_vtx_idx_expected, dn_face+1             );
  MPI_CHECK_EQ_C_ARRAY(0, dface_vtx    , dface_vtx_expected    , dface_vtx_idx[dn_face]);

  // printf(" dface_vtx_idx[dn_face]::%i\n",  dface_vtx_idx[dn_face]);
  // PDM_log_trace_array_long(dface_cell, 2*dn_face, "dface_cell:: ");
  // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1, "dface_vtx_idx:: ");
  // PDM_log_trace_array_long(dface_vtx, dface_vtx_idx[dn_face], "dface_vtx:: ");

  PDM_DMesh_nodal_free(dmesh_nodal_id);
}
