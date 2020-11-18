#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_logging.h"

// n_vtx = 18
// double coord_x[n_vtx] = {0.000499536, 0.000499536, -1.65387e-06, -1.65387e-06, 0.000499536, 0.000499536, -1.65387e-06, -1.65387e-06, 0.00100073, 0.00100073, 0.00100073, 0.00100073, 0.000499536, -1.65387e-06, 0.000499536, -1.65387e-06, 0.00100073, 0.00100073};
// double coord_y[n_vtx] = {0.000498807, -2.38663e-06, -2.38663e-06, 0.000498807, 0.000501645, 4.5126e-07, 4.5126e-07, 0.000501645, 0.000498807, -2.38663e-06, 0.000501645, 4.5126e-07, 0.001, 0.001, 0.00100284, 0.00100284, 0.001, 0.00100284};
// double coord_z[n_vtx] = {0.000999549, 0.000999549, 0.000999549, 0.000999549, 0., 0., 0., 0., 0.000999549, 0.000999549, 0., 0., 0.000999549, 0.000999549, 0., 0., 0.000999549, 0.};

// Cas simple EMMA : /stck2/stck2.3/bmaugars/dev/dev-Tools/maia/unit_tests_case/EMMA/cube_simple/Cube_ANSA_hexa_separated.cgns

MPI_TEST_CASE("[PDM_dmesh_nodal_to_dmesh] decomposes hexa ",1) {

  const PDM_g_num_t n_vtx            = 18;
  const PDM_g_num_t n_cell           = 4;
  const int         n_hexa_section_1 = 4;
  const int         n_quad_section_1 = 16;
  PDM_g_num_t connec_hexa_1[32] = {1,2,3,4,5,6,7,8,
                                   9,10,2,1,11,12,6,5,
                                   13,1,4,14,15,5,8,16,
                                   17,9,1,13,18,11,5,15};

  PDM_g_num_t connec_quad_1[64] = {6, 5, 8, 7,
                                   12, 11, 5, 6,
                                   5, 15, 16, 8,
                                   11, 18, 15, 5,
                                   15, 13, 14, 16,
                                   13, 15, 18, 17,
                                   1, 2, 3, 4,
                                   9, 10, 2, 1,
                                   13, 1, 4, 14,
                                   17, 9, 1, 13,
                                   2, 6, 7, 3,
                                   6, 2, 10, 12,
                                   8, 4, 3, 7,
                                   4, 8, 16, 14,
                                   9, 11, 12, 10,
                                   11, 9, 17, 18};

  // PDM_g_num_t elemt_group[16] = {5, 6, 7, 8,      // Bottom_1
  //                                9, 10,           // Left_1
  //                                11, 12, 13, 14,  // Top_1
  //                                15, 16,          // Right_1
  //                                17, 18,          // Inlet_1
  //                                19, 20};         // Outlet_1

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, n_cell, -1, -1);

  // The order of call is important for global numbering
  int hexa_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_MESH_NODAL_HEXA8);
  int quad_section_1 = PDM_DMesh_nodal_section_add(dmn, PDM_MESH_NODAL_QUAD4);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  hexa_section_1,
                                  n_hexa_section_1,
                                  connec_hexa_1);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  quad_section_1,
                                  n_quad_section_1,
                                  connec_quad_1);

  PDM_dmesh_nodal_generate_distribution(dmn);

  PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmntodm, PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE);

  // Plante car on a pas implementer la decomposition en edge des Hexa
  // PDM_dmesh_nodal_to_dmesh_compute(dmntodm, PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE);


  PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  PDM_DMesh_nodal_free(dmn, 0);
}

  // int n_face_elt_tot     = -1;
  // int n_sum_vtx_face_tot = -1;
  // PDM_dmesh_nodal_decompose_faces_get_size(dmn, &n_face_elt_tot, &n_sum_vtx_face_tot);

  // printf("n_face_elt_tot     = %i\n", n_face_elt_tot);
  // printf("n_sum_vtx_face_tot = %i\n", n_sum_vtx_face_tot);

  // CHECK( n_face_elt_tot     == 12 );
  // CHECK( n_sum_vtx_face_tot == 48 );

  // // PDM_g_num_t* delmt_face_cell    = (PDM_g_num_t*) malloc(  n_face_elt_tot     * sizeof(PDM_g_num_t));
  // // int*         dcell_face_vtx_idx = (int        *) malloc( (n_face_elt_tot +1) * sizeof(int        ));
  // // PDM_g_num_t* dcell_face_vtx     = (PDM_g_num_t*) malloc(  n_sum_vtx_face_tot * sizeof(PDM_g_num_t));

  // std::vector<int>         dcell_face_vtx_idx(n_face_elt_tot +1);
  // std::vector<PDM_g_num_t> dcell_face_vtx(n_sum_vtx_face_tot);
  // std::vector<PDM_g_num_t> delmt_face_cell(n_face_elt_tot);

  // dcell_face_vtx_idx[0] = 0;
  // PDM_dmesh_nodal_decompose_faces(dmn,
  //                                 dcell_face_vtx_idx.data(),
  //                                 dcell_face_vtx.data(),
  //                                 delmt_face_cell.data(),
  //                                 NULL, NULL);
  // // dcell_face_vtx_idx[0] = 0;
  // // PDM_dmesh_nodal_decompose_faces(dmn,
  // //                                 dcell_face_vtx_idx,
  // //                                 dcell_face_vtx,
  // //                                 delmt_face_cell,
  // //                                 NULL, NULL);

  // // PDM_log_trace_array_long(delmt_face_cell, n_face_elt_tot, "delmt_face_cell:: ");
  // // PDM_log_trace_array_int(dcell_face_vtx_idx, n_face_elt_tot+1, "dcell_face_vtx_idx:: ");
  // // PDM_log_trace_array_long(dcell_face_vtx, n_sum_vtx_face_tot, "dcell_face_vtx:: ");

  // std::vector<PDM_g_num_t> delmt_face_cell_expected    = {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2};
  // std::vector<int>         dcell_face_vtx_idx_expected = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48};
  // std::vector<PDM_g_num_t> dcell_face_vtx_expected     = {4, 5, 2, 1, 11, 10, 7, 8,
  //                                                         7, 10, 4, 1, 10, 11, 5, 4,
  //                                                         5, 11, 8, 2, 2, 8, 7, 1,
  //                                                         5, 6, 3, 2, 12, 11, 8, 9,
  //                                                         8, 11, 5, 2, 11, 12, 6, 5,
  //                                                         6, 12, 9, 3, 3, 9, 8, 2};

  // CHECK( delmt_face_cell    == delmt_face_cell_expected);
  // CHECK( dcell_face_vtx_idx == dcell_face_vtx_idx_expected);
  // CHECK( dcell_face_vtx     == dcell_face_vtx_expected);

  // PDM_DMesh_nodal_cell_face_compute(dmn);

  // // free(delmt_face_cell);
  // // free(dcell_face_vtx_idx);
  // // free(dcell_face_vtx);
  // PDM_g_num_t* dface_cell;
  // PDM_g_num_t* dface_vtx;
  // int*         dface_vtx_idx;
  // int dn_face = PDM_DMesh_nodal_face_cell_get(dmn, &dface_cell);
  // PDM_DMesh_nodal_face_vtx_get(dmn, &dface_vtx_idx, &dface_vtx);

  // CHECK( dn_face == 11);

  // PDM_g_num_t dface_cell_expected[22]    = {1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 0, 1, 0, 2, 0 };
  // int         dface_vtx_idx_expected[12] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44};
  // PDM_g_num_t dface_vtx_expected[44]     = {4, 5, 2, 1, 5, 6, 3, 2, 2, 8, 7, 1, 3, 9, 8, 2, 7, 10, 4, 1, 8, 11,
  //                                           5, 2, 6, 12, 9, 3, 10, 11, 5, 4, 11, 12, 6, 5, 11, 10, 7, 8, 12, 11, 8, 9};

  // MPI_CHECK_EQ_C_ARRAY(0, dface_cell   , dface_cell_expected   , 2*dn_face             );
  // MPI_CHECK_EQ_C_ARRAY(0, dface_vtx_idx, dface_vtx_idx_expected, dn_face+1             );
  // MPI_CHECK_EQ_C_ARRAY(0, dface_vtx    , dface_vtx_expected    , dface_vtx_idx[dn_face]);

  // // printf(" dface_vtx_idx[dn_face]::%i\n",  dface_vtx_idx[dn_face]);
  // // PDM_log_trace_array_long(dface_cell, 2*dn_face, "dface_cell:: ");
  // // PDM_log_trace_array_int(dface_vtx_idx, dn_face+1, "dface_vtx_idx:: ");
  // // PDM_log_trace_array_long(dface_vtx, dface_vtx_idx[dn_face], "dface_vtx:: ");
