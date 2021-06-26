#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"

/*
 *  Use case
 *
 *       3           4           2
 *       +-----------+-----------+
 *       |           |           |
 *       |           |           |
 *       |           |           |
 *       +-----------+-----------+
 *       5           1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */
MPI_TEST_CASE("[PDM_delmts_nodal_elmts_t] Constructor ",1) {
  // double dvtx_coord[27] = { 1. , 0. , 0.,
  //                           1. , 0.5, 0.,
  //                           1. , 1. , 0.,
  //                           1.5, 1. , 0.,
  //                           2. , 1. , 0.,
  //                           2. , 0.5, 0.,
  //                           2. , 0. , 0.,
  //                           1.5, 0. , 0.,
  //                           1.5, 0.5, 0.};
  // PDM_UNUSED(dvtx_coord);

  const PDM_g_num_t n_vtx            = 9;
  const PDM_g_num_t n_face           = 8;
  const int         n_tri_section_1  = 8;
  const int         n_bar_section_1  = 8;

  PDM_g_num_t connec_tri_1[24] = {6, 8, 9,
                                  9, 5, 6,
                                  2, 8, 1,
                                  9, 3, 4,
                                  6, 7, 8,
                                  9, 4, 5,
                                  2, 9, 8,
                                  9, 2, 3};

  PDM_g_num_t connec_bar_1[16] = {1, 2,
                                  2, 3,
                                  4, 5,
                                  3, 4,
                                  6, 7,
                                  5, 6,
                                  8, 1,
                                  7, 8};

  int n_group_elmt = 1;
  int dgroup_elmt_idx[2] = {0, 8};
  PDM_g_num_t dgroup_elmt[8] = {9, 10, 11, 12, 13, 14, 15, 16};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_DMesh_nodal_elmts_t* dmn_elmts = PDM_DMesh_nodal_elmts_create(pdm_comm, 2, n_face);

  // The order of call is important for global numbering
  int tri_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn_elmts, PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn_elmts, PDM_MESH_NODAL_BAR2);

  PDM_DMesh_nodal_elmts_section_std_set(dmn_elmts,
                                        tri_section_1,
                                        n_tri_section_1,
                                        connec_tri_1,
                                        PDM_OWNERSHIP_USER);
  PDM_DMesh_nodal_elmts_free(dmn_elmts);

  // PDM_DMesh_nodal_section_std_set(dmn,
  //                                 bar_section_1,
  //                                 n_bar_section_1,
  //                                 connec_bar_1,
  //                                 PDM_OWNERSHIP_USER);

  // PDM_DMesh_nodal_section_group_elmt_set(dmn, n_group_elmt, dgroup_elmt_idx, dgroup_elmt, PDM_OWNERSHIP_USER);

  // PDM_dmesh_nodal_generate_distribution(dmn);

  // PDM_dmesh_nodal_to_dmesh_t* dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);

  // PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

  // PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  // PDM_dmesh_nodal_to_dmesh_transform_to_coherent_dmesh(dmntodm, 2);

  // PDM_dmesh_t* dm;
  // PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dm);

  // // int         edge_face_idx_expected_p0[9] = {0, 1, 2, 3, 4, 5, 7, 9, 10 };
  // // PDM_g_num_t edge_face_expected_p0[32]    = {3, 0, 8, 0, 4, 0, 3, 0, 6, 0, 7, 3, 8, 7, 2, 0,
  // //                                             4, 8, 6, 4, 5, 0, 2, 6, 1, 5, 5, 0, 1, 2, 7, 1};

  // int dn_cell, dn_face, dn_vtx, dn_edge, n_bnd, n_join;
  // PDM_dmesh_dims_get(dm, &dn_cell, &dn_face, &dn_edge, &dn_vtx, &n_bnd, &n_join);

  // MPI_CHECK(0, dn_face == 8);

  // MPI_CHECK(0, dn_edge == 16);

  // PDM_g_num_t *edge_face;
  // int         *edge_face_idx;
  // PDM_dmesh_connectivity_get(dm, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
  //                            &edge_face, &edge_face_idx, PDM_OWNERSHIP_KEEP);

  // PDM_log_trace_array_long (edge_face, 2 * dn_edge, "edge_face:: ");
  // // MPI_CHECK_EQ_C_ARRAY(0, edge_face, edge_face_expected_p0, 2 * dn_edge);

  // PDM_dmesh_nodal_to_dmesh_free(dmntodm);
  // PDM_DMesh_nodal_free(dmn, 0);

}
