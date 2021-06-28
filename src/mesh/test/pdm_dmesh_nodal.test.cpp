#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

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
MPI_TEST_CASE("[PDM_delmts_nodal_elmts_t] Constructor",1) {
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
  const PDM_g_num_t n_ridge          = 8;
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
  PDM_dmesh_nodal_t* dmn = PDM_DMesh_nodal_create(pdm_comm, 3, n_vtx, -1, n_face, -1);

  PDM_DMesh_nodal_elmts_t* dmn_elmts_surf  = PDM_DMesh_nodal_elmts_create(pdm_comm, 2, n_face );
  PDM_DMesh_nodal_elmts_t* dmn_elmts_ridge = PDM_DMesh_nodal_elmts_create(pdm_comm, 1, n_ridge);

  // The order of call is important for global numbering
  int tri_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn_elmts_surf , PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_DMesh_nodal_elmts_section_add(dmn_elmts_ridge, PDM_MESH_NODAL_BAR2);

  PDM_DMesh_nodal_elmts_section_std_set(dmn_elmts_surf,
                                        tri_section_1,
                                        n_tri_section_1,
                                        connec_tri_1,
                                        PDM_OWNERSHIP_USER);

  PDM_DMesh_nodal_elmts_section_std_set(dmn_elmts_ridge,
                                        bar_section_1,
                                        n_bar_section_1,
                                        connec_bar_1,
                                        PDM_OWNERSHIP_USER);

  PDM_Mesh_nodal_add_desh_nodal_elmts(dmn, dmn_elmts_surf );
  PDM_Mesh_nodal_add_desh_nodal_elmts(dmn, dmn_elmts_ridge);
  /*
   * Generate the connectivity
   */
  //
  PDM_dmesh_nodal_generate_distribution2(dmn);
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, pdm_comm, PDM_OWNERSHIP_KEEP);
  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmn);
  PDM_dmesh_nodal_to_dmesh_compute2(dmn_to_dm,
                                    PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                    PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

  PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  /*
   *  Flow du partitionnement :
   *     - elmt_vtx_vol + elmt_vtx_surf --> face_ln_to_gn + dface_gnum (dans la surface ) + dface_to_surf_gnum (dans le volume)
   *
   */

  PDM_DMesh_nodal_elmts_free(dmn_elmts_surf);
  PDM_DMesh_nodal_elmts_free(dmn_elmts_ridge);

  PDM_DMesh_nodal_free(dmn, 0);
}
