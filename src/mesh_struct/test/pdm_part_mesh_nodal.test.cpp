
#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_generate_mesh.h"

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
MPI_TEST_CASE("[pdm_part_mesh_nodal_elmts] Constructor",1) {
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

  // const PDM_g_num_t n_vtx            = 9;
  // const PDM_g_num_t n_face           = 8;
  // const PDM_g_num_t n_ridge          = 8;
  const int         n_tri_section_1  = 8;
  const int         n_bar_section_1  = 8;

  int connec_tri_1[24] = {6, 8, 9,
                          9, 5, 6,
                          2, 8, 1,
                          9, 3, 4,
                          6, 7, 8,
                          9, 4, 5,
                          2, 9, 8,
                          9, 2, 3};

  int connec_bar_1[16] = {1, 2,
                          2, 3,
                          4, 5,
                          3, 4,
                          6, 7,
                          5, 6,
                          8, 1,
                          7, 8};

  int n_part = 1;
  // int n_group_elmt = 1;
  // int dgroup_elmt_idx[2] = {0, 8};
  // // PDM_g_num_t dgroup_elmt[8] = {9, 10, 11, 12, 13, 14, 15, 16};
  // PDM_g_num_t dgroup_elmt[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_part_mesh_nodal_t* pmn  = PDM_part_mesh_nodal_create(2, n_part, pdm_comm);

  PDM_part_mesh_nodal_elmts_t* pelmts_surf  = PDM_part_mesh_nodal_elmts_create(2, n_part, pdm_comm);
  PDM_part_mesh_nodal_elmts_t* pelmts_ridge = PDM_part_mesh_nodal_elmts_create(1, n_part, pdm_comm);

  int tri_section_1 = PDM_part_mesh_nodal_elmts_add(pelmts_surf , PDM_MESH_NODAL_TRIA3);
  int bar_section_1 = PDM_part_mesh_nodal_elmts_add(pelmts_ridge, PDM_MESH_NODAL_BAR2 );

  int i_part = 0;
  PDM_part_mesh_nodal_elmts_std_set(pelmts_surf,
                                    tri_section_1,
                                    i_part,
                                    n_tri_section_1,
                                    connec_tri_1,
                                    NULL,
                                    NULL,
                                    NULL,
                                    PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_elmts_std_set(pelmts_ridge,
                                    bar_section_1,
                                    i_part,
                                    n_bar_section_1,
                                    connec_bar_1,
                                    NULL,
                                    NULL,
                                    NULL,
                                    PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pelmts_surf );
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pelmts_ridge);

  // PDM_part_mesh_nodal_elmts_free(pelmts_surf);
  // PDM_part_mesh_nodal_elmts_free(pelmts_ridge);

  PDM_part_mesh_nodal_free(pmn);
}


MPI_TEST_CASE("[pdm_part_mesh_nodal] part_comm_graph from gnum", 2) {
  
  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
                                                                PDM_MESH_NODAL_TETRA4, 1, NULL,
                                                                0., 0., 0., // x/y/z min
                                                                1., 1., 1., // x/y/z length
                                                                3, 3, 3, // x/y/z n vertices
                                                                1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options

  PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_GEOMETRY_KIND_CORNER);
  std::vector<std::vector<int>> expected_graph = {{10, 1, 1, 1,
                                                   11, 1, 1, 2,
                                                   12, 1, 1, 3,
                                                   13, 1, 1, 4,
                                                   14, 1, 1, 5,
                                                   15, 1, 1, 6,
                                                   16, 1, 1, 7,
                                                   17, 1, 1, 8,
                                                   18, 1, 1, 9},
                                                  {1, 0, 1, 10,
                                                   2, 0, 1, 11,
                                                   3, 0, 1, 12,
                                                   4, 0, 1, 13,
                                                   5, 0, 1, 14,
                                                   6, 0, 1, 15,
                                                   7, 0, 1, 16,
                                                   8, 0, 1, 17,
                                                   9, 0, 1, 18,}};
  
  PDM_part_comm_graph_t *pcg_vtx = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn, PDM_GEOMETRY_KIND_CORNER, &pcg_vtx, PDM_OWNERSHIP_BAD_VALUE);

  int *computed_graph = NULL;
  int n_entity = PDM_part_comm_graph_entity_graph_get(pcg_vtx, 0, &computed_graph, PDM_OWNERSHIP_BAD_VALUE);
  assert(n_entity==9);
  for (int i_entity=0; i_entity<n_entity; ++i_entity) {
    assert(computed_graph[4*i_entity+0]==expected_graph[i_rank][4*i_entity+0]);
    assert(computed_graph[4*i_entity+1]==expected_graph[i_rank][4*i_entity+1]);
    assert(computed_graph[4*i_entity+2]==expected_graph[i_rank][4*i_entity+2]);
    assert(computed_graph[4*i_entity+3]==expected_graph[i_rank][4*i_entity+3]);
  }

  // > Ridge and surfacic pcg will be empty because no internal elements
  PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_GEOMETRY_KIND_RIDGE);
  PDM_part_comm_graph_t *pcg_ridge = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn, PDM_GEOMETRY_KIND_RIDGE, &pcg_ridge, PDM_OWNERSHIP_BAD_VALUE);
  n_entity = PDM_part_comm_graph_entity_graph_get(pcg_ridge, 0, &computed_graph, PDM_OWNERSHIP_BAD_VALUE);
  assert(n_entity==0);

  PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_GEOMETRY_KIND_SURFACIC);
  PDM_part_comm_graph_t *pcg_surfacic = NULL;
  PDM_part_mesh_nodal_part_comm_graph_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, &pcg_surfacic, PDM_OWNERSHIP_BAD_VALUE);
  n_entity = PDM_part_comm_graph_entity_graph_get(pcg_surfacic, 0, &computed_graph, PDM_OWNERSHIP_BAD_VALUE);
  assert(n_entity==0);


  PDM_part_mesh_nodal_free(pmn);
}
