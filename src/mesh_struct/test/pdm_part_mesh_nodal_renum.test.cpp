
#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_renum.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_logging.h"



MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_renum_cuthill", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  int n_vtx_seg = 4;
  int n_part    = 1;
  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_TRIA3,
                                                           1,
                                                           NULL,
                                                           -1,
                                                           -1,
                                                           0,
                                                           2,
                                                           2,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_part,
                                                           PDM_SPLIT_DUAL_WITH_PTSCOTCH);

  int visu = 0;
  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_in");
  }

  PDM_part_mesh_nodal_renum_cuthill(pmn);

  /* Some check */
  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);

  int *vtx_new_to_old = NULL;
  PDM_part_mesh_nodal_vtx_new_to_old_get(pmn,
                                         0,
                                         &vtx_new_to_old,
                                         PDM_OWNERSHIP_BAD_VALUE);

  int expected_vtx_new_to_old[16] = {0, 1, 4, 5, 8, 2, 9, 12, 6, 3, 13, 10, 7, 14, 11, 15};

  CHECK_EQ_C_ARRAY(vtx_new_to_old, expected_vtx_new_to_old, n_vtx);

  /* Some check */
  int n_surf_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, 0);

  int *surf_new_to_old = NULL;
  PDM_part_mesh_nodal_elmts_new_to_old_get(pmn,
                                           PDM_GEOMETRY_KIND_SURFACIC,
                                           0,
                                           &surf_new_to_old,
                                           PDM_OWNERSHIP_BAD_VALUE);

  int n_ridge_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_RIDGE, 0);

  int *ridge_new_to_old = NULL;
  PDM_part_mesh_nodal_elmts_new_to_old_get(pmn,
                                           PDM_GEOMETRY_KIND_RIDGE,
                                           0,
                                           &ridge_new_to_old,
                                           PDM_OWNERSHIP_BAD_VALUE);

  int n_corner_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_CORNER, 0);

  int *corner_new_to_old = NULL;
  PDM_part_mesh_nodal_elmts_new_to_old_get(pmn,
                                           PDM_GEOMETRY_KIND_CORNER,
                                           0,
                                           &corner_new_to_old,
                                           PDM_OWNERSHIP_BAD_VALUE);

  if(0 == 1) {
    PDM_log_trace_array_int(vtx_new_to_old   , n_vtx         , "vtx_new_to_old    :: ");
    PDM_log_trace_array_int(surf_new_to_old  , n_surf_elmts  , "surf_new_to_old   :: ");
    PDM_log_trace_array_int(ridge_new_to_old , n_ridge_elmts , "ridge_new_to_old  :: ");
    PDM_log_trace_array_int(corner_new_to_old, n_corner_elmts, "corner_new_to_old :: ");
  }

  int expected_surf_new_to_old  [18] = {0, 6, 2, 1, 12, 7, 8, 4, 3, 13, 9, 14, 10, 5, 15, 16, 11, 17};
  int expected_ridge_new_to_old [12] = {0, 6, 1, 7, 8, 2, 3, 9, 4, 10, 5, 11};
  int expected_corner_new_to_old[4 ] = {0, 2, 1, 3};

  CHECK_EQ_C_ARRAY(surf_new_to_old  , expected_surf_new_to_old  , n_surf_elmts  );
  CHECK_EQ_C_ARRAY(ridge_new_to_old , expected_ridge_new_to_old , n_ridge_elmts );
  CHECK_EQ_C_ARRAY(corner_new_to_old, expected_corner_new_to_old, n_corner_elmts);

  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_out");
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_RIDGE,
                                 "pmn_ridge_out");
  }

  PDM_part_mesh_nodal_free(pmn);
}



MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_renum_cuthill - mixed ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_nodal_from_file(pdm_comm,
                                                                 1,
                                                                 PDM_SPLIT_DUAL_WITH_HILBERT,
                                                                 PDM_MESH_DIR"mixed_elements_2d.mesh");

  int visu = 0;
  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_in");
  }

  PDM_part_mesh_nodal_renum_cuthill(pmn);

  /* Some check */
  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);

  int *vtx_new_to_old = NULL;
  PDM_part_mesh_nodal_vtx_new_to_old_get(pmn,
                                         0,
                                         &vtx_new_to_old,
                                         PDM_OWNERSHIP_BAD_VALUE);

  /* Some check */
  int n_surf_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, 0);

  int *surf_new_to_old = NULL;
  PDM_part_mesh_nodal_elmts_new_to_old_get(pmn,
                                           PDM_GEOMETRY_KIND_SURFACIC,
                                           0,
                                           &surf_new_to_old,
                                           PDM_OWNERSHIP_BAD_VALUE);

  int n_ridge_elmts = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_RIDGE, 0);

  int *ridge_new_to_old = NULL;
  PDM_part_mesh_nodal_elmts_new_to_old_get(pmn,
                                           PDM_GEOMETRY_KIND_RIDGE,
                                           0,
                                           &ridge_new_to_old,
                                           PDM_OWNERSHIP_BAD_VALUE);

  if(0 == 1) {
    PDM_log_trace_array_int(vtx_new_to_old   , n_vtx         , "vtx_new_to_old    :: ");
    PDM_log_trace_array_int(surf_new_to_old  , n_surf_elmts  , "surf_new_to_old   :: ");
    PDM_log_trace_array_int(ridge_new_to_old , n_ridge_elmts , "ridge_new_to_old  :: ");
  }

  int expected_vtx_new_to_old  [121] = {6, 36, 0, 41, 7, 35, 45, 8, 42, 34, 49, 9, 46, 43, 33, 53, 10, 50, 47, 44, 4, 57, 11, 61, 12, 65, 13, 69, 14, 73, 1, 37, 54, 51, 48, 32, 85, 84, 86, 83, 87, 82, 88, 81, 23, 2, 58, 62, 66, 70, 74, 38, 55, 52, 31, 89, 90, 91, 92, 22, 59, 63, 67, 71, 75, 39, 56, 30, 93, 94, 95, 96, 21, 60, 64, 68, 72, 76, 40, 29, 97, 98, 99, 100, 20, 28, 27, 26, 25, 24, 5, 101, 102, 103, 104, 19, 105, 109, 113, 77, 117, 106, 107, 108, 18, 110, 114, 78, 118, 111, 112, 17, 115, 79, 119, 116, 16, 120, 80, 15, 3};
  int expected_surf_new_to_old [150] = { 42, 92, 38, 43, 88, 34, 93, 39, 44, 84, 30, 89, 35, 94, 40, 45, 80, 26, 85, 31, 90, 36, 95, 41, 1, 76, 22, 72, 18, 68, 14, 64, 10, 60, 96, 46, 81, 27, 86, 32, 91, 37, 51, 2, 142, 143, 144, 145, 101, 77, 23, 73, 19, 69, 15, 65, 11, 61, 97, 47, 82, 28, 87, 33, 52, 3, 138, 139, 140, 141, 102, 78, 24, 74, 20, 70, 16, 66, 12, 62, 98, 48, 83, 29, 53, 4, 134, 135, 136, 137, 103, 79, 25, 75, 21, 71, 17, 67, 13, 63, 99, 49, 54, 5, 130, 131, 132, 133, 104, 55, 6, 56, 7, 57, 8, 9, 59, 58, 50, 0, 126, 127, 128, 129, 105, 122, 118, 114, 146, 110, 123, 124, 125, 106, 119, 115, 147, 111, 120, 121, 107, 116, 148, 112, 117, 108, 113, 149, 109, 100};
  int expected_ridge_new_to_old[50]  = {30, 28, 29, 31, 27, 32, 26, 33, 25, 20, 24, 48, 12, 23, 22, 21, 34, 36, 37, 47, 13, 14, 15, 3, 2, 4, 38, 46, 5, 39, 45, 6, 35, 44, 7, 43, 42, 41, 40, 49, 16, 8, 17, 9, 18, 10, 19, 11, 0, 1};

  CHECK_EQ_C_ARRAY(vtx_new_to_old   , expected_vtx_new_to_old   , n_vtx         );
  CHECK_EQ_C_ARRAY(surf_new_to_old  , expected_surf_new_to_old  , n_surf_elmts  );
  CHECK_EQ_C_ARRAY(ridge_new_to_old , expected_ridge_new_to_old , n_ridge_elmts );

  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_out");
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_RIDGE,
                                 "pmn_ridge_out");
  }

  PDM_part_mesh_nodal_free(pmn);

}


MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_reorder_vtx", 2) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  int n_vtx_seg = 4;
  int n_part    = 1;
  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_TRIA3,
                                                           1,
                                                           NULL,
                                                           -1,
                                                           -1,
                                                           0,
                                                           2,
                                                           2,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_part,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT);

  int visu = 0;
  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_in");
  }

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  PDM_g_num_t *vtx_gnum = PDM_part_mesh_nodal_vtx_g_num_get(pmn, 0, PDM_OWNERSHIP_BAD_VALUE);

  std::vector<std::vector<int>> new_to_old = {{5, 6, 7, 8, 9, 4, 1, 2, 3, 0, 10},
                                              {1, 0, 5, 3, 4, 2, 6, 7, 9, 8, 10}};
  int *perm[1] = {new_to_old[i_rank].data()};
  PDM_part_mesh_nodal_reorder_vtx(pmn, perm);

  /* Some check */
  if(0 == 1) {
    PDM_log_trace_array_long(vtx_gnum, n_vtx, "vtx_gnum ::");
  }

  PDM_g_num_t expected_vtx_gnum_p0[11] = {6, 7, 8, 9, 10, 5, 2, 3, 4, 1, 11};
  PDM_g_num_t expected_vtx_gnum_p1[11] = {7, 6, 11, 9, 10, 8, 12, 13, 15, 14, 16};

  MPI_CHECK_EQ_C_ARRAY(0, vtx_gnum, expected_vtx_gnum_p0, n_vtx);
  MPI_CHECK_EQ_C_ARRAY(1, vtx_gnum, expected_vtx_gnum_p1, n_vtx);

  PDM_part_comm_graph_t* pcg_vtx = NULL;
  PDM_part_mesh_nodal_part_comm_graph_vtx_get(pmn,
                                              &pcg_vtx,
                                              PDM_OWNERSHIP_BAD_VALUE);

  int *vtx_graph = NULL;
  int n_vtx_graph = PDM_part_comm_graph_entity_graph_get(pcg_vtx,
                                                         0,
                                                         &vtx_graph,
                                                         PDM_OWNERSHIP_BAD_VALUE);

  if(0 == 1) {
    PDM_log_trace_array_int(vtx_graph, 4 * n_vtx_graph, "vtx_graph ::");
  }

  int expected_vtx_graph_p0[24] = {1, 1, 1, 2, 2, 1, 1, 1, 3, 1, 1, 6, 4, 1, 1, 4, 5, 1, 1, 5, 11, 1, 1, 3};
  int expected_vtx_graph_p1[24] = {2, 0, 1, 1, 1, 0, 1, 2, 6, 0, 1, 3, 4, 0, 1, 4, 5, 0, 1, 5, 3, 0, 1, 11};

  MPI_CHECK_EQ_C_ARRAY(0, vtx_graph, expected_vtx_graph_p0, 4 * n_vtx_graph);
  MPI_CHECK_EQ_C_ARRAY(1, vtx_graph, expected_vtx_graph_p1, 4 * n_vtx_graph);

  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_out");
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_RIDGE,
                                 "pmn_ridge_out");
  }

  PDM_part_mesh_nodal_free(pmn);
}



MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_renum_cuthill - mixed - 3d ", 2) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_nodal_from_file(pdm_comm,
                                                                 1,
                                                                 PDM_SPLIT_DUAL_WITH_HILBERT,
                                                                 PDM_MESH_DIR"mixed_elements_3d.mesh");

  int visu = 0;
  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_VOLUMIC,
                                 "pmn_in");
  }

  PDM_part_mesh_nodal_renum_cuthill(pmn);

  if(visu == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_SURFACIC,
                                 "pmn_out");
    PDM_part_mesh_nodal_dump_vtk(pmn,
                                 PDM_GEOMETRY_KIND_RIDGE,
                                 "pmn_ridge_out");
  }

  PDM_part_mesh_nodal_free(pmn);

}
