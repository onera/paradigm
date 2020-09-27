#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"

// double coord_x[n_vtx] = {0., 1., 2., 0., 1., 2., 0., 1., 2., 0., 1., 2};
// double coord_y[n_vtx] = {0., 0., 0., 1., 1., 1., 0., 0., 0., 1., 1., 1};
// double coord_z[n_vtx] = {0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1};

MPI_TEST_CASE("decomposes hexa ",1) {
  const PDM_g_num_t n_vtx            = 10;
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

  PDM_DMesh_nodal_free(dmesh_nodal_id);
}
