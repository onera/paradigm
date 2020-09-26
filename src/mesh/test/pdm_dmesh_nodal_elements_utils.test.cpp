#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_elements_utils.h"


MPI_TEST_CASE("decomposes hexa ",1) {
  const PDM_g_num_t n_vtx            = 10;
  const PDM_g_num_t n_cell           = 10;
  const int         n_hexa_section_1 = 10;
  PDM_g_num_t connec_hexa_1[n_hexa_section_1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  int mesh_id = PDM_DMesh_nodal_create(test_comm, n_vtx, n_cell);

  int hexa_section_1 = PDM_DMesh_nodal_section_add(mesh_id, PDM_MESH_NODAL_HEXA8);

  // PDM_DMesh_nodal_section_std_set(mesh_id,
  //                                 hexa_section_1,
  //                                 n_hexa_section_1,
  //                                 connec_hexa_1);

  PDM_DMesh_nodal_free(mesh_id);
}
