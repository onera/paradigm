#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_geom.h"
#include "pdm_priv.h"


MPI_TEST_CASE("[pdm_part_mesh_geom] 1d - BAR_2 - 2p", 2) {
  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  
  PDM_part_mesh_t* pm = PDM_part_mesh_create(1, pdm_comm);

  // Line of points (0, 0.1, 0.25, 0.4, 0.5, 0.55, 0.8, 1)
  double coords_p0[15]  = {0,0,0,  0.1,0,0,  0.25,0,0,  0.4,0,0,  0.5,0,0,};
  double coords_p1[12] = {0.5,0,0,  0.55,0,0,  0.8,0,0,  1,0,0};
  PDM_g_num_t vtx_gnum_p0[5] = {1,2,3,4,5};
  PDM_g_num_t vtx_gnum_p1[4] = {5,6,7,8};
  
  int elt_vtx_p0[8] = {1,2, 2,3, 3,4, 4,5};
  int elt_vtx_p1[6] = {1,2, 2,3, 3,4};
  
  int n_vtx;
  double *coords;
  int *elt_vtx;
  PDM_g_num_t* vtx_gnum;
  if (i_rank == 0) {
    n_vtx = 5;
    coords = coords_p0;
    elt_vtx = elt_vtx_p0;
    vtx_gnum = vtx_gnum_p0;
  }
  if (i_rank == 1) {
    n_vtx = 4;
    coords = coords_p1;
    elt_vtx = elt_vtx_p1;
    vtx_gnum = vtx_gnum_p1;
  }
  
  PDM_part_mesh_n_entity_set(pm, 0, PDM_MESH_ENTITY_VTX, n_vtx);
  PDM_part_mesh_n_entity_set(pm, 0, PDM_MESH_ENTITY_EDGE, n_vtx-1);
  
  PDM_part_mesh_vtx_coord_set(pm, 0, coords, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(pm, 0, PDM_MESH_ENTITY_VTX, vtx_gnum, PDM_OWNERSHIP_USER);
  
  PDM_part_mesh_connectivity_set(pm, 0, PDM_CONNECTIVITY_TYPE_EDGE_VTX, elt_vtx, NULL, PDM_OWNERSHIP_USER);

  double **pdual_volume = NULL;
  PDM_part_mesh_dual_volume_compute(pm, PDM_TRUE, &pdual_volume);

  double expected_dual_volume_p0[5] = {0.1/2, 0.25/2, 0.3/2, 0.25/2, 0.15/2};
  double expected_dual_volume_p1[4] = {0.15/2, 0.3/2, 0.45/2, 0.2/2};
  
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
  }
  
  PDM_free(pdual_volume[0]);
  PDM_free(pdual_volume);
  PDM_part_mesh_free(pm);
}
