#include <ext/alloc_traits.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_sphere_vol_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_doctest.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"

MPI_TEST_CASE("[pdm_sphere_vol_gen] - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n               = 4;
  PDM_g_num_t n_layer         = 3;
  double      x_center        = 0;
  double      y_center        = 0;
  double      z_center        = 0;
  double      radius_interior = 1;
  double      radius_exterior = 2;
  double      geometric_ratio = 1.;
  int         visu            = 0;

  /*
   *  Generate distributed Icoball
   */
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_vol_hollow_gen_nodal(pdm_comm,
                                  n,
                                  n_layer,
                                  x_center,
                                  y_center,
                                  z_center,
                                  radius_interior,
                                  radius_exterior,
                                  geometric_ratio,
                                  &dmn);

  if (visu) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "hollow_volume_");

    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "hollow_surface_");
  }
  PDM_DMesh_nodal_free(dmn);

}


MPI_TEST_CASE("[pdm_sphere_vol_gen] - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t n               = 4;
  PDM_g_num_t n_layer         = 3;
  double      x_center        = 0;
  double      y_center        = 0;
  double      z_center        = 0;
  double      radius_interior = 1;
  double      radius_exterior = 2;
  double      geometric_ratio = 1.;
  int         visu            = 0;

  /*
   *  Generate distributed Icoball
   */
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_vol_hollow_gen_nodal(pdm_comm,
                                  n,
                                  n_layer,
                                  x_center,
                                  y_center,
                                  z_center,
                                  radius_interior,
                                  radius_exterior,
                                  geometric_ratio,
                                  &dmn);

  if (visu) {
    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_VOLUMIC,
                             "hollow_volume_");

    PDM_dmesh_nodal_dump_vtk(dmn,
                             PDM_GEOMETRY_KIND_SURFACIC,
                             "hollow_surface_");
  }
  PDM_DMesh_nodal_free(dmn);

}
