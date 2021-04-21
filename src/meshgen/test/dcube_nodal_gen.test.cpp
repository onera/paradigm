#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_dcube_nodal_gen.h"

MPI_TEST_CASE("[1p] dcube_nodal_gen HEXA ",1) {

  int n_vtx_seg = 3;

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_HEXA8,
                                                      PDM_OWNERSHIP_KEEP);

}
