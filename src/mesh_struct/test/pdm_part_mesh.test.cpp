
#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_part_mesh.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_priv.h"



MPI_TEST_CASE("[pdm_part_mesh] Constructor",1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        4,
                                                        4,
                                                        0,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build(dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);
  // free
  PDM_dcube_nodal_gen_free(dcube);

  int n_part = 1;
  PDM_multipart_t* mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                PDM_SPLIT_DUAL_WITH_HILBERT,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);
  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);
  PDM_DMesh_nodal_free(dmn);

  PDM_part_mesh_t* pm = NULL;
  PDM_multipart_get_part_mesh(mpart, 0, &pm, PDM_OWNERSHIP_USER);


  PDM_part_mesh_free(pm);

  PDM_multipart_free(mpart);

}
