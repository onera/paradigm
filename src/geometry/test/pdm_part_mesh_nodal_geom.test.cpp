#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_geom.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_logging.h"


MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - TRIA3 - simplices ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  int n_part = 1;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_TRIA3, 1, NULL,
                                                           0., 0., 0., // x/y/z min
                                                           1., 1.,     // x/y/z length
                                                           2, 2,       // x/y/z n vertices
                                                           1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options

  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, &pdual_volume);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - QUAD4 ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_QUAD4, 1, NULL,
                                                           0., 0., 0., // x/y/z min
                                                           1., 1.,     // x/y/z length
                                                           2, 2,       // x/y/z n vertices
                                                           1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options


  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, &pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 3d - TETRA4 - simplices ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
                                                                PDM_MESH_NODAL_TETRA4, 1, NULL,
                                                                0., 0., 0., // x/y/z min
                                                                1., 1., 1., // x/y/z length
                                                                3, 3, 3, // x/y/z n vertices
                                                                1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options



  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 3d - HEXA8 ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
                                                                PDM_MESH_NODAL_HEXA8, 1, NULL,
                                                                0., 0., 0., // x/y/z min
                                                                1., 1., 1., // x/y/z length
                                                                3, 3, 3, // x/y/z n vertices
                                                                1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options



  PDM_part_mesh_nodal_free(pmn);
}
