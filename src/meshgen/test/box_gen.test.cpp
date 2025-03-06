#include <ext/alloc_traits.h>
#include <vector>
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_box_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_doctest.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"


MPI_TEST_CASE("[pdm_box_gen] - cartesian - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  double      radius = 10.;

  int          n_box        = 0;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;
  PDM_box_gen_cartesian(pdm_comm,
                        5,                         // nx
                        5,                         // ny
                        5,                         // nz
                        -radius, -radius, -radius, // x,y,z_min
                        radius, radius, radius,    // x,y,z_max
                        &n_box,
                        &box_extents,
                        &box_ln_to_gn);

  CHECK(n_box == 64);

  free(box_extents);
  free(box_ln_to_gn);

}

MPI_TEST_CASE("[pdm_box_gen] - cartesian - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  double      radius = 10.;

  int          n_box        = 0;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;
  PDM_box_gen_cartesian(pdm_comm,
                        5,                         // nx
                        5,                         // ny
                        5,                         // nz
                        -radius, -radius, -radius, // x,y,z_min
                        radius, radius, radius,    // x,y,z_max
                        &n_box,
                        &box_extents,
                        &box_ln_to_gn);

  MPI_CHECK(0, n_box == 32);
  MPI_CHECK(1, n_box == 32);

  free(box_extents);
  free(box_ln_to_gn);

}




MPI_TEST_CASE("[pdm_box_gen] - random - 1p", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t gn_box = 10;
  double      radius = 10.;

  int          n_box        = 0;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;
  PDM_box_gen_random(pdm_comm,
                     0,                         // seed
                     0,                         // geometric_g_num
                     gn_box,                    // gn_box
                     0.5*radius,                // min_size
                     1.5*radius,                // max_size
                     -radius, -radius, -radius, // x,y,z_min
                     radius, radius, radius,    // x,y,z_max
                     &n_box,
                     &box_extents,
                     &box_ln_to_gn);

  CHECK(n_box == 10);

  free(box_extents);
  free(box_ln_to_gn);

}

MPI_TEST_CASE("[pdm_box_gen] - random - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  PDM_g_num_t gn_box = 10;
  double      radius = 10.;

  int          n_box        = 0;
  double      *box_extents  = NULL;
  PDM_g_num_t *box_ln_to_gn = NULL;
  PDM_box_gen_random(pdm_comm,
                     0,                         // seed
                     0,                         // geometric_g_num
                     gn_box,                    // gn_box
                     0.5*radius,                // min_size
                     1.5*radius,                // max_size
                     -radius, -radius, -radius, // x,y,z_min
                     radius, radius, radius,    // x,y,z_max
                     &n_box,
                     &box_extents,
                     &box_ln_to_gn);

  MPI_CHECK(0, n_box == 5);
  MPI_CHECK(1, n_box == 5);

  free(box_extents);
  free(box_ln_to_gn);

}
