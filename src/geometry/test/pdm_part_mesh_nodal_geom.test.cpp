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


MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - BAR_2 - 3p", 3) {
  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(1, 1, pdm_comm);

  // Line of points (0, 0.1, 0.25, 0.4, 0.5, 0.55, 0.8, 1)
  double coords_p0[9]  = {0,0,0,  0.1,0,0,  0.25,0,0};
  double coords_p1[12] = {0.25,0,0,  0.4,0,0,  0.5,0,0, 0.55,0,0};
  double coords_p2[9]  = {0.55,0,0,  0.8,0,0,  1,0,0};
  PDM_g_num_t vtx_gnum_p0[3] = {1,2,3};
  PDM_g_num_t vtx_gnum_p1[4] = {3,4,5,6};
  PDM_g_num_t vtx_gnum_p2[3] = {6,7,8};
  
  int elt_vtx_p0[4] = {1,2, 2,3};
  int elt_vtx_p1[6] = {1,2, 2,3, 3,4};
  int elt_vtx_p2[4] = {1,2, 2,3};
  PDM_g_num_t bar_gnum_p0[2] = {1,2};
  PDM_g_num_t bar_gnum_p1[3] = {3,4,5};
  PDM_g_num_t bar_gnum_p2[2] = {6,7};
  
  int n_vtx;
  double *coords;
  int *elt_vtx;
  PDM_g_num_t* vtx_gnum;
  PDM_g_num_t* bar_gnum;
  if (i_rank == 0) {
    n_vtx = 3;
    coords = coords_p0;
    elt_vtx = elt_vtx_p0;
    vtx_gnum = vtx_gnum_p0;
    bar_gnum = bar_gnum_p0;
  }
  if (i_rank == 1) {
    n_vtx = 4;
    coords = coords_p1;
    elt_vtx = elt_vtx_p1;
    vtx_gnum = vtx_gnum_p1;
    bar_gnum = bar_gnum_p1;
  }
  if (i_rank == 2) {
    n_vtx = 3;
    coords = coords_p2;
    elt_vtx = elt_vtx_p2;
    vtx_gnum = vtx_gnum_p2;
    bar_gnum = bar_gnum_p2;
  }
  
  PDM_part_mesh_nodal_coord_set(pmn, 0, n_vtx, coords, PDM_OWNERSHIP_USER);
  PDM_part_mesh_nodal_vtx_gnum_set(pmn, 0, vtx_gnum, PDM_OWNERSHIP_USER);
  
  int id_bar = PDM_part_mesh_nodal_section_add(pmn, PDM_MESH_NODAL_BAR2); 
  PDM_part_mesh_nodal_section_std_set(pmn, id_bar, 0, n_vtx-1, elt_vtx, bar_gnum,
                                      NULL, NULL, PDM_OWNERSHIP_USER);

  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);
  
  double expected_dual_volume_p0[3] = {0.1/2, 0.25/2, 0.3/2};
  double expected_dual_volume_p1[4] = {0.3/2, 0.25/2, 0.15/2, 0.3/2};
  double expected_dual_volume_p2[3] = {0.3/2, 0.45/2, 0.2/2};
  
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
    MPI_CHECK(2, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p2[i_vtx]).epsilon(0.01));
  }
  
  PDM_free(pdual_volume[0]);
  PDM_free(pdual_volume);
  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - TRIA3 - simplices - 1p", 1) {

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
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume[4] = {1.6666666666666666e-01, 3.3333333333333331e-01, 3.3333333333333331e-01, 1.6666666666666666e-01};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    CHECK(pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume[i_vtx]).epsilon(0.01));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}


MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - TRIA3 - simplices - 2p", 2) {

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

  // PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "surf_");

  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume_p0[3] = {1.6666666666666666e-01, 3.3333333333333331e-01, 3.3333333333333331e-01};
  double expected_dual_volume_p1[3] = {3.3333333333333331e-01, 3.3333333333333331e-01, 1.6666666666666666e-01};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
  }

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
  int n_part = 1;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_QUAD4, 1, NULL,
                                                           0., 0., 0., // x/y/z min
                                                           1., 1.,     // x/y/z length
                                                           3, 3,       // x/y/z n vertices
                                                           1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options


  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "surf_");
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume_p0[9] = {6.2500000000000000e-02, 1.2500000000000000e-01, 6.2500000000000000e-02,
                                       1.2500000000000000e-01, 2.5000000000000000e-01, 1.2500000000000000e-01,
                                       6.2500000000000000e-02, 1.2500000000000000e-01, 6.2500000000000000e-02};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}


MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 2d - QUAD4 - 2p", 2) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  int n_part = 1;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_QUAD4, 1, NULL,
                                                           0., 0., 0., // x/y/z min
                                                           1., 1.,     // x/y/z length
                                                           3, 3,       // x/y/z n vertices
                                                           1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options


  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume_p0[6] = {6.2500000000000000e-02, 1.2500000000000000e-01, 6.2500000000000000e-02,
                                       1.2500000000000000e-01, 2.5000000000000000e-01, 1.2500000000000000e-01};
  double expected_dual_volume_p1[6] = {1.2500000000000000e-01, 2.5000000000000000e-01, 1.2500000000000000e-01,
                                       6.2500000000000000e-02, 1.2500000000000000e-01, 6.2500000000000000e-02};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 3d - TETRA4 - simplices ", 2) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  int n_part = 1;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
                                                                PDM_MESH_NODAL_TETRA4, 1, NULL,
                                                                0., 0., 0., // x/y/z min
                                                                1., 1., 1., // x/y/z length
                                                                3, 2, 2, // x/y/z n vertices
                                                                1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options



  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume_p0[9] = {2.0833333333333332e-02, 2.0833333333333331e-01, 2.0833333333333332e-02,
                                       1.0416666666666666e-01, 4.1666666666666664e-02, 1.0416666666666666e-01,
                                       1.0416666666666666e-01, 1.0416666666666666e-01, 2.0833333333333331e-01};
  double expected_dual_volume_p1[10] = {2.0833333333333331e-01, 1.0416666666666666e-01, 4.1666666666666664e-02,
                                        1.0416666666666666e-01, 1.0416666666666666e-01, 4.1666666666666664e-02,
                                        1.0416666666666666e-01, 2.0833333333333332e-02, 2.0833333333333331e-01, 2.0833333333333332e-02};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);

  PDM_part_mesh_nodal_free(pmn);
}

MPI_TEST_CASE("[pdm_part_mesh_nodal_geom] 3d - HEXA8 - 2p", 2) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  int n_part = 1;

  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(pdm_comm,
                                                                PDM_MESH_NODAL_HEXA8, 1, NULL,
                                                                0., 0., 0., // x/y/z min
                                                                1., 1., 1., // x/y/z length
                                                                3, 3, 2, // x/y/z n vertices
                                                                1, PDM_SPLIT_DUAL_WITH_HILBERT); // part options
  double **pdual_volume = NULL;
  PDM_part_mesh_nodal_dual_volume_compute(pmn, PDM_TRUE, &pdual_volume);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      PDM_log_trace_array_double(pdual_volume[i_part], n_vtx, "pdual_volume ::");
    }
  }

  double expected_dual_volume_p0[12] = {3.1249999999999997e-02, 6.2500000000000000e-02, 3.1249999999999997e-02,
                                       6.2499999999999993e-02, 1.2500000000000000e-01, 6.2499999999999993e-02,
                                       3.1249999999999997e-02, 6.2500000000000000e-02, 3.1249999999999997e-02,
                                       6.2499999999999993e-02, 1.2500000000000000e-01, 6.2499999999999993e-02};
  double expected_dual_volume_p1[12] = {6.2499999999999993e-02, 1.2500000000000000e-01, 6.2499999999999993e-02,
                                        3.1249999999999997e-02, 6.2500000000000000e-02, 3.1249999999999997e-02,
                                        6.2499999999999993e-02, 1.2500000000000000e-01, 6.2499999999999993e-02,
                                        3.1249999999999997e-02, 6.2500000000000000e-02, 3.1249999999999997e-02};

  int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn, 0);
  for(int i_vtx = 0; i_vtx < n_vtx; ++i_vtx ) {
    MPI_CHECK(0, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p0[i_vtx]).epsilon(0.01));
    MPI_CHECK(1, pdual_volume[0][i_vtx] == doctest::Approx(expected_dual_volume_p1[i_vtx]).epsilon(0.01));
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pdual_volume[i_part]);
  }
  PDM_free(pdual_volume);


  PDM_part_mesh_nodal_free(pmn);
}
