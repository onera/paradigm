#include <stddef.h>
#include <vector>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_part_geom.h"


MPI_TEST_CASE("[pdm_part_geom] vtx normals from edges", 2) {

  /**
   *                  2  (1)  4
   *                  x-------x
   *                  |
   *                  | (3)
   *                  |
   *   rank0: x-------x 3
   *          1  (2)
   *
   *                  x 1
   *                  |
   *                  | (2)
   *                  |
   *   rank1: x-------x
   *          3  (1)  2
   *
   * Note:
   *   - edge 2 and 3 of rank 0 are the same that edge 1 and 2 of rank 1 (respectively)
   *   - this conf test the owner tricks (can occur in non-manifold conf)
   */

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);


  int n_part = 1;
  int dimension = 1;

  // > Mesh definition
  int n_edge = (i_rank==0) ? 3 : 2;
  int n_vtx  = (i_rank==0) ? 4 : 3;
  std::vector<std::vector<int>> vedge_vtx_idx = {{0,2,4,6},
                                                 {0,2,4  }};
  std::vector<std::vector<int>> vedge_vtx = {{2,4, 1,3, 3,2},
                                             {3,2, 2,1     }};
  std::vector<std::vector<double>> vcoord = {{0.,0.,0.,
                                              1.,1.,0.,
                                              1.,0.,0.,
                                              2.,1.,0.},
                                             {1.,1.,0.,
                                              1.,0.,0.,
                                              0.,0.,0.}};
  int    *edge_vtx_idx = vedge_vtx_idx[i_rank].data();
  int    *edge_vtx     = vedge_vtx    [i_rank].data();
  double *coord        = vcoord       [i_rank].data();

  // > Part comm graph definition
  int n_edge_graph = 2;
  std::vector<std::vector<int>> vedge_graph = {{3,1,1,2, 2,1,1,1},
                                               {1,0,1,2, 2,0,1,3}};
  int *edge_graph = vedge_graph[i_rank].data();
  PDM_part_comm_graph_t *pcg_edge = PDM_part_comm_graph_create(
    n_part,
   &n_edge_graph,
   &edge_graph,
    PDM_OWNERSHIP_USER,
    comm
  );

  int n_vtx_graph = 3;
  std::vector<std::vector<int>> vvtx_graph = {{1,1,1,3, 2,1,1,1, 3,1,1,2},
                                              {1,0,1,2, 2,0,1,3, 3,0,1,1}};
  int *vtx_graph = vvtx_graph[i_rank].data();
  PDM_part_comm_graph_t *pcg_vtx = PDM_part_comm_graph_create(
    n_part,
   &n_vtx_graph,
   &vtx_graph,
    PDM_OWNERSHIP_USER,
    comm
  );


  SUBCASE("full") {

    int n_selected_edge = n_edge;
    int n_selected_vtx  = n_vtx;
    double **vtx_normal = NULL;
    PDM_part_geom_vtx_normal_compute(
      comm,
      n_part,
      dimension,
     &n_selected_edge,
      NULL,
     &edge_vtx_idx,
     &edge_vtx,
      NULL,
      pcg_edge,
     &n_selected_vtx,
      NULL,
     &coord,
      pcg_vtx,
     &vtx_normal
    );

    std::vector<std::vector<double>> expected_vtx_normal = {{0.      , -1.      , 0.,
                                                             0.707107, -0.707107, 0.,
                                                             0.707107, -0.707107, 0.,
                                                             0.      , -1.      , 0.},
                                                            {0.707107, -0.707107, 0.,
                                                             0.707107, -0.707107, 0.,
                                                             0.      , -1.      , 0.}};

    for (int i_vtx=0; i_vtx<n_selected_vtx; ++i_vtx) {
      CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx  ]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+1]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+2]).epsilon(0.00001));
    }

    free(vtx_normal[0]);
    free(vtx_normal);
  }


  SUBCASE("partial edge") {

    int n_selected_edge = 1;
    int n_selected_vtx  = n_vtx;
    std::vector<std::vector<int>> vselected_edge = {{3},
                                                    {2}}; // what if it is not same edge selected ?
    int *selected_edge = vselected_edge[i_rank].data();

    double **vtx_normal = NULL;
    PDM_part_geom_vtx_normal_compute(
      comm,
      n_part,
      dimension,
     &n_selected_edge,
     &selected_edge,
     &edge_vtx_idx,
     &edge_vtx,
      NULL,
      pcg_edge,
     &n_selected_vtx,
      NULL,
     &coord,
      pcg_vtx,
     &vtx_normal
    );

    std::vector<std::vector<double>> expected_vtx_normal = {{0., 0., 0.,
                                                             1., 0., 0.,
                                                             1., 0., 0.,
                                                             0., 0., 0.},
                                                            {1., 0., 0.,
                                                             1., 0., 0.,
                                                             0., 0., 0.}};

    for (int i_vtx=0; i_vtx<n_selected_vtx; ++i_vtx) {
      CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx  ]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+1]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+2]).epsilon(0.00001));
    }

    free(vtx_normal[0]);
    free(vtx_normal);
  }


  SUBCASE("partial vtx") {

    int n_selected_edge = n_edge;
    int n_selected_vtx  = 2;
    std::vector<std::vector<int>> vselected_vtx = {{2,1},
                                                   {1,3}}; // what if it is not same vtx selected ?
    int *selected_vtx = vselected_vtx[i_rank].data();

    double **vtx_normal = NULL;
    PDM_part_geom_vtx_normal_compute(
      comm,
      n_part,
      dimension,
     &n_selected_edge,
      NULL,
     &edge_vtx_idx,
     &edge_vtx,
      NULL,
      pcg_edge,
     &n_selected_vtx,
     &selected_vtx,
     &coord,
      pcg_vtx,
     &vtx_normal
    );

    std::vector<std::vector<double>> expected_vtx_normal = {{0.707107, -0.707107, 0.,
                                                             0.      , -1.      , 0.},
                                                            {0.707107, -0.707107, 0.,
                                                             0.      , -1.      , 0.}};

    for (int i_vtx=0; i_vtx<n_selected_vtx; ++i_vtx) {
      CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx  ]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+1]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+2]).epsilon(0.00001));
    }

    free(vtx_normal[0]);
    free(vtx_normal);
  }


  SUBCASE("partial edge and vtx") {

    int n_selected_edge = 1;
    int n_selected_vtx  = 2;
    std::vector<std::vector<int>> vselected_edge = {{3},
                                                    {2}}; // what if it is not same edge selected ?
    std::vector<std::vector<int>> vselected_vtx = {{2,1},
                                                   {1,3}}; // what if it is not same vtx selected ?
    int *selected_edge = vselected_edge[i_rank].data();
    int *selected_vtx  = vselected_vtx [i_rank].data();

    double **vtx_normal = NULL;
    PDM_part_geom_vtx_normal_compute(
      comm,
      n_part,
      dimension,
     &n_selected_edge,
     &selected_edge,
     &edge_vtx_idx,
     &edge_vtx,
      NULL,
      pcg_edge,
     &n_selected_vtx,
     &selected_vtx,
     &coord,
      pcg_vtx,
     &vtx_normal
    );

    /**
     * Result can be strange because some referenced vertices won't have
     * edge contribution.
     * Is this ok, or should we check selected_edge/selected_vtx coherence ?
     */
    std::vector<std::vector<double>> expected_vtx_normal = {{1., 0., 0.,
                                                             0., 0., 0.},
                                                            {1., 0., 0.,
                                                             0., 0., 0.}};

    for (int i_vtx=0; i_vtx<n_selected_vtx; ++i_vtx) {
      CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx  ]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+1]).epsilon(0.00001));
      CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[i_rank][3*i_vtx+2]).epsilon(0.00001));
    }

    free(vtx_normal[0]);
    free(vtx_normal);
  }


  PDM_part_comm_graph_free(pcg_edge);
  PDM_part_comm_graph_free(pcg_vtx);
}


MPI_TEST_CASE("[pdm_part_geom] vtx normals from faces", 1) {

  /**
   *  2 x-------x 4
   *    |       |
   *    |       |
   *  1 x-------x 3
   *
   *   in XY plane, but rotated by 45° around Z
   *
   */

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_part = 1;
  int dimension = 2;

  // > Mesh definition
  int    n_face = 1;
  int    n_vtx  = 4;
  std::vector<int   > vface_vtx_idx = {0,4};
  std::vector<int   > vface_vtx     = {1,3,4,2};
  std::vector<double> vcoord        = {0., 0., 0.,
                                       1., 1., 0.,
                                       0., 0., 1.,
                                       1., 1., 1.};

  // int    face_vtx_idx[2] = {0,4};
  int    *face_vtx_idx = vface_vtx_idx.data();
  int    *face_vtx     = vface_vtx    .data();
  double *coord        = vcoord       .data();

  // > Part comm graph definition
  int  n_face_graph = 0;
  int *face_graph   = NULL;
  PDM_part_comm_graph_t *pcg_face = PDM_part_comm_graph_create(
    n_part,
   &n_face_graph,
   &face_graph,
    PDM_OWNERSHIP_USER,
    comm
  );

  int  n_vtx_graph = 0;
  int *vtx_graph   = NULL;
  PDM_part_comm_graph_t *pcg_vtx = PDM_part_comm_graph_create(
    n_part,
   &n_vtx_graph,
   &vtx_graph,
    PDM_OWNERSHIP_USER,
    comm
  );


  int n_selected_face = n_face;
  int n_selected_vtx  = n_vtx;
  double **vtx_normal = NULL;
  PDM_part_geom_vtx_normal_compute(
    comm,
    n_part,
    dimension,
    &n_selected_face,
    NULL,
   &face_vtx_idx,
   &face_vtx,
    NULL,
    pcg_face,
   &n_selected_vtx,
    NULL,
   &coord,
    pcg_vtx,
   &vtx_normal
  );

  std::vector<double> expected_vtx_normal = {-0.707107, 0.707107, 0.,
                                             -0.707107, 0.707107, 0.,
                                             -0.707107, 0.707107, 0.,
                                             -0.707107, 0.707107, 0.};

  for (int i_vtx=0; i_vtx<n_selected_vtx; ++i_vtx) {
    CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[3*i_vtx  ]).epsilon(0.00001));
    CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[3*i_vtx+1]).epsilon(0.00001));
    CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[3*i_vtx+2]).epsilon(0.00001));
  }

  free(vtx_normal[0]);
  free(vtx_normal);

  PDM_part_comm_graph_free(pcg_face);
  PDM_part_comm_graph_free(pcg_vtx);
}




MPI_TEST_CASE("[pdm_part_geom] vtx normals from edge normal", 1) {

  /**
   *
   *  x-------x
   *  1  (1)  2
   *
   */

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_part = 1;
  int dimension = 1;

  // > Mesh definition
  int n_edge = 1;
  int n_vtx  = 2;
  std::vector<int   > vedge_vtx_idx = {0,2};
  std::vector<int   > vedge_vtx     = {2,1};
  std::vector<double> vedge_normal  = {0.,1.,0.};
  std::vector<double> vcoord        = {0.,0.,0., 1.,0.,0.};
  int    *edge_vtx_idx = vedge_vtx_idx.data();
  int    *edge_vtx     = vedge_vtx    .data();
  double *edge_normal  = vedge_normal .data();
  double *coord        = vcoord       .data();

  // > Part comm graph definition
  int n_edge_graph = 0;
  int *edge_graph  = NULL;
  PDM_part_comm_graph_t *pcg_edge = PDM_part_comm_graph_create(
    n_part,
   &n_edge_graph,
   &edge_graph,
    PDM_OWNERSHIP_USER,
    comm
  );

  int n_vtx_graph = 0;
  int *vtx_graph  = NULL;
  PDM_part_comm_graph_t *pcg_vtx = PDM_part_comm_graph_create(
    n_part,
   &n_vtx_graph,
   &vtx_graph,
    PDM_OWNERSHIP_USER,
    comm
  );

  double **vtx_normal = NULL;
  PDM_part_geom_vtx_normal_compute(
    comm,
    n_part,
    dimension,
   &n_edge,
    NULL,
   &edge_vtx_idx,
   &edge_vtx,
   &edge_normal,
    pcg_edge,
   &n_vtx,
    NULL,
   &coord,
    pcg_vtx,
   &vtx_normal
  );

  std::vector<double> expected_vtx_normal = {0.,0.,1., 0.,0.,1.};

  for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
    CHECK(vtx_normal[0][3*i_vtx  ] == doctest::Approx(expected_vtx_normal[3*i_vtx  ]).epsilon(0.00001));
    CHECK(vtx_normal[0][3*i_vtx+1] == doctest::Approx(expected_vtx_normal[3*i_vtx+1]).epsilon(0.00001));
    CHECK(vtx_normal[0][3*i_vtx+2] == doctest::Approx(expected_vtx_normal[3*i_vtx+2]).epsilon(0.00001));
  }

  free(vtx_normal[0]);
  free(vtx_normal);

  PDM_part_comm_graph_free(pcg_edge);
  PDM_part_comm_graph_free(pcg_vtx);
}