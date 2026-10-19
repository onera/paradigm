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
   *   - edge 2 and 3 of rank 0 are the same as edge 1 and 2 of rank 1 (respectively)
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
     &n_vtx,
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
                                                    {2}};
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
     &n_vtx,
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
                                                   {1,3}};
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
     &n_vtx,
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
    int n_selected_vtx  = 1;
    std::vector<std::vector<int>> vselected_edge = {{3},
                                                    {2}};
    std::vector<std::vector<int>> vselected_vtx = {{2,1},
                                                   {1,3}};
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
     &n_vtx,
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
   &n_vtx,
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
   &n_vtx,
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


MPI_TEST_CASE("[pdm_part_geom] edge center", 1) {

  /**
   *          3 x (2,2)
   *            |
   *            | (2)
   *            |
   *  1 x-------x 2
   *   (0,0) (1) (2,0)
   */

  int n_part = 1;
  int n_edge = 2;

  std::vector<int>    vedge_vtx = {1,2, 2,3};
  std::vector<double> vcoord    = {0.,0.,0.,
                                   2.,0.,0.,
                                   2.,2.,0.};
  int    *edge_vtx = vedge_vtx.data();
  double *coord    = vcoord   .data();

  double **edge_center = NULL;
  PDM_part_geom_edge_center(n_part,
                           &n_edge,
                            NULL,          /* no extract */
                           &edge_vtx,
                           &coord,
                           &edge_center);

  std::vector<double> expected_center = {1.,0.,0.,    /* middle of edge 1 */
                                         2.,1.,0.};   /* middle of edge 2 */

  for (int i = 0; i < 3*n_edge; ++i) {
    CHECK(edge_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(edge_center[0]);
  PDM_free(edge_center);
}



MPI_TEST_CASE("[pdm_part_geom] edge center - extract_lnum", 1) {

  /**
   *          3 x (2,2)
   *            |
   *            | (2)    <- selected
   *            |
   *  1 x-------x 2
   *   (0,0) (1) (2,0)
   */
  int n_part    = 1;
  int n_extract = 1;

  std::vector<int>    vedge_vtx     = {1,2, 2,3};
  std::vector<int>    vextract_lnum = {2};
  std::vector<double> vcoord        = {0.,0.,0.,
                                       2.,0.,0.,
                                       2.,2.,0.};
  int    *edge_vtx     = vedge_vtx    .data();
  int    *extract_lnum = vextract_lnum.data();
  double *coord        = vcoord       .data();

  double **edge_center = NULL;
  PDM_part_geom_edge_center(n_part,
                           &n_extract,
                           &extract_lnum,
                           &edge_vtx,
                           &coord,
                           &edge_center);

  /* The output is indexed by the SELECTION, not by the edge number: size 3*n_extract,
   * edge 2 lands at position 0. */
  std::vector<double> expected_center = {2.,1.,0.};

  for (int i = 0; i < 3*n_extract; ++i) {
    CHECK(edge_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(edge_center[0]);
  PDM_free(edge_center);
}



MPI_TEST_CASE("[pdm_part_geom] face center from vtx", 1) {

  /**
   *   3 x (0,3)                       7 x-------x 6
   *     | \                             |       |
   *     |   \     (1)              (2)  |       |
   *     |     \                         |       |
   *   1 x-------x 2                    4 x-------x 5
   *    (0,0)   (3,0)                   (10,0)   (12,0)   top at y = 2
   *
   *  A triangle and a quadrangle: the face counts differ, so an index mixed up
   *  between the two shows up in the normalization.
   */


  int n_part = 1;
  int n_face = 2;

  std::vector<int>    vface_vtx_idx = {0,3,7};
  std::vector<int>    vface_vtx     = {1,2,3, 4,5,6,7};
  std::vector<double> vcoord        = { 0.,0.,0.,
                                        3.,0.,0.,
                                        0.,3.,0.,
                                       10.,0.,0.,
                                       12.,0.,0.,
                                       12.,2.,0.,
                                       10.,2.,0.};
  int    *face_vtx_idx = vface_vtx_idx.data();
  int    *face_vtx     = vface_vtx    .data();
  double *coord        = vcoord       .data();

  double **face_center = NULL;
  PDM_part_geom_face_center_from_vtx(n_part,
                                    &n_face,
                                     NULL,          /* no extract */
                                    &face_vtx_idx,
                                    &face_vtx,
                                    &coord,
                                    &face_center);

  std::vector<double> expected_center = { 1.,1.,0.,    /* triangle   */
                                         11.,1.,0.};   /* quadrangle */

  for (int i = 0; i < 3*n_face; ++i) {
    CHECK(face_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(face_center[0]);
  PDM_free(face_center);
}



MPI_TEST_CASE("[pdm_part_geom] face center from vtx - extract_lnum", 1) {

  /**
   *  Same two faces, selected in REVERSE order: extract_lnum = {2,1}
   */
  int n_part    = 1;
  int n_extract = 2;

  std::vector<int>    vface_vtx_idx = {0,3,7};
  std::vector<int>    vface_vtx     = {1,2,3, 4,5,6,7};
  std::vector<int>    vextract_lnum = {2,1};
  std::vector<double> vcoord        = { 0.,0.,0.,
                                        3.,0.,0.,
                                        0.,3.,0.,
                                       10.,0.,0.,
                                       12.,0.,0.,
                                       12.,2.,0.,
                                       10.,2.,0.};
  int    *face_vtx_idx = vface_vtx_idx.data();
  int    *face_vtx     = vface_vtx    .data();
  int    *extract_lnum = vextract_lnum.data();
  double *coord        = vcoord       .data();

  double **face_center = NULL;
  PDM_part_geom_face_center_from_vtx(n_part,
                                    &n_extract,
                                    &extract_lnum,
                                    &face_vtx_idx,
                                    &face_vtx,
                                    &coord,
                                    &face_center);

  /* Selection order, not face order: the quadrangle comes first.
   *
   * This is the convention _mesh_repart_2d gets wrong on its HILBERT path, where it
   * compacts the result a second time with a full-mesh index. */
  std::vector<double> expected_center = {11.,1.,0.,    /* face 2 */
                                          1.,1.,0.};   /* face 1 */

  for (int i = 0; i < 3*n_extract; ++i) {
    CHECK(face_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(face_center[0]);
  PDM_free(face_center);
}


MPI_TEST_CASE("[pdm_part_geom] face center from edge", 1) {

  /**
   *  The same two faces as the from_vtx case, described through their edges, so the
   *  expected values are the same ones. On a closed face each vertex belongs to
   *  exactly two edges, so the mean of the edge midpoints equals the mean of the
   *  vertices -- the two functions cannot legitimately disagree.
   *
   *  Edges:  1 = (1,2)  2 = (2,3)  3 = (3,1)
   *          4 = (4,5)  5 = (5,6)  6 = (6,7)  7 = (7,4)
   *
   *  Signs are scattered in face_edge: the function takes PDM_ABS.
   */
  int n_part = 1;
  int n_face = 2;

  std::vector<int>    vface_edge_idx = {0,3,7};
  std::vector<int>    vface_edge     = {1,2,-3, 4,5,-6,-7};
  std::vector<int>    vedge_vtx      = {1,2, 2,3, 3,1, 4,5, 5,6, 6,7, 7,4};
  std::vector<double> vcoord         = { 0.,0.,0.,
                                         3.,0.,0.,
                                         0.,3.,0.,
                                        10.,0.,0.,
                                        12.,0.,0.,
                                        12.,2.,0.,
                                        10.,2.,0.};
  int    *face_edge_idx = vface_edge_idx.data();
  int    *face_edge     = vface_edge    .data();
  int    *edge_vtx      = vedge_vtx     .data();
  double *coord         = vcoord        .data();

  double **face_center = NULL;
  PDM_part_geom_face_center_from_edge(n_part,
                                     &n_face,
                                      NULL,          /* no extract */
                                     &face_edge_idx,
                                     &face_edge,
                                     &edge_vtx,
                                     &coord,
                                     &face_center);

  std::vector<double> expected_center = { 1.,1.,0.,    /* triangle   */
                                         11.,1.,0.};   /* quadrangle */

  for (int i = 0; i < 3*n_face; ++i) {
    CHECK(face_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(face_center[0]);
  PDM_free(face_center);
}


MPI_TEST_CASE("[pdm_part_geom] face center from edge - extract_lnum", 1) {

  /**
   *  Same two faces, selected in REVERSE order: extract_lnum = {2,1}
   */
  int n_part    = 1;
  int n_extract = 2;

  std::vector<int>    vface_edge_idx = {0,3,7};
  std::vector<int>    vface_edge     = {1,2,-3, 4,5,-6,-7};
  std::vector<int>    vedge_vtx      = {1,2, 2,3, 3,1, 4,5, 5,6, 6,7, 7,4};
  std::vector<int>    vextract_lnum  = {2,1};
  std::vector<double> vcoord         = { 0.,0.,0.,
                                         3.,0.,0.,
                                         0.,3.,0.,
                                        10.,0.,0.,
                                        12.,0.,0.,
                                        12.,2.,0.,
                                        10.,2.,0.};
  int    *face_edge_idx = vface_edge_idx.data();
  int    *face_edge     = vface_edge    .data();
  int    *edge_vtx      = vedge_vtx     .data();
  int    *extract_lnum  = vextract_lnum .data();
  double *coord         = vcoord        .data();

  double **face_center = NULL;
  PDM_part_geom_face_center_from_edge(n_part,
                                     &n_extract,
                                     &extract_lnum,
                                     &face_edge_idx,
                                     &face_edge,
                                     &edge_vtx,
                                     &coord,
                                     &face_center);

  std::vector<double> expected_center = {11.,1.,0.,    /* face 2 */
                                          1.,1.,0.};   /* face 1 */

  for (int i = 0; i < 3*n_extract; ++i) {
    CHECK(face_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(face_center[0]);
  PDM_free(face_center);
}



MPI_TEST_CASE("[pdm_part_geom] cell center from vtx", 1) {

  /**
   *  cell 1 : tetrahedron, vertices 1..4 at (0,0,0) (3,0,0) (0,3,0) (0,0,3)
   *           4 triangular faces, centers (1,1,0) (1,0,1) (0,1,1) (1,1,1)
   *           mean -> (3/4, 3/4, 3/4)
   *
   *  cell 2 : unit cube, vertices 5..12, x in [10,11], y and z in [0,1]
   *           6 quadrangular faces
   *           mean -> the cube center, (10.5, 0.5, 0.5)
   */
  int n_part = 1;
  int n_cell = 2;

  std::vector<int> vcell_face_idx = {0,4,10};
  std::vector<int> vcell_face     = {1,-2,3,-4,  -5,6,-7,8,-9,10};

  std::vector<int> vface_vtx_idx  = {0,3,6,9,12,16,20,24,28,32,36};
  std::vector<int> vface_vtx      = {1,2,3,  1,2,4,  1,3,4,  2,3,4,
                                     5, 6, 7, 8,
                                     9,10,11,12,
                                     5, 6,10, 9,
                                     8, 7,11,12,
                                     5, 8,12, 9,
                                     6, 7,11,10};

  std::vector<double> vcoord = { 0.,0.,0.,
                                 3.,0.,0.,
                                 0.,3.,0.,
                                 0.,0.,3.,
                                10.,0.,0.,
                                11.,0.,0.,
                                11.,1.,0.,
                                10.,1.,0.,
                                10.,0.,1.,
                                11.,0.,1.,
                                11.,1.,1.,
                                10.,1.,1.};

  int    *cell_face_idx = vcell_face_idx.data();
  int    *cell_face     = vcell_face    .data();
  int    *face_vtx_idx  = vface_vtx_idx .data();
  int    *face_vtx      = vface_vtx     .data();
  double *coord         = vcoord        .data();

  double **cell_center = NULL;
  PDM_part_geom_cell_center(n_part,
                           &n_cell,
                            NULL,              /* no extract */
                           &cell_face_idx,
                           &cell_face,
                            NULL,
                            NULL,        /* no face_edge */
                           &face_vtx_idx,
                           &face_vtx,
                            NULL,              /* no edge_vtx  */
                           &coord,
                           &cell_center);

  std::vector<double> expected_center = { 0.75,0.75,0.75,    /* tetrahedron */
                                         10.5 ,0.5 ,0.5 };   /* cube center */

  for (int i = 0; i < 3*n_cell; ++i) {
    CHECK(cell_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(cell_center[0]);
  PDM_free(cell_center);
}


MPI_TEST_CASE("[pdm_part_geom] cell center from vtx - extract_lnum", 1) {

  /**
   *  Same two cells, selecting the CUBE only: extract_lnum = {2}.
   */
  int n_part    = 1;
  int n_extract = 1;

  std::vector<int> vcell_face_idx = {0,4,10};
  std::vector<int> vcell_face     = {1,-2,3,-4,  -5,6,-7,8,-9,10};
  std::vector<int> vextract_lnum  = {2};

  std::vector<int> vface_vtx_idx  = {0,3,6,9,12,16,20,24,28,32,36};
  std::vector<int> vface_vtx      = {1,2,3,  1,2,4,  1,3,4,  2,3,4,
                                     5, 6, 7, 8,
                                     9,10,11,12,
                                     5, 6,10, 9,
                                     8, 7,11,12,
                                     5, 8,12, 9,
                                     6, 7,11,10};

  std::vector<double> vcoord = { 0.,0.,0.,
                                 3.,0.,0.,
                                 0.,3.,0.,
                                 0.,0.,3.,
                                10.,0.,0.,
                                11.,0.,0.,
                                11.,1.,0.,
                                10.,1.,0.,
                                10.,0.,1.,
                                11.,0.,1.,
                                11.,1.,1.,
                                10.,1.,1.};

  int    *cell_face_idx = vcell_face_idx.data();
  int    *cell_face     = vcell_face    .data();
  int    *extract_lnum  = vextract_lnum .data();
  int    *face_vtx_idx  = vface_vtx_idx .data();
  int    *face_vtx      = vface_vtx     .data();
  double *coord         = vcoord        .data();

  double **cell_center = NULL;
  PDM_part_geom_cell_center(n_part,
                           &n_extract,
                           &extract_lnum,
                           &cell_face_idx,
                           &cell_face,
                            NULL,
                            NULL,
                           &face_vtx_idx,
                           &face_vtx,
                            NULL,
                           &coord,
                           &cell_center);

  std::vector<double> expected_center = {10.5,0.5,0.5};   /* cube center */

  for (int i = 0; i < 3*n_extract; ++i) {
    CHECK(cell_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(cell_center[0]);
  PDM_free(cell_center);
}


MPI_TEST_CASE("[pdm_part_geom] cell center from edge", 1) {

  /**
   *  The same two cells, described through face->edge and edge->vtx, so the expected
   *  values are the same ones as the from_vtx case. On a closed face each vertex
   *  belongs to exactly two edges, so the mean of the edge midpoints equals the mean
   *  of the vertices -- the two branches cannot legitimately disagree.
   *
   *  Edges:  1..6   tetrahedron: (1,2) (2,3) (3,1) (1,4) (2,4) (3,4)
   *          7..10  cube bottom: (5,6) (6,7) (7,8) (8,5)
   *         11..14  cube top   : (9,10) (10,11) (11,12) (12,9)
   *         15..18  cube sides : (5,9) (6,10) (7,11) (8,12)
   */


  int n_part = 1;
  int n_cell = 2;

  std::vector<int> vcell_face_idx = {0,4,10};
  std::vector<int> vcell_face     = {1,-2,3,-4,  -5,6,-7,8,-9,10};

  std::vector<int> vface_edge_idx = {0,3,6,9,12,16,20,24,28,32,36};
  std::vector<int> vface_edge     = {  1,  2,  3,            /* f1  */
                                       1,  5, -4,            /* f2  */
                                       3,  6, -4,            /* f3  */
                                       2,  6, -5,            /* f4  */
                                       7,  8,  9, 10,        /* f5  */
                                      11, 12, 13, 14,        /* f6  */
                                       7, 16,-11,-15,        /* f7  */
                                      -9, 17, 13, 18,        /* f8  */
                                      10, 18,-14,-15,        /* f9  */
                                       8, 17,-12,-16};       /* f10 */

  std::vector<int> vedge_vtx = {1,2,  2,3,  3,1,  1,4,  2,4,  3,4,
                                5, 6,   6, 7,   7, 8,   8, 5,
                                9,10,  10,11,  11,12,  12, 9,
                                5, 9,   6,10,   7,11,   8,12};

  std::vector<double> vcoord = { 0.,0.,0.,
                                 3.,0.,0.,
                                 0.,3.,0.,
                                 0.,0.,3.,
                                10.,0.,0.,
                                11.,0.,0.,
                                11.,1.,0.,
                                10.,1.,0.,
                                10.,0.,1.,
                                11.,0.,1.,
                                11.,1.,1.,
                                10.,1.,1.};

  int    *cell_face_idx = vcell_face_idx.data();
  int    *cell_face     = vcell_face    .data();
  int    *face_edge_idx = vface_edge_idx.data();
  int    *face_edge     = vface_edge    .data();
  int    *edge_vtx      = vedge_vtx     .data();
  double *coord         = vcoord        .data();

  double **cell_center = NULL;
  PDM_part_geom_cell_center(n_part,
                           &n_cell,
                            NULL,              /* no extract */
                           &cell_face_idx,
                           &cell_face,
                           &face_edge_idx,
                           &face_edge,
                            NULL,
                            NULL,        /* no face_vtx */
                           &edge_vtx,
                           &coord,
                           &cell_center);

  std::vector<double> expected_center = { 0.75,0.75,0.75,    /* same as from_vtx */
                                         10.5 ,0.5 ,0.5 };

  for (int i = 0; i < 3*n_cell; ++i) {
    CHECK(cell_center[0][i] == doctest::Approx(expected_center[i]).epsilon(1e-12));
  }

  PDM_free(cell_center[0]);
  PDM_free(cell_center);
}
