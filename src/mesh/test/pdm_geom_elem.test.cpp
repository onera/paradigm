#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_doctest.h"
#include "pdm_geom_elem.h"
#include "pdm_logging.h"
#include <vector>



MPI_TEST_CASE("PDM_geom_elem_edge_upwind_and_downwind_2d", 1) {

  std::vector<PDM_g_num_t> face_ln_to_gn = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<int        > face_edge_idx = {0,3,6,9,12,15,18,21,24};
  std::vector<int        > face_edge     = {1,2,4,-4,5,7,-5,3,6,-6,8,9,-7,10,11,-11,12,14,-12,-9,13,-13,15,16};

  std::vector<int        > edge_vtx      = {1,2,4,1,2,3,2,4,2,5,3,5,5,4,3,6,6,5,7,4,5,7,5,8,6,8,8,7,6,9,9,8};
  std::vector<int        > vtx_face_idx  = {0,1,4,6,9,15,18,20,23,24};
  std::vector<int        > vtx_face      = {1,1,2,3,3,4,1,2,5,2,3,4,5,6,7,4,7,8,5,6,6,7,8,8};
  std::vector<double     > vtx_coords    = {3.4018771715470957e-02,-1.0561707318090696e-02,0.0000000000000000e+00,
                                            5.2830992237586061e-01,2.9844003347607329e-02,0.0000000000000000e+00,
                                            1.0411647357936784e+00,-3.0244863070661605e-02,0.0000000000000000e+00,
                                            -1.6477724428511097e-02,5.2682295948119040e-01,0.0000000000000000e+00,
                                            4.7777747108031876e-01,5.0539699557954310e-01,0.0000000000000000e+00,
                                            9.9773970518621602e-01,5.1288709247619246e-01,0.0000000000000000e+00,
                                            -1.3521552720815667e-02,1.0013400910195616e+00,0.0000000000000000e+00,
                                            5.4522297251747132e-01,1.0416195068003700e+00,0.0000000000000000e+00,
                                            1.0135711727959902e+00,1.0217296929432682e+00,0.0000000000000000e+00};

  int n_vtx = vtx_face_idx.size()-1;

  int i_plane;
  SUBCASE("XY") {
    i_plane = 0;
  }

  SUBCASE("YZ") {
    i_plane = 1;
    for (int i = 0; i < n_vtx; i++) {
      double x = vtx_coords[3*i  ];
      double y = vtx_coords[3*i+1];
      double z = vtx_coords[3*i+2];
      vtx_coords[3*i  ] = z;
      vtx_coords[3*i+1] = x;
      vtx_coords[3*i+2] = y;
    }
  }

  SUBCASE("ZX") {
    i_plane = 2;
    for (int i = 0; i < n_vtx; i++) {
      double x = vtx_coords[3*i  ];
      double y = vtx_coords[3*i+1];
      double z = vtx_coords[3*i+2];
      vtx_coords[3*i  ] = y;
      vtx_coords[3*i+1] = z;
      vtx_coords[3*i+2] = x;
    }
  }


  int n_edge = edge_vtx.size()/2;

  int    *upwind_face_out    = NULL;
  int    *downwind_face_out  = NULL;
  int    *upwind_edge_out    = NULL;
  int    *downwind_edge_out  = NULL;
  double *upwind_point_out   = NULL;
  double *downwind_point_out = NULL;
  PDM_geom_elem_edge_upwind_and_downwind_2d(i_plane,
                                            face_ln_to_gn.data(),
                                            NULL,
                                            face_edge_idx.data(),
                                            face_edge.data(),
                                            n_edge,
                                            edge_vtx.data(),
                                            vtx_face_idx.data(),
                                            vtx_face.data(),
                                            vtx_coords.data(),
                                            &upwind_face_out,
                                            &downwind_face_out,
                                            &upwind_edge_out,
                                            &downwind_edge_out,
                                            &upwind_point_out,
                                            &downwind_point_out);


  int upwind_face_out_expected  [16] = {2,-1,-1,-1,5,4,-1,7,1,-1,-1,-1,-1,-1,-1,-1};
  int downwind_face_out_expected[16] = {-1,-1,0,-1,-1,-1,3,-1,-1,-1,2,1,-1,-1,3,-1};
  int upwind_edge_out_expected  [16] = {5,-1,-1,-1,13,9,-1,15,3,-1,-1,-1,-1,-1,-1,-1};
  int downwind_edge_out_expected[16] = {-1,-1,1,-1,-1,-1,7,-1,-1,-1,2,3,-1,-1,5,-1};

  CHECK_EQ_C_ARRAY(upwind_face_out  , upwind_face_out_expected  , n_edge);
  CHECK_EQ_C_ARRAY(downwind_face_out, downwind_face_out_expected, n_edge);
  CHECK_EQ_C_ARRAY(upwind_edge_out  , upwind_edge_out_expected  , n_edge);
  CHECK_EQ_C_ARRAY(downwind_edge_out, downwind_edge_out_expected, n_edge);

  PDM_free(upwind_face_out   );
  PDM_free(downwind_face_out );
  PDM_free(upwind_edge_out   );
  PDM_free(downwind_edge_out );
  PDM_free(upwind_point_out  );
  PDM_free(downwind_point_out);

}
