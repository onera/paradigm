
#include <stdlib.h>
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_convex.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_priv.h"



TEST_CASE("[pdm_convex] - PDM_points_inside_convex_hull 2d") {

  const int dim   = 2;
  const int n_src = 8;
  const int n_tgt = 3;

  double src_coord[n_src*3] = {
     0.0,  0.0, 0,
    12.0,  2.0, 0,
     4.0,  4.0, 0,
    14.0,  8.0, 0,
    10.0,  6.0, 0,
     4.0, 10.0, 0,
     8.0,  2.0, 0,
     0.0,  6.0, 0
  };


  double tgt_coord[n_tgt*3] = {
     2.0,  2.0, 0, // inside convex hull
    18.0,  2.0, 0, // outside source bounding box
    10.0, 10.0, 0  // inside  source bounding box but outside convex hull
  };

  int tgt_status[n_tgt];
  PDM_points_inside_convex_hull(dim,
                                n_src,
                                src_coord,
                                n_tgt,
                                tgt_coord,
                                tgt_status);

  int expected_tgt_status[n_tgt] = {1, 0, 0};
  CHECK_EQ_C_ARRAY(tgt_status, expected_tgt_status, n_tgt);
}



TEST_CASE("[pdm_convex] - PDM_points_inside_convex_hull 3d") {

  const int dim   = 3;
  const int n_src = 12;
  const int n_tgt = 3;

  double src_coord[n_src*3] = {
     2.0, 4.0, -1.0,
     1.0, 2.0,  2.0,
    -1.0, 0.0,  1.0,
     1.0, 1.0,  1.0,
     0.0, 1.0,  0.0,
    -1.0, 3.0,  0.0,
    -1.0, 1.0,  2.0,
     3.0, 1.0, -1.0,
    -1.0, 1.0, -1.0,
     2.0, 2.0,  0.0,
     1.0, 2.0, -1.0,
     1.0, 1.0,  0.0
  };


  double tgt_coord[n_tgt*3] = {
    0.0, 1.0, 0.0, // inside convex hull
    4.0, 1.0, 1.0, // outside source bounding box
    1.0, 3.0, 1.0  // inside  source bounding box but outside convex hull
  };

  int tgt_status[n_tgt];
  PDM_points_inside_convex_hull(dim,
                                n_src,
                                src_coord,
                                n_tgt,
                                tgt_coord,
                                tgt_status);

  int expected_tgt_status[n_tgt] = {1, 0, 0};
  CHECK_EQ_C_ARRAY(tgt_status, expected_tgt_status, n_tgt);
}