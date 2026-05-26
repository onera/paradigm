#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mesh_intersection_surf_surf_atomic.h"
#include "pdm_priv.h"


TEST_CASE("Mesh intersection surf-surf atomic") {

  std::vector<double> coord    [2];
  std::vector<int>    face_vtx [2];
  std::vector<int>    face_edge[2];
  std::vector<int>    edge_vtx [2];

  for (int i_poly = 0; i_poly < 2; i_poly++) {
    face_vtx [i_poly].clear();
    face_edge[i_poly].clear();
    edge_vtx [i_poly].clear();
  }

  double exp_area_ab;
  std::vector<double> exp_center_ab;


  SUBCASE("Quad / triangle") {
    printf("Quad / triangle\n");
    // First polygon
    coord[0] = {
      0,    12./13, 0,
      1./7, 12./13, 0,
      1./7, 13./13, 0,
      0,    13./13, 0
    };
    face_vtx [0] = {1, 2, 3, 4};
    face_edge[0] = {-1, 3, 2, 4};
    edge_vtx [0] = {
      2, 1,
      3, 4,
      2, 3,
      4, 1
    };

    // Second polygon
    coord[1] = {
      0,      9./10,  0,
      1./10,  9./10,  0,
      0,     10./10, 0
    };
    face_vtx [1] = {1, 2, 3};
    face_edge[1] = {1, 2, 3};
    edge_vtx [1] = {
      1, 2,
      2, 3,
      3, 1
    };

    // Expected result
    exp_area_ab   = 1./(2*13*13);
    exp_center_ab = {1./(13*3), (1 + 2*12./13)/3., 0.};
  }




  SUBCASE("Convex / convex") {
    printf("Convex / convex\n");
    // First polygon
    coord[0] = {
      7, 1, 0,
      6, 5, 0,
      1, 5, 0,
      3, 2, 0
    };
    face_vtx [0] = {1, 2, 3, 4};
    face_edge[0] = {-1, 3, 2, 4};
    edge_vtx [0] = {
      2, 1,
      3, 4,
      2, 3,
      4, 1
    };

    // Second polygon
    coord[1] = {
      3, 5, 0,
      3, 1, 0,
      1, 3, 0,
      7, 3, 0,
      5, 5, 0
    };
    face_vtx [1] = {1, 3, 2, 4, 5};
    face_edge[1] = {-1, 2, 3, 4, -5};
    edge_vtx [1] = {
      3, 1,
      3, 2,
      2, 4,
      4, 5,
      1, 5
    };

    // Expected result
    exp_area_ab   = 10.911111111111111;
    exp_center_ab = {4.177129063890774, 3.428618842875462, 0.};
  }




  SUBCASE("Convex / convex flipped") {
    printf("Convex / convex flipped\n");
    // First polygon
    coord[0] = {
      7, 1, 0,
      6, 5, 0,
      1, 5, 0,
      3, 2, 0
    };
    face_vtx [0] = {1, 2, 3, 4};
    face_edge[0] = {-1, 2, 3, 4};
    edge_vtx [0] = {
      2, 1,
      2, 3,
      3, 4,
      4, 1
    };

    // Second polygon
    coord[1] = {
      3, 5, 0,
      3, 1, 0,
      1, 3, 0,
      7, 3, 0,
      5, 5, 0
    };
    face_vtx [1] = {5, 4, 2, 3, 1};
    face_edge[1] = {1, -2, -3, -4, 5};
    edge_vtx [1] = {
      3, 1,
      3, 2,
      2, 4,
      4, 5,
      1, 5
    };

    // Expected result
    exp_area_ab   = -10.911111111111111;
    exp_center_ab = {4.177129063890774, 3.428618842875462, 0.};
  }




  SUBCASE("Non-convex / non-convex") {
    printf("Non-convex / non-convex\n");
    // First polygon
    coord[0] = {
      0, 1, 0,
      4, 2, 0,
      8, 0, 0,
      6, 4, 0,
      3, 5, 0
    };
    face_vtx [0] = {1, 2, 3, 4, 5};
    face_edge[0] = {1, 2, -3, -4, 5};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      4, 3,
      5, 4,
      5, 1
    };


    // Second polygon
    coord[1] = {
      3, 0, 0,
      8, 3, 0,
      5, 3, 0,
      1, 4, 0
    };
    face_vtx [1] = {1, 2, 3, 4};
    face_edge[1] = {1, -2, 3, 4};
    edge_vtx [1] = {
      1, 2,
      3, 2,
      3, 4,
      4, 1
    };

    // Expected result
    exp_area_ab   = 7.3737782685151121;
    exp_center_ab = {3.9757700808959036, 2.5251499884012016, 0.};
  }




  SUBCASE("Non-convex / convex") {
    printf("Non-convex / convex\n");
    // First polygon
    coord[0] = {
      0.5, 0.9, 0,
      0.2, 0.8, 0,
      0.1, 0.2, 0,
      0.9, 0.1, 0,
      1.0, 0.5, 0,
      0.4, 0.4, 0
    };
    face_vtx [0] = {1, 2, 3, 4, 5, 6};
    face_edge[0] = {1, 2, 3, 4, 5, 6};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 4,
      4, 5,
      5, 6,
      6, 1
    };


    // Second polygon
    coord[1] = {
      0.6, 1.1, 0,
      0.5, 0.5, 0,
      1.3, 0.6, 0,
      1.1, 0.9, 0
    };
    face_vtx [1] = {1, 2, 3, 4};
    face_edge[1] = {1, 2, 3, 4};
    edge_vtx [1] = {
      1, 2,
      2, 3,
      3, 4,
      4, 1
    };

    // Expected result
    exp_area_ab   = 0.;
    exp_center_ab = {0, 0, 0};
  }




  SUBCASE("Anisotropic triangles") {
    printf("Anisotropic triangles\n");
    // First polygon
    coord[0] = {
      1.76479,  0.149932, 0,
      0.064595, 0.95737,  0,
      1.95723,  0.047524, 0
    };
    face_vtx [0] = {1, 2, 3};
    face_edge[0] = {1, 2, 3};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 1
    };


    // Second polygon
    coord[1] = {
      0.631587, 0.655398, 0,
      1.80397,  0.15946,  0,
      0.192275, 0.846864, 0
    };
    face_vtx [1] = {1, 2, 3};
    face_edge[1] = {1, -2, 3};
    edge_vtx [1] = {
      1, 2,
      3, 2,
      3, 1
    };

    // Expected result
    exp_area_ab   = 0.00026652356360320085;
    exp_center_ab = {1.17907160607039, 0.42492320217419655, 0.};
  }




  SUBCASE("Star / star") {
    printf("Star / star\n");
    // First polygon
    coord[0] = {
      -1.0, -1.0, 0,
      -0.1, -0.2, 0,
       0.0, -1.0, 0,
       0.1, -0.2, 0,
       1.0, -1.0, 0,
       0.2, -0.1, 0,
       1.0,  0.0, 0,
       0.2,  0.1, 0,
       1.0,  1.0, 0,
       0.1,  0.2, 0,
       0.0,  1.0, 0,
      -0.1,  0.2, 0,
      -1.0,  1.0, 0,
      -0.2,  0.1, 0,
      -1.0,  0.0, 0,
      -0.2, -0.1, 0
    };
    face_vtx [0] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    face_edge[0] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 4,
      4, 5,
      5, 6,
      6, 7,
      7, 8,
      8, 9,
      9, 10,
      10, 11,
      11, 12,
      12, 13,
      13, 14,
      14, 15,
      15, 16,
      16, 1
    };


    // Second polygon
    coord[1] = {
       0.3,  0.0, 0,
       0.8,  0.6, 0,
       0.1,  0.3, 0,
      -0.3,  1.0, 0,
      -0.3,  0.2, 0,
      -1.0,  0.0, 0,
      -0.3, -0.2, 0,
      -0.3, -1.0, 0,
       0.1, -0.3, 0,
       0.8, -0.6, 0,
    };
    face_vtx [1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    face_edge[1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    edge_vtx [1] = {
      1, 2,
      2, 3,
      3, 4,
      4, 5,
      5, 6,
      6, 7,
      7, 8,
      8, 9,
      9, 10,
      10, 1
    };

    // Expected result
    exp_area_ab   = 0.4837148477164416;
    exp_center_ab = {-0.031173442189319722, -4.7816697080004124e-18, 0.};
  }




  SUBCASE("Hole / complementary") {
    printf("Hole / complementary\n");
    // First polygon
    coord[0] = {
      0.9, 0.1, 0,
      1.3, 0.7, 0,
      0.6, 1.2, 0,
      0.2, 0.6, 0,
      0.8, 0.4, 0,
      0.5, 0.6, 0,
      0.7, 0.9, 0,
      1.0, 0.6, 0
    };
    face_edge[0] = {1, 2, 3, 4, 5, 6, 7, 8};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 4,
      4, 1,
      5, 6,
      6, 7,
      7, 8,
      8, 5
    };


    // Second polygon
    coord[1] = {
      0.8, 0.4, 0,
      0.5, 0.6, 0,
      0.7, 0.9, 0,
      1.0, 0.6, 0
    };
    face_vtx [1] = {1, 2, 3, 4};
    face_edge[1] = {1, -2, 3, 4};
    edge_vtx [1] = {
      1, 2,
      3, 2,
      3, 4,
      4, 1
    };

    // Expected result
    exp_area_ab   = 0.;
    exp_center_ab = {0, 0, 0};
  }




  SUBCASE("Hole / triangle") {
    printf("Hole / triangle\n");
    // First polygon
    coord[0] = {
      0.9, 0.1, 0,
      1.3, 0.7, 0,
      0.6, 1.2, 0,
      0.2, 0.6, 0,
      0.8, 0.4, 0,
      0.5, 0.6, 0,
      0.7, 0.9, 0,
      1.0, 0.6, 0
    };
    face_edge[0] = {1, 2, 3, 4, 5, 6, 7, 8};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 4,
      4, 1,
      5, 6,
      6, 7,
      7, 8,
      8, 5
    };


    // Second polygon
    coord[1] = {
      1.1, 1.3, 0,
      0.3, 0.2, 0,
      0.7, 0.1, 0
    };
    face_vtx [1] = {1, 2, 3};
    face_edge[1] = {3, 1, 2};
    edge_vtx [1] = {
      3, 1,
      1, 2,
      2, 3
    };

    // Expected result
    exp_area_ab   = 0.0903333884160952;
    exp_center_ab = {0.725973010020932, 0.5631531345202446, 0};
  }




  SUBCASE("With 3D rotation") {
    printf("With 3D rotation\n");
    // First polygon
    coord[0] = {
      -0.5231886086180316, -0.0264123129141687, -0.2560782488593127,
       0.2095440742422240,  0.3871745186455422,  0.0786585854281783,
       0.4659181198269654, -0.2503669164800132,  0.2454316865223308,
      -0.1145409775821326, -0.5535584587883637, -0.0212931246739639
    };
    face_vtx [0] = {1, 2, 3, 4};
    face_edge[0] = {1, 2, 3, 4};
    edge_vtx [0] = {
      1, 2,
      2, 3,
      3, 4,
      4, 1
    };

    // Second polygon
    coord[1] = {
      -0.1138694005946041,  0.3879655249113231, -0.0807249099628307,
       0.5043223046835190,  0.2479938976826703,  0.2327187996506965,
       0.1802372528591624, -0.6927390797512355,  0.1327670895485542,
      -0.2193130147032446, -0.7471457181111347, -0.0606224173880530
    };
    face_vtx [1] = {1, 2, 3, 4};
    face_edge[1] = {1, 2, 3, 4};
    edge_vtx [1] = {
      1, 2,
      2, 3,
      3, 4,
      4, 1
    };

    // Expected result
    exp_area_ab   = 0.4012326802507837;
    exp_center_ab = {0.0793335379212110, -0.0979580292581317, 0.0453023845172782};
  }



  PDM_mesh_intersection_surf_surf_polygon_t poly[2];
  std::vector<double> coord_noise[2];

  srand(0);
  double max_error_area   = 0.;
  double max_error_center = 0.;

  for (int i_mode0 = 0; i_mode0 < 2; i_mode0++) {

    for (int i_mode1 = 0; i_mode1 < 2; i_mode1++) {

      int from_face_vtx[2] = {i_mode0, i_mode1};

      // Connectivity
      for (int i_poly = 0; i_poly < 2; i_poly++) {
        if (from_face_vtx[i_poly] && face_vtx[i_poly].size() > 0) {
          poly[i_poly].n_edge    = face_vtx[i_poly].size();
          poly[i_poly].face_vtx  = (int *) face_vtx[i_poly].data();
          poly[i_poly].face_edge = NULL;
          poly[i_poly].edge_vtx  = NULL;
        }
        else {
          poly[i_poly].n_edge    = face_edge[i_poly].size();
          poly[i_poly].face_vtx  = NULL;
          poly[i_poly].face_edge = (int *) face_edge[i_poly].data();
          poly[i_poly].edge_vtx  = (int *) edge_vtx [i_poly].data();
        }
      }

      // Coordinates
      for (int i_noise = 0; i_noise < 1000; i_noise++) {

        for (int i_poly = 0; i_poly < 2; i_poly++) {
          coord_noise[i_poly].resize(coord[i_poly].size());
          for (size_t i = 0; i < coord[i_poly].size(); i++) {
            coord_noise[i_poly][i] = coord[i_poly][i];
            if (i_noise > 0) {
              // Add small random perturbations
              double noise = rand() / ((double) RAND_MAX) - 0.5;
              coord_noise[i_poly][i] += noise * 1e-15;
            }
          }

          poly[i_poly].coord = (double *) coord_noise[i_poly].data();
        }


        // Compute intersection
        double area_ab;
        double center_ab[3];
        PDM_mesh_intersection_surf_surf_atomic_compute(&poly[0],
                                                       &poly[1],
                                                       &area_ab,
                                                       center_ab);

        double error_area = abs(area_ab - exp_area_ab);
        max_error_area = std::max(max_error_area, error_area);

        if (exp_area_ab > 0) {
          double error_center = 0;
          for (int i = 0; i < 3; i++) {
            error_center += (center_ab[i] - exp_center_ab[i]) * (center_ab[i] - exp_center_ab[i]);
          }
          max_error_center = std::max(max_error_center, sqrt(error_center));
        }
      }
    }
  }

  // Check result
  printf("  max error area   = %e\n",   max_error_area);
  printf("  max error center = %e\n\n", max_error_center);
  CHECK(max_error_area   < 1e-12);
  CHECK(max_error_center < 1e-12);
}
