#include "doctest/extensions/doctest_mpi.h"

#include <limits.h>
#include <float.h>
#include <math.h>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_priv.h"

#include "pdm_ho_basis.h"


/*
  Compute weights for different configuration
 */

static const double tol = 1e-14;

MPI_TEST_CASE("[pdm_ho_basis] edge",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_BAR2;
    const int order = 1;
    const int n_nodes = order + 1;
    double uvw[1] = {0.45};
    double weights_expected[n_nodes*n_pts] = {0.55, 0.45};    
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_BARHO;
    const int order = 2;
    const int n_nodes = order + 1;
    double uvw[1] = {0.45};
    double weights_expected[n_nodes*n_pts] = {0.055, 0.99, -0.045}; 
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_BARHO;
    const int order = 3;
    const int n_nodes = order + 1;
    double uvw[1] = {0.45};
    double weights_expected[n_nodes*n_pts] = {-1001./16000., 11583./16000., 6237./16000., -819./16000.};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 4") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_BARHO;
    const int order = 4;
    const int n_nodes = order + 1;
    double uvw[1] = {0.45};
    double weights_expected[n_nodes*n_pts] = {-11./625., 99./625., 594./625., -66./625., 9./625.};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 5") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_BARHO;
    const int order = 5;
    const int n_nodes = order + 1;
    double uvw[1] = {0.45};
    double weights_expected[n_nodes*n_pts] = {77./8192., -693/8192., 3465./4096., 1155./4096., -495./8192., 63./8192.};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] tria",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TRIA3;
    const int order = 1;
    const int n_nodes = (order + 2) * (order + 1) / 2;
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = {0.2, 0.45, 
                                              0.35};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TRIAHO;
    const int order = 2;
    const int n_nodes = (order + 2) * (order + 1) / 2;
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = { -3./25.,  9./25., -9./200., 
                                                7./25., 63./100., 
                                              -21./200.};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TRIAHO;
    const int order = 3;
    const int n_nodes = (order + 2) * (order + 1) / 2;
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = {   7./125.  ,  -81./500.  ,  567./4000. , -819./16000., 
                                               -63./500.  , 1701./2000. , 3969./16000.,
                                                63./4000. ,  567./16000., 
                                              -133./16000.};
    weights = (double*) malloc(sizeof(double) * n_nodes*n_pts);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 4") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TRIAHO;
    const int order = 4;
    const int n_nodes = (order + 2) * (order + 1) / 2;
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = {-0.0176,  0.0576, -0.0576, -0.0384,  0.0144, 
                                               0.0448, -0.2016,  0.8064, -0.0672, 
                                              -0.0224,  0.4032,  0.2016, 
                                              -0.0448, -0.1008, 
                                               0.0224};
    weights = (double*) malloc(sizeof(double) * n_nodes*n_pts);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] quad",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_QUAD4;
    const int order = 1;
    const int n_nodes = (order + 1) * (order + 1);
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = {0.3575, 0.2925, 
                                              0.1575, 0.1925};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_QUADHO;
    const int order = 2;
    const int n_nodes = (order + 1) * (order + 1);
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = { 0.010725,  0.193050, -0.008775, 
                                               0.050050,  0.900900, -0.040950, 
                                              -0.005775, -0.103950,  0.004725};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        // printf("weights : %.12f\n", weights[n_nodes*p+n]);
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_QUADHO;
    const int order = 3;
    const int n_nodes = (order + 1) * (order + 1);
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = { 0.0009658085937500, -0.0111757851562500, -0.0060177304687500,  0.0007902070312500, 
                                              -0.0608459414062500,  0.7040744648437497,  0.3791170195312500, -0.0497830429687500, 
                                              -0.0032024179687500,  0.0370565507812500,  0.0199535273437500, -0.0026201601562500, 
                                               0.0005200507812500, -0.0060177304687500, -0.0032403164062500,  0.0004254960937500};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 4") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_QUADHO;
    const int order = 4;
    const int n_nodes = (order + 1) * (order + 1);
    double uvw[2] = {0.45, 0.35};
    double weights_expected[n_nodes*n_pts] = { 0.00073216, -0.00658944, -0.03953664,  0.00439296, -0.00059904, 
                                              -0.01025024,  0.09225216,  0.55351296, -0.06150144,  0.00838656, 
                                              -0.01025024,  0.09225216,  0.55351296, -0.06150144,  0.00838656, 
                                               0.00256256, -0.02306304, -0.13837824,  0.01537536, -0.00209664, 
                                              -0.00039424,  0.00354816,  0.02128896, -0.00236544,  0.00032256};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] tetra",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TETRA4;
    const int order = 1;
    const int n_nodes = (order+1) * (order+2) * (order+3) / 6;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {0.05, 0.45, 0.35, 
                                              0.15};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TETRAHO;
    const int order = 2;
    const int n_nodes = (order+1) * (order+2) * (order+3) / 6;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {-0.045, 0.090, -0.045, 0.070, 0.630, -0.105, 
                                               0.030, 0.270,  0.210, 
                                              -0.105};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_TETRAHO;
    const int order = 3;
    const int n_nodes = (order+1) * (order+2) * (order+3) / 6;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = { 0.0393125, -0.0860625,  0.0354375, -0.0511875, -0.0669375, 0.2126250, 0.2480625, 0.0039375, 0.0354375, -0.0083125, 
                                              -0.0286875,  0.0911250,  0.1063125,  0.0708750,  0.6378750, 0.0118125, 
                                              -0.0185625, -0.1670625, -0.1299375, 
                                               0.0639375};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] pyra",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_PYRAMID5;
    const int order = 1;
    const int n_nodes = 5;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {0.2352941176470588, 0.2647058823529412, 0.1852941176470588, 0.1647058823529412, 
                                              0.1500000000000000};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_PYRAMIDHO;
    const int order = 2;
    const int n_nodes = 14;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {-0.0373702422145329, 0.0747404844290658, -0.0373702422145329, -0.0193771626297578, 0.6975778546712803, 0.0217993079584775, -0.0232525951557093, -0.0523183391003460, -0.0294290657439446, 
                                               0.1411764705882353, 0.1588235294117647,  0.0988235294117647,  0.1111764705882353, 
                                              -0.1050000000000000};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] prism",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_PRISM6;
    const int order = 1;
    const int n_nodes = (order+1) * (order+1) * (order+2) /2;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {0.1700, 0.3825, 0.2975, 
                                              0.0300, 0.0675, 0.0525};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_PRISMHO;
    const int order = 2;
    const int n_nodes = (order+1) * (order+1) * (order+2) /2;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {-0.071400,  0.214200, -0.026775,  0.166600,  0.374850, -0.062475, 
                                              -0.061200,  0.183600, -0.022950,  0.142800,  0.321300, -0.053550,
                                               0.012600, -0.037800,  0.004725, -0.029400, -0.066150,  0.011025};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_PRISMHO;
    const int order = 3;
    const int n_nodes = (order+1) * (order+1) * (order+2) /2;
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = { 0.02028950000000, -0.05869462500000,  0.05135779687500, -0.01854587109375, -0.04565137500000,  0.30814678125000,  0.08987614453125,  0.00570642187500,  0.01283944921875, -0.00301172265625, 
                                               0.04980150000000, -0.14406862500000,  0.12606004687500, -0.04552168359375, -0.11205337500000,  0.75636028125000,  0.22060508203125,  0.01400667187500,  0.03151501171875, -0.00739241015625, 
                                              -0.01767150000000,  0.05112112500000, -0.04473098437500,  0.01615285546875,  0.03976087500000, -0.26838590625000, -0.07827922265625, -0.00497010937500, -0.01118274609375,  0.00262311328125, 
                                               0.00358050000000, -0.01035787500000,  0.00906314062500, -0.00327280078125, -0.00805612500000,  0.05437884375000,  0.01586049609375,  0.00100701562500,  0.00226578515625, -0.00053148046875};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}

MPI_TEST_CASE("[PDM_ho_basis] hexa",1) {
  const int n_pts = 1;
  double *weights = NULL;

  SUBCASE("order 1") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_HEXA8;
    const int order = 1;
    const int n_nodes = (order+1) * (order+1) * (order+1);
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = {0.303875, 0.248625, 0.133875, 0.163625, 
                                              0.053625, 0.043875, 0.023625, 0.028875};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 2") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_HEXAHO;
    const int order = 2;
    const int n_nodes = (order+1) * (order+1) * (order+1);
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = { 0.006381375,  0.114864750, -0.005221125,  0.029779750,  0.536035500, -0.024365250, -0.003436125, -0.061850250,  0.002811375, 
                                               0.005469750,  0.098455500, -0.004475250,  0.025525500,  0.459459000, -0.020884500, -0.002945250, -0.053014500,  0.002409750, 
                                              -0.001126125, -0.020270250,  0.000921375, -0.005255250, -0.094594500,  0.004299750,  0.000606375,  0.010914750, -0.000496125};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }

  SUBCASE("order 3") {
    const PDM_Mesh_nodal_elt_t  type = PDM_MESH_NODAL_HEXAHO;
    const int order = 3;
    const int n_nodes = (order+1) * (order+1) * (order+1);
    double uvw[3] = {0.45, 0.35, 0.15};
    double weights_expected[n_nodes*n_pts] = { 0.0003499245261230, -0.0040491266594238, -0.0021802989704590,  0.0002863018850098, -0.0220452451457519,  0.2550949795437009,  0.1373588351389159, -0.0180370187556152, -0.0011602760603027,  0.0134260515549316,  0.0072294123757324, -0.0009493167766113,  0.0001884208986816, -0.0021802989704590, -0.0011740071379395,  0.0001541625534668, 
                                               0.0008589056550293, -0.0099387654367676, -0.0053516429274902,  0.0007027409904785, -0.0541110562668457,  0.6261422225163570,  0.3371535044318847, -0.0442726824001465, -0.0028479503298340,  0.0329548538166504,  0.0177449212858887, -0.0023301411789551,  0.0004624876604004, -0.0053516429274902, -0.0028816538840332,  0.0003783989948730, 
                                              -0.0003047729743652,  0.0035266587033691,  0.0018989700710449, -0.0002493597062988,  0.0192006973850098, -0.2221794983122557, -0.1196351144758300,  0.0157096614968262,  0.0010105630202637, -0.0116936578059082, -0.0062965849724121,  0.0008268242893066, -0.0001641085246582,  0.0018989700710449,  0.0010225223459473, -0.0001342706110840, 
                                               0.0000617513869629, -0.0007145517634277, -0.0003847586418457,  0.0000505238620605, -0.0038903373786621,  0.0450167610959472,  0.0242397944362793, -0.0031830033098145, -0.0002047545988770,  0.0023693032155762, 0.0012757786545410,  -0.0001675264899902,  0.0000332507468262, -0.0003847586418457, -0.0002071777302246,  0.0000272051564941};
    PDM_malloc(weights, n_nodes*n_pts, double);

    PDM_ho_basis(type, order, n_nodes, n_pts, uvw, weights);

    for (int p=0; p<n_pts; p++){
      double weights_summed = 0;
      for (int n=0; n<n_nodes; n++){
        CHECK(fabs(weights[n_nodes*p+n] - weights_expected[n_nodes*p+n]) < tol);
        weights_summed +=weights[n_nodes*p+n];
      }
      CHECK(fabs(weights_summed - 1.) < tol);
    }
  }
  PDM_free(weights);
}