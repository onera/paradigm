#include "doctest/extensions/doctest_mpi.h"
#include "pdm_doctest.h"
#include "float.h"
#include "pdm.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "math.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"

#include "pdm_mesh_location.h"


static double tol = 1e-10;

MPI_TEST_CASE("[pdm_mesh_location] - 2D", 1) {
  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_l_num_t n_vtx;
  double *vtx_coord=NULL;
  double fact;

  PDM_l_num_t n_face;
  PDM_l_num_t *face_vtx_idx = NULL;
  PDM_l_num_t *face_vtx = NULL;

  PDM_g_num_t *face_ln_to_gn = NULL;
  PDM_g_num_t *vtx_ln_to_gn = NULL;

  PDM_l_num_t n_pts;
  double *pts_coord = NULL;
  PDM_g_num_t *gnum = NULL;
  double *expected_weights = NULL;

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_method_set(ml, PDM_MESH_LOCATION_LOCATE_ALL_TGT);
  PDM_mesh_location_tolerance_set(ml, 1e-16);




  SUBCASE("tria") {

    printf("TRIA\n");
    n_vtx = 3;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-6; 

    // 1er point
    vtx_coord[ 0] = 0.0;
    vtx_coord[ 1] = 0.0;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] = 1.0*fact;
    vtx_coord[ 4] = 0.0;
    vtx_coord[ 5] = 0.0;
    // 2e point
    vtx_coord[ 6] = 1.0;
    vtx_coord[ 7] = 1.0;
    vtx_coord[ 8] = 0.0;


    n_face = 1;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0] = 0;
    face_vtx_idx[1] = 3;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[0] = 1;
    face_vtx[1] = 2;
    face_vtx[2] = 3;

    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    face_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0]  = 1;
    vtx_ln_to_gn[1]  = 2;
    vtx_ln_to_gn[2]  = 3;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = fact/2. + 0.5*(1.-fact/2.);
    pts_coord[1] = 0.5;
    pts_coord[2] = 0.0;
    
    PDM_malloc(expected_weights, n_vtx*n_pts, double);
    expected_weights[0] = 0.25;
    expected_weights[1] = 0.25;
    expected_weights[2] = 0.5;

  }

  SUBCASE("quad") {
    printf("QUAD\n");
    n_vtx = 4;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-5;  

    // 1er point
    vtx_coord[ 0] = 0.0;
    vtx_coord[ 1] = 0.0;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] = 0.5*(1.0+fact);
    vtx_coord[ 4] = 0.5*(1.0-fact);
    vtx_coord[ 5] = 0.0;
    // 3e point
    vtx_coord[ 6] = 1.0;
    vtx_coord[ 7] = 1.0;
    vtx_coord[ 8] = 0.0;
    // 4e point
    vtx_coord[ 9] = 0.5*(1.0-fact);
    vtx_coord[10] = 0.5*(1.0+fact);
    vtx_coord[11] = 0.0;


    n_face = 1;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0] = 0;
    face_vtx_idx[1] = 4;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[0] = 1;
    face_vtx[1] = 2;
    face_vtx[2] = 3;
    face_vtx[3] = 4;

    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    face_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0]  = 1;
    vtx_ln_to_gn[1]  = 2;
    vtx_ln_to_gn[2]  = 3;
    vtx_ln_to_gn[3]  = 4;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.65 + 0.5*fact;
    pts_coord[1] = 0.65;
    pts_coord[2] = 0.0;

    PDM_malloc(expected_weights, n_vtx*n_pts, double);
    expected_weights[0] = 5.9998250e-02;
    expected_weights[1] = 5.3999925e-01;
    expected_weights[2] = 3.6000325e-01;
    expected_weights[3] = 3.9999250e-02;    

  }

  SUBCASE("polygon convexe") {
    printf("POLY 1\n");
    n_vtx = 5;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-4;

    // 1er point
    vtx_coord[ 0] = 1.0;
    vtx_coord[ 1] = 0.0;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] = (cos(72.*PDM_PI/180.)-1.0)*fact + 1.0;
    vtx_coord[ 4] = sin(72.*PDM_PI/180.)*fact;
    vtx_coord[ 5] = 0.0;
    // 3e point
    vtx_coord[ 6] = cos(2*72.*PDM_PI/180.);
    vtx_coord[ 7] = sin(2*72.*PDM_PI/180.)*fact;
    vtx_coord[ 8] = 0.0;
    // 4e point
    vtx_coord[ 9] = cos(3*72.*PDM_PI/180.);
    vtx_coord[10] = sin(3*72.*PDM_PI/180.)*fact;
    vtx_coord[11] = 0.0;
    // 5e point
    vtx_coord[12] = (cos(4*72.*PDM_PI/180.)-1.0)*fact + 1.0;
    vtx_coord[13] = sin(4*72.*PDM_PI/180.)*fact;
    vtx_coord[14] = 0.0;


    n_face = 1;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0] = 0;
    face_vtx_idx[1] = 5;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[0] = 1;
    face_vtx[1] = 2;
    face_vtx[2] = 3;
    face_vtx[3] = 4;
    face_vtx[4] = 5;

    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    face_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;


    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 1.0-fact;
    pts_coord[1] = 0.1*fact;
    pts_coord[2] = 0.0;

    PDM_malloc(expected_weights, n_vtx*n_pts, double);
    expected_weights[0] = 2.527255770e-01;
    expected_weights[1] = 4.261966835e-01;
    expected_weights[2] = 1.380675184e-05;
    expected_weights[3] = 1.292956678e-05;
    expected_weights[4] = 3.210510032e-01;




  }
  SUBCASE("polygon concave"){
    printf("POLY 2\n");
    n_vtx = 8;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-6;

    // 1er point
    vtx_coord[ 0] = -1.0;
    vtx_coord[ 1] =  0.0;
    vtx_coord[ 2] =  0.0;
    // 2e point
    vtx_coord[ 3] = -fact;
    vtx_coord[ 4] =  0.0;
    vtx_coord[ 5] =  0.0;
    // 3e point
    vtx_coord[ 6] = -fact;
    vtx_coord[ 7] =  0.5;
    vtx_coord[ 8] =  0.0;
    // 4e point
    vtx_coord[ 9] =  fact;
    vtx_coord[10] =  0.5;
    vtx_coord[11] =  0.0;
    // 5e point
    vtx_coord[12] =  fact;
    vtx_coord[13] =  0.0;
    vtx_coord[14] =  0.0;
    // 6e point
    vtx_coord[15] =  1.0;
    vtx_coord[16] =  0.0;
    vtx_coord[17] =  0.0;
    // 7e point
    vtx_coord[18] =  1.0;
    vtx_coord[19] =  1.0;
    vtx_coord[20] =  0.0;
    // 8e point
    vtx_coord[21] = -1.0;
    vtx_coord[22] =  1.0;
    vtx_coord[23] =  0.0;

    n_face = 1;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0] = 0;
    face_vtx_idx[1] = 8;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[0] = 1;
    face_vtx[1] = 2;
    face_vtx[2] = 3;
    face_vtx[3] = 4;
    face_vtx[4] = 5;
    face_vtx[5] = 6;
    face_vtx[6] = 7;
    face_vtx[7] = 8;

    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    face_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;
    vtx_ln_to_gn[6] = 7;
    vtx_ln_to_gn[7] = 8;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    for (int i = 0; i<n_pts; i++){
      gnum[i] = i+1;
      pts_coord[3*i]   = -5.0*fact;
      pts_coord[3*i+1] =  0.4;
      pts_coord[3*i+2] =  0.0;
    }

    PDM_malloc(expected_weights, n_vtx*n_pts, double);
    expected_weights[0] =  6.5346563819e-06;
    expected_weights[1] =  5.9997677649e-01;
    expected_weights[2] =  2.3998664835e+00;
    expected_weights[3] = -1.5999109888e+00;
    expected_weights[4] = -3.9996759295e-01;
    expected_weights[5] =  6.5344852282e-06;
    expected_weights[6] =  1.1126286778e-05;
    expected_weights[7] =  1.1126393782e-05;
  }

  PDM_mesh_location_mesh_n_part_set(ml, 1);


  PDM_mesh_location_nodal_part_set_2d(ml,
                                      0,
                                      n_face,
                                      face_vtx_idx,
                                      face_vtx,
                                      face_ln_to_gn,
                                      n_vtx,
                                      vtx_coord,
                                      vtx_ln_to_gn);


  PDM_mesh_location_n_part_cloud_set(ml, 0, 1);

  PDM_mesh_location_cloud_set(ml,
                              0,
                              0,
                              n_pts,
                              pts_coord,
                              gnum);

  PDM_mesh_location_compute(ml);

  int *elt_pts_inside_idx = NULL;
  PDM_g_num_t *points_gnum = NULL;
  double *points_coords = NULL;
  double *points_uvw = NULL;
  int *points_weights_idx = NULL;
  double *points_weights = NULL;
  double *points_dist2 = NULL;
  double *points_projected_coords = NULL;

  PDM_mesh_location_points_in_elt_get(ml,
                                      0,
                                      0,
                                      &elt_pts_inside_idx,
                                      &points_gnum,
                                      &points_coords,
                                      &points_uvw,
                                      &points_weights_idx,
                                      &points_weights,
                                      &points_dist2,
                                      &points_projected_coords);



  int located = PDM_mesh_location_n_located_get(ml, 0, 0);

  for (int i=0; i<n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<3*n_pts; i++){
    CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  }


  PDM_mesh_location_free(ml);
  PDM_free(vtx_coord);
  PDM_free(face_vtx_idx);
  PDM_free(face_vtx);
  PDM_free(face_ln_to_gn);
  PDM_free(vtx_ln_to_gn);
  PDM_free(pts_coord);
  PDM_free(gnum);
  PDM_free(expected_weights);
}


MPI_TEST_CASE("[pdm_mesh_location] - 3D nodal", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_l_num_t n_vtx;
  double *vtx_coord=NULL;
  double fact;

  PDM_l_num_t n_cell;
  PDM_l_num_t *cell_vtx_idx = NULL;
  PDM_l_num_t *cell_vtx = NULL;

  PDM_g_num_t *cell_ln_to_gn = NULL;
  PDM_g_num_t *vtx_ln_to_gn = NULL;

  PDM_l_num_t n_pts;
  double *pts_coord = NULL;
  PDM_g_num_t *gnum = NULL;
  double *expected_weights = NULL;

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_method_set(ml, PDM_MESH_LOCATION_LOCATE_ALL_TGT);
  PDM_mesh_location_tolerance_set(ml, 1e-16);


  SUBCASE("tetra"){
    printf("TETRA\n");
    n_vtx = 4;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-2;

    // 1er point
    vtx_coord[ 0] = 0.0;
    vtx_coord[ 1] = -1.0*fact;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] =  2.0;
    vtx_coord[ 4] =  0.0;
    vtx_coord[ 5] = -1.0*fact;
    // 3e point
    vtx_coord[ 6] = 0.0;
    vtx_coord[ 7] = 1.0*fact;
    vtx_coord[ 8] = 0.0;
    // 4e point
    vtx_coord[ 9] = 2.0;
    vtx_coord[10] = 0.0;
    vtx_coord[11] = 1.0*fact;

    n_cell = 1;
    PDM_malloc(cell_vtx_idx, n_cell+1, PDM_l_num_t);
    cell_vtx_idx[0] = 0;
    cell_vtx_idx[1] = 4;
    PDM_malloc(cell_vtx, cell_vtx_idx[n_cell], PDM_l_num_t);
    cell_vtx[0] = 1;
    cell_vtx[1] = 2;
    cell_vtx[2] = 3;
    cell_vtx[3] = 4;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;

    n_pts = 2;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);
    gnum[0] = 1;
    gnum[1] = 2;

    pts_coord[0] = 0.45;
    pts_coord[1] = 0.05*fact;
    pts_coord[2] = -0.03*fact;

    pts_coord[3] = 0.75;
    pts_coord[4] = 0.25*fact;
    pts_coord[5] = 1.0*fact;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 0.3625;
    expected_weights[1] = 0.1275;
    expected_weights[2] = 0.4125;
    expected_weights[3] = 0.0975;
    expected_weights[4] = 0.18749218769507753;
    expected_weights[5] = 0.0;
    expected_weights[6] = 0.43749218769563264;
    expected_weights[7] = 0.37501562460928983;

  }

  SUBCASE("hexa"){
    printf("HEXA\n");
    n_vtx = 8;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0;  

    // 1er point
    vtx_coord[ 0] =  0.0;
    vtx_coord[ 1] = -1.0*fact;
    vtx_coord[ 2] = -1.0*fact;
    // 2e point
    vtx_coord[ 3] =  1.0;
    vtx_coord[ 4] =  0.0*fact;
    vtx_coord[ 5] = -1.0*fact;
    // 3e point
    vtx_coord[ 6] = -1.0;
    vtx_coord[ 7] =  0.0*fact;
    vtx_coord[ 8] = -1.0*fact;
    // 4e point
    vtx_coord[ 9] =  0.0;
    vtx_coord[10] =  1.0*fact;
    vtx_coord[11] = -1.0*fact;
    // 5e point
    vtx_coord[12] =  0.0;
    vtx_coord[13] = -1.0*fact;
    vtx_coord[14] =  1.0*fact;
    // 6e point
    vtx_coord[15] = 1.0;
    vtx_coord[16] = 0.0*fact;
    vtx_coord[17] = 1.0*fact;
    // 7e point
    vtx_coord[18] = -1.0;
    vtx_coord[19] = 0.0*fact;
    vtx_coord[20] = 1.0*fact;
    // 8e point
    vtx_coord[21] = 0.0;
    vtx_coord[22] =  1.0*fact;
    vtx_coord[23] =  1.0*fact;

    n_cell = 1;
    PDM_malloc(cell_vtx_idx, n_cell+1, PDM_l_num_t);
    cell_vtx_idx[0] = 0;
    cell_vtx_idx[1] = 8;
    PDM_malloc(cell_vtx, cell_vtx_idx[n_cell], PDM_l_num_t);
    cell_vtx[0] = 1;
    cell_vtx[1] = 2;
    cell_vtx[2] = 4;
    cell_vtx[3] = 3;
    cell_vtx[4] = 5;
    cell_vtx[5] = 6;
    cell_vtx[6] = 8;
    cell_vtx[7] = 7;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;
    vtx_ln_to_gn[6] = 7;
    vtx_ln_to_gn[7] = 8;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.35;
    pts_coord[1] = 0.45;
    pts_coord[2] = 0.15;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 1.91250e-02;
    expected_weights[1] = 1.72125e-01;
    expected_weights[2] = 2.10375e-01;
    expected_weights[3] = 2.33750e-02;
    expected_weights[4] = 2.58750e-02;
    expected_weights[5] = 2.32875e-01;
    expected_weights[6] = 2.84625e-01;
    expected_weights[7] = 3.16250e-02;

  }


  SUBCASE("pyra"){
    printf("PYRA\n");

    n_vtx = 5;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0e-11;

    // 1er point
    vtx_coord[ 0] = 0.0;
    vtx_coord[ 1] = 0.0;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] = 1.0;
    vtx_coord[ 4] = 0.0;
    vtx_coord[ 5] = 0.0;
    // 3e point
    vtx_coord[ 6] = 0.0;
    vtx_coord[ 7] = 1.0;
    vtx_coord[ 8] = -fact;
    // 4e point
    vtx_coord[ 9] = 1.0;
    vtx_coord[10] = 1.0;
    vtx_coord[11] = -fact;
    // 5e point
    vtx_coord[12] = 0.5;
    vtx_coord[13] = 0.5;
    vtx_coord[14] = fact;

    n_cell = 1;
    PDM_malloc(cell_vtx_idx, n_cell+1, PDM_l_num_t);
    cell_vtx_idx[0] = 0;
    cell_vtx_idx[1] = 5;
    PDM_malloc(cell_vtx, cell_vtx_idx[n_cell], PDM_l_num_t);
    cell_vtx[0] = 1;
    cell_vtx[1] = 2;
    cell_vtx[2] = 4;
    cell_vtx[3] = 3;
    cell_vtx[4] = 5;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 4;
    vtx_ln_to_gn[3] = 3;
    vtx_ln_to_gn[4] = 5;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.45 ;
    pts_coord[1] = 0.35;
    pts_coord[2] = 0.05*fact;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 0.293560606;
    expected_weights[1] = 0.2231060606;
    expected_weights[2] = 0.093560606;
    expected_weights[3] = 0.1231060606;
    expected_weights[4] = 0.2666666666;

  }


  SUBCASE("prism"){
    printf("PRISME\n");

    n_vtx = 6;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    // double fact = 1.0e-1;

    // 1er point
    vtx_coord[ 0] = 0.0;
    vtx_coord[ 1] = 0.0;
    vtx_coord[ 2] = 0.0;
    // 2e point
    vtx_coord[ 3] = 1.0;
    vtx_coord[ 4] = 0.0;
    vtx_coord[ 5] = 0.0;
    // 3e point
    vtx_coord[ 6] = 0.0;
    vtx_coord[ 7] = 1.0;
    vtx_coord[ 8] = 0.0;
    // 4e point
    vtx_coord[ 9] = 0.0;
    vtx_coord[10] = 0.0;
    vtx_coord[11] = 1.0;
    // 5e point
    vtx_coord[12] = 1.0;
    vtx_coord[13] = 0.0;
    vtx_coord[14] = 1.0;
    // 6e point
    vtx_coord[15] = 0.0;
    vtx_coord[16] = 1.0;
    vtx_coord[17] = 1.0;


    n_cell = 1;
    PDM_malloc(cell_vtx_idx, n_cell+1, PDM_l_num_t);
    cell_vtx_idx[0] = 0;
    cell_vtx_idx[1] = 6;
    PDM_malloc(cell_vtx, cell_vtx_idx[n_cell], PDM_l_num_t);
    cell_vtx[0] = 1;
    cell_vtx[1] = 2;
    cell_vtx[2] = 3;
    cell_vtx[3] = 4;
    cell_vtx[4] = 5;
    cell_vtx[5] = 6;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 4;
    vtx_ln_to_gn[3] = 3;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.45;
    pts_coord[1] = 0.35;
    pts_coord[2] = 0.15 ;


    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 0.17;
    expected_weights[1] = 0.3825;
    expected_weights[2] = 0.2975;
    expected_weights[3] = 0.03;
    expected_weights[4] = 0.0675;
    expected_weights[5] = 0.0525;

  }

  PDM_mesh_location_mesh_n_part_set(ml, 1);


  PDM_mesh_location_nodal_part_set(ml,
                                      0,
                                      n_cell,
                                      cell_vtx_idx,
                                      cell_vtx,
                                      cell_ln_to_gn,
                                      n_vtx,
                                      vtx_coord,
                                      vtx_ln_to_gn);


  PDM_mesh_location_n_part_cloud_set(ml, 0, 1);

  PDM_mesh_location_cloud_set(ml,
                              0,
                              0,
                              n_pts,
                              pts_coord,
                              gnum);

  PDM_mesh_location_compute(ml);

  int *elt_pts_inside_idx = NULL;
  PDM_g_num_t *points_gnum = NULL;
  double *points_coords = NULL;
  double *points_uvw = NULL;
  int *points_weights_idx = NULL;
  double *points_weights = NULL;
  double *points_dist2 = NULL;
  double *points_projected_coords = NULL;

  PDM_mesh_location_points_in_elt_get(ml,
                                      0,
                                      0,
                                      &elt_pts_inside_idx,
                                      &points_gnum,
                                      &points_coords,
                                      &points_uvw,
                                      &points_weights_idx,
                                      &points_weights,
                                      &points_dist2,
                                      &points_projected_coords);



  int located = PDM_mesh_location_n_located_get(ml, 0, 0);

  for (int i=0; i<n_pts*n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<3*n_pts; i++){
    CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  }


  PDM_mesh_location_free(ml);
  PDM_free(vtx_coord);
  PDM_free(cell_vtx_idx);
  PDM_free(cell_vtx);
  PDM_free(cell_ln_to_gn);
  PDM_free(vtx_ln_to_gn);
  PDM_free(pts_coord);
  PDM_free(gnum);
  PDM_free(expected_weights);

}






MPI_TEST_CASE("[pdm_mesh_location] - 3D", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_l_num_t n_vtx;
  double *vtx_coord=NULL;
  double fact;

  PDM_l_num_t n_cell;
  PDM_l_num_t *cell_face_idx = NULL;
  PDM_l_num_t *cell_face = NULL;

  PDM_l_num_t n_face;
  PDM_l_num_t *face_vtx_idx = NULL;
  PDM_l_num_t *face_vtx = NULL;

  PDM_g_num_t *cell_ln_to_gn = NULL;
  PDM_g_num_t *face_ln_to_gn = NULL;
  PDM_g_num_t *vtx_ln_to_gn = NULL;

  PDM_l_num_t n_pts;
  double *pts_coord = NULL;
  PDM_g_num_t *gnum = NULL;
  double *expected_weights = NULL;

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_mesh_location_method_set(ml, PDM_MESH_LOCATION_LOCATE_ALL_TGT);
  PDM_mesh_location_tolerance_set(ml, 1e-16);

  SUBCASE("hexa"){
    printf("HEXA\n");
    n_vtx = 8;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0;  

    // 1er point
    vtx_coord[ 0] =  0.0;
    vtx_coord[ 1] = -1.0*fact;
    vtx_coord[ 2] = -1.0*fact;
    // 2e point
    vtx_coord[ 3] =  1.0;
    vtx_coord[ 4] =  0.0*fact;
    vtx_coord[ 5] = -1.0*fact;
    // 3e point
    vtx_coord[ 6] = -1.0;
    vtx_coord[ 7] =  0.0*fact;
    vtx_coord[ 8] = -1.0*fact;
    // 4e point
    vtx_coord[ 9] =  0.0;
    vtx_coord[10] =  1.0*fact;
    vtx_coord[11] = -1.0*fact;
    // 5e point
    vtx_coord[12] =  0.0;
    vtx_coord[13] = -1.0*fact;
    vtx_coord[14] =  1.0*fact;
    // 6e point
    vtx_coord[15] = 1.0;
    vtx_coord[16] = 0.0*fact;
    vtx_coord[17] = 1.0*fact;
    // 7e point
    vtx_coord[18] = -1.0;
    vtx_coord[19] = 0.0*fact;
    vtx_coord[20] = 1.0*fact;
    // 8e point
    vtx_coord[21] = 0.0;
    vtx_coord[22] =  1.0*fact;
    vtx_coord[23] =  1.0*fact;

    n_cell = 1;
    PDM_malloc(cell_face_idx, n_cell+1, PDM_l_num_t);
    cell_face_idx[0] = 0;
    cell_face_idx[1] = 6;
    PDM_malloc(cell_face, cell_face_idx[n_cell], PDM_l_num_t);
    cell_face[0] = 1;
    cell_face[1] = 2;
    cell_face[2] = 3;
    cell_face[3] = 4;
    cell_face[4] = 5;
    cell_face[5] = 6;


    n_face = 6;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0] =  0;
    face_vtx_idx[1] =  4;
    face_vtx_idx[2] =  8;
    face_vtx_idx[3] = 12;
    face_vtx_idx[4] = 16;
    face_vtx_idx[5] = 20;
    face_vtx_idx[6] = 24;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[ 0] = 1;
    face_vtx[ 1] = 2;
    face_vtx[ 2] = 4;
    face_vtx[ 3] = 3;
    face_vtx[ 4] = 1;
    face_vtx[ 5] = 5;
    face_vtx[ 6] = 6;
    face_vtx[ 7] = 2;
    face_vtx[ 8] = 2;
    face_vtx[ 9] = 6;
    face_vtx[10] = 8;
    face_vtx[11] = 4;
    face_vtx[12] = 4;
    face_vtx[13] = 8;
    face_vtx[14] = 7;
    face_vtx[15] = 3;
    face_vtx[16] = 3;
    face_vtx[17] = 7;
    face_vtx[18] = 1;
    face_vtx[19] = 5;
    face_vtx[20] = 5;
    face_vtx[21] = 7;
    face_vtx[22] = 8;
    face_vtx[23] = 6;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    face_ln_to_gn[0] = 1;
    face_ln_to_gn[1] = 2;
    face_ln_to_gn[2] = 3;
    face_ln_to_gn[3] = 4;
    face_ln_to_gn[4] = 5;
    face_ln_to_gn[5] = 6;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;
    vtx_ln_to_gn[6] = 7;
    vtx_ln_to_gn[7] = 8;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.35;
    pts_coord[1] = 0.45;
    pts_coord[2] = 0.15;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 1.91250e-02;
    expected_weights[1] = 1.72125e-01;
    expected_weights[2] = 2.10375e-01;
    expected_weights[3] = 2.33750e-02;
    expected_weights[4] = 2.58750e-02;
    expected_weights[5] = 2.32875e-01;
    expected_weights[6] = 2.84625e-01;
    expected_weights[7] = 3.16250e-02;

  }

  SUBCASE("polyedre convexe"){
    printf("POLY 1\n");

    n_vtx = 8;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0;  

    // 1er point
    vtx_coord[ 0] = 0.0*fact;
    vtx_coord[ 1] = 0.0*fact;
    vtx_coord[ 2] = 0.0*fact;
    // 2e point
    vtx_coord[ 3] = 1.0*fact;
    vtx_coord[ 4] = 0.0*fact;
    vtx_coord[ 5] = 0.0*fact;
    // 3e point
    vtx_coord[ 6] = 1.0*fact;
    vtx_coord[ 7] = 1.0*fact;
    vtx_coord[ 8] = 0.0*fact;
    // 4e point
    vtx_coord[ 9] = 0.0*fact;
    vtx_coord[10] = 1.0*fact;
    vtx_coord[11] = 0.0*fact;
    // 5e point
    vtx_coord[12] = 0.5*fact;
    vtx_coord[13] = 0.0*fact;
    vtx_coord[14] = 1.0*fact;
    // 6e point
    vtx_coord[15] = 1.0*fact;
    vtx_coord[16] = 0.5*fact;
    vtx_coord[17] = 1.0*fact;
    // 7e point
    vtx_coord[18] = 0.5*fact;
    vtx_coord[19] = 1.0*fact;
    vtx_coord[20] = 1.0*fact;
    // 8e point
    vtx_coord[21] = 0.0*fact;
    vtx_coord[22] = 0.5*fact;
    vtx_coord[23] = 1.0*fact;

    n_cell = 1;
    PDM_malloc(cell_face_idx, n_cell+1, PDM_l_num_t);
    cell_face_idx[0] = 0;
    cell_face_idx[1] = 10;
    PDM_malloc(cell_face, cell_face_idx[n_cell], PDM_l_num_t);
    cell_face[0] = 1;
    cell_face[1] = 2;
    cell_face[2] = 3;
    cell_face[3] = 4;
    cell_face[4] = 5;
    cell_face[5] = 6;
    cell_face[6] = 7;
    cell_face[7] = 8;
    cell_face[8] = 9;
    cell_face[9] = 10;


    n_face = 10;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0]  =  0;
    face_vtx_idx[1]  =  4;
    face_vtx_idx[2]  =  8;
    face_vtx_idx[3]  = 11;
    face_vtx_idx[4]  = 14;
    face_vtx_idx[5]  = 17;
    face_vtx_idx[6]  = 20;
    face_vtx_idx[7]  = 23;
    face_vtx_idx[8]  = 26;
    face_vtx_idx[9]  = 29;
    face_vtx_idx[10] = 32;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[ 0] = 1;
    face_vtx[ 1] = 4;
    face_vtx[ 2] = 3;
    face_vtx[ 3] = 2;
    face_vtx[ 4] = 5;
    face_vtx[ 5] = 6;
    face_vtx[ 6] = 7;
    face_vtx[ 7] = 8;
    face_vtx[ 8] = 1;
    face_vtx[ 9] = 2;
    face_vtx[10] = 5;
    face_vtx[11] = 2;
    face_vtx[12] = 3;
    face_vtx[13] = 6;
    face_vtx[14] = 3;
    face_vtx[15] = 4;
    face_vtx[16] = 7;
    face_vtx[17] = 4;
    face_vtx[18] = 1;
    face_vtx[19] = 8;
    face_vtx[20] = 1;
    face_vtx[21] = 5;
    face_vtx[22] = 8;
    face_vtx[23] = 2;
    face_vtx[24] = 6;
    face_vtx[25] = 5;
    face_vtx[26] = 3;
    face_vtx[27] = 7;
    face_vtx[28] = 6;
    face_vtx[29] = 4;
    face_vtx[30] = 8;
    face_vtx[31] = 7;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    face_ln_to_gn[0] = 1;
    face_ln_to_gn[1] = 2;
    face_ln_to_gn[2] = 3;
    face_ln_to_gn[3] = 4;
    face_ln_to_gn[4] = 5;
    face_ln_to_gn[5] = 6;
    face_ln_to_gn[6] = 7;
    face_ln_to_gn[7] = 8;
    face_ln_to_gn[8] = 9;
    face_ln_to_gn[9] = 10;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;
    vtx_ln_to_gn[6] = 7;
    vtx_ln_to_gn[7] = 8;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.35;
    pts_coord[1] = 0.25;
    pts_coord[2] = 0.9;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 5.852188877e-02;
    expected_weights[1] = 1.846817191e-02;
    expected_weights[2] = 9.851165651e-03;
    expected_weights[3] = 1.315877366e-02;
    expected_weights[4] = 4.887459349e-01;
    expected_weights[5] = 5.594466685e-02;
    expected_weights[6] = 4.272605627e-02;
    expected_weights[7] = 3.125833420e-01;
  }


  SUBCASE("polyedre concave"){
    printf("POLY 2\n");

    n_vtx = 8;
    PDM_malloc(vtx_coord, 3*n_vtx, double);
    fact = 1.0;  

    // 1er point
    vtx_coord[ 0] = 0.0*fact;
    vtx_coord[ 1] = 0.0*fact;
    vtx_coord[ 2] = 0.0*fact;
    // 2e point
    vtx_coord[ 3] = 1.0*fact;
    vtx_coord[ 4] = 0.0*fact;
    vtx_coord[ 5] = 0.0*fact;
    // 3e point
    vtx_coord[ 6] = 1.0*fact;
    vtx_coord[ 7] = 1.0*fact;
    vtx_coord[ 8] = 0.0*fact;
    // 4e point
    vtx_coord[ 9] = 0.0*fact;
    vtx_coord[10] = 1.0*fact;
    vtx_coord[11] = 0.0*fact;
    // 5e point
    vtx_coord[12] = 0.5*fact;
    vtx_coord[13] = 0.0*fact;
    vtx_coord[14] = 1.0*fact;
    // 6e point
    vtx_coord[15] = 1.0*fact;
    vtx_coord[16] = 0.5*fact;
    vtx_coord[17] = 1.0*fact;
    // 7e point
    vtx_coord[18] = 0.25*fact;
    vtx_coord[19] = 0.25*fact;
    vtx_coord[20] = 1.0*fact;
    // 8e point
    vtx_coord[21] = 0.0*fact;
    vtx_coord[22] = 0.5*fact;
    vtx_coord[23] = 1.0*fact;

    n_cell = 1;
    PDM_malloc(cell_face_idx, n_cell+1, PDM_l_num_t);
    cell_face_idx[0] = 0;
    cell_face_idx[1] = 10;
    PDM_malloc(cell_face, cell_face_idx[n_cell], PDM_l_num_t);
    cell_face[0] = 1;
    cell_face[1] = 2;
    cell_face[2] = 3;
    cell_face[3] = 4;
    cell_face[4] = 5;
    cell_face[5] = 6;
    cell_face[6] = 7;
    cell_face[7] = 8;
    cell_face[8] = 9;
    cell_face[9] = 10;


    n_face = 10;
    PDM_malloc(face_vtx_idx, n_face+1, PDM_l_num_t);
    face_vtx_idx[0]  =  0;
    face_vtx_idx[1]  =  4;
    face_vtx_idx[2]  =  8;
    face_vtx_idx[3]  = 11;
    face_vtx_idx[4]  = 14;
    face_vtx_idx[5]  = 17;
    face_vtx_idx[6]  = 20;
    face_vtx_idx[7]  = 23;
    face_vtx_idx[8]  = 26;
    face_vtx_idx[9]  = 29;
    face_vtx_idx[10] = 32;
    PDM_malloc(face_vtx, face_vtx_idx[n_face], PDM_l_num_t);
    face_vtx[ 0] = 1;
    face_vtx[ 1] = 4;
    face_vtx[ 2] = 3;
    face_vtx[ 3] = 2;
    face_vtx[ 4] = 5;
    face_vtx[ 5] = 6;
    face_vtx[ 6] = 7;
    face_vtx[ 7] = 8;
    face_vtx[ 8] = 1;
    face_vtx[ 9] = 2;
    face_vtx[10] = 5;
    face_vtx[11] = 2;
    face_vtx[12] = 3;
    face_vtx[13] = 6;
    face_vtx[14] = 3;
    face_vtx[15] = 4;
    face_vtx[16] = 7;
    face_vtx[17] = 4;
    face_vtx[18] = 1;
    face_vtx[19] = 8;
    face_vtx[20] = 1;
    face_vtx[21] = 5;
    face_vtx[22] = 8;
    face_vtx[23] = 2;
    face_vtx[24] = 6;
    face_vtx[25] = 5;
    face_vtx[26] = 3;
    face_vtx[27] = 7;
    face_vtx[28] = 6;
    face_vtx[29] = 4;
    face_vtx[30] = 8;
    face_vtx[31] = 7;

    PDM_malloc(cell_ln_to_gn, n_cell, PDM_g_num_t);
    PDM_malloc(face_ln_to_gn, n_face, PDM_g_num_t);
    PDM_malloc(vtx_ln_to_gn, n_vtx, PDM_g_num_t);
    cell_ln_to_gn[0] = 1;
    face_ln_to_gn[0] = 1;
    face_ln_to_gn[1] = 2;
    face_ln_to_gn[2] = 3;
    face_ln_to_gn[3] = 4;
    face_ln_to_gn[4] = 5;
    face_ln_to_gn[5] = 6;
    face_ln_to_gn[6] = 7;
    face_ln_to_gn[7] = 8;
    face_ln_to_gn[8] = 9;
    face_ln_to_gn[9] = 10;
    vtx_ln_to_gn[0] = 1;
    vtx_ln_to_gn[1] = 2;
    vtx_ln_to_gn[2] = 3;
    vtx_ln_to_gn[3] = 4;
    vtx_ln_to_gn[4] = 5;
    vtx_ln_to_gn[5] = 6;
    vtx_ln_to_gn[6] = 7;
    vtx_ln_to_gn[7] = 8;

    n_pts = 1;
    PDM_malloc(pts_coord, 3*n_pts, double);
    PDM_malloc(gnum, n_pts, PDM_g_num_t);

    gnum[0] = 1;
    pts_coord[0] = 0.35;
    pts_coord[1] = 0.25;
    pts_coord[2] = 0.9;

    PDM_malloc(expected_weights, n_pts*n_vtx, double);
    expected_weights[0] = 4.4036112395e-02;
    expected_weights[1] = 1.3896791630e-02;
    expected_weights[2] = 3.2853012922e-02;
    expected_weights[3] = 9.2140830510e-03;
    expected_weights[4] = 2.7299871650e-01;
    expected_weights[5] = 6.1183099473e-02;
    expected_weights[6] = 4.2227095088e-01;
    expected_weights[7] = 1.4354723313e-01;

  }


  PDM_mesh_location_mesh_n_part_set(ml, 1);


  PDM_mesh_location_part_set(ml,
                             0,
                             n_cell,
                             cell_face_idx,
                             cell_face,
                             cell_ln_to_gn,
                             n_face,
                             face_vtx_idx,
                             face_vtx,
                             face_ln_to_gn,
                             n_vtx,
                             vtx_coord,
                             vtx_ln_to_gn);


  PDM_mesh_location_n_part_cloud_set(ml, 0, 1);

  PDM_mesh_location_cloud_set(ml,
                              0,
                              0,
                              n_pts,
                              pts_coord,
                              gnum);

  PDM_mesh_location_compute(ml);

  int *elt_pts_inside_idx = NULL;
  PDM_g_num_t *points_gnum = NULL;
  double *points_coords = NULL;
  double *points_uvw = NULL;
  int *points_weights_idx = NULL;
  double *points_weights = NULL;
  double *points_dist2 = NULL;
  double *points_projected_coords = NULL;

  PDM_mesh_location_points_in_elt_get(ml,
                                      0,
                                      0,
                                      &elt_pts_inside_idx,
                                      &points_gnum,
                                      &points_coords,
                                      &points_uvw,
                                      &points_weights_idx,
                                      &points_weights,
                                      &points_dist2,
                                      &points_projected_coords);



  int located = PDM_mesh_location_n_located_get(ml, 0, 0);

  for (int i=0; i<n_pts*n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<3*n_pts; i++){
    CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  }


  PDM_mesh_location_free(ml);
  PDM_free(vtx_coord);
  PDM_free(cell_face_idx);
  PDM_free(cell_face);
  PDM_free(face_vtx_idx);
  PDM_free(face_vtx);
  PDM_free(cell_ln_to_gn);
  PDM_free(face_ln_to_gn);
  PDM_free(vtx_ln_to_gn);
  PDM_free(pts_coord);
  PDM_free(gnum);
  PDM_free(expected_weights);

}
