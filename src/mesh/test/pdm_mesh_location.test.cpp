#include "doctest/extensions/doctest_mpi.h"
#include "float.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "math.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"

#include "pdm_mesh_location.h"




static double tol = 1e-10;

MPI_TEST_CASE("[pdm_mesh_location] - tria", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 3;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0e-6;  

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

  PDM_l_num_t n_face = 1;
  PDM_l_num_t face_vtx_idx[2] = {0, 3};
  PDM_l_num_t face_vtx[3] = {1, 2, 3};

  PDM_g_num_t face_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[3] = {1, 2, 3};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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
  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = fact/2. + 0.5*(1.-fact/2.);
  pts_coord[1] = 0.5;
  pts_coord[2] = 0.0;

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

  double expected_weights[3] = {0.25, 0.25, 0.5};

  for (int i=0; i<n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<3*n_pts; i++){
    CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - quad", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 4;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.5*1e-6;  

  // 1er point
  vtx_coord[ 0] = 0.0;
  vtx_coord[ 1] = 0.0;
  vtx_coord[ 2] = 0.0;
  // 2e point
  vtx_coord[ 3] = 1.0*fact;
  vtx_coord[ 4] = 0.0;
  vtx_coord[ 5] = 0.0;
  // 3e point
  vtx_coord[ 6] = 1.0;
  vtx_coord[ 7] = 1.0;
  vtx_coord[ 8] = 0.0;
  // 4e point
  vtx_coord[ 9] = 0.0;
  vtx_coord[10] = 1.0*fact;
  vtx_coord[11] = 0.0;

  PDM_l_num_t n_face = 1;
  PDM_l_num_t face_vtx_idx[2] = {0, 4};
  PDM_l_num_t face_vtx[4] = {1, 2, 3, 4};

  PDM_g_num_t face_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[4] = {1, 2, 3, 4};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.5*fact;
  pts_coord[1] = 0.5*fact;
  pts_coord[2] = 0.0;

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

  double expected_weights[4] = {0.33333325, 0.33333325, 2.5e-07, 0.33333325};

  for (int i=0; i<n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<3*n_pts; i++){
    CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - poly", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 6;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0;
  // double a = 1.0;//1.54919328;

  // 1er point
  vtx_coord[ 0] = 0.0;
  vtx_coord[ 1] = 0.0;
  vtx_coord[ 2] = 0.0;
  // 2e point
  vtx_coord[ 3] = 0.5;
  vtx_coord[ 4] = 0.5;
  vtx_coord[ 5] = 0.0;
  // 3e point
  vtx_coord[ 6] = 1.0;
  vtx_coord[ 7] = 0.0;
  vtx_coord[ 8] = 0.0;
  // 4e point
  vtx_coord[ 9] = 1.0;
  vtx_coord[10] = 1.0;
  vtx_coord[11] = 0.0;
  // 5e point
  vtx_coord[12] = 0.5;
  vtx_coord[13] = 1.5;
  vtx_coord[14] = 0.0;
  // 6e point
  vtx_coord[15] = 0.0;
  vtx_coord[16] = 1.0;
  vtx_coord[17] = 0.0;

  PDM_l_num_t n_face = 1;
  PDM_l_num_t face_vtx_idx[2] = {0, 6};
  PDM_l_num_t face_vtx[6] = {1, 2, 3, 4, 5, 6};

  PDM_g_num_t face_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[6] = {1, 2, 3, 4, 5, 6};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.5;
  pts_coord[1] = 1.0;
  pts_coord[2] = 0.0;

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


  printf("elt_pts_inside_idx : [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords : [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw : [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx : [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights : [%.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], 
                                                            points_weights[3], points_weights[4], points_weights[5]);
  printf("points_gnum : [%li]\n", points_gnum[0]);
  printf("points_dist2 : [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords : [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);

  // double expected_weights[n_vtx] = {0.33333325, 0.33333325, 2.5e-07, 0.33333325};
  // 
  // for (int i=0; i<n_vtx; i++){
  //   printf("diff poids = %.16e\n", fabs(points_weights[i] - expected_weights[i]));
  //   CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  // }
  // for (int i=0; i<3*n_pts; i++){
  //   CHECK(fabs(points_coords[i] - pts_coord[i]) < tol);
  // }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - tetra", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 4;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0e-2;  

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

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_vtx_idx[2] = {0,4};
  PDM_l_num_t cell_vtx[4] = {1, 2, 3, 4};

  PDM_g_num_t cell_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[4] = {1, 2, 3, 4};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 2;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[2] = {1, 2};

  pts_coord[0] = 0.45;
  pts_coord[1] = 0.05*fact;
  pts_coord[2] = -0.03*fact;

  pts_coord[3] = 0.75;
  pts_coord[4] = 0.25*fact;
  pts_coord[5] = 1.0*fact;

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


  double expected_weights[8] = {0.3625, 0.1275, 0.4125, 0.0975,
                                          0.18749218769507753, 0.0, 0.43749218769563264, 0.37501562460928983};

  for (int i=0; i<n_pts*n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<n_pts; i++){

    CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
    CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
    CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - hexa nodal", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 8;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0;  

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

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_vtx_idx[2] = {0, 8};
  PDM_l_num_t cell_vtx[8] = {1, 2, 4, 3, 5, 6, 8, 7};

  PDM_g_num_t cell_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 1.5;
  pts_coord[1] = 0.0;
  pts_coord[2] = -1.0;

  PDM_mesh_location_n_part_cloud_set(ml, 0, 1);

  PDM_mesh_location_cloud_set(ml,
                              0,
                              0,
                              n_pts,
                              pts_coord,
                              gnum);


  PDM_mesh_location_method_set(ml, PDM_MESH_LOCATION_LOCATE_ALL_TGT);


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


  printf("elt_pts_inside_idx p1: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], points_weights[3], 
                                                                                                           points_weights[4], points_weights[5], points_weights[6], points_weights[7]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);

  double recompute_point[3] = {0., 0., 0.};
  for (int i =0; i<n_vtx; i++){
    // printf("%.16e; %.16e; %.16e\n", vtx_coord[3*(cell_vtx[i]-1)], vtx_coord[3*(cell_vtx[i]-1)+1], vtx_coord[3*(cell_vtx[i]-1)+2]);
    // printf("%.16e\n", points_weights[i]);
    recompute_point[0] += points_weights[i]*vtx_coord[3*(cell_vtx[i]-1)];
    recompute_point[1] += points_weights[i]*vtx_coord[3*(cell_vtx[i]-1)+1];
    recompute_point[2] += points_weights[i]*vtx_coord[3*(cell_vtx[i]-1)+2];
  }

  printf("recomputed point : \n%.16e; %.16e; %.16e\n", recompute_point[0], recompute_point[1], recompute_point[2]);

  // printf("elt_pts_inside_idx p2: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  // printf("points_coords p2: [%.16e; %.16e; %.16e]\n", points_coords[3], points_coords[4], points_coords[5]);
  // printf("points_uvw p2: [%.16e; %.16e; %.16e]\n", points_uvw[3], points_uvw[4], points_uvw[5]);
  // printf("points_weights_idx p2: [%i; %i]\n", points_weights_idx[1], points_weights_idx[2]);
  // printf("points_weights p2: [%.16e; %.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e; %.16e]\n", points_weights[8], points_weights[9], points_weights[10], points_weights[11], 
  //                                                                                                          points_weights[12], points_weights[13], points_weights[14], points_weights[15]);
  // printf("points_gnum p2: [%li]\n", points_gnum[1]);
  // printf("points_dist2 p2: [%.16e]\n", points_dist2[1]);
  // printf("points_projected_coords p2: [%.16e; %.16e; %.16e]\n", points_projected_coords[3], points_projected_coords[4], points_projected_coords[5]);




  // double expected_weights[n_pts*n_vtx] = {0.303875, 0.248625, 0.133875, 0.163625,
  //                                         0.053625, 0.043875, 0.023625, 0.028875};

  // for (int i=0; i<n_pts*n_vtx; i++){
  //   CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  // }
  // for (int i=0; i<n_pts; i++){

  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  // }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - hexa", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 8;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0;  

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

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_face_idx[2] = {0, 6};
  PDM_l_num_t cell_face[6] = {1, 2, 3, 4, 5, 6};

  PDM_g_num_t cell_ln_to_gn[1] = {1};

  PDM_g_num_t n_face = 6;
  PDM_l_num_t face_vtx_idx[7] = {0, 4, 8, 12, 16, 20, 24};
  PDM_l_num_t face_vtx[24] = {1, 2, 4, 3,
                                                1, 5, 6, 2,
                                                2, 6, 8, 4,
                                                4, 8, 7, 3,
                                                3, 7, 1, 5,
                                                5, 7, 8, 6};

  PDM_g_num_t face_ln_to_gn[6] = {1, 2, 3, 4, 5, 6};
  PDM_g_num_t vtx_ln_to_gn[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.75;
  pts_coord[1] = 0.75;
  pts_coord[2] = 0.0;

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

  printf("elt_pts_inside_idx p1: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], points_weights[3], 
                                                                                                           points_weights[4], points_weights[5], points_weights[6], points_weights[7]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);


  // printf("elt_pts_inside_idx p2: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  // printf("points_coords p2: [%.16e; %.16e; %.16e]\n", points_coords[3], points_coords[4], points_coords[5]);
  // printf("points_uvw p2: [%.16e; %.16e; %.16e]\n", points_uvw[3], points_uvw[4], points_uvw[5]);
  // printf("points_weights_idx p2: [%i; %i]\n", points_weights_idx[1], points_weights_idx[2]);
  // printf("points_weights p2: [%.16e; %.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e; %.16e]\n", points_weights[8], points_weights[9], points_weights[10], points_weights[11], 
  //                                                                                                          points_weights[12], points_weights[13], points_weights[14], points_weights[15]);
  // printf("points_gnum p2: [%li]\n", points_gnum[1]);
  // printf("points_dist2 p2: [%.16e]\n", points_dist2[1]);
  // printf("points_projected_coords p2: [%.16e; %.16e; %.16e]\n", points_projected_coords[3], points_projected_coords[4], points_projected_coords[5]);


  for (int i =0; i<n_vtx; i++){
    printf("%.16e; %.16e; %.16e\n", vtx_coord[3*i], vtx_coord[3*i+1], vtx_coord[3*i+2]);
    printf("%.16e\n", points_weights[i]);
  }


  // printf("recomputed point : [%.16e, %.16e, %.16e]\n", points_weights[0]*vtx_coord[0]+points_weights[1]*vtx_coord[3]+points_weights[2]*vtx_coord[6]+points_weights[3]*vtx_coord[ 9]+points_weights[4]*vtx_coord[12]+points_weights[5]*vtx_coord[15]+points_weights[6]*vtx_coord[18]+points_weights[7]*vtx_coord[21],
                                                       // points_weights[0]*vtx_coord[1]+points_weights[1]*vtx_coord[4]+points_weights[2]*vtx_coord[7]+points_weights[3]*vtx_coord[10]+points_weights[4]*vtx_coord[13]+points_weights[5]*vtx_coord[16]+points_weights[6]*vtx_coord[19]+points_weights[7]*vtx_coord[22],
                                                       // points_weights[0]*vtx_coord[2]+points_weights[1]*vtx_coord[5]+points_weights[2]*vtx_coord[8]+points_weights[3]*vtx_coord[11]+points_weights[4]*vtx_coord[14]+points_weights[5]*vtx_coord[17]+points_weights[6]*vtx_coord[20]+points_weights[7]*vtx_coord[23]);



  // double expected_weights[n_pts*n_vtx] = {0.303875, 0.248625, 0.133875, 0.163625,
  //                                         0.053625, 0.043875, 0.023625, 0.028875};

  // for (int i=0; i<n_pts*n_vtx; i++){
  //   CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  // }
  // for (int i=0; i<n_pts; i++){

  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  // }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - pyra", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 5;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0e-11;

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

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_vtx_idx[2] = {0, 5};
  PDM_l_num_t cell_vtx[5] = {1, 2, 4, 3, 5};

  PDM_g_num_t cell_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[5] = {1, 2, 4, 3, 5};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.45 ;
  pts_coord[1] = 0.35;
  pts_coord[2] = 0.05*fact;

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

  printf("elt_pts_inside_idx p1: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], points_weights[3], points_weights[4]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);

  double expected_weights[5] = {0.293560606, 0.2231060606, 0.093560606, 0.1231060606, 0.2666666666};

  for (int i=0; i<n_pts*n_vtx; i++){
    CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  }
  for (int i=0; i<n_pts; i++){

    CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
    CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
    CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}

MPI_TEST_CASE("[pdm_mesh_location] - prism", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 6;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0e-1;

  // 1er point
  vtx_coord[ 0] = -0.1;
  vtx_coord[ 1] = 0.1;
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
  vtx_coord[ 9] = 0.1;
  vtx_coord[10] = 0.1;
  vtx_coord[11] = 1.0;
  // 5e point
  vtx_coord[12] = 1.0;
  vtx_coord[13] = 0.0;
  vtx_coord[14] = 1.0;
  // 6e point
  vtx_coord[15] = 0.0;
  vtx_coord[16] = 1.0;
  vtx_coord[17] = 1.0;

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_vtx_idx[2] = {0, 6};
  PDM_l_num_t cell_vtx[6] = {1, 2, 3, 4, 5, 6};

  PDM_g_num_t cell_ln_to_gn[1] = {1};
  PDM_g_num_t vtx_ln_to_gn[6] = {1, 2, 3, 4, 5, 6};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.0;//*(2-fact);
  pts_coord[1] = 0.0;//*(2-fact);
  pts_coord[2] = 0.5;

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

  printf("elt_pts_inside_idx p1: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e; %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], 
                                                                            points_weights[3], points_weights[4], points_weights[5]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);

  // double expected_weights[n_pts*n_vtx] = {0.293560606, 0.2231060606, 0.093560606, 0.1231060606, 0.2666666666};

  // for (int i=0; i<n_pts*n_vtx; i++){
  //   CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  // }
  // for (int i=0; i<n_pts; i++){

  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  // }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}






MPI_TEST_CASE("[pdm_mesh_location] - polyhedron", 1) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;

  int n_vtx = 8;
  double *vtx_coord = NULL;
  PDM_malloc(vtx_coord, 3*n_vtx, double);
  double fact = 1.0;  

  // 1er point
  vtx_coord[ 0] = 0.0*fact;
  vtx_coord[ 1] = 0.0*fact;
  vtx_coord[ 2] = 0.0*fact;
  // 2e point
  vtx_coord[ 3] = 2.0*fact;
  vtx_coord[ 4] = 0.0*fact;
  vtx_coord[ 5] = 0.0*fact;
  // 3e point
  vtx_coord[ 6] = 2.0*fact;
  vtx_coord[ 7] = 2.0*fact;
  vtx_coord[ 8] = 0.0*fact;
  // 4e point
  vtx_coord[ 9] = 0.0*fact;
  vtx_coord[10] = 2.0*fact;
  vtx_coord[11] = 0.0*fact;
  // 5e point
  vtx_coord[12] = 1.0*fact;
  vtx_coord[13] = 0.0*fact;
  vtx_coord[14] = 1.0*fact;
  // 6e point
  vtx_coord[15] = 2.0*fact;
  vtx_coord[16] = 1.0*fact;
  vtx_coord[17] = 1.0*fact;
  // 7e point
  vtx_coord[18] = 1.0*fact;
  vtx_coord[19] = 2.0*fact;
  vtx_coord[20] = 1.0*fact;
  // 8e point
  vtx_coord[21] = 0.0*fact;
  vtx_coord[22] = 1.0*fact;
  vtx_coord[23] = 1.0*fact;

  PDM_l_num_t n_cell = 1;
  PDM_l_num_t cell_face_idx[2] = {0, 10};
  PDM_l_num_t cell_face[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  PDM_g_num_t cell_ln_to_gn[1] = {1};

  PDM_g_num_t n_face = 10;
  PDM_l_num_t face_vtx_idx[11] = {0, 4, 8, 11, 14, 17, 20, 23, 26, 29, 32};
  PDM_l_num_t face_vtx[32] = {1, 4, 3, 2,
                              5, 6, 7, 8, 
                              1, 2, 5,
                              2, 3, 6,
                              3, 4, 7,
                              4, 1, 8,
                              1, 5, 8,
                              2, 6, 5,
                              3, 7, 6,
                              4, 8, 7};

  PDM_g_num_t face_ln_to_gn[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  PDM_g_num_t vtx_ln_to_gn[8] = {1, 2, 3, 4, 5, 6, 7, 8};

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
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

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.45;
  pts_coord[1] = 0.45;
  pts_coord[2] = 0.;

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

  printf("elt_pts_inside_idx p1: [%i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e; %.16e\n                  %.16e; %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], points_weights[3], 
                                                                                                           points_weights[4], points_weights[5], points_weights[6], points_weights[7]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);



  // double expected_weights[n_pts*n_vtx] = {0.303875, 0.248625, 0.133875, 0.163625,
  //                                         0.053625, 0.043875, 0.023625, 0.028875};

  // for (int i=0; i<n_pts*n_vtx; i++){
  //   CHECK(fabs(points_weights[i] - expected_weights[i]) < tol);
  // }
  // for (int i=0; i<n_pts; i++){

  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+0] - pts_coord[3*(points_gnum[i]-1)+0]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+1] - pts_coord[3*(points_gnum[i]-1)+1]) < tol);
  //   CHECK(fabs(points_coords[3*(gnum[i]-1)+2] - pts_coord[3*(points_gnum[i]-1)+2]) < tol);
  // }
  PDM_free(vtx_coord);
  PDM_free(pts_coord);
}







MPI_TEST_CASE("[pdm_mesh_location] - generated mesh", 1) {

  PDM_g_num_t n_vtx_seg               = 2;
  double      length                  = 1.;
  double      zero_x                  = 0.0;
  double      zero_y                  = 0.0;
  double      zero_z                  = 0.0;
  int         n_part                  = 1;

  /*
   * Generate a cube
   */
  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);;
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  int          *pn_vtx                 = NULL;
  int          *pn_edge                = NULL;
  int          *pn_face                = NULL;
  int          *pn_cell                = NULL;
  int          *pn_surface             = NULL;
  int          *pn_ridge               = NULL;
  double      **pvtx_coord             = NULL;
  int         **pedge_vtx              = NULL;
  int         **pface_edge_idx         = NULL;
  int         **pface_edge             = NULL;
  int         **pface_vtx              = NULL;
  int         **pcell_face_idx         = NULL;
  int         **pcell_face             = NULL;
  int         **psurface_face_idx      = NULL;
  int         **psurface_face          = NULL;
  int         **pridge_edge_idx        = NULL;
  int         **pridge_edge            = NULL;
  PDM_g_num_t **pvtx_ln_to_gn          = NULL;
  PDM_g_num_t **pedge_ln_to_gn         = NULL;
  PDM_g_num_t **pface_ln_to_gn         = NULL;
  PDM_g_num_t **pcell_ln_to_gn         = NULL;
  PDM_g_num_t **psurface_face_ln_to_gn = NULL;
  PDM_g_num_t **pridge_edge_ln_to_gn   = NULL;

  PDM_generate_mesh_parallelepiped_ngon(comm,
                                        PDM_MESH_NODAL_PRISM6,
                                        1,
                                        NULL,
                                        zero_x,
                                        zero_y,
                                        zero_z,
                                        length,
                                        length,
                                        length,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_part,
                                        PDM_SPLIT_DUAL_WITH_PTSCOTCH,
                                       &pn_vtx,
                                       &pn_edge,
                                       &pn_face,
                                       &pn_cell,
                                       &pvtx_coord,
                                       &pedge_vtx,
                                       &pface_edge_idx,
                                       &pface_edge,
                                       &pface_vtx,
                                       &pcell_face_idx,
                                       &pcell_face,
                                       &pvtx_ln_to_gn,
                                       &pedge_ln_to_gn,
                                       &pface_ln_to_gn,
                                       &pcell_ln_to_gn,
                                       &pn_surface,
                                       &psurface_face_idx,
                                       &psurface_face,
                                       &psurface_face_ln_to_gn,
                                       &pn_ridge,
                                       &pridge_edge_idx,
                                       &pridge_edge,
                                       &pridge_edge_ln_to_gn);

  PDM_mesh_location_t *ml = PDM_mesh_location_create(1,
                                                     comm,
                                                     PDM_OWNERSHIP_KEEP);
  PDM_mesh_location_mesh_n_part_set(ml, 1);

  PDM_l_num_t pface_vtx_idx[10] = {0, 3, 6, 10, 14, 17, 21, 25, 28, 32};

  PDM_mesh_location_part_set(ml,
                             0,
                             pn_cell[0],
                             pcell_face_idx[0],
                             pcell_face[0],
                             pcell_ln_to_gn[0],
                             pn_face[0],
                             pface_vtx_idx,
                             pface_vtx[0],
                             pface_ln_to_gn[0],
                             pn_vtx[0],
                             pvtx_coord[0],
                             pvtx_ln_to_gn[0]);

  int n_pts = 1;
  double *pts_coord = NULL;
  PDM_malloc(pts_coord, 3*n_pts, double);
  PDM_g_num_t gnum[1] = {1};

  pts_coord[0] = 0.66;
  pts_coord[1] = 0.8;
  pts_coord[2] = 0.15;

  // pts_coord[3] = 1.0;
  // pts_coord[4] = 1.0;
  // pts_coord[5] = 0.0;

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


  printf("elt_pts_inside_idx p1: [%i, %i, %i]\n", elt_pts_inside_idx[0], elt_pts_inside_idx[1], elt_pts_inside_idx[2]);
  printf("points_coords p1: [%.16e; %.16e; %.16e]\n", points_coords[0], points_coords[1], points_coords[2]);
  printf("points_uvw p1: [%.16e; %.16e; %.16e]\n", points_uvw[0], points_uvw[1], points_uvw[2]);
  printf("points_weights_idx p1: [%i; %i]\n", points_weights_idx[0], points_weights_idx[1]);
  printf("points_weights p1: [%.16e; %.16e; %.16e;\n                  %.16e; %.16e; %.16e]\n", points_weights[0], points_weights[1], points_weights[2], points_weights[3], 
                                                                                                           points_weights[4], points_weights[5]);
  printf("points_gnum p1: [%li]\n", points_gnum[0]);
  printf("points_dist2 p1: [%.16e]\n", points_dist2[0]);
  printf("points_projected_coords p1: [%.16e; %.16e; %.16e]\n", points_projected_coords[0], points_projected_coords[1], points_projected_coords[2]);
  
  PDM_free(pts_coord);
}