#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_dbbtree.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_multipart.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_plane.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const int    verbose = 1;
static const int    vtk     = 1;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -h               This message.\n\n");


  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp (argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else
      _usage (EXIT_FAILURE);
    i++;
  }
}

/* Create a 4 plane volume */

static void
_create_volume_4planes
(
double  *edge,
double  *direction_pt,
double   theta,
double   eps,
double **n_in,
double **pt_plane_in
)
{
  *n_in        = malloc(sizeof(double) * 12);
  *pt_plane_in = malloc(sizeof(double) * 12);

  double *n = *n_in;
  double *pt_plane = *pt_plane_in;

  // Determine eps translation planes
  // B--{eps}--G
  double CB[3] = {edge[3]-edge[6], edge[4]-edge[7], edge[5]-edge[8]};
  pt_plane[0] = edge[3] + (1+eps) * CB[0];
  pt_plane[1] = edge[4] + (1+eps) * CB[1];
  pt_plane[2] = edge[5] + (1+eps) * CB[2];
  n[0] = -CB[0];
  n[1] = -CB[1];
  n[2] = -CB[2];

  // A--{eps}--H
  double CA[3] = {edge[3]-edge[0], edge[4]-edge[1], edge[5]-edge[2]};
  pt_plane[3] = edge[3] + (1+eps) * CA[0];
  pt_plane[4] = edge[4] + (1+eps) * CA[1];
  pt_plane[5] = edge[5] + (1+eps) * CA[2];
  n[3] = -CA[0];
  n[4] = -CA[1];
  n[5] = -CA[2];

  // Determine theta angle planes E---D---F
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);

  double AB[3] = {edge[6]-edge[0], edge[7]-edge[1], edge[8]-edge[2]};
  double inverse_module_AB = 1 / PDM_MODULE(AB);
  double AB_normalised[3] = {AB[0] * inverse_module_AB, AB[1] * inverse_module_AB, AB[2] * inverse_module_AB};

  pt_plane[6]  = (cos_theta + (1 - cos_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[1] * (1 - cos_theta) - AB_normalised[2] * sin_theta) * direction_pt[1];
  pt_plane[6] += (AB_normalised[0] * AB_normalised[2] * (1 - cos_theta) + AB_normalised[1] * sin_theta) * direction_pt[2];
  pt_plane[7]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_theta) + AB_normalised[2] * sin_theta) * direction_pt[0];
  pt_plane[7] += (cos_theta + (1 - cos_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[7] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_theta) - AB_normalised[0] * sin_theta) * direction_pt[2];
  pt_plane[8]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_theta) - AB_normalised[1] * sin_theta) * direction_pt[0];
  pt_plane[8] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_theta) + AB_normalised[0] * sin_theta) * direction_pt[1];
  pt_plane[8] += (cos_theta + (1 - cos_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  double prod_vect[3];
  double CE[3] = {pt_plane[6]-edge[3], pt_plane[7]-edge[4], pt_plane[8]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CE);
  double ED[3] = {direction_pt[0] - pt_plane[6], direction_pt[1] - pt_plane[7], direction_pt[2] - pt_plane[8]};
  double sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, ED));
  n[6] = sign * prod_vect[0];
  n[7] = sign * prod_vect[1];
  n[8] = sign * prod_vect[2];

  double cos_minus_theta = cos(-theta);
  double sin_minus_theta = sin(-theta);

  pt_plane[9]   = (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[0] * AB_normalised[0]) * direction_pt[0];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[1] * (1 - cos_minus_theta) - AB_normalised[2] * sin_minus_theta) * direction_pt[1];
  pt_plane[9]  += (AB_normalised[0] * AB_normalised[2] * (1 - cos_minus_theta) + AB_normalised[1] * sin_minus_theta) * direction_pt[2];
  pt_plane[10]  = (AB_normalised[1] * AB_normalised[0] * (1 - cos_minus_theta) + AB_normalised[2] * sin_minus_theta) * direction_pt[0];
  pt_plane[10] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[1] * AB_normalised[1]) * direction_pt[1];
  pt_plane[10] += (AB_normalised[1] * AB_normalised[2] * (1 - cos_minus_theta) - AB_normalised[0] * sin_minus_theta) * direction_pt[2];
  pt_plane[11]  = (AB_normalised[2] * AB_normalised[0] * (1 - cos_minus_theta) - AB_normalised[1] * sin_minus_theta) * direction_pt[0];
  pt_plane[11] += (AB_normalised[2] * AB_normalised[1] * (1 - cos_minus_theta) + AB_normalised[0] * sin_minus_theta) * direction_pt[1];
  pt_plane[11] += (cos_minus_theta + (1 - cos_minus_theta) * AB_normalised[2] * AB_normalised[2]) * direction_pt[2];
  double CF[3] = {pt_plane[9]-edge[4], pt_plane[10]-edge[4], pt_plane[11]-edge[5]};
  PDM_CROSS_PRODUCT(prod_vect, CA, CF);
  double FD[3] = {direction_pt[0] - pt_plane[9], direction_pt[1] - pt_plane[10], direction_pt[2] -  pt_plane[11]};
  sign = PDM_SIGN(PDM_DOT_PRODUCT(prod_vect, FD));
  n[9]  = sign * prod_vect[0];
  n[10] = sign * prod_vect[1];
  n[11] = sign * prod_vect[2];
}

static void
_volume_4planes_to_4triangles
(
double       *edge,
double       *direction_pt,
double       *n,
double       *pt_plane,
double      **vtx_coord_in,
PDM_g_num_t **vtx_g_num_in,
int         **face_vtx_in
)
{
  *vtx_coord_in = malloc(sizeof(double) * 12 * 3);
  *vtx_g_num_in = malloc(sizeof(PDM_g_num_t) * 12);
  *face_vtx_in  = malloc(sizeof(int) * 12);

  double *vtx_coord = *vtx_coord_in;
  PDM_g_num_t *vtx_g_num = *vtx_g_num_in;
  int *face_vtx = *face_vtx_in;

  for (int i = 0; i < 12; i ++) {
    vtx_g_num[i] = i + 1;
    face_vtx[i]  = i +1;
  }

  double x[3];
  double origin[3];
  double normal[3];
  double proj_x[3];

  // double cross_prod[3];
  // PDM_CROSS_PRODUCT(cross_prod, n, n + 3);
  // log_trace("cross_prod = %f %f %f\n", cross_prod[0], cross_prod[1], cross_prod[2]);

  // log_trace("G = %f %f %f\n", pt_plane[0], pt_plane[1], pt_plane[2]);
  // log_trace("H = %f %f %f\n", pt_plane[3 + 0], pt_plane[3 + 1], pt_plane[3 + 2]);

  // translation plane 1
  // project E on 1
  for (int k = 0; k < 3; k++) {
    origin[k] = pt_plane[k]; // G
    normal[k] = n[k];
    x[k] = pt_plane[6 + k]; // E
    vtx_coord[k] = pt_plane[k]; // G
  }
  PDM_plane_projection2(x, origin, normal, proj_x);
  // project D on 1
  for (int k = 0; k < 3; k++) {
    x[k] = direction_pt[k]; // D
    vtx_coord[3 + k] = proj_x[k]; // proj_E
  }
  PDM_plane_projection2(x, origin, normal, proj_x);
  for (int k = 0; k < 3; k++) {
    vtx_coord[6 + k] = proj_x[k]; // proj_D
  }

  // translation plane 2
  // project E on 2
  for (int k = 0; k < 3; k++) {
    origin[k] = pt_plane[3 + k]; // H
    normal[k] = n[3 + k];
    x[k] = pt_plane[6 + k]; // E
    vtx_coord[9 + k] = pt_plane[3 + k]; // H
  }
  PDM_plane_projection2(x, origin, normal, proj_x);
  // project D on 2
  for (int k = 0; k < 3; k++) {
    x[k] = direction_pt[k]; // D
    vtx_coord[12 + k] = proj_x[k]; // proj_E
  }
  PDM_plane_projection2(x, origin, normal, proj_x);
  for (int k = 0; k < 3; k++) {
    vtx_coord[15 + k] = proj_x[k]; // proj_D
  }

  // double DDproj[3];
  // double CG[3];
  // for (int k = 0; k < 3; k++) {
  //   DDproj[k] = vtx_coord[6 + k] - direction_pt[k];
  //   CG[k]     = pt_plane[k] - edge[3 + k];
  // }

  // log_trace("DDproj = %f %f %f\n", DDproj[0], DDproj[1], DDproj[2]);
  // log_trace("CG     = %f %f %f\n", CG[0], CG[1], CG[2]);
  // log_trace("DDproj - CG = %f %f %f\n", CG[0] - DDproj[0], CG[1] - DDproj[1], CG[2] - DDproj[2]);

  // rotation plane 3
  for (int k = 0; k < 3; k++) {
    vtx_coord[18 + 3*0 + k] = edge[k]; // A
    vtx_coord[18 + 3*1 + k] = edge[6 + k]; // B
    vtx_coord[18 + 3*2 + k] = pt_plane[6 + k]; // E
  }

  // rotation plane 4
  for (int k = 0; k < 3; k++) {
    vtx_coord[27 + 3*0 + k] = edge[k]; // A
    vtx_coord[27 + 3*1 + k] = edge[6 + k]; // B
    vtx_coord[27 + 3*2 + k] = pt_plane[9 + k]; // E
  }

}



static void
_create_volume_4planes_tata
(
double  *edge,
double  *direction_pt,
double   theta,
double   eps,
double **n_in,
double **pt_plane_in
)
{
  *n_in        = malloc(sizeof(double) * 12);
  *pt_plane_in = malloc(sizeof(double) * 12);

  double *n = *n_in;
  double *pt_plane = *pt_plane_in;

  double u[3], v[3], w[3];
  for (int i = 0; i < 3; i++) {
    u[i] = direction_pt[i] - edge[3+i];
    w[i] = edge[6+i] - edge[i];
  }


  double mu = PDM_MODULE(u);
  double mw = PDM_MODULE(w);

  for (int i = 0; i < 3; i++) {
    u[i] /= mu;
    w[i] /= mw;
  }


  PDM_CROSS_PRODUCT(v, w, u);

  if (PDM_ABS(PDM_DOT_PRODUCT(u, w)) > 1e-13) {
    PDM_error(__FILE__, __LINE__, 0,"!!! u.w = %e\n", PDM_DOT_PRODUCT(u, w));
  }

  double c = cos(theta);
  double s = sin(theta);


  for (int i = 0; i < 3; i++) {

    pt_plane[i] = edge[i] - mw*eps*w[i];
    n[i] = w[i];

    pt_plane[3+i] = edge[6+i] + mw*eps*w[i];
    n[3+i] = -w[i];

    pt_plane[6+i] = edge[3+i] + c*u[i] + s*v[i];
    n[6+i] = s*u[i] - c*v[i];

    pt_plane[9+i] = edge[3+i] + c*u[i] - s*v[i];
    n[9+i] = s*u[i] + c*v[i];

  }

}

static void
_vtk_write_volume
(
 const char *filename,
 double     *plane_point,
 double     *plane_normal
 )
{
  const double scale = 1.;

  double u[3], v[3], w[3];
  for (int i = 0; i < 3; i++) {
    u[i] = plane_point[3+i] - plane_point[i];
  }

  PDM_CROSS_PRODUCT(v, u, plane_normal + 6);
  PDM_CROSS_PRODUCT(w, plane_normal + 9, u);

  double vtx_coord[18];
  memcpy(vtx_coord, plane_point, sizeof(double) * 6);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      vtx_coord[6  + 3*i + j] = plane_point[j] + i*u[j] + v[j];
      vtx_coord[12 + 3*i + j] = plane_point[j] + i*u[j] + w[j];
    }
  }

  int face_vtx[14] = {
    1, 5, 3,
    2, 4, 6,
    1, 3, 4, 2,
    1, 2, 6, 5
  };

  int face_vtx_idx[5] = {0, 3, 6, 10, 14};

  PDM_vtk_write_polydata(filename,
                         6,
                         vtx_coord,
                         NULL,
                         4,
                         face_vtx_idx,
                         face_vtx,
                         NULL,
                         NULL);
}

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int n_part = 1;

  _read_args (argc,
              argv);

  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_UNUSED(verbose);
  PDM_UNUSED(vtk);

  /* Create surface meshes */

  const double x_center = 0;
  const double y_center = 0;
  const double z_center = 0;

  // Background spherical mesh

  const int    back_nu     = 100;
  const int    back_nv     = 100;
  const double back_radius = 10;

  double      *d_back_vtx_coord    = NULL;
  int         *d_back_face_vtx_idx = NULL;
  PDM_g_num_t *d_back_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx    = NULL;
  PDM_g_num_t *back_distrib_face   = NULL;

  PDM_sphere_surf_gen(comm,
                      back_nu,
                      back_nv,
                      x_center,
                      y_center,
                      z_center,
                      back_radius,
                      &d_back_vtx_coord,
                      &d_back_face_vtx_idx,
                      &d_back_face_vtx,
                      &back_distrib_vtx,
                      &back_distrib_face);

  // "Volume" spherical mesh

  const int    vol_nu     = 10;
  const int    vol_nv     = 10;
  const double vol_radius = 10;

  PDM_dmesh_nodal_t *vol_dmn;

  PDM_sphere_surf_gen_nodal(comm,
                            vol_nu,
                            vol_nv,
                            x_center,
                            y_center,
                            z_center,
                            vol_radius,
                            &vol_dmn);

  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int n_zone                   = 1;
  int *n_part_zones            = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0]              = n_part;

  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart, 0, vol_dmn);

  PDM_multipart_run_ppart(mpart);

  free(n_part_zones);

  // Vertices
  int     i_part          = 0;
  double *p_vol_vtx_coord = NULL;
  int p_vol_n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &p_vol_vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *vol_vtx_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                  &vol_vtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Edges
  int *p_vol_edge_vtx     = NULL;
  int *p_vol_edge_vtx_idx = NULL;
  int  p_vol_n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                          0,
                                                          i_part,
                                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                          &p_vol_edge_vtx,
                                                          &p_vol_edge_vtx_idx,
                                                          PDM_OWNERSHIP_KEEP);

  if (p_vol_edge_vtx_idx != NULL) free(p_vol_edge_vtx_idx);

  PDM_g_num_t *vol_edge_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &vol_edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  // Faces
  int *p_vol_face_edge     = NULL;
  int *p_vol_face_edge_idx = NULL;
  int p_vol_n_face = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                         &p_vol_face_edge,
                                                         &p_vol_face_edge_idx,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_g_num_t *vol_face_ln_to_gn = NULL;
  PDM_multipart_part_ln_to_gn_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &vol_face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

  int *p_vol_face_vtx = NULL;
  PDM_compute_face_vtx_from_face_and_edge(p_vol_n_face,
                                          p_vol_face_edge_idx,
                                          p_vol_face_edge,
                                          p_vol_edge_vtx,
                                          &p_vol_face_vtx);

  // Get coordinates of the background mesh associated to the local face distribution

  int dn_back_face = back_distrib_face[i_rank+1] - back_distrib_face[i_rank];

  PDM_g_num_t *d_back_face_ln_to_gn = (PDM_g_num_t * ) malloc(dn_back_face * sizeof(PDM_g_num_t));
  for (int i = 0; i < dn_back_face; ++i) {
    d_back_face_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }

  int dn_back_vtx = back_distrib_vtx[i_rank+1] - back_distrib_vtx[i_rank];

  PDM_g_num_t *d_back_vtx_ln_to_gn = (PDM_g_num_t * ) malloc(dn_back_vtx * sizeof(PDM_g_num_t));
  for (int i = 0; i < dn_back_vtx; ++i) {
    d_back_vtx_ln_to_gn[i] = back_distrib_face[i_rank] + i + 1;
  }

  PDM_g_num_t *p_back_vtx_ln_to_gn = NULL;
  int         *p_back_face_vtx_idx = NULL;
  int         *p_back_face_vtx     = NULL;
  int          p_back_n_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           back_distrib_face,
                                                           d_back_face_vtx_idx,
                                                           d_back_face_vtx,
                                                           dn_back_face,
                                  (const PDM_g_num_t *)    d_back_face_ln_to_gn,
                                                           &p_back_n_vtx,
                                                           &p_back_vtx_ln_to_gn,
                                                           &p_back_face_vtx_idx,
                                                           &p_back_face_vtx);

  int      *n_part_p_back_n_vtx            = malloc(sizeof(int) * n_part);
  n_part_p_back_n_vtx[0]                   = p_back_n_vtx;
  double  **n_part_p_back_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        back_distrib_vtx,
                                        d_back_vtx_coord,
                                        n_part_p_back_n_vtx,
                 (const PDM_g_num_t **) &p_back_vtx_ln_to_gn,
                                        &n_part_p_back_vtx_coord);

  double *p_back_vtx_coord = n_part_p_back_vtx_coord[0];

  // Create the extents faces as a partition and get associated coords
  double      *background_box_extents = malloc(sizeof(double)      * dn_back_face * 6);
  int          eps                    = 1.0e-6;
  for (int iface = 0; iface < dn_back_face; iface++) {
    double *tmp_extents = background_box_extents + 6*iface;
    for (int k = 0; k < 3; k++) {
      tmp_extents[k]     =  1.0e15;
      tmp_extents[3 + k] = -1.0e15;
    }
    for (int ivtx = p_back_face_vtx_idx[iface]; ivtx < p_back_face_vtx_idx[iface+1]; ivtx++) {
      int vtx_id = p_back_face_vtx[ivtx]-1;
      double *tmp_coord = p_back_vtx_coord + 3*vtx_id;
      for (int k = 0; k < 3; k++) {
        tmp_extents[k]     =  PDM_MIN(tmp_coord[k], tmp_extents[k]);
        tmp_extents[3 + k] =  PDM_MAX(tmp_coord[k], tmp_extents[3 + k]);
      }
    } // end loop on vertices of iface
    for (int k = 0; k < 3; k++) {
      if (PDM_ABS(tmp_extents[k] - tmp_extents[3 + k]) < 1.0e-15) { //
        tmp_extents[k]     -= eps;
        tmp_extents[3 + k] += eps;
      }
    }
  } // end loop on background faces

  if (vtk) {
    char filename3[999];
    sprintf(filename3, "dbbtree_boxes_%d.vtk", i_rank);
    PDM_vtk_write_boxes(filename3,
                        dn_back_face,
                        background_box_extents,
                        d_back_face_ln_to_gn);
  }

  // Create dbbtree from surface mesh boxes

  const int dim = 3;
  double l_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                         -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i = 0; i < dn_back_face; i++) {
    for (int k = 0; k < 3; k++) {
      l_extents[k]     = PDM_MIN(l_extents[k],     background_box_extents[6*i + k]);
      l_extents[k + 3] = PDM_MAX(l_extents[k + 3], background_box_extents[6*i + k + 3]);
    }
  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX(max_range, g_extents[i+3] - g_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    g_extents[i]   -= max_range * 1.1e-3;
    g_extents[i+3] += max_range * 1.0e-3;
  }

  log_trace("g_extents = %f %f %f %f %f %f\n", g_extents[0], g_extents[1], g_extents[2], g_extents[3], g_extents[4], g_extents[5]);

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create(comm, dim, g_extents);

  PDM_box_set_t *box_set = PDM_dbbtree_boxes_set(dbbt,
                                                 n_part,
                                                 &dn_back_face,
                               (const double **) &background_box_extents,
                          (const PDM_g_num_t **) &d_back_face_ln_to_gn);

  // Create volumes using normals (first "volume" mesh then background mesh)

  int    *edge_face_normal_stride = PDM_array_zeros_int(p_vol_n_edge);
  double *edge_face_normal        = malloc(sizeof(double) * 3 * p_vol_n_edge * 2);

  double *face_center             = malloc(sizeof(double) * 3 * p_vol_n_face);
  double *face_normal             = malloc(sizeof(double) * 3 * p_vol_n_face);

  for (int iface = 0; iface < p_vol_n_face; iface++) {
    // compute normal vector
    int *tmp_face_edge = p_vol_face_edge + p_vol_face_edge_idx[iface];
    // with face_edge edge_vtx
    int edge_id1 = PDM_ABS(tmp_face_edge[0]) - 1;
    int vtx_id1, vtx_id2;
    if (tmp_face_edge[0] > 0) {
      vtx_id1 = p_vol_edge_vtx[2*edge_id1]-1;
      vtx_id2 = p_vol_edge_vtx[2*edge_id1+1]-1;
    } else {
      vtx_id1 = p_vol_edge_vtx[2*edge_id1+1]-1;
      vtx_id2 = p_vol_edge_vtx[2*edge_id1]-1;
    }
    int edge_id2 = PDM_ABS(tmp_face_edge[1]) - 1;
    int vtx_id3;
    if (p_vol_edge_vtx[2*edge_id2]-1 != vtx_id1 && p_vol_edge_vtx[2*edge_id2]-1 != vtx_id2) {
      vtx_id3 = p_vol_edge_vtx[2*edge_id2]-1;
    } else {
      vtx_id3 = p_vol_edge_vtx[2*edge_id2+1]-1;
    }

    double vect[6];
    for (int i = 0; i < 3; i++) {
      vect[i]     = p_vol_vtx_coord[3*vtx_id2 + i] - p_vol_vtx_coord[3*vtx_id1 + i];
      vect[3 + i] = p_vol_vtx_coord[3*vtx_id3 + i] - p_vol_vtx_coord[3*vtx_id1 + i];
    }
    int normal[3];
    PDM_CROSS_PRODUCT(normal, vect, vect + 3);

    for (int i = 0; i < 3; i++) {
      face_center[3*iface + i] = 1./3 * (p_vol_vtx_coord[3*vtx_id1 + i] +  p_vol_vtx_coord[3*vtx_id2 + i] + p_vol_vtx_coord[3*vtx_id3 + i]);
      face_normal[3*iface + i] = normal[i];
    }


    // fill in edge_face_normal
    for (int iedge = 0; iedge < 3; iedge++) {
      int edge_id = PDM_ABS(tmp_face_edge[iedge])-1;
      for (int i = 0; i < 3; i++) {
        edge_face_normal[6*edge_id + 3*edge_face_normal_stride[edge_id] + i] = normal[i];
      }
      edge_face_normal_stride[edge_id]++;
    } // end loop on all edges of iface
  } // end loop on volume faces

  // Remove empty cell in edge_face_normal
  int *edge_face_normal_idx = PDM_array_new_idx_from_sizes_int(edge_face_normal_stride, p_vol_n_edge);
  int idx_read  = 0;
  int idx_write = 0;
  for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
    idx_write = edge_face_normal_idx[iedge];
    idx_read  = 2*iedge;
    // if (idx_read != idx_write) {
      for (int inormal = 0; inormal < edge_face_normal_stride[iedge]; inormal++) {
        for (int i = 0; i < 3; i++) {
          edge_face_normal[3*idx_write + 3*inormal + i] = edge_face_normal[3*idx_read + 3*inormal + i];
        }
      }
    // }
  } // end loop on edges

  if (verbose) {
    PDM_log_trace_array_double(edge_face_normal, edge_face_normal_idx[p_vol_n_edge], "edge_face_normal");
  }

  // Get normal contributions for other procs

  // Part to Block

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &vol_edge_ln_to_gn,
                                                      NULL,
                                                      &p_vol_n_edge,
                                                      n_part,
                                                      comm);

  int    *block_stride = NULL;
  double *block_data   = NULL;
  PDM_part_to_block_exch(ptb,
                         3*sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &edge_face_normal_stride,
               (void **) &edge_face_normal,
                         &block_stride,
               (void **) &block_data);

  int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *block_g_num = PDM_part_to_block_block_gnum_get(ptb);

  // Block to Part

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_g_num,
                                                                        n_elt_block,
                                                 (const PDM_g_num_t **) &vol_edge_ln_to_gn,
                                                                        &p_vol_n_edge,
                                                                        n_part,
                                                                        comm);
  int two = 2;
  PDM_block_to_part_exch_in_place(btp,
                                  3*sizeof(double),
                                  PDM_STRIDE_CST_INTERLACED,
                                  &two,
                         (void *) block_data,
                                  NULL,
                        (void **) &edge_face_normal);

  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb);

  // compute volume angle

  // only for vtk
  double *middle_pt_coord = malloc(sizeof(double) * 3 * p_vol_n_edge);
  double *direction_vect  = malloc(sizeof(double) * 3 * p_vol_n_edge);

  double  theta_min         = 1.0e-1; // WARNING: randomly chosen value
  double  eps2              = 1.0e-1; // WARNING: randomly chosen value
  double *edge_normal       = malloc(sizeof(double) * p_vol_n_edge * 3 * 4);
  double *edge_pt_plane     = malloc(sizeof(double) * p_vol_n_edge * 3 * 4);
  for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
    double  theta;
    double  direction_pt[3];
    double *direction = malloc(sizeof(double) * 3);
    for (int i = 0; i < 3; i++) {
      direction[i] = 0.5 * (edge_face_normal[6*iedge + i] + edge_face_normal[6*iedge + 3 + i]);
    }
    double dot_prod = PDM_DOT_PRODUCT(edge_face_normal + 6*iedge, edge_face_normal + 6*iedge + 3);
    double module1 = PDM_MODULE(edge_face_normal + 6*iedge);
    double module2 = PDM_MODULE(edge_face_normal + 6*iedge + 3);
    theta = acos(dot_prod / (module1 * module2)); // make sure between -1 and 1
    theta += theta_min;
    if (theta > 3) { // WARNING: shouldn't be bigger than PI
      theta = 3;
    }
    if (theta < 0) {
      log_trace("theta is negative!!!\n");
      theta = PDM_ABS(theta);
    }



    int *tmp_edge_vtx = p_vol_edge_vtx + 2*iedge;
    int vtx_id1 = tmp_edge_vtx[0] - 1;
    int vtx_id2 = tmp_edge_vtx[1] - 1;
    double edge_vector[3];
    double edge[9];
    for (int i = 0; i < 3; i++) {
      edge[i]     = p_vol_vtx_coord[3*vtx_id1 + i];
      edge[6 + i] = p_vol_vtx_coord[3*vtx_id2 + i];
      edge[3 + i]     = 0.5 * (p_vol_vtx_coord[3*vtx_id1 + i] + p_vol_vtx_coord[3*vtx_id2 + i]);
      edge_vector[i] = edge[6 + i] - edge[i];
    }

    double fix = PDM_DOT_PRODUCT(edge_vector, direction) / PDM_DOT_PRODUCT(edge_vector, edge_vector);
    for (int i = 0; i < 3; i++) {
      direction[i] -= fix*edge_vector[i];
      direction_pt[i] = edge[3 + i] + 2 * direction[i]; // WARNING: 2 is a random distance choice

      // only vtk
      middle_pt_coord[3*iedge + i] = edge[3 + i];
      direction_vect[3*iedge + i]  = direction_pt[i];
    }

    double *normal   = edge_normal + 12*iedge;
    double *pt_plane = edge_pt_plane + 12*iedge;

    _create_volume_4planes_tata(edge,
                           direction_pt,
                           theta,
                           eps2,
                           &normal,
                           &pt_plane);

    memcpy(edge_normal   + 12*iedge, normal,   sizeof(double) * 12);
    memcpy(edge_pt_plane + 12*iedge, pt_plane, sizeof(double) * 12);

    if (vtk) {

      // TO DO: output volumes in a more elegant way
      double      *vtx_coord = NULL;
      PDM_g_num_t *vtx_g_num = NULL;
      int         *face_vtx  = NULL;
      _volume_4planes_to_4triangles(edge,
                                    direction_pt,
                                    normal,
                                    pt_plane,
                                    &vtx_coord,
                                    &vtx_g_num,
                                    &face_vtx);

      char filename4[999];
      sprintf(filename4, "planes_of_edge_id_%ld.vtk", vol_edge_ln_to_gn[iedge]);
      PDM_vtk_write_std_elements(filename4,
                                 12,
                                 vtx_coord,
                                 vtx_g_num,
                                 PDM_MESH_NODAL_TRIA3,
                                 4,
                                 face_vtx,
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);

      sprintf(filename4, "volume_of_edge_id_%ld.vtk", vol_edge_ln_to_gn[iedge]);
      _vtk_write_volume(filename4,
                        pt_plane,
                        normal);

    const char *normal_name = "normal";

    char filename38[999];
    sprintf(filename38, "normal_of_edge_id_%ld.vtk", vol_edge_ln_to_gn[iedge]);
    PDM_vtk_write_point_cloud_with_field(filename38,
                                         4,
                                         pt_plane,
                                         NULL,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &normal,
                                         0,
                                         NULL,
                                         NULL);

    }


  } // end loop on edges

  if (vtk) {

    const char *direction_name = "direction";

    char filename5[999];
    sprintf(filename5, "direction_%d.vtk", i_rank);
    PDM_vtk_write_point_cloud_with_field(filename5,
                                         p_vol_n_edge,
                                         middle_pt_coord,
                                         vol_edge_ln_to_gn,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &direction_name,
                       (const double **) &direction_vect,
                                         0,
                                         NULL,
                                         NULL);

    const char *normal_name = "face_normal";

    char filename6[999];
    sprintf(filename6, "face_normal_%d.vtk", i_rank);
    PDM_vtk_write_point_cloud_with_field(filename6,
                                         p_vol_n_face,
                                         face_center,
                                         vol_face_ln_to_gn,
                                         NULL,
                                         0,
                                         NULL,
                                         NULL,
                                         1,
                         (const char **) &normal_name,
                       (const double **) &face_normal,
                                         0,
                                         NULL,
                                         NULL);
  }

  // Call PDM_dbbtree_volumes_intersect_boxes

  int *volume_plane_idx = PDM_array_new_idx_from_const_stride_int(4, p_vol_n_edge);

  int         *volume_boxes_idx   = NULL;
  PDM_g_num_t *volume_boxes_g_num = NULL;
  PDM_dbbtree_volumes_intersect_boxes(dbbt,
                                      p_vol_n_edge,
                                      vol_edge_ln_to_gn,
                                      volume_plane_idx,
                                      edge_normal,
                                      edge_pt_plane,
                                      &volume_boxes_idx,
                                      &volume_boxes_g_num);

  if (verbose) {
    log_trace("VOLUME-BOX INTERSECTION\n");
    for (int iedge = 0; iedge < p_vol_n_edge; iedge++) {
      log_trace("--> volume "PDM_FMT_G_NUM" is intersected by ", vol_edge_ln_to_gn[iedge]);
      log_trace(" %d boxes ", volume_boxes_idx[iedge+1] - volume_boxes_idx[iedge]);
      for (int i = volume_boxes_idx[iedge]; i < volume_boxes_idx[iedge+1]; i++) {
        log_trace("%d ", volume_boxes_g_num[i]);
      }
      log_trace("\n");
    }
  }

  // VTK output of surface mesh with tagged elements for each volume mesh edges

  int n_boxes = volume_boxes_idx[p_vol_n_edge];

  PDM_g_num_t *p_box_volume_g_num  = malloc(sizeof(PDM_g_num_t) * n_boxes);

  for (int ivol = 0; ivol < p_vol_n_edge; ivol++) {
    for (int ibox = volume_boxes_idx[ivol]; ibox < volume_boxes_idx[ivol+1]; ibox++) {
      p_box_volume_g_num[ibox] = vol_edge_ln_to_gn[ivol];
    } // end loop on ivol boxes
  } // end loop on volumes

  // Part to Block

  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                                                    1.,
                                                                    &volume_boxes_g_num,
                                                                    back_distrib_face,
                                                                    &n_boxes,
                                                                    n_part,
                                                                    comm);

  int         *d_box_volume_stride = NULL;
  PDM_g_num_t *d_box_volume_g_num  = NULL;
  int         *part_stride         = PDM_array_const_int(n_boxes, 1);
  PDM_part_to_block_exch(ptb2,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &part_stride,
               (void **) &p_box_volume_g_num,
                         &d_box_volume_stride,
               (void **) &d_box_volume_g_num);

  PDM_g_num_t* distrib = PDM_part_to_block_adapt_partial_block_to_block(ptb2,
                                                                        &d_box_volume_stride,
                                                                        back_distrib_face[n_rank]);
  free(distrib);

  int dn_block_face =  PDM_part_to_block_n_elt_block_get(ptb2);
  log_trace("dn_block_face = %d and dn_back_face = %d\n", dn_block_face, dn_back_face);

  free(part_stride);
  PDM_part_to_block_free(ptb2);

  int *d_box_volume_idx = PDM_array_new_idx_from_sizes_int(d_box_volume_stride, dn_back_face);

  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
  PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, vol_dmn);
  PDM_dmesh_nodal_generate_distribution(vol_dmn);

  PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                   PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
  PDM_dmesh_t* vol_dmesh = NULL;
  PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &vol_dmesh);

  PDM_g_num_t* dedge_distrib = NULL;
  PDM_dmesh_distrib_get(vol_dmesh, PDM_MESH_ENTITY_EDGE, &dedge_distrib);

  int total_n_edges = dedge_distrib[n_rank];

  PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);

  int         **volume       = malloc(sizeof(int  *) * total_n_edges);
  const char  **volume_names = malloc(sizeof(char *) * total_n_edges);

  for (int ivol = 0; ivol < total_n_edges; ivol++) {
    volume[ivol] = PDM_array_zeros_int(dn_back_face);
    volume_names[ivol] = malloc(sizeof(char) * 99);
    sprintf((char * restrict) volume_names[ivol], "edge_%d.vtk", ivol+1);

  }

  log_trace("total_n_edges = %d\n", total_n_edges);
  log_trace("dn_back_face = %d\n", dn_back_face);

  for (int ibox = 0; ibox < dn_back_face; ibox++) {
    PDM_g_num_t box_gnum = d_back_face_ln_to_gn[ibox];
    for (int ivol = d_box_volume_idx[ibox]; ivol < d_box_volume_idx[ibox+1]; ivol++) {
      PDM_g_num_t vol_gnum = d_box_volume_g_num[ivol];
      volume[vol_gnum-1][ibox] = 1;
    }
  }

  if (vtk) {

    char filename1[999];
    sprintf(filename1, "background_mesh_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename1,
                               p_back_n_vtx,
                               p_back_vtx_coord,
                               p_back_vtx_ln_to_gn,
                               2, // PDM_MESH_NODAL_TRIA3
                               dn_back_face,
                               p_back_face_vtx,
                               d_back_face_ln_to_gn,
                               total_n_edges,
                               volume_names,
                (const int **) volume);

    char filename2[999];
    sprintf(filename2, "volume_mesh_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename2,
                               p_vol_n_vtx,
                               p_vol_vtx_coord,
                               vol_vtx_ln_to_gn,
                               2, // PDM_MESH_NODAL_TRIA3
                               p_vol_n_face,
                               p_vol_face_vtx,
                               vol_face_ln_to_gn,
                               0,
                               NULL,
                               NULL);

  }

  PDM_MPI_Finalize ();
}
