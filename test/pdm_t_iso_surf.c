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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_distrib.h"
#include "pdm_multipart.h"
#include "pdm_iso_surface.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *n_part,
           int                   *part_method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static void
_eval_sphere
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  *f = x * x + y * y + z * z - 0.125;

  *df_dx = 2*x;
  *df_dy = 2*y;
  *df_dz = 2*z;
}


static void
_eval_cylinder
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  const double axis[3]   = {0.6, 0.2, 0.3};
  const double center[3] = {0., 0., 0.};
  const double radius    = 0.3;

  double dx = x - center[0];
  double dy = y - center[1];
  double dz = z - center[2];

  double dot = dx*axis[0] + dy*axis[1] + dz*axis[2];
  dot /= (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);

  dx -= dot*axis[0];
  dy -= dot*axis[1];
  dz -= dot*axis[2];

  *f = dx*dx + dy*dy + dz*dz - radius*radius;

  *df_dx = 2*dx;
  *df_dy = 2*dy;
  *df_dz = 2*dz;
}


// length = 2.5
static const int it_max = 6;
static const int power  = 8;
static void
_eval_mandelbulb
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  PDM_UNUSED(df_dx);
  PDM_UNUSED(df_dy);
  PDM_UNUSED(df_dz);

  double dr = 1.;
  double r;
  double theta;
  double phi;

  double xk = x;
  double yk = y;
  double zk = z;

  for (int it = 0; it < it_max; it++) {

    // convert to polar coordinates
    r = sqrt(xk*xk + yk*yk + zk*zk);

    if (r > 2.) {
      break;
    }

    phi   = atan2(z, sqrt(xk*xk + yk*yk));
    theta = atan2(yk, xk);
    dr = pow(r, power-1) * power * dr + 1.;

    if (r > 2) {
      *f = 100.;
      break;
    }

    // scale and rotate the point
    r = pow(r, power);
    theta *= power;
    phi   *= power;

    // convert back to cartesian coordinates
    xk = x + r * cos(phi) * cos(theta);
    yk = y + r * cos(phi) * sin(theta);
    zk = z + r * sin(phi);
  }

  if (PDM_ABS(dr) < 1e-16) {
    *f = 1e9 * PDM_SIGN(dr);
  }
  else {
    if (r < 1e-16) {
      *f = -1;
    } else {
      *f = 0.5*log(r)*r / dr;
    }
  }
}


// length = 2.5
static void
_eval_heart
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  double a = (2*x*x + y*y + z*z - 1);
  *f = a*a*a - (0.1*x*x + y*y)*z*z*z;

  *df_dx = -0.2*x*z*z*z + 12*x*a*a;
  *df_dy = -2*y*z*z*z + 6*y*a*a;
  *df_dz = -(0.3*x*x + 3*y*y)*z*z + 6*z*a*a;
}


// https://ddcampayo.wordpress.com/2016/03/29/taylor-green-vortex-sheet-reduced-units/
static const double U0 = 1.;
static const double L  = 1.;
static const double nu = 1.e-3;
static const double t0 = 1.;

static void
_eval_taylor_green_vortex
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  PDM_UNUSED(z);
  double k  = 2. * PDM_PI / L;
  double F  = exp(-2.*nu*k*k*t0);

  double q = F*2*pow(k,2)*pow(U0,2)*pow(sin(k*x),2)*pow(sin(k*y),2)+2*pow(k,2)*pow(U0,2)*pow(cos(k*x),2)*pow(cos(k*y),2);

  *f = q - 20. * pow(L/U0, 2);

  *df_dx = F*4*pow(k,3)*pow(U0,2)*sin(k*x)*pow(sin(k*y),2)*cos(k*x)-4*pow(k,3)*pow(U0,2)*sin(k*x)*cos(k*x)*pow(cos(k*y),2);
  *df_dy = F*4*pow(k,3)*pow(U0,2)*pow(sin(k*x),2)*sin(k*y)*cos(k*y)-4*pow(k,3)*pow(U0,2)*sin(k*y)*pow(cos(k*x),2)*cos(k*y);
  *df_dz = 0.;
}


// length = 5
static void
_eval_pretzel
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  double b = (x*x + 0.25*y*y - 1);
  double c = (0.25*x*x + y*y - 1);
  double a = b * c;

  *f = a*a + 0.5*z*z - 0.2;

  *df_dx = x*(4*c*c*b + c*b*b);
  *df_dy = y*(c*c*b + 4*c*b*b);
  *df_dz = z;
}


// length = 5
static void
_eval_mcmullen
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  *f = (1 + x*x)*(1 + y*y)*(1 + z*z) + 8*x*y*z - 2.;

  *df_dx = 8*y*z + 2*x*(y*y + 1)*(z*z + 1);
  *df_dy = 8*x*z + 2*y*(x*x + 1)*(z*z + 1);
  *df_dz = 8*x*y + 2*z*(x*x + 1)*(y*y + 1);
}


static void
_eval_chmutov6
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  #define T6(x) (-1 + (x)*(x)*(18 + (x)*(x)*(-48 + 32*(x)*(x))))
  #define dT6dx(x) ((x)*(36 + (x)*(-192 + 192*(x))))

  *f = T6(x) + T6(y) + T6(z);

  *df_dx = dT6dx(x);
  *df_dy = dT6dx(y);
  *df_dz = dT6dx(z);

  #undef T6
  #undef dT6dx
}



#define PERLIN_NOISE_N 8
static double perlin_fx[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static double perlin_fy[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static double perlin_fz[PERLIN_NOISE_N][PERLIN_NOISE_N][PERLIN_NOISE_N];
static const double perlin_step = 1. / (double) (PERLIN_NOISE_N - 1);

static void
_init_perlin_noise (void)
{
  double i_rand_max = 1./(double) RAND_MAX;

  for (int k = 0; k < PERLIN_NOISE_N; k++) {
    for (int j = 0; j < PERLIN_NOISE_N; j++) {
      for (int i = 0; i < PERLIN_NOISE_N; i++) {
        perlin_fx[i][j][k] = 2*rand()*i_rand_max - 1;
        perlin_fy[i][j][k] = 2*rand()*i_rand_max - 1;
        perlin_fz[i][j][k] = 2*rand()*i_rand_max - 1;
      }
    }
  }
}

static inline double
_interpolate_perlin_noise
(
const double a0,
const double a1,
const double x
 )
{
 return (a1 - a0) * (3.0 - x * 2.0) * x * x + a0;
  // return (a1 - a0) * ((x * (x * 6.0 - 15.0) + 10.0) * x * x * x) + a0;
}

static void
_eval_perlin_noise
(
 const double  x,
 const double  y,
 const double  z,
       double *f,
       double *df_dx,
       double *df_dy,
       double *df_dz
 )
{
  PDM_UNUSED(df_dx);
  PDM_UNUSED(df_dy);
  PDM_UNUSED(df_dz);

  int pi = (int) (x / perlin_step);
  int pj = (int) (y / perlin_step);
  int pk = (int) (z / perlin_step);

  pi = PDM_MIN(PDM_MAX(pi, 0), PERLIN_NOISE_N-2);
  pj = PDM_MIN(PDM_MAX(pj, 0), PERLIN_NOISE_N-2);
  pk = PDM_MIN(PDM_MAX(pk, 0), PERLIN_NOISE_N-2);

  double x0 = pi*perlin_step;
  double y0 = pj*perlin_step;
  double z0 = pk*perlin_step;

  double x1 = x0 + perlin_step;
  double y1 = y0 + perlin_step;
  double z1 = z0 + perlin_step;

  double sx = (x - x0) / perlin_step;
  double sy = (y - y0) / perlin_step;
  double sz = (z - z0) / perlin_step;


  double n0, n1, ix0, ix1;

  n0 = (x - x0)*perlin_fx[pi][pj][pk]   + (y - y0)*perlin_fy[pi][pj][pk]   + (z - z0)*perlin_fz[pi][pj][pk];
  n1 = (x - x1)*perlin_fx[pi+1][pj][pk] + (y - y0)*perlin_fy[pi+1][pj][pk] + (z - z0)*perlin_fz[pi+1][pj][pk];
  ix0 = _interpolate_perlin_noise(n0, n1, sx);

  n0 = (x - x0)*perlin_fx[pi][pj+1][pk]   + (y - y1)*perlin_fy[pi][pj+1][pk]   + (z - z0)*perlin_fz[pi][pj][pk];
  n1 = (x - x1)*perlin_fx[pi+1][pj+1][pk] + (y - y1)*perlin_fy[pi+1][pj+1][pk] + (z - z0)*perlin_fz[pi+1][pj+1][pk];
  ix1 = _interpolate_perlin_noise(n0, n1, sx);

  double f0 = _interpolate_perlin_noise(ix0, ix1, sy);


  n0 = (x - x0)*perlin_fx[pi][pj][pk+1]   + (y - y0)*perlin_fy[pi][pj][pk+1]   + (z - z1)*perlin_fz[pi][pj][pk+1];
  n1 = (x - x1)*perlin_fx[pi+1][pj][pk+1] + (y - y0)*perlin_fy[pi+1][pj][pk+1] + (z - z1)*perlin_fz[pi+1][pj][pk+1];
  ix0 = _interpolate_perlin_noise(n0, n1, sx);

  n0 = (x - x0)*perlin_fx[pi][pj+1][pk+1]   + (y - y1)*perlin_fy[pi][pj+1][pk+1]   + (z - z1)*perlin_fz[pi][pj][pk+1];
  n1 = (x - x1)*perlin_fx[pi+1][pj+1][pk+1] + (y - y1)*perlin_fy[pi+1][pj+1][pk+1] + (z - z1)*perlin_fz[pi+1][pj+1][pk+1];
  ix1 = _interpolate_perlin_noise(n0, n1, sx);

  double f1 = _interpolate_perlin_noise(ix0, ix1, sy);

  *f = _interpolate_perlin_noise(f0, f1, sz);
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t          n_vtx_seg = 10;
  double               length    = 1.;
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TETRA4;
  //  5 -> tetra
  //  6 -> pyramid
  //  7 -> prism
  //  8 -> hexa
  //  9 -> poly3d
  int                  n_part    = -1;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &elt_type,
             &n_part,
     (int *) &part_method);
  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 3);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  void (*eval_field_and_gradient) (const double, const double, const double,
                                   double *,
                                   double *, double *, double *) = NULL;



  PDM_multipart_t *mpart = NULL;

  int n_zone = 1;
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;
  if (n_part > 0) {
    mpart = PDM_multipart_create(n_zone,
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
  }

  /*
   *  Create distributed cube
   */

  PDM_g_num_t *distrib_cell = NULL;
  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *vtx_distrib  = NULL;

  PDM_dcube_nodal_t          *dcube   = NULL;
  PDM_dmesh_nodal_to_dmesh_t *dmntodm = NULL;

  PDM_g_num_t  ng_cell = 0;
  PDM_g_num_t  ng_face = 0;
  PDM_g_num_t  ng_vtx  = 0;
  int          dn_cell = 0;
  int          dn_face = 0;
  int          dn_edge = 0;
  int          dn_vtx  = 0;
  double      *dvtx_coord     = NULL;
  int         *dcell_face_idx = NULL;
  PDM_g_num_t *dcell_face     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int          n_face_group = 0;
  PDM_g_num_t *dface_cell      = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    dcube = PDM_dcube_nodal_gen_create (comm,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        length,
                                        0,//-0.5*length,//-0.1,
                                        0,//-0.5*length,//-0.2,
                                        0,//-0.5*length,//0.07,
                                        elt_type,
                                        1,
                                        PDM_OWNERSHIP_KEEP);

    PDM_dcube_nodal_gen_build (dcube);

    PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

    if (0) {
      double noise = 0.2 * length / (double) (n_vtx_seg - 1);
      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 3; j++) {
          dvtx_coord[3*i+j] += noise * (rand() / (double) RAND_MAX - 0.5);
        }
        double y = dvtx_coord[3*i+1];
        dvtx_coord[3*i] += 0.15*length*sin(7*y);
        double x = dvtx_coord[3*i];
        dvtx_coord[3*i+2] += 0.2*length*cos(8*x);
      }
    }

    if(0 == 1) {
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    }

    dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

    PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

    PDM_dmesh_t* dmesh = NULL;
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);

    dn_cell = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                         &dcell_face,
                                         &dcell_face_idx,
                                         PDM_OWNERSHIP_KEEP);

    dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         &dface_edge,
                                         &dface_edge_idx,
                                         PDM_OWNERSHIP_KEEP);

    dn_edge = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         &dedge_vtx,
                                         &dedge_vtx_idx,
                                         PDM_OWNERSHIP_KEEP);

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &distrib_edge);
    assert(distrib_edge != NULL);

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);
    assert(distrib_face != NULL);

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_CELL, &distrib_cell);
    assert(distrib_cell != NULL);


    if (n_part > 0) {
      PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);
    }
  }

  else {
    // polyvol_gen
    PDM_poly_vol_gen (comm,
                      -0.1,//-0.5,
                      0.1,//-0.5,
                      0.13,//-0.5,
                      length,
                      length,
                      length,
                      n_vtx_seg,
                      n_vtx_seg,
                      n_vtx_seg,
                      1,
                      0,
                      &ng_cell,
                      &ng_face,
                      &ng_vtx,
                      &n_face_group,
                      &dn_cell,
                      &dn_face,
                      &dn_vtx,
                      &dcell_face_idx,
                      &dcell_face,
                      &dface_cell,
                      &dface_edge_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);

    printf("TO DO: build edges\n");
    abort();
  }


  if (n_part > 0) {
    PDM_multipart_run_ppart(mpart);
  }

  // eval_field_and_gradient = &_eval_cylinder;
  _init_perlin_noise();
  eval_field_and_gradient = &_eval_perlin_noise;

  PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

  PDM_iso_surface_eval_field_and_gradient_set(isos,
                                                eval_field_and_gradient);

  double  *dfield          = NULL;
  double  *dgradient_field = NULL;
  double **pfield          = NULL;
  double **pgradient_field = NULL;

  if (n_part > 0) {

    pfield          = (double **) malloc(sizeof(double *) * n_part);
    pgradient_field = (double **) malloc(sizeof(double *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      // printf("vtx: %d %p\n", n_vtx, (void *) vtx_coord);

      pfield[i_part]          = (double *) malloc(sizeof(double) * n_vtx);
      pgradient_field[i_part] = (double *) malloc(sizeof(double) * n_vtx * 3);

      for (int i = 0; i < n_vtx; i++) {
        eval_field_and_gradient(vtx_coord[3*i  ],
                                vtx_coord[3*i+1],
                                vtx_coord[3*i+2],
                                &pfield[i_part][i],
                                &pgradient_field[i_part][3*i],
                                &pgradient_field[i_part][3*i+1],
                                &pgradient_field[i_part][3*i+2]);
      }

      int *cell_face_idx = NULL;
      int *cell_face     = NULL;
      int *face_edge_idx = NULL;
      int *face_edge     = NULL;
      int *edge_vtx_idx  = NULL;
      int *edge_vtx      = NULL;

      int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                       &cell_face,
                                                       &cell_face_idx,
                                                       PDM_OWNERSHIP_KEEP);
      // printf("cell_face: %d %p %p\n", n_cell, (void *) cell_face_idx, (void *) cell_face);

      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                       &face_edge,
                                                       &face_edge_idx,
                                                       PDM_OWNERSHIP_KEEP);
      // printf("face_edge: %d %p %p\n", n_face, (void *) face_edge_idx, (void *) face_edge);

      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       0,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);
      // printf("edge_vtx: %d %p %p\n", n_edge, (void *) face_edge_idx, (void *) face_edge);

      PDM_g_num_t *cell_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      // printf("cell_ln_to_gn : %p\n", (void *) cell_ln_to_gn);

      PDM_g_num_t *face_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      // printf("face_ln_to_gn : %p\n", (void *) face_ln_to_gn);

      PDM_g_num_t *edge_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      // printf("edge_ln_to_gn : %p\n", (void *) edge_ln_to_gn);

      PDM_g_num_t *vtx_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      // printf("vtx_ln_to_gn : %p\n", (void *) vtx_ln_to_gn);



      PDM_iso_surface_part_set (isos,
                                i_part,
                                n_cell,
                                n_face,
                                n_edge,
                                n_vtx,
                                cell_face_idx,
                                cell_face,
                                face_edge_idx,
                                face_edge,
                                edge_vtx,
                                NULL,//face_vtx_idx,
                                NULL,//face_vtx,
                                cell_ln_to_gn,
                                face_ln_to_gn,
                                edge_ln_to_gn,
                                vtx_ln_to_gn,
                                vtx_coord);

      PDM_iso_surface_part_field_set (isos,
                                      i_part,
                                      pfield[i_part]);

      // PDM_iso_surface_part_gradient_field_set (isos,
      //                                          i_part,
      //                                          pgradient_field[i_part]);
    }

  }

  else {


    // Compute dfield and gradient field
    dfield          = (double *) malloc(     dn_vtx * sizeof(double));
    dgradient_field = (double *) malloc( 3 * dn_vtx * sizeof(double));

    for (int i = 0; i < dn_vtx; i++) {
      eval_field_and_gradient(dvtx_coord[3*i  ],
                              dvtx_coord[3*i+1],
                              dvtx_coord[3*i+2],
                              &dfield[i],
                              &dgradient_field[3*i],
                              &dgradient_field[3*i+1],
                              &dgradient_field[3*i+2]);
    }

    if (1 && n_rank == 1) {
      FILE *f = fopen("mandelbuld_sdf.vtk", "w");
      fprintf(f, "# vtk DataFile Version 2.0\n");
      fprintf(f, "mesh\n");
      fprintf(f, "ASCII\n");
      fprintf(f, "DATASET STRUCTURED_GRID\n");
      fprintf(f, "DIMENSIONS "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n", n_vtx_seg, n_vtx_seg, n_vtx_seg);
      fprintf(f, "POINTS %d double\n", dn_vtx);
      for(int i = 0; i < dn_vtx; ++i) {
        fprintf(f, "%f %f %f\n", dvtx_coord[3*i  ], dvtx_coord[3*i+1], dvtx_coord[3*i+2]);
      }
      fprintf(f, "POINT_DATA %d\n", dn_vtx);
      fprintf(f, "SCALARS sdf double 1\n");
      fprintf(f, "LOOKUP_TABLE default\n");
      for(int i = 0; i < dn_vtx; ++i) {
        fprintf(f, "%f\n", dfield[i]);
      }
      fclose(f);
    // abort();
    }


    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      dcell_face,
                                      dcell_face_idx);
    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      dface_edge,
                                      dface_edge_idx);
    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      dedge_vtx,
                                      dedge_vtx_idx);

    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_CELL  , distrib_cell);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_FACE  , distrib_face);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_EDGE  , distrib_edge);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, vtx_distrib);


    PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
    PDM_iso_surface_dfield_set     (isos, dfield         );
    // PDM_iso_surface_dgrad_field_set(isos, dgradient_field);
  }

  PDM_iso_surface_compute(isos);

  // char name[999];
  // sprintf(name, "iso_surface_%dproc", n_rank);
  // PDM_iso_surface_write(isos, name);

  PDM_iso_surface_free(isos);

  if (n_part > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      free(pfield[i_part]);
      free(pgradient_field[i_part]);
    }
    free(pfield);
    free(pgradient_field);

    PDM_multipart_free(mpart);
  }
  else {
    free(dfield);
    free(dgradient_field);
  }

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    PDM_dmesh_nodal_to_dmesh_free(dmntodm);
    PDM_dcube_nodal_gen_free(dcube);
  } else {
    free(distrib_cell);
    free(distrib_face);
    free(distrib_edge);
    free(vtx_distrib);
    free(dvtx_coord);
    free(dcell_face_idx);
    free(dcell_face);
    free(dface_edge_idx);
    free(dface_edge);
    free(dedge_vtx_idx);
    free(dedge_vtx);
    free(dface_cell);
    free(dface_vtx);
    free(dface_group_idx);
    free(dface_group);
  }


  double min_elaps_create;
  double max_elaps_create;
  double min_cpu_create;
  double max_cpu_create;
  double min_elaps_create2;
  double max_elaps_create2;
  double min_cpu_create2;
  double max_cpu_create2;
  double min_elaps_exch;
  double max_elaps_exch;
  double min_cpu_exch;
  double max_cpu_exch;

  PDM_part_to_block_global_timer_get (comm,
                                      &min_elaps_create,
                                      &max_elaps_create,
                                      &min_cpu_create,
                                      &max_cpu_create,
                                      &min_elaps_create2,
                                      &max_elaps_create2,
                                      &min_cpu_create2,
                                      &max_cpu_create2,
                                      &min_elaps_exch,
                                      &max_elaps_exch,
                                      &min_cpu_exch,
                                      &max_cpu_exch);

  if (i_rank == 0) {
    printf("Global time in PDM_part_to_block : \n");
    printf("   - min max elaps create  : %12.5e %12.5e\n", min_elaps_create, max_elaps_create);
    printf("   - min max elaps create2 : %12.5e %12.5e\n", min_elaps_create2, max_elaps_create2);
    printf("   - min max elaps exch    : %12.5e %12.5e\n", min_elaps_exch, max_elaps_exch);
    fflush(stdout);
  }


  PDM_block_to_part_global_timer_get (comm,
                                      &min_elaps_create,
                                      &max_elaps_create,
                                      &min_cpu_create,
                                      &max_cpu_create,
                                      &min_elaps_exch,
                                      &max_elaps_exch,
                                      &min_cpu_exch,
                                      &max_cpu_exch);

  if (i_rank == 0) {
    printf("Global time in pdm_block_to_part : \n");
    printf("   - min max elaps create  : %12.5e %12.5e\n", min_elaps_create, max_elaps_create);
    printf("   - min max elaps exch    : %12.5e %12.5e\n", min_elaps_exch, max_elaps_exch);
    fflush(stdout);
  }


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

}
