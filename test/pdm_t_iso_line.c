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
#include "pdm_poly_surf_gen.h"
#include "pdm_distrib.h"
#include "pdm_iso_surface.h"
#include "pdm_multipart.h"

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


static
void
_dump_dmesh
(
 PDM_MPI_Comm  comm,
 const int     dn_face,
 const int     dn_vtx,
 int          *dface_vtx_idx,
 PDM_g_num_t  *dface_vtx,
 double       *dvtx_coord,
 double       *dfield,
 double       *dgradient_field
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *distrib_face = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t *distrib_vtx  = PDM_compute_entity_distribution(comm, dn_vtx);

  PDM_g_num_t *dface_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face);
  for (int i = 0; i < dn_face; i++) {
    dface_ln_to_gn[i] = distrib_face[i_rank] + i + 1;
  }

  int          pn_vtx        = 0;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;
  int         *pface_vtx_idx = NULL;
  int         *pface_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_face,
                                                           dface_vtx_idx,
                                                           dface_vtx,
                                                           dn_face,
                                     (const PDM_g_num_t *) dface_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pface_vtx_idx,
                                                           &pface_vtx);

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_vtx,
                                        dvtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free (tmp_pvtx_coord);

  double** tmp_pfield = NULL;
  PDM_part_dfield_to_pfield(comm,
                            1,
                            sizeof(double),
                            distrib_vtx,
       (unsigned char    *) dfield,
                            &pn_vtx,
    (const PDM_g_num_t **)  &pvtx_ln_to_gn,
       (unsigned char ***)  &tmp_pfield);
  double *pfield = tmp_pfield[0];
  free(tmp_pfield);

  double** tmp_pgradient_field = NULL;
  PDM_part_dfield_to_pfield(comm,
                            1,
                            3*sizeof(double),
                            distrib_vtx,
       (unsigned char    *) dgradient_field,
                            &pn_vtx,
    (const PDM_g_num_t **)  &pvtx_ln_to_gn,
       (unsigned char ***)  &tmp_pgradient_field);
  double *pgradient_field = tmp_pgradient_field[0];
  free(tmp_pgradient_field);


  char filename[999];
  sprintf(filename, "mesh_%2.2d.vtk", i_rank);
  // PDM_vtk_write_polydata(filename,
  //                        pn_vtx,
  //                        pvtx_coord,
  //                        pvtx_ln_to_gn,
  //                        dn_face,
  //                        pface_vtx_idx,
  //                        pface_vtx,
  //                        dface_ln_to_gn,
  //                        NULL);
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", pn_vtx);
  for (int i = 0; i < pn_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", pvtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", dn_face, dn_face + pface_vtx_idx[dn_face]);
  for (int i = 0; i < dn_face; i++) {
    fprintf(f, "%d", pface_vtx_idx[i+1] - pface_vtx_idx[i]);
    for (int j = pface_vtx_idx[i]; j < pface_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", pface_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POINT_DATA %d\n", pn_vtx);
  // fprintf(f, "SCALARS field double 1\n");
  // fprintf(f, "LOOKUP_TABLE default\n");
  // for (int i = 0; i < pn_vtx; i++) {
  //   fprintf(f, "%f\n", pfield[i]);
  // }
  // fprintf(f, "VECTORS gradient double\n");
  // for (int i = 0; i < pn_vtx; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     fprintf(f, "%.20lf ", pgradient_field[3*i+j]);
  //   }
  //   fprintf(f, "\n");
  // }
  fprintf(f, "FIELD pts_field 2\n");
  fprintf(f, "field 1 %d double\n", pn_vtx);
  for (int i = 0; i < pn_vtx; i++) {
    fprintf(f, "%f\n", pfield[i]);
  }
  fprintf(f, "gradient 3 %d double\n", pn_vtx);
  for (int i = 0; i < pn_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", pgradient_field[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_DATA %d\n", dn_face);
  fprintf(f, "SCALARS face_gnum long 1\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (int i = 0; i < dn_face; i++) {
    fprintf(f, PDM_FMT_G_NUM"\n", dface_ln_to_gn[i]);
  }

  fclose(f);


  free(pvtx_ln_to_gn);
  free(pface_vtx_idx);
  free(pface_vtx);
  free(pvtx_coord);
  free(pfield);
  free(pgradient_field);
  free(dface_ln_to_gn);
  free(distrib_face);
  free(distrib_vtx);
}



static void
_eval_circle
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

  double x1 = 0.0;
  double y1 = 0.0;
  double r1 = 0.25;

  *f = (x-x1)*(x-x1) + (y-y1)*(y-y1) - r1*r1;

  *df_dx = 2*(x-x1);
  *df_dy = 2*(y-y1);
  *df_dz = 0.;
}


static void
_eval_mickey
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

  double x1 =  0.23;
  double y1 =  0.28-0.07;
  double r1 =  0.17;

  double x2 = -0.23;
  double y2 =  0.28-0.07;
  double r2 =  0.17;

  double x3 = 0.;
  double y3 = 0.-0.07;
  double r3 = 0.32;

  double f1 = (x-x1)*(x-x1) + (y-y1)*(y-y1) - r1*r1;
  double f2 = (x-x2)*(x-x2) + (y-y2)*(y-y2) - r2*r2;
  double f3 = (x-x3)*(x-x3) + (y-y3)*(y-y3) - r3*r3;

  *f = PDM_MIN(PDM_MIN(f1, f2), f3);

  double df_dx1 = 2*(x-x1);
  double df_dy1 = 2*(y-y1);
  double df_dx2 = 2*(x-x2);
  double df_dy2 = 2*(y-y2);
  double df_dx3 = 2*(x-x3);
  double df_dy3 = 2*(y-y3);

  if (f1 < f2) {
    if (f1 < f3) {
      *df_dx = df_dx1;
      *df_dy = df_dy1;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  } else {
    if (f2 < f3) {
      *df_dx = df_dx2;
      *df_dy = df_dy2;
    } else {
      *df_dx = df_dx3;
      *df_dy = df_dy3;
    }
  }

  *df_dz = 0.;
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
  double k = 2. * PDM_PI / L;
  double F = exp(-2.*nu*k*k*t0);

  double q = F*2*pow(k,2)*pow(U0,2)*pow(sin(k*x),2)*pow(sin(k*y),2)+2*pow(k,2)*pow(U0,2)*pow(cos(k*x),2)*pow(cos(k*y),2);

  *f = q - 20. * pow(L/U0, 2);

  *df_dx = F*4*pow(k,3)*pow(U0,2)*sin(k*x)*pow(sin(k*y),2)*cos(k*x)-4*pow(k,3)*pow(U0,2)*sin(k*x)*cos(k*x)*pow(cos(k*y),2);
  *df_dy = F*4*pow(k,3)*pow(U0,2)*pow(sin(k*x),2)*sin(k*y)*cos(k*y)-4*pow(k,3)*pow(U0,2)*sin(k*y)*pow(cos(k*x),2)*cos(k*y);
  *df_dz = 0.;
}



static const int it_max = 40;
static void
_eval_mandelbrot
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
  PDM_UNUSED(df_dx);
  PDM_UNUSED(df_dy);
  PDM_UNUSED(df_dz);

  double _x = x - 0.5;
  double _y = y;

  double xk = 0;
  double yk = 0;

  double zxk = 1.;
  double zyk = 0.;
  double rk;

  int it;
  for (it = 0; it < it_max; it++) {
    double xk_new = xk*xk - yk*yk + _x;
    double yk_new = 2*xk*yk       + _y;
    xk = xk_new;
    yk = yk_new;
    rk = sqrt(xk*xk + yk*yk);

    double zxk_new = 2*(xk*zxk - yk*zyk) + 1;
    double zyk_new = 2*(zxk*yk + xk*zyk) + 1;

    zxk = zxk_new;
    zyk = zyk_new;

    if (rk > 2.) {
      break;
    }
  }

  double mag_dz = sqrt(zxk*zxk + zyk*zyk);
  *f = rk * log(PDM_MAX(1e-9, rk)) / PDM_MAX(1e-9, mag_dz);
  // if (it < it_max-1) {
  //   *f = 1.;
  // } else {
  //   *f = -1.;
  // }
}


static const int it_max2 = 6;
static const int power  = 8;
static void
_eval_mandelbulb_slice
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
  PDM_UNUSED(z);

  double dr = 1.;
  double r;
  double theta;
  double phi;

  double xk = x;
  double yk = y;
  double zk = z;

  for (int it = 0; it < it_max2; it++) {

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
  double a = (2*z*z + x*x + y*y - 1);
  *f = a*a*a - (0.1*z*z + x*x)*y*y*y;

  *df_dz = -0.2*z*y*y*y + 12*z*a*a;
  *df_dx = -2*x*y*y*y + 6*x*a*a;
  *df_dy = -(0.3*z*z + 3*x*x)*y*y + 6*y*a*a;
}


#define PERLIN_NOISE_N 8
static double perlin_fx[PERLIN_NOISE_N][PERLIN_NOISE_N];
static double perlin_fy[PERLIN_NOISE_N][PERLIN_NOISE_N];
static const double perlin_step = 1. / (double) (PERLIN_NOISE_N - 1);

static void
_init_perlin_noise (void)
{
  double i_rand_max = 1./(double) RAND_MAX;

  for (int j = 0; j < PERLIN_NOISE_N; j++) {
    for (int i = 0; i < PERLIN_NOISE_N; i++) {
      perlin_fx[i][j] = 2*rand()*i_rand_max - 1;
      perlin_fy[i][j] = 2*rand()*i_rand_max - 1;
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
  PDM_UNUSED(z);
  PDM_UNUSED(df_dx);
  PDM_UNUSED(df_dy);
  PDM_UNUSED(df_dz);

  int pi = (int) (x / perlin_step);
  int pj = (int) (y / perlin_step);

  pi = PDM_MIN(PDM_MAX(pi, 0), PERLIN_NOISE_N-2);
  pj = PDM_MIN(PDM_MAX(pj, 0), PERLIN_NOISE_N-2);

  double x0 = pi*perlin_step;
  double y0 = pj*perlin_step;

  double x1 = x0 + perlin_step;
  double y1 = y0 + perlin_step;

  double sx = (x - x0) / perlin_step;
  double sy = (y - y0) / perlin_step;

  double n0, n1;
  n0 = (x - x0)*perlin_fx[pi][pj]   + (y - y0)*perlin_fy[pi][pj];
  n1 = (x - x1)*perlin_fx[pi+1][pj] + (y - y0)*perlin_fy[pi+1][pj];
  double ix0 = _interpolate_perlin_noise(n0, n1, sx);

  n0 = (x - x0)*perlin_fx[pi][pj+1]   + (y - y1)*perlin_fy[pi][pj+1];
  n1 = (x - x1)*perlin_fx[pi+1][pj+1] + (y - y1)*perlin_fy[pi+1][pj+1];
  double ix1 = _interpolate_perlin_noise(n0, n1, sx);

  *f = _interpolate_perlin_noise(ix0, ix1, sy);
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
  PDM_Mesh_nodal_elt_t elt_type  = PDM_MESH_NODAL_TRIA3;
  int                  n_part    = -1;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method   = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif
  //  2 -> tria
  //  3 -> quad
  //  4 -> poly2d

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

  assert(PDM_Mesh_nodal_elt_dim_get(elt_type) == 2);

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

  eval_field_and_gradient = &_eval_mickey;
  // _init_perlin_noise();
  // eval_field_and_gradient = &_eval_perlin_noise;


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

  PDM_g_num_t *distrib_face = NULL;
  PDM_g_num_t *distrib_edge = NULL;
  PDM_g_num_t *vtx_distrib  = NULL;

  PDM_dcube_nodal_t          *dcube = NULL;
  PDM_dmesh_nodal_to_dmesh_t *dmntodm = NULL;

  PDM_g_num_t     ng_face         = 0;
  PDM_g_num_t     ng_vtx          = 0;
  PDM_g_num_t     ng_edge         = 0;
  int             dn_vtx          = 0;
  int             dn_face         = 0;
  int             dn_edge         = 0;
  int             n_edge_group    = 0;
  double         *dvtx_coord      = NULL;
  int            *dface_edge_idx  = NULL;
  PDM_g_num_t    *dface_vtx       = NULL;
  PDM_g_num_t    *dface_edge      = NULL;
  int            *dedge_vtx_idx   = NULL;
  PDM_g_num_t    *dedge_vtx       = NULL;
  PDM_g_num_t    *dedge_face      = NULL;
  int            *dedge_group_idx = NULL;
  PDM_g_num_t    *dedge_group     = NULL;


  if (elt_type < PDM_MESH_NODAL_POLY_2D) {

    dcube = PDM_dcube_nodal_gen_create (comm,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        length,
                                        -0.5*length,// -0.35,
                                        -0.5*length,// -0.3,
                                        0.,
                                        elt_type,
                                        1,
                                        PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build (dcube);


    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
    PDM_dmesh_nodal_generate_distribution(dmn);

    vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];


    if (0) {
      double noise = 0.2 * length / (double) (n_vtx_seg - 1);
      for (int i = 0; i < dn_vtx; i++) {
        for (int j = 0; j < 2; j++) {
          dvtx_coord[3*i+j] += noise * (rand() / (double) RAND_MAX - 0.5);
        }
      }
    }

    if(1 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
      PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    }

    dmntodm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);

    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmntodm, 0, dmn);

    PDM_dmesh_nodal_to_dmesh_set_post_treat_result(dmntodm, 1);

    PDM_dmesh_nodal_to_dmesh_compute(dmntodm,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                     PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);

    PDM_dmesh_t* dmesh = NULL;
    PDM_dmesh_nodal_to_dmesh_get_dmesh(dmntodm, 0, &dmesh);


    dn_face = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         &dface_edge,
                                         &dface_edge_idx,
                                         PDM_OWNERSHIP_KEEP);

    dn_edge  = PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &dedge_vtx,
                                          &dedge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

    if(0 == 1) {
      PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dn_face, "dface_edge ::");
      PDM_log_trace_connectivity_long(dedge_vtx_idx , dedge_vtx , dn_edge, "dedge_vtx  ::");
    }

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE  , &distrib_edge);
    assert(distrib_edge != NULL);

    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE  , &distrib_face);
    assert(distrib_face != NULL);
  }

  else {

    PDM_poly_surf_gen (comm,
                       -0.5*length,
                       0.5*length,
                       -0.5*length,
                       0.5*length,
                       1,
                       0,
                       n_vtx_seg,
                       n_vtx_seg,
                       &ng_face,
                       &ng_vtx,
                       &ng_edge,
                       &dn_vtx,
                       &dvtx_coord,
                       &dn_face,
                       &dface_edge_idx,
                       &dface_vtx,
                       &dface_edge,
                       &dn_edge,
                       &dedge_vtx,
                       &dedge_face,
                       &n_edge_group,
                       &dedge_group_idx,
                       &dedge_group);

    // PDM_log_trace_connectivity_long(dface_edge_idx,
    //                                 dface_edge,
    //                                 dn_face,
    //                                 "dface_edge : ");
    // PDM_log_trace_connectivity_long(dface_edge_idx,
    //                                 dface_vtx,
    //                                 dn_face,
    //                                 "dface_vtx : ");

    // _dump_dmesh (comm,
    //              dn_face,
    //              dn_vtx,
    //              dface_edge_idx,
    //              dface_vtx,
    //              dvtx_coord);

    distrib_face = PDM_compute_entity_distribution(comm, dn_face);
    distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
    vtx_distrib  = PDM_compute_entity_distribution(comm, dn_vtx);

    dedge_vtx_idx = (int *) malloc(sizeof(int) * (dn_edge + 1));
    for (int i = 0; i <= dn_edge; i++) {
      dedge_vtx_idx[i] = 2*i;
    }

    if (n_part > 0) {
      int n_bound = n_edge_group;
      int n_join = 0;

      int *djoins_ids     = (int *) malloc(n_join * sizeof(int));
      int *dedge_bnd_idx  = (int *) malloc((n_bound + 1) * sizeof(int));
      int *dedge_join_idx = (int *) malloc((n_join  + 1) * sizeof(int));
      dedge_bnd_idx[0] = 0;
      dedge_join_idx[0] = 0;

      // First pass to count and allocate
      int i_bnd = 1;
      for (int igroup = 0; igroup < n_edge_group; igroup++) {
        int group_size = dedge_group_idx[igroup+1] - dedge_group_idx[igroup];
        dedge_bnd_idx[i_bnd++] = group_size;
      }
      for (int i = 0; i < n_bound; i++) {
        dedge_bnd_idx[i+1] = dedge_bnd_idx[i+1] + dedge_bnd_idx[i];
      }

      // Second pass to copy
      PDM_g_num_t *dedge_bnd  = (PDM_g_num_t *) malloc(dedge_bnd_idx[n_bound] * sizeof(PDM_g_num_t));
      PDM_g_num_t *dedge_join = (PDM_g_num_t *) malloc(dedge_join_idx[n_join] * sizeof(PDM_g_num_t));

      i_bnd = 0;
      for (int igroup = 0; igroup < n_edge_group; igroup++) {
        for (int i = dedge_group_idx[igroup]; i < dedge_group_idx[igroup+1]; i++) {
          dedge_bnd[i_bnd++] = dedge_group[i];
        }
      }
      PDM_dmesh_t *dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                             dn_face,
                                             dn_edge,
                                             -1,
                                             dn_vtx,
                                             n_bound,
                                             n_join,
                                             comm);

      PDM_dmesh_set (dmesh,
                     dvtx_coord,
                     dedge_vtx_idx,
                     dedge_vtx,
                     dedge_face,
                     dedge_bnd_idx,
                     dedge_bnd,
                     djoins_ids,
                     dedge_join_idx,
                     dedge_join);

      PDM_multipart_register_block (mpart, 0, dmesh);

      /* Connection between zones */
      int n_total_joins = 0;
      int *join_to_opposite = (int *) malloc(n_total_joins*sizeof(int));
      PDM_multipart_register_joins (mpart, n_total_joins, join_to_opposite);

      /* Run */
      PDM_multipart_run_ppart (mpart);
    }
  }




  PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(2, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

  PDM_iso_surface_eval_field_and_gradient_set(isos,
                                              eval_field_and_gradient);


  double  *dfield          = NULL;
  double  *dgradient_field = NULL;
  double **pfield          = NULL;
  double **pgradient_field = NULL;

  if (n_part > 0) {
    assert(elt_type == PDM_MESH_NODAL_POLY_2D);

    pfield          = (double **) malloc(sizeof(double *) * n_part);
    pgradient_field = (double **) malloc(sizeof(double *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);
      printf("vtx: %d %p\n", n_vtx, (void *) vtx_coord);

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



      int *face_edge_idx;
      int *face_edge;
      int *edge_vtx_idx;
      int *edge_vtx;

      int n_proc, tn_part;
      int _n_vtx, n_bounds, n_joins, n_part_joins;
      int sface_edge, sedge_vtx, sedge_bound, sedge_join;
      int  n_section;
      int* n_elt;

      int n_face, n_edge;
      PDM_multipart_part_dim_get(mpart, 0, i_part, &n_section, &n_elt,
                                 &n_face, &n_edge, &n_part_joins, &_n_vtx, &n_proc, &tn_part,
                                 &sface_edge, &sedge_vtx, &sedge_bound, &n_bounds, &sedge_join, &n_joins);

      double       *_vtx;
      int          *_edge_face;
      int          *edge_bound_idx, *edge_bound, *edge_join_idx, *edge_join;
      int          *edge_part_bound_proc_idx, *edge_part_bound_part_idx, *edge_part_bound;
      PDM_g_num_t  *_face_ln_to_gn, *edge_ln_to_gn, *_vtx_ln_to_gn, *edge_bound_ln_to_gn, *edge_join_ln_to_gn;
      int          *face_tag, *edge_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart, 0, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &face_tag, &face_edge_idx, &face_edge, &_face_ln_to_gn,
                                 &edge_tag, &_edge_face, &edge_vtx_idx, &edge_vtx, &edge_ln_to_gn,
                                 &edge_part_bound_proc_idx, &edge_part_bound_part_idx, &edge_part_bound,
                                 &vtx_tag, &_vtx, &_vtx_ln_to_gn, &edge_bound_idx, &edge_bound,
                                 &edge_bound_ln_to_gn, &edge_join_idx, &edge_join, &edge_join_ln_to_gn);

      // int n_face = PDM_multipart_part_connectivity_get(mpart,
      //                                                  0,
      //                                                  i_part,
      //                                                  PDM_CONNECTIVITY_TYPE_FACE_EDGE,
      //                                                  &face_edge,
      //                                                  &face_edge_idx,
      //                                                  PDM_OWNERSHIP_KEEP);
      printf("face_edge: %d %p %p\n", n_face, (void *) face_edge_idx, (void *) face_edge);

      // int n_edge = PDM_multipart_part_connectivity_get(mpart,
      //                                                  0,
      //                                                  i_part,
      //                                                  PDM_CONNECTIVITY_TYPE_EDGE_VTX,
      //                                                  &edge_vtx,
      //                                                  &edge_vtx_idx,
      //                                                  PDM_OWNERSHIP_KEEP);
      printf("edge_vtx: %d %p %p\n", n_edge, (void *) face_edge_idx, (void *) face_edge);

      PDM_g_num_t *face_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      printf("face_ln_to_gn : %p\n", (void *) face_ln_to_gn);

      // PDM_g_num_t *edge_ln_to_gn;
      // PDM_multipart_part_ln_to_gn_get(mpart,
      //                                 0,
      //                                 i_part,
      //                                 PDM_MESH_ENTITY_EDGE,
      //                                 &edge_ln_to_gn,
      //                                 PDM_OWNERSHIP_KEEP);
      printf("edge_ln_to_gn : %p\n", (void *) edge_ln_to_gn);

      PDM_g_num_t *vtx_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      printf("vtx_ln_to_gn : %p\n", (void *) vtx_ln_to_gn);


      PDM_iso_surface_part_set (isos,
                                i_part,
                                0,//n_cell,
                                n_face,
                                n_edge,
                                n_vtx,
                                NULL,//cell_face_idx,
                                NULL,//cell_face,
                                face_edge_idx,
                                face_edge,
                                edge_vtx,
                                NULL,//face_vtx_idx,
                                NULL,//face_vtx,
                                NULL,//cell_ln_to_gn,
                                face_ln_to_gn,
                                edge_ln_to_gn,
                                vtx_ln_to_gn,
                                vtx_coord);

      PDM_iso_surface_part_field_set (isos,
                                      i_part,
                                      pfield[i_part]);

      PDM_iso_surface_part_gradient_field_set (isos,
                                               i_part,
                                               pgradient_field[i_part]);
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

    if (elt_type == PDM_MESH_NODAL_POLY_2D) {
      _dump_dmesh (comm,
                   dn_face,
                   dn_vtx,
                   dface_edge_idx,
                   dface_vtx,
                   dvtx_coord,
                   dfield,
                   dgradient_field);
    }

    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      dface_edge,
                                      dface_edge_idx);
    PDM_iso_surface_dconnectivity_set(isos,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      dedge_vtx,
                                      dedge_vtx_idx);

    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_FACE  , distrib_face);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_EDGE  , distrib_edge);
    PDM_iso_surface_distrib_set(isos, PDM_MESH_ENTITY_VERTEX, vtx_distrib);



    PDM_iso_surface_dvtx_coord_set (isos, dvtx_coord     );
    PDM_iso_surface_dfield_set     (isos, dfield         );
    PDM_iso_surface_dgrad_field_set(isos, dgradient_field);
  }





  PDM_iso_surface_compute(isos);

  char name[999];
  sprintf(name, "iso_line_%dproc", n_rank);
  PDM_iso_surface_write(isos, name);


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


  if (elt_type < PDM_MESH_NODAL_POLY_2D) {
    PDM_dmesh_nodal_to_dmesh_free(dmntodm);
    PDM_dcube_nodal_gen_free(dcube);
  } else {
    free(distrib_face);
    free(distrib_edge);
    free(vtx_distrib);
    free(dvtx_coord);
    free(dface_edge_idx);
    free(dface_vtx);
    free(dface_edge);
    free(dedge_vtx_idx);
    free(dedge_vtx);
    free(dedge_face);
    free(dedge_group_idx);
    free(dedge_group);
  }

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
