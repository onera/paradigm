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
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_priv.h"
#include "pdm_writer.h"

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





static void
_eval_smiley
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
  double x0 = 0.;
  double y0 = 0.;
  double z0 = 0.;
  double r0 = 0.48;

  // left eye
  double x1 = -0.18;
  double y1 = 0.32;
  double z1 = 0.18;
  double r1 = 0.11;

  // right eye
  double x2 = -x1;
  double y2 = y1;
  double z2 = z1;
  double r2 = r1;

  // tear drop
  double x3 = -0.37;
  double y3 = 0.3;
  double z3 = 0.36;
  double a3 = 0.055;
  double b3 = 0.055;
  double c3 = 0.13;
  double ft = (1. + 1.*(z-z3))*((x-x3)*(x-x3)/(a3*a3) + (y-y3)*(y-y3)/(b3*b3)) + (z-z3)*(z-z3)/(c3*c3) - 1;
  double dft_dx = 2*(1. + 1.*(z-z3))*(x-x3)/(a3*a3);
  double dft_dy = 2*(1. + 1.*(z-z3))*(y-y3)/(b3*b3);
  double dft_dz = (x-x3)*(x-x3)/(a3*a3) + (y-y3)*(y-y3)/(b3*b3) + 2*(z-z3)/(c3*c3);

  // mouth
  double fm1 = -0.2*y + z + 0.18;
  double fm2 = -(y + z) + 0.13;
  double fm, dfm_dx, dfm_dy, dfm_dz;
  if (fm1 > fm2) {
    fm = -fm1;
    dfm_dx = 0.;
    dfm_dy = 0.;
    dfm_dz = -1.;
  } else {
    fm = -fm2;
    dfm_dx = 0.;
    dfm_dy = 1.;
    dfm_dz = 1.;
  }

  double F[4] = {
    (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) - r0*r0,
    -((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1) - r1*r1),
    -((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2) - r2*r2),
    fm
  };

  double dF_dx[4] = {
    2*(x-x0),
    -2*(x-x1),
    -2*(x-x2),
    dfm_dx
  };

  double dF_dy[4] = {
    2*(y-y0),
    -2*(y-y1),
    -2*(y-y2),
    dfm_dy
  };

  double dF_dz[4] = {
    2*(z-z0),
    -2*(z-z1),
    -2*(z-z2),
    dfm_dz
  };

  *f = -HUGE_VAL;
  for (int i = 0; i < 4; i++) {
    if (F[i] > (*f)) {
      *f     = F[i];
      *df_dx = dF_dx[i];
      *df_dy = dF_dy[i];
      *df_dz = dF_dz[i];
    }
  }

  if (*f > ft) {
    *f = ft;
    *df_dx = dft_dx;
    *df_dy = dft_dy;
    *df_dz = dft_dz;
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
  *f = x * x + y * y + z * z - 0.26;

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
static const int it_max = 12;
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


static void
_eval_helicoid
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
  const double h = 1.;
  double t = tan(z/h);

  *f = x*t - y;

  *df_dx = t;
  *df_dy = 1;
  *df_dz = (1 + t*t)*x/h;
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

static void
_update_perlin_noise (const double dt)
{
  const double c = cos(dt);
  const double s = sin(dt);

  for (int k = 0; k < PERLIN_NOISE_N; k++) {
    for (int j = 0; j < PERLIN_NOISE_N; j++) {
      for (int i = 0; i < PERLIN_NOISE_N; i++) {
        double x = perlin_fx[i][j][k];
        double y = perlin_fy[i][j][k];
        double z = perlin_fz[i][j][k];

        perlin_fx[i][j][k] =  c*x + s*y;
        perlin_fy[i][j][k] = -s*x + c*y;
        perlin_fz[i][j][k] = z;
      }
    }
  }
}


static double _smooth_step (const double x)
{
  return x*(3*x - 2*x*x);
}

static double _smooth_step_deriv (const double x)
{
  return 6*x*(1 - x);
}

static double _lerp (const double a, const double b, const double x)
{
  return (1 - x)*a + x*b;
}

static double _dot (
                    const double x, const double y, const double z,
                    const int i, const int j, const int k
                    )
{
  return
  x*perlin_fx[i][j][k] +
  y*perlin_fy[i][j][k] +
  z*perlin_fz[i][j][k];
}

static double
_eval_perlin
(
 const double  x,
 const double  y,
 const double  z
 )
 {
  // #define _smooth_step (x) ((x)*(3*(x) - 2*(x)*(x)))
  // #define _lerp (a,b,x) ((1 - (x))*(a) + (x)*(b))
  // #define _dot (x,y,z,i,j,k) ((x)*perlin_fx[(i)][(j)][(k)] + (y)*perlin_fy[(i)][(j)][(k)] + (z)*perlin_fz[(i)][(j)][(k)])

  int i = (int) (x / perlin_step);
  int j = (int) (y / perlin_step);
  int k = (int) (z / perlin_step);

  i = PDM_MIN(PDM_MAX(i, 0), PERLIN_NOISE_N-2);
  j = PDM_MIN(PDM_MAX(j, 0), PERLIN_NOISE_N-2);
  k = PDM_MIN(PDM_MAX(k, 0), PERLIN_NOISE_N-2);

  double tx = x - i*perlin_step;
  double ty = y - j*perlin_step;
  double tz = z - k*perlin_step;

  double u = _smooth_step(tx/perlin_step);
  double v = _smooth_step(ty/perlin_step);
  double w = _smooth_step(tz/perlin_step);

  double ux = _smooth_step_deriv(tx/perlin_step)/perlin_step;
  double vy = _smooth_step_deriv(ty/perlin_step)/perlin_step;
  double wz = _smooth_step_deriv(tz/perlin_step)/perlin_step;
  PDM_UNUSED(ux);
  PDM_UNUSED(vy);
  PDM_UNUSED(wz);


  double x0 = tx;
  double y0 = ty;
  double z0 = tz;
  double x1 = tx - perlin_step;
  double y1 = ty - perlin_step;
  double z1 = tz - perlin_step;


  double d000 = _dot(x0, y0, z0, i,   j,   k  );
  double d100 = _dot(x1, y0, z0, i+1, j,   k  );
  double d010 = _dot(x0, y1, z0, i,   j+1, k  );
  double d110 = _dot(x1, y1, z0, i+1, j+1, k  );
  double d001 = _dot(x0, y0, z1, i,   j,   k+1);
  double d101 = _dot(x1, y0, z1, i+1, j,   k+1);
  double d011 = _dot(x0, y1, z1, i,   j+1, k+1);
  double d111 = _dot(x1, y1, z1, i+1, j+1, k+1);

  double a = _lerp(d000, d100, u);
  double b = _lerp(d010, d110, u);
  double c = _lerp(d001, d101, u);
  double d = _lerp(d011, d111, u);

  double e = _lerp(a, b, v);
  double f = _lerp(c, d, v);

  double val = _lerp(e, f, w) - 0.041;

  // #undef _smooth_step
  // #undef _lerp
  // #undef _dot

  return val;
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
  const double h = 1e-4;

  *f = _eval_perlin(x, y, z);
  *df_dx = (_eval_perlin(x + h, y,     z    ) - (*f)) / h;
  *df_dy = (_eval_perlin(x,     y + h, z    ) - (*f)) / h;
  *df_dz = (_eval_perlin(x,     y,     z + h) - (*f)) / h;
}



static void
_generate_dedges
(
 PDM_MPI_Comm  comm,
 int           dn_face,
 int          *dface_vtx_idx,
 PDM_g_num_t  *dface_vtx,
 int           dn_vtx,
 int          *dn_edge,
 PDM_g_num_t **dedge_vtx,
 int         **dface_edge_idx,
 PDM_g_num_t **dface_edge
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t *face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  PDM_g_num_t *vtx_distribution = PDM_compute_entity_distribution(comm, dn_vtx);

  /*
   * Generate edge numbering
   */
  int n_edge_elt_tot = dface_vtx_idx[dn_face];
  PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
  int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
  PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

  int n_elmt_current = 0;
  int n_edge_current = 0;
  tmp_dface_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(dn_face,
                              &n_elmt_current,
                              &n_edge_current,
                              face_distribution[i_rank],
                              -1,
                              dface_vtx,
                              dface_vtx_idx,
                              tmp_dface_edge_vtx_idx,
                              tmp_dface_edge_vtx,
                              tmp_dface_edge,
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_elt_tot);
  free(tmp_parent_elmt_pos);

  // int  dn_edge = -1;
  PDM_g_num_t  *dedge_distrib;
  int          *dedge_vtx_idx;
  // PDM_g_num_t  *dedge_vtx;
  int          *dedge_face_idx;
  PDM_g_num_t  *dedge_face;

  PDM_generate_entitiy_connectivity_raw(comm,
                                        vtx_distribution[n_rank],
                                        n_edge_elt_tot,
                                        tmp_dface_edge,
                                        tmp_dface_edge_vtx_idx,
                                        tmp_dface_edge_vtx,
                                        dn_edge,
                                        &dedge_distrib,
                                        &dedge_vtx_idx,
                                        dedge_vtx,
                                        &dedge_face_idx,
                                        &dedge_face);

  /* Make ascending connectivity */
  // int          *dface_edge_idx;
  // PDM_g_num_t  *dface_edge;
  // PDM_log_trace_array_long(face_distribution, n_rank+1, "face_distribution::");
  // PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distrib::");
  // PDM_log_trace_array_int(dedge_vtx_idx, (*dn_edge)+1, "dedge_vtx_idx::");
  // PDM_log_trace_array_long(*dedge_vtx, dedge_vtx_idx[(*dn_edge)], "dedge_vtx::");
  // PDM_log_trace_array_int(dedge_face_idx, (*dn_edge), "dedge_face_idx::");
  // PDM_log_trace_array_long(dedge_face, dedge_face_idx[(*dn_edge)], "dedge_face::");

  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              face_distribution,
                              dedge_face_idx,
                              dedge_face,
                              1,
                              dface_edge_idx,
                              dface_edge);

 free(dedge_distrib);
 free(dedge_vtx_idx);
 free(dedge_face_idx);
 free(dedge_face);
 free(face_distribution);
 free(vtx_distribution);
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
_dump_dmesh
(
 int         dn_cell,
 int         dcell_face_idx[],
 PDM_g_num_t dcell_face[],
 int         dn_face,
 int         dface_vtx_idx[],
 PDM_g_num_t dface_vtx[],
 int         dn_vtx,
 double      dvtx_coord[],
 double      dfield[],
 double      dgradient[]
 )
{
  PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                          PDM_WRITER_FMT_ASCII,
                                          PDM_WRITER_TOPO_CONSTANTE,
                                          PDM_WRITER_OFF,
                                          "test_isosurf3d",
                                          "isosurf3d",
                                          PDM_MPI_COMM_WORLD,
                                          PDM_IO_ACCES_MPI_SIMPLE,
                                          1.,
                                          NULL);

  int id_geom = PDM_writer_geom_create(id_cs,
                                       "isosurf3d_geom",
                                       PDM_WRITER_OFF,
                                       PDM_WRITER_OFF,
                                       1);

  int id_var_field = PDM_writer_var_create(id_cs,
                                           PDM_WRITER_ON,
                                           PDM_WRITER_VAR_SCALAIRE,
                                           PDM_WRITER_VAR_SOMMETS,
                                           "field");

  int id_var_gradient;

  if (dgradient != NULL) {
    id_var_gradient = PDM_writer_var_create(id_cs,
                                            PDM_WRITER_ON,
                                            PDM_WRITER_VAR_VECTEUR,
                                            PDM_WRITER_VAR_SOMMETS,
                                            "gradient");
  }

  PDM_writer_step_beg(id_cs, 0.);

  PDM_g_num_t *cell_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_cell);
  for (int i = 0; i < dn_cell; i++) {
    cell_ln_to_gn[i] = i + 1;
  }

  PDM_g_num_t *face_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_face);
  for (int i = 0; i < dn_face; i++) {
    face_ln_to_gn[i] = i + 1;
  }

  PDM_g_num_t *vtx_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_vtx);
  for (int i = 0; i < dn_vtx; i++) {
    vtx_ln_to_gn[i] = i + 1;
  }


  int *cell_face = (int *) malloc(sizeof(int) * dcell_face_idx[dn_cell]);
  for (int i = 0; i < dcell_face_idx[dn_cell]; i++) {
    cell_face[i] = (int) dcell_face[i];
  }

  int *face_vtx = (int *) malloc(sizeof(int) * dface_vtx_idx[dn_face]);
  for (int i = 0; i < dface_vtx_idx[dn_face]; i++) {
    face_vtx[i] = (int) dface_vtx[i];
  }



  int *face_vtxNb  = (int *) malloc(sizeof(int) * dn_face);
  int *cell_faceNb = (int *) malloc(sizeof(int) * dn_cell);

  for (int i = 0; i < dn_cell; i++) {
    cell_faceNb[i] = dcell_face_idx[i+1] - dcell_face_idx[i];
  }

  for (int i = 0; i < dn_face; i++) {
    face_vtxNb[i] = dface_vtx_idx[i+1] - dface_vtx_idx[i];
  }

  PDM_writer_geom_coord_set(id_cs,
                            id_geom,
                            0,
                            dn_vtx,
                            dvtx_coord,
                            vtx_ln_to_gn);

  PDM_writer_geom_cell3d_cellface_add (id_cs,
                                       id_geom,
                                       0,
                                       dn_cell,
                                       dn_face,
                                       dface_vtx_idx,
                                       face_vtxNb,
                                       face_vtx,
                                       dcell_face_idx,
                                       cell_faceNb,
                                       cell_face,
                                       cell_ln_to_gn);


  PDM_writer_geom_write(id_cs,
                        id_geom);


  PDM_real_t *val_field    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * dn_vtx);
  PDM_real_t *val_gradient = (PDM_real_t *) malloc(sizeof(PDM_real_t) * dn_vtx * 3);

  for (int i = 0; i < dn_vtx; i++) {
    val_field[i] = (PDM_real_t) dfield[i];
    if (dgradient != NULL) {
      for (int j = 0; j < 3; j++) {
        val_gradient[3*i+j] = (PDM_real_t) dgradient[3*i+j];
      }
    }
  }


  PDM_writer_var_set(id_cs,
                     id_var_field,
                     id_geom,
                     0,
                     val_field);

  PDM_writer_var_write(id_cs,
                       id_var_field);


  if (dgradient != NULL) {
    PDM_writer_var_set(id_cs,
                       id_var_gradient,
                       id_geom,
                       0,
                       val_gradient);

    PDM_writer_var_write(id_cs,
                         id_var_gradient);
  }


  PDM_writer_step_end(id_cs);

  PDM_writer_free(id_cs);

  free(val_field);
  free(val_gradient);
  free(face_vtxNb);
  free(cell_faceNb);

  free(cell_ln_to_gn);
  free(face_ln_to_gn);
  free(vtx_ln_to_gn);
  free(cell_face);
  free(face_vtx);
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
  PDM_dmesh_t                *dmesh   = NULL;

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
  int          n_face_group    = 0;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  if (elt_type < PDM_MESH_NODAL_POLY_3D) {
    dcube = PDM_dcube_nodal_gen_create (comm,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        length,
                                        -0.5*length,//-0.1,
                                        -0.5*length,//-0.2,
                                        -0.5*length,//0.07,
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
                      0.,//-0.5*length,//-0.1,
                      0.,//-0.5*length,//-0.2,
                      0.,//-0.5*length,//0.07,
                      length,
                      length,
                      length,
                      n_vtx_seg,
                      n_vtx_seg,
                      n_vtx_seg,
                      1, // randomize
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
                      &dface_vtx_idx,
                      &dface_vtx,
                      &dvtx_coord,
                      &dface_group_idx,
                      &dface_group);

    _generate_dedges(comm,
                     dn_face,
                     dface_vtx_idx,
                     dface_vtx,
                     dn_vtx,
                     &dn_edge,
                     &dedge_vtx,
                     &dface_edge_idx,
                     &dface_edge);

    distrib_cell = PDM_compute_entity_distribution(comm, dn_cell);
    distrib_face = PDM_compute_entity_distribution(comm, dn_face);
    distrib_edge = PDM_compute_entity_distribution(comm, dn_edge);
    vtx_distrib  = PDM_compute_entity_distribution(comm, dn_vtx);


    if (n_part > 0) {
      /* Generate dmesh */
      int n_join = 0;
      dmesh = PDM_dmesh_create (PDM_OWNERSHIP_KEEP,
                                dn_cell,
                                dn_face,
                                dn_edge,
                                dn_vtx,
                                n_face_group,
                                n_join,
                                comm);

      int *djoins_ids = malloc (sizeof(int) * n_join);
      int *dface_join_idx = malloc (sizeof(int) * (n_join + 1));
      dface_join_idx[0] = 0;
      PDM_g_num_t *dface_join = malloc (sizeof(PDM_g_num_t) * dface_join_idx[n_join]);

      PDM_dmesh_set (dmesh,
                     dvtx_coord,
                     dface_vtx_idx,
                     dface_vtx,
                     dface_cell,
                     dface_group_idx,
                     dface_group,
                     djoins_ids,
                     dface_join_idx,
                     dface_join);

      //-->>
      dmesh->dn_edge    = dn_edge;
      dmesh->_dedge_vtx = dedge_vtx;
      //<<--

      PDM_multipart_register_block (mpart, 0, dmesh);

      /* Connection between zones */
      int n_total_joins = 0;
      int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
      PDM_multipart_register_joins (mpart, n_total_joins, join_to_opposite);
    }
  }


  if (n_part > 0) {
    PDM_multipart_run_ppart(mpart);
  }

<<<<<<< HEAD
  if(0 == 1) {
    _init_perlin_noise();
    _update_perlin_noise(1.);
  }
  if(0 == 1) {
    eval_field_and_gradient = &_eval_helicoid;
    eval_field_and_gradient = &_eval_chmutov6;
    eval_field_and_gradient = &_eval_mcmullen;
    eval_field_and_gradient = &_eval_pretzel;
    eval_field_and_gradient = &_eval_taylor_green_vortex;
    eval_field_and_gradient = &_eval_heart;
    eval_field_and_gradient = &_eval_cylinder;
    eval_field_and_gradient = &_eval_sphere;
    eval_field_and_gradient = &_eval_smiley;
    eval_field_and_gradient = &_eval_mandelbulb;
    eval_field_and_gradient = &_eval_perlin_noise;
  }


  eval_field_and_gradient = &_eval_mandelbulb;
=======
  // _init_perlin_noise();
  eval_field_and_gradient = &_eval_sphere;
>>>>>>> 19dec9e8 ([pdm_iso_surface] add visu for polyhedral mesh)

  PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_FIELD, 1, PDM_OWNERSHIP_KEEP, comm);
  // PDM_iso_surface_t* isos = PDM_iso_surface_create(3, PDM_ISO_SURFACE_KIND_PLANE, 1, PDM_OWNERSHIP_KEEP, comm);

  PDM_iso_surface_plane_equation_set(isos, 1., 0., 0., -0);
  // PDM_iso_surface_plane_equation_set(isos, 1., 0.5, 0.25, -0.0234);

  PDM_iso_surface_eval_field_and_gradient_set(isos, eval_field_and_gradient);

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

      // int *cell_face_idx = NULL;
      // int *cell_face     = NULL;
      int *face_edge_idx = NULL;
      int *face_edge     = NULL;
      int *edge_vtx_idx  = NULL;
      int *edge_vtx      = NULL;

      int n_cell;
      int n_face;
      // int n_edge;
      // int n_face_part_bound;
      // int n_vtx;
      int n_proc;
      int n_t_part;
      // int s_cell_face;
      // int s_face_vtx;
      // int s_face_group;
      // int n_edge_group2;

      int n_bounds, n_joins, n_part_joins;
      int scell_face, sface_vtx, sface_bound, sface_join;
      int n_section;
      int *n_elt;

      int         *cell_tag;
      int         *cell_face_idx;
      int         *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int         *face_tag;
      int         *face_cell;
      int         *face_vtx_idx;
      int         *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int         *face_part_bound_proc_idx;
      int         *face_part_bound_part_idx;
      int         *face_part_bound;
      int         *vtx_tag;
      // double      *vtx_coord;
      PDM_g_num_t *vtx_ln_to_gn;
      // int         *face_group_idx;
      // int         *face_group;
      // PDM_g_num_t *face_group_ln_to_gn;
      PDM_g_num_t *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int         *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_dim_get (mpart,
                                  0,
                                  i_part,
                                  &n_section,
                                  &n_elt,
                                  &n_cell,
                                  &n_face,
                                  &n_part_joins,
                                  &n_vtx,
                                  &n_proc,
                                  &n_t_part,
                                  &scell_face,
                                  &sface_vtx,
                                  &sface_bound,
                                  &n_bounds,
                                  &sface_join,
                                  &n_joins);

       PDM_multipart_part_val_get (mpart,
                                    0,
                                    i_part,
                                    &elt_vtx_idx,
                                    &elt_vtx,
                                    &elt_section_ln_to_gn,
                                    &cell_tag,
                                    &cell_face_idx,
                                    &cell_face,
                                    &cell_ln_to_gn,
                                    &face_tag,
                                    &face_cell,
                                    &face_vtx_idx,
                                    &face_vtx,
                                    &face_ln_to_gn,
                                    &face_part_bound_proc_idx,
                                    &face_part_bound_part_idx,
                                    &face_part_bound,
                                    &vtx_tag,
                                    &vtx_coord,
                                    &vtx_ln_to_gn,
                                    &face_bound_idx,
                                    &face_bound,
                                    &face_bound_ln_to_gn,
                                    &face_join_idx,
                                    &face_join,
                                    &face_join_ln_to_gn);

      // int n_cell = PDM_multipart_part_connectivity_get(mpart,
      //                                                  0,
      //                                                  i_part,
      //                                                  PDM_CONNECTIVITY_TYPE_CELL_FACE,
      //                                                  &cell_face,
      //                                                  &cell_face_idx,
      //                                                  PDM_OWNERSHIP_KEEP);
      // printf("cell_face: %d %p %p\n", n_cell, (void *) cell_face_idx, (void *) cell_face);

      PDM_multipart_part_connectivity_get(mpart,
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

      // PDM_g_num_t *cell_ln_to_gn;
      // PDM_multipart_part_ln_to_gn_get(mpart,
      //                                 0,
      //                                 i_part,
      //                                 PDM_MESH_ENTITY_CELL,
      //                                 &cell_ln_to_gn,
      //                                 PDM_OWNERSHIP_KEEP);
      // printf("cell_ln_to_gn : %p\n", (void *) cell_ln_to_gn);

      // PDM_g_num_t *face_ln_to_gn;
      // PDM_multipart_part_ln_to_gn_get(mpart,
      //                                 0,
      //                                 i_part,
      //                                 PDM_MESH_ENTITY_FACE,
      //                                 &face_ln_to_gn,
      //                                 PDM_OWNERSHIP_KEEP);
      // printf("face_ln_to_gn : %p\n", (void *) face_ln_to_gn);

      PDM_g_num_t *edge_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);
      // printf("edge_ln_to_gn : %p\n", (void *) edge_ln_to_gn);

      // PDM_g_num_t *vtx_ln_to_gn;
      // PDM_multipart_part_ln_to_gn_get(mpart,
      //                                 0,
      //                                 i_part,
      //                                 PDM_MESH_ENTITY_VERTEX,
      //                                 &vtx_ln_to_gn,
      //                                 PDM_OWNERSHIP_KEEP);
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


      _dump_dmesh (dn_cell,
                   dcell_face_idx,
                   dcell_face,
                   dn_face,
                   dface_vtx_idx,
                   dface_vtx,
                   dn_vtx,
                   dvtx_coord,
                   dfield,
                   dgradient_field);
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
    PDM_iso_surface_dgrad_field_set(isos, dgradient_field);
  }

  PDM_MPI_Barrier(comm);
  PDM_iso_surface_compute(isos);
  PDM_MPI_Barrier(comm);

  char name[999];
  sprintf(name, "iso_surface_%dproc", n_rank);
  PDM_iso_surface_write(isos, name);

  PDM_iso_surface_dump_times(isos);

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
    if (dedge_vtx_idx != NULL) free(dedge_vtx_idx);
    free(dface_vtx_idx);
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
