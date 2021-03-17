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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_octree.h"
#include "pdm_para_octree.h"
#include "pdm_timer.h"
#include "pdm_part.h"
#include "pdm_mpi_node_first_rank.h"

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
     "  -n      <level>  Number of vertices.\n\n"
     "  -l      <level>  Rocket length.\n\n"
     "  -t      <level>  Number of Target points (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_faceSeg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   nTgt       Number of Target points
 * \param [inout]   n_part     Number of partitions par process
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nv,
           double        *length,
           double        *domain_size,
           PDM_g_num_t   *nTgt,
           int           *n_max_per_leaf,
           int           *n_methods,
           int           *on_ground,
           int           *n_proc_data_src)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-nv") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _nv = atol(argv[i]);
        *nv = (PDM_g_num_t) _nv;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
      }
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *domain_size = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-mpl") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_max_per_leaf = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_methods = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-ground") == 0) {
      *on_ground = 1;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_proc_data_src = atoi(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double _rand01(void) {
  return (double) rand() / (double) RAND_MAX;
}

static void
_gen_cloud_random
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_pts,
 const double        origin[3],
 const double        length,
 int                *_n_pts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  // Define distribution
  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  distrib[0] = 0;
  PDM_g_num_t step = n_pts / n_rank;
  PDM_g_num_t remainder = n_pts % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] = step;
    const int i1 = i - 1;
    if (i1 < remainder) {
      distrib[i]++;
    }
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] += distrib[i-1];
  }

  PDM_g_num_t dn_pts = distrib[i_rank+1] - distrib[i_rank];
  *_n_pts = (int) dn_pts;

  *g_num = malloc (sizeof(PDM_g_num_t) * dn_pts);
  *coord = malloc (sizeof(double)      * dn_pts * 3);
  for (int i = 0; i < *_n_pts; i++) {
    (*g_num)[i] = 1 + i + distrib[i_rank];
    for (int j = 0; j < 3; j++) {
      (*coord)[3*i+j] = origin[j] + length * _rand01();
    }
  }

  free (distrib);
}


static void
_gen_rocket
(
 PDM_MPI_Comm     comm,
 PDM_g_num_t      nv,
 double           length,
 double           radius,
 PDM_g_num_t     *ng_face,
 PDM_g_num_t     *ng_vtx,
 PDM_g_num_t     *ng_edge,
 int             *dn_vtx,
 double         **dvtx_coord,
 int             *dn_face,
 int            **dface_vtx_idx,
 PDM_g_num_t    **dface_vtx,
 PDM_g_num_t    **dface_edge,
 int             *dn_edge,
 PDM_g_num_t    **dedge_vtx,
 PDM_g_num_t    **dedge_face,
 int             *n_edge_group,
 int            **dedge_group_idx,
 PDM_g_num_t    **dedge_group
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  double length_cap = length / 14.;
  double length_cyl = length - length_cap;

  double C = sqrt(length_cap*length_cap + radius*radius);
  PDM_g_num_t _nu = (int) (0.5 * (nv - 1) * (PDM_PI*radius/3. + C) / length_cyl);
  _nu = PDM_MAX (_nu, 1);

  PDM_g_num_t nu = 6*_nu;

  PDM_g_num_t ng_vtx_cyl    = nu * (nv - 1);
  PDM_g_num_t ng_edge_cyl_u = nu * (nv - 1);
  PDM_g_num_t ng_edge_cyl_v = nu * (nv - 1);
  PDM_g_num_t ng_face_cyl   = nu * (nv - 1);

  PDM_g_num_t nv_cap      = _nu;
  PDM_g_num_t ng_vtx_cap  = 1 + 3*nv_cap*(nv_cap + 1);
  PDM_g_num_t ng_face_cap = 6 * nv_cap * nv_cap;
  PDM_g_num_t ng_edge_cap = ng_face_cap + ng_vtx_cap - 1;

  PDM_g_num_t ng_edge_lim = nu;

  *ng_vtx  = ng_vtx_cyl + ng_vtx_cap;
  *ng_edge = ng_edge_cyl_u + ng_edge_cyl_v + ng_edge_cap;
  *ng_face = ng_face_cyl + ng_face_cap;
  *n_edge_group = 1;

  if (i_rank == 0) {
    printf("ng_vtx = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM"\n", *ng_vtx, *ng_face);
  }

  /* Define distributions */
  PDM_g_num_t *distrib_vtx  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_edge_lim = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  distrib_vtx[0]      = 0;
  distrib_edge[0]     = 0;
  distrib_face[0]     = 0;
  distrib_edge_lim[0] = 0;

  PDM_g_num_t step_vtx       = *ng_vtx / n_rank;
  PDM_g_num_t remainder_vtx  = *ng_vtx % n_rank;

  PDM_g_num_t step_edge      = *ng_edge / n_rank;
  PDM_g_num_t remainder_edge = *ng_edge % n_rank;

  PDM_g_num_t step_face      = *ng_face / n_rank;
  PDM_g_num_t remainder_face = *ng_face % n_rank;

  PDM_g_num_t step_edge_lim      = ng_edge_lim / n_rank;
  PDM_g_num_t remainder_edge_lim = ng_edge_lim % n_rank;

  for (int i = 0; i < n_rank; i++) {
    distrib_vtx[i+1] = distrib_vtx[i] + step_vtx;
    if (i < remainder_vtx) {
      distrib_vtx[i+1]++;
    }

    distrib_edge[i+1] = distrib_edge[i] + step_edge;
    if (i < remainder_edge) {
      distrib_edge[i+1]++;
    }

    distrib_face[i+1] = distrib_face[i] + step_face;
    if (i < remainder_face) {
      distrib_face[i+1]++;
    }

    distrib_edge_lim[i+1] = distrib_edge_lim[i] + step_edge_lim;
    if (i < remainder_edge_lim) {
      distrib_edge_lim[i+1]++;
    }
  }
  *dn_vtx  = (int) distrib_vtx[i_rank+1]  - distrib_vtx[i_rank];
  *dn_edge = (int) distrib_edge[i_rank+1] - distrib_edge[i_rank];
  *dn_face = (int) distrib_face[i_rank+1] - distrib_face[i_rank];
  int dn_edge_lim = (int) distrib_edge_lim[i_rank+1] - distrib_edge_lim[i_rank];


  /*
   *  Vertices
   */
  *dvtx_coord = malloc (sizeof(double) * (*dn_vtx) * 3);
  double *_dvtx_coord = *dvtx_coord;

  double step_u = 2.*PDM_PI / (double) nu;
  double step_v = length_cyl / (double) (nv - 1);

  PDM_g_num_t b_vtx_v;
  PDM_g_num_t r_vtx_v;
  PDM_g_num_t ivtx = 0;

  /* Cylinder */
  if (distrib_vtx[i_rank] < ng_vtx_cyl) {
    b_vtx_v = distrib_vtx[i_rank] / nu;
    r_vtx_v = distrib_vtx[i_rank] % nu;

    for (PDM_g_num_t j = b_vtx_v; j < nv-1; j++) {

      PDM_g_num_t _b_vtx_u = 0;
      if (j == b_vtx_v) {
        _b_vtx_u = r_vtx_v;
      }

      double v = j * step_v;

      for (PDM_g_num_t i = _b_vtx_u; i < nu; i++) {
        double u = i * step_u;
        _dvtx_coord[3*ivtx    ] = radius * cos(u);
        _dvtx_coord[3*ivtx + 1] = radius * sin(u);
        _dvtx_coord[3*ivtx + 2] = v;
        ivtx++;
        if (ivtx == *dn_vtx) break;
      }
      if (ivtx == *dn_vtx) break;
    }
  }

  /* Cap */
  if (ivtx < *dn_vtx) {
    PDM_g_num_t gvtx = ng_vtx_cyl;

    double d = length_cap / 3.;
    double ta = tan(PDM_PI * 60./180.);
    double z2 = length_cap - sqrt(d*d / (1. + ta*ta));
    double r2 = (length_cap - z2) * ta;

    for (PDM_g_num_t j = nv_cap; j > 0; j--) {
      float t = (float) j / (float) nv_cap;

      double b3 = (1. - t) * (1. - t) * (1. - t);
      double b2 = 3. * (1. - t) * (1. - t) * t;
      double b1 = 3. * (1. - t) * t * t;
      double b0 = t * t * t;

      double r = (b0 + b1)*radius + b2*r2;
      double z = length_cyl + b1*d + b2*z2 + b3*length_cap;

      PDM_g_num_t __nu = 6*j;
      step_u = 2.*PDM_PI / (double) __nu;

      for (PDM_g_num_t i = 0; i < __nu; i++) {
        if (gvtx >= distrib_vtx[i_rank]) {
          double u = i * step_u;
          _dvtx_coord[3*ivtx    ] = r * cos(u);
          _dvtx_coord[3*ivtx + 1] = r * sin(u);
          _dvtx_coord[3*ivtx + 2] = z;
          ivtx++;
          if (ivtx == *dn_vtx) break;
        }
        gvtx++;

      }
      if (ivtx == *dn_vtx) break;
    }

    if (ivtx < *dn_vtx) {
      _dvtx_coord[3*ivtx    ] = 0.;
      _dvtx_coord[3*ivtx + 1] = 0.;
      _dvtx_coord[3*ivtx + 2] = length;
      ivtx++;
    }
  }
  free (distrib_vtx);


  /*
   *  Edges
   */
  *dedge_vtx  = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  *dedge_face = malloc (sizeof(PDM_g_num_t ) * (*dn_edge) * 2);
  PDM_g_num_t  *_dedge_vtx = *dedge_vtx;
  PDM_g_num_t  *_dedge_face = *dedge_face;

  PDM_g_num_t iedg = 0;
  PDM_g_num_t ifac;

  /* Cylinder - horizontal */
  if (distrib_edge[i_rank] < ng_edge_cyl_u) {
    const PDM_g_num_t b_edge_uv = distrib_edge[i_rank] / nu;
    const PDM_g_num_t r_edge_uv = distrib_edge[i_rank] % nu;

    for (PDM_g_num_t j = b_edge_uv; j < nv-1; j++) {

      PDM_g_num_t _b_edge_uu = 0;
      if (j == b_edge_uv) {
        _b_edge_uu = r_edge_uv;
      }

      for (PDM_g_num_t i = _b_edge_uu; i < nu; i++) {
        _dedge_vtx[2*iedg    ] = 1 + i        + nu*j;
        _dedge_vtx[2*iedg + 1] = 1 + (i+1)%nu + nu*j;

        _dedge_face[2*iedg] = 1 + i + nu*j;
        if (j == 0) {
          _dedge_face[2*iedg + 1] = 0;
        } else {
          _dedge_face[2*iedg + 1] = 1 + i + nu*(j-1);
        }
        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }

  /* Cylinder - vertical */
  if (iedg < *dn_edge && distrib_edge[i_rank] <= ng_edge_cyl_u + ng_edge_cyl_v) {
    const PDM_g_num_t b_edge_vv = (distrib_edge[i_rank] + iedg - ng_edge_cyl_u) / nu;
    const PDM_g_num_t r_edge_vv = (distrib_edge[i_rank] + iedg - ng_edge_cyl_u) % nu;

    for (PDM_g_num_t j = b_edge_vv; j < nv-1; j++) {

      PDM_g_num_t _b_edge_vu = 0;
      if (j == b_edge_vv) {
        _b_edge_vu = r_edge_vv;
      }

      for (PDM_g_num_t i = _b_edge_vu; i < nu; i++) {
        _dedge_vtx[2*iedg    ] = 1 + i + nu*j;
        _dedge_vtx[2*iedg + 1] = 1 + i + nu*(j+1);

        _dedge_face[2*iedg    ] = 1 + (i+nu-1)%nu + nu*j;
        _dedge_face[2*iedg + 1] = 1 + i           + nu*j;
        iedg++;
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;
    }
  }

  /* Cap */
  if (iedg < *dn_edge) {
    ivtx = ng_vtx_cyl + 1;
    ifac = ng_face_cyl + 1;
    PDM_g_num_t gedg = ng_edge_cyl_u + ng_edge_cyl_v;
    for (PDM_g_num_t k = nv_cap; k > 0; k--) {
      PDM_g_num_t n_vtx_k = 6*k;
      PDM_g_num_t n_vtx_km = PDM_MAX (1, 6*(k-1));

      PDM_g_num_t n_tri_k = 6*(2*k - 1);
      PDM_g_num_t n_tri_kp = 6*(2*k + 1);

      for (int j = 0; j < 6; j++) {
        for (PDM_g_num_t i = 0; i < k; i++) {
          PDM_g_num_t v0 = ivtx + k*j + i;
          PDM_g_num_t q = 2*i + (2*k-1)*j;

          v0 = k*j + i;
          if (gedg >= distrib_edge[i_rank]) {
            _dedge_vtx[2*iedg    ] = ivtx + v0;
            _dedge_vtx[2*iedg + 1] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

            _dedge_face[2*iedg    ] = ifac + (q + n_tri_k - 1)%n_tri_k;
            _dedge_face[2*iedg + 1] = ifac + q;
            iedg++;
            if (iedg == *dn_edge) break;
          }
          gedg++;

          if (gedg >= distrib_edge[i_rank]) {
            _dedge_vtx[2*iedg    ] = ivtx + v0;
            _dedge_vtx[2*iedg + 1] = ivtx + (v0+1)%n_vtx_k;

            _dedge_face[2*iedg] = ifac + q;
            if (k == nv_cap) {
              _dedge_face[2*iedg + 1] = ng_face_cyl - nu + 1 + i + k*j;
            } else {
              _dedge_face[2*iedg + 1] = ifac - n_tri_kp + 1 + 2*i + (2*k+1)*j;
            }
            iedg++;
            if (iedg == *dn_edge) break;
          }
          gedg++;

          if (i < k-1) {
            if (gedg >= distrib_edge[i_rank]) {
              _dedge_vtx[2*iedg    ] = ivtx + (v0+1)%n_vtx_k;
              _dedge_vtx[2*iedg + 1] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

              _dedge_face[2*iedg    ] = ifac + q;
              _dedge_face[2*iedg + 1] = ifac + q + 1;
              iedg++;
              if (iedg == *dn_edge) break;
            }
            gedg++;
          }
        }
        if (iedg == *dn_edge) break;
      }
      if (iedg == *dn_edge) break;

      ivtx += n_vtx_k;
      ifac += n_tri_k;
    }
  }
  free (distrib_edge);


  /* Edge groups */
  *dedge_group_idx = malloc (sizeof(int) * (*n_edge_group + 1));
  int *_dedge_group_idx = *dedge_group_idx;
  _dedge_group_idx[0] = 0;
  _dedge_group_idx[1] = 0;

  *dedge_group = malloc (sizeof(PDM_g_num_t) * dn_edge_lim);
  PDM_g_num_t *_dedge_group = *dedge_group;

  for (PDM_g_num_t i = distrib_edge_lim[i_rank]; i < distrib_edge_lim[i_rank+1]; i++) {
    _dedge_group[_dedge_group_idx[1]++] = 1 + i;
  }
  free (distrib_edge_lim);

  /*
   *  Faces
   */
  *dface_vtx_idx = malloc (sizeof(int) * (*dn_face + 1));
  int *_dface_vtx_idx = *dface_vtx_idx;
  _dface_vtx_idx[0] = 0;

  int n_quad;
  if (distrib_face[i_rank+1] <= ng_face_cyl) {
    n_quad = *dn_face;
  } else {
    n_quad = (int) PDM_MAX (0, ng_face_cyl - distrib_face[i_rank]);
  }

  for (int i = 0; i < n_quad; i++) {
    _dface_vtx_idx[i+1] = 4 + _dface_vtx_idx[i];
  }
  for (int i = n_quad; i < *dn_face; i++) {
    _dface_vtx_idx[i+1] = 3 + _dface_vtx_idx[i];
  }


  *dface_vtx  = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  *dface_edge = malloc (sizeof(PDM_g_num_t) * _dface_vtx_idx[*dn_face]);
  PDM_g_num_t *_dface_vtx = *dface_vtx;
  PDM_g_num_t *_dface_edge = *dface_edge;


  ifac = 0;
  /* Cylinder */
  if (n_quad > 0) {
    const PDM_g_num_t b_face_v = distrib_face[i_rank] / nu;
    const PDM_g_num_t r_face_v = distrib_face[i_rank] % nu;

    for (PDM_g_num_t j = b_face_v; j < nv-1; j++) {

      PDM_g_num_t _b_face_u = 0;
      if (j == b_face_v) {
        _b_face_u = r_face_v;
      }

      for (PDM_g_num_t i = _b_face_u; i < nu; i++) {

        PDM_g_num_t ip = (i+1)%nu;

        _dface_vtx[4*ifac    ] = 1 + i  + nu*j;
        _dface_vtx[4*ifac + 1] = 1 + ip + nu*j;
        _dface_vtx[4*ifac + 2] = 1 + ip + nu*(j+1);
        _dface_vtx[4*ifac + 3] = 1 + i  + nu*(j+1);

        _dface_edge[4*ifac    ] = 1 + i  + nu*j;
        _dface_edge[4*ifac + 1] = 1 + ng_edge_cyl_u + ip + nu*j;
        if (j == nv-2) {
          PDM_g_num_t _j = i / nv_cap;
          PDM_g_num_t _i = i % nv_cap;
          _dface_edge[4*ifac + 2] = -(ng_edge_cyl_u + ng_edge_cyl_v + 3*_i + (3*nv_cap-1)*_j + 2);
        } else {
          _dface_edge[4*ifac + 2] = -(1 + i  + nu*(j+1));
        }
        _dface_edge[4*ifac + 3] = -(1 + ng_edge_cyl_u + i + nu*j);

        ifac++;
        if (ifac == n_quad) break;
      }
      if (ifac == n_quad) break;
    }
  }

  /* Cap */
  if (ifac < *dn_face) {
    ivtx = ng_vtx_cyl + 1;
    iedg = ng_edge_cyl_u + ng_edge_cyl_v + 1;
    PDM_g_num_t gfac = ng_face_cyl;

    for (PDM_g_num_t k = nv_cap; k > 0; k--) {
      PDM_g_num_t n_vtx_k = 6*k;
      PDM_g_num_t n_vtx_km = PDM_MAX (1, n_vtx_k - 6);
      PDM_g_num_t n_edge_k = 18*k - 6;

      for (int j = 0; j < 6; j++) {
        for (PDM_g_num_t i = 0; i < k; i++) {
          if (gfac >= distrib_face[i_rank]) {
            _dface_vtx[n_quad + 3*ifac    ] = ivtx + k*j + i;
            _dface_vtx[n_quad + 3*ifac + 1] = ivtx + (k*j + i + 1)%n_vtx_k;
            _dface_vtx[n_quad + 3*ifac + 2] = ivtx + n_vtx_k + ((k-1)*j + i)%n_vtx_km;

            PDM_g_num_t e = 3*i + (3*k-1)*j;
            _dface_edge[n_quad + 3*ifac    ] = iedg + e + 1;
            _dface_edge[n_quad + 3*ifac + 1] = iedg + (e+2)%n_edge_k;
            _dface_edge[n_quad + 3*ifac + 2] = iedg + e;
            ifac++;
            if (ifac == *dn_face) break;
          }
          gfac++;

          if (i < k-1) {
            if (gfac >= distrib_face[i_rank]) {
              _dface_vtx[n_quad + 3*ifac    ] = ivtx + k*j + i + 1;
              _dface_vtx[n_quad + 3*ifac + 1] = ivtx + n_vtx_k + ((k-1)*j + i + 1)%n_vtx_km;
              _dface_vtx[n_quad + 3*ifac + 2] = ivtx + n_vtx_k + (k-1)*j + i;

              PDM_g_num_t e = 3*i + (3*k-1)*j;
              _dface_edge[n_quad + 3*ifac    ] = iedg + (e+3)%n_edge_k;
              _dface_edge[n_quad + 3*ifac + 1] = iedg + n_edge_k + 1 + 3*i + (3*k-4)*j;
              _dface_edge[n_quad + 3*ifac + 2] = iedg + (e+2)%n_edge_k;
              ifac++;
              if (ifac == *dn_face) break;
            }
            gfac++;
          }
        }
        if (ifac == *dn_face) break;
      }
      if (ifac == *dn_face) break;
      ivtx += n_vtx_k;
      iedg += n_edge_k;
    }
  }
  free (distrib_face);
}




static void
_gen_src_point_cloud
(
 const int           active_rank,
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   nv,
 const double        length,
 int                *npts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{

  if (active_rank) {
    /*
     *  Distributed mesh generation
     */
    PDM_g_num_t     ng_face;
    PDM_g_num_t     ng_vtx;
    PDM_g_num_t     ng_edge;
    int             dn_vtx;
    double         *dvtx_coord      = NULL;
    int             dn_face;
    int            *dface_vtx_idx   = NULL;
    PDM_g_num_t    *dface_vtx       = NULL;
    PDM_g_num_t    *dface_edge      = NULL;
    int             dn_edge;
    PDM_g_num_t    *dedge_vtx       = NULL;
    PDM_g_num_t    *dedge_face      = NULL;
    int             n_edge_group;
    int            *dedge_group_idx = NULL;
    PDM_g_num_t    *dedge_group     = NULL;

    double radius = 0.05 * length;

    _gen_rocket (comm,
                 nv,
                 length,
                 radius,
                 &ng_face,
                 &ng_vtx,
                 &ng_edge,
                 &dn_vtx,
                 &dvtx_coord,
                 &dn_face,
                 &dface_vtx_idx,
                 &dface_vtx,
                 &dface_edge,
                 &dn_edge,
                 &dedge_vtx,
                 &dedge_face,
                 &n_edge_group,
                 &dedge_group_idx,
                 &dedge_group);

    /*
     *  Create mesh partitions
     */
    int have_dface_part = 0;

    int *dface_part    = (int *) malloc (sizeof(int) * dn_face);
    int *dedge_vtx_idx = (int *) malloc (sizeof(int) * (dn_edge + 1));

    dedge_vtx_idx[0] = 0;
    for (int i = 0; i < dn_edge; i++) {
      dedge_vtx_idx[i+1] = 2 + dedge_vtx_idx[i];
    }

    int n_part = 1;
    /*
     *  Split mesh
     */
    int ppart_id;

#ifdef PDM_HAVE_PARMETIS
    PDM_part_split_t part_method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
    PDM_part_split_t part_method  = PDM_PART_SPLIT_PTSCOTCH;
#else
    PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;
#endif
#endif

    int n_property_face = 0;
    int *renum_properties_face = NULL;
    int n_property_edge = 0;
    int *renum_properties_edge = NULL;

    PDM_part_create (&ppart_id,
                     comm,
                     part_method,
                     "PDM_PART_RENUM_CELL_NONE",
                     "PDM_PART_RENUM_FACE_NONE",
                     n_property_face,
                     renum_properties_face,
                     n_property_edge,
                     renum_properties_edge,
                     n_part,
                     dn_face,
                     dn_edge,
                     dn_vtx,
                     n_edge_group,
                     NULL,
                     NULL,
                     NULL,
                     NULL,
                     have_dface_part,
                     dface_part,
                     dedge_face,
                     dedge_vtx_idx,
                     dedge_vtx,
                     NULL,
                     dvtx_coord,
                     NULL,
                     dedge_group_idx,
                     dedge_group);

    free (dface_part);

    free (dvtx_coord);
    free (dface_vtx_idx);
    free (dface_vtx);
    free (dface_edge);
    free (dedge_vtx_idx);
    free (dedge_vtx);
    free (dedge_face);
    free (dedge_group_idx);
    free (dedge_group);




    int _nFace;
    int _nEdge;
    int _nVtx;
    int _nEdgePartBound;
    int _nProc;
    int _nTPart;
    int _sFaceEdge;
    int _sEdgeVtx;
    int _sEdgeGroup;
    int _nEdgeGroup2;

    PDM_part_part_dim_get (ppart_id,
                           0,
                           &_nFace,
                           &_nEdge,
                           &_nEdgePartBound,
                           &_nVtx,
                           &_nProc,
                           &_nTPart,
                           &_sFaceEdge,
                           &_sEdgeVtx,
                           &_sEdgeGroup,
                           &_nEdgeGroup2);

    int         *_faceTag;
    int         *_faceEdgeIdx;
    int         *_faceEdge;
    PDM_g_num_t *_faceLNToGN;
    int         *_edgeTag;
    int         *_edgeFace;
    int         *_edgeVtxIdx;
    int         *_edgeVtx;
    PDM_g_num_t *_edgeLNToGN;
    int         *_edgePartBoundProcIdx;
    int         *_edgePartBoundPartIdx;
    int         *_edgePartBound;
    int         *_vtxTag;
    double      *_vtx;
    PDM_g_num_t *_vtxLNToGN;
    int         *_edgeGroupIdx;
    int         *_edgeGroup;
    PDM_g_num_t *_edgeGroupLNToGN;

    PDM_part_part_val_get (ppart_id,
                           0,
                           &_faceTag,
                           &_faceEdgeIdx,
                           &_faceEdge,
                           &_faceLNToGN,
                           &_edgeTag,
                           &_edgeFace,
                           &_edgeVtxIdx,
                           &_edgeVtx,
                           &_edgeLNToGN,
                           &_edgePartBoundProcIdx,
                           &_edgePartBoundPartIdx,
                           &_edgePartBound,
                           &_vtxTag,
                           &_vtx,
                           &_vtxLNToGN,
                           &_edgeGroupIdx,
                           &_edgeGroup,
                           &_edgeGroupLNToGN);

    *npts = _nVtx;
    *g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _nVtx);
    *coord = (double *)      malloc(sizeof(double)      * _nVtx * 3);

    memcpy (*g_num, _vtxLNToGN, sizeof(PDM_g_num_t) * _nVtx);
    memcpy (*coord, _vtx,       sizeof(double)      * _nVtx * 3);

    PDM_part_free (ppart_id);
  }

  else {
    *npts = 0;
    *g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (*npts));
    *coord = (double *)      malloc (sizeof(double)      * (*npts) * 3);
  }
}






static void
_closest_point_par
(
 PDM_MPI_Comm         comm,
 const int            n_max_per_leaf,
 const int            local_search_fun,
 const int            n_part_src,
 const int           *n_src,
 const double       **src_coord,
 const PDM_g_num_t  **src_g_num,
 const int            n_part_tgt,
 const int           *n_tgt,
 const double       **tgt_coord,
 const PDM_g_num_t  **tgt_g_num,
 PDM_g_num_t       ***closest_point_g_num,
 double            ***closest_point_dist2
 )
{
  /* Build parallel octree */
  const int depth_max = 31;

  int octree_id = PDM_para_octree_create (n_part_src,
                                          depth_max,
                                          n_max_per_leaf,
                                          0,//use_neighbours,
                                          comm);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_para_octree_point_cloud_set (octree_id,
                                     i_part,
                                     n_src[i_part],
                                     src_coord[i_part],
                                     src_g_num[i_part]);
  }

#if 0
  /* Compute global extents of source and target point clouds */
  double local_min[3] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double local_max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i_part = 0; i_part < n_part_src; i_part++) {
    const double *x = src_coord[i_part];
    for (int i = 0; i < n_src[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    const double *x = tgt_coord[i_part];
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  double global_extents[6];
  PDM_MPI_Allreduce(local_min, global_extents,     3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(local_max, global_extents + 3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_para_octree_build (octree_id, global_extents);
#else
  PDM_para_octree_build (octree_id, NULL);
#endif

  //PDM_para_octree_dump (octree_id);
  PDM_para_octree_dump_times (octree_id);




  /* Concatenate partitions */
  int _n_tgt = 0;

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    _n_tgt += n_tgt[i_part];
  }

  double      *_tgt_coord = malloc (sizeof(double)      * _n_tgt * 3);
  PDM_g_num_t *_tgt_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  PDM_g_num_t *_closest_src_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  double      *_closest_src_dist2 = malloc (sizeof(double)      * _n_tgt);

  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        _tgt_coord[_n_tgt + 3*i + j] = tgt_coord[i_part][3*i + j];
      }
      _tgt_g_num[_n_tgt + i] = tgt_g_num[i_part][i];
    }
    _n_tgt += n_tgt[i_part];
  }


  /* Search closest source points */
  PDM_para_octree_single_closest_point (octree_id,
                                        (const _local_search_fun_t) local_search_fun,
                                        _n_tgt,
                                        _tgt_coord,
                                        _tgt_g_num,
                                        _closest_src_g_num,
                                        _closest_src_dist2);

  /* Restore partitions */
  free (_tgt_coord);
  free (_tgt_g_num);

  *closest_point_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *closest_point_dist2 = (double **) malloc (sizeof(double *) * n_part_tgt);
  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    (*closest_point_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_tgt[i_part]);
    (*closest_point_dist2)[i_part] = (double *) malloc (sizeof(double) * n_tgt[i_part]);

    for (int i = 0; i < n_tgt[i_part]; i++) {
      (*closest_point_g_num)[i_part][i] = _closest_src_g_num[_n_tgt + i];
      (*closest_point_dist2)[i_part][i] = _closest_src_dist2[_n_tgt + i];
    }
    _n_tgt += n_tgt[i_part];
  }
  free (_closest_src_g_num);
  free (_closest_src_dist2);

  /* Free parallel octree */
  PDM_para_octree_free (octree_id);
}



static void
_closest_point_seq
(
 PDM_MPI_Comm         comm,
 const int            n_part_src,
 const int           *n_src,
 const double       **src_coord,
 const PDM_g_num_t  **src_g_num,
 const int            n_part_tgt,
 const int           *n_tgt,
 const double       **tgt_coord,
 const PDM_g_num_t  **tgt_g_num,
 PDM_g_num_t       ***closest_point_g_num,
 double            ***closest_point_dist2
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  const double tolerance = 1e-4;
  const int depth_max = 31;
  const int points_in_leaf_max = 4;

  int octree_id = PDM_octree_create (n_part_src,
                                     depth_max,
                                     points_in_leaf_max,
                                     tolerance,
                                     comm);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_octree_point_cloud_set (octree_id,
                                i_part,
                                n_src[i_part],
                                src_coord[i_part],
                                src_g_num[i_part]);
  }

  /* Build octree */
  PDM_timer_t *timer = PDM_timer_create ();
  double t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  PDM_octree_build (octree_id);

  PDM_timer_hang_on (timer);
  double t_end = PDM_timer_elapsed (timer) - t_begin;
  t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  double t_max;
  PDM_MPI_Reduce (&t_end, &t_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, 0, PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf ("Octree build elapsed time = %.3gs\n", t_max);
  }

  /* Concatenate partitions */
  int _n_tgt = 0;

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    _n_tgt += n_tgt[i_part];
  }

  double      *_tgt_coord = malloc (sizeof(double)      * _n_tgt * 3);
  PDM_g_num_t *_tgt_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  PDM_g_num_t *_closest_src_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  double      *_closest_src_dist2 = malloc (sizeof(double)      * _n_tgt);

  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        _tgt_coord[_n_tgt + 3*i + j] = tgt_coord[i_part][3*i + j];
      }
      _tgt_g_num[_n_tgt + i] = tgt_g_num[i_part][i];
    }
    _n_tgt += n_tgt[i_part];
  }

  /* Search closest source points */
  PDM_timer_hang_on (timer);
  t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  PDM_octree_closest_point (octree_id,
                            _n_tgt,
                            _tgt_coord,
                            _tgt_g_num,
                            _closest_src_g_num,
                            _closest_src_dist2);

  PDM_timer_hang_on (timer);
  t_end = PDM_timer_elapsed (timer) - t_begin;
  PDM_timer_resume (timer);

  PDM_MPI_Reduce (&t_end, &t_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, 0, PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf ("Octree closest point elapsed time = %.3gs\n", t_max);
  }
  PDM_timer_free (timer);

  /* Restore partitions */
  free (_tgt_coord);
  free (_tgt_g_num);

  *closest_point_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *closest_point_dist2 = (double **) malloc (sizeof(double *) * n_part_tgt);
  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    (*closest_point_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_tgt[i_part]);
    (*closest_point_dist2)[i_part] = (double *) malloc (sizeof(double) * n_tgt[i_part]);

    for (int i = 0; i < n_tgt[i_part]; i++) {
      (*closest_point_g_num)[i_part][i] = _closest_src_g_num[_n_tgt + i];
      (*closest_point_dist2)[i_part][i] = _closest_src_dist2[_n_tgt + i];
    }
    _n_tgt += n_tgt[i_part];
  }
  free (_closest_src_g_num);
  free (_closest_src_dist2);

  /* Free octree */
  PDM_octree_free (octree_id);
}



static void
_write_point_cloud
(
 const char        *filename,
 const char        *header,
 const int          n_pts,
 const double       coord[],
 const PDM_g_num_t  g_num[],
 const double       field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (field != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS field double\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, "%f\n", field[i]);
    }
  }

  fclose(f);
}



/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init (&argc, &argv);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  srand (time(NULL) + i_rank);


  /*
   *  Set default values
   */
  PDM_g_num_t nv              = 10;
  double      length          = 1.;
  double      domain_size     = 3.;
  PDM_g_num_t n_tgt           = 10;
  int         n_max_per_leaf  = 10;
  int         n_methods       = 5;
  int         on_ground       = 0;
  int         n_proc_data_src = -1;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nv,
              &length,
              &domain_size,
              &n_tgt,
              &n_max_per_leaf,
              &n_methods,
              &on_ground,
              &n_proc_data_src);


  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  double      *tgt_coord = NULL;
  PDM_g_num_t *tgt_g_num = NULL;
  int _n_tgt;

  double xmin = -0.5*domain_size;
  double origin[3] = {xmin, xmin, xmin};
  _gen_cloud_random (PDM_MPI_COMM_WORLD,
                     n_tgt,
                     origin,
                     domain_size,
                     &_n_tgt,
                     &tgt_g_num,
                     &tgt_coord);


  /*
   *  Define the source point cloud
   */
  PDM_MPI_Comm src_comm = PDM_MPI_COMM_WORLD;
  int active_rank_src = 1;
  if (n_proc_data_src > 0 && n_proc_data_src < n_rank) {
    int rank_in_node = PDM_io_mpi_node_rank (PDM_MPI_COMM_WORLD);

    int n_nodes = 0;
    int i_node = -1;
    int master_rank = 0;
    if (rank_in_node == 0) {
      master_rank = 1;
    }

    int *rank_in_nodes = malloc(sizeof(int) * n_rank);

    PDM_MPI_Allreduce (&master_rank, &n_nodes, 1, PDM_MPI_INT, PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
    PDM_MPI_Allgather (&rank_in_node, 1, PDM_MPI_INT, rank_in_nodes, 1, PDM_MPI_INT, PDM_MPI_COMM_WORLD);

    active_rank_src = 0;

    for (int i = 0; i < i_rank; i++) {
      if (rank_in_nodes[i] == 0) {
        i_node += 1;
      }
    }

    if (n_proc_data_src <= n_nodes) {
      if (i_node < n_proc_data_src && rank_in_node == 0) {
        active_rank_src = 1;
      }
    }

    else {

      if (rank_in_node < (n_proc_data_src / n_nodes)) {
        active_rank_src = 1;
      }
      if ((rank_in_node == (n_proc_data_src / n_nodes)) &&
          (i_node < (n_proc_data_src % n_nodes))) {
        active_rank_src = 1;
      }

    }

    PDM_MPI_Comm_split (PDM_MPI_COMM_WORLD, active_rank_src, i_rank, &src_comm);

    free (rank_in_nodes);
  }


  int n_part_src = 1;
  int _n_src;
  double *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;

  _gen_src_point_cloud (active_rank_src,
                        src_comm,
                        nv,
                        length,
                        &_n_src,
                        &src_g_num,
                        &src_coord);


  PDM_g_num_t n_src_local, n_src_global;
  n_src_local = _n_src;
  PDM_MPI_Allreduce (&n_src_local, &n_src_global, 1, PDM__PDM_MPI_G_NUM,
                     PDM_MPI_SUM, PDM_MPI_COMM_WORLD);
  PDM_g_num_t src_g_num_min = 1;
  PDM_g_num_t src_g_num_max = n_src_global;

  double z_shift;
  if (on_ground) {
    z_shift = 0.5*domain_size;
  } else {
    z_shift = 0.5*length;
  }

  for (int i = 0; i < _n_src; i++) {
    src_coord[3*i+2] -= z_shift;
  }


  if (0) {
    char filename[999];

    sprintf(filename, "tgt_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "tgt",
                        _n_tgt,
                        tgt_coord,
                        tgt_g_num,
                        NULL);

    sprintf(filename, "src_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "src",
                        _n_src,
                        src_coord,
                        src_g_num,
                        NULL);
  }


  /*
   *  Compare methods
   */
  double *elapsed = malloc (sizeof(double) * n_methods);
  double t_begin, t_end;
  PDM_timer_t *timer = PDM_timer_create ();

  PDM_g_num_t **closest_point_g_num = NULL;
  double      **closest_point_dist2 = NULL;

  PDM_g_num_t *true_closest_point_g_num = NULL;
  double      *true_closest_point_dist2 = NULL;

  double elapsed_min = HUGE_VAL;
  for (int method = 0; method < n_methods; method++) {

    PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

    if (i_rank == 0) {
      printf("\n\nMethod %d\n", method);
    }

    if (method > 0) PDM_timer_hang_on (timer);
    t_begin = PDM_timer_elapsed (timer);
    PDM_timer_resume (timer);

    /* Compute closest points */
    if (method == 0) {
      _closest_point_seq (PDM_MPI_COMM_WORLD,
                          n_part_src,
                          (const int *) &_n_src,
                          (const double **) &src_coord,
                          (const PDM_g_num_t **) &src_g_num,
                          n_part_tgt,
                          (const int *) &_n_tgt,
                          (const double **) &tgt_coord,
                          (const PDM_g_num_t **) &tgt_g_num,
                          &closest_point_g_num,
                          &closest_point_dist2);
    }

    else if (method > 0) {
      _closest_point_par (PDM_MPI_COMM_WORLD,
                          n_max_per_leaf,
                          method - 1,
                          n_part_src,
                          (const int *) &_n_src,
                          (const double **) &src_coord,
                          (const PDM_g_num_t **) &src_g_num,
                          n_part_tgt,
                          (const int *) &_n_tgt,
                          (const double **) &tgt_coord,
                          (const PDM_g_num_t **) &tgt_g_num,
                          &closest_point_g_num,
                          &closest_point_dist2);
    }

    PDM_timer_hang_on (timer);
    t_end = PDM_timer_elapsed (timer) - t_begin;
    PDM_timer_resume (timer);

    PDM_MPI_Reduce (&t_end,
                    elapsed + method, 1,
                    PDM_MPI_DOUBLE,
                    PDM_MPI_MAX,
                    0,
                    PDM_MPI_COMM_WORLD);

    if (i_rank == 0) {
      printf ("\nTotal elapsed time = %.3gs\n", elapsed[method]);
      if (elapsed[method] < elapsed_min) {
        elapsed_min = elapsed[method];
      }
    }

    if (0) {
      // Check g num validity
      PDM_g_num_t gmin = 99999999;
      PDM_g_num_t gmax = -gmin;
      double dmin = HUGE_VAL;
      double dmax = -HUGE_VAL;
      for (int i = 0; i < _n_tgt; i++) {
        gmin = PDM_MIN (gmin, closest_point_g_num[0][i]);
        gmax = PDM_MAX (gmax, closest_point_g_num[0][i]);
        dmin = PDM_MIN (dmin, closest_point_dist2[0][i]);
        dmax = PDM_MAX (dmax, closest_point_dist2[0][i]);

        if (closest_point_g_num[0][i] < src_g_num_min) {
          printf("[%d] !!! pt ("PDM_FMT_G_NUM"): closest_point_g_num = "PDM_FMT_G_NUM" < min\n",
                 i_rank, tgt_g_num[i], closest_point_g_num[0][i]);
        }
        else if (closest_point_g_num[0][i] > src_g_num_max) {
          printf("[%d] !!! pt ("PDM_FMT_G_NUM"): closest_point_g_num = "PDM_FMT_G_NUM" > max\n",
                 i_rank, tgt_g_num[i], closest_point_g_num[0][i]);
        }
      }

      printf("[%d] gmin/gmax = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n", i_rank, gmin, gmax);
      printf("[%d] dmin/dmax = %f / %f\n", i_rank, dmin, dmax);
    }


    if (0) {
      char filename[999];
      double *_dist = malloc (sizeof(double) * _n_tgt);
      for (int i = 0; i < _n_tgt; i++) {
          _dist[i] = sqrt (closest_point_dist2[0][i]);
      }

      sprintf(filename, "check_dist_%d_%3.3d.vtk", method, i_rank);
      _write_point_cloud (filename,
                          "tgt",
                          _n_tgt,
                          tgt_coord,
                          NULL,
                          _dist);

      sprintf(filename, "check_gnum_%d_%3.3d.vtk", method, i_rank);
      _write_point_cloud (filename,
                          "tgt",
                          _n_tgt,
                          tgt_coord,
                          closest_point_g_num[0],
                          NULL);

      free (_dist);
    }

    if (method == 0) {
      true_closest_point_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
      true_closest_point_dist2 = malloc (sizeof(double)      * _n_tgt);
      memcpy (true_closest_point_g_num, closest_point_g_num[0], sizeof(PDM_g_num_t) * _n_tgt);
      memcpy (true_closest_point_dist2, closest_point_dist2[0], sizeof(double)      * _n_tgt);
    } else {
      if (true_closest_point_g_num != NULL) {
        for (int i = 0; i < _n_tgt; i++) {
          if (closest_point_g_num[0][i] != true_closest_point_g_num[i]) {
            double d0 = sqrt(closest_point_dist2[0][i]);
            double d1 = sqrt(true_closest_point_dist2[i]);
            printf("[%d] error point ("PDM_FMT_G_NUM"): ("PDM_FMT_G_NUM") / ("PDM_FMT_G_NUM"), %f / %f (rel. err. %e)\n",
                   i_rank, tgt_g_num[i], closest_point_g_num[0][i], true_closest_point_g_num[i], d0, d1, d0/d1 - 1.);
          }
        }
      }
    }



    for (int i = 0; i < n_part_tgt; i++) {
      free (closest_point_g_num[i]);
      free (closest_point_dist2[i]);
    }
    free (closest_point_g_num);
    free (closest_point_dist2);
  }




  PDM_timer_free (timer);


  /*
   *  Summary
   */
  if (i_rank == 0) {
    printf("\n\n\n");
    for (int method = 0; method < n_methods; method++) {
      printf ("method %d: elapsed = %.4fs, relative to min = %.3f\n",
              method, elapsed[method], elapsed[method] / elapsed_min);
    }
  }


  /*
   *  Finalize
   */
  free (elapsed);
  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
