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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_timer.h"
#include "pdm_part.h"
#include "pdm_geom_elem.h"
#include "pdm_dcube_gen.h"
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
           int           *grid,
           int           *n_max_per_leaf,
           int           *on_ground,
           int           *n_proc_data_src,
           int           *post)
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
    else if (strcmp(argv[i], "-grid") == 0) {
      *grid = 1;
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
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
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
  PDM_g_num_t *distrib;
  PDM_malloc(distrib, n_rank + 1, PDM_g_num_t);
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

  PDM_malloc(*g_num, dn_pts    , PDM_g_num_t);
  PDM_malloc(*coord, dn_pts * 3, double     );
  for (int i = 0; i < *_n_pts; i++) {
    (*g_num)[i] = 1 + i + distrib[i_rank];
    for (int j = 0; j < 3; j++) {
      (*coord)[3*i+j] = origin[j] + length * _rand01();
    }
  }

  PDM_free(distrib);
}



static void
_gen_cloud_grid
(
 PDM_MPI_Comm   comm,
 int            n_part,
 PDM_g_num_t    n_vtx_seg,
 double         origin[3],
 double         length,
 int          **n_pts,
 PDM_g_num_t ***pts_g_num,
 double      ***pts_coord
 )
{
  int n_rank, i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_dcube_t *dcube = PDM_dcube_gen_init (comm,
                                           n_vtx_seg,
                                           length,
                                           origin[0],
                                           origin[1],
                                           origin[2],
                                           PDM_OWNERSHIP_KEEP);

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  int          dface_vtx_l;
  int          dface_group_l;

  PDM_dcube_gen_dim_get (dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);


  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  PDM_dcube_gen_data_get (dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  /*
   *  Create mesh partitions
   */

  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

  // int ppart_id = 0;
  int have_dcell_part = 0;

  int *dcell_part = NULL;
  PDM_malloc(dcell_part, dn_cell, int);

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;
  PDM_part_t *ppart = PDM_part_create (PDM_MPI_COMM_WORLD,
                                       part_method,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       "PDM_PART_RENUM_FACE_NONE",
                                       n_property_cell,
                                       renum_properties_cell,
                                       n_property_face,
                                       renum_properties_face,
                                       n_part,
                                       dn_cell,
                                       dn_face,
                                       dn_vtx,
                                       n_face_group,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       have_dcell_part,
                                       dcell_part,
                                       dface_cell,
                                       dface_vtx_idx,
                                       dface_vtx,
                                       NULL,
                                       dvtx_coord,
                                       NULL,
                                       dface_group_idx,
                                       dface_group);

  PDM_free(dcell_part);


  PDM_malloc(*n_pts    , n_part, int          );
  PDM_malloc(*pts_g_num, n_part, PDM_g_num_t *);
  PDM_malloc(*pts_coord, n_part, double      *);

  for (int i_part = 0; i_part < n_part; i_part++) {
    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get (ppart,
                           i_part,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_t_part,
                           &s_cell_face,
                           &s_face_vtx,
                           &s_face_group,
                           &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart,
                           i_part,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    (*n_pts)[i_part] = n_cell;
    PDM_malloc((*pts_g_num)[i_part], n_cell, PDM_g_num_t);
    for (int i = 0; i < n_cell; i++) {
      (*pts_g_num)[i_part][i] = cell_ln_to_gn[i];
    }


    const int is_oriented = 0;
    PDM_malloc((*pts_coord)[i_part], n_cell * 3, double);
    double *cell_volume = NULL;
    PDM_malloc(cell_volume, n_cell, double);
    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume,
                                        (*pts_coord)[i_part],
                                        NULL,
                                        NULL);
    PDM_free(cell_volume);
  }

  PDM_part_free (ppart);
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
  PDM_g_num_t *distrib_vtx      = NULL;
  PDM_g_num_t *distrib_edge     = NULL;
  PDM_g_num_t *distrib_face     = NULL;
  PDM_g_num_t *distrib_edge_lim = NULL;
  PDM_malloc(distrib_vtx     , n_rank + 1, PDM_g_num_t);
  PDM_malloc(distrib_edge    , n_rank + 1, PDM_g_num_t);
  PDM_malloc(distrib_face    , n_rank + 1, PDM_g_num_t);
  PDM_malloc(distrib_edge_lim, n_rank + 1, PDM_g_num_t);
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
  PDM_malloc(*dvtx_coord, (*dn_vtx) * 3, double);
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
  PDM_free(distrib_vtx);


  /*
   *  Edges
   */
  PDM_malloc(*dedge_vtx , (*dn_edge) * 2, PDM_g_num_t);
  PDM_malloc(*dedge_face, (*dn_edge) * 2, PDM_g_num_t);
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
  PDM_free(distrib_edge);


  /* Edge groups */
  PDM_malloc(*dedge_group_idx, *n_edge_group + 1, int);
  int *_dedge_group_idx = *dedge_group_idx;
  _dedge_group_idx[0] = 0;
  _dedge_group_idx[1] = 0;

  PDM_malloc(*dedge_group, dn_edge_lim, PDM_g_num_t);
  PDM_g_num_t *_dedge_group = *dedge_group;

  for (PDM_g_num_t i = distrib_edge_lim[i_rank]; i < distrib_edge_lim[i_rank+1]; i++) {
    _dedge_group[_dedge_group_idx[1]++] = 1 + i;
  }
  PDM_free(distrib_edge_lim);

  /*
   *  Faces
   */
  PDM_malloc(*dface_vtx_idx, *dn_face + 1, int);
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


  PDM_malloc(*dface_vtx , _dface_vtx_idx[*dn_face], PDM_g_num_t);
  PDM_malloc(*dface_edge, _dface_vtx_idx[*dn_face], PDM_g_num_t);
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
  PDM_free(distrib_face);
}


static void
_get_connectivity
(
 PDM_part_t    *ppart,
 int            n_part,
 int          **n_face,
 int         ***face_edge_idx,
 int         ***face_edge,
 int         ***face_vtx_idx,
 int         ***face_vtx,
 PDM_g_num_t ***face_ln_to_gn,
 int          **n_edge,
 int         ***edge_vtx_idx,
 int         ***edge_vtx,
 int          **n_vtx,
 double      ***vtx_coord,
 PDM_g_num_t ***vtx_ln_to_gn
 )
{
  PDM_malloc(*n_face       , n_part, int          );
  PDM_malloc(*face_edge_idx, n_part, int         *);
  PDM_malloc(*face_edge    , n_part, int         *);
  PDM_malloc(*face_vtx_idx , n_part, int         *);
  PDM_malloc(*face_vtx     , n_part, int         *);
  PDM_malloc(*face_ln_to_gn, n_part, PDM_g_num_t *);

  PDM_malloc(*n_edge,       n_part, int  );
  PDM_malloc(*edge_vtx_idx, n_part, int *);
  PDM_malloc(*edge_vtx,     n_part, int *);

  PDM_malloc(*n_vtx,         n_part, int          );
  PDM_malloc(*vtx_coord,     n_part, double      *);
  PDM_malloc(*vtx_ln_to_gn,  n_part, PDM_g_num_t *);


  for (int ipart = 0; ipart < n_part; ipart++) {

    int _n_face;
    int _n_edge;
    int _n_edge_part_bound;
    int _n_vtx;
    int _n_proc;
    int _n_t_part;
    int _sFace_edge;
    int _s_edge_vtx;
    int _s_edge_group;
    int _n_edge_group2;

    PDM_part_part_dim_get (ppart,
                           ipart,
                           &_n_face,
                           &_n_edge,
                           &_n_edge_part_bound,
                           &_n_vtx,
                           &_n_proc,
                           &_n_t_part,
                           &_sFace_edge,
                           &_s_edge_vtx,
                           &_s_edge_group,
                           &_n_edge_group2);

    int         *_faceTag;
    int         *_face_edge_idx;
    int         *_face_edge;
    PDM_g_num_t *_face_ln_to_gn;
    int         *_edge_tag;
    int         *_edge_face;
    int         *_edge_vtx_idx;
    int         *_edge_vtx;
    PDM_g_num_t *_edge_ln_to_gn;
    int         *_edge_part_bound_proc_idx;
    int         *_edge_part_bound_part_idx;
    int         *_edge_part_bound;
    int         *_vtx_tag;
    double      *_vtx;
    PDM_g_num_t *_vtx_ln_to_gn;
    int         *_edge_group_idx;
    int         *_edge_group;
    PDM_g_num_t *_edge_groupLNToGN;

    PDM_part_part_val_get (ppart,
                           ipart,
                           &_faceTag,
                           &_face_edge_idx,
                           &_face_edge,
                           &_face_ln_to_gn,
                           &_edge_tag,
                           &_edge_face,
                           &_edge_vtx_idx,
                           &_edge_vtx,
                           &_edge_ln_to_gn,
                           &_edge_part_bound_proc_idx,
                           &_edge_part_bound_part_idx,
                           &_edge_part_bound,
                           &_vtx_tag,
                           &_vtx,
                           &_vtx_ln_to_gn,
                           &_edge_group_idx,
                           &_edge_group,
                           &_edge_groupLNToGN);

    /*for (int i = 0; i < _n_face; i++) {
      printf("face ("PDM_FMT_G_NUM"), edges =", _face_ln_to_gn[i]);
      for (int j = _face_edge_idx[i]; j < _face_edge_idx[i+1]; j++) {
      printf(" ("PDM_FMT_G_NUM")", _edge_ln_to_gn[PDM_ABS(_face_edge[j])-1]);
      }
      printf("\n");
      }*/


    /* Faces */
    (*n_face)[ipart] = _n_face;
    PDM_malloc((*face_edge_idx)[ipart], _n_face + 1, int        );
    PDM_malloc((*face_edge    )[ipart], _sFace_edge, int        );
    PDM_malloc((*face_vtx_idx )[ipart], _n_face + 1, int        );
    PDM_malloc((*face_vtx     )[ipart], _sFace_edge, int        );
    PDM_malloc((*face_ln_to_gn)[ipart], _n_face    , PDM_g_num_t);

    memcpy ((*face_edge_idx)[ipart], _face_edge_idx, (_n_face + 1) * sizeof(int        ));
    memcpy ((*face_edge    )[ipart], _face_edge    , _sFace_edge   * sizeof(int        ));
    memcpy ((*face_vtx_idx )[ipart], _face_edge_idx, (_n_face + 1) * sizeof(int        ));
    memcpy ((*face_ln_to_gn)[ipart], _face_ln_to_gn, _n_face       * sizeof(PDM_g_num_t));

    /* Edges */
    (*n_edge)[ipart] = _n_edge;
    PDM_malloc((*edge_vtx_idx)[ipart], (_n_edge + 1), int);
    PDM_malloc((*edge_vtx    )[ipart], _s_edge_vtx  , int);

    memcpy ((*edge_vtx_idx)[ipart], _edge_vtx_idx, (_n_edge + 1) * sizeof(int));
    memcpy ((*edge_vtx    )[ipart], _edge_vtx    , _s_edge_vtx   * sizeof(int));

    /* Vertices */
    (*n_vtx)[ipart] = _n_vtx;
    PDM_malloc((*vtx_coord   )[ipart], (3 * _n_vtx), double     );
    PDM_malloc((*vtx_ln_to_gn)[ipart], _n_vtx      , PDM_g_num_t);

    memcpy ((*vtx_coord   )[ipart], _vtx         , 3 *_n_vtx * sizeof(double     ));
    memcpy ((*vtx_ln_to_gn)[ipart], _vtx_ln_to_gn,    _n_vtx * sizeof(PDM_g_num_t));

    /* Compute face-vtx connectivity */
    int *_face_vtx = (*face_vtx)[ipart];

    int *vtx_edge_idx;
    PDM_malloc(vtx_edge_idx, _n_vtx + 1, int);

    for (int i = 0; i < _n_vtx + 1; i++) {
      vtx_edge_idx[i] = 0;
    }

    for (int i = 0; i < _n_edge; i++) {
      int ivtx1 = _edge_vtx[2*i];
      int ivtx2 = _edge_vtx[2*i + 1];

      vtx_edge_idx[ivtx1] += 1;
      vtx_edge_idx[ivtx2] += 1;
    }

    for (int i = 1; i < _n_vtx + 1; i++) {
      vtx_edge_idx[i] = vtx_edge_idx[i] + vtx_edge_idx[i-1];
    }

    int *vtx_edge   = NULL;
    int *vtx_edge_n = NULL;
    PDM_malloc(vtx_edge  , vtx_edge_idx[_n_vtx], int);
    PDM_malloc(vtx_edge_n, _n_vtx              , int);
    for (int i = 0; i < _n_vtx; i++) {
      vtx_edge_n[i] = 0;
    }

    for (int i = 0; i < _n_edge; i++) {
      int ivtx1 = _edge_vtx[2*i] - 1;
      int ivtx2 = _edge_vtx[2*i + 1] - 1;
      int iedge = i + 1;

      vtx_edge[vtx_edge_idx[ivtx1] + vtx_edge_n[ivtx1]] = iedge;
      vtx_edge[vtx_edge_idx[ivtx2] + vtx_edge_n[ivtx2]] = iedge;
      vtx_edge_n[ivtx1] += 1;
      vtx_edge_n[ivtx2] += 1;
    }
    PDM_free(vtx_edge_n);

    for (int i = 0; i < _n_face; i++) {
      int idx = _face_edge_idx[i];
      int __n_edge = _face_edge_idx[i+1] - idx;
      int *_edges = _face_edge + idx;
      int *_vertices = _face_vtx + idx;

      int edge_cur = _edges[0];
      int vtx_deb =  _edge_vtx[2*(edge_cur - 1)];
      _vertices[0] = vtx_deb;
      int vtx_cur =  _edge_vtx[2*(edge_cur - 1) + 1];
      int idx_vtx = 0;

      while (vtx_deb != vtx_cur) {
        _vertices[++idx_vtx] = vtx_cur;
        int find_vtx = 0;

        for (int j = vtx_edge_idx[vtx_cur - 1]; j <  vtx_edge_idx[vtx_cur]; j++) {
          for (int k = 0; k < __n_edge; k++) {
            if ((_edges[k] == vtx_edge[j]) && (_edges[k] != edge_cur)) {
              edge_cur = _edges[k];
              if (_edge_vtx[2*(_edges[k]-1)] == vtx_cur) {
                vtx_cur = _edge_vtx[2*(_edges[k]-1) + 1];
              }
              else {
                vtx_cur = _edge_vtx[2*(_edges[k]-1)];
              }
              find_vtx = 1;
              break;
            }
          }
          if (find_vtx)
            break;
        }
        if (!find_vtx) {
          printf("error face ("PDM_FMT_G_NUM"), vtx tmp:\n", _face_ln_to_gn[i]);
          for (int l = 0; l < idx_vtx; l++) {
            printf("  %d ("PDM_FMT_G_NUM")\n", _vertices[l], _vtx_ln_to_gn[_vertices[l]-1]);
          }
          printf("\n");
          PDM_error(__FILE__, __LINE__, 0,"Error to compute vtx_edge !!!!\n");
          abort();
        }
      }
      /*printf("face ("PDM_FMT_G_NUM"), vtx :", _face_ln_to_gn[i]);
        for (int l = 0; l < __n_edge; l++) {
        printf(" ("PDM_FMT_G_NUM")", _vtx_ln_to_gn[_vertices[l]-1]);
        }
        printf("\n");*/
    }

    PDM_free(vtx_edge);
    PDM_free(vtx_edge_idx);

  }
}




static void
_gen_src_mesh
(
 const int            active_rank,
 PDM_MPI_Comm         comm,
 int                  n_part,
 const PDM_g_num_t    nv,
 const double         length,
 int                **n_vtx,
 PDM_g_num_t       ***vtx_g_num,
 double            ***vtx_coord,
 int                **n_face,
 PDM_g_num_t       ***face_g_num,
 int               ***face_vtx_idx,
 int               ***face_vtx
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

    int *dface_part    = NULL;
    int *dedge_vtx_idx = NULL;
    PDM_malloc(dface_part   , dn_face    , int);
    PDM_malloc(dedge_vtx_idx, dn_edge + 1, int);

    dedge_vtx_idx[0] = 0;
    for (int i = 0; i < dn_edge; i++) {
      dedge_vtx_idx[i+1] = 2 + dedge_vtx_idx[i];
    }


    /*
     *  Split mesh
     */
    PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;

    int n_property_face = 0;
    int *renum_properties_face = NULL;
    int n_property_edge = 0;
    int *renum_properties_edge = NULL;

    PDM_part_t *ppart = PDM_part_create (comm,
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

    PDM_free(dface_part);
    PDM_free(dvtx_coord);
    PDM_free(dface_vtx_idx);
    PDM_free(dface_vtx);
    PDM_free(dface_edge);
    PDM_free(dedge_vtx_idx);
    PDM_free(dedge_vtx);
    PDM_free(dedge_face);
    PDM_free(dedge_group_idx);
    PDM_free(dedge_group);


    int **face_edge_idx = NULL;
    int **face_edge     = NULL;
    int  *n_edge        = NULL;
    int **edge_vtx_idx  = NULL;
    int **edge_vtx      = NULL;
    _get_connectivity (ppart,
                       n_part,
                       n_face,
                       &face_edge_idx,
                       &face_edge,
                       face_vtx_idx,
                       face_vtx,
                       face_g_num,
                       &n_edge,
                       &edge_vtx_idx,
                       &edge_vtx,
                       n_vtx,
                       vtx_coord,
                       vtx_g_num);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(face_edge_idx[i_part]);
      PDM_free(face_edge[i_part]);
      PDM_free(edge_vtx_idx[i_part]);
      PDM_free(edge_vtx[i_part]);
    }
    PDM_free(n_edge);
    PDM_free(face_edge_idx);
    PDM_free(face_edge);
    PDM_free(edge_vtx_idx);
    PDM_free(edge_vtx);

    PDM_part_free (ppart);
  }

  else {
    PDM_malloc(*n_vtx       , n_part, int          );
    PDM_malloc(*vtx_g_num   , n_part, PDM_g_num_t *);
    PDM_malloc(*vtx_coord   , n_part, double      *);
    PDM_malloc(*n_face      , n_part, int          );
    PDM_malloc(*face_vtx_idx, n_part, int         *);
    PDM_malloc(*face_vtx    , n_part, int         *);
    PDM_malloc(*face_g_num  , n_part, PDM_g_num_t *);

    for (int i_part = 0; i_part < n_part; i_part++) {
      (*n_vtx)[i_part] = 0;
      PDM_malloc((*vtx_g_num)[i_part], (*n_vtx)[i_part]    , PDM_g_num_t);
      PDM_malloc((*vtx_coord)[i_part], (*n_vtx)[i_part] * 3, double     );

      (*n_face)[i_part] = 0;
      PDM_malloc((*face_g_num  )[i_part], (*n_face)[i_part]     , PDM_g_num_t);
      PDM_malloc((*face_vtx_idx)[i_part],((*n_face)[i_part] + 1), int        );
      (*face_vtx_idx)[i_part][0] = 0;
      PDM_malloc((*face_vtx)[i_part], (*face_vtx_idx)[i_part][(*n_face)[i_part]], int);
    }
  }
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



static void
_write_polydata
(
 const char        *filename,
 const char        *header,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[]
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
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
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
  PDM_g_num_t ng_tgt          = 10;
  int         grid            = 0;
  int         n_max_per_leaf  = 10;
  int         on_ground       = 0;
  int         n_proc_data_src = -1;
  int         post            = 0;
  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &nv,
              &length,
              &domain_size,
              &ng_tgt,
              &grid,
              &n_max_per_leaf,
              &on_ground,
              &n_proc_data_src,
              &post);


  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  int          *n_tgt     = NULL;
  double      **tgt_coord = NULL;
  PDM_g_num_t **tgt_g_num = NULL;

  double xmin = -0.5*domain_size;
  double origin[3] = {xmin, xmin, xmin};
  if (grid) {
    _gen_cloud_grid (PDM_MPI_COMM_WORLD,
                     n_part_tgt,
                     ng_tgt,
                     origin,
                     domain_size,
                     &n_tgt,
                     &tgt_g_num,
                     &tgt_coord);
  }

  else {
    n_part_tgt = 1;
    PDM_malloc(n_tgt, n_part_tgt, int);
    PDM_malloc(tgt_g_num, n_part_tgt, PDM_g_num_t *);
    PDM_malloc(tgt_coord, n_part_tgt, double      *);

    _gen_cloud_random (PDM_MPI_COMM_WORLD,
                       ng_tgt,
                       origin,
                       domain_size,
                       n_tgt,
                       tgt_g_num,
                       tgt_coord);
  }

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

    int *rank_in_nodes;
    PDM_malloc(rank_in_nodes, n_rank, int);

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

    PDM_free(rank_in_nodes);
  }


  int n_part_src = 1;
  int          *n_vtx     = NULL;
  double      **vtx_coord = NULL;
  PDM_g_num_t **vtx_g_num = NULL;
  int          *n_face       = NULL;
  int         **face_vtx_idx = NULL;
  int         **face_vtx     = NULL;
  PDM_g_num_t **face_g_num   = NULL;

  _gen_src_mesh (active_rank_src,
                 src_comm,
                 n_part_src,
                 nv,
                 length,
                 &n_vtx,
                 &vtx_g_num,
                 &vtx_coord,
                 &n_face,
                 &face_g_num,
                 &face_vtx_idx,
                 &face_vtx);


  PDM_g_num_t n_g_face_loc = 0;
  PDM_g_num_t n_g_vtx_loc = 0;

  PDM_g_num_t n_g_face = 0;
  PDM_g_num_t n_g_vtx = 0;

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    for (int i = 0; i < n_face[i_part]; i++) {
      n_g_face_loc = PDM_MAX (n_g_face_loc, face_g_num[i_part][i]);
    }

    for (int i = 0; i < n_vtx[i_part]; i++) {
      n_g_vtx_loc = PDM_MAX (n_g_vtx_loc, vtx_g_num[i_part][i]);
    }
  }
  PDM_MPI_Allreduce (&n_g_face_loc, &n_g_face, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);

  PDM_MPI_Allreduce (&n_g_vtx_loc, &n_g_vtx, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX,
                     PDM_MPI_COMM_WORLD);


  double z_shift;
  if (on_ground) {
    z_shift = 0.5*domain_size;
  } else {
    z_shift = 0.5*length;
  }

  for (int i_part = 0; i_part <  n_part_src; i_part++) {
    for (int i = 0; i < n_vtx[i_part]; i++) {
      vtx_coord[i_part][3*i+2] -= z_shift;
    }
  }


  if (post) {
    char filename[999];

    for (int i_part = 0; i_part < n_part_tgt; i_part++) {
      sprintf(filename, "tgt_%3.3d.vtk", n_part_tgt*i_rank + i_part);
      _write_point_cloud (filename,
                          "tgt",
                          n_tgt[i_part],
                          tgt_coord[i_part],
                          tgt_g_num[i_part],
                          NULL);
    }

    for (int i_part = 0; i_part < n_part_src; i_part++) {
      sprintf(filename, "src_mesh_%3.3d.vtk", n_part_src*i_rank + i_part);
      _write_polydata (filename,
                       "src_mesh",
                       n_vtx[i_part],
                       vtx_coord[i_part],
                       vtx_g_num[i_part],
                       n_face[i_part],
                       face_vtx_idx[i_part],
                       face_vtx[i_part],
                       face_g_num[i_part]);
    }
  }



  int n_point_cloud = 1;
  PDM_dist_cloud_surf_t *id_dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                               n_point_cloud,
                                                               PDM_MPI_COMM_WORLD,
                                                               PDM_OWNERSHIP_KEEP);

  PDM_dist_cloud_surf_surf_mesh_global_data_set (id_dist,
                                                 n_part_src);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_dist_cloud_surf_surf_mesh_part_set (id_dist,
                                            i_part,
                                            n_face[i_part],
                                            face_vtx_idx[i_part],
                                            face_vtx[i_part],
                                            face_g_num[i_part],
                                            n_vtx[i_part],
                                            vtx_coord[i_part],
                                            vtx_g_num[i_part]);
  }


  PDM_dist_cloud_surf_n_part_cloud_set (id_dist,
                                        0,
                                        n_part_tgt);

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    PDM_dist_cloud_surf_cloud_set (id_dist,
                                   0,
                                   i_part,
                                   n_tgt[i_part],
                                   tgt_coord[i_part],
                                   tgt_g_num[i_part]);
  }

  /* Compute distance */
  // PDM_dist_cloud_surf_compute (id_dist);
  PDM_dist_cloud_surf_compute (id_dist);

  PDM_dist_cloud_surf_dump_times(id_dist);



  PDM_dist_cloud_surf_free (id_dist);

  /*
   *  Finalize
   */
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    PDM_free(tgt_coord[i_part]);
    PDM_free(tgt_g_num[i_part]);
  }
  PDM_free(n_tgt);
  PDM_free(tgt_coord);
  PDM_free(tgt_g_num);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_free(vtx_coord[i_part]);
    PDM_free(vtx_g_num[i_part]);
    PDM_free(face_vtx_idx[i_part]);
    PDM_free(face_vtx[i_part]);
    PDM_free(face_g_num[i_part]);
  }
  PDM_free(n_vtx);
  PDM_free(vtx_coord);
  PDM_free(vtx_g_num);
  PDM_free(n_face);
  PDM_free(face_vtx_idx);
  PDM_free(face_vtx);
  PDM_free(face_g_num);


  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
