/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_iso_surface.h"
#include "pdm_iso_surface_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_polygon.h"
#include "pdm_plane.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static
inline
double
_plane_field
(
  double x,
  double y,
  double z,
  double *plane_equation
)
{
  return plane_equation[0] * x + plane_equation[1] * y + plane_equation[2] * z + plane_equation[3];
}

static
inline
void
_plane_gradient_field
(
  double  x,
  double  y,
  double  z,
  double *plane_equation,
  double *gradx,
  double *grady,
  double *gradz
)
{
  PDM_UNUSED(x);
  PDM_UNUSED(y);
  PDM_UNUSED(z);
  *gradx = plane_equation[0];
  *grady = plane_equation[1];
  *gradz = plane_equation[2];
}

static void
_dump_vectors
(
 const char        *filename,
 const int          n_pts,
 const double       pts_coord[],
 const double       vector[],
 const PDM_g_num_t  pts_g_num[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", pts_coord[3*i+j]);
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

  fprintf(f, "POINT_DATA %d\n", n_pts);
  fprintf(f, "VECTORS vector double\n");
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vector[3*i+j]);
    }
    fprintf(f, "\n");
  }

  if (pts_g_num != NULL) {
    fprintf(f, "FIELD pts_field 1\n");
    fprintf(f, "g_num 1 %d long\n", n_pts);
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, "1 "PDM_FMT_G_NUM"\n", pts_g_num[i]);
    }
  }

  fclose(f);
}


#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
extern void dgesvd_(const char *jobu,
                    const char *jobvt,
                    int        *m,
                    int        *n,
                    double     *a,
                    int        *lda,
                    double     *s,
                    double     *u,
                    int        *ldu,
                    double     *vt,
                    int        *ldvt,
                    double     *work,
                    int        *lwork,
                    int        *info);
#endif


static void
_solve_least_square
(
 const int  m,
 const int  n,
 double    *A,
 double    *U,
 double    *S,
 double    *V,
 double    *b,
 double    *x
 )
{
  if (1) {
    log_trace("A =\n");
    for (int i = 0; i < m; i++) {
      log_trace("[");
      for (int j = 0; j < n; j++) {
        log_trace("%20.16f ", A[i+m*j]);
      }
      log_trace("]\n");
    }

    log_trace("b = [");
    for (int i = 0; i < m; i++) {
      log_trace("%20.16f ", b[i]);
    }
    log_trace("]\n");
  }

  /* Compute SVD of A */
#if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
  int info = 0;
  int n_row = m;
  int n_col = n;
  int lwork = 100;
  double work[100];
  dgesvd_("S",
          "S",
          &n_row,
          &n_col,
          A,
          &n_row,
          S,
          U,
          &n_row,
          V,
          &n_row,
          work,
          &lwork,
          &info);
#else
    printf("Error : LAPACK or MKL are mandatory, recompile with them. \n");
    abort();
#endif

  log_trace("U =\n");
  for (int i = 0; i < m; i++) {
    log_trace("[");
    for (int j = 0; j < n; j++) {
      log_trace("%20.16f ", U[i+m*j]);
    }
    log_trace("]\n");
  }

  log_trace("S = [");
  for (int i = 0; i < n; i++) {
    log_trace("%20.16f ", S[i]);
  }
  log_trace("]\n");

  log_trace("V =\n");
  for (int i = 0; i < m; i++) {
    log_trace("[");
    for (int j = 0; j < n; j++) {
      log_trace("%20.16f ", V[i+m*j]);
    }
    log_trace("]\n");
  }

  const double tol = 1.e-1;
  double tol_s = tol * PDM_ABS(S[0]);

  /* Compute y := S^{-1} U^T b */
  double y[n];

  log_trace("y = [");
  for (int i = 0; i < n; i++) {

    y[i] = 0.;

    if (PDM_ABS(S[i]) > tol_s) {
      for (int j = 0; j < m; j++) {
        y[i] += U[j + m*i] * b[j];
      }

      y[i] /= S[i];
    }
    log_trace("%f ", y[i]);
  }
  log_trace("]\n");

  /* Compute x := V^T y */
  log_trace("x = [");
  for (int i = 0; i < n; i++) {
    x[i] = 0.;

    for (int j = 0; j < m; j++) {
      x[i] += V[j + m*i] * y[j];
    }
  log_trace("%f ", x[i]);
  }
  log_trace("]\n");
}



static
void
_compute_edge_isosurf_intersection
(
 PDM_iso_surface_t  *isos,
 const int           n_edge,
 int                *pedge_vtx,
 double             *pvtx_coord,
 double             *pfield,
 double             *pgradient_field,
 int                *edge_tag,
 double             *edge_coord,
 double             *edge_gradient
 )
{
  for (int i = 0; i < n_edge; i++) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    double z1 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];
    double z2 = pvtx_coord[3*i_vtx2+2];

    double val1 = 0.;
    double val2 = 0.;
    if(isos->iso_kind == PDM_ISO_SURFACE_KIND_PLANE) {
      val1 = _plane_field(x1, y1, z1, isos->plane_equation);
      val2 = _plane_field(x2, y2, z2, isos->plane_equation);
    } else {
      val1 = pfield[i_vtx1];
      val2 = pfield[i_vtx2];
    }

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if (sgn1 != sgn2) {
      edge_tag[i] = sgn1;

      double grad1[3], grad2[3];

      if(isos->iso_kind == PDM_ISO_SURFACE_KIND_PLANE) {
        _plane_gradient_field(x1, y1, z1, isos->plane_equation, &grad1[0], &grad1[1], &grad1[2]);
        _plane_gradient_field(x2, y2, z2, isos->plane_equation, &grad2[0], &grad2[1], &grad2[2]);
      } else {
        if (pgradient_field != NULL) {
          grad1[0] = pgradient_field[3*i_vtx1  ];
          grad1[1] = pgradient_field[3*i_vtx1+1];
          grad1[2] = pgradient_field[3*i_vtx1+2];

          grad2[0] = pgradient_field[3*i_vtx2  ];
          grad2[1] = pgradient_field[3*i_vtx2+1];
          grad2[2] = pgradient_field[3*i_vtx2+2];
        }
      }

      // Linear interpolation
      double t = val1 / (val1 - val2);

      // Cubic (Hermite) interpolation
      if (pgradient_field != NULL) {
        double vec[3] = {x2 - x1, y2 - y1, z2 - z1};
        double m0 = PDM_DOT_PRODUCT(vec, grad1);
        double m1 = PDM_DOT_PRODUCT(vec, grad2);

        // Find a root of a_3*t^3 + a_2*t^2 + a_1*t + a_0 (for 0 <= t <= 1), with
        double a_3 = 2*(val1 - val2) + m0 + m1;
        double a_2 = 3*(val2 - val1) - 2*m0 - m1;
        double a_1 = m0;
        double a_0 = val1;
        log_trace("coefs = %f %f %f %f\n", a_3, a_2, a_1, a_0);

        double s = t;
        int stat = 0;
        log_trace("Newton:\n");
        for (int iter = 0; iter < 5; iter++) {
          double val = a_0 + s*(a_1 + s*(a_2 + s*a_3));
          log_trace("  it %d, s = %f, |val| = %e\n", iter, s, PDM_ABS(val));

          if (PDM_ABS(val) < 1e-6) {
            // Converged
            stat = 1;
            break;
          }
          double dval_ds = a_1 + s*(2*a_2 + s*3*a_3);

          if (PDM_ABS(dval_ds) < 1e-12) {
            // Singular derivative
            stat = -1;
            break;
          } else {
            // Apply Newton step
            s -= val / dval_ds;
            if (s < 0. || s > 1.) {
              // Iterate outside of valid bounds
              stat = -2;
              break;
            }
          }
        }

        if (stat == 1) {
          // Newton converged successfully
          t = s;
        }
      }

      edge_coord[3*i  ] = (1. - t)*x1 + t*x2;
      edge_coord[3*i+1] = (1. - t)*y1 + t*y2;
      edge_coord[3*i+2] = (1. - t)*z1 + t*z2;

      if (pgradient_field != NULL) {
        double gx = (1. - t)*grad1[0] + t*grad2[0];
        double gy = (1. - t)*grad1[1] + t*grad2[1];
        double gz = (1. - t)*grad1[2] + t*grad2[2];
        double mag = sqrt(gx*gx + gy*gy + gz*gz);
        if (mag > 1e-12) {
          mag = 1./ mag;
          gx *= mag;
          gy *= mag;
          gz *= mag;
        }

      edge_gradient[3*i  ] = gx;
      edge_gradient[3*i+1] = gy;
      edge_gradient[3*i+2] = gz;
      }
    }

    else {
      // Unnecessary
      edge_coord[3*i  ] = 0.;
      edge_coord[3*i+1] = 0.;
      edge_coord[3*i+2] = 0.;

      if (pgradient_field != NULL) {
        edge_gradient[3*i  ] = 0.;
        edge_gradient[3*i+1] = 0.;
        edge_gradient[3*i+2] = 0.;
      }
    }

  }
}




static
void
_dump_iso_line_dist
(
 PDM_MPI_Comm  comm,
 int           isoline_n_vtx,
 double       *isoline_vtx_coord,
 PDM_g_num_t  *isoline_vtx_parent_face,
 PDM_g_num_t  *isoline_vtx_ln_to_gn,
 int           isoline_n_edge,
 int          *isoline_edge_vtx_idx,
 PDM_g_num_t  *isoline_edge_vtx,
 PDM_g_num_t  *isoline_edge_ln_to_gn
 )
{
  PDM_UNUSED(isoline_vtx_parent_face);

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  char filename[999];
  sprintf(filename, "isoline_%2.2d.vtk", i_rank);

  // check field values at isoline vertices
  for (int i = 0; i < isoline_n_vtx; i++) {
    double x = isoline_vtx_coord[3*i  ];
    double y = isoline_vtx_coord[3*i+1];
    double z = isoline_vtx_coord[3*i+2];

    double v = x*x + y*y - 0.125;
    log_trace("<d> isoline_vtx "PDM_FMT_G_NUM", coords %f %f %f, value = %f\n",
              isoline_vtx_ln_to_gn[i], x, y, z, v);
  }


  PDM_g_num_t *distrib_edge = PDM_compute_entity_distribution(comm,
                                                              isoline_n_edge);

  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      &isoline_vtx_ln_to_gn,
                                                      NULL,
                                                      &isoline_n_vtx,
                                                      1,
                                                      comm);

  double *dvtx_coord = NULL;
  PDM_part_to_block_exch (ptb,
                          3*sizeof(double),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &isoline_vtx_coord,
                          NULL,
               (void **) &dvtx_coord);

  PDM_g_num_t *distrib_vtx = PDM_part_to_block_distrib_index_get(ptb);

  int          pn_vtx        = 0;
  PDM_g_num_t *pvtx_ln_to_gn = NULL;
  int         *pedge_vtx_idx = NULL;
  int         *pedge_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_edge,
                                                           isoline_edge_vtx_idx,
                                                           isoline_edge_vtx,
                                                           isoline_n_edge,
                                     (const PDM_g_num_t *) isoline_edge_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pedge_vtx_idx,
                                                           &pedge_vtx);

  PDM_log_trace_connectivity_int(pedge_vtx_idx, pedge_vtx, isoline_n_edge, "pedge_vtx : ");

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

  log_trace("\n\n");
  for (int i = 0; i < pn_vtx; i++) {
    double x = pvtx_coord[3*i  ];
    double y = pvtx_coord[3*i+1];
    double z = pvtx_coord[3*i+2];

    double v = x*x + y*y - 0.125;
    log_trace("<p> %d isoline_vtx "PDM_FMT_G_NUM", coords %f %f %f, value = %f\n",
              i+1, pvtx_ln_to_gn[i], x, y, z, v);
  }


  free(distrib_edge);
  PDM_part_to_block_free(ptb);





  // assert(n_rank == 1);
  // int *pedge_vtx = (int *) malloc(sizeof(int) * isoline_edge_vtx_idx[isoline_n_edge]);
  // for (int i = 0; i < isoline_edge_vtx_idx[isoline_n_edge]; i++) {
  //   pedge_vtx[i] = (int) isoline_edge_vtx[i];
  // }

  PDM_vtk_write_std_elements (filename,
                              pn_vtx,
                              pvtx_coord,
                              pvtx_ln_to_gn,
                              PDM_MESH_NODAL_BAR2,
                              isoline_n_edge,
                              pedge_vtx,
                              isoline_edge_ln_to_gn,
                              0,
                              NULL,
                              NULL);


  double *pedge_center = (double *) malloc(sizeof(double) * isoline_n_edge * 3);
  double *pedge_normal = (double *) malloc(sizeof(double) * isoline_n_edge * 3);
  for (int i = 0; i < isoline_n_edge; i++) {
    for (int k = 0; k < 3; k++) {
      pedge_center[3*i + k] = 0.;
      pedge_normal[3*i + k] = 0.;
    }


    for (int j = 0; j < 2; j++) {
      int ivtx = pedge_vtx[2*i+j] - 1;
      for (int k = 0; k < 3; k++) {
        pedge_center[3*i + k] += 0.5 * pvtx_coord[3*ivtx + k];
      }

      int s = (j == 0 ? -1 : 1);
      pedge_normal[3*i  ] -= s * pvtx_coord[3*ivtx+1];
      pedge_normal[3*i+1] += s * pvtx_coord[3*ivtx  ];
    }

  }

  sprintf(filename, "isoline_normal_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 isoline_n_edge,
                 pedge_center,
                 pedge_normal,
                 NULL);
  free(pedge_center);
  free(pedge_normal);

  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(pvtx_coord);


}


static
void
_iso_line_dist
(
 PDM_iso_surface_t  *isos,
 int                 n_face,
 int                 n_edge,
 int                 n_vtx,
 int                *pface_edge_idx,
 int                *pface_edge,
 int                *pedge_vtx,
 PDM_g_num_t        *pface_ln_to_gn,
 PDM_g_num_t        *pedge_ln_to_gn,
 PDM_g_num_t        *pvtx_ln_to_gn,
 double             *pvtx_coord,
 double             *pfield,
 double             *pgradient_field,
 int                *isoline_n_vtx,
 double            **isoline_vtx_coord,
 PDM_g_num_t       **isoline_vtx_parent_face,
 PDM_g_num_t       **isoline_vtx_ln_to_gn,
 int                *isoline_n_edge,
 int               **isoline_edge_vtx_idx,
 PDM_g_num_t       **isoline_edge_vtx,
 PDM_g_num_t       **isoline_edge_ln_to_gn
)
{
  PDM_MPI_Comm comm = isos->comm;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  char filename[999];
  sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            pvtx_coord,
                            pvtx_ln_to_gn,
                            NULL);

  /*
   *  Tag edges that cross the iso-line,
   *  compute the intersection point
   *  and the gradient at that point
   */
  int    *edge_tag      = PDM_array_zeros_int(n_edge);
  double *edge_coord    = (double *) malloc(sizeof(double) * n_edge * 3);
  double *edge_gradient = (double *) malloc(sizeof(double) * n_edge * 3);

  _compute_edge_isosurf_intersection(isos,
                                     n_edge,
                                     pedge_vtx,
                                     pvtx_coord,
                                     pfield,
                                     pgradient_field,
                                     edge_tag,
                                     edge_coord,
                                     edge_gradient);


  sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 n_edge,
                 edge_coord,
                 edge_gradient,
                 NULL);


  int n_face_edge_max = 0;
  for (int i = 0; i < n_face; i++) {
    n_face_edge_max = PDM_MAX(n_face_edge_max,
                              pface_edge_idx[i+1] - pface_edge_idx[i]);
  }

  int *pface_edge_inter_idx = (int *) malloc(sizeof(int) * (n_face + 1));
  pface_edge_inter_idx[0] = 0;
  PDM_g_num_t *pface_edge_inter_g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * pface_edge_idx[n_face]);

  const int dim = 2;
  double *mat = (double *) malloc (sizeof(double) * n_face_edge_max * dim);
  double *rhs = (double *) malloc (sizeof(double) * n_face_edge_max);

  double *S = (double *) malloc(sizeof(double) * dim);
  double *U = (double *) malloc(sizeof(double) * n_face_edge_max * dim);
  double *V = (double *) malloc(sizeof(double) * n_face_edge_max * dim);

  *isoline_n_vtx = n_face;
  *isoline_vtx_coord       = (double *)      malloc(sizeof(double)      * n_face * 3);
  *isoline_vtx_parent_face = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face );

  double *poly_coord = (double *) malloc(sizeof(double) * n_face_edge_max * 3);

  for (int i = 0; i < n_face; i++) {

    (*isoline_vtx_parent_face)[i] = pface_ln_to_gn[i];

    pface_edge_inter_idx[i+1] = pface_edge_inter_idx[i];

    int n_tagged_edge = 0;
    for (int j = pface_edge_idx[i]; j < pface_edge_idx[i+1]; j++) {
      int iedge = PDM_ABS(pface_edge[j]) - 1;
      if (edge_tag[iedge] != 0) {
        pface_edge_inter_g_num[pface_edge_inter_idx[i+1]++] = pedge_ln_to_gn[iedge];
        n_tagged_edge++;
      }
    }

    int k = 0;
    for (int j = pface_edge_idx[i]; j < pface_edge_idx[i+1]; j++) {
      int iedge = PDM_ABS(pface_edge[j]) - 1;

      if (edge_tag[iedge] != 0) {

        mat[k              ] = edge_gradient[3*iedge  ];
        mat[k+n_tagged_edge] = edge_gradient[3*iedge+1];

        rhs[k] = PDM_DOT_PRODUCT(edge_gradient + 3*iedge, edge_coord + 3*iedge);

        k++;
      }
    }

    assert(k >= 2);

    // log_trace("A =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", mat[j], mat[j+n_tagged_edge]);
    // }

    // log_trace("b = [%20.16f, %20.16f]\n", rhs[0], rhs[1]);

    /* Compute SVD of mat */
// #if defined(PDM_HAVE_MKL) || defined(PDM_HAVE_LAPACK)
//     int info = 0;
//     int n_row = n_tagged_edge;
//     int n_col = 2;
//     int lwork = 15;//>= MAX(1,5*MIN(M,N))
//     double work[15];
//     dgesvd_("S",
//             "S",
//             &n_row,
//             &n_col,
//             mat,
//             &n_row,
//             S,
//             U,
//             &n_row,
//             V,
//             &n_row,
//             work,
//             &lwork,
//             &info);
// #else
//     printf("Error : LAPACK or MKL are mandatory, recompile with them. \n");
//     abort();
// #endif

    // log_trace("S = [%20.16f, %20.16f]\n", S[0], S[1]);
    // log_trace("U =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", U[j], U[j+n_tagged_edge]);
    // }
    // log_trace("V  =\n");
    // for (int j = 0; j < n_tagged_edge; j++) {
    //   log_trace("[%20.16f, %20.16f]\n", V[j], V[j+n_tagged_edge]);
    // }

    /* Solve for iso-line vertex coordinates */
    double *sol = *isoline_vtx_coord + 3*i;

    _solve_least_square (n_tagged_edge,
                           dim,
                           mat,
                           U,
                           S,
                           V,
                           rhs,
                           sol);
    sol[2] = 0.;

    if (1) {
      /* Check if the solution is inside the polygon */
      // Get poly_coord
      int cur_vtx, next_vtx;
      int cur_edge = pface_edge[pface_edge_idx[i]];
      if (cur_edge < 0) {
        cur_edge = -cur_edge - 1;
        cur_vtx  = pedge_vtx[2*cur_edge+1];
        next_vtx = pedge_vtx[2*cur_edge  ];
      } else {
        cur_edge = cur_edge - 1;
        cur_vtx  = pedge_vtx[2*cur_edge  ];
        next_vtx = pedge_vtx[2*cur_edge+1];
      }

      for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
        memcpy(poly_coord + 3*ivtx, pvtx_coord + 3*(cur_vtx - 1), sizeof(double) * 3);

        for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
          cur_edge = pface_edge[iedg];
          int vtx1, vtx2;
          if (cur_edge < 0) {
            cur_edge = -cur_edge - 1;
            vtx1 = pedge_vtx[2*cur_edge+1];
            vtx2 = pedge_vtx[2*cur_edge  ];
          } else {
            cur_edge = cur_edge - 1;
            vtx1 = pedge_vtx[2*cur_edge  ];
            vtx2 = pedge_vtx[2*cur_edge+1];
          }

          if (vtx1 == next_vtx) {
            cur_vtx  = next_vtx;
            next_vtx = vtx2;
            break;
          }
        }
      }

      double poly_bounds[6] = {
        HUGE_VAL, -HUGE_VAL,
        HUGE_VAL, -HUGE_VAL,
        HUGE_VAL, -HUGE_VAL
      };
      double poly_normal[3];
      for (int j = 0; j < pface_edge_idx[i+1] - pface_edge_idx[i]; j++) {
        for (int l = 0; l < 3; l++) {
          double x = poly_coord[3*j + l];
          poly_bounds[2*l  ] = PDM_MIN(poly_bounds[2*l  ], x);
          poly_bounds[2*l+1] = PDM_MAX(poly_bounds[2*l+1], x);
        }
      }

      PDM_plane_normal(pface_edge_idx[i+1] - pface_edge_idx[i],
                       poly_coord,
                       poly_normal);

      PDM_polygon_status_t stat = PDM_polygon_point_in(sol,
                                                       pface_edge_idx[i+1] - pface_edge_idx[i],
                                                       poly_coord,
                                                       poly_bounds,
                                                       poly_normal);

      if (stat == PDM_POLYGON_OUTSIDE) {
        for (int l = 0; l < 3; l++) {
          sol[l] = 0.;
        }

        for (int j = pface_edge_idx[i]; j < pface_edge_idx[i+1]; j++) {
          int iedge = PDM_ABS(pface_edge[j]) - 1;
          if (edge_tag[iedge] != 0) {
            for (int l = 0; l < 3; l++) {
              sol[l] += edge_coord[3*iedge + l];
            }
          }
        }

        double normalization = 1. / (double) (n_tagged_edge);
        for (int l = 0; l < 3; l++) {
          sol[l] *= normalization;
        }
      }
    }

    // log_trace("sol = [%20.16f, %20.16f]\n", sol[0], sol[1]);
  }
  free(mat);
  free(rhs);
  free(S);
  free(U);
  free(V);
  free(edge_gradient);

  pface_edge_inter_g_num = realloc(pface_edge_inter_g_num, sizeof(PDM_g_num_t) * pface_edge_inter_idx[n_face]);

  PDM_log_trace_connectivity_long(pface_edge_inter_idx, pface_edge_inter_g_num, n_face, "pface_edge_inter : ");


  sprintf(filename, "isoline_vtx_coord_no_nbd_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_face,
                            *isoline_vtx_coord,
                            pface_ln_to_gn,
                            NULL);





  int n_face_edge_inter = pface_edge_inter_idx[n_face];

  PDM_g_num_t lmax_edge_ln_to_gn = 0;
  for (int i = 0; i < n_edge; i++) {
    lmax_edge_ln_to_gn = PDM_MAX(lmax_edge_ln_to_gn,
                                 pedge_ln_to_gn[i]);
  }

  PDM_g_num_t gmax_edge_ln_to_gn;
  PDM_MPI_Allreduce(&lmax_edge_ln_to_gn, &gmax_edge_ln_to_gn, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);


  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      &pface_edge_inter_g_num,
                                                      NULL,
                                                      &n_face_edge_inter,
                                                      1,
                                                      comm);
  free(pface_edge_inter_g_num);

  PDM_g_num_t *pface_edge_inter_face_g_num = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face_edge_inter);
  for (int i = 0; i < n_face; i++) {
    for (int j = pface_edge_inter_idx[i]; j < pface_edge_inter_idx[i+1]; j++) {
      pface_edge_inter_face_g_num[j] = pface_ln_to_gn[i];
    }
  }

  free(pface_edge_inter_idx);
  int *part_stride = PDM_array_const_int(n_face_edge_inter, 1);

  *isoline_n_edge = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distrib_edge = PDM_part_to_block_distrib_index_get(ptb);
  distrib_edge[n_rank] = gmax_edge_ln_to_gn;

  PDM_part_to_block_t *ptb2 = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                    PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                    1.,
                                                                    &pedge_ln_to_gn,
                                                                    distrib_edge,
                                                                    &n_edge,
                                                                    1,
                                                                    comm);

  /* Exchange edge_coord */
  int    *isoline_dedge_coord_n = NULL;
  double *isoline_dedge_coord   = NULL;
  int *edge_stride = PDM_array_const_int(n_edge, 1);
  PDM_part_to_block_exch (ptb2,
                          3*sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &edge_stride,
                (void **) &edge_coord,
                          &isoline_dedge_coord_n,
                (void **) &isoline_dedge_coord);

  /* Exchange edge_tag */
  int *dedge_tag_n = NULL;
  int *dedge_tag   = NULL;
  PDM_part_to_block_exch (ptb2,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &edge_stride,
                (void **) &edge_tag,
                          &dedge_tag_n,
                (void **) &dedge_tag);

  PDM_g_num_t *ptb_gnum  = PDM_part_to_block_block_gnum_get (ptb);
  PDM_g_num_t *ptb2_gnum = PDM_part_to_block_block_gnum_get (ptb2);
  int ptb2_n = PDM_part_to_block_n_elt_block_get(ptb2);

  PDM_log_trace_array_long(ptb2_gnum, ptb2_n, "ptb2_gnum : ");
  PDM_log_trace_array_int(isoline_dedge_coord_n, ptb2_n, "isoline_dedge_coord_n : ");
  PDM_log_trace_array_int(dedge_tag_n, ptb2_n, "dedge_tag_n : ");
  PDM_log_trace_array_int(dedge_tag, ptb2_n, "dedge_tag : ");

  free (edge_stride);
  free (edge_coord);
  free (edge_tag);
  free (isoline_dedge_coord_n);


  int *isoline_dedge_vtx_n = NULL;
  PDM_g_num_t *tmp_isoline_edge_vtx = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &part_stride,
                (void **) &pface_edge_inter_face_g_num,
                          &isoline_dedge_vtx_n,
                (void **) &tmp_isoline_edge_vtx);

  PDM_log_trace_array_long(ptb_gnum, *isoline_n_edge, "ptb_gnum : ");
  PDM_log_trace_array_int(isoline_dedge_vtx_n, *isoline_n_edge, "isoline_dedge_vtx_n : ");

  *isoline_edge_vtx = malloc(sizeof(PDM_g_num_t) * (*isoline_n_edge) * 2);

  int isoline_n_bnd_vtx = 0;
  for (int i = 0; i < *isoline_n_edge; i++) {
    if (isoline_dedge_vtx_n[i] == 1) {
      /* Boundary edge */
      isoline_n_bnd_vtx++;
    }
  }

  PDM_g_num_t *distrib_face    = PDM_compute_entity_distribution(comm,
                                                                 n_face);
  PDM_g_num_t *distrib_bnd_vtx = PDM_compute_entity_distribution(comm,
                                                                 isoline_n_bnd_vtx);

  *isoline_n_vtx += isoline_n_bnd_vtx;
  *isoline_vtx_parent_face = realloc(*isoline_vtx_parent_face,
                                      sizeof(PDM_g_num_t) * (*isoline_n_vtx));
  *isoline_vtx_coord       = realloc(*isoline_vtx_coord,
                                      sizeof(double) * (*isoline_n_vtx) * 3);

  *isoline_vtx_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*isoline_n_vtx));
  memcpy(*isoline_vtx_ln_to_gn, pface_ln_to_gn, sizeof(PDM_g_num_t) * n_face);

  isoline_n_bnd_vtx = 0;
  int idx = 0;
  for (int i = 0; i < *isoline_n_edge; i++) {

    for (int j = 0; j < isoline_dedge_vtx_n[i]; j++) {
      (*isoline_edge_vtx)[2*i+j] = tmp_isoline_edge_vtx[idx+j];
    }

    if (isoline_dedge_vtx_n[i] == 1) {
      /* Boundary edge */
      (*isoline_edge_vtx)[2*i+1] = distrib_face[n_rank] + distrib_bnd_vtx[i_rank] + isoline_n_bnd_vtx + 1;
      log_trace("bnd_vtx %d : g_num = "PDM_FMT_G_NUM"\n", isoline_n_bnd_vtx, (*isoline_edge_vtx)[2*i+1]);

      int idx2 = n_face + isoline_n_bnd_vtx;

      int pos = PDM_binary_search_long(ptb_gnum[i],
                                       ptb2_gnum,
                                       ptb2_n);
      assert(pos >= 0);
      memcpy((*isoline_vtx_coord) + 3*idx2, isoline_dedge_coord + 3*pos, sizeof(double) * 3);

      (*isoline_vtx_ln_to_gn)[idx2] = (*isoline_edge_vtx)[2*i+1];

      if (dedge_tag[pos] < 0) {
        PDM_g_num_t tmp = (*isoline_edge_vtx)[2*i];
        (*isoline_edge_vtx)[2*i] = (*isoline_edge_vtx)[2*i+1];
        (*isoline_edge_vtx)[2*i+1] = tmp;
      }

      (*isoline_vtx_parent_face)[idx2] = tmp_isoline_edge_vtx[idx];
      isoline_n_bnd_vtx++;
    }

    idx += isoline_dedge_vtx_n[i];

  }
  free(isoline_dedge_coord);
  PDM_part_to_block_free(ptb2);


  sprintf(filename, "isoline_vtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            *isoline_n_vtx,
                            *isoline_vtx_coord,
                            *isoline_vtx_ln_to_gn,
                            NULL);


  if (0) {
    double *pface_edge_inter_face_coord = (double *) malloc(sizeof(double) * n_face_edge_inter * 3);
    for (int i = 0; i < n_face; i++) {
      for (int j = pface_edge_inter_idx[i]; j < pface_edge_inter_idx[i+1]; j++) {
        for (int k = 0; k < 3; k++) {
          pface_edge_inter_face_coord[3*j+k] = (*isoline_vtx_coord)[3*i+k];
        }
      }
    }

    int    *isoline_dedge_vtx_coord_n = NULL;
    double *isoline_dedge_vtx_coord   = NULL;
    PDM_part_to_block_exch (ptb,
                            3*sizeof(double),
                            PDM_STRIDE_VAR_INTERLACED,
                            1,
                            &part_stride,
                  (void **) &pface_edge_inter_face_coord,
                            &isoline_dedge_vtx_coord_n,
                  (void **) &isoline_dedge_vtx_coord);

    PDM_log_trace_array_int(isoline_dedge_vtx_coord_n, *isoline_n_edge, "isoline_dedge_vtx_coord_n : ");

    sprintf(filename, "isoline_dedge_%2.2d.vtk", i_rank);
    PDM_vtk_write_lines (filename,
                         *isoline_n_edge,
                         isoline_dedge_vtx_coord,
                         NULL,
                         NULL);
    free(isoline_dedge_vtx_coord_n);
    free(isoline_dedge_vtx_coord  );
    free(pface_edge_inter_face_coord);
  }
  free(pface_edge_inter_face_g_num);

  free(part_stride);

  *isoline_edge_vtx_idx = (int *) malloc(sizeof(int) * (*isoline_n_edge + 1));
  for (int i = 0; i <= *isoline_n_edge; i++) {
    (*isoline_edge_vtx_idx)[i] = 2*i;
  }
  free(isoline_dedge_vtx_n);

  PDM_log_trace_connectivity_long(*isoline_edge_vtx_idx, *isoline_edge_vtx, *isoline_n_edge, "isoline_edge_vtx : ");

  PDM_g_num_t *distrib_isoline_edge = PDM_compute_entity_distribution(comm,
                                                                      *isoline_n_edge);

  (*isoline_edge_ln_to_gn) = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*isoline_n_edge));
  for (int i = 0; i < *isoline_n_edge; i++) {
    (*isoline_edge_ln_to_gn)[i] = distrib_isoline_edge[i_rank] + i + 1;
  }
  free(distrib_isoline_edge);

  PDM_part_to_block_free(ptb);

  free(distrib_face);
  free(distrib_bnd_vtx);
}


static void
_compute_face_vtx
(
 const int   n_face,
 int        *pface_edge_idx,
 int        *pface_edge,
 int        *pedge_vtx,
 int       **pface_vtx
 )
{
  *pface_vtx = (int *) malloc(sizeof(int) * pface_edge_idx[n_face]);

  for (int i = 0; i < n_face; i++) {

    int *_pface_vtx = *pface_vtx + pface_edge_idx[i];

    int cur_vtx, next_vtx;
    int cur_edge = pface_edge[pface_edge_idx[i]];
    if (cur_edge < 0) {
      cur_edge = -cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge+1];
      next_vtx = pedge_vtx[2*cur_edge  ];
    } else {
      cur_edge = cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge  ];
      next_vtx = pedge_vtx[2*cur_edge+1];
    }

    for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
      _pface_vtx[ivtx] = cur_vtx;

      for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
        cur_edge = pface_edge[iedg];
        int vtx1, vtx2;
        if (cur_edge < 0) {
          cur_edge = -cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge+1];
          vtx2 = pedge_vtx[2*cur_edge  ];
        } else {
          cur_edge = cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge  ];
          vtx2 = pedge_vtx[2*cur_edge+1];
        }

        if (vtx1 == next_vtx) {
          cur_vtx  = next_vtx;
          next_vtx = vtx2;
          break;
        }
      }
    }

  }
}


static
void
_iso_surf_dist
(
  PDM_iso_surface_t  *isos,
  int                 n_cell,
  int                 n_face,
  int                 n_edge,
  int                 n_vtx,
  int                *pcell_face_idx,
  int                *pcell_face,
  int                *pface_edge_idx,
  int                *pface_edge,
  int                *pedge_vtx,
  PDM_g_num_t        *pcell_ln_to_gn,
  PDM_g_num_t        *pface_ln_to_gn,
  PDM_g_num_t        *pedge_ln_to_gn,
  PDM_g_num_t        *pvtx_ln_to_gn,
  double             *pvtx_coord,
  double             *pfield,
  double             *pgradient_field,
  int                *isoline_n_vtx,
  double            **isoline_vtx_coord,
  PDM_g_num_t       **isoline_vtx_parent_cell,
  PDM_g_num_t       **isoline_vtx_ln_to_gn,
  int                *isoline_n_face,
  int               **isoline_face_vtx_idx,
  PDM_g_num_t       **isoline_face_vtx,
  PDM_g_num_t       **isoline_face_ln_to_gn
)
{
  const int use_gradient = (pgradient_field != NULL);

  PDM_MPI_Comm comm = isos->comm;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  char filename[999];
  sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            pvtx_coord,
                            pvtx_ln_to_gn,
                            NULL);

  /*
   *  Tag edges that cross the iso-line,
   *  compute the intersection point
   *  and the gradient at that point
   */
  int    *edge_tag      = PDM_array_zeros_int(n_edge);
  double *edge_coord    = (double *) malloc(sizeof(double) * n_edge * 3);
  double *edge_gradient = (double *) malloc(sizeof(double) * n_edge * 3);

  _compute_edge_isosurf_intersection(isos,
                                     n_edge,
                                     pedge_vtx,
                                     pvtx_coord,
                                     pfield,
                                     pgradient_field,
                                     edge_tag,
                                     edge_coord,
                                     edge_gradient);


  sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 n_edge,
                 edge_coord,
                 edge_gradient,
                 NULL);

  int n_face_edge_max = 0;
  for (int i = 0; i < n_face; i++) {
    n_face_edge_max = PDM_MAX(n_face_edge_max,
                              pface_edge_idx[i+1] - pface_edge_idx[i]);
  }

  int n_cell_edge_max = 0;
  for (int i = 0; i < n_cell; i++) {
    int n_cell_edge = 0;

    for (int j = pcell_face_idx[i]; j < pcell_face_idx[i+1]; j++) {
      int iface = PDM_ABS(pcell_face[j]) - 1;
      n_cell_edge += pface_edge_idx[iface+1] - pface_edge_idx[iface];
    }

    n_cell_edge_max = PDM_MAX(n_cell_edge_max, n_cell_edge);
  }
  n_cell_edge_max /= 2;


  int *is_treated = PDM_array_zeros_int(n_face);
  int *face_tag = PDM_array_zeros_int(n_face);

  const int dim = 3;
  double *mat_cell = (double *) malloc (sizeof(double) * n_cell_edge_max * dim);
  double *rhs_cell = (double *) malloc (sizeof(double) * n_cell_edge_max);
  double *mat_face = (double *) malloc (sizeof(double) * n_face_edge_max * dim);
  double *rhs_face = (double *) malloc (sizeof(double) * n_face_edge_max);

  double *S_cell = (double *) malloc(sizeof(double) * dim);
  double *U_cell = (double *) malloc(sizeof(double) * n_cell_edge_max * dim);
  double *V_cell = (double *) malloc(sizeof(double) * n_cell_edge_max * dim);
  double *S_face = (double *) malloc(sizeof(double) * dim);
  double *U_face = (double *) malloc(sizeof(double) * n_face_edge_max * dim);
  double *V_face = (double *) malloc(sizeof(double) * n_face_edge_max * dim);

  int *is_visited    = PDM_array_zeros_int(n_edge);
  int *visited_edges = (int *) malloc(sizeof(int) * n_edge);

  double *face_coord = (double *) malloc(sizeof(double) * n_face * 3);
  double *cell_coord = (double *) malloc(sizeof(double) * n_cell * 3);

  double face_center[3], face_normal[3];
  double *all_face_center = (double *) malloc(sizeof(double) * n_face * 3);
  double *all_face_normal = (double *) malloc(sizeof(double) * n_face * 3);

  // tmp -->>
  for (int i = 0; i < 3*n_cell; i++) {
    cell_coord[i] = 0.;
  }
  for (int i = 0; i < 3*n_face; i++) {
    face_coord[i] = 0.;
    all_face_center[i] = 0.;
    all_face_normal[i] = 0.;
  }
  //<<--

  for (int icell = 0; icell < n_cell; icell++) {

    int n_cell_edge = 0;
    double _cell_coord[3] = {0.};

    //-->>
    for (int i = 0; i < n_cell_edge_max * dim; i++) {
      mat_cell[i] = 123456789.;
    }
    for (int i = 0; i < n_cell_edge_max ; i++) {
     rhs_cell[i] = 987654321.;
    }
    //<<--

    // First loop to count number of tagged edges
    for (int jface = pcell_face_idx[icell]; jface < pcell_face_idx[icell+1]; jface++) {

      int iface = PDM_ABS(pcell_face[jface]) - 1;

      for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {

        int iedge = PDM_ABS(pface_edge[jedge]) - 1;

        if (edge_tag[iedge] != 0) {
          n_cell_edge++;
        }
      } // end of loop on edges of current face

    } // end of loop on faces of current cell


    log_trace("_n_cell_edge = %d\n", n_cell_edge);
    n_cell_edge /= 2;

    // Second loop to fill matrices
    int i_cell_edge = 0;
    for (int jface = pcell_face_idx[icell]; jface < pcell_face_idx[icell+1]; jface++) {

      int iface = PDM_ABS(pcell_face[jface]) - 1;
      int n_face_edge = 0;
      double _face_coord[3] = {0.};

      if (is_treated[iface] == 0) {
        for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {
          int iedge = PDM_ABS(pface_edge[jedge]) - 1;

          if (edge_tag[iedge] != 0) {
            n_face_edge++;
          }
        }

        n_face_edge++;
      }

      int i_face_edge = 0;
      for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {

        int iedge = PDM_ABS(pface_edge[jedge]) - 1;

        if (edge_tag[iedge] != 0) {

          double _rhs = PDM_DOT_PRODUCT(edge_gradient + 3*iedge, edge_coord + 3*iedge);

          if (is_visited[iedge] == 0) {
            if (use_gradient) {
              for (int l = 0; l < dim; l++) {
                mat_cell[i_cell_edge + l*n_cell_edge] = edge_gradient[3*iedge + l];
              }
              rhs_cell[i_cell_edge] = _rhs;
            }
            // else {
              for (int l = 0; l < dim; l++) {
                _cell_coord[l] += edge_coord[3*iedge + l];
              }
            // }

            visited_edges[i_cell_edge] = iedge;
            is_visited[iedge] = 1;
            i_cell_edge++;
          }

          if (is_treated[iface] == 0) {
            if (use_gradient) {
              for (int l = 0; l < dim; l++) {
                mat_face[i_face_edge + l*n_face_edge] = edge_gradient[3*iedge + l];
              }
              rhs_face[i_face_edge] = _rhs;
            }
            // else {
              for (int l = 0; l < dim; l++) {
                _face_coord[l] += edge_coord[3*iedge + l];
              }
            // }
            i_face_edge++;
          }

        }
      } // end of loop on edges of current face

      if (is_treated[iface] == 0) {

        // compute face center
        face_center[0] = face_center[1] = face_center[2] = 0.;
        for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {
          int iedge = pface_edge[jedge];
          int ivtx1;

          if (iedge < 0) {
            iedge = -iedge - 1;
            ivtx1 = pedge_vtx[2*iedge+1] - 1;
          } else {
            iedge = iedge - 1;
            ivtx1 = pedge_vtx[2*iedge  ] - 1;
          }

          for (int l = 0; l < dim; l++) {
            face_center[l] += pvtx_coord[3*ivtx1 + l];
          }
        }

        double normalization = 1. / (double) (pface_edge_idx[iface+1] - pface_edge_idx[iface]);
        for (int l = 0; l < dim; l++) {
          face_center[l] *= normalization;
        }

        // compute face normal
        face_normal[0] = face_normal[1] = face_normal[2] = 0.;
        for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {
          int iedge = pface_edge[jedge];
          int ivtx1, ivtx2;

          if (iedge < 0) {
            iedge = -iedge - 1;
            ivtx1 = pedge_vtx[2*iedge+1] - 1;
            ivtx2 = pedge_vtx[2*iedge  ] - 1;
          } else {
            iedge = iedge - 1;
            ivtx1 = pedge_vtx[2*iedge  ] - 1;
            ivtx2 = pedge_vtx[2*iedge+1] - 1;
          }

          double vec1[3] = {
            pvtx_coord[3*ivtx1  ] - face_center[0],
            pvtx_coord[3*ivtx1+1] - face_center[1],
            pvtx_coord[3*ivtx1+2] - face_center[2],
          };

          double vec2[3] = {
            pvtx_coord[3*ivtx2  ] - face_center[0],
            pvtx_coord[3*ivtx2+1] - face_center[1],
            pvtx_coord[3*ivtx2+2] - face_center[2],
          };

          double vec1xvec2[3];
          PDM_CROSS_PRODUCT(vec1xvec2, vec1, vec2);

          for (int l = 0; l < dim; l++) {
            face_normal[l] += vec1xvec2[l];
          }
        }
        double mag = PDM_MODULE(face_normal);
        if (PDM_ABS(mag) > 1.e-16) {
          mag = 1./ mag;
          for (int l = 0; l < dim; l++) {
            face_normal[l] *= mag;
          }
        }

        memcpy(all_face_center + 3*iface, face_center, sizeof(double)*3);
        memcpy(all_face_normal + 3*iface, face_normal, sizeof(double)*3);

        for (int l = 0; l < dim; l++) {
          mat_face[i_face_edge + l*n_face_edge] = face_normal[l];
        }
        rhs_face[i_face_edge] = PDM_DOT_PRODUCT(face_center, face_normal);
        i_face_edge++;

        log_trace("face %d ("PDM_FMT_G_NUM"), i_face_edge = %d\n",
                  iface, pface_ln_to_gn[iface], i_face_edge);
        if (1) {
          log_trace("  face_edge_tag :");
          for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {
            int iedge = PDM_ABS(pface_edge[jedge]) - 1;
            log_trace(" %d", edge_tag[iedge]);
          }
          log_trace("\n");

          log_trace("  face_center = %f %f %f\n",
                    face_center[0], face_center[1], face_center[2]);
          log_trace("  face_normal = %f %f %f\n",
                    face_normal[0], face_normal[1], face_normal[2]);
        }

        if (i_face_edge > 1) {
          if (use_gradient) {
            _solve_least_square (i_face_edge,
                                 dim,
                                 mat_face,
                                 U_face,
                                 S_face,
                                 V_face,
                                 rhs_face,
                                 face_coord + 3*iface);

            if (PDM_ABS(S_face[dim-1]) < 1e-3 * PDM_ABS(S_face[0])) {
              // rank-deficient matrix
              normalization = 1. / (double) (i_face_edge - 1);
              for (int l = 0; l < dim; l++) {
                face_coord[3*iface+l] = _face_coord[l] * normalization;
              }
            }

          } else {
            normalization = 1. / (double) (i_face_edge - 1);
            for (int l = 0; l < dim; l++) {
              face_coord[3*iface+l] = _face_coord[l] * normalization;
            }
          }

          face_tag[iface] = 1;

          log_trace("face_coord = %f %f %f\n",
              face_coord[3*iface], face_coord[3*iface+1], face_coord[3*iface+2]);
        }
        else {
          face_tag[iface] = 0;
        }

        is_treated[iface] = 1;
      }



    } // end of loop on faces of current cell

    log_trace("cell %d ("PDM_FMT_G_NUM"), i_cell_edge = %d\n",
              icell, pcell_ln_to_gn[icell], i_cell_edge);
    assert(i_cell_edge >= 3);
    if (use_gradient) {
      _solve_least_square (i_cell_edge,
                           dim,
                           mat_cell,
                           U_cell,
                           S_cell,
                           V_cell,
                           rhs_cell,
                           cell_coord + 3*icell);

      // Beware of NaNs
      int singular = 0;
      for (int l = 0; l < dim; l++) {
        if (!(cell_coord[3*icell+l] >= 0) &&
            !(cell_coord[3*icell+l] < 0)) {
          singular = 1;
          break;
        }
      }

      if (singular) {
        double normalization = 1. / (double) i_cell_edge;
        for (int l = 0; l < dim; l++) {
          cell_coord[3*icell+l] = _cell_coord[l] * normalization;
        }
      }

    } else {
      double normalization = 1. / (double) i_cell_edge;
      for (int l = 0; l < dim; l++) {
        cell_coord[3*icell+l] = _cell_coord[l] * normalization;
      }
    }
    log_trace("cell_coord = %f %f %f\n",
              cell_coord[3*icell], cell_coord[3*icell+1], cell_coord[3*icell+2]);

    for (int i = 0; i < i_cell_edge; i++) {
      is_visited[visited_edges[i]] = 0;
    }

  } // end of loop on cells


  if (1) {
    int *pface_vtx = NULL;
    _compute_face_vtx(n_face,
                      pface_edge_idx,
                      pface_edge,
                      pedge_vtx,
                      &pface_vtx);

    sprintf(filename, "out_equi_faces_%2.2d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           pvtx_coord,
                           pvtx_ln_to_gn,
                           n_face,
                           pface_edge_idx,
                           pface_vtx,
                           pface_ln_to_gn,
                           face_tag);
    free(pface_vtx);
  }

  sprintf(filename, "face_normal_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 n_face,
                 all_face_center,
                 all_face_normal,
                 NULL);
  free(all_face_center);
  free(all_face_normal);

  sprintf(filename, "face_intersection_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_face,
                            face_coord,
                            NULL,//pface_ln_to_gn,
                            face_tag);

  sprintf(filename, "cell_intersection_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            n_cell,
                            cell_coord,
                            pcell_ln_to_gn,
                            NULL);

  free (is_visited);
  free (visited_edges);

  free (mat_cell);
  free (rhs_cell);
  free (mat_face);
  free (rhs_face);

  free (is_treated);
  free (face_coord);
  free (cell_coord);
  free (face_tag);
}

static
void
_iso_surface_dist
(
  PDM_iso_surface_t        *isos
)
{
  PDM_UNUSED(isos);

  /*
   * Select gnum that contains iso-surface
   */
  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  char filename[999];
  sprintf(filename, "gradient_field_%2.2d.vtk", i_rank);
  _dump_vectors (filename,
                 (int) isos->distrib_vtx[i_rank+1] - isos->distrib_vtx[i_rank],
                 isos->dvtx_coord,
                 isos->dgradient_field,
                 NULL);

  assert(isos->distrib_edge != NULL);
  int dn_edge = isos->distrib_edge[i_rank+1] - isos->distrib_edge[i_rank];

  PDM_g_num_t* edge_ln_to_gn = (PDM_g_num_t * ) malloc( dn_edge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_edge; ++i) {
    edge_ln_to_gn[i] = isos->distrib_edge[i_rank] + i + 1;
  }


  int          pn_vtx           = 0;
  PDM_g_num_t *pvtx_ln_to_gn    = NULL;
  int         *pedge_vtx_idx    = NULL;
  int         *pedge_vtx        = NULL;

  int* dedge_vtx_idx = malloc( (dn_edge + 1) * sizeof(int));
  dedge_vtx_idx[0] = 0;
  for(int i = 0; i < dn_edge; ++i) {
    dedge_vtx_idx[i+1] = dedge_vtx_idx[i] + 2;
  }

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           isos->distrib_edge,
                                                           dedge_vtx_idx,
                                                           isos->dedge_vtx,
                                                           dn_edge,
                                     (const PDM_g_num_t *) edge_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pedge_vtx_idx,
                                                           &pedge_vtx);
  if(0 == 1) {
    PDM_log_trace_connectivity_int(pedge_vtx_idx, pedge_vtx, dn_edge, "pedge_vtx");
  }
  free(pedge_vtx_idx);
  free(edge_ln_to_gn);

  PDM_block_to_part_t* btp_vtx = PDM_block_to_part_create(isos->distrib_vtx,
                                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                                          &pn_vtx,
                                                          1,
                                                          isos->comm);

  int cst_stride = 1;
  double **tmp_pvtx_coord = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
            (void *  )   isos->dvtx_coord,
            (int  ***)   NULL,
            (void ***)  &tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);

  double* pfield = NULL;
  if(isos->iso_kind == PDM_ISO_SURFACE_KIND_FIELD) {
    double **tmp_pfield = NULL;
    PDM_block_to_part_exch(btp_vtx,
                           sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &cst_stride,
              (void *  )   isos->dfield,
              (int  ***)   NULL,
              (void ***)  &tmp_pfield);
    pfield = tmp_pfield[0];
    free(tmp_pfield);
  }

  PDM_block_to_part_free(btp_vtx);
  btp_vtx = NULL;


  /*
   *  Loop on edge to tag all edge
   */
  int    *dedge_tag    = (int    * ) malloc(    dn_edge * sizeof(int   ));
  double *dedge_center = (double * ) malloc(3 * dn_edge * sizeof(double));
  for(int i = 0; i < dn_edge; ++i) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    dedge_tag[i] = 0;

    // Besoin des coordonnés si call back
    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    double z1 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];
    double z2 = pvtx_coord[3*i_vtx1+2];

    double val1 = 0;
    double val2 = 0;
    if(isos->iso_kind == PDM_ISO_SURFACE_KIND_PLANE) {
      val1 = _plane_field(x1, y1, z1, isos->plane_equation);
      val2 = _plane_field(x2, y2, z2, isos->plane_equation);
    } else {
      val1 = pfield[i_vtx1];
      val2 = pfield[i_vtx2];
    }

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if(sgn1 * sgn2 < 0) {
      dedge_tag[i] = 1;
    }

    dedge_center[3*i  ] = 0.5 * (x1 + x2);
    dedge_center[3*i+1] = 0.5 * (y1 + y2);
    dedge_center[3*i+2] = 0.5 * (z1 + z2);
  }

  free(pvtx_ln_to_gn);
  free(pedge_vtx);
  free(pvtx_coord);
  free(pfield);

  PDM_g_num_t *dentity_edge     = NULL;
  int         *dentity_edge_idx = NULL;
  PDM_g_num_t *distrib_entity   = NULL;
  int          dn_entity        = -1;

  if(isos->dim == 3) {
    PDM_deduce_combine_connectivity(isos->comm,
                                    isos->distrib_cell,
                                    isos->distrib_face,
                                    isos->dcell_face_idx,
                                    isos->dcell_face,
                                    isos->dface_edge_idx,
                                    isos->dface_edge,
                                    1,
                                    &dentity_edge_idx,
                                    &dentity_edge);
    distrib_entity = isos->distrib_cell;
  } else {
    dentity_edge     = isos->dface_edge;
    dentity_edge_idx = isos->dface_edge_idx;
    distrib_entity   = isos->distrib_face;
  }

  dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  /*
   *  Deduce all entity concerns by the iso surface (could be optimize)
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(isos->distrib_edge,
                               (const PDM_g_num_t **) &dentity_edge,
                                                      &dentity_edge_idx[dn_entity],
                                                      1,
                                                      isos->comm);

  int strid_one = 1;
  int **tmp_dentity_edge_tag = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_tag,
            (int  ***)   NULL,
            (void ***)  &tmp_dentity_edge_tag);
  int *dentity_edge_tag = tmp_dentity_edge_tag[0];
  free(tmp_dentity_edge_tag);
  free(dedge_tag);

  double **tmp_dentity_edge_center = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_center,
            (int  ***)   NULL,
            (void ***)  &tmp_dentity_edge_center);
  double *dentity_edge_center = tmp_dentity_edge_center[0];
  free(tmp_dentity_edge_center);
  free(dedge_center);

  int         *dentity_tag            = malloc(     dn_entity * sizeof(int        ));
  PDM_g_num_t *entity_to_extract_gnum = malloc(     dn_entity * sizeof(PDM_g_num_t));
  double      *dentity_center         = malloc( 3 * dn_entity * sizeof(double     ));
  int  n_entity_tag = 0;
  int idx_write   = 0;
  for(int i = 0; i < dn_entity; ++i) {
    dentity_tag[i] = 0;

    for(int idx_entity = dentity_edge_idx[i]; idx_entity < dentity_edge_idx[i+1]; ++idx_entity) {
      if(dentity_edge_tag[idx_entity] == 1) {
        dentity_tag[i] = 1;
        entity_to_extract_gnum[n_entity_tag++] = distrib_entity[i_rank] + i + 1;
        break;
      }
    }

    if(dentity_tag[i] == 1) {
      dentity_center[3*idx_write  ] = 0.;
      dentity_center[3*idx_write+1] = 0.;
      dentity_center[3*idx_write+2] = 0.;

      double inv = 1./((double) (dentity_edge_idx[i+1] - dentity_edge_idx[i]));
      for(int idx_entity = dentity_edge_idx[i]; idx_entity < dentity_edge_idx[i+1]; ++idx_entity) {
        dentity_center[3*idx_write  ] += dentity_edge_center[3*idx_entity  ];
        dentity_center[3*idx_write+1] += dentity_edge_center[3*idx_entity+1];
        dentity_center[3*idx_write+2] += dentity_edge_center[3*idx_entity+2];
      }
      dentity_center[3*idx_write  ] = dentity_center[3*idx_write  ] * inv;
      dentity_center[3*idx_write+1] = dentity_center[3*idx_write+1] * inv;
      dentity_center[3*idx_write+2] = dentity_center[3*idx_write+2] * inv;

      idx_write++;
    }
  }
  free(dentity_edge_center);
  PDM_block_to_part_free(btp);

  free(dentity_tag);

  free(dentity_edge_tag);


  /*
   * Rebuild partition that contains entity and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_entity_tag, dentity_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_entity_gnum = PDM_gnum_get(gnum_equi, 0);

  if(1 == 1) {
    // char filename[999];
    sprintf(filename, "out_iso_surf_equi_entity_coord_%2.2d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_entity_tag,
                              dentity_center,
                              child_equi_entity_gnum,
                              NULL);
  }

  PDM_gnum_free(gnum_equi);
  free(dentity_center);

  /*
   * Equilibrage avec le part_to_block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &child_equi_entity_gnum,
                                                       NULL,
                                                       &n_entity_tag,
                                                       1,
                                                       isos->comm);

  int n_entity_equi = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_entity_equi_child_g_num = PDM_part_to_block_block_gnum_get (ptb);

  PDM_g_num_t *block_entity_equi_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &entity_to_extract_gnum,
                          NULL,
               (void **) &block_entity_equi_parent_g_num);

  free(entity_to_extract_gnum);

  /*
   * A refléchir pour le 3D
   */
  int          pn_face_equi               = 0;
  PDM_g_num_t *pequi_parent_face_ln_to_gn = NULL;
  int         *pequi_cell_face_idx        = NULL;
  int         *pequi_cell_face            = NULL;
  PDM_g_num_t *pequi_face_ln_to_gn        = NULL;

  if(isos->dim == 3) {
    PDM_g_num_t *block_cell_equi_parent_g_num = block_entity_equi_parent_g_num;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                             isos->distrib_cell,
                                                             isos->dcell_face_idx,
                                                             isos->dcell_face,
                                                             n_entity_equi,
                                       (const PDM_g_num_t *) block_cell_equi_parent_g_num,
                                                             &pn_face_equi,
                                                             &pequi_parent_face_ln_to_gn,
                                                             &pequi_cell_face_idx,
                                                             &pequi_cell_face);

    PDM_gen_gnum_t* gnum_face = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
    PDM_gnum_set_from_parents(gnum_face, 0, pn_face_equi, pequi_parent_face_ln_to_gn);
    PDM_gnum_compute(gnum_face);
    pequi_face_ln_to_gn = PDM_gnum_get(gnum_face, 0);
    PDM_gnum_free(gnum_face);

  } else {
    pn_face_equi               = n_entity_equi;
    pequi_parent_face_ln_to_gn = block_entity_equi_parent_g_num;
  }

  int          pn_edge_equi               = 0;
  PDM_g_num_t *pequi_parent_edge_ln_to_gn = NULL;
  int         *pequi_face_edge_idx        = NULL;
  int         *pequi_face_edge            = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           isos->distrib_face,
                                                           isos->dface_edge_idx,
                                                           isos->dface_edge,
                                                           pn_face_equi,
                                     (const PDM_g_num_t *) pequi_parent_face_ln_to_gn,
                                                           &pn_edge_equi,
                                                           &pequi_parent_edge_ln_to_gn,
                                                           &pequi_face_edge_idx,
                                                           &pequi_face_edge);

  PDM_gen_gnum_t* gnum_edge = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_edge, 0, pn_edge_equi, pequi_parent_edge_ln_to_gn);
  PDM_gnum_compute(gnum_edge);
  PDM_g_num_t* pequi_edge_ln_to_gn = PDM_gnum_get(gnum_edge, 0);
  PDM_gnum_free(gnum_edge);

  int          pn_vtx_equi               = 0;
  PDM_g_num_t *pequi_parent_vtx_ln_to_gn = NULL;
  int         *pequi_edge_vtx_idx        = NULL;
  int         *pequi_edge_vtx            = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           isos->distrib_edge,
                                                           dedge_vtx_idx,
                                                           isos->dedge_vtx,
                                                           pn_edge_equi,
                                     (const PDM_g_num_t *) pequi_parent_edge_ln_to_gn,
                                                           &pn_vtx_equi,
                                                           &pequi_parent_vtx_ln_to_gn,
                                                           &pequi_edge_vtx_idx,
                                                           &pequi_edge_vtx);
  free(dedge_vtx_idx);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_vtx, 0, pn_vtx_equi, pequi_parent_vtx_ln_to_gn);
  PDM_gnum_compute(gnum_vtx);
  PDM_g_num_t* pequi_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx, 0);
  PDM_gnum_free(gnum_vtx);

  assert(btp_vtx == NULL);
  btp_vtx = PDM_block_to_part_create(isos->distrib_vtx,
              (const PDM_g_num_t **) &pequi_parent_vtx_ln_to_gn,
                                     &pn_vtx_equi,
                                     1,
                                     isos->comm);

  double **tmp_pequi_vtx_coord = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
            (void *  )   isos->dvtx_coord,
            (int  ***)   NULL,
            (void ***)   &tmp_pequi_vtx_coord);
  double* pequi_vtx_coord = tmp_pequi_vtx_coord[0];
  free(tmp_pequi_vtx_coord);

  double* pequi_field          = NULL;
  double* pequi_gradient_field = NULL;
  if(isos->iso_kind == PDM_ISO_SURFACE_KIND_FIELD) {
    double **tmp_pequi_pfield = NULL;
    PDM_block_to_part_exch(btp_vtx,
                           sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &cst_stride,
              (void *  )   isos->dfield,
              (int  ***)   NULL,
              (void ***)   &tmp_pequi_pfield);
    pequi_field = tmp_pequi_pfield[0];
    free(tmp_pequi_pfield);


    double **tmp_pequi_pgradient_field = NULL;
    PDM_block_to_part_exch(btp_vtx,
                           3 * sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &cst_stride,
              (void *  )   isos->dgradient_field,
              (int  ***)   NULL,
              (void ***)   &tmp_pequi_pgradient_field);
    pequi_gradient_field = tmp_pequi_pgradient_field[0];
    free(tmp_pequi_pgradient_field);
  }
  PDM_block_to_part_free(btp_vtx);



  int          isoline_n_vtx           = 0;
  double      *isoline_vtx_coord       = NULL;
  PDM_g_num_t *isoline_vtx_parent_face = NULL;
  PDM_g_num_t *isoline_vtx_ln_to_gn    = NULL;
  int          isoline_n_edge          = 0;
  int         *isoline_edge_vtx_idx    = NULL;
  PDM_g_num_t *isoline_edge_vtx        = NULL;
  PDM_g_num_t *isoline_edge_ln_to_gn   = NULL;
  if(isos->dim == 2) {
    _iso_line_dist(isos,
                   n_entity_equi,
                   pn_edge_equi,
                   pn_vtx_equi,
                   pequi_face_edge_idx,
                   pequi_face_edge,
                   pequi_edge_vtx,
                   block_entity_equi_child_g_num,
                   pequi_edge_ln_to_gn,
                   pequi_vtx_ln_to_gn,
                   pequi_vtx_coord,
                   pequi_field,
                   pequi_gradient_field,
                   &isoline_n_vtx,
                   &isoline_vtx_coord,
                   &isoline_vtx_parent_face,
                   &isoline_vtx_ln_to_gn,
                   &isoline_n_edge,
                   &isoline_edge_vtx_idx,
                   &isoline_edge_vtx,
                   &isoline_edge_ln_to_gn);

    if (1) {
      _dump_iso_line_dist(isos->comm,
                          isoline_n_vtx,
                          isoline_vtx_coord,
                          isoline_vtx_parent_face,
                          isoline_vtx_ln_to_gn,
                          isoline_n_edge,
                          isoline_edge_vtx_idx,
                          isoline_edge_vtx,
                          isoline_edge_ln_to_gn);
    }
  } else {
    _iso_surf_dist(isos,
                   n_entity_equi,
                   pn_face_equi,
                   pn_edge_equi,
                   pn_vtx_equi,
                   pequi_cell_face_idx,
                   pequi_cell_face,
                   pequi_face_edge_idx,
                   pequi_face_edge,
                   pequi_edge_vtx,
                   block_entity_equi_child_g_num,
                   pequi_face_ln_to_gn,
                   pequi_edge_ln_to_gn,
                   pequi_vtx_ln_to_gn,
                   pequi_vtx_coord,
                   pequi_field,
                   pequi_gradient_field,
                   &isoline_n_vtx,
                   &isoline_vtx_coord,
                   &isoline_vtx_parent_face,
                   &isoline_vtx_ln_to_gn,
                   &isoline_n_edge,
                   &isoline_edge_vtx_idx,
                   &isoline_edge_vtx,
                   &isoline_edge_ln_to_gn);
  }

  if(isos->dim == 3) {
    free(pequi_face_ln_to_gn);
  }

  free(isoline_vtx_coord      );
  free(isoline_vtx_parent_face);
  free(isoline_vtx_ln_to_gn   );
  free(isoline_edge_vtx_idx   );
  free(isoline_edge_vtx       );
  free(isoline_edge_ln_to_gn  );

  if(isos->iso_kind == PDM_ISO_SURFACE_KIND_FIELD) {
    free(pequi_field);
    free(pequi_gradient_field);
  }


  free(pequi_edge_ln_to_gn);
  free(pequi_vtx_ln_to_gn);
  free(pequi_vtx_coord);
  free(pequi_parent_vtx_ln_to_gn);
  free(pequi_edge_vtx_idx);
  free(pequi_edge_vtx);
  free(pequi_parent_edge_ln_to_gn);
  free(pequi_face_edge_idx);
  free(pequi_face_edge);
  free(block_entity_equi_parent_g_num);
  PDM_part_to_block_free(ptb);

  free(child_equi_entity_gnum);

  if(isos->dim == 3) {
    free(dentity_edge);
    free(dentity_edge_idx);
  }

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_iso_surface_t*
PDM_iso_surface_create
(
 const int                    dim,
       PDM_iso_surface_kind_t iso_kind,
 const int                    n_part,
       PDM_ownership_t        ownership,
       PDM_MPI_Comm           comm
)
{
  PDM_iso_surface_t *isos = (PDM_iso_surface_t *) malloc(sizeof(PDM_iso_surface_t));

  isos->dim       = dim;
  isos->iso_kind  = iso_kind;
  isos->n_part    = n_part;
  isos->ownership = ownership;
  isos->comm      = comm;
  isos->is_dist   = -1;

  isos->n_cell         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_face         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_edge         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_vtx          = (int          *) malloc(n_part * sizeof(int          ));

  isos->pcell_face     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pcell_face_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pedge_vtx      = (int         **) malloc(n_part * sizeof(int         *));
  isos->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));

  isos->pface_vtx_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_vtx     = (int         **) malloc(n_part * sizeof(int         *));

  isos->pvtx_coord      = (double **) malloc(n_part * sizeof(double *));
  isos->pfield          = (double **) malloc(n_part * sizeof(double *));
  isos->pgradient_field = (double **) malloc(n_part * sizeof(double *));

  return isos;
}

void
PDM_iso_surface_compute
(
  PDM_iso_surface_t        *isos
)
{
  if(isos->is_dist == 0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_iso_surface_compute Not implemented with partition layout\n");
  } else {
    _iso_surface_dist(isos);
  }



}

// See with Eric et Bastien : par type ou une fonction avec 1000 arguments ?
void
PDM_iso_surface_part_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  isos->n_cell        [i_part] = n_cell;
  isos->n_face        [i_part] = n_face;
  isos->n_edge        [i_part] = n_edge;
  isos->n_vtx         [i_part] = n_vtx;
  isos->pcell_face    [i_part] = cell_face;
  isos->pcell_face_idx[i_part] = cell_face_idx;
  isos->pface_edge    [i_part] = face_edge;
  isos->pface_edge_idx[i_part] = face_edge_idx;
  isos->pedge_vtx     [i_part] = edge_vtx;
  isos->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  isos->face_ln_to_gn [i_part] = face_ln_to_gn;
  isos->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  isos->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  isos->pface_vtx_idx [i_part] = face_vtx_idx;
  isos->pface_vtx     [i_part] = face_vtx;
  isos->pvtx_coord    [i_part] = vtx_coord;
}


void
PDM_iso_surface_part_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *field
)
{
  isos->is_dist = 0;
  isos->pfield[i_part] = field;
}

void
PDM_iso_surface_part_gradient_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *gradient_field
)
{
  isos->is_dist = 0;
  isos->pgradient_field[i_part] = gradient_field;
}

void
PDM_iso_surface_dconnectivity_set
(
  PDM_iso_surface_t        *isos,
  PDM_connectivity_type_t   connectivity_type,
  PDM_g_num_t              *dconnect,
  int                      *dconnect_idx
)
{
  isos->is_dist = 1;
  switch (connectivity_type) {
   case PDM_CONNECTIVITY_TYPE_CELL_FACE:
     isos->dcell_face     = dconnect;
     isos->dcell_face_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
     isos->dface_edge     = dconnect;
     isos->dface_edge_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_VTX:
     isos->dface_vtx     = dconnect;
     isos->dface_vtx_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
     isos->dedge_vtx     = dconnect;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type for iso_surface %d\n", connectivity_type);
    break;
   }

}

void
PDM_iso_surface_distrib_set
(
  PDM_iso_surface_t        *isos,
  PDM_mesh_entities_t       entity_type,
  PDM_g_num_t              *distrib_entity
)
{
  isos->is_dist = 1;
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     isos->distrib_cell = distrib_entity;
     break;
   case PDM_MESH_ENTITY_FACE:
     isos->distrib_face = distrib_entity;
     break;
   case PDM_MESH_ENTITY_EDGE:
     isos->distrib_edge = distrib_entity;
     break;
   case PDM_MESH_ENTITY_VERTEX:
     isos->distrib_vtx = distrib_entity;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid entity_type for iso_surface %d\n", entity_type);
    break;
   }
}


void
PDM_iso_surface_dvtx_coord_set
(
  PDM_iso_surface_t *isos,
  double            *dvtx_coord
)
{
  isos->is_dist     = 1;
  isos->dvtx_coord  = dvtx_coord;
}

void
PDM_iso_surface_dfield_set
(
  PDM_iso_surface_t *isos,
  double            *dfield
)
{
  isos->is_dist = 1;
  isos->dfield  = dfield;
}

void
PDM_iso_surface_dgrad_field_set
(
  PDM_iso_surface_t *isos,
  double            *dgrad_field
)
{
  isos->is_dist         = 1;
  isos->dgradient_field = dgrad_field;
}


void
PDM_iso_surface_plane_equation_set
(
  PDM_iso_surface_t        *isos,
  double                    a,
  double                    b,
  double                    c,
  double                    d
)
{
  isos->plane_equation[0] = a;
  isos->plane_equation[1] = b;
  isos->plane_equation[2] = c;
  isos->plane_equation[3] = d;
}


void
PDM_iso_surface_free
(
  PDM_iso_surface_t        *isos
)
{
  free(isos->n_cell        );
  free(isos->n_face        );
  free(isos->n_edge        );
  free(isos->n_vtx         );
  free(isos->pcell_face    );
  free(isos->pcell_face_idx);
  // Si pface_edge a été calculé il faut le free
  free(isos->pface_edge    );
  free(isos->pface_edge_idx);

  free(isos->pface_vtx     );
  free(isos->pface_vtx_idx );
  free(isos->pedge_vtx     );
  free(isos->cell_ln_to_gn );
  free(isos->face_ln_to_gn );
  free(isos->edge_ln_to_gn );
  free(isos->vtx_ln_to_gn  );

  free(isos->pvtx_coord     );
  free(isos->pfield         );
  free(isos->pgradient_field);

  free(isos);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
