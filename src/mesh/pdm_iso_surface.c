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
#include "pdm_extract_part.h"
#include "pdm_part_to_part.h"
#include "pdm_gnum.h"
#include "pdm_unique.h"
#include "pdm_writer.h"
#include "pdm_multipart.h"

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
  const int dbg = 1;

  if (dbg) {
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

  if (dbg) {
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
  }

  const double tol = 1.e-2;
  double tol_s = tol * PDM_ABS(S[0]);
  // log_trace("cond = %e\n", S[0]/S[n-1]);
  // log_trace("S truncated = [");
  for (int i = 0; i < n; i++) {
    double si = S[i];
    if (PDM_ABS(si) <= tol_s) {
      si = 0;
    }
    // log_trace("%20.16f ", si);
  }
  // log_trace("]\n");

  /* Compute y := S^{-1} U^T b */
  double y[n];

  for (int i = 0; i < n; i++) {

    y[i] = 0.;

    if (PDM_ABS(S[i]) > tol_s) {
      for (int j = 0; j < m; j++) {
        y[i] += U[j + m*i] * b[j];
      }

      y[i] /= S[i];
    } else {
      S[i] = 0.;
    }
  }

  if (dbg) {
    log_trace("y = [");
    for (int i = 0; i < n; i++) {
      log_trace("%f ", y[i]);
    }
    log_trace("]\n");
  }

  /* Compute x := V^T y */
  for (int i = 0; i < n; i++) {
    x[i] = 0.;

    for (int j = 0; j < m; j++) {
      x[i] += V[j + m*i] * y[j];
    }
  }

  if (dbg) {
    log_trace("x = [");
    for (int i = 0; i < n; i++) {
      log_trace("%f ", x[i]);
    }
    log_trace("]\n");
  }
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
        // log_trace("coefs = %f %f %f %f\n", a_3, a_2, a_1, a_0);

        double s = t;
        int stat = 0;
        // log_trace("Newton:\n");
        for (int iter = 0; iter < 5; iter++) {
          double val = a_0 + s*(a_1 + s*(a_2 + s*a_3));
          // log_trace("  it %d, s = %f, |val| = %e\n", iter, s, PDM_ABS(val));

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
        double gx, gy, gz;
        if (isos->eval_field_and_gradient != NULL) {
          // Exact gradient
          double f;
          isos->eval_field_and_gradient(edge_coord[3*i  ],
                                        edge_coord[3*i+1],
                                        edge_coord[3*i+2],
                                        &f,
                                        &gx, &gy, &gz);
        } else {
          gx = (1. - t)*grad1[0] + t*grad2[0];
          gy = (1. - t)*grad1[1] + t*grad2[1];
          gz = (1. - t)*grad1[2] + t*grad2[2];
        }

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
PDM_polygon_status_t
_check_point_in_polygon
(
 const double  point[],
 const int     iface,
 const int     pface_edge_idx[],
 const int     pface_edge[],
 const int     pedge_vtx[],
 const double  pvtx_coord[],
 double       *poly_coord
 )
{
  int n_edge = pface_edge_idx[iface+1] - pface_edge_idx[iface];
  // Get poly_coord
  int cur_vtx, next_vtx;
  int cur_edge = pface_edge[pface_edge_idx[iface]];
  if (cur_edge < 0) {
    cur_edge = -cur_edge - 1;
    cur_vtx  = pedge_vtx[2*cur_edge+1];
    next_vtx = pedge_vtx[2*cur_edge  ];
  } else {
    cur_edge = cur_edge - 1;
    cur_vtx  = pedge_vtx[2*cur_edge  ];
    next_vtx = pedge_vtx[2*cur_edge+1];
  }

  for (int ivtx = 0; ivtx < n_edge; ivtx++) {
    memcpy(poly_coord + 3*ivtx, pvtx_coord + 3*(cur_vtx - 1), sizeof(double) * 3);

    for (int iedg = pface_edge_idx[iface]; iedg < pface_edge_idx[iface+1]; iedg++) {
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
  for (int j = 0; j < n_edge; j++) {
    for (int l = 0; l < 3; l++) {
      double x = poly_coord[3*j + l];
      poly_bounds[2*l  ] = PDM_MIN(poly_bounds[2*l  ], x);
      poly_bounds[2*l+1] = PDM_MAX(poly_bounds[2*l+1], x);
    }
  }

  PDM_plane_normal(n_edge,
                   poly_coord,
                   poly_normal);

  PDM_polygon_status_t stat = PDM_polygon_point_in_new(point,
                                                       n_edge,
                                                       poly_coord,
                                                       poly_bounds,
                                                       poly_normal);

  return stat;
}

static
PDM_polygon_status_t
_check_point_in_polygon2
(
 const double  point[],
 const int     iface,
 const int     pface_vtx_idx[],
 const int     pface_vtx[],
 const double  pvtx_coord[],
 double       *poly_coord
 )
 {
  int n_vtx = pface_vtx_idx[iface+1] - pface_vtx_idx[iface];
  const int *fv = pface_vtx + pface_vtx_idx[iface];

  for (int i = 0; i < n_vtx; i++) {
    int ivtx = fv[i] - 1;

    memcpy(poly_coord + 3*i,
           pvtx_coord + 3*ivtx,
           sizeof(double) * 3);
  }


  double poly_bounds[6] = {
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL
  };
  double poly_normal[3];
  for (int j = 0; j < n_vtx; j++) {
    for (int l = 0; l < 3; l++) {
      double x = poly_coord[3*j + l];
      poly_bounds[2*l  ] = PDM_MIN(poly_bounds[2*l  ], x);
      poly_bounds[2*l+1] = PDM_MAX(poly_bounds[2*l+1], x);
    }
  }

  PDM_plane_normal(n_vtx,
                   poly_coord,
                   poly_normal);

  PDM_polygon_status_t stat = PDM_polygon_point_in_new(point,
                                                       n_vtx,
                                                       poly_coord,
                                                       poly_bounds,
                                                       poly_normal);

  return stat;
}


// static void
// _compute_face_vtx
// (
//  const int   n_face,
//  int        *pface_edge_idx,
//  int        *pface_edge,
//  int        *pedge_vtx,
//  int       **pface_vtx
//  )
// {
//   int dbg = 0;
//   *pface_vtx = (int *) malloc(sizeof(int) * pface_edge_idx[n_face]);

//   for (int i = 0; i < n_face; i++) {

//     if (dbg) {
//       log_trace("\nFace %d\n", i);
//       for (int idx_edge = pface_edge_idx[i]; idx_edge < pface_edge_idx[i+1]; idx_edge++) {
//         int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;
//         log_trace("  edge %d: %d %d\n",
//                   pface_edge[idx_edge],
//                   pedge_vtx[2*iedge], pedge_vtx[2*iedge+1]);
//       }
//     }
//     int *_pface_vtx = *pface_vtx + pface_edge_idx[i];

//     int cur_vtx, next_vtx;
//     int cur_edge = pface_edge[pface_edge_idx[i]];
//     if (cur_edge < 0) {
//       cur_edge = -cur_edge - 1;
//       cur_vtx  = pedge_vtx[2*cur_edge+1];
//       next_vtx = pedge_vtx[2*cur_edge  ];
//     } else {
//       cur_edge = cur_edge - 1;
//       cur_vtx  = pedge_vtx[2*cur_edge  ];
//       next_vtx = pedge_vtx[2*cur_edge+1];
//     }

//     for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
//       _pface_vtx[ivtx] = cur_vtx;

//       for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
//         cur_edge = pface_edge[iedg];
//         int vtx1, vtx2;
//         if (cur_edge < 0) {
//           cur_edge = -cur_edge - 1;
//           vtx1 = pedge_vtx[2*cur_edge+1];
//           vtx2 = pedge_vtx[2*cur_edge  ];
//         } else {
//           cur_edge = cur_edge - 1;
//           vtx1 = pedge_vtx[2*cur_edge  ];
//           vtx2 = pedge_vtx[2*cur_edge+1];
//         }

//         if (vtx1 == next_vtx) {
//           cur_vtx  = next_vtx;
//           next_vtx = vtx2;
//           break;
//         }
//       }
//     }

//     if (dbg) {
//       log_trace("  face_vtx = ");
//       for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
//         log_trace("%d ", _pface_vtx[ivtx]);
//       }
//       log_trace("\n");
//     }

//   }
// }

static void
_compute_face_vtx
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  int dbg = 0;

  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *edge_tag = PDM_array_zeros_int(n_edge);

  for (int iface = 0; iface < n_face; iface++) {
    int *_face_vtx  = *face_vtx  + face_edge_idx[iface];
    int *_face_edge =  face_edge + face_edge_idx[iface];

    if (dbg) {
      log_trace("\nFace %d\n", iface);
      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }

    int _n_edge = face_edge_idx[iface+1] - face_edge_idx[iface];
    // first edge
    int iedge = PDM_ABS(_face_edge[0]) - 1;
    edge_tag[iedge] = 1;
    _face_vtx[0] = edge_vtx[2*iedge  ];
    _face_vtx[1] = edge_vtx[2*iedge+1];

    for (int i = 2; i < _n_edge; i++) {

      for (int j = 1; j < _n_edge; j++) {
        iedge = PDM_ABS(_face_edge[j]) - 1;

        if (edge_tag[iedge]) {
          continue;
        }

        if (edge_vtx[2*iedge] == _face_vtx[i-1]) {
          _face_vtx[i] = edge_vtx[2*iedge+1];
          edge_tag[iedge] = 1;
          break;
        }
        else if (edge_vtx[2*iedge+1] == _face_vtx[i-1]) {
          _face_vtx[i] = edge_vtx[2*iedge];
          edge_tag[iedge] = 1;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[iface+1] - face_edge_idx[iface]; ivtx++) {
        log_trace("%d ", _face_vtx[ivtx]);
      }
      log_trace("\n");
    }

    // reset tags
    for (int i = 0; i < _n_edge; i++) {
      iedge = PDM_ABS(_face_edge[i]) - 1;
      edge_tag[iedge] = 0;
    }
  }
  free(edge_tag);
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
 double             *pgradient_field
)
{
  const int use_gradient = (pgradient_field != NULL);

  PDM_MPI_Comm comm = isos->comm;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  char filename[999];
  if (isos->debug) {
    sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_vtx,
                              pvtx_coord,
                              pvtx_ln_to_gn,
                              NULL);
  }

  int *pface_vtx = NULL;
  _compute_face_vtx(n_face,
                    pface_edge_idx,
                    pface_edge,
                    pedge_vtx,
                    &pface_vtx);
  if (1) {
    // PDM_log_trace_connectivity_int(pface_edge_idx,
    //                                pface_edge,
    //                                n_face,
    //                                "face_edge : ");

    // PDM_log_trace_connectivity_int(pface_edge_idx,
    //                                pface_vtx,
    //                                n_face,
    //                                "pface_vtx : ");

    sprintf(filename, "extract_face_%2.2d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           pvtx_coord,
                           pvtx_ln_to_gn,
                           n_face,
                           pface_edge_idx,
                           pface_vtx,
                           pface_ln_to_gn,
                           NULL);
  }

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

  if (1) {//isos->debug) {
    sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
    _dump_vectors (filename,
                   n_edge,
                   edge_coord,
                   edge_gradient,
                   NULL);
  }

  int n_face_edge_max = 0;
  for (int i = 0; i < n_face; i++) {
    n_face_edge_max = PDM_MAX(n_face_edge_max,
                              pface_edge_idx[i+1] - pface_edge_idx[i]);
  }


  const int dim = 2;
  double *mat = (double *) malloc (sizeof(double) * n_face_edge_max * dim);
  double *rhs = (double *) malloc (sizeof(double) * n_face_edge_max);

  double *S = (double *) malloc(sizeof(double) * dim);
  double *U = (double *) malloc(sizeof(double) * n_face_edge_max * dim);
  double *V = (double *) malloc(sizeof(double) * n_face_edge_max * dim);

  double *face_coord = (double *) malloc(sizeof(double) * n_face * 3);
  double *poly_coord = (double *) malloc(sizeof(double) * n_face_edge_max * 3);

  for (int iface = 0; iface < n_face; iface++) {

    int n_face_edge = 0;
    double _face_coord[3] = {0.};

    // First loop to count number of tagged edges
    for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {

      int iedge = PDM_ABS(pface_edge[jedge]) - 1;

      if (edge_tag[iedge] != 0) {
        n_face_edge++;
      }
    } // end of loop on edges of current face


    // Second loop to fill the matrix
    int i_face_edge = 0;
    for (int jedge = pface_edge_idx[iface]; jedge < pface_edge_idx[iface+1]; jedge++) {

      int iedge = PDM_ABS(pface_edge[jedge]) - 1;

      if (edge_tag[iedge] != 0) {

        for (int l = 0; l < dim; l++) {
          _face_coord[l] += edge_coord[3*iedge + l];
        }

        if (use_gradient) {
          for (int l = 0; l < dim; l++) {
            mat[i_face_edge + l*n_face_edge] = edge_gradient[3*iedge + l];
          }
          rhs[i_face_edge] = PDM_DOT_PRODUCT(edge_gradient + 3*iedge, edge_coord + 3*iedge);
          i_face_edge++;
        }

      }
    } // end of loop on edges of current face

    assert (n_face_edge > 1);

    double normalization = 1. / (double) n_face_edge;
    for (int l = 0; l < dim; l++) {
      _face_coord[l] *= normalization;
    }

    if (use_gradient) {
      if (1) {//isos->debug) {
        log_trace("\nFace "PDM_FMT_G_NUM":\n", pface_ln_to_gn[iface]);
      }
      _solve_least_square (n_face_edge,
                           dim,
                           mat,
                           U,
                           S,
                           V,
                           rhs,
                           face_coord + 3*iface);
      face_coord[3*iface+2] = 0.;

      // PDM_polygon_status_t stat = _check_point_in_polygon(face_coord + 3*iface,
      //                                                     iface,
      //                                                     pface_edge_idx,
      //                                                     pface_edge,
      //                                                     pedge_vtx,
      //                                                     pvtx_coord,
      //                                                     poly_coord);
      PDM_polygon_status_t stat = _check_point_in_polygon2(face_coord + 3*iface,
                                                           iface,
                                                           pface_edge_idx,
                                                           pface_vtx,
                                                           pvtx_coord,
                                                           poly_coord);

      if (1) {//isos->debug) {
        log_trace("stat = %d\n", (int) stat);
      }
      if (stat != PDM_POLYGON_INSIDE) {
        memcpy(face_coord + 3*iface, _face_coord, sizeof(double) * 3);
      }
    }

    else {
      memcpy(face_coord + 3*iface, _face_coord, sizeof(double) * 3);
    }

  } // end of loop on faces
  free(mat);
  free(rhs);
  free(S);
  free(U);
  free(V);
  free(edge_gradient);

  free(pface_vtx);


  /*
   *  Build connectivity and global ids
   */
  int n_tagged_edge = 0;
  int *i_tagged_edge = PDM_array_const_int(n_edge, -1);

  PDM_g_num_t *tagged_edge_parent_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_edge);

  for (int i = 0; i < n_edge; i++) {
    if (edge_tag[i] != 0) {
      tagged_edge_parent_ln_to_gn[n_tagged_edge] = pedge_ln_to_gn[i];
      i_tagged_edge[i] = n_tagged_edge++;
    }
  }
  tagged_edge_parent_ln_to_gn = realloc(tagged_edge_parent_ln_to_gn,
                                        sizeof(PDM_g_num_t) * n_tagged_edge);


  PDM_gen_gnum_t *gen_gnum_edge = PDM_gnum_create(3,
                                                  1,
                                                  PDM_FALSE,
                                                  1.,
                                                  isos->comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gen_gnum_edge,
                            0,
                            n_tagged_edge,
                            tagged_edge_parent_ln_to_gn);

  PDM_gnum_compute(gen_gnum_edge);

  PDM_g_num_t *tagged_edge_ln_to_gn = PDM_gnum_get(gen_gnum_edge, 0);


  PDM_g_num_t gn_face;
  if (i_rank == n_rank-1) {
    gn_face = pface_ln_to_gn[n_face-1];
  }
  PDM_MPI_Bcast(&gn_face, 1, PDM__PDM_MPI_G_NUM, n_rank-1, comm);


  isos->isosurf_n_vtx = n_face + n_tagged_edge;
  isos->isosurf_vtx_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * isos->isosurf_n_vtx);
  memcpy(isos->isosurf_vtx_ln_to_gn, pface_ln_to_gn, sizeof(PDM_g_num_t) * n_face);

  isos->isosurf_vtx_coord = (double *) malloc(sizeof(double) * isos->isosurf_n_vtx * 3);
  memcpy(isos->isosurf_vtx_coord, face_coord, sizeof(double) * n_face * 3);
  free(face_coord);

  int i_vtx = n_face;
  for (int i = 0; i < n_edge; i++) {
    if (edge_tag[i] != 0) {
      isos->isosurf_vtx_ln_to_gn[i_vtx] = gn_face + tagged_edge_ln_to_gn[i_tagged_edge[i]];
      memcpy(isos->isosurf_vtx_coord + 3*i_vtx,
             edge_coord + 3*i,
             sizeof(double) * 3);
      i_vtx++;
    }
  }
  free(edge_coord);


  isos->isosurf_n_edge = 0;
  for (int iface = 0; iface < n_face; iface++) {

    for (int idx_edge = pface_edge_idx[iface]; idx_edge < pface_edge_idx[iface+1]; idx_edge++) {
      int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;

      if (edge_tag[iedge] != 0) {
        isos->isosurf_n_edge++;
      }
    }

  }

  isos->isosurf_edge_vtx = (int *) malloc(sizeof(int) * isos->isosurf_n_edge * 2);

  PDM_g_num_t *distrib_isosurf_edge = PDM_compute_entity_distribution(comm,
                                                                      isos->isosurf_n_edge);

  isos->isosurf_edge_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * isos->isosurf_n_edge);

  isos->isosurf_n_edge = 0;
  for (int iface = 0; iface < n_face; iface++) {

    for (int idx_edge = pface_edge_idx[iface]; idx_edge < pface_edge_idx[iface+1]; idx_edge++) {
      int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;

      if (edge_tag[iedge] == 0) {
        continue;
      }

      int sgn_edge = PDM_SIGN(pface_edge[idx_edge]);
      int ipt_edge = n_face + i_tagged_edge[iedge] + 1;

      int sgn = sgn_edge * edge_tag[iedge];

      if (sgn < 0) {
        isos->isosurf_edge_vtx[2*isos->isosurf_n_edge  ] = iface + 1;
        isos->isosurf_edge_vtx[2*isos->isosurf_n_edge+1] = ipt_edge;
      } else {
        isos->isosurf_edge_vtx[2*isos->isosurf_n_edge  ] = ipt_edge;
        isos->isosurf_edge_vtx[2*isos->isosurf_n_edge+1] = iface + 1;
      }

      isos->isosurf_edge_ln_to_gn[isos->isosurf_n_edge] = distrib_isosurf_edge[i_rank] + isos->isosurf_n_edge + 1;
      isos->isosurf_n_edge++;

    }

  }
  free(edge_tag);
  free(distrib_isosurf_edge);
  free(i_tagged_edge);

  free(tagged_edge_parent_ln_to_gn);
  PDM_gnum_free(gen_gnum_edge);
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
  double             *pgradient_field
)
{
  const int use_gradient = (pgradient_field != NULL);

  PDM_MPI_Comm comm = isos->comm;

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  char filename[999];
  if (isos->debug) {
    sprintf(filename, "out_equi_vtx_coord_%2.2d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_vtx,
                              pvtx_coord,
                              pvtx_ln_to_gn,
                              NULL);
  }

  /*
   *  Tag edges that cross the iso-surface,
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

  if (1) {//isos->debug) {
    sprintf(filename, "edge_intersection_%2.2d.vtk", i_rank);
    _dump_vectors (filename,
                   n_edge,
                   edge_coord,
                   edge_gradient,
                   NULL);
  }

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


  int *is_treated_face = PDM_array_zeros_int(n_face);
  int *face_tag = PDM_array_zeros_int(n_face);

  const int dim = 3;
  double *mat_cell = (double *) malloc (sizeof(double) * n_cell_edge_max * dim);
  double *rhs_cell = (double *) malloc (sizeof(double) * n_cell_edge_max);
  double *mat_face = (double *) malloc (sizeof(double) * (n_face_edge_max+1) * dim);
  double *rhs_face = (double *) malloc (sizeof(double) * (n_face_edge_max+1));

  double *S_cell = (double *) malloc(sizeof(double) * dim);
  double *U_cell = (double *) malloc(sizeof(double) * n_cell_edge_max * dim);
  double *V_cell = (double *) malloc(sizeof(double) * n_cell_edge_max * dim);
  double *S_face = (double *) malloc(sizeof(double) * dim);
  double *U_face = (double *) malloc(sizeof(double) * (n_face_edge_max+1) * dim);
  double *V_face = (double *) malloc(sizeof(double) * (n_face_edge_max+1) * dim);

  int *is_visited_edge = PDM_array_zeros_int(n_edge);
  int *visited_edges = (int *) malloc(sizeof(int) * n_edge);

  double *face_coord = (double *) malloc(sizeof(double) * n_face * 3);
  double *cell_coord = (double *) malloc(sizeof(double) * n_cell * 3);

  double normalization;
  double face_center[3], face_normal[3];
  double *poly_coord = (double *) malloc(sizeof(double) * n_face_edge_max * 3);
  double *all_face_center = NULL;
  double *all_face_normal = NULL;
  if (isos->debug) {
    all_face_center = (double *) malloc(sizeof(double) * n_face * 3);
    all_face_normal = (double *) malloc(sizeof(double) * n_face * 3);
  }

  // tmp -->>
  for (int i = 0; i < 3*n_cell; i++) {
    cell_coord[i] = 0.;
  }
  for (int i = 0; i < 3*n_face; i++) {
    face_coord[i] = 0.;
  }

  if (isos->debug) {
    for (int i = 0; i < 3*n_face; i++) {
      all_face_center[i] = 0.;
      all_face_normal[i] = 0.;
    }
  }
  //<<--

  for (int icell = 0; icell < n_cell; icell++) {

    int n_cell_edge = 0;
    double _cell_coord[3] = {0.};

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


    // log_trace("_n_cell_edge = %d\n", n_cell_edge);
    n_cell_edge /= 2;

    // Second loop to fill matrices
    int i_cell_edge = 0;
    for (int jface = pcell_face_idx[icell]; jface < pcell_face_idx[icell+1]; jface++) {

      int iface = PDM_ABS(pcell_face[jface]) - 1;
      int n_face_edge = 0;
      double _face_coord[3] = {0.};

      if (is_treated_face[iface] == 0) {
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

          double _rhs = 0.;
          if (use_gradient) {
            _rhs = PDM_DOT_PRODUCT(edge_gradient + 3*iedge, edge_coord + 3*iedge);
          }

          if (is_visited_edge[iedge] == 0) {
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
            is_visited_edge[iedge] = 1;
            i_cell_edge++;
          }

          if (is_treated_face[iface] == 0) {
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

      if (is_treated_face[iface] == 0) {

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

        normalization = 1. / (double) (pface_edge_idx[iface+1] - pface_edge_idx[iface]);
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

        if (isos->debug) {
          memcpy(all_face_center + 3*iface, face_center, sizeof(double)*3);
          memcpy(all_face_normal + 3*iface, face_normal, sizeof(double)*3);
        }

        for (int l = 0; l < dim; l++) {
          mat_face[i_face_edge + l*n_face_edge] = face_normal[l];
        }
        rhs_face[i_face_edge] = PDM_DOT_PRODUCT(face_center, face_normal);
        i_face_edge++;

        // log_trace("face %d ("PDM_FMT_G_NUM"), i_face_edge = %d\n",
                  // iface, pface_ln_to_gn[iface], i_face_edge);
        if (0) {
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
          normalization = 1. / (double) (i_face_edge - 1);
          for (int l = 0; l < dim; l++) {
            _face_coord[l] *= normalization;
          }

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
              memcpy(face_coord + 3*iface, _face_coord, sizeof(double) * 3);
            }

            else {
              PDM_polygon_status_t stat = _check_point_in_polygon(face_coord + 3*iface,
                                                                  iface,
                                                                  pface_edge_idx,
                                                                  pface_edge,
                                                                  pedge_vtx,
                                                                  pvtx_coord,
                                                                  poly_coord);
              if (stat != PDM_POLYGON_INSIDE) {
                memcpy(face_coord + 3*iface, _face_coord, sizeof(double) * 3);
              }
            }

          } else {
            memcpy(face_coord + 3*iface, _face_coord, sizeof(double) * 3);
          }

          face_tag[iface] = 1;

          // log_trace("face_coord = %f %f %f\n",
          //     face_coord[3*iface], face_coord[3*iface+1], face_coord[3*iface+2]);
        }
        else {
          face_tag[iface] = 0;
        }

        is_treated_face[iface] = 1;
      }



    } // end of loop on faces of current cell

    normalization = 1. / (double) i_cell_edge;
    for (int l = 0; l < dim; l++) {
      _cell_coord[l] *= normalization;
    }

    // log_trace("cell %d ("PDM_FMT_G_NUM"), i_cell_edge = %d\n",
    //           icell, pcell_ln_to_gn[icell], i_cell_edge);
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

      int mat_rank;
      for (mat_rank = 3; mat_rank >= 0; mat_rank--) {
        if (PDM_ABS(S_cell[mat_rank-1]) > 1e-9) {
          break;
        }
      }

      // Beware of NaNs
      for (int l = 0; l < dim; l++) {
        if (!(cell_coord[3*icell+l] >= 0) &&
            !(cell_coord[3*icell+l] < 0)) {
          mat_rank = -1;
        }
      }

      if (mat_rank == 2) {
        // Ridge
        double origin[3] = {cell_coord[3*icell], cell_coord[3*icell+1], cell_coord[3*icell+2]};
        double direction[3] = {V_cell[2], V_cell[i_cell_edge+2], V_cell[2*i_cell_edge+2]};

        // we have a line passing through 'origin' and directed by 'direction'
        // find the orthogonal projection of '_cell_coord' on that line
        double t =
        direction[0]*(_cell_coord[0] - origin[0]) +
        direction[1]*(_cell_coord[1] - origin[1]) +
        direction[2]*(_cell_coord[2] - origin[2]);

        cell_coord[3*icell  ] = origin[0] + t*direction[0];
        cell_coord[3*icell+1] = origin[1] + t*direction[1];
        cell_coord[3*icell+2] = origin[2] + t*direction[2];
      }

      else if (mat_rank <= 1) {
        // Plane
        memcpy(cell_coord + 3*icell, _cell_coord, sizeof(double) * 3);
      }

      /* Compute location in cell */
      //...
      memcpy(cell_coord + 3*icell, _cell_coord, sizeof(double) * 3);//!!!

    } else {
      memcpy(cell_coord + 3*icell, _cell_coord, sizeof(double) * 3);
    }
    // log_trace("cell_coord = %f %f %f\n",
    //           cell_coord[3*icell], cell_coord[3*icell+1], cell_coord[3*icell+2]);

    for (int i = 0; i < i_cell_edge; i++) {
      is_visited_edge[visited_edges[i]] = 0;
    }

  } // end of loop on cells






  // Visu
  if (isos->debug) {
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

  if (isos->debug) {
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
  }

  free (is_visited_edge);
  free (visited_edges);

  free (mat_cell);
  free (rhs_cell);
  free (mat_face);
  free (rhs_face);

  free (S_cell);
  free (U_cell);
  free (V_cell);
  free (S_face);
  free (U_face);
  free (V_face);

  free(poly_coord);


  /*
   *  Build connectivity and global ids
   */
  int n_tagged_edge = 0;
  int n_tagged_face = 0;
  int *i_tagged_edge = PDM_array_const_int(n_edge, -1);
  int *i_tagged_face = PDM_array_const_int(n_face, -1);

  PDM_g_num_t *tagged_edge_parent_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_edge);

  for (int i = 0; i < n_edge; i++) {
    if (edge_tag[i] != 0) {
      tagged_edge_parent_ln_to_gn[n_tagged_edge] = pedge_ln_to_gn[i];
      i_tagged_edge[i] = n_tagged_edge++;
    }
  }
  tagged_edge_parent_ln_to_gn = realloc(tagged_edge_parent_ln_to_gn,
                                        sizeof(PDM_g_num_t) * n_tagged_edge);


  PDM_g_num_t *tagged_face_parent_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face);

  for (int i = 0; i < n_face; i++) {
    if (face_tag[i] != 0) {
      tagged_face_parent_ln_to_gn[n_tagged_face] = pface_ln_to_gn[i];
      i_tagged_face[i] = n_tagged_face++;
    }
  }
  tagged_face_parent_ln_to_gn = realloc(tagged_face_parent_ln_to_gn,
                                        sizeof(PDM_g_num_t) * n_tagged_face);


  PDM_gen_gnum_t *gen_gnum_face = PDM_gnum_create(3,
                                                  1,
                                                  PDM_FALSE,
                                                  1.,
                                                  isos->comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gen_gnum_face,
                            0,
                            n_tagged_face,
                            tagged_face_parent_ln_to_gn);

  PDM_gnum_compute(gen_gnum_face);

  PDM_g_num_t *tagged_face_ln_to_gn = PDM_gnum_get(gen_gnum_face, 0);


  PDM_gen_gnum_t *gen_gnum_edge = PDM_gnum_create(3,
                                                  1,
                                                  PDM_FALSE,
                                                  1.,
                                                  isos->comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gen_gnum_edge,
                            0,
                            n_tagged_edge,
                            tagged_edge_parent_ln_to_gn);

  PDM_gnum_compute(gen_gnum_edge);

  PDM_g_num_t *tagged_edge_ln_to_gn = PDM_gnum_get(gen_gnum_edge, 0);

  PDM_g_num_t lmax_cell_ln_to_gn = 0;
  for (int i = 0; i < n_cell; i++) {
    lmax_cell_ln_to_gn = PDM_MAX(lmax_cell_ln_to_gn, pcell_ln_to_gn[i]);
  }

  PDM_g_num_t gn_cell;
  PDM_MPI_Allreduce(&lmax_cell_ln_to_gn, &gn_cell, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, isos->comm);

  PDM_g_num_t lmax_tagged_face_ln_to_gn = 0;
  for (int i = 0; i < n_tagged_face; i++) {
    lmax_tagged_face_ln_to_gn = PDM_MAX(lmax_tagged_face_ln_to_gn, tagged_face_ln_to_gn[i]);
  }

  PDM_g_num_t gmax_tagged_face_ln_to_gn;
  PDM_MPI_Allreduce(&lmax_tagged_face_ln_to_gn, &gmax_tagged_face_ln_to_gn, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, isos->comm);


  isos->isosurf_n_vtx = n_cell + n_tagged_face + n_tagged_edge;
  isos->isosurf_vtx_coord = (double *) malloc(sizeof(double) * isos->isosurf_n_vtx * 3);
  memcpy(isos->isosurf_vtx_coord, cell_coord, sizeof(double) * n_cell * 3);

  isos->isosurf_vtx_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * isos->isosurf_n_vtx);
  memcpy(isos->isosurf_vtx_ln_to_gn, pcell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);

  int i_vtx = n_cell;
  for (int i = 0; i < n_face; i++) {
    if (face_tag[i] != 0) {
      isos->isosurf_vtx_ln_to_gn[i_vtx] = gn_cell + tagged_face_ln_to_gn[i_tagged_face[i]];
      memcpy(isos->isosurf_vtx_coord + 3*i_vtx,
             face_coord + 3*i,
             sizeof(double) * 3);
      i_vtx++;
    }
  }

  for (int i = 0; i < n_edge; i++) {
    if (edge_tag[i] != 0) {
      isos->isosurf_vtx_ln_to_gn[i_vtx] = gn_cell + gmax_tagged_face_ln_to_gn + tagged_edge_ln_to_gn[i_tagged_edge[i]];
      memcpy(isos->isosurf_vtx_coord + 3*i_vtx,
             edge_coord + 3*i,
             sizeof(double) * 3);
      i_vtx++;
    }
  }

  int _isosurf_n_tri = 0;
  for (int icell = 0; icell < n_cell; icell++) {
    for (int idx_face = pcell_face_idx[icell]; idx_face < pcell_face_idx[icell+1]; idx_face++) {
      int iface = PDM_ABS(pcell_face[idx_face]) - 1;

      if (face_tag[iface] == 0) {
        continue;
      }

      for (int idx_edge = pface_edge_idx[iface]; idx_edge < pface_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;

        if (edge_tag[iedge] != 0) {
          _isosurf_n_tri++;
        }
      }
    }
  }

  isos->isosurf_n_face = _isosurf_n_tri;
  isos->isosurf_face_vtx = (int *) malloc(sizeof(int) * _isosurf_n_tri * 3);

  PDM_g_num_t *distrib_isosurf_face = PDM_compute_entity_distribution(comm,
                                                                      _isosurf_n_tri);

  isos->isosurf_face_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * _isosurf_n_tri);

  _isosurf_n_tri = 0;
  for (int icell = 0; icell < n_cell; icell++) {
    for (int idx_face = pcell_face_idx[icell]; idx_face < pcell_face_idx[icell+1]; idx_face++) {
      int iface = PDM_ABS(pcell_face[idx_face]) - 1;

      if (face_tag[iface] == 0) {
        continue;
      }
      int sgn_face = PDM_SIGN(pcell_face[idx_face]);
      int ipt_face = n_cell + i_tagged_face[iface] + 1;

      for (int idx_edge = pface_edge_idx[iface]; idx_edge < pface_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(pface_edge[idx_edge]) - 1;

        if (edge_tag[iedge] != 0) {

          int sgn_edge = PDM_SIGN(pface_edge[idx_edge]);
          int ipt_edge = n_cell + n_tagged_face + i_tagged_edge[iedge] + 1;

          int sgn = sgn_edge * sgn_face * edge_tag[iedge];

          isos->isosurf_face_vtx[3*_isosurf_n_tri] = icell + 1;
          if (sgn < 0) {
            isos->isosurf_face_vtx[3*_isosurf_n_tri+1] = ipt_face;
            isos->isosurf_face_vtx[3*_isosurf_n_tri+2] = ipt_edge;
          } else {
            isos->isosurf_face_vtx[3*_isosurf_n_tri+1] = ipt_edge;
            isos->isosurf_face_vtx[3*_isosurf_n_tri+2] = ipt_face;
          }

          isos->isosurf_face_ln_to_gn[_isosurf_n_tri] = distrib_isosurf_face[i_rank] + _isosurf_n_tri + 1;
          _isosurf_n_tri++;

        }
      }
    }
  }
  free(edge_tag     );
  free(edge_coord   );
  free(edge_gradient);

  free(distrib_isosurf_face);
  free(tagged_face_parent_ln_to_gn);
  free(tagged_edge_parent_ln_to_gn);
  PDM_gnum_free(gen_gnum_face);
  PDM_gnum_free(gen_gnum_edge);


  isos->isosurf_face_vtx_idx = (int *) malloc(sizeof(int) * (isos->isosurf_n_face + 1));
  for (int i = 0; i <= isos->isosurf_n_face; i++) {
    isos->isosurf_face_vtx_idx[i] = 3*i;
  }

  if (isos->debug) {
    sprintf(filename, "isosurf_tri_%2.2d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               isos->isosurf_n_vtx,
                               isos->isosurf_vtx_coord,
                               isos->isosurf_vtx_ln_to_gn,
                               PDM_MESH_NODAL_TRIA3,
                               isos->isosurf_n_face,
                               isos->isosurf_face_vtx,
                               isos->isosurf_face_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }

  free(i_tagged_edge);
  free(i_tagged_face);

  free (is_treated_face);
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
  double t1, t2;
  double delta_t;
  double delta_max;
  double delta_min;

  /*
   * Select gnum that contains iso-surface
   */
  int n_rank;
  int i_rank;
  PDM_MPI_Comm_size(isos->comm, &n_rank);
  PDM_MPI_Comm_rank(isos->comm, &i_rank);

  if (isos->debug) {
    char filename[999];
    sprintf(filename, "gradient_field_%2.2d.vtk", i_rank);
    _dump_vectors (filename,
                   (int) isos->distrib_vtx[i_rank+1] - isos->distrib_vtx[i_rank],
                   isos->dvtx_coord,
                   isos->dgradient_field,
                   NULL);
  }

  t1 = PDM_MPI_Wtime();

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
  if (0 && isos->debug) {
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

  t2 = PDM_MPI_Wtime();
  delta_t = t2 - t1;
  PDM_MPI_Allreduce (&delta_t, &delta_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, isos->comm);
  PDM_MPI_Allreduce (&delta_t, &delta_min, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, isos->comm);
  if(i_rank == 0) {
    printf("PDM_iso_surface : scan edges : min/max = %12.5e / %12.5e\n", delta_min, delta_max);
  }

  t1 = PDM_MPI_Wtime();

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
  // int         *unique_dentity_edge_order = NULL;
  // PDM_g_num_t *unique_dentity_edge       = NULL;
  // int dn_edge_sort = PDM_unique_long_with_distrib(isos->comm,
  //                                                 dentity_edge,
  //                                                 isos->distrib_edge,
  //                                                 dentity_edge_idx[dn_entity],
  //                                                 &unique_dentity_edge_order,
  //                                                 &unique_dentity_edge);



  PDM_block_to_part_t* btp = PDM_block_to_part_create(isos->distrib_edge,
                               (const PDM_g_num_t **) &dentity_edge,
                                                      &dentity_edge_idx[dn_entity],
                                                      1,
                                                      isos->comm);

  // PDM_block_to_part_t* btp = PDM_block_to_part_create(isos->distrib_edge,
  //                              (const PDM_g_num_t **) &unique_dentity_edge,
  //                                                     &dn_edge_sort,
  //                                                     1,
  //                                                     isos->comm);

  // free(unique_dentity_edge);

  if(isos->dim == 3) {
    free(dentity_edge);
  }

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

      // int idx_read = unique_dentity_edge_order[idx_entity];
      int idx_read = idx_entity;
      if(dentity_edge_tag[idx_read] == 1) {
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
        // int idx_read = unique_dentity_edge_order[idx_entity];
        int idx_read = idx_entity;
        dentity_center[3*idx_write  ] += dentity_edge_center[3*idx_read  ];
        dentity_center[3*idx_write+1] += dentity_edge_center[3*idx_read+1];
        dentity_center[3*idx_write+2] += dentity_edge_center[3*idx_read+2];
      }
      dentity_center[3*idx_write  ] = dentity_center[3*idx_write  ] * inv;
      dentity_center[3*idx_write+1] = dentity_center[3*idx_write+1] * inv;
      dentity_center[3*idx_write+2] = dentity_center[3*idx_write+2] * inv;

      idx_write++;
    }
  }
  free(dentity_edge_center);
  // free(unique_dentity_edge_order);
  PDM_block_to_part_free(btp);

  free(dentity_tag);

  if(isos->dim == 3) {
    free(dentity_edge_idx);
  }
  free(dentity_edge_tag);


  /*
   * Rebuild partition that contains entity and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_entity_tag, dentity_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_entity_gnum = PDM_gnum_get(gnum_equi, 0);

  if (isos->debug) {
    char filename[999];
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


    if (isos->dgradient_field != NULL) {
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
  }
  PDM_block_to_part_free(btp_vtx);


  t2 = PDM_MPI_Wtime();
  delta_t = t2 - t1;
  PDM_MPI_Allreduce (&delta_t, &delta_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, isos->comm);
  PDM_MPI_Allreduce (&delta_t, &delta_min, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, isos->comm);
  if(i_rank == 0) {
    printf("PDM_iso_surface : extraction : min/max = %12.5e / %12.5e\n", delta_min, delta_max);
  }

  t1 = PDM_MPI_Wtime();



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
                   pequi_gradient_field);
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
                   pequi_gradient_field);
  }

  t2 = PDM_MPI_Wtime();
  delta_t = t2 - t1;
  PDM_MPI_Allreduce (&delta_t, &delta_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, isos->comm);
  PDM_MPI_Allreduce (&delta_t, &delta_min, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, isos->comm);
  if(i_rank == 0) {
    printf("PDM_iso_surface : build isosurface : min/max = %12.5e / %12.5e\n", delta_min, delta_max);
  }

  t1 = PDM_MPI_Wtime();

  if(isos->dim == 3) {
    free(pequi_parent_face_ln_to_gn);
    free(pequi_face_ln_to_gn);
    free(pequi_cell_face_idx);
    free(pequi_cell_face);
  }


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

}



static
void
_iso_surface_part
(
  PDM_iso_surface_t        *isos
)
{
  /*
   *  TO DO: 2d case
   */

  /*
   *  Scan edges to find those which intersect the iso-surface
   */
  int **edge_tag = (int **) malloc(sizeof(int *) * isos->n_part);

  for (int i_part = 0; i_part < isos->n_part; i_part++) {

    edge_tag[i_part] = (int *) malloc(sizeof(int) * isos->n_edge[i_part]);

    for (int i = 0; i < isos->n_edge[i_part]; i++) {

      int i_vtx1 = isos->pedge_vtx[i_part][2*i  ]-1;
      int i_vtx2 = isos->pedge_vtx[i_part][2*i+1]-1;

      edge_tag[i_part][i] = 0;

      double x1 = isos->pvtx_coord[i_part][3*i_vtx1  ];
      double y1 = isos->pvtx_coord[i_part][3*i_vtx1+1];
      double z1 = isos->pvtx_coord[i_part][3*i_vtx1+2];

      double x2 = isos->pvtx_coord[i_part][3*i_vtx2  ];
      double y2 = isos->pvtx_coord[i_part][3*i_vtx2+1];
      double z2 = isos->pvtx_coord[i_part][3*i_vtx1+2];

      double val1 = 0;
      double val2 = 0;
      if (isos->iso_kind == PDM_ISO_SURFACE_KIND_PLANE) {
        val1 = _plane_field(x1, y1, z1, isos->plane_equation);
        val2 = _plane_field(x2, y2, z2, isos->plane_equation);
      }
      else {
        val1 = isos->pfield[i_part][i_vtx1];
        val2 = isos->pfield[i_part][i_vtx2];
      }

      int sgn1 = PDM_SIGN(val1);
      int sgn2 = PDM_SIGN(val2);

      if (sgn1 * sgn2 < 0) {
        edge_tag[i_part][i] = 1;
      }

    }
  }


  /*
   *  Scan cells/faces to find those which intersect the iso-surface
   */
  int  *n_extract_elt = (int *)  malloc(sizeof(int)   * isos->n_part);
  int **extract_elt   = (int **) malloc(sizeof(int *) * isos->n_part);

  if (isos->dim == 3) {
    for (int i_part = 0; i_part < isos->n_part; i_part++) {

      n_extract_elt[i_part] = 0;
      extract_elt[i_part] = (int *) malloc(sizeof(int) * isos->n_cell[i_part]);
      for (int icell = 0; icell < isos->n_cell[i_part]; icell++) {

        int is_tagged = 0;

        for (int idx_face = isos->pcell_face_idx[i_part][icell]; idx_face < isos->pcell_face_idx[i_part][icell+1]; idx_face++) {

          int iface = PDM_ABS(isos->pcell_face[i_part][idx_face]) - 1;

          for (int idx_edge = isos->pface_edge_idx[i_part][iface]; idx_edge < isos->pface_edge_idx[i_part][iface+1]; idx_edge++) {

            int iedge = PDM_ABS(isos->pface_edge[i_part][idx_edge]) - 1;

            if (edge_tag[i_part][iedge] != 0) {
              is_tagged = 1;
              extract_elt[i_part][n_extract_elt[i_part]++] = icell;
              break;
            }

          }

          if (is_tagged) break;
        }

      }

      extract_elt[i_part] = realloc(extract_elt[i_part], sizeof(int) * n_extract_elt[i_part]);
    }
  }

  else if (isos->dim == 2) {

    for (int i_part = 0; i_part < isos->n_part; i_part++) {

      n_extract_elt[i_part] = 0;
      extract_elt[i_part] = (int *) malloc(sizeof(int) * isos->n_face[i_part]);
      for (int iface = 0; iface < isos->n_face[i_part]; iface++) {

        for (int idx_edge = isos->pface_edge_idx[i_part][iface]; idx_edge < isos->pface_edge_idx[i_part][iface+1]; idx_edge++) {

          int iedge = PDM_ABS(isos->pface_edge[i_part][idx_edge]) - 1;

          if (edge_tag[i_part][iedge] != 0) {
            extract_elt[i_part][n_extract_elt[i_part]++] = iface;
            break;
          }

        }

      }

      extract_elt[i_part] = realloc(extract_elt[i_part], sizeof(int) * n_extract_elt[i_part]);
    }
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "invalid dimension %d\n", isos->dim);
  }


  /*
   *  Extract 'partition'
   */

  int n_part_out = 1;
  PDM_extract_part_t *extrp = PDM_extract_part_create(isos->dim,//3,
                                                      isos->n_part,
                                                      n_part_out,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_OWNERSHIP_KEEP,
                                                      isos->comm);

  for (int i_part = 0; i_part < isos->n_part; ++i_part) {

    PDM_extract_part_part_set(extrp,
                              i_part,
                              isos->n_cell[i_part],
                              isos->n_face[i_part],
                              isos->n_edge[i_part],
                              isos->n_vtx[i_part],
                              isos->pcell_face_idx[i_part],
                              isos->pcell_face[i_part],
                              isos->pface_edge_idx[i_part],
                              isos->pface_edge[i_part],
                              isos->pedge_vtx[i_part],
                              NULL,//isos->pface_vtx_idx[i_part],
                              NULL,//isos->pface_vtx[i_part],
                              isos->cell_ln_to_gn[i_part],
                              isos->face_ln_to_gn[i_part],
                              isos->edge_ln_to_gn[i_part],
                              isos->vtx_ln_to_gn[i_part],
                              isos->pvtx_coord[i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       n_extract_elt[i_part],
                                       extract_elt[i_part]);
  }

  PDM_extract_part_compute(extrp);

  int          n_cell = 0;
  int          n_face = 0;
  int          n_edge = 0;
  int          n_vtx  = 0;
  int         *pcell_face_idx  = NULL;
  int         *pcell_face      = NULL;
  int         *pface_edge_idx  = NULL;
  int         *pface_edge      = NULL;
  int         *pedge_vtx       = NULL;
  PDM_g_num_t *pcell_ln_to_gn  = NULL;
  PDM_g_num_t *pface_ln_to_gn  = NULL;
  PDM_g_num_t *pedge_ln_to_gn  = NULL;
  PDM_g_num_t *pvtx_ln_to_gn   = NULL;
  double      *pvtx_coord      = NULL;
  double      *pfield          = NULL;
  double      *pgradient_field = NULL;
  PDM_g_num_t *pcell_parent_ln_to_gn = NULL;
  PDM_g_num_t *pface_parent_ln_to_gn = NULL;
  PDM_g_num_t *pedge_parent_ln_to_gn = NULL;
  PDM_g_num_t *pvtx_parent_ln_to_gn  = NULL;


  int *tmp_n;
  PDM_extract_part_n_entity_get(extrp,
                                PDM_MESH_ENTITY_CELL,
                                &tmp_n);
  n_cell = tmp_n[0];
  free(tmp_n);

  PDM_extract_part_n_entity_get(extrp,
                                PDM_MESH_ENTITY_FACE,
                                &tmp_n);
  n_face = tmp_n[0];
  free(tmp_n);

  PDM_extract_part_n_entity_get(extrp,
                                PDM_MESH_ENTITY_EDGE,
                                &tmp_n);
  n_edge = tmp_n[0];
  free(tmp_n);

  PDM_extract_part_n_entity_get(extrp,
                                PDM_MESH_ENTITY_VERTEX,
                                &tmp_n);
  n_vtx = tmp_n[0];
  free(tmp_n);

  /* Connectivities */
  int **tmp_cell_face_idx;
  int **tmp_cell_face;
  PDM_extract_part_connectivity_get(extrp,
                                    PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                    &tmp_cell_face,
                                    &tmp_cell_face_idx,
                                    PDM_OWNERSHIP_KEEP);
  pcell_face_idx = tmp_cell_face_idx[0];
  pcell_face     = tmp_cell_face[0];
  // free(tmp_cell_face_idx);
  // free(tmp_cell_face);


  int **tmp_face_edge_idx;
  int **tmp_face_edge;
  PDM_extract_part_connectivity_get(extrp,
                                    PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                    &tmp_face_edge,
                                    &tmp_face_edge_idx,
                                    PDM_OWNERSHIP_KEEP);
  pface_edge_idx = tmp_face_edge_idx[0];
  pface_edge     = tmp_face_edge[0];
  // free(tmp_face_edge_idx);
  // free(tmp_face_edge);


  int **tmp_edge_vtx_idx;
  int **tmp_edge_vtx;
  PDM_extract_part_connectivity_get(extrp,
                                    PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                    &tmp_edge_vtx,
                                    &tmp_edge_vtx_idx,
                                    PDM_OWNERSHIP_KEEP);
  pedge_vtx = tmp_edge_vtx[0];
  // free(tmp_edge_vtx);


  /* Coordinates */
  double **tmp_vtx_coord;
  PDM_extract_part_vtx_coord_get(extrp,
                                 &tmp_vtx_coord,
                                 PDM_OWNERSHIP_KEEP);
  pvtx_coord = tmp_vtx_coord[0];
  // free(tmp_vtx_coord);


  /* (Child) Global ids */
  PDM_g_num_t **tmp_cell_ln_to_gn;
  PDM_extract_part_ln_to_gn_get(extrp,
                                PDM_MESH_ENTITY_CELL,
                                &tmp_cell_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
  pcell_ln_to_gn = tmp_cell_ln_to_gn[0];
  // free(tmp_cell_ln_to_gn);


  PDM_g_num_t **tmp_face_ln_to_gn;
  PDM_extract_part_ln_to_gn_get(extrp,
                                PDM_MESH_ENTITY_FACE,
                                &tmp_face_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
  pface_ln_to_gn = tmp_face_ln_to_gn[0];
  // free(tmp_face_ln_to_gn);


  PDM_g_num_t **tmp_edge_ln_to_gn;
  PDM_extract_part_ln_to_gn_get(extrp,
                                PDM_MESH_ENTITY_EDGE,
                                &tmp_edge_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
  pedge_ln_to_gn = tmp_edge_ln_to_gn[0];
  // free(tmp_edge_ln_to_gn);


  PDM_g_num_t **tmp_vtx_ln_to_gn;
  PDM_extract_part_ln_to_gn_get(extrp,
                                PDM_MESH_ENTITY_VERTEX,
                                &tmp_vtx_ln_to_gn,
                                PDM_OWNERSHIP_KEEP);
  pvtx_ln_to_gn = tmp_vtx_ln_to_gn[0];
  // free(tmp_vtx_ln_to_gn);



  /* Parent Global ids */
  PDM_g_num_t **tmp_cell_parent_ln_to_gn;
  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       PDM_MESH_ENTITY_CELL,
                                       &tmp_cell_parent_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
  pcell_parent_ln_to_gn = tmp_cell_parent_ln_to_gn[0];
  // free(tmp_cell_parent_ln_to_gn);


  PDM_g_num_t **tmp_face_parent_ln_to_gn;
  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       PDM_MESH_ENTITY_FACE,
                                       &tmp_face_parent_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
  pface_parent_ln_to_gn = tmp_face_parent_ln_to_gn[0];
  // free(tmp_face_parent_ln_to_gn);


  PDM_g_num_t **tmp_edge_parent_ln_to_gn;
  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       PDM_MESH_ENTITY_EDGE,
                                       &tmp_edge_parent_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
  pedge_parent_ln_to_gn = tmp_edge_parent_ln_to_gn[0];
  // free(tmp_edge_parent_ln_to_gn);


  PDM_g_num_t **tmp_vtx_parent_ln_to_gn;
  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       PDM_MESH_ENTITY_VERTEX,
                                       &tmp_vtx_parent_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
  pvtx_parent_ln_to_gn = tmp_vtx_parent_ln_to_gn[0];
  // free(tmp_vtx_parent_ln_to_gn);


  /* Field and gradient */
  int *child_to_parent_idx = (int *) malloc(sizeof(int) * (n_vtx + 1));
  child_to_parent_idx[0] = 0;
  for (int i = 0; i < n_vtx; i++) {
    child_to_parent_idx[i+1] = child_to_parent_idx[i] + 1;
  }
  PDM_part_to_part_t *ptp = PDM_part_to_part_create((const PDM_g_num_t **) &pvtx_ln_to_gn,
                                                    &n_vtx,
                                                    1,
                                                    (const PDM_g_num_t **) isos->vtx_ln_to_gn,
                                                    isos->n_vtx,
                                                    isos->n_part,
                                                    (const int **) &child_to_parent_idx,
                                                    (const PDM_g_num_t **) &pvtx_parent_ln_to_gn,
                                                    isos->comm);


  int request1;
  double **tmp_pfield;
  PDM_part_to_part_reverse_iexch(ptp,
                                 PDM_MPI_COMM_KIND_P2P,
                                 PDM_STRIDE_CST_INTERLACED,
                                 PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                 1,
                                 sizeof(double),
                                 NULL,
                 (const void **) isos->pfield,
                                 NULL,
                     (void ***) &tmp_pfield,
                                &request1);

  int request2;
  double **tmp_pgradient_field;
  if (isos->pgradient_field[0] != NULL) {
    PDM_part_to_part_reverse_iexch(ptp,
                                   PDM_MPI_COMM_KIND_P2P,
                                   PDM_STRIDE_CST_INTERLACED,
                                   PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                   1,
                                   3*sizeof(double),
                                   NULL,
                                   (const void **) isos->pgradient_field,
                                   NULL,
                                   (void ***) &tmp_pgradient_field,
                                   &request2);
  }

  PDM_part_to_part_reverse_iexch_wait(ptp, request1);
  pfield = tmp_pfield[0];
  free(tmp_pfield);

  if (isos->pgradient_field[0] != NULL) {
    PDM_part_to_part_reverse_iexch_wait(ptp, request2);
    pgradient_field = tmp_pgradient_field[0];
    free(tmp_pgradient_field);
  }

  PDM_part_to_part_free(ptp);
  free(child_to_parent_idx);



  // PDM_extract_part_free(extrp);

  for (int i_part = 0; i_part < isos->n_part; i_part++) {
    free(edge_tag[i_part]);
    free(extract_elt[i_part]);
  }
  free(edge_tag);
  free(extract_elt);
  free(n_extract_elt);



  /*
   *  Compute iso-surface
   */
  if (isos->dim == 3) {
    _iso_surf_dist(isos,
                   n_cell,
                   n_face,
                   n_edge,
                   n_vtx,
                   pcell_face_idx,
                   pcell_face,
                   pface_edge_idx,
                   pface_edge,
                   pedge_vtx,
                   pcell_ln_to_gn,
                   pface_ln_to_gn,
                   pedge_ln_to_gn,
                   pvtx_ln_to_gn,
                   pvtx_coord,
                   pfield,
                   pgradient_field);
  }
  else if (isos->dim == 2) {
    _iso_line_dist(isos,
                   n_face,
                   n_edge,
                   n_vtx,
                   pface_edge_idx,
                   pface_edge,
                   pedge_vtx,
                   pface_ln_to_gn,
                   pedge_ln_to_gn,
                   pvtx_ln_to_gn,
                   pvtx_coord,
                   pfield,
                   pgradient_field);
  }


  /*
   *  Free memory
   */
  // free(pcell_face_idx);
  // free(pcell_face);
  // free(pface_edge_idx);
  // free(pface_edge);
  // free(pedge_vtx);
  // free(pcell_ln_to_gn);
  // free(pface_ln_to_gn);
  // free(pedge_ln_to_gn);
  // free(pvtx_ln_to_gn);
  // free(pvtx_coord);
  PDM_extract_part_free(extrp);
  free(pfield);
  if (pgradient_field != NULL) {
    free(pgradient_field);
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

  for (int i = 0; i < n_part; i++) {
    isos->pgradient_field[i] = NULL;
  }
  isos->dgradient_field = NULL;


  isos->isosurf_n_vtx         = 0;
  isos->isosurf_n_edge        = 0;
  isos->isosurf_n_face        = 0;
  isos->isosurf_face_vtx_idx  = NULL;
  isos->isosurf_face_vtx      = NULL;
  isos->isosurf_edge_vtx_idx  = NULL;
  isos->isosurf_edge_vtx      = NULL;
  isos->isosurf_vtx_coord     = NULL;
  isos->isosurf_face_ln_to_gn = NULL;
  isos->isosurf_vtx_ln_to_gn  = NULL;

  isos->debug = 0;

  isos->eval_field_and_gradient = NULL;

  return isos;
}

void
PDM_iso_surface_compute
(
  PDM_iso_surface_t        *isos
)
{
  double t1 = PDM_MPI_Wtime();

  if(isos->is_dist == 0) {
    _iso_surface_part(isos);
  } else {
    _iso_surface_dist(isos);
  }

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;
  double delta_max;
  double delta_min;

  PDM_MPI_Allreduce (&delta_t, &delta_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, isos->comm);
  PDM_MPI_Allreduce (&delta_t, &delta_min, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, isos->comm);

  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  if(i_rank == 0) {
    printf("PDM_iso_surface : min/max = %12.5e / %12.5e\n", delta_min, delta_max);
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

  if (isos->ownership == PDM_OWNERSHIP_KEEP) {

    if (isos->isosurf_face_vtx_idx  != NULL) {
      free(isos->isosurf_face_vtx_idx);
    }
    if (isos->isosurf_face_vtx      != NULL) {
      free(isos->isosurf_face_vtx);
    }
    if (isos->isosurf_edge_vtx_idx  != NULL) {
      free(isos->isosurf_edge_vtx_idx);
    }
    if (isos->isosurf_edge_vtx      != NULL) {
      free(isos->isosurf_edge_vtx);
    }
    if (isos->isosurf_vtx_coord     != NULL) {
      free(isos->isosurf_vtx_coord);
    }
    if (isos->isosurf_face_ln_to_gn != NULL) {
      free(isos->isosurf_face_ln_to_gn);
    }
    if (isos->isosurf_vtx_ln_to_gn  != NULL) {
      free(isos->isosurf_vtx_ln_to_gn);
    }
  }

  free(isos);
}




void
PDM_iso_surface_write
(
 PDM_iso_surface_t  *isos,
 const char         *name
 )
 {
  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);


  char _name[999];

  PDM_writer_t *cs = PDM_writer_create("Ensight",
                                       PDM_WRITER_FMT_ASCII,
                                       PDM_WRITER_TOPO_CONSTANTE,
                                       PDM_WRITER_OFF,
                                       name,
                                       name,
                                       PDM_MPI_COMM_WORLD,
                                       PDM_IO_ACCES_MPI_SIMPLE,
                                       1.,
                                       NULL);

  sprintf(_name, "%s_geom", name);
  int id_geom = PDM_writer_geom_create(cs,
                                       _name,
                                       PDM_WRITER_OFF,
                                       PDM_WRITER_OFF,
                                       1);

  int id_var_part = PDM_writer_var_create(cs,
                                          PDM_WRITER_ON,
                                          PDM_WRITER_VAR_SCALAIRE,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "num_part");

  int id_var_vtx_gnum = PDM_writer_var_create(cs,
                                              PDM_WRITER_ON,
                                              PDM_WRITER_VAR_SCALAIRE,
                                              PDM_WRITER_VAR_SOMMETS,
                                              "vtx_g_num");

  int id_var_elt_gnum = PDM_writer_var_create(cs,
                                              PDM_WRITER_ON,
                                              PDM_WRITER_VAR_SCALAIRE,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "elt_g_num");

  PDM_writer_step_beg(cs, 0.);


  /*
   *  Write geometry
   */
  PDM_writer_geom_coord_set(cs,
                            id_geom,
                            0,
                            isos->isosurf_n_vtx,
                            isos->isosurf_vtx_coord,
                            isos->isosurf_vtx_ln_to_gn);

  int n_elt = 0;
  PDM_g_num_t *elt_ln_to_gn = NULL;
  if (isos->dim == 3) {
    n_elt = isos->isosurf_n_face;
    elt_ln_to_gn = isos->isosurf_face_ln_to_gn;

    int id_block = PDM_writer_geom_bloc_add(cs,
                                            id_geom,
                                            PDM_WRITER_ON,
                                            PDM_WRITER_POLY_2D);


    PDM_writer_geom_bloc_poly2d_set(cs,
                                    id_geom,
                                    id_block,
                                    0,
                                    isos->isosurf_n_face,
                                    isos->isosurf_face_vtx_idx,
                                    isos->isosurf_face_vtx,
                                    isos->isosurf_face_ln_to_gn);
  }

  else {
    n_elt = isos->isosurf_n_edge;
    elt_ln_to_gn = isos->isosurf_edge_ln_to_gn;

    int id_block = PDM_writer_geom_bloc_add(cs,
                                            id_geom,
                                            PDM_WRITER_ON,
                                            PDM_WRITER_BAR2);

    PDM_writer_geom_bloc_std_set(cs,
                                 id_geom,
                                 id_block,
                                 0,
                                 isos->isosurf_n_edge,
                                 isos->isosurf_edge_vtx,
                                 isos->isosurf_edge_ln_to_gn);
  }

  PDM_writer_geom_write(cs,
                          id_geom);


  /*
   *  Write variables
   */
  PDM_real_t *val_part     = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_elt);
  PDM_real_t *val_elt_gnum = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_elt);
  PDM_real_t *val_vtx_gnum = (PDM_real_t *) malloc(sizeof(PDM_real_t) * isos->isosurf_n_vtx);

  for (int i = 0; i < n_elt; i++) {
    val_part[i]     = (PDM_real_t) i_rank;
    val_elt_gnum[i] = (PDM_real_t) elt_ln_to_gn[i];
  }

  for (int i = 0; i < isos->isosurf_n_vtx; i++) {
    val_vtx_gnum[i] = (PDM_real_t) isos->isosurf_vtx_ln_to_gn[i];
  }

  PDM_writer_var_set(cs,
                     id_var_part,
                     id_geom,
                     0,
                     (const PDM_real_t *) val_part);
  PDM_writer_var_write(cs,
                       id_var_part);

  PDM_writer_var_set(cs,
                     id_var_elt_gnum,
                     id_geom,
                     0,
                     (const PDM_real_t *) val_elt_gnum);
  PDM_writer_var_write(cs,
                       id_var_elt_gnum);

  PDM_writer_var_set(cs,
                     id_var_vtx_gnum,
                     id_geom,
                     0,
                     (const PDM_real_t *) val_vtx_gnum);
  PDM_writer_var_write(cs,
                       id_var_vtx_gnum);



  PDM_writer_step_end(cs);

  // PDM_writer_free(cs);

  free(val_part);
  free(val_elt_gnum);
  free(val_vtx_gnum);
 }



void
PDM_iso_surface_eval_field_and_gradient_set
(
 PDM_iso_surface_t *isos,
 void (*eval_field_and_gradient) (const double, const double, const double,
                                   double *,
                                   double *, double *, double *)
 )
{
  assert(isos != NULL);
  isos->eval_field_and_gradient = eval_field_and_gradient;
}


void
PDM_iso_surface_surface_get
(
 PDM_iso_surface_t  *isos,
 int                *n_vtx,
 int                *n_elt,
 int               **elt_vtx_idx,
 int               **elt_vtx,
 double            **vtx_coord,
 PDM_g_num_t       **elt_ln_to_gn,
 PDM_g_num_t       **vtx_ln_to_gn
 )
{
  assert(isos != NULL);

  *n_vtx        = isos->isosurf_n_vtx;
  *vtx_coord    = isos->isosurf_vtx_coord;
  *vtx_ln_to_gn = isos->isosurf_vtx_ln_to_gn;

  if (isos->dim == 2) {
    *n_elt        = isos->isosurf_n_edge;
    *elt_vtx      = isos->isosurf_edge_vtx;
    *elt_ln_to_gn = isos->isosurf_edge_ln_to_gn;
  }

  else {
    *n_elt        = isos->isosurf_n_face;
    *elt_vtx_idx  = isos->isosurf_face_vtx_idx;
    *elt_vtx      = isos->isosurf_face_vtx;
    *elt_ln_to_gn = isos->isosurf_face_ln_to_gn;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
