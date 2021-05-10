#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_error.h"
#include "pdm_plane.h"
#include "pdm_surf_mesh.h"
#include "pdm_triangulate.h"
#include "pdm_polygon.h"
#include "pdm_edges_intersect.h"

static void
_read_args(int     argc,
           char  **argv,
           char  **filename,
           int    *n_pts_seg)
{
  int i = 1;

  while (i < argc) {
    if (strcmp(argv[i], "-f") == 0) {
      i++;
      *filename = argv[i];
    }

    if (strcmp(argv[i], "-n") == 0) {
      i++;
      *n_pts_seg = atoi(argv[i]);
    }

    i++;
  }
}


static int
_compute_normals
(
 const int      n_face,
 const int     *face_vtx_idx,
 const int     *face_vtx,
 const int      n_vtx,
 const double  *vtx_coord,
 double       **vtx_normal
 )
{
  *vtx_normal = malloc (sizeof(double) * n_vtx * 3);
  for (int i = 0; i < n_vtx*3; i++) {
    (*vtx_normal)[i] = 0.;
  }

  for (int j = 0; j < n_face; j++) {
    double face_normal[3] = {0., 0., 0.};
    double subNorm[3], v1[3], v2[3];

    int _n_vtx = face_vtx_idx[j+1] - face_vtx_idx[j];
    int ideb = face_vtx_idx[j];

    int idx0 = face_vtx[ideb] - 1;
    for (int k = 1; k < _n_vtx - 1; k++) {
      int idx1 = face_vtx[ideb + k] - 1;
      int idx2 = face_vtx[ideb + (k+1) % _n_vtx] - 1;
      for (int k1 = 0; k1 < 3; k1++) {
        v1[k1] = vtx_coord[3*idx1+k1] - vtx_coord[3*idx0+k1];
        v2[k1] = vtx_coord[3*idx2+k1] - vtx_coord[3*idx0+k1];
      }
      PDM_CROSS_PRODUCT (subNorm, v1, v2);
      for (int k1 = 0; k1 < 3; k1++) {
        face_normal[k1] += subNorm[k1];
      }
    }

    for (int i = face_vtx_idx[j]; i < face_vtx_idx[j+1]; i++) {
      int ivtx = face_vtx[i] - 1;
      for (int k = 0; k < 3; k++) {
        (*vtx_normal)[3*ivtx+k] += face_normal[k];
      }
    }
  }


  // Normalize
  int singular = 0;
  double *vn = *vtx_normal;
  for (int i = 0; i < n_vtx; i++) {
    double mag = PDM_DOT_PRODUCT (vn, vn);
    if (mag > 0.) {
      mag = 1. / sqrt(mag);
      for (int j = 0; j < 3; j++) {
        vn[j] *= mag;
      }
    } else {
      singular = 1;
    }
    vn += 3;
  }

  return singular;
}





static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

static inline double
_determinant_3x3
(
 const double a[3],
 const double b[3],
 const double c[3]
 )
{
  return a[0] * (b[1]*c[2] - b[2]*c[1])
    +    a[1] * (b[2]*c[0] - b[0]*c[2])
    +    a[2] * (b[0]*c[1] - b[1]*c[0]);
}

/*---------------------------------------------------------------------------
 * Solve the equation "A.x = b" with Cramer's rule.
 *
 * parameters:
 *   A[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   PDM_FALSE if matrix is singular, PDM_TRUE otherwise
 *----------------------------------------------------------------------------*/

static int
_solve_3x3(double  A[3][3],
           double  b[3],
           double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det = A[0][0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
    -   A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    +   A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  }

  else {
    det_inv = 1./det;
  }

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(A[1][1]*A[2][2] - A[2][1]*A[1][2])
          - b[1]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
          + b[2]*(A[0][1]*A[1][2] - A[1][1]*A[0][2])) * det_inv;

  x1 = (  A[0][0]*(b[1]*A[2][2] - b[2]*A[1][2])
          - A[1][0]*(b[0]*A[2][2] - b[2]*A[0][2])
          + A[2][0]*(b[0]*A[1][2] - b[1]*A[0][2])) * det_inv;

  x2 = (  A[0][0]*(A[1][1]*b[2] - A[2][1]*b[1])
          - A[1][0]*(A[0][1]*b[2] - A[2][1]*b[0])
          + A[2][0]*(A[0][1]*b[1] - A[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;

  return PDM_TRUE;
}


/**
 * Project a point q on a triangle with projection direction linearly interpolated from those
 * given at the triangle's vertices
 **/

static int
_projection_on_triangle_BtoA
(
 const double  q[3],
 const double  x0[3],
 const double  x1[3],
 const double  x2[3],
 const double  d0[3],
 const double  d1[3],
 const double  d2[3],
 double       *u,
 double       *v,
 double       *g,
 const int     verbose
 )
{
  const double told = 1e-12;
  const double tolr = 1e-6;
  const int max_iter = 10;

  /* Initial iterate */
  double _u = 1./3.;
  double _v = _u;
  double _g = 0.;

  /* Newton-Raphson iteration */
  double r[3];
  double J[3][3];
  double duvg[3] = {0., 0., 0.};
  int converged = 0;

  for (int iter = 0; iter < max_iter; iter++) {

    double w = 1. - _u - _v;

    // Jacobian matrix
    for (int i = 0; i < 3; i++) {
      J[i][0] = x0[i] - x2[i] + _g*(d2[i] - d0[i]);
      J[i][1] = x1[i] - x2[i] + _g*(d2[i] - d1[i]);
      J[i][2] = -(_u*d0[i] + _v*d1[i] + w*d2[i]);
    }

    // Right-hand-side term
    for (int i = 0; i < 3; i++) {
      r[i] = q[i] - (_u*x0[i] + _v*x1[i] + w*x2[i] + _g*J[i][2]);
    }
    double resr = PDM_DOT_PRODUCT (r, r);

    // Solve for Newton step
    int stat = _solve_3x3 (J, r, duvg);

    if (stat == 0) {
      // singular Jacobian matrix
      break;
    }

    double resd = PDM_MAX (PDM_ABS(duvg[0]), PDM_ABS(duvg[1]));
    //resd = PDM_MAX (PDM_ABS(duvg[2]), resd);
    if (verbose) printf("  it %d, u = %f, v = %f, resr = %f, resd = %f\n", iter, _u, _v, resr, resd);

    // update iterate
    _u += duvg[0];
    _v += duvg[1];
    _g += duvg[2];

    // Check for termination
    if (resd < told) {
      if (resr < tolr*tolr) {//?
        converged = 1;
      }
      break;
    }
  }

  *u = _u;
  *v = _v;
  *g = _g;

  return converged;
}


static int
_projection_on_triangle_AtoB
(
 const double  p[3],
 const double  n[3],
 const double  x0[3],
 const double  x1[3],
 const double  x2[3],
 double       *u,
 double       *v,
 double       *g
 )
{
  double A[3][3];
  double b[3];

  for (int i = 0; i < 3; i++) {
    A[i][0] = x0[i] - x2[i];
    A[i][1] = x1[i] - x2[i];
    A[i][2] = -n[i];
    b[i]    = p[i] - x2[i];
  }

  double uvg[3];
  int stat = _solve_3x3 (A, b, uvg);

  *u = uvg[0];
  *v = uvg[1];
  *g = uvg[2];

  return stat;
}


static void
_projection_on_polygon_BtoA
(
 const int             n_pts,
 const double         *pts_coord,
 const int             n_vtx,
 const double         *vtx_coord,
 const double         *vtx_normal,
 double               *proj_coord,
 PDM_polygon_status_t *pts_status,
 const int             verbose
 )
{
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx);
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx - 2)*3);

  int n_tri = PDM_triangulate_polygon (3,
                                       n_vtx,
                                       vtx_coord,
                                       NULL,
                                       NULL,
                                       PDM_TRIANGULATE_ELT_DEF,
                                       tri_vtx,
                                       state);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    pts_status[ipt] = PDM_POLYGON_OUTSIDE;
  }

  /* Loop on sub-triangles */
  for (int itri = 0; itri < n_tri; itri++) {

    int i0 = 3 * (tri_vtx[3*itri]     - 1);
    int i1 = 3 * (tri_vtx[3*itri  +1] - 1);
    int i2 = 3 * (tri_vtx[3*itri + 2] - 1);

    /* Loop on points */
    for (int ipt = 0; ipt < n_pts; ipt++) {

      if (pts_status[ipt] == PDM_POLYGON_INSIDE) continue;
      double u, v, g;

      int converged = _projection_on_triangle_BtoA (pts_coord + 3*ipt,
                                                    vtx_coord + i0,
                                                    vtx_coord + i1,
                                                    vtx_coord + i2,
                                                    vtx_normal + i0,
                                                    vtx_normal + i1,
                                                    vtx_normal + i2,
                                                    &u,
                                                    &v,
                                                    &g,
                                                    verbose);
      if (verbose) printf("tri %d, pt %d: converged? %d, u = %f, v = %f\n", itri, ipt, converged, u, v);

      if (converged) {
        double w = 1. - u - v;

        if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && w >= 0) {
          pts_status[ipt] = PDM_POLYGON_INSIDE;

          for (int i = 0; i < 3; i++) {
            proj_coord[3*ipt+i] = u*vtx_coord[i0+i] + v*vtx_coord[i1+i] + w*vtx_coord[i2+i];
          }

        }
      }
    }// End of loop on points

  } // End of loop on sub-triangles


  state = PDM_triangulate_state_destroy (state);
  free (tri_vtx);
}


static void
_projection_on_polygon_AtoB
(
 const int             n_pts,
 const double         *pts_coord,
 const double         *pts_normal,
 const int             n_vtx,
 const double         *vtx_coord,
 double               *proj_coord,
 PDM_polygon_status_t *pts_status
 )
{
  PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx);
  PDM_l_num_t *tri_vtx = malloc (sizeof(PDM_l_num_t) * (n_vtx - 2)*3);

  int n_tri = PDM_triangulate_polygon (3,
                                       n_vtx,
                                       vtx_coord,
                                       NULL,
                                       NULL,
                                       PDM_TRIANGULATE_ELT_DEF,
                                       tri_vtx,
                                       state);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    pts_status[ipt] = PDM_POLYGON_OUTSIDE;
  }

  /* Loop on sub-triangles */
  for (int itri = 0; itri < n_tri; itri++) {

    int i0 = 3 * (tri_vtx[3*itri]     - 1);
    int i1 = 3 * (tri_vtx[3*itri  +1] - 1);
    int i2 = 3 * (tri_vtx[3*itri + 2] - 1);

    /* Loop on points */
    for (int ipt = 0; ipt < n_pts; ipt++) {

      if (pts_status[ipt] == PDM_POLYGON_INSIDE) continue;
      double u, v, g;

      int converged = _projection_on_triangle_AtoB (pts_coord + 3*ipt,
                                                    pts_normal + 3*ipt,
                                                    vtx_coord + i0,
                                                    vtx_coord + i1,
                                                    vtx_coord + i2,
                                                    &u,
                                                    &v,
                                                    &g);

      if (converged) {
        double w = 1. - u - v;

        if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && w >= 0) {
          pts_status[ipt] = PDM_POLYGON_INSIDE;

          for (int i = 0; i < 3; i++) {
            proj_coord[3*ipt+i] = u*vtx_coord[i0+i] + v*vtx_coord[i1+i] + w*vtx_coord[i2+i];
          }

        }
      }
    }// End of loop on points

  } // End of loop on sub-triangles


  state = PDM_triangulate_state_destroy (state);
  free (tri_vtx);
}











static void
_read_polygonal_mesh
(
 const char  *filename,
 int         *n_face,
 int         *n_vtx,
 double     **vtx_coord,
 int        **face_vtx_idx,
 int        **face_vtx
 )
{
  FILE *f = fopen(filename, "r");

  if (f == NULL) {
    printf("Could not read %s\n", filename);
    return;
  }

  fscanf(f, "%d %d", n_vtx, n_face);

  // Read vtx_coord
  *vtx_coord = malloc (sizeof(double) * (*n_vtx) * 3);
  for (int i = 0; i < *n_vtx; i++) {
    fscanf(f, "%lf %lf %lf", *vtx_coord + 3*i, *vtx_coord + 3*i+1, *vtx_coord + 3*i+2);
  }

  // Read face_vtx_idx
  *face_vtx_idx = malloc (sizeof(int) * (*n_face + 1));
  for (int i = 0; i <= *n_face; i++) {
    fscanf(f, "%d", *face_vtx_idx + i);
  }

  // Read face_vtx
  *face_vtx = malloc (sizeof(int) * (*face_vtx_idx)[*n_face]);
  for (int i = 0; i < (*face_vtx_idx)[*n_face]; i++) {
    fscanf(f, "%d", *face_vtx + i);
  }

  fclose(f);
}

static void
_read_edge
(
 const char *filename,
 double      x0[3],
 double      x1[3]
 )
{
  FILE *f = fopen(filename, "r");

  if (f == NULL) {
    printf("Could not read %s\n", filename);
    return;
  }

  fscanf(f, "%lf %lf %lf", x0, x0+1, x0+2);
  fscanf(f, "%lf %lf %lf", x1, x1+1, x1+2);

  fclose(f);
}




static void
_oriented_planar_grid
(
 const int      n_vtx,
 const double  *vtx_coord,
 const int      n_pts_seg,
 double       **pts_coord,
 int           *n_pts
 )
{
  /* Barycenter */
  double barycenter[3];
  PDM_plane_barycenter (n_vtx,
                        vtx_coord,
                        barycenter);

  /* Normal */
  double normal[3];
  PDM_plane_normal (n_vtx,
                    vtx_coord,
                    normal);

  /* Orthgonal frame for median plane */
  double max_l = 0.;
  double vec1[3];
  for (int i = 0; i < n_vtx; i++) {
    double v[3] = {vtx_coord[3*i]   - barycenter[0],
                   vtx_coord[3*i+1] - barycenter[1],
                   vtx_coord[3*i+2] - barycenter[2]};
    double dot = PDM_DOT_PRODUCT (v, normal);
    v[0] -= dot * normal[0];
    v[1] -= dot * normal[1];
    v[2] -= dot * normal[2];

    double l = PDM_MODULE (v);
    if (l > max_l) {
      vec1[0] = v[0] / l;
      vec1[1] = v[1] / l;
      vec1[2] = v[2] / l;
    }
  }

  double vec2[3];
  PDM_CROSS_PRODUCT (vec2, normal, vec1);


  /* Bounds in median plane frame */
  double umin = DBL_MAX;
  double umax = -umin;
  double vmin = umin;
  double vmax = -vmin;
  for (int i = 0; i < n_vtx; i++) {
    double p[3] = {vtx_coord[3*i]   - barycenter[0],
                   vtx_coord[3*i+1] - barycenter[1],
                   vtx_coord[3*i+2] - barycenter[2]};
    double u = PDM_DOT_PRODUCT (vec1, p);
    double v = PDM_DOT_PRODUCT (vec2, p);

    umin = PDM_MIN (umin, u);
    umax = PDM_MAX (umax, u);
    vmin = PDM_MIN (vmin, v);
    vmax = PDM_MAX (vmax, v);
  }

  double mrg = 1e-2;
  double umrg = mrg * (umax - umin);
  double vmrg = mrg * (vmax - vmin);

  umin -= umrg;
  umax += umrg;
  vmin -= vmrg;
  vmax += vmrg;

  /* Grid */
  *n_pts = n_pts_seg * n_pts_seg;
  *pts_coord = malloc (sizeof(double) * (*n_pts) * 3);

  double step_u = (umax - umin) / (n_pts_seg - 1);
  double step_v = (vmax - vmin) / (n_pts_seg - 1);

  double *_p = *pts_coord;
  for (int j = 0; j < n_pts_seg; j++) {
    double v = vmin + j * step_v;

    for (int i = 0; i < n_pts_seg; i++) {
      double u = umin + i * step_u;

      for (int k = 0; k < 3; k++) {
        _p[k] = barycenter[k] + u * vec1[k] + v * vec2[k];
      }
      _p += 3;

    }
  }
}



static void
_write_vtk_polydata
(
 const char   *filename,
 const int     n_vtx,
 const double *vtx_coord,
 const double *vtx_normal,
 const int     n_face,
 const int    *face_vtx_idx,
 const int    *face_vtx,
 const int    *face_num
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "polydata\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  const double *_pt = vtx_coord;
  for (int i = 0; i < n_vtx; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d ", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, "%d ", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }

  if (vtx_normal != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "VECTORS vtx_normal double\n");
    _pt = vtx_normal;
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
      _pt += 3;
    }
  }

  fprintf(f, "CELL_DATA %d\n", n_face);
  fprintf(f, "SCALARS lnum int 1\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  if (face_num == NULL) {
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d\n", i+1);
    }
  }
  else {
    for (int i = 0; i < n_face; i++) {
      fprintf(f, "%d\n", face_num[i]);
    }
  }

  fclose(f);
}


static void
_write_vtk_structured_grid
(
 const char   *filename,
 const int     n_i,
 const int     n_j,
 const int     n_pts,
 const double *pts_coord,
 const int    *pts_scalar
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "polydata\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET STRUCTURED_GRID\n");

  fprintf(f, "DIMENSIONS %d %d %d\n", n_i, n_j, 1);
  fprintf(f, "POINTS %d double\n", n_pts);
  const double *_pt = pts_coord;
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }

  if (pts_scalar != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_pts);
    fprintf(f, "SCALARS lnum int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, "%d\n", pts_scalar[i]);
    }
  }

  fclose(f);
}






int main(int argc, char *argv[])
{

  char *filename = NULL;
  int n_pts_seg = 10;

  _read_args (argc,
              argv,
              &filename,
              &n_pts_seg);

  int n_face, n_vtx;
  double *vtx_coord = NULL;
  int *face_vtx_idx = NULL;
  int *face_vtx = NULL;

  _read_polygonal_mesh (filename,
                        &n_face,
                        &n_vtx,
                        &vtx_coord,
                        &face_vtx_idx,
                        &face_vtx);

  double *vtx_normal = NULL;
  int singular = _compute_normals (n_face,
                                   face_vtx_idx,
                                   face_vtx,
                                   n_vtx,
                                   vtx_coord,
                                   &vtx_normal);
  if (singular) {
    PDM_error (__FILE__, __LINE__, 0,"Error singular vertex normals\n");
  }


  _write_vtk_polydata ("polygonal_mesh.vtk",
                       n_vtx,
                       vtx_coord,
                       vtx_normal,
                       n_face,
                       face_vtx_idx,
                       face_vtx,
                       NULL);




  int n_tri = 0;
  int max_n_vtx_poly = 0;
  for (int ipoly = 0; ipoly < n_face; ipoly++) {
    int n_vtx_poly = face_vtx_idx[ipoly+1] - face_vtx_idx[ipoly];
    n_tri += (n_vtx_poly - 2)*3;
    max_n_vtx_poly = PDM_MAX (max_n_vtx_poly, n_vtx_poly);
  }

  PDM_triangulate_state_t *state = PDM_triangulate_state_create (max_n_vtx_poly);

  int *tri_vtx_idx = malloc (sizeof(int) * (n_tri+1));
  int *tri_vtx = malloc (sizeof(int) * n_tri * 3);
  int *tri_poly = malloc (sizeof(int) * n_tri);
  tri_vtx_idx[0] = 0;
  n_tri = 0;
  for (int ipoly = 0; ipoly < n_face; ipoly++) {
    int n_vtx_poly = face_vtx_idx[ipoly+1] - face_vtx_idx[ipoly];

    int _n_tri = PDM_triangulate_polygon (3,
                                          n_vtx_poly,
                                          vtx_coord,
                                          NULL,
                                          face_vtx + face_vtx_idx[ipoly],
                                          PDM_TRIANGULATE_MESH_DEF,
                                          tri_vtx + tri_vtx_idx[n_tri],
                                          state);

    for (int i = 0; i < _n_tri; i++) {
      tri_vtx_idx[n_tri+1] = tri_vtx_idx[n_tri] + 3;
      tri_poly[n_tri] = ipoly + 1;
      n_tri++;
    }
  }

  _write_vtk_polydata ("polygonal_mesh_triangulation.vtk",
                       n_vtx,
                       vtx_coord,
                       vtx_normal,
                       n_tri,
                       tri_vtx_idx,
                       tri_vtx,
                       tri_poly);





  int n_pts;
  double *pts_coord = NULL;
  _oriented_planar_grid (n_vtx,
                         vtx_coord,
                         n_pts_seg,
                         &pts_coord,
                         &n_pts);

  double *pts_proj = malloc (sizeof(double) * n_pts * 3);
  int    *pts_poly = malloc (sizeof(int) * n_pts);

  for (int ipt = 0; ipt < n_pts; ipt++) {
    pts_poly[ipt] = -1;
  }

  double *poly_coord = malloc (sizeof(double) * max_n_vtx_poly * 3);
  double *poly_normal = malloc (sizeof(double) * max_n_vtx_poly * 3);

  PDM_polygon_status_t *pts_status = malloc (sizeof(PDM_polygon_status_t) * n_pts);

  for (int ipoly = 3; ipoly < 4; ipoly++) {//n_face; ipoly++) {
    int n_vtx_poly = face_vtx_idx[ipoly+1] - face_vtx_idx[ipoly];
    int *poly_vtx = face_vtx + face_vtx_idx[ipoly];

    for (int i = 0; i < n_vtx_poly; i++) {
      int ivtx = poly_vtx[i] - 1;
      for (int j = 0; j < 3; j++) {
        poly_coord[3*i+j] = vtx_coord[3*ivtx+j];
        poly_normal[3*i+j] = vtx_normal[3*ivtx+j];
      }
    }

    _projection_on_polygon_BtoA (n_pts,
                                 pts_coord,
                                 n_vtx_poly,
                                 poly_coord,
                                 poly_normal,
                                 pts_proj,
                                 pts_status,
                                 0);

    for (int ipt = 0; ipt < n_pts; ipt++) {
      if (pts_poly[ipt] >= 0) continue;

      if (pts_status[ipt] == PDM_POLYGON_INSIDE) {
        pts_poly[ipt] = ipoly + 1;
      }
    }

  }

  for (int ipt = 0; ipt < n_pts; ipt++) {
    if (pts_poly[ipt] <= 0) {
      for (int j = 0; j < 3; j++) {
        pts_proj[3*ipt+j] = pts_coord[3*ipt+j];
      }
    }
  }




  _write_vtk_structured_grid ("polygonal_projection.vtk",
                              n_pts_seg,
                              n_pts_seg,
                              n_pts,
                              pts_coord,
                              pts_poly);


  _write_vtk_structured_grid ("polygonal_projection_coord.vtk",
                              n_pts_seg,
                              n_pts_seg,
                              n_pts,
                              pts_proj,
                              pts_poly);


  FILE *f = fopen("polygonal_projection_lines.vtk", "w");
  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "proj_lines\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 2*n_pts);
  double *_pt = pts_coord;
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }
  _pt = pts_proj;
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 3*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "2 %d %d\n", i, n_pts+i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "3\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_pts);
  fprintf(f, "SCALARS lnum int 1\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "%d\n", pts_poly[i]);
  }

  fclose(f);










  /*****************************************************/

  double b0[3], b1[3];
  _read_edge ("edge.dat",
              b0,
              b1);



int ifaceA = 3;
  int iA0 = face_vtx[face_vtx_idx[ifaceA]] - 1;
  int iA1 = face_vtx[face_vtx_idx[ifaceA]+1] - 1;

  double *a0 = vtx_coord + 3*iA0;
  double *a1 = vtx_coord + 3*iA1;
  double *n0 = vtx_normal + 3*iA0;
  double *n1 = vtx_normal + 3*iA1;

  double uA, uB;
  PDM_line_intersect_t tIntersect = _line_intersection_projection (a0,
                                                                   a1,
                                                                   n0,
                                                                   n1,
                                                                   b0,
                                                                   b1,
                                                                   &uA,
                                                                   &uB);

  printf("tIntersect = %d, uA = %f, uB = %f\n", (int) tIntersect, uA, uB);

  double xA[3] = {(1 - uA)*a0[0] + uA*a1[0],
                  (1 - uA)*a0[1] + uA*a1[1],
                  (1 - uA)*a0[2] + uA*a1[2]};

  double xB[3] = {(1 - uB)*b0[0] + uB*b1[0],
                  (1 - uB)*b0[1] + uB*b1[1],
                  (1 - uB)*b0[2] + uB*b1[2]};


  f = fopen("edge_intersection.vtk", "w");
  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "intersection\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS 2 double\n");
  fprintf(f, "%lf %lf %lf\n", xA[0], xA[1], xA[2]);
  fprintf(f, "%lf %lf %lf\n", xB[0], xB[1], xB[2]);

  fprintf(f, "CELLS 2 4\n1 0\n1 1\n");

  fprintf(f, "CELL_TYPES 2\n1\n1\n");

  fprintf(f, "POINT_DATA 2\n");
  fprintf(f, "SCALARS imesh int 1\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  fprintf(f, "1\n2\n");

  fclose(f);








  const int n_sample = 1000;
  double *t_edge = malloc (sizeof(double) * n_sample);
  double *edgeB_coord = malloc (sizeof(double) * n_sample * 3);
  double *proj_edgeB_coord = malloc (sizeof(double) * n_sample * 3);
  int    *proj_edgeB_poly = malloc (sizeof(int) * n_sample);

  for (int i = 0; i < n_sample; i++) {
    proj_edgeB_poly[i] = -1;
    t_edge[i] = (double) i / (double) (n_sample - 1);
  }

  int l = 0;
  int r = n_sample;
  while (l < r - 1) {
    int m = l + (r - l) / 2;

    if (t_edge[m] < uB) {
      l = m;
    } else {
      r = m;
    }
  }

  t_edge[l] = uB;

  for (int i = 0; i < n_sample; i++) {
    double t = t_edge[i];
    double s = 1. - t;

    for (int j = 0; j < 3; j++) {
      edgeB_coord[3*i+j] = s*b0[j] + t*b1[j];
    }
  }



  PDM_polygon_status_t *proj_edgeB_status = malloc (sizeof(PDM_polygon_status_t) * n_sample);

  for (int ipoly = 0; ipoly < n_face; ipoly++) {
    int n_vtx_poly = face_vtx_idx[ipoly+1] - face_vtx_idx[ipoly];
    int *poly_vtx = face_vtx + face_vtx_idx[ipoly];

    for (int i = 0; i < n_vtx_poly; i++) {
      int ivtx = poly_vtx[i] - 1;
      for (int j = 0; j < 3; j++) {
        poly_coord[3*i+j] = vtx_coord[3*ivtx+j];
        poly_normal[3*i+j] = vtx_normal[3*ivtx+j];
      }
    }

    //printf("\npolygon %d\n", ipoly);
    _projection_on_polygon_BtoA (n_sample,
                                 edgeB_coord,
                                 n_vtx_poly,
                                 poly_coord,
                                 poly_normal,
                                 proj_edgeB_coord,
                                 proj_edgeB_status,
                                 0);

    for (int ipt = 0; ipt < n_sample; ipt++) {
      if (proj_edgeB_poly[ipt] >= 0) continue;

      if (proj_edgeB_status[ipt] == PDM_POLYGON_INSIDE) {
        proj_edgeB_poly[ipt] = ipoly + 1;
      }
    }

  }

  for (int ipt = 0; ipt < n_sample; ipt++) {
    if (proj_edgeB_poly[ipt] <= 0) {
      for (int j = 0; j < 3; j++) {
        proj_edgeB_coord[3*ipt+j] = edgeB_coord[3*ipt+j];
      }
    }
  }


  f = fopen("edgeB.vtk", "w");
  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "edgeB\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_sample);
  _pt = edgeB_coord;
  for (int i = 0; i < n_sample; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }

  fprintf(f, "CELLS %d %d\n", n_sample-1, 3*(n_sample - 1));
  for (int i = 0; i < n_sample-1; i++) {
    fprintf(f, "2 %d %d\n", i, i+1);
  }

  fprintf(f, "CELL_TYPES %d\n", n_sample-1);
  for (int i = 0; i < n_sample-1; i++) {
    fprintf(f, "3\n");
  }

  fclose(f);




  f = fopen("projection_edgeB.vtk", "w");
  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "proj_edgeB\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_sample);
  _pt = proj_edgeB_coord;
  for (int i = 0; i < n_sample; i++) {
    fprintf(f, "%lf %lf %lf\n", _pt[0], _pt[1], _pt[2]);
    _pt += 3;
  }

  fprintf(f, "CELLS %d %d\n", n_sample-1, 3*(n_sample - 1));
  for (int i = 0; i < n_sample-1; i++) {
    fprintf(f, "2 %d %d\n", i, i+1);
  }

  fprintf(f, "CELL_TYPES %d\n", n_sample-1);
  for (int i = 0; i < n_sample-1; i++) {
    fprintf(f, "3\n");
  }

  fprintf(f, "POINT_DATA %d\n", n_sample);
  fprintf(f, "SCALARS lnum int 1\n");
  fprintf(f, "LOOKUP_TABLE default\n");
  for (int i = 0; i < n_sample; i++) {
    fprintf(f, "%d\n", proj_edgeB_poly[i]);
  }

  fclose(f);






  return 0;
}
