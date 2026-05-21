/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_intersection_surf_surf_atomic.h"
#include "pdm.h"
#include "pdm_error.h"
#include "pdm_line.h"
#include "pdm_logging.h"
#include "pdm_priv.h"

/*============================================================================
 * Global variables
 *============================================================================*/

// static int dbg_enabled = 0;

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static int
_solve_quadratic
(
 const double a,
 const double b,
 const double c,
       double solutions[2]
 )
{
  if (a == 0) {
    // Linear equation
    if (b == 0) {
      if (c == 0) {
        // trivial equation 0 = 0
        return -1;
      }
      else {
        // no solution
        return 0;
      }
    }
    else {
      solutions[0] = -c/b;
      return 1;
    }
  }

  else {
    // True quadratic equation
    double ia = 1./a;
    double _b = b*ia;
    double _c = c*ia;

    if (_c == 0) {
      solutions[0] = PDM_MIN(0, -_b);
      solutions[1] = PDM_MAX(0, -_b);
      return 2;
    }
    else {
      double d = _b*_b - 4.*_c;
      if (d < 0) {
        // no real solution
        return 0;
      }
      else {
        // two real solutions (possibly one double solution)
        d = sqrt(d);
        double x1, x2;
        if (_b < 0) {
          x1 = 0.5*(-_b + d);
        }
        else {
          x1 = 0.5*(-_b - d);
        }
        x2 = _c/x1;
        solutions[0] = PDM_MIN(x1, x2);
        solutions[1] = PDM_MAX(x1, x2);
        return 2;
      }
    }
  }
}
PDM_GCC_SUPPRESS_WARNING_POP

static PDM_line_intersect_t
_line_intersection_projection
(
 const double  a0[3],
 const double  a1[3],
 const double  b0[3],
 const double  b1[3],
 const double  d0[3],
 const double  d1[3],
       double *ta,
       double *tb
)
{
  double a0a1[3] = {a1[0] - a0[0],
                    a1[1] - a0[1],
                    a1[2] - a0[2]};

  double b0b1[3] = {b1[0] - b0[0],
                    b1[1] - b0[1],
                    b1[2] - b0[2]};

  double d0d1[3] = {d1[0] - d0[0],
                    d1[1] - d0[1],
                    d1[2] - d0[2]};

  double a0b0[3] = {b0[0] - a0[0],
                    b0[1] - a0[1],
                    b0[2] - a0[2]};

  double a0a1_x_d0d1[3];
  PDM_CROSS_PRODUCT(a0a1_x_d0d1, a0a1, d0d1);

  double a0a1_x_d0[3];
  PDM_CROSS_PRODUCT(a0a1_x_d0, a0a1, d0);

  double c2 = PDM_DOT_PRODUCT(a0a1_x_d0d1, b0b1);
  double c1 = PDM_DOT_PRODUCT(a0a1_x_d0, b0b1) + PDM_DOT_PRODUCT(a0a1_x_d0d1, a0b0);
  double c0 = PDM_DOT_PRODUCT(a0a1_x_d0, a0b0);


  double solutions[2];
  int n_solutions = _solve_quadratic(c2, c1, c0, solutions);


  if (n_solutions < 0) {
    return PDM_LINE_INTERSECT_ON_LINE;
  }
  else if (n_solutions == 0) {
    return PDM_LINE_INTERSECT_NO;
  }
  else if (n_solutions == 1) {
    *tb = solutions[0];
    if (*tb < 0 || *tb > 1) {
      return PDM_LINE_INTERSECT_NO;
    }
  }
  else {
    if (solutions[0] >= 0 && solutions[0] <= 1) {
      if (solutions[1] >= 0 && solutions[0] <= 1) {
        // two valid solutions, which one do we choose??
        return PDM_LINE_INTERSECT_UNDEF;
      }
      else {
        *tb = solutions[0];
      }
    }
    else if (solutions[1] >= 0 && solutions[0] <= 1) {
      *tb = solutions[1];
    }
    else {
      return PDM_LINE_INTERSECT_NO;
    }
  }



  double d[3];
  for (int i = 0; i < 3; i++) {
    d[i] = (1-(*tb))*d0[i] + (*tb)*d1[i];
  }

  double l[3];
  PDM_CROSS_PRODUCT(l, b0b1, d);

  double denom = PDM_DOT_PRODUCT(l, a0a1);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (denom == 0) {
    return PDM_LINE_INTERSECT_UNDEF;
  }
  else {
    *ta = PDM_DOT_PRODUCT(l, a0b0) / denom;
    if (*ta >= 0 && *ta <= 1) {
      return PDM_LINE_INTERSECT_YES;
    }
    else {
      return PDM_LINE_INTERSECT_NO;
    }
  }
  PDM_GCC_SUPPRESS_WARNING_POP

}







static double
_line_signed_distance
(
 const double coord[2],
 const int    iline
 )
{
  switch (iline) {
  case 0: // x = 0
    return coord[0];
  case 1: // x = 1
    return 1. - coord[0];
  case 2: // y = 0
    return coord[1];
  case 3: // x + y = 1
    return 1. - coord[0] - coord[1];
  default:
    log_error("_line_signed_distance: wrong line number %d\n", iline);
  }
  return 0;
}


static void
_line_vtxA_id
(
 const int  iline,
       int *vtxA_id0,
       int *vtxA_id1
 )
{
  switch (iline) {
  case 0: // x = 0
    *vtxA_id0 = 2;
    *vtxA_id1 = 0;
    break;
  case 1: // x = 1
    *vtxA_id0 = -1;
    *vtxA_id1 = -1;
    break;
  case 2: // y = 0
    *vtxA_id0 = 0;
    *vtxA_id1 = 1;
    break;
  case 3: // x + y = 1
    *vtxA_id0 = 1;
    *vtxA_id1 = 2;
    break;
  default:
    log_error("_line_vtxA_id: wrong line number %d\n", iline);
  }
}



static inline void
_project_on_line3
(
       double *coord,
 const int     n
 )
{
  for (int i = 0; i < n; i++) {
    coord[2*i+1] = 1 - coord[2*i];
  }
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static int
_grandy2d
(
 double *coord,
 double  triaA_coord[9],
 double  edgeB_coord[6],
 double  edgeB_normal[6]
 )
{
  int dbg_enabled = 1;
  PDM_UNUSED(coord);
  PDM_UNUSED(triaA_coord);
  PDM_UNUSED(edgeB_coord);
  PDM_UNUSED(edgeB_normal);

  double uvA[6] = {
    0, 0,
    1, 0,
    0, 1
  };

  double tmin = 0;
  double tmax = 1;

  /* Check if the initial segment is inside the unit triangle */
  int inside = 1;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (coord[2*i+j] < 0 || coord[2*i+j] > 1) {
        inside = 0;
        break;
      }
    }
    if (!inside) {
      break;
    }

    if (1 - coord[2*i] - coord[2*i+1] < 0) {
        inside = 0;
        break;
      }
  }

  if (inside) {
    if (dbg_enabled) {
      log_trace("initial segment in triangle\n");
    }
    return 2;
  }


  /* Intersect/clip with all 4 lines */
  for (int iline = 0; iline < 4; iline++) {

    double f0 = _line_signed_distance(&coord[0], iline);
    double f1 = _line_signed_distance(&coord[2], iline);


    if (f0 == 0) {

      if (f1 > 0) {
        // 0 on, rest in => all in
        continue;
      }
      else {
        // 0 on, rest on or out
        if (iline == 3) {
          _project_on_line3(&coord[2], 1);
          return 2;
        }
        else {
          return 0;
        }
      }

    } // end if f0 == 0

    else if (f1 == 0) {

      if (f0 < 0) {
        // 1 on, rest out => all out
        if (iline == 3) {
          _project_on_line3(&coord[0], 1);
          return 2;
        }
        else {
          return 0;
        }
      }
      else {
        // 1 on, rest in => all in
        continue;
      }

    }

    else if (f0*f1 < 0) {
      // intersection
      int vtxA_id0 = -1;
      int vtxA_id1 = -1;
      _line_vtxA_id(iline, &vtxA_id0, &vtxA_id1);

      double tA, tB;
      double inter_coord[2];
      PDM_line_intersect_t stat = PDM_LINE_INTERSECT_UNDEF;
      if (edgeB_normal != NULL && vtxA_id0 >= 0 && vtxA_id1 >= 0) {
        stat = _line_intersection_projection(&triaA_coord[3*vtxA_id0],
                                             &triaA_coord[3*vtxA_id1],
                                             &edgeB_coord[0],
                                             &edgeB_coord[3],
                                             &edgeB_normal[0],
                                             &edgeB_normal[3],
                                             &tA,
                                             &tB);
        if (stat == PDM_LINE_INTERSECT_YES) {
          // intersection coord from tA or tB?????
          assert(tB >= tmin);
          assert(tB <= tmax);
          for (int i = 0; i < 2; i++) {
            inter_coord[i] = (1-tA)*uvA[2*vtxA_id0+i] + tA*uvA[2*vtxA_id1+i];
          }
        }
      }

      if (stat != PDM_LINE_INTERSECT_YES) {
        tB = f0 / (f0 - f1);
        for (int i = 0; i < 2; i++) {
          inter_coord[i] = (1-tB)*coord[i] + tB*coord[2+i];
        }

        tB = tmin + (tmax-tmin)*tB;
      }


      if (iline == 3) {
        memcpy(&coord[4], &coord[2],   sizeof(double)*2);
        memcpy(&coord[2], inter_coord, sizeof(double)*2);
        if (f0 < 0) {
          _project_on_line3(&coord[0], 1);
        }
        else {
          _project_on_line3(&coord[4], 1);
        }
        return 3;
      }
      else {
        if (f0 < 0) {
          tmin = tB;
          memcpy(&coord[0], inter_coord, sizeof(double)*2);
        }
        else {
          tmax = tB;
          memcpy(&coord[2], inter_coord, sizeof(double)*2);
        }
      }

    }

    else if (f0 < 0) {
      // no intersection, all out
      if (iline == 3) {
        _project_on_line3(coord, 2);
        return 2;
      }
      else {
        return 0;
      }
    }

    else {
      // no intersection, all in
      continue;
    }

  } // End of loop on lines

  return 0;
}
PDM_GCC_SUPPRESS_WARNING_POP


static inline void
_get_edge_vtx
(
  PDM_mesh_intersection_surf_surf_polygon_t  *poly,
  int                                         i,
  int                                        *i_vtx0,
  int                                        *i_vtx1,
  double                                    **p0,
  double                                    **p1
)
{
  if (poly->face_vtx != NULL) {
    *i_vtx0 = poly->face_vtx[ i                ] - 1;
    *i_vtx1 = poly->face_vtx[(i+1)%poly->n_edge] - 1;
  }
  else {
    int i_edge = poly->face_edge[i];
    if (i_edge < 0) {
      i_edge = -i_edge - 1;
      *i_vtx0 = poly->edge_vtx[2*i_edge+1] - 1;
      *i_vtx1 = poly->edge_vtx[2*i_edge  ] - 1;
    }
    else {
      i_edge = i_edge - 1;
      *i_vtx0 = poly->edge_vtx[2*i_edge  ] - 1;
      *i_vtx1 = poly->edge_vtx[2*i_edge+1] - 1;
    }
  }

  *p0 = &poly->coord[3*(*i_vtx0)];
  *p1 = &poly->coord[3*(*i_vtx1)];
}


static int
_clip_segment
(
  double *uv0,
  double *uv1,
  double *uv2
)
{
  double f0, f1, t, s;

  int n_comp = 0;

  // Clip segment to keep intersection with the positive quadrant (u >= 0, v >= 0)
  for (int i = 0; i < 2; i++) {

    f0 = uv0[i];
    f1 = uv1[i];

    if (f0*f1 < 0) {
      t = f0/(f0 - f1);

      if (f0 < 0) {
        // keep segment [t, 1]
        uv0[0] = (1-t)*uv0[0] + t*uv1[0];
        uv0[1] = (1-t)*uv0[1] + t*uv1[1];
      }
      else { // f1 < 0
        // keep segment [0, t]
        uv1[0] = (1-t)*uv0[0] + t*uv1[0];
        uv1[1] = (1-t)*uv0[1] + t*uv1[1];
      }
    }
    else if (f0 < 0) { // f1 < 0 as well
      // The segment does not intersect the positive quadrant
      return n_comp;
    }
  }

  // Find intersection with diagonal (u + v = 1)
  f0 = 1 - uv0[0] - uv0[1];
  f1 = 1 - uv1[0] - uv1[1];
  if (f0*f1 < 0) {
    n_comp = 2;
    t = f0/(f0 - f1);
    uv2[0] = (1-t)*uv0[0] + t*uv1[0];
    uv2[1] = (1-t)*uv0[1] + t*uv1[1];
  }
  else {
    // no intersection
    n_comp = 1;
  }

  // Project outer part on diagonal
  if (f0 < 0) {
    s = 1./(1 - f0);
    uv0[0] *= s;
    uv0[1] *= s;
  }
  if (f1 < 0) {
    s = 1./(1 - f1);
    uv1[0] *= s;
    uv1[1] *= s;
  }

  return n_comp;
}


#define ONE_THRID 0.3333333333333333
static inline void
_geom_comp
(
  const double *uv0,
  const double *uv1,
        double *area,
        double *center
)
{
  *area = uv0[0]*uv1[1] - uv0[1]*uv1[0];
  center[0] = ONE_THRID*(uv0[0] + uv1[0]);
  center[1] = ONE_THRID*(uv0[1] + uv1[1]);
}
#undef ONE_THIRD


static inline void
_vector_ab
(
        double ab[3],
  const double a[3],
  const double b[3]
)
{
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

double PDM_mesh_intersection_surf_surf_atomic_compute
(
 double triaA_coord[9],
 double edgeB_coord[6],
 double edgeB_normal[6]
 )
{
  int dbg_enabled = 0;

  double mat[3][3];
  double rhs[3] = {0, 0, 0};

  for (int i = 0; i < 3; i++) {
    mat[i][0] = triaA_coord[3+i] - triaA_coord[i];
    mat[i][1] = triaA_coord[6+i] - triaA_coord[i];
  }

  double uv[6];
  for (int i = 0; i < 2; i++) {

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
    if (edgeB_normal == NULL) {
      // projection = barycentric coordinates in A
      double m00 = 0;
      double m01 = 0;
      double m11 = 0;
      for (int j = 0; j < 3; j++) {
        m00 += mat[j][0] * mat[j][0];
        m01 += mat[j][0] * mat[j][1];
        m11 += mat[j][1] * mat[j][1];
        rhs[0] += mat[j][0] * (edgeB_coord[3*i+j] - triaA_coord[i]);
        rhs[1] += mat[j][1] * (edgeB_coord[3*i+j] - triaA_coord[i]);
      }

      double det = m00*m11 - m01*m01;
PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
      if (det == 0) {
        PDM_error(__FILE__, __LINE__, 0,
                  "cannot project (degenerate triangle A)\n");
      }
PDM_GCC_SUPPRESS_WARNING_POP

      double idet = 1./det;

      uv[2*i  ] = (rhs[0]*m11 - rhs[1]*m01) * idet;
      uv[2*i+1] = (rhs[1]*m00 - rhs[0]*m01) * idet;
    }

    else {

      for (int j = 0; j < 3; j++) {
        mat[j][2] = -edgeB_normal[3*i+j];
        rhs[j] = edgeB_coord[3*i+j] - triaA_coord[i];
      }

      double det = mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
      -            mat[1][0]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
      +            mat[2][0]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]);

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
      if (det == 0) {
        PDM_error(__FILE__, __LINE__, 0,
                  "cannot project : mat =\n"
                  "%f %f %f\n"
                  "%f %f %f\n"
                  "%f %f %f\n",
                  mat[0][0], mat[0][1], mat[0][2],
                  mat[1][0], mat[1][1], mat[1][2],
                  mat[2][0], mat[2][1], mat[2][2]);
      }
PDM_GCC_SUPPRESS_WARNING_POP

      double idet = 1./det;

      uv[2*i    ] = idet *
      (  rhs[0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
       - rhs[1]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
       + rhs[2]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]));

      uv[2*i + 1] = idet *
      (  mat[0][0]*(rhs[1]*mat[2][2] - rhs[2]*mat[1][2])
       - mat[1][0]*(rhs[0]*mat[2][2] - rhs[2]*mat[0][2])
       + mat[2][0]*(rhs[0]*mat[1][2] - rhs[1]*mat[0][2]));
    }

  }

  /* Clipping */
  if (dbg_enabled) {
    log_trace("initial segment :\n");
    for (int i = 0; i < 2; i++) {
      log_trace("%20.16f %20.16f\n", uv[2*i], uv[2*i+1]);
    }
  }
  int n_vtx = _grandy2d(uv,
                        triaA_coord,
                        edgeB_coord,
                        edgeB_normal);

  if (dbg_enabled) {
    log_trace("clipped :\n");
    for (int i = 0; i < n_vtx; i++) {
      log_trace("%20.16f %20.16f\n", uv[2*i], uv[2*i+1]);
    }
  }

  /* Compute column area */
  double area = 0;
  for (int i = 0; i < n_vtx-1; i++) {
    area += 0.5 * (uv[2*i+1] + uv[2*i+3]) * (uv[2*i] - uv[2*i+2]);
  }

  return area;
}



void PDM_mesh_intersection_surf_surf_atomic_compute2 // temporary name
(
  PDM_mesh_intersection_surf_surf_polygon_t *poly_a,
  PDM_mesh_intersection_surf_surf_polygon_t *poly_b,
  double                                    *area_ab,
  double                                    *center_ab
)
{
  *area_ab = 0.;
  for (int i = 0; i < 3; i++) {
    center_ab[i] = 0.;
  }

  int i_vtx_ref = -1;

  double *pi, *pj;

  double ai[3], aj[3], bi[3], bj[3];
  double uvi[2], uvj[2], uvk[2];

  for (int ia = 0; ia < poly_a->n_edge; ia++) {

    int i_vtx_a, j_vtx_a;
    _get_edge_vtx(poly_a, ia, &i_vtx_a, &j_vtx_a, &pi, &pj);

    i_vtx_ref = (i_vtx_ref < 0) ? i_vtx_a : i_vtx_ref;

    if ((i_vtx_a == i_vtx_ref) || (j_vtx_a == i_vtx_ref)) {
      continue;
    }

    _vector_ab(ai, &poly_a->coord[3*i_vtx_ref], pi);
    _vector_ab(aj, &poly_a->coord[3*i_vtx_ref], pj);

    double aiai = PDM_DOT_PRODUCT(ai, ai);
    double aiaj = PDM_DOT_PRODUCT(ai, aj);
    double ajaj = PDM_DOT_PRODUCT(aj, aj);

    double det_aij = aiai*ajaj - aiaj*aiaj;

    if (det_aij <= 0) {
      // vertices i_vtx_ref, i_vtx_a and j_vtx_a are collinear => ignore this subtriangle
      continue;
    }

    double inv_det_aij = 1./det_aij;

    double normal_aij[3];
    PDM_CROSS_PRODUCT(normal_aij, ai, aj);

    double area_aij = 0.5*PDM_MODULE(normal_aij);

    for (int ib = 0; ib < poly_b->n_edge; ib++) {

      int i_vtx_b, j_vtx_b;
      _get_edge_vtx(poly_b, ib, &i_vtx_b, &j_vtx_b, &pi, &pj);

      _vector_ab(bi, &poly_a->coord[3*i_vtx_ref], pi);
      _vector_ab(bj, &poly_a->coord[3*i_vtx_ref], pj);


      double biai = PDM_DOT_PRODUCT(bi, ai);
      double biaj = PDM_DOT_PRODUCT(bi, aj);
      double bjai = PDM_DOT_PRODUCT(bj, ai);
      double bjaj = PDM_DOT_PRODUCT(bj, aj);

      uvi[0] = biai*ajaj - biaj*aiaj;
      uvj[0] = bjai*ajaj - bjaj*aiaj;

      if ((uvi[0] <= 0) && (uvj[0] <= 0)) {
        continue;
      }

      uvi[1] = aiai*biaj - aiaj*biai;
      uvj[1] = aiai*bjaj - aiaj*bjai;

      if ((uvi[1] <= 0) && (uvj[1] <= 0)) {
        continue;
      }

      uvi[0] *= inv_det_aij;
      uvi[1] *= inv_det_aij;
      uvj[0] *= inv_det_aij;
      uvj[1] *= inv_det_aij;


      int n_comp = _clip_segment(uvi, uvj, uvk);

      if (n_comp == 0) {
        continue;
      }

      double area_uv;
      double center_uv[2];
      if (n_comp == 1) {
        _geom_comp(uvi, uvj, &area_uv, center_uv);
      }
      else { // n_comp == 2
        double area_uv_ik, area_uv_kj;
        double center_uv_ik[2], center_uv_kj[2];
        _geom_comp(uvi, uvk, &area_uv_ik, center_uv_ik);
        _geom_comp(uvk, uvj, &area_uv_kj, center_uv_kj);
        area_uv = area_uv_ik + area_uv_kj;
PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
        if (area_uv != 0) {
          double inv_area_uv = 1./area_uv;
          for (int l = 0; l < 2; l++) {
            center_uv[l] = (area_uv_ik*center_uv_ik[l] + area_uv_kj*center_uv_kj[l]) * inv_area_uv;
          }
        }
PDM_GCC_SUPPRESS_WARNING_POP
      }

      area_uv *= area_aij;

      *area_ab += area_uv;
      for (int l = 0; l < 3; l++) {
        center_ab[l] += area_uv * (
          poly_a->coord[3*i_vtx_ref+l] +
          center_uv[0]*ai[l] +
          center_uv[1]*aj[l]
        );
      }

    } // End loop on edges B
  } // End loop on edges A

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (*area_ab != 0) {
    double inv_area_ab = 1./(*area_ab);
    for (int l = 0; l < 3; l++) {
      center_ab[l] *= inv_area_ab;
    }
  }
PDM_GCC_SUPPRESS_WARNING_POP
}


