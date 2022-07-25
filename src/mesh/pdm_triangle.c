/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_triangle.h"
#include "pdm_line.h"
#include "pdm_plane.h"
#include "pdm_logging.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/* Computes the determinant of a 3x3 matrix defined by its columns */

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

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Evaluates the position in a triangle
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  pts      Triangle vertices coordinates
 * \param [out] closest  Closest Point in Triangle or NULL
 * \param [out] min_dist2 Square of the distance
 * \param [out] weights  Vertices weights or NULL
 *
 * \return      \ref PDM_TRIANGLE_INSIDE or \ref PDM_TRIANGLE_OUTSIDE
 *              if the projected is in the triangle or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */


PDM_triangle_status_t
PDM_triangle_evaluate_position
(
 const double  x[3],
 const double  pts[9],
       double *closest_point,
       double *min_dist2,
       double *weights
)
{
  double *pt1, *pt2, *pt3;
  double n[3], fabsn;
  double rhs[2], c1[2], c2[2];
  double det;
  double maxComponent;
  int idx=0, indices[2];
  double dist_to_point, dist_to_line1, dist_to_line2;
  double *closest, closest_point1[3], closest_point2[3], cp[3];
  double pcoords[3];

  pcoords[2] = 0.0;

  double weights_local[3];
  double *_weights = weights_local;
  if (weights != NULL) {
    _weights = weights;
  }

  /*
   * Get normal for triangle, only the normal direction is needed, i.e. the
   * normal need not be normalized (unit length)
   */

  PDM_plane_normal (3, pts, n);

  pt1 = (double *) pts;
  pt2 = (double *) pts + 3;
  pt3 = (double *) pts + 6;

  /*
   * Project point to plane
   */

  PDM_plane_projection2 (x, pt1, n, cp);

  /*
   * Construct matrices.  Since we have over determined system, need to find
   * which 2 out of 3 equations to use to develop equations. (Any 2 should
   * work since we've projected point to plane.)
   */

  maxComponent = 0.0;
  for (int i = 0; i < 3; i++) {

    /*
     * Trying to avoid an expensive call to fabs()
     */

    if (n[i] < 0) {
      fabsn = -n[i];
    }
    else {
      fabsn = n[i];
    }
    if (fabsn > maxComponent) {
      maxComponent = fabsn;
      idx = i;
    }
  }

  for (int j = 0, i = 0; i < 3; i++) {
    if (i != idx) {
      indices[j++] = i;
    }
  }

  for (int i = 0; i < 2; i++) {
    rhs[i] = cp[indices[i]] - pt3[indices[i]];
    c1[i] = pt1[indices[i]] - pt3[indices[i]];
    c2[i] = pt2[indices[i]] - pt3[indices[i]];
  }

PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

  if ((det = PDM_DETERMINANT2X2(c1,c2)) == 0.0) {
    pcoords[0] = 0.0;
    pcoords[1] = 0.0;
    return PDM_TRIANGLE_DEGENERATED;
  }

PDM_GCC_SUPPRESS_WARNING_POP

  pcoords[0] = PDM_DETERMINANT2X2(rhs,c2) / det;
  pcoords[1] = PDM_DETERMINANT2X2(c1,rhs) / det;

  /*
   * Okay, now find closest point to element
   */

  _weights[0] = 1 - (pcoords[0] + pcoords[1]);
  _weights[1] = pcoords[0];
  _weights[2] = pcoords[1];

  if ( _weights[0] >= 0.0 && _weights[0] <= 1.0 &&
       _weights[1] >= 0.0 && _weights[1] <= 1.0 &&
       _weights[2] >= 0.0 && _weights[2] <= 1.0 ) {

    /*
     * Projection distance
     */

    if (closest_point) {
      double v_cp_x[3];
      for (int i = 0; i < 3; i++) {
        v_cp_x[i] = cp[i] - x[i];
      }

      *min_dist2 = PDM_DOT_PRODUCT(v_cp_x, v_cp_x);
      closest_point[0] = cp[0];
      closest_point[1] = cp[1];
      closest_point[2] = cp[2];
    }
    return PDM_TRIANGLE_INSIDE;
  }

  else {
    double t1, t2;
    if (closest_point) {
      if (_weights[1] < 0.0 && _weights[2] < 0.0) {
        double v_pt3_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt3_x[i] = pt3[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT (v_pt3_x, v_pt3_x);
        dist_to_line1 = PDM_line_distance (x, pt1, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt3, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt3;
          _weights[0] = 1.;
          _weights[1] = 0.;
          _weights[2] = 0.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = t1;
          _weights[1] = 1. - t1;
          _weights[2] = 0.;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
          if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 1. - t2;
          _weights[1] = 0.;
          _weights[2] = t2;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if (_weights[2] < 0.0 && _weights[0] < 0.0) {
        double v_pt1_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt1_x[i] = pt1[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT(v_pt1_x, v_pt1_x);
        dist_to_line1 = PDM_line_distance (x, pt1, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt1, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt1;
          _weights[0] = 0.;
          _weights[1] = 1.;
          _weights[2] = 0.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = t1;
          _weights[1] = 1. - t1;
          _weights[2] = 0.;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
           if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = 1. - t2;
          _weights[2] = t2;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if ( _weights[1] < 0.0 && _weights[0] < 0.0 ) {
        double v_pt2_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt2_x[i] = pt2[i] - x[i];
        }
        dist_to_point = PDM_DOT_PRODUCT (v_pt2_x, v_pt2_x);
        dist_to_line1 = PDM_line_distance (x, pt2, pt3, &t1, closest_point1);
        dist_to_line2 = PDM_line_distance (x, pt1, pt2, &t2, closest_point2);
        if (dist_to_point < dist_to_line1) {
          *min_dist2 = dist_to_point;
          closest = pt2;
          _weights[0] = 0.;
          _weights[1] = 0.;
          _weights[2] = 1.;
        }
        else {
          *min_dist2 = dist_to_line1;
          closest = closest_point1;
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = t1;
          _weights[1] = 0.;
          _weights[2] = 1. - t1;
        }
        if (dist_to_line2 < *min_dist2) {
          *min_dist2 = dist_to_line2;
          closest = closest_point2;
           if (t2 < 0.) {
            t2 = 0.;
          } else if (t2 > 1.) {
            t2 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = 1. - t2;
          _weights[2] = t2;
        }
        for (int i = 0; i < 3; i++) {
          closest_point[i] = closest[i];
        }

      }
      else if (_weights[0] < 0.0) {
        *min_dist2 = PDM_line_distance (x, pt1, pt2, &t1, closest_point);
        if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = 0.;
          _weights[1] = 1. - t1;
          _weights[2] = t1;
      }
      else if (_weights[1] < 0.0) {
          *min_dist2 = PDM_line_distance (x, pt2, pt3, &t1, closest_point);
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = t1;
          _weights[1] = 0;
          _weights[2] = 1. - t1;
        }
      else if (_weights[2] < 0.0) {
          *min_dist2 = PDM_line_distance (x, pt1, pt3, &t1, closest_point);
          if (t1 < 0.) {
            t1 = 0.;
          } else if (t1 > 1.) {
            t1 = 1.;
          }
          _weights[0] = t1;
          _weights[1] = 1. - t1;
          _weights[2] = 0.;
        }

        }
    return PDM_TRIANGLE_OUTSIDE;
  }
}

/**
 * \brief Computes the intersection between a line and a triangle
 *
 * \param [in]   line        Points of the line
 * \param [in]   tria_coord  Points of the triangle
 * \param [out]  ip          Intersection point
 *
 */

PDM_triangle_status_t
PDM_triangle_line_intersection
(
const double line[6],
const double tria_coord[9],
      double ip[3]
)
{

  PDM_plane_line_intersection(line, tria_coord, ip);

  PDM_triangle_status_t found = PDM_TRIANGLE_OUTSIDE;

  /* Check if intersection point is in the triangle */
  double c1[3];
  double c2[3];
  double c3[3];

  double det1, det2, det3;

  // ABD
  c1[0] = tria_coord[0]; c1[1] = tria_coord[3]; c1[2] = ip[0];
  c2[0] = tria_coord[1]; c2[1] = tria_coord[4]; c2[2] = ip[1];
  c3[0] = tria_coord[2]; c3[1] = tria_coord[5]; c3[2] = ip[2];

  det1 = _determinant_3x3(c1, c2, c3);

  // DBC
  c1[0] = ip[0]; c1[1] = tria_coord[3]; c1[2] = tria_coord[6];
  c2[0] = ip[1]; c2[1] = tria_coord[4]; c2[2] = tria_coord[7];
  c3[0] = ip[2]; c3[1] = tria_coord[5]; c3[2] = tria_coord[8];

  det2 = _determinant_3x3(c1, c2, c3);

  // ADC
  c1[0] = tria_coord[0]; c1[1] = ip[0]; c1[2] = tria_coord[6];
  c2[0] = tria_coord[1]; c2[1] = ip[1]; c2[2] = tria_coord[7];
  c3[0] = tria_coord[2]; c3[1] = ip[2]; c3[2] = tria_coord[8];

  det3 = _determinant_3x3(c1, c2, c3);

  if (det1 > 0 && det2 > 0 && det3 > 0) {

    found = PDM_TRIANGLE_INSIDE;

  } // end if point is inside triangle

  return found;

}


/**
 * \brief Computes triangle barycenter
 *
 * \param [in]   pts     Triangle vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_triangle_compute_barycenter
(
 const double pts[9],
       double bary[3]
)
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < 3; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= 3;
  }
}



/**
 * \brief Computes closest point on a triangle
 *
 * \param [in]   p              Point coordinates
 * \param [in]   v              Triangle vertices coordinates
 * \param [out]  closest_point  Closest point coordinates
 * \param [out]  min_dist2      Squared distance from point to triangle
 * \param [out]  weights        Barycentric coordinates of closest point (or NULL)
 *
 */

PDM_triangle_status_t
PDM_triangle_closest_point
(
 const double  x[3],
 const double  v[9],
 double       *closest_point,
 double       *min_dist2,
 double       *weights
 )
{
  const double *v0 = v;
  const double *v1 = v + 3;
  const double *v2 = v + 6;

  // double v20[3] = {v0[0] - v2[0],
  //                  v0[1] - v2[1],
  //                  v0[2] - v2[2]};

  // double v21[3] = {v1[0] - v2[0],
  //                  v1[1] - v2[1],
  //                  v1[2] - v2[2]};

  // double v2x[3] = {x[0] - v2[0],
  //                  x[1] - v2[1],
  //                  x[2] - v2[2]};

  // double a = PDM_DOT_PRODUCT (v20, v20);
  // double b = PDM_DOT_PRODUCT (v20, v21);
  // double c = PDM_DOT_PRODUCT (v21, v21);
  double e1[3] = {v1[0] - v0[0],
                  v1[1] - v0[1],
                  v1[2] - v0[2]};
  double e2[3] = {v2[0] - v0[0],
                  v2[1] - v0[1],
                  v2[2] - v0[2]};
  double e[3]  = {x[0] - v0[0],
                  x[1] - v0[1],
                  x[2] - v0[2]};

  double a = PDM_DOT_PRODUCT(e1, e1);
  double b = PDM_DOT_PRODUCT(e1, e2);
  double c = PDM_DOT_PRODUCT(e2, e2);

  double det = a*c - b*b;

  if (det < 1e-14) {
    return PDM_TRIANGLE_DEGENERATED;
  }

  // double r = PDM_DOT_PRODUCT (v20, v2x);
  // double s = PDM_DOT_PRODUCT (v21, v2x);
  double r = PDM_DOT_PRODUCT(e1, e);
  double s = PDM_DOT_PRODUCT(e2, e);


  /* Solve for weights of orthogonal projection of point on triangle's plane */
  double weights_local[3];
  double *_weights = weights_local;
  if (weights != NULL) {
    _weights = weights;
  }
  _weights[1] = (r*c - s*b) / det;
  _weights[2] = (s*a - r*b) / det;
  log_trace("  u = %f, v = %f\n", _weights[1], _weights[2]);
  _weights[0] = 1. - _weights[1] - _weights[2];

  /* Projection inside triangle (= closest point) */
  if ( _weights[0] >= 0.0 && _weights[0] <= 1.0 &&
       _weights[1] >= 0.0 && _weights[1] <= 1.0 &&
       _weights[2] >= 0.0 && _weights[2] <= 1.0 ) {

    *min_dist2 = 0.;
    for (int idim = 0; idim < 3; idim++) {
      // closest_point[idim] = v[6 + idim] + _weights[1]*v20[idim] + _weights[2]*v21[idim];
      closest_point[idim] = v[idim] + _weights[1]*e1[idim] + _weights[2]*e2[idim];
      double delta = x[idim] - closest_point[idim];
      *min_dist2 += delta * delta;
    }

    return PDM_TRIANGLE_INSIDE;
  }

  /* Projection inside triangle --> find closest point on triangle's boundary */
  else {

    double t01, t12, t20, d01, d12, d20, c01[3], c12[3], c20[3];

    int i, j, k;
    double t;
    d01 = PDM_line_distance (x, v0, v1, &t01, c01);
    d12 = PDM_line_distance (x, v1, v2, &t12, c12);
    d20 = PDM_line_distance (x, v2, v0, &t20, c20);
    log_trace("d01/12/20 = %f / %f / %f\n", d01, d12, d20);

    if (d01 <= d12 && d01 <= d20) {
      // i = 0;
      // j = 1;
      // k = 2;
      i = 2; j = 0; k = 1;
      *min_dist2 = d01;
      t = t01;

      closest_point[0] = c01[0];
      closest_point[1] = c01[1];
      closest_point[2] = c01[2];
    }

    else if (d12 <= d01 && d12 <= d20) {
      // i = 1;
      // j = 2;
      // k = 0;
      i = 0; j = 1; k = 2;
      *min_dist2 = d12;
      t = t12;

      closest_point[0] = c12[0];
      closest_point[1] = c12[1];
      closest_point[2] = c12[2];
    }

    else {
      // i = 2;
      // j = 0;
      // k = 1;
      i = 1; j = 0; k = 2;
      *min_dist2 = d20;
      t = t20;

      closest_point[0] = c20[0];
      closest_point[1] = c20[1];
      closest_point[2] = c20[2];
    }

    if (t < 0.) {
      t = 0.;
    } else if (t > 1.) {
      t = 1.;
    }

    // _weights[i] = 0.;
    // _weights[j] = 1. - t;
    // _weights[k] = t;

    return PDM_TRIANGLE_OUTSIDE;
  }

}



/**
 * \brief Computes the center and radius of a triangle's circumcircle
 *
 * \param [in]   vtx_coord  Triangle vertices coordinates
 * \param [out]  center     Circumcircle center
 * \param [out]  radius     Circumcircle radius
 *
 */

void
PDM_triangle_circumcircle
(
 const double  vtx_coord[9],
 double        center[3],
 double       *radius
 )
{
  double a[3] = {vtx_coord[0] - vtx_coord[6],
                 vtx_coord[1] - vtx_coord[7],
                 vtx_coord[2] - vtx_coord[8]};

  double b[3] = {vtx_coord[3] - vtx_coord[6],
                 vtx_coord[4] - vtx_coord[7],
                 vtx_coord[5] - vtx_coord[8]};

  double ab[3] = {b[0] - a[0],
                  b[1] - a[1],
                  b[2] - a[2]};


  double axb[3];
  PDM_CROSS_PRODUCT (axb, a, b);

  double a2 = PDM_DOT_PRODUCT (a, a);
  double b2 = PDM_DOT_PRODUCT (b, b);
  double ab2 = PDM_DOT_PRODUCT (ab, ab);
  double axb2 = PDM_DOT_PRODUCT (axb, axb);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

    if (axb2 == 0) {
      center[0] = vtx_coord[6];
      center[1] = vtx_coord[7];
      center[2] = vtx_coord[8];
      *radius = 0.;
    }

    else {
      axb2 = 0.25 / axb2;

      a[0] = a2*b[0] - b2*a[0];
      a[1] = a2*b[1] - b2*a[1];
      a[2] = a2*b[2] - b2*a[2];
      PDM_CROSS_PRODUCT (b, a, axb);

      center[0] = vtx_coord[6] + 2. * b[0] * axb2;
      center[1] = vtx_coord[7] + 2. * b[1] * axb2;
      center[2] = vtx_coord[8] + 2. * b[2] * axb2;
      *radius = sqrt(a2 * b2 * ab2 * axb2);
    }

  PDM_GCC_SUPPRESS_WARNING_POP
    }



#ifdef __cplusplus
}
#endif /* __cplusplus */
