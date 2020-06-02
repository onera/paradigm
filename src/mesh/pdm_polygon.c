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
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_polygon.h"

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
 * Type
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/

/**
 * \brief Return a random value in [min, max]
 *
 * \param [in]  min  Minimum
 * \param [in]  max  Maximum
 *
 * \return      Random value
 */

static double
_randomVal
(
 const double min,
 const double max
)
{
  double resultat = ((double)rand())/((double)RAND_MAX);
  resultat = min + resultat * (max - min);

  return resultat;
}


/*=============================================================================
 * Public function definition
 *============================================================================*/

/**
 * \brief Get bounds
 *
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 *
 * \return      Bounds
 *
 */


double *
PDM_bounds_get
(
 const int     numPts,
 const double *pts
)
{
  double *bounds = malloc (sizeof(double) * 6);

  bounds[0] = DBL_MAX;
  bounds[1] = -DBL_MAX;
  bounds[2] = DBL_MAX;
  bounds[3] = -DBL_MAX;
  bounds[4] = DBL_MAX;
  bounds[5] = -DBL_MAX;

  for (int isom = 0; isom < numPts; isom++) {
    for (int l = 0; l < 3; l++) {
      double coord = pts[3*isom + l];
      if (bounds[2*l] > coord) {
        bounds[2*l] = coord;
      }
      if (bounds[2*l+1] < coord) {
        bounds[2*l+1] = coord;
      }
    }
  }

  return bounds;
}


/**
 * \brief Evaluates the position in a polygon
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 * \param [out] closest  Closest Point in Polygon or NULL
 * \param [out] minDist2 Square of the distance
 *
 * \return      \ref PDM_POLYGON_INSIDE or \ref PDM_POLYGON_OUTSIDE
 *              if the projected is in the polygon or not
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_evaluate_position
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double        closestPoint[3],
 double        *minDist2
 )
{
  double p0[3];
  double p10[3];
  double l10;
  double p20[3];
  double l20;
  double n[3];
  double cp[3];
  double ray[3];
  double bary[3];

  double pts_p[3*numPts];

  /*
   * Average plane
   */

  PDM_plane_normal (numPts, pts, n);

  PDM_polygon_compute_barycenter (numPts, pts, bary);

  for (int k = 0; k < numPts; k++) {
    double *pt = (double *) pts + 3*k;
    double *pt_p = pts_p + 3*k;
    PDM_plane_projection2 (pt, bary, n, pt_p);
  }

  double *_pts_p = pts_p;

  PDM_polygon_parameterize (numPts, _pts_p, p0, p10, &l10, p20, &l20, n);

  PDM_plane_projection (x,p0,n,cp);

  for (int i = 0; i < 3; i++) {
    ray[i] = cp[i] - p0[i];
  }

  double pcoords[3];

  pcoords[0] = PDM_DOT_PRODUCT(ray,p10) / (l10*l10);
  pcoords[1] = PDM_DOT_PRODUCT(ray,p20) / (l20*l20);
  pcoords[2] = 0.0;

  double bounds[6] = {DBL_MAX, -DBL_MAX,
                      DBL_MAX, -DBL_MAX,
                      DBL_MAX, -DBL_MAX};

  for (int isom = 0; isom < numPts; isom++) {
    for (int l = 0; l < 3; l++) {
      double coord = _pts_p[3*isom + l];
      if (bounds[2*l] > coord) {
        bounds[2*l] = coord;
      }
      if (bounds[2*l+1] < coord) {
        bounds[2*l+1] = coord;
      }
    }
  }

  if (pcoords[0] >= 0.0 && pcoords[0] <= 1.0 &&
      pcoords[1] >= 0.0 && pcoords[1] <= 1.0 &&
      (PDM_polygon_point_in (cp, numPts, _pts_p,
                             bounds, n) == PDM_POLYGON_INSIDE) ) {
    if (closestPoint) {
      closestPoint[0] = cp[0];
      closestPoint[1] = cp[1];
      closestPoint[2] = cp[2];
      double v[3] = {x[0] - closestPoint[0],
                     x[1] - closestPoint[1],
                     x[2] - closestPoint[2]};

      *minDist2 = PDM_DOT_PRODUCT (v, v);
    }
    return PDM_POLYGON_INSIDE;
  }

  /*
   * If here, point is outside of polygon, so need to find distance to boundary
   */

  else {
    double t, dist2;
    double closest[3];
    double *pt1, *pt2;

    if (closestPoint) {
      *minDist2 = DBL_MAX;
      for (int i=0; i<numPts; i++) {
        pt1 = (double *) pts + 3 * i;
        pt2 = (double *) pts + 3 * ((i+1)%numPts);
        dist2 = PDM_line_distance (x, pt1, pt2, &t, closest);
        if ( dist2 < *minDist2 ) {
          closestPoint[0] = closest[0];
          closestPoint[1] = closest[1];
          closestPoint[2] = closest[2];
          *minDist2 = dist2;
        }
      }
    }
    return PDM_POLYGON_OUTSIDE;
  }
}


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref PDM_TRUE except for a triangle
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_bool_t
PDM_polygon_parameterize
(
 const int     numPts,
 const double *pts,
 double       *p0,
 double       *p10,
 double       *l10,
 double       *p20,
 double       *l20,
 double       *n
)
{
  double s, t, p[3], p1[3], p2[3], sbounds[2], tbounds[2];

  if (numPts < 3) {
    return PDM_FALSE;
  }

  /*
   *  This is a two pass process: first create a p' coordinate system
   *  that is then adjusted to insure that the polygon points are all in
   *  the range 0<=s,t<=1.  The p' system is defined by the polygon normal,
   *  first vertex and the first edge.
   */

  PDM_plane_normal (numPts, pts, n);

  double x1[3];
  x1[0] = pts[0];
  x1[1] = pts[1];
  x1[2] = pts[2];

  double x2[3];
  x2[0] = pts[3];
  x2[1] = pts[3+1];
  x2[2] = pts[3+2];

  for (int i = 0; i < 3; i++) {
    p0[i] = x1[i];
    p10[i] = x2[i] - x1[i];
  }

  PDM_CROSS_PRODUCT(p20,n,p10);

  /*
   * Determine lengths of edges
   */

  if ( ((*l10)= PDM_DOT_PRODUCT(p10,p10)) == 0.0
    || ((*l20)= PDM_DOT_PRODUCT(p20,p20)) == 0.0 ) {
    return PDM_FALSE;
  }

  /*
   *  Now evalute all polygon points to determine min/max parametric
   *  coordinate values.
   *
   * first vertex has (s,t) = (0,0)
   */

  sbounds[0] = 0.0;
  sbounds[1] = 0.0;

  tbounds[0] = 0.0;
  tbounds[1] = 0.0;

  for (int i = 1; i < numPts; i++) {
    x1[0] = pts[3*i];
    x1[1] = pts[3*i+1];
    x1[2] = pts[3*i+2];
    for (int j = 0; j < 3; j++) {
      p[j] = x1[j] - p0[j];
    }

    s = PDM_DOT_PRODUCT (p,p10) / (*l10);
    t = PDM_DOT_PRODUCT (p,p20) / (*l20);

    sbounds[0] = (s<sbounds[0]?s:sbounds[0]);
    sbounds[1] = (s>sbounds[1]?s:sbounds[1]);
    tbounds[0] = (t<tbounds[0]?t:tbounds[0]);
    tbounds[1] = (t>tbounds[1]?t:tbounds[1]);
  }

  /*
   * Re-evaluate coordinate system
   */

  for (int i = 0; i < 3; i++) {
    p1[i] = p0[i] + sbounds[1]*p10[i] + tbounds[0]*p20[i];
    p2[i] = p0[i] + sbounds[0]*p10[i] + tbounds[1]*p20[i];
    p0[i] = p0[i] + sbounds[0]*p10[i] + tbounds[0]*p20[i];
    p10[i] = p1[i] - p0[i];
    p20[i] = p2[i] - p0[i];
  }

  (*l10) = PDM_MODULE(p10);
  (*l20) = PDM_MODULE(p20);

  return PDM_TRUE;
}


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */


#define _TOL_POLY 1.e-05
#define _POLYGON_INTERSECTION 2
#define _POLYGON_ON_LINE 3

#define _POLYGON_CERTAIN 1
#define _POLYGON_UNCERTAIN 0
#define _POLYGON_RAY_TOL 1.e-03 //Tolerance for ray firing
#define _POLYGON_MAX_ITER 10    //Maximum iterations for ray-firing
#define _POLYGON_VOTE_THRESHOLD 2


PDM_polygon_status_t
PDM_polygon_point_in
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double       *bounds,
 double       *n
)
{
  double *x1, *x2, xray[3], u, v;
  double rayMag, mag=1, ray[3];
  int rayOK, status;
  int iterNumber;
  int maxComp, comps[2];
  int deltaVotes;

  /*
   * Do a quick bounds check
   */

  if ( x[0] < bounds[0] || x[0] > bounds[1] ||
       x[1] < bounds[2] || x[1] > bounds[3] ||
       x[2] < bounds[4] || x[2] > bounds[5]) {

    return PDM_POLYGON_OUTSIDE;
  }

  /*
   * Define a ray to fire. The ray is a random ray normal to the
   *  normal of the face. The length of the ray is a function of the
   *  size of the face bounding box.
   */

  for (int i = 0; i < 3; i++) {
    ray[i] = ( bounds[2*i+1] - bounds[2*i] )*1.1 +
          fabs((bounds[2*i+1] + bounds[2*i])/2.0 - x[i]);
  }

  if ((rayMag = PDM_MODULE(ray)) < 1.e-15 ) {
    return PDM_POLYGON_OUTSIDE;
  }

  /* Get the maximum component of the normal. */

  if (fabs(n[0]) > fabs(n[1])) {
    if (fabs(n[0]) > fabs(n[2])) {
      maxComp  = 0;
      comps[0] = 1;
      comps[1] = 2;
    }
    else {
      maxComp  = 2;
      comps[0] = 0;
      comps[1] = 1;
    }
  }
  else {
    if (fabs(n[1]) > fabs(n[2])) {
      maxComp  = 1;
      comps[0] = 0;
      comps[1] = 2;
    }
    else {
      maxComp = 2;
      comps[0] = 0;
      comps[1] = 1;
    }
  }

  /* Check that max component is non-zero */

  if (fabs(n[maxComp]) < 1.e-15) {
    return PDM_POLYGON_DEGENERATED;
  }

  /*
   * Enough information has been acquired to determine the random ray.
   * Random rays are generated until one is satisfactory (i.e.,
   * produces a ray of non-zero magnitude).  Also, since more than one
   * ray may need to be fired, the ray-firing occurs in a large loop.
   *
   * The variable iterNumber counts the number of iterations and is
   * limited by the defined variable _POLYGON_MAX_ITER.
   *
   * The variable deltaVotes keeps track of the number of votes for
   * "in" versus "out" of the face.  When delta_vote > 0, more votes
   * have counted for "in" than "out".  When delta_vote < 0, more votes
   * have counted for "out" than "in".  When the delta_vote exceeds or
   * equals the defined variable _POLYGON_VOTE_THRESHOLD, than the
   * appropriate "in" or "out" status is returned.
   */

  for (deltaVotes = 0, iterNumber = 1;
       (iterNumber < _POLYGON_MAX_ITER)
         && (PDM_ABS(deltaVotes) < _POLYGON_VOTE_THRESHOLD);
       iterNumber++) {

    /*
     * Generate ray
     */

    for (rayOK = PDM_FALSE; rayOK == PDM_FALSE; ) {
      ray[comps[0]] = _randomVal (-rayMag, rayMag);
      ray[comps[1]] = _randomVal (-rayMag, rayMag);
      ray[maxComp] = -(n[comps[0]]*ray[comps[0]] +
                       n[comps[1]]*ray[comps[1]]) / n[maxComp];
      if ( (mag = PDM_MODULE(ray)) > rayMag * _TOL_POLY ) {
        rayOK = PDM_TRUE;
      }
    }

    /*
     * The ray must be appropriately sized.
     */

    for (int i = 0; i < 3; i++) {
      xray[i] = x[i] + (rayMag/mag)*ray[i];
    }

    /*
     * The ray may now be fired against all the edges
     */


    int testResult = _POLYGON_CERTAIN;
    int numInts = 0;
    for (int i = 0; i < numPts; i++) {
      x1 = (double *) pts + 3*i;
      x2 = (double *) pts + 3*((i+1)%numPts);

      /*
       * Fire the ray and compute the number of intersections.  Be careful
       * of degenerate cases (e.g., ray intersects at vertex).
       */

      if ((status = PDM_line_intersection (x,xray,x1,x2, &u,&v)) == PDM_LINE_INTERSECT_YES) {

        /*
         * This test checks for vertex and edge intersections
         * For example
         *  Vertex intersection
         *    (u=0 v=0), (u=0 v=1), (u=1 v=0), (u=1 v=0)
         *  Edge intersection
         *    (u=0 v!=0 v!=1), (u=1 v!=0 v!=1)
         *    (u!=0 u!=1 v=0), (u!=0 u!=1 v=1)
         */

        if ( (0. < u) && (u < 1.0) &&
             (0. < v) && (v < 1.0) ) {
        /* if ( (_POLYGON_RAY_TOL < u) && (u < 1.0 - _POLYGON_RAY_TOL) && */
        /*      (_POLYGON_RAY_TOL < v) && (v < 1.0 - _POLYGON_RAY_TOL) ) { */
          numInts++;
        }
        else {
          testResult = _POLYGON_UNCERTAIN;
        }
      }

      else if ( status == _POLYGON_ON_LINE ) {
        testResult = _POLYGON_UNCERTAIN;
      }

    }

    if ( testResult == _POLYGON_CERTAIN ) {
      if ( numInts % 2 == 0) {
        --deltaVotes;
      }
      else {
        ++deltaVotes;
      }
    }
  } /* try another ray */

  /*
   * If the number of intersections is odd, the point is in the polygon.
   */

  if (deltaVotes <= 0) {
    return PDM_POLYGON_OUTSIDE;
  }
  else {
    return PDM_POLYGON_INSIDE;
  }
}

#undef _TOL_POLY
#undef _POLYGON_INTERSECTION
#undef _POLYGON_ON_LINE

#undef _POLYGON_CERTAIN
#undef _POLYGON_UNCERTAIN
#undef _POLYGON_RAY_TOL
#undef _POLYGON_MAX_ITER
#undef _POLYGON_VOTE_THRESHOLD


/**
 * \brief Computes polygon barycenter
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  bary    Barycenter
 *
 *
 */


void
PDM_polygon_compute_barycenter
(
const int numPts,
const double *pts,
double bary[3]
 )
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < numPts; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= numPts;
  }
}


//--------------------
void
PDM_polygon_orthonormal_basis
(
 const int     n_vtx,
 const double *vtx_xyz,
 double       *tangent_u,
 double       *tangent_v,
 double       *normal
 )
{
  const double eps = 1e-15;

  /* Compute polygon normal */
  PDM_plane_normal (n_vtx, vtx_xyz, normal);

  /* Find one suitable tangent direction */
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    double mag = 0;
    int jvtx = (ivtx + 1) % n_vtx;

    for (int idim = 0; idim < 3; idim++) {
      tangent_u[idim] = vtx_xyz[3*jvtx + idim] - vtx_xyz[3*ivtx + idim];
      mag += tangent_u[idim] * tangent_u[idim];
    }

    if (mag > eps) {
      mag = 1. / sqrt(mag);
      for (int idim = 0; idim < 3; idim++) {
        tangent_u[idim] *= mag;
      }
    }
  }

  /* Complete orthonormal basis with second tangent direction */
  PDM_CROSS_PRODUCT (tangent_v, normal, tangent_u);
}



void PDM_polygon_compute_uv_coordinates
(
 const int     n_pts,
 const double *xyz,
 const double  orig_xyz[3],
 const double  tangent_u[3],
 const double  tangent_v[3],
 double       *uv
 )
{
  for (int ipt = 0; ipt < n_pts; ipt++) {
    const double *_xyz = xyz + 3 * ipt;
    double *_uv  = uv + 2 * ipt;

    _uv[0] = 0;
    _uv[1] = 0;

    for (int idim = 0; idim < 3; idim++) {
      double d = _xyz[idim] - orig_xyz[idim];

      _uv[0] += d * tangent_u[idim];
      _uv[1] += d * tangent_v[idim];
    }
  }
}


PDM_polygon_status_t PDM_polygon_point_in_2d
(
 const double  uv[2],
 const int     n_vtx,
 const double *vtx_uv,
 double       *bounds
 )
{

  /*
   * Do a quick bounds check
   */
  if (bounds == NULL) {
    bounds = malloc (sizeof(double) * 4);

    bounds[0] =  DBL_MAX;
    bounds[1] = -DBL_MAX;
    bounds[2] =  DBL_MAX;
    bounds[3] = -DBL_MAX;

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 2; idim++) {
        bounds[2*idim]     = PDM_MIN (bounds[2*idim],     vtx_uv[2*ivtx + idim]);
        bounds[2*idim + 1] = PDM_MAX (bounds[2*idim + 1], vtx_uv[2*ivtx + idim]);
      }
    }
  }

  if ( uv[0] < bounds[0] || uv[0] > bounds[1] ||
       uv[1] < bounds[2] || uv[1] > bounds[3]) {
    return PDM_POLYGON_OUTSIDE;
  }


  /*
   * Define a random ray
   */
  double scale = sqrt( (bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
                       (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) );

  double ray[2], ray_mag;
  do {
    ray[0] = _randomVal(0., 1.);
    ray[1] = _randomVal(0., 1.);

    ray_mag = sqrt (ray[0] * ray[0] + ray[1] * ray[1]);
  } while (ray_mag < 1.e-15);

  scale /= ray_mag;
  ray[0] = uv[0] + scale * ray[0];
  ray[1] = uv[1] + scale * ray[1];

  /*
   * Count intersections between the ray and the edges of the polygon
   */
  int count = 0;

  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    int jvtx = (ivtx + 1) % n_vtx;

    double t1, t2;
    PDM_line_intersect_t stat = PDM_line_intersection_2d (uv,
                                                          ray,
                                                          vtx_uv + 2*ivtx,
                                                          vtx_uv + 2*jvtx,
                                                          &t1,
                                                          &t2);
    if (stat == PDM_LINE_INTERSECT_YES) {
      count++;
    }

  }

  if (count % 2 == 0) {
    return PDM_POLYGON_OUTSIDE;
  } else {
    return PDM_POLYGON_INSIDE;
  }
}




#if 0
PDM_polygon_status_t PDM_polygon_point_in_3d
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 double       *normal[3],
 double       *bounds
 )
{
  if (normal == NULL) {
    normal = malloc (sizeof(double) * 3);

    PDM_plane_normal (n_vtx,
                      vtx_xyz,
                      normal);
  }


  /*
   * Do a quick bounds check
   */
  if (bounds == NULL) {
    bounds = malloc (sizeof(double) * 6);

    bounds[0] =  DBL_MAX;
    bounds[1] = -DBL_MAX;
    bounds[2] =  DBL_MAX;
    bounds[3] = -DBL_MAX;
    bounds[4] =  DBL_MAX;
    bounds[5] = -DBL_MAX;

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 3; idim++) {
        bounds[2*idim]     = PDM_MIN (bounds[2*idim],     vtx_xyz[3*ivtx + idim]);
        bounds[2*idim + 1] = PDM_MAX (bounds[2*idim + 1], vtx_xyz[3*ivtx + idim]);
      }
    }
  }

  if ( xyz[0] < bounds[0] || xyz[0] > bounds[1] ||
       xyz[1] < bounds[2] || xyz[1] > bounds[3] ||
       xyz[2] < bounds[4] || xyz[2] > bounds[5]) {

    return PDM_POLYGON_OUTSIDE;
  }


  /*
   * Define a random ray
   */
  double scale = sqrt( (bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
                       (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
                       (bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));

  double ray[3], ray_mag;
  do {
    double ray_n = 0.0;
    for (int idim = 0; idim < 3; idim++) {
      ray[idim] = _randomVal(0., 1.);
      ray_n += ray[idim] * normal[idim];
    }

    double ray_mag = 0.0;
    for (int idim = 0; idim < 3; idim++) {
      ray[idim] -= ray_n * normal[idim];
      ray_mag += ray[idim] * ray[idim];
    }

  } while (ray_mag < 1e-15);

  scale /= sqrt (ray_mag);
  for (int idim = 0; idim < 3; idim++) {
    ray[idim] = xyz[idim] + scale * ray[idim];
  }

  /* ... */
}
#endif


/*
 *
 * http://geomalgorithms.com/a03-_inclusion.html
 */

PDM_polygon_status_t PDM_polygon_point_in_2d_wn
(
 const double  xy[2],
 const int     n_vtx,
 const double *vtx_xy,
 double       *bounds
 )
{
  double eps = 1e-12;

  /*
   * Do a quick bounds check
   */
  if (bounds == NULL) {
    bounds = malloc (sizeof(double) * 4);

    bounds[0] =  DBL_MAX;
    bounds[1] = -DBL_MAX;
    bounds[2] =  DBL_MAX;
    bounds[3] = -DBL_MAX;

    for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
      for (int idim = 0; idim < 2; idim++) {
        bounds[2*idim]     = PDM_MIN (bounds[2*idim],     vtx_xy[2*ivtx + idim]);
        bounds[2*idim + 1] = PDM_MAX (bounds[2*idim + 1], vtx_xy[2*ivtx + idim]);
      }
    }
  }

  if ( xy[0] < bounds[0] || xy[0] > bounds[1] ||
       xy[1] < bounds[2] || xy[1] > bounds[3]) {
    return PDM_POLYGON_OUTSIDE;
  }

  double scale = 0;
  for (int idim = 0; idim < 2; idim++) {
    scale = PDM_MAX (scale, bounds[2*idim + 1] - bounds[2*idim]);
  }

  if (scale > 1.0) {
    eps *= scale;
  }

  double eps2 = eps * eps;


  /* Compute winding number */
  int wn = 0;

  for (int i = 0; i < n_vtx; i++) {
    int ip = (i+1) % n_vtx;
    const double *p0 = vtx_xy + 2*i;
    const double *p1 = vtx_xy + 2*ip;

    double dx = xy[0] - p0[0];
    double dy = xy[1] - p0[1];

    //-->>
#if 1
    if (dx*dx + dy*dy < eps2) {
      /* Point coincident with vertex i */
      return PDM_POLYGON_INSIDE;
    }
#endif
    //<<--

    if (dy >= 0) {//

      if (p1[1] > xy[1]) {
        /* Upward crossing */
        double ex = p1[0] - p0[0];
        double ey = p1[1] - p0[1];

        double s = ex*dy - ey*dx;
        if (s > eps) {//if (s > 0) {
          wn++;
        } else if (s > -eps) {
          double denom = ex*ex + ey*ey;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along line (i, ip) */
            double t = (ex*dx + ey*dy) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }
      }

    } else {

      if (p1[1] < xy[1]) {
        /* Downward crossing */
        double ex = p1[0] - p0[0];
        double ey = p1[1] - p0[1];

        double s = ex*dy - ey*dx;
        if (s < -eps) {//if (s < 0) {
          wn--;
        } else if (s < eps) {
          double denom = ex*ex + ey*ey;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along line (i, ip) */
            double t = (ex*dx + ey*dy) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }

      }

    }
  } /* Loop over edges */

  /* Result */
  if (wn == 0) {
    return PDM_POLYGON_OUTSIDE;
  } else {
    return PDM_POLYGON_INSIDE;
  }

}







PDM_polygon_status_t PDM_polygon_point_in_3d_wn
(
 const double  xyz[2],
 const int     n_vtx,
 const double *vtx_xyz,
 const double  normal[3],
 const double  characteristic_length,
 double       *bounds
 )
{
  const double eps_base = 1e-12;

  const double eps = eps_base * characteristic_length;
  const double eps2 = eps * eps;

  /*
   * Do a quick bounds check
   */
  if (bounds != NULL) {
    if ( xyz[0] < bounds[0] || xyz[0] > bounds[1] ||
         xyz[1] < bounds[2] || xyz[1] > bounds[3] ||
         xyz[2] < bounds[4] || xyz[2] > bounds[5]) {
      return PDM_POLYGON_OUTSIDE;
    }
  }

  double max_n = 0.0;
  int dim_max = -1;

  for (int idim = 0; idim < 3; idim++) {
    if (-normal[idim] > max_n) {
      max_n = -normal[idim];
      dim_max = idim;
    } else if (normal[idim] > max_n) {
      max_n = normal[idim];
      dim_max = idim;
    }
  }

  if (max_n < eps) {
    return PDM_POLYGON_DEGENERATED;
  }

  /* 'Project' onto cartesian plane most orthogonal to normal */
  int dim1 = (dim_max + 1) % 3;
  int dim2 = (dim_max + 2) % 3;

  /* Compute winding number */
  int wn = 0;

  double u = xyz[dim1];
  double v = xyz[dim2];

  double u0 = vtx_xyz[dim1];
  double v0 = vtx_xyz[dim2];

  double u1, v1;
  double du, dv, eu, ev, s, denom, t;

  for (int i = 0; i < n_vtx; i++) {
    int ip = (i+1) % n_vtx;

    u1 = vtx_xyz[3*ip + dim1];
    v1 = vtx_xyz[3*ip + dim2];

    du = u - u0;
    dv = v - v0;

#if 1
    if (du*du + dv*dv < eps2) {
      /* Point coincident with vertex i */
      return PDM_POLYGON_INSIDE;
    }
#endif

    if (dv >= 0) {

      if (v1 > v) {
        /* Upward crossing */
        eu = u1 - u0;
        ev = v1 - v0;

        s = eu*dv - ev*du;
        if (s > eps) {
          wn++;
        } else if (s > -eps) {
          denom = eu*eu + ev*ev;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along line (i, ip) */
            t = (eu*du + ev*dv) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }
      }

    } else {

      if (v1 < v) {
        /* Downward crossing */
        eu = u1 - u0;
        ev = v1 - v0;

        s = eu*dv - ev*du;
        if (s < -eps) {//if (s < 0) {
          wn--;
        } else if (s < eps) {
          denom = eu*eu + ev*ev;
          if (denom < eps2) {
            /* Vertices i and ip are coincident */
            return PDM_POLYGON_DEGENERATED;
          } else {
            /* Compute parameter along line (i, ip) */
            t = (eu*du + ev*dv) / denom;
            if (t > -eps && t < 1 + eps) {
              return PDM_POLYGON_INSIDE;
            }
          }
        }

      }

    }

    u0 = u1;
    v0 = v1;

  } /* Loop over edges */

  /* Result */
  if (wn == 0) {
    return PDM_POLYGON_OUTSIDE;
  } else {
    return PDM_POLYGON_INSIDE;
  }

}






#ifdef __cplusplus
}
#endif /* __cplusplus */
