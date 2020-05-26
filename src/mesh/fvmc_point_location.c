/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2007-2009  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*#include "fvmc_config.h"*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/* #include <bftc_error.h> */
/* #include <bftc_mem.h> */
/* #include <bftc_printf.h> */

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
/* #include "fvmc_config_defs.h" */
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
/* #include "fvmc_ho_basis.h" */
/* #include "fvmc_ho_location.h" */
/* #include "fvmc_nodal_priv.h" */
#include "fvmc_triangulate.h"
#include "pdm_plane.h"
#include "pdm_triangle.h"

#include "pdm_mesh_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_point_location.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-30;       /* Minimum denominator */


static int idebug = 0;
/*============================================================================
 * Private function definitions
 *============================================================================*/



/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static PDM_bool_t
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  PDM_bool_t retval = PDM_TRUE;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
           || (extents_2[i] > extents_1[i + dim])) {
      retval = PDM_FALSE;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static PDM_bool_t
_within_extents(int                 dim,
                const fvmc_coord_t  coords[],
                const double        extents[])
{
  int i;
  PDM_bool_t retval = PDM_TRUE;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
           || (coords[i] > extents[i + dim])) {
      retval = PDM_FALSE;
      break;
    }
  }

  return retval;
}


/*----------------------------------------------------------------------------
 * Updates the location[] and distance[] arrays associated with a set
 * of 1d points for points that are in a given element extent, based only
 * on this extent only (temporary, unoptimzed location).
 *
 * parameters:
 *   elt_num         <-- number of element corresponding to extents
 *   extents         <-> extents associated with element:
 *                       x_min, x_max (size: 2)
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *   location        <-> number of element containing or closest to each
 *                       point (size: n_points)
 *   distance        <-> distance from point to element indicated by
 *                       location[]: < 0 if unlocated, 0 - 1 if inside,
 *                       > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_by_extents_1d(fvmc_lnum_t         elt_num,
                      const double        extents[],
                      fvmc_lnum_t         n_points,
                      const fvmc_coord_t  point_coords[],
                      fvmc_lnum_t         location[],
                      float               distance[])
{
  fvmc_lnum_t  i;

  /* For now, we base a minimal location test on the element extents */
  /* The behavior is quadradic, nothing is optimized yet */

  for (i = 0; i < n_points; i++) {

    double elt_coord_max = -1;
    double elt_coord = -1;

    double cur_coord = point_coords[i];

    elt_coord =   (cur_coord - 0.5*(extents[1] + extents[0]))
      / (            0.5*(extents[1] - extents[0]));

    elt_coord = PDM_ABS(elt_coord);

    if (elt_coord > elt_coord_max)
      elt_coord_max = elt_coord;

    if (  (distance[i] < 0 && elt_coord_max < 1)
          || elt_coord_max < distance[i]) {

      location[i] = elt_num;
      distance[i] = (float) elt_coord_max;

    }

  }

}


/*----------------------------------------------------------------------------
 * Locate points on a 3d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_3d(fvmc_lnum_t         elt_num,
                   const fvmc_lnum_t   element_vertex_num[],
                   const fvmc_lnum_t  *parent_vertex_num,
                   const fvmc_coord_t  vertex_coords[],
                   const fvmc_coord_t  point_coords[],
                   fvmc_lnum_t         n_points_in_extents,
                   const fvmc_lnum_t   points_in_extents[],
                   double              tolerance,
                   fvmc_lnum_t         location[],
                   float               distance[])
{
  fvmc_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[3], v[3];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  if (parent_vertex_num == NULL) {
    coord_idx_0 = element_vertex_num[0] - 1;
    coord_idx_1 = element_vertex_num[1] - 1;
  }
  else {
    coord_idx_0 = parent_vertex_num[element_vertex_num[0] - 1] - 1;
    coord_idx_1 = parent_vertex_num[element_vertex_num[1] - 1] - 1;
  }

  /* Calculate edge vector and length */

  for (j = 0; j < 3; j++)
    u[j] =   vertex_coords[(coord_idx_1*3) + j]
      - vertex_coords[(coord_idx_0*3) + j];

  len2 = PDM_DOT_PRODUCT(u, u);

  if (len2 < _epsilon_denom){
    PDM_printf("warning _locate_on_edge_3d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", len2, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 3; j++)
      v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_0*3) + j];

    uv = PDM_DOT_PRODUCT(u, v);

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 3; j++)
        v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_1*3) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 3; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = PDM_DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] =  (float) sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}


/*----------------------------------------------------------------------------
 * Locate points on a 2d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_2d(fvmc_lnum_t         elt_num,
                   const fvmc_lnum_t   element_vertex_num[],
                   const fvmc_lnum_t  *parent_vertex_num,
                   const fvmc_coord_t  vertex_coords[],
                   const fvmc_coord_t  point_coords[],
                   fvmc_lnum_t         n_points_in_extents,
                   const fvmc_lnum_t   points_in_extents[],
                   double              tolerance,
                   fvmc_lnum_t         location[],
                   float               distance[])
{
  fvmc_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[2], v[2];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  if (parent_vertex_num == NULL) {
    coord_idx_0 = element_vertex_num[0] - 1;
    coord_idx_1 = element_vertex_num[1] - 1;
  }
  else {
    coord_idx_0 = parent_vertex_num[element_vertex_num[0] - 1] - 1;
    coord_idx_1 = parent_vertex_num[element_vertex_num[1] - 1] - 1;
  }

  /* Calculate edge vector and length */

  for (j = 0; j < 2; j++)
    u[j] =   vertex_coords[(coord_idx_1*2) + j]
      - vertex_coords[(coord_idx_0*2) + j];

  len2 = PDM_DOT_PRODUCT_2D(u, u);

  if (len2 < _epsilon_denom){
    PDM_printf("warning _locate_on_edge_2d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", len2, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 2; j++)
      v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_0*2) + j];

    uv = u[0]*v[0] + u[1]*v[1];

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 2; j++)
        v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_1*2) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 2; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = PDM_DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = (float) sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}


/*----------------------------------------------------------------------------
 * Locate points in a given set of 3d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 3d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_3d(fvmc_lnum_t         elt_num,
                        int                 n_triangles,
                        const fvmc_lnum_t   triangle_vertices[],
                        const fvmc_lnum_t  *parent_vertex_num,
                        const fvmc_coord_t  vertex_coords[],
                        const fvmc_coord_t  point_coords[],
                        fvmc_lnum_t         n_points_in_extents,
                        const fvmc_lnum_t   points_in_extents[],
                        const double        tolerance,
                        fvmc_lnum_t         location[],
                        float               distance[])
{

  fvmc_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  const int _order = 1;

  const int n_vtx_tria = (_order+1)*(_order+2)/2;

  double u[3], v[3], w[3];
  double uu, vv, ww, tmp_max;
  double epsilon2, dist2, vertex_dist2;

  double tolerance2 = tolerance*tolerance;

  double coords[9];
  double *pt1 = coords;
  double *pt2 = coords + 3;
  double *pt3 = coords + 6;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    if (parent_vertex_num == NULL) {
      coord_idx_0 = triangle_vertices[tria_id*n_vtx_tria]     - 1;
      coord_idx_1 = triangle_vertices[tria_id*n_vtx_tria + 1] - 1;
      coord_idx_2 = triangle_vertices[tria_id*n_vtx_tria + 2] - 1;
    }
    else {
      coord_idx_0 = parent_vertex_num[triangle_vertices[tria_id*n_vtx_tria]    - 1] - 1;
      coord_idx_1 = parent_vertex_num[triangle_vertices[tria_id*n_vtx_tria+ 1] - 1] - 1;
      coord_idx_2 = parent_vertex_num[triangle_vertices[tria_id*n_vtx_tria+ 2] - 1] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 3; j++) {
      coords[j]   = vertex_coords[(coord_idx_0*3) + j];
      coords[3+j] = vertex_coords[(coord_idx_1*3) + j];
      coords[6+j] = vertex_coords[(coord_idx_2*3) + j];
    }

    for (j = 0; j < 3; j++) {
      u[j] = pt1[j] - pt2[j];
      v[j] = pt1[j] - pt3[j];
      w[j] = pt2[j] - pt3[j];
      uu = PDM_DOT_PRODUCT(u, u);
      vv = PDM_DOT_PRODUCT(v, v);
      ww = PDM_DOT_PRODUCT(w, w);
    }

    /* epsilon2 is based on maximum edge length (squared) */

    tmp_max = PDM_MAX(uu, vv);
    tmp_max = PDM_MAX(ww, tmp_max);

    if (tolerance < 0.)
      epsilon2 = HUGE_VAL;
    else
      epsilon2 = tmp_max * tolerance2;

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      vertex_dist2 = distance[i]*distance[i];

      double *x = (double *) point_coords + 3*i;
      double closestPoint[3];
      /* double closestPointpcoords[3]; */

      double closestPointweights[3];

      /* int error = fvmc_triangle_evaluate_Position (x, coords, closestPoint, closestPointpcoords, */
      /*                                              &dist2, closestPointweights); */
      int error = PDM_triangle_evaluate_position (x,
                                                  coords,
                                                  closestPoint,
                                                  &dist2,
                                                  closestPointweights);

      if (error == -1) {
        continue;
      }

      if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
        location[i] = elt_num;
        distance[i] = (float) sqrt(dist2);
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}


/*----------------------------------------------------------------------------
 * Locate points in a given set of 2d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 2d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_2d(fvmc_lnum_t         elt_num,
                        int                 n_triangles,
                        const fvmc_lnum_t   triangle_vertices[],
                        const fvmc_lnum_t  *parent_vertex_num,
                        const fvmc_coord_t  vertex_coords[],
                        const fvmc_coord_t  point_coords[],
                        fvmc_lnum_t         n_points_in_extents,
                        const fvmc_lnum_t   points_in_extents[],
                        double              tolerance,
                        fvmc_lnum_t         location[],
                        float               distance[])
{
  fvmc_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[2], u[2], v[2], shapef[3];
  double uu, vv, uv, ut, vt, det;
  double dist, max_dist, isop_0, isop_1;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    if (parent_vertex_num == NULL) {
      coord_idx_0 = triangle_vertices[tria_id*3]     - 1;
      coord_idx_1 = triangle_vertices[tria_id*3 + 1] - 1;
      coord_idx_2 = triangle_vertices[tria_id*3 + 2] - 1;
    }
    else {
      coord_idx_0 = parent_vertex_num[triangle_vertices[tria_id*3]    - 1] - 1;
      coord_idx_1 = parent_vertex_num[triangle_vertices[tria_id*3+ 1] - 1] - 1;
      coord_idx_2 = parent_vertex_num[triangle_vertices[tria_id*3+ 2] - 1] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 2; j++) {
      u[j] = - vertex_coords[(coord_idx_0*2) + j]
        + vertex_coords[(coord_idx_1*2) + j];
      v[j] = - vertex_coords[(coord_idx_0*2) + j]
        + vertex_coords[(coord_idx_2*2) + j];
    }

    uu = PDM_DOT_PRODUCT_2D(u, u);
    vv = PDM_DOT_PRODUCT_2D(v, v);
    uv = PDM_DOT_PRODUCT_2D(u, v);

    det = (uu*vv - uv*uv);

    if (det < _epsilon_denom) {
      PDM_printf("warning _locate_on_triangles_2d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", det, _epsilon_denom);
      PDM_printf_flush();
      continue;
    }

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      /* Calculation of the barycenter coordinates for the projected node */

      for (j = 0; j < 2; j++)
        t[j] = - vertex_coords[(coord_idx_0*2) + j]
          + point_coords[i*2 + j];

      ut = u[0]*t[0] + u[1]*t[1];
      vt = v[0]*t[0] + v[1]*t[1];

      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;

      shapef[0] = 1. - isop_0 - isop_1;
      shapef[1] =      isop_0;
      shapef[2] =               isop_1;

      max_dist = -1.0;

      for (j = 0; j < 3; j++){

        dist = 2.*PDM_ABS(shapef[j] - 0.5);

        if (max_dist < dist)
          max_dist = dist;
      }

      if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
             && (max_dist < distance[i] || distance[i] < 0)) {
        location[i] = elt_num;
        distance[i] = (float) max_dist;
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron whose coordinates are pre-computed,
 * updating the location[] and distance[] arrays associated with a set
 * of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   unable_degenerated  <-- unable degenerated tetrahedra < 0
 *   tetra_coords[]      <-- tetrahedra vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_tetra(fvmc_lnum_t         elt_num,
                 fvmc_coord_t        tetra_coords[4][3],
                 const fvmc_coord_t  point_coords[],
                 fvmc_lnum_t         n_points_in_extents,
                 const fvmc_lnum_t   points_in_extents[],
                 double              tolerance,
                 fvmc_lnum_t         location[],
                 float               distance[])
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3], shapef[4];

  for (i = 0; i < 3; i++) {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }


  vol6 =  fabs( v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
                - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
                + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));

  if (vol6 < _epsilon_denom){
    PDM_printf("warning _locate_in_tetra : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  for (k = 0; k < n_points_in_extents; k++) {

    i = points_in_extents[k];

    t00  =   point_coords[i*3]     - tetra_coords[0][0];
    t10  =   point_coords[i*3 + 1] - tetra_coords[0][1];
    t20  =   point_coords[i*3 + 2] - tetra_coords[0][2];

    t01  = - tetra_coords[0][0] + tetra_coords[1][0];
    t02  = - tetra_coords[0][0] + tetra_coords[2][0];
    t03  = - tetra_coords[0][0] + tetra_coords[3][0];

    t11  = - tetra_coords[0][1] + tetra_coords[1][1];
    t12  = - tetra_coords[0][1] + tetra_coords[2][1];
    t13  = - tetra_coords[0][1] + tetra_coords[3][1];

    t21  = - tetra_coords[0][2] + tetra_coords[1][2];
    t22  = - tetra_coords[0][2] + tetra_coords[2][2];
    t23  = - tetra_coords[0][2] + tetra_coords[3][2];

    isop_0 = (  t00 * (t12*t23 - t13*t22)
                - t10 * (t02*t23 - t22*t03)
                + t20 * (t02*t13 - t12*t03)) / vol6;
    isop_1 = (- t00 * (t11*t23 - t13*t21)
              + t10 * (t01*t23 - t21*t03)
              - t20 * (t01*t13 - t03*t11)) / vol6;
    isop_2 = (  t00 * (t11*t22 - t21*t12)
                - t10 * (t01*t22 - t21*t02)
                + t20 * (t01*t12 - t11*t02)) / vol6;

    shapef[0] = 1. - isop_0 - isop_1 - isop_2;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;
    shapef[3] =                        isop_2;

    max_dist = -1.0;

    //TODO : faire un calcul fin de la distance

    for (j = 0; j < 4; j++){

      dist = 2.*PDM_ABS(shapef[j] - 0.5);

      if (max_dist < dist)
        max_dist = dist;
    }

    if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
           && (max_dist < distance[i] || distance[i] < 0)) {
      location[i] = elt_num;
      distance[i] = (float) max_dist;
    }

  }

}


/*---------------------------------------------------------------------------
 * Solve the equation "matrix.x = b" with Cramer's rule.
 *
 * parameters:
 *   m[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   1 if matrix is singular, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
    - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
    + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if (PDM_ABS(det) < _epsilon_denom)
    return 1;
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
          - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
          + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
          - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
          + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
          - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
          + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}


/*----------------------------------------------------------------------------
 * Compute 3d shape functions and their derivatives given element
 * parametric coordinates.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
 *----------------------------------------------------------------------------*/

static void
_compute_shapef_3d(fvmc_element_t  elt_type,
                   const double    uvw[3],
                   double          shapef[8],
                   double          deriv[8][3])

{

  switch (elt_type) {

  case FVMC_CELL_HEXA:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
    }

    break;

  case FVMC_CELL_PRISM:

    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }

    break;

  case FVMC_CELL_PYRAM:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] =  0.0;
      deriv[4][1] =  0.0;
      deriv[4][2] =  1.0;
    }

    break;

  default:
    PDM_error(__FILE__, __LINE__, 0,
              "_compute_shapef: unhandled element type %s\n",/* _("_compute_shapef: unhandled element type %s\n"), */
              fvmc_element_type_name[elt_type]);

  }

}


/*----------------------------------------------------------------------------
 * Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type            <-- type of element
 *   point_coords        <-- point coordinates
 *   vertex_coords[]     <-- pointer to element vertex coordinates
 *   tolerance           <-- location tolerance factor
 *   uvw[]               --> parametric coordinates of point in element
 *----------------------------------------------------------------------------*/
static int
_compute_uvw(fvmc_element_t      elt_type,
             const fvmc_coord_t  point_coords[],
             double              vertex_coords[8][3],
             double              tolerance,
             double              uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  const int order = 1;

  n_elt_vertices = fvmc_nodal_n_vertices_element(elt_type, order);

  assert(   elt_type == FVMC_CELL_HEXA
            || elt_type == FVMC_CELL_PRISM
            || elt_type == FVMC_CELL_PYRAM);

  /* Use Newton-method to determine parametric coordinates and shape function */

  for (i = 0; i < 3; i++)
    uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        a[i][j] = 0.0;
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x))
      return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance))
      return 1;

  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Locate points in a given 3d cell (other than tetrahedra or polyhedra,
 * handlesd elsewhere), updating the location[] and distance[] arrays
 * associated with a set of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   elt_type            <-- type of element
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords[]     <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_cell_3d(fvmc_lnum_t         elt_num,
                   fvmc_element_t      elt_type,
                   const fvmc_lnum_t   element_vertex_num[],
                   const fvmc_lnum_t  *parent_vertex_num,
                   const fvmc_coord_t  vertex_coords[],
                   const fvmc_coord_t  point_coords[],
                   fvmc_lnum_t         n_points_in_extents,
                   const fvmc_lnum_t   points_in_extents[],
                   double              tolerance,
                   fvmc_lnum_t         location[],
                   float               distance[])
{

  int i, j, k, n_vertices;
  fvmc_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, shapef[8],max_dist;
  double  _vertex_coords[8][3];

  const int order = 1;

  n_vertices = fvmc_nodal_n_vertices_element(elt_type, order);

  /* Initialize local element coordinates copy */

  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) {

    if (parent_vertex_num == NULL)
      coord_idx = element_vertex_num[vertex_id] -1;
    else
      coord_idx = parent_vertex_num[element_vertex_num[vertex_id] - 1] - 1;

    for (j = 0; j < 3; j++)
      _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];

  }

  /* Shape functions may be computed directly with tetrahedra */

  if (elt_type == FVMC_CELL_TETRA)

    _locate_in_tetra(elt_num,
                     _vertex_coords,
                     point_coords,
                     n_points_in_extents,
                     points_in_extents,
                     tolerance,
                     location,
                     distance);

  /* For cell shapes other than tetrahedra, find shape functions iteratively */

  else {

    for (k = 0; k < n_points_in_extents; k++) {

      i = points_in_extents[k];


      /* Check vertices (To not have a problem with pyramids) */

      int onVtx = 0;
      for (int k1 = 0; k1 < n_vertices; k1++) {
        const double *_pt =  point_coords + 3*i;
        double v[3] = {_vertex_coords[k1][0] - _pt[0],
                       _vertex_coords[k1][1] - _pt[1],
                       _vertex_coords[k1][2] - _pt[2]};

        double _dist = PDM_MODULE(v);

        if (_dist < 1e-6 * tolerance) {
          location[i] = elt_num;
          distance[i] = 0.;
          onVtx = 1;
          break;
        }
      }

      if (!onVtx) {

        if (_compute_uvw(elt_type,
                         point_coords + 3*i,
                         _vertex_coords,
                         tolerance,
                         uvw)) {

          max_dist = -1.0;

          /* For hexahedra, no need to compute shape functions, as
             the 3 parametric coordinates are simpler to use */

          if (elt_type == FVMC_CELL_HEXA) {

            for (j = 0; j < 3; j++){

              dist = 2.*PDM_ABS(uvw[j] - 0.5);

              if (max_dist < dist)
                max_dist = dist;
            }

          }

          /* For pyramids ands prisms, we need to compute shape functions */

          else {

            _compute_shapef_3d(elt_type, uvw, shapef, NULL);

            for (j = 0; j < n_vertices; j++){

              dist = 2.*PDM_ABS(shapef[j] - 0.5);

              if (max_dist < dist)
                max_dist = dist;
            }

          }

          /* For all element types, update location and distance arrays */

          if ((   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
                  && (max_dist < distance[i] || distance[i] < 0)) || (location[i] == -1)) {
            location[i] = elt_num;
            distance[i] = (float) max_dist;
          }

        }

        else {

          if (location[i] == -1) {
            location[i] = elt_num;
            distance[i] = 1.e12; // Pour les pyramides pb de convergence
          }

        }

      }

    } /* End of loop on points in extents */

  }

}



/*----------------------------------------------------------------------------
 * Find elements in a given polygonal section closest to 3d points: updates
 * the location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate (size: n_points)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, absolute distance
 *                         to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polygons_section_closest_3d(/* const fvmc_nodal_section_t   *this_section, */
                             const fvmc_lnum_t   n_elements,
                             const fvmc_lnum_t  *parent_element_num,
                             const fvmc_lnum_t  *connectivity_idx,
                             const fvmc_lnum_t  *connectivity,
                             const fvmc_lnum_t            *parent_vertex_num,
                             const fvmc_coord_t            vertex_coords[],
                             fvmc_lnum_t                   base_element_num,
                             const fvmc_coord_t            point_coords[],
                             fvmc_lnum_t                   n_point_ids,
                             const fvmc_lnum_t             point_id[],
                             fvmc_lnum_t                   location[],
                             float                        distance[])
{
  fvmc_lnum_t  i, n_vertices, vertex_id, elt_num;
  int n_triangles;

  int n_vertices_max = 0;

  fvmc_lnum_t *triangle_vertices = NULL;
  fvmc_triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (n_elements == 0)
    return;

  /* Counting loop on elements */

  for (i = 0; i < n_elements; i++) {

    n_vertices = (  connectivity_idx[i + 1]
                  - connectivity_idx[i]);

    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;

  }

  if (n_vertices_max < 3)
    return;

  triangle_vertices = malloc (sizeof(int) * (n_vertices_max-2)*3);
  state = fvmc_triangulate_state_create(n_vertices_max);

  /* Main loop on elements */

  for (i = 0; i < n_elements; i++) {

    if (base_element_num < 0) {
      if (parent_element_num != NULL)
        elt_num = parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    /* Triangulate polygon */

    n_vertices = (  connectivity_idx[i + 1]
                  - connectivity_idx[i]);
    vertex_id = connectivity_idx[i];

    //TODO: Correction provisoire bug triangulation si que des triangles dans le bloc polygons

    if (n_vertices > 4)  {

      n_triangles = fvmc_triangulate_polygon(3,
                                             n_vertices,
                                             vertex_coords,
                                             parent_vertex_num,
                                             (connectivity + vertex_id),
                                             FVMC_TRIANGULATE_MESH_DEF,
                                             triangle_vertices,
                                             state);
    }

    else if (n_vertices == 4) {


      n_triangles = fvmc_triangulate_quadrangle(3,
                                                vertex_coords,
                                                parent_vertex_num,
                                                (connectivity + vertex_id),
                                                triangle_vertices);

    }

    else {
      n_triangles = 1;

      fvmc_lnum_t *ptCur = (fvmc_lnum_t *) connectivity + vertex_id;

      triangle_vertices[0] = ptCur[0];
      triangle_vertices[1] = ptCur[1];
      triangle_vertices[2] = ptCur[2];

    }

    /* Locate on triangulated polygon */

    _locate_on_triangles_3d(elt_num,
                            n_triangles,
                            triangle_vertices,
                            parent_vertex_num,
                            vertex_coords,
                            point_coords,
                            n_point_ids,
                            point_id,
                            -1.,
                            location,
                            distance);

  } /* End of loop on elements */

  free (triangle_vertices);
  state = fvmc_triangulate_state_destroy(state);
}


/*----------------------------------------------------------------------------
 * Find elements in a given section closest to 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a line element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_closest_2d(/* const fvmc_nodal_section_t  *this_section, */
                          const fvmc_lnum_t     n_elements,
                          const fvmc_element_t  element_type,
                          const int             entity_dim,
                          const fvmc_lnum_t    *parent_element_num,
                          const fvmc_lnum_t    *connectivity,
                          const int             stride,
                          const fvmc_lnum_t           *parent_vertex_num,
                          const fvmc_coord_t           vertex_coords[],
                          fvmc_lnum_t                  base_element_num,
                          const fvmc_coord_t           point_coords[],
                          fvmc_lnum_t                  n_point_ids,
                          const fvmc_lnum_t            point_id[],
                          fvmc_lnum_t                  location[],
                          float                       distance[])
{
  fvmc_lnum_t  i, elt_num;

  /* Return immediately if nothing to do for this rank */

  if (   n_elements == 0
      || entity_dim != 1)
    return;

  assert(element_type == FVMC_EDGE);

  /* Main loop on elements */

  for (i = 0; i < n_elements; i++) {

    if (base_element_num < 0) {
      if (parent_element_num != NULL)
        elt_num = parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    /* Locate on edge */

    _locate_on_edge_2d(elt_num,
                       connectivity + i*stride,
                       parent_vertex_num,
                       vertex_coords,
                       point_coords,
                       n_point_ids,
                       point_id,
                       -1.0,
                       location,
                       distance);

  } /* End of loop on elements */
}



/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute distance to polygons
 *
 * parameters:
 *   dim               <-- dimension
 *   n_poly            <-- number of polygon
 *   connectivity_idx  <-- polygon connectivity index
 *   connectivity      <-- polygon connectivity
 *   vertex_coords     <-- polygon connectivity
 *   n_points          <-- polygon connectivity
 *   point_coords      <-- polygon connectivity
 *   distance          --> 3d surf : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *                         2d or 3d : distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_dist_closest_polygon(const int            dim,
                                const fvmc_lnum_t    n_poly,
                                const fvmc_lnum_t    connectivity_idx[],
                                const fvmc_lnum_t    connectivity[],
                                const fvmc_coord_t   vertex_coords[],
                                const fvmc_lnum_t    n_points,
                                const fvmc_lnum_t    point_ids[],
                                const fvmc_coord_t   point_coords[],
                                fvmc_lnum_t          location[],
                                float                distance[])
{
  int order = -1;
  /* fvmc_nodal_section_t* section = fvmc_nodal_section_create(FVMC_FACE_POLY, order); */
  /* section->entity_dim        = dim; */
  /* section->n_elements        = n_poly; */
  /* section->type              = FVMC_FACE_POLY; */
  /* section->connectivity_size = connectivity_idx[n_poly]; */
  /* section->stride            = 0; */
  /* section->vertex_index      = connectivity_idx; */
  /* section->vertex_num        = connectivity; */
  /* section->parent_element_num   = NULL; */
  const int base_element_num = 1;

  if (dim == 3)

    _polygons_section_closest_3d(/* section ,*/
                                 n_poly,
                                 NULL,
                                 connectivity_idx,
                                 connectivity,
                                 NULL,
                                 vertex_coords,
                                 base_element_num,
                                 point_coords,
                                 n_points,
                                 point_ids,
                                 location,
                                 distance);

  else

    _nodal_section_closest_2d(/* section, */
                              n_poly,
                              FVMC_FACE_POLY,
                              dim,
                              NULL,
                              connectivity,
                              0,
                              NULL,
                              vertex_coords,
                              base_element_num,
                              point_coords,
                              n_points,
                              point_ids,
                              location,
                              distance);


  /* fvmc_nodal_section_destroy(section); */

}






#if 0
void
PDM_nodal_block_locate_3d
(
 const int          mesh_nodal_id,
 const int          id_block,
 const int          id_part,
 const int         *pts_in_extents_idx,
 const PDM_g_num_t *pts_in_extents_g_num,
 const double      *pts_in_extents_coords,
 double            *distance
 )
{
  /*int n_parts = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);

  int id_elt = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                id_block,
                                                ipart);

    for (int ielt = 0; ielt < n_elt; ielt++) {
      //...

      id_elt++;
    }
    }*/

  PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, mesh_nodal_id);

  if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    int _id_block = id_block - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_Mesh_nodal_block_poly3d_t *block =
      (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, _id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    //...

  } // Polyhedra block

  else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
  } // Polygons block

  else {

    int _id_block = id_block;

    PDM_Mesh_nodal_block_std_t *block =
      (PDM_Mesh_nodal_block_std_t *) PDM_Handles_get (mesh->blocks_std, _id_block);

    if (block == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }

    PDM_l_num_t n_elt = block->n_elt[id_part];

    fvmc_element_t elt_type;
    switch (block->t_elt) {
      case PDM_MESH_NODAL_POINT:
        PDM_error (__FILE__, __LINE__, 0, "Invalid FMVC element type\n");
        break;
      case PDM_MESH_NODAL_BAR2:
        elt_type = FVMC_EDGE;
        break;
      case PDM_MESH_NODAL_TRIA3:
        elt_type = FVMC_FACE;
        break;
      case PDM_MESH_NODAL_QUAD4:
        elt_type = FVMC_FACE_QUAD;
        break;
      case PDM_MESH_NODAL_TETRA4:
        elt_type = FVMC_CELL_TETRA;
        break;
      case PDM_MESH_NODAL_PYRAMID5:
        elt_type = FVMC_CELL_PYRAM;
        break;
      case PDM_MESH_NODAL_PRISM6:
        elt_type = FVMC_CELL_PRISM;
        break;
      case PDM_MESH_NODAL_HEXA8:
        elt_type = FVMC_CELL_HEXA;
        break;
      default:
        printf("Unknown standard element type %d\n", block->t_elt);
        assert (1 == 0);
      }

    double *vtx_coords = (double *) PDM_Mesh_nodal_vertices_get (mesh_nodal_id, id_part);
    PDM_l_num_t cell_vtx = block->_connec[id_part];


    /*_locate_in_cell_3d (n_elt,
      elt_type,
      cell_vtx,
      const fvmc_lnum_t  *parent_vertex_num,
      vtx_coords,
      const fvmc_coord_t  point_coords[],
      fvmc_lnum_t         n_points_in_extents,
      const fvmc_lnum_t   points_in_extents[],
      double              tolerance,
      fvmc_lnum_t         location[],
      float               distance[]);*/

  } // Standard elements block

}
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */
