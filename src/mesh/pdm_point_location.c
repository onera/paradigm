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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_plane.h"
#include "pdm_line.h"
#include "pdm_triangle.h"
#include "pdm_polygon.h"
#include "pdm_triangulate.h"
#include "pdm_mesh_nodal.h"
//#include "pdm_mesh_nodal_priv.h"
#include "pdm_mean_values.h"

#include "pdm_point_location.h"

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
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static double _epsilon_denom = 1.e-30;       /* Minimum denominator */

/*=============================================================================
 * Private function definition
 *============================================================================*/

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
    - A[1][0]*(A[0][1]*A[2][2] - A[2][1]*A[0][2])
    + A[2][0]*(A[0][1]*A[1][2] - A[1][1]*A[0][2]);

  if (PDM_ABS(det) < _epsilon_denom) {
    printf("_solve_3x3: det = %e\n", det);
    return PDM_FALSE;
  } else {
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

  x[0] = x0; x[1] = x1; x[2] = x2;

  return PDM_TRUE;
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
_compute_shapef_3d
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               uvw[3],
 double                     shapef[8],
 double                     deriv[8][3]
 )
{
  switch (elt_type) {
  case PDM_MESH_NODAL_HEXA8:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] =        uvw[0]  * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] =        uvw[0]  *        uvw[1]  * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) *        uvw[1]  * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) *        uvw[2];
    shapef[5] =        uvw[0]  * (1.0 - uvw[1]) *        uvw[2];
    shapef[6] =        uvw[0]  *        uvw[1]  *        uvw[2];
    shapef[7] = (1.0 - uvw[0]) *        uvw[1]  *        uvw[2];

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

  case PDM_MESH_NODAL_PRISM6:

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

  case PDM_MESH_NODAL_PYRAMID5:

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

    //--->>
  case PDM_MESH_NODAL_QUAD4:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]);
    shapef[1] =        uvw[0]  * (1.0 - uvw[1]);
    shapef[2] =        uvw[0]  *        uvw[1] ;
    shapef[3] = (1.0 - uvw[0]) *        uvw[1] ;

    if (deriv != NULL) {
      deriv[0][0] = uvw[1] - 1.0;
      deriv[0][1] = uvw[0] - 1.0;
      deriv[1][0] = 1.0 - uvw[1];
      deriv[1][1] = -uvw[0];
      deriv[2][0] = uvw[1];
      deriv[2][1] = uvw[0];
      deriv[3][0] = -uvw[1];
      deriv[3][1] = 1.0 - uvw[0];
    }

    break;
    //<<---

  default:
    /*PDM_error(__FILE__, __LINE__, 0,
      "_compute_shapef: unhandled element type %s\n",
      fvmc_element_type_name[elt_type]);*/
    PDM_error (__FILE__, __LINE__, 0, "Wrong element type\n");

  }

}

//-->>
#if 1
void PDM_point_location_distance
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const int                  n_pts,
 const double               uvw[],
 double                     shapef[],
 double                     distance[]
 )
{
  int i, j;
  const int order = 1;
  int n_vtx = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  for (i = 0; i < n_pts; i++) {
    double *s = shapef + n_vtx * i;

    /* Shape functions */
    if (elt_type == PDM_MESH_NODAL_TETRA4) {
      s[0] = 1.0 - uvw[3*i] - uvw[3*i+1] - uvw[3*i+2];
      s[1] =       uvw[3*i];
      s[2] =                  uvw[3*i+1];
      s[3] =                               uvw[3*i+2];
    }

    else {
      _compute_shapef_3d (elt_type,
                          uvw + 3*i,
                          s,
                          NULL);
    }


    /* "Distance" */
    double max_dist2 = 0., dist2;
    if (elt_type == PDM_MESH_NODAL_HEXA8) {
      for (j = 0; j < 3; j++) {
        dist2 = 2. * PDM_ABS (uvw[3*i+j] - 0.5);
        max_dist2 = PDM_MAX (max_dist2, dist2);
      }
    }

    else {
      for (j = 0; j < n_vtx; j++) {
        dist2 = 2. * PDM_ABS (s[j] - 0.5);
        max_dist2 = PDM_MAX (max_dist2, dist2);
      }
    }

    distance[i] = sqrt(max_dist2);
  }
}
#endif
//<<--

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

static PDM_bool_t
_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[8][3],
 const double               tolerance,
 double                     uvw[3]
 )
{
  int i, j, n_elt_vertices, iter;
  const int max_iter = 30;
  const double tolerance2 = tolerance * tolerance;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  /* Get number of vertices */
  const int order = 1;
  n_elt_vertices = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

  assert (elt_type == PDM_MESH_NODAL_PYRAMID5 ||
          elt_type == PDM_MESH_NODAL_PRISM6   ||
          elt_type == PDM_MESH_NODAL_HEXA8);

  /* Use Newton-method to determine parametric coordinates and shape function */
  for (i = 0; i < 3; i++) {
    uvw[i] = 0.5;
  }

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d (elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        a[i][j] = 0.0;
      }
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

    if (_solve_3x3(a, b, x) == PDM_FALSE) {
      printf("_compute_uvw: singular matrix\n");
      return PDM_FALSE;
    }

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= tolerance2) {
      return PDM_TRUE;
    }

  }

  return PDM_FALSE;
}



/*----------------------------------------------------------------------------
 * Locate points on a 3d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   edge_vtx           <-- ids of edge vertices
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   vtx_coord          <-- pointer to vertex coordinates
 *   pts_coord          <-- xyz-coordinates of points to locate
 *   n_pts              <-- number of points to locate
 *   distance           --> distance from point to edge (size: n_pts)
 *   bar_coord          --> barycentric coordinates of points projections on edge (size: 2*n_pts)
 *----------------------------------------------------------------------------*/
static void
_locate_on_edge_3d
(
 const PDM_l_num_t  edge_vtx[2],
 const PDM_l_num_t *parent_vertex_num,
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 const PDM_g_num_t  pts_g_num[],//debug only
 float              distance[],
 double             bar_coord[]
 )
{
  const int DEBUG = 0;

  int idim, ipt;
  PDM_l_num_t id0, id1;

  /* vertex index of current edge */
  id0 = edge_vtx[0] - 1;
  id1 = edge_vtx[1] - 1;

  if (parent_vertex_num != NULL) {
    id0 = parent_vertex_num[id0] - 1;
    id1 = parent_vertex_num[id1] - 1;
  }

  /* Calculate edge vector and length */
  double u[3], uu;
  for (idim = 0; idim < 3; idim++) {
    u[idim] = vtx_coord[3*id1 + idim] - vtx_coord[3*id0 + idim];
  }

  uu = PDM_DOT_PRODUCT (u, u);

  if (uu < _epsilon_denom){
    PDM_printf("warning _locate_on_edge_3d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", uu, _epsilon_denom);
    PDM_printf_flush();
    return;
  }

  else {
    uu = 1. / uu;
  }

  /* Loop on points to locate on edge */
  double v[3], uv, t;
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 2 * ipt;

    if (DEBUG) {
      printf("pt %d (%ld) : (%f %f %f)\n",
             ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);
    }

    /* Calculate linear coordinates of projection of point on edge axis */
    for (idim = 0; idim < 3; idim++) {
      v[idim] = _pt[idim] - vtx_coord[3*id0 + idim];
    }

    uv = PDM_DOT_PRODUCT (u, v);

    t = uv * uu;


    /* Set v to be the vector from the point to the closest point on
       the segment (if t < 0, v is already that vector) */
    if (t >= 1.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] = _pt[idim] - vtx_coord[3*id1 + idim];
      }
    }

    else if (t > 0.) {
      for (idim = 0; idim < 3; idim++) {
        v[idim] -= t * u[idim];
      }
    }

    /* Distance between point to locate and its projection */
    distance[ipt] = (float) sqrt(PDM_DOT_PRODUCT (v, v));

    /* Barycentric coordinates */
    _bc[0] = 1. - t;
    _bc[1] = t;

    if (DEBUG) {
      double sum = 0;
      printf(" bc = ");
      for (int ivtx = 0; ivtx < 2; ivtx++) {
        printf("%f ", _bc[ivtx]);
        sum += _bc[ivtx];
      }
      printf("  (sum = %f)\n", sum);
      printf("dist = %f\n\n", distance[ipt]);
    }

  } // End of loop on points

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
_locate_on_triangles_3d
(
 const int          n_tri,
 const PDM_l_num_t  tri_vtx[],
 const PDM_l_num_t *parent_vertex_num,
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 const PDM_g_num_t  pts_g_num[],//debug only
 const double       tolerance,
 int                location[],
 float              distance[],
 double             bar_coord[]
 )
{
  const int DEBUG = 0;
  const int CHECK = 0;

  int ipt, itri, idim, ivtx;
  PDM_l_num_t id0, id1, id2;

  double u[3], v[3], w[3];
  double uu, vv, ww;
  double epsilon2;

  double tolerance2 = tolerance * tolerance;

  double tri_coord[9];
  double *vtx0 = tri_coord;
  double *vtx1 = tri_coord + 3;
  double *vtx2 = tri_coord + 6;

  /* const int _order = 1;
     const int n_vtx_tria = (_order+1)*(_order+2)/2;*/
  const int n_vtx_tria = 3;

  /* Initialize distance of points to locate */
  for (ipt = 0; ipt < n_pts; ipt++) {
    distance[ipt] = HUGE_VAL;
  }

  if (location != NULL) {
    for (ipt = 0; ipt < n_pts; ipt++) {
      location[ipt] = -1;
    }
  }

  /* Loop on triangles */
  for (itri = 0; itri < n_tri; itri++) {

    /* vertex index of current triangle */
    id0 = tri_vtx[itri*n_vtx_tria]     - 1;
    id1 = tri_vtx[itri*n_vtx_tria + 1] - 1;
    id2 = tri_vtx[itri*n_vtx_tria + 2] - 1;

    if (parent_vertex_num != NULL) {
      id0 = parent_vertex_num[id0] - 1;
      id1 = parent_vertex_num[id1] - 1;
      id2 = parent_vertex_num[id2] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */
    for (idim = 0; idim < 3; idim++) {
      tri_coord[    idim] = vtx_coord[3*id0 + idim];
      tri_coord[3 + idim] = vtx_coord[3*id1 + idim];
      tri_coord[6 + idim] = vtx_coord[3*id2 + idim];
    }

    for (idim = 0; idim < 3; idim++) {
      u[idim] = vtx0[idim] - vtx1[idim];
      v[idim] = vtx0[idim] - vtx2[idim];
      w[idim] = vtx1[idim] - vtx2[idim];
    }
    uu = PDM_DOT_PRODUCT(u, u);
    vv = PDM_DOT_PRODUCT(v, v);
    ww = PDM_DOT_PRODUCT(w, w);

    /* epsilon2 is based on maximum edge length (squared) */
    if (tolerance < 0.) {
      epsilon2 = HUGE_VAL;
    } else {
      epsilon2 = PDM_MAX (uu, vv);
      epsilon2 = PDM_MAX (epsilon2, ww);
      epsilon2 = tolerance2 * epsilon2;
    }

    /* Loop on points to locate */
    for (ipt = 0; ipt < n_pts; ipt++) {
      const double *_pt = pts_coord + 3 * ipt;
      double       *_bc = bar_coord + 3 * ipt;

      if (DEBUG) {
        printf("\npt %d (%ld) : (%f %f %f)\n",
               ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);

        printf("tri_coord =\n");
        for (ivtx = 0; ivtx < 3; ivtx++) {
          for (idim = 0; idim < 3; idim++) {
            printf("%f ", tri_coord[3*ivtx + idim]);
          }
          printf("\n");
        }
      }

      double closest_point[3];
      double dist2;
      double weights[3];
      PDM_triangle_status_t stat = PDM_triangle_evaluate_position (_pt,
                                                                   tri_coord,
                                                                   closest_point,
                                                                   &dist2,
                                                                   weights);
      if (DEBUG) {
        printf("  closest point : (%f %f %f), at dist = %f, with weights %f %f %f \n",
               closest_point[0], closest_point[1], closest_point[2],
               sqrt(dist2), weights[0], weights[1], weights[2]);
      }

      if (stat == PDM_TRIANGLE_DEGENERATED) {
        continue;
      }

      if (dist2 < distance[ipt]*distance[ipt]) {
        if (bar_coord != NULL) {
          _bc[0] = weights[1];
          _bc[1] = weights[2];
          _bc[2] = weights[0];
        }

        distance[ipt] = (float) sqrt(dist2);

        if (location != NULL) {
          location[ipt] = itri;
        }
      }

      if (DEBUG && bar_coord != NULL) {
        double sum = 0;
        printf(" bc = ");
        for (ivtx = 0; ivtx < 3; ivtx++) {
          printf("%f ", _bc[ivtx]);
          sum += _bc[ivtx];
        }
        printf("  (sum = %f)\n", sum);
        printf("dist = %f\n", distance[ipt]);
      }

      if (CHECK && bar_coord != NULL) {
        double err[3] = {_pt[0], _pt[1], _pt[2]};
        for (ivtx = 0; ivtx < 3; ivtx++) {
          for (idim = 0; idim < 3; idim++) {
            err[idim] -= _bc[ivtx] * tri_coord[3*ivtx + idim] ;
          }
        }

        printf("pt %d (%ld) : dist = %f\t ; linear precision = %f\n",
               ipt, pts_g_num[ipt], distance[ipt], PDM_MODULE(err));
      }

    } // End of loop on points

  } // End of loop on triangles
}







static void
_locate_on_quad_3d
(
 const PDM_l_num_t  quad_vtx[4],
 const PDM_l_num_t *parent_vertex_num,
 const double       vtx_coord[],
 const int          n_pts,
 const double       pts_coord[],
 const PDM_g_num_t  pts_g_num[],//debug only
 const double       tolerance,
 float              distance[],
 double             bar_coord[]
 )
{
  const int DEBUG = 0;
  const int CHECK = 0;

  int ipt, ivtx, id_vtx, idim;

  double quad_coord[12];
  for (ivtx = 0; ivtx < 4; ivtx++) {
    id_vtx = quad_vtx[ivtx] - 1;

    if (parent_vertex_num != NULL) {
      id_vtx = parent_vertex_num[id_vtx] - 1;
    }

    for (idim = 0; idim < 3; idim++) {
      quad_coord[3*ivtx + idim] = vtx_coord[3*id_vtx + idim];
    }
  }

  if (DEBUG) {
    printf("\n\nquad_coord =\n");
    for (ivtx = 0; ivtx < 4; ivtx++) {
      printf("  %f %f %f\n", quad_coord[3*ivtx], quad_coord[3*ivtx+1], quad_coord[3*ivtx+2]);
    }
  }

  PDM_mean_values_polygon_compute3 (n_pts,
                                    pts_coord,
                                    4,
                                    quad_coord,
                                    bar_coord);


  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 3 * ipt;
    double *_bc = bar_coord + 4 * ipt;

    if (DEBUG) {
      printf("\npt %d (%ld) : (%f %f %f)\n",
             ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);
    }

    double v_cp_p[3] = {_pt[0], _pt[1], _pt[2]};
    for (ivtx = 0; ivtx < 4; ivtx++) {
      for (idim = 0; idim < 3; idim++) {
        v_cp_p[idim] -= _bc[ivtx] * quad_coord[3*ivtx+idim];
      }
    }
    double dist2 = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
    distance[ipt] = (float) sqrt(dist2);

    if (DEBUG) {
      double sum = 0;
      printf(" bc = ");
      for (ivtx = 0; ivtx < 4; ivtx++) {
        printf("%f ", _bc[ivtx]);
        sum += _bc[ivtx];
      }
      printf("  (sum = %f)\n", sum);
      printf("dist = %f\n", distance[ipt]);
    }

    if (CHECK) {
      double err[3] = {_pt[0], _pt[1], _pt[2]};
      for (ivtx = 0; ivtx < 4; ivtx++) {
        for (idim = 0; idim < 3; idim++) {
          err[idim] -= _bc[ivtx] * quad_coord[3*ivtx+idim] ;
        }
      }

      printf("pt %d (%ld) : dist = %f\t ; linear precision = %f\n",
             ipt, pts_g_num[ipt], distance[ipt], PDM_MODULE(err));
    }
  }
}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron whose coordinates are pre-computed,
 * updating the location[] and distance[] arrays associated with a set
 * of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   tetra_coords[]      <-- tetrahedra vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points            <-- number of points in element extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/
static void
_locate_in_tetra
(
 const double        tetra_coord[4][3],
 const int           n_pts,
 const double        pts_coord[],
 const PDM_g_num_t   pts_g_num[],//debug only
 float              *distance,
 double             *bar_coord
 )
{
  const int DEBUG = 0;
  const int CHECK = 0;
  int ivtx, idim, ipt, i, j, k;

  double v[3][3];
  for (ivtx = 0; ivtx < 3; ivtx++) {
    for (idim = 0; idim < 3; idim++) {
      v[ivtx][idim] = tetra_coord[ivtx+1][idim] - tetra_coord[0][idim];
    }
  }

  double vol6 = v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
    v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
    v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]);

  PDM_bool_t flip = PDM_FALSE;
  if (vol6 > 0) {
    flip = PDM_TRUE;
  } else {
    vol6 = -vol6;
  }

  if (vol6 < _epsilon_denom){
    PDM_printf("warning _locate_in_tetra : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    PDM_printf_flush();
    return;
  }
  vol6 = 1./vol6;

  double r[3][3];
  for (i = 0; i < 3; i++) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    PDM_CROSS_PRODUCT (r[i], v[j], v[k]);

    for (idim = 0; idim < 3; idim++) {
      r[i][idim] *= vol6;
      if (flip == PDM_TRUE) {
        r[i][idim] = -r[i][idim];
      }
    }
  }

  /* Loop on points */
  double vp0[3];
  for (ipt = 0; ipt < n_pts; ipt++) {
    const double *_pt = pts_coord + 3 * ipt;
    double       *_bc = bar_coord + 4 * ipt;

    if (DEBUG) {
      printf("pt %d (%ld) : (%f %f %f)\n",
             ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);
    }

    for (idim = 0; idim < 3; idim++) {
      vp0[idim] = tetra_coord[0][idim] - _pt[idim];
    }

    /* Compute barycentric coordinates of current point in tetrahedron */
    _bc[0] = 1.;
    for (i = 0; i < 3; i++) {
      _bc[i+1] = PDM_DOT_PRODUCT (vp0, r[i]);
      _bc[0] -= _bc[i+1];
    }

    /* Compute 'distance' from current point to tetrahedron */
    double max_dist = -1.;
    for (ivtx = 0; ivtx < 4; ivtx++) {
      double dist = PDM_ABS (_bc[ivtx] - 0.5);

      if (dist > max_dist) {
        max_dist = dist;
      }
    }

    distance[ipt] = (float) (2 * max_dist);

    /* Point strictly inside cell */
    if (distance[ipt] < 1) {
      distance[ipt] = 0.;
    }

    if (DEBUG) {
      double sum = 0;
      printf(" bc = ");
      for (ivtx = 0; ivtx < 4; ivtx++) {
        printf("%f ", _bc[ivtx]);
        sum += _bc[ivtx];
      }
      printf("  (sum = %f)\n", sum);
      printf("dist = %f\n\n", distance[ipt]);
    }

    if (CHECK) {
      double err[3] = {_pt[0], _pt[1], _pt[2]};
      for (ivtx = 0; ivtx < 4; ivtx++) {
        for (idim = 0; idim < 3; idim++) {
          err[idim] -= _bc[ivtx] * tetra_coord[ivtx][idim] ;
        }
      }

      printf("pt %d (%ld) : dist = %f\t ; linear precision = %f\n",
             ipt, pts_g_num[ipt], distance[ipt], PDM_MODULE(err));
    }

  } // End of loop on points

}

#if 1
//-->>
void
PDM_locate_points_in_tetra (const double       vtx_xyz[12],
                            const int          n_pts,
                            const double       pts_xyz[],
                            const PDM_g_num_t  pts_g_num[],
                            float             *distance,
                            double            *bary_coords)
{
  //const double tolerance = 1e-6;

  double tetra_coords[4][3];
  for (int ivtx = 0; ivtx < 4; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      tetra_coords[ivtx][idim] = vtx_xyz[3*ivtx + idim];
    }
  }

  _locate_in_tetra (tetra_coords,
                    n_pts,
                    pts_xyz,
                    pts_g_num,
                    distance,
                    bary_coords);
}
//<<--
#endif


static void
_locate_in_cell_3d
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const PDM_l_num_t           cell_vtx[],
 const PDM_l_num_t          *parent_vertex_num,//?
 const double                vtx_coord[],
 const int                   n_pts,
 const double                pts_coord[],
 const PDM_g_num_t           pts_g_num[],//debug only
 const double                tolerance,
 float                      *distance,
 double                     *bar_coord
 )
{
  const int DEBUG = 0;
  const int CHECK = 0;

  double uvw[3];
  double _vtx_coord[8][3];
  double dist, max_dist;

  const int order = 1;
  const int n_vtx = PDM_Mesh_nodal_n_vertices_element (elt_type, order);


  /* Initialize local element coordinates copy */
  PDM_l_num_t id_vtx;
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {

    if (parent_vertex_num == NULL) {
      id_vtx = cell_vtx[ivtx] - 1;
    } else {
      id_vtx = parent_vertex_num[cell_vtx[ivtx] - 1] - 1;
    }

    for (int idim = 0; idim < 3; idim++) {
      _vtx_coord[ivtx][idim] = vtx_coord[3*id_vtx + idim];
    }
  }


  /* Tetrahedron */
  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    /* Shape functions may be computed directly */
    _locate_in_tetra (_vtx_coord,
                      n_pts,
                      pts_coord,
                      pts_g_num,//debug only
                      distance,
                      bar_coord);
  }

  /* Other cell types */
  else {
    /* Shape functions must be computed iteratively */
    for (int ipt = 0; ipt < n_pts; ipt++) {

      const double *_pt = pts_coord + 3 * ipt;
      if (DEBUG) {
        printf("pt %d (%ld) : (%f %f %f)\n",
               ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);
      }

      double *_bc = bar_coord + n_vtx * ipt;

      /* Check vertices (To not have a problem with pyramids) */
      PDM_bool_t on_vtx = PDM_FALSE;
      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {

        double v[3] = {_vtx_coord[ivtx][0] - _pt[0],
                       _vtx_coord[ivtx][1] - _pt[1],
                       _vtx_coord[ivtx][2] - _pt[2]};

        double _dist = PDM_MODULE (v);

        if (_dist < 1e-6 * tolerance) {
          distance[ipt] = 0.;
          for (int i = 0; i < n_vtx; i++) {
            _bc[i] = 0.;
          }
          _bc[ivtx] = 1.;

          on_vtx = PDM_TRUE;
          break;
        }

      } // Loop on vertices


      if (on_vtx == PDM_FALSE) {

        /* Compute parametric coordinates */
        if (_compute_uvw(elt_type,
                         _pt,
                         _vtx_coord,
                         tolerance,
                         uvw) == PDM_TRUE) {

          _compute_shapef_3d (elt_type, uvw, _bc, NULL);

          max_dist = -1.0;
          if (elt_type == PDM_MESH_NODAL_HEXA8) {
            for (int i = 0; i < 3; i++){
              dist = PDM_ABS (uvw[i] - 0.5);

              if (max_dist < dist) {
                max_dist = dist;
              }
            }
          } else {
            for (int ivtx = 0; ivtx < n_vtx; ivtx++){
              dist = PDM_ABS (_bc[ivtx] - 0.5);

              if (max_dist < dist) {
                max_dist = dist;
              }
            }
          }

          distance[ipt] = (float) 2. * max_dist;
        }


        /* Failed to compute parametric coordinates */
        else {
          distance[ipt] = 1.e12;
          for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
            _bc[ivtx] = 0.;
          }

        }
      }

      /* Point strictly inside cell */
      if (distance[ipt] < 1) {
        distance[ipt] = 0.;
      }

      if (DEBUG) {
        double sum = 0;
        printf(" bc = ");
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          printf("%f ", _bc[ivtx]);
          sum += _bc[ivtx];
        }
        printf("  (sum = %f)\n", sum);
        printf("dist = %f\n\n", distance[ipt]);
      }

      if (CHECK) {
        double err[3] = {_pt[0], _pt[1], _pt[2]};
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          for (int idim = 0; idim < 3; idim++) {
            err[idim] -= _bc[ivtx] * _vtx_coord[ivtx][idim] ;
          }
        }

        printf("pt %d (%ld) : dist = %f\t ; linear precision = %f\n",
               ipt, pts_g_num[ipt], distance[ipt], PDM_MODULE(err));
      }

    } // End of loop on points

  }

}


static void
_std_block_locate_3d
(
 const int           mesh_nodal_id,
 const int           id_block,
 const int          *pts_idx,
 const double       *pts_coord,
 const PDM_g_num_t  *pts_g_num,//
 double              tolerance,
 float              *distance,
 // double             *projected_coord,// high-order
 double             *bar_coord
 )
{
  const int DEBUG = 0;

  int n_parts = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);
  int n_elt = 0;
  int n_pts = 0;

  PDM_l_num_t *element_vertex_num = NULL;
  const PDM_l_num_t          *parent_vertex_num = NULL;//???

  const int order = 1;
  PDM_Mesh_nodal_elt_t elt_type = PDM_Mesh_nodal_block_type_get (mesh_nodal_id,
                                                                 id_block);

  int entity_dim = 0;
  if (elt_type == PDM_MESH_NODAL_POINT) {
    entity_dim = 0;
  }
  else if (elt_type == PDM_MESH_NODAL_BAR2) {
    entity_dim = 1;
  }
  else if (elt_type == PDM_MESH_NODAL_TRIA3 ||
           elt_type == PDM_MESH_NODAL_QUAD4) {
    entity_dim = 2;
  }
  else if (elt_type == PDM_MESH_NODAL_TETRA4   ||
           elt_type == PDM_MESH_NODAL_PYRAMID5 ||
           elt_type == PDM_MESH_NODAL_PRISM6   ||
           elt_type == PDM_MESH_NODAL_HEXA8) {
    entity_dim = 3;
  }
  else {
    PDM_error (__FILE__, __LINE__, 0, "Wrong standard element type\n");
  }


  int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (elt_type,
                                                     order);


  /* Loop on parts */
  int idx = -1;
  int ipt = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {

    if (DEBUG) {
      printf("\n\n------\nid_block = %d, ipart = %d\n",
             id_block, ipart);
    }

    n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                            id_block,
                                            ipart);

    const double *vertex_coord = PDM_Mesh_nodal_vertices_get (mesh_nodal_id,
                                                              ipart);
    PDM_Mesh_nodal_block_std_get (mesh_nodal_id,
                                  id_block,
                                  ipart,
                                  &element_vertex_num);

    //-->>//if (DEBUG)
    PDM_g_num_t *_gnum = PDM_Mesh_nodal_g_num_get (mesh_nodal_id,
                                                   id_block,
                                                   ipart);

    //<<--

    /* Loop on elements */
    for (int ielt = 0; ielt < n_elt; ielt++) {
      idx++;

      n_pts = pts_idx[idx + 1] - pts_idx[idx];
      //-->>
      if (DEBUG) {
        printf("\n\nelt (%ld) : idx = %d, pts_idx = %d, n_pts = %d\n",
               _gnum[ielt],
               idx,
               pts_idx[idx],
               n_pts);
      }
      //<<--


      if (n_pts < 1) {
        continue;
      }

      /* Volume elements */
      if (entity_dim == 3) {

        /* First-order elements */
        if (order == 1) {
          _locate_in_cell_3d (elt_type,
                              element_vertex_num + ielt * n_vtx_elt,
                              parent_vertex_num,
                              vertex_coord,
                              n_pts,
                              pts_coord + pts_idx[idx] * 3,
                              pts_g_num + pts_idx[idx],
                              tolerance,
                              distance + ipt,
                              bar_coord + ipt * n_vtx_elt);
        }

        /* High-order elements */
        else {
          assert (0 == 1);
        }

      }


      /* Surface elements */
      else if (entity_dim == 2) {

        /* First-order elements */
        if (order == 1) {
          /*    Quadrangles */
          if (elt_type == PDM_MESH_NODAL_QUAD4) {
            _locate_on_quad_3d (element_vertex_num + ielt * n_vtx_elt,
                                parent_vertex_num,
                                vertex_coord,
                                n_pts,
                                pts_coord + pts_idx[idx] * 3,
                                pts_g_num + pts_idx[idx],//debug
                                tolerance,
                                distance + ipt,
                                bar_coord + ipt * n_vtx_elt);
          }

          /*    Triangles */
          else {
            assert (elt_type == PDM_MESH_NODAL_TRIA3);

            _locate_on_triangles_3d (1,
                                     element_vertex_num + ielt * n_vtx_elt,
                                     parent_vertex_num,
                                     vertex_coord,
                                     n_pts,
                                     pts_coord + pts_idx[idx] * 3,
                                     pts_g_num + pts_idx[idx],//debug
                                     tolerance,
                                     NULL,
                                     distance + ipt,
                                     bar_coord + ipt * n_vtx_elt);
          }

        }

        /* High-order elements */
        else {
          assert (0 == 1);
        }

      }


      /* Linear elements */
      else if (entity_dim == 1) {

        /* First-order elements */
        if (order == 1) {
          _locate_on_edge_3d (element_vertex_num + ielt * n_vtx_elt,
                              parent_vertex_num,
                              vertex_coord,
                              n_pts,
                              pts_coord + pts_idx[idx] * 3,
                              pts_g_num + pts_idx[idx],//debug
                              distance + ipt,
                              bar_coord + ipt * n_vtx_elt);
        }

        /* High-order elements */
        else {
          assert (0 == 1);
        }

      }


      /* Points */
      else {
        const double *vtx_coord = vertex_coord + 3*(element_vertex_num[ielt] - 1);

        if (DEBUG) {
          printf(" elt (%ld) is a single point (%f %f %f)\n",
                 _gnum[ielt], vtx_coord[0], vtx_coord[1], vtx_coord[2]);
        }

        for (int i = 0; i < n_pts; i++) {
          int id_pt = pts_idx[idx] + i;

          double dist2 = 0.;
          for (int idim = 0; idim < 3; idim++) {
            double delta = pts_coord[3*id_pt + idim] - vtx_coord[idim];
            dist2 += delta * delta;
          }
          distance[ipt] = (float) sqrt(dist2);
          bar_coord[ipt] = 1.;

          if (DEBUG) {
            printf("\tpt (%ld) (%f %f %f) distance = %f\n",
                   pts_g_num[id_pt],
                   pts_coord[3*id_pt], pts_coord[3*id_pt+1], pts_coord[3*id_pt+2],
                   distance[ipt]);
          }
        }

      }


      ipt += n_pts;
    } // End of loop on elements

  } // End of loop on parts

}


void
_poly2d_block_locate
(
 const int          mesh_nodal_id,
 const int          id_block,
 const int         *pts_idx,
 const double      *pts_coord,
 const PDM_g_num_t *pts_g_num,//debug
 float             *distance,
 int               *bar_coord_idx,
 double            *bar_coord
 )
{
  const int DEBUG = 0;
  const int CHECK = 0;

  int n_elt, n_pts, n_vtx;

  int n_parts = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);

  PDM_l_num_t *connec_idx = NULL;
  PDM_l_num_t *connec = NULL;
  const PDM_l_num_t          *parent_vertex_num = NULL;//???


  int n_vtx_max = 10;
  double *poly_coord = malloc (sizeof(double) * n_vtx_max * 3);


  /* Loop on parts */
  int idx = -1;
  int ipt = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {

    n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                            id_block,
                                            ipart);

    const double *vertex_coord = PDM_Mesh_nodal_vertices_get (mesh_nodal_id,
                                                              ipart);

    PDM_Mesh_nodal_block_poly2d_get (mesh_nodal_id,
                                     id_block,
                                     ipart,
                                     &connec_idx,
                                     &connec);

    /* Loop on elements */
    for (int ielt = 0; ielt < n_elt; ielt++) {
      idx++;

      n_pts = pts_idx[idx + 1] - pts_idx[idx];
      n_vtx = connec_idx[ielt + 1] - connec_idx[ielt];

      if (n_vtx_max < n_vtx) {
        n_vtx_max = PDM_MAX (2*n_vtx_max, n_vtx);
        poly_coord = realloc (poly_coord, sizeof(double) * n_vtx_max * 3);
      }

      for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
        int id_vtx = connec[connec_idx[ielt] + ivtx] - 1;

        if (parent_vertex_num != NULL) {
          id_vtx = parent_vertex_num[id_vtx];
        }

        for (int idim = 0; idim < 3; idim++) {
          poly_coord[3*ivtx + idim] = vertex_coord[3*id_vtx + idim];
        }
      }
      if (DEBUG) {
        printf("\n\npoly_coord =\n");
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          printf("  %f %f %f\n", poly_coord[3*ivtx], poly_coord[3*ivtx+1], poly_coord[3*ivtx+2]);
        }
      }

      PDM_mean_value_coordinates_polygon_3d (n_vtx,
                                             poly_coord,
                                             n_pts,
                                             pts_coord + pts_idx[idx] * 3,
                                             bar_coord + bar_coord_idx[ipt]);

      for (int i = 0; i < n_pts; i++) {
        int id_pt = pts_idx[idx] + i;
        const double *_pt = pts_coord + 3 * id_pt;
        double       *_bc = bar_coord + bar_coord_idx[ipt];

        if (DEBUG) {
          printf("\npt %d (%ld) : (%f %f %f)\n",
                 ipt, pts_g_num[id_pt], _pt[0], _pt[1], _pt[2]);
        }

        double v_cp_p[3] = {_pt[0], _pt[1], _pt[2]};
        for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
          for (int idim = 0; idim < 3; idim++) {
            v_cp_p[idim] -= _bc[ivtx] * poly_coord[3*ivtx+idim];
          }
        }
        double dist2 = PDM_DOT_PRODUCT (v_cp_p, v_cp_p);
        distance[ipt] = (float) sqrt(dist2);

        if (DEBUG) {
          double sum = 0;
          printf(" bc = ");
          for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
            printf("%f ", _bc[ivtx]);
            sum += _bc[ivtx];
          }
          printf("  (sum = %f)\n", sum);
          printf("dist = %f\n", distance[ipt]);
        }

        if (CHECK) {
          double err[3] = {_pt[0], _pt[1], _pt[2]};
          for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
            for (int idim = 0; idim < 3; idim++) {
              err[idim] -= _bc[ivtx] * poly_coord[3*ivtx+idim] ;
            }
          }

          printf("pt %d (%ld) : dist = %f\t ; linear precision = %f\n",
                 ipt, pts_g_num[id_pt], distance[ipt], PDM_MODULE(err));
        }

        ipt++;
      }

    } //End loop on elements

  } // End loop on parts
}



void
_poly3d_block_locate
(
 const int          mesh_nodal_id,
 const int          id_block,
 const int         *pts_idx,
 const double      *pts_coords,
 const PDM_g_num_t *pts_g_num,//debug
 double             tolerance,//?
 float             *distance,
 int               *bar_coords_idx,
 double            *bar_coords
 )
{
  const int DEBUG = 1;
  const PDM_l_num_t          *parent_vertex_num = NULL;//???

  const double _eps_loc = 1e-5;

  double *solid_angle = NULL;
  int s_solid_angle = 0;

  float *min_dist = NULL;
  int s_min_dist = 0;

  /* Loop on parts */
  int n_parts = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);

  int ielt = 0;
  int ipt = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {
    int n_cell = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                 id_block,
                                                 ipart);

    const double *vertex_coords = PDM_Mesh_nodal_vertices_get (mesh_nodal_id,
                                                               ipart);

    PDM_l_num_t  n_face;
    PDM_l_num_t *face_vtx_idx  = NULL;
    PDM_l_num_t *face_vtx      = NULL;
    PDM_l_num_t *cell_face_idx = NULL;
    PDM_l_num_t *cell_face     = NULL;
    PDM_Mesh_nodal_block_poly3d_get (mesh_nodal_id,
                                     id_block,
                                     ipart,
                                     &n_face,
                                     &face_vtx_idx,
                                     &face_vtx,
                                     &cell_face_idx,
                                     &cell_face);

    /* Counting loop on faces */
    PDM_l_num_t n_vtx_face = 0;
    PDM_l_num_t n_vtx_face_max = 0;

    for (int iface = 0; iface < n_face; iface++) {
      n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];
      n_vtx_face_max = PDM_MAX (n_vtx_face, n_vtx_face_max);
    }


    PDM_l_num_t *triangle_vertices = malloc (sizeof(PDM_l_num_t) * (n_vtx_face_max-2)*3);
    PDM_triangulate_state_t *state = PDM_triangulate_state_create (n_vtx_face_max);

    /* Loop on elements */
    for (int icell = 0; icell < n_cell; icell++) {
      int n_pts = pts_idx[icell+1] - pts_idx[icell];

      if (n_pts < 1) {
        continue;
      }

      if (solid_angle == NULL) {
        s_solid_angle = n_pts;
        solid_angle = malloc (sizeof(double) * s_solid_angle);
      } else {
        if (s_solid_angle < n_pts) {
          s_solid_angle = n_pts;
          solid_angle = realloc (solid_angle, sizeof(double) * s_solid_angle);
        }
      }

      if (min_dist == NULL) {
        s_min_dist = n_pts;
        min_dist = malloc (sizeof(float) * s_min_dist);
      } else {
        if (s_min_dist < n_pts) {
          s_min_dist = n_pts;
          min_dist = realloc (min_dist, sizeof(float) * s_min_dist);
        }
      }

      for (int i = 0; i < n_pts; i++) {
        solid_angle[i] = 0.;
        min_dist[i] = FLT_MAX;
      }

      /* Loop on element faces */
      for (int j = cell_face_idx[icell]; j < cell_face_idx[icell+1]; j++) {
        PDM_l_num_t iface = PDM_ABS(cell_face[j]) - 1;

        n_vtx_face = face_vtx_idx[iface+1] - face_vtx_idx[iface];

        const PDM_l_num_t *_vertex_num = face_vtx + face_vtx_idx[iface];

        PDM_l_num_t n_triangles;
        /* Triangular face */
        if (n_vtx_face == 3) {
          n_triangles = 1;
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            triangle_vertices[ivtx] = _vertex_num[ivtx];
          }
        }

        /* Quadrilateral face */
        else if (n_vtx_face == 4) {
          n_triangles = PDM_triangulate_quadrangle (3,
                                                    vertex_coords,
                                                    parent_vertex_num,
                                                    _vertex_num,
                                                    triangle_vertices);
        }

        /* Polygonal face */
        else {
          n_triangles = PDM_triangulate_polygon(3,
                                                n_vtx_face,
                                                vertex_coords,
                                                parent_vertex_num,
                                                _vertex_num,
                                                PDM_TRIANGULATE_MESH_DEF,
                                                triangle_vertices,
                                                state);
        }

        /* Loop on face triangles so as to loop on tetrahedra
           built by joining face triangles and pseudo-center */
        for (int itri = 0; itri < n_triangles; itri++) {

          PDM_l_num_t coord_id[3];
          coord_id[0] = triangle_vertices[itri*3    ] - 1;
          coord_id[1] = triangle_vertices[itri*3 + 1] - 1;
          coord_id[2] = triangle_vertices[itri*3 + 2] - 1;

          if (parent_vertex_num != NULL) {
            coord_id[0] = parent_vertex_num[coord_id[0]] - 1;
            coord_id[1] = parent_vertex_num[coord_id[1]] - 1;
            coord_id[2] = parent_vertex_num[coord_id[2]] - 1;
          }

          if (cell_face[j] < 0) {
            PDM_l_num_t tmp = coord_id[0];
            coord_id[0] = coord_id[2];
            coord_id[2] = tmp;
          }

          double closest[3];
          double tria_coords[9];
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            for (int idim = 0; idim < 3; idim++) {
              tria_coords[3*ivtx + idim] = vertex_coords[3*coord_id[ivtx] + idim];
            }
          }

          /* Loop on points (compute solid angle and min_dist) */
          for (int i = 0; i < n_pts; i++) {
            const double *_point_coords = pts_coords + 3 * (pts_idx[ielt] + i);

            double minDist2;
            double closestPointweights[3];

            PDM_triangle_status_t error = PDM_triangle_evaluate_position ((double *) _point_coords,
                                                                          tria_coords,
                                                                          closest,
                                                                          &minDist2,
                                                                          closestPointweights);
            if (error == PDM_TRIANGLE_DEGENERATED) {
              continue;
            }

            float dist = (float) sqrt (minDist2);

            if (min_dist[i] > dist) {
              min_dist[i] = dist;
            }



            double v_pt_to_a[3], v_pt_to_b[3], v_pt_to_c[3];
            for (int idim = 0; idim < 3; idim++) {
              v_pt_to_a[idim] = vertex_coords[3*coord_id[0] + idim] - _point_coords[idim];
              v_pt_to_b[idim] = vertex_coords[3*coord_id[1] + idim] - _point_coords[idim];
              v_pt_to_c[idim] = vertex_coords[3*coord_id[2] + idim] - _point_coords[idim];
            }

            double det_abc = _determinant_3x3 (v_pt_to_a, v_pt_to_b, v_pt_to_c);

            double n_a = PDM_MODULE (v_pt_to_a);
            double n_b = PDM_MODULE (v_pt_to_b);
            double n_c = PDM_MODULE (v_pt_to_c);

            double dot_ab = PDM_DOT_PRODUCT (v_pt_to_a, v_pt_to_b);
            double dot_ac = PDM_DOT_PRODUCT (v_pt_to_a, v_pt_to_c);
            double dot_bc = PDM_DOT_PRODUCT (v_pt_to_b, v_pt_to_c);

            double denom = n_a * n_b * n_c
              + dot_ab * n_c
              + dot_ac * n_b
              + dot_bc * n_a;

            double half_angle = atan2(det_abc, denom);

            if ((half_angle < 0.) && (det_abc > 0)) {
              half_angle = 2*M_PI - half_angle;
            }
            else if ((half_angle > 0.) && (det_abc < 0)){
              half_angle = -2*M_PI + half_angle;
            }

            solid_angle[i] += 2 * half_angle;
          } // End of loop on points

        } // End of loop on face triangles

      } // End of loop on element faces


      /* Loop on points */
      for (int i = 0; i < n_pts; i++) {
#if 1
        distance[ipt] = 1.e12;
#else
        /* Point inside polyhedron */
        if (fabs(solid_angle[i]) >= _eps_loc) {
          distance[ipt] = 1.e-3;
          // compute barycentric coordinates ...
        }

        /* Point outside polyhedron */
        else {
          distance[ipt] = 1 + min_dist[i];
          // compute barycentric coordinates ...
        }
#endif

        //-->>
        for (int j = bar_coords_idx[ipt]; j < bar_coords_idx[ipt+1]; j++) {
          bar_coords[j] = -1.;
        }
        //<<--

        ipt++;
      } // End of loop on points

    } // End of loop on elements

    ielt += n_cell;

    free (triangle_vertices);
    state = PDM_triangulate_state_destroy (state);
  } // End of loop on parts

  if (solid_angle != NULL) {
    free (solid_angle);
  }

  if (min_dist != NULL) {
    free (min_dist);
  }
}


/*============================================================================
 * Public function definitions
 *============================================================================*/
void
PDM_point_location_nodal
(
 const int            mesh_nodal_id,
 const int            n_pts,
 const int           *pts_idx,
 const double        *pts_coords,
 const PDM_g_num_t   *pts_g_num,//debug only
 const double         tolerance,
 int                  base_element_num,//?
 float              **distance,
 double             **projected_coords,//high-order
 int                **bar_coords_idx,
 double             **bar_coords
 )
{
  const int DEBUG = 0;

  const int order = 1;

  int idx = 0;
  int n_elt = 0;
  int n_vtx_elt = 0;
  int n_blocks = PDM_Mesh_nodal_n_blocks_get (mesh_nodal_id);
  int n_parts  = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);
  int *blocks_id = PDM_Mesh_nodal_blocks_id_get (mesh_nodal_id);

  *distance         = malloc (sizeof(float) * n_pts);
  //*projected_coords = malloc (sizeof(double) * n_pts * 3);//high-order
  *bar_coords_idx   = malloc (sizeof(int) * (n_pts+1));


  /* First loop to allocate */
  (*bar_coords_idx)[0] = 0;
  int ipt = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
      for (int ipart = 0; ipart < n_parts; ipart++) {
        n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                id_block,
                                                ipart);

        PDM_l_num_t *cellvtx_idx = NULL;
        PDM_l_num_t *cellvtx = NULL;
        PDM_Mesh_nodal_block_poly3d_cell_vtx_connect_get (mesh_nodal_id,
                                                          id_block,
                                                          ipart,
                                                          &cellvtx_idx,
                                                          &cellvtx);
        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_elt = cellvtx_idx[ielt+1] - cellvtx_idx[ielt];

          for (int i = pts_idx[idx]; i < pts_idx[idx+1]; i++) {
            (*bar_coords_idx)[ipt + 1] = (*bar_coords_idx)[ipt] + n_vtx_elt;
            ipt++;
          }
          idx++;
        }
      }
    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
      for (int ipart = 0; ipart < n_parts; ipart++) {
        n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                id_block,
                                                ipart);

        PDM_l_num_t *connec_idx = NULL;
        PDM_l_num_t *connec     = NULL;
        PDM_Mesh_nodal_block_poly2d_get (mesh_nodal_id,
                                         id_block,
                                         ipart,
                                         &connec_idx,
                                         &connec);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_elt = connec_idx[ielt + 1] - connec_idx[ielt];

          for (int i = pts_idx[idx]; i < pts_idx[idx+1]; i++) {
            (*bar_coords_idx)[ipt + 1] = (*bar_coords_idx)[ipt] + n_vtx_elt;
            ipt++;
          }
          idx++;
        }
      }
    }

    /* Standard elements */
    else {
      PDM_Mesh_nodal_elt_t elt_type = PDM_Mesh_nodal_block_type_get (mesh_nodal_id,
                                                                     id_block);
      n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (elt_type, order);

      for (int ipart = 0; ipart < n_parts; ipart++) {
        n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                id_block,
                                                ipart);
        for (int ielt = 0; ielt < n_elt; ielt++) {
          for (int i = pts_idx[idx]; i < pts_idx[idx+1]; i++) {
            (*bar_coords_idx)[ipt + 1] = (*bar_coords_idx)[ipt] + n_vtx_elt;
            ipt++;
          }
          idx++;
        }
      }

    }
  }



  *bar_coords = malloc (sizeof(double) * ((*bar_coords_idx)[n_pts]));


  /* Loop over nodal blocks */
  idx = 0;
  for (int iblock = 0; iblock < n_blocks; iblock++) {
    int id_block = blocks_id[iblock];

    /* Polyhedra */
    if (id_block >= PDM_BLOCK_ID_BLOCK_POLY3D) {
      _poly3d_block_locate (mesh_nodal_id,
                            id_block,
                            pts_idx + idx,
                            pts_coords,
                            pts_g_num,//debug
                            tolerance,
                            *distance + pts_idx[idx],
                            *bar_coords_idx + pts_idx[idx],
                            *bar_coords);
    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
      _poly2d_block_locate (mesh_nodal_id,
                            id_block,
                            pts_idx + idx,
                            pts_coords,
                            pts_g_num,//debug
                            *distance + pts_idx[idx],
                            *bar_coords_idx + pts_idx[idx],
                            *bar_coords);
    }

    /* Standard elements */
    else {
      if (DEBUG) {
        printf("nodal block %d (id = %d) --> _std_block_locate_3d (type = %d), idx = %d\n",
               iblock, id_block,
               PDM_Mesh_nodal_block_type_get(mesh_nodal_id, id_block),
               idx);
      }
      _std_block_locate_3d (mesh_nodal_id,
                            id_block,
                            pts_idx + idx,
                            pts_coords,
                            pts_g_num,//debug
                            tolerance,
                            *distance + pts_idx[idx],
                            //*projected_coords,
                            //*bar_coords_idx + idx,
                            *bar_coords + (*bar_coords_idx)[pts_idx[idx]]);

    }

    for (int ipart = 0; ipart < n_parts; ipart++) {
      n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                              id_block,
                                              ipart);
      idx += n_elt;
    }
  }

}
