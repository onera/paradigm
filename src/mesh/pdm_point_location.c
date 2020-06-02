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

#if 0
void
_poly3d_block_locate
(
 const int                            mesh_nodal_id,
 const int                            id_block,
 const int                           *pts_in_extents_idx,
 const PDM_g_num_t                   *pts_in_extents_g_num,
 const double                        *pts_in_extents_coords,
 double                               tolerance,
 double                             **pts_in_extents_distance,
 int                                **pts_in_extents_bar_coord_idx,
 double                             **pts_in_extents_bar_coord
 )
{
  /*PDM_Mesh_nodal_t *mesh = (PDM_Mesh_nodal_t *) PDM_Handles_get (mesh_handles, mesh_nodal_id);

    if (mesh == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
    }

    PDM_Mesh_nodal_block_poly3d_t *block =
    (PDM_Mesh_nodal_block_poly3d_t *) PDM_Handles_get (mesh->blocks_poly3d, id_block);

    if (block == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard block identifier\n");
    }*/

  /* Loop over parts */
  int n_part = PDM_Mesh_nodal_n_part_get (mesh_nodal_id);
  for (int ipart = 0; ipart < n_part; ipart++) {
    // Get vertex coordinates
    const double *vtx_coords = PDM_Mesh_nodal_vertices_get (mesh_nodal_id, ipart);


  }
}
#endif



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
 const double        tetra_coords[4][3],
 const double        point_coords[],
 const int           n_points,
 const double        tolerance,
 float              *distance,
 double             *bary_coords
 )
{
  //const double upper_bound = 0.5 + tolerance;

  double v[3][3];
  for (int ivtx = 0; ivtx < 3; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      v[ivtx][idim] = tetra_coords[ivtx+1][idim] - tetra_coords[0][idim];
    }
  }

  //double vol6 = PDM_ABS (_determinant_3x3 (v[0], v[1], v[2]));
  double vol6 = fabs (v[0][0] * (v[1][1]*v[2][2] - v[1][2]*v[2][1]) +
                      v[0][1] * (v[1][2]*v[2][0] - v[1][0]*v[2][2]) +
                      v[0][2] * (v[1][0]*v[2][1] - v[1][1]*v[2][0]));
  if (vol6 < _epsilon_denom){
    PDM_printf("warning _locate_in_tetra : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    PDM_printf_flush();
    return;
  }
  vol6 = 1./vol6;

  double r[3][3];
  for (int i = 0; i < 3; i++) {
    int j = (i + 1) % 3;
    int k = (i + 2) % 3;

    PDM_CROSS_PRODUCT (r[i], v[j], v[k]);

    for (int idim = 0; idim < 3; idim++) {
      r[i][idim] *= vol6;
    }
  }

  /* Loop over points */
  double vp0[3];
  for (int ipt = 0; ipt < n_points; ipt++) {
    const double *p = point_coords + 3*ipt;
    double *bco = bary_coords + 4*ipt;

    for (int idim = 0; idim < 3; idim++) {
      vp0[idim] = tetra_coords[0][idim] - p[idim];
    }

    /* Compute barycentric coordinates of current point in tetrahedron */
    bco[0] = 1.;
    for (int i = 0; i < 3; i++) {
      bco[i+1] = PDM_DOT_PRODUCT (vp0, r[i]);
      bco[0] -= bco[i+1];
    }

    /* Compute 'distance' from current point to tetrahedron */
    double max_dist = -1.;
    for (int ivtx = 0; ivtx < 4; ivtx++) {
      double dist = PDM_ABS (bco[ivtx] - 0.5);

      if (dist > max_dist) {
        max_dist = dist;
      }
    }

    /* ... */
    //if (max_dist > -0.25 && max_dist < upper_bound)
    distance[ipt] = (float) (2 * max_dist);

  }

}

//-->>
void
PDM_locate_points_in_tetra (const double  vtx_xyz[12],
                            const int     n_pts,
                            const double  pts_xyz[],
                            float        *distance,
                            double       *bary_coords)
{
  const double tolerance = 1e-6;

  double tetra_coords[4][3];
  for (int ivtx = 0; ivtx < 4; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      tetra_coords[ivtx][idim] = vtx_xyz[3*ivtx + idim];
    }
  }

  _locate_in_tetra (tetra_coords,
                    pts_xyz,
                    n_pts,
                    tolerance,
                    distance,
                    bary_coords);
}
//<<--

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

  default:
    /*PDM_error(__FILE__, __LINE__, 0,
      "_compute_shapef: unhandled element type %s\n",
      fvmc_element_type_name[elt_type]);*/
    PDM_error (__FILE__, __LINE__, 0, "Wrong element type\n");

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

static PDM_bool_t
_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[],
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





static void
_locate_in_cell_3d
(
 const PDM_Mesh_nodal_elt_t  elt_type,
 const PDM_l_num_t           element_vertex_num[],
 const PDM_l_num_t          *parent_vertex_num,
 const double                vertex_coords[],
 const int                   n_pts,
 const double                pts_coords[],
 const PDM_g_num_t           pts_g_num[],//
 const double                tolerance,
 float                      *distance,
 double                     *bar_coords
 )
{
  const int DEBUG = 0;
  PDM_l_num_t coord_idx;

  double uvw[3], dist, shapef[8], max_dist;
  double _vertex_coords[8][3];

  const int order = 1;
  const double upper_bound = 1. + 2.*tolerance;

  const int n_vertices = PDM_Mesh_nodal_n_vertices_element (elt_type, order);





  /* Initialize local element coordinates copy */
  //printf("vtx = \n");
  for (int ivtx = 0; ivtx < n_vertices; ivtx++) {

    if (parent_vertex_num == NULL) {
      coord_idx = element_vertex_num[ivtx] - 1;
    } else {
      coord_idx = parent_vertex_num[element_vertex_num[ivtx] - 1] - 1;
    }
    //printf("%d ", coord_idx);

    for (int idim = 0; idim < 3; idim++) {
      _vertex_coords[ivtx][idim] = vertex_coords[coord_idx * 3 + idim];
      //printf("%f ", _vertex_coords[ivtx][idim]);
    }
    //printf("\n");

  }

  /*printf("%f < x < %f, %f < y < %f, %f < z = %f\n",
         _vertex_coords[3][0], _vertex_coords[5][0],
         _vertex_coords[3][1], _vertex_coords[5][1],
         _vertex_coords[3][2], _vertex_coords[5][2]);*/


  if (elt_type == PDM_MESH_NODAL_TETRA4) {
    /* Shape functions may be computed directly with tetrahedra */
    _locate_in_tetra (_vertex_coords,
                      pts_coords,
                      n_pts,
                      tolerance,
                      distance,
                      bar_coords);

  } else {
    /* For cell shapes other than tetrahedra, find shape functions iteratively */
    for (int ipt = 0; ipt < n_pts; ipt++) {
      /* Check vertices (To not have a problem with pyramids) */
      const double *_pt = pts_coords + 3 * ipt;
      if (DEBUG) {
        printf("pt %d (%ld) : (%f %f %f)\n",
                        ipt, pts_g_num[ipt], _pt[0], _pt[1], _pt[2]);
      }
      double *_bco = bar_coords + n_vertices * ipt;
      PDM_bool_t on_vtx = PDM_FALSE;

      for (int ivtx = 0; ivtx < n_vertices; ivtx++) {

        double v[3] = {_vertex_coords[ivtx][0] - _pt[0],
                       _vertex_coords[ivtx][1] - _pt[1],
                       _vertex_coords[ivtx][2] - _pt[2]};

        double _dist = PDM_MODULE (v);

        if (_dist < 1e-6 * tolerance) {
          distance[ipt] = 0.;
          for (int i = 0; i < n_vertices; i++) {
            _bco[i] = 0.;
          }
          _bco[ivtx] = 1.;

          on_vtx = PDM_TRUE;
          break;
        }
      }

      if (on_vtx == PDM_FALSE) {

        if (_compute_uvw(elt_type,
                         _pt,
                         _vertex_coords,
                         tolerance,
                         uvw) == PDM_TRUE) {
          max_dist = -1.0;

          /* For hexahedra, no need to compute shape functions, as
             the 3 parametric coordinates are simpler to use */

          if (elt_type == PDM_MESH_NODAL_HEXA8) {

            for (int i = 0; i < 3; i++){
              dist = 2.*PDM_ABS(uvw[i] - 0.5);

              if (max_dist < dist) {
                max_dist = dist;
              }
            }

            //-->>
            /* Barycentric coordinates */
            _compute_shapef_3d (elt_type,
                                uvw,
                                _bco,
                                NULL);
            if (DEBUG) {
              double sum = 0;
              printf(" bco = ");
              for (int ivtx = 0; ivtx < n_vertices; ivtx++) {
                printf("%f ", _bco[ivtx]);
                sum += _bco[ivtx];
              }
              printf("  (sum = %f)\n", sum);
            }
            //<<--

          }

          /* For pyramids ands prisms, we need to compute shape functions */

          else {

            _compute_shapef_3d (elt_type, uvw, shapef, NULL);

            for (int ivtx = 0; ivtx < n_vertices; ivtx++){
              dist = 2.*PDM_ABS(shapef[ivtx] - 0.5);

              if (max_dist < dist) {
                max_dist = dist;
              }

              _bco[ivtx] = shapef[ivtx];
            }

          }

          /* For all element types, update location and distance arrays */

          /*if (max_dist > -0.5 && max_dist < upper_bound) {
            distance[ipt] = (float) max_dist;
            printf("dist = %f\n\n", distance[ipt]);
          } else {
            distance[ipt] = HUGE_VAL;
            }*/
          distance[ipt] = (float) max_dist;
          if (DEBUG) {
            printf("dist = %f\n\n", distance[ipt]);
          }

        } else {

          distance[ipt] = 1.e12;
          //bar_coords?
          for (int ivtx = 0; ivtx < n_vertices; ivtx++) {
            _bco[ivtx] = -1.;
          }

        }

      }

    }
  }
}







static void
_std_block_locate_3d
(
 const int           mesh_nodal_id,
 const int           id_block,
 const int          *pts_idx,
 const double       *pts_coords,
 const PDM_g_num_t  *pts_g_num,//
 double              tolerance,
 float              *distance,
 // double             *projected_coords,// high-order
 int                *bar_coords_idx,
 double             *bar_coords
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

  } else if (elt_type == PDM_MESH_NODAL_BAR2) {
    entity_dim = 1;

  } else if (elt_type == PDM_MESH_NODAL_TRIA3 ||
             elt_type == PDM_MESH_NODAL_QUAD4) {
    entity_dim = 2;

  } else if (elt_type == PDM_MESH_NODAL_TETRA4   ||
             elt_type == PDM_MESH_NODAL_PYRAMID5 ||
             elt_type == PDM_MESH_NODAL_PRISM6   ||
             elt_type == PDM_MESH_NODAL_HEXA8) {
    entity_dim = 3;

  } else {
    PDM_error (__FILE__, __LINE__, 0, "Wrong standard element type\n");
  }


  int n_vtx = PDM_Mesh_nodal_n_vertices_element (elt_type,
                                                 order);


  /* Loop over parts */
  int idx = -1;
  int ipt = 0;
  for (int ipart = 0; ipart < n_parts; ipart++) {

    n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                            id_block,
                                            ipart);

    const double *vertex_coords = PDM_Mesh_nodal_vertices_get (mesh_nodal_id,
                                                               ipart);
    PDM_Mesh_nodal_block_std_get (mesh_nodal_id,
                                  id_block,
                                  ipart,
                                  &element_vertex_num);

    /* Loop over elements */
    for (int ielt = 0; ielt < n_elt; ielt++) {
      idx++;

      n_pts = pts_idx[idx + 1] - pts_idx[idx];

      if (n_pts < 1) {
        continue;
      }

      /* Volume elements */
      if (entity_dim == 3) {

        if (order == 1) {
          /* First-order elements */
          if (DEBUG) {
            printf("id_block = %d, ipart = %d, ielt = %d --> _locate_in_cell_3d (ipt = %d)\n",
                   id_block, ipart, ielt, ipt);
          }
          _locate_in_cell_3d (elt_type,
                              element_vertex_num + ielt * n_vtx,
                              parent_vertex_num,
                              vertex_coords,
                              n_pts,
                              pts_coords + ipt * 3,
                              pts_g_num + ipt,
                              tolerance,
                              distance + ipt,
                              bar_coords + ipt * n_vtx);

        } else {
          /* High-order elements */
          assert (0 == 1);
        }

      }


      /* Surface elements */
      else if (entity_dim == 2) {

        if (order == 1) {
          /* First-order elements */
          /*
            _locate_on_face_3d...
          */
          assert (0 == 1);

        } else {
          /* High-order elements */
          assert (0 == 1);
        }

      }


      /* Linear elements */
      else if (entity_dim == 1) {

        if (order == 1) {
          /* First-order elements */
          /*
            _locate_on_edge_3d...
          */
          assert (0 == 1);

        } else {
          /* High-order elements */
          assert (0 == 1);
        }

      }


      /* Points */
      else { // entity_dim == 0
        /*...*/
        assert (0 == 1);

      }


      ipt += n_pts;
    } // Loop over elements

  } // Loop over parts

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
 const PDM_g_num_t   *pts_g_num,//
 const double         tolerance,
 int                  base_element_num,//?
 float              **distance,
 double             **projected_coords,//high-order
 int                **bar_coords_idx,
 double             **bar_coords
 )
{
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

        PDM_l_num_t  n_face;
        PDM_l_num_t *facvtx_idx = NULL;
        PDM_l_num_t *facvtx = NULL;
        PDM_l_num_t *cellfac_idx = NULL;
        PDM_l_num_t *cellfac = NULL;
        PDM_Mesh_nodal_block_poly3d_get (mesh_nodal_id,
                                         id_block,
                                         ipart,
                                         &n_face,
                                         &facvtx_idx,
                                         &facvtx,
                                         &cellfac_idx,
                                         &cellfac);

        for (int ielt = 0; ielt < n_elt; ielt++) {
          n_vtx_elt = 0;
          for (int i = cellfac_idx[ielt]; i < cellfac_idx[ielt+1]; i++) {
            int iface = cellfac[i] - 1;
            n_vtx_elt += facvtx_idx[iface + 1] - facvtx_idx[iface];
          }
          n_vtx_elt /= 2;

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
      /*_poly3d_block_locate (mesh_nodal_id,
        id_block,
        ...);*/
      assert (0);

    }

    /* Polygons */
    else if (id_block >= PDM_BLOCK_ID_BLOCK_POLY2D) {
      assert (0);

    }

    /* Standard elements */
    else {
      _std_block_locate_3d (mesh_nodal_id,
                            id_block,
                            pts_idx + idx,
                            pts_coords,
                            pts_g_num,
                            tolerance,
                            *distance,
                            //*projected_coords,
                            *bar_coords_idx + idx,
                            *bar_coords);

    }

    for (int ipart = 0; ipart < n_parts; ipart++) {
      n_elt = PDM_Mesh_nodal_block_n_elt_get (mesh_nodal_id,
                                                  id_block,
                                                  ipart);
      idx += n_elt;
    }
  }

}
