/*============================================================================
 * Functions about high order meshes
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2018       ONERA

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_vtk.h"
#include "pdm_ho_location.h"
#include "pdm_box_tree.h"
#include "pdm_mesh_nodal.h"
#include "pdm_ho_seg_intersect.h"


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Optimization
 *----------------------------------------------------------------------------*/

// double *
// _scalar_product_function
// (
// double *dir_vector,
// double *inserted_point,
// double *x,
// int     mesh_dimension,
// )
// {
//   double *f;

//   if (mesh_dimension == 2) {
//     f[0] = (((dir_vector[0] - dir_vector[3])*(inserted_point[0] - x[0]))/((dir_vector[0] - dir_vector[3])*(dir_vector[0] - dir_vector[3]))) -1;
//     f[1] = (((dir_vector[1] - dir_vector[4])*(inserted_point[1] - x[1]))/((dir_vector[1] - dir_vector[4])*(dir_vector[1] - dir_vector[4]))) -1;
//   } // end 2D case

//   if (mesh_dimension == 3) {
//     f[0] = (((dir_vector[0] - dir_vector[3])*(inserted_point[0] - x[0]))/((dir_vector[0] - dir_vector[3])*(dir_vector[0] - dir_vector[3]))) -1;
//     f[1] = (((dir_vector[1] - dir_vector[4])*(inserted_point[1] - x[1]))/((dir_vector[1] - dir_vector[4])*(dir_vector[1] - dir_vector[4]))) -1;
//     f[2] = (((dir_vector[2] - dir_vector[5])*(inserted_point[2] - x[2]))/((dir_vector[2] - dir_vector[5])*(dir_vector[2] - dir_vector[5]))) -1;
//   } // end 3D case

//   return f;
// }

// double *
// _scalar_product_function_derivative
// (
// double *dir_vector,
// int     mesh_dimension,
// )
// {
//   double *df;

//   if (mesh_dimension == 2) {
//     df[0] = 1 / (dir_vector[0] - dir_vector[3]);
//     df[1] = 1 / (dir_vector[1] - dir_vector[4]);
//   } // end 2D case

//   if (mesh_dimension == 3) {
//     df[0] = 1 / (dir_vector[0] - dir_vector[3]);
//     df[1] = 1 / (dir_vector[1] - dir_vector[4]);
//     df[2] = 1 / (dir_vector[2] - dir_vector[5]);
//   } // end 3D case

//   return df;
// }

/* Apply Newton method */
static void
_newton_method
(
int     n_boxes,
double  *box_initial_points,
double  *f,
double  *df,
double  *norm,
double   eps,
double **line_box_intersection_point)
{


  *line_box_intersection_point = malloc(sizeof(double) * 3 * n_boxes);

  double *x[3];
  double *x_in[3];

  for (int ibox = 0; ibox < n_boxes; ibox++) {

    x[0] = box_initial_points[3*ibox];
    x[1] = box_initial_points[3*ibox + 1];
    x[2] = box_initial_points[3*ibox + 2];

    x_in[0] = HUGE_VAL;
    x_in[1] = HUGE_VAL;
    x_in[2] = HUGE_VAL;

    while (norm(x_in, x) > eps) {

      x_in[0] = x[0];
      x_in[1] = x[1];
      x_in[2] = x[2];

      x[0] = x[0] - f(x[0])/df(x[0]);
      x[1] = x[1] - f(x[1])/df(x[1]);
      x[2] = x[2] - f(x[2])/df(x[2]);

    } // end eps criterion satisfied

    (*line_box_intersection_point)[3*ibox] = x[0];
    (*line_box_intersection_point)[3*ibox+1] = x[1];
    (*line_box_intersection_point)[3*ibox+2] = x[2];

  }

}

/*----------------------------------------------------------------------------
 *  Get bounding boxes
 *----------------------------------------------------------------------------*/

/* Get Bezier coordinates from Lagrange coordinates for a bar (copied from pdm_t_dcube_nodal_gen.c) */
static void
_lagrange_to_bezier_bar
(
 const int  order,
 double    *lag,
 double    *bez
)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_BARHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }
  }

  else if (order == 3) {

    double f833 = 5. / 6.;
    double f333 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f833*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f333*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f333*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f833*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }

}

/* Get Bezier coordinates from Lagrange coordinates for a triangle (copied from pdm_t_dcube_nodal_gen.c) */
static void
_lagrange_to_bezier_tria
(
 const int  order,
 double    *lag,
 double    *bez
)
{
   int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get(PDM_MESH_NODAL_TRIAHO, order);

  if (order == 1) {
    memcpy (bez, lag, sizeof(double) * n_nodes * 3);
  }

  else if (order == 2) {
    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -0.5*lag[j] + 2*lag[3+j] - 0.5*lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = lag[6+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = -0.5*lag[j] + 2*lag[9+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -0.5*lag[6+j] + 2*lag[12+j] - 0.5*lag[15+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = lag[15+j];
    }
  }

  else if (order == 3) {
    double f5_6 = 5. / 6.;
    double f1_3 = 1. / 3.;

    for (int j = 0; j < 3; j++) {
      bez[j] = lag[j];
    }

    for (int j = 0; j < 3; j++) {
      bez[3+j] = -f5_6*lag[j] + 3*lag[3+j] - 1.5*lag[6+j] + f1_3*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[6+j] = f1_3*lag[j] - 1.5*lag[3+j] + 3*lag[6+j] - f5_6*lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[9+j] = lag[9+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[12+j] = -f5_6*lag[j] + 3*lag[12+j] - 1.5*lag[21+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[15+j] = f1_3*lag[j] - 0.75*lag[3+j] - 0.75*lag[6+j] + f1_3*lag[9+j] - 0.75*lag[12+j] + 4.5*lag[15+j] - 0.75*lag[18+j] - 0.75*lag[21+j] - 0.75*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[18+j] = -f5_6*lag[9+j] + 3*lag[18+j] - 1.5*lag[24+j] + f1_3*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[21+j] = f1_3*lag[j] - 1.5*lag[12+j] + 3*lag[21+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[24+j] = f1_3*lag[9+j] - 1.5*lag[18+j] + 3*lag[24+j] - f5_6*lag[27+j];
    }

    for (int j = 0; j < 3; j++) {
      bez[27+j] = lag[27+j];
    }
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for order > 3\n");
  }

}

/* Get bezier bounding box (copied from pdm_t_dcube_nodal_gen.c) */
static void
_bezier_bounding_boxes
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 PDM_g_num_t                *pvtx_ln_to_gn,
 int                        *pelt_vtx_idx,
 int                        *pelt_vtx,
 int                         n_elt,
 int                         pn_vtx,
 double                     *lagrange_coord,
 double                    **extents
)
{
  double *bezier_coord   = malloc (sizeof(double) * n_nodes * 3);
  int idx = 0;

  for (int i = 0; i < n_elt; i++) {
    double *_min = (*extents) + 6*i;
    double *_max = _min + 3;

    for (int j = 0; j < 3; j++) {
      _min[j] =  1e30;
      _max[j] = -1e30;
    }

    if (t_elt == PDM_MESH_NODAL_BAR2 ||
        t_elt == PDM_MESH_NODAL_BARHO) {
      _lagrange_to_bezier_bar (order, lagrange_coord, bezier_coord);
    }
    else if (t_elt == PDM_MESH_NODAL_TRIA3 ||
             t_elt == PDM_MESH_NODAL_TRIAHO) {
      _lagrange_to_bezier_tria (order, lagrange_coord, bezier_coord);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Only implemented for other elements in pdm_t_dcube_nodal_gen.c\n");
    }

    for (int k = 0; k < n_nodes; k++) {
      for (int j = 0; j < 3; j++) {
        _min[j] = _MIN(_min[j], bezier_coord[3*k + j]);
        _max[j] = _MAX(_max[j], bezier_coord[3*k + j]);
      }
    }
  }
}

/*----------------------------------------------------------------------------
 *  Intersect bounding boxes
 *----------------------------------------------------------------------------*/

/* Get the intersection between a face and a line */
static int
_face_line_intersect
(
double                     *face_points,
double                     *line_points,
double                    **intersection_point
)
{
  /*
   * WARNING: - A needs to be "smaller"  than B
   */

  int out = 0;
  int *x_int[3];
  int x_min, y_min, z_min, x_max, y_max, z_max;
  int line0_side, line1_side;

  double *face_coord0[3];
  double *face_coord1[3];
  double *face_coord2[3];

  double *line_coord0[3];
  double *line_coord1[3];

  double *n[3];
  double *face01[3];
  double *face02[3];
  double *face0line0[3];
  double *face0line1[3];

  face_coord0[0] = face_points[0];
  face_coord0[1] = face_points[1];
  face_coord0[2] = face_points[2];

  face_coord1[0] = face_points[3];
  face_coord1[1] = face_points[4];
  face_coord1[2] = face_points[5];

  face_coord2[0] = face_points[6];
  face_coord2[1] = face_points[7];
  face_coord2[2] = face_points[8];

  line_coord0[0] = line_points[0];
  line_coord0[1] = line_points[1];
  line_coord0[2] = line_points[2];

  line_coord1[0] = line_points[3];
  line_coord1[1] = line_points[4];
  line_coord1[2] = line_points[5];

  face01[0] = face_coord1[0] - face_coord0[0];
  face01[1] = face_coord1[1] - face_coord0[1];
  face01[2] = face_coord1[2] - face_coord0[2];

  face02[0] = face_coord2[0] - face_coord0[0];
  face02[1] = face_coord2[1] - face_coord0[1];
  face02[2] = face_coord2[2] - face_coord0[2];

  face0line0[0] = line_coord0[0] - face_coord0[0];
  face0line0[1] = line_coord0[1] - face_coord0[1];
  face0line0[2] = line_coord0[2] - face_coord0[2];

  face0line1[0] = line_coord1[0] - face_coord0[0];
  face0line1[1] = line_coord1[1] - face_coord0[1];
  face0line1[2] = line_coord1[2] - face_coord0[2];

  PDM_CROSS_PRODUCT(n, face01, face02);

  int d = -(n[0]*face_coord0[0] + n[1]*face_coord0[1] + n[2]*face_coord0[2]);

  line0_side = PDM_DOT_PRODUCT(n, line_coord0);
  line1_side = PDM_DOT_PRODUCT(n, line_coord1);

  if (line0_side * line1_side < 0) {
    *intersection_point = malloc(sizeof(double) * 3);
    out = 1;
    int t = ((n[0] + 1) * face_coord0[0] + (n[1] + 1) *  face_coord0[1] +  (n[2] + 1) *  face_coord0[2] + d);
    t -= (line_coord0[0] + line_coord0[1] + line_coord0[2]);
    t /= ((line_coord0[0] - line_coord1[0]) + (line_coord0[1] - line_coord1[1]) + (line_coord0[2] - line_coord1[2]));

    (*intersection_point)[0] = line_coord0[0] + t * (line_coord0[0] - line_coord1[0]);
    (*intersection_point)[1] = line_coord0[1] + t * (line_coord0[1] - line_coord1[1]);
    (*intersection_point)[2] = line_coord0[2] + t * (line_coord0[2] - line_coord1[2]);
  }
  return out;
}

/* Get points at which line intersects bounding box */
static void
_ho_bounding_box_line_intersect_points_get
(
 const int        n_line,                   // number of direction lines for the set of point to project
 double          *line_coords,              // extrema coordinates of those lines
 int             *line_boxes_idx,           // boxes associated to a given line
 double          *extents,                  // extents of the background mesh element bounding boxes
 double         **box_line_intersect_points // for each box the intersection points with a given line (2 per box)
)
{
  *box_line_intersect_points     = malloc(sizeof(double) * (*line_boxes_idx)[n_line] * 6);
  double *_box_line_intersect_points = *box_line_intersect_points;

  int count_plane_intersect;
  double x_min, y_min, z_min, x_max, y_max, z_max;

  int *face_points[9];
  int *line_points[6];


  for (int iline = 0; iline < n_line; iline++) {
    for (int ibox = line_boxes_idx[iline]; ibox < line_boxes_idx[iline+1]; ibox++) {

      x_min = extents[6*ibox];
      y_min = extents[6*ibox+1];
      z_min = extents[6*ibox+2];
      x_max = extents[6*ibox+3];
      y_max = extents[6*ibox+4];
      z_max = extents[6*ibox+5];

      line_point = &line_coord[6*ibox];

      // x_min side

      face_points[0] = x_min;
      face_points[1] = y_min;
      face_points[2] = z_min;
      face_points[3] = x_min;
      face_points[4] = y_min;
      face_points[5] = z_max;
      face_points[6] = x_max;
      face_points[7] = y_min;
      face_points[8] = z_min;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

      // y_min side

      face_points[0] = x_min;
      face_points[1] = y_min;
      face_points[2] = z_min;
      face_points[3] = x_min;
      face_points[4] = y_max;
      face_points[5] = z_min;
      face_points[6] = x_min;
      face_points[7] = y_min;
      face_points[8] = z_max;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

      // z_min side

      face_points[0] = x_min;
      face_points[1] = y_min;
      face_points[2] = z_min;
      face_points[3] = x_min;
      face_points[4] = y_min;
      face_points[5] = z_max;
      face_points[6] = x_min;
      face_points[7] = y_max;
      face_points[8] = z_min;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

      // x_max side

      face_points[0] = x_min;
      face_points[1] = y_max;
      face_points[2] = z_min;
      face_points[3] = x_min;
      face_points[4] = y_max;
      face_points[5] = z_max;
      face_points[6] = x_max;
      face_points[7] = y_max;
      face_points[8] = z_min;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

      // y_max side

      face_points[0] = x_max;
      face_points[1] = y_min;
      face_points[2] = z_min;
      face_points[3] = x_max;
      face_points[4] = y_min;
      face_points[5] = z_max;
      face_points[6] = x_max;
      face_points[7] = y_max;
      face_points[8] = z_min;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

      // z_max side

      face_points[0] = x_min;
      face_points[1] = y_min;
      face_points[2] = z_max;
      face_points[3] = x_max;
      face_points[4] = y_min;
      face_points[5] = z_max;
      face_points[6] = x_min;
      face_points[7] = y_max;
      face_points[8] = z_max;

      double **intersection_point = NULL;

      int out _face_line_intersect(face_coord, line_coord, intersection_point);

      if (out == 1) {

        if ((x_min <= x_int[0] && x_max >= x_int[0] ) && (y_min <= x_int[1] && y_max >= x_int[1] ) && (z_min <= x_int[2] && z_max >= x_int[2] )) {
          _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = (*intersection_point)[0];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = (*intersection_point)[1];
          _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = (*intersection_point)[2];
          count_plane_intersect++;
        } // end if on the face

        free(*intersection_point);
        free(intersection_point);

      } // end if there is an intersection

    } // end loop on boxes crossed by iline
  } // end loop on lines

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Get bezier bounding boxes of ho elements
 *
 * \param [in]   type              Element type
 * \param [out]  t                 Parametric coordinates of the projection on the segment
 *
 */

void
PDM_ho_seg_intersect_boxes_get
(
)
{
// int **line_boxes_idx  = NULL;
// int **line_boxes_lnum = NULL;

// PDM_box_tree_intersect_lines_boxes(bt,
//                                    i_copied_rank,
//                                    n_line,
//                                    line_coord,
//                                    line_boxes_idx,
//                                    line_boxes_lnum);
}

/**
 * \brief Determine intersection between P1 element and segment from line intersecting ho face bounding box
 *
 * \param [in]   t_elt                           Element type
 * \param [in]   n_back_face_to_intersect        Number of background mesh faces to intersect
 * \param [in]   back_face_to_intersect_ln_to_gn Set of background mesh faces to intersect
 * \param [in]   back_face_vtx                   Background face -> vertex connectivity
 * \param [in]   back_vtx_coord                  Coordinates of vertices of the background mesh
 * \param [in]   back_face_line_idx              Index of background face -> line connectivity
 * \param [in]   back_face_line                  Background face -> line connectivity
 * \param [in]   line_boxes_idx                  Index of line -> box connectivity
 * \param [in]   box_line_intersect_points       Line -> box intersection points connectivity
 * \param [out]  newton_initial_point            Initial point for Newthon Method
 *
 */

void
PDM_ho_seg_intersect_P1_line
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 int                         n_back_face_to_intersect,
 PDM_g_num_t                *back_face_to_intersect_ln_to_gn,
 PDM_g_num_t                *back_face_vtx,
 double                     *back_vtx_coord,
 int                        *back_face_line_idx,
 PDM_g_num_t                *back_face_line,
 int                        *line_boxes_idx,
 double                     *box_line_intersect_points,
 double                    **newton_initial_point
 )
{

  /*
   * WARNING: - Make sure intersection is in triangle
   *          - If not implement moving from one triangle to the other
   */

  switch (t_elt)
  {
  case PDM_MESH_NODAL_TRIA3:

    *newton_initial_point = malloc(sizeof(double) * 3 * line_boxes_idx[back_face_line_idx[n_back_face_to_intersect]]);

    // Get intersection between line and P1 approximation (to get starting point)

    int face_id;
    int line_id;
    int vtx_id1, vtx_id2, vtx_id3;

    int *face_coord[9];
    int *line_coord[6];

    for (int iface = 0; iface < n_back_face_to_intersect; iface++) {

      face_id = PDM_ABS(back_face_to_intersect_ln_to_gn[iface])-1;
      vtx_id1 = back_face_vtx[3*face_id]-1;
      vtx_id2 = back_face_vtx[3*face_id+1]-1;
      vtx_id3 = back_face_vtx[3*face_id+2]-1;
      line_id = back_face_line[back_face_line_idx[face_id]];

      face_coord[0] = back_vtx_coord[3*vtx_id1];
      face_coord[1] = back_vtx_coord[3*vtx_id1+1];
      face_coord[2] = back_vtx_coord[3*vtx_id1+2];
      face_coord[3] = back_vtx_coord[3*vtx_id2];
      face_coord[4] = back_vtx_coord[3*vtx_id2+1];
      face_coord[5] = back_vtx_coord[3*vtx_id2+2];
      face_coord[6] = back_vtx_coord[3*vtx_id3];
      face_coord[7] = back_vtx_coord[3*vtx_id3+1];
      face_coord[8] = back_vtx_coord[3*vtx_id3+2];

      for (int ibox = line_boxes_idx[line_id]; ibox < line_boxes_idx[line_id+1]; ibox++ ) {

        line_coord[0] = box_line_intersect_points[6*ibox];
        line_coord[1] = box_line_intersect_points[6*ibox+1];
        line_coord[2] = box_line_intersect_points[6*ibox+2];
        line_coord[3] = box_line_intersect_points[6*ibox+3];
        line_coord[4] = box_line_intersect_points[6*ibox+4];
        line_coord[5] = box_line_intersect_points[6*ibox+5];

        double **intersection_point = NULL;

        int out _face_line_intersect(face_coord, line_coord, intersection_point);

        if (out == 1) {

          (*newton_initial_point)[3*ibox    ] = (**intersection_point)[0];
          (*newton_initial_point)[3*ibox + 1] = (**intersection_point)[1];
          (*newton_initial_point)[3*ibox + 2] = (**intersection_point)[2];

        } // end if there is an intersection

      } // end loop on boxes
    } // end loop on faces

    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Not implement for elements other than TRIA3\n");
    break;
  }
}

/**
 * \brief Compute intersection between line and ho element using Newton Method
 *
 * \param [in]   type                            Element type
 * \param [in]   n_line                          Number of lines
 * \param [in]   line_coords                     Coordinates of lines
 * \param [in]   n_back_face_to_intersect        Number of background mesh faces to intersect
 * \param [in]   back_face_to_intersect_ln_to_gn Set of background mesh faces to intersect
 * \param [in]   back_face_vtx                   Background face -> vertex connectivity
 * \param [in]   back_vtx_coord                  Coordinates of vertices of the background mesh
 * \param [in]   back_face_line_idx              Index of background face -> line connectivity
 * \param [in]   back_face_line                  Background face -> line connectivity
 * \param [out]  line_box_intersection_point     Line->Box->Point connectivity
 *
 */

void
PDM_ho_seg_intersect_compute
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   n_line,
 double                     *line_coords,
 int                         n_back_face_to_intersect,
 PDM_g_num_t                *back_face_to_intersect_ln_to_gn,
 PDM_g_num_t                *back_face_vtx,
 double                     *back_vtx_coord,
 int                        *back_face_line_idx,
 PDM_g_num_t                *back_face_line,
 double                    **line_box_intersection_point
)
{
  // Get boxes from background faces to intersect

  /*
   * Out of PDM_ho_seg_intersect_boxes_get : - line_boxes_idx
   *                                         - extents (equivalent to line_boxes_extents_coords)
   */

  // Get the intersection between the lines and those boxes

  double **box_line_intersect_points = NULL;

  _ho_bounding_box_line_intersect_points_get(n_line,
                                             line_coords,
                                             line_boxes_idx,
                                             extents,
                                             box_line_intersect_points);

  // Get points at which line intersects box (to have line parametrization)

  double **newton_initial_point = NULL;

  PDM_ho_seg_intersect_P1_line(t_elt,
                               n_back_face_to_intersect,
                               back_face_to_intersect_ln_to_gn,
                               back_face_vtx,
                               back_vtx_coord,
                               back_face_line_idx,
                               back_face_line,
                               line_boxes_idx,
                               *box_line_intersect_points,
                               **newton_initial_point);

  // Operate intersection

  _newton_method(box_initial_points,
                 box_ho_element_function,
                 box_ho_element_function_derivative,
                 line_box_intersection_point);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
