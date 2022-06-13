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
 * Type definitions
 *============================================================================*/

#define _MIN(a,b) ((a) < (b) ? (a) : (b))
#define _MAX(a,b) ((a) > (b) ? (a) : (b))

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

// double *
// _optimum_scalar_product
// (
// double *x,
// double *dir_vector,
// double *inserted_point,
// int     mesh_dimension
// )
// {
//   double *x_in;
//   double *f;
//   double *df;
//   while (norm(x_in,x) > eps) {
//     x_in = x;
//     f    = _scalar_product_function(dir_vector, inserted_point, x, mesh_dimension);
//     df   = _scalar_product_function_derivative(dir_vector, mesh_dimension);

//     if (mesh_dimension == 2) {
//       x[0] = x[0] - f[0]/df[0];
//       x[1] = x[1] - f[1]/df[1];
//     } // end 2D case

//     if (mesh_dimension == 3) {
//       x[0] = x[0] - f[0]/df[0];
//       x[1] = x[1] - f[1]/df[1];
//       x[2] = x[2] - f[2]/df[2];
//     } // end 3D case
//   }
// }

/*----------------------------------------------------------------------------
 *  Bounding boxes
 *----------------------------------------------------------------------------*/

/* Get Bezier coordinates from Lagrange coordinates for a bar */
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

/* Get Bezier coordinates from Lagrange coordinates for a triangle */
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

/* Get bezier bounding box */
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

/* Get points at which line intersects bounding box */
static void
_ho_bounding_box_line_intersect_points_get
(
 const int        n_line,
 double          *line_coords,
 int             *line_boxes_idx,
 double          *extents,
 int            **box_line_intersect_points_idx
)
{
  *box_line_intersect_points     = malloc(sizeof(double) * (*line_boxes_idx)[n_line] * 6);
  double *_box_line_intersect_points = *box_line_intersect_points;

  double *n[3];
  double *side_seg_start[3];
  double *side_seg_end[3];

  double start_side;
  double end_side;

  int count_plane_intersect;
  double t;
  double u,v,w;
  double x_min, y_min, z_min, x_max, y_max, z_max;
  double x_a, y_a, z_a, x_b, y_b, z_b;


  for (int iline = 0; iline < n_line; iline++) {
    for (int ibox = line_boxes_idx[iline]; ibox < line_boxes_idx[iline+1]; ibox++) {

      x_min = extents[6*ibox];
      y_min = extents[6*ibox+1];
      z_min = extents[6*ibox+2];
      x_max = extents[6*ibox+3];
      y_max = extents[6*ibox+4];
      z_max = extents[6*ibox+5];

      x_a = line_coord[6*ibox];
      y_a = line_coord[6*ibox+1];
      z_a = line_coord[6*ibox+2];
      x_b = line_coord[6*ibox+3];
      y_b = line_coord[6*ibox+4];
      z_b = line_coord[6*ibox+5];

      u = x_b - x_a;
      v = y_b - y_a;
      w = z_b - z_a;

      side_seg_start[0] = x_min - x_a;
      side_seg_start[1] = y_min - y_a;
      side_seg_start[2] = z_min - z_a;

      side_seg_end[0] = x_b - x_min;
      side_seg_end[1] = y_b - y_min;
      side_seg_end[2] = z_b - z_min;

      count_plane_intersect = 0;

      /*
       * WARNING: - A needs to be "smaller"  than B
       *          - I might need to replace 1 by value on the face for the normal
       *          - Check if intersection really is on the face before adding it
       */

      // x_min side

      n[0] = 0;
      n[1] = -1;
      n[2] = 0;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_min + z_min + y_min) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

      // y_min side

      n[0] = -1;
      n[1] = 0;
      n[2] = 0;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_min + z_min + y_min) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

      // z_min side

      n[0] = 0;
      n[1] = 0;
      n[2] = -1;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_min + z_min + y_min) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

      // x_max side

      n[0] = 0;
      n[1] = 1;
      n[2] = 0;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_min + z_min + y_max) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

      // y_max side

      n[0] = 1;
      n[1] = 0;
      n[2] = 0;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_max + z_min + y_min) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

      // z_max side

      n[0] = 0;
      n[1] = 0;
      n[2] = 1;

      start_side = PDM_DOT_PRODUCT(n, side_seg_start);
      end_side = PDM_DOT_PRODUCT(n, side_seg_end);

      if (start_side * end_side < 0) {

        t = ((x_min + z_max + y_min) + (x_a + y_a + z_a)) / (u + v + w);

        _box_line_intersect_points[ibox+ (3*count_plane_intersect)]   = x_a + t * u;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+1)] = y_a + t * v;
        _box_line_intersect_points[ibox+ (3*count_plane_intersect+2)] = z_a + t * w;
        count_plane_intersect++;

      } // == 0 point on face, > 0 no intersection

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
 * \brief Determine intersection between a segment and a high order element
 *
 * \param [in]   type              Element type
 * \param [in]   order             Element order
 * \param [in]   n_nodes           Number of nodes
 * \param [in]   nodes_coords      Coordinates of the nodes (size = 3 * \ref n_nodes)
 * \param [in]   seg_coords        Coordinates of the point to locate (size = 3 * 2)
 * \param [out]  projected_coords  Coordinates of the projection on the element (size = 3)
 * \param [out]  uvw               Parametric coordinates of the projection on the element
 * \param [out]  t                 Parametric coordinates of the projection on the segment
 *
 */

void
PDM_ho_seg_intersect_element
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 const double               *nodes_coords,
 const double               *seg_coords,
 const int                   n_face,
 PDM_g_num_t                *face_vtx,
 double                     *vtx_coord,
 double                     *newton_initial_point
 )
{
  switch (t_elt)
  {
  case PDM_MESH_NODAL_TRIA3:

  // Get intersection between line and P1 approximation (to get starting point)

  int vtx_id1, vtx_id2, vtx_id3;

  double *start[3];
  double *end[3];

  double *u = NULL;
  double *v = NULL;

  double *n[3];
  double *AB[3];
  double *AC[3];
  double *side_seg_start[3];
  double *side_seg_end[3];

  double d, t;

  PDM_line_intersect_t found;

  for (int iface = 0; iface < n_face; iface++-) {

    vtx_id1 = PDM_ABS(face_vtx[3*iface])  -1;
    vtx_id2 = PDM_ABS(face_vtx[3*iface+1])-1;
    vtx_id3 = PDM_ABS(face_vtx[3*iface+2])-1;

    AB[0] = vtx_coord[3*vtx_id1] - vtx_coord[3*vtx_id2];
    AB[1] = vtx_coord[3*vtx_id1+1] - vtx_coord[3*vtx_id2+1];
    AB[2] = vtx_coord[3*vtx_id1+2] - vtx_coord[3*vtx_id2+2];

    AC[0] = vtx_coord[3*vtx_id1] - vtx_coord[3*vtx_id3];
    AC[1] = vtx_coord[3*vtx_id1+1] - vtx_coord[3*vtx_id3+1];
    AC[2] = vtx_coord[3*vtx_id1+2] - vtx_coord[3*vtx_id3+2];

    PDM_CROSS_PRODUCT(n, AB, AC);

    d = -(n[0] * vtx_coord[3*vtx_id1] + n[1] * vtx_coord[3*vtx_id1+1] + n[2] * vtx_coord[3*vtx_id1+2]);

    side_seg_start[0] = vtx_coord[3*vtx_id1]   - box_in[0];
    side_seg_start[1] = vtx_coord[3*vtx_id1+1] - box_in[1];
    side_seg_start[2] = vtx_coord[3*vtx_id1+2] - box_in[2];

    side_seg_end[0] = vtx_coord[3*vtx_id1]   - box_out[0];
    side_seg_end[1] = vtx_coord[3*vtx_id1+1] - box_out[1];
    side_seg_end[2] = vtx_coord[3*vtx_id1+2] - box_out[2];

    start_side = PDM_DOT_PRODUCT(n, side_seg_start);
    end_side = PDM_DOT_PRODUCT(n, side_seg_end);

    if (start_side * end_side < 0) {

      t = (vtx_coord[3*vtx_id1] + vtx_coord[3*vtx_id1+1] + vtx_coord[3*vtx_id1+2] + n[0]*vtx_coord[3*vtx_id1] + n[1]*vtx_coord[3*vtx_id1+1] + n[2]*vtx_coord[3*vtx_id1+2] + d);
      t -= (box_in[0] + box_in[1] + box_in[2]);
      t /= ((box_in[0] - box_out[0]) + (box_in[1] - box_out[1]) + (box_in[2] - box_out[2]));

      newton_initial_point[3*iface    ] = box_in[0] + t * (box_in[0] - box_out[0]);
      newton_initial_point[3*iface + 1] = box_in[1] + t * (box_in[1] - box_out[1]);
      newton_initial_point[3*iface + 2] = box_in[2] + t * (box_in[2] - box_out[2]);

    } // == 0 point on face, > 0 no intersection

  } // end loop on faces

    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Not implement for elements other than TRIA3\n");
    break;
  }
}

/**
 * \brief Get all intersections for a set op points
 *
 * \param [in]   type              Element type
 * \param [out]  t                 Parametric coordinates of the projection on the segment
 *
 */

void
PDM_ho_seg_intersect_compute
(
)
{
  // Get mesh boxes

   // ---> PDM_ho_seg_intersect_boxes_get

  // Get points at which line intersects box (to have line parametrization)

  // ---> PDM_ho_seg_intersect_element

  // Operate intersection

  // ---> Newton method

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
