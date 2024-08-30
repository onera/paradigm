/*
  This file is part of the CWIPI library.

  Copyright (C) 2012  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_priv.h"
#include "pdm_triangle.h"
#include "pdm_triangulate.h"
#include "pdm_polygon.h"
#include "pdm_plane.h"
#include "pdm_geom_elem.h"
#include "pdm_hash_tab.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_order.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

enum {false, true};

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double GEOM_EPS_MIN  = 1e-30; /*!< Minimum value allowed for geometric computation */

static const double GEOM_EPS_VOL  = 1e-9; /*!< Constant value used to compute geomtric epsilon for volume */

static const double GEOM_EPS_SURF = 1e-9; /*!< Constant value used to compute geomtric epsilon for surface */

static const double GEOM_EPS_DIST = 1e-9; /*!< Minimum distance between two vertices */

enum {
  NOT_DEFINE,
  IN_STACK,
  UNCHANGED_CYCLE,
  CHANGED_CYCLE,
};

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static PDM_polygon_status_t
_intersect_ray_triangulated_face
(
       PDM_triangulate_state_t *tri_state,
       int                     *tri_vtx,
       double                  *vtx_coord,
 const int                      face_vtx_n,
 const int                     *face_vtx,
 const double                  *ray_origin,
 const double                  *ray_direction,
       double                  *intersection_coord
 )
{
  /* Triangulate current face */
  int n_tri;
  double tri_coord[9];

  if (face_vtx_n == 3) {
    // Triangle
    n_tri = 1;
    for (int ivtx = 0; ivtx < 3; ivtx++) {
      tri_vtx[ivtx] = face_vtx[ivtx];
    }
  }
  else if (face_vtx_n == 4) {
    // Quadrangle
    n_tri = PDM_triangulate_quadrangle(3,
                                       vtx_coord,
                                       NULL,
                                       face_vtx,
                                       tri_vtx);
  }
  else {
    // Polygon
    n_tri = PDM_triangulate_polygon(3,
                                    face_vtx_n,
                                    vtx_coord,
                                    NULL,
                                    face_vtx,
                                    PDM_TRIANGULATE_MESH_DEF,
                                    tri_vtx,
                                    tri_state);
  }

  /* Intersect sub-triangles with edge ray */
  /* If multiple intersections, pick the one closest from the ray's origin */
  double t = HUGE_VAL;
  PDM_polygon_status_t stat = PDM_POLYGON_OUTSIDE;
  for (int itri = 0; itri < n_tri; itri++) {

    int *_tri_vtx = tri_vtx + 3*itri;

    for (int i = 0; i < 3; i++) {
      int _vtx_id = _tri_vtx[i] - 1;
      memcpy(tri_coord + 3*i,
             vtx_coord + 3*_vtx_id,
             sizeof(double) * 3);
    }

    double _t;
    double ip[3];
    PDM_triangle_status_t _stat = PDM_triangle_ray_intersection(ray_origin,
                                                                ray_direction,
                                                                tri_coord,
                                                                ip,
                                                                &_t,
                                                                NULL);

    if (_stat == PDM_TRIANGLE_INSIDE && _t >= 0.) {
      stat = PDM_POLYGON_INSIDE;
      if (_t < t) {
        t = _t;
      }
    }

  } // End of loop on subtriangles

  if (stat == PDM_POLYGON_INSIDE) {
    intersection_coord[0] = ray_origin[0] + t*ray_direction[0];
    intersection_coord[1] = ray_origin[1] + t*ray_direction[1];
    intersection_coord[2] = ray_origin[2] + t*ray_direction[2];
  }

  return stat;
}


static PDM_polygon_status_t
_intersect_ray_face
(
       double                  *face_center,
       double                  *face_normal,
       double                  *face_coord,
       double                  *vtx_coord,
 const int                      face_vtx_n,
 const int                     *face_vtx,
 const double                  *ray_origin,
 const double                  *ray_direction,
       double                  *intersection_coord
 )
{
  double face_bound[6] = {
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL,
    HUGE_VAL, -HUGE_VAL
  };
  for (int i = 0; i < face_vtx_n; i++) {
    int vtx_id = face_vtx[i] - 1;
    double *vc = vtx_coord + 3*vtx_id;
    for (int j = 0; j < 3; j++) {
      face_coord[3*i+j] = vc[j];
      face_bound[2*j  ] = PDM_MIN(face_bound[2*j  ], vc[j]);
      face_bound[2*j+1] = PDM_MAX(face_bound[2*j+1], vc[j]);
    }
  }

  // Inflate the face's bounding box
  double d = 0.;
  for (int j = 0; j < 3; j++) {
    d += (face_bound[2*j+1] - face_bound[2*j])*(face_bound[2*j+1] - face_bound[2*j]);
  }
  d = 0.1*sqrt(d);
  for (int j = 0; j < 3; j++) {
    face_bound[2*j  ] -= d;
    face_bound[2*j+1] += d;
  }

  double t;
  return PDM_polygon_ray_intersection(ray_origin,
                                      ray_direction,
                                      face_vtx_n,
                                      face_coord,
                                      face_center,
                                      face_normal,
                                      face_bound,
                                      intersection_coord,
                                      &t,
                                      NULL);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *  \brief Compute a dynamic geometric epsilon from a characteristic length
 *
 *    @param [in]  characteristic_length  Characteristic length
 *    @param [in]  consEpsilon           Constant part
 *    @return                            Geometric epsilon
 */

double
PDM_geom_elem_geometric_epsilon
(
 const double characteristic_length,
 const double const_epsilon
 )
{
  return PDM_MAX(const_epsilon * characteristic_length, GEOM_EPS_MIN);
}

/**
 *  \brief Triangle surface vector
 *
 *  @param [in]  n_triangle      Number of triangles
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] surface_vector  Surface Vector
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tria_surface_vector
(
 const int     n_triangle,
 const int    *connectivity,
 const double *coords,
 double       *surface_vector,
 double       *characteristic_length,
 int          *is_degenerated
 )
{
  for (int itri = 0; itri < n_triangle; itri++) {

    const int i = connectivity[3*itri    ] - 1;
    const int j = connectivity[3*itri + 1] - 1;
    const int k = connectivity[3*itri + 2] - 1;

    double v1[3];
    double v2[3];
    double *surface_vectorTri = surface_vector + 3*itri;

    v1[0] = coords[3*j    ] - coords[3*i    ];
    v1[1] = coords[3*j + 1] - coords[3*i + 1];
    v1[2] = coords[3*j + 2] - coords[3*i + 2];

    v2[0] = coords[3*k    ] - coords[3*i    ];
    v2[1] = coords[3*k + 1] - coords[3*i + 1];
    v2[2] = coords[3*k + 2] - coords[3*i + 2];

    if (characteristic_length != NULL) {

      double v3[3];

      v3[0] = coords[3*k    ] - coords[3*j    ];
      v3[1] = coords[3*k + 1] - coords[3*j + 1];
      v3[2] = coords[3*k + 2] - coords[3*j + 2];

      double normV1 = PDM_MODULE(v1);
      double normV2 = PDM_MODULE(v2);
      double normV3 = PDM_MODULE(v3);

      characteristic_length[itri] = PDM_MIN(PDM_MIN(normV1, normV2), normV3);

    }

    PDM_CROSS_PRODUCT(surface_vectorTri, v1, v2);

    surface_vectorTri[0] *= 0.5;
    surface_vectorTri[1] *= 0.5;
    surface_vectorTri[2] *= 0.5;

    if ((characteristic_length != NULL) && (is_degenerated != NULL)) {

      double normsurface_vectorTri = PDM_MODULE(surface_vectorTri);
      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristic_length[itri], GEOM_EPS_SURF);
      is_degenerated[itri] = 0;
      if (normsurface_vectorTri <= eps_loc)
        is_degenerated[itri] = 1;

    }
  }
}

/**
 *  \brief Triangle area
 *
 *  @param [in]  n_triangle      Number of triangles
 *  @param [in]  surface_vector         surface_vector vectors
 *  @param [out] area           Area
 */

void
PDM_geom_elem_tria_area
(
 const int     n_triangle,
 const double *surface_vector,
 double *area
 )

{
  for (int itri = 0; itri < n_triangle; itri++) {
    const double *surface_vectorTri = surface_vector + 3*itri;
    area[itri] = PDM_MODULE(surface_vectorTri);
  }
}

/**
 *  \brief Triangle center
 *
 *  @param [in]  n_triangle      Number of triangles
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] center         center
 */

void
PDM_geom_elem_tria_center
(
 const int     n_triangle,
 const int    *connectivity,
 const double *coords,
       double *center
)
{
  for (int itri = 0; itri < n_triangle; itri++) {

    const int i = connectivity[3*itri    ] - 1;
    const int j = connectivity[3*itri + 1] - 1;
    const int k = connectivity[3*itri + 2] - 1;

    double *center_tri = center + 3*itri;

    center_tri[0] = (coords[3*i    ] + coords[3*j    ] + coords[3*k    ])/3;
    center_tri[1] = (coords[3*i + 1] + coords[3*j + 1] + coords[3*k + 1])/3;
    center_tri[2] = (coords[3*i + 2] + coords[3*j + 2] + coords[3*k + 2])/3;
  }
}


/**
 *  \brief Tetrahedra oriented volume
 *
 *  @param [in]  n_tetrahedra           Number of tetrahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  coords                Vertice coordinates
 *  @param [out] volume                Volume
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tetra_oriented_volume
(
 const int     n_tetrahedra,
 const int    *connectivity,
 const double *coords,
 double       *volume,
 double       *characteristic_length,
 int          *is_degenerated
 )
{
  for (int itet = 0; itet < n_tetrahedra; itet++) {

    const int *tetra_vtx = connectivity + 4*itet;

    const int i1 = tetra_vtx[0] - 1;
    const int i2 = tetra_vtx[1] - 1;
    const int i3 = tetra_vtx[2] - 1;
    const int i4 = tetra_vtx[3] - 1;

    const int n_triangle = 1;
    double surface_vector[3];
    double face_center[3];
    double fC_i4[3];

    PDM_geom_elem_tria_surface_vector(n_triangle,
                                      tetra_vtx,
                                      coords,
                                      surface_vector,
                                      NULL,
                                      NULL);

    if (characteristic_length != NULL) {

      double vectV1V2[3] =
        {coords[3*i2    ] - coords[3*i1    ],
         coords[3*i2 + 1] - coords[3*i1 + 1],
         coords[3*i2 + 2] - coords[3*i1 + 2]};

      double vectV1V3[3] =
        {coords[3*i3    ] - coords[3*i1    ],
         coords[3*i3 + 1] - coords[3*i1 + 1],
         coords[3*i3 + 2] - coords[3*i1 + 2]};

      double vectV2V3[3] =
        {coords[3*i3    ] - coords[3*i2    ],
         coords[3*i3 + 1] - coords[3*i2 + 1],
         coords[3*i3 + 2] - coords[3*i2 + 2]};

      double vectV1V4[3] =
        {coords[3*i4    ] - coords[3*i1    ],
         coords[3*i4 + 1] - coords[3*i1 + 1],
         coords[3*i4 + 2] - coords[3*i1 + 2]};

      double normV1V2 = PDM_MODULE(vectV1V2);
      double normV1V3 = PDM_MODULE(vectV1V3);
      double normV2V3 = PDM_MODULE(vectV2V3);
      double normV1V4 = PDM_MODULE(vectV1V4);

      characteristic_length[itet] = PDM_MIN(normV1V2,                    normV1V3);
      characteristic_length[itet] = PDM_MIN(characteristic_length[itet],  normV2V3);
      characteristic_length[itet] = PDM_MIN(characteristic_length[itet],  normV1V2);
      characteristic_length[itet] = PDM_MIN(characteristic_length[itet],  normV1V3);
      characteristic_length[itet] = PDM_MIN(characteristic_length[itet],  normV1V4);

    }

    PDM_geom_elem_tria_center (n_triangle,
                               connectivity,
                               coords,
                               face_center);

    fC_i4[0] = coords[3*i4    ] - face_center[0];
    fC_i4[1] = coords[3*i4 + 1] - face_center[1];
    fC_i4[2] = coords[3*i4 + 2] - face_center[2];

    volume[itet] = 1./3. * PDM_DOT_PRODUCT(fC_i4, surface_vector);

    if ((characteristic_length != NULL) && (is_degenerated != NULL)) {

      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristic_length[itet], GEOM_EPS_VOL);
      is_degenerated[itet] = 0;
      if (volume[itet]  <= eps_loc)
        is_degenerated[itet] = 1;

    }
  }
}

/**
 *  \brief Tetrahedra center
 *
 *  @param [in]  n_tetrahedra    Number of tetrahedra
 *  @param [in]  connectivity   Connectivity
 *  @param [in]  coords         Vertice coordinates
 *  @param [out] center         center
 */

void
PDM_geom_elem_tetra_center
(
 const int     n_tetrahedra,
 const int    *connectivity,
 const double *coords,
 double *center
 )

{
  for (int itet = 0; itet < n_tetrahedra; itet++) {

    const int i1 = connectivity[4*itet    ] - 1;
    const int i2 = connectivity[4*itet + 1] - 1;
    const int i3 = connectivity[4*itet + 2] - 1;
    const int i4 = connectivity[4*itet + 3] - 1;

    double *centerTet = center + 3*itet;

    centerTet[0] = (coords[3*i1    ] + coords[3*i2    ]
                  + coords[3*i3    ] + coords[3*i4    ])/4;
    centerTet[1] = (coords[3*i1 + 1] + coords[3*i2 + 1]
                  + coords[3*i3 + 1] + coords[3*i4 + 1])/4;
    centerTet[2] = (coords[3*i1 + 2] + coords[3*i2 + 2]
                  + coords[3*i3 + 2] + coords[3*i4 + 2])/4;
  }
}


/**
 *  \brief Tetrahedra Faces
 *
 *  @param [in]  n_tetrahedra       Number of tetrahedra
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] face_connectivity  Face connectivity
 */

void
PDM_geom_elem_tetra_faces
(
 const int     n_tetrahedra,
 const int     orientation,
 const int    *connectivity,
 int          *face_connectivity_idx,
 int          *face_connectivity
 )
{
  const int n_face       = 4;         /* 4 triangles */
  const int n_vertex     = n_face * 3; /* 4 * 3 vertices */
  const int n_vertex_elmt  = 4;         /* 4 vertices */

  face_connectivity_idx[0] = 0;

  for (int ielt = 0; ielt < n_tetrahedra; ielt++) {

    for (int iface = 0; iface < n_face; iface++)
      face_connectivity_idx[ielt * n_face + iface + 1] =
        face_connectivity_idx[ielt * n_face + iface] + 3;

    if (orientation == 0) {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 2];

      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt + 1];

      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 2];

    }

    else {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 1];

    }
  }
}


/**
 *  \brief HexahedraFaces
 *
 *  @param [in]  n_hexahedra        Number of hexahedra
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] face_connectivity  Face connectivity
 */

void
PDM_geom_elem_hexa_faces
(
 const int     n_hexahedra,
 const int     orientation,
 const int    *connectivity,
 int          *face_connectivity_idx,
 int          *face_connectivity
 )
{
  const int n_face        = 6;          /* 6 quadrangles  */
  const int n_vertex      = n_face * 4; /* 6 * 4 vertices */
  const int n_vertex_elmt = 8;          /* 8 vertices     */

  face_connectivity_idx[0] = 0;

  for (int ielt = 0; ielt < n_hexahedra; ielt++) {

    for (int iface = 0; iface < n_face; iface++) {
      face_connectivity_idx[ielt * n_face + iface + 1] = face_connectivity_idx[ielt * n_face + iface] + 4;
    }

    if (orientation == 0) {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 7];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 6];

      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 7];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt + 6];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt + 7];

      face_connectivity[n_vertex * ielt + 16] = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 17] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 18] = connectivity[n_vertex_elmt * ielt + 6];
      face_connectivity[n_vertex * ielt + 19] = connectivity[n_vertex_elmt * ielt + 2];

      face_connectivity[n_vertex * ielt + 20] = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 21] = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 22] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 23] = connectivity[n_vertex_elmt * ielt + 1];

    }

    else {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 6];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt + 7];
      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 5];

      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 7];
      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt + 7];
      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt + 6];
      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 16] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 17] = connectivity[n_vertex_elmt * ielt + 6];
      face_connectivity[n_vertex * ielt + 18] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 19] = connectivity[n_vertex_elmt * ielt + 1];

      face_connectivity[n_vertex * ielt + 20] = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 21] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 22] = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 23] = connectivity[n_vertex_elmt * ielt + 0];

    }
  }
}


/**
 *  \brief Prism Faces
 *
 *  @param [in]  n_prism            Number of Prism
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] face_connectivity  Face connectivity
 */

void
PDM_geom_elem_prism_faces
(
 const int     n_prism,
 const int     orientation,
 const int    *connectivity,
 int          *face_connectivity_idx,
 int          *face_connectivity
 )
{
  const int n_face        = 5;         /* 3 quadrangles + 2 triangles */
  const int n_vertex      = 3*4 + 2*3; /* */
  const int n_vertex_elmt = 6;         /* 6 vertices */

  face_connectivity_idx[0] = 0;

  for (int ielt = 0; ielt < n_prism; ielt++) {

    face_connectivity_idx[ielt * n_face + 1] = face_connectivity_idx[ielt * n_face    ] + 3;
    face_connectivity_idx[ielt * n_face + 2] = face_connectivity_idx[ielt * n_face + 1] + 3;
    face_connectivity_idx[ielt * n_face + 3] = face_connectivity_idx[ielt * n_face + 2] + 4;
    face_connectivity_idx[ielt * n_face + 4] = face_connectivity_idx[ielt * n_face + 3] + 4;
    face_connectivity_idx[ielt * n_face + 5] = face_connectivity_idx[ielt * n_face + 4] + 4;

    if (orientation == 0) {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 2];

      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 5];

      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 16] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 17] = connectivity[n_vertex_elmt * ielt + 3];

    }

    else {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 2];

      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt + 1];

      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt + 5];
      face_connectivity[n_vertex * ielt + 16] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 17] = connectivity[n_vertex_elmt * ielt    ];

    }
  }
}

/**
 *  \brief Pyramid Faces
 *
 *  @param [in]  n_pyramid          Number of pyramid
 *  @param [in]  orientation       Surface vector oriented towards inside cell (0) or outside (1)
 *  @param [in]  connectivity      Connectivity
 *  @param [out] face_connectivity  Face connectivity
 */

void
PDM_geom_elem_pyramid_faces
(
 const int   n_pyramid,
 const int   orientation,
 const int  *connectivity,
 int        *face_connectivity_idx,
 int        *face_connectivity
 )
{
  const int n_face        = 5;         /* 1 quadrangle + 4 triangles */
  const int n_vertex      = 1*4 + 4*3; /* */
  const int n_vertex_elmt = 5;         /* 5 vertices */

  face_connectivity_idx[0] = 0;

  for (int ielt = 0; ielt < n_pyramid; ielt++) {

    face_connectivity_idx[ielt * n_face + 1] = face_connectivity_idx[ielt * n_face    ] + 4;
    face_connectivity_idx[ielt * n_face + 2] = face_connectivity_idx[ielt * n_face + 1] + 3;
    face_connectivity_idx[ielt * n_face + 3] = face_connectivity_idx[ielt * n_face + 2] + 3;
    face_connectivity_idx[ielt * n_face + 4] = face_connectivity_idx[ielt * n_face + 3] + 3;
    face_connectivity_idx[ielt * n_face + 5] = face_connectivity_idx[ielt * n_face + 4] + 3;

    if (orientation == 0) {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt + 4];

      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt + 4];

    }

    else {

      face_connectivity[n_vertex * ielt + 0]  = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 1]  = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 2]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 3]  = connectivity[n_vertex_elmt * ielt    ];

      face_connectivity[n_vertex * ielt + 4]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 5]  = connectivity[n_vertex_elmt * ielt    ];
      face_connectivity[n_vertex * ielt + 6]  = connectivity[n_vertex_elmt * ielt + 1];

      face_connectivity[n_vertex * ielt + 7]  = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 8]  = connectivity[n_vertex_elmt * ielt + 1];
      face_connectivity[n_vertex * ielt + 9]  = connectivity[n_vertex_elmt * ielt + 2];

      face_connectivity[n_vertex * ielt + 10] = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 11] = connectivity[n_vertex_elmt * ielt + 2];
      face_connectivity[n_vertex * ielt + 12] = connectivity[n_vertex_elmt * ielt + 3];

      face_connectivity[n_vertex * ielt + 13] = connectivity[n_vertex_elmt * ielt + 4];
      face_connectivity[n_vertex * ielt + 14] = connectivity[n_vertex_elmt * ielt + 3];
      face_connectivity[n_vertex * ielt + 15] = connectivity[n_vertex_elmt * ielt    ];

    }
  }
}


/**
 *  \brief Edges properties
 *
 *  @param [in]  n_edges                Number of edges
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] length                Length
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_edges_properties
(
 const int     n_edges,
 const int    *connectivity,
 const double *coords,
 double       *length,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
 )
{

  for (int iedge = 0; iedge < n_edges; iedge++) {

    const int *edge_vtx = connectivity + 2*iedge;

    const int i1 = edge_vtx[0] - 1;
    const int i2 = edge_vtx[1] - 1;

    double *centerEdge = center + 3*iedge;

    for (int i = 0; i < 3; i++)
      centerEdge[i] = 0.5 * (coords[3*i1 + i] + coords[3*i2 + i]);

    length[iedge] = sqrt((coords[3*i2 + 0] - coords[3*i1 + 0])
                         * (coords[3*i2 + 0] - coords[3*i1 + 0])
                         + (coords[3*i2 + 1] - coords[3*i1 + 1])
                         * (coords[3*i2 + 1] - coords[3*i1 + 1])
                         + (coords[3*i2 + 2] - coords[3*i1 + 2])
                         * (coords[3*i2 + 2] - coords[3*i1 + 2]));

    if (characteristic_length != NULL)
      characteristic_length[iedge] = length[iedge];

    if (is_degenerated != NULL) {
      is_degenerated[iedge] = 0;
      if (length[iedge] < GEOM_EPS_DIST)
        is_degenerated[iedge] = 1;
    }
  }
}


/**
 *  \brief Triangle properties
 *
 *  @param [in]  n_triangle             Number of triangles
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surface_vector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tria_properties
(
 const int     n_triangle,
 const int    *connectivity,
 const double *coords,
 double       *surface_vector,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
 )
{
  PDM_geom_elem_tria_surface_vector (n_triangle,
                                     connectivity,
                                     coords,
                                     surface_vector,
                                     characteristic_length,
                                     is_degenerated);

  PDM_geom_elem_tria_center (n_triangle,
                             connectivity,
                             coords,
                             center);

}


/**
 *  \brief Tetrahedra properties
 *
 *  @param [in]  n_tetrahedra           Number of tetrahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tetra_properties
(
 const int     n_tetrahedra,
 const int    *connectivity,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
 )

{
  PDM_geom_elem_tetra_oriented_volume (n_tetrahedra,
                                       connectivity,
                                       coords,
                                       volume,
                                       characteristic_length,
                                       is_degenerated);

  PDM_geom_elem_tetra_center (n_tetrahedra,
                              connectivity,
                              coords,
                              center);
}


/**
 * \brief Quadrangle properties
 *
 *  @param [in]  n_triangle             Number of quadrangles
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surface_vector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 *
 *  @return                     The status of properties computation convergence
 */

int
PDM_geom_elem_quad_properties
(
 const int     n_quadrangle,
 const int    *connectivity,
 const double *coords,
       double *surface_vector,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
)
{
  int *_connectivity_index = NULL;
  PDM_malloc(_connectivity_index, n_quadrangle + 1, int);
  int convergence;

  _connectivity_index[0] = 0;
  for (int i = 1; i < n_quadrangle + 1; i++) {
    _connectivity_index[i] = _connectivity_index[i-1] + 4;
  }

  convergence = PDM_geom_elem_polygon_properties(n_quadrangle,
                                                 _connectivity_index,
                                                 connectivity,
                                                 coords,
                                                 surface_vector,
                                                 center,
                                                 characteristic_length,
                                                 is_degenerated);

  PDM_free(_connectivity_index);

  return convergence;
}

/**
 * \brief Compute the barycentric coordinates of a set of points inside
          their belonging polygons.
 *
 *  @param [in]  n_points               Number of points
 *  @param [in]  pts_locations          Numbering of the belonging polygons inside the _connectivity_index
 *  @param [in]  _connectivity_index     Mesh connectivity Index
 *  @param [in]  connectivity          Mesh connectivity
 *  @param [in]  coords                Mesh coordinates
 *  @param [out] bar_coords_idx        Pointer to the barycentric coordinates index
 *  @param [out] bar_coords_idx        Pointer to the barycentric coordinates
 *
 *  @return                     The status of properties computation convergence
 */

int
PDM_geom_elem_compute_polygon_barycentric_coordinates
(
  const int           n_points,
  const int          *pts_locations,
  const double       *pts_coords,
  const int          *_connectivity_index,
  const int          *connectivity,
  const double       *coords,
        int         **bar_coords_idx,
        double      **bar_coords
)
{

  int convergence = 1;
  double local_pts[3];

  /* Tableaux locaux */

  const double eps_base = 1e-10;
  double* coords_sommets = NULL;
  double* s = NULL;
  double* dist = NULL;
  double* aire = NULL;
  double* scal_prod = NULL;

  PDM_malloc(*bar_coords_idx, n_points+1, int);
  int* _bar_coords_index = *bar_coords_idx;

  int prev_n_sommets = 0;
  int n_sommets = 0;

  _bar_coords_index[0] = 0;
  /* Boucle sur les points distants */
  for (int ipoint =  0; ipoint < n_points; ipoint++ ) {

    /* Initialisation - Copie locale */

    int is_on_edge = 0;
    int isVertex = 0;
    int ielt = pts_locations[ipoint] - 1;
    prev_n_sommets = n_sommets;
    n_sommets =  _connectivity_index[ielt+1] - _connectivity_index[ielt];

    local_pts[0] = pts_coords[3*ipoint    ];
    local_pts[1] = pts_coords[3*ipoint + 1];
    local_pts[2] = pts_coords[3*ipoint + 2];

    if (ipoint == 0) {
      PDM_malloc(coords_sommets, 3 * n_sommets, double);
      PDM_malloc(s             , 3 * n_sommets, double);
      PDM_malloc(dist          ,     n_sommets, double);
      PDM_malloc(aire          ,     n_sommets, double);
      PDM_malloc(scal_prod     ,     n_sommets, double);
    }
    else {
      if (prev_n_sommets < n_sommets) {
        PDM_realloc(coords_sommets, coords_sommets , 3 * n_sommets, double);
        PDM_realloc(s             , s              , 3 * n_sommets, double);
        PDM_realloc(dist          , dist           ,     n_sommets, double);
        PDM_realloc(aire          , aire           ,     n_sommets, double);
        PDM_realloc(scal_prod     , scal_prod      ,     n_sommets, double);
      }
    }

    for (int isom = 0; isom < n_sommets; isom++) {
      coords_sommets[3*isom]   =
        coords[3*(connectivity[_connectivity_index[ielt]+isom]-1)];

      coords_sommets[3*isom+1] =
        coords[3*(connectivity[_connectivity_index[ielt]+isom]-1)+1];

      coords_sommets[3*isom+2] =
        coords[3*(connectivity[_connectivity_index[ielt]+isom]-1)+2];
    }

    /* Projection sur un plan moyen */

    double bary[3];

    PDM_polygon_compute_barycenter (n_sommets, &(coords_sommets[0]), bary);

    double n[3]   = {0, 0, 1};
    double p0 [3] = {0 ,0, 0};
    double p10[3] = {0 ,0, 0};
    double l10[3] = {0 ,0, 0};
    double p20[3] = {0 ,0, 0};
    double l20[3] = {0 ,0, 0};

    /*Compute polygon normal*/
    PDM_polygon_parameterize (n_sommets, &(coords_sommets[0]),p0,p10,l10,p20,l20, n);

    PDM_plane_projection2 (local_pts, bary, n, local_pts);

    for (int isom = 0; isom < n_sommets; isom++) {

      double *pt1 = &(coords_sommets[0]) + 3 *isom;
      PDM_plane_projection2 (pt1, bary, n, pt1);

    }

    double bounds[6] = {DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX,
                        DBL_MAX, -DBL_MAX};

    for (int isom = 0; isom < n_sommets; isom++) {
      bounds[0] = PDM_MIN(bounds[0], coords_sommets[3*isom]);
      bounds[1] = PDM_MAX(bounds[1], coords_sommets[3*isom]);

      bounds[2] = PDM_MIN(bounds[2], coords_sommets[3*isom + 1]);
      bounds[3] = PDM_MAX(bounds[3], coords_sommets[3*isom + 1]);

      bounds[4] = PDM_MIN(bounds[4], coords_sommets[3*isom + 2]);
      bounds[5] = PDM_MAX(bounds[5], coords_sommets[3*isom + 2]);
    }


    /* Verification que le point est dans l'element */
    double closest[3];
    double dist_min = DBL_MAX;

    PDM_polygon_status_t position_inout
      = PDM_polygon_evaluate_position(local_pts, n_sommets, &(coords_sommets[0]), closest, &dist_min);

    if (position_inout == PDM_POLYGON_OUTSIDE) {
      local_pts[0] = closest[0];
      local_pts[1] = closest[1];
      local_pts[2] = closest[2];
    }

    /* Calcul des coordonnnees barycentriques */

    double min_dist = DBL_MAX;
    for (int isom = 0; isom < n_sommets; isom++) {

      int inext = (isom + 1) % n_sommets;
      double *vect = &s[0] + 3*isom;
      double l_edge;
      vect[0] = coords_sommets[3*inext]   - coords_sommets[3*isom];
      vect[1] = coords_sommets[3*inext+1] - coords_sommets[3*isom+1];
      vect[2] = coords_sommets[3*inext+2] - coords_sommets[3*isom+2];
      l_edge  = PDM_MODULE (vect);
      min_dist = PDM_MIN(l_edge, min_dist);
    }
    double eps = PDM_MAX(min_dist * eps_base, 1.e-30);

    for (int isom = 0; isom < n_sommets; isom++) {

      double *vect = &s[0] + 3*isom;
      vect[0] = coords_sommets[3*isom]   - local_pts[0];
      vect[1] = coords_sommets[3*isom+1] - local_pts[1];
      vect[2] = coords_sommets[3*isom+2] - local_pts[2];
      dist[isom] = PDM_MODULE (vect);
    }

    int current_vertex;
    for (int isom = 0; isom < n_sommets; isom++) {
      int inext = (isom + 1) % n_sommets;
      double *vect1 = &s[0] + 3 * isom;
      double *vect2 = &s[0] + 3 * inext;
      double pvect[3];

      scal_prod[isom] = PDM_DOT_PRODUCT (vect1, vect2);
      PDM_CROSS_PRODUCT(pvect, vect1, vect2);

      double sign = PDM_DOT_PRODUCT (pvect, n);
      aire[isom] = PDM_MODULE(pvect);

      if (sign < 0) {
        aire[isom] = -aire[isom];
      }

      if (dist[isom] <= eps) {

        isVertex = 1;
        current_vertex = isom;
        break;
      }

      else if ((fabs(aire[isom]) <= eps)  && (scal_prod[isom] < 0)) {

        is_on_edge = 1;
        current_vertex = isom;
        break;

      }

    }

    _bar_coords_index[ipoint+1] = _bar_coords_index[ipoint] + n_sommets;
    //Vector/Pointer containing Barycentric coordinates

    PDM_realloc(*bar_coords ,*bar_coords ,(_bar_coords_index[ipoint+1]) ,double);

    double *_bar_coords = *bar_coords;

    double* _local_bary_coords  = &(_bar_coords[ _bar_coords_index[ipoint] ]);

    /* Le point distant est un sommet */

    if (isVertex) {
      for (int isom = 0; isom < n_sommets; isom++)
        _local_bary_coords[isom] = 0.;
      _local_bary_coords[current_vertex] = 1.;
    }
    else if (is_on_edge) {
      /* Le point distant est sur arete */
      for (int isom = 0; isom < n_sommets; isom++)
        _local_bary_coords[isom] = 0.;

      int nextPoint = (current_vertex + 1) % n_sommets;

      _local_bary_coords[current_vertex] =
        dist[nextPoint]     / (dist[nextPoint]+dist[current_vertex]);
      _local_bary_coords[nextPoint]     =
        dist[current_vertex] / (dist[nextPoint]+dist[current_vertex]);
    }
    else {
      /* Cas general */
      double sigma = 0;
      for (int isom = 0; isom < n_sommets; isom++) {
        double coef = 0.;
        int previousVertex = (isom - 1 + n_sommets) % n_sommets;
        int nextVertex = (isom + 1) % n_sommets;

        if (fabs(aire[previousVertex]) > eps)
          coef += (dist[previousVertex] - scal_prod[previousVertex]/dist[isom]) / aire[previousVertex];
        if (fabs(aire[isom]) > eps)
          coef += (dist[nextVertex]     - scal_prod[isom]/dist[isom])           / aire[isom];
        sigma += coef;
        _local_bary_coords[isom] = coef;
      }

      if (PDM_ABS(sigma) >= eps ) {
        for (int isom = 0; isom < n_sommets; isom++) {
          _local_bary_coords[isom] /= sigma;
        }
      }
      else {
        double abs_sigma = fabs(sigma);
        printf("Warning : Mise à NAN %f %f\n", abs_sigma,  eps);
        for (int isom = 0; isom < n_sommets; isom++) {
          _local_bary_coords[isom] = NAN;
        }
      }

      /* Check Result */
      for (int isom = 0; isom <  n_sommets; isom++) {
        if ( isnan(_local_bary_coords[isom])     ||
                   _local_bary_coords[isom] < 0. ||
                   _local_bary_coords[isom] > 1. ) {

          convergence = 0;
  /*        double dist_min = DBL_MAX;
          int k_min = 0;
          double t_min;

          for (int k = 0; k < n_sommets; k++) {
            _local_bary_coords[k] = 0.0;
          }

          for (int k = 0; k < n_sommets; k++) {
            double *p1 = &(coords_sommets[3 * k]);
            double *p2 = &(coords_sommets[3 * ((k+1) % n_sommets)]);
            double closest[3];
            double t;

            double dist2 = fvmc_distance_to_line (local_pts,
                                                 p1,
                                                 p2,
                                                 &t,
                                                 closest);
            if (dist2 < dist_min) {
              t_min = t;
              k_min = k;
            }
          }

          _local_bary_coords[k_min] = 1 - t_min;
          _local_bary_coords[(k_min + 1) % n_sommets] = t_min;
*/
          break;

        }

      }

    }

    if (0 == 1) {
      if ((n_points == 1) && (pts_locations[0] == 1)) {

        PDM_printf("coord %i %i :", ipoint+1, ielt+1);
        PDM_printf(" %12.5e %12.5e %12.5e", pts_coords[3*ipoint],
                    pts_coords[3*ipoint+1],
                    pts_coords[3*ipoint+2] );
        PDM_printf("\n");

        PDM_printf("coo b %i :", ipoint+1);
        for (int isom = 0; isom < n_sommets; isom++) {
          PDM_printf(" %f", _local_bary_coords[isom]);
        }
        PDM_printf("\n");
      }
    }
  }

  PDM_free(coords_sommets);
  PDM_free(s);
  PDM_free(aire);
  PDM_free(dist);
  PDM_free(scal_prod);

  return convergence;
}


/**
 *  \brief Polygon properties
 *
 *  @param [in]  nPolygon              Number of polygon
 *  @param [in]  _connectivity_index     Connectivity Index
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surface_vector         Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 *
 *  @return                        The status of properties computation convergence
 */

int
PDM_geom_elem_polygon_properties
(
 const int     nPolygon,
 const int    *_connectivity_index,
 const int    *connectivity,
 const double *coords,
 double       *surface_vector,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
)
{

  int convergence = 1;

  const double dispMin = 1e-9; /* Minimum displacement */
  const double big = 1e30;     /* Big value */

  const int nIterMax = 100;    /* Maximum iteration number */

  for (int ifac = 0; ifac < nPolygon; ifac++) {

    /*
     * Define local pointer
     * --------------------
     */

    const int n_verticesFace = _connectivity_index[ifac + 1]
      - _connectivity_index[ifac    ];

    const int *connectivityFace = connectivity + _connectivity_index[ifac];

    double *surface_vectorFace = surface_vector + 3*ifac; /* Face surface vector */

    double *centerFace = center + 3*ifac; /* Face Center */

    if (characteristic_length != NULL) {
      characteristic_length[ifac] = big;
    }

    /*
     * Initialization
     * --------------
     */

    int nIter = 0;       /* Number of iteration */

    for (int i = 0; i < 3; i++) {
      surface_vectorFace[i] = 0;
      centerFace[i]        = 0;
    }

    /* Initialize face center to the barycenter */

    for (int ivert = 0; ivert < n_verticesFace; ivert++) {
      const int vert = connectivityFace[ivert] - 1;
      for (int i = 0; i < 3; i++)
        centerFace[i] += coords[3*vert + i];
    }

    for (int i = 0; i < 3; i++)
      centerFace[i] /= n_verticesFace;

    /*
     * Compute cell center and surface vector
     * --------------------------------------
     */

    while (1) {

      double displacement[3] = {0, 0, 0}; /* Cell center displacement */

      nIter += 1;

      for (int i = 0; i < 3; i++)
        surface_vectorFace[i] = 0;

      double areaFace = 0; /* Face area */

      for (int ivert = 0; ivert < n_verticesFace; ivert++) {

        const int vert1 = connectivityFace[ivert] - 1;
        const int vert2 = connectivityFace[(ivert + 1) % n_verticesFace] - 1;

        /* Edge center */

        double edgeCenter[3] = {0.5 * (coords[3*vert1    ] + coords[3*vert2    ]),
                                0.5 * (coords[3*vert1 + 1] + coords[3*vert2 + 1]),
                                0.5 * (coords[3*vert1 + 2] + coords[3*vert2 + 2])};

        /* Edge vector */

        double vectV1V2[3] = {coords[3*vert2    ] - coords[3*vert1    ],
                              coords[3*vert2 + 1] - coords[3*vert1 + 1],
                              coords[3*vert2 + 2] - coords[3*vert1 + 2]};

        if (characteristic_length != NULL) {
          double norm_V1V2 = PDM_MODULE(vectV1V2);
          characteristic_length[ifac] = PDM_MIN(characteristic_length[ifac], norm_V1V2);
        }

        /* Vector face center -> edge center */

        double vectFECenter[3] = {edgeCenter[0] - centerFace[0],
                                  edgeCenter[1] - centerFace[1],
                                  edgeCenter[2] - centerFace[2]};

        /* Compute area of the triangle (face center, vert1, vert2) */

        double surface_vectorTria[3];
        PDM_CROSS_PRODUCT(surface_vectorTria, vectFECenter, vectV1V2);

        for (int i = 0; i < 3; i++)
          surface_vectorTria[i] *= 0.5;

        const double areaTri = PDM_MODULE(surface_vectorTria);

        areaFace += areaTri;
        for (int i = 0; i < 3; i++) {
          surface_vectorFace[i] += surface_vectorTria[i];
          displacement[i] += areaTri * vectFECenter[i];
        }

        areaFace += areaTri;

      }

      double denomAreaFace = 1. / PDM_MAX(fabs(areaFace), GEOM_EPS_MIN);

      for (int i = 0; i < 3; i++) {
        displacement[i] = 2./3. * denomAreaFace * displacement[i];
        centerFace[i] += displacement[i];
      }

      /*
       * Check convergence
       */

      const double normDisp = PDM_MODULE(displacement);

      if (normDisp < dispMin) {
        break;
      }

      /*
       * Check Number of iteration
       */

      else if (nIterMax < nIter) {
        convergence = false;
        break;
      }
    } /* while (1) */

    if ((characteristic_length != NULL) && (is_degenerated != NULL)) {

      double normsurface_vector = PDM_MODULE(surface_vectorFace);
      double eps_loc = PDM_geom_elem_geometric_epsilon(characteristic_length[ifac], GEOM_EPS_SURF);
      is_degenerated[ifac] = 0;
      if (normsurface_vector <= eps_loc)
        is_degenerated[ifac] = 1;
    }
  } /* for (int ifac = 0; ifac < nPolygon; ifac++) */

  if (0 == 1) {
    PDM_printf( "surface_vector : ");
    for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
      PDM_printf( "%12.5e ",surface_vector[ifac]);
    }
    PDM_printf( "\n");

    PDM_printf( "center : ");
    for (int ifac = 0; ifac < 3*nPolygon; ifac++) {
      PDM_printf( "%12.5e ",center[ifac]);
    }
    PDM_printf( "\n");
  }
  return convergence;
}


/**
 *  \brief Hexahedra properties
 *
 *  @param [in]  n_hexahedra            Number of hexahedra
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_hexa_properties
(
 const int     n_hexahedra,
 const int    *connectivity,
 const int     n_vertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
)
{

  const int orientation = 1; /*  Surface vector oriented towards inside cell outside */

  const int n_quadrangle     = 6;
  const int n_triangle       = 0;
  const int n_hexahedraFaces = n_quadrangle + n_triangle;
  const int n_faces          = n_hexahedraFaces * n_hexahedra;

  int *face_connectivity          = NULL;
  int *face_connectivity_idx      = NULL;
  int *cell_face_connectivity_idx = NULL;
  int *cell_face_connectivity     = NULL;
  PDM_malloc(face_connectivity         , (n_quadrangle*4 + n_triangle*3) * n_hexahedra, int);
  PDM_malloc(face_connectivity_idx     , n_faces + 1                                  , int);
  PDM_malloc(cell_face_connectivity_idx, n_hexahedra + 1                              , int);
  PDM_malloc(cell_face_connectivity    , n_faces                                      , int);

  /*
   * Get hexahedra faces
   */

  PDM_geom_elem_hexa_faces (n_hexahedra,
                            orientation,
                            connectivity,
                            face_connectivity_idx,
                            face_connectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < n_faces; i++) {
    cell_face_connectivity[i] = i + 1;
  }

  cell_face_connectivity_idx[0] = 0;
  for (int i = 1; i < n_hexahedra + 1; i++)
    cell_face_connectivity_idx[i] = cell_face_connectivity_idx[i-1] + n_hexahedraFaces;

  /*
   * Compute Volume and center
   */

  PDM_geom_elem_polyhedra_properties (1,
                                      n_hexahedra,
                                      n_faces,
                                      face_connectivity_idx,
                                      face_connectivity,
                                      cell_face_connectivity_idx,
                                      cell_face_connectivity,
                                      n_vertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristic_length,
                                      is_degenerated);

  /*
   * Free
   */
  PDM_free(face_connectivity         );
  PDM_free(face_connectivity_idx     );
  PDM_free(cell_face_connectivity    );
  PDM_free(cell_face_connectivity_idx);

}


/**
 *  \brief Prism properties
 *
 *  @param [in]  n_prism                Number of prism
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_prism_properties
(
 const int     n_prism,
 const int    *connectivity,
 const int     n_vertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
)
{

  const int orientation = 1; //  Surface vector oriented towards inside cell outside

  const int n_quadrangle  = 3;
  const int n_triangle    = 2;
  const int n_prism_faces = n_quadrangle + n_triangle;
  const int n_faces       = n_prism_faces * n_prism;

  int *face_connectivity          = NULL;
  int *face_connectivity_idx      = NULL;
  int *cell_face_connectivity_idx = NULL;
  int *cell_face_connectivity     = NULL;
  PDM_malloc(face_connectivity         , (n_quadrangle*4 + n_triangle*3) * n_prism, int);
  PDM_malloc(face_connectivity_idx     , n_faces + 1                              , int);
  PDM_malloc(cell_face_connectivity_idx, n_prism + 1                              , int);
  PDM_malloc(cell_face_connectivity    , n_faces                                  , int);

  /*
   * Get prism faces
   */
  PDM_geom_elem_prism_faces (n_prism,
                             orientation,
                             connectivity,
                             face_connectivity_idx,
                             face_connectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < n_faces; i++) {
    cell_face_connectivity[i] = i + 1;
  }

  cell_face_connectivity_idx[0] = 0;
  for (int i = 1; i < n_prism + 1; i++) {
    cell_face_connectivity_idx[i] = cell_face_connectivity_idx[i-1] + n_prism_faces;
  }

  /*
   * Compute Volume and center
   */
  PDM_geom_elem_polyhedra_properties (1,
                                      n_prism,
                                      n_faces,
                                      face_connectivity_idx,
                                      face_connectivity,
                                      cell_face_connectivity_idx,
                                      cell_face_connectivity,
                                      n_vertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristic_length,
                                      is_degenerated);


  /*
   * Free
   */
  PDM_free(face_connectivity         );
  PDM_free(face_connectivity_idx     );
  PDM_free(cell_face_connectivity    );
  PDM_free(cell_face_connectivity_idx);

}


/**
 *  \brief Pyramid properties
 *
 *  @param [in]  n_pyramid              Number of pyramid
 *  @param [in]  connectivity           Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                 Vertices coordinates
 *  @param [out] volume                 Volume
 *  @param [out] center                 Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_pyramid_properties
(
 const int      n_pyramid,
 const int     *connectivity,
 const int      n_vertices,
 const double  *coords,
       double  *volume,
       double  *center,
       double  *characteristic_length,
       int     *is_degenerated
)
{
  const int orientation = 1; /*  Surface vector oriented towards inside cell outside */

  const int n_quadrangle    = 1;
  const int n_triangle      = 4;
  const int n_pyramid_faces = n_quadrangle + n_triangle;
  const int n_faces         = n_pyramid_faces * n_pyramid;

  int *face_connectivity          = NULL;
  int *face_connectivity_idx      = NULL;
  int *cell_face_connectivity_idx = NULL;
  int *cell_face_connectivity     = NULL;
  PDM_malloc(face_connectivity         , (n_quadrangle*4 + n_triangle*3) * n_pyramid, int);
  PDM_malloc(face_connectivity_idx     , n_faces + 1                                , int);
  PDM_malloc(cell_face_connectivity_idx, n_pyramid + 1                              , int);
  PDM_malloc(cell_face_connectivity    , n_faces                                    , int);

  /*
   * Get pyramid faces
   */
  PDM_geom_elem_pyramid_faces (n_pyramid,
                               orientation,
                               connectivity,
                               face_connectivity_idx,
                               face_connectivity);

  /*
   * Define cell to face connectivity
   */

  for (int i = 0; i < n_faces; i++)
    cell_face_connectivity[i] = i + 1;

  cell_face_connectivity_idx[0] = 0;
  for (int i = 1; i < n_pyramid + 1; i++)
    cell_face_connectivity_idx[i] = cell_face_connectivity_idx[i-1] + n_pyramid_faces;

  /*
   * Compute Volume and center
   */
  PDM_geom_elem_polyhedra_properties (1,
                                      n_pyramid,
                                      n_faces,
                                      face_connectivity_idx,
                                      face_connectivity,
                                      cell_face_connectivity_idx,
                                      cell_face_connectivity,
                                      n_vertices,
                                      coords,
                                      volume,
                                      center,
                                      characteristic_length,
                                      is_degenerated);

  /*
   * Free
   */
  PDM_free(face_connectivity);
  PDM_free(face_connectivity_idx);
  PDM_free(cell_face_connectivity);
  PDM_free(cell_face_connectivity_idx);

}


/**
 *  \brief Polyhedra properties
 *
 *  Compute cellCenter and volume. Reorient cell_face_connectivity
 *
 *  @param [in]  is_oriented                 1 if cell_face_connectivity is already oriented, 0 otherwise
 *  @param [in]  n_polyhedra                 Number of polyhedra
 *  @param [in]  n_face                      Number of faces
 *  @param [in]  face_connectivity_idx        Face connectivity index
 *  @param [in]  face_connectivity           Face connectivity
 *  @param [in,out]  cell_face_connectivity_idx  Cell to face connectivity index
 *  @param [in,out]  cell_face_connectivity     Cell to face connectivity
 *  @param [in]  n_vertices                  Number of vertices
 *  @param [in]  coords                     Vertices coordinates
 *  @param [out] volume                     Volume
 *  @param [out] center                     Center
 *  @param [out] characteristic_length       Characteristic length (active if != NULL)
 *  @param [out] is_degenerated              Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_polyhedra_properties
(
 const int     is_oriented,
 const int     n_polyhedra,
 const int     n_face,
 const int    *face_connectivity_idx,
 const int    *face_connectivity,
 const int    *cell_face_connectivity_idx,
       int    *cell_face_connectivity,
 const int     n_vertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristic_length,
 int          *is_degenerated
)
{
  const double big = 1e30;
  int convergence = 1;
  int *colorVertice;
  PDM_malloc(colorVertice,n_vertices,int);

  int  lpolyhedra_vertices = 24;
  int *polyhedra_vertices;
  PDM_malloc(polyhedra_vertices,lpolyhedra_vertices,int); //First allocation

  for (int i = 0; i < n_vertices; i++){
    colorVertice[i] = false;
  }

  /*
   * Compute face properties
   */

  double *surface_vector = NULL;
  double *face_center    = NULL;
  PDM_malloc(surface_vector, 3 *n_face, double);
  PDM_malloc(face_center   , 3 *n_face, double);

  int convergenceFace = PDM_geom_elem_polygon_properties (n_face,
                                                          face_connectivity_idx,
                                                          face_connectivity,
                                                          coords,
                                                          surface_vector,
                                                          face_center,
                                                          NULL,
                                                          NULL);

  if (0 == 1) {

    PDM_printf( "face_connectivity : \n");
    for (int ipoly = 0; ipoly < n_face; ipoly++) {
      PDM_printf( "  - face %i : ", ipoly+1);
      for (int j = face_connectivity_idx[ipoly]; j < face_connectivity_idx[ipoly+1]; j++) {
        PDM_printf( "%i ",face_connectivity[j]);
      }
      PDM_printf( "\n");
    }

    PDM_printf( "surface_vector : ");
    for (int ipoly = 0; ipoly < 3 * n_face; ipoly++) {
      PDM_printf( "%12.5e ",surface_vector[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "face_center : ");
    for (int ipoly = 0; ipoly < 3 * n_face; ipoly++) {
      PDM_printf( "%12.5e ",face_center[ipoly]);
    }
    PDM_printf( "\n");
  }

  /*
   * Loop on polyhedra
   */

  int keyMax = 2 * n_vertices;

  PDM_hash_tab_t *hash_orient = NULL;
  int n_key_poly = 0;
  int s_key_poly = 10;

  int *key_poly = NULL;

  double volume_t = 0;
  int max_n_poly_face = 0;

  int *stack = NULL;
  int *tag_face = NULL;

  if (!is_oriented) {

    hash_orient = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, &keyMax);

    PDM_malloc(key_poly, s_key_poly, int);

    for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {
      const int poly_idx   = cell_face_connectivity_idx[ipoly];
      const int n_poly_face = cell_face_connectivity_idx[ipoly + 1] - poly_idx;
      max_n_poly_face = PDM_MAX (max_n_poly_face, n_poly_face);
    }

    PDM_malloc(stack   , max_n_poly_face, int);
    PDM_malloc(tag_face, max_n_poly_face, int);

  }

  for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {

    double *poly_center = center + 3*ipoly;
    double  disp[3];

    const int poly_idx   = cell_face_connectivity_idx[ipoly];
    const int n_poly_face = cell_face_connectivity_idx[ipoly + 1] - poly_idx;

    if (characteristic_length != NULL)
      characteristic_length[ipoly] = big;

    /*
     * Intialize cell center to the barycenter
     */

    n_key_poly = 0;

    /* Search polyhedra vertices */

    int n_polyhedra_vertices = 0;

    for (int iface = 0; iface < n_poly_face; iface++) {

      if (!is_oriented) {
        tag_face[iface] = NOT_DEFINE;
      }

      const int face          = abs(cell_face_connectivity[poly_idx + iface]) - 1;
      if (!is_oriented) {
        cell_face_connectivity[poly_idx + iface] = face+1;
      }
      const int face_idx       = face_connectivity_idx[face];
      const int n_face_vertices = face_connectivity_idx[face+1] - face_idx;

      for (int ivert = 0; ivert < n_face_vertices; ivert++) {
        const int vertex = face_connectivity[face_idx + ivert] - 1;

        if (!is_oriented) {
          const int inext = (ivert + 1) % n_face_vertices;
          const int vertex_next = face_connectivity[face_idx + inext] - 1;
          const int key = vertex + vertex_next;

          if (n_key_poly >= s_key_poly) {
            s_key_poly *= 2;
            PDM_realloc(key_poly, key_poly, s_key_poly, int);
          }

          key_poly[n_key_poly++] = key;

          int *edge;
          PDM_malloc(edge,3,int);
          edge[0] = vertex;
          edge[1] = vertex_next;
          edge[2] = iface;

          PDM_hash_tab_data_add (hash_orient, (void *) &key, edge);
        }

        if (!colorVertice[vertex]) {
          colorVertice[vertex] = 1;

          if (n_polyhedra_vertices >= lpolyhedra_vertices) {
            lpolyhedra_vertices *= 2;
            PDM_realloc(polyhedra_vertices ,polyhedra_vertices , lpolyhedra_vertices,int);
          }
          polyhedra_vertices[n_polyhedra_vertices++] = vertex;
        }
      }
    }

    if (!is_oriented) {
      int n_stack = -1;
      stack[++n_stack] = 0;
      tag_face[0] = IN_STACK;

      while (n_stack >= 0) {

        int iFace = stack[n_stack--];

        const int face          = abs(cell_face_connectivity[poly_idx + iFace]) - 1;
        const int face_idx       = face_connectivity_idx[face];
        const int n_face_vertices = face_connectivity_idx[face+1] - face_idx;

        for (int ivert = 0; ivert < n_face_vertices; ivert++) {
          const int inext = (ivert + 1) % n_face_vertices;

          const int vertex = face_connectivity[face_idx + ivert] - 1;
          const int vertex_next = face_connectivity[face_idx + inext] - 1;
          int key = vertex + vertex_next;

          int n_data =          PDM_hash_tab_n_data_get(hash_orient, &key);
          int **data = (int **) PDM_hash_tab_data_get  (hash_orient, &key);

          int j_current_edge = -1;
          for (int j = 0; j < n_data; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int is_same_edge    = (vertex == _edge[0]) && (vertex_next == _edge[1]);
              int is_same_face    = iFace == _edge[2];
              if (is_same_edge && is_same_face) {
                j_current_edge = j;
                break;
              }
            }
          }

          assert (j_current_edge > -1);

          for (int j = 0; j < n_data; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int is_inverse_edge = (vertex == _edge[1]) && (vertex_next == _edge[0]);
              int is_same_edge    = (vertex == _edge[0]) && (vertex_next == _edge[1]);
              int is_same_face    = iFace == _edge[2];

              int neighbour = _edge[2];

              if (!is_same_face) {
                if (is_same_edge || is_inverse_edge) {

                  if (tag_face[iFace] < UNCHANGED_CYCLE) {

                    if (tag_face[neighbour] >= UNCHANGED_CYCLE) {
                      if (tag_face[neighbour] == UNCHANGED_CYCLE) {
                        if (is_same_edge) {
                          tag_face[iFace] = CHANGED_CYCLE;
                        }
                        else  {
                          tag_face[iFace] = UNCHANGED_CYCLE;
                        }
                      }
                      else {
                        if (is_same_edge) {
                          tag_face[iFace] = UNCHANGED_CYCLE;
                        }
                        else  {
                          tag_face[iFace] = CHANGED_CYCLE;
                        }
                      }

                      PDM_free(data[j]);
                      PDM_free(data[j_current_edge]);
                      data[j] = NULL;
                      data[j_current_edge] = NULL;

                    }
                  }

                  if (tag_face[neighbour] == NOT_DEFINE) {
                    stack[++n_stack] = neighbour;
                    tag_face[neighbour] = IN_STACK;
                  }

                  break;

                }
              }
            }
          }
        }

        if (tag_face[iFace] == IN_STACK) {
          tag_face[iFace] = UNCHANGED_CYCLE;
        }
      }


      for (int iface = 0; iface < n_poly_face; iface++) {

        if (tag_face[iface] == CHANGED_CYCLE) {
          cell_face_connectivity[poly_idx + iface] = -cell_face_connectivity[poly_idx + iface];
        }
      }

    }

    /* Compute barycenter */

    for (int i = 0; i < 3; i++)
      poly_center[i] = 0.;

    for (int j = 0; j < n_polyhedra_vertices; j++) {
      const int vertex = polyhedra_vertices[j];
      colorVertice[vertex] = false;
      for (int i = 0; i < 3; i++)
        poly_center[i] += coords[3*vertex + i];
    }

    for (int i = 0; i < 3; i++)
      poly_center[i] /= n_polyhedra_vertices;

    n_polyhedra_vertices = 0;

    /*
     * Intialize volume
     */

    volume[ipoly] = 0.;

    /*
     * Intialize cell center displacement
     */

    for (int i = 0; i < 3; i++)
      disp[i] = 0.;

    /*
     * Loop on faces
     */

    for (int iface = 0; iface < n_poly_face; iface++) {

      const int face          = abs(cell_face_connectivity[poly_idx + iface]) - 1;
      const int direction     = (cell_face_connectivity[poly_idx + iface] < 0) ? -1 : 1;

      const int face_idx        = face_connectivity_idx[face];
      const int n_face_vertices = face_connectivity_idx[face+1] - face_idx;

      /*
       * Loop on vertices
       */

      for (int ivert = 0; ivert < n_face_vertices; ivert++) {

        int vert1;
        int vert2;
        int ivert1 = ivert;

        if (direction > 0) {
          vert1 = face_connectivity[face_idx + ivert1] - 1;
          vert2 = face_connectivity[face_idx + (ivert1 + 1) % n_face_vertices] - 1;
        }
        else {
          ivert1 = n_face_vertices - 1 - ivert;
          vert1 = face_connectivity[face_idx + (ivert1 + 1) % n_face_vertices] - 1;
          vert2 = face_connectivity[face_idx + ivert1] - 1;
        }

        if (characteristic_length != NULL)  {

          /* Vector vert1 -> vert2 */

          const double vectV1V2[3] =
            {coords[3*vert2    ] - coords[3*vert1    ],
             coords[3*vert2 + 1] - coords[3*vert1 + 1],
             coords[3*vert2 + 2] - coords[3*vert1 + 2]};

          double normV1V2 = PDM_MODULE (vectV1V2);

          characteristic_length[ipoly] = PDM_MIN(characteristic_length[ipoly], normV1V2);

        }

        /* Vector face center -> vert1 */

        double vectFCV1[3] =
          {coords[3*vert1    ] - face_center[3 * face    ],
           coords[3*vert1 + 1] - face_center[3 * face + 1],
           coords[3*vert1 + 2] - face_center[3 * face + 2]};

        /* Vector face center -> vert2 */

        double vectFCV2[3] =
          {coords[3*vert2    ] - face_center[3 * face    ],
           coords[3*vert2 + 1] - face_center[3 * face + 1],
           coords[3*vert2 + 2] - face_center[3 * face + 2]};

        /* Vector cell center -> face center */

        double vectCCFC[3] =
          {face_center[3 * face    ] - poly_center[0],
           face_center[3 * face + 1] - poly_center[1],
           face_center[3 * face + 2] - poly_center[2]};

        double surface_vectorTri[3];
        PDM_CROSS_PRODUCT (surface_vectorTri, vectFCV1, vectFCV2);

        for (int i = 0; i < 3; i++)
          surface_vectorTri[i] *= 0.5;

        /* Oriented volume */

        double volumeTet = 1./3 * PDM_DOT_PRODUCT (surface_vectorTri, vectCCFC);

        volume[ipoly] += volumeTet;

        for (int i = 0; i < 3; i++)
          disp[i] = disp[i] + volumeTet * vectCCFC[i];

      }

    }

    int signeVol = (volume[ipoly] < 0.) ? -1 : 1 ;

    if (signeVol == -1) {

      if (!is_oriented) {
        volume[ipoly] = -volume[ipoly];
        for (int iface = 0; iface < n_poly_face; iface++) {
          cell_face_connectivity[poly_idx + iface] = - cell_face_connectivity[poly_idx + iface];
        }
      }

      else {
        PDM_printf( "Warning polyhedraProperties : volume < 0 for polyhedron '%i' (%f)\n",
                    ipoly + 1, volume[ipoly]);
      }
    }

    double denomVol = 1 / PDM_MAX(fabs(volume[ipoly]), GEOM_EPS_MIN);

    for (int i = 0; i < 3; i++)
      poly_center[i] =  poly_center[i] + signeVol * denomVol * disp[i];


    volume_t += volume[ipoly];

    /*
     * Check convergence
     */

    if (!convergenceFace)
      convergence = false;

    if ((characteristic_length != NULL) && (is_degenerated != NULL)) {

      double eps_loc = PDM_geom_elem_geometric_epsilon (characteristic_length[ipoly], GEOM_EPS_VOL);
      is_degenerated[ipoly] = 0;
      if (fabs(volume[ipoly]) <= eps_loc)
        is_degenerated[ipoly] = 1;
    }

    if (!is_oriented) {
      PDM_hash_tab_purge (hash_orient, PDM_TRUE);
    }

  }
  PDM_UNUSED(volume_t);

  PDM_free(polyhedra_vertices);

  if (!is_oriented) {
    PDM_free(key_poly);
    PDM_free(tag_face);
    PDM_free(stack);
    PDM_hash_tab_free (hash_orient);
  }

  if (0 == 1) {

    PDM_printf( "surface_vector : ");
    for (int ipoly = 0; ipoly < 3 * n_face; ipoly++) {
      PDM_printf( "%12.5e ",surface_vector[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "face_center : ");
    for (int ipoly = 0; ipoly < 3 * n_face; ipoly++) {
      PDM_printf( "%12.5e ",face_center[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "isDegenrated : ");
    for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {
      PDM_printf( "%i ",is_degenerated[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "characteristic_length  : ");
    for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {
      PDM_printf( "%12.5e ",characteristic_length[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "volume  : ");
    for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {
      PDM_printf( "%12.5e ",volume[ipoly]);
    }
    PDM_printf( "\n");

    PDM_printf( "center  : ");
    for (int ipoly = 0; ipoly < 3 * n_polyhedra; ipoly++) {
      PDM_printf( "%12.5e ",center[ipoly]);
    }
    PDM_printf( "\n");

  }

  PDM_free(surface_vector);
  PDM_free(face_center);
  PDM_free(colorVertice);

  if (!convergence)
    PDM_printf( "Warning polyhedraProperties : some polyhedra faces are not planar\n");

}



void
PDM_geom_elem_polyhedra_properties_triangulated
(
 const int     is_oriented,
 const int     n_polyhedra,
 const int     n_face,
 const int    *face_connectivity_idx,
 const int    *face_connectivity,
 const int    *cell_face_connectivity_idx,
       int    *cell_face_connectivity,
 const int     n_vertices,
 const double *coords,
       double *volume,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
)
{
  PDM_UNUSED(n_vertices);
  PDM_UNUSED(characteristic_length);

  /**
   * TO DO :
   *  - compute characteristic_length
   *  - handle (!is_oriented) case
   */

  /*
   * Triangulate faces
   */
  int max_face_vtx_n = 0;
  int *face_tria_idx;
  PDM_malloc(face_tria_idx, n_face + 1, int);
  face_tria_idx[0] = 0;
  for (int iface = 0; iface < n_face; iface++) {
    int face_vtx_n = face_connectivity_idx[iface+1] - face_connectivity_idx[iface];
    max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    int face_tria_n = face_vtx_n - 2;
    face_tria_idx[iface+1] = face_tria_idx[iface] + face_tria_n;
  }

  PDM_triangulate_state_t *tri_state = PDM_triangulate_state_create(max_face_vtx_n);

  int *tria_vtx;
  PDM_malloc(tria_vtx, face_tria_idx[n_face] * 3, int);

  for (int iface = 0; iface < n_face; iface++) {

    const int *_face_vtx = face_connectivity + face_connectivity_idx[iface];
    int face_vtx_n = face_connectivity_idx[iface+1] - face_connectivity_idx[iface];

    int *_tria_vtx = tria_vtx + 3*face_tria_idx[iface];

    int n_tria;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tria = 1;
      memcpy(_tria_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tria = PDM_triangulate_quadrangle(3,
                                          coords,
                                          NULL,
                                          _face_vtx,
                                          _tria_vtx);
    }
    else {
      /* Polygonal face */
      n_tria = PDM_triangulate_polygon(3,
                                       face_vtx_n,
                                       coords,
                                       NULL,
                                       _face_vtx,
                                       PDM_TRIANGULATE_MESH_DEF,
                                       _tria_vtx,
                                       tri_state);
    }

    assert(n_tria == face_tria_idx[iface+1] - face_tria_idx[iface]);
  }
  PDM_triangulate_state_destroy(tri_state);

  assert(is_oriented == 1);


  for (int ipoly = 0; ipoly < n_polyhedra; ipoly++) {

    double *poly_center = center + 3*ipoly;
    const int poly_idx   = cell_face_connectivity_idx[ipoly];
    // const int n_poly_face = cell_face_connectivity_idx[ipoly + 1] - poly_idx;

    int ref_face = PDM_ABS(cell_face_connectivity[poly_idx]) - 1;
    int ref_point = face_connectivity[face_connectivity_idx[ref_face]] - 1;

    int connec[4];
    connec[0] = ref_point+1;
    for (int i = 0; i < 3; i++) {
      poly_center[i] = 0;
    }
    volume[ipoly] = 0;


    for (int iface = cell_face_connectivity_idx[ipoly]; iface < cell_face_connectivity_idx[ipoly+1]; iface++) {

      int face_id = PDM_ABS (cell_face_connectivity[iface]) - 1;
      int sign    = PDM_SIGN(cell_face_connectivity[iface]);

      for (int itria = face_tria_idx[face_id]; itria < face_tria_idx[face_id+1]; itria++) {

        for (int i = 0; i < 3; i++) {
          connec[i+1] = tria_vtx[3*itria + i];
        }

        double tetra_center[3];
        double tetra_volume;
        PDM_geom_elem_tetra_properties(1,
                                       connec,
                                       coords,
                                       &tetra_volume,
                                       tetra_center,
                                       NULL,
                                       NULL);

        tetra_volume *= sign;

        volume[ipoly] += tetra_volume;
        for (int i = 0; i < 3; i++) {
          poly_center[i] += tetra_volume * tetra_center[i];
        }

      }

    }

    if (PDM_ABS(volume[ipoly]) < 1e-15) {
      if (is_degenerated != NULL) {
        is_degenerated[ipoly] = 1;
      }
    }
    else {
      double ivol = 1./volume[ipoly];
      for (int i = 0; i < 3; i++) {
        poly_center[i] *= ivol;
      }
    }

  }
  //PDM_free(surface_vector);
  PDM_free(face_tria_idx);
  PDM_free(tria_vtx);
}






/**
 *  \brief Compute downwind and updind elemt of all edges (or -1 if not found )
 *
 *  If the face centers and normals are not provided, the faces are triangulated
 *
 *  @param [in]  n_face               Number of faces
 *  @param [in]  n_edge               Number of edges
 *  @param [in]  cell_face_idx        Index for cell-face connectivity
 *  @param [in]  cell_face            Cell-face connectivity
 *  @param [in]  face_vtx_idx         Index for face-vertex connectivity
 *  @param [in]  face_vtx             Face-vertex connectivity
 *  @param [in]  vtx_cell_idx         Index for vertex-cell connectivity
 *  @param [in]  vtx_cell             Vertex-cell connectivity
 *  @param [in]  edge_vtx             Edge-vertex connectivity
 *  @param [in]  vtx_coord            Vertex coordinates (size = 3*n_vtx)
 *  @param [in]  face_center          Face center (or NULL)
 *  @param [in]  face_normal          Face normal vectors (or NULL, need not be normalized)
 *  @param [out] upwind_cell_out      Cell number corresponding of upwind cell (or -1)   (size =   n_edge)
 *  @param [out] downwind_cell_out    Cell number corresponding of downwind cell (or -1) (size =   n_edge)
 *  @param [out] upwind_face_out      Face number corresponding of upwind face (or -1)   (size =   n_edge)
 *  @param [out] downwind_face_out    Face number corresponding of downwind face (or -1) (size =   n_edge)
 *  @param [out] upwind_point_out     Coordinates of upwind point                        (size = 3*n_edge)
 *  @param [out] downwind_point_out   Coordinates of downwind point                      (size = 3*n_edge)
 *
 */
void
PDM_geom_elem_edge_upwind_and_downwind
(
 int         n_face,
 int         n_edge,
 PDM_g_num_t *cell_ln_to_gn,
 int         *cell_face_idx,
 int         *cell_face,
 int         *face_vtx_idx,
 int         *face_vtx,
 int         *vtx_cell_idx,
 int         *vtx_cell,
 int         *edge_vtx,
 double      *vtx_coord,
 double      *face_center,
 double      *face_normal,
 int        **upwind_cell_out,
 int        **downwind_cell_out,
 int        **upwind_face_out,
 int        **downwind_face_out,
 double     **upwind_point_out,
 double     **downwind_point_out
)
{
  /* Allocate stuff */
  PDM_malloc(*upwind_cell_out   ,     n_edge, int   );
  PDM_malloc(*downwind_cell_out ,     n_edge, int   );
  PDM_malloc(*upwind_face_out   ,     n_edge, int   );
  PDM_malloc(*downwind_face_out ,     n_edge, int   );
  PDM_malloc(*upwind_point_out  , 3 * n_edge, double);
  PDM_malloc(*downwind_point_out, 3 * n_edge, double);

  int    *upwind_cell    = *upwind_cell_out;
  int    *downwind_cell  = *downwind_cell_out;
  int    *upwind_face    = *upwind_face_out;
  int    *downwind_face  = *downwind_face_out;
  double *upwind_point   = *upwind_point_out;
  double *downwind_point = *downwind_point_out;

  for(int i = 0; i < n_edge; ++i) {
    upwind_cell  [i] = -1;
    downwind_cell[i] = -1;
    upwind_face  [i] = -1;
    downwind_face[i] = -1;
    upwind_point  [3*i  ] = DBL_MAX;
    upwind_point  [3*i+1] = DBL_MAX;
    upwind_point  [3*i+2] = DBL_MAX;
    downwind_point[3*i  ] = DBL_MAX;
    downwind_point[3*i+1] = DBL_MAX;
    downwind_point[3*i+2] = DBL_MAX;
  }


  PDM_triangulate_state_t *tri_state  = NULL;
  int                     *tri_vtx    = NULL;
  double                  *poly_coord = NULL;

  int triangulate_faces = (face_center == NULL ||
                           face_normal == NULL);

  int max_face_vtx_n = 0;
  for (int iface = 0; iface < n_face; iface++) {
    max_face_vtx_n = PDM_MAX(max_face_vtx_n,
                             face_vtx_idx[iface+1] - face_vtx_idx[iface]);
  }

  if (triangulate_faces) {
    if (max_face_vtx_n > 4) {
      tri_state = PDM_triangulate_state_create(max_face_vtx_n);
    }

    PDM_malloc(tri_vtx, (max_face_vtx_n - 2)*3, int);
  }
  else {
    PDM_malloc(poly_coord, max_face_vtx_n * 3, double);
  }

  int *is_visited_face = PDM_array_zeros_int(n_face);
  int *visited_face;
  PDM_malloc(visited_face, n_face, int);
  int n_visited_face = 0;

  for (int iedge = 0; iedge < n_edge; iedge++) {

    /* Initialize */
    upwind_face  [iedge] = -1;
    downwind_face[iedge] = -1;

    /* Reset visited faces */
    for (int i = 0; i< n_visited_face; i++) {
      is_visited_face[visited_face[i]] = 0;
    }
    n_visited_face = 0;

    /* Set up edge-ray */
    int ivtx1 = edge_vtx[2*iedge  ] - 1;
    int ivtx2 = edge_vtx[2*iedge+1] - 1;

    double ray_direction[3] = {
      vtx_coord[3*ivtx1  ] - vtx_coord[3*ivtx2  ],
      vtx_coord[3*ivtx1+1] - vtx_coord[3*ivtx2+1],
      vtx_coord[3*ivtx1+2] - vtx_coord[3*ivtx2+2],
    };

    double ray_origin[3] = {
      0.5*(vtx_coord[3*ivtx1  ] + vtx_coord[3*ivtx2  ]),
      0.5*(vtx_coord[3*ivtx1+1] + vtx_coord[3*ivtx2+1]),
      0.5*(vtx_coord[3*ivtx1+2] + vtx_coord[3*ivtx2+2]),
    };

    int found[2] = {0, 0};

    for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {

      // /* Reset visited faces */
      // for (int i = 0; i< n_visited_face; i++) {
      //   is_visited_face[visited_face[i]] = 0;
      // }
      // n_visited_face = 0;

      int vtx_id  = edge_vtx[2*iedge + idx_vtx] - 1;

      int n_cell_per_vtx = (vtx_cell_idx[vtx_id+1]-vtx_cell_idx[vtx_id]);
      int *order;
      PDM_malloc(order, n_cell_per_vtx, int);

      if(cell_ln_to_gn != NULL) {
        PDM_g_num_t *vtx_cell_gnum;
        PDM_malloc(vtx_cell_gnum, n_cell_per_vtx, PDM_g_num_t);
        int idx_write = 0;
        for (int idx_cell = vtx_cell_idx[vtx_id]; idx_cell < vtx_cell_idx[vtx_id+1]; idx_cell++) {
          int cell_id = PDM_ABS(vtx_cell[idx_cell]) - 1;
          vtx_cell_gnum[idx_write++] = cell_ln_to_gn[cell_id];
        }
        PDM_order_gnum_s(vtx_cell_gnum, 1, order, n_cell_per_vtx);
        PDM_free(vtx_cell_gnum);
      } else {
        for(int i = 0; i < n_cell_per_vtx; ++i) {
          order[i] = i;
        }
      }

      for (int idx_cell0 = vtx_cell_idx[vtx_id]; idx_cell0 < vtx_cell_idx[vtx_id+1]; idx_cell0++) {
        int idx_cell = order[idx_cell0-vtx_cell_idx[vtx_id]]+vtx_cell_idx[vtx_id];

        int cell_id = PDM_ABS(vtx_cell[idx_cell]) - 1;

        /* Check if one face has the 2 vertex, it is an adjacent cell to the edge */
        int has_face_edge = 0;
        for (int idx_face = cell_face_idx[cell_id]; idx_face < cell_face_idx[cell_id+1]; idx_face++) {

          int face_id = PDM_ABS(cell_face[idx_face]) - 1;
          for (int i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {
            if (face_vtx[i] - 1 == ivtx1) {
              for (int j = face_vtx_idx[face_id]; j < face_vtx_idx[face_id+1]; j++) {
                if ((j != i) && (face_vtx[j] - 1 == ivtx2)) {
                  has_face_edge = 1;
                  break;
                }
              }
            }
          }
        }

        for (int idx_face = cell_face_idx[cell_id]; idx_face < cell_face_idx[cell_id+1]; idx_face++) {

          int face_id = PDM_ABS(cell_face[idx_face]) - 1;

          if (has_face_edge == 1) {
            is_visited_face[face_id] = 1;
            visited_face[n_visited_face++] = face_id;
          }

          if (is_visited_face[face_id]) {
            continue;
          }

          is_visited_face[face_id] = 1;
          visited_face[n_visited_face++] = face_id;

          /* Skip face if incident to current vertex */
          int has_current_vtx = 0;
          for (int i = face_vtx_idx[face_id]; i < face_vtx_idx[face_id+1]; i++) {
            if (face_vtx[i] - 1 == vtx_id) {
              has_current_vtx = 1;
              break;
            }
          }

          if (has_current_vtx) {
            continue;
          }

          int stat;
          int *_face_vtx = face_vtx + face_vtx_idx[face_id];
          int face_vtx_n = face_vtx_idx[face_id+1] - face_vtx_idx[face_id];
          double intersection_coord[3];
          if (triangulate_faces) {
            stat = _intersect_ray_triangulated_face(tri_state,
                                                    tri_vtx,
                                                    vtx_coord,
                                                    face_vtx_n,
                                                    _face_vtx,
                                                    ray_origin,
                                                    ray_direction,
                                                    intersection_coord);
          }
          else {
            stat = _intersect_ray_face(&face_center[3*face_id],
                                       &face_normal[3*face_id],
                                       poly_coord,
                                       vtx_coord,
                                       face_vtx_n,
                                       _face_vtx,
                                       ray_origin,
                                       ray_direction,
                                       intersection_coord);
          }

          if (stat == PDM_POLYGON_INSIDE) {

            found[idx_vtx] = 1;

            if (idx_vtx == 1) {
              upwind_cell[iedge] = cell_id;
              upwind_face[iedge] = face_id;
              memcpy(upwind_point + 3*iedge,
                     intersection_coord,
                     sizeof(double) * 3);
              // if (iedge == iedge_test) {
              //   log_trace("FOUND[1]: iedge:%d, face_id:%d, cell_id:%d\n", iedge, face_id, cell_id);
              // }
            }
            else {
              downwind_cell[iedge] = cell_id;
              downwind_face[iedge] = face_id;
              memcpy(downwind_point + 3*iedge,
                     intersection_coord,
                     sizeof(double) * 3);
              // if (iedge == iedge_test) {
              //   log_trace("FOUND[2]: iedge:%d, face_id:%d, cell_id:%d\n", iedge, face_id, cell_id);
              // }
            }

            break;
          }

        } // End of loop on current cell's faces

        if (found[idx_vtx]) break;

      } // End of loop on current vertex's cells

      PDM_free(order);

      /* Reverse ray direction */
      ray_direction[0] = -ray_direction[0];
      ray_direction[1] = -ray_direction[1];
      ray_direction[2] = -ray_direction[2];

    } // End of loop on current edge's vertices

  } // End of loop on edges

  PDM_free(is_visited_face);
  PDM_free(visited_face);

  if (tri_state != NULL) {
    PDM_triangulate_state_destroy(tri_state);
  }

  if (tri_vtx != NULL) {
    PDM_free(tri_vtx);
  }

  if (poly_coord != NULL) {
    PDM_free(poly_coord);
  }
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
