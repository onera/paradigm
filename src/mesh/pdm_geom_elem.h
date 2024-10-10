/*
 * \file
 */

#ifndef __PDM_GEOM_ELEM_H__
#define __PDM_GEOM_ELEM_H__
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


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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
 * Public function prototypes
 *============================================================================*/

/**
 *  \brief Compute a dynamic geometric epsilon from a characteristic length
 *
 *    @param [in]  characteristic_length  Characteristic length
 *    @param [in]  const_epsilon          Constant part
 *    @return                             Geometric epsilon
 */

double
PDM_geom_elem_geometric_epsilon
(
 const double characteristic_length,
 const double const_epsilon
);

/**
 *  \brief Triangle surface vector
 *
 *  @param [in]  n_triangle             Number of triangles
 *  @param [in]  connectivity           Connectivity
 *  @param [in]  coords                 Vertice coordinates
 *  @param [out] surface_vector         Surface Vector
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
);

/**
 *  \brief Triangle area
 *
 *  @param [in]  n_triangle      Number of triangles
 *  @param [in]  surface_vector  surface_vector vectors
 *  @param [out] area            Area
 */

void
PDM_geom_elem_tria_area
(
 const int     n_triangle,
 const double *surface_vector,
       double *area
);


/**
 *  \brief Triangle center
 *
 *  @param [in]  n_triangle    Number of triangles
 *  @param [in]  connectivity  Connectivity
 *  @param [in]  coords        Vertice coordinates
 *  @param [out] center        center
 */

void
PDM_geom_elem_tria_center
(
 const int     n_triangle,
 const int    *connectivity,
 const double *coords,
       double *center
);


/**
 *  \brief Tetrahedra oriented volume
 *
 *  @param [in]  n_tetrahedra           Number of tetrahedra
 *  @param [in]  connectivity           Connectivity
 *  @param [in]  coords                 Vertice coordinates
 *  @param [out] volume                 Volume
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated         Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tetra_oriented_volume
(
 const int     n_tetrahedra,
 const int    *connectivity,
 const double *coords,
       double *volume,
       double *characteristic_length,
       int    *is_degenerated
);

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
);

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
       int    *face_connectivity_idx,
       int    *face_connectivity
);


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
       int    *face_connectivity_idx,
       int    *face_connectivity
);


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
       int    *face_connectivity_idx,
       int    *face_connectivity
);

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
       int  *face_connectivity_idx,
       int  *face_connectivity
);


/**
 *  \brief Edges properties
 *
 *  @param [in]  nEdges                Number of edges
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
 const int     nEdges,
 const int    *connectivity,
 const double *coords,
       double *length,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


/**
 *  \brief Triangle properties
 *
 *  @param [in]  n_triangle            Number of triangles
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices            Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] surface_vector        Surface vector
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated        Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_tria_properties
(
 const int     n_triangle,
 const int    *connectivity,
 const double *coords,
       double *surface_vector,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


/**
 * \brief Quadrangle properties
 *
 *  @param [in]  n_quadrangle             Number of quadrangles
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
 const int     *connectivity,
 const double  *coords,
       double  *surface_vector,
       double  *center,
       double  *characteristic_length,
       int     *is_degenerated
);


/**
 * \brief Compute the barycentric coordinates of a set of points inside
          their belonging polygons.
 *
 *  @param [in]  nPoints               Number of points
 *  @param [in]  ptsLocations          Numbering of the belonging polygons inside the connectivityIndex
 *  @param [in]  connectivityIndex     Mesh connectivity Index
 *  @param [in]  connectivity          Mesh connectivity
 *  @param [in]  coords                Mesh coordinates
 *  @param [out] barCoordsIndex        Pointer to the barycentric coordinates index
 *  @param [out] barCoordsIndex        Pointer to the barycentric coordinates
 *
 *  @return                     The status of properties computation convergence
 */

int
PDM_geom_elem_compute_polygon_barycentric_coordinates
(
 const int           n_points,
 const int          *pts_locations,
 const double       *pts_coords,
 const int          *connectivityIndex,
 const int          *connectivity,
 const double       *coords,
       int         **barCoordsIndex,
       double      **barCoords
);

/**
 *  \brief Polygon properties
 *
 *  @param [in]  n_polygon              Number of polygon
 *  @param [in]  connectivityIndex     Connectivity Index
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
 const int      n_polygon,
 const int     *connectivityIndex,
 const int     *connectivity,
 const double  *coords,
       double  *surface_vector,
       double  *center,
       double  *characteristic_length,
       int     *is_degenerated
);


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
       double *volume,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


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
       double *volume,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


/**
 *  \brief Prism properties
 *
 *  @param [in]  n_prism               Number of prism
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated        Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_prism_properties
(
 const int     n_prism,
 const int    *connectivity,
 const int     n_vertices,
 const double *coords,
       double *volume,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


/**
 *  \brief Pyramid properties
 *
 *  @param [in]  n_pyramid             Number of pyramid
 *  @param [in]  connectivity          Connectivity
 *  @param [in]  n_vertices             Number of vertices
 *  @param [in]  coords                Vertices coordinates
 *  @param [out] volume                Volume
 *  @param [out] center                Center
 *  @param [out] characteristic_length  Characteristic length (active if != NULL)
 *  @param [out] is_degenerated        Degenerated edge indicator (active if != NULL)
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
);


/**
 *  \brief Polyhedra properties
 *
 *  @param [in]  nPolyhedra                   Number of polyhedra
 *  @param [in]  n_face                       Number of faces
 *  @param [in]  face_connectivity_idx        Face connectivity index
 *  @param [in]  face_connectivity            Face connectivity
 *  @param [in]  cellToFace_connectivity_idx  Cell to face connectivity index
 *  @param [in]  cellToFace_connectivity      Cell to face connectivity
 *  @param [in]  n_vertices                    Number of vertices
 *  @param [in]  coords                       Vertices coordinates
 *  @param [out] volume                       Volume
 *  @param [out] center                       Center
 *  @param [out] characteristic_length         Characteristic length (active if != NULL)
 *  @param [out] is_degenerated                Degenerated edge indicator (active if != NULL)
 */

void
PDM_geom_elem_polyhedra_properties
(
 const  int    isOriented,
 const int     nPolyhedra,
 const int     n_face,
 const int    *face_connectivity_idx,
 const int    *face_connectivity,
 const int    *cellToFace_connectivity_idx,
       int    *cellToFace_connectivity,
 const int     n_vertices,
 const double *coords,
 double       *volume,
 double       *center,
 double       *characteristic_length,
 int         *is_degenerated
);


void
PDM_geom_elem_polyhedra_properties_triangulated
(
 const int     isOriented,
 const int     nPolyhedra,
 const int     n_face,
 const int    *face_connectivity_idx,
 const int    *face_connectivity,
 const int    *cellToFace_connectivity_idx,
       int    *cellToFace_connectivity,
 const int     n_vertices,
 const double *coords,
       double *volume,
       double *center,
       double *characteristic_length,
       int    *is_degenerated
);


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
 int          n_face,
 int          n_edge,
 PDM_g_num_t  *cell_ln_to_gn,
 int          *cell_face_idx,
 int          *cell_face,
 int          *face_vtx_idx,
 int          *face_vtx,
 int          *vtx_cell_idx,
 int          *vtx_cell,
 int          *edge_vtx,
 double       *vtx_coord,
 double       *face_center,
 double       *face_normal,
 int         **upwind_cell_out,
 int         **downwind_cell_out,
 int         **upwind_face_out,
 int         **downwind_face_out,
 double      **upwind_point_out,
 double      **downwind_point_out
);


/**
 *  \brief Compute downwind and upwind element of all edges (or -1 if not found )
 *
 *  @param [in]  i_plane              Cartesian plane (XY : 0, YZ : 1, ZX : 2)
 *  @param [in]  face_ln_to_gn        Face global IDs (optional, used to guarantee deterministic results)
 *  @param [in]  face_center          Face centers (optional, used to guarantee deterministic results)
 *  @param [in]  face_edge_idx        Index for face-edge connectivity
 *  @param [in]  face_edge            Face-edge connectivity
 *  @param [in]  n_edge               Number of edges
 *  @param [in]  edge_vtx             Edge-vertex connectivity
 *  @param [in]  vtx_face_idx         Index for vertex-face connectivity
 *  @param [in]  vtx_face             Vertex-face connectivity
 *  @param [in]  vtx_coord            Vertex coordinates (size = 3*n_vtx)
 *  @param [out] upwind_face_out      Face number corresponding of upwind face   (or -1) (size =     \p n_edge)
 *  @param [out] downwind_face_out    Face number corresponding of downwind face (or -1) (size =     \p n_edge)
 *  @param [out] upwind_edge_out      Edge number corresponding of upwind edge   (or -1) (size =     \p n_edge)
 *  @param [out] downwind_edge_out    Edge number corresponding of downwind edge (or -1) (size =     \p n_edge)
 *  @param [out] upwind_point_out     Coordinates of upwind point                        (size = 3 * \p n_edge)
 *  @param [out] downwind_point_out   Coordinates of downwind point                      (size = 3 * \p n_edge)
 *
 */
void
PDM_geom_elem_edge_upwind_and_downwind_2d
(
 int          i_plane,
 PDM_g_num_t *face_ln_to_gn,
 double      *face_center,
 int         *face_edge_idx,
 int         *face_edge,
 int          n_edge,
 int         *edge_vtx,
 int         *vtx_face_idx,
 int         *vtx_face,
 double      *vtx_coord,
 int        **upwind_face_out,
 int        **downwind_face_out,
 int        **upwind_edge_out,
 int        **downwind_edge_out,
 double     **upwind_point_out,
 double     **downwind_point_out
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ELEM_GEOM__ */
