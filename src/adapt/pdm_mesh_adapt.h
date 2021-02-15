#ifndef PDM_MESH_ADAPT_H
#define PDM_MESH_ADAPT_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \enum  PDM_Mesh_adapt_method_t
 * \brief Mesh adaptation method
 *
 */

typedef enum {

  PDM_MESH_ADAPT_REMESHING  = 0,   /*!< Remeshing */
  PDM_MESH_ADAPT_REFINMENT  = 1,   /*!< Refinment */
  PDM_MESH_ADAPT_METHOD_N   = 2    /*!< Number of methods */

} PDM_Mesh_adapt_method_t;


/**
 * \enum  PDM_Mesh_adapt_method_t
 * \brief Mesh adaptation method
 *
 */

typedef enum {

  PDM_MESH_ADAPT_ADAPTCELLS  = 0, /*!< Nuga - Adapt cells */
  PDM_MESH_ADAPT_FEFLO       = 1, /*!< Feflo.a */
  PDM_MESH_ADAPT_PARMMG      = 2, /*!< Feflo.a */
  PDM_MESH_ADAPT_TREEPART    = 3, /*!< Feflo.a */
  PDM_MESH_ADAPT_TOOL_N      = 5  /*!< Number of tools */

} PDM_Mesh_adapt_tool_t;


/**
 * \struct PDM_Mesh_adapt_t
 * \brief  Mesh adaptation workflow
 *
 */

typedef struct _PDM_Mesh_adapt_t PDM_Mesh_adapt_t;


/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*
 *
 * General functions
 *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a parallel mesh adaptation workflow
 *
 * \param [in] comm   PDM_MPI communicator
 * \param [in] method Mesh adaptation method
 *                       - PDM_MESH_ADAPT_REMESHING
 *                       - PDM_MESH_ADAPT_REFINMENT
 * \param [in] tool   Mesh adaptation tool
 *                       - PDM_MESH_ADAPT_ADAPTCELLS
 *                       - PDM_MESH_ADAPT_FEFLO
 *                       - PDM_MESH_ADAPT_PARMMG
 *                       - PDM_MESH_ADAPT_TREEPART
 *                       - ...
 * \param [in] order  Mesh order
 * \param [in] n_part Number of local mesh partition
 * \param [in] owner  Results owner
 *                       - PDM_OWNERSHIP_KEEP                 : Paradigm will
 *                       - PDM_OWNERSHIP_USER                 : Ownership is given
 *                       - PDM_OWNERSHIP_UNGET_RESULT_IS_FREE : Free all memory that
 *                                                              not be getted by user
 *
 * \return     Identifier
 *
 */

PDM_Mesh_adapt_t *
PDM_Mesh_adapt_create
(
 const PDM_MPI_Comm            comm,
 const PDM_Mesh_adapt_method_t method,
 const PDM_Mesh_adapt_tool_t   tool,
 const int                     order,
 const int                     n_part,
 const PDM_ownership_t         owner
);


/**
 *
 * \brief Free a parallel mesh adaptation workflow according to the selected
 *        ownership property
 *
 * \param [in]  ma    Mesh adaptation workfow
 *
 */

void
PDM_Mesh_adapt_free
(
 PDM_Mesh_adapt_t *ma
);


/*----------------------------------------------------------------------------*
 *
 * Functions about source mesh definition
 *   - Volume mesh    : PDM_Mesh_adapt_src_block*
 *   - Boundary mesh  : PDM_Mesh_adapt_boundary__src_block*
 *   - Face group (To describe boundary conditions) : PDM_Mesh_adapt_face_group*
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set source mesh vertices.
 *
 * \param [in]  ma           Mesh adaptation workfow
 * \param [in]  i_part       Current partition
 * \param [in]  n_pts        Number of points
 * \param [in]  coord        Coordinates (size = 3 * \p n_pts)
 * \param [in]  g_num        Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_vtx_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         n_pts,
 double            coord[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Add a connectivity block to the source mesh.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  block_type       Block type
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_src_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the source mesh.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
 *
 *   \code
 *             x 4
 *            /|\
 *           / | \
 *          /  |  \
 *       1 x- -|- -x 3
 *          \  |  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
 *
 *   \code
 *              5 x
 *               /|\
 *              //| \
 *             // |  \
 *          4 x/--|---x 3
 *           //   |  /
 *          //    | /
 *       1 x-------x 2
 *   \endcode
 *
 *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
 *
 *   \code
 *       4 x-------x 6
 *         |\     /|
 *         | \   / |
 *       1 x- \-/ -x 3
 *          \ 5x  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
 *
 *   \code
 *          8 x-------x 7
 *           /|      /|
 *          / |     / |
 *       5 x-------x6 |
 *         | 4x----|--x 3
 *         | /     | /
 *         |/      |/
 *       1 x-------x 2
 *   \endcode
 *
 * \param [in]  ma           Mesh adaptation workfow
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 * \return                   Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_src_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workfow
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 * \return                           Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_src_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the source mesh.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  g_num            Global element number (or NULL)
 *
 * \return                       Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_src_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 PDM_g_num_t       g_num[]
);

/**
 *
 * \brief Set the connectivity of a polyhedron block of the source mesh.
 *
 * Connectivity is supposed to be oriented. Connectivity can be oriented by calling
 * PDM_cellface_orient
 *
 * \param [in]  ma                Mesh adaptation workfow
 * \param [in]  i_part            Partition identifier
 * \param [in]  i_block           Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [in]  face_vtx          Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  face_g_num        Face global element number (or NULL)
 * \param [in]  cell_face_idx     Polyhedron to face index (or NULL)
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Polyhedron to face connectivity (or NULL)
 *                                The connectivity is oriented :
 *                                  - > 0 if outgoing normal,
 *                                  - < 0 otherwise
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  face_cell         Face to polyhedron connectivity (or NULL)
 *                                  - left value  : outgoing normal,
 *                                  - right value : incoming normal
 *                                (size = 2 * \p n_faces)
 * \param [in]  cell_g_num        Cell global element number (or NULL)
 *
 * \return                        Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_src_block_c_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 const int         n_faces,
 int               face_vtx_idx[],
 int               face_vtx[],
 PDM_g_num_t       face_g_num[],
 int               cell_face_idx[],
 int               cell_face[],
 int               face_cell[],
 PDM_g_num_t       cell_g_num[]
);

/**
 * \brief Add a connectivity block to the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  block_type       Block type
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_boundary_src_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the boundary source mesh.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *
 * \param [in]  ma           Mesh adaptation workfow
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 * \return                   Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_boundary_src_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the boundary source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workfow
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 * \return                           Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_boundary_src_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  g_num            Global element number (or NULL)
 *
 * \return                       Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_boundary_src_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the number of face groups
 *
 * Face groups are used to define boundary conditions.
 *
 * \param [in]  ma        Mesh adaptation workfow
 * \param [in]  n_group   Number of face groups
 *
 */

void
PDM_Mesh_adapt_face_group_n_set
(
 PDM_Mesh_adapt_t *ma,
 const int        n_group
);


/**
 * \brief Set a face group
 *
 * A face group is used to define a boundary condittion. As a face can be
 * contained in several groups, the notion of group  is more general than
 * the notion of boudary condition
 *
 * \param [in]  ma        Mesh adaptation workfow
 * \param [in]  i_group   Group identifier
 * \param [in]  n_face    Number of faces in the group
 * \param [in]  faces     List of faces
 * \param [in]  g_num     Global element number in the group (or NULL)
 *
 */

void
PDM_Mesh_adapt_face_group_set
(
 PDM_Mesh_adapt_t *ma,
 const int        i_group,
 const int        n_face,
 int              faces[],
 PDM_g_num_t      g_num[]
);

/**
 * \brief Finalize the source mesh definition.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided (MPI collective communications)
 *
 * \param [in]  ma                Mesh adaptation workfow
 *
 */

void
PDM_Mesh_adapt_src_finalize
(
 PDM_Mesh_adapt_t *ma
);

/*----------------------------------------------------------------------------*
 *
 * Functions about the geometric representation of the boundary
 * - set geometry reprsentation of the boundary (Boundary source mesh (default),
 *                                               P1 mesh, P2 mesh, P3 mesh, STL,
 *                                               CAD iges, CAD step, ...)
 * - set ridge and corners on meshes
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Add a connectivity block to the mesh that describes the geometric
 *        representation of the boundary.
 *
 *  This function is called only if \ref PDM_MESH_ADAPT_GEOM_REPR_ANOOTHER_MESH
 *  is selected.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  block_type       Block type
 *
 * \return block identifier (i_block)
 */

int
PDM_Mesh_adapt_geom_repr_block_add
(
 PDM_Mesh_adapt_t           *ma,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Set a standard block of the the mesh that describes the geometric
 *        representation of the boundary.
 *
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *
 * \param [in]  ma           Mesh adaptation workfow
 * \param [in]  i_part       Partition identifier
 * \param [in]  i_block      Block identifier
 * \param [in]  n_elts       Number of elements
 * \param [in]  elt_vtx      Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  g_num        Global element number (or NULL) (size = n_elts)
 *
 * \return                   Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_geom_repr_block_std_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               elt_vtx[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set a generic high order block of the the mesh that describes
 *        the geometric representation of the boundary.the boundary source mesh.
 *
 * \param [in]  ma                   Mesh adaptation workfow
 * \param [in]  i_part               Partition identifier
 * \param [in]  i_block              Block identifier
 * \param [in]  local_node_location  Node location in the (u, v, w) grid
 *                                   (size = 3 * n_nodes_elt)
 * \param [in]  n_elts               Number of elements
 * \param [in]  elt_node             Connectivity (size = n_nodes_elt * n_elts)
 * \param [in]  g_num                Global element number (or NULL)
 *
 * \return                           Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_geom_repr_block_ho_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 int               local_node_location[],
 const int         n_elts,
 int               elt_node[],
 PDM_g_num_t       g_num[]
);


/**
 * \brief Set the connectivity of a polygon block of the boundary source mesh.
 *
 * \param [in]  ma               Mesh adaptation workfow
 * \param [in]  i_part           Partition identifier
 * \param [in]  i_block          Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  face_vtx_idx     Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  face_vtx         face_vtx (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  g_num            Global element number (or NULL)
 *
 * \return                       Local element number in the mesh
 *
 */

int *
PDM_Mesh_adapt_geom_repr_block_f_poly_set
(
 PDM_Mesh_adapt_t *ma,
 const int         i_part,
 const int         i_block,
 const int         n_elts,
 int               face_vtx_idx[],
 int               face_vtx[],
 PDM_g_num_t       g_num[]
);



/*----------------------------------------------------------------------------*
 *
 * Functions about geomtric criteria
 *   - Volume mesh
 *   - Boundary faces ???
 *
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*
 *
 * Functions about target mesh
 *   - Volume mesh
 *   - Boundary faces ???
 *
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*
 *
 * Functions about field transfer
 *   - Defined on the Volume (cell-center/node/user points())
 *   - Boundary faces (cell-center/node/user points())
 *   - Interpolation
 *
 *----------------------------------------------------------------------------*/


#ifdef	__cplusplus
}
#endif

#endif // PDM_CLOSEST_POINTS_H
