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

#ifdef	__cplusplus
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

int
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
 * \param [in]  ima      Identifier
 *
 */

void
PDM_Mesh_adapt_free
(
 const int  ima
);


/*----------------------------------------------------------------------------*
 *
 * Functions about source mesh definition
 *   * Volume mesh
 *   * Boundary faces ???
 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Set source mesh vertices.
 *
 * \param [in]  ima          Identifier
 * \param [in]  i_part       Current partition
 * \param [in]  n_pts        Number of points
 * \param [in]  coord        Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num   Global element number (or NULL)
 *
 */

void
PDM_Mesh_adapt_src_vtx_set
(
 const int    ima,
 const int    i_part,
 const int    n_pts,
 double       coord[],
 PDM_g_num_t  global_num[]
);


/**
 * \brief Add a connectivity block to the source mesh.
 *
 * \param [in]  ima              Identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier (ibl)
 */

int
PDM_Mesh_adapt_src_block_add
(
 const int                   ima,
 const PDM_Mesh_nodal_elt_t  block_type
);


/**
 * \brief Finalize the source mesh definition.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided (MPI collective communications)
 *
 * \param [in]  ima          Identifier
 *
 */

void
PDM_Mesh_adapt_src_finalize
(
 const int  ima
);


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
