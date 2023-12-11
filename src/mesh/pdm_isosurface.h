#ifndef PDM_ISOSURFACE_H
#define PDM_ISOSURFACE_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

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
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_dmesh.h"
#include "pdm_part_to_part.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

  struct _pdm_isosurface_t {

  };

typedef struct _pdm_isosurface_t PDM_isosurface_t;


typedef void (*PDM_isosurface_field_function_t)
(
 const double  x,
 const double  y,
 const double  z,
 double       *value
);


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// Doit-on gérer multi-domaines?

/**
 * \brief Create a \ref PDM_isosurface_t instance
 *
 * \param [in]  comm            PDM MPI communicator
 * \param [in]  mesh_dimension  Dimension of source mesh (2 or 3)
 *
 * \return Pointer to a new \ref PDM_isosurface_t instance
 *
 */

PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm           comm,
 int                    mesh_dimension
);


/* --- Input mesh definition --- */

/* Partitioned */

// Ngon

/**
 *
 * \brief Set number of partitions
 *
 * \param [in]  isos    \ref PDM_isosurface_t instance
 * \param [in]  n_part  Number of partitions
 *
 */

void
PDM_isosurface_n_part_set
(
 PDM_isosurface_t *isos,
 int               n_part
);


/**
 *
 * \brief Set connectivity
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  connect            Connectivity (1-based)
 *
 */

void
PDM_isosurface_connectivity_set
(
 PDM_isosurface_t        *isos,
 int                      i_part,
 PDM_connectivity_type_t  connectivity_type,
 int                     *connect_idx,
 int                     *connect
);


/**
 *
 * \brief Set vertex coordinates
 *
 * \param [in]  isos       \ref PDM_isosurface_t instance
 * \param [in]  i_part     Partition identifier
 * \param [in]  vtx_coord  Vertex coordinates (size = 3 * *n_vtx*)
 *
 */

void
PDM_isosurface_vtx_coord_set
(
 PDM_isosurface_t *isos,
 int               i_part,
 double           *vtx_coord
);


/**
 *
 * \brief Set global ids
 *
 * \param [in]  isos         \ref PDM_isosurface_t instance
 * \param [in]  i_part       Partition identifier
 * \param [in]  mesh_entity  Type of mesh entity
 * \param [in]  n_entity     Local number of entities
 * \param [in]  ln_to_gn     Global ids (size = \p n_entity)
 *
 */

void
PDM_isosurface_ln_to_gn_set
(
 PDM_isosurface_t    *isos,
 int                  i_part,
 PDM_mesh_entities_t  mesh_entity,
 int                  n_entity,
 PDM_g_num_t         *ln_to_gn
);


/**
 *
 * \brief Set group description
 *
 * \param [in]  isos                   \ref PDM_isosurface_t instance
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 * \param [in]  group_entity_idx       Index for group→entity connectivity (size = \p n_group)
 * \param [in]  group_entity           Group→entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  group_entity_ln_to_gn  Group→entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 *
 */

void
PDM_isosurface_group_set
(
 PDM_isosurface_t    *isos,
 int                  i_part,
 PDM_mesh_entities_t  entity_type,
 int                  n_group,
 int                 *group_entity_idx,
 int                 *group_entity,
 PDM_g_num_t         *group_entity_ln_to_gn
);


/**
 *
 * \brief Set partitioned mesh
 *
 * \param [in]  isos   \ref PDM_isosurface_t instance
 * \param [in]  pmesh  \ref PDM_part_mesh_t instance
 *
 */

void
PDM_isosurface_part_mesh_set
(
 PDM_isosurface_t *isos,
 PDM_part_mesh_t  *pmesh
);



// Nodal

/**
 *
 * \brief Set nodal mesh
 *
 * \param [in]  isos  \ref PDM_isosurface_t instance
 * \param [in]  pmn   \p PDM_part_mesh_nodal_t instance
 *
 */

void
PDM_isosurface_mesh_nodal_set
(
 PDM_isosurface_t      *isos,
 PDM_part_mesh_nodal_t *pmn
);


/* Block-distributed */

// Ngon

/**
 *
 * \brief Set block-distributed connectivity
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  connectivity_type  Type of connectivity
 * \param [in]  connect_idx        Index for connectivity (can be \p NULL for \p PDM_CONNECTIVITY_TYPE_EDGE_VTX)
 * \param [in]  connect            Connectivity (1-based *global* ids)
 *
 */

void
PDM_isosurface_dconnectivity_set
(
 PDM_isosurface_t        *isos,
 PDM_connectivity_type_t  connectivity_type,
 int                     *dconnect_idx,
 PDM_g_num_t             *dconnect
);


/**
 *
 * \brief Set block-distributed vertex coordinates
 *
 * \param [in]  isos        \ref PDM_isosurface_t instance
 * \param [in]  dvtx_coord  Vertex coordinates (size = 3 * *dn_vtx*)
 *
 */

void
PDM_isosurface_dvtx_coord_set
(
 PDM_isosurface_t *isos,
 double           *dvtx_coord
);


/**
 *
 * \brief Set entity block distribution index
 *
 * \param [in]  isos         \ref PDM_isosurface_t instance
 * \param [in]  entity_type  Type of mesh entity
 * \param [in]  distrib      Block-distribution (size = *n_rank* + 1)
 *
 */

void
PDM_isosurface_distrib_set
(
 PDM_isosurface_t    *isos,
 PDM_mesh_entities_t  entity_type,
 PDM_g_num_t         *distrib
);


/**
 *
 * \brief Set block-distributed group description
 *
 * \param [in]  isos                   \ref PDM_isosurface_t instance
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 * \param [in]  dgroup_entity_idx      Index for group→entity connectivity (size = \p n_group)
 * \param [in]  dgroup_entity          Group→entity connectivity (1-based global ids, size = \p dgroup_entity_idx[\p n_group])
 *
 */

void
PDM_isosurface_dgroup_set
(
 PDM_isosurface_t    *isos,
 PDM_mesh_entities_t  entity_type,
 int                  n_group,
 int                 *dgroup_entity_idx,
 PDM_g_num_t         *dgroup_entity
);


/**
 *
 * \brief Set distributed mesh
 *
 * \param [in]  isos   \ref PDM_isosurface_t instance
 * \param [in]  dmesh  \ref PDM_dmesh_t instance
 *
 */

void
PDM_isosurface_dmesh_set
(
 PDM_isosurface_t *isos,
 PDM_dmesh_t      *dmesh
);


// Nodal

/**
 *
 * \brief Set block-distributed nodal mesh
 *
 * \param [in]  isos  \ref PDM_isosurface_t instance
 * \param [in]  dmn   \ref PDM_dmesh_nodal_t instance
 *
 */

void
PDM_isosurface_dmesh_nodal_set
(
 PDM_isosurface_t  *isos,
 PDM_dmesh_nodal_t *dmn
);


/* --- Iso-surface settings --- */

/**
 *
 * \brief Add a requested set of iso-surfaces
 *
 * \param [in]  isos         \ref PDM_isosurface_t instance
 * \param [in]  kind         Iso-surface kind (discrete field, slice equation or function pointer)
 * \param [in]  n_isovalues  Number of iso-values to capture
 * \param [in]  isovalues    Iso-values to capture (size = \p n_isovalues)
 *
 * \return Iso-surface identifier
 *
 */
// elt_type, extract_kind, part_method?

int PDM_isosurface_add
(
 PDM_isosurface_t       *isos,
 PDM_iso_surface_kind_t  kind,
 int                     n_isovalues,
 double                 *isovalues
 );


/**
 *
 * \brief Set source field equation
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  coeff          Equation coefficients
 *
 * - \ref PDM_ISO_SURFACE_KIND_PLANE (3 coefficients):
 *   \f$\phi(x,y,z) = \texttt{coeff[0]} \cdot x + \texttt{coeff[1]} \cdot y + \texttt{coeff[2]} \cdot z\f$
 * - \ref PDM_ISO_SURFACE_KIND_SPHERE (4 coefficients):
 *   \f$\phi(x,y,z) = (x - \texttt{coeff[0]})^2 + (y - \texttt{coeff[1]})^2 + (z - \texttt{coeff[2]})^2 - \texttt{coeff[3]}^2\f$
 * - \ref PDM_ISO_SURFACE_KIND_ELLIPSE (7 coefficients):
 *   \f$\phi(x,y,z) = \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[6]}^2\f$
 * - \ref PDM_ISO_SURFACE_KIND_QUADRIC (10 coefficients):
 *   \f$\phi(x,y,z) = \texttt{coeff[6]} \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \texttt{coeff[7]} \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \texttt{coeff[8]} \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[9]}^2\f$
 *
 */

void PDM_isosurface_equation_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *coeff
);


/**
 *
 * \brief Set source field function
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  func           Function pointer
 *
 */
// ajouter PDM_isosurface_gradient_function_set ?

void PDM_isosurface_field_function_set
(
 PDM_isosurface_t                *isos,
 int                              id_isosurface,
 PDM_isosurface_field_function_t  func
 );


/**
 *
 * \brief Set field values
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  i_part         Partition identifier
 * \param [in]  field          Field values (size = *n_vtx*)
 *
 */
// on restreint à un seul field par id_isosurface?

void PDM_isosurface_field_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 int               i_part,
 double           *field
);


/**
 *
 * \brief Set gradient values
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  i_part         Partition identifier
 * \param [in]  gradient       Gradient values (size = 3 * *n_vtx*)
 *
 */

void PDM_isosurface_gradient_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 int               i_part,
 double           *gradient
);


/**
 *
 * \brief Set block-distributed field values
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  dfield         Field values (size = *dn_vtx*)
 *
 */
// Doit-on gérer cell-centered?

void PDM_isosurface_dfield_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *dfield
);


/**
 *
 * \brief Set block-distributed gradient values
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  dgradient      Gradient values (size = 3 * *dn_vtx*)
 *
 */

void PDM_isosurface_dgradient_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *dgradient
);


//---->>>> DANS CREATE OU ADD?
/**
 *
 * \brief Set the iso-surface mesh element type
 *
 * \param [in]  isos      \ref PDM_isosurface_t instance
 * \param [in]  elt_type  Desired element type for iso-surface mesh
 *
 * \note Admissible values for \p elt_type are:
 *   - \p PDM_MESH_NODAL_TRIA3
 *   - \p PDM_MESH_NODAL_QUAD4
 *   - \p PDM_MESH_NODAL_POLY_2D
 *
 * \note Maybe ajouter un "contouring_method_set"?
 *       - Dual contouring: tria ou quad
 *       - Marching tetra/cube/poly3d: poly2d ou tria
 *
 * --> déplacer dans "create"?
 */

void PDM_isosurface_elt_type_set
(
 PDM_isosurface_t     *isos,
 PDM_Mesh_nodal_elt_t  elt_type
);


/**
 *
 * \brief Set the iso-surface mesh redistribution kind
 *
 * \param [in]  isos          \ref PDM_isosurface_t instance
 * \param [in]  extract_kind  Redistribution kind
 *
 * \note Admissible values for \p extract_kind are:
 *   - \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE: the iso-surface is evenly redistributed (Default kind)
 *   - \ref PDM_EXTRACT_PART_KIND_LOCAL: the iso-surface is not redistributed (same partitioning as the input mesh)
 *
 *  --> déplacer dans "create"?
 */

void
PDM_isosurface_extract_kind_set
(
 PDM_isosurface_t        *isos,
 PDM_extract_part_kind_t  extract_kind
);


/**
 *
 * \brief Set the iso-surface mesh partitioning method
 *
 * \param [in]  isos        \ref PDM_isosurface_t instance
 * \param [in]  part_method  Partitioning method
 *
 * \note Not used if \ref PDM_EXTRACT_PART_KIND_LOCAL is passed to \ref PDM_isosurface_extract_kind_set
 *
 * --> déplacer dans "create"?
 */

void
PDM_isosurface_part_method_set
(
 PDM_isosurface_t *isos,
 PDM_split_dual_t  part_method
 );
//<<<<----


/**
 *
 * \brief Clear the constructed iso-surface meshes
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier (if < 0, all iso-surfaces are reset)
 *
 */

void
PDM_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
);


/* --- Compute --- */

/**
 *
 * \brief Compute the iso-surface mesh for all requested iso-values
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier (if < 0, all iso-surfaces are computed)
 *
 */

void
PDM_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
);


/**
 *
 * \brief Dump elapsed and CPU times
 *
 * \param [in]  isos  \ref PDM_isosurface_t instance
 *
 */

void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
);


/* --- Outputs --- */

// Partitioned

/**
 *
 * \brief Get iso-surface mesh connectivity for a given iso-value
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
 * \param [in]  i_isovalue         Iso-value identifier
 * \param [in]  i_part             Partition identifier
 * \param [in]  connectivity_type  Connectivity type
 * \param [out] connect_idx        Connectivity index
 * \param [out] connect            Connectivity
 * \param [in]  ownership          Ownership
 *
 * \return Number of leading entities
 *
 */
// si on ne gère que du nodal, peut-on retourner un connec_idx NULL? (si on est pas en POLY2D)
// connectivity_type n'a que deux valeurs admissibles:
//  FACE_VTX
//  EDGE_VTX (iso-line ou bien edges sur des groupes de surface)

int
PDM_isosurface_connectivity_get
(
 PDM_isosurface_t         *isos,
 int                       id_isosurface,
 int                       i_isovalue,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                     **connect_idx,
 int                     **connect,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Get coordinates of iso-surface vertices
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  i_isovalue     Iso-value identifier
 * \param [in]  i_part         Partition identifier
 * \param [out] vtx_coord      Vertex coordinates (size = 3 * *n_vtx*)
 * \param [in]  ownership      Ownership
 *
 * \return  Number of vertices
 *
 */

int
PDM_isosurface_vtx_coord_get
(
 PDM_isosurface_t  *isos,
 int                id_isosurface,
 int                i_isovalue,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
);


/**
 *
 * \brief Get global ids of iso-surface entities
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  i_isovalue     Iso-value identifier
 * \param [in]  i_part         Partition identifier
 * \param [in]  entity_type    Entity type
 * \param [out] ln_to_gn       Global ids
 * \param [in]  ownership      Ownership
 *
 * \return  Number of entities
 *
 */

int
PDM_isosurface_ln_to_gn_get
(
 PDM_isosurface_t     *isos,
 int                   id_isosurface,
 int                   i_isovalue,
 int                   i_part,
 PDM_mesh_entities_t   entity_type,
 PDM_g_num_t         **ln_to_gn,
 PDM_ownership_t       ownership
);


// Get groups?

// Sorties en pmesh/pmesh_nodal ?

// Block-distributed

/**
 *
 * \brief Get iso-surface block-distributed mesh connectivity for a given iso-value
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
 * \param [in]  i_isovalue         Iso-value identifier
 * \param [in]  connectivity_type  Connectivity type
 * \param [out] dconnect_idx       Connectivity index
 * \param [out] dconnect           Connectivity
 * \param [in]  ownership          Ownership
 *
 * \return Number of leading entities
 *
 */

int
PDM_isosurface_dconnectivity_get
(
 PDM_isosurface_t         *isos,
 int                       id_isosurface,
 int                       i_isovalue,
 PDM_connectivity_type_t   connectivity_type,
 int                     **dconnect_idx,
 PDM_g_num_t             **dconnect,
 PDM_ownership_t           ownership
);


/**
 *
 * \brief Get block-distributed coordinates of iso-surface vertices
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  i_isovalue     Iso-value identifier
 * \param [out] dvtx_coord     Vertex coordinates (size = 3 * *dn_vtx*)
 * \param [in]  ownership      Ownership
 *
 * \return  Number of vertices
 *
 */

int
PDM_isosurface_dvtx_coord_get
(
 PDM_isosurface_t  *isos,
 int                id_isosurface,
 int                i_isovalue,
 double           **dvtx_coord,
 PDM_ownership_t    ownership
);

// Communication graphs

// enable/disable construction of ptps?

/**
 * \brief Get \ref PDM_part_to_part_t instance to exchange data
 * between source mesh entities and iso-surface entities.
 *
 * \param [in]  isos         \ref PDM_isosurface_t instance
 * \param [in]  i_isovalue   Iso-value identifier
 * \param [in]  entity_type  Entity type
 * \param [out] ptp          Pointer to \ref PDM_part_to_part_t instance
 * \param [in ] ownership    Ownership for \p ptp
 *
 */

// PDM_MESH_ENTITY_VTX:  iso_vtx  → src_vtx  (gérer trace des ridges si mesh_dimension == 2 ?)
// PDM_MESH_ENTITY_EDGE: iso_edge → src_face (only group faces if mesh_dimension == 3)
// PDM_MESH_ENTITY_FACE: iso_face → src_cell

void PDM_isosurface_part_to_part_get
(
 PDM_isosurface_t     *isos,
 int                   i_isovalue,
 PDM_mesh_entities_t   entity_type,
 PDM_part_to_part_t  **ptp,
 PDM_ownership_t       ownership
 );


/**
 *
 * \brief Free a \ref PDM_isosurface_t instance
 *
 * \param [inout] isos  \ref PDM_isosurface_t instance
 */

void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
);



#ifdef  __cplusplus
}
#endif

#endif // PDM_ISOSURFACE_H
