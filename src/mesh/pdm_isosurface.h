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

typedef struct _pdm_isosurface_t PDM_isosurface_t;


/**
 * \brief User-defined field function.
 *
 * \param [in]  x      X-coordinate
 * \param [in]  y      Y-coordinate
 * \param [in]  z      Z-coordinate
 * \param [out] value  Field value at point (x, y, z)
 *
 */
typedef void (*_pdm_isosurface_field_function_t)
(
 const double  x,
 const double  y,
 const double  z,
 double       *value
);
typedef _pdm_isosurface_field_function_t PDM_isosurface_field_function_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


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
  PDM_MPI_Comm             comm,
  int                      mesh_dimension
);



/**
 *
 * \brief Set isosurface tolerance. May improve resulting mesh quality.
 *
 * \param [in]  isos       \ref PDM_isosurface_t instance
 * \param [in]  tolerance  Field tolerance (default at 0)
 *
 */
void
PDM_isosurface_set_tolerance
(
  PDM_isosurface_t *isos,
  double            tolerance
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
 * \param [in]  n_entity           Local number of leading entities
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
  int                      n_entity,
  int                     *connect_idx,
  int                     *connect
);


/**
 *
 * \brief Set vertex coordinates
 *
 * \param [in]  isos       \ref PDM_isosurface_t instance
 * \param [in]  i_part     Partition identifier
 * \param [in]  n_vtx      Local number of vertices
 * \param [in]  vtx_coord  Vertex coordinates (size = 3 * \p n_vtx)
 *
 */

void
PDM_isosurface_vtx_coord_set
(
  PDM_isosurface_t *isos,
  int               i_part,
  int               n_vtx,
  double           *vtx_coord
);


/**
 *
 * \brief Set global ids
 *
 * \param [in]  isos         \ref PDM_isosurface_t instance
 * \param [in]  i_part       Partition identifier
 * \param [in]  entity_type  Type of mesh entity
 * \param [in]  ln_to_gn     Global ids
 *
 */

void
PDM_isosurface_ln_to_gn_set
(
  PDM_isosurface_t    *isos,
  int                  i_part,
  PDM_mesh_entities_t  entity_type,
  PDM_g_num_t         *ln_to_gn
);


/**
 *
 * \brief Set number of groups
 *
 * \param [in]  isos                   \ref PDM_isosurface_t instance
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  n_group                Number of groups
 *
 */

void
PDM_isosurface_n_group_set
(
  PDM_isosurface_t    *isos,
  PDM_mesh_entities_t  entity_type,
  int                  n_group
);


/**
 *
 * \brief Set group description
 *
 * \param [in]  isos                   \ref PDM_isosurface_t instance
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Type of mesh entity
 * \param [in]  group_entity_idx       Index for group→entity connectivity (size = \p n_group + 1)
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
 * \param [in]  pmn   \ref PDM_part_mesh_nodal_t instance
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
 * \param [in]  dgroup_entity_idx      Index for group→entity connectivity (size = \p n_group + 1)
 * \param [in]  dgroup_entity          Group→entity connectivity (1-based global ids, size = \p dgroup_entity_idx[\p n_group])
 *
 */

void
PDM_isosurface_dgroup_set
(
  PDM_isosurface_t    *isos,
  PDM_mesh_entities_t  entity_type,
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

int
PDM_isosurface_add
(
  PDM_isosurface_t       *isos,
  PDM_iso_surface_kind_t  kind,
  int                     n_isovalues,
  double                 *isovalues
);


/**
 *
 * \brief Reset isovalues of given isosurface.
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  n_isovalues    Number of iso-values to capture
 * \param [in]  isovalues      Iso-values to capture (size = \p n_isovalues)
 *
 */

void
PDM_isosurface_set_isovalues
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  int               n_isovalues,
  double           *isovalues
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
 * 
 * 
 * - \ref PDM_ISO_SURFACE_KIND_SPHERE (4 coefficients):
 *   \f$\phi(x,y,z) = (x - \texttt{coeff[0]})^2 + (y - \texttt{coeff[1]})^2 + (z - \texttt{coeff[2]})^2 - \texttt{coeff[3]}^2\f$
 * 
 * 
 * - \ref PDM_ISO_SURFACE_KIND_ELLIPSE (7 coefficients):
 *   \f$\phi(x,y,z) = \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[6]}^2\f$
 * 
 * 
 * - \ref PDM_ISO_SURFACE_KIND_QUADRIC (10 coefficients):
 *   \f$\phi(x,y,z) = \texttt{coeff[6]} \left(\frac{x - \texttt{coeff[0]}}{\texttt{coeff[3]}}\right)^2 + \texttt{coeff[7]} \left(\frac{y - \texttt{coeff[1]}}{\texttt{coeff[4]}}\right)^2 + \texttt{coeff[8]} \left(\frac{z - \texttt{coeff[2]}}{\texttt{coeff[5]}}\right)^2 - \texttt{coeff[9]}^2\f$
 *
 */

void
PDM_isosurface_equation_set
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

void
PDM_isosurface_field_function_set
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

void
PDM_isosurface_field_set
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  int               i_part,
  double           *field
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

void
PDM_isosurface_dfield_set
(
  PDM_isosurface_t *isos,
  int               id_isosurface,
  double           *dfield
);


/**
 * \brief Set the isosurface redistribution options
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  extract_kind    Redistribution kind
 * \param [in]  part_method     Partitioning method (only used if \p extract_kind is set to PDM_EXTRACT_PART_KIND_REEQUILIBRATE)
 *
 * \note Admissible values for \p extract_kind are:
 *   - \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE: the iso-surface is evenly redistributed (Default kind)
 *   - \ref PDM_EXTRACT_PART_KIND_LOCAL: the iso-surface is not redistributed (same partitioning as the input mesh)
 *
 */

void
PDM_isosurface_redistribution_set
(
  PDM_isosurface_t        *isos,
  PDM_extract_part_kind_t  extract_kind,
  PDM_split_dual_t         part_method
);


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


/**
 * \brief Set the number of partitions in the isosurface mesh (Optional).
 *
 * \param [in]  isos        \ref PDM_isosurface_t instance
 * \param [in]  n_part_out  Number of partitions
 *
 * \warning This function must be called prior to \ref PDM_isosurface_compute.
 *
 * \note By default, the number of partitions in the isosurface mesh is set to
 *  - 1 in \ref PDM_EXTRACT_PART_KIND_REEQUILIBRATE mode
 *  - the number of partitions in the source mesh in \ref PDM_EXTRACT_PART_KIND_LOCAL mode (mandatory)
 */

void
PDM_isosurface_n_part_out_set
(
  PDM_isosurface_t *isos,
  int               n_part_out
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

/*
 * TODO:
 *   - PDM_isosurface_n_part_get ???
 */



/**
 *
 * \brief Get iso-surface mesh connectivity
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
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
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  PDM_g_num_t         **ln_to_gn,
  PDM_ownership_t       ownership
);


// Groups

/**
 *
 * \brief Set group description
 *
 * \param [in]  isos                   \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface          Iso-surface identifier
 * \param [in]  i_part                 Partition identifier
 * \param [in]  entity_type            Entity type
 * \param [out] group_entity_idx       Index for group→entity connectivity (size = \p n_group + 1)
 * \param [out] group_entity           Group→entity connectivity (1-based local ids, size = \p group_entity_idx[\p n_group])
 * \param [out] group_entity_ln_to_gn  Group→entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  ownership              Ownership
 *
 * \return Number of groups
 *
 */

int
PDM_isosurface_group_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **group_entity_idx,
  int                 **group_entity,
  PDM_g_num_t         **group_entity_ln_to_gn,
  PDM_ownership_t       ownership
);


// Sorties en part_mesh_nodal ?

// Block-distributed

/**
 *
 * \brief Get iso-surface block-distributed mesh connectivity
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
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
  PDM_connectivity_type_t   connectivity_type,
  int                     **dconnect_idx,
  PDM_g_num_t             **dconnect,
  PDM_ownership_t           ownership
);


/**
 *
 * \brief Get iso-surface entities parent gnum
 *
 * \param [in]  isos          \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface Iso-surface identifier
 * \param [in]  entity_type   Entity type
 * \param [out] dparent_idx   Parent index
 * \param [out] dparent       Parent gnum
 * \param [in]  ownership     Ownership
 *
 * \return Number of leading entities
 *
 */

int
PDM_isosurface_parent_gnum_get
(
  PDM_isosurface_t     *isos,
  int                   id_iso,
  PDM_mesh_entities_t   entity_type,
  int                 **parent_idx,
  PDM_g_num_t         **parent_gnum,
  PDM_ownership_t       ownership
);


/**
 *
 * \brief Get iso-surface parent weight for iso vertices
 *
 * \param [in]  isos                \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface       Iso-surface identifier
 * \param [out] dvtx_parent_idx     Index for parent weights
 * \param [out] dvtx_parent_weight  Parent weight
 * \param [in]  ownership           Ownership
 *
 * \return  Number of iso-surface vertices
 *
 */

int
PDM_isosurface_dvtx_parent_weight_get
(
  PDM_isosurface_t     *isos,
  int                   id_iso,
  int                 **dvtx_parent_idx,
  double              **dvtx_parent_weight,
  PDM_ownership_t       ownership
);

/**
 *
 * \brief TODO
 *
 */
int
PDM_isosurface_dvtx_protocol_get
(
  PDM_isosurface_t         *isos
);




/**
 *
 * \brief Get coordinates of block-distributed iso-surface vertices
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
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
  double           **dvtx_coord,
  PDM_ownership_t    ownership
);


/**
 *
 * \brief Get block distribution
 *
 * \param [in]  isos          \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface Iso-surface identifier
 * \param [in]  entity_type   Entity type
 * \param [out] distribution  Entity distribution
 *
 * \return Number of groups
 *
 */

void
PDM_isosurface_distrib_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  PDM_mesh_entities_t   entity_type,
  PDM_g_num_t         **distribution
);


// Groups

/**
 *
 * \brief Get block-distributed group description
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
 * \param [in]  entity_type        Entity type
 * \param [out] dgroup_entity_idx  Index for group→entity connectivity (size = \p n_group + 1)
 * \param [out] dgroup_entity      Group→entity connectivity (group-specific global ids, size = \p group_entity_idx[\p n_group])
 * \param [in]  ownership          Ownership
 *
 * \return Number of groups
 *
 */

int
PDM_isosurface_dgroup_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  PDM_mesh_entities_t   entity_type,
  int                 **dgroup_entity_idx,
  PDM_g_num_t         **dgroup_entity,
  PDM_ownership_t       ownership
);


/**
 * \brief Get isovalue
 *
 * \param [in]  isos                 \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface        Iso-surface identifier
 * \param [in]  i_part               Partition identifier
 * \param [in]  entity_type          Entity type
 * \param [in]  isovalue_entity_idx  Index for isovalue→entity connectivity (size = \p n_isovalue + 1)
 * \param [in]  ownership            Ownership
 *
 * \return Number of isovalues
 *
 * \warning comment on fait en block-distribué?
 *
 */

int
PDM_isosurface_isovalue_entity_idx_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **isovalue_entity_idx,
  PDM_ownership_t       ownership
);


/**
 *
 * \brief Get local parents of iso-surface entities
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
 * \param [in]  i_part             Partition identifier
 * \param [in]  entity_type        Entity type
 * \param [out] entity_parent_idx  Index for isosurface entity → parent connectivity (size = \p n_entity + 1)
 * \param [out] entity_parent      Isosurface entity → parent connectivity (size = \p entity_parent_idx[\p n_entity])
 * \param [in]  ownership          Ownership
 *
 * \return  Number of entities
 *
 * \note This function can only be called if \p extract_kind has been set to \ref PDM_EXTRACT_PART_KIND_LOCAL
 * in \ref PDM_isosurface_redistribution_set.
 * The nature of parent entities depends on \p entity_type as follows:
 *   - PDM_MESH_ENTITY_VTX  : parents are vertices
 *   - PDM_MESH_ENTITY_EDGE : parents are faces
 *   - PDM_MESH_ENTITY_FACE : parents are cells
 *
 * \warning comment on fait en block-distribué?
 *
 */

int
PDM_isosurface_local_parent_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  int                   i_part,
  PDM_mesh_entities_t   entity_type,
  int                 **entity_parent_idx,
  int                 **entity_parent,
  PDM_ownership_t       ownership
);


/**
 *
 * \brief Get interpolation weights of iso-surface vertices
 *
 * \param [in]  isos               \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface      Iso-surface identifier
 * \param [in]  i_part             Partition identifier
 * \param [out] vtx_parent_idx     Index for interpolation weights
 * \param [out] vtx_parent_weight  Interpolation weights
 * \param [in]  ownership          Ownership
 *
 * \warning These weights are only computed if the construction
 * of the vertex Part-to-Part has been enabled (see \ref PDM_isosurface_enable_part_to_part).
 *
 * \return  Number of iso-surface vertices
 *
 */

int
PDM_isosurface_vtx_parent_weight_get
(
  PDM_isosurface_t  *isos,
  int                id_isosurface,
  int                i_part,
  int              **vtx_parent_idx,
  double           **vtx_parent_weight,
  PDM_ownership_t    ownership
);

// Communication graphs

/**
 * \brief Enable construction of a communication graph between source mesh entities and iso-surface entities.
 *
 * \param [in]  isos              \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface     Iso-surface identifier
 * \param [in]  entity_type       Entity type
 * \param [in]  unify_parent_info Get all parent over all procs (not implemented)
 *
 * \warning This function must be called prior to \ref PDM_isosurface_compute
 *
 */

void
PDM_isosurface_enable_part_to_part
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
  PDM_mesh_entities_t   entity_type,
  int                   unify_parent_info
);

/**
 * \brief Get \ref PDM_part_to_part_t instance to exchange data
 * between source mesh entities and iso-surface entities.
 *
 * \param [in]  isos           \ref PDM_isosurface_t instance
 * \param [in]  id_isosurface  Iso-surface identifier
 * \param [in]  entity_type    Entity type
 * \param [out] ptp            Pointer to \ref PDM_part_to_part_t instance
 * \param [in ] ownership      Ownership for \p ptp
 *
 */

// PDM_MESH_ENTITY_VTX:  iso_vtx  → src_vtx  (gérer trace des ridges si mesh_dimension == 2 ?)
// PDM_MESH_ENTITY_EDGE: iso_edge → src_face (only group faces if mesh_dimension == 3)
// PDM_MESH_ENTITY_FACE: iso_face → src_cell

void
PDM_isosurface_part_to_part_get
(
  PDM_isosurface_t     *isos,
  int                   id_isosurface,
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



/*
 * ==========
 * Algorithms
 */

void
PDM_isosurface_marching_algo
(
  PDM_isosurface_t        *isos,
  int                      id_iso
);


void
PDM_isosurface_ngon_algo
(
  PDM_isosurface_t        *isos,
  int                      id_iso
);

/*
 * End Algorithms
 * ==============
 */



#ifdef  __cplusplus
}
#endif

#endif // PDM_ISOSURFACE_H
