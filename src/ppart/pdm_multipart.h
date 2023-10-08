/*
 * \file
 */

#ifndef __PDM_MULTIPART_H__
#define __PDM_MULTIPART_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_domain_interface.h"

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

typedef struct _pdm_multipart_t PDM_multipart_t;

/**
 * \enum PDM_part_size_t
 * \brief Use homogeneous or heterogeneous partition sizes (only for ParMetis method)
 */
typedef enum {
  PDM_PART_SIZE_HOMOGENEOUS   = 1, /*!< All requested partition have the same size */
  PDM_PART_SIZE_HETEROGENEOUS = 2, /*!< Each requested partition can have a portion (within 0. and 1.) of the mesh */
} PDM_part_size_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure. This method allows to split multiple zone
 *
 * \param [in]   n_zone           Number of zones in the original mesh
 * \param [in]   n_part           Number of partition per proc in each zone
 * \param [in]   merge_blocks     Merge or not the zones before splitting
 * \param [in]   split_method     Choice of library used to split the mesh
 * \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
 * \param [in]   part_weight      Weight (in %) of each partition in heterogeneous case if \ref PDM_part_size_t is set at PDM_PART_SIZE_HETEROGENEOUS
 * \param [in]   comm             PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_multipart_t object
 */

PDM_multipart_t *
PDM_multipart_create
(
 const int              n_zone,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_split_dual_t split_method,
 const PDM_part_size_t  part_size_method,
 const double          *part_fraction,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
);


/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh          Pointer on \ref PDM_dmesh_t containaing all distributed connectivities
 */

void PDM_multipart_register_block
(
 PDM_multipart_t   *multipart,
 const int          zone_id,
       PDM_dmesh_t *dmesh
);

/**
 *
 * \brief Set distributed mesh data for the input zone. The mesh is describe by nodal connectiviy
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_nodal    Pointer on \ref PDM_dmesh_nodal_t
 */

void PDM_multipart_register_dmesh_nodal
(
 PDM_multipart_t         *multipart,
 const int                zone_id,
       PDM_dmesh_nodal_t *dmesh_nodal
);

/**
 * \brief Set block
 *
 * \param [in]   multipart              Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone                 Id of zone
 * \param [in]   dn_cell                Number of distributed cells
 * \param [in]   dn_face                Number of distributed faces
 * \param [in]   dn_vtx                 Number of distributed vertices
 * \param [in]   n_face_group           Number of face groups
 * \param [in]   dcell_face_idx         Distributed cell face connectivity index or NULL
 *                                      (size : dn_cell + 1, numbering : 0 to n-1)
 * \param [in]   dcell_face             Distributed cell face connectivity or NULL
 *                                      (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 * \param [in]   dface_cell             Distributed face cell connectivity or NULL
 *                                      (size : 2 * dn_face, numbering : 1 to n)
 * \param [in]   dface_vtx_idx          Distributed face to vertex connectivity index
 *                                      (size : dn_face + 1, numbering : 0 to n-1)
 * \param [in]   dface_vtx              Distributed face to vertex connectivity
 *                                      (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 * \param [in]   dvtx_coord             Distributed vertex coordinates
 *                                      (size : 3*dn_vtx)
 * \param [in]   dface_group_idx        Index of distributed faces list of each group
 *                                      (size = n_face_group + 1) or NULL
 * \param [in]   dface_group            Distributed faces list of each group
 *                                      (size = dface_group[dface_group_idx[n_face_group]], numbering : 1 to n)
 *                                      or NULL
 *
 */
void
PDM_multipart_block_set
(
 PDM_multipart_t             *multipart,
 const int                    i_zone,
 const int                    dn_cell,
 const int                    dn_face,
 const int                    dn_vtx,
 const int                    n_face_group,
 const int                   *dcell_face_idx,
 const PDM_g_num_t           *dcell_face,
 const PDM_g_num_t           *dface_cell,
 const int                   *dface_vtx_idx,
 const PDM_g_num_t           *dface_vtx,
 const double                *dvtx_coord,
 const int                   *dface_group_idx,
 const PDM_g_num_t           *dface_group
);

/**
 *
 * \brief Set connecting data between all the zones
 *
 * \param [in]   multipart         Pointer to \ref PDM_multipart_t object
 * \param [in]   n_total_joins     Total number of interfaces
 * \param [in]   join_to_opposite  For each global join id, give the global id
 *                                   of the opposite join (size = n_total_joins)
 *
 * \note Join global id numbering must start at 0 and be continuous.
 */

void PDM_multipart_register_joins
(
 PDM_multipart_t *multipart,
 const int        n_total_joins,
 const int       *join_to_opposite
);

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
 * \param [in]   renum_cell_method     Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cacheblocking method :
 *                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
                                        n_vect_face, split_method]
 * \param [in]   renum_face_method     Choice of renumbering method for faces
 *
 */
void PDM_multipart_set_reordering_options
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const char      *renum_cell_method,
 const int       *renum_cell_properties,
 const char      *renum_face_method
);

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
 * \param [in]   renum_vtx_method      Choice of renumbering method for vertices
 *
 */

void PDM_multipart_set_reordering_options_vtx
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const char      *renum_vtx_method
);


/**
 *
 * \brief Construct the partitioned meshes on every zones
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 */
void
PDM_multipart_run_ppart
(
 PDM_multipart_t *multipart
);


/**
 * \brief Retreive the partitionned mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [out] pmesh                 Partitionned mesh
 * \param [in]  ownership             Who is responsible to free retreived data ?
 *
 */

void
PDM_multipart_get_part_mesh_nodal
(
      PDM_multipart_t        *multipart,
const int                     i_zone,
      PDM_part_mesh_nodal_t **pmesh_nodal,
      PDM_ownership_t         ownership
);

/**
 * \brief Retreive the partitionned mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [out] pmesh                 Partitionned mesh
 *
 */

// void
// PDM_multipart_get_part_mesh
// (
//        PDM_multipart_t  *multipart,
//  const int               i_zone,
//        PDM_part_mesh_t **pmesh
// );

void
PDM_multipart_domain_interface_shared_set
(
  PDM_multipart_t        *multipart,
  PDM_domain_interface_t *ditrf
);

/**
 *
 * \brief Returns the dimensions of a given partition
 */
void
PDM_multipart_part_dim_get
(
PDM_multipart_t *multipart,
const int        i_zone,
const int        i_part,
      int       *n_cell,
      int       *n_face,
      int       *n_face_part_bound,
      int       *n_vtx,
      int       *n_proc,
      int       *n_total_part,
      int       *s_cell_face,
      int       *s_face_vtx,
      int       *s_face_bound,
      int       *n_bound_groups
);


/**
 *
 * \brief Returns the connexion graph between partition for the request \ref PDM_bound_type_t
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_zone                Id of part
 * \param [in]  bound_type            Bound type
 * \param [out] ppart_bound_proc_idx  Partitioning boundary entities block distribution from processus (size = n_proc + 1)
 * \param [out] ppart_bound_part_idx  Partitioning boundary entities block distribution from partition (size = n_total_part + 1)
 * \param [out] ppart_bound           Partitioning boundary entities (size = 4 * n_entity_part_bound)
 * \param [in]  ownership             Choice of ownership of the resulting arrays \ref PDM_ownership_t
 */
void
PDM_multipart_part_graph_comm_get
(
 PDM_multipart_t    *multipart,
 const int           i_zone,
 const int           i_part,
 PDM_bound_type_t    bound_type,
 int               **ppart_bound_proc_idx,
 int               **ppart_bound_part_idx,
 int               **ppart_bound,
 PDM_ownership_t     ownership
);

/**
 *
 * \brief Returns the data arrays of a given partition
 *
 * \deprecated Use \ref PDM_multipart_part_connectivity_get instead
 */
void
PDM_multipart_part_val_get
(
PDM_multipart_t     *multipart,
const int            i_zone,
const int            i_part,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_bound_idx,
      int          **face_bound,
      PDM_g_num_t  **face_bound_ln_to_gn
);


/**
 *
 * \brief Returns the total number of part among all process
 */
int
PDM_multipart_part_tn_part_get
(
PDM_multipart_t                *multipart,
const int                       i_zone
);

/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  connectivity_type     Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect               Connectivity array (size = connect_idx[n_entity] )
 * \param [in]  connect_idx           Connectivity index (size = n_entity+1 )
 */
int
PDM_multipart_part_connectivity_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      PDM_connectivity_type_t   connectivity_type,
      int                     **connect,
      int                     **connect_idx,
      PDM_ownership_t           ownership
);


/**
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_zone                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t
 */
int
PDM_multipart_part_n_entity_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type
);

/**
 *
 * \brief Return size of leading connectivity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_ln_to_gn       Entity local numbering to global numbering (size = n_entity, numbering : 1 to n)
 * \param [in]  ownership             Ownership for entity_ln_to_gn ( \ref PDM_ownership_t )
 */
int
PDM_multipart_part_ln_to_gn_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      PDM_g_num_t         **entity_ln_to_gn,
      PDM_ownership_t       ownership
);


/**
 *
 * \brief Return number of entity on current partition ( n_entity )
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  entity_type           Entity kind \ref PDM_mesh_entities_t)
 * \param [out] entity_color          Entity color (only for specific renumbering option )
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
int
PDM_multipart_partition_color_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      int                 **entity_color,
      PDM_ownership_t       ownership
);

/**
 *
 * \brief Get array containing hyperplane color
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  hyperplane_color      Hyperplane color
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_hyperplane_color_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **hyperplane_color,
      PDM_ownership_t   ownership
);

/**
 *
 * \brief Get array containing thread color - Only if specific reordering (in paradigma plugins)
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  thread_color          Thread color
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_thread_color_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **thread_color,
      PDM_ownership_t   ownership
);


/**
 *
 * \brief Get array containing vtx_ghost_information, usefull to have a priority on vertex between 2 partitions
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  i_part                Id of part
 * \param [in]  vtx_ghost_information Integer that give the current priority of vertices on current partitions
 * \param [in]  ownership             Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_multipart_part_ghost_infomation_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **vtx_ghost_information,
      PDM_ownership_t   ownership
);


/**
 *
 * \brief Return times for a given zone
 * (NOT IMPLEMENTED)
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone         Id of current zone
 * \param [out]  elapsed        Elapsed time
 * \param [out]  cpu            CPU time
 * \param [out]  cpu_user       User CPU time
 * \param [out]  cpu_sys        System CPU time
 *
 */
void
PDM_multipart_time_get
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 double         **elapsed,
 double         **cpu,
 double         **cpu_user,
 double         **cpu_sys
);


/**
 *
 * \brief Free the structure
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 */
void
PDM_multipart_free
(
 PDM_multipart_t *multipart
);

/**
 *
 * \brief Get the vertex coordinates on current i_zone, i_part partition and return number of vertices
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone         Id of current zone
 * \param [in]   i_part         Id of part
 * \param [out]  vtx_coord      Vertex coordinate (size = 3 * n_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 *
 */
int
PDM_multipart_part_vtx_coord_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      double                  **vtx_coord,
      PDM_ownership_t           ownership
);


/**
 *
 * \brief Get the bound description for the entity
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone         Id of current zone
 * \param [in]   i_part         Id of part
 * \param [in]   bound_type     Bound type \ref PDM_bound_type_t
 * \param [out]  n_bound        Number of bound for bound_type
 * \param [out]  bound_idx      Entity group index (size = n_bound )
 * \param [out]  bound          Entity id for each group (size = bound_idx[n_bound])
 * \param [out]  bound_ln_to_gn Entity global numbering for each group (size = bound_idx[n_bound])
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 *
 */
void PDM_multipart_bound_get
(
 PDM_multipart_t   *multipart,
 const int          i_zone,
 const int          i_part,
 PDM_bound_type_t   bound_type,
 int               *n_bound,
 int              **bound_idx,
 int              **bound,
 PDM_g_num_t      **bound_ln_to_gn,
 PDM_ownership_t    ownership
);


/**
 *
 * \brief Return statistics
 *
 * \param [in]   ppart                          Pointer to \ref PDM_part object
 * \param [out]  cells_average                  average of cells number
 * \param [out]  cells_median                   median of cells number
 * \param [out]  cells_std_deviation            standard deviation of cells number
 * \param [out]  cells_min                      minimum of cells nummber
 * \param [out]  cells_max                      maximum of cells nummber
 * \param [out]  bound_part_faces_average       average of partitioning boundary faces
 * \param [out]  bound_part_faces_median        median of partitioning boundary faces
 * \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 * \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 * \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 *
 */
void
PDM_multipart_stat_get
(
 PDM_multipart_t  *multipart,
 int               i_zone,
 int              *cells_average,
 int              *cells_median,
 double           *cells_std_deviation,
 int              *cells_min,
 int              *cells_max,
 int              *bound_part_faces_average,
 int              *bound_part_faces_median,
 double           *bound_part_faces_std_deviation,
 int              *bound_part_faces_min,
 int              *bound_part_faces_max,
 int              *bound_part_faces_sum
);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MULTIPART_H__ */
