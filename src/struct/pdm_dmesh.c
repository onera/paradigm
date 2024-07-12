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

/*============================================================================
 * Interface structure to represent a distributed mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_dmesh_priv.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_dmesh.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a distributed mesh structure
 *
 * \param [in]   dn_cell   Number of distributed cells
 * \param [in]   dn_face   Number of distributed faces
 * \param [in]   dn_edge   Number of distributed edges
 * \param [in]   dn_vtx    Number of distributed vertices
 * \param [in]   comm      PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_dmesh_t object
 */
PDM_dmesh_t*
PDM_dmesh_create
(
       PDM_ownership_t owner,
 const int             dn_cell,
 const int             dn_face,
 const int             dn_edge,
 const int             dn_vtx,
       PDM_MPI_Comm    comm
)
{
  PDM_dmesh_t *dmesh = (PDM_dmesh_t *) malloc(sizeof(PDM_dmesh_t));

  dmesh->comm              = comm;
  dmesh->owner             = owner;
  dmesh->results_is_getted = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t) );

  dmesh->dn_cell           = dn_cell;
  dmesh->dn_face           = dn_face;
  dmesh->dn_edge           = dn_edge;
  dmesh->dn_vtx            = dn_vtx;

  dmesh->n_g_cell          = 0;
  dmesh->n_g_face          = 0;
  dmesh->n_g_edge          = 0;
  dmesh->n_g_vtx           = 0;

  PDM_g_num_t _dn_cell = dmesh->dn_cell;
  PDM_g_num_t _dn_face = dmesh->dn_face;
  PDM_g_num_t _dn_edge = dmesh->dn_edge;
  PDM_g_num_t _dn_vtx  = dmesh->dn_vtx;

  PDM_MPI_Allreduce(&_dn_cell, &dmesh->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_face, &dmesh->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_edge, &dmesh->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(&_dn_vtx , &dmesh->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  dmesh->cell_distrib      = NULL;
  dmesh->face_distrib      = NULL;
  dmesh->edge_distrib      = NULL;
  dmesh->vtx_distrib       = NULL;

  dmesh->_dvtx_coord       = NULL;
  dmesh->is_owner_vtx_coord  = PDM_TRUE;

  dmesh->dconnectivity         = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dconnectivity_idx     = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_connectivity = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    dmesh->is_owner_connectivity[i] = PDM_FALSE;
    dmesh->dconnectivity        [i] = NULL;
    dmesh->dconnectivity_idx    [i] = NULL;
  }

  dmesh->dbound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t *) );
  dmesh->dbound_idx      = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         *) );
  dmesh->is_owner_bound  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
    dmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    dmesh->is_owner_bound[i] = PDM_FALSE;
    dmesh->dbound        [i] = NULL;
    dmesh->dbound_idx    [i] = NULL;
  }

  dmesh->is_computed_g_extents = PDM_FALSE;

  return dmesh;
}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh    Pointer to \ref PDM_dmesh_t object
 * \param [out]   dn_cell  Number of distributed cells
 * \param [out]   dn_face  Number of distributed faces
 * \param [out]   dn_edge  Number of distributed edges
 * \param [out]   dn_vtx   Number of distributed vertices
 */
void
PDM_dmesh_dims_get
(
 PDM_dmesh_t *dmesh,
 int         *dn_cell,
 int         *dn_face,
 int         *dn_edge,
 int         *dn_vtx
)
{
  *dn_cell = dmesh->dn_cell;
  *dn_face = dmesh->dn_face;
  *dn_edge = dmesh->dn_edge;
  *dn_vtx  = dmesh->dn_vtx;
}


/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]    entity_type       Kind of entity
 * \param [in]    dn_cell           Number of distributed cells
 */
void
PDM_dmesh_dn_entity_set
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type,
 int                  dn_entity
)
{
  if(entity_type == PDM_MESH_ENTITY_CELL) {
    dmesh->dn_cell = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_FACE) {
    dmesh->dn_face = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_EDGE) {
    dmesh->dn_edge = dn_entity;
  } else if(entity_type == PDM_MESH_ENTITY_VTX) {
    dmesh->dn_vtx = dn_entity;
  }
}

/**
 *
 * \brief Get the dimensions of the distributed mesh
 *
 * \param [in]    dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]    entity_type       Kind of mesh entity \ref PDM_mesh_entities_t
 */
int
PDM_dmesh_dn_entity_get
(
 PDM_dmesh_t         *dmesh,
 PDM_mesh_entities_t  entity_type
)
{
  if(entity_type == PDM_MESH_ENTITY_CELL) {
    return dmesh->dn_cell;
  } else if(entity_type == PDM_MESH_ENTITY_FACE) {
    return dmesh->dn_face;
  } else if(entity_type == PDM_MESH_ENTITY_EDGE) {
    return dmesh->dn_edge;
  } else if(entity_type == PDM_MESH_ENTITY_VTX) {
    return dmesh->dn_vtx;
  } else {
    return -1;
  }
}

/**
 *
 * \brief Get the distributed coordinates array
 *
 * \param [in]   dmesh          Pointer to \ref PDM_dmesh_t object
 * \param [out]  dvtx_coord     Vertex coordinate (size = 3 * dn_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_dmesh_vtx_coord_get
(
 PDM_dmesh_t      *dmesh,
 double          **dvtx_coord,
 PDM_ownership_t   ownership
)
{
  *dvtx_coord      = dmesh->_dvtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

/**
 *
 * \brief Set the distributed coordinates array
 *
 * \param [in]   dmesh          Pointer to \ref PDM_dmesh_t object
 * \param [out]  dvtx_coord     Vertex coordinate (size = 3 * dn_vtx)
 * \param [in]   ownership      Ownership for color ( \ref PDM_ownership_t )
 */
void
PDM_dmesh_vtx_coord_set
(
 PDM_dmesh_t      *dmesh,
 double           *dvtx_coord,
 PDM_ownership_t   ownership
)
{
  dmesh->_dvtx_coord = dvtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}


/**
 *
 * \brief Set the distributed connectivity array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  connectivity_type Connectivity kind \ref PDM_connectivity_type_t
 * \param [in]  connect           Connectivity array (size = connect_idx[n_entity] )
 * \param [in]  connect_idx       Connectivity index (size = n_entity+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 */
void
PDM_dmesh_connectivity_set
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t              *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
)
{
  assert(dmesh != NULL);
  dmesh->dconnectivity    [connectivity_type] = connect;
  dmesh->dconnectivity_idx[connectivity_type] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}

/**
 *
 * \brief Get the distributed connectivity array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  connectivity_type Connectivity kind \ref PDM_connectivity_type_t
 * \param [out] connect           Connectivity array (size = connect_idx[n_entity] )
 * \param [out] connect_idx       Connectivity index (size = n_entity+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 * \return Number of element of entity kind
 */
int
PDM_dmesh_connectivity_get
(
 PDM_dmesh_t              *dmesh,
 PDM_connectivity_type_t   connectivity_type,
 PDM_g_num_t             **connect,
 int                     **connect_idx,
 PDM_ownership_t           ownership
)
{
  assert(dmesh != NULL);

  *connect     = dmesh->dconnectivity    [connectivity_type];
  *connect_idx = dmesh->dconnectivity_idx[connectivity_type];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }

  int dn_entity = -1;
  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    dn_entity = dmesh->dn_cell;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    dn_entity = dmesh->dn_face;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    dn_entity = dmesh->dn_edge;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    dn_entity = dmesh->dn_vtx;
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_ELMT_VTX )
  {
    dn_entity = -1;
  }

  return dn_entity;
}



/**
 *
 * \brief Get the distributed connectivity bound array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  bound_type        Connectivity kind \ref PDM_bound_type_t
 * \param [out] connect           Connectivity array (size = connect_idx[n_bound] )
 * \param [out] connect_idx       Connectivity index (size = n_bound+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 * \return Number of group for the requested entity (n_bound)
 */
int
PDM_dmesh_bound_get
(
 PDM_dmesh_t       *dmesh,
 PDM_bound_type_t   bound_type,
 PDM_g_num_t      **connect,
 int              **connect_idx,
 PDM_ownership_t    ownership
)
{
  assert(dmesh != NULL);

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }

  // assert(dmesh->dbound[bound_type] != NULL);

  *connect     = dmesh->dbound    [bound_type];
  *connect_idx = dmesh->dbound_idx[bound_type];

  return dmesh->n_group_bnd[bound_type];
}


/**
 *
 * \brief Get the distribution of requested entity
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [out] entity_type       entity_type       Kind of mesh entity \ref PDM_mesh_entities_t
 * \param [out] distrib           Distribution array (size = n_rank+1, numbering start at 0)
 * \return Number of process on this distribution ( n_rank )
 */
int
PDM_dmesh_distrib_get
(
 PDM_dmesh_t              *dmesh,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **distrib
)
{
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     *distrib = dmesh->cell_distrib;
     break;
   case PDM_MESH_ENTITY_FACE:
     *distrib = dmesh->face_distrib;
     break;
   case PDM_MESH_ENTITY_EDGE:
     *distrib = dmesh->edge_distrib;
     break;
   case PDM_MESH_ENTITY_VTX:
     *distrib = dmesh->vtx_distrib;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "PDM_dmesh_distrib_get invalid entity_type %d\n", entity_type);
    break;
   }
   int n_rank;
   PDM_MPI_Comm_size(dmesh->comm, &n_rank);
   return n_rank;
}


/**
 *
 * \brief Free
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 *
 */
void
PDM_dmesh_free
(
 PDM_dmesh_t         *dmesh
)
{
  if (dmesh == NULL) {
    return;
  }
  dmesh->dn_cell           = 0;
  dmesh->dn_face           = 0;
  dmesh->dn_edge           = 0;
  dmesh->dn_vtx            = 0;

  if(dmesh->is_owner_vtx_coord ==  PDM_TRUE) {
    if(dmesh->_dvtx_coord != NULL) {
     PDM_free(dmesh->_dvtx_coord);
    }
  }
  dmesh->_dvtx_coord       = NULL;

  // On doit gérer les cas ou la structure est partagé en python et auquel cas
  // On est owner des resultats et il faut free le reste
  // Donc il faut un is_getted + is_owner pour s'en sortir

  if(( dmesh->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmesh->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE)){
    for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {

      if(dmesh->is_owner_connectivity[i] == PDM_TRUE) {

        if(dmesh->dconnectivity[i] != NULL){
         PDM_free(dmesh->dconnectivity[i]);
        }
        if(dmesh->dconnectivity_idx[i] != NULL){
         PDM_free(dmesh->dconnectivity_idx[i]);
        }
        dmesh->dconnectivity    [i] = NULL;
        dmesh->dconnectivity_idx[i] = NULL;

      }
    }

    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {

      if(dmesh->is_owner_bound[i] == PDM_TRUE) {

        //printf(" dmesh_free :: %i \n", i);
        if(dmesh->dbound[i] != NULL) {
         PDM_free(dmesh->dbound[i]);
        }
        if(dmesh->dbound_idx[i] != NULL){
         PDM_free(dmesh->dbound_idx[i]);
        }
        dmesh->dbound    [i] = NULL;
        dmesh->dbound_idx[i] = NULL;

      }
    }
  }

 PDM_free(dmesh->results_is_getted    );
 PDM_free(dmesh->dconnectivity        );
 PDM_free(dmesh->dconnectivity_idx    );
 PDM_free(dmesh->is_owner_connectivity);

 PDM_free(dmesh->dbound        );
 PDM_free(dmesh->dbound_idx    );
 PDM_free(dmesh->is_owner_bound);

  /* This result is never getted so we can free them */
  if(dmesh->cell_distrib != NULL) {
   PDM_free(dmesh->cell_distrib);
    dmesh->cell_distrib = NULL;
  }

  if(dmesh->face_distrib != NULL) {
   PDM_free(dmesh->face_distrib);
    dmesh->face_distrib = NULL;
  }

  if(dmesh->edge_distrib != NULL) {
   PDM_free(dmesh->edge_distrib);
    dmesh->edge_distrib = NULL;
  }

  if(dmesh->vtx_distrib != NULL) {
   PDM_free(dmesh->vtx_distrib);
    dmesh->vtx_distrib = NULL;
  }

 PDM_free(dmesh);
}



/**
 *
 * \brief Compute the bounding box extend of current distributed mesh
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \return Extents of current mesh (6 components Xmin, Ymin, Zmin, Xmax, Ymax, Zmax )
 *
 */
const double *
PDM_dmesh_global_extents_get
(
 PDM_dmesh_t         *dmesh
 )
{
  if (dmesh->is_computed_g_extents == PDM_FALSE) {

    double l_min[3] = { HUGE_VAL};
    double l_max[3] = {-HUGE_VAL};

    for (int i = 0; i < dmesh->dn_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        double x = dmesh->_dvtx_coord[3*i + j];
        l_min[j] = PDM_MIN(l_min[j], x);
        l_max[j] = PDM_MAX(l_max[j], x);
      }
    }

    PDM_MPI_Allreduce(l_min, dmesh->g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, dmesh->comm);
    PDM_MPI_Allreduce(l_max, dmesh->g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, dmesh->comm);

    dmesh->is_computed_g_extents = PDM_TRUE;
  }

  return dmesh->g_extents;
}


/**
 *
 * \brief Get the distributed connectivity bound array
 *
 * \param [in]  dmesh             Pointer to \ref PDM_dmesh_t object
 * \param [in]  bound_type        Connectivity kind \ref PDM_bound_type_t
 * \param [in]  n_bound           Number of bound for current entity
 * \param [in]  connect           Connectivity array (size = connect_idx[n_bound] )
 * \param [in]  connect_idx       Connectivity index (size = n_bound+1 )
 * \param [in]  ownership         Choice of ownership of the input arrays \ref PDM_ownership_t
 */
void
PDM_dmesh_bound_set
(
 PDM_dmesh_t      *dmesh,
 PDM_bound_type_t  bound_type,
 int               n_bound,
 PDM_g_num_t      *connect,
 int              *connect_idx,
 PDM_ownership_t   ownership
)
{
  assert(dmesh != NULL);

  dmesh->n_group_bnd[bound_type] = n_bound;
  dmesh->dbound     [bound_type] = connect;
  dmesh->dbound_idx [bound_type] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    dmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if(ownership == PDM_OWNERSHIP_KEEP) {
    dmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_dmesh_find_topological_ridges
(
  PDM_MPI_Comm   comm,
  PDM_g_num_t   *distrib_face,
  int           *dface_vtx_idx,
  PDM_g_num_t   *dface_vtx,
  int            n_group_face,
  int           *dgroup_face_idx,
  PDM_g_num_t   *dgroup_face,
  PDM_g_num_t  **out_distrib_ridge,
  PDM_g_num_t  **out_dridge_vtx,
  int           *out_n_group_ridge,
  int          **out_dgroup_edge_idx,
  PDM_g_num_t  **out_dgroup_edge,
  int          **out_dridge_face_group_idx,
  int          **out_dridge_face_group
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Computation of edges
   */
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  /*
   * Extract faces only link to a group
   */
  PDM_g_num_t         *distrib_extract_face     = NULL;
  PDM_g_num_t         *dparent_face             = NULL;
  int                 *dextract_face_vtx_idx    = NULL;
  PDM_g_num_t         *dextract_face_vtx        = NULL;
  PDM_block_to_part_t *btp_face_to_extract_face = NULL;
  PDM_g_num_t         *distrib_extract_vtx      = NULL;
  PDM_g_num_t         *dparent_vtx              = NULL;

  PDM_dconnectivity_to_extract_dconnectivity(comm,
                                             dgroup_face_idx[n_group_face],
                                             dgroup_face,
                                             distrib_face,
                                             dface_vtx_idx,
                                             dface_vtx,
                                             &distrib_extract_face,
                                             &dparent_face,
                                             &dextract_face_vtx_idx,
                                             &dextract_face_vtx,
                                             &btp_face_to_extract_face,
                                             &distrib_extract_vtx,
                                             &dparent_vtx);

  int dn_face_extract = distrib_extract_face[i_rank+1] - distrib_extract_face[i_rank];

  int n_edge_elt_tot = dextract_face_vtx_idx[dn_face_extract];
  PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
  int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
  int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
  PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

  int n_elmt_current = 0;
  int n_edge_current = 0;
  tmp_dface_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(dn_face_extract,
                              &n_elmt_current,
                              &n_edge_current,
                              distrib_extract_face[i_rank],
                              -1,
                              dextract_face_vtx,
                              dextract_face_vtx_idx,
                              tmp_dface_edge_vtx_idx,
                              tmp_dface_edge_vtx,
                              tmp_dface_edge,
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_elt_tot);
 PDM_free(tmp_parent_elmt_pos);

  int dn_edge = 0;
  PDM_g_num_t  *edge_distrib   = NULL;
  int          *dedge_vtx_idx  = NULL;
  PDM_g_num_t  *dedge_vtx      = NULL;
  int          *dedge_face_idx = NULL;
  PDM_g_num_t  *dedge_face     = NULL;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        distrib_extract_vtx[n_rank],
                                        n_edge_elt_tot,
                                        tmp_dface_edge,
                                        tmp_dface_edge_vtx_idx,
                                        tmp_dface_edge_vtx,
                                        &dn_edge,
                                        &edge_distrib,
                                        &dedge_vtx_idx,
                                        &dedge_vtx,
                                        &dedge_face_idx,
                                        &dedge_face);

 PDM_free(dparent_face         );
 PDM_free(dextract_face_vtx_idx);
 PDM_free(dextract_face_vtx    );

  if(0 == 1) {
    PDM_log_trace_array_long(edge_distrib, n_rank+1               , "edge_distrib::");
    PDM_log_trace_array_long(dedge_vtx   , dedge_vtx_idx [dn_edge], "dedge_vtx::"   );
    PDM_log_trace_array_long(dedge_face  , dedge_face_idx[dn_edge], "dedge_face::"  );
    PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face  , dn_edge, "dedge_face::"  );
  }


  PDM_part_to_block_t* ptb = PDM_part_to_block_create_from_distrib(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                   PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                   1.,
                                                                   &dgroup_face,
                                                                   distrib_face,
                                                                   &dgroup_face_idx[n_group_face],
                                                                   1,
                                                                   comm);

  int *pface_group   = malloc(dgroup_face_idx[n_group_face] * sizeof(int));
  int *pface_group_n = malloc(dgroup_face_idx[n_group_face] * sizeof(int));
  for(int i_group = 0; i_group < n_group_face; ++i_group) {
    for(int idx_face = dgroup_face_idx[i_group]; idx_face < dgroup_face_idx[i_group+1]; ++idx_face) {
      pface_group  [idx_face] = (i_group+1);
      pface_group_n[idx_face] = 1;
    }
  }

  if(0 == 1) {
    PDM_log_trace_array_long(dgroup_face  , dgroup_face_idx[n_group_face], "dgroup_face   ::");
    PDM_log_trace_array_int (pface_group_n, dgroup_face_idx[n_group_face], "pface_group_n ::");
  }

  int *dface_group   = NULL;
  int *tmp_dface_group_n = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                        &pface_group_n,
             (void **)  &pface_group,
                        &tmp_dface_group_n,
             (void **)  &dface_group);

 PDM_free(pface_group_n);
 PDM_free(pface_group);

  int n_elt_face = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t* blk_gnum_face = PDM_part_to_block_block_gnum_get(ptb);

  int *dface_group_n = PDM_array_zeros_int(dn_face);

  for(int i1 = 0; i1 < n_elt_face; ++i1) {
    int i = (int) (blk_gnum_face[i1] - distrib_face[i_rank] - 1);
    dface_group_n[i] = tmp_dface_group_n[i1];
  }
 PDM_free(tmp_dface_group_n);

  PDM_part_to_block_free(ptb);

  /*
   * Go to extract_face frame
   */
  int **tmp_dextract_face_group_n = NULL;
  int **tmp_dextract_face_group   = NULL;
  PDM_block_to_part_exch(btp_face_to_extract_face,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         dface_group_n,
                         dface_group,
             (int ***)   &tmp_dextract_face_group_n,
             (void ***)  &tmp_dextract_face_group);
  int *dextract_face_group_n = tmp_dextract_face_group_n[0];
  int *dextract_face_group   = tmp_dextract_face_group  [0];
 PDM_free(tmp_dextract_face_group_n);
 PDM_free(tmp_dextract_face_group  );

 PDM_free(dface_group_n);
 PDM_free(dface_group);
 PDM_free(dextract_face_group_n);
  PDM_block_to_part_free(btp_face_to_extract_face);

  /* Go back to edge */
  int dn_edge_twice = dedge_face_idx[dn_edge];
  PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_extract_face,
                            (const PDM_g_num_t **)    &dedge_face,
                                                      &dn_edge_twice,
                                                      1,
                                                      comm);

  int **dedge_face_group_tmp = NULL;
  int stride_one = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dextract_face_group,
                         NULL,
             (void ***) &dedge_face_group_tmp);
  int *dedge_face_group = dedge_face_group_tmp[0];
 PDM_free(dedge_face_group_tmp);
 PDM_free(dextract_face_group);
  PDM_block_to_part_free(btp);

 PDM_free(distrib_extract_face );

  /*
   * Prepare gnum
   */
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             1,
                                             PDM_TRUE,
                                             1e-4,
                                             comm,
                                             PDM_OWNERSHIP_USER);


  int n_max_nuplet = 0;
  for(int i = 0; i < dn_edge; ++i) {
    n_max_nuplet = PDM_MAX(n_max_nuplet, dedge_face_idx[i+1] - dedge_face_idx[i]);
  }

  PDM_g_num_t *edge_group   = NULL;
  int         *pridge_edge  = malloc( dn_edge       * sizeof(int        ));

  /* For each ridge keep the link with the face group associated */
  int         *pridge_face_group_idx  = malloc( (dn_edge+1)            * sizeof(int        ));
  int         *pridge_face_group      = malloc( n_max_nuplet * dn_edge * sizeof(int        ));

  int          dn_ridge     = 0;

  PDM_gnum_set_parents_nuplet(gen_gnum, n_max_nuplet);

  PDM_g_num_t *edge_doublet = malloc( n_max_nuplet * dn_edge * sizeof(PDM_g_num_t));

  int *group_list = malloc(n_max_nuplet * sizeof(int));
  int idx_write = 0;
  pridge_face_group_idx[0] = 0;
  for(int i = 0; i < dn_edge; ++i) {

    int beg     = dedge_face_idx[i];
    int n_strid = dedge_face_idx[i+1] - beg;

    for(int k = 0; k < n_max_nuplet; ++k) {
      group_list[k] = 10000000;
    }

    for(int k = 0; k < n_strid; ++k) {
      group_list[k] = dedge_face_group[beg+k];
    }

    // PDM_sort_int(group_list, NULL, n_strid);
    int n_unique = PDM_inplace_unique(group_list, 0, n_strid-1);

    if(n_strid != 1 && n_unique == 1) {
      continue;
    }

    for(int k = 0; k < n_unique; ++k) {
      edge_doublet[n_max_nuplet*idx_write+k] = group_list[k];
    }
    for(int k = n_unique; k < n_max_nuplet; ++k) {
      edge_doublet[n_max_nuplet*idx_write+k] = -1;
    }


    pridge_face_group_idx[dn_ridge+1] = pridge_face_group_idx[dn_ridge];
    for(int k = 0; k < n_unique; ++k) {
      pridge_face_group[pridge_face_group_idx[dn_ridge+1]++] = group_list[k];
    }

    // Cas 1 : 1 neihbor -> A ridge
    if(n_strid == 1) {
      pridge_edge[dn_ridge++] = i;
    } else if(n_unique > 1) {
      pridge_edge[dn_ridge++] = i;
    }
    idx_write++;

    // Dans le cas particulier manifold et pas de truc tordu
    // int igroup1 = dedge_face_group[2*i  ];
    // int igroup2 = dedge_face_group[2*i+1];

    // if(igroup1 != igroup2) {
    //   edge_doublet[2*idx_write  ] = PDM_MIN(igroup1, igroup2);
    //   edge_doublet[2*idx_write+1] = PDM_MAX(igroup1, igroup2);
    //   pridge_edge[dn_ridge++] = i;
    //   idx_write++;
    // }

  }
 PDM_free(dedge_face_group);
 PDM_free(group_list);

  PDM_realloc(pridge_face_group_idx ,pridge_face_group_idx ,                    (dn_ridge+1) ,int);
  PDM_realloc(pridge_face_group     ,pridge_face_group     , pridge_face_group_idx[dn_ridge] ,int);

  PDM_gnum_set_from_parents(gen_gnum, 0, dn_ridge, edge_doublet);
  // PDM_gnum_set_parents_nuplet(gen_gnum, 2);

  PDM_gnum_compute(gen_gnum);
  edge_group = PDM_gnum_get(gen_gnum, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(edge_group, dn_ridge, "edge_group ::");
  }

  int _n_group_ridge = 0;
  for(int i = 0; i < dn_ridge; ++i) {
    _n_group_ridge = PDM_MAX(_n_group_ridge, edge_group[i]);
  }
  int n_group_ridge = 0;
  PDM_MPI_Allreduce(&_n_group_ridge, &n_group_ridge, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

 PDM_free(edge_doublet);
  /*
   * Hook edge
   */
  PDM_g_num_t *dridge_vtx = malloc(2 * dn_ridge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_ridge; ++i) {
    int i_edge = pridge_edge[i];
    dridge_vtx[2*i  ] = dedge_vtx[2*i_edge  ];
    dridge_vtx[2*i+1] = dedge_vtx[2*i_edge+1];
  }
 PDM_free(pridge_edge);
  // PDM_log_trace_array_long(dridge_vtx, 2 * dn_ridge, "dridge_vtx ::");

  /*
   * Re-création des groupes
   */
  int *dgroup_edge_n = malloc(n_group_ridge * sizeof(int));

  for(int i = 0; i < n_group_ridge; ++i) {
    dgroup_edge_n[i] = 0;
  }

  for(int i = 0; i < dn_ridge; ++i) {
    dgroup_edge_n[edge_group[i]-1]++;
  }

  if(0 == 1) {
    PDM_log_trace_array_int(dgroup_edge_n, n_group_ridge, "dgroup_edge_n ::");
  }

  int *dgroup_edge_idx = malloc((n_group_ridge+1) * sizeof(int));
  dgroup_edge_idx[0] = 0;
  for(int i = 0; i < n_group_ridge; ++i) {
    dgroup_edge_idx[i+1] = dgroup_edge_idx[i] + dgroup_edge_n[i];
    dgroup_edge_n[i] = 0;
  }

  PDM_g_num_t* dgroup_edge = malloc(dgroup_edge_idx[n_group_ridge] * sizeof(PDM_g_num_t));

  PDM_g_num_t* distrib_ridge = PDM_compute_entity_distribution(comm, dn_ridge);

  for(int i = 0; i < dn_ridge; ++i) {
    int i_group = edge_group[i]-1;
    int idx = dgroup_edge_idx[i_group] + dgroup_edge_n[i_group]++;
    dgroup_edge[idx] = distrib_ridge[i_rank] + i + 1;
  }

  if(0 == 1) {
    PDM_log_trace_array_int(dgroup_edge_idx, n_group_ridge+1, "dgroup_edge_idx ::");
    PDM_log_trace_connectivity_long(dgroup_edge_idx, dgroup_edge, n_group_ridge, "dgroup_edge :: ");
  }

 PDM_free(edge_group);
  PDM_gnum_free(gen_gnum);

  /*
   * Update ridge gnum
   */
  int dn_ridge_twice = 2 * dn_ridge;
  PDM_block_to_part_t* btp_vtx = PDM_block_to_part_create(distrib_extract_vtx,
                                   (const PDM_g_num_t **) &dridge_vtx,
                                                          &dn_ridge_twice,
                                                          1,
                                                          comm);

  PDM_g_num_t **tmp_dridge_vtx_parent = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         &stride_one,
                         dparent_vtx,
                         NULL,
             (void ***) &tmp_dridge_vtx_parent);
  PDM_g_num_t *dridge_vtx_parent = tmp_dridge_vtx_parent[0];
 PDM_free(tmp_dridge_vtx_parent);


 PDM_free(dparent_vtx);
 PDM_free(distrib_extract_vtx);
 PDM_free(dridge_vtx);

  PDM_block_to_part_free(btp_vtx);


  /* Fix output */
  *out_distrib_ridge   = distrib_ridge;
  *out_dridge_vtx      = dridge_vtx_parent;
  *out_n_group_ridge   = n_group_ridge;
  *out_dgroup_edge     = dgroup_edge;
  *out_dgroup_edge_idx = dgroup_edge_idx;

  *out_dridge_face_group_idx = pridge_face_group_idx;
  *out_dridge_face_group     = pridge_face_group;

  /*
   * Free
   */
 PDM_free(dgroup_edge_n );
 PDM_free(edge_distrib  );
 PDM_free(dedge_vtx_idx );
 PDM_free(dedge_vtx     );
 PDM_free(dedge_face_idx);
 PDM_free(dedge_face    );

}





#ifdef __cplusplus
}
#endif /* __cplusplus */
