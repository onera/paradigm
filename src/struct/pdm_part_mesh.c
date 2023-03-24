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
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_dmesh.h"

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
 * \param [in]   dn_cell             Number of distributed cells
 * \param [in]   dn_face             Number of distributed faces
 * \param [in]   dn_vtx              Number of distributed vertices
 * \param [in]   n_bnd               Number of boundaries
 * \param [in]   n_join              Number of interfaces with other zones
 *
 * \return     Identifier
 */

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
)
{
  PDM_part_mesh_t *pmesh = (PDM_part_mesh_t *) malloc(sizeof(PDM_part_mesh_t));

  pmesh->n_part = n_part;
  pmesh->comm   = comm;

  int tn_part;
  PDM_MPI_Allreduce(&pmesh->n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmesh->tn_part = tn_part;

  pmesh->pconnectivity          = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pconnectivity_idx      = (int         *** ) malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(int         **) );
  pmesh->pentity_ln_to_gn       = (PDM_g_num_t *** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_g_num_t **) );
  pmesh->pentity_color          = (int         *** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(int         **) );

  pmesh->pn_entity              = (int          ** ) malloc( PDM_MESH_ENTITY_MAX       * sizeof(int          *) );

  pmesh->is_owner_connectivity  = malloc( PDM_CONNECTIVITY_TYPE_MAX * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_ln_to_gn      = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );
  pmesh->is_owner_color         = malloc( PDM_MESH_ENTITY_MAX       * sizeof(PDM_bool_t   ) );

  for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
    pmesh->is_owner_connectivity[i] = PDM_FALSE;
    pmesh->pconnectivity        [i] = NULL;
    pmesh->pconnectivity_idx    [i] = NULL;
  }

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    pmesh->is_owner_ln_to_gn[i] = PDM_FALSE;
    pmesh->is_owner_color   [i] = PDM_FALSE;
    pmesh->pentity_ln_to_gn [i] = NULL;
    pmesh->pentity_color    [i] = NULL;
    pmesh->pn_entity        [i] = NULL;
  }

  pmesh->pn_bound                = malloc( PDM_BOUND_TYPE_MAX * sizeof(int          **) );
  pmesh->pbound                  = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         ***) );
  pmesh->pbound_ln_to_gn         = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t ***) );
  pmesh->is_owner_bound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t     ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i ) {
    pmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_bound         [i] = PDM_FALSE;

    pmesh->pn_bound       [i] = malloc( pmesh->n_part * sizeof(int          *) );
    pmesh->pbound         [i] = malloc( pmesh->n_part * sizeof(int         **) );
    pmesh->pbound_ln_to_gn[i] = malloc( pmesh->n_part * sizeof(PDM_g_num_t **) );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pmesh->pn_bound               [i][i_part] = NULL;
      pmesh->pbound                 [i][i_part] = NULL;
      pmesh->pbound_ln_to_gn        [i][i_part] = NULL;
    }
  }

  pmesh->pconcat_bound_idx       = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->pconcat_bound           = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->pconcat_bound_ln_to_gn  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_g_num_t **) );
  pmesh->is_owner_concat_bound   = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );
  pmesh->is_compute_concat_bound = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );
  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_concat_bound         [i] = PDM_FALSE;
    pmesh->is_compute_concat_bound       [i] = PDM_FALSE;

    pmesh->pconcat_bound_idx     [i] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->pconcat_bound         [i] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->pconcat_bound_ln_to_gn[i] = malloc( pmesh->n_part * sizeof(PDM_g_num_t *) );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pmesh->pconcat_bound_idx     [i][i_part] = NULL;
      pmesh->pconcat_bound         [i][i_part] = NULL;
      pmesh->pconcat_bound_ln_to_gn[i][i_part] = NULL;
    }
  }


  pmesh->vtx_coords = malloc(pmesh->n_part * sizeof(double));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pmesh->vtx_coords[i_part] = NULL;
  }
  pmesh->is_owner_vtx_coord = PDM_FALSE;

  pmesh->ppart_bound_proc_idx = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->ppart_bound_part_idx = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->ppart_bound          = malloc( PDM_BOUND_TYPE_MAX * sizeof(int         **) );
  pmesh->is_owner_part_bound  = malloc( PDM_BOUND_TYPE_MAX * sizeof(PDM_bool_t    ) );

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->is_owner_part_bound [i] = PDM_FALSE;
    pmesh->ppart_bound_proc_idx[i] = NULL;
    pmesh->ppart_bound_part_idx[i] = NULL;
    pmesh->ppart_bound         [i] = NULL;
  }

  return pmesh;
}


void
PDM_part_mesh_n_entity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                       pn_entity
)
{
  if(pmesh->pn_entity[entity_type] == NULL) {
    pmesh->pn_entity[entity_type] = malloc(pmesh->n_part * sizeof(int));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pn_entity[entity_type][i] = 0;
    }
  }
  pmesh->pn_entity[entity_type][i_part] = pn_entity;
}


int
PDM_part_mesh_n_entity_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type
)
{
  if(pmesh->pn_entity[entity_type] != NULL) {
    return pmesh->pn_entity[entity_type][i_part];
  } else {
    return 0;
  }

}


void
PDM_part_mesh_entity_ln_to_gn_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t              *pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn [entity_type] == NULL) {
    pmesh->pentity_ln_to_gn [entity_type] = malloc(pmesh->n_part * sizeof(PDM_g_num_t *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_ln_to_gn [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_ln_to_gn [entity_type][i_part] = pentity_ln_to_gn;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_entity_ln_to_gn_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 PDM_g_num_t             **pentity_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_ln_to_gn[entity_type] != NULL) {
    *pentity_ln_to_gn = pmesh->pentity_ln_to_gn[entity_type][i_part];
  } else {
    *pentity_ln_to_gn = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_entity_color_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                      *pentity_color,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_color [entity_type] == NULL) {
    pmesh->pentity_color [entity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_color [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_color [entity_type][i_part] = pentity_color;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_color[entity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_color[entity_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_entity_color_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_mesh_entities_t       entity_type,
 int                     **pentity_color,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pentity_color[entity_type] != NULL) {
    *pentity_color = pmesh->pentity_color[entity_type][i_part];
  } else {
    *pentity_color = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_color[entity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_color[entity_type] = PDM_TRUE;
  }
}



void
PDM_part_mesh_connectivity_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_connectivity_type_t   connectivity_type,
 int                      *connect,
 int                      *connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity [connectivity_type] == NULL) {
    pmesh->pconnectivity    [connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    pmesh->pconnectivity_idx[connectivity_type] = malloc(pmesh->n_part * sizeof(int *));
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pconnectivity    [connectivity_type][i] = NULL;
      pmesh->pconnectivity_idx[connectivity_type][i] = NULL;
    }
  }
  pmesh->pconnectivity        [connectivity_type][i_part] = connect;
  pmesh->pconnectivity_idx    [connectivity_type][i_part] = connect_idx;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_set
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double            *vtx_coord,
 PDM_ownership_t    ownership
)
{
  pmesh->vtx_coords[i_part] = vtx_coord;
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_vtx_coord_get
(
 PDM_part_mesh_t   *pmesh,
 int                i_part,
 double           **vtx_coord,
 PDM_ownership_t    ownership
)
{
  if(pmesh->vtx_coords[i_part] != NULL) {
    *vtx_coord = pmesh->vtx_coords[i_part];
  } else {
    *vtx_coord = NULL;
  }
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_vtx_coord = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_vtx_coord = PDM_TRUE;
  }
}

void
PDM_part_mesh_connectivity_get
(
 PDM_part_mesh_t           *pmesh,
 int                        i_part,
 PDM_connectivity_type_t    connectivity_type,
 int                      **connect,
 int                      **connect_idx,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pconnectivity    [connectivity_type] != NULL && pmesh->pconnectivity    [connectivity_type][i_part] != NULL) {
    *connect     = pmesh->pconnectivity    [connectivity_type][i_part];
    *connect_idx = pmesh->pconnectivity_idx[connectivity_type][i_part];
  } else {
    *connect     = NULL;
    *connect_idx = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_connectivity[connectivity_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_n_bound_set
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type,
 int                       n_bound
)
{
  assert(pmesh->n_group_bnd[bound_type] == 0);
  pmesh->n_group_bnd[bound_type] = n_bound;

  for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
    pmesh->pn_bound       [bound_type][i_part] = malloc( n_bound * sizeof(int          ) );
    pmesh->pbound         [bound_type][i_part] = malloc( n_bound * sizeof(int         *) );
    pmesh->pbound_ln_to_gn[bound_type][i_part] = malloc( n_bound * sizeof(PDM_g_num_t *) );
    for(int i_group = 0; i_group < n_bound; ++i_group) {
      pmesh->pn_bound       [bound_type][i_part][i_group] = 0;
      pmesh->pbound         [bound_type][i_part][i_group] = NULL;
      pmesh->pbound_ln_to_gn[bound_type][i_part][i_group] = NULL;
    }
  }
}

int
PDM_part_mesh_n_bound_get
(
 PDM_part_mesh_t          *pmesh,
 PDM_bound_type_t          bound_type
)
{
  return pmesh->n_group_bnd[bound_type];
}

int
PDM_part_mesh_tn_part_get
(
 PDM_part_mesh_t          *pmesh
)
{
  return pmesh->tn_part;
}

void
PDM_part_mesh_bound_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                       pn_bound,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  pmesh->pn_bound       [bound_type][i_part][i_group] = pn_bound;
  pmesh->pbound         [bound_type][i_part][i_group] = pbound;
  pmesh->pbound_ln_to_gn[bound_type][i_part][i_group] = pbound_ln_to_gn;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 int                       i_group,
 PDM_bound_type_t          bound_type,
 int                      *pn_bound,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  if(pmesh->pbound[bound_type] != NULL) {
    *pn_bound        = pmesh->pn_bound       [bound_type][i_part][i_group];
    *pbound          = pmesh->pbound         [bound_type][i_part][i_group];
    *pbound_ln_to_gn = pmesh->pbound_ln_to_gn[bound_type][i_part][i_group];
  } else {
    *pn_bound        = 0;
    *pbound          = NULL;
    *pbound_ln_to_gn = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_bound[bound_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_bound_concat_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                       n_bound,
 int                      *pbound_idx,
 int                      *pbound,
 PDM_g_num_t              *pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{
  pmesh->n_group_bnd           [bound_type]         = n_bound;
  pmesh->pconcat_bound_idx     [bound_type][i_part] = pbound_idx;
  pmesh->pconcat_bound         [bound_type][i_part] = pbound;
  pmesh->pconcat_bound_ln_to_gn[bound_type][i_part] = pbound_ln_to_gn;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_FALSE;
    pmesh->is_compute_concat_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_TRUE;
    pmesh->is_compute_concat_bound[bound_type] = PDM_TRUE;
  }

}


void
PDM_part_mesh_bound_concat_compute
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 PDM_ownership_t           ownership
)
{
  if(pmesh->is_compute_concat_bound[bound_type] == PDM_FALSE) {
    int n_group = pmesh->n_group_bnd[bound_type];
    assert(pmesh->pconcat_bound_idx[bound_type][i_part] == NULL);
    pmesh->pconcat_bound_idx[bound_type][i_part] = malloc( (n_group+1) * sizeof(int));
    int *_pconcat_bound_idx = pmesh->pconcat_bound_idx[bound_type][i_part];

    _pconcat_bound_idx[0] = 0;
    for(int i_group = 0; i_group < n_group; ++i_group) {
      _pconcat_bound_idx[i_group+1] = _pconcat_bound_idx[i_group] + pmesh->pn_bound[bound_type][i_part][i_group];
    }

    pmesh->pconcat_bound         [bound_type][i_part] = malloc(_pconcat_bound_idx[n_group] * sizeof(int        ));
    pmesh->pconcat_bound_ln_to_gn[bound_type][i_part] = malloc(_pconcat_bound_idx[n_group] * sizeof(PDM_g_num_t));
    int         *_pconcat_bound          = pmesh->pconcat_bound         [bound_type][i_part];
    PDM_g_num_t *_pconcat_bound_ln_to_gn = pmesh->pconcat_bound_ln_to_gn[bound_type][i_part];

    for(int i_group = 0; i_group < n_group; ++i_group) {
      int idx_write = _pconcat_bound_idx[i_group];
      for(int i = 0; i < pmesh->pn_bound[bound_type][i_part][i_group]; ++i) {
        _pconcat_bound         [idx_write+i] = pmesh->pbound         [bound_type][i_part][i_group][i];
        _pconcat_bound_ln_to_gn[idx_write+i] = pmesh->pbound_ln_to_gn[bound_type][i_part][i_group][i];
      }
    }
    pmesh->is_compute_concat_bound[bound_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_bound_concat_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **pbound_idx,
 int                     **pbound,
 PDM_g_num_t             **pbound_ln_to_gn,
 PDM_ownership_t           ownership
)
{

  *pbound_idx      = pmesh->pconcat_bound_idx     [bound_type][i_part];
  *pbound          = pmesh->pconcat_bound         [bound_type][i_part];
  *pbound_ln_to_gn = pmesh->pconcat_bound_ln_to_gn[bound_type][i_part];

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_concat_bound  [bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_part_graph_comm_set
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                      *ppart_bound_proc_idx,
 int                      *ppart_bound_part_idx,
 int                      *ppart_bound,
 PDM_ownership_t           ownership
)
{
  if(pmesh->ppart_bound[bound_type] == NULL) {

    pmesh->ppart_bound_proc_idx[bound_type] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->ppart_bound_part_idx[bound_type] = malloc( pmesh->n_part * sizeof(int         *) );
    pmesh->ppart_bound         [bound_type] = malloc( pmesh->n_part * sizeof(PDM_g_num_t *) );
    for(int j = 0; j < pmesh->n_part; ++j) {
      pmesh->ppart_bound_proc_idx[bound_type][j] = NULL;
      pmesh->ppart_bound_part_idx[bound_type][j] = NULL;
      pmesh->ppart_bound         [bound_type][j] = NULL;
    }
  }

  pmesh->ppart_bound_proc_idx[bound_type][i_part] = ppart_bound_proc_idx;
  pmesh->ppart_bound_part_idx[bound_type][i_part] = ppart_bound_part_idx;
  pmesh->ppart_bound         [bound_type][i_part] = ppart_bound;

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_part_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_part_bound[bound_type] = PDM_TRUE;
  }
}


void
PDM_part_mesh_part_graph_comm_get
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type,
 int                     **ppart_bound_proc_idx,
 int                     **ppart_bound_part_idx,
 int                     **ppart_bound,
 PDM_ownership_t           ownership
)
{
  if(pmesh->ppart_bound[bound_type] != NULL) {
    *ppart_bound_proc_idx = pmesh->ppart_bound_proc_idx[bound_type][i_part];
    *ppart_bound_part_idx = pmesh->ppart_bound_part_idx[bound_type][i_part];
    *ppart_bound          = pmesh->ppart_bound         [bound_type][i_part];
  } else {
    *ppart_bound_proc_idx = NULL;
    *ppart_bound_part_idx = NULL;
    *ppart_bound          = NULL;
  }

  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    pmesh->is_owner_part_bound[bound_type] = PDM_FALSE;
  } else {
    pmesh->is_owner_part_bound[bound_type] = PDM_TRUE;
  }
}

void
PDM_part_mesh_free
(
 PDM_part_mesh_t        *pmesh
)
{
  if (pmesh != NULL) {
    /* Free connectivity */
    for(int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; ++i) {
      if(pmesh->is_owner_connectivity[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pconnectivity[i][i_part] != NULL) {
            free(pmesh->pconnectivity[i][i_part]);
          }
          if(pmesh->pconnectivity_idx[i][i_part] != NULL) {
            free(pmesh->pconnectivity_idx[i][i_part]);
          }
        }
      }

      if(pmesh->pconnectivity[i] != NULL) {
        free(pmesh->pconnectivity[i]);
        pmesh->pconnectivity[i] = NULL;
      }

      if(pmesh->pconnectivity_idx[i] != NULL) {
        free(pmesh->pconnectivity_idx[i]);
        pmesh->pconnectivity_idx[i] = NULL;
      }
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_ln_to_gn[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pentity_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pentity_ln_to_gn[i] != NULL) {
        free(pmesh->pentity_ln_to_gn[i]);
        pmesh->pentity_ln_to_gn[i] = NULL;
      }
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_color[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_color[i][i_part] != NULL) {
            free(pmesh->pentity_color[i][i_part]);
          }
        }
      }

      if(pmesh->pentity_color[i] != NULL) {
        free(pmesh->pentity_color[i]);
        pmesh->pentity_color[i] = NULL;
      }
    }

    /* Free vtx__coord */
    if(pmesh->is_owner_vtx_coord == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        free(pmesh->vtx_coords[i_part]);
        pmesh->vtx_coords[i_part] = NULL;
      }
    }
    free(pmesh->vtx_coords);

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          for(int i_group = 0; i_group < pmesh->n_group_bnd[i]; ++i_group) {
            if(pmesh->pbound[i][i_part][i_group] != NULL) {
              free(pmesh->pbound[i][i_part][i_group]);
            }
            if(pmesh->pbound_ln_to_gn[i][i_part][i_group] != NULL) {
              free(pmesh->pbound_ln_to_gn[i][i_part][i_group]);
            }
          }

          if(pmesh->pn_bound[i][i_part] != NULL) {
            free(pmesh->pn_bound[i][i_part]);
          }
          if(pmesh->pbound[i][i_part] != NULL) {
            free(pmesh->pbound[i][i_part]);
          }
          if(pmesh->pbound_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pbound_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pn_bound[i] != NULL) {
        free(pmesh->pn_bound[i]);
        pmesh->pn_bound[i] = NULL;
      }

      if(pmesh->pbound[i] != NULL) {
        free(pmesh->pbound[i]);
        pmesh->pbound[i] = NULL;
      }

      if(pmesh->pbound_ln_to_gn[i] != NULL) {
        free(pmesh->pbound_ln_to_gn[i]);
        pmesh->pbound_ln_to_gn[i] = NULL;
      }
    }


    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_concat_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pconcat_bound_idx[i][i_part] != NULL) {
            free(pmesh->pconcat_bound_idx[i][i_part]);
          }
          if(pmesh->pconcat_bound[i][i_part] != NULL) {
            free(pmesh->pconcat_bound[i][i_part]);
          }
          if(pmesh->pconcat_bound_ln_to_gn[i][i_part] != NULL) {
            free(pmesh->pconcat_bound_ln_to_gn[i][i_part]);
          }
        }
      }

      if(pmesh->pconcat_bound_idx[i] != NULL) {
        free(pmesh->pconcat_bound_idx[i]);
        pmesh->pconcat_bound_idx[i] = NULL;
      }

      if(pmesh->pconcat_bound[i] != NULL) {
        free(pmesh->pconcat_bound[i]);
        pmesh->pconcat_bound[i] = NULL;
      }

      if(pmesh->pconcat_bound_ln_to_gn[i] != NULL) {
        free(pmesh->pconcat_bound_ln_to_gn[i]);
        pmesh->pconcat_bound_ln_to_gn[i] = NULL;
      }
    }

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_part_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->ppart_bound_proc_idx[i][i_part] != NULL) {
            free(pmesh->ppart_bound_proc_idx[i][i_part]);
          }
          if(pmesh->ppart_bound_part_idx[i][i_part] != NULL) {
            free(pmesh->ppart_bound_part_idx[i][i_part]);
          }
          if(pmesh->ppart_bound[i][i_part] != NULL) {
            free(pmesh->ppart_bound[i][i_part]);
          }
        }

        if(pmesh->ppart_bound_proc_idx[i] != NULL) {
          free(pmesh->ppart_bound_proc_idx[i]);
          pmesh->ppart_bound_proc_idx[i] = NULL;
        }

        if(pmesh->ppart_bound_part_idx[i] != NULL) {
          free(pmesh->ppart_bound_part_idx[i]);
          pmesh->ppart_bound_part_idx[i] = NULL;
        }

        if(pmesh->ppart_bound[i] != NULL) {
          free(pmesh->ppart_bound[i]);
          pmesh->ppart_bound[i] = NULL;
        }
      }
    }


    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->pn_entity[i] !=NULL){
        free(pmesh->pn_entity[i]);
      }
    }

    free(pmesh->pn_entity);
    free(pmesh->pconnectivity);
    free(pmesh->pconnectivity_idx);
    free(pmesh->pentity_ln_to_gn);
    free(pmesh->pentity_color);
    free(pmesh->is_owner_connectivity);
    free(pmesh->is_owner_ln_to_gn    );
    free(pmesh->is_owner_color       );

    free(pmesh->pn_bound       );
    free(pmesh->pbound         );
    free(pmesh->pbound_ln_to_gn);
    free(pmesh->is_owner_bound );

    free(pmesh->pconcat_bound         );
    free(pmesh->pconcat_bound_ln_to_gn);
    free(pmesh->pconcat_bound_idx);
    free(pmesh->is_owner_concat_bound);
    free(pmesh->is_compute_concat_bound);

    free(pmesh->ppart_bound_proc_idx);
    free(pmesh->ppart_bound_part_idx);
    free(pmesh->ppart_bound         );
    free(pmesh->is_owner_part_bound );

    free(pmesh);
    pmesh = NULL;
  }
}

