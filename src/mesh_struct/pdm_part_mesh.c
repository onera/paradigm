/*============================================================================
 * Interface structure to represent a distributed mesh
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_extract_part.h"
#include "pdm_gnum.h"
#include "pdm_io.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_part_mesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_priv.h"
#include "pdm_writer.h"


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
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

/* Local bound extractions for \ref PDM_part_mesh_dump_ensight */

static void
_build_extract_part_bound
(
 PDM_part_mesh_t     *pmesh,
 PDM_bound_type_t     bound_type,
 PDM_extract_part_t **extrp
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(pmesh->comm, &i_rank);

  int dim = 0;
  PDM_mesh_entities_t entity_type;
  switch (bound_type) {
    case PDM_BOUND_TYPE_VTX:
      dim = 0;
      entity_type = PDM_MESH_ENTITY_VTX;
      break;

    case PDM_BOUND_TYPE_EDGE:
      dim = 1;
      entity_type = PDM_MESH_ENTITY_EDGE;
      break;

    case PDM_BOUND_TYPE_FACE:
      dim = 2;
      entity_type = PDM_MESH_ENTITY_FACE;
      break;

    case PDM_BOUND_TYPE_CELL:
      dim = 3;
      entity_type = PDM_MESH_ENTITY_CELL;
      break;

    default:
      PDM_error("Invalid bound_type %d", bound_type);
  }

  PDM_UNUSED(entity_type);

  int   n_group          = pmesh->n_group_bnd[bound_type];
  int **group_entity_idx = pmesh->pconcat_bound_idx[bound_type];
  int **group_entity     = pmesh->pconcat_bound    [bound_type];

  *extrp = PDM_extract_part_create(dim,
                                   pmesh->n_part,
                                   pmesh->n_part,
                                   PDM_EXTRACT_PART_KIND_LOCAL,
                                   PDM_SPLIT_DUAL_WITH_HILBERT, // unused
                                   PDM_TRUE,
                                   PDM_OWNERSHIP_KEEP,
                                   pmesh->comm);

  /* Set parts */
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {

    PDM_extract_part_selected_lnum_set(*extrp,
                                       i_part,
                                       group_entity_idx[i_part][n_group],
                                       group_entity[i_part],
                                       PDM_OWNERSHIP_USER);

    int          n_cell        = 0;
    int          n_face        = 0;
    int          n_edge        = 0;
    int         *cell_face_idx = NULL;
    int         *cell_face     = NULL;
    int         *face_edge_idx = NULL;
    int         *face_edge     = NULL;
    int         *face_vtx_idx  = NULL;
    int         *face_vtx      = NULL;
    int         *edge_vtx      = NULL;
    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;

    if (dim >= 3 && pmesh->pn_entity[PDM_MESH_ENTITY_CELL] != NULL) {
      n_cell        = pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
      cell_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE] != NULL) {
        cell_face_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part];
        cell_face     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part];
      }
    }

    if (pmesh->pn_entity[PDM_MESH_ENTITY_FACE] != NULL) {
      n_face        = pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part];
      face_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE] != NULL) {
        face_edge_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
        face_edge     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
      }

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX] != NULL) {
        face_vtx_idx = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
        face_vtx     = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
      }
    }

    if (dim >= 1 && pmesh->pn_entity[PDM_MESH_ENTITY_EDGE] != NULL) {
      n_edge        = pmesh->pn_entity[PDM_MESH_ENTITY_EDGE][i_part];
      edge_ln_to_gn = pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_EDGE][i_part];

      if (pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] != NULL) {
        edge_vtx = pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part];
      }
    }

    PDM_extract_part_part_set(*extrp,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              pmesh->pn_entity[PDM_MESH_ENTITY_VTX][i_part],
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part],
                              pmesh->vtx_coords[i_part]);

  }

  PDM_extract_part_compute(*extrp);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_part_mesh_t*
PDM_part_mesh_create
(
 const int             n_part,
       PDM_MPI_Comm    comm
)
{
  PDM_part_mesh_t *pmesh;
  PDM_malloc(pmesh, 1, PDM_part_mesh_t);
  memset(pmesh, 0, sizeof(PDM_part_mesh_t));

  pmesh->n_part = n_part;
  pmesh->comm   = comm;

  int tn_part;
  PDM_MPI_Allreduce(&pmesh->n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmesh->tn_part = tn_part;

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

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->n_group_bnd[i] = 0;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    pmesh->pn_bound       [i] = NULL;
    pmesh->pbound         [i] = NULL;
    pmesh->pbound_ln_to_gn[i] = NULL;
    pmesh->is_owner_bound [i] = NULL;
  }

  for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
    PDM_malloc(pmesh->pconcat_bound_idx      [i], pmesh->n_part, int         *);
    PDM_malloc(pmesh->pconcat_bound          [i], pmesh->n_part, int         *);
    PDM_malloc(pmesh->pconcat_bound_ln_to_gn [i], pmesh->n_part, PDM_g_num_t *);
    PDM_malloc(pmesh->is_compute_concat_bound[i], pmesh->n_part, PDM_bool_t   );
    PDM_malloc(pmesh->is_owner_concat_bound  [i], pmesh->n_part, PDM_bool_t   );
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pmesh->pconcat_bound_idx      [i][i_part] = NULL;
      pmesh->pconcat_bound          [i][i_part] = NULL;
      pmesh->pconcat_bound_ln_to_gn [i][i_part] = NULL;
      pmesh->is_compute_concat_bound[i][i_part] = PDM_FALSE;
      pmesh->is_owner_concat_bound  [i][i_part] = PDM_FALSE;
    }
  }

  PDM_malloc(pmesh->vtx_coords, pmesh->n_part, double *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pmesh->vtx_coords[i_part] = NULL;
  }
  pmesh->is_owner_vtx_coord = PDM_FALSE;

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
    PDM_malloc(pmesh->pn_entity[entity_type], pmesh->n_part, int);
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
  if(pmesh->pentity_ln_to_gn[entity_type] == NULL) {
    PDM_malloc(pmesh->pentity_ln_to_gn[entity_type], pmesh->n_part, PDM_g_num_t *);
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_ln_to_gn [entity_type][i] = NULL;
    }
  }
  pmesh->pentity_ln_to_gn[entity_type][i_part] = pentity_ln_to_gn;
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
  if(pmesh->pentity_color[entity_type] == NULL) {
    PDM_malloc(pmesh->pentity_color [entity_type], pmesh->n_part, int *);
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pentity_color[entity_type][i] = NULL;
    }
  }
  pmesh->pentity_color[entity_type][i_part] = pentity_color;
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
    PDM_malloc(pmesh->pconnectivity    [connectivity_type], pmesh->n_part, int *);
    PDM_malloc(pmesh->pconnectivity_idx[connectivity_type], pmesh->n_part, int *);
    for(int i = 0; i < pmesh->n_part; ++i) {
      pmesh->pconnectivity    [connectivity_type][i] = NULL;
      pmesh->pconnectivity_idx[connectivity_type][i] = NULL;
    }
  }
  pmesh->pconnectivity    [connectivity_type][i_part] = connect;
  pmesh->pconnectivity_idx[connectivity_type][i_part] = connect_idx;

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
  if(pmesh->pconnectivity[connectivity_type] != NULL && pmesh->pconnectivity[connectivity_type][i_part] != NULL) {
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

  PDM_malloc(pmesh->pn_bound       [bound_type], pmesh->n_part, int          *);
  PDM_malloc(pmesh->pbound         [bound_type], pmesh->n_part, int         **);
  PDM_malloc(pmesh->pbound_ln_to_gn[bound_type], pmesh->n_part, PDM_g_num_t **);
  PDM_malloc(pmesh->is_owner_bound [bound_type], pmesh->n_part, PDM_bool_t   *);

  for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {

    PDM_malloc(pmesh->pn_bound       [bound_type][i_part], n_bound, int          );
    PDM_malloc(pmesh->pbound         [bound_type][i_part], n_bound, int         *);
    PDM_malloc(pmesh->pbound_ln_to_gn[bound_type][i_part], n_bound, PDM_g_num_t *);
    PDM_malloc(pmesh->is_owner_bound [bound_type][i_part], n_bound, PDM_bool_t   );

    for (int i_group = 0; i_group < n_bound; i_group++) {
      pmesh->pn_bound       [bound_type][i_part][i_group] = 0;
      pmesh->pbound         [bound_type][i_part][i_group] = NULL;
      pmesh->pbound_ln_to_gn[bound_type][i_part][i_group] = NULL;
      pmesh->is_owner_bound [bound_type][i_part][i_group] = PDM_FALSE;
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

int
PDM_part_mesh_n_part_get
(
 PDM_part_mesh_t          *pmesh
)
{
  return pmesh->n_part;
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
    pmesh->is_owner_bound[bound_type][i_part][i_group] = PDM_FALSE;
  } else {
    pmesh->is_owner_bound[bound_type][i_part][i_group] = PDM_TRUE;
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
    pmesh->is_owner_bound[bound_type][i_part][i_group] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_bound[bound_type][i_part][i_group] = PDM_TRUE;
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
    pmesh->is_owner_concat_bound  [bound_type][i_part] = PDM_FALSE;
    pmesh->is_compute_concat_bound[bound_type][i_part] = PDM_FALSE;
  } else {
    pmesh->is_owner_concat_bound  [bound_type][i_part] = PDM_TRUE;
    pmesh->is_compute_concat_bound[bound_type][i_part] = PDM_TRUE;
  }

}


void
PDM_part_mesh_bound_concat_compute
(
 PDM_part_mesh_t          *pmesh,
 int                       i_part,
 PDM_bound_type_t          bound_type
)
{
  if(pmesh->is_compute_concat_bound[bound_type][i_part] == PDM_FALSE) {
    int n_group = pmesh->n_group_bnd[bound_type];
    assert(pmesh->pconcat_bound_idx[bound_type][i_part] == NULL);
    PDM_malloc(pmesh->pconcat_bound_idx[bound_type][i_part], n_group+1, int);
    int *_pconcat_bound_idx = pmesh->pconcat_bound_idx[bound_type][i_part];

    _pconcat_bound_idx[0] = 0;
    for(int i_group = 0; i_group < n_group; ++i_group) {
      _pconcat_bound_idx[i_group+1] = _pconcat_bound_idx[i_group] + pmesh->pn_bound[bound_type][i_part][i_group];
    }

    PDM_malloc(pmesh->pconcat_bound         [bound_type][i_part], _pconcat_bound_idx[n_group], int        );
    PDM_malloc(pmesh->pconcat_bound_ln_to_gn[bound_type][i_part], _pconcat_bound_idx[n_group], PDM_g_num_t);
    int         *_pconcat_bound          = pmesh->pconcat_bound         [bound_type][i_part];
    PDM_g_num_t *_pconcat_bound_ln_to_gn = pmesh->pconcat_bound_ln_to_gn[bound_type][i_part];

    for(int i_group = 0; i_group < n_group; ++i_group) {
      int idx_write = _pconcat_bound_idx[i_group];
      for(int i = 0; i < pmesh->pn_bound[bound_type][i_part][i_group]; ++i) {
        _pconcat_bound         [idx_write+i] = pmesh->pbound         [bound_type][i_part][i_group][i];
        _pconcat_bound_ln_to_gn[idx_write+i] = pmesh->pbound_ln_to_gn[bound_type][i_part][i_group][i];
      }
    }
    pmesh->is_compute_concat_bound[bound_type][i_part] = PDM_TRUE;
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
    pmesh->is_owner_concat_bound  [bound_type][i_part] = PDM_FALSE;
  } else if (ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->is_owner_concat_bound  [bound_type][i_part] = PDM_TRUE;
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

    PDM_malloc(pmesh->ppart_bound_proc_idx[bound_type], pmesh->n_part, int *);
    PDM_malloc(pmesh->ppart_bound_part_idx[bound_type], pmesh->n_part, int *);
    PDM_malloc(pmesh->ppart_bound[bound_type], pmesh->n_part, int *);
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
PDM_part_mesh_part_comm_graph_set
(
  PDM_part_mesh_t       *pmesh,
  PDM_part_comm_graph_t *pcg,
  PDM_mesh_entities_t    entity_type,
  PDM_ownership_t        ownership
)
{
  pmesh->pcg[entity_type] = pcg;
  if (ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_KEEP) {
    pmesh->pcg_ownership[entity_type] = ownership;
  } else {
    PDM_error("Invalid ownership (got %d, must be %d or %d)", ownership, PDM_OWNERSHIP_KEEP, PDM_OWNERSHIP_USER);
  }
}

void
PDM_part_mesh_part_comm_graph_get
(
  PDM_part_mesh_t        *pmesh,
  PDM_mesh_entities_t     entity_type,
  PDM_part_comm_graph_t **pcg,
  PDM_ownership_t         ownership
)
{
  *pcg = pmesh->pcg[entity_type];
  if (ownership != PDM_OWNERSHIP_BAD_VALUE) {
    pmesh->pcg_ownership[entity_type] = ownership;
  }
}


void
PDM_part_mesh_part_comm_graph_compute_from_gnum
(
  PDM_part_mesh_t      *pmesh,
  PDM_mesh_entities_t   entity_type
)
{
  if (pmesh->pcg[entity_type] != NULL) {
    PDM_error("Invalid pointer : pmesh->pcg[entity_type=%d]!=NULL", entity_type);
  }

  int n_rank;
  PDM_MPI_Comm_size(pmesh->comm, &n_rank);

  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(pmesh->comm, pmesh->n_part);

  int  *n_entity_part_bound        = NULL;
  int **entity_proc_bound_idx      = NULL;
  int **entity_part_bound_idx      = NULL;
  int **entity_part_bound          = NULL;
  int **entity_part_bound_priority = NULL;
  PDM_part_generate_entity_graph_comm(pmesh->comm,
                                      part_distribution,
                                      NULL,
                                      pmesh->n_part,
                                      pmesh->pn_entity       [entity_type],
               (const PDM_g_num_t **) pmesh->pentity_ln_to_gn[entity_type],
                                      NULL,
                                      &entity_proc_bound_idx,
                                      &entity_part_bound_idx,
                                      &entity_part_bound,
                                      &entity_part_bound_priority);

  PDM_malloc(n_entity_part_bound, pmesh->n_part, int);
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
    n_entity_part_bound[i_part] = entity_part_bound_idx[i_part][part_distribution[n_rank]];

    PDM_free(entity_proc_bound_idx     [i_part]);
    PDM_free(entity_part_bound_idx     [i_part]);
    PDM_free(entity_part_bound_priority[i_part]);
  }
  PDM_free(entity_proc_bound_idx);
  PDM_free(entity_part_bound_idx);
  PDM_free(entity_part_bound_priority);
  PDM_free(part_distribution);

  // Build part comm graph
  pmesh->pcg[entity_type] = PDM_part_comm_graph_create(pmesh->n_part,
                                                       n_entity_part_bound,
                                                       entity_part_bound,
                                                       PDM_OWNERSHIP_KEEP,
                                                       pmesh->comm);
  pmesh->pcg_ownership[entity_type] = PDM_OWNERSHIP_KEEP;

  PDM_free(n_entity_part_bound);
  PDM_free(entity_part_bound);

}


void
PDM_part_mesh_gnum_compute_from_part_comm_graph
(
  PDM_part_mesh_t      *pmesh,
  PDM_mesh_entities_t   entity_type
)
{
  if (pmesh->pcg[entity_type] == NULL) {
    PDM_error("Invalid pmesh->pcg[entity_type=%d] == NULL but it's required to compute", entity_type);
  }

  if(pmesh->is_owner_ln_to_gn[entity_type] == PDM_TRUE) {
    for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
      if(pmesh->pentity_ln_to_gn[entity_type] != NULL) {
        PDM_free(pmesh->pentity_ln_to_gn[entity_type][i_part]);
      }
    }
  }
  PDM_free(pmesh->pentity_ln_to_gn[entity_type]);

  PDM_malloc(pmesh->pentity_ln_to_gn[entity_type], pmesh->n_part, PDM_g_num_t *);

  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3,
                                             pmesh->n_part,
                                             PDM_TRUE,
                                             1.e-6,
                                             pmesh->comm,
                                             PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_part_comm_graph(gen_gnum,
                                    pmesh->pn_entity[entity_type],
                                    pmesh->pcg[entity_type]);
  PDM_gnum_compute(gen_gnum);
  for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
    pmesh->pentity_ln_to_gn[entity_type][i_part] = PDM_gnum_get(gen_gnum, i_part);
  }
  pmesh->is_owner_ln_to_gn[entity_type] = PDM_TRUE;

  PDM_gnum_free(gen_gnum);
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
          if(pmesh->pconnectivity[i] != NULL) {
            PDM_free(pmesh->pconnectivity[i][i_part]);
          }
          if(pmesh->pconnectivity_idx[i] != NULL) {
            PDM_free(pmesh->pconnectivity_idx[i][i_part]);
          }
        }
      }
      PDM_free(pmesh->pconnectivity    [i]);
      PDM_free(pmesh->pconnectivity_idx[i]);
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_ln_to_gn[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_ln_to_gn[i] != NULL) {
            PDM_free(pmesh->pentity_ln_to_gn[i][i_part]);
          }
        }
      }
      PDM_free(pmesh->pentity_ln_to_gn[i]);
    }

    /* Free ln_to_gn */
    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->is_owner_color[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->pentity_color[i] != NULL) {
            PDM_free(pmesh->pentity_color[i][i_part]);
          }
        }
      }
      PDM_free(pmesh->pentity_color[i]);
    }

    /* Free vtx__coord */
    if(pmesh->is_owner_vtx_coord == PDM_TRUE) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        PDM_free(pmesh->vtx_coords[i_part]);
      }
    }
    PDM_free(pmesh->vtx_coords);

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if (pmesh->pn_bound[i] == NULL) {
        continue;
      }
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        for(int i_group = 0; i_group < pmesh->n_group_bnd[i]; ++i_group) {
          if(pmesh->is_owner_bound[i][i_part][i_group] == PDM_TRUE) {
            PDM_free(pmesh->pbound         [i][i_part][i_group]);
            PDM_free(pmesh->pbound_ln_to_gn[i][i_part][i_group]);
          }
        }
        PDM_free(pmesh->pn_bound       [i][i_part]);
        PDM_free(pmesh->pbound         [i][i_part]);
        PDM_free(pmesh->pbound_ln_to_gn[i][i_part]);
        PDM_free(pmesh->is_owner_bound [i][i_part]);
      }
      PDM_free(pmesh->pn_bound       [i]);
      PDM_free(pmesh->pbound         [i]);
      PDM_free(pmesh->pbound_ln_to_gn[i]);
      PDM_free(pmesh->is_owner_bound [i]);
    }

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
        if(pmesh->is_owner_concat_bound[i][i_part] == PDM_TRUE) {
          PDM_free(pmesh->pconcat_bound_idx     [i][i_part]);
          PDM_free(pmesh->pconcat_bound         [i][i_part]);
          PDM_free(pmesh->pconcat_bound_ln_to_gn[i][i_part]);
        }
      }
      PDM_free(pmesh->pconcat_bound_idx      [i]);
      PDM_free(pmesh->pconcat_bound          [i]);
      PDM_free(pmesh->pconcat_bound_ln_to_gn [i]);
      PDM_free(pmesh->is_owner_concat_bound  [i]);
      PDM_free(pmesh->is_compute_concat_bound[i]);
    }

    /* Free group */
    for(int i = 0; i < PDM_BOUND_TYPE_MAX; ++i) {
      if(pmesh->is_owner_part_bound[i] == PDM_TRUE) {
        for(int i_part = 0; i_part < pmesh->n_part; ++i_part) {
          if(pmesh->ppart_bound_proc_idx[i] != NULL) {
            PDM_free(pmesh->ppart_bound_proc_idx[i][i_part]);
          }
          if(pmesh->ppart_bound_part_idx[i] != NULL) {
            PDM_free(pmesh->ppart_bound_part_idx[i][i_part]);
          }
          if(pmesh->ppart_bound[i] != NULL) {
            PDM_free(pmesh->ppart_bound[i][i_part]);
          }
        }
      }
      PDM_free(pmesh->ppart_bound_proc_idx[i]);
      PDM_free(pmesh->ppart_bound_part_idx[i]);
      PDM_free(pmesh->ppart_bound         [i]);
    }

    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      if(pmesh->pcg_ownership[i] == PDM_OWNERSHIP_KEEP) {
        PDM_part_comm_graph_free(pmesh->pcg[i]);
      }
    }

    for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
      PDM_free(pmesh->pn_entity[i]);
    }

    PDM_free(pmesh);
  }
}

void
PDM_part_mesh_dump_ensight
(
 PDM_part_mesh_t *pmesh,
 const char      *directory,
 const char      *name,
 PDM_bool_t       export_bounds
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(pmesh->comm, &i_rank);
  PDM_MPI_Comm_size(pmesh->comm, &n_rank);


  int part_offset = 0;
  PDM_MPI_Scan(&pmesh->n_part, &part_offset, 1, PDM_MPI_INT, PDM_MPI_SUM, pmesh->comm);


  /* Compute mesh highest dimension */
  int mesh_dimension = 2;
  int tn_cell = 0;
  if (pmesh->pn_entity[PDM_MESH_ENTITY_CELL] != NULL) {
    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      tn_cell += pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
    }
  }

  int max_tn_cell = 0;
  PDM_MPI_Allreduce(&tn_cell, &max_tn_cell, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, pmesh->comm);

  if (max_tn_cell > 0) {
    mesh_dimension = 3;
  }

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        directory,
                                        name,
                                        pmesh->comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  /* Define main geometry */
  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       pmesh->n_part);

  /* Define bound geometries and perform extractions */
  int id_geom_bound[PDM_BOUND_TYPE_MAX];
  PDM_extract_part_t *extrp[PDM_BOUND_TYPE_MAX];

  // Idea: move bound_type_name to pdm.h and generalize to
  // PDM_mesh_entities_t, PDM_geometry_kind_t, PDM_Mesh_nodal_elt_t, etc ?
  const char *bound_type_name[PDM_BOUND_TYPE_MAX];
  bound_type_name[PDM_BOUND_TYPE_VTX ] = "corners";
  bound_type_name[PDM_BOUND_TYPE_EDGE] = "ridges";
  bound_type_name[PDM_BOUND_TYPE_FACE] = "surfaces";
  bound_type_name[PDM_BOUND_TYPE_CELL] = "volumes";
  bound_type_name[PDM_BOUND_TYPE_ELMT] = "???";

  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
    id_geom_bound[bound_type] = -1;
    extrp        [bound_type] = NULL;

    if (bound_type == PDM_BOUND_TYPE_VTX) {
      continue;
    }

    // log_trace("%d %s\n", pmesh->n_group_bnd[bound_type], bound_type_name[bound_type]);

    if (export_bounds && pmesh->n_group_bnd[bound_type] > 0) {
      char geom_name[999];
      sprintf(geom_name, "%s_%s", name, bound_type_name[bound_type]);
      id_geom_bound[bound_type] = PDM_writer_geom_create(wrt,
                                                         geom_name,
                                                         pmesh->n_part);

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        if (!pmesh->is_compute_concat_bound[bound_type][i_part]) {
          pmesh->is_owner_concat_bound[bound_type][i_part] = PDM_TRUE;
          PDM_part_mesh_bound_concat_compute(pmesh, i_part, (PDM_bound_type_t) bound_type);
        }
      }

      _build_extract_part_bound(pmesh,
             (PDM_bound_type_t) bound_type,
                                &extrp[bound_type]);
    }
  }

  /* Define variables */
  int id_var_num_part = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "num_part");

  int id_var_bound_id = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "bound_id");

  int id_var_bound_type = PDM_writer_var_create(wrt,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAR,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "bound_type");

  PDM_writer_step_beg(wrt, 0.);

  /* Set coordinates */
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              i_part,
                              pmesh->pn_entity[PDM_MESH_ENTITY_VTX][i_part],
                              pmesh->vtx_coords[i_part],
                              pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_VTX][i_part],
                              PDM_OWNERSHIP_USER);

    /* Bounds */
    for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
      if (id_geom_bound[bound_type] >= 0) {

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        double      *vtx_coord    = NULL;
        int n_vtx = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);
        PDM_extract_part_vtx_coord_get(extrp[bound_type],
                                       i_part,
                                       &vtx_coord,
                                       PDM_OWNERSHIP_KEEP);


        PDM_writer_geom_coord_set(wrt,
                                  id_geom_bound[bound_type],
                                  i_part,
                                  n_vtx,
                                  vtx_coord,
                                  vtx_ln_to_gn,
                                  PDM_OWNERSHIP_USER);

      }
    }
  }

  /* Set elements */
  int has_face_vtx = (pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] != NULL);
  if (!has_face_vtx) {
    assert(pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] != NULL);
  }

  int **pface_vtx_idx = NULL;
  int **pface_vtx     = NULL;
  PDM_malloc(pface_vtx_idx, pmesh->n_part, int *);
  PDM_malloc(pface_vtx    , pmesh->n_part, int *);

  if (mesh_dimension == 3) {
    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      if (has_face_vtx) {
        pface_vtx_idx[i_part] = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
        pface_vtx    [i_part] = pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part];
      } else {
        // Compute face->vtx if missing
        pface_vtx_idx[i_part] = pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part];
        PDM_compute_face_vtx_from_face_and_edge(pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                                pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                                pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                                pmesh->pconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX] [i_part],
                                                &pface_vtx[i_part]);
      }
    }

    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      PDM_writer_geom_cell3d_cellface_add(wrt,
                                          id_geom,
                                          i_part,
                                          pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part],
                                          pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                          pface_vtx_idx[i_part],
                                          NULL,
                                          pface_vtx[i_part],
                                          pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part],
                                          NULL,
                                          pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE][i_part],
                                          pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_CELL][i_part]);
    }
  } else {

    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      if (has_face_vtx) {
        PDM_writer_geom_faces_facesom_add(wrt,
                                          id_geom,
                                          i_part,
                                          pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                          pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part],
                                          NULL,
                                          pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX][i_part],
                                          pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part]);
      } else {
        PDM_writer_geom_cell2d_cellface_add(wrt,
                                            id_geom,
                                            i_part,
                                            pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part],
                                            pmesh->pn_entity[PDM_MESH_ENTITY_EDGE][i_part],
                                            NULL,
                                            NULL,
                                            pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX][i_part],
                                            pmesh->pconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                            NULL,
                                            pmesh->pconnectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE][i_part],
                                            pmesh->pentity_ln_to_gn[PDM_MESH_ENTITY_FACE][i_part]);
      }
    }
  }

  /* Set bound elements */
  int **vtx_vtx = NULL;
  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {

    if (id_geom_bound[bound_type] < 0) continue;

    int id_block = -1;

    if (bound_type == PDM_BOUND_TYPE_CELL) {

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        PDM_g_num_t *cell_ln_to_gn = NULL;
        int n_cell = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_CELL,
                                                   &cell_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int *cell_face     = NULL;
        int *cell_face_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                          &cell_face,
                                          &cell_face_idx,
                                          PDM_OWNERSHIP_KEEP);

        int *face_vtx     = NULL;
        int *face_vtx_idx = NULL;
        int n_face = PDM_extract_part_connectivity_get(extrp[bound_type],
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                       &face_vtx,
                                                       &face_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);

        int *_face_vtx_idx = face_vtx_idx;

        if (face_vtx_idx == NULL) {
          int *face_edge     = NULL;
          n_face = PDM_extract_part_connectivity_get(extrp[bound_type],
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge,
                                                     &_face_vtx_idx,
                                                     PDM_OWNERSHIP_KEEP);

          int *edge_vtx     = NULL;
          int *edge_vtx_idx = NULL;
          PDM_extract_part_connectivity_get(extrp[bound_type],
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &edge_vtx,
                                            &edge_vtx_idx,
                                            PDM_OWNERSHIP_KEEP);

          PDM_compute_face_vtx_from_face_and_edge(n_face,
                                                  _face_vtx_idx,
                                                  face_edge,
                                                  edge_vtx,
                                                  &face_vtx);
        }

        PDM_writer_geom_cell3d_cellface_add(wrt,
                                            id_geom_bound[bound_type],
                                            i_part,
                                            n_cell,
                                            n_face,
                                            _face_vtx_idx,
                                            NULL,
                                            face_vtx,
                                            cell_face_idx,
                                            NULL,
                                            cell_face,
                                            cell_ln_to_gn);

        if (face_vtx_idx == NULL) {
          PDM_free(face_vtx);
        }
      }

    }
    else if (bound_type == PDM_BOUND_TYPE_FACE) {

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        PDM_g_num_t *face_ln_to_gn = NULL;
        int n_face = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_FACE,
                                                   &face_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        int n_edge = PDM_extract_part_connectivity_get(extrp[bound_type],
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);

        int *face_vtx     = NULL;
        int *face_vtx_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          &face_vtx,
                                          &face_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

        int *face_edge     = NULL;
        int *face_edge_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          &face_edge,
                                          &face_edge_idx,
                                          PDM_OWNERSHIP_KEEP);

        if (face_vtx_idx != NULL) {
          PDM_writer_geom_faces_facesom_add(wrt,
                                            id_geom_bound[bound_type],
                                            i_part,
                                            n_face,
                                            face_vtx_idx,
                                            NULL,
                                            face_vtx,
                                            face_ln_to_gn);
        } else {
          assert(face_edge_idx != NULL);
          PDM_writer_geom_cell2d_cellface_add(wrt,
                                              id_geom_bound[bound_type],
                                              i_part,
                                              n_face,
                                              n_edge,
                                              NULL,
                                              NULL,
                                              edge_vtx,
                                              face_edge_idx,
                                              NULL,
                                              face_edge,
                                              face_ln_to_gn);
        }
      }
    }
    else if (bound_type == PDM_BOUND_TYPE_EDGE) {
      id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom_bound[bound_type],
                                          PDM_WRITER_BAR2,
                                          PDM_OWNERSHIP_USER);

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
        PDM_g_num_t *edge_ln_to_gn = NULL;
        int n_edge = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                   i_part,
                                                   PDM_MESH_ENTITY_EDGE,
                                                   &edge_ln_to_gn,
                                                   PDM_OWNERSHIP_KEEP);

        int *edge_vtx     = NULL;
        int *edge_vtx_idx = NULL;
        PDM_extract_part_connectivity_get(extrp[bound_type],
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx,
                                          &edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);

        PDM_writer_geom_bloc_std_set(wrt,
                                     id_geom_bound[bound_type],
                                     id_block,
                                     i_part,
                                     n_edge,
                                     edge_vtx,
                                     edge_ln_to_gn);
      }
    }
    else if (bound_type == PDM_BOUND_TYPE_VTX) {

      id_block = PDM_writer_geom_bloc_add(wrt,
                                          id_geom_bound[bound_type],
                                          PDM_WRITER_POINT,
                                          PDM_OWNERSHIP_USER);

      PDM_malloc(vtx_vtx, pmesh->n_part, int *);

      for (int i_part = 0; i_part < pmesh->n_part; i_part++) {

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        int n_vtx = PDM_extract_part_ln_to_gn_get(extrp[bound_type],
                                                  i_part,
                                                  PDM_MESH_ENTITY_VTX,
                                                  &vtx_ln_to_gn,
                                                  PDM_OWNERSHIP_KEEP);

        PDM_malloc(vtx_vtx[i_part],n_vtx,int);
        for (int i = 0; i < n_vtx; i++) {
          vtx_vtx[i_part][i] = i + 1;
        }

        PDM_writer_geom_bloc_std_set(wrt,
                                     id_geom_bound[bound_type],
                                     id_block,
                                     i_part,
                                     n_vtx,
                                     vtx_vtx[i_part],
                                     vtx_ln_to_gn);
      }
    }
  }

  /* Write geometries */
  PDM_writer_geom_write(wrt, id_geom);

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    if (id_geom_bound[i] >= 0) {
      PDM_writer_geom_write(wrt, id_geom_bound[i]);
    }
  }

  /* Set variables */
  PDM_real_t **val_num_part   = NULL;
  PDM_real_t **val_bound_id   = NULL;
  PDM_real_t **val_bound_type = NULL;
  PDM_malloc(val_num_part,   pmesh->n_part, PDM_real_t *);
  PDM_malloc(val_bound_id,   pmesh->n_part, PDM_real_t *);
  PDM_malloc(val_bound_type, pmesh->n_part, PDM_real_t *);
  for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
    int n_entity = 0;
    if (mesh_dimension == 2) {
      n_entity = pmesh->pn_entity[PDM_MESH_ENTITY_FACE][i_part];
    } else {
      n_entity = pmesh->pn_entity[PDM_MESH_ENTITY_CELL][i_part];
    }

    PDM_malloc(val_num_part  [i_part], n_entity, PDM_real_t);
    PDM_malloc(val_bound_id  [i_part], n_entity, PDM_real_t);
    PDM_malloc(val_bound_type[i_part], n_entity, PDM_real_t);
    for (int i = 0; i < n_entity; i++) {
      val_num_part  [i_part][i] = part_offset + i_part;
      val_bound_id  [i_part][i] = 0;
      val_bound_type[i_part][i] = -1;
    }

    PDM_writer_var_set(wrt,
                       id_var_num_part,
                       id_geom,
                       i_part,
                       val_num_part[i_part]);
    PDM_free(val_num_part[i_part]);

    PDM_writer_var_set(wrt,
                       id_var_bound_id,
                       id_geom,
                       i_part,
                       val_bound_id[i_part]);
    PDM_free(val_bound_id[i_part]);

    PDM_writer_var_set(wrt,
                       id_var_bound_type,
                       id_geom,
                       i_part,
                       val_bound_type[i_part]);
    PDM_free(val_bound_type[i_part]);
  }


  for (int bound_type = 0; bound_type < PDM_BOUND_TYPE_MAX; bound_type++) {
    if (id_geom_bound[bound_type] < 0) {
      continue;
    }

    for (int i_part = 0; i_part < pmesh->n_part; i_part++) {
      int n_entity = 0;
      if (bound_type == PDM_BOUND_TYPE_CELL) {
        n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                 i_part,
                                                 PDM_MESH_ENTITY_CELL);
      }
      else if (bound_type == PDM_BOUND_TYPE_FACE) {
        n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                 i_part,
                                                 PDM_MESH_ENTITY_FACE);
      }
      else if (bound_type == PDM_BOUND_TYPE_EDGE) {
        n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                 i_part,
                                                 PDM_MESH_ENTITY_EDGE);
      }
      else if (bound_type == PDM_BOUND_TYPE_VTX) {
        n_entity = PDM_extract_part_n_entity_get(extrp[bound_type],
                                                 i_part,
                                                 PDM_MESH_ENTITY_VTX);
      }

      PDM_malloc(val_num_part  [i_part], n_entity, PDM_real_t);
      PDM_malloc(val_bound_type[i_part], n_entity, PDM_real_t);
      for (int i_entity = 0; i_entity < n_entity; i_entity++) {
        val_num_part  [i_part][i_entity] = part_offset + i_part;
        val_bound_type[i_part][i_entity] = bound_type;
      }

      PDM_malloc(val_bound_id[i_part], n_entity, PDM_real_t);
      for (int bound_id = 0; bound_id < pmesh->n_group_bnd[bound_type]; bound_id++) {
        for (int k = pmesh->pconcat_bound_idx[bound_type][i_part][bound_id]; k < pmesh->pconcat_bound_idx[bound_type][i_part][bound_id+1]; k++) {
          val_bound_id[i_part][k] = bound_id+1;
        }
      }

      PDM_writer_var_set(wrt,
                          id_var_num_part,
                          id_geom_bound[bound_type],
                          i_part,
                          val_num_part[i_part]);
      PDM_free(val_num_part[i_part]);

      PDM_writer_var_set(wrt,
                          id_var_bound_id,
                          id_geom_bound[bound_type],
                          i_part,
                          val_bound_id[i_part]);
      PDM_free(val_bound_id[i_part]);

      PDM_writer_var_set(wrt,
                          id_var_bound_type,
                          id_geom_bound[bound_type],
                          i_part,
                          val_bound_type[i_part]);
      PDM_free(val_bound_type[i_part]);
    } // End loop on parts
  } // End loop on bound types


  /* Write variables */
  PDM_writer_var_write(wrt,
                       id_var_num_part);

  PDM_writer_var_write(wrt,
                       id_var_bound_id);

  PDM_writer_var_write(wrt,
                       id_var_bound_type);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);


  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    PDM_extract_part_free(extrp[i]);
  }

  if (mesh_dimension == 3 && !has_face_vtx) {
    for (int i = 0; i < pmesh->n_part; i++) {
      PDM_free(pface_vtx[i]);
    }
  }
  PDM_free(pface_vtx_idx );
  PDM_free(pface_vtx     );
  PDM_free(val_num_part  );
  PDM_free(val_bound_id  );
  PDM_free(val_bound_type);

  if (vtx_vtx != NULL) {
    for (int i = 0; i < pmesh->n_part; i++) {
      PDM_free(vtx_vtx[i]);
    }
    PDM_free(vtx_vtx);
  }
}
