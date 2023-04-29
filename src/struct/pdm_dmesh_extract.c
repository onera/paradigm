/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_part_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_dmesh_extract.h"
#include "pdm_dmesh_extract_priv.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_gnum_location.h"
#include "pdm_unique.h"
#include "pdm_dmesh_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
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

static
void
_dmesh_extract_3d
(
 PDM_dmesh_extract_t *dme
)
{
  PDM_UNUSED(dme);

}


static
void
_dmesh_extract_2d
(
 PDM_dmesh_extract_t *dme
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dme->comm, &i_rank);

  int from_face_edge = 0;
  int from_face_vtx  = 0;

  PDM_g_num_t *dface_vtx     = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_vtx_idx != NULL){
    from_face_vtx = 1;
  }

  // face edge
  PDM_g_num_t *dface_edge     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_dmesh_connectivity_get(dme->dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_edge_idx != NULL) {
    from_face_edge = 1;
  }

  int dn_face = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_FACE  );
  int dn_edge = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_EDGE  );
  int dn_vtx  = PDM_dmesh_dn_entity_get(dme->dmesh, PDM_MESH_ENTITY_VERTEX);

  PDM_g_num_t* distrib_face = PDM_compute_entity_distribution(dme->comm, dn_face);
  PDM_g_num_t* distrib_vtx  = PDM_compute_entity_distribution(dme->comm, dn_vtx );
  PDM_g_num_t* distrib_edge = NULL;

  if(from_face_edge == 1) {

    PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                               dme->n_selected,
                                               dme->selected_gnum,
                                               distrib_face,
                                               dface_edge_idx,
                                               dface_edge,
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_FACE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_FACE],
                                               &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                               &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_EDGE],
                                               &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_EDGE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_EDGE]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;

    // edge_vtx
    distrib_edge = PDM_compute_entity_distribution(dme->comm, dn_edge);
    PDM_g_num_t *dedge_vtx     = NULL;
    int         *dedge_vtx_idx = NULL;
    PDM_dmesh_connectivity_get(dme->dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
    int *_dedge_vtx_idx = NULL;
    if(dedge_vtx_idx == NULL)  {
      _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
        _dedge_vtx_idx[i_edge] = 2*i_edge;
      }
    } else {
      _dedge_vtx_idx = dedge_vtx_idx;
    }
    dme->dmesh_extract->dn_edge = dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_EDGE][i_rank];
    PDM_dconnectivity_to_extract_dconnectivity_block(dme->comm,
                                                     dme->dmesh_extract->dn_edge,
                                                     dme->parent_extract_gnum[PDM_MESH_ENTITY_EDGE],
                                                     distrib_edge,
                                                     _dedge_vtx_idx,
                                                     dedge_vtx,
                                                     &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_EDGE_VTX],
                                                     &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_EDGE],
                                                     &dme->distrib_extract                 [PDM_MESH_ENTITY_VERTEX],
                                                     &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VERTEX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = PDM_TRUE;
    free(dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX]);
    dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_EDGE_VTX] = NULL;
    if(dedge_vtx_idx == NULL)  {
      free(_dedge_vtx_idx);
    }

  } else if(from_face_vtx == 1) {
    PDM_dconnectivity_to_extract_dconnectivity(dme->comm,
                                               dme->n_selected,
                                               dme->selected_gnum,
                                               distrib_face,
                                               dface_vtx_idx,
                                               dface_vtx,
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_FACE],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_FACE],
                                               &dme->dmesh_extract->dconnectivity_idx[PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                               &dme->dmesh_extract->dconnectivity    [PDM_CONNECTIVITY_TYPE_FACE_VTX],
                                               &dme->btp_entity_to_extract_entity    [PDM_MESH_ENTITY_FACE],
                                               &dme->distrib_extract                 [PDM_MESH_ENTITY_VERTEX],
                                               &dme->parent_extract_gnum             [PDM_MESH_ENTITY_VERTEX]);

    dme->dmesh_extract->is_owner_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX] = PDM_TRUE;
    dme->dmesh_extract->dn_face = dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_FACE][i_rank];
  }

  dme->dmesh_extract->dn_vtx = dme->distrib_extract[PDM_MESH_ENTITY_VERTEX][i_rank+1] - dme->distrib_extract[PDM_MESH_ENTITY_VERTEX][i_rank];
  dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VERTEX] = PDM_block_to_part_create(distrib_vtx,
                                                                (const PDM_g_num_t **) &dme->parent_extract_gnum[PDM_MESH_ENTITY_VERTEX],
                                                                                       &dme->dmesh_extract->dn_vtx,
                                                                                       1,
                                                                                       dme->comm);


  free(distrib_face);
  if(distrib_edge != NULL) {
    free(distrib_edge);
  }
  free(distrib_vtx );
}


static
void
_dmesh_extract_1d
(
 PDM_dmesh_extract_t *dme
)
{
  PDM_UNUSED(dme);

}

static
void
_dmesh_extract_0d
(
 PDM_dmesh_extract_t *dme
)
{
  PDM_UNUSED(dme);

}



/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_extract_t*
PDM_dmesh_extract_create
(
 const int                     dim,
       PDM_MPI_Comm            comm
)
{
  PDM_dmesh_extract_t *dme = (PDM_dmesh_extract_t *) malloc(sizeof(PDM_dmesh_extract_t));

  dme->dim  = dim;
  dme->comm = comm;

  // Utilisation privé
  dme->dmesh         = PDM_dmesh_create(PDM_OWNERSHIP_KEEP, 0, 0, 0, 0, comm);
  dme->dmesh_extract = NULL;
  dme->dmesh_extract_ownership = PDM_OWNERSHIP_KEEP;

  for(int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    dme->btp_entity_to_extract_entity [i] = NULL;
    dme->btp_ownership                [i] = PDM_OWNERSHIP_KEEP;

    dme->distrib_extract              [i] = NULL;
    dme->parent_extract_gnum          [i] = NULL;

    dme->distrib_extract_ownership    [i] = PDM_OWNERSHIP_KEEP;
    dme->parent_extract_gnum_ownership[i] = PDM_OWNERSHIP_KEEP;
  }

  return dme;
}

void
PDM_dmesh_extract_compute
(
 PDM_dmesh_extract_t *dme
)
{
  // Synchronize dmesh
  PDM_g_num_t _dn_cell = dme->dmesh->dn_cell;
  PDM_g_num_t _dn_face = dme->dmesh->dn_face;
  PDM_g_num_t _dn_edge = dme->dmesh->dn_edge;
  PDM_g_num_t _dn_vtx  = dme->dmesh->dn_vtx;

  PDM_MPI_Allreduce(&_dn_cell, &dme->dmesh->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_face, &dme->dmesh->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_edge, &dme->dmesh->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_vtx , &dme->dmesh->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);

  dme->dmesh_extract = PDM_dmesh_create(PDM_OWNERSHIP_KEEP, 0, 0, 0, 0, dme->comm);
  dme->dmesh_extract->_dvtx_coord = NULL;

  if(dme->dim == 3) {
    _dmesh_extract_3d(dme);
  } else if(dme->dim == 2) {
    _dmesh_extract_2d(dme);
  } else if(dme->dim == 1) {
    _dmesh_extract_1d(dme);
  } else {
    _dmesh_extract_0d(dme);
  }

  /* Exchange coordinates */
  if(dme->dmesh->_dvtx_coord != NULL) {
    double** tmp_dvtx_coord   = NULL;
    int stride_one = 1;
    PDM_block_to_part_exch(dme->btp_entity_to_extract_entity[PDM_MESH_ENTITY_VERTEX],
                           3 * sizeof(double),
                           PDM_STRIDE_CST_INTERLACED,
                           &stride_one,
               (void *  )  dme->dmesh->_dvtx_coord,
                           NULL,
               (void ***)  &tmp_dvtx_coord);
    PDM_dmesh_vtx_coord_set(dme->dmesh_extract,
                            tmp_dvtx_coord[0],
                            PDM_OWNERSHIP_KEEP);
    free(tmp_dvtx_coord);
  }

  PDM_g_num_t _dn_extract_cell = dme->dmesh_extract->dn_cell;
  PDM_g_num_t _dn_extract_face = dme->dmesh_extract->dn_face;
  PDM_g_num_t _dn_extract_edge = dme->dmesh_extract->dn_edge;
  PDM_g_num_t _dn_extract_vtx  = dme->dmesh_extract->dn_vtx;

  PDM_MPI_Allreduce(&_dn_extract_cell, &dme->dmesh_extract->n_g_cell, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_extract_face, &dme->dmesh_extract->n_g_face, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_extract_edge, &dme->dmesh_extract->n_g_edge, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);
  PDM_MPI_Allreduce(&_dn_extract_vtx , &dme->dmesh_extract->n_g_vtx , 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, dme->comm);

}


void
PDM_dmesh_extract_selected_gnum_set
(
 PDM_dmesh_extract_t *dme,
 PDM_mesh_entities_t  entity_type,
 int                  n_selected,
 PDM_g_num_t         *selected_gnum
)
{
  PDM_UNUSED(entity_type);
  dme->n_selected    = n_selected;
  dme->selected_gnum = selected_gnum;

}

void
PDM_dmesh_extract_dn_entity_set
(
 PDM_dmesh_extract_t *dme,
 PDM_mesh_entities_t  entity_type,
 int                  dn_entity
)
{
  PDM_dmesh_dn_entity_set(dme->dmesh, entity_type, dn_entity);
}


void
PDM_dmesh_extract_vtx_coord_set
(
 PDM_dmesh_extract_t *dme,
 double              *dvtx_coord
)
{
  PDM_dmesh_vtx_coord_set(dme->dmesh, dvtx_coord, PDM_OWNERSHIP_USER);
}

void
PDM_dmesh_extract_dmesh_bound_set
(
 PDM_dmesh_extract_t *dme,
 PDM_bound_type_t     bound_type,
 int                  n_bound,
 PDM_g_num_t         *connect,
 int                 *connect_idx
)
{
  PDM_dmesh_bound_set(dme->dmesh, bound_type, n_bound, connect, connect_idx, PDM_OWNERSHIP_USER);
}


void
PDM_dmesh_extract_dconnectivity_set
(
       PDM_dmesh_extract_t     *dme,
       PDM_connectivity_type_t  connectivity_type,
       PDM_g_num_t             *dconnect,
       int                     *dconnect_idx
)
{
  PDM_dmesh_connectivity_set(dme->dmesh,
                             connectivity_type,
                             dconnect,
                             dconnect_idx,
                             PDM_OWNERSHIP_USER);
}


void
PDM_dmesh_extract_dmesh_get
(
 PDM_dmesh_extract_t     *dme,
 PDM_dmesh_t            **dmesh_extract,
 PDM_ownership_t          ownership
)
{
  *dmesh_extract = dme->dmesh_extract;
  dme->dmesh_extract_ownership = ownership;

}

void
PDM_dmesh_extract_free
(
  PDM_dmesh_extract_t  *dme
)
{


  PDM_dmesh_free(dme->dmesh);

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {
    if (dme->btp_ownership[i] == PDM_OWNERSHIP_KEEP && dme->btp_entity_to_extract_entity[i] != NULL) {
      PDM_block_to_part_free(dme->btp_entity_to_extract_entity[i]);
    }
  }


  for (int i = 0; i < PDM_MESH_ENTITY_MAX; ++i) {

    if(dme->distrib_extract_ownership[i] == PDM_OWNERSHIP_KEEP && dme->distrib_extract[i] != NULL) {
      free(dme->distrib_extract[i]);
    }

    if(dme->parent_extract_gnum_ownership[i] == PDM_OWNERSHIP_KEEP && dme->parent_extract_gnum[i] != NULL) {
      free(dme->parent_extract_gnum[i]);
    }

  }

  if(dme->dmesh_extract_ownership == PDM_OWNERSHIP_KEEP && dme->dmesh_extract != NULL) {
    PDM_dmesh_free(dme->dmesh_extract);
  }

  free(dme);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
