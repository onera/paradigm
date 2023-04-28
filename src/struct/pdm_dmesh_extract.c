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
  PDM_UNUSED(dme);

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

  dme->dim                   = dim;
  dme->comm                  = comm;

  // Utilisation privé
  dme->dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP, 0, 0, 0, 0, comm);


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

  if(dme->dim == 3) {
    _dmesh_extract_3d(dme);
  } else if(dme->dim == 2) {
    _dmesh_extract_2d(dme);
  } else if(dme->dim == 1) {
    _dmesh_extract_1d(dme);
  } else {
    _dmesh_extract_0d(dme);
  }



}


void
PDM_dmesh_extract_selected_gnum_set
(
 PDM_dmesh_extract_t *dme,
 PDM_mesh_entities_t  entity_type,
 int                  dn_extract_entity,
 PDM_g_num_t         *extract_gnum
)
{

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
PDM_dmesh_extract_free
(
  PDM_dmesh_extract_t  *dme
)
{


  PDM_dmesh_free(dme->dmesh);


  free(dme);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
