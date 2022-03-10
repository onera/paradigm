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
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_iso_surface.h"
#include "pdm_iso_surface_priv.h"


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
_iso_surface_dist
(
  PDM_iso_surface_t        *isos
)
{

  PDM_UNUSED(isos);

}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_iso_surface_t*
PDM_iso_surface_create
(
 const int             dim,
 const int             n_part,
       PDM_ownership_t ownership,
       PDM_MPI_Comm    comm
)
{
  PDM_iso_surface_t *isos = (PDM_iso_surface_t *) malloc(sizeof(PDM_iso_surface_t));

  isos->dim       = dim;
  isos->n_part    = n_part;
  isos->ownership = ownership;
  isos->comm      = comm;
  isos->is_dist   = -1;

  isos->n_cell         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_face         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_edge         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_vtx          = (int          *) malloc(n_part * sizeof(int          ));

  isos->pcell_face     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pcell_face_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pedge_vtx      = (int         **) malloc(n_part * sizeof(int         *));
  isos->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));

  isos->pface_vtx_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_vtx     = (int         **) malloc(n_part * sizeof(int         *));

  isos->pvtx_coord      = (double **) malloc(n_part * sizeof(double *));
  isos->pfield          = (double **) malloc(n_part * sizeof(double *));
  isos->pgradient_field = (double **) malloc(n_part * sizeof(double *));

  return isos;
}

void
PDM_iso_surface_compute
(
  PDM_iso_surface_t        *isos
)
{
  if(isos->is_dist == 0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_iso_surface_compute Not implemented with partition layout\n");
  } else {
    _iso_surface_dist(isos);
  }



}

// See with Eric et Bastien : par type ou une fonction avec 1000 arguments ?
void
PDM_iso_surface_part_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
)
{
  isos->n_cell        [i_part] = n_cell;
  isos->n_face        [i_part] = n_face;
  isos->n_edge        [i_part] = n_edge;
  isos->n_vtx         [i_part] = n_vtx;
  isos->pcell_face    [i_part] = cell_face;
  isos->pcell_face_idx[i_part] = cell_face_idx;
  isos->pface_edge    [i_part] = face_edge;
  isos->pface_edge_idx[i_part] = face_edge_idx;
  isos->pedge_vtx     [i_part] = edge_vtx;
  isos->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  isos->face_ln_to_gn [i_part] = face_ln_to_gn;
  isos->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  isos->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  isos->pface_vtx_idx [i_part] = face_vtx_idx;
  isos->pface_vtx     [i_part] = face_vtx;
  isos->pvtx_coord    [i_part] = vtx_coord;
}


void
PDM_iso_surface_part_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *field
)
{
  isos->is_dist = 0;
  isos->pfield[i_part] = field;
}

void
PDM_iso_surface_part_gradient_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *gradient_field
)
{
  isos->is_dist = 0;
  isos->pgradient_field[i_part] = gradient_field;
}

void
PDM_iso_surface_dconnectivity_set
(
  PDM_iso_surface_t        *isos,
  PDM_connectivity_type_t   connectivity_type,
  PDM_g_num_t              *dconnect,
  int                      *dconnect_idx
)
{
  isos->is_dist = 1;
  switch (connectivity_type) {
   case PDM_CONNECTIVITY_TYPE_CELL_FACE:
     isos->dcell_face     = dconnect;
     isos->dcell_face_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
     isos->dface_edge     = dconnect;
     isos->dface_edge_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_VTX:
     isos->dface_vtx     = dconnect;
     isos->dface_vtx_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
     isos->dedge_vtx     = dconnect;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type for iso_surface %d\n", connectivity_type);
    break;
   }

}

void
PDM_iso_surface_distrib_set
(
  PDM_iso_surface_t        *isos,
  PDM_mesh_entities_t       entity_type,
  PDM_g_num_t              *distrib_entity
)
{
  isos->is_dist = 1;
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     isos->distrib_cell = distrib_entity;
     break;
   case PDM_MESH_ENTITY_FACE:
     isos->distrib_face = distrib_entity;
     break;
   case PDM_MESH_ENTITY_EDGE:
     isos->distrib_edge = distrib_entity;
     break;
   case PDM_MESH_ENTITY_VERTEX:
     isos->distrib_vtx = distrib_entity;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid entity_type for iso_surface %d\n", entity_type);
    break;
   }
}


void
PDM_iso_surface_dvtx_coord_set
(
  PDM_iso_surface_t *isos,
  double            *dvtx_coord
)
{
  isos->is_dist     = 1;
  isos->dvtx_coord  = dvtx_coord;
}

void
PDM_iso_surface_dfield_set
(
  PDM_iso_surface_t *isos,
  double            *dfield
)
{
  isos->is_dist = 1;
  isos->dfield  = dfield;
}

void
PDM_iso_surface_dgrad_field_set
(
  PDM_iso_surface_t *isos,
  double            *dgrad_field
)
{
  isos->is_dist         = 1;
  isos->dgradient_field = dgrad_field;
}

void
PDM_iso_surface_free
(
  PDM_iso_surface_t        *isos
)
{
  free(isos->n_cell        );
  free(isos->n_face        );
  free(isos->n_edge        );
  free(isos->n_vtx         );
  free(isos->pcell_face    );
  free(isos->pcell_face_idx);
  // Si pface_edge a été calculé il faut le free
  free(isos->pface_edge    );
  free(isos->pface_edge_idx);

  free(isos->pface_vtx     );
  free(isos->pface_vtx_idx );
  free(isos->pedge_vtx     );
  free(isos->cell_ln_to_gn );
  free(isos->face_ln_to_gn );
  free(isos->edge_ln_to_gn );
  free(isos->vtx_ln_to_gn  );

  free(isos->pvtx_coord     );
  free(isos->pfield         );
  free(isos->pgradient_field);

  free(isos);
}
