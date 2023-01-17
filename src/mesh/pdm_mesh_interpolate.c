
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_interpolate_priv.h"
#include "pdm_mesh_interpolate.h"

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
 * \brief Create a structure that compute a global mean
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_mesh_interpolate object
 */
PDM_mesh_interpolate_t*
PDM_mesh_interpolate_create
(
 const int            n_domain,
 const int           *n_part,
 const int           *n_group,
 const int            interp_kind,
 const PDM_MPI_Comm   comm
)
{
  PDM_mesh_interpolate_t* mi = malloc(sizeof(PDM_mesh_interpolate_t));

  mi->n_domain    = n_domain;
  mi->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < mi->n_domain; ++i) {
    mi->n_part[i] = n_part[i];
  }
  mi->comm        = comm;

  mi->n_part_idx    = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->n_part_g_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->parts = malloc(n_domain * sizeof(_part_t *));


  mi->n_part_idx  [0] = 0;
  mi->n_part_g_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    mi->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
    mi->n_part_idx[i_domain+1] = mi->n_part_idx[i_domain] + mi->n_part[i_domain];

    int n_part_l = n_part[i_domain];
    int n_part_g = -100;
    PDM_MPI_Allreduce(&n_part_l, &n_part_g, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    mi->n_part_g_idx[i_domain+1] = mi->n_part_idx[i_domain] + n_part_g;

  }

  /* Graph comm  */
  mi->entity_part_bound_proc_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->entity_part_bound_part_idx = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->entity_part_bound          = malloc( PDM_MESH_ENTITY_MAX * sizeof(int ***));
  mi->graph_comm_is_defined      = malloc( PDM_MESH_ENTITY_MAX * sizeof(int *  ));

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    mi->entity_part_bound_proc_idx[i_kind] = malloc(n_domain * sizeof(int **));
    mi->entity_part_bound_part_idx[i_kind] = malloc(n_domain * sizeof(int **));
    mi->entity_part_bound         [i_kind] = malloc(n_domain * sizeof(int **));
    mi->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
      mi->entity_part_bound_proc_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      mi->entity_part_bound_part_idx[i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      mi->entity_part_bound         [i_kind][i_domain] = malloc(n_part[i_domain] * sizeof(int *));
      for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
        mi->entity_part_bound_proc_idx[i_kind][i_domain][i_part] = NULL;
        mi->entity_part_bound_part_idx[i_kind][i_domain][i_part] = NULL;
        mi->entity_part_bound         [i_kind][i_domain][i_part] = NULL;
      }
    }
  }


  return mi;
}

void
PDM_mesh_interpolate_compute
(
  PDM_mesh_interpolate_t *part_ext
)
{


  /*
   * Compute graph comm from gnum with PDM_part_generate_entity_graph_comm is not provided
   */

  /* Deduce graph with all graphe inside same domain and between domain */
  /* Compute weight */
  /* Create protocol only between join - Distant neighbor */
}



void
PDM_mesh_interpolate_part_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
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
  mi->parts[i_domain][i_part].n_cell            = n_cell;
  mi->parts[i_domain][i_part].n_face            = n_face;
  mi->parts[i_domain][i_part].n_edge            = n_edge;
  mi->parts[i_domain][i_part].n_vtx             = n_vtx;

  mi->parts[i_domain][i_part].cell_face_idx = cell_face_idx;
  mi->parts[i_domain][i_part].cell_face     = cell_face;

  mi->parts[i_domain][i_part].face_edge_idx = face_edge_idx;
  mi->parts[i_domain][i_part].face_edge     = face_edge;

  mi->parts[i_domain][i_part].face_vtx_idx  = face_vtx_idx;
  mi->parts[i_domain][i_part].face_vtx      = face_vtx;

  mi->parts[i_domain][i_part].edge_vtx      = edge_vtx;

  mi->parts[i_domain][i_part].cell_ln_to_gn = cell_ln_to_gn;
  mi->parts[i_domain][i_part].face_ln_to_gn = face_ln_to_gn;
  mi->parts[i_domain][i_part].edge_ln_to_gn = edge_ln_to_gn;
  mi->parts[i_domain][i_part].vtx_ln_to_gn  = vtx_ln_to_gn;

  mi->parts[i_domain][i_part].vtx = vtx_coord;
}


void
PDM_mesh_interpolate_graph_comm_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                      *entity_part_bound_proc_idx,
  int                      *entity_part_bound_part_idx,
  int                      *entity_part_bound
)
{
  mi->entity_part_bound_proc_idx[mesh_entity][i_domain][i_part] = entity_part_bound_proc_idx;
  mi->entity_part_bound_part_idx[mesh_entity][i_domain][i_part] = entity_part_bound_part_idx;
  mi->entity_part_bound         [mesh_entity][i_domain][i_part] = entity_part_bound;
  mi->graph_comm_is_defined[mesh_entity] = 1;
}


void
PDM_mesh_interpolate_part_domain_interface_shared_set
(
  PDM_mesh_interpolate_t      *mi,
  PDM_part_domain_interface_t *pdi
)
{
  mi->pdi = pdi;
}


void
PDM_mesh_interpolate_free
(
 PDM_mesh_interpolate_t *mi
)
{

  for(int i_kind = 0; i_kind < PDM_MESH_ENTITY_MAX; ++i_kind) {
    mi->graph_comm_is_defined     [i_kind] = 0;
    for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
      free(mi->entity_part_bound_proc_idx[i_kind][i_domain]);
      free(mi->entity_part_bound_part_idx[i_kind][i_domain]);
      free(mi->entity_part_bound         [i_kind][i_domain]);
    }
    free(mi->entity_part_bound_proc_idx[i_kind]);
    free(mi->entity_part_bound_part_idx[i_kind]);
    free(mi->entity_part_bound         [i_kind]);
  }
  free(mi->entity_part_bound_proc_idx);
  free(mi->entity_part_bound_part_idx);
  free(mi->entity_part_bound         );
  free(mi->graph_comm_is_defined     );


  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    free(mi->parts[i_domain]);
  }
  free(mi->parts);
  free(mi->n_part_idx);
  free(mi->n_part_g_idx);
  free(mi->n_part);


  free(mi);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
