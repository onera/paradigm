/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_unique.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"

#include "pdm_part_domain_interface.h"
#include "pdm_part_domain_interface_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/


/*============================================================================
 * Public function definitions
 *============================================================================*/
PDM_part_domain_interface_t*
PDM_part_domain_interface_create
(
const int                   n_interface,
const int                   n_domain,
const int                   n_part,
PDM_domain_interface_mult_t multidomain_interface,
PDM_ownership_t             ownership,
PDM_MPI_Comm                comm
)
{
  PDM_part_domain_interface_t *dom_intrf = (PDM_part_domain_interface_t *) malloc (sizeof(PDM_part_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  dom_intrf->n_part            = n_part;
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_pn_vtx        = (int          **) malloc(n_part * sizeof(int          *));
  dom_intrf->interface_vtx_ln_to_gn  = (PDM_g_num_t ***) malloc(n_part * sizeof(PDM_g_num_t **));
  dom_intrf->interface_ids_vtx       = (int         ***) malloc(n_part * sizeof(int         **));
  dom_intrf->interface_dom_vtx       = (int         ***) malloc(n_part * sizeof(int         **));

  dom_intrf->interface_pn_face       = (int         ** ) malloc(n_part * sizeof(int          *));
  dom_intrf->interface_face_ln_to_gn = (PDM_g_num_t ***) malloc(n_part * sizeof(PDM_g_num_t **));
  dom_intrf->interface_ids_face      = (int         ***) malloc(n_part * sizeof(int         **));
  dom_intrf->interface_dom_face      = (int         ***) malloc(n_part * sizeof(int         **));

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++)
    dom_intrf->is_result[i] = 0;

  return dom_intrf;
}

void
PDM_part_domain_interface_set
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind,
 int                           i_part,
 int                          *interface_pn,
 PDM_g_num_t                 **interface_ln_to_gn,
 int                         **interface_ids,
 int                         **interface_dom
)
{
  assert(i_part < dom_intrf->n_part);
  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_pn_vtx      [i_part] = interface_pn;
    dom_intrf->interface_vtx_ln_to_gn[i_part] = interface_ln_to_gn;
    dom_intrf->interface_ids_vtx     [i_part] = interface_ids;
    dom_intrf->interface_dom_vtx     [i_part] = interface_dom;
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_pn_face      [i_part] = interface_pn;
    dom_intrf->interface_face_ln_to_gn[i_part] = interface_ln_to_gn;
    dom_intrf->interface_ids_face     [i_part] = interface_ids;
    dom_intrf->interface_dom_face     [i_part] = interface_dom;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}



void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
)
{

  free(dom_intrf->interface_pn_vtx       );
  free(dom_intrf->interface_vtx_ln_to_gn );
  free(dom_intrf->interface_ids_vtx      );
  free(dom_intrf->interface_dom_vtx      );

  free(dom_intrf->interface_pn_face      );
  free(dom_intrf->interface_face_ln_to_gn);
  free(dom_intrf->interface_ids_face     );
  free(dom_intrf->interface_dom_face     );


  free(dom_intrf);
}
