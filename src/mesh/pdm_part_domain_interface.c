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
#include "pdm_distant_neighbor.h"

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
const int                    n_interface,
const int                    n_domain,
const int                   *n_part,
PDM_domain_interface_mult_t  multidomain_interface,
PDM_ownership_t              ownership,
PDM_MPI_Comm                 comm
)
{
  PDM_part_domain_interface_t *dom_intrf = (PDM_part_domain_interface_t *) malloc (sizeof(PDM_part_domain_interface_t));
  dom_intrf->n_interface       = n_interface;
  dom_intrf->n_domain          = n_domain;
  dom_intrf->n_part            = malloc(n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->n_part[i_domain] = n_part[i_domain];
  }
  dom_intrf->multidomain_intrf = multidomain_interface;
  dom_intrf->ownership         = ownership;
  dom_intrf->comm              = comm;

  dom_intrf->interface_pn_vtx        = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_vtx_ln_to_gn  = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_vtx_idx   = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_vtx       = (int         ****) malloc(n_domain * sizeof(int         ***));

  dom_intrf->interface_pn_face       = (int          ***) malloc(n_domain * sizeof(int          **));
  dom_intrf->interface_face_ln_to_gn = (PDM_g_num_t ****) malloc(n_domain * sizeof(PDM_g_num_t ***));
  dom_intrf->interface_sgn_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face      = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_ids_face_idx  = (int         ****) malloc(n_domain * sizeof(int         ***));
  dom_intrf->interface_dom_face      = (int         ****) malloc(n_domain * sizeof(int         ***));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain ) {
    dom_intrf->interface_pn_vtx       [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_vtx_ln_to_gn [i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_vtx_idx  [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_vtx      [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    dom_intrf->interface_pn_face      [i_domain] = (int          **) malloc(n_part[i_domain] * sizeof(int          *));
    dom_intrf->interface_face_ln_to_gn[i_domain] = (PDM_g_num_t ***) malloc(n_part[i_domain] * sizeof(PDM_g_num_t **));
    dom_intrf->interface_sgn_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_ids_face_idx [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));
    dom_intrf->interface_dom_face     [i_domain] = (int         ***) malloc(n_part[i_domain] * sizeof(int         **));

    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
      dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
      dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
      dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
      dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
      dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
      dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
      dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
      dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
    }
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    dom_intrf->is_result[i] = 0;
  }

  dom_intrf->translation_vect   = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_direction = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_center    = (double **) malloc(n_interface * sizeof(double *));
  dom_intrf->rotation_angle     = (double  *) malloc(n_interface * sizeof(double  ));

  for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
    dom_intrf->translation_vect  [i_interface] = NULL;
    dom_intrf->rotation_direction[i_interface] = NULL;
    dom_intrf->rotation_center   [i_interface] = NULL;
    dom_intrf->rotation_angle    [i_interface] = 0.;
  }

  dom_intrf->interface_describe_by_face = 0;
  dom_intrf->interface_describe_by_vtx  = 0;

  return dom_intrf;
}

void
PDM_part_domain_interface_set
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind,
 int                           i_domain,
 int                           i_part,
 int                          *interface_pn,
 PDM_g_num_t                 **interface_ln_to_gn,
 int                         **interface_sgn,
 int                         **interface_ids,
 int                         **interface_ids_idx,
 int                         **interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    dom_intrf->interface_pn_vtx      [i_domain][i_part] = interface_pn;
    dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part] = interface_ln_to_gn;
    dom_intrf->interface_sgn_vtx     [i_domain][i_part] = interface_sgn;
    dom_intrf->interface_ids_vtx     [i_domain][i_part] = interface_ids;
    dom_intrf->interface_ids_vtx_idx [i_domain][i_part] = interface_ids_idx;
    dom_intrf->interface_dom_vtx     [i_domain][i_part] = interface_dom;
    dom_intrf->interface_describe_by_vtx  = 1;
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    dom_intrf->interface_pn_face      [i_domain][i_part] = interface_pn;
    dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = interface_ln_to_gn;
    dom_intrf->interface_sgn_face     [i_domain][i_part] = interface_sgn;
    dom_intrf->interface_ids_face     [i_domain][i_part] = interface_ids;
    dom_intrf->interface_ids_face_idx [i_domain][i_part] = interface_ids_idx;
    dom_intrf->interface_dom_face     [i_domain][i_part] = interface_dom;
    dom_intrf->interface_describe_by_face  = 1;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

void
PDM_part_domain_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf,
 PDM_bound_type_t               interface_kind,
 int                            i_domain,
 int                            i_part,
 int                          **interface_pn,
 PDM_g_num_t                 ***interface_ln_to_gn,
 int                         ***interface_sgn,
 int                         ***interface_ids,
 int                         ***interface_ids_idx,
 int                         ***interface_dom
)
{
  assert(i_part < dom_intrf->n_part[i_domain]);

  assert (dom_intrf != NULL);
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    *interface_pn       = dom_intrf->interface_pn_vtx      [i_domain][i_part];
    *interface_ln_to_gn = dom_intrf->interface_vtx_ln_to_gn[i_domain][i_part];
    *interface_sgn      = dom_intrf->interface_sgn_vtx     [i_domain][i_part];
    *interface_ids      = dom_intrf->interface_ids_vtx     [i_domain][i_part];
    *interface_ids_idx  = dom_intrf->interface_ids_vtx_idx [i_domain][i_part];
    *interface_dom      = dom_intrf->interface_dom_vtx     [i_domain][i_part];
  }
  else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    *interface_pn       = dom_intrf->interface_pn_face      [i_domain][i_part];
    *interface_ln_to_gn = dom_intrf->interface_face_ln_to_gn[i_domain][i_part];
    *interface_sgn      = dom_intrf->interface_sgn_face     [i_domain][i_part];
    *interface_ids      = dom_intrf->interface_ids_face     [i_domain][i_part];
    *interface_ids_idx  = dom_intrf->interface_ids_face_idx [i_domain][i_part];
    *interface_dom      = dom_intrf->interface_dom_face     [i_domain][i_part];
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Kind of interface not supported\n");
  }
}

int
PDM_part_domain_interface_n_interface_get
(
 PDM_part_domain_interface_t   *dom_intrf
)
{
  return dom_intrf->n_interface;
}


int
PDM_part_domain_interface_exist_get
(
 PDM_part_domain_interface_t  *dom_intrf,
 PDM_bound_type_t              interface_kind
)
{
  if (interface_kind == PDM_BOUND_TYPE_VTX) {
    return dom_intrf->interface_describe_by_vtx;
  } else if (interface_kind == PDM_BOUND_TYPE_FACE) {
    return dom_intrf->interface_describe_by_face;
  }
  return 0;
}

void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
)
{


  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {

    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      if(dom_intrf->ownership == PDM_OWNERSHIP_KEEP) {
        if(dom_intrf->interface_pn_vtx       [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_vtx       [i_domain][i_part]);
          dom_intrf->interface_pn_vtx       [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part]);
          dom_intrf->interface_vtx_ln_to_gn [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_vtx      [i_domain][i_part]);
          dom_intrf->interface_sgn_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx      [i_domain][i_part]);
          dom_intrf->interface_ids_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_vtx_idx  [i_domain][i_part]);
          dom_intrf->interface_ids_vtx_idx  [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_vtx      [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_vtx      [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_vtx      [i_domain][i_part]);
          dom_intrf->interface_dom_vtx      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_pn_face      [i_domain][i_part] != NULL) {
          free(dom_intrf->interface_pn_face      [i_domain][i_part]);
          dom_intrf->interface_pn_face      [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_face_ln_to_gn[i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_face_ln_to_gn[i_domain][i_part]);
          dom_intrf->interface_face_ln_to_gn[i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_sgn_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_sgn_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_sgn_face     [i_domain][i_part]);
          dom_intrf->interface_sgn_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face     [i_domain][i_part]);
          dom_intrf->interface_ids_face     [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_ids_face_idx [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_ids_face_idx [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_ids_face_idx [i_domain][i_part]);
          dom_intrf->interface_ids_face_idx [i_domain][i_part] = NULL;
        };
        if(dom_intrf->interface_dom_face     [i_domain][i_part] != NULL) {
          for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface){
            free(dom_intrf->interface_dom_face     [i_domain][i_part][i_interface]);
          }
          free(dom_intrf->interface_dom_face     [i_domain][i_part]);
          dom_intrf->interface_dom_face     [i_domain][i_part] = NULL;
        };
      }
    }

    free(dom_intrf->interface_pn_vtx       [i_domain]);
    free(dom_intrf->interface_vtx_ln_to_gn [i_domain]);
    free(dom_intrf->interface_sgn_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx      [i_domain]);
    free(dom_intrf->interface_ids_vtx_idx  [i_domain]);
    free(dom_intrf->interface_dom_vtx      [i_domain]);

    free(dom_intrf->interface_pn_face      [i_domain]);
    free(dom_intrf->interface_face_ln_to_gn[i_domain]);
    free(dom_intrf->interface_sgn_face     [i_domain]);
    free(dom_intrf->interface_ids_face     [i_domain]);
    free(dom_intrf->interface_ids_face_idx [i_domain]);
    free(dom_intrf->interface_dom_face     [i_domain]);

  }
  free(dom_intrf->n_part);

  free(dom_intrf->interface_pn_vtx       );
  free(dom_intrf->interface_vtx_ln_to_gn );
  free(dom_intrf->interface_sgn_vtx      );
  free(dom_intrf->interface_ids_vtx      );
  free(dom_intrf->interface_ids_vtx_idx  );
  free(dom_intrf->interface_dom_vtx      );

  free(dom_intrf->interface_pn_face      );
  free(dom_intrf->interface_face_ln_to_gn);
  free(dom_intrf->interface_sgn_face     );
  free(dom_intrf->interface_ids_face     );
  free(dom_intrf->interface_ids_face_idx );
  free(dom_intrf->interface_dom_face     );

  for(int i_interface = 0; i_interface < dom_intrf->n_interface; ++i_interface) {
    if(dom_intrf->translation_vect[i_interface]   != NULL) {
      free(dom_intrf->translation_vect[i_interface]);
      dom_intrf->translation_vect[i_interface] = NULL;
    }
    if(dom_intrf->rotation_direction[i_interface]   != NULL) {
      free(dom_intrf->rotation_direction[i_interface]);
      dom_intrf->rotation_direction[i_interface] = NULL;
    }
    if(dom_intrf->rotation_center[i_interface]   != NULL) {
      free(dom_intrf->rotation_center[i_interface]);
      dom_intrf->rotation_center[i_interface] = NULL;
    }
  }

  free(dom_intrf->translation_vect  );
  free(dom_intrf->rotation_direction);
  free(dom_intrf->rotation_center   );
  free(dom_intrf->rotation_angle    );

  free(dom_intrf);
}




void
PDM_part_domain_interface_translation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
  const double                       *vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->translation_vect[i_interface] == NULL);

  dom_intrf->translation_vect[i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->translation_vect[i_interface][i] = vect[i];
  }

}

void
PDM_part_domain_interface_rotation_set
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
  const double                       *direction,
  const double                       *center,
  const double                        angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  assert(dom_intrf->rotation_direction[i_interface] == NULL);
  assert(dom_intrf->rotation_center   [i_interface] == NULL);

  dom_intrf->rotation_direction[i_interface] = (double *) malloc( 3 * sizeof(double));
  dom_intrf->rotation_center   [i_interface] = (double *) malloc( 3 * sizeof(double));

  for(int i = 0; i < 3; ++i) {
    dom_intrf->rotation_direction[i_interface][i] = direction[i];
    dom_intrf->rotation_center   [i_interface][i] = center   [i];
  }
  dom_intrf->rotation_angle[i_interface] = angle;
}



void
PDM_part_domain_interface_translation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
        int                           i_interface,
        double                      **vect
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->translation_vect[i_interface] != NULL){

    *vect = (double *) malloc( 3 * sizeof(double));
    double* _vect = *vect;

    for(int i = 0; i < 3; ++i) {
      _vect[i] = dom_intrf->translation_vect[i_interface][i];
    }
  } else {
    *vect = NULL;
  }
}

void
PDM_part_domain_interface_rotation_get
(
        PDM_part_domain_interface_t  *dom_intrf,
  const int                           i_interface,
        double                      **direction,
        double                      **center,
        double                       *angle
)
{
  assert(i_interface < dom_intrf->n_interface);
  if(dom_intrf->rotation_direction[i_interface] != NULL) {
    assert(dom_intrf->rotation_center   [i_interface] != NULL);

    *direction = (double *) malloc( 3 * sizeof(double));
    *center    = (double *) malloc( 3 * sizeof(double));
    double *_direction = *direction;
    double *_center    = *center   ;

    for(int i = 0; i < 3; ++i) {
      _direction[i] = dom_intrf->rotation_direction[i_interface][i];
      _center   [i] = dom_intrf->rotation_center   [i_interface][i];
    }
    *angle = dom_intrf->rotation_angle[i_interface];
  } else {
    *direction = NULL;
    *center    = NULL;
    *angle     = 0;
  }
}


void
PDM_part_domain_interface_as_graph
(
  PDM_part_domain_interface_t    *dom_intrf,
  PDM_bound_type_t                interface_kind,
  int                           **n_entity,
  PDM_g_num_t                  ***entity_ln_to_gn
)
{
  int i_rank;
  PDM_MPI_Comm_rank(dom_intrf->comm, &i_rank);

  assert(dom_intrf->n_domain == 1); // TODO --> shift of gnum AND part_id

  int* n_tot_part_by_domain = (int *) malloc( dom_intrf->n_domain * sizeof(int));
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain) {

    n_tot_part_by_domain[i_domain] = -1;
    int n_part_loc = dom_intrf->n_part[i_domain];
    PDM_MPI_Allreduce(&n_part_loc, &n_tot_part_by_domain[i_domain], 1, PDM_MPI_INT, PDM_MPI_SUM, dom_intrf->comm);
  }

  /* Si multidomain on fait un shift et tt roule */
  int n_part_loc_all_domain = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain) {
    n_part_loc_all_domain += dom_intrf->n_part[i_domain];
  }

  int n_interface = PDM_part_domain_interface_n_interface_get(dom_intrf);

  int **neighbor_n         = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_idx       = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_desc      = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_interface = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int  *n_entity_bound     = (int  *) malloc( n_part_loc_all_domain * sizeof(int  ) );

  int **neighbor_opp_n    = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_opp_idx  = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );
  int **neighbor_opp_desc = (int **) malloc( n_part_loc_all_domain * sizeof(int *) );

  /*
   * Loop over all interfaces to create distant neighbor structure
   */

  int shift_part   = 0;
  int shift_part_g = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {
    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      n_entity_bound[i_part+shift_part] = n_entity[i_domain][i_part];
      neighbor_idx  [i_part+shift_part] = (int *) malloc( (n_entity_bound[i_part+shift_part]+1) * sizeof(int) );
      neighbor_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_n   = neighbor_n    [i_part+shift_part];
      int* _neighbor_idx = neighbor_idx  [i_part+shift_part];

      neighbor_opp_idx  [i_part+shift_part] = (int *) malloc( (n_entity_bound[i_part+shift_part]+1) * sizeof(int) );
      neighbor_opp_n    [i_part+shift_part] = PDM_array_zeros_int(n_entity_bound[i_part+shift_part]);

      int* _neighbor_opp_n   = neighbor_opp_n    [i_part+shift_part];
      int* _neighbor_opp_idx = neighbor_opp_idx  [i_part+shift_part];

      int           *interface_pn       = NULL;
      PDM_g_num_t  **interface_ln_to_gn = NULL;
      int          **interface_sgn      = NULL;
      int          **interface_ids      = NULL;
      int          **interface_ids_idx  = NULL;
      int          **interface_dom      = NULL;
      PDM_part_domain_interface_get(dom_intrf,
                                    interface_kind,
                                    i_domain,
                                    i_part,
                                    &interface_pn,
                                    &interface_ln_to_gn,
                                    &interface_sgn,
                                    &interface_ids,
                                    &interface_ids_idx,
                                    &interface_dom);

      /*
       * First step : Count interface to add in distant neighbor due to connectivity betwenn domain
       */
      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {

        // log_trace("-------------------------------- i_interface = %i  -------------------------------- \n", i_interface);
        // PDM_log_trace_array_int(interface_sgn[i_interface], interface_pn[i_interface], "interface_sgn :: ");

        for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

          // Search the first in list that is in current part/proc
          // int i_proc_cur   = -1;
          // int i_part_cur   = -1;
          int i_entity_cur = -1;
          int found        = 0;
          int idx_current  = -1;
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
              // i_proc_cur   = i_proc_opp;
              // i_part_cur   = i_part_opp;
              i_entity_cur = i_entity_opp;
              idx_current  = j;
              found = 1;
            }
          }

          if(!found) {
            continue;
          }

          // Il manque une notion de direction sinon on sait pas dans quelle sens va le raccord

          assert(found == 1);

          // log_trace("i_proc_cur = %i | i_part_cur = %i | i_entity_cur = %i | sgn = %i \n", i_proc_cur, i_part_cur, i_entity_cur, interface_sgn[i_interface][idx_entity]);

          // Only add the opposite part of the graph
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            // int i_proc_opp   = interface_ids[i_interface][3*j  ];
            // int i_part_opp   = interface_ids[i_interface][3*j+1];
            // int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(idx_current != j) {
              // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
              _neighbor_n[i_entity_cur] += 1;
            } else {
              _neighbor_opp_n[i_entity_cur] += 1;
            }
          }
        }
      }

      /* Compute index */
      _neighbor_idx[0] = 0;
      _neighbor_opp_idx[0] = 0;
      for(int i_entity = 0; i_entity < n_entity_bound[i_part+shift_part]; ++i_entity) {
        _neighbor_idx    [i_entity+1] = _neighbor_idx    [i_entity] + _neighbor_n    [i_entity];
        _neighbor_opp_idx[i_entity+1] = _neighbor_opp_idx[i_entity] + _neighbor_opp_n[i_entity];
        _neighbor_n    [i_entity] = 0;
        _neighbor_opp_n[i_entity] = 0;
      }

      neighbor_desc     [i_part+shift_part] = (int *) malloc( 3 * _neighbor_idx    [n_entity_bound[i_part+shift_part]] * sizeof(int) );
      neighbor_opp_desc [i_part+shift_part] = (int *) malloc( 4 * _neighbor_opp_idx[n_entity_bound[i_part+shift_part]] * sizeof(int) );
      neighbor_interface[i_part+shift_part] = (int *) malloc(     _neighbor_idx[n_entity_bound[i_part+shift_part]] * sizeof(int) );
      int* _neighbor_desc      = neighbor_desc     [i_part+shift_part];
      int* _neighbor_opp_desc  = neighbor_opp_desc [i_part+shift_part];
      int* _neighbor_interface = neighbor_interface[i_part+shift_part];

      for(int i_interface = 0; i_interface < n_interface; ++i_interface) {
        for(int idx_entity = 0; idx_entity < interface_pn[i_interface]; ++idx_entity) {

          // Search the first in list that is in current part/proc
          // int i_proc_cur   = -1;
          // int i_part_cur   = -1;
          int i_entity_cur = -1;
          int found        = 0;
          int idx_current  = -1;
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(i_proc_opp == i_rank && i_part_opp == i_part && found == 0) {
              // i_proc_cur   = i_proc_opp;
              // i_part_cur   = i_part_opp;
              i_entity_cur = i_entity_opp;
              idx_current  = j;
              found = 1;
            }
          }

          if(!found) {
            continue;
          }

          assert(found == 1);

          // Only add the opposite part of the graph
          for(int j = interface_ids_idx[i_interface][idx_entity]; j < interface_ids_idx[i_interface][idx_entity+1]; ++j) {
            int i_proc_opp   = interface_ids[i_interface][3*j  ];
            int i_part_opp   = interface_ids[i_interface][3*j+1];
            int i_entity_opp = interface_ids[i_interface][3*j+2];

            if(idx_current != j) {
              // log_trace("\t i_proc_opp = %i | i_part_opp = %i | i_entity_opp = %i \n", i_proc_opp, i_part_opp, i_entity_opp);
              int idx_write = _neighbor_idx[i_entity_cur] + _neighbor_n[i_entity_cur]++;
              _neighbor_desc[3*idx_write  ] = i_proc_opp;              // i_proc_opp;
              _neighbor_desc[3*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
              _neighbor_desc[3*idx_write+2] = i_entity_opp;            // i_entity_opp
              _neighbor_interface[idx_write] = (i_interface+1) * interface_sgn[i_interface][idx_entity];
              // _neighbor_interface[idx_write] = (i_interface+1);
            } else {
              int idx_write = _neighbor_opp_idx[i_entity_cur] + _neighbor_opp_n[i_entity_cur]++;
              _neighbor_opp_desc[4*idx_write  ] = i_proc_opp;              // i_proc_opp;
              _neighbor_opp_desc[4*idx_write+1] = i_part_opp+shift_part_g; // i_part_opp
              _neighbor_opp_desc[4*idx_write+2] = i_entity_opp;            // i_entity_opp
              _neighbor_opp_desc[4*idx_write+3] = - (i_interface+1) * interface_sgn[i_interface][idx_entity];
            }
          }
        }
      }


      if(1 == 1) {
        printf(" ------------------------------------- _neighbor graph : \n");
        for(int i = 0; i < n_entity_bound[i_part+shift_part]; ++i) {
          printf("i_elmt [%i] = ", i);
          for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
              printf(" (%i, %i, %i, %i) ", _neighbor_desc[3*idx], _neighbor_desc[3*idx+1], _neighbor_desc[3*idx+2], _neighbor_interface[idx]);
          }
          printf("\n");
        }
        printf(" ------------------------------------- _neighbor graph END \n");
      }

    }
    shift_part   += dom_intrf->n_part[i_domain];
    shift_part_g += n_tot_part_by_domain[i_domain];
  }



  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(dom_intrf->comm,
                                                           n_part_loc_all_domain,
                                                           n_entity_bound,
                                                           neighbor_idx,
                                                           neighbor_desc);

  int **next_neighbor_opp_n = NULL;
  int **next_neighbor_opp = NULL;
  PDM_distant_neighbor_exch(dn,
                            4 * sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            neighbor_opp_n,
                  (void **) neighbor_opp_desc,
                           &next_neighbor_opp_n,
                 (void ***)&next_neighbor_opp);


  shift_part   = 0;
  shift_part_g = 0;
  for(int i_domain = 0; i_domain < dom_intrf->n_domain; ++i_domain ) {
    for(int i_part = 0; i_part < dom_intrf->n_part[i_domain]; ++i_part) {

      int   n_elmt = n_entity_bound[i_part+shift_part];
      int  *_neighbor_idx   = neighbor_idx  [i_part+shift_part];
      int  *_next_neighbor_opp   = next_neighbor_opp  [i_part+shift_part];
      int  *_next_neighbor_opp_n = next_neighbor_opp_n[i_part+shift_part];

      log_trace("n_elmt = %i \n", n_elmt);
      int n_data = 0;
      int idx_read = 0;
      for(int i = 0; i < n_elmt; ++i) {

        printf("i_elmt [%i] = ", i);
        for(int idx = _neighbor_idx[i]; idx < _neighbor_idx[i+1]; ++idx) {
          for(int k = 0; k < _next_neighbor_opp_n[idx]; ++k) {
            printf(" (%i, %i, %i, %i) ",
                   _next_neighbor_opp[4*idx_read],
                   _next_neighbor_opp[4*idx_read+1],
                   _next_neighbor_opp[4*idx_read+2],
                   _next_neighbor_opp[4*idx_read+3]);
            idx_read += 1;
          }
          n_data += _next_neighbor_opp_n[idx];
        }
        printf("\n");
      }

      PDM_log_trace_array_int(_next_neighbor_opp, 4 * n_data, "next_neighbor_opp :: ");

    }
    shift_part   += dom_intrf->n_part[i_domain];
    shift_part_g += n_tot_part_by_domain[i_domain];
  }

  PDM_distant_neighbor_free(dn);



  free(n_tot_part_by_domain);

}
