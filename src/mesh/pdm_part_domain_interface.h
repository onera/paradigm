#ifndef __PDM_PART_DOMAIN_INTERFACE_H__
#define __PDM_PART_DOMAIN_INTERFACE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_domain_interface.h"

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
 * Type
 *============================================================================*/


typedef struct _pdm_part_domain_interface_t PDM_part_domain_interface_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
);


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
);




void
PDM_part_domain_interface_free
(
 PDM_part_domain_interface_t  *dom_intrf
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_DOMAIN_INTERFACE_H__ */
