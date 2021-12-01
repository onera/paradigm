#ifndef __PDM_DOMAIN_INTERFACE_H__
#define __PDM_DOMAIN_INTERFACE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * From some connected domains, generated the connectivity information for the
 * vertices from the connectivity information of the faces.
 *
 * The input face-based interface connectivity should be specified through :
 *  - Some global data about interfaces : n_interface is the number of interfaces
 *  - Interface description : for each interface I, interfaces_size[I] is the number of face pairs
 *    in the interface; interface_face_ids[I] is the list of (face, opposite face) in the corresponding
 *    domains numbering (size = 2*interfaces_size[I]); interface_domains_ids is the list of (domain, opposite
 *    domain) of the faces (size = 2*interfaces_size[I])
 *  - Some mesh data : n_zone is the number of domain, dn_vtx, dn_face and dface_vtx_idx are standard
 *    topological data for each domain
 * Ouput are the array vtx_interface_size, interface_vtx_ids, interface_vtx_dom_ids who have the same meaning than
 * interface_size, interface_face_ids, interface_domains_ids.
 * Output arrays must be allocated at size n_interface ; inner level will be allocated within the function
*/
void PDM_domain_interface_face_to_vertex
(
 int            n_interface,             /* Total number of interfaces */
 int           *interfaces_size,         /* Number of face pairs in each interface */
 PDM_g_num_t  **interface_face_ids,      /* For each interface, list of pairs face,face_opp */
 int          **interface_domains_ids,   /* For each interface, list of domains dom,dom_opp */
 int            n_zone,                  /* Number of zones */
 PDM_g_num_t   *dn_vtx,                  /* Number of vertex in each zone (distributed) */
 PDM_g_num_t   *dn_face,                 /* Number of face in each zone (distributed) */
 int          **dface_vtx_idx,           /* Face->vertex connectivity for each domain */
 PDM_g_num_t  **dface_vtx,
 int           *vtx_interface_size,      /* [OUT] Number of vtx pairs in each interface */
 PDM_g_num_t  **interface_vtx_ids,       /* [OUT] For each interface, list of pairs vtx,vtx_opp */
 int          **interface_vtx_dom_ids,   /* [OUT] For each interface, list of domains dom,dom_opp */
 PDM_MPI_Comm   comm                     /* Mpi comunicator */
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOMAIN_INTERFACE_H__ */
