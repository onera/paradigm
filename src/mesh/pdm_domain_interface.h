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
 *  - Some global data about interfaces : n_interface is the number of interfaces, (each
 *    interface is counted twice), group_join_to_join_opp is an array indicating the number
 *    of the opposite interface, group_join_to_zone_cur and group_join_to_zone_opp gives the
 *    id of domain and opposite domain for each interface
 *  - The faces belonging to interfaces in a distributed flat jagged array : dface_join_idx is the number of
 *    faces in each interface, dface_join is the global id of theese faces in the domain numbering and
 *    dface_join_opp is the id of the matching faces in the opposite domain numbering
 *  - Some mesh data : n_zone is the number of domain, dn_vtx, dn_face and dface_vtx_idx are standard
 *    topological data for each domain
 * Ouput are the array dvtx_group_idx, dvtx_group, dvtx_group_opp who have the same meaning than
 * dface_join_idx, dface_join and dface_join_opp but for vertices.
 * Output arrays are allocated.
*/
void PDM_domain_interface_face_to_vertex
(
 int            n_interface,             /* Total number of interfaces */
 int           *group_join_to_join_opp,  /* Link between interfaces (size=n_interface) */
 int           *group_join_to_zone_cur,  /* Id of domain for each interface (size=n_interface) */
 int           *group_join_to_zone_opp,  /* Id of opposite domain for each interface (size=n_interface) */
 int           *dface_join_idx,          /* Idx array of all faces belong to interfaces (size=n_interface) */
 PDM_g_num_t   *dface_join,              /* Faces belong to interfaces (size=dface_join_idx[n_interface]) */
 PDM_g_num_t   *dface_join_opp,          /* Opposite face id of the faces belonging to interface (same size) */
 int            n_zone,                  /* Number of zones */
 PDM_g_num_t   *dn_vtx,                  /* Number of vertex in each zone (distributed) */
 PDM_g_num_t   *dn_face,                 /* Number of face in each zone (distributed) */
 int          **dface_vtx_idx,           /* Face->vertex connectivity for each domain */
 PDM_g_num_t  **dface_vtx,
 int          **dvtx_group_idx,          /* [OUT] Idx array of all vertex belonging to interfaces */
 PDM_g_num_t  **dvtx_group,              /* [OUT] Vertices belong to interfaces (size=dvtx_group_idx[n_interface]) */
 PDM_g_num_t  **dvtx_group_opp,          /* [OUT] Opposite vertex id for theses vertices (same size) */
 PDM_MPI_Comm   comm                     /* Mpi comunicator */
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DOMAIN_INTERFACE_H__ */
