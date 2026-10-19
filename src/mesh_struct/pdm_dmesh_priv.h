#ifndef __PDM_DMESH_PRIV_H__
#define __PDM_DMESH_PRIV_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_dmesh_t
 * \brief  Define a distributed mesh. Arrays are shared: this structure does not
 *         holds the data.
 *
 */

struct _pdm_dmesh_t
{
  PDM_MPI_Comm comm;

  int         dn_cell;                          /*!< Number of distributed cells         */
  int         dn_face;                          /*!< Number of distributed faces         */
  int         dn_edge;                          /*!< Number of distributed edges         */
  int         dn_vtx;                           /*!< Number of distributed vertices      */

  PDM_g_num_t n_g_cell;                         /*!< Number of distributed cells         */
  PDM_g_num_t n_g_face;                         /*!< Number of distributed faces         */
  PDM_g_num_t n_g_edge;                         /*!< Number of distributed edges         */
  PDM_g_num_t n_g_vtx;                          /*!< Number of distributed vertices      */

  PDM_g_num_t *cell_distrib;
  PDM_g_num_t *face_distrib;
  PDM_g_num_t *edge_distrib;
  PDM_g_num_t *vtx_distrib;

  double      *_dvtx_coord;                     /*!< Coordinates of ditributed vertices
                                                  (size = 3 * dn_vtx)                    */
  PDM_ownership_t owner_vtx_coord;

  int           n_group_bnd[PDM_BOUND_TYPE_MAX]; /*!< Number of group by elememnt type            */
  PDM_g_num_t **dconnectivity;                   /* Array of connectivty (size = PDM_CONNECTIVITY_TYPE_MAX) */
  int         **dconnectivity_idx;               /* Array of connectivty_idx if any (size = PDM_CONNECTIVITY_TYPE_MAX) */

  PDM_ownership_t owner_connectivity[PDM_CONNECTIVITY_TYPE_MAX];

  PDM_g_num_t **dbound;                   /* Array of connectivty (size = PDM_CONNECTIVITY_TYPE_MAX) */
  int         **dbound_idx;               /* Array of connectivty_idx if any (size = PDM_CONNECTIVITY_TYPE_MAX) */

  PDM_ownership_t owner_bound[PDM_BOUND_TYPE_MAX];

  PDM_bool_t   is_computed_g_extents;
  double       g_extents[6];

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DMESH_PRIV_H__ */
