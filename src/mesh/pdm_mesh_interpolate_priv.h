#ifndef __PDM_MESH_INTERPOLATE_PRIV_H__
#define __PDM_MESH_INTERPOLATE_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

// #include "pdm_multipart.h"
#include "pdm_part_priv.h"
#include "pdm_part_domain_interface.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_mesh_interpolate_t
 * \brief  Distributed cube
 *
 * _pdm_mesh_interpolate_t define a distributed mesh of a cube
 *
 */

struct _pdm_mesh_interpolate_t {
  PDM_MPI_Comm      comm;            /*!< MPI communicator                          */

  int               n_domain;
  int              *n_group;
  int              *n_part;
  int              *n_part_idx;
  int              *n_part_g_idx;

  PDM_part_domain_interface_t  *pdi;
  int                           n_interface;

  _part_t  **parts;

  /* Graphe communication inside domain */
  int ****entity_part_bound_proc_idx;
  int ****entity_part_bound_part_idx;
  int ****entity_part_bound;
  int    *graph_comm_is_defined;

  int ****group_entity_idx; // (i_kind, i_domain, i_part)
  int ****group_entity;
  int    *group_is_defined;

};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_INTERPOLATE_PRIV_H__ */
