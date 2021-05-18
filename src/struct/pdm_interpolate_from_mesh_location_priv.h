#ifndef __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__
#define __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__

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

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _points_in_element_t
 * \brief
 *
 */

typedef struct {

  int           n_part;
  int          *n_elts; // Aredescendre dans la structure parente
  int         **pts_inside_idx;
  PDM_g_num_t **gnum;
  double      **coords;
  double      **uvw;
  double      **projected_coords;
  int         **weights_idx;
  double      **weights;
  double      **dist2;

} _points_in_element_t;

/**
 * \struct _pdm_interpolate_from_mesh_location_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

struct _pdm_interpolate_from_mesh_location_t {
  PDM_MPI_Comm   comm;                /*!< MSG communicator */
  int            n_rank;              /*!< Communicator size */
  int            i_rank;              /*!< Current rank in comm */

  int            n_part_src;          /*!< Number of partitions */
  int           *n_cell;              /*!< Number of elements for each partition */
  int           *n_face;              /*!< Number of elements for each partition */
  int           *n_vtx;               /*!< Number of elements for each partition */

  int          **cell_face_idx;
  int          **cell_face;
  PDM_g_num_t  **cell_ln_to_gn;

  int          **face_vtx_idx;
  int          **face_vtx;
  PDM_g_num_t  **face_ln_to_gn;

  PDM_g_num_t  **vtx_ln_to_gn;
  double       **coords;

  int            n_cloud_target;

  _points_in_element_t *points_in_elements; /*!< System CPU time */


} _pdm_interpolate_from_mesh_location_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__ */
