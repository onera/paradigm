/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_interpolate_from_mesh_location_priv.h"
#include "pdm_interpolate_from_mesh_location.h"
#include "pdm_timer.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute interpolation from mesh_location information
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

PDM_interpolate_from_mesh_location_t*
PDM_interpolate_from_mesh_location_create
(
 const int                    n_part_src,
 const int                    n_cloud_target,
       PDM_interpolate_kind_t interp_kind,
 const PDM_MPI_Comm           comm
)
{
  PDM_interpolate_from_mesh_location_t *interp_from_ml = (PDM_interpolate_from_mesh_location_t *) malloc (sizeof(PDM_interpolate_from_mesh_location_t));

  interp_from_ml->comm = comm;

  PDM_MPI_Comm_size (comm, &interp_from_ml->n_rank);
  PDM_MPI_Comm_rank (comm, &interp_from_ml->i_rank);

  interp_from_ml->n_part_src     = n_part_src;
  interp_from_ml->n_cloud_target = n_cloud_target;

  interp_from_ml->points_in_elements = malloc (sizeof(_points_in_element_t) * interp_from_ml->n_cloud_target);

  for (int icloud = 0; icloud < interp_from_ml->n_cloud_target; icloud++) {

    _points_in_element_t *pcloud = interp_from_ml->points_in_elements + icloud;

    pcloud->pts_inside_idx   = (int         **) malloc( n_part_src * sizeof(int         *));
    pcloud->gnum             = (PDM_g_num_t **) malloc( n_part_src * sizeof(PDM_g_num_t *));
    pcloud->uvw              = (double      **) malloc( n_part_src * sizeof(double      *));
    pcloud->coords           = (double      **) malloc( n_part_src * sizeof(double      *));
    pcloud->projected_coords = (double      **) malloc( n_part_src * sizeof(double      *));
    pcloud->weights_idx      = (int         **) malloc( n_part_src * sizeof(int         *));
    pcloud->weights          = (double      **) malloc( n_part_src * sizeof(double      *));
    pcloud->dist2            = (double      **) malloc( n_part_src * sizeof(double      *));

  }

  PDM_UNUSED(interp_kind);

  return interp_from_ml;
}

void
PDM_mesh_location_part_free
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
)
{

  for (int icloud = 0; icloud < interp_from_ml->n_cloud_target; icloud++) {

    _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + icloud;

    for (int i_part = 0; i_part < _points_in_elements->n_part; ++i_part) {

      free (_points_in_elements->pts_inside_idx[i_part]);
      free (_points_in_elements->gnum[i_part]);
      free (_points_in_elements->uvw[i_part]);
      free (_points_in_elements->coords[i_part]);
      free (_points_in_elements->projected_coords[i_part]);
      free (_points_in_elements->weights_idx[i_part]);
      free (_points_in_elements->weights[i_part]);
      free (_points_in_elements->dist2[i_part]);

    }
  }

  free(interp_from_ml->points_in_elements);
  free(interp_from_ml);
}


void
PDM_mesh_location_part_compute
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml
)
{
  PDM_UNUSED(interp_from_ml);

  /*
   * On prepare les échanges pour faire le part_to_part
   */

}

void
PDM_mesh_location_part_exch
(
 PDM_interpolate_from_mesh_location_t   *interp_from_ml,
 size_t                                  s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);


}


void
PDM_mesh_location_part_send
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 size_t                                 s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);
}



void
PDM_mesh_location_part_recv
(
 PDM_interpolate_from_mesh_location_t  *interp_from_ml,
 size_t                                 s_data,
 void                                  **part_data_in,
 void                                 ***cloud_data_out
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(s_data);
  PDM_UNUSED(part_data_in);
  PDM_UNUSED(cloud_data_out);
}

/**
 *
 * \brief Set a part of a mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of cells
 * \param [in]   cell_face_idx Index in the cell -> face connectivity
 * \param [in]   cell_face     cell -> face connectivity
 * \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */
void
PDM_interpolate_from_mesh_location_part_set
(
       PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                                   i_part,
 const int                                   n_cell,
 const int                                  *cell_face_idx,
 const int                                  *cell_face,
 const PDM_g_num_t                          *cell_ln_to_gn,
 const int                                   n_face,
 const int                                  *face_vtx_idx,
 const int                                  *face_vtx,
 const PDM_g_num_t                          *face_ln_to_gn,
 const int                                   n_vtx,
 const double                               *coords,
 const PDM_g_num_t                          *vtx_ln_to_gn
)
{
  PDM_UNUSED(interp_from_ml);
  PDM_UNUSED(i_part);
  PDM_UNUSED(n_cell);
  PDM_UNUSED(cell_face_idx);
  PDM_UNUSED(cell_face);
  PDM_UNUSED(cell_ln_to_gn);
  PDM_UNUSED(n_face);
  PDM_UNUSED(face_vtx_idx);
  PDM_UNUSED(face_vtx);
  PDM_UNUSED(face_ln_to_gn);
  PDM_UNUSED(n_vtx);
  PDM_UNUSED(coords);
  PDM_UNUSED(vtx_ln_to_gn);

}


/**
 *
 * \brief Set point list located in elements
 *
 * \param [in]   id                      Identifier
 * \param [in]   i_part                  Index of partition of the mesh
 * \param [in]   i_point_cloud           Index of cloud
 * \param [in]   elt_pts_inside_idx      Points index (size = n_elt + 1)
 * \param [in]   points_gnum             Points global number
 * \param [in]   points_coords           Points coordinates
 * \param [in]   points_uvw              Points parametric coordinates in elements
 * \param [in]   points_weights_idx      Interpolation weights index (size = elt_pts_inside_idx[n_elt] + 1)
 * \param [in]   points_weights          Interpolation weights
 * \param [in]   points_dist2            Distance element-points (dist < 0 if the point is inside)
 * \param [in]   points_projected_coords Point projection on element if the point is outside
 *
 */
void
PDM_interpolate_from_mesh_location_points_in_elt_set
(
 PDM_interpolate_from_mesh_location_t *interp_from_ml,
 const int                             i_part,
 const int                             i_point_cloud,
 int                                  *elt_pts_inside_idx,
 PDM_g_num_t                          *points_gnum,
 double                               *points_coords,
 double                               *points_uvw,
 int                                  *points_weights_idx,
 double                               *points_weights,
 double                               *points_dist2,
 double                               *points_projected_coords
)
{
  assert (i_point_cloud < interp_from_ml->n_cloud_target);

  assert (interp_from_ml->points_in_elements != NULL);
  assert (i_part < interp_from_ml->points_in_elements->n_part);

  _points_in_element_t *_points_in_elements = interp_from_ml->points_in_elements + i_point_cloud;

  _points_in_elements->pts_inside_idx[i_part]   = elt_pts_inside_idx;
  _points_in_elements->gnum[i_part]             = points_gnum;
  _points_in_elements->coords[i_part]           = points_coords;
  _points_in_elements->uvw[i_part]              = points_uvw;
  _points_in_elements->weights_idx[i_part]      = points_weights_idx;
  _points_in_elements->weights[i_part]          = points_weights;
  _points_in_elements->dist2[i_part]            = points_dist2;
  _points_in_elements->projected_coords[i_part] = points_projected_coords;

}



#ifdef  __cplusplus
}
#endif
