/*
 * \file
 */

#ifndef PDM_DIST_CLOUD_SURF_H
#define PDM_DIST_CLOUD_SURF_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_surf_mesh.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_dist_cloud_surf_t PDM_dist_cloud_surf_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nature    Nature of the mesh
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 * \param [in]   owner          Ownership of \ref PDM_dist_cloud_surf_t
 *
 * \return     Pointer to \ref PDM_dist_cloud_surf object
 */


PDM_dist_cloud_surf_t*
PDM_dist_cloud_surf_create
(
 const PDM_mesh_nature_t mesh_nature,
 const int               n_point_cloud,
 const PDM_MPI_Comm      comm,
 const PDM_ownership_t   owner
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_dist_cloud_surf_n_part_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    n_part
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_dist_cloud_surf_cloud_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
 const int                    n_points,
       double                *coords,
       PDM_g_num_t           *gnum
);



/**
 *
 * \brief Set the mesh nodal
 *
 * \param [in]   dist           Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   mesh_nodal_id  Mesh nodal Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_nodal_mesh_set
(
 PDM_dist_cloud_surf_t *dist,
 PDM_part_mesh_nodal_t *mesh_nodal
 // PDM_Mesh_nodal_t      *mesh_nodal
);

/**
 *
 * \brief Map a surface mesh
 *
 * \param [in]   dist       Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   surf_mesh  Surface mesh pointer
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_map
(
 PDM_dist_cloud_surf_t *dist,
 PDM_surf_mesh_t       *surf_mesh
);


/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   dist           Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_dist_cloud_surf_surf_mesh_global_data_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    n_part
);


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   dist          Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_part        Partition to define
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
PDM_dist_cloud_surf_surf_mesh_part_set
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_part,
 const int                    n_face,
 const int                   *face_vtx_idx,
 const int                   *face_vtx,
 const PDM_g_num_t           *face_ln_to_gn,
 const int                    n_vtx,
 const double                *coords,
 const PDM_g_num_t           *vtx_ln_to_gn
);


/**
 *
 * \brief Compute distance
 *
 * \param [in]   dist  Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_compute
(
 PDM_dist_cloud_surf_t *dist
);


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   dist                  Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud         Current cloud
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  closest_elt_distance  Distance
 * \param [out]  closest_elt_projected Projected point coordinates
 * \param [out]  closest_elt_gnum      Global number of the closest element
 *
 */

void
PDM_dist_cloud_surf_get
(
       PDM_dist_cloud_surf_t  *dist,
 const int                     i_point_cloud,
 const int                     i_part,
       double                **distance,
       double                **projected,
       PDM_g_num_t           **closest_elt_gnum
);


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_dist_cloud_surf_free
(
 PDM_dist_cloud_surf_t  *dist
);


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  dist     Pointer to \ref PDM_dist_cloud_surf object
 *
 */

void
PDM_dist_cloud_surf_dump_times
(
 PDM_dist_cloud_surf_t  *dist
);



/**
 *
 * \brief Get the dimension of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  n_points        Number of points
 *
 */

void
PDM_dist_cloud_surf_cloud_dim_get
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud,
 const int                    i_part,
       int                   *n_points
);


/**
 *
 * \brief Get the dimension of a point cloud
 *
 * \param [in]   dist            Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [out]  n_part          Number of partition
 *
 */

int
PDM_dist_cloud_surf_cloud_n_part_get
(
       PDM_dist_cloud_surf_t *dist,
 const int                    i_point_cloud
);



/**
 *
 * \brief Distribute data from surf to cloud
 *
 * \param [in]   dist              Pointer to \ref PDM_dist_cloud_surf object
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   stride            Stride
 * \param [in]   surf_data         Data over the surface
 * \param [out]  cloud_data        Data over the cloud
 *
 */

void
PDM_dist_cloud_surf_distri_data
(
       PDM_dist_cloud_surf_t  *dist,
 const int                     i_point_cloud,
 const int                     stride,
 const void                  **surf_data,
       void                 ***cloud_data
);


#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_DIST_H
