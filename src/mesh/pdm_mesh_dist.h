#ifndef PDM_MESH_DIST_H
#define PDM_MESH_DIST_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

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

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_create
(
 const int mesh_nodal_id,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_dist_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number 
 *
 */

void
PDM_mesh_dist_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
);


/**
 *
 * \brief Set a point cloud with initial distance
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   initial_dist    Initial distance  
 * \param [in]   coords          Point coordinates
 *
 */

/* void */
/* PDM_mesh_dist_cloud_with_initial_set */
/* ( */
/*  const int          id, */
/*  const int          i_point_cloud, */
/*  const int          i_part, */
/*  const int          n_points, */
/*  const double      *initial_dist, */
/*  const double      *coords */
/* ); */


/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

/* void */
/* PDM_mesh_dist_normal_set */
/* ( */
/*  const int          id, */
/*  const int          i_part, */
/*  const double      *normal */
/* ); */

/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

/* void */
/* PDM_mesh_dist_center_set */
/* ( */
/*  const int          id, */
/*  const int          i_part, */
/*  const double      *center */
/* ); */


/**
 *
 * \brief Process merge points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_dist_process
(
 const int id
);


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                Identifier
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   i_part            Index of partition of the cloud
 * \param [out]  distance          Distance
 * \param [out]  projected         Projected point coordinates
 * \param [out]  closest_elt_rank  Closest element rank
 * \param [out]  closest_elt_part  Closest element partition
 * \param [out]  closest_elt_l_num Local number of the closest element
 * \param [out]  closest_elt_g_num Global number of the closest element
 *
 */

void
PDM_mesh_dist_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       double      **distance,
       double      **projected,
       int         **closest_elt_rank,
       int         **closest_elt_part,
       int         **closest_elt_lnum,
       PDM_g_num_t **closest_elt_gnum
 );




/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed. Otherwise, results are kept. 
 *
 * \return     Identifier
 */

void
PDM_mesh_dist_free
(
 const int id,
 const int partial
 );


#ifdef	__cplusplus
}
#endif

#endif // PDM_MESH_DIST_H
