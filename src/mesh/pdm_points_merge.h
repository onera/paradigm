/*
 * File:   pdm_points_merge.h
 * Author: equemera
 *
 * Created on November 16, 2017, 3:36 PM
 */

#ifndef PDM_POINTS_MERGE_H
#define	PDM_POINTS_MERGE_H

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

typedef struct _pdm_points_merge_t PDM_points_merge_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a points merge structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier
 */

PDM_points_merge_t*
PDM_points_merge_create
(
 const int             n_point_cloud,
 const double          tolerance,
 const PDM_MPI_Comm    comm,
 const PDM_ownership_t owner
);


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_points_merge_free
(
 PDM_points_merge_t* pm
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id             Identifier
 * \param [in]   i_point_cloud  Index of point cloud
 * \param [in]   n_points       Number of points
 * \param [in]   coords         Point coordinates
 * \param [in]   char_length    Characteristic length (or NULL)
 *
 */

void
PDM_points_merge_cloud_set
(
       PDM_points_merge_t *pm,
 const int                 i_point_cloud,
 const int                 n_points,
 const double             *coords,
 const double             *char_length
);


/**
 *
 * \brief Process merge points
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_points_merge_process
(
 PDM_points_merge_t *pm
);


/**
 *
 * \brief Get candidates to merge for each point
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Current cloud
 * \param [out]  candidates_idx  Indexes of candidate for each current cloud point
 *                               (size = number of points in the current cloud + 1)
 * \param [out]  candidates_desc Candidates description (process,
 *                                                       cloud in the process,
 *                                                       point in the cloud)
 *
 */

void
PDM_points_merge_candidates_get
(
       PDM_points_merge_t  *pm,
 const int                  i_point_cloud,
       int                **candidates_idx,
       int                **candidates_desc
);

/**
 *
 * \brief Get size of the resulting array
 *
 * \param [in]   id                Identifier
 * \param [in]   i_point_cloud     Current cloud
 * \param [out]  n_point_cloud     Number of points in the current cloud
 * \param [out]  n_candidates_desc Size of candidates_desc = candidates_idx[n_point_cloud+1]
 *
 */
void
PDM_points_merge_candidates_size_get
(
       PDM_points_merge_t *pm,
 const int                 i_point_cloud,
       int                *n_point_cloud,
       int                *n_candidates_desc
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_POINTS_MERGE_H */

