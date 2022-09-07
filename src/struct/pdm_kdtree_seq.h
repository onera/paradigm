#ifndef PDM_KDTREE_SEQ_H
#define PDM_KDTREE_SEQ_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_kdtree_seq_t PDM_kdtree_seq_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create an kdtree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_kdtree_seq object
 */

PDM_kdtree_seq_t *
PDM_kdtree_seq_create
(
 const int    n_point_cloud,
 const int    depth_max,
 const int    points_in_leaf_max,
 const double tolerance
);


/**
 *
 * \brief Free an kdtree structure
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_free
(
 PDM_kdtree_seq_t *kdtree
);



/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 *
 */


void
PDM_kdtree_seq_point_cloud_set
(
 PDM_kdtree_seq_t *kdtree,
 const int         i_point_cloud,
 const int         n_points,
 const double     *coords
);



/**
 *
 * \brief Build kdtree
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_build
(
 PDM_kdtree_seq_t *kdtree
);



/**
 *
 * \brief Write kdtree nodes in a VTK file
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [in]   filename               Output file name
 *
 */

void PDM_kdtree_seq_write_nodes
(
 PDM_kdtree_seq_t *kdtree,
 const char       *filename
 );




void
PDM_kdtree_seq_points_inside_ball
(
 const PDM_kdtree_seq_t  *kdtree,
 const int                n_pts,
 double                  *pts_coord,
 double                  *ball_radius2,
 int                    **pts_inside_ball_idx,
 int                    **pts_inside_ball_l_num,
 double                 **pts_inside_ball_dist2
 );



void
PDM_kdtree_seq_extract_extent
(
  PDM_kdtree_seq_t  *kdtree,
  int                root_id,
  int                n_depth,
  int               *n_box,
  double           **box_extents
);

#ifdef  __cplusplus
}
#endif

#endif  /* PDM_KDTREE_SEQ_H */

