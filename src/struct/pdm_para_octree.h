#ifndef PDM_PARA_OCTREE_H
#define	PDM_PARA_OCTREE_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
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

/**
 * \enum pdm_para_octree_child_t
 * \brief Names of 8 children of a node
 *
 */
#ifndef PDM_OCTREE_H
typedef enum {
  PDM_BOTTOM,
  PDM_UP,
  PDM_SOUTH,
  PDM_NORTH,
  PDM_WEST,
  PDM_EAST,
  PDM_N_DIRECTION,
} PDM_para_octree_direction_t;
#endif
/**
 * \enum PDM_para_octree_child_t
 * \brief Names of 8 children of a node
 *
 * If the type is BOX_TREE_NODE, the ordering of children is defined as follows,
 *  using notation B: bottom, U: up, E: east, W: west, S: south,  N: north.
 *
 *  octant:   0: BSW, 1: BSE, 2: BNW, 3: BNE, 4: USW, 5: USE, 6: UNW, 7: UNE
 *  quadrant: 0:  SW, 1:  SE, 2:  NW, 3:  NE
 *  segment:  0:   W, 1:   E
 *
 */

typedef enum {
  PDM_BSW,
  PDM_BSE,
  PDM_BNW,
  PDM_BNE,
  PDM_USW,
  PDM_USE,
  PDM_UNW,
  PDM_UNE,
} PDM_para_octree_child_t;

/*  */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier
 */

int
PDM_para_octree_create
(
 const int n_point_cloud,
 const int depth_max,
 const int points_in_leaf_max,
 const int build_leaf_neighbours,
 const PDM_MPI_Comm comm
);

/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free
(
 const int          id
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 * \param [in]   g_num              Point global number or NULL
 *
 */


void
PDM_para_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords,
 const PDM_g_num_t *g_num
);


/**
 *
 * \brief Build octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_build
(
 const int          id,
 double            *global_extents
);


/**
 *
 * \brief Dump octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_dump
(
 const int          id
);


/**
 *
 * \brief Get extents
 *
 * \param [in]   id                 Identifier
 *
 * \return     Extents
 *
 */

double *
PDM_para_octree_extents_get
(
 const int          id
);


/**
 *
 * \brief Dump octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_dump
(
 const int          id
 );


/**
 *
 * Look for closest points stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_closest_points       Number of closest points to find
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest points in octree global number
 * \param [out]  closest_octree_pt_dist Closest points in octree distance
 *
 */

void
PDM_para_octree_closest_points
(
const int    id,
const int    n_closest_points,
const int    n_pts,
double      *pts,
PDM_g_num_t *pts_g_num,
PDM_g_num_t *closest_octree_pt_g_num,
double      *closest_octree_pt_dist2
);


/**
 *
 * Look for single closest point stored inside an octree
 *
 * \param [in]   id                     Identifier
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [in]   pts_g_num              Point global numbers
 * \param [out]  closest_octree_pt_id   Closest points in octree global number
 * \param [out]  closest_octree_pt_dist Closest points in octree distance
 *
 */

void
PDM_para_octree_single_closest_point
(
const int    id,
const int    n_pts,
double      *pts,
PDM_g_num_t *pts_g_num,
PDM_g_num_t *closest_octree_pt_g_num,
double      *closest_octree_pt_dist2
);

/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_para_octree_dump_times
(
 const int id
 );


/**
 *
 * Get points located inside a set of boxes
 *
 * \param [in]   id                     Octree identifier
 * \param [in]   n_boxes                Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [in]   box_g_num              Global numbers of boxes
 * \param [out]  pts_in_box_idx         Index of points located in boxes
 * \param [out]  pts_in_box_g_num       Global numbers of points located in boxes
 * \param [out]  pts_in_box_coord       Coordinates of points located in boxes
 *
 */

void
PDM_para_octree_points_inside_boxes
(
 const int           octree_id,
 const int           n_boxes,
 const double       *box_extents,
 const PDM_g_num_t  *box_g_num,
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 );

void
PDM_para_octree_points_inside_boxes2
(
 const int           octree_id,
 const int           n_boxes,
 const double       *box_extents,
 const PDM_g_num_t  *box_g_num,
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 );

/**
 *
 * \brief Copy octree data of some ranks into all ranks
 *
 * \param [in]   id                 Identifier
 * \param [in]   n_copied_ranks     Number of ranks to copy
 * \param [in]   copied_ranks       Array of ranks to copy
 *
 */

void
PDM_para_octree_copy_ranks
(
 const int  id,
 const int  n_copied_ranks,
 const int *copied_ranks
 );

/**
 *
 * \brief Free copied data in an octree structure
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_free_copies
(
 const int          id
 );






#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PARA_OCTREE_H */
