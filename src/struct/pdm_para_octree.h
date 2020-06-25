#ifndef PDM_PARA_OCTREE_H
#define	PDM_PARA_OCTREE_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_handles.h"
#include "pdm_morton.h"
#include "pdm_timer.h"

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

#define NTIMER2 11

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \enum pdm_para_octree_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_BOTTOM,
  PDM_UP,
  PDM_SOUTH,
  PDM_NORTH,
  PDM_WEST,
  PDM_EAST,
} PDM_para_octree_direction_t;

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

/**
 * \struct _l_octant_t
 * \brief  Define a list of octants
 *
 */

typedef struct  {

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */

  PDM_morton_code_t *codes;        /*!< Morton codes */

  int  *n_points;          /*!< Number of points in octant*/
  int  *range;             /*!< Start index of point list for each octant */

  int   *neighbour_idx;
  int   *neighbours;               /*!< rank + id_node size = 2 * n_nodes */
  int   dim;

} _l_octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /*!< Translation for the normalization */
  double      d[3];           /*!< Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t        t_n_points;     /*!< total number of points */
  int                n_points;       /*!< Number of points in each cloud */
  double            *points;         /*!< Point coordinates */
  int               *points_icloud;  /*!< Point cloud */
  PDM_g_num_t       *points_gnum;    /*!< Point global number */
  PDM_morton_code_t *points_code;    /*!< Morton codes */

  PDM_morton_code_t *rank_octants_index;
  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

  int  n_part_boundary_elt;    /*!< Number of partitioning boundary element */
  int *part_boundary_elt_idx; /*!< Index for part_boundary_elt (size=\ref n_part_boundary_elt + 1 */
  int *part_boundary_elt;     /*!< Partitioning boundary elements description (proc number + element number) */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER2]; /*!< Elapsed time */

  double times_cpu[NTIMER2];     /*!< CPU time */

  double times_cpu_u[NTIMER2];  /*!< User CPU time */

  double times_cpu_s[NTIMER2];  /*!< System CPU time */

  int neighboursToBuild;

  int  n_connected;
  int *connected_idx;

} _octree_t;

typedef _octree_t PDM_octree_t;

/*  */

/*============================================================================
 * Global variable
 *============================================================================*/

 //extern MANAGED PDM_Handles_t *_octrees;

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
PDM_para_octree_closest_point
(
const int    id,
const int    n_closest_points,
const int    n_pts,
double      *pts,
PDM_g_num_t *pts_g_num,
PDM_g_num_t *closest_octree_pt_g_num,
double      *closest_octree_pt_dist2
);


void
PDM_para_octree_closest_point2
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
 * \brief  transfert _octrees var to gpu through PDM_closest_points_compute
 *
 */

PDM_octree_t *
PDM_para_octree_octrees_transfert
(
  int id
 );


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PARA_OCTREE_H */
