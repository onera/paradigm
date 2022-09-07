/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_octree_seq.h"
#include "pdm_octree_seq_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

static const double _eps_default = 1.e-12;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

static void
_l_octant_set
(
 PDM_octree_seq_t *octree
 )
{
  octree->octants = malloc(sizeof(_l_octant_t));

  _l_octant_t *octants = octree->octants;

  octants->ancestor_id          = malloc(sizeof(int                   ) * octree->n_nodes);
  octants->is_leaf              = malloc(sizeof(int                   ) * octree->n_nodes);
  octants->location_in_ancestor = malloc(sizeof(PDM_octree_seq_child_t) * octree->n_nodes);
  octants->depth                = malloc(sizeof(int                   ) * octree->n_nodes);
  octants->children_id          = malloc(sizeof(int                   ) * octree->n_nodes * 8);
  octants->range                = malloc(sizeof(int                   ) * octree->n_nodes * 2);
  octants->idx                  = malloc(sizeof(int                   ) * octree->n_nodes * 9);
  octants->n_points             = malloc(sizeof(int                   ) * octree->n_nodes);
  octants->extents              = malloc(sizeof(double                ) * octree->n_nodes * 6);

  for (int i = 0; i < octree->n_nodes; i++) {
    _octant_t *node = octree->nodes + i;

    octants->ancestor_id[i] = node->ancestor_id;
    octants->is_leaf[i] = node->is_leaf;
    octants->location_in_ancestor[i] = node->location_in_ancestor;
    octants->depth[i] = node->depth;
    for (int j = 0; j < 8; j++) {
      octants->children_id[8*i+j] = node->children_id[8*i+j];
    }
    for (int j = 0; j < 2; j++) {
      octants->range[2*i+j] = node->range[2*i+j];
    }
    for (int j = 0; j < 9; j++) {
      octants->idx[9*i+j] = node->idx[9*i+j];
    }
    octants->n_points[i] = node->n_points;
    for (int j = 0; j < 6; j++) {
      octants->extents[6*i+j] = node->extents[6*i+j];
    }
  }

}

static void
_l_octant_free
(
 PDM_octree_seq_t *octree
 )
{
  if (octree->octants != NULL) {

    free(octree->octants->ancestor_id);
    free(octree->octants->is_leaf);
    free(octree->octants->location_in_ancestor);
    free(octree->octants->depth);
    free(octree->octants->children_id);
    free(octree->octants->range);
    free(octree->octants->idx);
    free(octree->octants->n_points);
    free(octree->octants->extents);

    free(octree->octants);

    octree->octants = NULL;
  }
}


/**
 *
 * \brief Compute distance to a box
 *
 * \param [in]   dim        Dimension
 * \param [in]   extents    Box extents
 * \param [in]   coords     Point coords
 * \param [out]  min_dist2  Square of minimum distance
 * \param [out]  max_dist2  Sqaure of maximum distance
 *
 * \return 1 if point is in the box, 0 otherwise
 *
 */

inline static int
_box_dist2
(
const int              dim,
const double          *restrict extents,
const double          *restrict coords,
double                *restrict min_dist2
//double                *restrict max_dist2
)
{

  int inbox = 0;

  double _min_dist2 = 0.;
  //  double _max_dist2 = 0.;
  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = extents[i];
    double xmax = extents[i+dim];

    if (x > xmax) {
      double diff_max = x - xmax;
      //double diff_min = x - xmin;
      _min_dist2 += diff_max * diff_max;
      //      _max_dist2 += diff_min * diff_min;
    }

    else if (x < xmin) {
      //double diff_max = x - xmax;
      double diff_min = x - xmin;
      _min_dist2 += diff_min * diff_min;
      //      _max_dist2 += diff_max * diff_max;
    }

    else {
      //double diff_max = x - xmax;
      //double diff_min = x - xmin;
      inbox += 1;
      // _max_dist2 += PDM_MAX (diff_min * diff_min,
      //                       diff_max * diff_max);
    }

  }

  *min_dist2 = _min_dist2;
//  *max_dist2 = _max_dist2;
  return inbox == dim;

}

static int
_intersect_node_box_explicit
(
 int           dim,
 const double *node_extents,
 const double *box_extents,
 int          *inside
)
{
  *inside = 1;

  for (int i = 0; i < dim; i++) {
    if (node_extents[i]   > box_extents[i+3] ||
        node_extents[i+3] < box_extents[i]) {
      return 0;
    }
    else if (node_extents[i]  < box_extents[i] ||
             node_extents[i+3] > box_extents[i+3]) {
      *inside = 0;
    }
  }

  return 1;
}

/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppart_id        ppart identifier
 *
 */

// static _octree_seq_t *
// _get_from_id
// (
//  int  id
// )
// {
//   _octree_seq_t *octree = (_octree_seq_t *) PDM_Handles_get (_octrees, id);

//   if (octree == NULL) {
//     PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
//   }

//   return octree;
// }

/**
 *
 * \brief Build a local octree's leaves.
 *
 * \param[in]  ancestor_id          Ancestor identifier
 * \param[in]  location_in_ancestor Location in ancestor
 * \param[in]  depth                Depth in the tree
 * \param[in]  extents              Extents associated with node:
 *                                  x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 * \param[in]  point_coords         Point coordinates
 * \param[in]  point_ids_tmp        Temporary point indexes
 * \param[in]  pos_tmp              Temporary point position in octree
 * \param[inout]  octree            Current octree structure
 * \param[inout]  point_range       Start and past-the end index in point_idx
 *                                  for current node (size: 2)
 */

static void
_build_octree_seq_leaves(const int                      ancestor_id,
                         const PDM_octree_seq_child_t   location_in_ancestor,
                         const int                      depth,
                         const double                   extents[],
                         const double                 **point_coords,
                         int                           *point_icloud_tmp,
                         int                           *point_ids_tmp,
                         PDM_octree_seq_t              *octree,
                         int                            point_range[2])
{
  int i, j, k, _n_nodes, _n_points, tmp_size;

  int count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];
  _octant_t  *_node;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  _n_nodes = octree->n_nodes;
  tmp_size = octree->n_nodes;

  /* Resize octree if necessary */

  if (octree->n_nodes >= octree->n_nodes_max) {
    if (octree->n_nodes == 0) {
      octree->n_nodes = 1;
      octree->n_nodes_max = 8;
    }
    octree->n_nodes_max *= 2;
    octree->nodes = realloc (octree->nodes,
                             octree->n_nodes_max * sizeof(_octant_t));
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  if (depth < octree->depth_max && _n_points > octree->points_in_leaf_max) {
    /* Extents center */

    for (j = 0; j < 3; j++) {
      mid[j]= (extents[j] + extents[j + 3]) * 0.5;
    }


    /* Count points in each octant */

    for (i = point_range[0]; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      count[k] += 1;
    }

    /* Build index */

    idx[0] = 0;
    for (j = 0; j < 8; j++) {
      idx[j+1] = idx[j] + count[j];
    }

    for (j = 0; j < 8; j++) {
      count[j] = 0;
    }

    for (i = point_range[0], j = 0; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      point_icloud_tmp[idx[k] + count[k]] = octree->point_icloud[i];
      point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
      count[k] += 1;
    }

    /* Check if this subdivision is static
       and check coordinates to find multi point */

    for (i = point_range[0], j = 0; i < point_range[1]; i++, j++) {
      octree->point_icloud[i] = point_icloud_tmp[j];
      octree->point_ids[i] = point_ids_tmp[j];
    }

    for (i = 0; i < 9; i++)
      idx[i] = point_range[0] + idx[i];

    /* Build leaves recursively */

    for (i = 0; i < 8; i++) {

      if ((idx[i+1] - idx[i]) > 0) {

        tmp_size++;

        octant_id[i] = tmp_size;

        if (i < 4) {
          sub_extents[0] = extents[0];
          sub_extents[3] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
          sub_extents[3] = extents[3];
        }
        /* 1.0e-12 term in assert() used to allow for
           truncation error in for xmin = xmax case */
        assert(sub_extents[0] < sub_extents[3] + 1.0e-14);

        if (i%4 < 2) {
          sub_extents[1] = extents[1];
          sub_extents[4] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
          sub_extents[4] = extents[4];
        }
        assert(sub_extents[1] < sub_extents[4] + 1.0e-14);

        if (i%2 < 1) {
          sub_extents[2] = extents[2];
          sub_extents[5] = mid[2];
        }
        else {
          sub_extents[2] = mid[2];
          sub_extents[5] = extents[5];
        }
        assert(sub_extents[2] < sub_extents[5] + 1.0e-14);

        octree->n_nodes = tmp_size;

        _build_octree_seq_leaves(_n_nodes,
                             (PDM_octree_seq_child_t) i,
                             depth+1,
                             sub_extents,
                             point_coords,
                             point_icloud_tmp,
                             point_ids_tmp,
                             octree,
                             idx + i);

        tmp_size = octree->n_nodes;
      }

    }

  }
  /* Finalize node */


  _node = octree->nodes + _n_nodes;

  _node->range[0] = point_range[0];
  _node->range[1] = point_range[1];

  for (i = 0; i < 9; i++) {
    _node->idx[i] = idx[i];
  }

  for (i = 0; i < 6; i++) {
    _node->extents[i] = extents[i];
  }

  for (i = 0; i < 8; i++) {
    _node->children_id[i] = octant_id[i];
  }

  _node->is_leaf = (_node->children_id[0] == -1) &&
                (_node->children_id[1] == -1) &&
                (_node->children_id[2] == -1) &&
                (_node->children_id[3] == -1) &&
                (_node->children_id[4] == -1) &&
                (_node->children_id[5] == -1) &&
                (_node->children_id[6] == -1) &&
                (_node->children_id[7] == -1);


  _node->ancestor_id = ancestor_id;
  _node->depth = depth;

  _node->n_points = _n_points;
  _node->location_in_ancestor = location_in_ancestor;
}




static void
_build_octree_seq_leaves2(const int                      ancestor_id,
                          const PDM_octree_seq_child_t   location_in_ancestor,
                          const int                      depth,
                          const double                   extents[],
                          const double                 **point_coords,
                          int                           *point_icloud_tmp,
                          int                           *point_ids_tmp,
                          PDM_octree_seq_t              *octree,
                          int                            point_range[2])
{
  int i, j, k, _n_nodes, _n_points, tmp_size;

  int count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];

  _l_octant_t *octants = octree->octants;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  _n_nodes = octree->n_nodes2;
  tmp_size = octree->n_nodes2;

  /* Resize octree if necessary */

  if (octree->n_nodes2 >= octree->n_nodes2_max) {
    if (octree->n_nodes2 == 0) {
      octree->n_nodes2 = 1;
      octree->n_nodes2_max = 8;
    }
    octree->n_nodes2_max *= 2;

    octants->ancestor_id          = realloc(octants->ancestor_id,          sizeof(int                   ) * octree->n_nodes2_max);
    octants->is_leaf              = realloc(octants->is_leaf,              sizeof(int                   ) * octree->n_nodes2_max);
    octants->location_in_ancestor = realloc(octants->location_in_ancestor, sizeof(PDM_octree_seq_child_t) * octree->n_nodes2_max);
    octants->depth                = realloc(octants->depth,                sizeof(int                   ) * octree->n_nodes2_max);
    octants->children_id          = realloc(octants->children_id,          sizeof(int                   ) * octree->n_nodes2_max * 8);
    octants->range                = realloc(octants->range,                sizeof(int                   ) * octree->n_nodes2_max * 2);
    octants->idx                  = realloc(octants->idx,                  sizeof(int                   ) * octree->n_nodes2_max * 9);
    octants->n_points             = realloc(octants->n_points,             sizeof(int                   ) * octree->n_nodes2_max);
    octants->extents              = realloc(octants->extents,              sizeof(double                ) * octree->n_nodes2_max * 6);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  if (depth < octree->depth_max && _n_points > octree->points_in_leaf_max) {
    /* Extents center */

    for (j = 0; j < 3; j++) {
      mid[j]= (extents[j] + extents[j + 3]) * 0.5;
    }


    /* Count points in each octant */

    for (i = point_range[0]; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      count[k] += 1;
    }

    /* Build index */

    idx[0] = 0;
    for (j = 0; j < 8; j++) {
      idx[j+1] = idx[j] + count[j];
    }

    for (j = 0; j < 8; j++) {
      count[j] = 0;
    }

    for (i = point_range[0], j = 0; i < point_range[1]; i++) {

      for (j = 0, k = 0; j < 3; j++) {
        if (point_coords[octree->point_icloud[i]][octree->point_ids[i]*3 + j] > mid[j])
          k += octant_mask[j];
      }

      point_icloud_tmp[idx[k] + count[k]] = octree->point_icloud[i];
      point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
      count[k] += 1;
    }

    /* Check if this subdivision is static
       and check coordinates to find multi point */

    for (i = point_range[0], j = 0; i < point_range[1]; i++, j++) {
      octree->point_icloud[i] = point_icloud_tmp[j];
      octree->point_ids[i] = point_ids_tmp[j];
    }

    for (i = 0; i < 9; i++)
      idx[i] = point_range[0] + idx[i];

    /* Build leaves recursively */

    for (i = 0; i < 8; i++) {

      if ((idx[i+1] - idx[i]) > 0) {

        tmp_size++;

        octant_id[i] = tmp_size;

        if (i < 4) {
          sub_extents[0] = extents[0];
          sub_extents[3] = mid[0];
        }
        else {
          sub_extents[0] = mid[0];
          sub_extents[3] = extents[3];
        }
        /* 1.0e-12 term in assert() used to allow for
           truncation error in for xmin = xmax case */
        assert(sub_extents[0] < sub_extents[3] + 1.0e-14);

        if (i%4 < 2) {
          sub_extents[1] = extents[1];
          sub_extents[4] = mid[1];
        }
        else {
          sub_extents[1] = mid[1];
          sub_extents[4] = extents[4];
        }
        assert(sub_extents[1] < sub_extents[4] + 1.0e-14);

        if (i%2 < 1) {
          sub_extents[2] = extents[2];
          sub_extents[5] = mid[2];
        }
        else {
          sub_extents[2] = mid[2];
          sub_extents[5] = extents[5];
        }
        assert(sub_extents[2] < sub_extents[5] + 1.0e-14);

        octree->n_nodes2 = tmp_size;

        _build_octree_seq_leaves2(_n_nodes,
                                  (PDM_octree_seq_child_t) i,
                                  depth+1,
                                  sub_extents,
                                  point_coords,
                                  point_icloud_tmp,
                                  point_ids_tmp,
                                  octree,
                                  idx + i);

        tmp_size = octree->n_nodes2;
      }

    }

  }
  /* Finalize node */



  for (i = 0; i < 2; i++) {
    octants->range[2*_n_nodes + i] = point_range[i];
  }

  for (i = 0; i < 9; i++) {
    octants->idx[9*_n_nodes + i] = idx[i];
  }

  for (i = 0; i < 6; i++) {
    octants->extents[6*_n_nodes + i] = extents[i];
  }

  for (i = 0; i < 8; i++) {
    octants->children_id[8*_n_nodes + i] = octant_id[i];
  }

  octants->is_leaf[_n_nodes] =
  (octants->children_id[8*_n_nodes + 0] == -1) &&
  (octants->children_id[8*_n_nodes + 1] == -1) &&
  (octants->children_id[8*_n_nodes + 2] == -1) &&
  (octants->children_id[8*_n_nodes + 3] == -1) &&
  (octants->children_id[8*_n_nodes + 4] == -1) &&
  (octants->children_id[8*_n_nodes + 5] == -1) &&
  (octants->children_id[8*_n_nodes + 6] == -1) &&
  (octants->children_id[8*_n_nodes + 7] == -1);

  octants->ancestor_id[_n_nodes] = ancestor_id;
  octants->depth[_n_nodes]       = depth;

  octants->n_points[_n_nodes]             = _n_points;
  octants->location_in_ancestor[_n_nodes] = location_in_ancestor;
}


/**
 *
 * \brief   Compute extents of a point set
 *
 *  \param [in] dim         Space dimension of points to locate_3d
 *  \param [in] n_points    Number of points to locate
 *  \param [in] point_index optional indirection array to point_coords
 *                          (1 to n_points numbering)
 *  \param [in] point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 */

static void
_point_extents(const int     dim,
               const int     n_points,
               const int     point_index[],
               const double  point_coords[],
               double        extents[])
{
  int i;
  int j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }

    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }
}


/**
 *
 * \brief Build an octree
 *
 * \param[in]  octree    Current octree
 * .
 */

static void
_build_octree
(
PDM_octree_seq_t *octree
)
{

  int point_range[2];

  /* Initialization */

  octree->n_nodes = 0;
  octree->n_nodes_max = 0;
  octree->n_nodes2 = 0;
  octree->n_nodes2_max = 0;
  octree->nodes = NULL;


  octree->octants = malloc(sizeof(_l_octant_t));
  octree->octants->ancestor_id          = NULL;
  octree->octants->is_leaf              = NULL;
  octree->octants->location_in_ancestor = NULL;
  octree->octants->depth                = NULL;
  octree->octants->children_id          = NULL;
  octree->octants->range                = NULL;
  octree->octants->idx                  = NULL;
  octree->octants->n_points             = NULL;
  octree->octants->extents              = NULL;

  for (int i = 0; i < octree->n_point_clouds; i++) {
    octree->t_n_points += octree->n_points[i];
  };

  octree->point_ids = malloc (sizeof(int) * octree->t_n_points);
  octree->point_icloud = malloc (sizeof(int) * octree->t_n_points);

  int cpt = 0;
  for (int i = 0; i < octree->n_point_clouds; i++) {

    int n_points = octree->n_points[i];
    double extents[6];

    if (n_points > 0) {

      _point_extents(3,
                     n_points,
                     NULL,
                     octree->point_clouds[i],
                     extents);

      for (int i1 = 0; i1 < 3; i1++) {
        octree->extents[i1] = PDM_MIN (extents[i1], octree->extents[i1]);
        octree->extents[i1 + 3] = PDM_MAX (extents[i1 + 3], octree->extents[i1 + 3]);
      }

      for (int j = 0; j < n_points; j++) {
        octree->point_ids[cpt] = j;
        octree->point_icloud[cpt] = i;
        cpt +=1;
      }
    }
  }

  double delta = -1;
  for (int i = 0; i < 3; i++) {
    delta = PDM_MAX (octree->tolerance * (octree->extents[i + 3] - octree->extents[i]),
                     delta);
  }
  delta = PDM_MAX (delta,_eps_default);

  for (int i = 0; i < 3; i++) {
    octree->extents[i] += -1.01*delta;
    octree->extents[i + 3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = octree->t_n_points;

  int *point_ids_tmp = malloc (sizeof(int) * octree->t_n_points);
  int *point_icloud_tmp = malloc (sizeof(int) * octree->t_n_points);

  _build_octree_seq_leaves(-1,
                       (PDM_octree_seq_child_t) 0,
                       -1,
                       octree->extents,
                       (const double **) octree->point_clouds,
                       point_icloud_tmp,
                       point_ids_tmp,
                       octree,
                       point_range);

  octree->n_nodes +=1;

  free (point_ids_tmp);
  free (point_icloud_tmp);


  point_range[0] = 0;
  point_range[1] = octree->t_n_points;

  point_ids_tmp    = malloc (sizeof(int) * octree->t_n_points);
  point_icloud_tmp = malloc (sizeof(int) * octree->t_n_points);

  _build_octree_seq_leaves2(-1,
                       (PDM_octree_seq_child_t) 0,
                       -1,
                       octree->extents,
                       (const double **) octree->point_clouds,
                       point_icloud_tmp,
                       point_ids_tmp,
                       octree,
                       point_range);

  octree->n_nodes2 +=1;

  assert(octree->n_nodes ==octree->n_nodes2);

  free (point_ids_tmp);
  free (point_icloud_tmp);

}

/*=============================================================================
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
 *
 * \return     Pointer to \ref PDM_octree_seq object
 */

PDM_octree_seq_t *
PDM_octree_seq_create
(
 const int    n_point_cloud,
 const int    depth_max,
 const int    points_in_leaf_max,
 const double tolerance
)
{
  PDM_octree_seq_t *octree = (PDM_octree_seq_t *) malloc(sizeof(PDM_octree_seq_t));

  octree->n_point_clouds = n_point_cloud;
  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;
  octree->tolerance = tolerance;

  octree->n_nodes = 0;
  octree->n_nodes_max = 0;

  octree->n_points = malloc (sizeof(double) * n_point_cloud);
  octree->point_clouds = malloc (sizeof(double *) * n_point_cloud);
  for (int i = 0; i < n_point_cloud; i++) {
    octree->n_points[i] = 0;
    octree->point_clouds[i] = NULL;
  }

  octree->point_icloud = NULL;
  octree->point_ids = NULL;
  octree->nodes = NULL;
  octree->t_n_points = 0;
  for (int i = 0; i < 3; i++) {
    octree->extents[i]     =  HUGE_VAL;
    octree->extents[i + 3] = -HUGE_VAL;
  }

  octree->octants = NULL;

  return (PDM_octree_seq_t *) octree;
}


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 */

void
PDM_octree_seq_free
(
 PDM_octree_seq_t *octree
)
{

  free (octree->n_points);
  free (octree->point_clouds);
  free (octree->point_ids);
  free (octree->nodes);
  free (octree->point_icloud);

  _l_octant_free(octree);

  free (octree);
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 *
 */


void
PDM_octree_seq_point_cloud_set
(
 PDM_octree_seq_t *octree,
 const int         i_point_cloud,
 const int         n_points,
 const double     *coords
)
{

  octree->n_points[i_point_cloud] = n_points;
  octree->point_clouds[i_point_cloud] = coords;
}


/**
 *
 * \brief Build octree
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 */

void
PDM_octree_seq_build
(
 PDM_octree_seq_t *octree
)
{

  if (octree->nodes == NULL) {
    _build_octree (octree);

    // _l_octant_set(_octree);
  }

}


/**
 *
 * \brief Get root node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 * \return     Root node identifier (-1 if octree is not built)
 *
 */

int
PDM_octree_seq_root_node_id_get
(
 PDM_octree_seq_t *octree
)
{

  if (octree->nodes == NULL) {
    return -1;
  }
  else {
    return 0;
  }
}

/**
 *
 * \brief Get extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 * \return     Extents
 *
 */

double *
PDM_octree_seq_extents_get
(
 PDM_octree_seq_t *octree
)
{

  return octree->extents;
}


/**
 *
 * \brief Get ancestor node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return     Ancestor node identifier
 *
 */

int
PDM_octree_seq_ancestor_node_id_get
(
 PDM_octree_seq_t *octree,
 const int         node_id
)
{

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].ancestor_id;
}


/**
 *
 * \brief Get node extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return     Extents
 *
 */

const double *
PDM_octree_seq_node_extents_get
(
 PDM_octree_seq_t *octree,
 const int         node_id
)
{

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].extents;
}


/**
 *
 * \brief Get children of a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [in]   child              Children
 *
 * \return     Children node id
 *
 */

int
PDM_octree_seq_children_get
(
 PDM_octree_seq_t             *octree,
 const int                     node_id,
 const PDM_octree_seq_child_t  child
)
{

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].children_id[child];
}


/**
 *
 * \brief Get Neighbor of node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [in]   direction          Neighbor direction
 *
 * \return     Neighbor node id (-1 if no neighbor)
 *
 */

int
PDM_octree_seq_neighbor_get
(
 PDM_octree_seq_t                 *octree,
 const int                         node_id,
 const PDM_octree_seq_direction_t  direction
)
{

  assert (node_id < octree->n_nodes);

  int neighbor_id = -1;

  _octant_t *octant = &(octree->nodes[node_id]);

  int isSameDiretion = 1;
  while (isSameDiretion && (octant->ancestor_id != 0)) {

    switch (direction) {
    case PDM_OCTREE_SEQ_NADIR:
    case PDM_OCTREE_SEQ_ZENITH:
      isSameDiretion = octant->location_in_ancestor%2  == direction;
      break;
    case PDM_OCTREE_SEQ_WEST:
    case PDM_OCTREE_SEQ_EAST:
      isSameDiretion = ((octant->location_in_ancestor%4 < 2) + 2) == direction;
      break;
    case PDM_OCTREE_SEQ_NORTH:
    case PDM_OCTREE_SEQ_SOUTH:
      isSameDiretion = ((octant->location_in_ancestor < 4) + 4) == direction;
      break;
    }

   octant = &(octree->nodes[octant->ancestor_id]);
  }

  if (octant->ancestor_id != 0) {
    if (direction < 2) {
      neighbor_id = (octant->location_in_ancestor/2) * 2 +
                    (octant->location_in_ancestor + 1)%2;
    }
    if (direction < 4) {
      neighbor_id = (octant->location_in_ancestor/4) * 4 +
                    (octant->location_in_ancestor + 2)%4;
    }
    else {
      neighbor_id = (octant->location_in_ancestor + 4)%8;
    }
  }

  return neighbor_id;
}

/**
 *
 * \brief Get the number of point inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return   Number of points
 *
 */

int
PDM_octree_seq_n_points_get
(
 PDM_octree_seq_t        *octree,
 const int                node_id
)
{

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].range[1] - octree->nodes[node_id].range[0];

}


/**
 *
 * \brief Get indexes of points inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [out]  point_clouds_id    Point clouds number
 *                                  (size = Number of points inside the node)
 * \param [out]  point_indexes      Point indexes
 *                                  (size = Number of points inside the node)
 *
 */

void
PDM_octree_seq_points_get
(
 PDM_octree_seq_t        *octree,
 const int                node_id,
 int                    **point_clouds_id,
 int                    **point_indexes
)
{

  assert (node_id < octree->n_nodes);

  *point_clouds_id = octree->point_icloud + octree->nodes[node_id].range[0];

  *point_indexes = octree->point_ids + octree->nodes[node_id].range[0];
}

/**
 *
 * \brief Get indexes of points inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [out]  point_clouds_id    Point clouds number
 *                                  (size = Number of points inside the node)
 * \param [out]  point_indexes      Point indexes
 *                                  (size = Number of points inside the node)
 *
 */

void
PDM_octree_seq_extract_extent
(
  PDM_octree_seq_t  *octree,
  int                root_id,
  int                n_depth,
  int               *n_box,
  double           **box_extents
)
{
  _l_octant_t *octants = octree->octants;

  int n_children   = 8;
  int s_pt_stack   = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int *stack_id    = malloc (s_pt_stack * sizeof(int              ));
  int *stack_depth  = malloc (s_pt_stack * sizeof(int              ));

  // int n_extract_max = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int *id_to_extract = malloc( octree->n_nodes * sizeof(int));

  int n_extract = 0;
  int pos_stack = 0;
  stack_id   [pos_stack] = root_id;
  stack_depth[pos_stack] = 0;
  pos_stack++;
  while(pos_stack > 0) {

    /* Inspect node */
    --pos_stack;
    int node_id = stack_id   [pos_stack];
    int depth   = stack_depth[pos_stack];

    if(octants->is_leaf[node_id] || depth == n_depth) {
      if(octants->n_points[node_id] > 0) {
        id_to_extract[n_extract++] = node_id;
      }
    } else {
      for (int i = 0; i < n_children; i++) {
        int child_id = octants->children_id[8*node_id+i];
        if (child_id < 0) {
          continue;
        }

        if(depth < n_depth) {
          stack_id   [pos_stack] = child_id;
          stack_depth[pos_stack] = depth + 1;
          pos_stack++;
        }
      }
    }
  }
  free(stack_id);
  free(stack_depth);

  double* _extents = malloc(n_extract * 6 * sizeof(double));
  for(int i = 0; i < n_extract; ++i) {
    int node_id = id_to_extract[i];
    for(int k = 0; k < 6; ++k) {
      _extents[6*i+k] = octants->extents[6*node_id+k];
    }
  }

  *n_box       = n_extract;
  *box_extents = _extents;

  free(id_to_extract);

}

/**
 *
 * \brief Is it a leaf
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return   1 or 0
 *
 */

int
PDM_octree_seq_leaf_is
(
 PDM_octree_seq_t        *octree,
 const int                node_id
)
{

  assert (node_id < octree->n_nodes);

  return octree->nodes[node_id].is_leaf;
}


/**
 *
 * \brief Look for closest points stored inside an octree
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree_seq object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_octree_pt_id   Closest point in octree index
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_octree_seq_closest_point
(
PDM_octree_seq_t *octree,
const int         n_pts,
double           *pts,
int              *closest_octree_pt_id,
double           *closest_octree_pt_dist2
)
{

  const int n_children = 8;


  int s_pt_stack = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int sort_child[n_children];
  double dist_child[n_children];
  int inbox_child[n_children];

  int *stack = malloc ((sizeof(int)) * s_pt_stack);
  int *inbox_stack = malloc ((sizeof(int)) * s_pt_stack);
  double *min_dist2_stack = malloc ((sizeof(double)) * s_pt_stack);

  int dim = 3;
  /* double ptprintc[3] = {5.83333e-01,  6.25000e-01,  5.00000e-01}; */
  /* double ptprintc2[3] = {5.83333e-01,  1.,  5.00000e-01}; */
  //int count = 0;
  for (int i = 0; i < n_pts; i++) {

    int pos_stack = 0;
    const double *_pt = pts + dim * i;

      /* double t1 = (_pt[0] - ptprintc[0]) * (_pt[0] - ptprintc[0]); */
      /* double t2 = (_pt[1] - ptprintc[1]) * (_pt[1] - ptprintc[1]); */
      /* double t3 = (_pt[2] - ptprintc[2]) * (_pt[2] - ptprintc[2]); */
      /* t1=1; */
      /* if ((t1+t2+t3) < 1e-6) */
      /* printf ("\n ******** pt %d : %12.5e, %12.5e, %12.5e\n", i, _pt[0], _pt[1], _pt[2]); */

    /* Init stack */

    closest_octree_pt_id[2*i] = -1;
    closest_octree_pt_id[2*i+1] = -1;
    closest_octree_pt_dist2[i] = HUGE_VAL;

    stack[pos_stack] = 0; /* push root in th stack */

    double _min_dist2;
    _octant_t *root_node = &(octree->nodes[0]);
    int inbox1 =  _box_dist2 (dim,
                              root_node->extents,
                              _pt,
                              &_min_dist2);

    inbox_stack[pos_stack] = inbox1;
    min_dist2_stack[pos_stack] = _min_dist2;
    pos_stack++;

    while (pos_stack > 0) {
      //count++;
      int id_curr_node = stack[--pos_stack];
      _octant_t *curr_node = &(octree->nodes[id_curr_node]);

      double min_dist2 = min_dist2_stack[pos_stack];
      int inbox =  inbox_stack[pos_stack];

      /* int inbox =  _box_dist2 (dim, */
      /*                          curr_node->extents, */
      /*                          _pt, */
      /*                          &min_dist2); */

      if ((min_dist2 <= closest_octree_pt_dist2[i]) || (inbox == 1)) {

        if (!curr_node->is_leaf) {

          /* Sort children and store them into the stack */

          const int *_child_ids = curr_node->children_id;

          for (int j = 0; j < n_children; j++) {
            dist_child[j] = HUGE_VAL;
          }

          int n_selec = 0;
          for (int j = 0; j < n_children; j++) {


            int child_id = _child_ids[j];


            int child_inbox = 0;

            if (child_id != -1) {

              _octant_t *child_node = &(octree->nodes[child_id]);
              double child_min_dist2;

              child_inbox = _box_dist2 (dim,
                                        child_node->extents,
                                        _pt,
                                        &child_min_dist2);

              int i1 = 0;
              for (i1 = n_selec;
                   (i1 > 0) && (dist_child[i1-1] > child_min_dist2) ; i1--) {
                dist_child[i1] = dist_child[i1-1];
                sort_child[i1] = sort_child[i1-1];
                inbox_child[i1] = inbox_child[i1-1];
              }

              sort_child[i1] = child_id;
              dist_child[i1] = child_min_dist2;
              inbox_child[i1] = child_inbox;

              n_selec += 1;

            }
          }

          for (int j = 0; j < n_selec; j++) {
            int j1 = n_selec- 1 - j;
            int child_id = sort_child[j1];
            if (child_id != -1) {
              _octant_t *child_node = &(octree->nodes[child_id]);
              if ((dist_child[j1] < closest_octree_pt_dist2[i]) &&
                  (child_node->n_points > 0)) {

                min_dist2_stack[pos_stack] = dist_child[j1];
                inbox_stack[pos_stack] = inbox_child[j1];

                stack[pos_stack++] = child_id; /* push root in th stack */
              }
            }
          }
        }

        else {

          int *point_clouds_id = octree->point_icloud + curr_node->range[0];
          int *point_indexes = octree->point_ids + curr_node->range[0];

          for (int j = 0; j < curr_node->n_points; j++) {

            double point_dist2 = 0;
            const double *_coords = octree->point_clouds[point_clouds_id[j]]
                                    + dim * point_indexes[j];

            for (int k = 0; k < dim; k++) {
              point_dist2 += (_coords[k] - _pt[k]) *
                             (_coords[k] - _pt[k]);
            }

            if (point_dist2 < closest_octree_pt_dist2[i]) {
              closest_octree_pt_id[2*i] = point_clouds_id[j];
              closest_octree_pt_id[2*i+1] = point_indexes[j];
              closest_octree_pt_dist2[i] = point_dist2;
            }
          }
        }
      }
    }
      /*     if ((t1+t2+t3) < 1e-6) */
      /* printf ("\n ******** fin point ******************\n"); */

  }

  //if (n_pts != 0) printf("n nodes per point = %d\n", count / n_pts);

  free (inbox_stack);
  free (min_dist2_stack);
  free (stack);

}



/**
 *
 * \brief Get points located inside a set of boxes
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree_seq object
 * \param [in]   n_box                  Number of boxes
 * \param [in]   box_extents            Extents of boxes
 * \param [out]  pts_idx                Index of points located in boxes
 * \param [out]  pts_l_num              Local ids of points located in boxes
 *
 */

void
PDM_octree_seq_inside_boxes
(
       PDM_octree_seq_t   *octree,
 const int                 n_box,
 const double              box_extents[],
       int               **pts_idx,
       int               **pts_l_num
)
{

  *pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_pts_idx = *pts_idx;
  _pts_idx[0] = 0;

  if (n_box < 1) {
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }

  int dim = 3;
  const int n_children = 8;


  int s_pt_stack = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int *stack_id  = malloc (s_pt_stack * sizeof(int              ));

  int node_inside_box;
  int intersect;

  int tmp_size = 4 * n_box;
  *pts_l_num = malloc (sizeof(int) * tmp_size);
  int *_pts_l_num = *pts_l_num;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int dbg_enabled = 0; //(box_g_num[ibox] == 2793384);//

    _pts_idx[ibox+1] = _pts_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;
    const double *box_min      = box_extents + 6*ibox;
    const double *box_max      = box_min + 3;

    _octant_t *root_node = &(octree->nodes[0]);
    intersect = _intersect_node_box_explicit (3,
                                              root_node->extents,
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      continue;
    }


    if (node_inside_box) {
      /* The box must contain all points */
      if (dbg_enabled) {
        printf("    add pts with lnum %d through %d\n", root_node->range[0], root_node->range[0] + root_node->n_points);
      }
      int new_size = _pts_idx[ibox+1] + root_node->n_points;

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
        _pts_l_num = *pts_l_num;

      }

      for (int j = 0; j < root_node->n_points; j++) {
        _pts_l_num[_pts_idx[ibox+1]++] = root_node->range[0] + j;
      }
      continue;
    } /* End node_inside_box */


    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {
      int node_id = stack_id[--pos_stack];
      _octant_t *curr_node = &(octree->nodes[node_id]);
      if (dbg_enabled) {
        printf("  node %d, range=%d/%d, n_points=%d, leaf_id=%d\n",
               node_id,
               curr_node->range[0], curr_node->range[1],
               curr_node->n_points,
               curr_node->is_leaf);
      }

      /* is leaf */
      if(curr_node->is_leaf == 1) {

        int *point_clouds_id = octree->point_icloud + curr_node->range[0];
        int *point_indexes   = octree->point_ids    + curr_node->range[0];

        for (int i = 0; i < curr_node->n_points; i++) {
          int ipt = curr_node->range[0] + i;
          const double      *_pt       = octree->point_clouds[point_clouds_id[i]] + dim * point_indexes[i];
          // const PDM_g_num_t *_pt_g_num = octree->point_gnum  [point_clouds_id[i]] +       point_indexes[i];

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _pts_idx[ibox+1] + 1);
              *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
              _pts_l_num = *pts_l_num;
            }

            _pts_l_num[_pts_idx[ibox+1]++] = ipt;
            // if (dbg_enabled) {
            //   printf("    add point %d ("PDM_FMT_G_NUM")\n", ipt, _pt_g_num[ipt]);
            // }
          }
        }
      } else { /* Internal nodes */

        const int *_child_ids = curr_node->children_id;
        for (int i = 0; i < n_children; i++) {
          int child_id = _child_ids[i];
          if (child_id < 0) {
            continue;
          }
          _octant_t *child_node = &(octree->nodes[child_id]);

          if (dbg_enabled) {
            printf("    child %d: id=%d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                   i,
                   child_id,
                   child_node->range[0],
                   child_node->range[1],
                   child_node->n_points,
                   child_node->is_leaf);
            printf("    pts_extents = %f %f %f %f %f %f\n",
                   child_node->extents[0],
                   child_node->extents[1],
                   child_node->extents[2],
                   child_node->extents[3],
                   child_node->extents[4],
                   child_node->extents[5]);
          }

          intersect = _intersect_node_box_explicit (3,
                                                    child_node->extents,
                                                    _box_extents,
                                                    &node_inside_box);

          if (dbg_enabled) {
            printf("    intersect = %d\n", intersect);
          }

          if (intersect) {
            if (node_inside_box) {
              /* The box must contain all points */
              if (dbg_enabled) {
                printf("    add pts with lnum %d through %d\n", child_node->range[0], child_node->range[1]);
              }

              int new_size = _pts_idx[ibox+1] + child_node->n_points;

              if (tmp_size <= new_size) {
                tmp_size = PDM_MAX (2*tmp_size, new_size);
                *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
                _pts_l_num = *pts_l_num;
              }

              for (int j = 0; j < child_node->n_points; j++) {
                _pts_l_num[_pts_idx[ibox+1]++] = child_node->range[0] + j;
              }
            }

            else {
              /* Push child in stack */
              stack_id[pos_stack++] = child_id;
            }
          }
        } // End of loop on children
      }

    } /* End While */
  } /* End boxe loop */

  free (stack_id);
  *pts_l_num = realloc (*pts_l_num, sizeof(int) * _pts_idx[n_box]);
}



void
PDM_octree_seq_inside_boxes2
(
       PDM_octree_seq_t   *octree,
 const int                 n_box,
 const double              box_extents[],
       int               **pts_idx,
       int               **pts_l_num
)
{
  *pts_idx = malloc (sizeof(int) * (n_box + 1));
  int *_pts_idx = *pts_idx;
  _pts_idx[0] = 0;

  if (n_box < 1) {
    *pts_l_num = malloc (sizeof(int) * _pts_idx[n_box]);
    return;
  }

  int dim = 3;
  const int n_children = 8;


  _l_octant_t *octants = octree->octants;

  int s_pt_stack = ((n_children - 1) * (octree->depth_max - 1) + n_children);
  int *stack_id  = malloc (s_pt_stack * sizeof(int              ));

  int node_inside_box;
  int intersect;

  int tmp_size = 4 * n_box;
  *pts_l_num = malloc (sizeof(int) * tmp_size);
  int *_pts_l_num = *pts_l_num;

  for (int ibox = 0; ibox < n_box; ibox++) {
    int dbg_enabled = 0; //(box_g_num[ibox] == 2793384);//
    if (dbg_enabled) {
      log_trace("box %d\n", ibox);
    }

    _pts_idx[ibox+1] = _pts_idx[ibox];

    const double *_box_extents = box_extents + 6*ibox;
    const double *box_min      = box_extents + 6*ibox;
    const double *box_max      = box_min + 3;

    intersect = _intersect_node_box_explicit (3,
                                              &octants->extents[0],
                                              _box_extents,
                                              &node_inside_box);

    if (!intersect) {
      continue;
    }


    if (node_inside_box) {
      /* The box must contain all points */
      if (dbg_enabled) {
        log_trace("    add pts with lnum %d through %d\n", octants->range[0], octants->range[1] + octants->n_points[0]);
      }
      int new_size = _pts_idx[ibox+1] + octants->n_points[0];

      if (tmp_size <= new_size) {
        tmp_size = PDM_MAX (2*tmp_size, new_size);
        *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
        _pts_l_num = *pts_l_num;

      }

      for (int j = 0; j < octants->n_points[0]; j++) {
        _pts_l_num[_pts_idx[ibox+1]++] = octants->range[0] + j;
      }
      continue;
    } /* End node_inside_box */


    /* Push root in stack */
    int pos_stack = 0;
    stack_id[pos_stack++] = 0;

    while (pos_stack > 0) {
      int node_id = stack_id[--pos_stack];

      if (dbg_enabled) {
        log_trace("  node %d, range=%d/%d, n_points=%d, leaf_id=%d\n",
               node_id,
               octants->range[2*node_id], octants->range[2*node_id+1],
               octants->n_points[node_id],
               octants->is_leaf[node_id]);
      }

      /* is leaf */
      if(octants->is_leaf[node_id]) {

        int *point_clouds_id = octree->point_icloud + octants->range[2*node_id];
        int *point_indexes   = octree->point_ids    + octants->range[2*node_id];

        for (int i = 0; i < octants->n_points[node_id]; i++) {
          int ipt = octants->range[2*node_id] + i;
          const double      *_pt       = octree->point_clouds[point_clouds_id[i]] + dim * point_indexes[i];

          int pt_inside_box = 1;
          for (int idim = 0; idim < 3; idim++) {
            if (_pt[idim] < box_min[idim] || _pt[idim] > box_max[idim]) {
              pt_inside_box = 0;
              break;
            }
          }

          if (pt_inside_box) {
            if (_pts_idx[ibox+1] >= tmp_size) {
              tmp_size = PDM_MAX (2*tmp_size, _pts_idx[ibox+1] + 1);
              *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
              _pts_l_num = *pts_l_num;
            }

            _pts_l_num[_pts_idx[ibox+1]++] = ipt;
          }
        }
      } else { /* Internal nodes */

        const int *_child_ids = octants->children_id + 8*node_id;
        for (int i = 0; i < n_children; i++) {
          int child_id = _child_ids[i];
          if (child_id < 0) {
            continue;
          }

          if (dbg_enabled) {
            log_trace("    child %d: id=%d, range=%d/%d, n_points=%d, leaf_id=%d\n",
                   i,
                   child_id,
                   octants->range[2*child_id+0],
                   octants->range[2*child_id+1],
                   octants->n_points[child_id],
                   octants->is_leaf[child_id]);
            log_trace("    pts_extents = %f %f %f %f %f %f\n",
                   octants->extents[6*child_id+0],
                   octants->extents[6*child_id+1],
                   octants->extents[6*child_id+2],
                   octants->extents[6*child_id+3],
                   octants->extents[6*child_id+4],
                   octants->extents[6*child_id+5]);
          }

          intersect = _intersect_node_box_explicit (3,
                                                    // child_node->extents,
                                                    octants->extents + 6*child_id,
                                                    _box_extents,
                                                    &node_inside_box);

          if (dbg_enabled) {
            log_trace("    intersect = %d\n", intersect);
          }

          if (intersect) {
            if (node_inside_box) {
              /* The box must contain all points */
              if (dbg_enabled) {
                log_trace("    add pts with lnum %d through %d\n", octants->range[2*child_id+0], octants->range[2*child_id+1]);
              }

              int new_size = _pts_idx[ibox+1] + octants->n_points[child_id];

              if (tmp_size <= new_size) {
                tmp_size = PDM_MAX (2*tmp_size, new_size);
                *pts_l_num = realloc (*pts_l_num, sizeof(int) * tmp_size);
                _pts_l_num = *pts_l_num;
              }

              for (int j = 0; j < octants->n_points[child_id]; j++) {
                _pts_l_num[_pts_idx[ibox+1]++] = octants->range[2*child_id] + j;
              }
            }

            else {
              /* Push child in stack */
              stack_id[pos_stack++] = child_id;
            }
          }
        } // End of loop on children
      }

    } /* End While */
  } /* End boxe loop */

  free (stack_id);
  *pts_l_num = realloc (*pts_l_num, sizeof(int) * _pts_idx[n_box]);
}



/**
 *
 * \brief Write octants in a VTK file
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree_seq object
 * \param [in]   filename               Output file name
 *
 */

void PDM_octree_seq_write_octants
(
 PDM_octree_seq_t *octree,
 const char       *filename
 )
{

  // count leaves
  int n_leaves = 0;
  for (int i = 0; i < octree->n_nodes; i++) {
    if (octree->nodes[i].is_leaf) n_leaves++;
  }

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "octree_seq\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*octree->n_nodes);//n_leaves);
  for (int inode = 0; inode < octree->n_nodes; inode++) {
    if (1) {//octree->nodes[inode].is_leaf) {
      double *ext = octree->nodes[inode].extents;
      for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
          for (int i = 0; i < 2; i++) {
            int ii = (1-j)*i + j*(1-i);
            fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
          }
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", octree->n_nodes, 9*octree->n_nodes);//n_leaves, 9*n_leaves);
  int ileaf = 0;
  for (int inode = 0; inode < octree->n_nodes; inode++) {
    if (1) {//octree->nodes[inode].is_leaf) {
      fprintf(f, "8 ");
      for (int j = 0; j < 8; j++) {
        fprintf(f, "%d ", 8*ileaf+j);
      }
      fprintf(f, "\n");
      ileaf++;
    }
  }

  fprintf(f, "CELL_TYPES %d\n", octree->n_nodes);//octree->n_nodesn_leaves);
  for (int i = 0; i < octree->n_nodes; i++) {//n_leaves; i++) {
    fprintf(f, "%d\n", 12);
  }

  fclose(f);
}



void PDM_octree_seq_write_octants2
(
 PDM_octree_seq_t *octree,
 const char       *filename
 )
{

  _l_octant_t *octants = octree->octants;

  // count leaves
  int n_leaves = 0;
  for (int i = 0; i < octree->n_nodes; i++) {
    if (octants->is_leaf[i]) n_leaves++;
  }

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "octree_seq\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*octree->n_nodes);//n_leaves);
  for (int inode = 0; inode < octree->n_nodes; inode++) {
    if (1) {//octants->is_leaf[inode]) {
      double *ext = octants->extents + 6*inode;
      for (int k = 0; k < 2; k++) {
        for (int j = 0; j < 2; j++) {
          for (int i = 0; i < 2; i++) {
            int ii = (1-j)*i + j*(1-i);
            fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
          }
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", octree->n_nodes, 9*octree->n_nodes);//n_leaves, 9*n_leaves);
  int ileaf = 0;
  for (int inode = 0; inode < octree->n_nodes; inode++) {
    if (1) {//octants->is_leaf[inode]) {
      fprintf(f, "8 ");
      for (int j = 0; j < 8; j++) {
        fprintf(f, "%d ", 8*ileaf+j);
      }
      fprintf(f, "\n");
      ileaf++;
    }
  }

  fprintf(f, "CELL_TYPES %d\n", octree->n_nodes);//octree->n_nodesn_leaves);
  for (int i = 0; i < octree->n_nodes; i++) {//n_leaves; i++) {
    fprintf(f, "%d\n", 12);
  }

  fclose(f);
}




PDM_octree_seq_shm_t*
PDM_octree_make_shared
(
  PDM_octree_seq_t* local_octree,
  PDM_MPI_Comm      comm_shared
)
{
  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_octree_seq_shm_t* shm_octree = malloc(sizeof(PDM_octree_seq_shm_t));

  shm_octree->comm_shared = comm_shared;
  shm_octree->octrees     = malloc(n_rank_in_shm * sizeof(PDM_octree_seq_t));

  /*
   * Exchange size
   */
  int s_shm_data_in_rank[3] = {0};
  s_shm_data_in_rank[0] = local_octree->n_nodes;
  s_shm_data_in_rank[1] = local_octree->n_nodes_max;
  s_shm_data_in_rank[2] = local_octree->t_n_points;
  int *s_shm_data_in_all_nodes = malloc(3 * n_rank_in_shm * sizeof(int));

  PDM_MPI_Allgather(s_shm_data_in_rank     , 3, PDM_MPI_INT,
                    s_shm_data_in_all_nodes, 3, PDM_MPI_INT, comm_shared);

  int *octants_n_nodes_idx = malloc((n_rank_in_shm+1) * sizeof(int));
  int *n_nodes_max_idx     = malloc((n_rank_in_shm+1) * sizeof(int));
  int *t_n_points_idx      = malloc((n_rank_in_shm+1) * sizeof(int));
  octants_n_nodes_idx[0] = 0;
  n_nodes_max_idx    [0] = 0;
  t_n_points_idx     [0] = 0;
  for(int i = 0; i < n_rank_in_shm; ++i) {
    octants_n_nodes_idx[i+1] = octants_n_nodes_idx[i] + s_shm_data_in_all_nodes[3*i  ];
    n_nodes_max_idx    [i+1] = n_nodes_max_idx    [i] + s_shm_data_in_all_nodes[3*i+1];
    t_n_points_idx     [i+1] = t_n_points_idx     [i] + s_shm_data_in_all_nodes[3*i+2];
  }

  int n_nodes_shared_tot     = octants_n_nodes_idx[n_rank_in_shm];
  int n_nodes_max_shared_tot = n_nodes_max_idx    [n_rank_in_shm];
  int t_n_points_shared_tot  = t_n_points_idx     [n_rank_in_shm];

  // octants->ancestor_id          = realloc(octants->ancestor_id,          sizeof(int                   ) * octree->n_nodes2_max);
  // octants->is_leaf              = realloc(octants->is_leaf,              sizeof(int                   ) * octree->n_nodes2_max);
  // octants->location_in_ancestor = realloc(octants->location_in_ancestor, sizeof(PDM_octree_seq_child_t) * octree->n_nodes2_max);
  // octants->depth                = realloc(octants->depth,                sizeof(int                   ) * octree->n_nodes2_max);
  // octants->children_id          = realloc(octants->children_id,          sizeof(int                   ) * octree->n_nodes2_max * 8);
  // octants->range                = realloc(octants->range,                sizeof(int                   ) * octree->n_nodes2_max * 2);
  // octants->idx                  = realloc(octants->idx,                  sizeof(int                   ) * octree->n_nodes2_max * 9);
  // octants->n_points             = realloc(octants->n_points,             sizeof(int                   ) * octree->n_nodes2_max);
  // octants->extents              = realloc(octants->extents,              sizeof(double                ) * octree->n_nodes2_max * 6);

  free(s_shm_data_in_all_nodes);
  free(octants_n_nodes_idx);
  free(n_nodes_max_idx    );
  free(t_n_points_idx     );

  return shm_octree;
}

void
PDM_octree_seq_shm_free
(
 PDM_octree_seq_shm_t* shm_octree
)
{

  free(shm_octree->octrees);
  free(shm_octree);
}


#ifdef	__cplusplus
}
#endif
