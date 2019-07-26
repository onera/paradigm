
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_handles.h"
#include "pdm_para_octree.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _l_octant_t
 * \brief  Define a list of octants
 *
 */

/* typedef struct  { */

/*   PDM_morton_code_t code; /\*!< morton code *\/ */

/*   int  n_points;          /\*!< Number of points in octant*\/ */
/*   int  range;             /\*!< Start index of point list for each octant *\/ */
/*   int  is_leaf;           /\*!< IS a leaf >*\/ */

/* } _octant_t; */


/**
 * \struct _octant_t
 * \brief  Define an octant
 *
 */

typedef struct  {

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */

  PDM_morton_code_t *codes;        /*!< Morton codes */

  int  *n_points;          /*!< Number of points in octant*/
  int  *range;             /*!< Start index of point list for each octant */
  int  *is_leaf;           /*!< IS a leaf >*/

  int   *neighbor_idx;
  int   *neighbors;               /*!< rank + id_node size = 2 * n_nodes */
  int   dim;

} _l_octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 *
 */

typedef struct  {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /* Translation for the normalization */
  double      d[3];           /* Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t    t_n_points;         /*!< total number of points */
  int            n_points;           /*!< Number of points in each cloud */
  double *points;                    /*!< Point coordinates */
  int *points_icloud;                /*!< Point cloud */
  PDM_g_num_t *points_gnum;          /*!< Point global number */
  PDM_morton_code_t  *points_code;   /*!< Morton codes */

  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

} _octree_t;



/**
 * \struct _neighbors_tmp_t
 * \brief  Define a temporary neighbor structure
 *
 */


typedef struct  {

  int n_neighbor[6];     /*!< Number of neighbors in the arrays  */
  int s_neighbor[6];     /*!< Size of arrays */
  int *neighbors[6];     /*!< Arrays */

} _neighbors_tmp_t;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees    = NULL;

static const double _eps_default  = 1.e-12;

static const int max_morton_level = 31;

/*============================================================================
 * Private function definitions
 *============================================================================*/



/**
 *
 * \brief Neighbor
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static PDM_para_octree_direction_t
_inv_direction
(
 PDM_para_octree_direction_t direc
)
{
  if (direc == PDM_BOTTOM) {
    return PDM_UP;
  }
  else if (direc == PDM_UP) {
    return PDM_BOTTOM;
  }
  else if (direc == PDM_SOUTH) {
    return PDM_NORTH;
  }
  else if (direc == PDM_NORTH) {
    return PDM_SOUTH;
  }
  else if (direc == PDM_WEST) {
    return PDM_EAST;
  }
  else if (direc == PDM_EAST) {
    return PDM_WEST;
  }
  else {
    abort();
  }

}

/**
 *
 * \brief Neighbor
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static PDM_morton_code_t *
_neighbour
(
 PDM_morton_code_t code,
 PDM_para_octree_direction_t direction
)
{
  const int dim = direction / 2;
  const int _direction = 2 * (direction % 2) - 1;

  PDM_morton_code_t *neighbour = NULL;

  if (((_direction > 0) && (code.X[dim] < (code.L - 1))) ||
      ((_direction < 0) && (code.X[dim] > 0))) {

    neighbour = malloc(sizeof(PDM_morton_code_t));

    neighbour->L = code.L;
    neighbour->X[0] = code.X[0];
    neighbour->X[1] = code.X[1];
    neighbour->X[2] = code.X[2];

    neighbour->X[dim] = code.X[dim] + _direction;
  }

  return neighbour;
}


/**
 *
 * \brief Free octants
 *
 * \param [inout]   octants     Octants
 *
 * \return NULL
 *
 */

static _l_octant_t *
_octants_free
(
 _l_octant_t *octants
)
{
  octants->n_nodes_max = 0;
  octants->n_nodes     = 0;

  if (octants->codes != NULL) {
    free (octants->codes);
  }

  if (octants->n_points != NULL) {
    free (octants->n_points);
  }

  if (octants->is_leaf != NULL) {
    free (octants->is_leaf);
  }

  if (octants->range != NULL) {
    free (octants->range);
  }

  if (octants->neighbor_idx != NULL) {
    free (octants->neighbor_idx);
  }

  if (octants->neighbors != NULL) {
    free (octants->neighbors);
  }

  free(octants);
  return NULL;
}


/**
 *
 * \brief Initialize list of octants
 *
 * \param [inout]   octants     Octants
 * \param [in]      octant_dim  Dimension of an octant
 * \param [in]      init_size   Initial size of octants
 *
 */

static void
_octants_init
(
 _l_octant_t *octants,
 const int   octant_dim,
 const int   init_size
)
{
  octants->n_nodes_max = init_size;
  octants->n_nodes     = 0;

  octants->codes    = malloc (sizeof(PDM_morton_code_t) * octants->n_nodes_max);
  octants->n_points = malloc (sizeof(int) * octants->n_nodes_max);
  octants->range    = malloc (sizeof(int) * (octants->n_nodes_max+1));
  octants->is_leaf  = malloc (sizeof(int) * octants->n_nodes_max);

  octants->neighbor_idx = NULL;
  octants->neighbors    = NULL;
  octants->dim = octant_dim;
}


/**
 *
 * \brief Check size of the size of a list of octants
 *
 * \param [in]   octants       Octants
 * \param [in]   n_free_node   Number of required fre nodes
 *
 */

static void
_octants_check_alloc
(
 _l_octant_t *octants,
 const int n_free_node
)
{
  if (octants->n_nodes + n_free_node > octants->n_nodes_max) {

    octants->n_nodes_max *= 2;

    octants->codes    = realloc (octants->codes,
                                 sizeof(PDM_morton_code_t) * octants->n_nodes_max);
    octants->n_points = realloc (octants->n_points,
                                 sizeof(int) * octants->n_nodes_max);
    octants->range = realloc (octants->range,
                              sizeof(int) * (octants->n_nodes_max+1));
    octants->is_leaf = realloc (octants->is_leaf,
                                sizeof(int) * octants->n_nodes_max);

    octants->neighbor_idx = NULL;
    octants->neighbors    = NULL;

  }
}


/**
 *
 * \brief Push back a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_back
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range,
 const int is_leaf
)
{

  _octants_check_alloc (octants, 1);

  const int idx = octants->n_nodes;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->is_leaf[idx] = is_leaf;

  octants->n_nodes += 1;

}


/**
 *
 * \brief Push front a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_push_front
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range,
 const int is_leaf
)
{

  _octants_check_alloc (octants, 1);

  for (int i = octants->n_nodes; i > 0; i--) {

    PDM_morton_copy (octants->codes[i - 1], octants->codes + i);

    octants->n_points[i] =  octants->n_points[i-1];

    octants->range[i] = octants->range[i-1];

    octants->is_leaf[i] = octants->is_leaf[i-1];

  }

  const int idx = 0;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->is_leaf[idx] = is_leaf;

  octants->n_nodes += 1;

}


/**
 *
 * \brief Push back a octant to a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_replace_node_by_child
(
 _l_octant_t *octants,
 const int points_in_leaf_max,
 const int node_id,
 const int n_stored_points,
 const PDM_morton_code_t  *stored_points_code,
 _neighbors_tmp_t **neighbors_tmp
)
{
  const int dim = 3;
  const int n_child = 8;
  const int n_direction = 6;

  _neighbors_tmp_t *_neighbors_tmp = *neighbors_tmp;

  if ((octants->codes[node_id].L < max_morton_level) &&
      !octants->is_leaf[node_id]) {

    _octants_check_alloc (octants, 8);

    if (octants->n_nodes >= octants->n_nodes_max) {
      int pre_nodes_max = octants->n_nodes_max;
      octants->n_nodes_max *= 2;

      *neighbors_tmp =
        realloc (neighbors_tmp, sizeof(_neighbors_tmp_t) * octants->n_nodes_max);
      _neighbors_tmp = *neighbors_tmp;

      for (int i = pre_nodes_max; i < octants->n_nodes_max; i++) {
        for (int j = 0; j < n_direction; j++) {
          _neighbors_tmp[i].n_neighbor[j] = 0;
          _neighbors_tmp[i].s_neighbor[j] = 1;
          _neighbors_tmp[i].neighbors[j] =
            malloc (sizeof(int) * _neighbors_tmp[i].s_neighbor[j]);
        }
      }
    }

    _neighbors_tmp_t cp_neighbors;

    for (int i = 0; i < n_direction; i++) {
      cp_neighbors.s_neighbor[i] = _neighbors_tmp[node_id].s_neighbor[i];
      cp_neighbors.n_neighbor[i] = _neighbors_tmp[node_id].n_neighbor[i];
      cp_neighbors.neighbors[i] = malloc(sizeof(int) *  cp_neighbors.s_neighbor[i]);
      for (int j = 0; j < _neighbors_tmp[node_id].n_neighbor[i]; j++) {
        cp_neighbors.neighbors[i][j] =  _neighbors_tmp[node_id].neighbors[i][j];
      }
    }

    PDM_morton_code_t children[n_child];
    PDM_morton_get_children(dim,
                            octants->codes[node_id],
                            children);

    const int step = dim - 1;
    for (int i = octants->n_nodes - 1; i > node_id; i--) {
      PDM_morton_copy (octants->codes[i], octants->codes + step + i);
      octants->n_points[step+i] =  octants->n_points[i];
      octants->range[step+i] = octants->range[i];
      octants->is_leaf[step+i] = octants->is_leaf[i];
      for (int j = 0; j < n_direction; j++) {
        _neighbors_tmp[step+i].n_neighbor[j] = _neighbors_tmp[i].n_neighbor[j];
        _neighbors_tmp[step+i].s_neighbor[j] = _neighbors_tmp[i].s_neighbor[j];
        _neighbors_tmp[step+i].neighbors[j] = _neighbors_tmp[i].neighbors[j];
      }
    }

    const int n_points_node = octants->n_points[node_id];
    const int range_node = octants->range[node_id];

    int range_children[n_child];
    int n_points_children[n_child];

    for (int i = 0; i < dim; i++) {
      range_children[i] =
        PDM_morton_binary_search (n_points_node,
                                  children[i],
                                  (PDM_morton_code_t *) stored_points_code + range_node);
    }

    for (int i = 0; i < dim - 1; i++) {
      n_points_children[i] = range_children[i+1] - range_children[i];
    }

    n_points_children[dim-1] = n_points_node - range_children[dim-1];

    int k = 0;
    for (int i = node_id; i < node_id + n_child; i++) {
      PDM_morton_copy (children[k], octants->codes + step + i);
      octants->n_points[i] = n_points_children[k];
      octants->range[i] = range_children[k];
      if (n_points_children[k] <= points_in_leaf_max) {
        octants->is_leaf[i] = 1;
      }
      else {
        octants->is_leaf[i] = 0;
      }
      k += 1;
    }

    /* Inter children neighborhood */

    int id_1 = node_id;
    int id_2 = node_id + 1;
    int dir_1 = PDM_WEST;
    int dir_2 = PDM_EAST;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id;
    id_2 = node_id + 2;
    dir_1 = PDM_NORTH;
    dir_2 = PDM_SOUTH;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id;
    id_2 = node_id + 4;
    dir_1 = PDM_UP;
    dir_2 = PDM_BOTTOM;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 1;
    id_2 = node_id + 3;
    dir_1 = PDM_NORTH;
    dir_2 = PDM_SOUTH;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 1;
    id_2 = node_id + 5;
    dir_1 = PDM_UP;
    dir_2 = PDM_BOTTOM;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 2;
    id_2 = node_id + 3;
    dir_1 = PDM_EAST;
    dir_2 = PDM_WEST;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 2;
    id_2 = node_id + 6;
    dir_1 = PDM_UP;
    dir_2 = PDM_BOTTOM;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 3;
    id_2 = node_id + 7;
    dir_1 = PDM_UP;
    dir_2 = PDM_BOTTOM;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 4;
    id_2 = node_id + 5;
    dir_1 = PDM_EAST;
    dir_2 = PDM_WEST;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 4;
    id_2 = node_id + 6;
    dir_1 = PDM_NORTH;
    dir_2 = PDM_SOUTH;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 5;
    id_2 = node_id + 7;
    dir_1 = PDM_NORTH;
    dir_2 = PDM_SOUTH;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    id_1 = node_id + 6;
    id_2 = node_id + 7;
    dir_1 = PDM_EAST;
    dir_2 = PDM_WEST;
    _neighbors_tmp[id_1].neighbors[dir_1][_neighbors_tmp[id_1].n_neighbor[dir_1]] = id_2;
    _neighbors_tmp[id_1].n_neighbor[dir_1] += 1;
    _neighbors_tmp[id_2].neighbors[dir_2][_neighbors_tmp[id_2].n_neighbor[dir_2]] = id_1;
    _neighbors_tmp[id_2].n_neighbor[dir_2] += 1;

    /* Extern neighborhood */

    /* PDM_para_octree_direction_t dir_children[8][3] = {{PDM_BOTTOM, PDM_SOUTH, PDM_WEST}, */
    /*                                                   {PDM_BOTTOM, PDM_SOUTH, PDM_EAST}, */
    /*                                                   {PDM_BOTTOM, PDM_NORTH, PDM_WEST}, */
    /*                                                   {PDM_BOTTOM, PDM_NORTH, PDM_EAST}, */
    /*                                                   {PDM_UP, PDM_SOUTH, PDM_WEST}, */
    /*                                                   {PDM_UP, PDM_SOUTH, PDM_EAST}, */
    /*                                                   {PDM_UP, PDM_NORTH, PDM_WEST}, */
    /*                                                   {PDM_UP, PDM_NORTH, PDM_EAST}}; */

    const int dir_children[6][4] = {{0, 1, 2, 3},
                                    {4, 5, 6, 7},
                                    {0, 1, 4, 5},
                                    {2, 3, 6, 7},
                                    {0, 2, 4, 6},
                                    {1, 3, 5, 7}};

    int sel_children[4] = {0, 0, 0, 0};

    for (int i = 0; i < n_direction; i++) {

       for (int j = 0; j < cp_neighbors.n_neighbor[i]; j++) {
        int neighbor_id = cp_neighbors.neighbors[i][j];

        int n_neighbor_child = 0;
        for (k = 0; k < 4; k++) {
          int child_id = node_id + dir_children[i][k];

          PDM_morton_code_t *neighbour_code = _neighbour (children[dir_children[i][k]], i);

          if (neighbor_id >= 0) {

            PDM_morton_compare_t pmc =
              PDM_morton_compare(dim,  *neighbour_code, octants->codes[neighbor_id]);

            if ((pmc == PDM_MORTON_SAME_ANCHOR) || (pmc == PDM_MORTON_EQUAL_ID)) {
              sel_children[n_neighbor_child++] = k;

              _neighbors_tmp[child_id].neighbors[i][_neighbors_tmp[child_id].n_neighbor[i]] =
                neighbor_id;
              _neighbors_tmp[child_id].n_neighbor[i] += 1;
            }

          }
          else {
            _neighbors_tmp[child_id].neighbors[i][_neighbors_tmp[child_id].n_neighbor[i]] =
              neighbor_id;
            _neighbors_tmp[child_id].n_neighbor[i] += 1;
          }

        }

        if (n_neighbor_child > 0) {

          for (k = 0; k < 4; k++) {

            while ((_neighbors_tmp[neighbor_id].n_neighbor[i] + (n_neighbor_child - 1))
                   >= _neighbors_tmp[neighbor_id].s_neighbor[i]) {
              _neighbors_tmp[neighbor_id].s_neighbor[i] *= 2;
              _neighbors_tmp[neighbor_id].neighbors[i] =
                realloc (_neighbors_tmp[neighbor_id].neighbors[i],
                         sizeof(int) *  _neighbors_tmp[neighbor_id].s_neighbor[i]);
            }

            int idx = -1;
            for (int k1 = 0; k1 < _neighbors_tmp[neighbor_id].n_neighbor[i]; k1++) {
              if (_neighbors_tmp[neighbor_id].neighbors[i][k1] == node_id) {
                idx = k1;
                break;
              }
            }

            assert (idx != -1);
            _neighbors_tmp[neighbor_id].neighbors[i][idx] = node_id + sel_children[0];

            for (int k1 = 1; k1 < n_neighbor_child; k1++) {

              _neighbors_tmp[neighbor_id].neighbors[i][_neighbors_tmp[neighbor_id].n_neighbor[i]++] =
                node_id + sel_children[k1];
            }

          }
        }
      }
    }

    octants->n_nodes += n_child - 1;

    for (int i = node_id; i < node_id + n_child; i++) {
      _octants_replace_node_by_child (octants,
                                      points_in_leaf_max,
                                      i,
                                      n_stored_points,
                                      stored_points_code,
                                      neighbors_tmp);
    }
  }
}


/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _octree_t *
_get_from_id
(
 int  id
)
{
  _octree_t *octree = (_octree_t *) PDM_Handles_get (_octrees, id);

  if (octree == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return octree;
}


/**
 *
 * \brief Build minimal octree between two octants
 *
 * \param [in]     octree    Current octree
 * \param [in]     code      Morton code
 * \param [inout]  extents   Extents associated to the Morton code
 *
 */

static void
_extents
(
 _octree_t *octree,
 PDM_morton_code_t code,
 double    extents[]
)
{
  for (int i = 0; i < octree->dim; i++) {
    extents[i] =
      ((double) code.X[i]/((double) pow(2,code.L)))* octree->d[i] + octree->s[i];
    extents[octree->dim + i] =
      (((double) code.X[i] + 1)/((double) pow(2,code.L))) * octree->d[i] + octree->s[i];
  }
}

/**
 *
 * \brief Remove duplicates
 *
 * \param [in]  octants List of octants
 *
 * \return octants without duplicates
 *
 */

static _l_octant_t *
_remove_duplicates
(
 _l_octant_t *octants
)
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  PDM_morton_code_t prev_code;

  prev_code.L = -1;
  prev_code.X[0] = 0;
  prev_code.X[1] = 0;
  prev_code.X[2] = 0;

  for (int i = 0; i < octants->n_nodes; i++) {

    if (_codes[i].L == prev_code.L) {
      if ((prev_code.X[0] == _codes[i].X[0]) &&
          (prev_code.X[1] == _codes[i].X[1]) &&
          (prev_code.X[2] == _codes[i].X[2])) {

        break;
      }
    }

    prev_code.L    = _codes[i].L;
    prev_code.X[0] = _codes[i].X[0];
    prev_code.X[1] = _codes[i].X[1];
    prev_code.X[2] = _codes[i].X[2];

    _octants_push_back (r_octants,
                        _codes[i],
                        0,
                        0,
                        0);

  }

  return r_octants;
}

/**
 *
 * \brief Removing overlaps from a sorted lis of octants
 *
 * \param [inout]  octants A lis of octants
 *
 */

static _l_octant_t *
_linearize
(
 _l_octant_t *octants
)
{
  PDM_morton_code_t *_codes = octants->codes;
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (r_octants, dim, octants->n_nodes);

  for (int i = 0; i < octants->n_nodes - 1; i++) {

    if (!PDM_morton_ancestor_is (_codes[i], _codes[i+1])) {
      _octants_push_back (r_octants,
                          _codes[i],
                          0,
                          0,
                          0);
    }

  }

  _octants_push_back (r_octants,
                      _codes[octants->n_nodes-1],
                      0,
                      0,
                      0);

  return r_octants;
}


/**
 *
 * \brief Constructing a minimal linear octree between two octants
 *
 * \param [in]  a     Morton code a
 * \param [in]  b     Morton code b
 *
 * \return octants The minimal linear octree between a and b
 *
 */

static _l_octant_t *
_complete_region
(
 PDM_morton_code_t a,
 PDM_morton_code_t b
)
{
  _l_octant_t *w_octants = malloc(sizeof(_l_octant_t));
  _l_octant_t *r_octants = malloc(sizeof(_l_octant_t));

  const int dim = 3;

  _octants_init (w_octants, dim, 4);
  _octants_init (r_octants, dim, 4);

  assert (PDM_morton_a_gt_b (b, a));

  PDM_morton_code_t nca;
  PDM_morton_nearest_common_ancestor (a, b, &nca);

  PDM_morton_code_t children[8];
  PDM_morton_get_children(dim,
                          nca,
                          children);

  for (int i = 0; i < 8; i++) {
    _octants_push_back (w_octants,
                        children[i],
                        0,
                        0,
                        0);
  }

  int i1 = 0;
  while (i1 < w_octants->n_nodes) {

    PDM_morton_code_t *_code = w_octants->codes + i1;

    if (PDM_morton_a_gt_b (*_code, a) &&
        PDM_morton_a_gt_b (b, *_code) &&
        !PDM_morton_ancestor_is (*_code, b)) {

      _octants_push_back (r_octants,
                          *_code,
                          0,
                          0,
                          0);
    }

    else if (PDM_morton_ancestor_is (*_code, a) || PDM_morton_ancestor_is (*_code, b)) {

      PDM_morton_get_children(dim,
                              nca,
                              children);

      for (int i = 0; i < 8; i++) {
        _octants_push_back (w_octants,
                            children[i],
                            0,
                            0,
                            0);
      }

    }

    i1 += 1;
  }

  _octants_free (w_octants);

  return r_octants;
}



/**
 *
 * \brief Redistribute octants
 *
 * \param [inout]  L             Distributed list of octants
 * \param [in]     morton_index  Morton index
 * \param [in]     comm          MPI communicator
 *
 */

static void
_distribute_octants
(
 _l_octant_t       *L,
 PDM_morton_code_t *morton_index,
 PDM_MPI_Comm       comm
)
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int *send_count = malloc(sizeof(int) * n_ranks);
  size_t *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  size_t *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  int irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {
    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search (n_ranks - (irank + 1),
                                             L->codes[i],
                                             morton_index + irank + 1);
    }
    send_count[irank] += L->dim + 1;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < L->n_nodes; i++) {

    if (PDM_morton_a_ge_b (L->codes[i], morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                            L->codes[i],
                                            morton_index + irank + 1);
    }

    int shift = send_shift[irank] + send_count[irank];
    send_codes[shift++] = L->codes[i].L;

    for (int j = 0; j < L->dim; j++) {
      send_codes[shift++] = L->codes[i].X[j];
    }

    send_count[irank] += L->dim + 1;
  }

  PDM_morton_int_t * recv_codes = malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv_l(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                      recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                      comm);

  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  /* - tri des codes recus */

  _octants_free (L);
  _octants_init (L, L->dim, recv_shift[n_ranks]/(L->dim + 1));

  int idx = 0;
  for (int i = 0; i < recv_shift[n_ranks]; i++) {
    PDM_morton_code_t _code;
    _code.L = recv_codes[idx++];
    for (int j = 0; j < L->dim; j++) {
     _code.X[j] = recv_codes[idx++];
    }
    _octants_push_back (L,
                        _code,
                        0,
                        0,
                        0);
  }

  free (recv_shift);
  free (recv_codes);

  PDM_morton_local_sort (L->n_nodes, L->codes);

}


/**
 *
 * \brief Constructing a complete linear octree from partial set of octants
 *
 * \param [in]  L     Distributed list of octants
 * \param [in]  comm  MPI Communicator
 *
 * \return octants The complete linear octree
 *
 */

static _l_octant_t *
_complete_octree
(
 _l_octant_t *L,
 PDM_MPI_Comm comm
)
{
  const int dim = 3;

  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_size (comm, &rank);

  /* Remove duplicates */

  _l_octant_t *L1 = _remove_duplicates (L);

  /* Linearize */

  _l_octant_t *L2 = _linearize (L1);

  _octants_free (L1);

  PDM_morton_code_t *L2_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));

  int *order = malloc (sizeof(int) * L2->n_nodes);
  int *weight = malloc (sizeof(int) * L2->n_nodes);

  for (int i = 0; i < L2->n_nodes; i++) {
    weight[i] = 0;
    order[i] = i;
  }

  int max_level = -1;
  for (int i = 0; i < L2->n_nodes; i++) {
    max_level = PDM_MAX (L2->codes[i].L, max_level);
  }

  int max_max_level;
  PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, comm);

  PDM_morton_build_rank_index(dim,
                              max_level,
                              L2->n_nodes,
                              L2->codes,
                              weight,
                              order,
                              L2_morton_index,
                              comm);

  free(weight);
  free(order);

  _distribute_octants (L2, L2_morton_index, comm);

  free (L2_morton_index);

  if (rank == 0) {
    PDM_morton_code_t root_DFD;

    root_DFD.L = max_morton_level;
    root_DFD.X[0] = 0;
    root_DFD.X[1] = 0;
    root_DFD.X[2] = 0;

    PDM_morton_code_t FINA;

    PDM_morton_nearest_common_ancestor (root_DFD,
                                        L2->codes[0],
                                        &FINA);

    PDM_morton_code_t child[8];
    PDM_morton_get_children(dim,
                            FINA,
                            child);

    _octants_push_front (L2,
                         child[0],
                         0,
                         0,
                         0);
  }

  if (rank == n_ranks - 1) {
    PDM_morton_code_t root_DLD;

    root_DLD.L = max_morton_level;
    root_DLD.X[0] = 1 << max_morton_level;
    root_DLD.X[1] = 1 << max_morton_level;
    root_DLD.X[2] = 1 << max_morton_level;

    PDM_morton_code_t FINA;
    PDM_morton_nearest_common_ancestor (root_DLD,
                                        L2->codes[0],
                                        &FINA);

    PDM_morton_code_t child[8];
    PDM_morton_get_children(dim,
                            FINA,
                            child);

    _octants_push_back (L2,
                        child[7],
                        0,
                        0,
                        0);

  }

  unsigned int sbuff[4];
  unsigned int rbuff[4];
  PDM_MPI_Request srequest;
  PDM_MPI_Request rrequest;

  if (rank < n_ranks - 1) {

    PDM_MPI_Irecv ((void *) rbuff,
                   4,
                   PDM_MPI_UNSIGNED,
                   rank+1,
                   0,
                   comm,
                   &rrequest);

  }

  if (rank > 0) {
    sbuff[0] = L2->codes[0].L;
    sbuff[1] = L2->codes[0].X[0];
    sbuff[2] = L2->codes[0].X[1];
    sbuff[3] = L2->codes[0].X[2];

    PDM_MPI_Issend ((void *) sbuff,
                    4,
                    PDM_MPI_UNSIGNED,
                    rank-1,
                    0,
                    comm,
                    &srequest);


  }

  if (rank < n_ranks - 1) {

    PDM_MPI_Wait (&rrequest);
    PDM_morton_code_t code;

    code.L = rbuff[0];
    code.X[0] = rbuff[1];
    code.X[1] = rbuff[2];
    code.X[2] = rbuff[3];

    _octants_push_back (L2,
                        code,
                        0,
                        0,
                        0);

  }

  if (rank > 0) {

    PDM_MPI_Wait (&srequest);

  }

  _l_octant_t *R = malloc(sizeof(_l_octant_t));
  _octants_init (R, dim, L2->n_nodes);

  for (int i = 0; i < L2->n_nodes - 1; i++) {
    _l_octant_t *A = _complete_region (L2->codes[i], L2->codes[i+1]);

    _octants_push_back (R,
                        L2->codes[i],
                        0,
                        0,
                        0);

    for (int j = 0; j < A->n_nodes - 1; j++) {
      _octants_push_back (R,
                          A->codes[j],
                          0,
                          0,
                          0);
    }
    _octants_free (A);
  }

  if (rank == n_ranks - 1) {
    _octants_push_back (R,
                        L2->codes[L2->n_nodes-1],
                        0,
                        0,
                        0);
  }

  _octants_free (L2);

  return R;
}



/**
 *
 * \brief Distribute points
 *
 * \param [in]   id                 Identifier
 *
 */

static void
_distribute_points
(
 int *n_points,
 double **points,
 int **points_icloud,
 PDM_g_num_t **points_gnum,
 PDM_morton_code_t **points_code,
 PDM_morton_code_t *morton_index,
 const PDM_MPI_Comm comm,
 const int dim,
 const int max_level,
 const double *global_extents
)
{
  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  int _n_points = *n_points;

  double *__points = *points;
  int *__points_icloud = *points_icloud;
  PDM_g_num_t *__points_gnum = *points_gnum;
  PDM_morton_code_t *__points_code = *points_code;

  int *c_rank = malloc (_n_points * sizeof(int));

  for (int i = 0; i < _n_points; i++) {
    size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                __points_code[i],
                                                morton_index);
    c_rank[i] = (int) _c_rank;
  }

  int *send_count = malloc (n_ranks * sizeof (int));
  int *recv_count = malloc (n_ranks * sizeof (int));
  int *send_shift = malloc ((n_ranks + 1) * sizeof (int));
  int *recv_shift = malloc ((n_ranks + 1) * sizeof (int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    send_count[c_rank[i]] += dim;
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  double *send_coords = malloc (send_shift[n_ranks] * sizeof(double));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++)
    send_count[rank_id] = 0;

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    for (int j = 0; j < dim; j++)
      send_coords[shift + j] = __points[i*dim + j];
    send_count[rank_id] += dim;
  }

  double *recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

  /* Exchange coords between processes */

  PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                    recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                    comm);

  free(send_coords);

  /* Build send and receive buffers */

  for (int rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
    send_shift[rank_id] = send_shift[rank_id]/dim;
    recv_shift[rank_id] = recv_shift[rank_id]/dim;
  }

  int *send_points_icloud = malloc (send_shift[n_ranks] * sizeof(int));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    recv_count[rank_id] = recv_count[rank_id]/dim;
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_icloud[shift] = __points_icloud[i];
    send_count[rank_id] += 1;
  }

  int *recv_points_icloud = malloc (recv_shift[n_ranks] * sizeof(int));

  /* Exchange points_icloud between processes */

  PDM_MPI_Alltoallv(send_points_icloud, send_count, send_shift, PDM_MPI_INT,
                    recv_points_icloud, recv_count, recv_shift, PDM_MPI_INT,
                    comm);

  free(send_points_icloud);


  /* Build send and receive buffers : points_gnum*/

  PDM_g_num_t *send_points_gnum =
    malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  for (int i = 0; i < _n_points; i++) {
    int rank_id = c_rank[i];
    int shift = send_shift[rank_id] + send_count[rank_id];
    send_points_gnum[shift] = __points_gnum[i];
    send_count[rank_id] += 1;
  }

  free (c_rank);

  PDM_g_num_t *recv_points_gnum =
    malloc (recv_shift[n_ranks] * sizeof(PDM_g_num_t));

  /* Exchange points_gnum between processes */

  PDM_MPI_Alltoallv(send_points_gnum, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                    recv_points_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                    comm);

  free(send_points_gnum);

  _n_points = recv_shift[n_ranks];

  free (send_count);
  free (recv_count);
  free (send_shift);
  free (recv_shift);

  __points = realloc (__points, sizeof(double) * _n_points);

  __points_icloud =
    realloc (__points_icloud, sizeof(int) * _n_points);

  __points_gnum =
    realloc (__points_gnum, sizeof(PDM_g_num_t) * _n_points);

  /* Re-encode points */

  __points_code = realloc (__points_code,
                                 sizeof(PDM_morton_code_t) * _n_points);

  double d[3];
  double s[3];

  PDM_morton_encode_coords(dim,
                           max_level,
                           global_extents,
                           _n_points,
                           __points,
                           __points_code,
                           d,
                           s);

  int *order = malloc (sizeof(int) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    order[i] = i;
  }

  PDM_morton_local_order(_n_points, __points_code, order);

  for (int i = 0; i < _n_points; i++) {
    __points_icloud[i] = recv_points_icloud[order[i]];
    __points_gnum[i] = recv_points_gnum[order[i]];
    for (int j = 0; j < dim; j++) {
      __points[dim*i+j] = recv_coords[dim*order[i]+j];
    }
  }

  free (recv_points_icloud);
  free (recv_points_gnum);
  free (recv_coords);

  PDM_morton_code_t *_points_code =
    malloc (sizeof(PDM_morton_code_t) * _n_points);

  for (int i = 0; i < _n_points; i++) {
    _points_code[i].L = __points_code[order[i]].L;
    _points_code[i].X[0] = __points_code[order[i]].X[0];
    _points_code[i].X[1] = __points_code[order[i]].X[1];
    _points_code[i].X[2] = __points_code[order[i]].X[2];
  }

  free (__points_code);
  free (order);

  *points_code = _points_code;

  *points = __points;
  *points_icloud = __points_icloud;
  *points_gnum = __points_gnum;

  *n_points = _n_points;
}


/**
 *
 * \brief Partitioning octants into large contiguous blocks. The list of octants
 *        is redistributed
 *
 * \param [in]  octant_list  a list of distributed octants,
 *                           octant_list is not redistributed at the end
 *
 * \return block_octants  A list of distributed blocks
 *
 */

static _l_octant_t *
_block_partition
(
 _l_octant_t *octant_list,
 const PDM_MPI_Comm comm,
 PDM_morton_code_t **G_morton_index
)
{

  /* Complete region */

  _l_octant_t *T = _complete_region (octant_list->codes[0],
                                     octant_list->codes[octant_list->n_nodes]);

  int max_level = -1;
  int min_level = 31;
  for (int i = 0; i < octant_list->n_nodes; i++) {
    max_level = PDM_MAX (octant_list->codes[i].L, max_level);
    min_level = PDM_MIN (octant_list->codes[i].L, min_level);
  }

  int max_max_level;
  PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, comm);

  /* Complete octree */

  _l_octant_t C;

  _octants_init (&C, octant_list->dim, octant_list->n_nodes);

  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (octant_list->codes[i].L <= min_level) {
      _octants_push_back (&C,
                          octant_list->codes[i],
                          octant_list->n_points[i],
                          octant_list->range[i],
                          octant_list->is_leaf[i]);
    }
  }

  _octants_free (&C);

  _l_octant_t *G = _complete_octree (&C, comm);

  _octants_free (T);

  /*
   * Compute weight
   */

  /* - exchange codes to ranks (weight per rank)*/

  int n_ranks;
  PDM_MPI_Comm_size(comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank(comm, &rank);

  int *code_buff = malloc (sizeof(int) * (octant_list->dim + 1));
  int *rank_buff = malloc (sizeof(int) * n_ranks * (octant_list->dim + 1));
  code_buff[0] = G->codes[0].L;

  for (int i = 0; i < octant_list->dim; i++) {
    code_buff[i+1] =  G->codes[0].X[i];
  }

  PDM_MPI_Allgather (code_buff, octant_list->dim + 1, PDM_MPI_INT,
                     rank_buff, octant_list->dim + 1, PDM_MPI_INT,
                     comm);

  PDM_morton_code_t *rank_codes = malloc (sizeof(PDM_morton_code_t) * n_ranks);

  for (int i = 0; i < n_ranks; i++) {
    rank_codes[i].L = rank_buff[(octant_list->dim + 1) * i];
    for (int j = 0; j < octant_list->dim; j++) {
      rank_codes[i].X[j] = rank_buff[(octant_list->dim + 1) * i + j];
    }
  }

  free (code_buff);

  int *send_count = malloc(sizeof(int) * n_ranks);
  int *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  int *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  int irank = 0;
  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }
    send_count[irank] += send_shift[irank] + (octant_list->dim + 1);
  }

  /* Exchange number of coords to send to each process */

  PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT,
                   recv_count, 1, PDM_MPI_INT, comm);

  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
    recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
  }

  /* Build send and receive buffers */

  PDM_morton_int_t *send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }

    int shift = send_shift[irank] + send_count[irank];
    send_codes[shift++] = octant_list->codes[i].L;

    for (int j = 0; j < octant_list->dim; j++) {
      send_codes[shift++] = octant_list->codes[i].X[j];
    }

    send_count[irank] += octant_list->dim + 1;
  }

  PDM_morton_int_t *recv_codes =
    malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */

  PDM_MPI_Alltoallv(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                    recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                    comm);


  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  const int n_recv_codes = recv_shift[n_ranks] / (1+octant_list->dim);

  free (recv_shift);

  int *weight = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i < G->n_nodes; i++) {
    weight[i] = 0;
  }

  /* - compute weight of each cell */

  const int _stride = octant_list->dim + 1;

  for (int i = 0; i < n_recv_codes; i++) {

    PDM_morton_code_t code;

    code.L = recv_codes[i*_stride];

    for (int j = 0; j < _stride-1; j++) {
      code.X[j] = recv_codes[i*_stride+j+1];
    }

    int G_node =  PDM_morton_binary_search(G->n_nodes,
                                           code,
                                           G->codes);

    weight[G_node] += octant_list->n_points[i];
  }

  free (recv_codes);

  /*
   * Load balancing G from weight
   */

  int *order = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i <  G->n_nodes; i++) {
    order[i] = i;
  }

  *G_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));
  PDM_morton_code_t *_G_morton_index = *G_morton_index;

  PDM_morton_build_rank_index (octant_list->dim,
                               max_max_level,
                               G->n_nodes,
                               G->codes,
                               weight,
                               order,
                               _G_morton_index,
                               comm);

  free (order);
  free (weight);

  _distribute_octants (G, _G_morton_index, comm);

  /*
   * Redistribute octant list from coarse load balancing
   */

  _distribute_octants (octant_list, _G_morton_index, comm);

  return G;

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
 const PDM_MPI_Comm comm
)
{

  if (_octrees == NULL) {
    _octrees = PDM_Handles_create (4);
  }

  _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

  int id = PDM_Handles_store (_octrees, octree);

  octree->dim = 3;

  for (int i = 0; i < octree->dim; i++) {
    octree->global_extents[i]   = -HUGE_VAL;
    octree->global_extents[octree->dim+i] =  HUGE_VAL;
    octree->s[i]         = 0.;
    octree->d[i]         = 0.;
  }

  octree->depth_max = depth_max;
  octree->points_in_leaf_max = points_in_leaf_max;

  octree->n_point_clouds = n_point_cloud;
  octree->t_n_points = 0;
  octree->n_points = 0;
  octree->points = NULL;
  octree->points_icloud = NULL;
  octree->points_gnum = NULL;
  octree->points_code = NULL;

  octree->octants = NULL;

  octree->comm = comm;

  return id;

}


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
)
{
  _octree_t *octree = _get_from_id (id);

  if (octree->points != NULL) {
    free (octree->points);
  }

  if (octree->points_icloud != NULL) {
    free (octree->points_icloud);
  }

  if (octree->points_gnum != NULL) {
    free (octree->points_gnum);
  }

  if (octree->points_code != NULL) {
    free (octree->points_code);
  }

  free (octree->octants);

  free (octree);

  PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

  const int n_octrees = PDM_Handles_n_get (_octrees);

  if (n_octrees == 0) {
    _octrees = PDM_Handles_free (_octrees);
  }
}


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
)
{
  _octree_t *octree = _get_from_id (id);

  const int idx = octree->n_points;

  octree->n_points += n_points;
  octree->points =
    realloc (octree->points, octree->n_points * sizeof(double) * octree->dim);
  octree->points_icloud =
    realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum =
    realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code =
    realloc (octree->points_code, octree->n_points * sizeof(PDM_morton_code_t));

  for (int i = 0; i < octree->dim * n_points; i++) {
    octree->points[octree->dim*idx + i] = coords[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_gnum[idx + i] = g_num[i];
  }

  for (int i = 0; i < n_points; i++) {
    octree->points_icloud[idx + i] = i_point_cloud;
  }

}



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
 const int  id
)
{
  _octree_t *octree = _get_from_id (id);

  const int dim = octree->dim;
  const int max_level = sizeof(PDM_morton_int_t)*8 - 1;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  /*
   * Get coord extents
   */

  PDM_morton_get_coord_extents(dim,
                               octree->n_points,
                               octree->points,
                               octree->global_extents,
                               octree->comm);

  /*
   * Encode coords
   */

  PDM_morton_encode_coords(dim,
                           max_level,
                           octree->global_extents,
                           octree->n_points,
                           octree->points,
                           octree->points_code,
                           octree->d,
                           octree->s);

  int *order = malloc (sizeof(int) * octree->n_points);

  for (int i = 0; i < octree->n_points; i++) {
    order[i] = i;
  }

  /**************************************
   *
   * Global order of codes and balancing
   *
   **************************************/

  PDM_morton_local_order (octree->n_points,
                          octree->points_code,
                          order);

  if (n_ranks > 1) {

    int *weight = malloc (sizeof(int) * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      weight[i] = 1;
    }

    PDM_morton_code_t *morton_index =
      malloc (sizeof(PDM_morton_code_t) * (n_ranks + 1));

    PDM_morton_build_rank_index(dim,
                                max_level,
                                octree->n_points,
                                octree->points_code,
                                weight,
                                order,
                                morton_index,
                                octree->comm);

    free (weight);
    free (order);

    /* distribute point from morton_index */

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        morton_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    free(morton_index);

  }

  else {

    int *_points_icloud = malloc (sizeof(int) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_icloud[i] =  octree->points_icloud[order[i]];
    }

    free (octree->points_icloud);
    octree->points_icloud = _points_icloud;

    PDM_g_num_t *_points_gnum = malloc (sizeof(PDM_g_num_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_gnum[i] =  octree->points_gnum[order[i]];
    }

    free (octree->points_gnum);
    octree->points_gnum = _points_gnum;

    PDM_morton_code_t *_points_code =
      malloc (sizeof(PDM_morton_code_t) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      _points_code[i].L = octree->points_code[order[i]].L;
      _points_code[i].X[0] = octree->points_code[order[i]].X[0];
      _points_code[i].X[1] = octree->points_code[order[i]].X[1];
      _points_code[i].X[2] = octree->points_code[order[i]].X[2];
    }

    free (octree->points_code);
    octree->points_code = _points_code;

    double *_points = malloc (sizeof(double) * dim * octree->n_points);
    for (int i = 0; i < octree->n_points; i++) {
      for (int j = 0; j < dim; j++) {
        _points[dim*i+j] = octree->points[dim*order[i]+j];
      }
    }
    free (octree->points_code);
    free (order);
  }

  PDM_morton_code_t *block_octants_index = NULL;
  if (n_ranks > 1) {

    /*************************************************************************
     *
     * Store points in the octants (leaves) at the maximum depth of the octree
     * to build
     *
     *************************************************************************/

    int chg_code = 1;
    _l_octant_t *point_octants = malloc(sizeof(_l_octant_t));

    int curr_node = -1;

    _octants_init (point_octants, octree->n_points, octree->dim);

    for (int i = 0; i < octree->n_points; i++) {

      PDM_morton_code_t _point_code;
      PDM_morton_copy (octree->points_code[i], &_point_code);

      PDM_morton_assign_level (&_point_code, octree->depth_max);

      if (curr_node != -1) {
        chg_code = !(PDM_morton_a_eq_b (point_octants->codes[curr_node],
                                        _point_code));
      }

      if (chg_code) {

        _octants_check_alloc (point_octants, 1);

        int idx = point_octants->n_nodes;

        curr_node = idx;

        PDM_morton_copy (octree->points_code[i], &(point_octants->codes[idx]));

        point_octants->is_leaf[idx] = 1;
        point_octants->n_points[idx] = 1;
        point_octants->range[idx] = i;

        point_octants->n_nodes += 1;
      }

      else {
        point_octants->n_points[curr_node] += 1;
      }

    }

    /*************************************************************************
     *
     * Block partition (algo 2 sundar)
     *
     *************************************************************************/

    octree->octants = _block_partition (point_octants,
                                        octree->comm,
                                        &block_octants_index);

    /*************************************************************************
     *
     * Redistribute points
     *
     *************************************************************************/

    _distribute_points (&octree->n_points,
                        &octree->points,
                        &octree->points_icloud,
                        &octree->points_gnum,
                        &octree->points_code,
                        block_octants_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    int iblock = 0;
    for (int i = 0; i < octree->n_points; i++) {
      if (iblock < (octree->octants->n_nodes - 1)) {
        if (PDM_morton_a_ge_b (octree->points_code[i], octree->octants->codes[iblock+1])) {

          iblock = iblock + 1 + PDM_morton_binary_search (octree->octants->n_nodes - (iblock + 1),
                                                          octree->points_code[i],
                                                          octree->octants->codes + iblock + 1);
        }
      }
      octree->octants->n_points[iblock] += 1;
    }

    octree->octants->range[0] = 0;

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      octree->octants->range[i+1] =
        octree->octants->range[i] +
        octree->octants->n_points[i];
    }


  }

  else {

    octree->octants = malloc(sizeof(_l_octant_t));

    _octants_init (octree->octants, octree->dim, octree->n_points);

    PDM_morton_code_t code;

    code.L = 0;
    code.X[0] = 0;
    code.X[1] = 0;
    code.X[2] = 0;

    int is_leaf = 0;
    if (octree->n_points <= octree->points_in_leaf_max) {
      is_leaf = 1;
    }

    _octants_push_back (octree->octants,
                        code,
                        octree->n_points,
                        0,
                        is_leaf);

    octree->octants->range[0] = 0;
    octree->octants->range[1] = octree->n_points;
    octree->octants->n_points[0] = octree->n_points;

  }

  /*************************************************************************
   *
   * Add child (while n points > N points max)
   *     - Build neighbor
   *
   *************************************************************************/

  const int n_direction = 6;

  _neighbors_tmp_t *neighbors_tmp = malloc(sizeof(_neighbors_tmp_t) * octree->octants->n_nodes_max);

  for (int i = 0; i < octree->octants->n_nodes_max; i++) {
    for (int j = 0; j < n_direction; j++) {
      neighbors_tmp[i].n_neighbor[j] = 0;
      neighbors_tmp[i].s_neighbor[j] = 1;
      neighbors_tmp[i].neighbors[j] = malloc (sizeof(int) * neighbors_tmp[i].s_neighbor[j]);
    }
  }

  if (n_ranks > 1) {

    int n_init_node = octree->octants->n_nodes;

    /* Initialize neighbor */

    for (int i = 0; i < octree->octants->n_nodes; i++) {

      for (PDM_para_octree_direction_t j = 0; j < n_direction; j++) {

        PDM_morton_code_t *neighbour_code = _neighbour (octree->octants->codes[i], j);
        PDM_para_octree_direction_t inv_dir = _inv_direction(j);

        if (neighbour_code != NULL) {

          int i_rank = PDM_morton_quantile_search (n_ranks,
                                                   *neighbour_code,
                                                   block_octants_index);

          if (i_rank == rank) {

            int first_octant = PDM_morton_binary_search (octree->octants->n_nodes,
                                                         *neighbour_code,
                                                         octree->octants->codes);

            if (first_octant > i) {

              PDM_morton_compare_t pmc =
                PDM_morton_compare(dim,  *neighbour_code, octree->octants->codes[first_octant]);

              if ((pmc == PDM_MORTON_SAME_ANCHOR) || pmc ==  (PDM_MORTON_EQUAL_ID)) {

                if (neighbors_tmp[i].n_neighbor[j] >= neighbors_tmp[i].s_neighbor[j]) {
                  neighbors_tmp[i].s_neighbor[j] *= 2;
                  neighbors_tmp[i].neighbors[j] =
                    realloc (neighbors_tmp[i].neighbors[j],
                             sizeof(int) *  neighbors_tmp[i].s_neighbor[j]);
                }
                neighbors_tmp[i].neighbors[j][neighbors_tmp[i].n_neighbor[j]] = first_octant;
                neighbors_tmp[i].n_neighbor[j]++;

                if (neighbors_tmp[first_octant].n_neighbor[inv_dir] >=
                    neighbors_tmp[first_octant].s_neighbor[inv_dir]) {

                  neighbors_tmp[first_octant].s_neighbor[inv_dir] *= 2;
                  neighbors_tmp[first_octant].neighbors[inv_dir] =
                    realloc (neighbors_tmp[first_octant].neighbors[inv_dir],
                             sizeof(int) *  neighbors_tmp[first_octant].s_neighbor[inv_dir]);
                }
                neighbors_tmp[first_octant].neighbors[inv_dir][neighbors_tmp[first_octant].n_neighbor[inv_dir]] = i;
                neighbors_tmp[first_octant].n_neighbor[inv_dir]++;

              }

              while (first_octant < n_ranks - 1
                     && PDM_morton_a_ge_b (*neighbour_code, octree->octants->codes[first_octant+1])) {

                first_octant++;

                if (first_octant > i) {
                  PDM_morton_compare_t pmc2 =
                    PDM_morton_compare(dim,  *neighbour_code, octree->octants->codes[first_octant]);

                  if ((pmc2 == PDM_MORTON_SAME_ANCHOR) || (pmc2  == PDM_MORTON_EQUAL_ID)) {
                    if (neighbors_tmp[i].n_neighbor[j] >= neighbors_tmp[i].s_neighbor[j]) {
                      neighbors_tmp[i].s_neighbor[j] *= 2;
                      neighbors_tmp[i].neighbors[j] =
                        realloc (neighbors_tmp[i].neighbors[j],
                                 sizeof(int) *  neighbors_tmp[i].s_neighbor[j]);
                    }
                    neighbors_tmp[i].neighbors[j][neighbors_tmp[i].n_neighbor[j]] = first_octant;
                    neighbors_tmp[i].n_neighbor[j]++;

                    if (neighbors_tmp[first_octant].n_neighbor[inv_dir] >=
                        neighbors_tmp[first_octant].s_neighbor[inv_dir]) {

                      neighbors_tmp[first_octant].s_neighbor[inv_dir] *= 2;
                      neighbors_tmp[first_octant].neighbors[inv_dir] =
                        realloc (neighbors_tmp[first_octant].neighbors[inv_dir],
                                 sizeof(int) *  neighbors_tmp[first_octant].s_neighbor[inv_dir]);
                    }
                    neighbors_tmp[first_octant].neighbors[inv_dir][neighbors_tmp[first_octant].n_neighbor[inv_dir]] = i;
                    neighbors_tmp[first_octant].n_neighbor[inv_dir]++;

                  }
                }
              }
            }
          }

          else {

            if (neighbors_tmp[i].n_neighbor[j] >= neighbors_tmp[i].s_neighbor[j]) {
              neighbors_tmp[i].s_neighbor[j] *= 2;
              neighbors_tmp[i].neighbors[j] =
                realloc (neighbors_tmp[i].neighbors[j],
                         sizeof(int) *  neighbors_tmp[i].s_neighbor[j]);
            }

            neighbors_tmp[i].neighbors[j][neighbors_tmp[i].n_neighbor[j]] = -(i_rank + 1);
            neighbors_tmp[i].n_neighbor[j]++;

          }
        }

        free (neighbour_code);

      }
    }

    while (n_init_node > 0) {
      int current_id = octree->octants->n_nodes - n_init_node;
      _octants_replace_node_by_child (octree->octants,
                                      octree->points_in_leaf_max,
                                      current_id,
                                      octree->n_points,
                                      octree->points_code,
                                      &neighbors_tmp);
      n_init_node += 1;
    }
  }

  else {

    _octants_replace_node_by_child (octree->octants,
                                    octree->points_in_leaf_max,
                                    0,
                                    octree->n_points,
                                    octree->points_code,
                                    &neighbors_tmp);

  }


  /*************************************************************************
   *
   * Copy temporary neighbours in the neighbor structure
   *
   *************************************************************************/

  octree->octants->neighbor_idx =
    malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  int idx = 0;
  octree->octants->neighbor_idx[0] = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbor_idx[idx+1] =
        octree->octants->neighbor_idx[idx] + neighbors_tmp[i].n_neighbor[j];
      idx += 1;
    }
  }

  octree->octants->neighbors =
    malloc(sizeof(int) *
           octree->octants->neighbor_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      for (int k = 0; k < neighbors_tmp[i].n_neighbor[j]; k++) {
        octree->octants->neighbors[idx++] = neighbors_tmp[i].neighbors[j][k];
      }
    }
  }

  /* Free temporary arrays */

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      if (neighbors_tmp[i].neighbors[j] != NULL) {
        free (neighbors_tmp[i].neighbors[j]);
      }
    }
  }

  free (neighbors_tmp);


  /*************************************************************************
   *
   * Build parallel partition boundary
   *
   *************************************************************************/

  if (n_ranks > 1) {

    const int n_quantile =  n_ranks * n_direction;

    int *neighbor_rank_n = malloc (sizeof(int) * n_quantile);
    int *neighbor_rank_idx = malloc (sizeof(int) * (n_quantile + 1));
    int n_neighbor_rank = 0;

    for (int i = 0; i < n_quantile; i++) {
      neighbor_rank_n[i] = 0;
    }

    /* Premiere boucle pour compter */

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = octree->octants->neighbor_idx[n_direction * i + j];
             k < octree->octants->neighbor_idx[n_direction * i + j + 1];
             k++) {
          if (octree->octants->neighbors[k] < 0) {
            neighbor_rank_n[(PDM_ABS(octree->octants->neighbors[k]) - 1)*n_direction +j]++;
          }
        }
      }
    }

    neighbor_rank_idx[0] = 0;
    for (int i = 0; i < n_quantile; i++) {
      neighbor_rank_idx[i+1] = neighbor_rank_idx[i] + neighbor_rank_n[i];
      neighbor_rank_n[i] = 0;
    }

    /* Allocation */

    int *neighbor_rank_node_id = malloc (sizeof(int) * neighbor_rank_idx[n_quantile]);
    PDM_morton_code_t *neighbor_rank_code = malloc (sizeof(PDM_morton_code_t) *
                                                neighbor_rank_idx[n_quantile]);

    /* Deuxieme boucle pour stocker avec tri suivant la direction */

    int max_node_dir = -1;
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = octree->octants->neighbor_idx[n_direction * i + j];
             k < octree->octants->neighbor_idx[n_direction * i + j + 1];
             k++) {
          max_node_dir = PDM_MAX (max_node_dir,
                                  octree->octants->neighbor_idx[n_direction * i + j + 1] -
                                  octree->octants->neighbor_idx[n_direction * i + j]);
          if (octree->octants->neighbors[k] < 0) {
            int index = (PDM_ABS(octree->octants->neighbors[k]) - 1) * n_direction + j;
            int index2 = neighbor_rank_idx[index] + neighbor_rank_n[index];

            neighbor_rank_node_id[index2] = i;
            PDM_morton_copy (octree->octants->codes[i],
                             neighbor_rank_code + index2);

            neighbor_rank_n[index]++;
          }
        }
      }
    }

    /* Tri des codes */

    order = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_id = malloc (sizeof(int) * max_node_dir);
    PDM_morton_code_t *tmp_code = malloc (sizeof(int) * max_node_dir);

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        PDM_morton_local_order (octree->octants->neighbor_idx[n_direction * i + j +1] -
                                octree->octants->neighbor_idx[n_direction * i + j],
                                neighbor_rank_code + octree->octants->neighbor_idx[n_direction * i + j],
                                order);
        int idx1 = 0;
        for (int k = octree->octants->neighbor_idx[n_direction * i + j];
             k < octree->octants->neighbor_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (neighbor_rank_code[k], tmp_code + idx1);
          tmp_node_id[idx1] = neighbor_rank_node_id[k];
          idx1 += 1;
        }

        idx1 = 0;
        for (int k = octree->octants->neighbor_idx[n_direction * i + j];
             k < octree->octants->neighbor_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (tmp_code[order[idx1]], neighbor_rank_code + k );
          neighbor_rank_node_id[k] = tmp_node_id[order[idx1]];
          idx1 += 1;
        }
      }
    }

    free (order);
    free (neighbor_rank_node_id);
    free (neighbor_rank_code);

    neighbor_rank_node_id = tmp_node_id;
    neighbor_rank_code = tmp_code;

    /* Envoi reception (Les donnees recues sont triees) */

    for (int i = 0; i < n_ranks; i++) {


    }

    // Au retour des vacances : Envoyer neighbor_rank_node_id[k] et neighbor_rank_code */


    /* Troisieme boucle pour stocker les infos distantes */

    // neighbor_rank[i_rank] += 1;

    /* - Pour chaque rang en regard : envoi des noeuds (direction id et des codes de morton) en frontiere de ce rang avec au prealable Tri en fonction de la direction puis chaque direction tri des codes de morton  */

    /* - A reception : Pour chaque noeud frontiere en regard du rang recu : recherche des noeuds distants frontiere puis stockage */


    /* Remplacement du num proc dans neighbor par num_part_bound */

    free (neighbor_rank_n);
    free (neighbor_rank_idx);
    free (neighbor_rank_node_id);
    free (neighbor_rank_code);

  }

}


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
 const int  id
)
{
 _octree_t *octree = _get_from_id (id);

 return octree->global_extents;
}


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
)
{
 _octree_t *octree = _get_from_id (id);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
