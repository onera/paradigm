
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
 * \struct _octant_t
 * \brief  Define an octant
 * 
 */

typedef struct  {

  PDM_morton_code_t code; /*!< morton code */

  int  n_points;          /*!< Number of points in octant*/
  int  range;             /*!< Start index of point list for each octant */
  int  is_leaf;           /*!< IS a leaf >*/
  
} _octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 * 
 */

typedef struct  {

  double  global_extents[6];            /*!< Extents of current process */ 
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

  int   n_nodes;                 /*!< Current number of nodes in octree */
  int   n_nodes_max;             /*!< Maximum number of nodes in octree */
  _octant_t   *nodes;            /*!< Array of octree nodes
                                       (size: n_nodes_max) */
  int   *neighbor_idx;
  int   *neighbors;               /*!< rank + id_node size = 2 * n_nodes */
  int   *ancestor;                /*!< rank + id_node size = n_nodes */
  int   *child;                /*!< rank + id_node size = 8 * n_nodes */
  double   *extents;                /*!< rank + id_node size = 6 * n_nodes */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */
  
} _octree_t;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees   = NULL;

static const double _eps_default = 1.e-12;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Check size of the octree and realloc it if necessary
 *
 * \param [in]   octree   octree to realloc
 *
 */

static void
_check_alloc
(
 _octree_t *octree
)
{
  if (octree->n_nodes >= octree->n_nodes_max) {
    
    if (octree->n_nodes_max == 0) {
      octree->n_nodes_max = octree->n_points;
    }
    else {
      octree->n_nodes_max *= 2;
    }
    
    octree->nodes    = realloc (octree->nodes,
                                sizeof(_octant_t) * octree->n_nodes_max);
    octree->ancestor = realloc (octree->ancestor,
                                sizeof(_octant_t) * octree->n_nodes_max);
    octree->child    = realloc (octree->child,
                                sizeof(_octant_t) * octree->dim * octree->n_nodes_max);
    octree->extents  = realloc (octree->extents,
                                sizeof(_octant_t) * octree->dim * 2 * octree->n_nodes_max);
    octree->neighbor_idx = NULL;
    octree->neighbors      = NULL;
    
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

  octree->n_nodes = 0;
  octree->n_nodes_max = -1;
  octree->nodes = NULL;
  
  octree->neighbor_idx = NULL;
  octree->neighbors = NULL;
  octree->ancestor = NULL;
  octree->child = NULL;
  octree->extents = NULL;
  
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

  if (octree->nodes != NULL) {
    free (octree->nodes);
  }
  
  if (octree->neighbor_idx != NULL) {
    free (octree->neighbor_idx);
  }
        
  if (octree->neighbors != NULL) {
    free (octree->neighbors);
  }
        
  if (octree->ancestor != NULL) {
    free (octree->ancestor);
  }
        
  if (octree->child != NULL) {
    free (octree->child);
  }
        
  if (octree->extents != NULL) {
    free (octree->extents);
  }

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
  octree->points = realloc (octree->points, octree->n_points * sizeof(double) * octree->dim);
  octree->points_icloud = realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum = realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code = realloc (octree->points_code, octree->n_points * sizeof(_octant_t));

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
    
    PDM_morton_code_t *morton_index = malloc (sizeof(PDM_morton_code_t) * (n_ranks + 1));
    
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
    
    int *c_rank = malloc (octree->n_points * sizeof(int));
    
    for (int i = 0; i < octree->n_points; i++) {
      size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                  octree->points_code[i],
                                                  morton_index);
      c_rank[i] = (int) _c_rank; 
    }
    
    free(morton_index);

    int *send_count = malloc (n_ranks * sizeof (int));
    int *recv_count = malloc (n_ranks * sizeof (int));
    int *send_shift = malloc ((n_ranks + 1) * sizeof (int));
    int *recv_shift = malloc ((n_ranks + 1) * sizeof (int));

    for (int rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (int i = 0; i < octree->n_points; i++)
      send_count[c_rank[i]] += octree->dim;

    /* Exchange number of coords to send to each process */

    PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, octree->comm);

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

    for (int i = 0; i < octree->n_points; i++) {
      int rank_id = c_rank[i];
      int shift = send_shift[rank_id] + send_count[rank_id];
      for (int j = 0; j < dim; j++)
        send_coords[shift + j] = octree->points[i*dim + j];
      send_count[rank_id] += dim;
    }

    double *recv_coords = malloc (recv_shift[n_ranks] * sizeof(double));

    /* Exchange coords between processes */

    PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                      recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                      octree->comm);

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

    for (int i = 0; i < octree->n_points; i++) {
      int rank_id = c_rank[i];
      int shift = send_shift[rank_id] + send_count[rank_id];
      send_points_icloud[shift] = octree->points_icloud[i];
      send_count[rank_id] += 1;
    }

    int *recv_points_icloud = malloc (recv_shift[n_ranks] * sizeof(int));

    /* Exchange points_icloud between processes */

    PDM_MPI_Alltoallv(send_points_icloud, send_count, send_shift, PDM_MPI_INT,
                      recv_points_icloud, recv_count, recv_shift, PDM_MPI_INT,
                      octree->comm);

    free(send_points_icloud);


    /* Build send and receive buffers : points_gnum*/

    PDM_g_num_t *send_points_gnum = malloc (send_shift[n_ranks] * sizeof(PDM_g_num_t));

    for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] = 0;
    }

    for (int i = 0; i < octree->n_points; i++) {
      int rank_id = c_rank[i];
      int shift = send_shift[rank_id] + send_count[rank_id];
      send_points_gnum[shift] = octree->points_gnum[i];
      send_count[rank_id] += 1;
    }

    PDM_g_num_t *recv_points_gnum = malloc (recv_shift[n_ranks] * sizeof(PDM_g_num_t));

    /* Exchange points_gnum between processes */

    PDM_MPI_Alltoallv(send_points_gnum, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                      recv_points_gnum, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                      octree->comm);

    free(send_points_gnum);

    octree->n_points = recv_shift[n_ranks];

    free (send_count);
    free (recv_count);
    free (send_shift);
    free (recv_shift);

    octree->points = realloc (octree->points, sizeof(double) * octree->n_points);

    octree->points_icloud = realloc (octree->points_icloud, sizeof(int) * octree->n_points);

    octree->points_gnum = realloc (octree->points_gnum, sizeof(PDM_g_num_t) * octree->n_points);

    /* Re-encode points */

    octree->points_code = realloc (octree->points_code,
                                   sizeof(PDM_morton_code_t) * octree->n_points);

    PDM_morton_encode_coords(dim,
                             max_level,
                             octree->global_extents,
                             octree->n_points,
                             octree->points,
                             octree->points_code,
                             octree->d,
                             octree->s);

    order = realloc (order, sizeof(int) * octree->n_points);

    for (int i = 0; i < octree->n_points; i++) {
      order[i] = i;
    }

    PDM_morton_local_order(octree->n_points, octree->points_code, order);

    for (int i = 0; i < octree->n_points; i++) {
      octree->points_icloud[i] = recv_points_icloud[order[i]];
      octree->points_gnum[i] = recv_points_gnum[order[i]];
      for (int j = 0; j < dim; j++) {
        octree->points[dim*i+j] = recv_coords[dim*order[i]+j];
      }
    }

    free (recv_points_icloud);
    free (recv_points_gnum);
    free (recv_coords);
    
    
    PDM_morton_code_t *_points_code = malloc (sizeof(PDM_morton_code_t) * octree->n_points);
    
    for (int i = 0; i < octree->n_points; i++) {
      _points_code[i].L = octree->points_code[order[i]].L;
      _points_code[i].X[0] = octree->points_code[order[i]].X[0];
      _points_code[i].X[1] = octree->points_code[order[i]].X[1];
      _points_code[i].X[2] = octree->points_code[order[i]].X[2];
    }

    free (octree->points_code);
    octree->points_code = _points_code;
    
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
    
    PDM_morton_code_t *_points_code = malloc (sizeof(PDM_morton_code_t) * octree->n_points);
    
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
  }
 
  free (order);
  
  /*************************************************************************
   *
   * Assign point morton code level to octree->depth_max
   *
   *************************************************************************/

  for (int i = 0; i < octree->n_points; i++) {
    PDM_morton_assign_level (octree->points_code[i], octree->depth_max);
  }

  /*************************************************************************
   *
   * Store points in the octants (leaves) at the maximum depth of the octree 
   * to build
   *
   *************************************************************************/

  int chg_code = 1;
  _octant_t *curr_node = NULL;

  _check_alloc (octree);

  for (int i = 0; i < octree->n_points; i++) {

    if (curr_node != NULL) {
      chg_code = !(PDM_morton_a_eq_b(curr_node->code,
                                     octree->points_code[i]));
    }
    
    if (chg_code) {

      _check_alloc (octree);

      int idx = octree->n_nodes;

      _octant_t *_node =  octree->nodes + idx;
      curr_node = _node;
      
      PDM_morton_copy (octree->points_code[i], &(_node->code));
      
      _node->is_leaf = 1;
      _node->n_points = 1;
      _node->range = i;
      
      octree->n_nodes += 1;
    }

    else {
      curr_node->n_points += 1;
    }

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

 return NULL;
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
