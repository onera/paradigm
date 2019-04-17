
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
  int   *ancestor;                /*!< rank + id_node size = n_nodes */
  int   *child;                /*!< rank + id_node size = 8 * n_nodes */
  int   dim;
  
} _l_octant_t;

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

  _l_octant_t *octants;       /*!< list of octants */

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

  if (octants->ancestor != NULL) {
    free (octants->ancestor);
  }

  if (octants->child != NULL) {
    free (octants->child);
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
  octants->range = malloc (sizeof(int) * octants->n_nodes_max);
  octants->is_leaf = malloc (sizeof(int) * octants->n_nodes_max);
  octants->ancestor = malloc (sizeof(int) * octants->n_nodes_max);
  octants->child    = malloc (sizeof(int) * octant_dim * octants->n_nodes_max);
  
  octants->neighbor_idx = NULL;
  octants->neighbors    = NULL;
  octants->dim = octant_dim;
}


/**
 *
 * \brief Check size of the size of a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_check_alloc
(
 _l_octant_t *octants
)
{
  if (octants->n_nodes >= octants->n_nodes_max) {
    
    octants->n_nodes_max *= 2;
    
    octants->codes    = realloc (octants->codes,
                                 sizeof(PDM_morton_code_t) * octants->n_nodes_max);
    octants->n_points = realloc (octants->n_points,
                                 sizeof(int) * octants->n_nodes_max);
    octants->range = realloc (octants->range,
                              sizeof(int) * octants->n_nodes_max);
    octants->is_leaf = realloc (octants->is_leaf,
                                sizeof(int) * octants->n_nodes_max);
    octants->ancestor = realloc (octants->ancestor,
                                 sizeof(int) * octants->n_nodes_max);
    octants->child    = realloc (octants->child,
                                 sizeof(int) * octants->dim * octants->n_nodes_max);

    octants->neighbor_idx = NULL;
    octants->neighbors    = NULL;
    
  }
}



/**
 *
 * \brief Check size of the size of a list of octants
 *
 * \param [in]   octants     Octants
 *
 */

static void
_octants_add
(
 _l_octant_t *octants,
 const PDM_morton_code_t code,
 const int n_points,
 const int range,
 const int is_leaf,
 const int ancestor,
 const int *child_id
)
{

  _octants_check_alloc (octants);
    
  const int idx = octants->n_nodes;
    
  PDM_morton_copy (code, octants->codes + idx);
  
  octants->n_points[idx] = n_points;
  
  octants->range[idx] = range;
                             
  octants->is_leaf[idx] = is_leaf;
    
  octants->ancestor[idx] = ancestor;

  if (child_id != NULL) {
    int idx2 = idx * octants->dim;
    for (int i = 0; i < octants->dim; i++) {
      octants->child[idx2+i] = child_id[i];
    }
  }


  octants->neighbor_idx = NULL;
  octants->neighbors    = NULL;
    
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
      ((double) code.X[i]/(double) code.L)* octree->d[i] + octree->s[i]; 
    extents[octree->dim + i] =
      (((double) code.X[i] + 1)/(double) code.L)* octree->d[i] + octree->s[i]; 
  }
}


/**
 *
 * \brief Removing overlaps from a sorted lis of octants
 *
 * \param [inout]  octants A lis of octants
 *
 */

void
_linearize
(
 _l_octant_t octants
)
{
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

_l_octant_t *
_complete_region
(
 PDM_morton_code_t a,
 PDM_morton_code_t b
)
{
  _l_octant_t *_octants = NULL;


  //TODO:
  
  return _octants;
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

_l_octant_t *
_complete_octree
(
 _l_octant_t *L,
 PDM_MPI_Comm comm
)
{
  _l_octant_t *_octants = NULL;

   /* _remove_duplicates (L, comm); */

   /*     PDM_morton_build_rank_index(dim, */
   /*                              max_level, */
   /*                              octree->n_points, */
   /*                              octree->points_code, */
   /*                              weight, */
   /*                              order, */
   /*                              morton_index, */
   /*                              octree->comm); */

  
  return _octants;
}


/**
 *
 * \brief Partitioning octants into large contiguous blocks. The list of octants
 *        is redistributed
 *
 * \param [inout]  octant_list  a list of octants,  
 *
 * \return block_octants A list of blocks
 *
 */

_l_octant_t *
_block_partition
(
 _l_octant_t *octant_list,
 const PDM_MPI_Comm comm
)
{
 
  /* Complete region */

  _l_octant_t *T = _complete_region (octant_list->codes[0],
                                     octant_list->codes[octant_list->n_nodes]);

  int max_level = -1;
  for (int i = 0; i < octant_list->n_nodes; i++) {
    max_level = PDM_MAX (octant_list->codes[i].L, max_level);
  }

  /* Complete octree */

  _l_octant_t C;
  
  _octants_init (&C, octant_list->dim, octant_list->n_nodes);
                 
  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (octant_list->codes[i].L >= max_level) {
      _octants_add (&C,
                    octant_list->codes[i],
                    octant_list->n_points[i],
                    octant_list->range[i],
                    octant_list->is_leaf[i],
                    octant_list->ancestor[i],
                    NULL);
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

  for (int i = 0; i < octant_list->n_nodes; i++) {
    int irank = PDM_morton_binary_search(n_ranks,
                                         octant_list->codes[i],
                                         rank_codes);
    //rank_counts[irank] += 1;
  }

  /* - compute weight of each cell (*/





  
  /* 
   * Load balancing G from weight 
   */

  /* PDM_morton_build_rank_index(dim, */
  /*                             max_level, */
  /*                             octree->n_points, */
  /*                             octree->points_code, */
  /*                             weight, */
  /*                             order, */
  /*                             morton_index, */
  /*                             octree->comm); */
     
  /* Redistribute octant_list : MPI_allgather about codes */

  /* 
   * Redistirbute octant list from coarse load balancing 
   */


  
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
  octree->points = realloc (octree->points, octree->n_points * sizeof(double) * octree->dim);
  octree->points_icloud = realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum = realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code = realloc (octree->points_code, octree->n_points * sizeof(PDM_morton_code_t));

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
  }
 
  free (order);
  
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

      _octants_check_alloc (point_octants);

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

  _l_octant_t *block_octants = _block_partition (point_octants, octree->comm);

  /*************************************************************************
   *
   * Add child (whilen points > N points max)
   *     - sundar does not keep ancestors (Keep ancestor for our application ?)
   *     - Build neighbor 
   *
   *************************************************************************/
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
