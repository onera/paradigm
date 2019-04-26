
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
  octants->range = malloc (sizeof(int) * (octants->n_nodes_max+1));
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
                              sizeof(int) * (octants->n_nodes_max+1));
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

static void
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

static _l_octant_t *
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

static _l_octant_t *
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

  /* 
   * Redistribute octant list from coarse load balancing 
   */

  send_count = malloc(sizeof(int) * n_ranks);
  send_shift = malloc(sizeof(int) * (n_ranks+1));

  recv_count = malloc(sizeof(int) * n_ranks);
  recv_shift = malloc(sizeof(int) * (n_ranks+1));
  
  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  irank = 0;
  for (int i = 0; i < G->n_nodes; i++) {
    if (PDM_morton_a_ge_b (G->codes[i], _G_morton_index[irank+1])) {

      irank += 1 + PDM_morton_binary_search (n_ranks - (irank + 1),
                                             G->codes[i],
                                             _G_morton_index + irank + 1);
    }
    //send_count[irank] += send_shift[irank] + (octant_list->dim + 1);
    send_count[irank] += octant_list->dim + 1;
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

  send_codes =
    malloc (send_shift[n_ranks] * sizeof(PDM_morton_int_t));

  for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
    send_count[rank_id] = 0;
  }

  irank = 0;
  for (int i = 0; i < G->n_nodes; i++) {

    if (PDM_morton_a_ge_b (G->codes[i], rank_codes[irank+1])) {
        
      irank += 1 + PDM_morton_binary_search(n_ranks - (irank + 1),
                                            G->codes[i],
                                            _G_morton_index + irank + 1);
    }

    int shift = send_shift[irank] + send_count[irank];
    send_codes[shift++] = G->codes[i].L;
    
    for (int j = 0; j < octant_list->dim; j++) {
      send_codes[shift++] = G->codes[i].X[j];
    }
    
    send_count[irank] += octant_list->dim + 1;
  }

  recv_codes = malloc (recv_shift[n_ranks] * sizeof(PDM_morton_int_t));

  /* - exchange codes between processes */
  
  PDM_MPI_Alltoallv(send_codes, send_count, send_shift, PDM_MPI_UNSIGNED,
                    recv_codes, recv_count, recv_shift, PDM_MPI_UNSIGNED,
                    comm);

  free (send_codes);
  free (send_count);
  free (send_shift);
  free (recv_count);

  /* - tri des codes recus */

  _octants_free (G);
  _octants_init (G, octant_list->dim, recv_shift[n_ranks]/(octant_list->dim + 1));

  int idx = 0;
  for (int i = 0; i < G->n_nodes; i++) {
    G->codes[i].L = recv_codes[idx++];
    for (int j = 0; j < octant_list->dim; j++) {
      G->codes[i].X[j] = recv_codes[idx++];
    }
  }
  
  PDM_morton_local_sort (G->n_nodes, G->codes);

  /* - pas de tri des donnees */
  
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
      if (iblock < octree->octants->n_nodes) {
        if (PDM_morton_a_ge_b (octree->points_code[i], octree->octants->codes[iblock+1])) {
        
          iblock += 1 + PDM_morton_binary_search (octree->octants->n_nodes - (iblock + 1),
                                                  octree->points_code[i],
                                                  octree->octants->codes + iblock + 1);
        }
      }
      octree->octants->n_points[iblock] += 1;
    }

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      octree->octants->range[i] = 0;
    }
    
    for (int i = 0; i < octree->octants->n_nodes - 1; i++) {
      octree->octants->range[i+1] =
        octree->octants->range[i] + octree->octants->n_points[iblock];
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
    
    _octants_add (octree->octants,
                  code,
                  octree->n_points,
                  0,
                  0,
                  0,
                  NULL);
    
  }
  
  /*************************************************************************
   *
   * Add child (while n points > N points max)
   *     - Build neighbor 
   *
   *************************************************************************/

  int _max_size = octree->octants->n_nodes_max;

  int **_neighbor_idx = malloc(sizeof(int*) * _max_size);

  int **_neighbor = malloc(sizeof(int*) * _max_size);

  for (int i = 0; i < _max_size; i++) {
    _neighbor_idx[i] = NULL;
    _neighbor[i] = NULL;
  }

  //TODO 26/04: Continuer a partir d'ici + faire algo 4, 3 et 7 de sundar. Puis faire algo de
  //  de localisation de points puis d'intersection avec une boite englobante.
  
  if (n_ranks > 1) {

    /*  on part des octants blocks et on subdivise */

    /* Recherche des voisins (processus ou octant local) */

    /* Raffinement automatique :
         - remplacement du noeud par ses enfants
         - mise a jour des voisins (remplancement du noeud par les 4 de la direction concernee */
    
    /* Echange de la table de voisinage (octants fantome à ajouter ?)
         - remplacement du noeud par ses enfants
         - mise a jour des voisins (remplancement du noeud par les 4 de la direction concernee */
    
  }
  
  else {
    /* on part de la racine et on subdivise */


    

  }

  /* Copy temporary neighbours in the octree structure*/

  const int n_direction = 1 << octree->dim;
  
  octree->octants->neighbor_idx =
    malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  octree->octants->neighbor_idx[0] = 0;
  int idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbor_idx[idx+1] =
        octree->octants->neighbor_idx[idx] + (_neighbor_idx[i][j+1] - _neighbor_idx[i][j]);
    }
  }

  octree->octants->neighbors =
    malloc(sizeof(int) * octree->octants->neighbor_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < _neighbor_idx[i][n_direction]; j++) {
      octree->octants->neighbors[idx++] = _neighbor[i][j];
    }
  }

  /* Free temporary arrays */

  for (int i = 0; i < _max_size; i++) {
    if (_neighbor_idx[i] != NULL) {
      free (_neighbor_idx[i]);
    }
    if (_neighbor[i] != NULL) {
      free (_neighbor[i]);
    }
  }

  free (_neighbor_idx);
  free (_neighbor);
  
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
