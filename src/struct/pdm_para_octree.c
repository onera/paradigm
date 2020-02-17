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
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_hash_tab.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

#define NTIMER 11

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  BUILD_ORDER_POINTS            = 1,
  BUILD_BLOCK_PARTITION         = 2,
  BUILD_LOCAL_NODES             = 3,
  BUILD_LOCAL_NEIGHBOURS_STEP1  = 4,
  BUILD_LOCAL_NEIGHBOURS_STEP2  = 5,
  BUILD_LOCAL_NEIGHBOURS_STEP3  = 6,
  BUILD_LOCAL_NEIGHBOURS        = 7,
  BUILD_DISTANT_NEIGHBOURS      = 8,
  BUILD_TOTAL                   = 9,
  END                           = 10,

} _ol_timer_step_t;


/**
 * \struct _heap_t
 * \brief  Heap used to recursively subdivide nodes
 *
 */

typedef struct  {

  int   top;                  /*!< Top of head  */
  int   size;                 /*!< Size of heap */
  PDM_morton_code_t *codes;   /*!< Morton codes */
  int *range;                 /*!< Points range */
  int *n_points;              /*!< Points number */
  int   max_top;

} _heap_t;


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

typedef struct  {

  double  global_extents[6];     /*!< Extents of current process */
  int     depth_max;             /*!< Maximum depth of the three */
  int     points_in_leaf_max;    /*!< Maximum number of points in a leaf */
  double      s[3];           /*!< Translation for the normalization */
  double      d[3];           /*!< Dilatation for the normalization */

  int     n_point_clouds;        /*!< Number of point cloud */

  PDM_g_num_t    t_n_points;         /*!< total number of points */
  int            n_points;           /*!< Number of points in each cloud */
  double *points;                    /*!< Point coordinates */
  int *points_icloud;                /*!< Point cloud */
  PDM_g_num_t *points_gnum;          /*!< Point global number */
  PDM_morton_code_t  *points_code;   /*!< Morton codes */

  PDM_morton_code_t *rank_octants_index;
  _l_octant_t *octants;       /*!< list of octants */

  PDM_MPI_Comm comm;           /*!< MPI communicator */
  int   dim;                     /*!< Dimension */

  int n_part_boundary_elt;    /*!< Number of partitioning boundary element */
  int *part_boundary_elt_idx; /*!< Index for part_boundary_elt (size=\ref n_part_boundary_elt + 1 */
  int *part_boundary_elt;     /*!< Partitioning boundary elements description (proc number + element number) */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

  int neighboursToBuild;

  int  n_connected;
  int *connected_idx;

} _octree_t;



/**
 * \struct _neighbours_tmp_t
 * \brief  Define a temporary neighbour structure
 *
 */


typedef struct  {

  int n_neighbour[6];     /*!< Number of neighbours in the arrays  */
  int s_neighbour[6];     /*!< Size of arrays */
  int *neighbours[6];     /*!< Arrays */

} _neighbours_tmp_t;




//-->>
typedef struct {

  int          size;
  int          count;
  int         *lnum;
  PDM_g_num_t *gnum;
  double      *priority;

} _min_heap_t;
//<<--


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_octrees    = NULL;

//static const double _eps_default  = 1.e-12;

static const PDM_morton_int_t max_morton_level = 15;
//static const int max_morton_level = 2;

/*============================================================================
 * Private function definitions
 *============================================================================*/

static _min_heap_t *
_min_heap_create
(
 const int size
 )
{
  _min_heap_t *h = malloc (sizeof(_min_heap_t));
  h->size        = size;
  h->count       = 0;
  h->lnum        = malloc (sizeof(int)         * size);
  h->gnum        = malloc (sizeof(PDM_g_num_t) * size);
  h->priority    = malloc (sizeof(double)      * size);

  for (int i = 0; i < h->size; i++)
    h->priority[i] = HUGE_VAL;

  return h;
}


static void
_min_heap_reset
(
 _min_heap_t *h
 )
{
  h->count = 0;
}


static void
_min_heap_free
(
 _min_heap_t *h
 )
{
  free (h->lnum);
  free (h->gnum);
  free (h->priority);
  free (h);
}



static void
_min_heap_swap
(
 _min_heap_t *h,
 const int    a,
 const int    b
 )
{
  int l          = h->lnum[a];
  PDM_g_num_t g  = h->gnum[a];
  double p       = h->priority[a];

  h->lnum[a]     = h->lnum[b];
  h->gnum[a]     = h->gnum[b];
  h->priority[a] = h->priority[b];

  h->lnum[b]     = l;
  h->gnum[b]     = g;
  h->priority[b] = p;
}



static void
_min_heap_push
(
 _min_heap_t       *h,
 const int          lnum,
 const PDM_g_num_t  gnum,
 const double       priority
 )
{
  int i = 0;
  /* make sure the heap is large enough to contain the new element */
  if (h->count >= h->size) {
    h->size = (h->size) ? 2*h->size : 10;
    h->lnum =     (int *)         realloc (h->lnum,     sizeof(int)         * h->size);
    h->gnum =     (PDM_g_num_t *) realloc (h->gnum,     sizeof(PDM_g_num_t) * h->size);
    h->priority = (double *)      realloc (h->priority, sizeof(double)      * h->size);
    for (i = h->count+1; i < h->size; i++)
      h->priority[i] = HUGE_VAL;
  }

  i = h->count;
  h->count++;
  h->lnum[i]     = lnum;
  h->gnum[i]     = gnum;
  h->priority[i] = priority;

  int parent = (i - 1)/2;
  while (i > 0 && h->priority[parent] > priority) {
    _min_heap_swap (h, parent, i);
    i = parent;
    parent = (i - 1)/2;
  }
}



static void
_min_heap_heapify
(
 _min_heap_t *h,
 const int    i
 )
{
  int l = 2*i+1; // left child node
  int r = 2*i+2; // right child node
  int s = i; // node with smallest priority

  if (l < h->count && h->priority[l] < h->priority[i])
    s = l;

  if (r < h->count && h->priority[r] < h->priority[s])
    s = r;

  if (s != i) {
    _min_heap_swap (h, i, s);
    _min_heap_heapify (h, s);
  }
}


static int
_min_heap_pop
(
 _min_heap_t *h,
 int         *lnum,
 PDM_g_num_t *gnum,
 double      *priority
 )
{
  if (h->count < 1)
    return 0;

  *lnum     = h->lnum[0];
  *gnum     = h->gnum[0];
  *priority = h->priority[0];

  h->count--;
  h->lnum[0]     = h->lnum[h->count];
  h->gnum[0]     = h->gnum[h->count];
  h->priority[0] = h->priority[h->count];

  _min_heap_heapify (h, 0);

  return 1;
}




inline static double
_octant_min_dist2
(
 const int          dim,
 PDM_morton_code_t  code,
 const double      *d,
 const double      *s,
 const double      *coords
 )
{
  double min_dist2 = 0., delta = 0.;
  double side = 1./pow(2, code.L);

  for (int i = 0; i < dim; i++) {
    double x = coords[i];
    double xmin = s[i] + d[i] * side * code.X[i];
    double xmax = xmin + d[i] * side;

    if (x > xmax) {
      delta = x - xmax;
      min_dist2 += delta * delta;
    } else if (x < xmin) {
      delta = x - xmin;
      min_dist2 += delta * delta;
    }
  }

  return min_dist2;
}





/**
 *
 * \brief Create a heap
 *
 * \param [in]  Size    Size of heap
 *
 * \return   a new heap
 *
 */

static _heap_t *
_heap_create
(
const int size
)
{
  _heap_t *heap = malloc(sizeof(_heap_t));
  heap->top = 0;
  heap->max_top = 0;
  heap->size = size;
  heap->codes = malloc(sizeof(PDM_morton_code_t) * size);
  heap->range =  malloc(sizeof(int) * size);
  heap->n_points =  malloc(sizeof(int) * size);
  return heap;
}


/**
 *
 * \brief Free a heap
 *
 * \param [in]  heap   Heap to free
 *
 * \return NULL
 *
 */

static _heap_t *
_heap_free
(
_heap_t *heap
)
{
  free (heap->codes);
  free (heap->range);
  free (heap->n_points);
  free (heap);
  return NULL;
}


/**
 *
 * \brief Push a new element in the heap
 *
 * \param [inout]  heap      Heap
 * \param [in]     code      Morton code
 * \param [in]     range     Range
 * \param [in]     n_points  Number of points
 *
 * \return  1 if pushed 0 otherwise
 *
 */

static int
_heap_push
(
 _heap_t *heap,
 const PDM_morton_code_t code,
 const int range,
 const int n_points
)
{
  if (heap->top >= heap->size) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (code, &(heap->codes[idx]));
  heap->range[idx] = range;
  heap->n_points[idx] = n_points;
  heap->top++;
  heap->max_top = PDM_MAX(heap->top, heap->max_top);
  return 1;
}


/**
 *
 * \brief Pull top element of the heap
 *
 * \param [inout]  heap      Heap
 * \param [out]    code      Morton code
 * \param [out]    range     Range
 * \param [out]    n_points  Number of points
 *
 * \return  1 if pulled 0 otherwise
 *
 */

static int
_heap_pull
(
 _heap_t *heap,
 PDM_morton_code_t *code,
 int *range,
 int *n_points
)
{
  heap->top--;
  if (heap->top < 0) {
    return 0;
  }
  int idx = heap->top;
  PDM_morton_copy (heap->codes[idx], code);
  *range = heap->range[idx];
  *n_points = heap->n_points[idx];
  return 1;
}

/**
 *
 * \brief Neighbour
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
 * \brief Neighbour
 *
 * \param [inout]   octants     Octants
 *
 * \return neighbour or NULL ()
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

  if (((_direction > 0) && (code.X[dim] < (pow(2,code.L) - 1))) ||
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

static void
_octants_purge
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

  if (octants->range != NULL) {
    free (octants->range);
  }

  if (octants->neighbour_idx != NULL) {
    free (octants->neighbour_idx);
  }

  if (octants->neighbours != NULL) {
    free (octants->neighbours);
  }
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

  _octants_purge (octants);

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

  octants->neighbour_idx = NULL;
  octants->neighbours    = NULL;
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

static int
_octants_check_alloc
(
 _l_octant_t *octants,
 const int n_free_node
)
{
  int is_realloc = 0;
  if (octants->n_nodes + n_free_node > octants->n_nodes_max) {

    octants->n_nodes_max *= 2;

    octants->codes    = realloc (octants->codes,
                                 sizeof(PDM_morton_code_t) * octants->n_nodes_max);
    octants->n_points = realloc (octants->n_points,
                                 sizeof(int) * octants->n_nodes_max);
    octants->range = realloc (octants->range,
                              sizeof(int) * (octants->n_nodes_max+1));
    octants->neighbour_idx = NULL;
    octants->neighbours    = NULL;

    is_realloc = 1;
  }
  return is_realloc;
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
 const int range
)
{

  _octants_check_alloc (octants, 1);

  const int idx = octants->n_nodes;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

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
 const int range
)
{

  _octants_check_alloc (octants, 1);

  for (int i = octants->n_nodes; i > 0; i--) {

    PDM_morton_copy (octants->codes[i - 1], octants->codes + i);

    octants->n_points[i] =  octants->n_points[i-1];

    octants->range[i] = octants->range[i-1];
  }

  const int idx = 0;

  PDM_morton_copy (code, octants->codes + idx);

  octants->n_points[idx] = n_points;

  octants->range[idx] = range;

  octants->n_nodes += 1;

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

/* static void */
/* g_extents */
/* ( */
/*  _octree_t *octree, */
/*  PDM_morton_code_t code, */
/*  double    extents[] */
/* ) */
/* { */
/*   for (int i = 0; i < octree->dim; i++) { */
/*     extents[i] = */
/*       ((double) code.X[i]/((double) pow(2,code.L)))* octree->d[i] + octree->s[i]; */
/*     extents[octree->dim + i] = */
/*       (((double) code.X[i] + 1)/((double) pow(2,code.L))) * octree->d[i] + octree->s[i]; */
/*   } */
/* } */

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

  if  (octants->n_nodes > 0) {
    for (int i = 0; i < octants->n_nodes - 1; i++) {

      if (!PDM_morton_ancestor_is (_codes[i], _codes[i+1])) {
        _octants_push_back (r_octants,
                            _codes[i],
                            0,
                            0);
      }

    }

    _octants_push_back (r_octants,
                        _codes[octants->n_nodes-1],
                        0,
                        0);
  }

  return r_octants;
}


/**
 *
 * \brief Constructing a minimal linear octree between two octants
 *
 * \param [in]  a     Morton code a
 * \param [in]  b     Morton code b
 *
 * \return octants The minimal linear octree between a and b or NULL if a >= b
 *
 */

static _l_octant_t *
_complete_region
(
 PDM_morton_code_t a,
 PDM_morton_code_t b
)
{
  const int dim = 3;

  _l_octant_t *r_octants = NULL;

  if (PDM_morton_a_gt_b (b, a)) {

    _l_octant_t *w_octants = malloc(sizeof(_l_octant_t));
    r_octants = malloc(sizeof(_l_octant_t));

    _octants_init (w_octants, dim, 4);
    _octants_init (r_octants, dim, 4);

    /* printf("_complete_region\n"); */

    /* printf("a_d\n"); */
    /* PDM_morton_dump(3, a); */
    /* printf("a_f\n"); */

    /* printf("b_d\n"); */
    /* PDM_morton_dump(3, b); */
    /* printf("b_f\n"); */

    PDM_morton_code_t nca;
    PDM_morton_nearest_common_ancestor (a, b, &nca);

    const int n_child = 8;

    int  size = PDM_morton_max_level * 8;
    _heap_t *heap = _heap_create (size);

    PDM_morton_code_t children[8];
    PDM_morton_get_children(dim,
                            nca,
                            children);

    for (int i = n_child - 1; i >= 0; i--) {
      int is_pushed = _heap_push (heap,
                                  children[i],
                                  0,
                                  0);
      if (!is_pushed) {
        printf ("Internal error PDM_para_octree 1 : heap is full\n");
        exit(1);
      }
    }

    PDM_morton_code_t code;
    int range;
    int n_points;

    while (_heap_pull (heap, &code, &range, &n_points)) {
      if (PDM_morton_a_gt_b (code, a) &&
          PDM_morton_a_gt_b (b, code) &&
          !PDM_morton_ancestor_is (code, b)) {

        _octants_push_back (r_octants,
                            code,
                            0,
                            0);
      }

      else if ((PDM_morton_ancestor_is (code, b) ||
                PDM_morton_ancestor_is (code, a)) &&
               !((code.X[0] == a.X[0]) &&
                 (code.X[1] == a.X[1]) &&
                 (code.X[2] == a.X[2]) &&
                 (code.L == a.L)) &&
               !((code.X[0] == b.X[0]) &&
                 (code.X[1] == b.X[1]) &&
                 (code.X[2] == b.X[2]) &&
                 (code.L == b.L))) {

        PDM_morton_get_children(dim,
                                code,
                                children);

        for (int i = n_child - 1; i >= 0; i--) {
          int is_pushed = _heap_push (heap,
                                      children[i],
                                      0,
                                      0);

          if (!is_pushed) {
            printf ("Internal error PDM_para_octree 2 : heap is full\n");
            exit(1);
          }
        }
      }
    }

    _octants_free (w_octants);

    _heap_free (heap);
  }

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
  size_t *send_shift = malloc(sizeof(size_t) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  size_t *recv_shift = malloc(sizeof(size_t) * (n_ranks+1));

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

  const int _dim = L->dim;
  _octants_purge (L);
  _octants_init (L, _dim, recv_shift[n_ranks]/(_dim + 1));

  int idx = 0;
  for (int i = 0; i < recv_shift[n_ranks]/4; i++) {
    PDM_morton_code_t _code;
    _code.L = recv_codes[idx++];
    for (int j = 0; j < L->dim; j++) {
     _code.X[j] = recv_codes[idx++];
    }
    _octants_push_back (L,
                        _code,
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
  PDM_MPI_Comm_rank (comm, &rank);

  /* Remove duplicates */

  _l_octant_t *L1 = _remove_duplicates (L);

  /* Linearize */

  _l_octant_t *L2 = _linearize (L1);

  _octants_free (L1);

  PDM_g_num_t _n_nodes_global = 0;
  PDM_g_num_t _n_nodes_local =  L2->n_nodes;

  PDM_MPI_Allreduce (&_n_nodes_local, &_n_nodes_global, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  _l_octant_t *R = malloc(sizeof(_l_octant_t));
  _octants_init (R, dim, PDM_MAX(L2->n_nodes, 1));

  if (_n_nodes_global > 0) {

    PDM_morton_code_t *L2_morton_index = malloc(sizeof(PDM_morton_code_t) * (n_ranks + 1));

    int *order = malloc (sizeof(int) * L2->n_nodes);
    int *weight = malloc (sizeof(int) * L2->n_nodes);

    for (int i = 0; i < L2->n_nodes; i++) {
      weight[i] = 1;
      order[i] = i;
    }

    PDM_morton_int_t max_level = 0;
    for (int i = 0; i < L2->n_nodes; i++) {
    max_level = PDM_MAX (L2->codes[i].L, max_level);
    }

    PDM_morton_int_t max_max_level;
    PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                      PDM_MPI_UNSIGNED, PDM_MPI_MAX, comm);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\n-------------\nL2 avant d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    PDM_morton_ordered_build_rank_index (dim,
                                         max_max_level,
                                         L2->n_nodes,
                                         L2->codes,
                                         weight,
                                         L2_morton_index,
                                         comm);

    /* printf("\nL2 morton index d\n"); */
    /* for (int i = 0; i < n_ranks + 1; i++) { */
    /*   PDM_morton_dump (3,  L2_morton_index[i]); */
    /* } */
    /* printf("L2 morton, index f\n"); */

    free(weight);
    free(order);

    _distribute_octants (L2, L2_morton_index, comm);

    free (L2_morton_index);

    /* printf("n L2 : %d %d\n" , rank, L2->n_nodes); */
    /* printf("\nL2 d\n"); */
    /* for (int i = 0; i < L2->n_nodes; i++) { */
    /*   PDM_morton_dump (3, L2->codes[i]); */
    /* } */
    /* printf("L2 f\n--------------\n"); */

  /* PDM_MPI_Barrier(comm); */
  /* exit(1); */

    int *rank_n_nodes = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Allgather(&L2->n_nodes, 1, PDM_MPI_INT,
                      rank_n_nodes, 1, PDM_MPI_INT,
                      comm);

    int first_rank = 0;
    while (first_rank < n_ranks-1
           && rank_n_nodes[first_rank] == 0) {
      first_rank++;
    }

    int last_rank = n_ranks-1;
    while (last_rank > 0
           && rank_n_nodes[last_rank] == 0) {
      last_rank--;
    }

    int next_rank = rank + 1;
    while (next_rank < n_ranks-1
           && rank_n_nodes[next_rank] == 0) {
      next_rank++;
    }

    int prev_rank = rank - 1;
    while (prev_rank > 0
           && rank_n_nodes[prev_rank] == 0) {
      prev_rank--;
    }

    /* printf ("[%d] first last next prev : %d %d %d %d\n", rank, */
    /*         first_rank, last_rank, next_rank, prev_rank); */

    if (rank == first_rank && rank_n_nodes[rank] > 0) {
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
                           0);
    }


    if (rank == last_rank && rank_n_nodes[rank] > 0) {
      PDM_morton_code_t root_DLD;

      root_DLD.L = max_morton_level;
      root_DLD.X[0] = (1u << max_morton_level) - 1u;
      root_DLD.X[1] = (1u << max_morton_level) - 1u;
      root_DLD.X[2] = (1u << max_morton_level) - 1u;

      PDM_morton_code_t FINA;
      PDM_morton_nearest_common_ancestor (root_DLD,
                                          L2->codes[L2->n_nodes -1],
                                          &FINA);

      PDM_morton_code_t child[8];
      PDM_morton_get_children(dim,
                              FINA,
                              child);

      _octants_push_back (L2,
                          child[7],
                          0,
                          0);

    }

    unsigned int sbuff[4];
    unsigned int rbuff[4];
    PDM_MPI_Request srequest;
    PDM_MPI_Request rrequest;

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Irecv ((void *) rbuff,
                     4,
                     PDM_MPI_UNSIGNED,
                     next_rank,
                     0,
                     comm,
                     &rrequest);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      assert (L2->n_nodes > 0);
      sbuff[0] = L2->codes[0].L;
      sbuff[1] = L2->codes[0].X[0];
      sbuff[2] = L2->codes[0].X[1];
      sbuff[3] = L2->codes[0].X[2];

      PDM_MPI_Issend ((void *) sbuff,
                      4,
                      PDM_MPI_UNSIGNED,
                      prev_rank,
                      0,
                      comm,
                      &srequest);


    }

    if (rank < last_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&rrequest);
      PDM_morton_code_t code;

      code.L = rbuff[0];
      code.X[0] = rbuff[1];
      code.X[1] = rbuff[2];
      code.X[2] = rbuff[3];

      _octants_push_back (L2,
                          code,
                          0,
                          0);

    }

    if (rank > first_rank && rank_n_nodes[rank] > 0) {

      PDM_MPI_Wait (&srequest);

    }

    for (int i = 0; i < L2->n_nodes - 1; i++) {
      _l_octant_t *A = _complete_region (L2->codes[i], L2->codes[i+1]);

      _octants_push_back (R,
                          L2->codes[i],
                          0,
                          0);

      if (A != NULL) {
        for (int j = 0; j < A->n_nodes; j++) {
          _octants_push_back (R,
                              A->codes[j],
                              0,
                              0);
        }
        _octants_free (A);
      }
    }

    if (rank == last_rank  && rank_n_nodes[rank] > 0) {
      _octants_push_back (R,
                          L2->codes[L2->n_nodes-1],
                          0,
                          0);
    }

    _octants_free (L2);

    free (rank_n_nodes);
  }

  else {
    if (rank == n_ranks - 1) {

      PDM_morton_code_t _code;
      _code.L = 0;
      _code.X[0] = 0;
      _code.X[1] = 0;
      _code.X[2] = 0;

      _octants_push_back (R,
                          _code,
                          0,
                          0);

    }

  }

  /* printf("fin complete_octree\n"); */
  /* fflush(stdout); */

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
 const PDM_morton_int_t max_level,
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

  __points = realloc (__points, sizeof(double)  * 3 * _n_points);

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
                           recv_coords,
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

  _l_octant_t *T = NULL;

  int max_level = -1;
  int min_level = 32;

  int comm_rank;
  PDM_MPI_Comm_rank(comm, &comm_rank);

  /* printf("\n_block_partition : octant_list %d : d\n", comm_rank); */
  /* if (octant_list!=NULL) { */
  /*   printf("\n_nodes : %d\n", octant_list->n_nodes); */
  /*   for (int i = 0; i < octant_list->n_nodes; i++) { */
  /*     PDM_morton_dump (3, octant_list->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("octant_list NULL\n"); */
  /* } */

  /* printf("_block_partition :octant_list f\n\n"); */

  if (octant_list->n_nodes > 1 ) {

    T = _complete_region (octant_list->codes[0],
                          octant_list->codes[octant_list->n_nodes-1]);

    if (T !=NULL) {
      for (int i = 0; i < T->n_nodes; i++) {
        max_level = PDM_MAX (T->codes[i].L, max_level);
        min_level = PDM_MIN (T->codes[i].L, min_level);
      }
    }
  }

  /* printf("\n_block_partition : complete_region %d :d\n", comm_rank); */
  /* if (T!=NULL) { */
  /*   printf("\n_nodes : %d\n", T->n_nodes); */
  /*   for (int i = 0; i < T->n_nodes; i++) { */
  /*     PDM_morton_dump (3, T->codes[i]); */
  /*   } */
  /* } */
  /* else { */
  /*   printf ("T NULL\n"); */
  /* } */
  /* printf("_block_partition : complete_region f\n\n"); */

  int max_max_level;
  PDM_MPI_Allreduce(&max_level, &max_max_level, 1,
                    PDM_MPI_INT, PDM_MPI_MAX, comm);

  /* Intialize C with coarse octants */

  _l_octant_t *C = malloc(sizeof(_l_octant_t));

  if (T != NULL) {

    _octants_init (C, T->dim, T->n_nodes);

    for (int i = 0; i < T->n_nodes; i++) {

      if (T->codes[i].L <= min_level) {
        _octants_push_back (C,
                            T->codes[i],
                            T->n_points[i],
                            T->range[i]);
      }
    }
  }

  else {
    _octants_init (C, octant_list->dim, 1);
  }

  /* Complete octree */

  /* printf("\n_block_partition : before complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < C->n_nodes; i++) { */
  /*   PDM_morton_dump (3, C->codes[i]); */
  /* } */
  /* printf("_block_partition : before complete_octree f\n\n"); */

  _l_octant_t *G = _complete_octree (C, comm);


  double vol = 0;
  for (int i = 0; i < G->n_nodes; i++) {
    vol += pow(1./pow(2, G->codes[i].L),3);
    G->range[i+1] =
      G->range[i] +
      G->n_points[i];
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
    printf("Erreur volume different de 1 apres complete_octree : %12.5e\n", total_vol);
    for (int i = 0; i < G->n_nodes; i++) {
      PDM_morton_dump (3, G->codes[i]);
    }
  }

  assert (PDM_ABS(total_vol - 1.) < 1e-15);


  /* printf("\n_block_partition : after complete_octree %d %d : d\n", comm_rank, C->n_nodes); */
  /* for (int i = 0; i < G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("_block_partition : after complete_octree f\n\n"); */

  /* PDM_MPI_Barrier (comm); */
  /* exit(1); */
  _octants_free (C);

  if (T != NULL) {
    _octants_free (T);
  }

  /*
   * Compute weight
   */

  /* - exchange codes to ranks (weight per rank)*/

  int n_ranks;
  PDM_MPI_Comm_size(comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank(comm, &rank);

  PDM_morton_int_t *code_buff = malloc (sizeof(PDM_morton_int_t) * (G->dim + 1));
  PDM_morton_int_t *rank_buff = malloc (sizeof(PDM_morton_int_t) * n_ranks * (G->dim + 1));
  int *n_nodes_rank = malloc (sizeof(int) * n_ranks);

  for (int i = 0; i < G->dim + 1; i++) {
    code_buff[i] = 0;
  }

  if ( G->n_nodes > 0) {
    code_buff[0] = G->codes[0].L;
    for (int i = 0; i < G->dim; i++) {
      code_buff[i+1] =  G->codes[0].X[i];
    }
  }

  PDM_MPI_Allgather (&(G->n_nodes), 1, PDM_MPI_INT,
                     n_nodes_rank,  1, PDM_MPI_INT,
                     comm);

  PDM_MPI_Allgather (code_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     rank_buff, G->dim + 1, PDM_MPI_UNSIGNED,
                     comm);

  int n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      n_active_ranks++;
    }
  }

  //assert (n_active_ranks > 0);

  PDM_morton_code_t *rank_codes = malloc (sizeof(PDM_morton_code_t) * n_active_ranks);
  int *active_ranks = malloc (sizeof(int) * n_active_ranks);

  n_active_ranks = 0;
  for (int i = 0; i < n_ranks; i++) {
    if (n_nodes_rank[i] > 0) {
      active_ranks[n_active_ranks] = i;
      rank_codes[n_active_ranks].L = rank_buff[(G->dim + 1) * i];
      for (int j = 0; j < G->dim; j++) {
        rank_codes[n_active_ranks].X[j] = rank_buff[(G->dim + 1) * i + j + 1];
      }
      n_active_ranks++;
    }
  }

  free (code_buff);
  free (rank_buff);

  int *send_count = malloc(sizeof(int) * n_ranks);
  int *send_shift = malloc(sizeof(int) * (n_ranks+1));

  int *recv_count = malloc(sizeof(int) * n_ranks);
  int *recv_shift = malloc(sizeof(int) * (n_ranks+1));

  int irank = 0;
  for (int i = 0; i < n_ranks; i++) {
    send_count[i] = 0;
  }

  /* printf("rank codes deb\n"); */
  /* for (int i = 0; i < n_active_ranks; i++) { */
  /*   PDM_morton_dump(3, rank_codes[i]); */
  /* } */
  /* printf("rank codes fin\n"); */

  for (int i = 0; i < octant_list->n_nodes; i++) {

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }
    send_count[active_ranks[irank]] += octant_list->dim + 2;
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

    if (irank < (n_active_ranks - 1)) {
      if (PDM_morton_a_ge_b (octant_list->codes[i], rank_codes[irank+1])) {

        irank += 1 + PDM_morton_binary_search(n_active_ranks - (irank + 1),
                                              octant_list->codes[i],
                                              rank_codes + irank + 1);
      }
    }

    int shift = send_shift[active_ranks[irank]] + send_count[active_ranks[irank]];

    assert(octant_list->n_points[i] >= 0);

    send_codes[shift++] = octant_list->codes[i].L;

    for (int j = 0; j < octant_list->dim; j++) {
      send_codes[shift++] = octant_list->codes[i].X[j];
    }

    send_codes[shift++] = (PDM_morton_int_t) octant_list->n_points[i];

    send_count[active_ranks[irank]] += octant_list->dim + 2;
  }

  free (rank_codes);

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

  const int _stride = octant_list->dim + 2;
  const int n_recv_codes = recv_shift[n_ranks] / _stride;

  free (recv_shift);

  int *weight = malloc (sizeof(int) * G->n_nodes);

  for (int i = 0; i < G->n_nodes; i++) {
    weight[i] = 0;
  }

  /* - compute weight of each cell */

  for (int i = 0; i < n_recv_codes; i++) {

    PDM_morton_code_t code;

    code.L = recv_codes[i*_stride];

    for (int j = 0; j < octant_list->dim; j++) {
      code.X[j] = recv_codes[i*_stride+j+1];
    }

    int G_node =  PDM_morton_binary_search(G->n_nodes,
                                           code,
                                           G->codes);

    weight[G_node] +=  recv_codes[i*_stride + 1 + octant_list->dim];
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

  PDM_morton_ordered_build_rank_index (octant_list->dim,
                                       max_max_level,
                                       G->n_nodes,
                                       G->codes,
                                       weight,
                                       _G_morton_index,
                                       comm);

  free (order);
  free (weight);

  /* printf("\nblock_partition avant d\n"); */
  /* for (int i = 0; i <  G->n_nodes; i++) { */
  /*   PDM_morton_dump (3, G->codes[i]); */
  /* } */
  /* printf("block_partition avant f\n\n"); */

  _distribute_octants (G, _G_morton_index, comm);

  /*
   * Redistribute octant list from coarse load balancing
   */

  _distribute_octants (octant_list, _G_morton_index, comm);

  free (active_ranks);
  free (n_nodes_rank);
  return G;

}



/**
 *
 * \brief Compute neighbours
 *
 * \param [inout]   octree
 * \param [in]      b_t_elapsed
 * \param [in]      b_t_cpu
 * \param [in]      b_t_cpu_u
 * \param [in]      b_t_cpu_s
 *
 */

static void
_compute_neighbours
(
 _octree_t *octree,
 double   b_t_elapsed,
 double   b_t_cpu,
 double   b_t_cpu_u,
 double   b_t_cpu_s
)
{
  double   e_t_elapsed;
  double   e_t_cpu;
  double   e_t_cpu_u;
  double   e_t_cpu_s;

  const int n_direction = 6;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  _neighbours_tmp_t *neighbours_tmp = malloc (sizeof(_neighbours_tmp_t) * octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j =  0; j < n_direction; j++) {
      neighbours_tmp[i].n_neighbour[j] = 0;
      neighbours_tmp[i].s_neighbour[j] = 1;
      neighbours_tmp[i].neighbours[j] = malloc (sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
      neighbours_tmp[i].neighbours[j][0] = 0;
    }
  }

  /* Boucle sur les noeuds : */
  size_t start_intersect, end_intersect;
  size_t start_intersect_quantile, end_intersect_quantile;

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 1; j < n_direction; j+=2) {
      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t)j);
      PDM_para_octree_direction_t inv_j =
        _inv_direction((PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (octree->rank_octants_index != NULL) {

          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         octree->rank_octants_index,
                                         &start_intersect_quantile,
                                         &end_intersect_quantile);
        }
        else {
          start_intersect_quantile = 0;
          end_intersect_quantile = 1;
        }

        for (int neighbour_rank = start_intersect_quantile;
             neighbour_rank < end_intersect_quantile; neighbour_rank++) {

          if (neighbour_rank == rank) {

            PDM_morton_list_intersect (octree->octants->n_nodes - (i+1),
                                       *neighbour_code,
                                       octree->octants->codes + i + 1,
                                       &start_intersect,
                                       &end_intersect);

            for (int k = start_intersect; k < end_intersect; k++) {
              int idx = k + i + 1;
              PDM_morton_code_t *neighbour_neighbour_code =
                _neighbour (octree->octants->codes[idx], inv_j);

              assert (neighbour_neighbour_code != NULL);

              if (PDM_morton_ancestor_is (octree->octants->codes[i], *neighbour_neighbour_code) ||
                  PDM_morton_ancestor_is (*neighbour_neighbour_code, octree->octants->codes[i])) {

                if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
                  neighbours_tmp[i].s_neighbour[j] *= 2;
                  neighbours_tmp[i].neighbours[j] =
                    realloc (neighbours_tmp[i].neighbours[j],
                             sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
                }
                neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = idx;

                if (neighbours_tmp[idx].n_neighbour[inv_j] >= neighbours_tmp[idx].s_neighbour[inv_j]) {
                  neighbours_tmp[idx].s_neighbour[inv_j] *= 2;
                  neighbours_tmp[idx].neighbours[inv_j] =
                    realloc (neighbours_tmp[idx].neighbours[inv_j],
                             sizeof(int) * neighbours_tmp[idx].s_neighbour[inv_j]);
                }
                neighbours_tmp[idx].neighbours[inv_j][neighbours_tmp[idx].n_neighbour[inv_j]++] = i;
              }

              free (neighbour_neighbour_code);
            }

          }

          else {
            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j+=2) {
      PDM_morton_code_t *neighbour_code =
        _neighbour (octree->octants->codes[i], (PDM_para_octree_direction_t) j);

      if (neighbour_code != NULL) {

        if (octree->rank_octants_index != NULL) {

          PDM_morton_quantile_intersect (n_ranks,
                                         *neighbour_code,
                                         octree->rank_octants_index,
                                         &start_intersect_quantile,
                                         &end_intersect_quantile);
        }
        else {
          start_intersect_quantile = 0;
          end_intersect_quantile = 1;
        }

        for (int neighbour_rank = start_intersect_quantile;
             neighbour_rank < end_intersect_quantile; neighbour_rank++) {

          if (neighbour_rank != rank) {
            if (neighbours_tmp[i].n_neighbour[j] >= neighbours_tmp[i].s_neighbour[j]) {
              neighbours_tmp[i].s_neighbour[j] *= 2;
              neighbours_tmp[i].neighbours[j] = realloc (neighbours_tmp[i].neighbours[j],
                                                         sizeof(int) * neighbours_tmp[i].s_neighbour[j]);
            }
            neighbours_tmp[i].neighbours[j][neighbours_tmp[i].n_neighbour[j]++] = - (neighbour_rank + 1);
          }
        }
        free (neighbour_code);
      }
    }
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Compute connected parts
   *
   *************************************************************************/
  //_compute_connected_parts (octree); // TIMER ??
  int s_connected = 3;
  octree->connected_idx = malloc (sizeof(int) * s_connected);
  octree->connected_idx[0] = 0;

  int *visited = malloc (sizeof(int) * octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++)
    visited[i] = 0;

  int *stack = malloc (sizeof(int) * octree->octants->n_nodes);
  int pos_stack = 0;

  int max = 0;
  while (max < octree->octants->n_nodes) {

    stack[pos_stack++] = max;
    visited[max] = 1;

    while (pos_stack > 0) {
      int i_node = stack[--pos_stack];

      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i_node].n_neighbour[j]; k++) {
          int i_ngb = neighbours_tmp[i_node].neighbours[j][k];
          if (i_ngb >= 0) {
            if (!visited[i_ngb]) {
              stack[pos_stack++] = i_ngb;
              visited[i_ngb] = 1;
              max = PDM_MAX (max, i_ngb);
            }
          }
        }
      }

    }

    max++;
    if (s_connected <= octree->n_connected) {
      s_connected *= 2;
      octree->connected_idx = realloc (octree->connected_idx, sizeof(int) * s_connected);
    }
    octree->connected_idx[++octree->n_connected] = max;
  }
  free (visited);
  free (stack);

  octree->connected_idx = realloc (octree->connected_idx, sizeof(int) * (octree->n_connected+1));


  /*************************************************************************
   *
   * Build parallel partition boundary
   *
   *************************************************************************/

  int FALSE_NEIGHBOUR = -1;
  if (n_ranks > 1) {
    const int n_quantile =  n_ranks * n_direction;

    int *neighbour_rank_n = malloc (sizeof(int) * n_quantile);
    int *neighbour_rank_idx = malloc (sizeof(int) * (n_quantile + 1));

    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_n[i] = 0;
    }

    /* Premiere boucle pour compter */

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
          if (neighbours_tmp[i].neighbours[j][k] < 0) {
            neighbour_rank_n[-(neighbours_tmp[i].neighbours[j][k] + 1)*n_direction +j]++;
          }
        }
      }
    }

    int max_node_dir = -1;
    neighbour_rank_idx[0] = 0;
    for (int i = 0; i < n_quantile; i++) {
      neighbour_rank_idx[i+1] = neighbour_rank_idx[i] + neighbour_rank_n[i];
      max_node_dir = PDM_MAX (max_node_dir, neighbour_rank_n[i]);
      neighbour_rank_n[i] = 0;
    }

    /* Allocation */

    int *neighbour_rank_node_id = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *neighbour_rank_code = malloc (sizeof(PDM_morton_code_t) *
                                                     neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_k = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);
    int *neighbour_rank_node_part = malloc (sizeof(int) * neighbour_rank_idx[n_quantile]);//

    /* Deuxieme boucle pour stocker avec tri suivant la direction */

    for (int i_part = 0; i_part < octree->n_connected; i_part++) {
      for (int i = octree->connected_idx[i_part]; i < octree->connected_idx[i_part+1]; i++) {
        for (int j = 0; j < n_direction; j++) {
          for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
            if (neighbours_tmp[i].neighbours[j][k] < 0) {
              int index = -(neighbours_tmp[i].neighbours[j][k] + 1)*n_direction +j;
              int index2 = neighbour_rank_idx[index] + neighbour_rank_n[index];

              neighbour_rank_node_id[index2] = i;
              PDM_morton_copy (octree->octants->codes[i],
                               neighbour_rank_code + index2);
              neighbour_rank_node_k[index2] = k;
              neighbour_rank_node_part[index2] = i_part;

              neighbour_rank_n[index]++;
            }
          }
        }
      }
    }


    /* Tri des codes pour chaque direction de chaque rang */
    int *order = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_id = malloc (sizeof(int) * max_node_dir);
    PDM_morton_code_t *tmp_code = malloc (sizeof(PDM_morton_code_t) * max_node_dir);
    int *tmp_node_k = malloc (sizeof(int) * max_node_dir);
    int *tmp_node_part = malloc (sizeof(int) * max_node_dir);

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        PDM_morton_local_order (neighbour_rank_n[n_direction * i + j],
                                neighbour_rank_code + neighbour_rank_idx[n_direction * i + j],
                                order);
        int idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (neighbour_rank_code[k], tmp_code + idx1);
          tmp_node_id[idx1]   = neighbour_rank_node_id[k];
          tmp_node_k[idx1]    = neighbour_rank_node_k[k];
          tmp_node_part[idx1] = neighbour_rank_node_part[k];
          idx1 += 1;
        }

        idx1 = 0;
        for (int k = neighbour_rank_idx[n_direction * i + j];
             k < neighbour_rank_idx[n_direction * i + j + 1];
             k++) {
          PDM_morton_copy (tmp_code[order[idx1]], neighbour_rank_code + k );
          neighbour_rank_node_id[k]   = tmp_node_id[order[idx1]];
          neighbour_rank_node_k[k]    = tmp_node_k[order[idx1]];
          neighbour_rank_node_part[k] = tmp_node_part[order[idx1]];
          idx1 += 1;
        }
      }
    }

    free (tmp_code);
    free (order);
    free (tmp_node_id);
    free (tmp_node_k);
    free (tmp_node_part);


    /* Envoi / reception (Les donnees recues sont triees) */

    int *recv_neighbour_rank_n = malloc (sizeof(int) * n_quantile);

    for (int i = 0; i < n_quantile; i++) {
      recv_neighbour_rank_n[i] = 0;
    }

    PDM_MPI_Request *recv_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);
    PDM_MPI_Request *send_request = malloc (sizeof(PDM_MPI_Request) * n_ranks);

    int *used_ranks = malloc (sizeof(int) * n_ranks);

    PDM_MPI_Alltoall (neighbour_rank_n, n_direction, PDM_MPI_INT,
                      recv_neighbour_rank_n, n_direction, PDM_MPI_INT,
                      octree->comm);

    int *recv_neighbour_rank_idx = malloc (sizeof(int) * (n_direction * n_ranks + 1));
    recv_neighbour_rank_idx[0] = 0;


    for (int i = 0; i <  n_direction * n_ranks; i++)
      recv_neighbour_rank_idx[i+1] = recv_neighbour_rank_idx[i] + recv_neighbour_rank_n[i];



    int *recv_neighbour_rank_node_id   = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    int *recv_neighbour_rank_node_part = malloc (sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
    PDM_morton_code_t *recv_neighbour_rank_code =
      malloc (sizeof(PDM_morton_code_t) * recv_neighbour_rank_idx[n_quantile]);


    unsigned int *_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * neighbour_rank_idx[n_quantile]);
    unsigned int *_recv_neighbour_rank_code =
      malloc (sizeof(unsigned int) * 4 * recv_neighbour_rank_idx[n_quantile]);

    int idx = 0;
    for (int i = 0; i < neighbour_rank_idx[n_quantile]; i++) {
      _neighbour_rank_code[idx++] = neighbour_rank_code[i].L;
      for (int j = 0; j < 3; j++) {
        _neighbour_rank_code[idx++] = neighbour_rank_code[i].X[j];
      }
    }

    int *rank_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));
    int *rank_recv_neighbour_rank_n = malloc (sizeof(int) * n_ranks);
    int *rank_recv_neighbour_rank_idx = malloc (sizeof(int) * (n_ranks + 1));

    rank_neighbour_rank_idx[0] = 0;
    rank_recv_neighbour_rank_idx[0] = 0;

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] = 0;
      rank_recv_neighbour_rank_n[i] = 0;
    }

    for (int i = 0; i < n_ranks; i++) {
      for (int j = 0; j < n_direction; j++) {
        rank_neighbour_rank_n[i] += neighbour_rank_n[i*n_direction+j];
        rank_recv_neighbour_rank_n[i] += recv_neighbour_rank_n[i*n_direction+j];
      }
      rank_neighbour_rank_idx[i+1] = rank_neighbour_rank_n[i] + rank_neighbour_rank_idx[i];
      rank_recv_neighbour_rank_idx[i+1] = rank_recv_neighbour_rank_n[i] + rank_recv_neighbour_rank_idx[i];
    }

    PDM_MPI_Alltoallv (neighbour_rank_node_id,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_id,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    PDM_MPI_Alltoallv (neighbour_rank_node_part,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_INT,
                       recv_neighbour_rank_node_part,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_INT,
                       octree->comm);

    for (int i = 0; i < n_ranks; i++) {
      rank_neighbour_rank_n[i] *= 4;
      rank_recv_neighbour_rank_n[i] *= 4;
      rank_neighbour_rank_idx[i+1] *= 4;
      rank_recv_neighbour_rank_idx[i+1] *= 4;
    }


    PDM_MPI_Alltoallv (_neighbour_rank_code,
                       rank_neighbour_rank_n,
                       rank_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       _recv_neighbour_rank_code,
                       rank_recv_neighbour_rank_n,
                       rank_recv_neighbour_rank_idx,
                       PDM_MPI_UNSIGNED,
                       octree->comm);


    free (_neighbour_rank_code);

    free (rank_neighbour_rank_n);
    free (rank_neighbour_rank_idx);
    free (rank_recv_neighbour_rank_n);
    free (rank_recv_neighbour_rank_idx);

    idx = 0;
    for (int i = 0; i < recv_neighbour_rank_idx[n_quantile]; i++) {
      recv_neighbour_rank_code[i].L = _recv_neighbour_rank_code[idx++];
      for (int j = 0; j < 3; j++) {
        recv_neighbour_rank_code[i].X[j] = _recv_neighbour_rank_code[idx++];
      }
    }
    free (_recv_neighbour_rank_code);

    free (recv_request);
    free (send_request);

    free (used_ranks);


    octree->n_part_boundary_elt = neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt_idx = malloc (sizeof(int) * (octree->n_part_boundary_elt + 1));

    int s_part_boundary_elt = 2 * 3 * neighbour_rank_idx[n_quantile];
    octree->part_boundary_elt = malloc (sizeof(int) * s_part_boundary_elt);

    int n_part_boundary_elt = 0;

    idx = 0;
    int idx_part_boundary_elt = 0;
    for (int i = 0; i <= octree->n_part_boundary_elt; i++)
      octree->part_boundary_elt_idx[i] = 0;


    FALSE_NEIGHBOUR = -(octree->n_part_boundary_elt + 1);

    for (int i = 0; i < octree->octants->n_nodes; i++) {
      for (int j = 0; j < n_direction; j++) {
        for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
          if (neighbours_tmp[i].neighbours[j][k] < 0)
            neighbours_tmp[i].neighbours[j][k] = FALSE_NEIGHBOUR;
        }
      }
    }


    /*int *intersect = malloc(sizeof(int) * recv_neighbour_rank_idx[n_quantile]);
      size_t n_intersect;*/

    for (int i = 0; i < n_ranks; i++ ) {
      for (PDM_para_octree_direction_t j = PDM_BOTTOM; j < n_direction; j++) {
        PDM_para_octree_direction_t inv_j = _inv_direction(j);

        int idx_recv = n_direction * i + inv_j;

        int idx_candidate = recv_neighbour_rank_idx[idx_recv];
        int n_candidate = recv_neighbour_rank_idx[idx_recv+1] - idx_candidate;

        if (n_candidate > 0) {

          for (int k = neighbour_rank_idx[i * n_direction + j];
                   k < neighbour_rank_idx[i * n_direction + j + 1]; k++) {
            PDM_morton_code_t *neighbour_code = _neighbour (neighbour_rank_code[k], j);

            PDM_morton_list_intersect (n_candidate,
                                       *neighbour_code,
                                       recv_neighbour_rank_code + idx_candidate,
                                       &start_intersect,
                                       &end_intersect);

            int n_intersect_neighbours = 0;

            if (end_intersect > start_intersect) {
              for (int k1 = start_intersect; k1 < end_intersect; k1++) {
                int k2 = idx_candidate + k1;

                PDM_morton_code_t *neighbour_neighbour_code =
                  _neighbour (recv_neighbour_rank_code[k2], inv_j);

                assert (neighbour_neighbour_code != NULL);

                if (PDM_morton_ancestor_is (neighbour_rank_code[k], *neighbour_neighbour_code) ||
                    PDM_morton_ancestor_is (*neighbour_neighbour_code, neighbour_rank_code[k])) {
                  n_intersect_neighbours++;

                  if ((s_part_boundary_elt - idx_part_boundary_elt) <= 3) {
                    s_part_boundary_elt *= 2;
                    octree->part_boundary_elt = realloc (octree->part_boundary_elt,
                                                         sizeof(int) * s_part_boundary_elt);
                  }
                  octree->part_boundary_elt_idx[n_part_boundary_elt+1]++;
                  octree->part_boundary_elt[idx_part_boundary_elt++] = i; // rank
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_id[k2]; // neighbour's local number in rank i
                  octree->part_boundary_elt[idx_part_boundary_elt++] = recv_neighbour_rank_node_part[k2]; // neighbour's part number in rank i
                }

                free (neighbour_neighbour_code);
              }
            }

            int k3 = neighbour_rank_node_k[k];
            int i2 = neighbour_rank_node_id[k];

            if (n_intersect_neighbours > 0) {
              neighbours_tmp[i2].neighbours[j][k3] = -(n_part_boundary_elt+1);

              assert (neighbours_tmp[i2].neighbours[j][k3] != FALSE_NEIGHBOUR);

              n_part_boundary_elt++;
            }
            else {
              neighbours_tmp[i2].neighbours[j][k3] = FALSE_NEIGHBOUR;
            }

            free (neighbour_code);
          }
        }
      }
    }

    //free (intersect);

    free (neighbour_rank_n);
    free (neighbour_rank_idx);
    free (neighbour_rank_node_id);
    free (neighbour_rank_node_k);
    free (neighbour_rank_node_part);
    free (neighbour_rank_code);

    free (recv_neighbour_rank_n);
    free (recv_neighbour_rank_idx);

    free (recv_neighbour_rank_node_id);
    free (recv_neighbour_rank_node_part);
    free (recv_neighbour_rank_code);

    octree->n_part_boundary_elt = n_part_boundary_elt;
  }

  for (int i = 0; i < octree->n_part_boundary_elt; i++)
    octree->part_boundary_elt_idx[i+1] += octree->part_boundary_elt_idx[i];



  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_DISTANT_NEIGHBOURS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_DISTANT_NEIGHBOURS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_DISTANT_NEIGHBOURS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);



  /*************************************************************************
   *
   * Copy temporary neighbours in the neighbour structure
   *
   *************************************************************************/

  octree->octants->neighbour_idx =
    malloc(sizeof(int) * (n_direction * octree->octants->n_nodes + 1));

  int idx = 0;
  octree->octants->neighbour_idx[0] = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      octree->octants->neighbour_idx[idx+1] =
        octree->octants->neighbour_idx[idx] + neighbours_tmp[i].n_neighbour[j];

      /* account for false distant neighbours */
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {
        if (neighbours_tmp[i].neighbours[j][k] == FALSE_NEIGHBOUR)
          octree->octants->neighbour_idx[idx+1]--;
      }

      idx += 1;
    }
  }

  octree->octants->neighbours =
    malloc(sizeof(int) *
           octree->octants->neighbour_idx[n_direction * octree->octants->n_nodes]);

  idx = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      for (int k = 0; k < neighbours_tmp[i].n_neighbour[j]; k++) {

        if (neighbours_tmp[i].neighbours[j][k] != FALSE_NEIGHBOUR)
          octree->octants->neighbours[idx++] = neighbours_tmp[i].neighbours[j][k];

      }
    }
  }

  /* Free temporary arrays */
  /* printf("sortie 2 neighbours_tmp debut\n"); */
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    for (int j = 0; j < n_direction; j++) {
      if (neighbours_tmp[i].neighbours[j] != NULL) {
        free (neighbours_tmp[i].neighbours[j]);
      }
    }
  }
  /* printf("sortie 2 neighbours_tmp fin\n"); */

  free (neighbours_tmp);




  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3]   += e_t_cpu_s - b_t_cpu_s;


  octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS] += octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_elapsed[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_u[BUILD_LOCAL_NEIGHBOURS_STEP3];

  octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS] += octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP1]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP2]
    + octree->times_cpu_s[BUILD_LOCAL_NEIGHBOURS_STEP3];

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);
}




static void
_check_neighbours_area
(
 const _octree_t *octree
 )
{
  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);


  _l_octant_t *octants = octree->octants;
  double *area = malloc (sizeof(double) * octants->n_nodes);

  int *rank_ngb_n = malloc (sizeof(double) * lComm);
  for (int i = 0; i < lComm; i++)
    rank_ngb_n[i] = 0;

  for (int i = 0; i < octants->n_nodes; i++) {
    PDM_morton_code_t code = octants->codes[i];

    /* check neighbours of current octant */
    area[i] = 0;
    for (int j = 0; j < 6; j++) {
      for (int k = octants->neighbour_idx[6*i+j];
           k < octants->neighbour_idx[6*i+j+1]; k++) {
        int ingb = octants->neighbours[k];

        if (ingb < 0) {
          // distant neighbour
          ingb = -(ingb + 1);
          for (int l = octree->part_boundary_elt_idx[ingb];
               l < octree->part_boundary_elt_idx[ingb+1]; l++) {
            int ngb_rank = octree->part_boundary_elt[3*l];
            rank_ngb_n[ngb_rank]++;
          }

        } else {
          // local neighbour
          PDM_morton_code_t ngb_code = octants->codes[ingb];
          double side = 1./pow(2, PDM_MAX(code.L, ngb_code.L));
          area[i] += side * side;
        }
      }
    }
  }

  // MPI communications to check distant neighbours
  int *rank_ngb_idx = NULL;
  int *rank_ngb_id_level = NULL;
  int *recv_rank_ngb_n = NULL;
  int *recv_rank_ngb_idx = NULL;
  int *recv_rank_ngb_id_level = NULL;

  if (lComm > 1) {
    rank_ngb_idx = malloc (sizeof(int) * (lComm + 1));
    rank_ngb_idx[0] = 0;
    for (int i = 0; i < lComm; i++) {
      rank_ngb_idx[i+1] = rank_ngb_idx[i] + 2*rank_ngb_n[i];
      rank_ngb_n[i] = 0;
    }
    rank_ngb_id_level = malloc (sizeof(int) * rank_ngb_idx[lComm]);

    for (int i = 0; i < octants->n_nodes; i++) {
      for (int j = 0; j < 6; j++) {
        for (int k = octants->neighbour_idx[6*i+j];
             k < octants->neighbour_idx[6*i+j+1]; k++) {
          int ingb = octants->neighbours[k];

          if (ingb < 0) {
            ingb = -(ingb + 1);

            for (int l = octree->part_boundary_elt_idx[ingb];
                 l < octree->part_boundary_elt_idx[ingb+1]; l++) {
              int ngb_rank = octree->part_boundary_elt[3*l];
              int ngb_id = octree->part_boundary_elt[3*l+1];
              int idx = rank_ngb_idx[ngb_rank] + rank_ngb_n[ngb_rank];
              rank_ngb_id_level[idx++] = ngb_id;
              rank_ngb_id_level[idx++] = (int) octants->codes[i].L;
              rank_ngb_n[ngb_rank] += 2;
            }
          }
        }
      }
    }

    recv_rank_ngb_n = malloc (sizeof(int) * lComm);
    PDM_MPI_Alltoall (rank_ngb_n, 1, PDM_MPI_INT,
                      recv_rank_ngb_n, 1, PDM_MPI_INT,
                      octree->comm);

    recv_rank_ngb_idx = malloc (sizeof(int) * (lComm + 1));
    recv_rank_ngb_idx[0] = 0;
    for (int i = 0; i < lComm; i++) {
      recv_rank_ngb_idx[i+1] = recv_rank_ngb_idx[i] + recv_rank_ngb_n[i];
    }

    recv_rank_ngb_id_level = malloc (sizeof(int) * recv_rank_ngb_idx[lComm]);
    PDM_MPI_Alltoallv (rank_ngb_id_level, rank_ngb_n, rank_ngb_idx, PDM_MPI_INT,
                       recv_rank_ngb_id_level, recv_rank_ngb_n, recv_rank_ngb_idx, PDM_MPI_INT,
                       octree->comm);

    free (rank_ngb_id_level);
    free (rank_ngb_n);
    free (rank_ngb_idx);


    for (int i = 0; i < lComm; i++) {
      recv_rank_ngb_n[i] /= 2;
      recv_rank_ngb_idx[i+1] /= 2;
    }
    for (int i = 0; i < lComm; i++) {
      for (int j = recv_rank_ngb_idx[i]; j < recv_rank_ngb_idx[i+1]; j++) {
        int id = recv_rank_ngb_id_level[2*j];
        PDM_morton_int_t level = (PDM_morton_int_t) recv_rank_ngb_id_level[2*j+1];
        double side = 1./pow(2, PDM_MAX (level, octants->codes[id].L));

        area[id] += side * side;
      }
    }

    free (recv_rank_ngb_id_level);
    free (recv_rank_ngb_n);
    free (recv_rank_ngb_idx);
  }



  for (int i = 0; i < octants->n_nodes; i++) {
    /* compute exact interior surface area of current octant */
    PDM_morton_code_t code = octants->codes[i];
    double side = 1./pow(2, code.L);
    int ndir = 0;
    for (int j = 0; j < 6; j++) {
      if (_neighbour (code, (PDM_para_octree_direction_t) j) != NULL)
        ndir++;
    }
    double exact_area = ndir * side * side;

    /* compare with actual area */
    /*printf("node #%d: exact area = %f, actual area = %f, err = %f\n",
           i,
           exact_area,
           area[i],
           PDM_ABS(area[i]/exact_area - 1));*/
    //-->>
    /*if (PDM_ABS(area[i]/exact_area - 1) > 1e-15) {
      printf("[%d] node %d, level %u, area = %f, exact = %f, relative area error = %f\n",
             myRank, i, octants->codes[i].L, area[i], exact_area, area[i]/exact_area - 1);
      printf("\tneighbours (rank, node_id):\n");
      for (int j = 0; j < 6; j++) {
        printf("\t\tdirection %d:", j);
        for (int k = octants->neighbour_idx[6*i+j];
             k < octants->neighbour_idx[6*i+j+1]; k++) {
          int ingb = octants->neighbours[k];
          if (ingb < 0) {
            ingb = -(ingb + 1);

            for (int l = octree->part_boundary_elt_idx[ingb];
                 l < octree->part_boundary_elt_idx[ingb+1]; l++) {
              int ngb_rank = octree->part_boundary_elt[2*l];
              int ngb_id = octree->part_boundary_elt[2*l+1];
              printf(" (%d, %d)", ngb_rank, ngb_id);
            }
          } else {
            printf(" (%d, %d)", myRank, ingb);
          }
        }
        printf("\n");
      }
      printf("\n");
      }*/
    //<<--
    assert (PDM_ABS(area[i]/exact_area - 1) < 1e-15);
  }

  free (area);
}




static void
_closest_points_local
(
 const _octree_t  *octree,
 const int         n_closest_points,
 const int         n_pts,
 double           *pts_coord,
 PDM_g_num_t      *pts_g_num,
 int              *start_leaf,
 double           *upper_bound_dist,
 PDM_g_num_t      *local_closest_src_gnum,
 double           *local_closest_src_dist,
 int              *send_count,
 int              *send_shift,
 int             **send_tgt_lnum,
 int             **send_ngb_id
 )
{
  const int DEBUG = 0;


  _l_octant_t *octants = octree->octants;

  int dim = octree->dim;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);


  int s_heap = octants->n_nodes;
  _min_heap_t *heap = _min_heap_create (s_heap);

  int *is_visited     = malloc (sizeof(int) * octants->n_nodes);
  int *visited_leaves = malloc (sizeof(int) * octants->n_nodes);// could be smaller...
  int visited_count = 0;
  for (int i = 0; i < octants->n_nodes; i++)
    is_visited[i] = 0;

  int *is_receiver    = malloc (sizeof(int) * lComm);
  int *receiver_ranks = malloc (sizeof(int) * lComm);
  int receiver_count = 0;
  for (int i = 0; i < lComm; i++)
    is_receiver[i] = 0;



  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < n_closest_points; j++) {
      local_closest_src_dist[n_closest_points*i + j] = upper_bound_dist[i];
      local_closest_src_gnum[n_closest_points*i + j] = -1; // USEFUL ONLY FOR DEBUG
    }
  }


  for (int i = 0; i < lComm; i++)
    send_count[i] = 0;


  int **send_tgt_leaf = malloc (sizeof(int *) * lComm);
  for (int i = 0; i < lComm; i++) {
    if (i != myRank)
      send_tgt_leaf[i] = malloc (sizeof(int) * n_pts * 2);// could be smaller and dynamically re-alloc'ed when necesseray
  }


  double min_dist;

  for (int i_tgt = 0; i_tgt < n_pts; i_tgt++) {

    for (int i = 0; i < visited_count; i++)
      is_visited[visited_leaves[i]] = 0;
    visited_count = 0;

    for (int i = 0; i < receiver_count; i++)
      is_receiver[receiver_ranks[i]] = 0;
    receiver_count = 0;


    /* Get current target point data */
    const double *_pt = pts_coord + i_tgt * dim;

    PDM_g_num_t *closest_src_gnum = local_closest_src_gnum + n_closest_points * i_tgt;
    double      *closest_src_dist = local_closest_src_dist + n_closest_points * i_tgt;

    //--->>>
    if (DEBUG) {
      printf("\t\tpoint %d (%ld) BEFORE\n", i_tgt, pts_g_num[i_tgt]);
      for (int i = 0; i < n_closest_points; i++) {
        printf("\t\t\t%ld\t%f\n", closest_src_gnum[i], closest_src_dist[i]);
      }
      printf("\n");
    }
    //<<<---


    /* Push start leaf in priority queue */
    _min_heap_reset (heap);

    int i_leaf = start_leaf[i_tgt];
    _min_heap_push (heap,
                    i_leaf,
                    0,
                    0.);

    is_visited[i_leaf] = 1;
    visited_leaves[visited_count++] = i_leaf;

    /* Visit octree leaves using a priority queue */
    PDM_g_num_t unused;
    while (_min_heap_pop (heap, &i_leaf, &unused, &min_dist)) {

      if (min_dist > closest_src_dist[n_closest_points - 1]) break;

      /* inspect source points inside popped leaf */
      for (int i = 0; i < octants->n_points[i_leaf]; i++) {

        // get source point coords and gnum
        int i_src = octants->range[i_leaf] + i;
        double *src_pt = octree->points + i_src * dim;
        PDM_g_num_t src_gnum = octree->points_gnum[i_src];

        // compute (squared) distance from target point
        double src_dist = 0;
        for (int j = 0; j < dim; j++) {
          double delta = _pt[j] - src_pt[j];
          src_dist += delta * delta;
        }

        // insertion sort
        if (src_dist < closest_src_dist[n_closest_points-1]) {
          int j = n_closest_points - 1;
          while (j > 0 && src_dist < closest_src_dist[j-1]) {
            closest_src_gnum[j] = closest_src_gnum[j-1];
            closest_src_dist[j] = closest_src_dist[j-1];
            j--;
          }

          closest_src_gnum[j] = src_gnum;
          closest_src_dist[j] = src_dist;
        }
      }



      /* inspect neighbours of popped leaf */
      for (int i = octants->neighbour_idx[6*i_leaf];
           i < octants->neighbour_idx[6*(i_leaf+1)]; i++) {
        int ngb_id = octants->neighbours[i];

        if (ngb_id < 0) {
          // distant neighbour
          ngb_id = -(ngb_id + 1);

          for (int j = octree->part_boundary_elt[ngb_id];
               j < octree->part_boundary_elt[ngb_id+1]; j++) {
            int ngb_rank = octree->part_boundary_elt[2*j];
            ngb_id       = octree->part_boundary_elt[2*j+1];

            // do not send the same point to the same rank more than once
            if (is_receiver[ngb_rank] == 0) {
              is_receiver[ngb_rank] = 1;
              receiver_ranks[receiver_count++] = ngb_rank;

              // store i_tgt and ngb_id
              send_tgt_leaf[ngb_rank][2*send_count[ngb_rank]] = i_tgt;
              send_tgt_leaf[ngb_rank][2*send_count[ngb_rank]+1] = ngb_id;

              send_count[ngb_rank]++;
            }
          }

        } else {
          // local neighbour
          if (is_visited[ngb_id] == 1) continue;

          is_visited[ngb_id] = 1;
          visited_leaves[visited_count++] = ngb_id;

          // compute min dist from target point to current neighbor leaf
          double ngb_min_dist = _octant_min_dist2 (dim,
                                                   octants->codes[ngb_id],
                                                   octree->d,
                                                   octree->s,
                                                   _pt);

          if (ngb_min_dist < closest_src_dist[n_closest_points - 1]) {
            // push current neighbour in priority queue
            _min_heap_push (heap,
                            ngb_id,
                            0,
                            ngb_min_dist);
          }

        }
      }
    }

    //--->>>
    if (DEBUG) {
      printf("\t\tpoint %d (%ld) AFTER\n", i_tgt, pts_g_num[i_tgt]);
      for (int i = 0; i < n_closest_points; i++) {
      printf("\t\t\t%ld\t%f\n", closest_src_gnum[i], closest_src_dist[i]);
    }
      printf("\n");
    }
    //<<<---

  }

  _min_heap_free (heap);
  free (is_visited);
  free (visited_leaves);
  free (is_receiver);
  free (receiver_ranks);


  send_shift[0] = 0;
  for (int i = 0; i < lComm; i++)
    send_shift[i+1] = send_shift[i] + send_count[i];

  *send_tgt_lnum = malloc (sizeof(int) * send_shift[lComm]);
  *send_ngb_id   = malloc (sizeof(int) * send_shift[lComm]);


  int idx = 0;
  for (int i = 0; i < lComm; i++) {
    for (int j = 0; j < send_count[i]; j++) {
      (*send_tgt_lnum)[idx] = send_tgt_leaf[i][2*j];
      (*send_ngb_id)[idx]   = send_tgt_leaf[i][2*j+1];
      idx++;
    }

    if (i != myRank)
      free (send_tgt_leaf[i]);
  }
  free (send_tgt_leaf);



  if (DEBUG) {
    printf("\n\nSend count =");
    for (int i = 0; i < lComm; i++)
      printf(" %d", send_count[i]);
    printf("\n");
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud          Number of point cloud
 * \param [in]   depth_max              Maximum depth
 * \param [in]   points_in_leaf_max     Maximum points in a leaf
 * \param [in]   build_leaf_neighbours  Build leaf nieghbours (1 = true)
 * \param [in]   comm                   MPI communicator
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

  octree->rank_octants_index = NULL;
  octree->octants = NULL;

  octree->n_part_boundary_elt = 0;
  octree->part_boundary_elt_idx = NULL;
  octree->part_boundary_elt = NULL;

  octree->neighboursToBuild = build_leaf_neighbours;

  octree->comm = comm;

  octree->n_connected = 0;
  octree->connected_idx = NULL;

  octree->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    octree->times_elapsed[i] = 0.;
    octree->times_cpu[i] = 0.;
    octree->times_cpu_u[i] = 0.;
    octree->times_cpu_s[i] = 0.;
  }

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

  if (octree->part_boundary_elt_idx != NULL) {
    free (octree->part_boundary_elt_idx);
  }

  if (octree->part_boundary_elt != NULL) {
    free (octree->part_boundary_elt);
  }

  if (octree->rank_octants_index != NULL) {
    free (octree->rank_octants_index);
  }

  if (octree->octants != NULL) {

    if (octree->octants->codes != NULL) {
      free (octree->octants->codes);
    }

    if (octree->octants->range != NULL) {
      free (octree->octants->range);
    }

    if (octree->octants->n_points != NULL) {
      free (octree->octants->n_points);
    }

    if (octree->octants->neighbour_idx != NULL) {
      free (octree->octants->neighbour_idx);
    }

    if (octree->octants->neighbours != NULL) {
      free (octree->octants->neighbours);
    }

    free (octree->octants);
  }

  if (octree->connected_idx != NULL) {
    free (octree->connected_idx);
  }

  PDM_timer_free (octree->timer);

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
  const PDM_morton_int_t max_level = PDM_morton_max_level;

  int n_ranks;
  PDM_MPI_Comm_size (octree->comm, &n_ranks);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  octree->times_elapsed[BEGIN] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BEGIN]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(octree->timer);

  b_t_elapsed = octree->times_elapsed[BEGIN];
  b_t_cpu     = octree->times_cpu[BEGIN];
  b_t_cpu_u   = octree->times_cpu_u[BEGIN];
  b_t_cpu_s   = octree->times_cpu_s[BEGIN];
  PDM_timer_resume(octree->timer);

  /*
   * Get coord extents
   */

  PDM_morton_get_coord_extents(dim,
                               octree->n_points,
                               octree->points,
                               octree->global_extents,
                               octree->comm);

  /*
   * Dilate extents
   */
  const double EPS_range  = 1.e-6;
  const double EPS_double = 1.e-12;
  for (int i = 0; i < dim; i++) {
    double range = octree->global_extents[i+dim] - octree->global_extents[i];
    double epsilon = PDM_MAX (EPS_double, range * EPS_range);
    octree->global_extents[i]     -= epsilon;
    octree->global_extents[i+dim] += epsilon;
  }

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

    free (octree->points);
    octree->points = _points;

    free (order);
  }


  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_ORDER_POINTS] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_ORDER_POINTS]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_ORDER_POINTS]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_ORDER_POINTS]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

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

    _octants_init (point_octants, octree->dim, octree->n_points);

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
                                        &octree->rank_octants_index);

    _octants_free (point_octants);

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
                        octree->rank_octants_index,
                        octree->comm,
                        octree->dim,
                        max_level,
                        octree->global_extents);

    int iblock = 0;
    for (int i = 0; i < octree->n_points; i++) {
      while (!PDM_morton_ancestor_is (octree->octants->codes[iblock],
                                      octree->points_code[i])) {
        iblock++;
      }
      assert (iblock < octree->octants->n_nodes);
      octree->octants->n_points[iblock] += 1;
    }

    octree->octants->range[0] = 0;

    double vol = 0;
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      vol += pow(1./pow(2, octree->octants->codes[i].L),3);
      octree->octants->range[i+1] =
        octree->octants->range[i] +
        octree->octants->n_points[i];
    }
    double total_vol;
    PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

    if ( (PDM_ABS(total_vol - 1.)>= 1e-15)) {
      printf("Erreur volume different de 1 : %12.5e\n", total_vol);
    }

    assert (PDM_ABS(total_vol - 1.) < 1e-15);

  }

  else {

    octree->octants = malloc(sizeof(_l_octant_t));

    _octants_init (octree->octants, octree->dim, octree->n_points);

    PDM_morton_code_t code;

    code.L = 0;
    code.X[0] = 0;
    code.X[1] = 0;
    code.X[2] = 0;

    _octants_push_back (octree->octants,
                        code,
                        octree->n_points,
                        0);

    octree->octants->range[0] = 0;
    octree->octants->range[1] = octree->n_points;
    octree->octants->n_points[0] = octree->n_points;

  }

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_BLOCK_PARTITION] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_BLOCK_PARTITION]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_BLOCK_PARTITION]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_BLOCK_PARTITION]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Build local octree
   *
   *************************************************************************/

  const int n_child = 8;
  //const int n_direction = 6;

  int  size = octree->depth_max * 8;
  _heap_t *heap = _heap_create (size);
  for (int i = octree->octants->n_nodes - 1; i >= 0; i--) {
    int is_pushed = _heap_push (heap,
                                octree->octants->codes[i],
                                octree->octants->range[i],
                                octree->octants->n_points[i]);
    if (!is_pushed) {
      printf ("Internal error PDM_para_octree 3 : heap is full\n");
      exit(1);
    }
  }

  PDM_morton_code_t code;
  int range;
  int n_points;

  octree->octants->n_nodes = 0;

  while (_heap_pull (heap, &code, &range, &n_points)) {

    /* Add children into the heap*/

    if ((code.L < max_morton_level) && (code.L < max_level) &&
        (n_points > octree->points_in_leaf_max)) {

      PDM_morton_code_t children[n_child];
      PDM_morton_get_children(dim,
                              code,
                              children);

      int range_children[n_child];
      int n_points_children[n_child];

      for (int i = 0; i < n_child; i++) {
        n_points_children[i] = 0;
      }

      int ichild = 0;
      for (int i = 0; i < n_points; i++) {
        assert ((range + i) < octree->n_points);
        if (!PDM_morton_ancestor_is(code, octree->points_code[range + i])) {
          printf("Erreur : n'est pas un ancetre !!!!!\n");

        }
        assert (PDM_morton_ancestor_is(code, octree->points_code[range + i]));
        while (!PDM_morton_ancestor_is (children[ichild], octree->points_code[range + i])) {
          ichild += 1;
        }
        assert (ichild < n_child);
        n_points_children[ichild] += 1;
      }

      range_children[0] = 0;
      for (int i = 0; i < n_child - 1; i++) {
        range_children[i+1] = range_children[i] + n_points_children[i];
      }

      for (int i = n_child - 1; i >= 0; i--) {
        int is_pushed = _heap_push (heap,
                                    children[i],
                                    range + range_children[i],
                                    n_points_children[i]);
        if (!is_pushed) {
          printf ("Internal error PDM_para_octree 4 : heap is full\n");
          exit(1);
        }
      }

    }

    /* Store the leaf */

    else {
      _octants_push_back (octree->octants, code, n_points, range);
    }

  }


  double vol = 0;
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    vol += pow(1./pow(2, octree->octants->codes[i].L),3);
  }
  double total_vol;
  PDM_MPI_Allreduce(&vol, &total_vol, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, octree->comm);

  assert (PDM_ABS(total_vol - 1.) < 1e-15);

  heap = _heap_free (heap);

  PDM_timer_hang_on(octree->timer);
  e_t_elapsed = PDM_timer_elapsed(octree->timer);
  e_t_cpu     = PDM_timer_cpu(octree->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(octree->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(octree->timer);

  octree->times_elapsed[BUILD_LOCAL_NODES] += e_t_elapsed - b_t_elapsed;
  octree->times_cpu[BUILD_LOCAL_NODES]     += e_t_cpu - b_t_cpu;
  octree->times_cpu_u[BUILD_LOCAL_NODES]   += e_t_cpu_u - b_t_cpu_u;
  octree->times_cpu_s[BUILD_LOCAL_NODES]   += e_t_cpu_s - b_t_cpu_s;

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(octree->timer);

  /*************************************************************************
   *
   * Neighbours
   *
   *************************************************************************/

  if (octree->neighboursToBuild) {
    _compute_neighbours (octree,
                         b_t_elapsed,
                         b_t_cpu,
                         b_t_cpu_u,
                         b_t_cpu_s);

#if 1
    _check_neighbours_area (octree);
#endif
  }

  PDM_timer_hang_on(octree->timer);
  octree->times_elapsed[BUILD_TOTAL] = PDM_timer_elapsed(octree->timer);
  octree->times_cpu[BUILD_TOTAL]     = PDM_timer_cpu(octree->timer);
  octree->times_cpu_u[BUILD_TOTAL]   = PDM_timer_cpu_user(octree->timer);
  octree->times_cpu_s[BUILD_TOTAL]   = PDM_timer_cpu_sys(octree->timer);
  PDM_timer_resume(octree->timer);

  octree->times_elapsed[END] = octree->times_elapsed[BUILD_TOTAL];
  octree->times_cpu[END]     = octree->times_cpu[BUILD_TOTAL];
  octree->times_cpu_u[END]   = octree->times_cpu_u[BUILD_TOTAL];
  octree->times_cpu_s[END]   = octree->times_cpu_s[BUILD_TOTAL];

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
 * \brief Dump octree
 *
 * \param [in]   id                 Identifier
 *
 */

void
PDM_para_octree_dump
(
 const int          id
)
{
  _octree_t *octree = _get_from_id (id);
  PDM_printf("PDM_dump_para_octree : %d\n",id);

  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  printf("  - global_extents :");
  for (int i = 0; i < 2*octree->dim; i++) {
    printf(" %12.5e", octree->global_extents[i]);
  }
  PDM_printf("\n");
  PDM_printf("  - depth_max : %d\n", octree->depth_max);
  PDM_printf("  - points_in_leaf_max : %d\n", octree->points_in_leaf_max);

  PDM_printf("  - s : %12.5e %12.5e %12.5e\n", octree->s[0], octree->s[1], octree->s[2]);
  PDM_printf("  - d : %12.5e %12.5e %12.5e\n", octree->d[0], octree->d[1], octree->d[2]);

  PDM_printf("  - n_point_clouds : %d\n", octree->n_point_clouds);
  PDM_printf("  - t_n_points : "PDM_FMT_G_NUM"\n", octree->t_n_points);
  PDM_printf("  - n_points : %d\n", octree->n_points);
  for (int i = 0; i < octree->n_points; i++) {
    PDM_printf("  %d "PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e / level %u code [%u, %u, %u]\n",
               i, octree->points_gnum[i],
               octree->points[3*i], octree->points[3*i+1], octree->points[3*i+2],
               octree->points_code[i].L,
               octree->points_code[i].X[0],
               octree->points_code[i].X[1],
               octree->points_code[i].X[2]);
  }
  PDM_printf("  - n_nodes : %d\n", octree->octants->n_nodes);
  for (int i = 0; i < octree->octants->n_nodes; i++) {
    PDM_printf("  %d : level %u code [%u, %u, %u], range %d, n_points %d\n",
               i,
               octree->octants->codes[i].L,
               octree->octants->codes[i].X[0],
               octree->octants->codes[i].X[1],
               octree->octants->codes[i].X[2],
               octree->octants->range[i],
               octree->octants->n_points[i]
               );
  }
  if (octree->neighboursToBuild) {
    for (int i = 0; i < octree->octants->n_nodes; i++) {
      PDM_printf("  %d : neighbors\n", i);
      for (int j = 0; j < 6; j++) {
        PDM_printf("    - direction %d : ", j);
        for (int k = octree->octants->neighbour_idx[6*i+j];
             k < octree->octants->neighbour_idx[6*i+j+1]; k++) {
          PDM_printf(" %d",octree->octants->neighbours[k]);
        }
        PDM_printf("\n");
      }
    }
  }
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

#define VERBOSE 0

  _octree_t *octree = _get_from_id (id);
  _l_octant_t *octants = octree->octants;

  int dim = octree->dim;
  int _n_closest_points = n_closest_points;
  int _n_pts = n_pts;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);


  /* Force target points inside octree extents /!\ /!\ /!\ -->> */
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < dim; j++) {
      pts[dim*i+j] = PDM_MAX (pts[dim*i+j], octree->global_extents[j]);
      pts[dim*i+j] = PDM_MIN (pts[dim*i+j], octree->global_extents[dim+j]);
    }
  }
  /* <<-- */


  /* Part-to-block create (only to get block distribution) */
  PDM_part_to_block_t *ptb1 = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        &pts_g_num,
                                                        NULL,
                                                        &_n_pts,
                                                        1,
                                                        octree->comm);

  PDM_g_num_t *block_distrib_idx1 = PDM_part_to_block_distrib_index_get (ptb1);
  const int n_pts_block1 = PDM_part_to_block_n_elt_block_get (ptb1);

  PDM_g_num_t *block_closest_src_gnum1 =
    malloc (sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points);
  double      *block_closest_src_dist1 =
    malloc (sizeof(double)      * n_pts_block1 * n_closest_points);

  //-->>
  for (int i = 0; i < n_pts_block1; i++) {
    for (int j = 0; j < n_closest_points; j++) {
      block_closest_src_dist1[n_closest_points*i + j] = HUGE_VAL;
      block_closest_src_gnum1[n_closest_points*i + j] = -1; // USEFUL ONLY FOR DEBUG
    }
  }
  //<<--

  /*************************************************************************
   *
   * Distribute the target points
   *
   *************************************************************************/
  /*   1) Encode the coordinates of every target point */
  PDM_morton_code_t *pts_code = malloc (sizeof(PDM_morton_code_t) * n_pts);
  double d[3], s[3];
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_pts,
                            pts,
                            pts_code,
                            d,
                            s);

  /*   2) Get Morton code of first leaf for each rank */
  unsigned int send_first_code[4];
  send_first_code[0] = octants->codes[0].L;
  for (int i = 0; i < 3; i++)
    send_first_code[i+1] = octants->codes[0].X[i];

  unsigned int *recv_first_code = malloc (sizeof(int) * 4 * lComm);
  PDM_MPI_Allgather (send_first_code, 4, PDM_MPI_UNSIGNED,
                     recv_first_code, 4, PDM_MPI_UNSIGNED,
                     octree->comm);

  PDM_morton_code_t *first_code = malloc (sizeof(PDM_morton_code_t) * lComm);
  for (int i = 0; i < lComm; i++) {
    first_code[i].L = recv_first_code[4*i];
    for (int j = 0; j < 3; j++)
      first_code[i].X[j] = recv_first_code[1+4*i+j];
  }
  free (recv_first_code);

  /*   3) Use binary search to associate each target point to the right process */
  int *send_count = malloc (sizeof(int) * lComm);
  int *recv_count = malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++)
    send_count[i] = 0;

  int *rank_pt = malloc (sizeof(int) * n_pts);
  for (int i = 0; i < n_pts; i++) {
    rank_pt[i] = PDM_morton_binary_search (lComm,
                                           pts_code[i],
                                           first_code);

    send_count[rank_pt[i]]++;
  }
  free (first_code);
  free (pts_code);

  /*   4) Exchange send/recv counts */
  PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                    recv_count, 1, PDM_MPI_INT,
                    octree->comm);

  int *send_shift = malloc (sizeof(int) * (lComm+1));
  int *recv_shift = malloc (sizeof(int) * (lComm+1));
  send_shift[0] = 0;
  recv_shift[0] = 0;
  for (int i = 0; i < lComm; i++) {
    send_shift[i+1] = send_shift[i] + send_count[i];
    recv_shift[i+1] = recv_shift[i] + recv_count[i];
    send_count[i] = 0;
  }

  /*   5) Fill send buffers */
  PDM_g_num_t *send_g_num = malloc (sizeof(PDM_g_num_t) * send_shift[lComm]);
  PDM_g_num_t *recv_g_num = malloc (sizeof(PDM_g_num_t) * recv_shift[lComm]);
  double      *send_coord = malloc (sizeof(double)      * send_shift[lComm]*dim);
  double      *recv_coord = malloc (sizeof(double)      * recv_shift[lComm]*dim);
  for (int i = 0; i < n_pts; i++) {
    int rank = rank_pt[i];
    int k = send_shift[rank] + send_count[rank];
    send_g_num[k] = pts_g_num[i];
    for (int j = 0; j < dim; j++)
      send_coord[dim*k+j] = pts[dim*i+j];

    send_count[rank]++;
  }
  free (rank_pt);

  /*   6) Send gnum buffer */
  PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                     recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                     octree->comm);

  /*   7) Send coord buffer */
  int n_recv_pts = recv_shift[lComm];
  for (int i = 0; i < lComm; i++) {
    send_count[i] *= dim;
    recv_count[i] *= dim;
    send_shift[i+1] *= dim;
    recv_shift[i+1] *= dim;
  }
  PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                     recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                     octree->comm);




  /* Encode the coordinates of the received target points */
  pts_code = malloc (sizeof(PDM_morton_code_t) * n_recv_pts);
  PDM_morton_encode_coords (dim,
                            PDM_morton_max_level,
                            octree->global_extents,
                            (size_t) n_recv_pts,
                            recv_coord,
                            pts_code,
                            d,
                            s);

  /* Find start leaf for each target point */
  int *start_leaf = malloc (sizeof(int) * n_recv_pts);
  for (int i = 0; i < n_recv_pts; i++) {
    start_leaf[i] = PDM_morton_binary_search (octants->n_nodes,
                                              pts_code[i],
                                              octants->codes);
  }
  free (pts_code);







  int s_processed_tgt_pts = PDM_MAX (128, 2 * n_recv_pts);
  int n_processed_tgt_pts = 0;
  PDM_g_num_t *processed_tgt_pts = malloc (sizeof(PDM_g_num_t) * s_processed_tgt_pts);


  PDM_g_num_t keyMax = 1000;//???
  PDM_hash_tab_t *ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_LONG,
                                            &keyMax);


  PDM_g_num_t *local_closest_src_gnum =
    malloc (sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  double      *local_closest_src_dist =
    malloc (sizeof(double)      * n_recv_pts * n_closest_points);

  double *upper_bound_dist = malloc (sizeof(double) * n_recv_pts);
  for (int i = 0; i < n_recv_pts; i++)
    upper_bound_dist[i] = HUGE_VAL;

  int *send_tgt_lnum = NULL;
  int *send_ngb_id   = NULL;
  int *block_stride2 = NULL;
  double      *block_closest_src_dist2 = NULL;
  PDM_g_num_t *block_closest_src_gnum2 = NULL;
  double      *block_upper_bound_dist = malloc (sizeof(double) * n_pts_block1);
  int one = 1;

  _min_heap_t *merge_heap = _min_heap_create (2); // used to merge part-to-block results

  double      *tmp_closest_src_dist = malloc (sizeof(double)      * n_closest_points);
  PDM_g_num_t *tmp_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);

  int iteration = 0;
  while (1) {
    //-->>
    iteration++;
#if VERBOSE
    printf("[%d] iteration %d\n", myRank, iteration);
#else
    if (myRank == 0)
      printf("iteration %d\n", iteration);
#endif
    if (iteration > 10*lComm) break; // emergency exit
    //<<--



    //--->>>
    /* Clean-up/Filter 'recv' data */
    /*   1) remove duplicate target points */
    /*   2) remove target points already processed in current rank */
    int tmp_n = 0;
    PDM_g_num_t *tmp_g_num       = malloc (sizeof(PDM_g_num_t) * n_recv_pts);
    double      *tmp_coord       = malloc (sizeof(double)      * n_recv_pts * 3);
    double      *tmp_upper_bound = malloc (sizeof(double)      * n_recv_pts);
    int         *tmp_start_leaf  = malloc (sizeof(int)         * n_recv_pts);

    if (s_processed_tgt_pts <= n_processed_tgt_pts + n_recv_pts) {
      while (s_processed_tgt_pts <= n_processed_tgt_pts + n_recv_pts)
        s_processed_tgt_pts *= 2;

      processed_tgt_pts = realloc (processed_tgt_pts, sizeof(PDM_g_num_t) * s_processed_tgt_pts);
    }

    for (int i = 0; i < n_recv_pts; i++) {
      PDM_g_num_t gnum = recv_g_num[i];
      PDM_g_num_t key = gnum % keyMax;
      int n_data = PDM_hash_tab_n_data_get (ht, &key);

      int found = 0;
      PDM_g_num_t **data = (PDM_g_num_t **) PDM_hash_tab_data_get (ht, &key);
      for (int k = 0; k < n_data; k++) {
        if (*(data[k]) == gnum) {
          found = 1;
          break;
        }
      }

      // check if current tgt point has already been processed in current rank
      int already_processed = 0;
      int lo = 0;
      int hi = n_processed_tgt_pts;

      if (n_processed_tgt_pts == 0) {
        processed_tgt_pts[0] = gnum;
        n_processed_tgt_pts++;

      } else {

        while (lo < hi - 1) {
          int mid = lo + (hi - lo)/2;
          if (gnum < processed_tgt_pts[mid])
            hi = mid;
          else
            lo = mid;
        }

        if (processed_tgt_pts[lo] == gnum) {
          already_processed = 1;
        } else {
          lo += 1;
          for (int j = n_processed_tgt_pts; j > lo; j--)
            processed_tgt_pts[j] = processed_tgt_pts[j-1];

          processed_tgt_pts[lo] = gnum;
          n_processed_tgt_pts++;
        }
      }


      if (!found && !already_processed) {
        tmp_g_num[tmp_n] = gnum;

        for (int j = 0; j < 3; j++)
          tmp_coord[3*tmp_n + j] = recv_coord[3*i+j];

        tmp_upper_bound[tmp_n] = upper_bound_dist[i];
        tmp_start_leaf[tmp_n] = start_leaf[i];

        tmp_n++;
      }
    }

    if (tmp_n < n_recv_pts) {
      printf("Filtering points: tmp_n / n_recv_pts = %d / %d\n", tmp_n, n_recv_pts);

      for (int i = 0; i < tmp_n; i++) {
        recv_g_num[i] = tmp_g_num[i];

        for (int j = 0; j < 3; j++)
          recv_coord[3*i+j] = tmp_coord[3*i+j];

        upper_bound_dist[i] = tmp_upper_bound[i];
        start_leaf[i] = tmp_start_leaf[i];
      }

      recv_g_num       = realloc (recv_g_num,       sizeof(PDM_g_num_t) * tmp_n);
      recv_coord       = realloc (recv_coord,       sizeof(double)      * tmp_n * 3);
      upper_bound_dist = realloc (upper_bound_dist, sizeof(double)      * tmp_n);
      start_leaf       = realloc (start_leaf,       sizeof(int)         * tmp_n);
      n_recv_pts = tmp_n;
    }

    free (tmp_g_num);
    free (tmp_coord);
    free (tmp_upper_bound);
    free (tmp_start_leaf);
    //<<<---

    if (myRank == 0)
      printf("Filtering OK\n");

    _closest_points_local (octree,
                           n_closest_points,
                           n_recv_pts,
                           recv_coord,
                           recv_g_num,
                           start_leaf,
                           upper_bound_dist,
                           local_closest_src_gnum,
                           local_closest_src_dist,
                           send_count,
                           send_shift,
                           &send_tgt_lnum,
                           &send_ngb_id);

    if (myRank == 0)
      printf("_closest_points_local OK\n");


    PDM_part_to_block_t *ptb2 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_MERGE,
                                                           1.,
                                                           &recv_g_num,
                                                           block_distrib_idx1,
                                                           &n_recv_pts,
                                                           1,
                                                           octree->comm);


    int n_pts_block2 = PDM_part_to_block_n_elt_block_get (ptb2);

    int *stride2 = malloc (sizeof(int) * n_recv_pts);

    /* Part-to-block exchanges */
    PDM_g_num_t *block_tgt_gnum2 = PDM_part_to_block_block_gnum_get (ptb2);

    for (int i = 0; i < n_recv_pts; i++) {
      stride2[i] = n_closest_points;
    }

    PDM_part_to_block_exch (ptb2,
                            sizeof(double),
                            PDM_STRIDE_VAR,
                            1,
                            &stride2,
                            (void **) &local_closest_src_dist,
                            &block_stride2,
                            (void **) &block_closest_src_dist2);
    free (block_stride2);

    PDM_part_to_block_exch (ptb2,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_VAR,
                            1,
                            &stride2,
                            (void **) &local_closest_src_gnum,
                            &block_stride2,
                            (void **) &block_closest_src_gnum2);
    free (stride2);


    /* Merge block data */
    int *block_idx2 = malloc (sizeof(int) * (n_pts_block2 + 1));
    block_idx2[0] = 0;
    for (int i = 0; i < n_pts_block2; i++) {
      block_idx2[i+1] = block_idx2[i] + block_stride2[i];
    }

#if VERBOSE
    printf("\n\n- - - Merge - - -");
#endif
    for (int i = 0; i < n_pts_block2; i++) {
      int id_block1 = (int) n_closest_points * (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]);

#if VERBOSE
      printf("\nblock2 i = %d (%ld) <%d>\n", i, block_tgt_gnum2[i],
             (int) (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]));
      printf("id_block1 = %d / %d\n", id_block1, n_pts_block1*n_closest_points);
#endif

      int n_procs_pt = block_stride2[i] / n_closest_points;
      int *idx_proc = malloc (sizeof(int) * (n_procs_pt + 1));

      _min_heap_reset (merge_heap);

      for (int j = 0; j < n_closest_points; j++) {
        tmp_closest_src_gnum[j] = block_closest_src_gnum1[id_block1 + j];
        tmp_closest_src_dist[j] = block_closest_src_dist1[id_block1 + j];
      }

#if VERBOSE
      printf("\tj = %d: (%ld, %f)\n",
             -1, tmp_closest_src_gnum[0], tmp_closest_src_dist[0]);
#endif
      _min_heap_push (merge_heap,
                      -1,
                      tmp_closest_src_gnum[0],
                      tmp_closest_src_dist[0]);
      idx_proc[0] = 0;

      for (int j = 0; j < n_procs_pt; j++) {
        int k = block_idx2[i] + j*n_closest_points;
#if VERBOSE
        printf("\tj =  %d: (%ld, %f)\n",
               j, block_closest_src_gnum2[k], block_closest_src_dist2[k]);
#endif

        _min_heap_push (merge_heap,
                        j,
                        block_closest_src_gnum2[k],
                        block_closest_src_dist2[k]);
        idx_proc[j+1] = 0;
      }

      int _proc;
      PDM_g_num_t _gnum;
      double _dist;
      for (int j = 0; j < n_closest_points; j++) {
        int popped = _min_heap_pop (merge_heap,
                                    &_proc,
                                    &_gnum,
                                    &_dist);
        assert (popped);

        if (_dist >= block_closest_src_dist1[id_block1 + j]) break;

        block_closest_src_dist1[id_block1 + j] = _dist;
        block_closest_src_gnum1[id_block1 + j] = _gnum;

#if VERBOSE
        printf("\t\t%d/%d: %ld, %f\n",
               j+1, n_closest_points,
               block_closest_src_gnum1[id_block1 + j], block_closest_src_dist1[id_block1 + j]);
#endif

        if (j >= n_closest_points - 1)
          break;

        idx_proc[_proc+1]++;

        if (_proc < 0) {
          _min_heap_push (merge_heap,
                          -1,
                          tmp_closest_src_gnum[idx_proc[0]],
                          tmp_closest_src_dist[idx_proc[0]]);
        } else {
          int k = block_idx2[i] + _proc*n_closest_points + idx_proc[_proc+1];
          _min_heap_push (merge_heap,
                          _proc,
                          block_closest_src_gnum2[k],
                          block_closest_src_dist2[k]);
        }

      }
      free (idx_proc);
    }
#if VERBOSE
    printf("- - - - - - - - -\n\n");
#endif

    free (block_idx2);
    free (block_stride2);
    free (block_closest_src_gnum2);
    free (block_closest_src_dist2);

    ptb2 = PDM_part_to_block_free (ptb2);

    //-->>
#if VERBOSE
    printf("\n--- Block ---\n");
    for (int i = 0; i < n_pts_block1; i++) {
      printf("%d (%ld)\n", i, block_distrib_idx1[myRank] + i + 1);
      for (int j = 0; j < n_closest_points; j++) {
        printf("\t%ld\t%f\n",
               block_closest_src_gnum1[n_closest_points*i + j],
               block_closest_src_dist1[n_closest_points*i + j]);
      }
    }
    printf("-------------\n\n");
#endif
    //<---


    /* Update upper_bound_dist */
    PDM_block_to_part_t *btp2 = PDM_block_to_part_create (block_distrib_idx1,
                                                          (const PDM_g_num_t **) &recv_g_num,
                                                          &n_recv_pts,
                                                          1,
                                                          octree->comm);

    for (int i = 0; i < n_pts_block1; i++) {
      block_upper_bound_dist[i] = block_closest_src_dist1[n_closest_points*(i+1) - 1];
    }

    PDM_block_to_part_exch (btp2,
                            sizeof(double),
                            PDM_STRIDE_CST,
                            &one,
                            block_upper_bound_dist,
                            NULL,
                            (void **) &upper_bound_dist);
    btp2 = PDM_block_to_part_free (btp2);



    /* Redistribute target points for next iteration */
    // send count
    PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
                      recv_count, 1, PDM_MPI_INT,
                      octree->comm);
    recv_shift[0] = 0;
    for (int i = 0; i < lComm; i++)
      recv_shift[i+1] = recv_shift[i] + recv_count[i];

    n_recv_pts = recv_shift[lComm];

    /* Termination criterion */
    int max_n_recv_pts = 0;
    PDM_MPI_Allreduce (&n_recv_pts, &max_n_recv_pts, 1,
                       PDM_MPI_INT, PDM_MPI_MAX, octree->comm);

    if (max_n_recv_pts == 0) {
      if (myRank == 0)
        printf("*** max_n_recv_pts = 0 --> DONE :-)\n");

      free (send_ngb_id);
      free (send_tgt_lnum);
      break;
    }

    send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[lComm]);
    for (int i = 0; i < send_shift[lComm]; i++) {
      send_g_num[i] = recv_g_num[send_tgt_lnum[i]];
    }

    /* send g_num buffer */
    recv_g_num = realloc (recv_g_num, sizeof(PDM_g_num_t) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                       recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                       octree->comm);


    /* send start_leaf buffer */
    start_leaf = realloc (start_leaf, sizeof(int) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_ngb_id, send_count, send_shift, PDM_MPI_INT,
                       start_leaf,  recv_count, recv_shift, PDM_MPI_INT,
                       octree->comm);
    free (send_ngb_id);


    /* send upper bound buffer */
    /* /!\ UPDATE upper_bound_dist before sending /!\ */
    double *send_upper_bound_dist = malloc (sizeof(double) * send_shift[lComm]);
    for (int i = 0; i < send_shift[lComm]; i++) {
      send_upper_bound_dist[i] = upper_bound_dist[send_tgt_lnum[i]];
    }
    upper_bound_dist = realloc (upper_bound_dist, sizeof(double) * n_recv_pts);
    PDM_MPI_Alltoallv (send_upper_bound_dist, send_count, send_shift, PDM_MPI_DOUBLE,
                       upper_bound_dist,      recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);
    free (send_upper_bound_dist);


    /* send coords buffer */
    send_coord = realloc (send_coord, sizeof(double) * send_shift[lComm] * dim);

    for (int i = 0; i < send_shift[lComm]; i++) {
      for (int j = 0; j < dim; j++)
        send_coord[dim*i + j] = recv_coord[dim*send_tgt_lnum[i] + j];
    }
    free (send_tgt_lnum);

    n_recv_pts = recv_shift[lComm];
    for (int i = 0; i < lComm; i++) {
      send_count[i]   *= dim;
      recv_count[i]   *= dim;
      send_shift[i+1] *= dim;
      recv_shift[i+1] *= dim;
    }
    recv_coord = realloc (recv_coord, sizeof(double) * recv_shift[lComm]);
    PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
                       recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
                       octree->comm);


    local_closest_src_gnum =
      realloc (local_closest_src_gnum, sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
    local_closest_src_dist =
      realloc (local_closest_src_dist, sizeof(double) * n_recv_pts * n_closest_points);

    PDM_hash_tab_purge (ht, PDM_FALSE);
  }
  PDM_hash_tab_free (ht);

  free (start_leaf);
  free (local_closest_src_gnum);
  free (local_closest_src_dist);

  free (send_coord);
  free (send_g_num);
  free (send_count);
  free (send_shift);

  free (recv_coord);
  free (recv_g_num);
  free (recv_count);
  free (recv_shift);

  _min_heap_free (merge_heap);
  free (tmp_closest_src_dist);
  free (tmp_closest_src_gnum);

  free (block_upper_bound_dist);

  free (processed_tgt_pts);




  /* Final Block-to-part exchanges */
  PDM_block_to_part_t *btp1 = PDM_block_to_part_create (block_distrib_idx1,
                                                        (const PDM_g_num_t **) &pts_g_num,
                                                        &n_pts,
                                                        1,
                                                        octree->comm);
  PDM_part_to_block_free (ptb1);

  PDM_block_to_part_exch (btp1,
                          sizeof(double),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_dist1,
                          NULL,
                          (void **) &closest_octree_pt_dist2);

  PDM_block_to_part_exch (btp1,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &_n_closest_points,
                          block_closest_src_gnum1,
                          NULL,
                          (void **) &closest_octree_pt_g_num);

  free (block_closest_src_dist1);
  free (block_closest_src_gnum1);

  btp1 = PDM_block_to_part_free (btp1);
}

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
)
{
  _octree_t *octree = _get_from_id (id);

  double t1 = octree->times_elapsed[END] - octree->times_elapsed[BEGIN];
  double t2 = octree->times_cpu[END] - octree->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (octree->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, octree->comm);

  int rank;
  PDM_MPI_Comm_rank (octree->comm, &rank);

  if (rank == 0) {

    PDM_printf( "PDM_para_octree timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "PDM_para_octree timer : build octree : total (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_TOTAL],
                t_cpu_max[BUILD_TOTAL]);
    PDM_printf( "PDM_para_octree timer : build octree : step order points (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_ORDER_POINTS],
                t_cpu_max[BUILD_ORDER_POINTS]);
    PDM_printf( "PDM_para_octree timer : build octree : step block partition (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_BLOCK_PARTITION],
                t_cpu_max[BUILD_BLOCK_PARTITION]);
    PDM_printf( "PDM_para_octree timer : build octree : step local nodes (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NODES],
                t_cpu_max[BUILD_LOCAL_NODES]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 1 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP1],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP1]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 2 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP2],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP2]);
    PDM_printf( "PDM_para_octree timer : build octree : step local neighbours - step 3 (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_LOCAL_NEIGHBOURS_STEP3],
                t_cpu_max[BUILD_LOCAL_NEIGHBOURS_STEP3]);
    PDM_printf( "PDM_para_octree timer : build octree : step distant neighbours (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[BUILD_DISTANT_NEIGHBOURS],
                t_cpu_max[BUILD_DISTANT_NEIGHBOURS]);

  }

}



#ifdef __cplusplus
}
#endif /* __cplusplus */
