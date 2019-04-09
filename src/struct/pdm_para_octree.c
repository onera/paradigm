
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

  int  n_points;       /*!< Number of points in octant*/
  int  range;         /*!< Start index of point list for each octant */
  int  is_leaf; /*!< IS a leaf >*/
  
} _octant_t;


/**
 * \struct _octree_t
 * \brief  Define an octree
 * 
 */

typedef struct  {

  double  extents[6];            /*!< Extents of current process */ 
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

  for (int i = 0; i < 3; i++) {
    octree->extents[i]   = -HUGE_VAL;
    octree->extents[3+i] =  HUGE_VAL;
    octree->s[i] = 0.;
    octree->d[i] = 0.;
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
  octree->points = realloc (octree->points, octree->n_points * sizeof(double) * 3);
  octree->points_icloud = realloc (octree->points_icloud, octree->n_points * sizeof(int));
  octree->points_gnum = realloc (octree->points_gnum, octree->n_points * sizeof(PDM_g_num_t));
  octree->points_code = realloc (octree->points_code, octree->n_points * sizeof(_octant_t));

  for (int i = 0; i < 3 * n_points; i++) {
    octree->points[3*idx + i] = coords[i];
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

   /* PDM_morton_get_coord_extents(int               dim, */
   /*                              size_t            n_coords, */
   /*                              const double      coords[], */
   /*                              double             g_extents[], */
   /*                              PDM_MPI_Comm          comm); */

 
/* void */
/* PDM_morton_encode_coords(int                dim, */
/*                          PDM_morton_int_t   level, */
/*                          const double   extents[], */
/*                          size_t             n_coords, */
/*                          const double   coords[], */
/*                          PDM_morton_code_t  m_code[]); */

 
/* void */
/* PDM_morton_local_order(int                n_codes, */
/*                        const PDM_morton_code_t  morton_codes[], */
/*                        int                order[]); */


/* void */
/* PDM_morton_local_sort(int          n_codes, */
/*                       PDM_morton_code_t  morton_codes[]); */

 
 /* Normalize */

 /* Encode points */

 /* Sort points local */


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
