/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_mesh_dist.h"
#include "pdm_mesh_nodal.h"
#include "pdm_handles.h"


/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int           n_part;
  int          *n_points;
  double      **coords;
  PDM_g_num_t **gnum;
  double      **dist;
  double      **proj;
  int         **closest_elt_rank;
  int         **closest_elt_part;
  int         **closest_elt_lnum;
  PDM_g_num_t **closest_elt_gnum;
  
} _points_cloud_t;

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 * 
 */

typedef struct {

  int  n_point_cloud;
  PDM_MPI_Comm comm;
  int  mesh_nodal_id;

  _points_cloud_t *points_cloud;
  
} _PDM_dist_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dists   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/
  
/**
 *
 * \brief Return ppart object from it identifier
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static _PDM_dist_t *
_get_from_id
(
 int  id
)
{
  _PDM_dist_t *dist = (_PDM_dist_t *) PDM_Handles_get (_dists, id);
    
  if (dist == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_mesh_dist error : Bad identifier\n");
  }

  return dist;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to compute distance to a mesh nodal
 *
 * \param [in]   mesh_nodal_id  Mesh nodal identifier
 * \param [in]   n_point_cloud  Number of point cloud
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_create
(
 const int mesh_nodal_id,
 const int n_point_cloud,
 const PDM_MPI_Comm comm
)
{
  if (_dists == NULL) {
    _dists = PDM_Handles_create (4);
  }

  _PDM_dist_t *dist = (_PDM_dist_t *) malloc(sizeof(_PDM_dist_t));

  int id = PDM_Handles_store (_dists, dist);

  dist->mesh_nodal_id = mesh_nodal_id;
  dist->n_point_cloud = n_point_cloud;
  dist->comm = comm;
  dist->points_cloud = (_points_cloud_t*) malloc (sizeof(_points_cloud_t) * n_point_cloud);

  for (int i = 0; i <  n_point_cloud; i++) {
    dist->points_cloud[i].n_part = -1;
    dist->points_cloud[i].n_points = NULL;
    dist->points_cloud[i].coords = NULL;
    dist->points_cloud[i].gnum = NULL;
    dist->points_cloud[i].dist = NULL;
    dist->points_cloud[i].proj = NULL;
    dist->points_cloud[i].closest_elt_rank = NULL;
    dist->points_cloud[i].closest_elt_part = NULL;
    dist->points_cloud[i].closest_elt_lnum = NULL;
    dist->points_cloud[i].closest_elt_gnum = NULL;
  }
  
  return id;
}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   n_part          Number of partitions
 *
 */

void
PDM_mesh_dist_n_part_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_part
)
{
 _PDM_dist_t *dist = _get_from_id (id);

 dist->points_cloud[i_point_cloud].n_part = n_part;
 dist->points_cloud[i_point_cloud].n_points =
   realloc(dist->points_cloud[i_point_cloud].n_points, n_part * sizeof(int));
 dist->points_cloud[i_point_cloud].coords =
   realloc(dist->points_cloud[i_point_cloud].coords, n_part * sizeof(double *));
 dist->points_cloud[i_point_cloud].gnum =
   realloc(dist->points_cloud[i_point_cloud].gnum, n_part * sizeof(PDM_g_num_t * ));
 dist->points_cloud[i_point_cloud].dist =
   realloc(dist->points_cloud[i_point_cloud].dist, n_part * sizeof(double *));
 dist->points_cloud[i_point_cloud].proj =
   realloc(dist->points_cloud[i_point_cloud].proj, n_part * sizeof(double *));
 dist->points_cloud[i_point_cloud].closest_elt_rank =
   realloc(dist->points_cloud[i_point_cloud].closest_elt_rank, n_part * sizeof(int *));
 dist->points_cloud[i_point_cloud].closest_elt_part =
   realloc(dist->points_cloud[i_point_cloud].closest_elt_part, n_part * sizeof(int *));
 dist->points_cloud[i_point_cloud].closest_elt_lnum =
   realloc(dist->points_cloud[i_point_cloud].closest_elt_lnum, n_part * sizeof(int *));
 dist->points_cloud[i_point_cloud].closest_elt_gnum =
   realloc(dist->points_cloud[i_point_cloud].closest_elt_gnum, n_part * sizeof(PDM_g_num_t * ));

 for (int i = 0; i < n_part; i++) {
   dist->points_cloud[i_point_cloud].n_points[i] = -1;
   dist->points_cloud[i_point_cloud].coords[i] = NULL;
   dist->points_cloud[i_point_cloud].gnum[i] = NULL;
   dist->points_cloud[i_point_cloud].dist[i] = NULL;
   dist->points_cloud[i_point_cloud].proj[i] = NULL;
   dist->points_cloud[i_point_cloud].closest_elt_rank[i] = NULL;
   dist->points_cloud[i_point_cloud].closest_elt_part[i] = NULL;
   dist->points_cloud[i_point_cloud].closest_elt_lnum[i] = NULL;
   dist->points_cloud[i_point_cloud].closest_elt_gnum[i] = NULL;
 }
 
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number 
 *
 */

void
PDM_mesh_dist_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
 _PDM_dist_t *dist = _get_from_id (id);

 dist->points_cloud[i_point_cloud].n_points[i_part] = n_points;
 dist->points_cloud[i_point_cloud].coords[i_part] = coords;
 dist->points_cloud[i_point_cloud].gnum[i_part] = gnum;
}


/**
 *
 * \brief Set a point cloud with initial distance
 *
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Index of point cloud
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   initial_dist    Initial distance  
 * \param [in]   coords          Point coordinates
 *
 */

/* void */
/* PDM_mesh_dist_cloud_with_initial_set */
/* ( */
/*  const int          id, */
/*  const int          i_point_cloud, */
/*  const int          i_part, */
/*  const int          n_points, */
/*  const double      *initial_dist, */
/*  const double      *coords */
/* ) */
/* { */
/* } */


/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

/* void */
/* PDM_mesh_dist_normal_set */
/* ( */
/*  const int          id, */
/*  const int          i_part, */
/*  const double      *normal */
/* ) */
/* { */
/* } */

  

/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

/* void */
/* PDM_mesh_dist_center_set */
/* ( */
/*  const int          id, */
/*  const int          i_part, */
/*  const double      *center */
/* ) */
/* { */
/* } */


/**
 *
 * \brief Process merge points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_mesh_dist_process
(
 const int id
)
{

  /* 
   * Construction octree distribue avec les sommets de la surface 
   * (ou centre face a voir) 
   */

  /*
   *  Pour chaque point recherche du sommet le plus proche 
   *  (initialisation du bbtree) 
   */

  /*
   *  Construction du dbbtree 
   */
  
  /* 
   * Pour chaque point determination des boites situee 
   * a une plus courte distance que le maxima produit par l'octree 
   */

  /* 
   * Repartition des sommets suivant la numerotation absolue en fonction 
   * (poids sur le nombre de candidats)
   *  Necessite de faire une block_to_part avec poids
   * 
   * Il faut envoyer les coordonnees des sommets de chaque triangle ou 
   * quadrangle ou polygone
   *
   */

  /* 
   * Calcul des distances pour chaque candidat 
   */

  /* 
   * Envoi du resultat selon la repartition initiale des points  
   */
  
}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                Identifier
 * \param [in]   i_point_cloud     Current cloud
 * \param [in]   i_part            Index of partition of the cloud
 * \param [out]  distance          Distance
 * \param [out]  projected         Projected point coordinates
 * \param [out]  closest_elt_rank  Closest element rank
 * \param [out]  closest_elt_part  Closest element partition
 * \param [out]  closest_elt_l_num Local number of the closest element
 * \param [out]  closest_elt_g_num Global number of the closest element
 *
 */

void
PDM_mesh_dist_get
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
       double      **distance,
       double      **projected,
       int         **closest_elt_rank,
       int         **closest_elt_part,
       int         **closest_elt_lnum,
       PDM_g_num_t **closest_elt_gnum
)
{
 _PDM_dist_t *dist = _get_from_id (id);

 *distance = dist->points_cloud[i_point_cloud].dist[i_part];
 *projected = dist->points_cloud[i_point_cloud].proj[i_part];
 *closest_elt_rank = dist->points_cloud[i_point_cloud].closest_elt_rank[i_part];
 *closest_elt_part = dist->points_cloud[i_point_cloud].closest_elt_part[i_part];
 *closest_elt_lnum = dist->points_cloud[i_point_cloud].closest_elt_lnum[i_part];
 *closest_elt_gnum = dist->points_cloud[i_point_cloud].closest_elt_gnum[i_part];
}


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed. Otherwise, results are kept. 
 *
 */

void
PDM_mesh_dist_free
(
 const int id,
 const int partial
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  if (!partial) {
    for (int i_point_cloud = 0; i_point_cloud < dist->n_point_cloud; i_point_cloud++) {
      for (int i = 0; i < (dist->points_cloud[i_point_cloud]).n_part; i++) {
        free (dist->points_cloud[i_point_cloud].dist[i]);
        free (dist->points_cloud[i_point_cloud].proj[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_rank[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_part[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_lnum[i]);
        free (dist->points_cloud[i_point_cloud].closest_elt_gnum[i]); 
      }
    }
  }
  
  for (int i_point_cloud = 0; i_point_cloud < dist->n_point_cloud; i_point_cloud++) {
    free (dist->points_cloud[i_point_cloud].n_points);
    free (dist->points_cloud[i_point_cloud].coords);
    free (dist->points_cloud[i_point_cloud].gnum);
    free (dist->points_cloud[i_point_cloud].dist);
    free (dist->points_cloud[i_point_cloud].proj);
    free (dist->points_cloud[i_point_cloud].closest_elt_rank);
    free (dist->points_cloud[i_point_cloud].closest_elt_part);
    free (dist->points_cloud[i_point_cloud].closest_elt_lnum);
    free (dist->points_cloud[i_point_cloud].closest_elt_gnum);
  }
  
  free (dist->points_cloud);
  
  free (dist);
  
  PDM_Handles_handle_free (_dists, id, PDM_FALSE);
  
  const int n_dists = PDM_Handles_n_get (_dists);
  
  if (n_dists == 0) {
    _dists = PDM_Handles_free (_dists);
  }

}

  
#ifdef	__cplusplus
}
#endif
