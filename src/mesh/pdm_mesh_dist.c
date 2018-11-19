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

  int  n_part;
  int  **coords;
  PDM_g_num_t **gnum;
  double **dist;
  
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
    dist->points_cloud[i].coords = NULL;
    dist->points_cloud[i].gnum = NULL;
    dist->points_cloud[i].dist = NULL;
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
 *
 */

void
PDM_mesh_dist_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
 const double      *coords
)
{
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

void
PDM_mesh_dist_cloud_with_initial_set
(
 const int          id,
 const int          i_point_cloud,
 const int          i_part,
 const int          n_points,
 const double      *initial_dist,
 const double      *coords
)
{
}


/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

void
PDM_mesh_dist_normal_set
(
 const int          id,
 const int          i_part,
 const double      *normal
)
{
}

  

/**
 *
 * \brief Set normal surface mesh
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   normal          Normal
 *
 */

void
PDM_mesh_dist_center_set
(
 const int          id,
 const int          i_part,
 const double      *center
)
{
}


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
 * \param [in]   id              Identifier
 * \param [in]   i_point_cloud   Current cloud
 * \param [in]   i_part          Index of partition
 * \param [out]  dist            Distance
 * \param [out]  proj            Projected point coordinates
 * \param [out]  closest_part    Closest partition
 * \param [out]  closest_elt     Closest element
 *
 */

void
PDM_mesh_dist_get
(
 const int       id,
 const int       i_point_cloud,
       double  **dist,
       double  **proj,
       int     **closest_part,
       int     **closest_elt
)
{
}

/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id  Identifier
 *
 * \return     Identifier
 */

int
PDM_mesh_dist_free
(
 const int id
)
{
  return 0;
}

  
#ifdef	__cplusplus
}
#endif
