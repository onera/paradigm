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
#include "pdm_error.h"
#include "pdm_handles.h"
#include "pdm_mpi.h"
#include "pdm_octree.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_points_merge.h"

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

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/**
 * \struct _point_merge_t
 * \brief  Define a point merge structures
 * 
 */


typedef struct  {

  PDM_MPI_Comm comm;             /*!< MPI communicator */
  double tolerance;              /*!< Relative geometric tolerance */
  int   n_point_clouds;          /*!< Number of point cloud */
  int   *n_points;               /*!< Number of points in each cloud */
  const double **point_clouds;   /*!< points cloud */
  const double **char_length;    /*!< Characteristic length of points (optionnal) */
  int   octree_id;               /*!< Octree identifier */
  int   **candidates_idx;        /*!< Candidates indexes for each cloud */
  int   **candidates_desc;       /*!< Candidates description for each cloud */

} _point_merge_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_ppms   = NULL;

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

static _point_merge_t *
_get_from_id
(
 int  id
)
{
  _point_merge_t *ppm = (_point_merge_t *) PDM_Handles_get (_ppms, id);
    
  if (ppm == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_octree error : Bad identifier\n");
  }

  return ppm;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a points merge structure   
 *
 * \param [in]   n_point_cloud      Number of point cloud 
 * \param [in]   tolerance          Relative geometric tolerance
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_points_merge_create
(
 const int n_point_cloud,
 const double tolerance, 
 const PDM_MPI_Comm comm
)
{
  if (_ppms == NULL) {
    _ppms = PDM_Handles_create (4);
  }

  _point_merge_t *ppm = (_point_merge_t *) malloc(sizeof(_point_merge_t));

  int id = PDM_Handles_store (_ppms, ppm);
  
  ppm->comm = comm;             
  ppm->tolerance = tolerance;               
  ppm->n_point_clouds = n_point_cloud;          
  ppm->n_points = malloc (sizeof(int) * n_point_cloud);              
  ppm->point_clouds = malloc (sizeof(double *) * n_point_cloud);   
  ppm->char_length = malloc (sizeof(double *) * n_point_cloud);   
  ppm->octree_id = -1;               
  ppm->candidates_idx = malloc (sizeof(int *) * n_point_cloud);
  ppm->candidates_desc = malloc (sizeof(int *) * n_point_cloud);
  
  for (int i = 0; i < n_point_cloud; i++) {
    ppm->candidates_idx[i] = NULL;
    ppm->point_clouds[i] = NULL;
    ppm->char_length[i] = NULL;
    ppm->candidates_desc[i] = NULL;
  }

  const int depth_max = 1000;
  const int points_in_leaf_max = 4;
  
  
  ppm->octree_id = PDM_octree_create (n_point_cloud, depth_max, 
                                      points_in_leaf_max, tolerance, comm);
  
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
PDM_points_merge_free
(
 const int          id
)
{
  _point_merge_t *ppm = _get_from_id (id);

  for (int i = 0; i < ppm->n_point_clouds; i++) {
    if (ppm->candidates_idx[i] != NULL) {
      free (ppm->candidates_idx[i]);
    }
    if (ppm->candidates_desc[i] != NULL) {
      free (ppm->candidates_desc[i]);
    }
  }  
  
  free (ppm->candidates_idx);
  free (ppm->candidates_desc);
  free (ppm->point_clouds);
  free (ppm->char_length);
  free (ppm->n_points);
  
  PDM_octree_free (ppm->octree_id);
  
  free (ppm);
  
  PDM_Handles_handle_free (_ppms, id, PDM_FALSE);

  const int n_ppm = PDM_Handles_n_get (_ppms);
  
  if (n_ppm == 0) {
    _ppms = PDM_Handles_free (_ppms);
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
 * \param [in]   char_length        Characteristic length (or NULL)
 * 
 */

void
PDM_points_merge_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords, 
 const double      *char_length 
)
{
  _point_merge_t *ppm = _get_from_id (id);
  
  ppm->char_length[i_point_cloud] = char_length;
  ppm->point_clouds[i_point_cloud] = coords;
  ppm->n_points[i_point_cloud] = n_points;
  
  PDM_octree_point_cloud_set (ppm->octree_id, i_point_cloud, n_points, coords);
  
}


/**
 *
 * \brief Process merge points  
 *
 * \param [in]   id                 Identifier 
 * 
 */

void
PDM_points_merge_process
(
 const int          id
)
{
  _point_merge_t *ppm = _get_from_id (id);

  PDM_octree_build (ppm->octree_id);

  // Envoi des points + char length en option sur les autres procs (test bounding box)
  
  // Récupération des points distants
  
  // Boucle sur les points pour les localiser dans les octree (en option prendre en compte la taille charecteristic)
  
  // On renvoie le resultat quand il y a un candidate
  
  // Merge et tri candidats par point
  
 
}


/**
 *
 * \brief Get candidates to merge for each point 
 *
 * \param [in]   id             Identifier 
 * \param [in]   i_point_cloud  Current cloud 
 * \param [out]  candidates_idx Indexes of candidate for each current cloud point 
 *                              (size = number of points in the current cloud + 1) 
 * \param [out]  candidates_desc Candidates description (process, 
 *                                                       cloud in the process, 
 *                                                       point in the cloud)
 *
 */

void
PDM_points_merge_candidates_get
(
 const int     id,
 const int     i_point_cloud,
 const int    **candidates_idx, 
 const int    **candidates_desc 
) 
{
  _point_merge_t *ppm = _get_from_id (id);

  PDM_octree_build (ppm->octree_id);

}
