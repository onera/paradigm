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
#include "pdm_octree_seq.h"

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
  int  depth_max;                /*!< Maximum depth of internal octrees */
  int  points_in_leaf_max;       /*!< Maximum number of point in a leaf 
                                  *   of internal octrees */
  int   *n_points;               /*!< Number of points in each cloud */
  int  max_n_points;             /*!< Maximum number of points in each cloud */
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

static const double _default_eps = 1e-9;

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


/**
 *
 * \brief Search a point 
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static int
_intersect_extents 
(
const double *first_extents,
const double *second_extents
)
{
  int intersect = 0;
  
  for (int i = 0; i < 3; i++) {
    if ((first_extents[3*i] >= second_extents[3*i]) &&
        (first_extents[3*i] <= second_extents[3*i+3])) {
      intersect = 1;
      break;
    }
    else if ((first_extents[3*i+3] >= second_extents[3*i]) &&
            (first_extents[3*i+3] <= second_extents[3*i+3])) {
      intersect = 1;
      break;
    }
    else if ((first_extents[3*i] <= second_extents[3*i]) &&
             (first_extents[3*i+3] >= second_extents[3*i+3])) {
      intersect = 1;
      break;
    }
  }
  
  return intersect;
}


/**
 *
 * \brief Search a point in local partitions 
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static void
_search_local_couple 
(
int **local_couple, 
int  *n_couple, 
int  *s_couple,
const int  point_cloud, 
const int  point_idx, 
const double *point_coords,
const double *point_box,
const int search_cloud,
const int associated_octree_id,        
const int associated_octree_node_id,
const double *associated_coords,
const double *associated_char_length,
const double tolerance        
)
{  
  int node_id = associated_octree_node_id;
  int octree_id = associated_octree_id;
  const double *coords = associated_coords;
  const double *char_length = associated_char_length;        
  
  if (PDM_octree_seq_leaf_is (octree_id, node_id)) {
    
    int *points_clouds_id;
    int *point_indexes;
    int n_candidates = PDM_octree_seq_n_points_get (octree_id, node_id);
    PDM_octree_seq_points_get (octree_id, node_id, 
                               &points_clouds_id, &point_indexes);
    
    for (int i = 0; i < n_candidates; i++) {
      double dist2;
      double tol;

      if (point_box == NULL) {
        const double *coords_candidate = associated_coords + 3 * point_indexes[i];
        double char_length_candidate = associated_char_length[point_indexes[i]];
        dist2 = (coords_candidate[0] - point_coords[0]) * 
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) * 
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) * 
                (coords_candidate[2] - point_coords[2]);
        double tol1 = char_length_candidate * tolerance * 
                      char_length_candidate * tolerance;

        double tol2 = (point_box[0] - point_coords[0]) * (point_box[0] - point_coords[0]) * 
                      tolerance * tolerance;

        tol = PDM_MIN (tol1, tol2);
      }
      else {
        const double *coords_candidate = associated_coords + 3 * point_indexes[i];
        dist2 = (coords_candidate[0] - point_coords[0]) * 
                (coords_candidate[0] - point_coords[0]) +
                (coords_candidate[1] - point_coords[1]) * 
                (coords_candidate[1] - point_coords[1]) +
                (coords_candidate[2] - point_coords[2]) * 
                (coords_candidate[2] - point_coords[2]);
        tol = _default_eps * tolerance *
              _default_eps * tolerance;
      }

      if (dist2 <= tol) {
        if (*n_couple >= *s_couple) {
          if (*s_couple == 0) {
            *s_couple = 4;
          }
          else {
            *s_couple *= 2;
          }
          *local_couple = realloc(*local_couple, sizeof(int) * (*s_couple) * 4);
        }

        int _n_couple = *n_couple;
        int *_local_couple = *local_couple;
        _local_couple[4*_n_couple]   = search_cloud;
        _local_couple[4*_n_couple+1] = point_indexes[i];
        _local_couple[4*_n_couple+2] = point_cloud;
        _local_couple[4*_n_couple+3] = point_idx;

      }
    }
  }
  
  else {
    for (int i = 0; i < 8; i++) {
      const int node_child = 
            PDM_octree_seq_children_get (octree_id, node_id,
                                         (PDM_octree_seq_child_t) i);
      if (node_child != -1) {
        if (point_box != NULL) {
          if (_intersect_extents (PDM_octree_node_extents_get (octree_id, node_child),
                                  point_box)) {

            _search_local_couple (local_couple, n_couple, s_couple, point_cloud, 
                            point_idx, point_coords, point_box, search_cloud,
                            node_child, octree_id, coords, char_length, tolerance);
          }
        }
        else {
          double _extents[6] = {point_coords[0] - _default_eps,
                                point_coords[1] - _default_eps,
                                point_coords[2] - _default_eps,
                                point_coords[0] + _default_eps,
                                point_coords[1] + _default_eps,
                                point_coords[2] + _default_eps};
          
          if (_intersect_extents (PDM_octree_node_extents_get (octree_id, node_child),
                                  _extents)) {

            _search_local_couple (local_couple, n_couple, s_couple, point_cloud, 
                            point_idx, point_coords, point_box, search_cloud,
                            node_child, octree_id, coords, char_length, tolerance);
          }
        }
      }
    }    
  }  
}


/**
 *
 * \brief Search a point in distant partitions 
 *
 * \param [in]   ppartId        ppart identifier
 *
 */

static void
_search_distant_couple 
(
int          **distant_couple, 
int           *n_couple, 
int           *s_couple,
const int      point_proc, 
const int      point_cloud, 
const int      point_idx, 
const double  *point_coords,
const double  *point_box,
const int      associated_octree_id,        
const int      associated_octree_node_id,
const double **associated_coords,
const double **associated_char_length,
const double   tolerance        
)
{
  
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

  ppm->depth_max = 1000;
  ppm->points_in_leaf_max = 4;
  
  ppm->octree_id = PDM_octree_create (n_point_cloud, ppm->depth_max, 
                                      ppm->points_in_leaf_max, tolerance, comm);
  
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
  
  int n_proc;
  PDM_MPI_Comm_size(ppm->comm , &n_proc);
  
  int *local_couple = NULL;
  int n_local_couple = 0;
  int s_local_couple = 0;
  
  /*
   * Local fusion
   */
  
  double *point_box = NULL;
  double _point_box[6];
  
  if (ppm->char_length != NULL) {
    point_box = _point_box;
  }
  
  ppm->max_n_points = 0;
  for (int i = 0; i < ppm->n_point_clouds; i++) {
    ppm->max_n_points = PDM_MAX (ppm->max_n_points, ppm->n_points[i]);
  }  
  for (int i = 0; i < ppm->n_point_clouds; i++) {
    const int octree_seq_id = PDM_octree_seq_create (1, 
                                                     ppm->depth_max, 
                                                     ppm->points_in_leaf_max, 
                                                     ppm->tolerance);
    
    PDM_octree_seq_build (octree_seq_id);
    
    const int root_id = PDM_octree_seq_root_node_id_get (octree_seq_id);
    
    const double *_char_length = NULL;
    if (ppm->char_length != NULL) {
      _char_length = ppm->char_length[i];
    }

    for (int j = i + 1; j < ppm->n_point_clouds; j++) {
      for (int k = 0; k < ppm->n_points[j]; k++) {
        const double *_coord = ppm->point_clouds[j] + 3 * k;
        if (point_box != NULL) {
          double char_length_point = ppm->char_length[j][k];
          double tolerance = ppm->tolerance;
          point_box[0] = _coord[0] - tolerance * char_length_point;
          point_box[1] = _coord[1] - tolerance * char_length_point;       
          point_box[2] = _coord[2] - tolerance * char_length_point;       
          point_box[3] = _coord[0] + tolerance * char_length_point;
          point_box[4] = _coord[1] + tolerance * char_length_point;       
          point_box[5] = _coord[2] + tolerance * char_length_point;
        }

        _search_local_couple (&local_couple, &n_local_couple, &s_local_couple,
                        j, k, _coord, point_box, i, octree_seq_id,
                        root_id, ppm->point_clouds[i], _char_length,
                        ppm->tolerance);

      }        
    }  

    PDM_octree_seq_free (octree_seq_id);
  }
  
  /*
   * Distant fusion
   *   - Send/recv points between processes candidates
   *   - Search points candidates in the octree
   *   - Update distant table couple
   *   - Check if the umber of couple is coherent between other processes 
   */

  const double *extents_proc = PDM_octree_processes_extents_get (ppm->octree_id);

  int s_tmp_store = sizeof(int) * ppm->max_n_points;
  int n_tmp_store = 0;

  int *tmp_store = malloc (sizeof(int) * s_tmp_store * 3);
  
  for (int i_cloud = 0; i_cloud < ppm->n_point_clouds; i_cloud++) {
    int n_points = ppm->n_points[i_cloud]; 
    const double *_coord = ppm->point_clouds[i_cloud]; 
    for (int i = 0; i < n_points; i++) {
      const double *__coord = _coord + 3 * i; 
      double box[8];
      
      if (ppm->char_length != NULL) {
        box[0] = __coord[0] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[1] = __coord[1] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[2] = __coord[2] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[3] = __coord[0] + ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[4] = __coord[1] + ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[5] = __coord[2] + ppm->char_length[i_cloud][i] * ppm->tolerance;
      }
      else {
        box[0] = __coord[0] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[1] = __coord[1] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[2] = __coord[2] - ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[3] = __coord[0] + ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[4] = __coord[1] + ppm->char_length[i_cloud][i] * ppm->tolerance;
        box[5] = __coord[2] + ppm->char_length[i_cloud][i] * ppm->tolerance;
      }

      for (int k = 0; k < n_proc ; k++) {
        const double *_extents_proc = extents_proc + 6;
      
        if (_intersect_extents(box, _extents_proc)) {
          if (n_tmp_store >= s_tmp_store) {
            s_tmp_store *= 2;
            tmp_store = realloc (tmp_store, sizeof(int) * s_tmp_store * 3);
          }
          tmp_store[3*n_tmp_store]   = k;
          tmp_store[3*n_tmp_store+1] = i_cloud;
          tmp_store[3*n_tmp_store+2] = i;
          n_tmp_store += 1;
        }
      }
    }
  }
  
  int *val_send_n = malloc(sizeof(int)*n_proc);
  
  for (int i = 0; i < n_proc; i++) {
    val_send_n[i] = 0;
  }
  
  for (int i = 0; i < n_tmp_store; i++) {
    val_send_n[tmp_store[3*n_tmp_store]]++;
  }

  int *val_recv_n = malloc (sizeof(int)*n_proc);
  PDM_MPI_Alltoall (val_send_n, 1, PDM_MPI_INT, val_recv_n, 1, PDM_MPI_INT, ppm->comm);
  
  // Envoi des points + char length en option sur les autres procs (test bounding box)

  int *val_send_idx = malloc (sizeof(int)*n_proc);
  int *val_recv_idx = malloc (sizeof(int)*n_proc);

  int _stride = 3 * 8 + 4 + 4; /* Coords + icloud + ipoint */
  if (ppm->char_length != NULL) {
    _stride += 8; /* char_length */
  }
  
  for (int i = 0; i < n_proc; i++) {
    val_send_n[i] *= _stride;
  }

  val_send_idx[0] = 0;
  for (int i = 1; i < n_proc; i++) {
    val_send_idx[i] = val_send_idx[i-1] + val_send_idx[i];
    val_recv_idx[i] = val_recv_idx[i-1] + val_recv_idx[i];
  }

  unsigned char *val_send = 
        malloc (sizeof(unsigned char) * (val_send_idx[n_proc-1] + val_send_n[n_proc-1]));
  unsigned char *val_recv = 
        malloc (sizeof(unsigned char) * (val_recv_idx[n_proc-1] + val_recv_n[n_proc-1]));
  
  unsigned char *_tmp_val_send = val_send;
  for (int i = 0; i < n_tmp_store; i++) {
    size_t idx = 0;
    int iproc   = tmp_store[3*n_tmp_store];
    int i_cloud = tmp_store[3*n_tmp_store+1];
    int i_point = tmp_store[3*n_tmp_store+2];

    double *_coord = (double *) ppm->point_clouds[i_cloud] + 3 * i_point;
    double *_tmp_val_double = (double *) (_tmp_val_send + val_send_idx[iproc]);
    
    _tmp_val_double[0] = _coord[0]; 
    _tmp_val_double[1] = _coord[1]; 
    _tmp_val_double[2] = _coord[2]; 

    idx += 24;
    
    if (ppm->char_length != NULL) {
      double _char_length = ppm->char_length[i_cloud][i_point];
      _tmp_val_double[3] = _char_length;       
      idx += 32;
    }
    
    int *_tmp_val_int = (int *) (_tmp_val_send + idx);

    _tmp_val_int[0] = i_cloud; 
    _tmp_val_int[1] = i_point;

    val_send_idx[iproc] += idx + 4 + 4;
    
  }
    
  PDM_MPI_Alltoallv(val_send, val_send_n, val_send_idx, PDM_MPI_UNSIGNED_CHAR,
                    val_recv, val_recv_n, val_recv_idx, PDM_MPI_UNSIGNED_CHAR,
                    ppm->comm);
  
  /*
   * Build candidates_idx and candidates_desc arrays
   * from distant_couples and local_couples arrays
   * 
   */

  int *distant_couple  = NULL;
  int n_distant_couple = 0;
  int s_distant_couple = 0;
  
  double distant_coord[3];
  double point_box[6];
  
  unsigned char *_tmp_recv = val_recv;
  for (int i = 0; i < n_proc; i++) {
    int _deb = val_recv_idx[i] / _stride;
    int _end = _deb + val_recv_n[i] / _stride;
  
    for (int j = _deb; j < _end; j++) {
      double distant_coord[0] = (double) _tmp_recv[0];
      _tmp_recv += 8;
      double distant_coord[1] = (double) _tmp_recv[0];
      _tmp_recv += 8;
      double distant_coord[2] = (double) _tmp_recv[0];
      _tmp_recv += 8;
      double _char_length = -1;
      if (ppm->char_length != NULL) {
        _char_length  = (double) _tmp_recv[0];
        _tmp_recv += 8;       
        point_box[0] = 1;
        point_box[1] = 1;
        point_box[2] = 1;
        point_box[3] = 1;
        point_box[4] = 1;
        point_box[5] = 1;
      }
      else {
        point_box[0] = 1;
        point_box[1] = 1;
        point_box[2] = 1;
        point_box[3] = 1;
        point_box[4] = 1;
        point_box[5] = 1;
      }
      int distant_cloud  = _tmp_recv[0];
      _tmp_recv += 4;
      int distant_point = _tmp_recv[0];
      _tmp_recv += 4;

      _search_distant_couple (&distant_couple, 
                              &n_distant_couple, 
                              &s_distant_couple,
                              i,
                              distant_cloud, 
                              distant_point,
                              distant_coord,
                              point_box,
                              associated_octree_id,        
                              associated_octree_node_id,
                              associated_coords,
                              associated_char_length,
                               tolerance);      







    }
  }
  
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
