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
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_mesh_dist.h"
//#include "pdm_mesh_nodal.h"
#include "pdm_surf_mesh.h"
#include "pdm_handles.h"
//#include "pdm_octree.h"
//#include "pdm_dbbtree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
//#include "pdm_triangle.h"
//#include "pdm_polygon.h"
#include "pdm_timer.h"
//#include "pdm_hash_tab.h"
#include "pdm_gnum.h"

#include "pdm_wall_dist.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 4

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_Queue_t
 * \brief  Define a queue
 *
 */

typedef struct {

  int beg;       /*!< Beginning of the queue */
  int end;       /*!< End of the queue */

  int size;      /*!< Size of the queue */
  int n_elt;    /*!< Number of elements */

  int *queue;    /*!< queue array */

} _PDM_queue_t;


/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                         = 0,
  FIRST_THICKNESS               = 1,
  PROPAGATION                   = 2,
  END                           = 3,

} _ol_timer_step_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {
  PDM_g_num_t          g_n_cell;      /*!< Global number of cells */
  PDM_g_num_t          g_n_face;      /*!< Global number of faces */
  PDM_g_num_t          g_n_vtx;       /*!< Global number of vertices */
  int                  n_part;        /*!< Number of partitions */
  int                 *n_cell;        /*!< Number of cells (size = \n_part) */
  int                 *n_face;        /*!< Number of faces (size = \n_part) */
  int                 *n_vtx;         /*!< Number of vertices (size = \n_part) */
  const int          **cell_face_idx; /*!< Cell -> face connectivity index */
  const int          **cell_face;     /*!< Cell -> face connectivity */
  const double       **cell_center;   /*!< Cell center or NULL */
  const PDM_g_num_t  **cell_ln_to_gn; /*!< Cell local numbering to global numbering */
  const int          **face_vtx_idx;  /*!< Face -> vtx connectivity index */
  const int          **face_vtx;      /*!< Face -> vtx connectivity */
  const int          **face_cell;     /*!< Face -> cell connectivity */
  const PDM_g_num_t  **face_ln_to_gn; /*!< Face local numbering to global numbering */
  const double       **coords;        /*!< Vertex coordinates */
  const PDM_g_num_t  **vtx_ln_to_gn;  /*!< Vertex local numbering to global numbering */

  double      **closest_elt_distance;  /*!< Distance to the closest element */
  double      **closest_elt_projected; /*!< Projected point on the closest element */
  PDM_g_num_t **closest_elt_gnum;      /*!< Global numbering of the closest element */

} _PDM_vol_mesh_t;


/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  PDM_MPI_Comm comm;  /*!< MPI communicator */

  PDM_surf_mesh_t *surf_mesh;  /*!< Surface mesh pointer */

  _PDM_vol_mesh_t *vol_mesh; /*!< Volume mesh pointer */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


} _PDM_dist_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dists   = NULL;

static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return a new queue
 *
 * \param [in]   ini_size      Initial size of the queue
 *
 */

static _PDM_queue_t *
_PDM_queue_init
(
 const int ini_size
)
{
  _PDM_queue_t *q = malloc (sizeof(_PDM_queue_t));

  q->size = ini_size;

  q->queue = malloc(sizeof(int) * ini_size);
  q->beg = 0;
  q->end = -1;
  q->n_elt = 0;

  return q;

}


/**
 *
 * \brief Free a new queue
 *
 * \param [in]   q      Queue
 *
 */

static _PDM_queue_t *
_PDM_queue_free
(
 _PDM_queue_t *q
)
{
  free (q->queue);
  free (q);
  return NULL;
}


/**
 *
 * \brief Update size of the queue
 *
 * \param [in]   q      Queue
 *
 */

static void
_PDM_queue_update_size
(
 _PDM_queue_t *q,
 const int new_size
)
{
  assert (new_size > q->size);

  q->queue = realloc (q->queue, sizeof(int) * new_size);

  const int _beg = q->beg;
  const int _end = q->end;

  assert (_beg != _end);

  const int step = new_size - q->size;

  if (_beg > _end) {

    for (int i = q->size - 1; i >= _beg; i--) {
      q->queue[i+step] = q->queue[i];
    }
    q->beg += step;

  }

}


/**
 *
 * \brief Push value into the queue
 *
 * \param [in]   q      Queue
 * \param [in]   val    Value
 *
 */

static void
_PDM_queue_push
(
 _PDM_queue_t *q,
 const int val
)
{

  if (q->n_elt >= q->size) {
    _PDM_queue_update_size (q, 2*q->size);
  }

  q->end = (q->end + 1) % q->size;

  q->queue[q->end] = val;

  q->n_elt += 1;

}


/**
 *
 * \brief Pull the firts value in the queue
 *
 * \param [in]   q      Queue
 * \param [out]  val    Value
 *
 * \return 1 if the queue is empty, 0 otherwise
 */

static int
_PDM_queue_pull
(
 _PDM_queue_t *q,
 int  *val
)
{

  int is_empty = 0;

  if (q->n_elt <= 0) {
    is_empty = 1;
  }

  else {
    *val = q->queue[q->beg];
    q->beg = (q->beg + 1) % q->size;
    q->n_elt += -1;
  }

  return is_empty;
}


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
 * \brief Create a structure to compute distance of cell centers of a volume mesh
 * to a surface mesh
 *
 * \param [in]   comm           MPI communicator
 *
 * \return     Identifier
 *
 */

int
PDM_wall_dist_create
(
 const PDM_MPI_Comm comm
)
{
  if (_dists == NULL) {
    _dists = PDM_Handles_create (4);
  }

  _PDM_dist_t *dist = (_PDM_dist_t *) malloc(sizeof(_PDM_dist_t));

  int id = PDM_Handles_store (_dists, dist);

  dist->surf_mesh = NULL;
  dist->vol_mesh = NULL;

  dist->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    dist->times_elapsed[i] = 0.;
    dist->times_cpu[i] = 0.;
    dist->times_cpu_u[i] = 0.;
    dist->times_cpu_s[i] = 0.;
  }

  return id;
}

void
PDM_wall_dist_create_cf
(
 const PDM_MPI_Fint comm,
 int *id
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_wall_dist_create (_comm);
}


/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_wall_dist_surf_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  if (dist->surf_mesh != NULL) {
    dist->surf_mesh = PDM_surf_mesh_free (dist->surf_mesh);
  }

  dist->surf_mesh =
    PDM_surf_mesh_create (n_g_face, n_g_vtx, n_part, dist->comm);
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_wall_dist_surf_mesh_part_set
(
 const int          id,
 const int          i_part,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{

  _PDM_dist_t *dist = _get_from_id (id);

  PDM_surf_mesh_part_input (dist->surf_mesh,
                            i_part,
                            n_face,
                            face_vtx_idx,
                            face_vtx,
                            face_ln_to_gn,
                            n_vtx,
                            coords,
                            vtx_ln_to_gn);
}

/**
 *
 * \brief Set global data of a surface mesh
 *
 * \param [in]   id             Identifier
 * \param [in]   n_g_face       Global number of cells
 * \param [in]   n_g_face       Global number of faces
 * \param [in]   n_g_vtx        Global number of vertices
 * \param [in]   n_part         Number of partition
 *
 */

void
PDM_wall_dist_vol_mesh_global_data_set
(
 const int         id,
 const PDM_g_num_t n_g_cell,
 const PDM_g_num_t n_g_face,
 const PDM_g_num_t n_g_vtx,
 const int         n_part
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  assert (dist->vol_mesh == NULL);

  dist->vol_mesh = malloc (sizeof(_PDM_vol_mesh_t));

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  _vol_mesh->g_n_cell = n_g_cell;
  _vol_mesh->g_n_face = n_g_face;
  _vol_mesh->g_n_vtx = n_g_vtx;
  _vol_mesh->n_part = n_part;

  _vol_mesh->n_cell = malloc (sizeof(int) * n_part);
  _vol_mesh->n_face = malloc (sizeof(int) * n_part);
  _vol_mesh->n_vtx = malloc (sizeof(int) * n_part);

  _vol_mesh->cell_face_idx = malloc (sizeof(int *) * n_part);
  _vol_mesh->cell_face = malloc (sizeof(int *) * n_part);
  _vol_mesh->cell_center = malloc (sizeof(double *) * n_part);
  _vol_mesh->cell_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part);

  _vol_mesh->face_vtx_idx = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_vtx = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_cell = malloc (sizeof(int *) * n_part);
  _vol_mesh->face_ln_to_gn = malloc (sizeof (PDM_g_num_t *) * n_part);

  _vol_mesh->coords = malloc (sizeof(double *) * n_part);
  _vol_mesh->vtx_ln_to_gn = malloc (sizeof(PDM_g_num_t *) * n_part);
}


/**
 *
 * \brief Set a part of a surface mesh
 *
 * \param [in]   id            Identifier
 * \param [in]   i_part        Partition to define
 * \param [in]   n_cell        Number of faces
 * \param [in]   cell_face_idx Index in the face -> vertex connectivity
 * \param [in]   cell_face     face -> vertex connectivity
 * \param [in]   cell_center   Cell center or NULL
 * \param [in]   cell_ln_to_gn Local cell numbering to global face numbering
 * \param [in]   n_face        Number of faces
 * \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
 * \param [in]   face_vtx      face -> vertex connectivity
 * \param [in]   face_cell     face -> cell   connectivity
 * \param [in]   face_ln_to_gn Local face numbering to global face numbering
 * \param [in]   n_vtx         Number of vertices
 * \param [in]   coords        Coordinates
 * \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
 *
 */

void
PDM_wall_dist_vol_mesh_part_set
(
 const int          id,
 const int          i_part,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const double      *cell_center,
 const PDM_g_num_t *cell_ln_to_gn,
 const int          n_face,
 const int         *face_vtx_idx,
 const int         *face_vtx,
 const int         *face_cell,
 const PDM_g_num_t *face_ln_to_gn,
 const int          n_vtx,
 const double      *coords,
 const PDM_g_num_t *vtx_ln_to_gn
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  assert (dist->vol_mesh == NULL);

  dist->vol_mesh = malloc (sizeof(_PDM_vol_mesh_t));

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  _vol_mesh->n_cell[i_part] = n_cell;
  _vol_mesh->cell_face_idx[i_part] = cell_face_idx;
  _vol_mesh->cell_face[i_part] = cell_face;
  _vol_mesh->cell_center[i_part] = cell_center;
  _vol_mesh->cell_ln_to_gn[i_part] = cell_ln_to_gn;

  _vol_mesh->n_face[i_part] = n_face;
  _vol_mesh->face_vtx_idx[i_part] = face_vtx_idx;
  _vol_mesh->face_vtx[i_part] = face_vtx;
  _vol_mesh->face_cell[i_part] = face_cell;
  _vol_mesh->face_ln_to_gn[i_part] = face_ln_to_gn;

  _vol_mesh->n_vtx[i_part] = n_vtx;
  _vol_mesh->coords[i_part] = coords;
  _vol_mesh->vtx_ln_to_gn[i_part] = vtx_ln_to_gn;

}


/**
 *
 * \brief Compute distance
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_wall_dist_compute
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;
  PDM_surf_mesh_t *_surf_mesh = dist->surf_mesh;

  double b_t_elapsed;
  double b_t_cpu;
  double b_t_cpu_u;
  double b_t_cpu_s;

  double e_t_elapsed;
  double e_t_cpu;
  double e_t_cpu_u;
  double e_t_cpu_s;

  dist->times_elapsed[BEGIN] = PDM_timer_elapsed(dist->timer);
  dist->times_cpu[BEGIN]     = PDM_timer_cpu(dist->timer);
  dist->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(dist->timer);
  dist->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  /* First step : Look for boundary cells in the volume mesh */

  PDM_timer_hang_on(dist->timer);
  b_t_elapsed = PDM_timer_elapsed(dist->timer);
  b_t_cpu     = PDM_timer_cpu(dist->timer);
  b_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  b_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);
  PDM_timer_resume(dist->timer);

  int _t_n_face = 0;
  for (int i = 0; _vol_mesh->n_part; i++) {
    _t_n_face += _vol_mesh->n_face[i];
  }

  int *bound_cell_idx = malloc (sizeof(int) * (_vol_mesh->n_part + 1));
  bound_cell_idx[0] = 0;

  int n_bound_cell = 0;
  int *bound_cell = malloc (sizeof(int) * _t_n_face);
  PDM_g_num_t *bound_parent_gnum = malloc (sizeof(PDM_g_num_t) * _t_n_face);
  double * bound_cell_center = malloc (sizeof(double) * 3 * _t_n_face);

  for (int i = 0; _vol_mesh->n_part; i++) {

    bound_cell_idx[i+1] = bound_cell_idx[i];

    const int *_face_cell = _vol_mesh->face_cell[i];
    const int _n_face = _vol_mesh->n_face[i];
    const double *_cell_center = _vol_mesh->cell_center[i];
    const PDM_g_num_t *_gnum = _vol_mesh->cell_ln_to_gn[i];

    for (int j = 0; j < _n_face; j++) {
      if (_face_cell[2*j+1] == 0) {
        int icell = _face_cell[2*j] - 1;
        bound_cell[n_bound_cell] = icell;
        bound_parent_gnum[n_bound_cell] = _gnum[icell];
        bound_cell_idx[i+1] += 1;
        for (int k = 0; k < 3; k++) {
          bound_cell_center[3*n_bound_cell+k] = _cell_center[3*icell+k];
        }
        n_bound_cell += 1;
      }
    }

  }

  bound_cell = realloc (bound_cell, sizeof(int) * n_bound_cell);
  bound_cell_center = realloc (bound_cell_center, sizeof(double) * 3 * n_bound_cell);
  bound_parent_gnum = realloc (bound_parent_gnum, sizeof(PDM_g_num_t) * n_bound_cell);

  /* Second step : Compute distance to the surface mesh for the cell centers of boundary cells
                   Call PDM_mesh_dist */

  int id_bound_dist = PDM_mesh_dist_create (PDM_MESH_NATURE_SURFACE_MESH,
                                            1,
                                            dist->comm);

  PDM_mesh_dist_n_part_cloud_set (id_bound_dist, 0, 1);

  int id_gnum = PDM_gnum_create (3, 1, 0, 0.001, dist->comm);

  PDM_gnum_set_from_parents (id_gnum, 0,  n_bound_cell, bound_parent_gnum);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *bound_gnum = PDM_gnum_get (id_gnum, 0);

  PDM_gnum_free (id_gnum, 1);

  PDM_mesh_dist_cloud_set (id_bound_dist,
                           0,
                           0,
                           n_bound_cell,
                           bound_cell_center,
                           bound_gnum);

  PDM_mesh_dist_surf_mesh_map (id_bound_dist, dist->surf_mesh);

  PDM_mesh_dist_compute (id_bound_dist);

  double      *closest_elt_distance_bound = NULL;
  double      *closest_elt_projected_bound = NULL;
  PDM_g_num_t *closest_elt_gnum_bound = NULL;

  PDM_mesh_dist_get (id_bound_dist,
                     0,
                     0,
                     &closest_elt_distance_bound,
                     &closest_elt_projected_bound,
                     &closest_elt_gnum_bound);

  PDM_mesh_dist_free (id_bound_dist, 1);

  free (bound_gnum);
  free (bound_cell_center);
  free (bound_parent_gnum);

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[FIRST_THICKNESS] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[FIRST_THICKNESS]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[FIRST_THICKNESS]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[FIRST_THICKNESS]   += e_t_cpu_s - b_t_cpu_s;

  /* Third step : Get vertices about closest faces :
     - Build connectivity with coordinates
     - Part to block faces */

  b_t_elapsed = e_t_elapsed;
  b_t_cpu     = e_t_cpu;
  b_t_cpu_u   = e_t_cpu_u;
  b_t_cpu_s   = e_t_cpu_s;

  PDM_timer_resume(dist->timer);

  int n_part_sm = PDM_surf_mesh_n_part_get (_surf_mesh);

  int *n_face_sm = malloc (sizeof(int) * n_part_sm);

  const PDM_g_num_t **gnum_sm = malloc (sizeof(PDM_g_num_t *) * n_part_sm);

  for (int i = 0; i < n_part_sm; i++) {
    n_face_sm[i] = PDM_surf_mesh_part_n_face_get (_surf_mesh, i);
    gnum_sm[i]   = PDM_surf_mesh_part_face_g_num_get (_surf_mesh, i);
  }

  PDM_part_to_block_t *ptb =
    PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                              1.,
                              (PDM_g_num_t **) gnum_sm,
                              NULL,
                              n_face_sm,
                              n_part_sm,
                              dist->comm);

  /* const int * */
/* PDM_surf_mesh_part_face_vtx_get */
/* ( */
/*  PDM_surf_mesh_t      *mesh, */
/*  int                   iPart */
/* ); */


/* const int * */
/* PDM_surf_mesh_part_face_vtx_idx_get */
/* ( */
/*  PDM_surf_mesh_t      *mesh, */
/*  int                   iPart */
/* ); */

/* const double * */
/* PDM_surf_mesh_part_vtx_get */
/* ( */
/*  PDM_surf_mesh_t      *mesh, */
/*  int                   iPart */
/* ); */


/* PDM_part_to_block_exch */
/* ( */
/*  PDM_part_to_block_t       *ptb, */
/*  size_t                     s_data, */
/*  PDM_stride_t               t_stride, */
/*  int                        cst_stride, */
/*  int                      **part_stride, */
/*  void                     **part_data, */
/*  int                      **block_stride, */
/*  void                     **block_data */
/* ); */


  /* Fourth step : Get connectivity with coordinates +
        Compute distance to the surface mesh for the other centers from the distance
        of the cell centers of boundary cells  */

  for (int i = 0; i < _vol_mesh->n_part; i++) {

    int n_bound_cell_part = bound_cell_idx[i+1] - bound_cell_idx[i];
    double      *_closest_elt_distance = _vol_mesh->closest_elt_distance[i];
    double      *_closest_elt_projected = _vol_mesh->closest_elt_projected[i];
    PDM_g_num_t *_closest_elt_gnum = _vol_mesh->closest_elt_gnum[i];
    const int _n_cell = _vol_mesh->n_cell[i];

    int *in_queue = malloc(sizeof(int) * _n_cell);

    for (int j = 0; j < _n_cell; j++) {
      _closest_elt_distance[j] = HUGE_VAL;
      _closest_elt_gnum[j] = -1;
      in_queue[j] = 0;
      for (int k = 0; k < 3; k++) {
        _closest_elt_projected[3*j+k] = HUGE_VAL;
      }
    }

    _PDM_queue_t * q = _PDM_queue_init(2 * n_bound_cell_part);

    // Faux : il faut ajouter les voisins dans la pile

    for (int j = 0; j < n_bound_cell_part; j++) {
      int idx = bound_cell_idx[i] + j;
      int icell = bound_cell[idx];
      in_queue[icell] = 1;
      _PDM_queue_push (q, bound_cell[idx]);
      _closest_elt_distance[icell] = closest_elt_distance_bound[idx];
      _closest_elt_gnum[icell] = closest_elt_gnum_bound[idx];
      for (int k = 0; k < 3; k++) {
        _closest_elt_projected[3*icell+k] = closest_elt_projected_bound[3*idx+k];
      }
    }

    int curr_cell;
    while (!_PDM_queue_pull(q, &curr_cell)) {
      in_queue[curr_cell] = 0;

      //TODO : calcul

    }

    free (in_queue);
    _PDM_queue_free (q);
  }

  PDM_part_to_block_free (ptb);

  free (n_face_sm);
  free (gnum_sm);

  free (closest_elt_distance_bound);
  free (closest_elt_projected_bound);
  free (closest_elt_gnum_bound);
  free (bound_cell);
  free (bound_cell_idx);

  PDM_timer_hang_on(dist->timer);
  e_t_elapsed = PDM_timer_elapsed(dist->timer);
  e_t_cpu     = PDM_timer_cpu(dist->timer);
  e_t_cpu_u   = PDM_timer_cpu_user(dist->timer);
  e_t_cpu_s   = PDM_timer_cpu_sys(dist->timer);

  dist->times_elapsed[PROPAGATION] += e_t_elapsed - b_t_elapsed;
  dist->times_cpu[PROPAGATION]     += e_t_cpu - b_t_cpu;
  dist->times_cpu_u[PROPAGATION]   += e_t_cpu_u - b_t_cpu_u;
  dist->times_cpu_s[PROPAGATION]   += e_t_cpu_s - b_t_cpu_s;
}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part                Index of partition of the cloud
 * \param [out]  closest_elt_distance  Distance
 * \param [out]  closest_elt_projected Projected point coordinates
 * \param [out]  closest_elt_g_num     Global number of the closest element
 *
 */

void
PDM_wall_dist_get
(
 const int          id,
 const int          i_part,
       double      **closest_elt_distance,
       double      **closest_elt_projected,
       PDM_g_num_t **closest_elt_gnum
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  assert(_vol_mesh->closest_elt_distance != NULL);

  *closest_elt_distance = _vol_mesh->closest_elt_distance[i_part];
  *closest_elt_projected = _vol_mesh->closest_elt_projected[i_part];
  *closest_elt_gnum = _vol_mesh->closest_elt_gnum[i_part];
}


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_wall_dist_free
(
 const int id,
 const int partial
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  _PDM_vol_mesh_t *_vol_mesh = dist->vol_mesh;

  if (_vol_mesh->closest_elt_distance != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_distance[i] != NULL) {
          free (_vol_mesh->closest_elt_distance[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_distance);
  }

  if (_vol_mesh->closest_elt_projected != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_projected[i] != NULL) {
          free (_vol_mesh->closest_elt_projected[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_projected);
  }

  if (_vol_mesh->closest_elt_gnum != NULL) {
    if (!partial) {
      for (int i = 0; i < _vol_mesh->n_part; i++) {
        if (_vol_mesh->closest_elt_gnum[i] != NULL) {
          free (_vol_mesh->closest_elt_gnum[i]);
        }
      }
    }
    free (_vol_mesh->closest_elt_gnum);
  }

  if (_vol_mesh->n_cell != NULL) {
    free (_vol_mesh->n_cell);
  }

  if (_vol_mesh->n_face != NULL) {
    free (_vol_mesh->n_face);
  }

  if (_vol_mesh->n_vtx != NULL) {
    free (_vol_mesh->n_vtx);
  }

  if (_vol_mesh->cell_face_idx != NULL) {
    free (_vol_mesh->cell_face_idx);
  }

  if (_vol_mesh->cell_face != NULL) {
    free (_vol_mesh->cell_face);
  }

  if (_vol_mesh->cell_center != NULL) {
    free (_vol_mesh->cell_center);
  }

  if (_vol_mesh->cell_ln_to_gn != NULL) {
    free (_vol_mesh->cell_ln_to_gn);
  }

  if (_vol_mesh->face_vtx_idx != NULL) {
    free (_vol_mesh->face_vtx_idx);
  }

  if (_vol_mesh->face_vtx != NULL) {
    free (_vol_mesh->face_vtx);
  }

  if (_vol_mesh->face_cell != NULL) {
    free (_vol_mesh->face_cell);
  }

  if (_vol_mesh->face_ln_to_gn != NULL) {
    free (_vol_mesh->face_ln_to_gn);
  }

  if (_vol_mesh->vtx_ln_to_gn != NULL) {
    free (_vol_mesh->vtx_ln_to_gn);
  }

  if (_vol_mesh->coords != NULL) {
    free (_vol_mesh->coords);
  }

  free (_vol_mesh);

  if (dist->surf_mesh != NULL) {
    PDM_surf_mesh_free (dist->surf_mesh);
  }

  PDM_Handles_handle_free (_dists, id, PDM_FALSE);

  const int n_dists = PDM_Handles_n_get (_dists);

  if (n_dists == 0) {
    _dists = PDM_Handles_free (_dists);
  }

}


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_wall_dist_dump_times
(
 const int id
)
{
  _PDM_dist_t *dist = _get_from_id (id);

  double t1 = dist->times_elapsed[END] - dist->times_elapsed[BEGIN];
  double t2 = dist->times_cpu[END] - dist->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_elaps_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_elapsed, t_elaps_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  double t_cpu_max[NTIMER];
  PDM_MPI_Allreduce (dist->times_cpu, t_cpu_max, NTIMER, PDM_MPI_DOUBLE, PDM_MPI_MAX, dist->comm);

  int rank;
  PDM_MPI_Comm_rank (dist->comm, &rank);

  if (rank == 0) {

    PDM_printf( "Cell center distance timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
    PDM_printf( "Cell center distance timer : Distance of the cell centers of the first thickness (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[FIRST_THICKNESS],
                t_cpu_max[FIRST_THICKNESS]);
    PDM_printf( "Cell center distance timer : Distance of the other cell centers (elapsed and cpu) :"
                " %12.5es %12.5es\n",
                t_elaps_max[PROPAGATION],
                t_cpu_max[PROPAGATION]);
  }

}

#ifdef	__cplusplus
}
#endif
