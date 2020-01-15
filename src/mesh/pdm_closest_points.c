/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
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
#include "pdm_handles.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_closest_points.h"

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


#define NTIMER 2

/*============================================================================
 * Type definitions
 *============================================================================*/


/**
 * \enum _timer_step_t
 *
 */

typedef enum {

  BEGIN    = 0,
  END      = 1,

} _timer_step_t;


/**
 * \struct _tgt_point_cloud_t
 * \brief  Target point cloud structure
 *
 */

typedef struct {

  int           n_part;            /*!< Number of partition */
  int          *n_points;          /*!< Number of points of each partition */
  double      **coords;            /*!< Point coordinates points of each partition */
  PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */
  PDM_g_num_t **closest_src_gnum;  /*!< Global numbering of the n_closest source points
                                        for each point of each partition  */
  double      **closest_src_dist; /*!< Distance to the n_closest source points
                                        for each point of each partition  */

} _tgt_point_cloud_t;


/**
 * \struct _src_point_cloud_t
 * \brief  Src point cloud structure
 *
 */

typedef struct {

  int           n_part;            /*!< Number of partition */
  int          *n_points;          /*!< Number of points of each partition */
  double      **coords;            /*!< Point coordinates points of each partition */
  PDM_g_num_t **gnum;              /*!< Point global numbering of each partition */

} _src_point_cloud_t;


/**
 * \struct _PDM_closest_t
 * \brief  Closest points structure
 *
 */

typedef struct {

  PDM_MPI_Comm comm;  /*!< MPI communicator */

  int n_closest;  /*!< Number of closest source points to find for each
                    target point  */

  _src_point_cloud_t *src_cloud; /*!< Source point cloud */

  _tgt_point_cloud_t *tgt_cloud; /*!< Target point cloud */

  PDM_timer_t *timer; /*!< Timer */

  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */


} _PDM_closest_t;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_closest_pts   = NULL;

static int idebug = 0;

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

static _PDM_closest_t *
_get_from_id
(
 int  id
)
{
  _PDM_closest_t *closest = (_PDM_closest_t *) PDM_Handles_get (_closest_pts, id);

  if (closest == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_closest_points error : Bad identifier\n");
  }

  return closest;
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure to look for the closest points of a point cloud
 * (target cloud) in an other point cloud (source cloud)
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_closest      Number of closest source points to find for each
 *                              target point
 *
 * \return     Identifier
 *
 */

int
PDM_closest_points_create
(
 const PDM_MPI_Comm comm,
 const int          n_closest
)
{
  if (_closest_pts == NULL) {
    _closest_pts = PDM_Handles_create (4);
  }

  _PDM_closest_t *closest = (_PDM_closest_t *) malloc(sizeof(_PDM_closest_t));

  int id = PDM_Handles_store (_closest_pts, closest);

  closest->comm = comm;
  closest->n_closest = n_closest;
  closest->src_cloud = NULL;
  closest->tgt_cloud = NULL;

  closest->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    closest->times_elapsed[i] = 0.;
    closest->times_cpu[i] = 0.;
    closest->times_cpu_u[i] = 0.;
    closest->times_cpu_s[i] = 0.;
  }

  return id;
}

void
PDM_closest_points_create_cf
(
 const PDM_MPI_Fint comm,
 const int          n_closest,
 int *id
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  *id = PDM_closest_points_create(_comm, n_closest);
}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id                Identifier
 * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
 * \param [in]   n_part_cloud_tgt  Number of partitions of the target cloud
 *
 */

void
PDM_closest_points_n_part_cloud_set
(
 const int  id,
 const int  n_part_cloud_src,
 const int  n_part_cloud_tgt
)
{
  _PDM_closest_t *cls = _get_from_id (id);
  assert(cls->src_cloud == NULL);
  assert(cls->tgt_cloud == NULL);

  cls->src_cloud = malloc (sizeof(_src_point_cloud_t));
  cls->tgt_cloud = malloc (sizeof(_tgt_point_cloud_t));

  cls->src_cloud->n_part = n_part_cloud_src;
  cls->src_cloud->coords = malloc (sizeof(double *) * n_part_cloud_src);
  cls->src_cloud->gnum = malloc (sizeof(int *) * n_part_cloud_src);
  cls->src_cloud->n_points = malloc (sizeof(int) * n_part_cloud_src);

  cls->tgt_cloud->n_part = n_part_cloud_tgt;
  cls->tgt_cloud->coords = malloc (sizeof(double *) * n_part_cloud_tgt);
  cls->tgt_cloud->gnum = malloc (sizeof(int *) * n_part_cloud_tgt);
  cls->tgt_cloud->n_points = malloc (sizeof(int) * n_part_cloud_tgt);
  cls->tgt_cloud->closest_src_gnum = NULL;
  cls->tgt_cloud->closest_src_dist = NULL;
}


/**
 *
 * \brief Set the target point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_tgt_cloud_set
(
 const int          id,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
  _PDM_closest_t *cls = _get_from_id (id);
  assert(cls->tgt_cloud == NULL);
  cls->tgt_cloud->coords[i_part] = coords;
  cls->tgt_cloud->gnum[i_part] = gnum;
}


/**
 *
 * \brief Set the source point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_src_cloud_set
(
 const int          id,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
)
{
  _PDM_closest_t *cls = _get_from_id (id);
  assert(cls->src_cloud == NULL);
  cls->src_cloud->coords[i_part] = coords;
  cls->src_cloud->gnum[i_part] = gnum;
}

/**
 *
 * \brief Look for closest points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_closest_points_compute
(
 const int id
)
{
  //  _PDM_closest_t *cls = _get_from_id (id);
  //TODO: PDM_closest_points_compute algorithm

  /* int */
  /*   PDM_para_octree_create */
  /*   ( */
  /*    const int n_point_cloud, (n_part_source) */
  /*    const int depth_max, */
  /*    const int points_in_leaf_max, */
  /*    const PDM_MPI_Comm comm */
  /*    ); */



/* void */
/* PDM_para_octree_point_cloud_set */
/* ( */
/*  const int          id, */
/*  const int          i_point_cloud, */
/*  const int          n_points, */
/*  const double      *coords, */
/*  const PDM_g_num_t *g_num */
/* ); */

/* void */
/* PDM_para_octree_build */
/* ( */
/*  const int          id */
/* ); */

/* void */
/* PDM_para_octree_closest_point */
/* ( */
/* const int    id, */
/* const int    n_closest_points, */
/* const int    n_pts, */
/* double      *pts, */
/* PDM_g_num_t *pts_g_num, */
/* PDM_g_num_t *closest_octree_pt_g_num, */
/* double      *closest_octree_pt_dist2 */
/* ); */

/* void */
/* PDM_para_octree_free */
/* ( */
/*  const int          id */
/* ); */


}


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part_tgt            Index of partition of the cloud
 * \param [out]  closest_src_g_num     Global number of the closest element (size = n_closest * n_tgt_points)
 * \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
 *
 */

void
PDM_closest_points_get
(
 const int        id,
 const int        i_part_tgt,
 PDM_g_num_t    **closest_src_gnum,
       double   **closest_src_distance
)
{
  _PDM_closest_t *cls = _get_from_id (id);

  assert (cls->tgt_cloud->closest_src_gnum != NULL);
  assert (cls->tgt_cloud->closest_src_dist != NULL);

  *closest_src_gnum = cls->tgt_cloud->closest_src_gnum[i_part_tgt];
  *closest_src_distance = cls->tgt_cloud->closest_src_dist[i_part_tgt];
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
PDM_closest_points_free
(
 const int id,
 const int partial
)
{
  _PDM_closest_t *cls = _get_from_id (id);

  if (!partial) {
    if (cls->tgt_cloud->closest_src_gnum != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_gnum[j] != NULL) {
          free (cls->tgt_cloud->closest_src_gnum[j]);
        }
      }
      free (cls->tgt_cloud->closest_src_gnum);
    }

    if (cls->tgt_cloud->closest_src_dist != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_dist[j] != NULL) {
          free (cls->tgt_cloud->closest_src_dist[j]);
        }
      }
      free (cls->tgt_cloud->closest_src_dist);
    }
  }

  if (cls->tgt_cloud->gnum != NULL) {
    free (cls->tgt_cloud->gnum);
  }
  if (cls->tgt_cloud->coords != NULL) {
    free (cls->tgt_cloud->coords);
  }

  if (cls->src_cloud->gnum != NULL) {
    free (cls->src_cloud->gnum);
  }
  if (cls->src_cloud->coords != NULL) {
    free (cls->src_cloud->coords);
  }

  PDM_timer_free(cls->timer);

  free (cls);

  PDM_Handles_handle_free (_closest_pts, id, PDM_FALSE);

  const int n_closest_pts = PDM_Handles_n_get (_closest_pts);

  if (n_closest_pts == 0) {
    _closest_pts = PDM_Handles_free (_closest_pts);
  }
}


/**
 *
 * \brief  Dump elapsed and CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_closest_points_dump_times
(
 const int id
)
{
  _PDM_closest_t *cls = _get_from_id (id);
  double t1 = cls->times_elapsed[END] - cls->times_elapsed[BEGIN];
  double t2 = cls->times_cpu[END] - cls->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  int rank;
  PDM_MPI_Comm_rank (cls->comm, &rank);

  if (rank == 0) {

    PDM_printf( "closest_points timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
  }
}

#ifdef	__cplusplus
}
#endif
#undef NTIMER
