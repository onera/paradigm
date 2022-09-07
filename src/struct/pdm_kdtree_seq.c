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
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_kdtree_seq.h"
#include "pdm_kdtree_seq_priv.h"

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

#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

static const double _eps_default = 1.e-12;
static const int    dbg_kdtree   = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * \brief   Evaluate a distribution array.
 *
 * \param [in]  n_ranges     Number of ranges in the distribution
 * \param [in]  distribution Number of elements associated to each range of the distribution
 * \param [in]  optim        Optimal count in each range
 *
 * \return  a fit associated to the distribution. If fit = 0, distribution is perfect.
 *
 */
static double
_evaluate_distribution(int          n_ranges,
                       int         *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}




static double
_median_point
(
 int      split_direction,
 double   extents_min,
 double   extents_max,
 int      point_range[2],
 int     *point_icloud,
 int     *point_ids,
 const double **pts_coord
 )
{
  double mid = 0.5*(extents_min + extents_max);

  int n_pts = point_range[1] - point_range[0];
  if (n_pts > 2) {
    double *x = malloc(sizeof(double) * n_pts);

    for (int j = 0; j < n_pts; j++) {
      int i = point_range[0] + j;

      x[j] = pts_coord[point_icloud[i]][point_ids[i]*3 + split_direction];
    }

    PDM_sort_double(x, NULL, n_pts);

    int h = n_pts / 2;
    // log_trace("n_pts = %d, h = %d\n", n_pts, h);
    mid = 0.5*(x[h] + x[h+1]);
    free(x);
  }

  return mid;
}



static void
_define_rank_distrib
(
 int      split_direction,
 int      point_range[2],
 int     *point_icloud,
 int     *point_ids,
 const double **pts_coord,
 int      n_sample,
 double  *sampling,
 double  *cfreq,
 int      distrib[2]
 )
{
  double i_npts = 1. / (double) (point_range[1] - point_range[0]);

  int l_distrib[n_sample];
  for (int i = 0; i < n_sample; i++) {
    l_distrib[i] = 0;
  }


  for (int i = point_range[0]; i < point_range[1]; i++) {
    double x = pts_coord[point_icloud[i]][point_ids[i]*3 + split_direction];

    int isample = PDM_binary_search_gap_double(x,
                                               sampling,
                                               n_sample+1);


    l_distrib[isample]++;
  }

  // PDM_log_trace_array_int(l_distrib, n_sample, "l_distrib : ");


  /* Define the cumulative frequency related to g_distribution */
  cfreq[0] = 0.;
  for (int id = 0; id < n_sample; id++) {
    cfreq[id+1] = cfreq[id] + l_distrib[id] * i_npts;
  }
  cfreq[n_sample] = 1.0;
  // PDM_log_trace_array_double(cfreq, n_sample+1, "cfreq : ");

  distrib[0] = 0;
  for (int id = 0; id < n_sample/2; id++) {
    distrib[0] += l_distrib[id];
  }
  distrib[1] = (point_range[1] - point_range[0]) - distrib[0];

  // PDM_log_trace_array_int(distrib, 2, "distrib : ");

}


static void
_update_sampling(int     n_sample,
                 double  c_freq[],
                 double *sampling[])
{
  int  i, j, next_id;
  double  target_freq, f_high, f_low, delta;
  double  s_low, s_high;

  // double new_sampling[n_sample+1];
  double *new_sampling = malloc(sizeof(double) * (n_sample+1));
  double *_sampling = *sampling;


  const double unit = 1/(double)n_sample;

  /* Compute new_sampling */
  new_sampling[0] = _sampling[0];
  next_id = 1;

  for (i = 0; i < n_sample; i++) {

    target_freq = (i+1)*unit;

    /* Find the next id such as c_freq[next_id] >= target_freq */

    for (j = next_id; j < n_sample + 1; j++) {
      if (c_freq[j] >= target_freq) {
        next_id = j;
        break;
      }
    }

    /* Find new s such as new_s is equal to target_freq by
       a linear interpolation */

    f_low = c_freq[next_id-1];
    f_high = c_freq[next_id];

    s_low = _sampling[next_id-1];
    s_high = _sampling[next_id];

    if (f_high - f_low > 0) {
      delta = (target_freq - f_low) * (s_high - s_low) / (f_high - f_low);
      new_sampling[i+1] = s_low + delta;
    }
    else /* f_high = f_low */
      new_sampling[i+1] = s_low + 0.5 * (s_low + s_high);

  } /* End of loop on samples */


  new_sampling[n_sample] = _sampling[n_sample];

  free(_sampling);

  /* Return pointers */
  *sampling = new_sampling;
}



static double
_approx_median_point
(
 const int n_sample,
 int      split_direction,
 double   extents_min,
 double   extents_max,
 int      point_range[2],
 int     *point_icloud,
 int     *point_ids,
 const double **pts_coord
 )
{
  double  fit, best_fit, optim;

  int n_pts = point_range[1] - point_range[0];

  // double sampling[n_sample+1];
  double *sampling = malloc(sizeof(double) * (n_sample+1));

   /* Define a naive sampling (uniform distribution) */
  double step = (extents_max - extents_min) / (double) n_sample;
  for (int i = 0; i <= n_sample; i++) {
    sampling[i] = extents_min + i*step;
  }
  sampling[n_sample] += 1e-3;

  // PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");


  int distrib[2];
  double cfreq[n_sample+1];
  _define_rank_distrib(split_direction,
                       point_range,
                       point_icloud,
                       point_ids,
                       pts_coord,
                       n_sample,
                       sampling,
                       cfreq,
                       distrib);

  optim = 0.5*n_pts;

  /* Initialize best choice */

  fit = _evaluate_distribution(2, distrib, optim);
  best_fit = fit;

  double best_sampling[n_sample+1];
  for (int i = 0; i < (n_sample + 1); i++) {
    best_sampling[i] = sampling[i];
  }

  /* Loop to get a better sampling array */

  // log_trace(">> loop\n");
  for (int n_iters = 0; (n_iters < 5 && fit > 0.10); n_iters++) {

    _update_sampling(n_sample, cfreq, &sampling);

    /* Compute the new distribution associated to the new sampling */

    _define_rank_distrib(split_direction,
                         point_range,
                         point_icloud,
                         point_ids,
                         pts_coord,
                         n_sample,
                         sampling,
                         cfreq,
                         distrib);

    fit = _evaluate_distribution(2, distrib, optim);
    // log_trace("n_iters = %d, fit = %f\n", n_iters, fit);
    // PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");

    /* Save the best sampling array and its fit */

    if (fit < best_fit) {

      best_fit = fit;
      for (int i = 0; i < (n_sample + 1); i++){
        best_sampling[i] = sampling[i];
      }
    }

  } /* End of while */
  free(sampling);

  // PDM_log_trace_array_double(best_sampling, n_sample+1, "best_sampling : ");


  double mid = best_sampling[n_sample/2];
  // log_trace("mid = %f\n", mid);

  return mid;
}




/**
 *
 * \brief Build a local kdtree's leaves.
 *
 * \param[in]  ancestor_id          Ancestor identifier
 * \param[in]  location_in_ancestor Location in ancestor
 * \param[in]  depth                Depth in the tree
 * \param[in]  extents              Extents associated with node:
 *                                  x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 * \param[in]  point_coords         Point coordinates
 * \param[in]  point_ids_tmp        Temporary point indexes
 * \param[in]  pos_tmp              Temporary point position in kdtree
 * \param[inout]  kdtree            Current kdtree structure
 * \param[inout]  point_range       Start and past-the end index in point_idx
 *                                  for current node (size: 2)
 */

static void
_build_kdtree_seq_leaves(const int                      ancestor_id,
                         // const PDM_kdtree_seq_child_t   location_in_ancestor,
                         const int                      split_direction,
                         const int                      depth,
                         const double                   extents[],
                         const double                 **point_coords,
                         int                           *point_icloud_tmp,
                         int                           *point_ids_tmp,
                         PDM_kdtree_seq_t              *kdtree,
                         int                            point_range[2])
{
  if (dbg_kdtree) {
    log_trace("\nnode_id = %d\n",
              kdtree->n_nodes);
    log_trace("ancestor_id = %d, split_direction = %d, depth = %d, point_range = %d/%d\n",
              ancestor_id, split_direction, depth, point_range[0], point_range[1]);
    log_trace("extents = %f %f %f  %f %f %f\n",
              extents[0], extents[1], extents[2], extents[3], extents[4], extents[5]);
  }
  _l_nodes_t *nodes = kdtree->nodes;

  int i, j, k, _n_nodes, _n_points, tmp_size;

  int count[2], idx[3], node_id[2];
  double mid, sub_extents[6];

  _n_nodes = kdtree->n_nodes;
  tmp_size = kdtree->n_nodes;

  /* Resize kdtree if necessary */

  if (kdtree->n_nodes >= kdtree->n_nodes_max) {
    if (kdtree->n_nodes == 0) {
      kdtree->n_nodes     = 1;
      kdtree->n_nodes_max = 8;//?
    }
    kdtree->n_nodes_max *= 2;

    nodes->ancestor_id          = realloc(nodes->ancestor_id,          sizeof(int                   ) * kdtree->n_nodes_max);
    nodes->is_leaf              = realloc(nodes->is_leaf,              sizeof(int                   ) * kdtree->n_nodes_max);
    // nodes->location_in_ancestor = realloc(nodes->location_in_ancestor, sizeof(PDM_kdtree_seq_child_t) * kdtree->n_nodes_max);
    nodes->depth                = realloc(nodes->depth,                sizeof(int                   ) * kdtree->n_nodes_max);
    nodes->children_id          = realloc(nodes->children_id,          sizeof(int                   ) * kdtree->n_nodes_max * 2);
    nodes->range                = realloc(nodes->range,                sizeof(int                   ) * kdtree->n_nodes_max * 2);
    nodes->idx                  = realloc(nodes->idx,                  sizeof(int                   ) * kdtree->n_nodes_max * 3);
    nodes->n_points             = realloc(nodes->n_points,             sizeof(int                   ) * kdtree->n_nodes_max);
    nodes->extents              = realloc(nodes->extents,              sizeof(double                ) * kdtree->n_nodes_max * 6);
  }


  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  for (j = 0; j < 2; j++) {
    count[j]   = 0;
    node_id[j] = -1;
  }

  if (depth < kdtree->depth_max && _n_points > kdtree->points_in_leaf_max) {

    /* Choose split direction */
    double max_range = -1.;
    int _split_direction = split_direction;

    if (_split_direction < 0) {
      for (int direction = 0; direction < 3; direction++) {
        double range = extents[3+direction] - extents[direction];
        if (range > max_range) {
          max_range        = range;
          _split_direction = direction;
        }
      }
    }

    /* Choose split point */
    if (dbg_kdtree) {
      log_trace("_split_direction = %d\n", _split_direction);
    }

    if (0) {
      mid = _median_point(_split_direction,
                          extents[_split_direction],
                          extents[3+_split_direction],
                          point_range,
                          kdtree->point_icloud,
                          kdtree->point_ids,
                          point_coords);

      // log_trace("true median = %f\n", mid);
    }
    else {
      mid = _approx_median_point(6,
                                 _split_direction,
                                 extents[_split_direction],
                                 extents[3+_split_direction],
                                 point_range,
                                 kdtree->point_icloud,
                                 kdtree->point_ids,
                                 point_coords);
    }
    if (dbg_kdtree) {
      log_trace("_split_direction = %d, mid = %f\n",
                _split_direction, mid);
    }


    /* Count points in each child node */

    for (i = point_range[0]; i < point_range[1]; i++) {

      k = 0;
      if (point_coords[kdtree->point_icloud[i]][kdtree->point_ids[i]*3 + _split_direction] > mid) {
        k = 1;
      }

      count[k] += 1;
    }

    if (dbg_kdtree) {
      log_trace("count = %d %d\n", count[0], count[1]);
    }


    /* Build index */

    idx[0] = 0;
    for (j = 0; j < 2; j++) {
      idx[j+1] = idx[j] + count[j];
    }

    for (j = 0; j < 2; j++) {
      count[j] = 0;
    }

    for (i = point_range[0], j = 0; i < point_range[1]; i++) {

      k = 0;
      if (point_coords[kdtree->point_icloud[i]][kdtree->point_ids[i]*3 + _split_direction] > mid) {
        k = 1;
      }

      point_icloud_tmp[idx[k] + count[k]] = kdtree->point_icloud[i];
      point_ids_tmp   [idx[k] + count[k]] = kdtree->point_ids[i];
      count[k] += 1;
    }


    /* Check if this subdivision is static
       and check coordinates to find multi point */

    for (i = point_range[0], j = 0; i < point_range[1]; i++, j++) {
      kdtree->point_icloud[i] = point_icloud_tmp[j];
      kdtree->point_ids[i]    = point_ids_tmp[j];
    }

    for (i = 0; i < 3; i++) {
      idx[i] = point_range[0] + idx[i];
    }

    if (dbg_kdtree) {
      log_trace("idx = %d %d %d\n", idx[0], idx[1], idx[2]);
    }

    /* Build leaves recursively */

    for (i = 0; i < 2; i++) {

      if ((idx[i+1] - idx[i]) > 0) {

        tmp_size++;

        node_id[i] = tmp_size;

        memcpy(sub_extents, extents, sizeof(double) * 6);
        if (i == 0) {
          sub_extents[3+_split_direction] = mid;
        }
        else {
          sub_extents[_split_direction]   = mid;
        }

        if (1) {
          // Tighter extents (fix contained points)
          for (int l = 0; l < 3; l++) {
            // if (l == _split_direction) continue;
            sub_extents[l  ] =  HUGE_VAL;
            sub_extents[l+3] = -HUGE_VAL;
          }
          for (int ipt = idx[i]; ipt < idx[i+1]; ipt++) {
            for (int l = 0; l < 3; l++) {
              // if (l == _split_direction) continue;
              double x = point_coords[kdtree->point_icloud[ipt]][kdtree->point_ids[ipt]*3 + l];
              sub_extents[l  ] = PDM_MIN(sub_extents[l  ], x);
              sub_extents[l+3] = PDM_MAX(sub_extents[l+3], x);
            }
          }
          for (int l = 0; l < 3; l++) {
            if (sub_extents[l+3] < sub_extents[l] + _eps_default) {
              sub_extents[l  ] -= 0.5*_eps_default;
              sub_extents[l+3] += 0.5*_eps_default;
            }
          }
        }

        if (dbg_kdtree) {
          log_trace("child %d, sub_extents = %f %f %f  %f %f %f\n",
                    i,
                    sub_extents[0], sub_extents[1], sub_extents[2], sub_extents[3], sub_extents[4], sub_extents[5]);
        }

        /* 1.0e-12 term in assert() used to allow for
           truncation error in for min = max case */
        assert(sub_extents[_split_direction] < sub_extents[_split_direction] + 1.0e-14);

        kdtree->n_nodes = tmp_size;

        int next_split_direction = split_direction;
        if (next_split_direction >= 0) {
          next_split_direction = (next_split_direction+1) % 3;
        }
        _build_kdtree_seq_leaves(_n_nodes,
                                  // (PDM_kdtree_seq_child_t) i,
                                 next_split_direction,
                                 depth+1,
                                 sub_extents,
                                 point_coords,
                                 point_icloud_tmp,
                                 point_ids_tmp,
                                 kdtree,
                                 idx + i);

        tmp_size = kdtree->n_nodes;
      }

    }

  }

  /* Finalize node */



  for (i = 0; i < 2; i++) {
    nodes->range[2*_n_nodes + i] = point_range[i];
  }

  for (i = 0; i < 3; i++) {
    nodes->idx[3*_n_nodes + i] = idx[i];
  }

  for (i = 0; i < 6; i++) {
    nodes->extents[6*_n_nodes + i] = extents[i];
  }

  for (i = 0; i < 2; i++) {
    nodes->children_id[2*_n_nodes + i] = node_id[i];
  }

  nodes->is_leaf[_n_nodes] =
  (nodes->children_id[2*_n_nodes + 0] == -1) &&
  (nodes->children_id[2*_n_nodes + 1] == -1);

  nodes->ancestor_id[_n_nodes] = ancestor_id;
  nodes->depth[_n_nodes]       = depth;

  nodes->n_points[_n_nodes] = _n_points;
  // nodes->location_in_ancestor[_n_nodes] = location_in_ancestor;
}




/**
 *
 * \brief   Compute extents of a point set
 *
 *  \param [in] dim         Space dimension of points to locate_3d
 *  \param [in] n_points    Number of points to locate
 *  \param [in] point_index optional indirection array to point_coords
 *                          (1 to n_points numbering)
 *  \param [in] point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *
 */

static void
_point_extents(const int     dim,
               const int     n_points,
               const int     point_index[],
               const double  point_coords[],
               double        extents[])
{
  int i;
  int j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }

    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }
}

/**
 *
 * \brief Build a kdtree
 *
 * \param[in]  kdtree    Current kdtree
 * .
 */

static void
_build_kdtree
(
 PDM_kdtree_seq_t *kdtree
)
{
  int point_range[2];

  /* Initialization */

  kdtree->n_nodes = 0;
  kdtree->n_nodes_max = 0;

  kdtree->nodes = malloc(sizeof(_l_nodes_t));
  kdtree->nodes->ancestor_id          = NULL;
  kdtree->nodes->is_leaf              = NULL;
  // kdtree->nodes->location_in_ancestor = NULL;
  kdtree->nodes->depth                = NULL;
  kdtree->nodes->children_id          = NULL;
  kdtree->nodes->range                = NULL;
  kdtree->nodes->idx                  = NULL;
  kdtree->nodes->n_points             = NULL;
  kdtree->nodes->extents              = NULL;


  for (int i = 0; i < kdtree->n_point_clouds; i++) {
    kdtree->t_n_points += kdtree->n_points[i];
  };

  kdtree->point_ids    = malloc (sizeof(int) * kdtree->t_n_points);
  kdtree->point_icloud = malloc (sizeof(int) * kdtree->t_n_points);

  int cpt = 0;
  for (int i = 0; i < kdtree->n_point_clouds; i++) {

    int n_points = kdtree->n_points[i];
    double extents[6];

    if (n_points > 0) {

      _point_extents(3,
                     n_points,
                     NULL,
                     kdtree->point_clouds[i],
                     extents);

      for (int i1 = 0; i1 < 3; i1++) {
        kdtree->extents[i1] = PDM_MIN (extents[i1], kdtree->extents[i1]);
        kdtree->extents[i1 + 3] = PDM_MAX (extents[i1 + 3], kdtree->extents[i1 + 3]);
      }

      for (int j = 0; j < n_points; j++) {
        kdtree->point_ids[cpt] = j;
        kdtree->point_icloud[cpt] = i;
        cpt +=1;
      }
    }
  }

  double delta = -1;
  for (int i = 0; i < 3; i++) {
    delta = PDM_MAX (kdtree->tolerance * (kdtree->extents[i + 3] - kdtree->extents[i]),
                     delta);
  }
  delta = PDM_MAX (delta,_eps_default);

  for (int i = 0; i < 3; i++) {
    kdtree->extents[i  ] += -1.01*delta;
    kdtree->extents[i+3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = kdtree->t_n_points;


  /* Build kd-tree recursively */
  int *point_ids_tmp    = malloc (sizeof(int) * kdtree->t_n_points);
  int *point_icloud_tmp = malloc (sizeof(int) * kdtree->t_n_points);

  int split_direction = -1;

  if (0) {
    /* Split dimension with greatest range */
    double max_range = -1.;
    for (int i = 0; i < 3; i++) {
      double range = kdtree->extents[3+i] - kdtree->extents[i];
      if (range > max_range) {
        max_range       = range;
        split_direction = i;
      }
    }
  }

  if (dbg_kdtree) {
    log_trace(">> _build_kdtree_seq_leaves\n");
  }
  _build_kdtree_seq_leaves(-1,
                           // location_in_ancestor,
                           split_direction,
                           -1,
                           kdtree->extents,
         (const double **) kdtree->point_clouds,
                           point_icloud_tmp,
                           point_ids_tmp,
                           kdtree,
                           point_range);

  kdtree->n_nodes += 1;

  if (dbg_kdtree) {
    // Dump kd-tree
    _l_nodes_t *nodes = kdtree->nodes;

    for (int i = 0; i < kdtree->n_nodes; i++) {
      if (1) {//nodes->is_leaf[i]) {
        log_trace("\nNode %d :", i);
        log_trace("  is_leaf = %d\n", nodes->is_leaf[i]);
        log_trace("  depth = %d\n", nodes->depth[i]);
        log_trace("  extents = %f %f %f  %f %f %f\n",
                  nodes->extents[6*i+0],
                  nodes->extents[6*i+1],
                  nodes->extents[6*i+2],
                  nodes->extents[6*i+3],
                  nodes->extents[6*i+4],
                  nodes->extents[6*i+5]);
        log_trace("  point_range = %d / %d\n", nodes->range[2*i+0], nodes->range[2*i+1]);
      }
    }
  }


  free (point_ids_tmp);
  free (point_icloud_tmp);
}


static void
_l_nodes_free
(
 PDM_kdtree_seq_t *kdtree
 )
{
  if (kdtree->nodes != NULL) {

    free(kdtree->nodes->ancestor_id);
    free(kdtree->nodes->is_leaf);
    // free(kdtree->nodes->location_in_ancestor);
    free(kdtree->nodes->depth);
    free(kdtree->nodes->children_id);
    free(kdtree->nodes->range);
    free(kdtree->nodes->idx);
    free(kdtree->nodes->n_points);
    free(kdtree->nodes->extents);

    free(kdtree->nodes);

    kdtree->nodes = NULL;
  }
}



inline static int
_box_dist2_min
(
 const int              dim,
 const double          *restrict extents,
 const double          *restrict coords,
 double                *restrict min_dist2
 )
{

  int inbox = 0;
  *min_dist2 = 0.;

  for (int i = 0; i < dim; i++) {
    if (coords[i] > extents[i+dim]) {
      double _min_dist2 = coords[i] - extents[dim+i];
      *min_dist2 += _min_dist2 * _min_dist2;
    }

    else if (coords[i] < extents[i]) {
      double _min_dist2 = coords[i] - extents[i];
      *min_dist2 += _min_dist2 * _min_dist2;
    }

    else {
      inbox += 1;
    }
  }

  return inbox == dim;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a kdtree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_kdtree_seq object
 */

PDM_kdtree_seq_t *
PDM_kdtree_seq_create
(
 const int    n_point_cloud,
 const int    depth_max,
 const int    points_in_leaf_max,
 const double tolerance
)
{
  PDM_kdtree_seq_t *kdtree = (PDM_kdtree_seq_t *) malloc(sizeof(PDM_kdtree_seq_t));

  kdtree->n_point_clouds     = n_point_cloud;
  kdtree->depth_max          = depth_max;
  kdtree->points_in_leaf_max = points_in_leaf_max;
  kdtree->tolerance          = tolerance;

  kdtree->n_nodes     = 0;
  kdtree->n_nodes_max = 0;

  kdtree->n_points     = malloc (sizeof(double  ) * n_point_cloud);
  kdtree->point_clouds = malloc (sizeof(double *) * n_point_cloud);
  for (int i = 0; i < n_point_cloud; i++) {
    kdtree->n_points[i]     = 0;
    kdtree->point_clouds[i] = NULL;
  }

  kdtree->point_icloud = NULL;
  kdtree->point_ids    = NULL;
  kdtree->nodes        = NULL;
  kdtree->t_n_points = 0;
  for (int i = 0; i < 3; i++) {
    kdtree->extents[i]     =  HUGE_VAL;
    kdtree->extents[i + 3] = -HUGE_VAL;
  }

  return kdtree;
}


/**
 *
 * \brief Free a kdtree structure
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_free
(
 PDM_kdtree_seq_t *kdtree
)
{
  free (kdtree->n_points);
  free (kdtree->point_clouds);
  free (kdtree->point_ids);
  free (kdtree->point_icloud);

  _l_nodes_free(kdtree);

  free (kdtree);
}



/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 *
 */


void
PDM_kdtree_seq_point_cloud_set
(
 PDM_kdtree_seq_t *kdtree,
 const int         i_point_cloud,
 const int         n_points,
 const double     *coords
)
{
  kdtree->n_points    [i_point_cloud] = n_points;
  kdtree->point_clouds[i_point_cloud] = coords;
}



/**
 *
 * \brief Build kdtree
 *
 * \param [in]   kdtree             Pointer to \ref PDM_kdtree_seq object
 *
 */

void
PDM_kdtree_seq_build
(
 PDM_kdtree_seq_t *kdtree
)
{
  if (kdtree->nodes == NULL) {
    _build_kdtree(kdtree);
  }

}



/**
 *
 * \brief Write kdtree nodes in a VTK file
 *
 * \param [in]   kdtree                 Pointer to \ref PDM_kdtree_seq object
 * \param [in]   filename               Output file name
 *
 */

void PDM_kdtree_seq_write_nodes
(
 PDM_kdtree_seq_t *kdtree,
 const char       *filename
 )
{
  _l_nodes_t *nodes = kdtree->nodes;

  double tol_visu = 1e-3;
  double _ext[6];

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "kdtree_seq\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*kdtree->n_nodes);
  for (int inode = 0; inode < kdtree->n_nodes; inode++) {
    double *ext = nodes->extents + 6*inode;

    if (1 && (kdtree->nodes->range[2*inode+1] - kdtree->nodes->range[2*inode] == 1)) {
      // Trick to visualize nodes with degenerate extents (single point)
      ext = _ext;
      for (int i = 0; i < 3; i++) {
        double x = nodes->extents[6*inode + i];
        double eps = tol_visu*(kdtree->extents[i+3] - kdtree->extents[i]);
        ext[i  ] = x - eps;
        ext[i+3] = x + eps;
      }
    }

    for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 2; i++) {
          int ii = (1-j)*i + j*(1-i);
          fprintf(f, "%f %f %f\n", ext[3*ii], ext[3*j+1], ext[3*k+2]);
        }
      }
    }
  }

  fprintf(f, "CELLS %d %d\n", kdtree->n_nodes, 9*kdtree->n_nodes);
  for (int inode = 0; inode < kdtree->n_nodes; inode++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*inode+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", kdtree->n_nodes);
  for (int i = 0; i < kdtree->n_nodes; i++) {
    fprintf(f, "%d\n", 12);
  }

  fprintf(f, "CELL_DATA %d\n", kdtree->n_nodes);

  fprintf(f, "FIELD node_field 3\n");
  fprintf(f, "depth 1 %d int\n", kdtree->n_nodes);
  for (int i = 0; i < kdtree->n_nodes; i++) {
    fprintf(f, "%d\n", kdtree->nodes->depth[i]);
  }
  fprintf(f, "is_leaf 1 %d int\n", kdtree->n_nodes);
  for (int i = 0; i < kdtree->n_nodes; i++) {
    fprintf(f, "%d\n", kdtree->nodes->is_leaf[i]);
  }
  fprintf(f, "n_pts 1 %d int\n", kdtree->n_nodes);
  for (int i = 0; i < kdtree->n_nodes; i++) {
    fprintf(f, "%d\n", kdtree->nodes->range[2*i+1] - kdtree->nodes->range[2*i]);
  }

  fclose(f);
}





/**
 *
 * \brief Look for points inside at set of balls
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  n_ball               Number of balls
 * \param [in]  ball_center          Center of balls (size = \ref n_ball * 3)
 * \param [in]  ball_radius2         Squared radius of balls (size = \ref n_ball)
 * \param [out] ball_pts_idx         Index for ball->points graph (size \ref n_ball + 1)
 * \param [out] ball_pts_l_num       Ball->points graph (cloud_id, point_id)
 * \param [out] ball_pts_dist2       Distance from points to ball centers
 *
 */

void
PDM_kdtree_seq_points_inside_ball
(
 const PDM_kdtree_seq_t  *kdtree,
 const int                n_ball,
 double                  *ball_center,
 double                  *ball_radius2,
 int                    **ball_pts_idx,
 int                    **ball_pts_l_num,
 double                 **ball_pts_dist2
 )
{
  const int n_children = 2;

  int s_pt_stack = ((n_children - 1) * (kdtree->depth_max - 1) + n_children);


  *ball_pts_idx = malloc(sizeof(int) * (n_ball + 1));
  int *pib_idx = *ball_pts_idx;
  pib_idx[0] = 0;

  int s_pib = 4*n_ball;
  *ball_pts_l_num = malloc(sizeof(int   ) * s_pib * 2);
  *ball_pts_dist2 = malloc(sizeof(double) * s_pib);

  int    *pib_l_num = *ball_pts_l_num;
  double *pib_dist2 = *ball_pts_dist2;


  _l_nodes_t *nodes = kdtree->nodes;


  int *stack = malloc(sizeof(int) * s_pt_stack);


  for (int iball = 0; iball < n_ball; iball++) {

    pib_idx[iball+1] = pib_idx[iball];

    double *_center  = ball_center + 3*iball;
    double  _radius2 = ball_radius2[iball];


    /* Start by root */
    int pos_stack = 0;
    double min_dist2;
    int inside_box = _box_dist2_min(3,
                                    &nodes->extents[0],
                                    _center,
                                    &min_dist2);

    if (inside_box || min_dist2 <= _radius2) {
      stack[pos_stack++] = 0;
    }


    while (pos_stack > 0) {

      int node_id = stack[--pos_stack];

      if (nodes->is_leaf[node_id]) {
        /* Leaf node */

        int *point_clouds_id = kdtree->point_icloud + nodes->range[2*node_id];
        int *point_indexes   = kdtree->point_ids    + nodes->range[2*node_id];

        for (int i = 0; i < nodes->n_points[node_id]; i++) {
          const double *_pt = kdtree->point_clouds[point_clouds_id[i]] + 3*point_indexes[i];

          double dist2 = 0.;
          for (int j = 0; j < 3; j++) {
            double delta = _pt[j] - _center[j];
            dist2 += delta*delta;
          }

          if (dist2 <= _radius2) {
            /* Check size and realloc if necessary */
            if (pib_idx[iball+1] >= s_pib) {
              s_pib *= 2;

              *ball_pts_l_num = realloc(ball_pts_l_num, sizeof(int   ) * s_pib * 2);
              *ball_pts_dist2 = realloc(ball_pts_dist2, sizeof(double) * s_pib);

              pib_l_num = *ball_pts_l_num;
              pib_dist2 = *ball_pts_dist2;
            }

            /* Add point */
            pib_l_num[2*pib_idx[iball+1]  ] = point_clouds_id[i];
            pib_l_num[2*pib_idx[iball+1]+1] = point_indexes[i];

            pib_dist2[pib_idx[iball+1]] = dist2;

            pib_idx[iball+1]++;
          }
        } // End of loop on current leaf's points
      }

      else {
        /* Internal node */
        for (int ichild = 0; ichild < n_children; ichild++) {

          int child_id = nodes->children_id[n_children*node_id + ichild];

          if (nodes->n_points[child_id] == 0) {
            continue;
          }

          inside_box = _box_dist2_min(3,
                                      &nodes->extents[6*child_id],
                                      _center,
                                      &min_dist2);

          if (inside_box || min_dist2 <= _radius2) {
            stack[pos_stack++] = child_id;
          }

        }

      }




    } // End of while loop


  } // End of loop on points
  free(stack);

  s_pib = pib_idx[n_ball];
  *ball_pts_l_num = realloc(*ball_pts_l_num, sizeof(int   ) * s_pib * 2);
  *ball_pts_dist2 = realloc(*ball_pts_dist2, sizeof(double) * s_pib);


}



/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_kdtree_seq_extract_extent
(
  PDM_kdtree_seq_t  *kdtree,
  int                root_id,
  int                n_depth,
  int               *n_node,
  double           **node_extents,
  int              **node_weight
)
{
  _l_nodes_t *nodes = kdtree->nodes;

  int n_children   = 2;
  int s_pt_stack   = ((n_children - 1) * (kdtree->depth_max - 1) + n_children);
  int *stack_id    = malloc (s_pt_stack * sizeof(int              ));
  int *stack_depth  = malloc (s_pt_stack * sizeof(int              ));

  // int n_extract_max = ((n_children - 1) * (kdtree->depth_max - 1) + n_children);
  int *id_to_extract = malloc( kdtree->n_nodes * sizeof(int));

  int n_extract = 0;
  int pos_stack = 0;
  stack_id   [pos_stack] = root_id;
  stack_depth[pos_stack] = 0;
  pos_stack++;
  while(pos_stack > 0) {

    /* Inspect node */
    --pos_stack;
    int node_id = stack_id   [pos_stack];
    int depth   = stack_depth[pos_stack];

    if(nodes->is_leaf[node_id] || depth == n_depth) {
      if(nodes->n_points[node_id] > 0) {
        id_to_extract[n_extract++] = node_id;
      }
    } else {
      for (int i = 0; i < n_children; i++) {
        int child_id = nodes->children_id[n_children*node_id+i];
        if (child_id < 0) {
          continue;
        }

        if(depth < n_depth) {
          stack_id   [pos_stack] = child_id;
          stack_depth[pos_stack] = depth + 1;
          pos_stack++;
        }
      }
    }
  }
  free(stack_id);
  free(stack_depth);

  double* _extents = malloc(n_extract * 6 * sizeof(double));
  int   * _n_pts   = malloc(n_extract *     sizeof(int   ));
  for(int i = 0; i < n_extract; ++i) {
    int node_id = id_to_extract[i];
    _n_pts[i] = nodes->n_points[node_id];
    for(int k = 0; k < 6; ++k) {
      _extents[6*i+k] = nodes->extents[6*node_id+k];
    }
  }

  *n_node       = n_extract;
  *node_extents = _extents;
  *node_weight  = _n_pts;

  free(id_to_extract);

}



#ifdef  __cplusplus
}
#endif
