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

#include "pdm_point_tree_seq.h"
#include "pdm_point_tree_seq_priv.h"

/*----------------------------------------------------------------------------*/


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
static const int    dbg_ptree   = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/


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
 int     split_direction,
 double  extents_min,
 double  extents_max,
 int     point_range[2],
 double *pts_coord
 )
{
  double mid = 0.5*(extents_min + extents_max);

  int n_pts = point_range[1] - point_range[0];
  if (n_pts > 2) {
    double *x = malloc(sizeof(double) * n_pts);

    for (int j = 0; j < n_pts; j++) {
      int i = point_range[0] + j;

      x[j] = pts_coord[3*i + split_direction];
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
       int     split_direction,
       int     point_range[2],
 const double *pts_coord,
       int     n_sample,
       double *sampling,
       double *cfreq,
       int     distrib[2]
 )
{
  double i_npts = 1. / (double) (point_range[1] - point_range[0]);

  int l_distrib[n_sample];
  for (int i = 0; i < n_sample; i++) {
    l_distrib[i] = 0;
  }


  for (int i = point_range[0]; i < point_range[1]; i++) {
    double x = pts_coord[3*i + split_direction];

    int isample = PDM_binary_search_gap_double(x,
                                               sampling,
                                               n_sample+1);

    if (isample >= n_sample || isample < 0) {
      PDM_log_trace_array_double(sampling, n_sample+1, "sampling : ");
      log_trace("!!! x = %f, isample = %d / %d\n", x, isample, n_sample);
    }

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
 const int     n_sample,
       int     split_direction,
       double  extents_min,
       double  extents_max,
       int     point_range[2],
 const double *pts_coord
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



static void
_build_point_tree_seq_leaves
(
 const int                         ancestor_id,
 const PDM_point_tree_seq_child_t  location_in_ancestor,
 const int                         depth,
 const double                      extents[],
       PDM_point_tree_seq_t       *ptree,
       int                         point_range[2],
       int                        *new_to_old
 )
{
  _l_nodes_t *nodes = ptree->nodes;

  if (dbg_ptree) {
    log_trace("\nnode_id = %d\n",
              ptree->n_nodes);
    log_trace("ancestor_id = %d, location_in_ancestor = %d, depth = %d, point_range = %d/%d\n",
              ancestor_id, (int) location_in_ancestor, depth, point_range[0], point_range[1]);
    log_trace("extents = %f %f %f  %f %f %f\n",
              extents[0], extents[1], extents[2], extents[3], extents[4], extents[5]);
  }

  int n_children = PDM_point_tree_n_children_get(ptree);

  /* Resize point_tree if necessary */
  int _n_nodes = ptree->n_nodes;
  int tmp_size = ptree->n_nodes;

  if (ptree->n_nodes >= ptree->n_nodes_max) {
    if (ptree->n_nodes == 0) {
      ptree->n_nodes     = 1;
      ptree->n_nodes_max = 8;
    }
    ptree->n_nodes_max *= 2;

    nodes->ancestor_id = realloc(nodes->ancestor_id, sizeof(int   ) * ptree->n_nodes_max);
    nodes->is_leaf     = realloc(nodes->is_leaf,     sizeof(int   ) * ptree->n_nodes_max);
    nodes->depth       = realloc(nodes->depth,       sizeof(int   ) * ptree->n_nodes_max);
    nodes->children_id = realloc(nodes->children_id, sizeof(int   ) * ptree->n_nodes_max * n_children);
    nodes->range       = realloc(nodes->range,       sizeof(int   ) * ptree->n_nodes_max * 2);
    nodes->idx         = realloc(nodes->idx,         sizeof(int   ) * ptree->n_nodes_max * (n_children+1));
    nodes->n_points    = realloc(nodes->n_points,    sizeof(int   ) * ptree->n_nodes_max);
    nodes->extents     = realloc(nodes->extents,     sizeof(double) * ptree->n_nodes_max * 6);
    nodes->location_in_ancestor = realloc(nodes->location_in_ancestor, sizeof(PDM_point_tree_seq_child_t) * ptree->n_nodes_max);
  }


  /* Number of points */
  int _n_points = point_range[1] - point_range[0];

  int idx[9];
  int child_id[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */

  int is_leaf = 1;
  if (depth < ptree->depth_max && _n_points > ptree->points_in_leaf_max) {

    /* Choose split direction */
    double max_range = -1.;
    int split_direction;

    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      for (int direction = 0; direction < 3; direction++) {
        double range = extents[3+direction] - extents[direction];
        if (range > max_range) {
          max_range        = range;
          split_direction = direction;
        }
      }
      if (dbg_ptree) {
        log_trace("split_direction = %d\n", split_direction);
      }
    }

    /* Choose split point */
    double mid[3];
    if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      if (1) {
        mid[0] = _approx_median_point(6,
                                      split_direction,
                                      extents[split_direction],
                                      extents[3+split_direction],
                                      point_range,
                                      ptree->_pts_coord);
      }
      else
      {
        mid[0] = _median_point(split_direction,
                               extents[split_direction],
                               extents[3+split_direction],
                               point_range,
                               ptree->_pts_coord);
      }
      if (dbg_ptree) {
        log_trace("mid = %f\n", mid[0]);
      }
    }
    else {
      for (int i = 0; i < 3; i++) {
        mid[i] = 0.5*(extents[i] + extents[3+i]);
      }
    }


    /* Count and reorder points in each child node */
    int count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }
      count[ichild]++;
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(count, n_children, "count : ");
    }


    /* Build index */
    idx[0] = 0;
    for (int ichild = 0; ichild < n_children; ichild++) {
      idx[ichild+1] = idx[ichild] + count[ichild];
      count[ichild] = 0;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      int ichild = 0;
      if (ptree->tree_type == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
        ichild = (ptree->_pts_coord[3*i + split_direction] > mid[0]);
      }
      else {
        for (int j = 0; j < 3; j++) {
          if (ptree->_pts_coord[3*i + j] > mid[j]) {
            ichild += octant_mask[j];
          }
        }
      }

      int old = ptree->new_to_old[i];
      int pos = point_range[0] + idx[ichild] + count[ichild];

      new_to_old[pos]         = old;
      ptree->old_to_new[old] = pos;

      count[ichild]++;
    }

    for (int i = point_range[0]; i < point_range[1]; i++) {
      ptree->new_to_old[i] = new_to_old[i];
    }
    if (dbg_ptree) {
      PDM_log_trace_array_int(ptree->new_to_old + point_range[0],
                              _n_points,
                              "new_to_old: ");
    }

    /* Reorder points */
    for (int i = point_range[0]; i < point_range[1]; i++) {
      memcpy(ptree->_pts_coord + 3*i,
             ptree->pts_coord  + 3*ptree->new_to_old[i],
             sizeof(double) * 3);
    }


    for (int i = 0; i <= n_children; i++) {
      idx[i] += point_range[0];
    }

    if (dbg_ptree) {
      PDM_log_trace_array_int(idx, n_children+1, "idx : ");
    }

    /* Build leaves recursively */
    double sub_extents[6];
    for (int ichild = 0; ichild < n_children; ichild++) {

      if (idx[ichild+1] <= idx[ichild]) {
        continue;
      }

      tmp_size++;

      child_id[ichild] = tmp_size;
      is_leaf = 0;

      // memcpy(sub_extents, extents, sizeof(double) * 6);
      // if (ichild == 0) {
      //   sub_extents[3+_split_direction] = mid;
      // }
      // else {
      //   sub_extents[_split_direction]   = mid;
      // }

      // Tight extents (fit contained points)
      for (int j = 0; j < 3; j++) {
        sub_extents[j  ] =  HUGE_VAL;
        sub_extents[j+3] = -HUGE_VAL;
      }
      for (int ipt = idx[ichild]; ipt < idx[ichild+1]; ipt++) {
        for (int j = 0; j < 3; j++) {
          double x = ptree->_pts_coord[3*ipt+j];
          sub_extents[j  ] = PDM_MIN(sub_extents[j  ], x);
          sub_extents[j+3] = PDM_MAX(sub_extents[j+3], x);
        }
      }
      for (int j = 0; j < 3; j++) {
        // if (sub_extents[j+3] < sub_extents[j] + _eps_default) {
        sub_extents[j  ] -= 0.5*_eps_default;
        sub_extents[j+3] += 0.5*_eps_default;
        // }
      }

      if (dbg_ptree) {
        log_trace("child %d, id %d, sub_extents = %f %f %f  %f %f %f\n",
                  ichild, child_id[ichild],
                  sub_extents[0], sub_extents[1], sub_extents[2],
                  sub_extents[3], sub_extents[4], sub_extents[5]);
      }

      ptree->n_nodes = tmp_size;

      _build_point_tree_seq_leaves(_n_nodes,
                                   (PDM_point_tree_seq_child_t) ichild,
                                   depth + 1,
                                   sub_extents,
                                   ptree,
                                   idx + ichild,
                                   new_to_old);

      tmp_size = ptree->n_nodes;
    }

  }

  /* Finalize node */
  for (int i = 0; i < 2; i++) {
    nodes->range[2*_n_nodes + i] = point_range[i];
  }

  for (int i = 0; i <= n_children; i++) {
    nodes->idx[(n_children+1)*_n_nodes + i] = idx[i];
  }

  for (int i = 0; i < 6; i++) {
    nodes->extents[6*_n_nodes + i] = extents[i];
  }

  for (int i = 0; i < n_children; i++) {
    nodes->children_id[n_children*_n_nodes + i] = child_id[i];
  }

  nodes->is_leaf[_n_nodes] = is_leaf;

  nodes->ancestor_id[_n_nodes] = ancestor_id;
  nodes->depth[_n_nodes]       = depth;

  nodes->n_points[_n_nodes] = _n_points;
  nodes->location_in_ancestor[_n_nodes] = location_in_ancestor;
}



/**
 *
 * \brief Build a point_tree
 *
 * \param[in]  ptree    Current point_tree
 * .
 */

static void
_build_point_tree
(
 PDM_point_tree_seq_t *ptree
)
{
  int point_range[2];

  /* Initialization */

  ptree->n_nodes     = 0;
  ptree->n_nodes_max = 0;

  ptree->nodes = malloc(sizeof(_l_nodes_t));
  ptree->nodes->ancestor_id          = NULL;
  ptree->nodes->is_leaf              = NULL;
  ptree->nodes->location_in_ancestor = NULL;
  ptree->nodes->depth                = NULL;
  ptree->nodes->children_id          = NULL;
  ptree->nodes->range                = NULL;
  ptree->nodes->idx                  = NULL;
  ptree->nodes->n_points             = NULL;
  ptree->nodes->extents              = NULL;

  for (int i = 0; i < ptree->n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      double x = ptree->pts_coord[3*i+j];
      ptree->extents[j  ] = PDM_MIN(ptree->extents[j  ], x);
      ptree->extents[j+3] = PDM_MAX(ptree->extents[j+3], x);
    }
  }

  ptree->new_to_old = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->new_to_old[i] = i;
  }

  ptree->old_to_new = malloc(sizeof(int) * ptree->n_pts);
  for (int i = 0; i < ptree->n_pts; i++) {
    ptree->old_to_new[i] = i;
  }

  ptree->_pts_coord = malloc(sizeof(double) * ptree->n_pts * 3);
  memcpy(ptree->_pts_coord, ptree->pts_coord, sizeof(double) * ptree->n_pts * 3);

  double delta = -1;
  for (int i = 0; i < 3; i++) {
    delta = PDM_MAX (ptree->tolerance * (ptree->extents[i + 3] - ptree->extents[i]),
                     delta);
  }
  delta = PDM_MAX (delta,_eps_default);

  for (int i = 0; i < 3; i++) {
    ptree->extents[i  ] += -1.01*delta;
    ptree->extents[i+3] += delta;;
  }

  point_range[0] = 0;
  point_range[1] = ptree->n_pts;


  /* Build point_tree recursively */

  if (dbg_ptree) {
    log_trace(">> _build_point_tree_seq_leaves\n");
  }
  int *tmp_new_to_old = malloc(sizeof(int) * ptree->n_pts);
  _build_point_tree_seq_leaves(-1,
                               (PDM_point_tree_seq_child_t) 0,
                               -1,
                               ptree->extents,
                               ptree,
                               point_range,
                               tmp_new_to_old);
  free(tmp_new_to_old);


  if (ptree->n_nodes > 1) {
    ptree->n_nodes += 1;
  }

  if (dbg_ptree) {
    // PDM_log_trace_array_int(ptree->old_to_new,
    //                         ptree->n_pts,
    //                         "old_to_new : ");
    // PDM_log_trace_array_int(ptree->new_to_old,
    //                         ptree->n_pts,
    //                         "new_to_old : ");
    for (int i = 0; i < ptree->n_pts; i++) {
      if (ptree->new_to_old[ptree->old_to_new[i]] != i) {
        log_trace("!!! point %d error with old_to_new_to_old\n", i);
      }
    }

    // Dump kd-tree
    _l_nodes_t *nodes = ptree->nodes;

    int n_children = PDM_point_tree_n_children_get(ptree);

    for (int i = 0; i < ptree->n_nodes; i++) {
      if (1) {//nodes->is_leaf[i]) {
        log_trace("\nNode %d :", i);
        log_trace("  depth = %d\n", nodes->depth[i]);
        log_trace("  is_leaf = %d\n", nodes->is_leaf[i]);
        // log_trace("  children_id = %d %d\n", nodes->children_id[2*i], nodes->children_id[2*i+1]);
        PDM_log_trace_array_int(nodes->children_id + n_children*i,
                                n_children,
                                "  children_id : ");
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

}




static void
_l_nodes_free
(
 PDM_point_tree_seq_t *ptree
 )
{
  if (ptree->nodes != NULL) {

    free(ptree->nodes->ancestor_id);
    free(ptree->nodes->is_leaf);
    free(ptree->nodes->depth);
    free(ptree->nodes->location_in_ancestor);
    free(ptree->nodes->children_id);
    free(ptree->nodes->range);
    free(ptree->nodes->idx);
    free(ptree->nodes->n_points);
    free(ptree->nodes->extents);

    free(ptree->nodes);

    ptree->nodes = NULL;
  }
}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a point_tree structure
 *
 * \param [in]   tree_type          Tree type
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_point_tree_seq object
 */

PDM_point_tree_seq_t *
PDM_point_tree_seq_create
(
 const PDM_doctree_local_tree_t tree_type,
 const int                      depth_max,
 const int                      points_in_leaf_max,
 const double                   tolerance
)
{
  PDM_point_tree_seq_t *ptree = (PDM_point_tree_seq_t *) malloc(sizeof(PDM_point_tree_seq_t));

  ptree->tree_type = tree_type;

  ptree->depth_max          = depth_max;
  ptree->points_in_leaf_max = points_in_leaf_max;
  ptree->tolerance          = tolerance;

  ptree->nodes = NULL;

  ptree->n_nodes     = 0;
  ptree->n_nodes_max = 0;

  for (int i = 0; i < 3; i++) {
    ptree->extents[i]     =  HUGE_VAL;
    ptree->extents[i + 3] = -HUGE_VAL;
  }

  ptree->n_pts = 0;
  ptree->pts_coord  = NULL;
  ptree->_pts_coord = NULL;
  ptree->new_to_old = NULL;

  return ptree;
}


/**
 *
 * \brief Free a point_tree structure
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_free
(
 PDM_point_tree_seq_t *ptree
)
{
  if (ptree->_pts_coord != NULL) {
    free(ptree->_pts_coord);
  }

  if (ptree->new_to_old != NULL) {
    free(ptree->new_to_old);
  }

  if (ptree->old_to_new != NULL) {
    free(ptree->old_to_new);
  }

  _l_nodes_free(ptree);

  free(ptree);
}


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 * \param [in]   n_pts             Number of points in cloud
 * \param [in]   pts_coord         Point coordinates
 *
 */

void
PDM_point_tree_seq_point_cloud_set
(
       PDM_point_tree_seq_t *ptree,
 const int                   n_pts,
 const double               *pts_coord
)
{
  ptree->n_pts     = n_pts;
  ptree->pts_coord = pts_coord;
}


/**
 *
 * \brief Build point_tree
 *
 * \param [in]   ptree             Pointer to \ref PDM_point_tree_seq object
 *
 */

void
PDM_point_tree_seq_build
(
 PDM_point_tree_seq_t *ptree
)
{
  if (ptree->nodes == NULL) {
    if (dbg_ptree) {
      log_trace(">> _build_point_tree\n");
    }
    _build_point_tree(ptree);
  }

}


/**
 *
 * \brief Write point_tree nodes in a VTK file
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   filename              Output file name
 *
 */

void PDM_point_tree_seq_write_nodes
(
       PDM_point_tree_seq_t *ptree,
 const char                 *filename
 )
{
  _l_nodes_t *nodes = ptree->nodes;

  double tol_visu = 1e-3;
  double _ext[6];

  // write VTK
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "ptree_seq\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 8*ptree->n_nodes);
  for (int inode = 0; inode < ptree->n_nodes; inode++) {
    double *ext = nodes->extents + 6*inode;

    if (1 && (ptree->nodes->range[2*inode+1] - ptree->nodes->range[2*inode] == 1)) {
      // Trick to visualize nodes with degenerate extents (single point)
      ext = _ext;
      for (int i = 0; i < 3; i++) {
        double x = nodes->extents[6*inode + i];
        double eps = tol_visu*(ptree->extents[i+3] - ptree->extents[i]);
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

  fprintf(f, "CELLS %d %d\n", ptree->n_nodes, 9*ptree->n_nodes);
  for (int inode = 0; inode < ptree->n_nodes; inode++) {
    fprintf(f, "8 ");
    for (int j = 0; j < 8; j++) {
      fprintf(f, "%d ", 8*inode+j);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELL_TYPES %d\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", 12);
  }

  fprintf(f, "CELL_DATA %d\n", ptree->n_nodes);

  fprintf(f, "FIELD node_field 3\n");
  fprintf(f, "depth 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->depth[i]);
  }
  fprintf(f, "is_leaf 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->is_leaf[i]);
  }
  fprintf(f, "n_pts 1 %d int\n", ptree->n_nodes);
  for (int i = 0; i < ptree->n_nodes; i++) {
    fprintf(f, "%d\n", ptree->nodes->range[2*i+1] - ptree->nodes->range[2*i]);
  }

  fclose(f);
}


/**
 *
 * \brief Get number of children per node in  point_tree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 *
 * \return   Number of children per node in point_tree
 */

int
PDM_point_tree_n_children_get
(
 PDM_point_tree_seq_t *ptree
 )
{
  switch (ptree->tree_type) {
    case PDM_DOCTREE_LOCAL_TREE_OCTREE:
    case PDM_DOCTREE_LOCAL_TREE_LINEAR_OCTREE:
    return 8;

    case PDM_DOCTREE_LOCAL_TREE_KDTREE:
    return 2;

    default:
    PDM_error(__FILE__, __LINE__, 0,
              "Invalid tree_type %d\n", (int) ptree->tree_type);
    break;
  }

  return -1;
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_ptree_seq object
 * \param [out]  new_to_old             New to old order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_new_to_old_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **new_to_old
)
{
  *new_to_old = ptree->new_to_old;
}


/**
 *
 * \brief Get point order in ptree
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  old_to_new             Old to new order of points in ptree
 *
 */

void
PDM_point_tree_seq_point_old_to_new_get
(
 PDM_point_tree_seq_t  *ptree,
 int                  **old_to_new
)
{
  *old_to_new = ptree->old_to_new;
}


/**
 *
 * \brief Get point coords in point_tree's order
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [out]  pts_coord             Point coordinates
 *
 */

void
PDM_point_tree_seq_sorted_points_get
(
 PDM_point_tree_seq_t  *ptree,
 double               **pts_coord
)
{
  *pts_coord = ptree->_pts_coord;
}


/**
 *
 * \brief Get point range of a point_tree node
 *
 * \param [in]   ptree                 Pointer to \ref PDM_point_tree_seq object
 * \param [in]   node_id               Current node ID (zero-based)
 * \param [out]  point_range           Point range of current node
 *
 * \return   Number of points inside current node
 *
 */

int
PDM_point_tree_seq_point_range_get
(
       PDM_point_tree_seq_t *ptree,
 const int                   node_id,
       int                  *point_range
)
{
  assert(node_id < ptree->n_nodes);
  for (int i = 0; i < 2; i++) {
    point_range[i] = ptree->nodes->range[2*node_id + i];
  }

  return point_range[1] - point_range[0];
}



/**
 *
 * \brief Get node extents of subtree of given depth and starting from given root
 *
 * \param [in]  kdtree               Pointer to \ref PDM_kdtree_seq object
 * \param [in]  root_id              ID of subtree root
 * \param [in]  n_depth              Depth of subtree
 * \param [out] n_node               Number of subtree nodes
 * \param [out] node_ids             IDs of subtree nodes
 * \param [out] node_extents         Extents of subtree nodes
 * \param [out] node_weight          Weights of subtree nodes
 *
 */

void
PDM_point_tree_seq_extract_nodes
(
  PDM_point_tree_seq_t  *ptree,
  int                    root_id,
  int                    n_depth,
  int                   *n_node,
  int                  **node_ids,
  double               **node_extents,
  int                  **node_weight
)
{
  _l_nodes_t *nodes = ptree->nodes;

  int n_children   = PDM_point_tree_n_children_get(ptree);
  int s_pt_stack   = ((n_children - 1) * (ptree->depth_max - 1) + n_children);
  int *stack_id    = malloc (s_pt_stack * sizeof(int));
  int *stack_depth = malloc (s_pt_stack * sizeof(int));

  int *id_to_extract = malloc(ptree->n_nodes * sizeof(int));

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
  *node_ids     = id_to_extract;
  *node_extents = _extents;
  *node_weight  = _n_pts;
}
