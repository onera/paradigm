
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
#include "pdm_mpi.cuh"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_cuda_error.cuh"
#include "pdm_cuda.cuh"
#include "pdm_gnum.h"
#include "pdm_morton.h"
#include "pdm_morton.cuh"
#include "pdm_handles.h"
#include "pdm_handles.cuh"
#include "pdm_para_octree.h"
#include "pdm_para_octree.cuh"
#include "pdm_timer.h"
#include "pdm_timer.cuh"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"


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

#define N_TIMER_KNN 12
typedef enum {
  KNN_BEGIN,
  KNN_INIT,
  KNN_MAIN,
  LOOP_BEGIN,
  LOOP_FILTER,
  LOOP_LOCAL_SEARCH,
  LOOP_PTB_EXCH,
  LOOP_MERGE_PTB,
  LOOP_PREP_NEXT,
  LOOP_END,
  KNN_BTP_EXCH,
  KNN_END,
} _timer_step_knn_t;


/*============================================================================
 * Global variable
 *============================================================================*/

//__device__ static PDM_Handles_t *d_octrees    = NULL;
//static const double _eps_default  = 1.e-12;

/*============================================================================
 * Private function definitions
 *============================================================================*/


// static void
// _closest_points_local
// (
//  const _octree_t  *octree,
//  const int         n_closest_points,
//  const int         n_pts,
//  double           *pts_coord,
//  PDM_g_num_t      *pts_g_num, // ONLY FOR DEBUG
//  int              *start_leaves,
//  int              *start_leaves_idx,
//  double           *upper_bound_dist,
//  PDM_g_num_t      *local_closest_src_gnum,
//  double           *local_closest_src_dist,
//  int              *send_count,
//  int              *send_shift,
//  int             **send_tgt_lnum,
//  int             **send_start_leaves,
//  int             **send_start_leaves_count,
//  int             **send_start_leaves_rank_shift
//  )
// {
//   const int DEBUG = 0;
//   const int CHECK_FACE_DIST = 1;
//   //const int CHECK_EMPTY_START = 0;
//   const int CHECK_CLOSEST = 1;
//   const double THRESHOLD_CLOSEST = 0.75; // check closest if min_dist > threshold * upper bound
//   const int CHECK_INTERSECT = 1;
//   const int NORMALIZE = 1;

//   const _l_octant_t *octants = octree->octants;
//   const int dim = octree->dim;

//   int myRank, lComm;
//   PDM_MPI_Comm_rank (octree->comm, &myRank);
//   PDM_MPI_Comm_size (octree->comm, &lComm);

//   for (int i = 0; i < lComm; i++) {
//     send_count[i] = 0;
//   }

//   //-->> DETAIL TIMERS
//   double timer_ls[N_TIMER_LS];
//   for (int i = 0; i < N_TIMER_LS; i++) {
//     timer_ls[i] = 0;
//   }
//   double b_timer, e_timer, start_main;

//   PDM_timer_hang_on(octree->timer);
//   timer_ls[LS_BEGIN] = PDM_timer_elapsed(octree->timer);
//   b_timer = timer_ls[LS_BEGIN];
//   PDM_timer_resume(octree->timer);
//   //<<--

//   //--->>>
//   int **tmp_send_tgt_lnum       = malloc (sizeof(int *) * lComm);
//   int **tmp_send_tgt_n_leaves   = malloc (sizeof(int *) * lComm);
//   int **tmp_send_start_leaves   = malloc (sizeof(int *) * lComm);
//   int  *s_tmp_send_start_leaves = malloc (sizeof(int) * lComm);
//   int  *n_tmp_send_start_leaves = malloc (sizeof(int) * lComm);

//   int **send_to_rank_leaves   = malloc (sizeof(int *) * lComm);
//   int  *n_send_to_rank_leaves = malloc (sizeof(int) * lComm);
//   int  *s_send_to_rank_leaves = malloc (sizeof(int) * lComm);
//   for (int i = 0; i < lComm; i++) {
//     s_send_to_rank_leaves[i] = 16;// ?
//     if (i != myRank) {
//       send_to_rank_leaves[i] = malloc (sizeof(int) * 2 * s_send_to_rank_leaves[i]);

//       tmp_send_tgt_lnum[i] = malloc (sizeof(int) * n_pts);// could be smaller and dynamically re-alloc'ed when necessary
//       tmp_send_tgt_n_leaves[i] = malloc (sizeof(int) * n_pts);// idem
//       s_tmp_send_start_leaves[i] = 16;//?
//       tmp_send_start_leaves[i] = malloc (sizeof(int) * s_tmp_send_start_leaves[i]);
//     } else {
//       s_tmp_send_start_leaves[i] = 0;
//     }
//     n_tmp_send_start_leaves[i] = 0;
//   }
//   //<<<---

//   int *is_visited_part = malloc (sizeof(int) * octree->n_connected);

//   int *is_visited     = malloc (sizeof(int) * octants->n_nodes);
//   int *visited_leaves = malloc (sizeof(int) * octants->n_nodes);// could be smaller...
//   int n_visited = 0;
//   for (int i = 0; i < octants->n_nodes; i++) {
//     is_visited[i] = 0;
//   }

//   /* Min heap used to sort start leaves */
//   _min_heap_t *start_heap = _min_heap_create (start_leaves_idx[n_pts]); // smaller size?

//   /* Min heap used to visit leaves from neighbour to neighbour */
//   _min_heap_t *leaf_heap = _min_heap_create (octants->n_nodes);

//   //-->> DETAIL TIMERS
//   PDM_timer_hang_on(octree->timer);
//   e_timer = PDM_timer_elapsed(octree->timer);
//   timer_ls[LS_INIT] = e_timer - b_timer;
//   PDM_timer_resume(octree->timer);
//   //<<--


//   /* Loop over target points */
//   for (int i_tgt = 0; i_tgt < n_pts; i_tgt++) {

//     //-->> DETAIL TIMERS
//     PDM_timer_hang_on(octree->timer);
//     b_timer = PDM_timer_elapsed(octree->timer);
//     PDM_timer_resume(octree->timer);
//     //<<--

//     /* Init */
//     for (int i = 0; i < octree->n_connected; i++) {
//       is_visited_part[i] = 0;
//     }

//     for (int i = 0; i < lComm; i++) {
//       n_send_to_rank_leaves[i] = 0;
//     }

//     for (int i = 0; i < n_visited; i++) {
//       is_visited[visited_leaves[i]] = 0;
//     }
//     n_visited = 0;

//     PDM_g_num_t *closest_src_gnum = local_closest_src_gnum + n_closest_points * i_tgt;
//     double      *closest_src_dist = local_closest_src_dist + n_closest_points * i_tgt;
//     for (int j = 0; j < n_closest_points; j++) {
//       closest_src_dist[j] = upper_bound_dist[i_tgt];
//       closest_src_gnum[j] = -1; // USEFUL ONLY FOR DEBUG
//     }

//     double *max_src_dist = closest_src_dist + n_closest_points - 1;


//     /* Get current target point data */
//     const double *_pt = pts_coord + i_tgt * dim;
//     double _ptn[3];
//     if (NORMALIZE) {
//       for (int i_dim = 0; i_dim < dim; i_dim++) {
//         _ptn[i_dim] = (_pt[i_dim] - octree->s[i_dim]) / octree->d[i_dim];
//       }
//     }

//     int n_start_leaves = start_leaves_idx[i_tgt+1] - start_leaves_idx[i_tgt];
//     if (DEBUG) {
//       printf("\n=== pt (%ld) (upper_bound_dist = %f) ===\nstart leaves (%d):\n",
//              pts_g_num[i_tgt], upper_bound_dist[i_tgt], n_start_leaves);
//     }

//     //-->> DETAIL TIMERS
//     PDM_timer_hang_on(octree->timer);
//     e_timer = PDM_timer_elapsed(octree->timer);
//     timer_ls[LS_POINT_INIT] += e_timer - b_timer;
//     b_timer = e_timer;
//     PDM_timer_resume(octree->timer);
//     //<<--

//     if (n_start_leaves < 1) {
//       continue; /* move on to next target point */
//     }

//     /* Sort start leaves in ascending order of min distance from tgt point */
//     _min_heap_reset (start_heap);
//     double min_start_dist = *max_src_dist;
//     for (int i_start = start_leaves_idx[i_tgt]; i_start < start_leaves_idx[i_tgt+1]; i_start++) {

//       int leaf_id = start_leaves[i_start];
//       int part_id = _get_octant_part_id (octree,
//                                          leaf_id);
//       is_visited_part[part_id] = 1;

//       double start_dist;
//       if (NORMALIZE) {
//         start_dist = _octant_min_dist2_normalized (dim,
//                                                    octants->codes[leaf_id],
//                                                    octree->d,
//                                                    _ptn);
//       } else {
//         start_dist = _octant_min_dist2 (dim,
//                                         octants->codes[leaf_id],
//                                         octree->d,
//                                         octree->s,
//                                         _pt);
//       }

//       min_start_dist = PDM_MIN (min_start_dist, start_dist);

//       if (start_dist < *max_src_dist) {
//         _min_heap_push (start_heap,
//                         leaf_id,
//                         0,
//                         start_dist);
//       }

//       if (DEBUG) {
//         int i_part = _get_octant_part_id (octree, leaf_id);
//         printf("\t%d (part %d): dist = %f\n", leaf_id, i_part, start_dist);
//       }
//     }
//     if (DEBUG) {
//       printf("============================\n");
//     }

//     //-->> DETAIL TIMERS
//     PDM_timer_hang_on(octree->timer);
//     e_timer = PDM_timer_elapsed(octree->timer);
//     timer_ls[LS_SORT_START_LEAVES] += e_timer - b_timer;
//     b_timer = e_timer;
//     PDM_timer_resume(octree->timer);
//     //<<--

//     /* Check whether start_heap is empty */
//     //if (CHECK_EMPTY_START) {
//     //if (start_heap->count == 0 && n_start_leaves > 0) {
//     if (CHECK_CLOSEST && n_start_leaves > 0) {
//       if (min_start_dist >= THRESHOLD_CLOSEST * (*max_src_dist)) {
//         /* Look for the closest octant (within the search radius), which would have been missed if we returned too early) */
//         double            box[2*dim];
//         PDM_morton_code_t box_corners[2];
//         double            width = sqrt (*max_src_dist);
//         double s[3], d[3];
//         /* Encode corners of search box */
//         for (int i_dim = 0; i_dim < dim; i_dim++) {
//           box[i_dim]       = PDM_MAX (octree->global_extents[i_dim],
//                                       _pt[i_dim] - width);
//           box[dim + i_dim] = PDM_MAX (octree->global_extents[dim + i_dim],
//                                       _pt[i_dim] + width);
//         }

//         PDM_morton_encode_coords (dim,
//                                   PDM_morton_max_level,
//                                   octree->global_extents,
//                                   (size_t) 2,
//                                   box,
//                                   box_corners,
//                                   d,
//                                   s);

//         /* Loop over visited parts */
//         for (int i_part = 0; i_part < octree->n_connected; i_part++) {

//           if (!is_visited_part[i_part]) {
//             continue;
//           }

//           /* Get range of octants that intersect the search box in the sense of Morton codes */
//           size_t l, r;
//           _maximal_intersecting_range (box_corners[0],
//                                        box_corners[1],
//                                        octants->codes + octree->connected_idx[i_part],
//                                        octree->connected_idx[i_part+1] - octree->connected_idx[i_part],
//                                        &l,
//                                        &r);

//           l += octree->connected_idx[i_part];
//           r += octree->connected_idx[i_part];

//           /* Narrow down this range as much as possible */
//           // BIGMIN/LITMAX? ...

//           /* Find closest leaf within this range */
//           int    closest_id = -1;
//           double min_dist = *max_src_dist;
//           double dist;
//           for (int i = l; i < r; i++) {

//             if (CHECK_INTERSECT) {
//               /* make sure leaf intersects search box so it's worth computing min distance */
//               int intersect = 1;
//               PDM_morton_code_t *code = octants->codes + i;

//               double min_octant[dim];
//               double max_octant[dim];
//               double side = 1. / pow (2., code->L);

//               for (int j = 0; j < dim; j++) {
//                 min_octant[j] = octree->s[j] + octree-> d[j] * side * code->X[j];
//                 max_octant[j] = min_octant[j] + octree->d[j] * side;

//                 if (max_octant[j] < box[j] || min_octant[j] > box[dim + j]) {
//                   intersect = 0;
//                   break;
//                 }
//               }

//               if (intersect) {
//                 dist = 0.;
//                 double delta;
//                 for (int j = 0; j < dim; j++) {
//                   double x = _pt[j];

//                   if (x > max_octant[j]) {
//                     delta = x - max_octant[j];
//                     dist += delta * delta;
//                   } else if (x < min_octant[j]) {
//                     delta = x - min_octant[j];
//                     dist += delta * delta;
//                   }
//                 }
//               } else {
//                 dist = HUGE_VAL;
//               }

//             } else {
//               if (NORMALIZE) {
//                 dist = _octant_min_dist2_normalized (dim,
//                                                      octants->codes[i],
//                                                      octree->d,
//                                                      _ptn);
//               } else {
//                 dist = _octant_min_dist2 (dim,
//                                           octants->codes[i],
//                                           octree->d,
//                                           octree->s,
//                                           _pt);
//               }
//             } // end if (CHECK_INTERSECT)

//             if (dist < min_dist) {
//               closest_id = i;
//               min_dist = dist;
//             }
//           }

//           if (DEBUG) {
//             printf("[%d] closest octant: %d, dist = %f / %f\n",
//                    myRank, closest_id, min_dist, *max_src_dist);
//           }

//           /* If the closest octant is within the search radius, push it in start_heap */
//           if (closest_id >= 0) {
//             _min_heap_push (start_heap,
//                             closest_id,
//                             0,
//                             min_dist);
//           }
//         } // end loop over visited parts
//       } // end if (min_start_dist >= THRESHOLD_CLOSEST * (*max_src_dist))
//     } // end if (CHECK_CLOSEST && n_start_leaves > 0)

//     //-->> DETAIL TIMERS
//     PDM_timer_hang_on(octree->timer);
//     e_timer = PDM_timer_elapsed(octree->timer);
//     timer_ls[LS_CHECK_CLOSEST] += e_timer - b_timer;
//     b_timer = e_timer;
//     PDM_timer_resume(octree->timer);
//     //<<--


//     /* Loop over (sorted) start leaves */
//     int start_id;
//     PDM_g_num_t unused;
//     double start_dist;
//     while (_min_heap_pop (start_heap, &start_id, &unused, &start_dist)) {

//       if (DEBUG) {
//         printf("tgt point (%ld) start leaf %d: dist = %f / %f  (is visited? %d)\n",
//                pts_g_num[i_tgt], start_id, start_dist, *max_src_dist, is_visited[start_id]);
//       }

//       if (start_dist >= *max_src_dist) {
//         break;
//       }

//       if (is_visited[start_id]) {
//         continue;
//       }

//       /* Push start leaf in priority queue */
//       _min_heap_reset (leaf_heap);
//       _min_heap_push (leaf_heap,
//                       start_id,
//                       0,
//                       start_dist);
//       is_visited[start_id] = 1;
//       visited_leaves[n_visited++] = start_id;

//       /* Visit octree leaves from neighbour to neighbour using priority queue */
//       int leaf_id;
//       double leaf_dist;
//       while (_min_heap_pop (leaf_heap, &leaf_id, &unused, &leaf_dist)) {

//         if (DEBUG) {
//           printf("tgt point (%ld) inspecting leaf %d: dist = %f / %f\n",
//                  pts_g_num[i_tgt], leaf_id, leaf_dist, *max_src_dist);
//         }

//         if (leaf_dist >= *max_src_dist) {
//           break;
//         }

//         /* inspect source points inside popped leaf */
//         for (int i = 0; i < octants->n_points[leaf_id]; i++) {

//           // get source point coords and gnum
//           int i_src = octants->range[leaf_id] + i;
//           double *src_pt = octree->points + i_src * dim;
//           PDM_g_num_t src_gnum = octree->points_gnum[i_src];

//           // compute (squared) distance from target point
//           double src_dist = 0;
//           for (int j = 0; j < dim; j++) {
//             double delta = _pt[j] - src_pt[j];
//             src_dist += delta * delta;
//           }

//           if (DEBUG && pts_g_num[i_tgt] == 38) {
//             printf("  src pt (%ld) [%f %f %f] at dist %f / %f\n",
//                    src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist);
//           }

//           // insertion sort
//           if (src_dist < *max_src_dist) {
//             int j = n_closest_points - 1;
//             while (j > 0 && src_dist < closest_src_dist[j-1]) {
//               closest_src_gnum[j] = closest_src_gnum[j-1];
//               closest_src_dist[j] = closest_src_dist[j-1];
//               j--;
//             }

//             closest_src_gnum[j] = src_gnum;
//             closest_src_dist[j] = src_dist;

//             if (DEBUG) {
//               printf("  src pt (%ld) [%f %f %f] at dist %f / %f --> insert at pos %d\n",
//                      src_gnum, src_pt[0], src_pt[1], src_pt[2], src_dist, *max_src_dist, j);
//             }
//           }
//         } // end loop over source points inside popped leaf


//         /* inspect neighbours of popped leaf */
//         double side = 1./pow(2., octants->codes[leaf_id].L);
//         int check_dist;
//         for (PDM_para_octree_direction_t dir = PDM_BOTTOM; dir < 6; dir++) {
//           if (CHECK_FACE_DIST) {
//             /* coarse test to discard neighbours that are clearly beyond the current search radius */
//             int i_dim = dir / 2;
//             int sign = dir % 2;

//             double x_face = octree->s[i_dim] + octree->d[i_dim] * side * (octants->codes[leaf_id].X[i_dim] + sign);

//             if (sign == 0) {
//               check_dist = _pt[i_dim] > x_face;
//             } else {
//               check_dist = _pt[i_dim] < x_face;
//             }

//             if (check_dist) {
//               double dist_face = (_pt[i_dim] -  x_face) * (_pt[i_dim] -  x_face);

//               if (dist_face >= *max_src_dist) {
//                 /* all neighbours in direction dir are at least
//                    at a distance from the tgt point equal to dist_face
//                    so they need not be inspected */
//                 continue;
//               }
//             }
//           }

//           /* inspect neighbours in direction dir */
//           for (int i = octants->neighbour_idx[6*leaf_id + dir];
//                i < octants->neighbour_idx[6*leaf_id + dir + 1]; i++) {

//             int ngb = octants->neighbours[i];

//             if (ngb < 0) {
//               // distant neighbour(s)
//               ngb = -(ngb + 1);

//               for (int j = octree->part_boundary_elt_idx[ngb];
//                    j < octree->part_boundary_elt_idx[ngb+1]; j++) {

//                 int ngb_rank = octree->part_boundary_elt[3*j];
//                 int ngb_leaf = octree->part_boundary_elt[3*j+1];
//                 int ngb_part = octree->part_boundary_elt[3*j+2];

//                 if (n_send_to_rank_leaves[ngb_rank] == 0) {
//                   tmp_send_tgt_lnum[ngb_rank][send_count[ngb_rank]] = i_tgt;
//                   tmp_send_tgt_n_leaves[ngb_rank][send_count[ngb_rank]] = 0;
//                   send_count[ngb_rank]++;
//                 }

//                 if (s_send_to_rank_leaves[ngb_rank] <= n_send_to_rank_leaves[ngb_rank]) {
//                   s_send_to_rank_leaves[ngb_rank] *= 2;
//                   send_to_rank_leaves[ngb_rank] = realloc (send_to_rank_leaves[ngb_rank],
//                                                            sizeof(int) * 2 * s_send_to_rank_leaves[ngb_rank]);
//                 }
//                 if (DEBUG) {
//                   printf("  Send pt (%ld) to rank %d, leaf %d (part %d)\n",
//                          pts_g_num[i_tgt], ngb_rank, ngb_leaf, ngb_part);
//                 }
//                 send_to_rank_leaves[ngb_rank][2*n_send_to_rank_leaves[ngb_rank]]   = ngb_leaf;
//                 send_to_rank_leaves[ngb_rank][2*n_send_to_rank_leaves[ngb_rank]+1] = ngb_part;
//                 n_send_to_rank_leaves[ngb_rank]++;

//                 tmp_send_tgt_n_leaves[ngb_rank][send_count[ngb_rank]-1]++;

//               } // end loop over distant neighbours (j)

//             } else {
//               // local neighbour
//               if (is_visited[ngb] == 1) continue;

//               is_visited[ngb] = 1;
//               visited_leaves[n_visited++] = ngb;

//               // compute min dist from target point to current neighbor leaf
//               double ngb_min_dist;
//               if (NORMALIZE) {
//                 ngb_min_dist = _octant_min_dist2_normalized (dim,
//                                                              octants->codes[ngb],
//                                                              octree->d,
//                                                              _ptn);
//               } else {
//                 ngb_min_dist  = _octant_min_dist2 (dim,
//                                                    octants->codes[ngb],
//                                                    octree->d,
//                                                    octree->s,
//                                                    _pt);
//               }

//               if (ngb_min_dist < *max_src_dist) {
//                 // push current neighbour in priority queue
//                 _min_heap_push (leaf_heap,
//                                 ngb,
//                                 0,
//                                 ngb_min_dist);
//               }

//             } // end if local/distant neighbour
//           } // end loop over neighbours in direction dir (i)
//         } // end loop over directions (dir)
//       } // end while (_min_heap_pop (leaf_heap))

//     } // end loop over (sorted) start leaves

//     for (int rank = 0; rank < lComm; rank++) {
//       if (rank == myRank) continue;

//       if (s_tmp_send_start_leaves[rank] <= n_tmp_send_start_leaves[rank] + n_send_to_rank_leaves[rank]) {
//         s_tmp_send_start_leaves[rank] = PDM_MAX (2 * s_tmp_send_start_leaves[rank],
//                                                  n_tmp_send_start_leaves[rank] + n_send_to_rank_leaves[rank]);
//         tmp_send_start_leaves[rank] = realloc (tmp_send_start_leaves[rank],
//                                                sizeof(int) * s_tmp_send_start_leaves[rank]);
//       }

//       for (int i = 0; i < n_send_to_rank_leaves[rank]; i++) {
//         tmp_send_start_leaves[rank][n_tmp_send_start_leaves[rank]++] = send_to_rank_leaves[rank][2*i];
//       }
//     }

//     //-->> DETAIL TIMERS
//     PDM_timer_hang_on(octree->timer);
//     e_timer = PDM_timer_elapsed(octree->timer);
//     timer_ls[LS_TRAVERSAL] += e_timer - b_timer;
//     PDM_timer_resume(octree->timer);
//     //<<--

//   } // end loop over target points (i_tgt)
//   free (is_visited);
//   free (visited_leaves);
//   _min_heap_free (start_heap);
//   _min_heap_free (leaf_heap);

//   //-->> DETAIL TIMERS
//   PDM_timer_hang_on(octree->timer);
//   b_timer = PDM_timer_elapsed(octree->timer);
//   PDM_timer_resume(octree->timer);
//   //<<--

//   free (n_send_to_rank_leaves);
//   free (s_send_to_rank_leaves);
//   for (int i = 0; i < lComm; i++) {
//     if (i != myRank) {
//       free (send_to_rank_leaves[i]);
//     }
//   }
//   free (send_to_rank_leaves);

//   /* Manage send_tgt_lnum buffer */
//   send_shift[0] = 0;
//   for (int i = 0; i < lComm; i++) {
//     send_shift[i+1] = send_shift[i] + send_count[i];
//   }

//   *send_tgt_lnum = malloc (sizeof(int) * send_shift[lComm]);
//   int idx = 0;
//   for (int i = 0; i < lComm; i++) {
//     for (int j = 0; j < send_count[i]; j++) {
//       (*send_tgt_lnum)[idx] = tmp_send_tgt_lnum[i][j];
//       idx++;
//     }

//     if (i != myRank) {
//       free (tmp_send_tgt_lnum[i]);
//     }
//   }
//   free (tmp_send_tgt_lnum);


//   /* Manage send_start_leaves buffer */
//   *send_start_leaves_rank_shift = malloc (sizeof(int) * (lComm+1));
//   *send_start_leaves_count = malloc (sizeof(int) * send_shift[lComm]);
//   (*send_start_leaves_rank_shift)[0] = 0;
//   idx = 0;
//   for (int i = 0; i < lComm; i++) {
//     (*send_start_leaves_rank_shift)[i+1] = (*send_start_leaves_rank_shift)[i] + n_tmp_send_start_leaves[i];
//     for (int j = 0; j < send_count[i]; j++) {
//       (*send_start_leaves_count)[idx++] = tmp_send_tgt_n_leaves[i][j];
//     }
//   }

//   *send_start_leaves = malloc (sizeof(int) * (*send_start_leaves_rank_shift)[lComm]);
//   idx = 0;
//   for (int i = 0; i < lComm; i++) {
//     for (int j = 0; j < n_tmp_send_start_leaves[i]; j++) {
//       (*send_start_leaves)[idx++] = tmp_send_start_leaves[i][j];
//     }
//   }



//   for (int i = 0; i < lComm; i++) {
//     if (i != myRank) {
//       free (tmp_send_start_leaves[i]);
//       free (tmp_send_tgt_n_leaves[i]);
//     }
//   }
//   free (tmp_send_tgt_n_leaves);
//   free (s_tmp_send_start_leaves);
//   free (n_tmp_send_start_leaves);
//   free (tmp_send_start_leaves);

//   //-->> DETAIL TIMERS
//   PDM_timer_hang_on(octree->timer);
//   e_timer = PDM_timer_elapsed(octree->timer);
//   timer_ls[LS_FINALIZE] = e_timer - b_timer;
//   timer_ls[LS_END] = e_timer;
//   PDM_timer_resume(octree->timer);

//   _dump_timer_knnls (timer_ls);
//   //<<--
// }

/*============================================================================
 * CUDA Kernels
 *============================================================================*/

/**
 *
 * Force target points inside octree extents
 *
 * \param [in]      n_pts            Number of points
 * \param [in]      dim              Octree dimension
 * \param [inout]   pts              Point Coordinates
 * \param [in]      octree           Octree structure
 *
 */

__global__
static
void
_force_target_points
(
 int n_pts,
 int dim,
 double *pts,
 const _octree_t *octree
 )
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if ((i >= n_pts) || (j >= dim))
  {
    return;
  }

  pts[dim*i+j] = PDM_MAX (pts[dim*i+j], octree->global_extents[j]);
  pts[dim*i+j] = PDM_MIN (pts[dim*i+j], octree->global_extents[dim+j]);

}

/**
 *
 * Vector initialization
 *
 * \param [in]      n_pts_block1 
 * \param [in]      n_closest_points 
 * \param [inout]   block_closest_src_gnum1 
 * \param [inout]   block_closest_src_dist1 
 *
 */

__global__
static
void
_closest_src_init
(
  int n_pts_block1,
  int n_closest_points,
  PDM_g_num_t *block_closest_src_gnum1,
  double *block_closest_src_dist1
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  if ((i >= n_pts_block1) || (j >= n_closest_points))
  {
    return;
  }


  block_closest_src_dist1[n_closest_points*i + j] = HUGE_VAL;
  block_closest_src_gnum1[n_closest_points*i + j] = -1; // USEFUL ONLY FOR DEBUG
}

/**
 *
 * Perform Morton code binary search
 *
 * \param [in]      n_pts
 * \param [inout]   rank_pt 
 * \param [in]      lComm 
 * \param [in]      pts_code 
 * \param [in]      octree
 * \param [inout]   send_count
 *
 */

__global__
static
void
_binary_search
(
  int n_pts,
  int *rank_pt,
  int lComm,
  PDM_morton_code_t *pts_code,
  const _octree_t *octree,
  int *send_count
)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= n_pts)
  {
    return;
  }

  int k = PDM_morton_binary_search_GPU (lComm,
                                        pts_code[i],
                                        octree->rank_octants_index);

  send_count[i + k*n_pts] = 1;
  rank_pt[i] = k;
}



/**
 *
 * Get the k-th row from d_send_count_tab, and copy it in d_temp
 *
 * \param [in]      d_send_count_tab
 * \param [out]      d-temp
 * \param [in]   k
 * \param [in]   length
 *
 */

__global__
static
void
_get_row
(
  int *d_send_count_tab,
  int *d_temp, 
  int k,
  int length
  )
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= length)
  {
    return;
  }

  d_temp[i] = d_send_count_tab[i + k*length];
}



/*=============================================================================
 * Public function definitions
 *============================================================================*/


// /**
//  *
//  * \brief Create an octree structure on a GPU
//  *
//  * \param [in]   n_point_cloud          Number of point cloud
//  * \param [in]   depth_max              Maximum depth
//  * \param [in]   points_in_leaf_max     Maximum points in a leaf
//  * \param [in]   build_leaf_neighbours  Build leaf nieghbours (1 = true)
//  * \param [in]   comm                   MPI communicator
//  *
//  * \return     Identifier
//  */


// int
// PDM_para_octree_create_GPU
// (
//  const int n_point_cloud,
//  const int depth_max,
//  const int points_in_leaf_max,
//  const int build_leaf_neighbours,
//  const PDM_MPI_Comm comm
//  )
// {
//   if (_octrees == NULL) {
//     _octrees = PDM_Handles_create_GPU (4);
//   }

//   printf("octrees size %d\n", _octrees->s_array);

//   _octree_t *octree = (_octree_t *) malloc(sizeof(_octree_t));

//   int id = PDM_Handles_store_GPU (_octrees, octree);

//   octree->dim = 3;
//   _octree_t *octree_test = _get_from_id2(id);
//   printf("octree test dim %d\n", octree_test->dim);


//   for (int i = 0; i < octree->dim; i++) {
//     octree->global_extents[i]   = -HUGE_VAL;
//     octree->global_extents[octree->dim+i] =  HUGE_VAL;
//     octree->s[i]         = 0.;
//     octree->d[i]         = 0.;
//   }

//   octree->depth_max = depth_max;
//   octree->points_in_leaf_max = points_in_leaf_max;

//   octree->n_point_clouds = n_point_cloud;
//   octree->t_n_points = 0;
//   octree->n_points = 0;
//   octree->points = NULL;
//   octree->points_icloud = NULL;
//   octree->points_gnum = NULL;
//   octree->points_code = NULL;

//   octree->rank_octants_index = NULL;
//   octree->octants = NULL;

//   octree->n_part_boundary_elt = 0;
//   octree->part_boundary_elt_idx = NULL;
//   octree->part_boundary_elt = NULL;

//   octree->neighboursToBuild = build_leaf_neighbours;

//   octree->comm = comm;

//   octree->n_connected = 0;
//   octree->connected_idx = NULL;

//   octree->timer = PDM_timer_create ();

//   for (int i = 0; i < NTIMER2; i++) {
//     octree->times_elapsed[i] = 0.;
//     octree->times_cpu[i] = 0.;
//     octree->times_cpu_u[i] = 0.;
//     octree->times_cpu_s[i] = 0.;
//   }

//   PDM_printf("Depth max from cpu : %d\n", octree->depth_max);

//   return id;

// }


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
PDM_para_octree_closest_point_GPU
(
 PDM_octree_t    *octree,
 const int    n_closest_points,
 const int    n_pts,
 double      *pts,
 PDM_g_num_t *pts_g_num,
 PDM_g_num_t *closest_octree_pt_g_num,
 double      *closest_octree_pt_dist2
 )
{
  //-->> DETAIL TIMERS
  double timer_knn[N_TIMER_KNN];
  for (int i = 0; i < N_TIMER_KNN; i++) {
    timer_knn[i] = 0;
  }
  double b_timer, e_timer, start_main;
  //<<--

  const int DEBUG = 0;
  const int DEBUG_FILTER = 0;
  const int DEBUG_MERGE = 0;

  const int COMPUTE_FIRST_UPPER_BOUND = 1;

  //printf("before get from id octree\n");
  //const _octree_t *octree = (_octree_t*)octrees->array[0];
  const _l_octant_t *octants = octree->octants;

  const int dim = octree->dim;
  int _n_closest_points = n_closest_points;
  int _n_pts = n_pts;

  int myRank, lComm;
  PDM_MPI_Comm_rank (octree->comm, &myRank);
  PDM_MPI_Comm_size (octree->comm, &lComm);

  //-->> DETAIL TIMERS
  PDM_timer_hang_on(octree->timer);
  timer_knn[KNN_BEGIN] = PDM_timer_elapsed(octree->timer);
  b_timer = timer_knn[KNN_BEGIN];
  PDM_timer_resume(octree->timer);
  //<<--

  //Preparing octree copy on GPU - allocation on GPU
  PDM_octree_t *d_octree = NULL;
  gpuErrchk(cudaMalloc(&d_octree, sizeof(PDM_octree_t)));
  double *d_pts = NULL;
  gpuErrchk(cudaMalloc(&d_pts, 3*n_pts*sizeof(double)));


  //copy octree data on GPU
  gpuErrchk(cudaMemcpy(d_octree, octree, sizeof(PDM_octree_t), cudaMemcpyHostToDevice));
  gpuErrchk(cudaMemcpy(d_pts, pts, 3*n_pts*sizeof(double), cudaMemcpyHostToDevice));

  /* /!\ /!\ /!\ Force target points inside octree extents /!\ /!\ /!\ -->> */
  dim3 n_threads(32, 32);
  dim3 n_blocks((n_pts + 31)/32, (dim + 31)/32);
  _force_target_points<<<n_blocks,n_threads>>>(n_pts, dim, d_pts, d_octree);



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

  PDM_g_num_t *block_closest_src_gnum1 = (PDM_g_num_t*)malloc (sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points);
  double      *block_closest_src_dist1 = (double*)malloc (sizeof(double)      * n_pts_block1 * n_closest_points);
  PDM_g_num_t *d_block_closest_src_gnum1 = NULL;
  gpuErrchk(cudaMalloc(&d_block_closest_src_gnum1, sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points));
  double *d_block_closest_src_dist1 = NULL;
  gpuErrchk(cudaMalloc(&d_block_closest_src_dist1, sizeof(double) * n_pts_block1 * n_closest_points));

  n_threads = {32, 32, 1};
  n_blocks = {(n_pts_block1 + 31)/32, (n_closest_points + 31)/32, 1};
  _closest_src_init<<<n_blocks,n_threads>>>(n_pts_block1, n_closest_points, d_block_closest_src_gnum1, d_block_closest_src_dist1);


  gpuErrchk(cudaMemcpy(block_closest_src_gnum1, d_block_closest_src_gnum1, sizeof(PDM_g_num_t) * n_pts_block1 * n_closest_points, cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(block_closest_src_dist1, d_block_closest_src_dist1, sizeof(double) * n_pts_block1 * n_closest_points, cudaMemcpyDeviceToHost));

  gpuErrchk(cudaFree(d_block_closest_src_gnum1));
  gpuErrchk(cudaFree(d_block_closest_src_dist1));

  /*************************************************************************
   *
   * Distribute the target points
   *
   *************************************************************************/
  /*   1) Encode the coordinates of every target point */
  PDM_morton_code_t *pts_code = (PDM_morton_code_t*)malloc (sizeof(PDM_morton_code_t) * n_pts);
  PDM_morton_code_t *d_pts_code = NULL;
  gpuErrchk(cudaMalloc(&d_pts_code, sizeof(PDM_morton_code_t) * n_pts));
  double d[3], s[3];
  PDM_morton_encode_coords_GPU (dim,
                                n_pts,
                                PDM_morton_max_level,
                                octree->global_extents,
                                (size_t) n_pts,
                                d_pts,
                                d_pts_code,
                                d,
                                s);
                          

  gpuErrchk(cudaMemcpy(pts, d_pts, 3 * n_pts * sizeof(double), cudaMemcpyDeviceToHost));

  gpuErrchk(cudaFree(d_pts));                        
                           

  /*   2) Use binary search to associate each target point to the appropriate process */
  int *send_count = (int*)malloc (sizeof(int) * lComm);
  int *recv_count = (int*)malloc (sizeof(int) * lComm);
  for (int i = 0; i < lComm; i++) {
    send_count[i] = 0;
  }

  int *rank_pt = (int*)malloc (sizeof(int) * n_pts);
  int *d_rank_pt = NULL;
  gpuErrchk(cudaMalloc(&d_rank_pt, sizeof(int) * n_pts));

  n_threads = {1024, 1, 1};
  n_blocks = {(n_pts + 1023)/1024, 1, 1};
  
  int *d_send_count_tab = NULL;
  gpuErrchk(cudaMalloc(&d_send_count_tab, sizeof(int) * n_pts * lComm));

  // cudaEvent_t start, stop;
  // gpuErrchk(cudaEventCreate(&start));
  // gpuErrchk(cudaEventCreate(&stop));

  // gpuErrchk(cudaEventRecord(start));
  _binary_search<<<n_blocks,n_threads>>>(n_pts, d_rank_pt, lComm, d_pts_code, d_octree, d_send_count_tab);
  //The output of this kernel is an array of size lComm*total_threads, values are 0 or 1
  //There is one row per thread (and one thread per point), the proc where we will send the point is the index of the 1 value
  //Concurrent writing is impossible, so we will reduce this array to get the send_count array by adding the values of each column (column index are proc rank)
  // gpuErrchk(cudaEventRecord(stop));
  // gpuErrchk(cudaEventSynchronize(stop));

  // float milliseconds = 0;
  // gpuErrchk(cudaEventElapsedTime(&milliseconds, start, stop));
  // printf("n points : %d\n", n_pts);
  // printf("kernel time : %f ms\n", milliseconds);


  gpuErrchk(cudaMemcpy(rank_pt, d_rank_pt, sizeof(int) * n_pts, cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(octree, d_octree, sizeof(PDM_octree_t), cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(pts_code, d_pts_code, sizeof(PDM_morton_code_t) * n_pts, cudaMemcpyDeviceToHost));


  gpuErrchk(cudaFree(d_octree));
  gpuErrchk(cudaFree(d_rank_pt));
  gpuErrchk(cudaFree(d_pts_code));


  //Reduce send_count tab from each thread block to one tab
  int *d_temp = NULL;
  gpuErrchk(cudaMalloc(&d_temp, sizeof(int) * n_pts));
  for (int k = 0; k < lComm; k++)
  {
    n_threads = {1024, 1, 1};
    n_blocks = {(n_pts + 1023)/1024, 1, 1};
 
    //Get the column (or row) for the rank k proc
    _get_row<<<n_blocks,n_threads>>>(d_send_count_tab, d_temp, k, n_pts);
    //Reduce the row once using d_temp
    int blockSize = n_threads.x*n_threads.y*n_threads.z;
    int shared_size = blockSize*sizeof(int);
    printf("blocksize : %d, shared size : %d\n", blockSize, shared_size);
    Reduce_kernel(n_threads, n_blocks, shared_size, blockSize, n_pts, d_temp, d_temp);
    cudaDeviceSynchronize();


    //When we need to use multiple blocks of threads to reduce, the output is an array of size n_block
    //containing the reduction for each block
    //To get the complete reduction, we add the results of each block
    if (n_blocks.x > 1)
    {
      int n_temp = n_blocks.x;
      n_blocks = {(n_blocks.x + 1023)/1024, 1, 1};
      Reduce_kernel(n_threads, n_blocks, shared_size, blockSize, n_temp, d_temp, d_temp);
      cudaDeviceSynchronize();
    }

    //Copy result on the CPU, in the send_count array
    int *send_count_k = (int*)malloc(sizeof(int) * n_pts);
    gpuErrchk(cudaMemcpy(send_count_k, d_temp, sizeof(int), cudaMemcpyDeviceToHost));
    send_count[k] = send_count_k[0];
    printf("send_count[%d] = %d\n", k, send_count[k]);
  }

  gpuErrchk(cudaFree(d_temp));
  gpuErrchk(cudaFree(d_send_count_tab));

  
  free (pts_code);

  // /*   3) Exchange send/recv counts */
  // PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
  //                   recv_count, 1, PDM_MPI_INT,
  //                   octree->comm);

  // int *send_shift = (int *)malloc (sizeof(int) * (lComm+1));
  // int *recv_shift = (int *)malloc (sizeof(int) * (lComm+1));
  // send_shift[0] = 0;
  // recv_shift[0] = 0;
  // for (int i = 0; i < lComm; i++) {
  //   send_shift[i+1] = send_shift[i] + send_count[i];
  //   recv_shift[i+1] = recv_shift[i] + recv_count[i];
  //   send_count[i] = 0;
  // }

  // /*   4) Fill send buffers */
  // PDM_g_num_t *send_g_num = (PDM_g_num_t *)malloc (sizeof(PDM_g_num_t) * send_shift[lComm]);
  // PDM_g_num_t *recv_g_num = (PDM_g_num_t *)malloc (sizeof(PDM_g_num_t) * recv_shift[lComm]);
  // double      *send_coord = (double *)malloc (sizeof(double)      * send_shift[lComm]*dim);
  // double      *recv_coord = (double *)malloc (sizeof(double)      * recv_shift[lComm]*dim);


  // for (int i = 0; i < n_pts; i++) {
  //   int rank = rank_pt[i];
  //   int k = send_shift[rank] + send_count[rank];
  //   send_g_num[k] = pts_g_num[i];
  //   for (int j = 0; j < dim; j++)
  //     send_coord[dim*k+j] = pts[dim*i+j];

  //   send_count[rank]++;
  // }

  // free (rank_pt);

  // /*   5) Send gnum buffer */
  // PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
  //                    recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
  //                    octree->comm);

  // /*   6) Send coord buffer */
  // int n_recv_pts = recv_shift[lComm];
  // for (int i = 0; i < lComm; i++) {
  //   send_count[i] *= dim;
  //   recv_count[i] *= dim;
  //   send_shift[i+1] *= dim;
  //   recv_shift[i+1] *= dim;
  // }
  // PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
  //                    recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
  //                    octree->comm);




  // PDM_g_num_t *local_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  // double      *local_closest_src_dist = malloc (sizeof(double)      * n_recv_pts * n_closest_points);
  // double      *upper_bound_dist       = malloc (sizeof(double)      * n_recv_pts);

  // /* Encode the coordinates of the received target points */
  // pts_code = malloc (sizeof(PDM_morton_code_t) * n_recv_pts);
  // PDM_morton_encode_coords (dim,
  //                           PDM_morton_max_level,
  //                           octree->global_extents,
  //                           (size_t) n_recv_pts,
  //                           recv_coord,
  //                           pts_code,
  //                           d,
  //                           s);

  // if (COMPUTE_FIRST_UPPER_BOUND && octree->n_points >= n_closest_points) {
  //   const double EPS_max_dist = 1.e-6;
  //   const int window_width = (int) ceil (0.5 * n_closest_points);
  //   int window_start, window_end;
  //   /* Inspect a window of src points around each tgt point on the Z-order curve */
  //   for (int i = 0; i < n_recv_pts; i++) {
  //     int pos = PDM_morton_binary_search (octree->n_points,
  //                                         pts_code[i],
  //                                         octree->points_code);

  //     if (pos < window_width) {
  //       window_start = 0;
  //       window_end   = n_closest_points;
  //     } else if (pos >= octree->n_points - window_width) {
  //       window_end   = octree->n_points;
  //       window_start = window_end - n_closest_points;
  //     } else {
  //       window_start = pos - window_width;
  //       window_end   = window_start + n_closest_points;
  //     }

  //     const double *_pt = recv_coord + i * dim;
  //     double max_dist = 0.;
  //     for (int i_src = window_start; i_src < window_end; i_src++) {
  //       double *src_pt = octree->points + i_src * dim;
  //       double src_dist = 0;
  //       for (int j = 0; j < dim; j++) {
  //         double delta = _pt[j] - src_pt[j];
  //         src_dist += delta * delta;
  //       }

  //       max_dist = PDM_MAX (max_dist, src_dist);
  //     }

  //     upper_bound_dist[i] = max_dist + EPS_max_dist * max_dist;
  //   }

  // } else {

  //   for (int i = 0; i < n_recv_pts; i++) {
  //     upper_bound_dist[i] = HUGE_VAL;
  //   }

  // }


  // /* Find start leaf for each target point */
  // int *start_leaves = malloc (sizeof(int) * n_recv_pts);
  // for (int i = 0; i < n_recv_pts; i++) {
  //   start_leaves[i] = PDM_morton_binary_search (octants->n_nodes,
  //                                               pts_code[i],
  //                                               octants->codes);
  // }
  // free (pts_code);

  // int *start_leaves_idx = malloc (sizeof(int) * (n_recv_pts+1));
  // start_leaves_idx[0] = 0;
  // for (int i = 0; i < n_recv_pts; i++) {
  //   start_leaves_idx[i+1] = start_leaves_idx[i] + 1;
  // }


  // /* Stuff used for filering data coming from concurrent ranks */
  // /*PDM_hash_tab_t *processed_tgt = malloc (sizeof(PDM_hash_tab_t) * octree->n_connected);
  //   PDM_g_num_t keyMax = n_recv_pts; //?
  //   for (int i_part = 0; i_part < octree->n_connected; i_part++) {
  //   processed_tgt[i_part] = PDM_hash_tab_create (PDM_HASH_TAB_KEY_LONG
  //   &keyMax);
  //   }*/
  // int *s_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  // int *n_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  // PDM_g_num_t **processed_tgt = malloc (sizeof(PDM_g_num_t*) * octree->n_connected);
  // for (int i = 0; i < octree->n_connected; i++) {
  //   s_processed_tgt[i] = PDM_MAX (128, 2 * n_recv_pts);//?
  //   n_processed_tgt[i] = 0;
  //   processed_tgt[i] = malloc (sizeof(PDM_g_num_t) * s_processed_tgt[i]);
  // }

  // PDM_g_num_t **new_processed_tgt   = malloc (sizeof(PDM_g_num_t *) * octree->n_connected);
  // int          *n_new_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  // int          *s_new_processed_tgt = malloc (sizeof(int) * octree->n_connected);
  // for (int i = 0; i < octree->n_connected; i++) {
  //   s_new_processed_tgt[i] = PDM_MAX (128, n_recv_pts);//?
  //   new_processed_tgt[i] = malloc (sizeof(PDM_g_num_t) * s_new_processed_tgt[i]);
  // }


  // /* Stuff used to merge results received from 'part-to-block'2 */
  // _min_heap_t *merge_heap = _min_heap_create (10 * n_closest_points); // used to merge part-to-block results

  // double      *tmp_closest_src_dist = malloc (sizeof(double)      * n_closest_points);
  // PDM_g_num_t *tmp_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);

  // int *block_stride2 = NULL;
  // double      *block_closest_src_dist2 = NULL;
  // PDM_g_num_t *block_closest_src_gnum2 = NULL;
  // double      *block_upper_bound_dist = malloc (sizeof(double) * n_pts_block1);
  // int one = 1;


  // /* Stuff used to redistribute target points */
  // int *send_tgt_lnum = NULL;
  // int *send_start_leaves = NULL;
  // int *send_start_leaves_count = NULL;
  // int *send_start_leaves_rank_shift = NULL;

  // //-->> DETAIL TIMERS
  // PDM_timer_hang_on(octree->timer);
  // e_timer = PDM_timer_elapsed(octree->timer);
  // timer_knn[KNN_INIT] = e_timer - b_timer;
  // start_main = e_timer;
  // PDM_timer_resume(octree->timer);
  // //<<--

  // // while loop...
  // int iteration = 0;
  // while (1) {

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   timer_knn[LOOP_BEGIN] = PDM_timer_elapsed(octree->timer);
  //   b_timer = timer_knn[LOOP_BEGIN];
  //   PDM_timer_resume(octree->timer);
  //   //<<--

  //   //-->>
  //   iteration++;
  //   if (DEBUG) {
  //     printf("\n\n\n[%d] iteration %d\n", myRank, iteration);
  //   } else {
  //     if (myRank == 0)
  //       printf("\n\n\niteration %d\n", iteration);
  //   }
  //   if (iteration > 10*lComm) break; // emergency exit
  //   //<<--

  //   /* Filter 'recv' data: */
  //   //--->>>
  //   if (DEBUG) {
  //     printf("\n\n++++ BEFORE FILTERING ++++\n");
  //     printf("[%d] processed_pts:\n", myRank);
  //     for (int i = 0; i < octree->n_connected; i++) {
  //       printf("\tpart %d:", i);
  //       for (int j = 0; j < n_processed_tgt[i]; j++) {
  //         printf(" %ld", processed_tgt[i][j]);
  //       }
  //       printf("\n");
  //     }


  //     printf("\n[%d] recv_g_num = [\n", myRank);
  //     for (int i = 0; i < lComm; i++) {
  //       printf("\t{rank %d:", i);
  //       for (int j = recv_shift[i]/dim; j < recv_shift[i+1]/dim; j++) {
  //         printf(" %ld", recv_g_num[j]);
  //       }
  //       printf(" }\n");
  //     }
  //     printf(" ]\n");


  //     printf("\n[%d] start_leaves = [\n", myRank);
  //     for (int i = 0; i < n_recv_pts; i++) {
  //       printf("  (%ld) : {", recv_g_num[i]);
  //       for (int j = start_leaves_idx[i]; j < start_leaves_idx[i+1]; j++) {
  //         printf(" %d", start_leaves[j]);
  //       }
  //       printf("}\n");
  //     }
  //     printf(" ]\n\n\n");
  //   }
  //   //<<<---

  //   if (iteration == 1) {
  //     // at this point, there is exactly one start leaf per tgt point
  //     for (int i = 0; i < n_recv_pts; i++) {
  //       // get part number
  //       int i_part = _get_octant_part_id (octree,
  //                                         start_leaves[i]);
  //       processed_tgt[i_part][n_processed_tgt[i_part]++] = recv_g_num[i];
  //     }

  //     // sort in ascending order
  //     for (int i_part = 0; i_part < octree->n_connected; i_part++) {
  //       PDM_sort_long (processed_tgt[i_part],
  //                      NULL,
  //                      n_processed_tgt[i_part]);
  //     }

  //   } else {
  //     // at this point, there are possibly
  //     //   1) duplicate tgt points (but with different start leaves)
  //     //   2) tgt points already processed in some parts of current rank
  //     int tmp_n = 0;
  //     PDM_g_num_t *tmp_g_num       = malloc (sizeof(PDM_g_num_t) * n_recv_pts);
  //     double      *tmp_coord       = malloc (sizeof(double)      * n_recv_pts * dim);
  //     double      *tmp_upper_bound = malloc (sizeof(double)      * n_recv_pts);

  //     int  *n_tmp_start_leaves = malloc (sizeof(int) * n_recv_pts);
  //     int  *s_tmp_start_leaves = malloc (sizeof(int) * n_recv_pts);
  //     int **tmp_start_leaves   = malloc (sizeof(int *) * n_recv_pts);


  //     for (int i_part = 0; i_part < octree->n_connected; i_part++) {
  //       n_new_processed_tgt[i_part] = 0;
  //     }

  //     for (int i_tgt = 0; i_tgt < n_recv_pts; i_tgt++) {
  //       PDM_g_num_t tgt_gnum = recv_g_num[i_tgt];
  //       if (DEBUG && DEBUG_FILTER) {
  //         printf("\nFilter: recv pt (%ld)\n", tgt_gnum);
  //       }

  //       // check whether that tgt point has already been added to tmp_* arrays
  //       int found_tmp = 0;
  //       int pos_tmp = _binary_search_long (tgt_gnum,
  //                                          tmp_g_num,
  //                                          tmp_n,
  //                                          &found_tmp);
  //       if (DEBUG && DEBUG_FILTER) {
  //         printf(" found_tmp = %d\tpos_tmp = %d\n", found_tmp, pos_tmp);
  //       }

  //       if (!found_tmp) {
  //         pos_tmp = tmp_n;
  //         // add tgt point to tmp_* arrays
  //         tmp_g_num[pos_tmp] = tgt_gnum;
  //         for (int i = 0; i < dim; i++) {
  //           tmp_coord[dim*pos_tmp + i] = recv_coord[dim*i_tgt + i];
  //         }

  //         tmp_upper_bound[pos_tmp] = upper_bound_dist[i_tgt];

  //         s_tmp_start_leaves[pos_tmp] = 2 * (start_leaves_idx[i_tgt+1] - start_leaves_idx[i_tgt]);
  //         tmp_start_leaves[pos_tmp] = malloc (sizeof(int) * s_tmp_start_leaves[pos_tmp]);
  //         n_tmp_start_leaves[pos_tmp] = 0;

  //         tmp_n++;
  //       }

  //       for (int i_start = start_leaves_idx[i_tgt];
  //            i_start < start_leaves_idx[i_tgt+1]; i_start++) {
  //         int leaf_id = start_leaves[i_start];
  //         int leaf_part = _get_octant_part_id (octree,
  //                                              leaf_id);
  //         if (DEBUG && DEBUG_FILTER) {
  //           printf("\tstart leaf id %d (part %d)\n", leaf_id, leaf_part);
  //         }

  //         // check whether that tgt point has already been processed in part #leaf_part
  //         int found = 0;
  //         _binary_search_long (tgt_gnum,
  //                              processed_tgt[leaf_part],
  //                              n_processed_tgt[leaf_part],
  //                              &found);
  //         if (DEBUG && DEBUG_FILTER) {
  //           printf("\t already processed? %d\n", found);
  //         }

  //         if (!found) {
  //           // check whether that start leaf has already been added to tmp_start_leaves[leaf_part]
  //           int found_leaf = 0;
  //           int pos_leaf = _binary_search (leaf_id,
  //                                          tmp_start_leaves[pos_tmp],
  //                                          n_tmp_start_leaves[pos_tmp],
  //                                          &found_leaf);
  //           if (DEBUG && DEBUG_FILTER) {
  //             printf("\t\tfound_leaf = %d, pos_leaf = %d\n", found_leaf, pos_leaf);
  //           }

  //           if (!found_leaf) {
  //             /* add start leaf to tmp_start_leaves[pos_tmp] */
  //             // realloc tmp_start_leaves[pos_tmp] if necessary
  //             if (s_tmp_start_leaves[pos_tmp] <= n_tmp_start_leaves[pos_tmp]) {
  //               s_tmp_start_leaves[pos_tmp] *= 2;
  //               tmp_start_leaves[pos_tmp] = realloc (tmp_start_leaves[pos_tmp],
  //                                                    sizeof(int) * s_tmp_start_leaves[pos_tmp]);
  //             }

  //             // insert-sort leaf_id in tmp_start_leaves[pos_tmp]
  //             for (int i = n_tmp_start_leaves[pos_tmp]; i > pos_leaf; i--) {
  //               tmp_start_leaves[pos_tmp][i] = tmp_start_leaves[pos_tmp][i-1];
  //             }
  //             tmp_start_leaves[pos_tmp][pos_leaf] = leaf_id;
  //             n_tmp_start_leaves[pos_tmp]++;



  //             // check whether that tgt point has already been added to new_processed_tgt[leaf_part]
  //             int found_new = 0;
  //             int pos_new = _binary_search_long (tgt_gnum,
  //                                                new_processed_tgt[leaf_part],
  //                                                n_new_processed_tgt[leaf_part],
  //                                                &found_new);
  //             if (DEBUG && DEBUG_FILTER) {
  //               printf("\t\tfound_new = %d, pos_new = %d\n", found_new, pos_new);
  //             }

  //             if (!found_new) {
  //               /* add tgt point to new_processed_tgt[leaf_part] */
  //               // realloc new_processed_tgt[leaf_part] if necessary
  //               if (s_new_processed_tgt[leaf_part] <= n_new_processed_tgt[leaf_part]) {
  //                 s_new_processed_tgt[leaf_part] *= 2;
  //                 new_processed_tgt[leaf_part] = realloc (new_processed_tgt[leaf_part],
  //                                                         sizeof(PDM_g_num_t) * s_new_processed_tgt[leaf_part]);
  //               }

  //               // insert-sort tgt_gnum in new_processed_tgt[leaf_part]
  //               for (int i = n_new_processed_tgt[leaf_part]; i > pos_new; i--) {
  //                 new_processed_tgt[leaf_part][i] = new_processed_tgt[leaf_part][i-1];
  //               }
  //               new_processed_tgt[leaf_part][pos_new] = tgt_gnum;
  //               n_new_processed_tgt[leaf_part]++;
  //             }

  //           } // end if (leaf_id not found in processed_tgt[leaf_part])
  //         } // end if (tgt_gnum not found in processed_tgt[leaf_part])
  //       } // end loop over start leaves (i_start)
  //     } // end loop over received tgt points (i_tgt)

  //     int k = 0;
  //     start_leaves_idx[0] = 0;
  //     for (int i = 0; i < tmp_n; i++) {
  //       if (n_tmp_start_leaves[i] > 0) {
  //         start_leaves_idx[k+1] = start_leaves_idx[k] + n_tmp_start_leaves[i];

  //         recv_g_num[k] = tmp_g_num[i];

  //         for (int j = 0; j < dim; j++) {
  //           recv_coord[dim*k+j] = tmp_coord[dim*i+j];
  //         }

  //         upper_bound_dist[k] = tmp_upper_bound[i];
  //         k++;
  //       }
  //     }

  //     if (k < n_recv_pts) {
  //       recv_g_num       = realloc (recv_g_num,       sizeof(PDM_g_num_t) * k);
  //       recv_coord       = realloc (recv_coord,       sizeof(double)      * k * dim);
  //       upper_bound_dist = realloc (upper_bound_dist, sizeof(double)      * k);
  //     }
  //     free (tmp_g_num);
  //     free (tmp_coord);
  //     free (tmp_upper_bound);

  //     /* manage start leaves */
  //     start_leaves = realloc (start_leaves, sizeof(int) * start_leaves_idx[k]);
  //     k = 0;
  //     int idx = 0;
  //     for (int i = 0; i < tmp_n; i++) {
  //       if (n_tmp_start_leaves[i] > 0) {
  //         for (int j = 0; j < n_tmp_start_leaves[i]; j++) {
  //           start_leaves[idx++] = tmp_start_leaves[i][j];
  //         }
  //         k++;
  //       }
  //     }

  //     for (int i = 0; i < tmp_n; i++) {
  //       free (tmp_start_leaves[i]);
  //     }
  //     free (tmp_start_leaves);
  //     free (s_tmp_start_leaves);
  //     free (n_tmp_start_leaves);

  //     n_recv_pts = k;


  //     /* merge new_processed_tgt into processed_tgt */
  //     for (int i_part = 0; i_part < octree->n_connected; i_part++) {
  //       // realloc processed_tgt[i_part] if necessary
  //       if (s_processed_tgt[i_part] <= n_processed_tgt[i_part] + n_new_processed_tgt[i_part]) {
  //         s_processed_tgt[i_part] = PDM_MAX (2 * s_processed_tgt[i_part],
  //                                            n_processed_tgt[i_part] + n_new_processed_tgt[i_part]);
  //         processed_tgt[i_part] = realloc (processed_tgt[i_part],
  //                                          sizeof(PDM_g_num_t) * s_processed_tgt[i_part]);
  //       }

  //       for (int i = 0; i < n_new_processed_tgt[i_part]; i++) {
  //         int found = 0;
  //         int pos = _binary_search_long (new_processed_tgt[i_part][i],
  //                                        processed_tgt[i_part],
  //                                        n_processed_tgt[i_part],
  //                                        &found);

  //         // insert-sort
  //         for (int j = n_processed_tgt[i_part]; j > pos; j--) {
  //           processed_tgt[i_part][j] = processed_tgt[i_part][j-1];
  //         }
  //         processed_tgt[i_part][pos] = new_processed_tgt[i_part][i];
  //         n_processed_tgt[i_part]++;
  //       }
  //     }

  //   } // end if/else (iteration == 1)


  //   //--->>>
  //   if (DEBUG) {
  //     printf("++++ AFTER FILTERING ++++\n");
  //     printf("[%d] processed_pts:\n", myRank);
  //     for (int i = 0; i < octree->n_connected; i++) {
  //       printf("\tpart %d:", i);
  //       for (int j = 0; j < n_processed_tgt[i]; j++) {
  //         printf(" %ld", processed_tgt[i][j]);
  //       }
  //       printf("\n");
  //     }


  //     printf("\n[%d] recv_g_num = [", myRank);
  //     for (int i = 0; i < n_recv_pts; i++) {
  //       printf(" %ld", recv_g_num[i]);
  //     }
  //     printf(" ]\n");


  //     printf("\n[%d] start_leaves = [\n", myRank);
  //     for (int i = 0; i < n_recv_pts; i++) {
  //       printf("  (%ld) : {", recv_g_num[i]);
  //       for (int j = start_leaves_idx[i]; j < start_leaves_idx[i+1]; j++) {
  //         printf(" %d", start_leaves[j]);
  //       }
  //       printf("}\n");
  //     }
  //     printf(" ]\n\n\n");
  //   }
  //   //<<<---

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   e_timer = PDM_timer_elapsed(octree->timer);
  //   timer_knn[LOOP_FILTER] = e_timer - b_timer;

  //   b_timer = PDM_timer_elapsed(octree->timer);
  //   PDM_timer_resume(octree->timer);
  //   //<<--

  //   /* Search closest src points in local octree */
  //   _closest_points_local (octree,
  //                          n_closest_points,
  //                          n_recv_pts,
  //                          recv_coord,
  //                          recv_g_num, // ONLY FOR DEBUG
  //                          start_leaves,
  //                          start_leaves_idx,
  //                          upper_bound_dist,
  //                          local_closest_src_gnum,
  //                          local_closest_src_dist,
  //                          send_count,
  //                          send_shift,
  //                          &send_tgt_lnum,
  //                          &send_start_leaves,
  //                          &send_start_leaves_count,
  //                          &send_start_leaves_rank_shift);

  //   printf("Exit program\n");
  //   exit(1); 

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   e_timer = PDM_timer_elapsed(octree->timer);
  //   timer_knn[LOOP_LOCAL_SEARCH] = e_timer - b_timer;
  //   b_timer = PDM_timer_elapsed(octree->timer);
  //   PDM_timer_resume(octree->timer);
  //   //<<--

  //   /* Part-to-block exchanges to merge results in block arrays */
  //   PDM_part_to_block_t *ptb2 = PDM_part_to_block_create2 (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
  //                                                          PDM_PART_TO_BLOCK_POST_MERGE,
  //                                                          1.,
  //                                                          &recv_g_num,
  //                                                          block_distrib_idx1,
  //                                                          &n_recv_pts,
  //                                                          1,
  //                                                          octree->comm);

  //   int n_pts_block2 = PDM_part_to_block_n_elt_block_get (ptb2);
  //   PDM_g_num_t *block_tgt_gnum2 = PDM_part_to_block_block_gnum_get (ptb2);

  //   int *stride2 = malloc (sizeof(int) * n_recv_pts);
  //   for (int i = 0; i < n_recv_pts; i++) {
  //     stride2[i] = n_closest_points;
  //   }

  //   PDM_part_to_block_exch (ptb2,
  //                           sizeof(double),
  //                           PDM_STRIDE_VAR,
  //                           1,
  //                           &stride2,
  //                           (void **) &local_closest_src_dist,
  //                           &block_stride2,
  //                           (void **) &block_closest_src_dist2);
  //   free (block_stride2);

  //   PDM_part_to_block_exch (ptb2,
  //                           sizeof(PDM_g_num_t),
  //                           PDM_STRIDE_VAR,
  //                           1,
  //                           &stride2,
  //                           (void **) &local_closest_src_gnum,
  //                           &block_stride2,
  //                           (void **) &block_closest_src_gnum2);
  //   free (stride2);

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   e_timer = PDM_timer_elapsed(octree->timer);
  //   timer_knn[LOOP_PTB_EXCH] = e_timer - b_timer;
  //   b_timer = PDM_timer_elapsed(octree->timer);
  //   PDM_timer_resume(octree->timer);
  //   //<<--

  //   /* Merge block data */
  //   if (DEBUG && DEBUG_MERGE) {
  //     printf("\n\n- - - Merge - - -");
  //   }
  //   int *block_idx2 = malloc (sizeof(int) * (n_pts_block2 + 1));
  //   block_idx2[0] = 0;
  //   for (int i = 0; i < n_pts_block2; i++) {
  //     block_idx2[i+1] = block_idx2[i] + block_stride2[i];
  //   }

  //   for (int i = 0; i < n_pts_block2; i++) {
  //     int id_block1 = (int) n_closest_points * (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]);

  //     if (DEBUG && DEBUG_MERGE) {
  //       /*printf("\nblock2 i = %d (%ld) <%d>\n", i, block_tgt_gnum2[i],
  //         (int) (block_tgt_gnum2[i]-1 - block_distrib_idx1[myRank]));
  //         printf("id_block1 = %d / %d\n", id_block1, n_pts_block1*n_closest_points);*/
  //       printf("\npoint (%ld)\n", block_tgt_gnum2[i]);
  //     }

  //     int n_procs_pt = block_stride2[i] / n_closest_points;

  //     int *idx_proc = malloc (sizeof(int) * (n_procs_pt + 1));

  //     _min_heap_reset (merge_heap);

  //     for (int j = 0; j < n_closest_points; j++) {
  //       tmp_closest_src_gnum[j] = block_closest_src_gnum1[id_block1 + j];
  //       tmp_closest_src_dist[j] = block_closest_src_dist1[id_block1 + j];
  //     }

  //     if (DEBUG && DEBUG_MERGE) {
  //       printf("\tj = %d: (%ld, %f)\n",
  //              -1, tmp_closest_src_gnum[0], tmp_closest_src_dist[0]);
  //     }

  //     _min_heap_push (merge_heap,
  //                     -1,
  //                     tmp_closest_src_gnum[0],
  //                     tmp_closest_src_dist[0]);
  //     idx_proc[0] = 0;


  //     for (int j = 0; j < n_procs_pt; j++) {
  //       int k = block_idx2[i] + j*n_closest_points;

  //       if (DEBUG && DEBUG_MERGE) {
  //         printf("\tj =  %d: (%ld, %f)\n",
  //                j, block_closest_src_gnum2[k], block_closest_src_dist2[k]);
  //       }

  //       _min_heap_push (merge_heap,
  //                       j,
  //                       block_closest_src_gnum2[k],
  //                       block_closest_src_dist2[k]);
  //       idx_proc[j+1] = 0;
  //     }


  //     int _proc;
  //     PDM_g_num_t _gnum;
  //     double _dist;
  //     for (int j = 0; j < n_closest_points; j++) {
  //       int popped = _min_heap_pop (merge_heap,
  //                                   &_proc,
  //                                   &_gnum,
  //                                   &_dist);
  //       assert (popped);
  //       if (DEBUG && DEBUG_MERGE) {
  //         printf("\tpopped: %ld %f \t / %f\n", _gnum, _dist, block_closest_src_dist1[id_block1 + j]);
  //       }

  //       if (_dist >= block_closest_src_dist1[id_block1 + n_closest_points - 1]) break;

  //       block_closest_src_dist1[id_block1 + j] = _dist;
  //       block_closest_src_gnum1[id_block1 + j] = _gnum;

  //       if (DEBUG && DEBUG_MERGE) {
  //         printf("\t\t%d/%d: %ld, %f\n",
  //                j+1, n_closest_points,
  //                block_closest_src_gnum1[id_block1 + j], block_closest_src_dist1[id_block1 + j]);
  //       }

  //       if (j >= n_closest_points - 1)
  //         break;

  //       idx_proc[_proc+1]++;

  //       if (_proc < 0) {
  //         _min_heap_push (merge_heap,
  //                         -1,
  //                         tmp_closest_src_gnum[idx_proc[0]],
  //                         tmp_closest_src_dist[idx_proc[0]]);
  //       } else {
  //         int k = block_idx2[i] + _proc*n_closest_points + idx_proc[_proc+1];
  //         _min_heap_push (merge_heap,
  //                         _proc,
  //                         block_closest_src_gnum2[k],
  //                         block_closest_src_dist2[k]);
  //       }

  //     }
  //     free (idx_proc);
  //   }
  //   if (DEBUG && DEBUG_MERGE) {
  //     printf("- - - - - - - - -\n\n");
  //   }

  //   free (block_idx2);
  //   free (block_stride2);
  //   free (block_closest_src_gnum2);
  //   free (block_closest_src_dist2);
  //   ptb2 = PDM_part_to_block_free (ptb2);
  //   // end merge

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   e_timer = PDM_timer_elapsed(octree->timer);
  //   timer_knn[LOOP_MERGE_PTB] = e_timer - b_timer;
  //   b_timer = PDM_timer_elapsed(octree->timer);
  //   PDM_timer_resume(octree->timer);
  //   //<<--


  //   /* Update upper_bound_dist */
  //   PDM_block_to_part_t *btp2 = PDM_block_to_part_create (block_distrib_idx1,
  //                                                         (const PDM_g_num_t **) &recv_g_num,
  //                                                         &n_recv_pts,
  //                                                         1,
  //                                                         octree->comm);

  //   for (int i = 0; i < n_pts_block1; i++) {
  //     block_upper_bound_dist[i] = block_closest_src_dist1[n_closest_points*(i+1) - 1];
  //   }

  //   PDM_block_to_part_exch (btp2,
  //                           sizeof(double),
  //                           PDM_STRIDE_CST,
  //                           &one,
  //                           block_upper_bound_dist,
  //                           NULL,
  //                           (void **) &upper_bound_dist);
  //   btp2 = PDM_block_to_part_free (btp2);


  //   /* Redistribute target points for next iteration */
  //   // send count
  //   PDM_MPI_Alltoall (send_count, 1, PDM_MPI_INT,
  //                     recv_count, 1, PDM_MPI_INT,
  //                     octree->comm);
  //   recv_shift[0] = 0;
  //   for (int i = 0; i < lComm; i++)
  //     recv_shift[i+1] = recv_shift[i] + recv_count[i];

  //   n_recv_pts = recv_shift[lComm];


  //   /* Termination criterion */
  //   int max_n_recv_pts = 0;
  //   PDM_MPI_Allreduce (&n_recv_pts, &max_n_recv_pts, 1,
  //                      PDM_MPI_INT, PDM_MPI_MAX, octree->comm);

  //   if (max_n_recv_pts == 0) {
  //     if (myRank == 0)
  //       printf("*** max_n_recv_pts = 0 --> DONE :-)\n");

  //     free (send_tgt_lnum);
  //     free (send_start_leaves);
  //     free (send_start_leaves_count);
  //     free (send_start_leaves_rank_shift);

  //     //-->> DETAIL TIMERS
  //     PDM_timer_hang_on(octree->timer);
  //     e_timer = PDM_timer_elapsed(octree->timer);
  //     timer_knn[LOOP_PREP_NEXT] = e_timer - b_timer;
  //     timer_knn[LOOP_END] = e_timer;
  //     PDM_timer_resume(octree->timer);

  //     _dump_timer_knn_loop (iteration, timer_knn);
  //     //<<--
  //     break;
  //   }


  //   // send g_num
  //   send_g_num = realloc (send_g_num, sizeof(PDM_g_num_t) * send_shift[lComm]);
  //   for (int i = 0; i < send_shift[lComm]; i++) {
  //     send_g_num[i] = recv_g_num[send_tgt_lnum[i]];
  //   }
  //   recv_g_num = realloc (recv_g_num, sizeof(PDM_g_num_t) * recv_shift[lComm]);
  //   PDM_MPI_Alltoallv (send_g_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
  //                      recv_g_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
  //                      octree->comm);

  //   // send upper bound dist
  //   double *send_upper_bound_dist = malloc (sizeof(double) * send_shift[lComm]);
  //   for (int i = 0; i < send_shift[lComm]; i++) {
  //     send_upper_bound_dist[i] = upper_bound_dist[send_tgt_lnum[i]];
  //   }
  //   upper_bound_dist = realloc (upper_bound_dist, sizeof(double) * n_recv_pts);
  //   PDM_MPI_Alltoallv (send_upper_bound_dist, send_count, send_shift, PDM_MPI_DOUBLE,
  //                      upper_bound_dist,      recv_count, recv_shift, PDM_MPI_DOUBLE,
  //                      octree->comm);
  //   free (send_upper_bound_dist);

  //   // send start leaves
  //   int *send_start_leaves_rank_count = malloc (sizeof(int) * lComm);
  //   int *recv_start_leaves_rank_count = malloc (sizeof(int) * lComm);
  //   for (int i = 0; i < lComm; i++) {
  //     send_start_leaves_rank_count[i] =
  //       send_start_leaves_rank_shift[i+1] - send_start_leaves_rank_shift[i];
  //   }


  //   PDM_MPI_Alltoall (send_start_leaves_rank_count, 1, PDM_MPI_INT,
  //                     recv_start_leaves_rank_count, 1, PDM_MPI_INT,
  //                     octree->comm);

  //   int *recv_start_leaves_rank_shift = malloc (sizeof(int) * (lComm+1));
  //   recv_start_leaves_rank_shift[0] = 0;
  //   for (int i = 0; i < lComm; i++) {
  //     recv_start_leaves_rank_shift[i+1] =
  //       recv_start_leaves_rank_shift[i] + recv_start_leaves_rank_count[i];
  //   }

  //   start_leaves = realloc (start_leaves, sizeof(int) * recv_start_leaves_rank_shift[lComm]);
  //   PDM_MPI_Alltoallv (send_start_leaves,
  //                      send_start_leaves_rank_count,
  //                      send_start_leaves_rank_shift,
  //                      PDM_MPI_INT,
  //                      start_leaves,
  //                      recv_start_leaves_rank_count,
  //                      recv_start_leaves_rank_shift,
  //                      PDM_MPI_INT,
  //                      octree->comm);
  //   free (send_start_leaves);
  //   free (send_start_leaves_rank_shift);


  //   int *recv_start_leaves_count = malloc (sizeof(int) * recv_shift[lComm]);
  //   PDM_MPI_Alltoallv (send_start_leaves_count,
  //                      send_count,
  //                      send_shift,
  //                      PDM_MPI_INT,
  //                      recv_start_leaves_count,
  //                      recv_count,
  //                      recv_shift,
  //                      PDM_MPI_INT,
  //                      octree->comm);

  //   start_leaves_idx = realloc (start_leaves_idx, sizeof(int) * (recv_shift[lComm]+1));
  //   start_leaves_idx[0] = 0;
  //   for (int i = 0; i < recv_shift[lComm]; i++) {
  //     start_leaves_idx[i+1] = start_leaves_idx[i] + recv_start_leaves_count[i];
  //   }


  //   free (send_start_leaves_rank_count);
  //   free (recv_start_leaves_rank_count);
  //   free (recv_start_leaves_rank_shift);

  //   free (send_start_leaves_count);
  //   free (recv_start_leaves_count);


  //   // send coords
  //   send_coord = realloc (send_coord, sizeof(double) * send_shift[lComm] * dim);

  //   for (int i = 0; i < send_shift[lComm]; i++) {
  //     for (int j = 0; j < dim; j++)
  //       send_coord[dim*i + j] = recv_coord[dim*send_tgt_lnum[i] + j];
  //   }
  //   free (send_tgt_lnum);

  //   n_recv_pts = recv_shift[lComm];
  //   for (int i = 0; i < lComm; i++) {
  //     send_count[i]   *= dim;
  //     recv_count[i]   *= dim;
  //     send_shift[i+1] *= dim;
  //     recv_shift[i+1] *= dim;
  //   }
  //   recv_coord = realloc (recv_coord, sizeof(double) * recv_shift[lComm]);
  //   PDM_MPI_Alltoallv (send_coord, send_count, send_shift, PDM_MPI_DOUBLE,
  //                      recv_coord, recv_count, recv_shift, PDM_MPI_DOUBLE,
  //                      octree->comm);


  //   local_closest_src_gnum = realloc (local_closest_src_gnum, sizeof(PDM_g_num_t) * n_recv_pts * n_closest_points);
  //   local_closest_src_dist = realloc (local_closest_src_dist, sizeof(double)      * n_recv_pts * n_closest_points);

  //   //-->> DETAIL TIMERS
  //   PDM_timer_hang_on(octree->timer);
  //   e_timer = PDM_timer_elapsed(octree->timer);
  //   timer_knn[LOOP_PREP_NEXT] = e_timer - b_timer;
  //   timer_knn[LOOP_END] = e_timer;
  //   PDM_timer_resume(octree->timer);
  //   _dump_timer_knn_loop (iteration, timer_knn);
  //   //<<--
  // } // end while loop

  // //-->> DETAIL TIMERS
  // PDM_timer_hang_on(octree->timer);
  // e_timer = PDM_timer_elapsed(octree->timer);
  // timer_knn[KNN_MAIN] = e_timer - start_main;
  // b_timer = e_timer;
  // PDM_timer_resume(octree->timer);
  // //<<--

  // /* Free stuff */
  // free (start_leaves);
  // free (start_leaves_idx);
  // free (local_closest_src_gnum);
  // free (local_closest_src_dist);

  // free (send_coord);
  // free (send_g_num);
  // free (send_count);
  // free (send_shift);

  // free (recv_coord);
  // free (recv_g_num);
  // free (recv_count);
  // free (recv_shift);

  // _min_heap_free (merge_heap);
  // free (tmp_closest_src_dist);
  // free (tmp_closest_src_gnum);

  // free (block_upper_bound_dist);


  // for (int i = 0; i < octree->n_connected; i++) {
  //   free (processed_tgt[i]);
  //   free (new_processed_tgt[i]);
  // }
  // free (s_processed_tgt);
  // free (n_processed_tgt);
  // free (processed_tgt);
  // free (s_new_processed_tgt);
  // free (n_new_processed_tgt);
  // free (new_processed_tgt);


  // /* Final Block-to-part exchanges */
  // PDM_block_to_part_t *btp1 = PDM_block_to_part_create (block_distrib_idx1,
  //                                                       (const PDM_g_num_t **) &pts_g_num,
  //                                                       &n_pts,
  //                                                       1,
  //                                                       octree->comm);
  // PDM_part_to_block_free (ptb1);

  // PDM_block_to_part_exch (btp1,
  //                         sizeof(double),
  //                         PDM_STRIDE_CST,
  //                         &_n_closest_points,
  //                         block_closest_src_dist1,
  //                         NULL,
  //                         (void **) &closest_octree_pt_dist2);

  // PDM_block_to_part_exch (btp1,
  //                         sizeof(PDM_g_num_t),
  //                         PDM_STRIDE_CST,
  //                         &_n_closest_points,
  //                         block_closest_src_gnum1,
  //                         NULL,
  //                         (void **) &closest_octree_pt_g_num);

  // free (block_closest_src_dist1);
  // free (block_closest_src_gnum1);

  // btp1 = PDM_block_to_part_free (btp1);

  // //-->> DETAIL TIMERS
  // PDM_timer_hang_on(octree->timer);
  // e_timer = PDM_timer_elapsed(octree->timer);
  // timer_knn[KNN_BTP_EXCH] = e_timer - b_timer;
  // timer_knn[KNN_END] = e_timer;
  // PDM_timer_resume(octree->timer);

  // _dump_timer_knn (timer_knn);
  // //<<--
}



// /**
//  *
//  * \brief Free an octree structure
//  *
//  * \param [in]   id                 Identifier
//  *
//  */

// __host__
// void
// PDM_para_octree_free_GPU
// (
//   const int          id
//   )
// {
//   _octree_t *octree = _get_from_id2 (id);

//   if (octree->points != NULL) {
//     gpuErrchk(cudaFree (octree->points));
//   }

//   if (octree->points_icloud != NULL) {
//     gpuErrchk(cudaFree (octree->points_icloud));
//   }

//   if (octree->points_gnum != NULL) {
//     gpuErrchk(cudaFree (octree->points_gnum));
//   }

//   if (octree->points_code != NULL) {
//     gpuErrchk(cudaFree (octree->points_code));
//   }

//   if (octree->part_boundary_elt_idx != NULL) {
//     gpuErrchk(cudaFree (octree->part_boundary_elt_idx));
//   }

//   if (octree->part_boundary_elt != NULL) {
//     gpuErrchk(cudaFree (octree->part_boundary_elt));
//   }

//   if (octree->rank_octants_index != NULL) {
//     gpuErrchk(cudaFree (octree->rank_octants_index));
//   }

//   if (octree->octants != NULL) {

//     if (octree->octants->codes != NULL) {
//       gpuErrchk(cudaFree (octree->octants->codes));
//     }

//     if (octree->octants->range != NULL) {
//       gpuErrchk(cudaFree (octree->octants->range));
//     }

//     if (octree->octants->n_points != NULL) {
//       gpuErrchk(cudaFree (octree->octants->n_points));
//     }

//     if (octree->octants->neighbour_idx != NULL) {
//       gpuErrchk(cudaFree (octree->octants->neighbour_idx));
//     }

//     if (octree->octants->neighbours != NULL) {
//       gpuErrchk(cudaFree (octree->octants->neighbours));
//     }

//     gpuErrchk(cudaFree (octree->octants));
//   }

//   if (octree->connected_idx != NULL) {
//     gpuErrchk(cudaFree (octree->connected_idx));
//   }

//   PDM_timer_free_GPU (octree->timer);

//   gpuErrchk(cudaFree (octree));

//   PDM_Handles_handle_free (_octrees, id, PDM_FALSE);

//   const int n_octrees = PDM_Handles_n_get (_octrees);

//   if (n_octrees == 0) {
//     _octrees = PDM_Handles_free (_octrees);
//   }
// }






















//test
__global__
void
print_from_gpu
(
 PDM_octree_t *octree
 )
{
  printf("Hello World! from thread [%d,%d] From device\n", threadIdx.x,blockIdx.x);
  //const _octree_t *octree = (_octree_t*)octrees->array[0];
  printf("octree %p\n", octree);
  printf("Depth max from GPU : %d\n", octree->depth_max);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
