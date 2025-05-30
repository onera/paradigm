/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

/*============================================================================
 * Main structure for an I/O numbering scheme associated with mesh entities
 * (such as cells, faces, and vertices);
 *
 * In parallel mode, such a scheme is important so as to redistribute
 * locally numbered entities on n processes to files written by p
 * processes, with p <= n.
 *
 * Only the case where p = 1 is presently implemented, so the numbering
 * scheme is simply based on entity's global labels.
 *
 * For p > 1, it would probably be necessary to extend the numbering
 * schemes so as to account for the fact that a given entity may have
 * a main index on its main associated domain, but may be present
 * as a ghost entity with another index on neighboring domains.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_mpi.h"
#include "pdm_points_merge.h"
#include "pdm_timer.h"
#include "pdm_distrib.h"
#include "pdm_logging.h"
#include "pdm_unique.h"
#include "pdm_order.h"
#include "pdm_part_comm_graph_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_gnum.h"
#include "pdm_gnum_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Use bubble sort on an expectedly short sequence of coordinates
 * to ensure lexicographical ordering.
 *
 *  \param[in]      dim        <-- spatial dimension
 *  \param[in]      start_id   <-- start id in array
 *  \param[in]      end_id     <-- past-the-end id in array
 *  \param[in]      coords     <-- pointer to entity coordinates (interlaced)
 *  \param[in, out] order      <-> ordering array base on Morton encoding, or
 *                                 lexicographical coordinate ordering for ties
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

inline static void
_reorder_coords_lexicographic
(
int                dim,
size_t             start_id,
size_t             end_id,
const double       coords[],
PDM_l_num_t        order[]
)
{
  size_t  i;
  bool g_swap;

  do {

    g_swap = false;

    // remove i-1 to avoid compiler warning since 0-1 in size_t is max_size_t
    for (i = start_id; i <= end_id; i++) {

      size_t j_prev = order[i], j = order[i+1];
      bool l_swap = false;

      if (dim == 3) {
        if (coords[j_prev*3] < coords[j*3])
          continue;
        else if (coords[j_prev*3] > coords[j*3])
          l_swap = true;
        else if (coords[j_prev*3 + 1] < coords[j*3 + 1])
          continue;
        else if (   coords[j_prev*3 + 1] > coords[j*3 + 1]
                 || coords[j_prev*3 + 2] > coords[j*3 + 2])
          l_swap = true;
      }
      else if (dim == 2) {
        if (coords[j_prev*2] < coords[j*2 + 1])
          continue;
        else if (   coords[j_prev*2]     > coords[j*2]
                 || coords[j_prev*2 + 1] > coords[j*2 + 1])
          l_swap = true;
      }
      else { /* if (dim == 1) */
        if (coords[j_prev] > coords[j])
          l_swap = true;
      }

      if (l_swap) {
        PDM_l_num_t o_save = order[i];
        order[i] = order[i+1];
        order[i+1] = o_save;
        g_swap = true;
      }
    }

  } while (g_swap);
}


/**
 *
 * \brief Creation of an I/O numbering structure based on coordinates.
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be different, but their order is undetermined.
 *
 *  \param [in]      dim        spatial dimension
 *  \param [in]      n_entities number of entities considered
 *  \param [in]      coords     pointer to entity coordinates (interlaced)
 *  \param [in]      m_code     Morton code associated with each entity
 *  \param [in, out] order      ordering array base on Morton encoding, or
 *                              lexicographical coordinate ordering for ties
 *
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

static void
_check_morton_ordering
(
int                      dim,
size_t                   n_entities,
const double             coords[],
const PDM_morton_code_t  m_code[],
PDM_l_num_t              order[]
)
{
  size_t  i_prev = 0, i = 1;

  if (n_entities == 0)
    return;

  /* Check ordering; if two entities have the same Morton codes,
     use lexicographical coordinates ordering to ensure the
     final order is deterministic. */

  for (i = 1; i < n_entities; i++) {

    size_t j_prev = order[i_prev], j = order[i];

    if (   m_code[j_prev].X[0] != m_code[j].X[0]
        || m_code[j_prev].X[1] != m_code[j].X[1]
        || m_code[j_prev].X[2] != m_code[j].X[2]) {

      /* If successive values have the same Morton code,
         order them lexicographically */
      if (i_prev < i - 1)
        _reorder_coords_lexicographic(dim, i_prev, i-1, coords, order);

    }
    i_prev = i;
  }

  if (i_prev < n_entities - 1)
    _reorder_coords_lexicographic(dim, i_prev, n_entities - 1, coords, order);
}


/**
 *
 * \brief Compute from coords
 *
 * \param [in]   _gnum          Current _pdm_gen_gnum_t structure
 *
 */

/*
  This function comes from "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2004-2009  EDF

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

static void
_gnum_from_coords_compute
(
 PDM_gen_gnum_t *gen_gnum
)
{
  double extents[6];
  PDM_l_num_t  *order = NULL;
  PDM_morton_code_t *m_code = NULL;

  PDM_MPI_Comm comm = gen_gnum->comm;

  const int level = sizeof(PDM_morton_int_t)*8 - 1;

  int n_ranks;
  PDM_MPI_Comm_size (comm, &n_ranks);

  /* Merge double points */

  PDM_points_merge_t* pts_merge = NULL;

  int n_entities = 0;

  int iproc;
  PDM_MPI_Comm_rank (comm, &iproc);

  if (gen_gnum->merge) {

    PDM_malloc(gen_gnum->index, gen_gnum->n_part, int *);

    pts_merge = PDM_points_merge_create (gen_gnum->n_part, gen_gnum->tolerance, gen_gnum->comm, PDM_OWNERSHIP_KEEP);

    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
      PDM_malloc(gen_gnum->index[i_part], gen_gnum->n_elts[i_part], int);
      PDM_points_merge_cloud_set (pts_merge, i_part, gen_gnum->n_elts[i_part],
                                  gen_gnum->coords[i_part], gen_gnum->char_length[i_part]);
      for (int i = 0; i < gen_gnum->n_elts[i_part]; i++) {
        gen_gnum->index[i_part][i] = 0;
      }
    }

//    PDM_timer_t *timer = PDM_timer_create();
//    PDM_timer_resume(timer);
    PDM_points_merge_process (pts_merge);
//    PDM_timer_hang_on(timer);
//    printf("Compute points merge %12.5es\n", PDM_timer_elapsed(timer));
//    PDM_timer_free(timer);

    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {

      int *candidates_idx;
      int *candidates_desc;

      PDM_points_merge_candidates_get (pts_merge, i_part, &candidates_idx, &candidates_desc);

      for (int i = 0; i < gen_gnum->n_elts[i_part]; i++) {
        for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
          int idx = j;
          int distant_proc = candidates_desc[3*idx    ];
          int distant_part = candidates_desc[3*idx + 1];
          int distant_pt = candidates_desc[3*idx + 2];

          if ((distant_proc < iproc) ||
              ((distant_proc == iproc) && (distant_part < i_part)) ||
              ((distant_proc == iproc) && (distant_part == i_part)
               && (distant_pt < i))) {
            gen_gnum->index[i_part][i] = -1;
          }
        }
        if (gen_gnum->index[i_part][i] == 0) {
          gen_gnum->index[i_part][i] = n_entities;
          n_entities++;
        }
      }
    }
  }
  else {
    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
      n_entities += gen_gnum->n_elts[i_part];
    }
  }

  double *coords;
  PDM_malloc(coords, gen_gnum->dim * n_entities, double);

  if (gen_gnum->merge) {
    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
      for (int i = 0; i <  gen_gnum->n_elts[i_part]; i++) {
        if (gen_gnum->index[i_part][i] != -1) {
          for (int k = 0; k < gen_gnum->dim; k++) {
            coords[gen_gnum->dim * gen_gnum->index[i_part][i] + k] =
            gen_gnum->coords[i_part][gen_gnum->dim * i + k];
          }
        }
      }
    }
  }
  else {
    int k = 0;
    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
      for (int j = 0; j < gen_gnum->dim * gen_gnum->n_elts[i_part]; j++) {
        coords[k++] = gen_gnum->coords[i_part][j];
      }
    }
  }

  /* Build Morton encoding and order it */

  PDM_morton_get_coord_extents(gen_gnum->dim, n_entities, coords, extents, comm);
  if (0 && iproc == 0) {
    printf("  _gnum_from_coords_compute : PDM_morton_get_coord_extents OK\n");
    fflush(stdout);
  }

  PDM_malloc(m_code, n_entities, PDM_morton_code_t);
  order = NULL;

  double d[3];
  double s[3];
  PDM_morton_encode_coords(gen_gnum->dim, level, extents, n_entities, coords, m_code, d, s);
  if (0 && iproc == 0) {
    printf("  _gnum_from_coords_compute : PDM_morton_encode_coords OK\n");
    fflush(stdout);
  }

  if (0 && iproc == 0) {
    printf("  _gnum_from_coords_compute : PDM_morton_local_order OK\n");
    fflush(stdout);
  }

  if (n_ranks > 1) {

    int rank_id;
//    PDM_l_num_t j;
    PDM_l_num_t shift;

    int n_block_ents = 0;
    PDM_g_num_t current_global_num = 0, global_num_shift = 0;

    int *c_rank = NULL;
    int *send_count = NULL, *send_shift = NULL;
    int *recv_count = NULL, *recv_shift = NULL;
    double *send_coords = NULL, *recv_coords = NULL;
    double *weight = NULL;
    PDM_g_num_t *block_global_num = NULL, *part_global_num = NULL;
    PDM_morton_code_t *morton_index = NULL;

    PDM_malloc(weight      , n_entities , double           );
    PDM_malloc(morton_index, n_ranks + 1, PDM_morton_code_t);

    for (int i = 0; i < n_entities; i++) {
      weight[i] = 1;
    }

    PDM_morton_build_rank_index(gen_gnum->dim,
                                level,
                                n_entities,
                                m_code,
                                weight,
                                NULL, //order
                                morton_index,
                                comm);
    if (0 && iproc == 0) {
      printf("  _gnum_from_coords_compute : PDM_morton_build_rank_index OK\n");
      fflush(stdout);
    }

    //PDM_free(order);
    PDM_free(weight);
    PDM_malloc(c_rank, n_entities, int);

    for (int i = 0; i < n_entities; i++) {
      size_t _c_rank = PDM_morton_quantile_search((size_t) n_ranks,
                                                   m_code[i],
                                                   morton_index);
      c_rank[i] = (int) _c_rank;
    }

    PDM_free(morton_index);
    PDM_free(m_code);

    /* Build send_buf, send_count and send_shift
       to build a rank to coords indexed list */

    PDM_malloc(send_count, n_ranks    , int);
    PDM_malloc(recv_count, n_ranks    , int);
    PDM_malloc(send_shift, n_ranks + 1, int);
    PDM_malloc(recv_shift, n_ranks + 1, int);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (int i = 0; i < n_entities; i++)
      send_count[c_rank[i]] += gen_gnum->dim;

    /* Exchange number of coords to send to each process */

    PDM_MPI_Alltoall(send_count, 1, PDM_MPI_INT, recv_count, 1, PDM_MPI_INT, comm);
    if (0 && iproc == 0) {
      printf("  _gnum_from_coords_compute : Alltoall OK\n");
      fflush(stdout);
    }

    send_shift[0] = 0;
    recv_shift[0] = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_shift[rank_id + 1] = send_shift[rank_id] + send_count[rank_id];
      recv_shift[rank_id + 1] = recv_shift[rank_id] + recv_count[rank_id];
    }

    /* Build send and receive buffers */

    PDM_malloc(send_coords, send_shift[n_ranks], double);

    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      send_count[rank_id] = 0;

    for (int i = 0; i < n_entities; i++) {
      rank_id = c_rank[i];
      shift = send_shift[rank_id] + send_count[rank_id];
      for (int j = 0; j < gen_gnum->dim; j++)
        send_coords[shift + j] = coords[i*gen_gnum->dim + j];
      send_count[rank_id] += gen_gnum->dim;
    }

    PDM_malloc(recv_coords, recv_shift[n_ranks], double);

    /* Exchange coords between processes */

    PDM_MPI_Alltoallv(send_coords, send_count, send_shift, PDM_MPI_DOUBLE,
                      recv_coords, recv_count, recv_shift, PDM_MPI_DOUBLE,
                      comm);
    if (0 && iproc == 0) {
      printf("  _gnum_from_coords_compute : Alltoallv OK\n");
      fflush(stdout);
    }

    PDM_free(send_coords);

    /* Now re-build Morton codes on block distribution */

    n_block_ents = recv_shift[n_ranks] / gen_gnum->dim;

    PDM_malloc(m_code, n_block_ents, PDM_morton_code_t);
    PDM_malloc(order , n_block_ents, PDM_l_num_t      );

    PDM_morton_encode_coords(gen_gnum->dim,
                             level,
                             extents,
                             n_block_ents,
                             recv_coords,
                             m_code,
                             d,
                             s);

    PDM_morton_local_order((int) n_block_ents, m_code, order);

    /* Check ordering; if two entities have the same Morton codes,
       use lexicographical coordinates ordering to ensure the
       final order is deterministic. */

    _check_morton_ordering(gen_gnum->dim, n_block_ents, recv_coords, m_code, order);

    /* Determine global order; requires ordering to loop through buffer by
       increasing number (slice blocks associated with each process are
       already sorted, but the whole "gathered" slice is not).
       We build an initial global order based on the initial global numbering,
       such that for each slice, the global number of an entity is equal to
       the cumulative number of sub-entities */

    PDM_free(m_code);
    PDM_free(recv_coords);
    PDM_malloc(block_global_num, n_block_ents, PDM_g_num_t);

    for (int i = 0; i < n_block_ents; i++) {
      block_global_num[order[i]] = (PDM_g_num_t) i + 1;
    }

    PDM_free(order);

    current_global_num = (PDM_g_num_t) n_block_ents;

    /* At this stage, block_global_num[] is valid for this process, and
       current_global_num indicates the total number of entities handled
       by this process; we must now shift global numberings on different
       processes by the cumulative total number of entities handled by
       each process */

    PDM_MPI_Scan(&current_global_num, &global_num_shift, 1, PDM__PDM_MPI_G_NUM,
                 PDM_MPI_SUM, comm);
    global_num_shift -= current_global_num;

    for (int i = 0; i < n_block_ents; i++)
      block_global_num[i] += global_num_shift;

    /* Return global order to all processors */

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] /= gen_gnum->dim;
      recv_count[rank_id] /= gen_gnum->dim;
      send_shift[rank_id] /= gen_gnum->dim;
      recv_shift[rank_id] /= gen_gnum->dim;
    }

    send_shift[n_ranks] /= gen_gnum->dim;

    PDM_malloc(part_global_num, send_shift[n_ranks], PDM_g_num_t);

    PDM_MPI_Alltoallv(block_global_num, recv_count, recv_shift, PDM__PDM_MPI_G_NUM,
                      part_global_num, send_count, send_shift, PDM__PDM_MPI_G_NUM,
                      comm);

    for (rank_id = 0; rank_id < n_ranks; rank_id++) {
      send_count[rank_id] = 0;
    }

    PDM_g_num_t _max_loc = -1;

    if (gen_gnum->merge) {

      /*
       * Define local points
       */

      int k = 0;
      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
        PDM_malloc(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
        for (int j1 = 0; j1 < gen_gnum->n_elts[i_part]; j1++) {
          gen_gnum->g_nums[i_part][j1] = -1;
        }
        for (int j1 = 0; j1 < gen_gnum->n_elts[i_part]; j1++) {
          if (gen_gnum->index[i_part][j1] != -1) {
            rank_id = c_rank[k++];
            shift = send_shift[rank_id] + send_count[rank_id];
            gen_gnum->g_nums[i_part][j1] = part_global_num[shift];
            _max_loc = PDM_MAX (_max_loc, part_global_num[shift]);
            send_count[rank_id] += 1;
          }
        }
      }

      int *send_count2 = NULL;
      int *recv_count2 = NULL;
      int *send_shift2 = NULL;
      int *recv_shift2 = NULL;
      PDM_malloc(send_count2, n_ranks    , int);
      PDM_malloc(recv_count2, n_ranks    , int);
      PDM_malloc(send_shift2, n_ranks + 1, int);
      PDM_malloc(recv_shift2, n_ranks + 1, int);

      /*
       * Count number of values to send
       */

      for (rank_id = 0; rank_id < n_ranks; rank_id++) {
        send_count2[rank_id] = 0;
      }

      for (rank_id = 0; rank_id < n_ranks + 1; rank_id++) {
        send_shift2[rank_id] = 0;
        recv_shift2[rank_id] = 0;
      }

      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {

        int *candidates_idx;
        int *candidates_desc;


        PDM_points_merge_candidates_get (pts_merge, i_part, &candidates_idx, &candidates_desc);

        for (int i = 0; i < gen_gnum->n_elts[i_part]; i++) {
          int update_proc = iproc;

          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
//            int distant_part = candidates_desc[3*idx + 1];
//            int distant_pt   = candidates_desc[3*idx + 2];
            update_proc = PDM_MIN(update_proc, distant_proc);
          }
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];
            if ((update_proc == iproc) && (distant_proc != iproc)) {
              send_count2[distant_proc] += 1;
              assert (gen_gnum->g_nums[i_part][i] != -1);
            }
            else if (  (update_proc == iproc)
                    && (distant_proc == iproc)
                    && (i_part <= distant_part)) {
              assert (gen_gnum->g_nums[i_part][i] != -1);
              gen_gnum->g_nums[distant_part][distant_pt] = gen_gnum->g_nums[i_part][i];
            }
          }
        }
      }

      PDM_MPI_Alltoall (send_count2, 1, PDM_MPI_INT,
                        recv_count2, 1, PDM_MPI_INT, comm);

      for (rank_id = 0; rank_id < n_ranks; rank_id++) {
        send_shift2[rank_id + 1] = send_shift2[rank_id] + 3 * send_count2[rank_id];
        send_count2[rank_id] = 0;

        recv_count2[rank_id] *= 3;
        recv_shift2[rank_id + 1] = recv_shift2[rank_id] + recv_count2[rank_id];
      }

      /*
       * Send values
       */

      PDM_g_num_t *send_buff = NULL;
      PDM_g_num_t *recv_buff = NULL;
      PDM_malloc(send_buff, send_shift2[n_ranks], PDM_g_num_t);
      PDM_malloc(recv_buff, recv_shift2[n_ranks], PDM_g_num_t);

      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {

        int *candidates_idx;
        int *candidates_desc;

        PDM_points_merge_candidates_get (pts_merge, i_part, &candidates_idx, &candidates_desc);

        for (int i = 0; i < gen_gnum->n_elts[i_part]; i++) {
          int update_proc = iproc;

          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
//            int distant_part = candidates_desc[3*idx + 1];
//            int distant_pt   = candidates_desc[3*idx + 2];
            update_proc = PDM_MIN (update_proc, distant_proc);
          }

          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];

            if ((iproc == update_proc) && (distant_proc != iproc)){
              int idx2 = send_shift2[distant_proc] + send_count2[distant_proc];
              send_buff[idx2]     = distant_part;
              send_buff[idx2 + 1] = distant_pt;
              send_buff[idx2 + 2] = gen_gnum->g_nums[i_part][i];;

              send_count2[distant_proc] += 3;
            }
          }
        }
      }

      /*
       * Send : distant_part, distant_pt, gnum
       */

      PDM_MPI_Alltoallv(send_buff, send_count2, send_shift2, PDM__PDM_MPI_G_NUM,
                        recv_buff, recv_count2, recv_shift2, PDM__PDM_MPI_G_NUM,
                        comm);

      /*
       * update gnum
       */

      k = 0;
      while (k < recv_shift2[n_ranks]) {
        int i_part        = (int) recv_buff[k++];
        int ipt          = (int) recv_buff[k++];
        PDM_g_num_t gnum =       recv_buff[k++];
         gen_gnum->g_nums[i_part][ipt] = gnum;
      }

      PDM_free(send_buff);
      PDM_free(recv_buff);
      PDM_free(send_count2);
      PDM_free(recv_count2);
      PDM_free(send_shift2);
      PDM_free(recv_shift2);

    }

    else {
      int k = 0;
      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
        PDM_malloc(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
        for (int j1 = 0; j1 < gen_gnum->n_elts[i_part]; j1++) {
          rank_id = c_rank[k++];
          shift = send_shift[rank_id] + send_count[rank_id];
          gen_gnum->g_nums[i_part][j1] = part_global_num[shift];
          _max_loc = PDM_MAX (_max_loc, part_global_num[shift]);
          send_count[rank_id] += 1;
        }
      }
    }

    /* Free memory */

    PDM_free(c_rank);

    PDM_free(block_global_num);
    PDM_free(part_global_num);

    PDM_free(send_count);
    PDM_free(recv_count);
    PDM_free(send_shift);
    PDM_free(recv_shift);

    /* Get final maximum global number value */

    PDM_MPI_Allreduce (&_max_loc,
                       &gen_gnum->n_g_elt,
                       1,
                       PDM__PDM_MPI_G_NUM,
                       PDM_MPI_MAX,
                       comm);

  }

  else if (n_ranks == 1) {

    PDM_malloc(order, n_entities, PDM_l_num_t);

    PDM_morton_local_order(n_entities, m_code, order);

    _check_morton_ordering(gen_gnum->dim, n_entities, coords, m_code, order);

    PDM_free(m_code);

    PDM_g_num_t *tmp_gnum = NULL;
    PDM_malloc(tmp_gnum, n_entities, PDM_g_num_t);

    for (int i = 0; i < n_entities; i++) {
      tmp_gnum[order[i]] = (PDM_g_num_t) i+1;
    }

    if (gen_gnum->merge) {

      int *_entities = NULL;
      PDM_malloc(_entities, 2 * n_entities, int);
      int k = 0;
      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
        PDM_malloc(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
        for (int j1 = 0; j1 < gen_gnum->n_elts[i_part]; j1++) {
          gen_gnum->g_nums[i_part][j1] = -1;
          if (gen_gnum->index[i_part][j1] != -1) {
            _entities[k++] = i_part;
            _entities[k++] = j1;
          }
        }
      }

      k = 0;
      for (int i = 0; i < n_entities; i++) {
        gen_gnum->g_nums[_entities[2*k]][_entities[2*k+1]] = tmp_gnum[k];
        k += 1;
      }

      PDM_free(_entities);

      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
        int *candidates_idx;
        int *candidates_desc;
        PDM_points_merge_candidates_get (pts_merge, i_part, &candidates_idx, &candidates_desc);

        for (int i = 0; i < gen_gnum->n_elts[i_part]; i++) {
          for (int j = candidates_idx[i]; j < candidates_idx[i+1]; j++) {
            int idx = j;
            int distant_proc = candidates_desc[3*idx    ];
            int distant_part = candidates_desc[3*idx + 1];
            int distant_pt   = candidates_desc[3*idx + 2];
            if ((iproc == distant_proc) && (i_part <= distant_part)) {
              assert (gen_gnum->g_nums[i_part][i] != -1);
              gen_gnum->g_nums[distant_part][distant_pt] = gen_gnum->g_nums[i_part][i];
            }
          }
        }

      }
    }

    else {

      int k = 0;
      for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
        PDM_malloc(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
        for (int j1 = 0; j1 <  gen_gnum->n_elts[i_part]; j1++) {
          gen_gnum->g_nums[i_part][j1] = tmp_gnum[k++];
        }
      }
    }

    PDM_free(order);
    PDM_free(tmp_gnum);

    gen_gnum->n_g_elt = n_entities;

  }

  if (gen_gnum->merge) {
    for (int i_part = 0; i_part < gen_gnum->n_part; i_part++) {
      PDM_free(gen_gnum->index[i_part]);
    }
    PDM_free(gen_gnum->index);
    PDM_points_merge_free (pts_merge);
  }

  PDM_free(coords);

}


/**
 *
 * \brief Compute from coords
 *
 * \param [in]   _gnum          Current _pdm_gen_gnum_t structure
 *
 */

static void
_gnum_from_parent_compute
(
 PDM_gen_gnum_t *gen_gnum
)
{
  int n_procs = 0;
  PDM_MPI_Comm_size(gen_gnum->comm, &n_procs);

  int i_proc = 0;
  PDM_MPI_Comm_rank(gen_gnum->comm,
                &i_proc);

  int *send_buff_n   = NULL;
  int *send_buff_idx = NULL;
  PDM_malloc(send_buff_n  , n_procs, int);
  PDM_malloc(send_buff_idx, n_procs, int);

  int *recv_buff_n   = NULL;
  int *recv_buff_idx = NULL;
  PDM_malloc(recv_buff_n  , n_procs, int);
  PDM_malloc(recv_buff_idx, n_procs, int);

  /* Calcul du nombre total d'elements du bloc */

  PDM_l_num_t n_elt_loc_total = 0;
  PDM_g_num_t l_max_parent = 0;

  for (int j = 0; j < gen_gnum->n_part; j++) {
    n_elt_loc_total += gen_gnum->n_elts[j];
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      l_max_parent = PDM_MAX (l_max_parent, gen_gnum->parent[j][k]);
    }
  }

  PDM_g_num_t max_parent = 0;
  PDM_MPI_Allreduce (&l_max_parent, &max_parent, 1,
                     PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, gen_gnum->comm);

  /* Comptage du nombre d'elements a envoyer a chaque processus */

  for (int j = 0; j < n_procs; j++) {
    send_buff_n[j]   = 0;
    send_buff_idx[j] = 0;
    recv_buff_n[j]   = 0;
    recv_buff_idx[j] = 0;
  }

  PDM_g_num_t *d_elt_proc = NULL;
  PDM_malloc(d_elt_proc, n_procs + 1, PDM_g_num_t);


  PDM_g_num_t div_entiere = max_parent / n_procs;
  PDM_g_num_t div_reste = max_parent % n_procs;

  d_elt_proc[0] = 1;
  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] =  div_entiere;
    if (i < div_reste) {
      d_elt_proc[i+1] += 1;
    }
  }

  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] += d_elt_proc[i];
  }

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long (gen_gnum->parent[j][k],
                                                         d_elt_proc,
                                                         n_procs + 1);
      send_buff_n[i_elt_proc] += 1;
    }
  }

  send_buff_idx[0] = 0;
  for (int j = 1; j < n_procs; j++) {
    send_buff_idx[j] = send_buff_idx[j-1] + send_buff_n[j-1];
  }

  /* Determination du nombre d'elements recu de chaque processus */

  PDM_MPI_Alltoall(send_buff_n,
               1,
               PDM_MPI_INT,
               recv_buff_n,
               1,
               PDM_MPI_INT,
               gen_gnum->comm);

  recv_buff_idx[0] = 0;
  for(int j = 1; j < n_procs; j++) {
    recv_buff_idx[j] = recv_buff_idx[j-1] + recv_buff_n[j-1];
  }

  /* Transmission des numeros absolus  */

  PDM_g_num_t _l_numabs_tmp = d_elt_proc[i_proc+1] - d_elt_proc[i_proc];
  int l_numabs_tmp = (int) _l_numabs_tmp;

  PDM_g_num_t *numabs_tmp         = NULL;
  PDM_g_num_t *n_elt_stocke_procs = NULL;
  PDM_malloc(numabs_tmp        , l_numabs_tmp, PDM_g_num_t);
  PDM_malloc(n_elt_stocke_procs, n_procs + 1 , PDM_g_num_t);

  PDM_g_num_t *send_buff_numabs = NULL;
  PDM_g_num_t *recv_buff_numabs = NULL;
  PDM_malloc(send_buff_numabs, n_elt_loc_total                                      , PDM_g_num_t);
  PDM_malloc(recv_buff_numabs, recv_buff_idx[n_procs - 1] + recv_buff_n[n_procs - 1], PDM_g_num_t);

  for (int j = 0; j < n_procs; j++) {
    send_buff_n[j] = 0;
  }

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(gen_gnum->parent[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      send_buff_numabs[send_buff_idx[i_elt_proc] + send_buff_n[i_elt_proc]] = gen_gnum->parent[j][k];
      send_buff_n[i_elt_proc] += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) send_buff_numabs,
                send_buff_n,
                send_buff_idx,
                PDM__PDM_MPI_G_NUM,
                (void *) recv_buff_numabs,
                recv_buff_n,
                recv_buff_idx,
                PDM__PDM_MPI_G_NUM,
                gen_gnum->comm);

  /* Echange du nombre d'elements stockes sur chaque processus */

  const PDM_g_num_t n_elt_stocke =
    (PDM_g_num_t) (recv_buff_idx[n_procs - 1] + recv_buff_n[n_procs - 1]);

  PDM_MPI_Allgather((void *) &n_elt_stocke,
                1,
                PDM__PDM_MPI_G_NUM,
                (void *) (n_elt_stocke_procs + 1),
                1,
                PDM__PDM_MPI_G_NUM,
                gen_gnum->comm);

  n_elt_stocke_procs[0] = 1;
  for (int j = 1; j < n_procs + 1; j++) {
    n_elt_stocke_procs[j] += n_elt_stocke_procs[j-1];
  }

  /* Stockage du resultat et determination de la nouvelle numerotation absolue
     independante du parallelisme */

  for (int j = 0; j < l_numabs_tmp; j++) {
    numabs_tmp[j] = 0;
  }

  for (int j = 0; j < n_procs; j++) {

    const int ideb = recv_buff_idx[j];
    const int ifin = recv_buff_idx[j] + recv_buff_n[j];

    for (int k = ideb; k < ifin; k++) {

      PDM_g_num_t _idx = recv_buff_numabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;
      assert((idx < l_numabs_tmp) && (idx >= 0));

      numabs_tmp[idx] = 1; /* On marque les elements */
    }
  }

  PDM_g_num_t cpt_elt_proc = 0;
  for (int j = 0; j < l_numabs_tmp; j++) {
    if (numabs_tmp[j] == 1) {
      cpt_elt_proc += 1;
    }
  }

  /* Mise a jour de n_elt_stocke_procs */

  PDM_MPI_Allgather((void *) &cpt_elt_proc,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (n_elt_stocke_procs + 1),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);

  n_elt_stocke_procs[0] = 1;
  for (int j = 1; j < n_procs + 1; j++) {
    n_elt_stocke_procs[j] += n_elt_stocke_procs[j-1];
  }

  /* On fournit une numerotation independante du parallelisme */

  cpt_elt_proc = 0;
  for (int j = 0; j < l_numabs_tmp; j++) {
    if (numabs_tmp[j] == 1) {

      numabs_tmp[j] = n_elt_stocke_procs[i_proc] + cpt_elt_proc;
      cpt_elt_proc += 1;
    }
  }

  /* On remplit le buffer de reception qui devient le buffer d'envoi
     Le buffer d'envoi devient lui le buffer de reception */

  cpt_elt_proc = 0;
  for (int j = 0; j < n_procs; j++) {

    const int ideb = recv_buff_idx[j];
    const int ifin = recv_buff_idx[j] + recv_buff_n[j];

    for (int k = ideb; k < ifin; k++) {

      PDM_g_num_t _idx = recv_buff_numabs[k] - d_elt_proc[i_proc];
      const int idx = (int) _idx;

      recv_buff_numabs[cpt_elt_proc] = numabs_tmp[idx];

      cpt_elt_proc += 1;
    }
  }

  PDM_MPI_Alltoallv((void *) recv_buff_numabs,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    (void *) send_buff_numabs,
                    send_buff_n,
                    send_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);

  /* On Stocke l'information recue */

  for (int j = 0; j < n_procs; j++) {
    send_buff_n[j] = 0;
  }

  for (int j = 0; j < gen_gnum->n_part; j++) {

    PDM_malloc(gen_gnum->g_nums[j], gen_gnum->n_elts[j], PDM_g_num_t);

    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(gen_gnum->parent[j][k],
                                                        d_elt_proc,
                                                        n_procs+1);
      gen_gnum->g_nums[j][k] = send_buff_numabs[send_buff_idx[i_elt_proc] + send_buff_n[i_elt_proc]];
      send_buff_n[i_elt_proc] += 1;
    }
  }

  /* Liberation memoire */

  PDM_free(send_buff_idx);
  PDM_free(send_buff_n);
  PDM_free(recv_buff_idx);
  PDM_free(recv_buff_n);
  PDM_free(send_buff_numabs);
  PDM_free(recv_buff_numabs);
  PDM_free(d_elt_proc);
  PDM_free(numabs_tmp);
  PDM_free(n_elt_stocke_procs);

}


static void
_gnum_from_parent_compute_opt
(
 PDM_gen_gnum_t *gen_gnum
)
{
  int i_rank = -1;
  int n_rank = -1;
  PDM_MPI_Comm_rank (gen_gnum->comm, &i_rank);
  PDM_MPI_Comm_size (gen_gnum->comm, &n_rank);

  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol          = 0.10;
  PDM_g_num_t* distrib = NULL;
  PDM_distrib_weight(    sampling_factor,
                         n_rank,
                         gen_gnum->n_part,
                         gen_gnum->n_elts,
  (const PDM_g_num_t **) gen_gnum->parent,
                         NULL,
                         n_iter_max,
                         tol,
                         gen_gnum->comm,
                         &distrib);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib, n_rank+1, "distrib_key :");
  }

  /*
   * Create send buffer
   */
  int *send_buff_n   = NULL;
  int *send_buff_idx = NULL;
  PDM_malloc(send_buff_n  , n_rank    , int);
  PDM_malloc(send_buff_idx, n_rank + 1, int);

  int *recv_buff_n   = NULL;
  int *recv_buff_idx = NULL;
  PDM_malloc(recv_buff_n  , n_rank    , int);
  PDM_malloc(recv_buff_idx, n_rank + 1, int);

  /* Comptage du nombre d'elements a envoyer a chaque processus */
  for (int j = 0; j < n_rank; j++) {
    send_buff_n  [j] = 0;
    send_buff_idx[j] = 0;
    recv_buff_n  [j] = 0;
    recv_buff_idx[j] = 0;
  }

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (gen_gnum->parent[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);
      // log_trace("t_rank = %i | key = %i \n", t_rank, key_ln_to_gn[j][k]);
      send_buff_n[t_rank] += 1;
    }
  }

  send_buff_idx[0] = 0;
  for (int j = 0; j < n_rank; j++) {
    send_buff_idx[j+1] = send_buff_idx[j] + send_buff_n[j];
  }

  /* Determination du nombre d'elements recu de chaque processus */
  PDM_MPI_Alltoall(send_buff_n,
                   1,
                   PDM_MPI_INT,
                   recv_buff_n,
                   1,
                   PDM_MPI_INT,
                   gen_gnum->comm);


  recv_buff_idx[0] = 0;
  for (int j = 0; j < n_rank; j++) {
    recv_buff_idx[j+1] = recv_buff_idx[j] + recv_buff_n[j];
    send_buff_n  [j]   = 0;
  }

  /*
   * Exchange key and val associate
   */
  PDM_g_num_t *send_gnum = NULL;
  PDM_malloc(send_gnum, send_buff_idx[n_rank], PDM_g_num_t);

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (gen_gnum->parent[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);

      int idx_write = send_buff_idx[t_rank] + send_buff_n[t_rank]++;
      send_gnum[idx_write] = gen_gnum->parent[j][k];
    }
  }

  PDM_g_num_t *recv_gnum;
  PDM_malloc(recv_gnum, recv_buff_idx[n_rank], PDM_g_num_t);
  PDM_MPI_Alltoallv((void *) send_gnum,
                    send_buff_n,
                    send_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    (void *) recv_gnum,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);


  /*
   *
   */
  int *unique_order = NULL;
  PDM_malloc(unique_order, recv_buff_idx[n_rank] ,int);
  int n_unique = PDM_inplace_unique_long2(recv_gnum, unique_order, 0, recv_buff_idx[n_rank]-1);

  PDM_g_num_t current_global_num = n_unique;
  PDM_g_num_t global_num_shift   = 0;
  PDM_MPI_Scan(&current_global_num, &global_num_shift, 1, PDM__PDM_MPI_G_NUM,
               PDM_MPI_SUM, gen_gnum->comm);
  global_num_shift -= current_global_num;

  for(int i = 0; i < recv_buff_idx[n_rank]; ++i) {
    recv_gnum[i] = (PDM_g_num_t) unique_order[i] + global_num_shift + 1;
  }

  PDM_free(unique_order);

  /*
   * Reverse exchange
   */
  PDM_MPI_Alltoallv((void *) recv_gnum,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    (void *) send_gnum,
                    send_buff_n,
                    send_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);

  // PDM_log_trace_array_long(send_gnum, send_buff_idx[n_rank], "send_gnum :");

  /* On Stocke l'information recue */
  for (int j = 0; j < n_rank; j++) {
    send_buff_n  [j]   = 0;
  }
  for (int j = 0; j < gen_gnum->n_part; j++) {

    PDM_malloc(gen_gnum->g_nums[j], gen_gnum->n_elts[j], PDM_g_num_t);

    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (gen_gnum->parent[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);

      int idx_read = send_buff_idx[t_rank] + send_buff_n[t_rank]++;
      gen_gnum->g_nums[j][k] = send_gnum[idx_read];
    }
  }

  PDM_free(send_gnum);
  PDM_free(recv_gnum);

  PDM_free(send_buff_n  );
  PDM_free(send_buff_idx);
  PDM_free(recv_buff_n  );
  PDM_free(recv_buff_idx);

  PDM_free(distrib);

}

static void
_gnum_from_parent_compute_nuplet
(
 PDM_gen_gnum_t *gen_gnum
)
{
  int i_rank = -1;
  int n_rank = -1;
  PDM_MPI_Comm_rank (gen_gnum->comm, &i_rank);
  PDM_MPI_Comm_size (gen_gnum->comm, &n_rank);

  /* Generate a keys */
  int nuplet = gen_gnum->nuplet;
  PDM_g_num_t **key_ln_to_gn = NULL;
  PDM_malloc(key_ln_to_gn, gen_gnum->n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    PDM_malloc(key_ln_to_gn[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
    for(int i = 0; i < gen_gnum->n_elts[i_part]; ++i) {
      /* La boucle suivante (ou l'imbrication des trois boucles) semble poser probleme
       * au compilo intel (icc) en optimisé.
       * On la remplace par l'utilisation d'une variable temporaire (cf ci-dessous)
      key_ln_to_gn[i_part][i] = 1;
      for(int k = 0; k < nuplet; ++k) {
        key_ln_to_gn[i_part][i] += PDM_ABS(gen_gnum->parent[i_part][nuplet * i + k]);
      }
      */
      PDM_g_num_t tmp = 1;
      for(int k = 0; k < nuplet; ++k) {
        tmp += PDM_ABS(gen_gnum->parent[i_part][nuplet * i + k]);
      }
      key_ln_to_gn[i_part][i] = tmp;
    }
  }

  int sampling_factor = 2;
  int n_iter_max      = 5;
  double tol          = 0.10;
  PDM_g_num_t* distrib = NULL;
  PDM_distrib_weight(    sampling_factor,
                         n_rank,
                         gen_gnum->n_part,
                         gen_gnum->n_elts,
  (const PDM_g_num_t **) key_ln_to_gn,
                         NULL,
                         n_iter_max,
                         tol,
                         gen_gnum->comm,
                         &distrib);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib, n_rank+1, "distrib_key :");
  }

  /*
   * Create send buffer
   */
  int *send_buff_n   = NULL;
  int *send_buff_idx = NULL;
  PDM_malloc(send_buff_n  , n_rank    , int);
  PDM_malloc(send_buff_idx, n_rank + 1, int);

  int *recv_buff_n   = NULL;
  int *recv_buff_idx = NULL;
  PDM_malloc(recv_buff_n  , n_rank    , int);
  PDM_malloc(recv_buff_idx, n_rank + 1, int);

  /* Comptage du nombre d'elements a envoyer a chaque processus */
  for (int j = 0; j < n_rank; j++) {
    send_buff_n  [j] = 0;
    send_buff_idx[j] = 0;
    recv_buff_n  [j] = 0;
    recv_buff_idx[j] = 0;
  }

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (key_ln_to_gn[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);
      // log_trace("t_rank = %i | key = %i \n", t_rank, key_ln_to_gn[j][k]);
      send_buff_n[t_rank] += 1;
    }
  }

  send_buff_idx[0] = 0;
  for (int j = 0; j < n_rank; j++) {
    send_buff_idx[j+1] = send_buff_idx[j] + send_buff_n[j];
  }

  /* Determination du nombre d'elements recu de chaque processus */
  PDM_MPI_Alltoall(send_buff_n,
                   1,
                   PDM_MPI_INT,
                   recv_buff_n,
                   1,
                   PDM_MPI_INT,
                   gen_gnum->comm);


  recv_buff_idx[0] = 0;
  for (int j = 0; j < n_rank; j++) {
    recv_buff_idx[j+1] = recv_buff_idx[j] + recv_buff_n[j];
    send_buff_n  [j]   = 0;
  }

  /*
   * Exchange key and val associate
   */
  PDM_g_num_t *send_key   = NULL;
  PDM_g_num_t *send_elmts = NULL;
  PDM_malloc(send_key  ,          send_buff_idx[n_rank], PDM_g_num_t);
  PDM_malloc(send_elmts, nuplet * send_buff_idx[n_rank], PDM_g_num_t);

  for (int j = 0; j < gen_gnum->n_part; j++) {
    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (key_ln_to_gn[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);

      int idx_write = send_buff_idx[t_rank] + send_buff_n[t_rank]++;
      send_key[idx_write] = key_ln_to_gn[j][k];
      for(int p = 0; p < nuplet; ++p) {
        send_elmts[nuplet*idx_write+p] = gen_gnum->parent[j][nuplet * k + p];
      }
    }
  }


  PDM_MPI_Datatype mpi_entity_type;
  int min_nuplet, max_nuplet;
  PDM_MPI_Allreduce(&nuplet, &min_nuplet, 1, PDM_MPI_INT, PDM_MPI_MIN, gen_gnum->comm);
  PDM_MPI_Allreduce(&nuplet, &max_nuplet, 1, PDM_MPI_INT, PDM_MPI_MAX, gen_gnum->comm);
  if (min_nuplet != nuplet && max_nuplet != nuplet) {
    PDM_error(__FILE__, __LINE__, 0, "Error : nuplet mismatch (min = %d, max = %d)\n", min_nuplet, max_nuplet);
  }
  PDM_MPI_Type_create_contiguous(nuplet, PDM__PDM_MPI_G_NUM, &mpi_entity_type);
  PDM_MPI_Type_commit(&mpi_entity_type);

  PDM_g_num_t *recv_key   = NULL;
  PDM_g_num_t *recv_elmts = NULL;
  PDM_malloc(recv_key  ,          recv_buff_idx[n_rank], PDM_g_num_t);
  PDM_malloc(recv_elmts, nuplet * recv_buff_idx[n_rank], PDM_g_num_t);
  PDM_MPI_Alltoallv((void *) send_key,
                    send_buff_n,
                    send_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    (void *) recv_key,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);

  PDM_MPI_Alltoallv((void *) send_elmts,
                    send_buff_n,
                    send_buff_idx,
                    mpi_entity_type,
                    (void *) recv_elmts,
                    recv_buff_n,
                    recv_buff_idx,
                    mpi_entity_type,
                    gen_gnum->comm);

  PDM_MPI_Type_free(&mpi_entity_type);

  /*
   * Sort incoming key
   */
  int n_recv_key = recv_buff_idx[n_rank];
  int *order = NULL;
  PDM_malloc(order, n_recv_key, int);
  PDM_order_gnum_s(recv_key, 1, order, n_recv_key);

  int n_conflit_to_solve = 0;
  PDM_g_num_t last_gnum = -1;

  int *key_conflict_idx;
  PDM_malloc(key_conflict_idx, n_recv_key+1, int);
  key_conflict_idx[0] = 0;
  for(int i = 0; i < n_recv_key; ++i) {
    if(recv_key[order[i]] != last_gnum){
      key_conflict_idx[n_conflit_to_solve+1] = key_conflict_idx[n_conflit_to_solve]+1;
      n_conflit_to_solve++;
      last_gnum = recv_key[order[i]];
    } else {
      key_conflict_idx[n_conflit_to_solve]++;
    }
  }

  int n_max_entity_per_key = 0;
  for(int i = 0; i < n_conflit_to_solve; ++i) {
    n_max_entity_per_key = PDM_MAX(n_max_entity_per_key, key_conflict_idx[i+1]-key_conflict_idx[i]);
  }

  /*
   * Solve conflict
   */
  if(0 == 1) {
    PDM_log_trace_array_int(key_conflict_idx, n_conflit_to_solve, "key_conflict_idx ::  ");
    for(int i = 0; i < n_conflit_to_solve; ++i) {
      log_trace(" ------ i = %i \n", i);
      for(int i_key = key_conflict_idx[i]; i_key < key_conflict_idx[i+1]; ++i_key) {
        int i_conflict = order[i_key];
        // int beg = recv_entity_vtx_idx[i_conflict];
        // int n_vtx_in_entity = recv_entity_vtx_idx[i_conflict+1] - beg;
        log_trace(" \t i_key = %i \n", recv_key[i_conflict]);
      }
    }
  }

  int         *already_treat   = NULL;
  int         *same_entity_idx = NULL;
  PDM_g_num_t *tmp_parent      = NULL;
  int         *order_parent    = NULL;
  PDM_malloc(already_treat  ,          n_max_entity_per_key   , int        );
  PDM_malloc(same_entity_idx,         (n_max_entity_per_key+1), int        );
  PDM_malloc(tmp_parent     , nuplet * n_max_entity_per_key   , PDM_g_num_t);
  PDM_malloc(order_parent   ,          n_max_entity_per_key   , int        );

  int i_abs_entity   = 0;
  for(int i = 0; i < n_conflit_to_solve; ++i) {

    int n_conflict_entitys = key_conflict_idx[i+1] - key_conflict_idx[i];
    for(int j = 0; j < n_conflict_entitys; ++j ) {
      already_treat[j] = -1;

      int i_conflict = order[key_conflict_idx[i]+j];
      int beg_elmt   = nuplet * i_conflict;

      for(int k = 0; k < nuplet; ++k) {
        tmp_parent[nuplet * j + k] = recv_elmts[beg_elmt + k];
      }
    }

    PDM_order_gnum_s(tmp_parent, nuplet, order_parent, n_conflict_entitys);

    // PDM_log_trace_array_int (order_parent, n_conflict_entitys, "order_parent  :" );
    // PDM_log_trace_array_long(tmp_parent, nuplet * n_conflict_entitys, "tmp_parent  :" );

    for(int idx_entity = 0; idx_entity < n_conflict_entitys; ++idx_entity) {
      int i_entity  = order[key_conflict_idx[i]+order_parent[idx_entity]];
      int beg_elmt1 = nuplet * i_entity;

      int idx_next_same_entity = 0;
      same_entity_idx[idx_next_same_entity++] = idx_entity;

      if(already_treat[idx_entity] != 1) {

        for(int idx_entity2 = 0; idx_entity2 < n_conflict_entitys; ++idx_entity2) {
          int i_entity_next = order[key_conflict_idx[i]+order_parent[idx_entity2]];
          int beg_elmt2 = nuplet * i_entity_next;

          if (i_entity_next == i_entity) {
            continue;
          }

          // printf("conflict : i_entity = %d, i_entity_next = %d...\n", i_entity, i_entity_next);
          if(already_treat[idx_entity2] == 1) {
            continue;
          }

          int is_same_entity = 1;
          for(int k = 0; k < nuplet; ++k) {
            if(recv_elmts[beg_elmt1 + k] != recv_elmts[beg_elmt2 + k]){
              is_same_entity = -1;
            }
          }

          if(is_same_entity == 1 ){
            same_entity_idx[idx_next_same_entity++] = idx_entity2;
          }
        }

        /* Conflict is solve save it */
        for(int k = 0; k < idx_next_same_entity; ++k) {
          int i_same_entity = same_entity_idx[k];
          int t_entity      = order[key_conflict_idx[i]+order_parent[i_same_entity]];
          recv_key[t_entity] = i_abs_entity+1;
          already_treat[i_same_entity] = 1;
        }
        i_abs_entity++;

      } /* End already_treat */
    }
  }

  PDM_free(already_treat  );
  PDM_free(same_entity_idx);
  PDM_free(tmp_parent     );
  PDM_free(order_parent   );

  PDM_free(order);
  PDM_free(key_conflict_idx);


  PDM_g_num_t current_global_num = i_abs_entity;
  PDM_g_num_t global_num_shift   = 0;
  PDM_MPI_Scan(&current_global_num, &global_num_shift, 1, PDM__PDM_MPI_G_NUM,
               PDM_MPI_SUM, gen_gnum->comm);
  global_num_shift -= current_global_num;

  for(int i = 0; i < recv_buff_idx[n_rank]; ++i) {
    recv_key[i] = recv_key[i] + global_num_shift;
  }

  /*
   * Reverse exchange
   */
  PDM_MPI_Alltoallv(recv_key,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    send_key,
                    send_buff_n,
                    send_buff_idx,
                    PDM__PDM_MPI_G_NUM,
                    gen_gnum->comm);

  /* On Stocke l'information recue */
  for (int j = 0; j < n_rank; j++) {
    send_buff_n  [j]   = 0;
  }
  for (int j = 0; j < gen_gnum->n_part; j++) {

    PDM_malloc(gen_gnum->g_nums[j], gen_gnum->n_elts[j], PDM_g_num_t);

    for (int k = 0; k < gen_gnum->n_elts[j]; k++) {
      const int t_rank = PDM_binary_search_gap_long (key_ln_to_gn[j][k]-1,
                                                     distrib,
                                                     n_rank + 1);

      int idx_read = send_buff_idx[t_rank] + send_buff_n[t_rank]++;
      gen_gnum->g_nums[j][k] = send_key[idx_read];
    }
  }

  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    PDM_free(key_ln_to_gn[i_part]);
  }
  PDM_free(key_ln_to_gn);
  PDM_free(send_key);
  PDM_free(send_elmts);
  PDM_free(recv_key);
  PDM_free(recv_elmts);

  PDM_free(send_buff_n  );
  PDM_free(send_buff_idx);
  PDM_free(recv_buff_n  );
  PDM_free(recv_buff_idx);

  PDM_free(distrib);
}

static
void
_gnum_from_comm_graph
(
 PDM_gen_gnum_t *gen_gnum,
 int             build_pcg
)
{
  PDM_part_comm_graph_t *pcg = NULL;
  if(build_pcg == 1) {
    pcg = PDM_part_comm_graph_create(gen_gnum->n_part,
                                     gen_gnum->pn_entity_graph,
                                     gen_gnum->pentity_graph,
                                     gen_gnum->comm);
  } else {
    pcg = gen_gnum->pcg;
  }

  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    PDM_malloc(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], PDM_g_num_t);
    for(int i = 0; i < gen_gnum->n_elts[i_part]; ++i) {
      gen_gnum->g_nums[i_part][i] = -1;
    }
  }

  /* Count interface entity */
  int *n_l_entity_owner_graph = malloc(gen_gnum->n_part * sizeof(int));
  int *n_l_entity_ghost_graph = malloc(gen_gnum->n_part * sizeof(int));
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {

    const int *lowner = PDM_part_comm_graph_owner_get(pcg, i_part);
    n_l_entity_owner_graph[i_part] = 0;
    n_l_entity_ghost_graph[i_part] = 0;
    int n_entity_graph = pcg->n_entity_graph[i_part];
    for(int i = 0; i < n_entity_graph; ++i) {
      int i_entity = pcg->pentity_graph[i_part][4*i  ]-1;
      if(gen_gnum->g_nums[i_part][i_entity] == -1) {
        if(lowner[i] == 1) {
          n_l_entity_owner_graph[i_part]++;
          gen_gnum->g_nums[i_part][i_entity] = -2;
        } else {
          n_l_entity_ghost_graph[i_part]++;
          gen_gnum->g_nums[i_part][i_entity] = -3;
        }
      }
    }
  }

  int *n_l_entity_interior = NULL;
  PDM_malloc(n_l_entity_interior, gen_gnum->n_part, int);
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    n_l_entity_interior[i_part] = gen_gnum->n_elts[i_part] - n_l_entity_owner_graph[i_part] - n_l_entity_ghost_graph[i_part];
  }

  /* Syncho rank */
  PDM_g_num_t l_shift = 0;
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    l_shift += n_l_entity_interior[i_part] + n_l_entity_owner_graph[i_part];
  }

  PDM_g_num_t g_shift = 0;
  PDM_MPI_Exscan(&l_shift, &g_shift, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, gen_gnum->comm);

  /* Generation ln_to_gn */
  PDM_g_num_t l_part_shift = 0;
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    int idx_write_interior = 0;
    int idx_write_graph    = 0;
    for(int i = 0; i < gen_gnum->n_elts[i_part]; ++i) {
      if(gen_gnum->g_nums[i_part][i] == -1) { // Interior
        gen_gnum->g_nums[i_part][i] = g_shift + idx_write_interior + 1;
        idx_write_interior++;
      } else if (gen_gnum->g_nums[i_part][i] == -2) {
        gen_gnum->g_nums[i_part][i] = g_shift + n_l_entity_interior[i_part] + idx_write_graph + 1;
        idx_write_graph++;
      }
    }
    l_part_shift += n_l_entity_owner_graph[i_part];
  }

  /* Echange */
  int         **send_gnum_n = NULL;
  PDM_g_num_t **send_gnum   = NULL;
  PDM_malloc(send_gnum_n, gen_gnum->n_part, int         *);
  PDM_malloc(send_gnum  , gen_gnum->n_part, PDM_g_num_t *);
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    int n_entity_graph = pcg->n_entity_graph[i_part];
    PDM_malloc(send_gnum  [i_part], n_entity_graph, PDM_g_num_t);
    PDM_malloc(send_gnum_n[i_part], n_entity_graph, int        );

    const int *lowner = PDM_part_comm_graph_owner_get(pcg, i_part);

    int idx_write = 0;
    for(int i = 0; i < n_entity_graph; ++i) {
      int i_entity = pcg->pentity_graph[i_part][4*i  ]-1;
      if(gen_gnum->g_nums[i_part][i_entity] > 0) {
        if(lowner[i] == 1) {
          send_gnum_n[i_part][i] = 1;
          send_gnum  [i_part][idx_write] = gen_gnum->g_nums[i_part][i_entity];
          idx_write++;
        } else {
          send_gnum_n[i_part][i] = 0;
        }
      } else {
        send_gnum_n[i_part][i] = 0;
      }
    }

    // PDM_log_trace_array_int (send_gnum_n[i_part], n_entity_graph, "send_gnum_n ::" );
    // PDM_log_trace_array_long(send_gnum  [i_part], idx_write, "send_gnum ::" );
  }

  PDM_g_num_t **recv_gnum   = NULL;
  int         **recv_gnum_n = NULL;
  PDM_part_comm_graph_exch(pcg,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           send_gnum_n,
              (void **)    send_gnum,
                           &recv_gnum_n,
              (void ***)   &recv_gnum);

  /*
   * Copy
   */
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    int n_entity_graph = pcg->n_entity_graph[i_part];

    if(0 == 1) {
      int n_recv = 0;
      for(int i = 0; i < n_entity_graph; ++i) {
        n_recv += recv_gnum_n[i_part][i];
      }
      PDM_log_trace_array_long(gen_gnum->g_nums[i_part], gen_gnum->n_elts[i_part], "_pentity_ln_to_gn ::");
      PDM_log_trace_array_long(recv_gnum       [i_part], n_recv                  , "recv_gnum         ::");
      PDM_log_trace_array_int (recv_gnum_n     [i_part], n_entity_graph          , "recv_gnum_n       ::");
    }

    int idx_read = 0;
    for(int i = 0; i < n_entity_graph; ++i) {
      int i_entity = pcg->pentity_graph[i_part][4*i  ]-1;
      if(recv_gnum_n[i_part][i] == 1) {
        // log_trace("i_entity = %i replace by recv_gnum[%i] = %i (before = %i) \n", i_entity, idx_read, recv_gnum[i_part][idx_read], _pentity_ln_to_gn[i_part][i_entity] );
        gen_gnum->g_nums[i_part][i_entity] = recv_gnum [i_part][idx_read];
        idx_read++;
      }
    }
  }


  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    PDM_free(send_gnum_n[i_part]);
    PDM_free(send_gnum  [i_part]);
    PDM_free(recv_gnum_n[i_part]);
    PDM_free(recv_gnum  [i_part]);
  }
  PDM_free(send_gnum_n);
  PDM_free(recv_gnum_n);
  PDM_free(send_gnum);
  PDM_free(recv_gnum);
  PDM_free(n_l_entity_interior);
  PDM_free(n_l_entity_owner_graph);
  PDM_free(n_l_entity_ghost_graph);


  if(build_pcg == 1) {
    PDM_part_comm_graph_free(pcg);
  }
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a global numbering structure
 *
 * \param [in]   dim          Spatial dimension
 * \param [in]   n_part       Number of local partitions
 * \param [in]   merge        Merge double points or not
 * \param [in]   tolerance    Geometric tolerance (used if merge double points is activated)
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to PDM_gen_gnum object
 */

PDM_gen_gnum_t *
PDM_gnum_create
(
 const int             dim,
 const int             n_part,
 const PDM_bool_t      merge,
 const double          tolerance,
 const PDM_MPI_Comm    comm,
 const PDM_ownership_t owner
)
{
  PDM_gen_gnum_t *gen_gnum;
  PDM_malloc(gen_gnum, 1, PDM_gen_gnum_t);

  gen_gnum->comm              = comm;
  gen_gnum->owner             = owner;
  gen_gnum->results_is_getted = PDM_FALSE;

  gen_gnum->n_part      = n_part;
  gen_gnum->dim         = dim;
  gen_gnum->nuplet      = 0;
  gen_gnum->merge       = merge;
  gen_gnum->tolerance   = tolerance;
  gen_gnum->n_g_elt     = -1;
  PDM_malloc(gen_gnum->g_nums, n_part, PDM_g_num_t * );
  gen_gnum->coords      = NULL;
  gen_gnum->char_length = NULL;
  gen_gnum->parent      = NULL;
  PDM_malloc(gen_gnum->n_elts, n_part, int);
  gen_gnum->index       = NULL;

  for (int i = 0; i < n_part; i++) {
    gen_gnum->g_nums[i] = NULL;
  }

  gen_gnum->pn_entity_graph = NULL;
  gen_gnum->pentity_graph   = NULL;
  gen_gnum->pcg             = NULL;

  return gen_gnum;

}


/**
 *
 * \brief Set from coordinates
 *
 * The ordering is based on a Morton code, and it is expected that
 * entities are unique (i.e. not duplicated on 2 or more ranks).
 * In the case that 2 entities have a same Morton code, their global
 * number will be determined by lexicographical ordering of coordinates.
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   coords       Coordinates (size = 3 * \ref n_elts)
 * \param [in]   char_length  Characteristic length (or NULL)
 *                            (used if merge double points is activated)
 *
 */

void
PDM_gnum_set_from_coords
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part,
 const int             n_elts,
 const double         *coords,
 const double         *char_length
)
{
  if (gen_gnum->coords == NULL) {
    PDM_malloc(gen_gnum->coords, gen_gnum->n_part, double * );
    for (int i = 0; i < gen_gnum->n_part; i++) {
      gen_gnum->coords[i_part] = NULL;
    }
  }

  if (gen_gnum->merge && gen_gnum->char_length == NULL) {
    PDM_malloc(gen_gnum->char_length, gen_gnum->n_part, double * );
    for (int i = 0; i < gen_gnum->n_part; i++) {
      gen_gnum->char_length[i_part] = NULL;
    }
  }
  gen_gnum->coords[i_part]      = (double *) coords;
  if (gen_gnum->merge) {
    gen_gnum->char_length[i_part] = (double *) char_length;
  }
  gen_gnum->n_elts[i_part]      = n_elts;
  gen_gnum->g_nums[i_part]      = NULL;

}


/**
 *
 * \brief Set Parent global numbering
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
 * \param [in]   i_part       Current partition
 * \param [in]   n_elts       Number of elements
 * \param [in]   parent_gnum  Parent global numbering (size = \ref n_elts)
 *
 */

void
PDM_gnum_set_from_parents
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part,
 const int             n_elts,
 const PDM_g_num_t    *parent_gnum
)
{

  if (gen_gnum->parent == NULL) {
    PDM_malloc(gen_gnum->parent, gen_gnum->n_part, PDM_g_num_t * );
    for (int i = 0; i < gen_gnum->n_part; i++) {
      gen_gnum->parent[i_part] = NULL;
    }
  }

  gen_gnum->parent[i_part] = (PDM_g_num_t *) parent_gnum;
  gen_gnum->n_elts[i_part] = n_elts;
  gen_gnum->g_nums[i_part] = NULL;

}

void
PDM_gnum_set_parents_nuplet
(
       PDM_gen_gnum_t  *gen_gnum,
 const int              nuplet
)
{
  gen_gnum->nuplet = nuplet;
}

void
PDM_gnum_set_from_part_comm_graph
(
       PDM_gen_gnum_t        *gen_gnum,
 const int                   *n_elts,
       PDM_part_comm_graph_t *pcg
)
{
  gen_gnum->pcg            = pcg;
  for(int i_part = 0; i_part < gen_gnum->n_part; ++i_part) {
    gen_gnum->n_elts[i_part] = n_elts[i_part];
  }
}

void
PDM_gnum_set_from_entity_graph
(
       PDM_gen_gnum_t  *gen_gnum,
 const int              i_part,
 const int              n_elts,
       int              pn_entity_graph,
       int             *pentity_graph
)
{

  if (gen_gnum->pn_entity_graph == NULL) {
    PDM_malloc(gen_gnum->pn_entity_graph, gen_gnum->n_part, int  );
    PDM_malloc(gen_gnum->pentity_graph  , gen_gnum->n_part, int *);
    for (int i = 0; i < gen_gnum->n_part; i++) {
      gen_gnum->pn_entity_graph[i_part] = 0;
      gen_gnum->pentity_graph  [i_part] = NULL;
    }
  }

  gen_gnum->n_elts         [i_part] = n_elts;
  gen_gnum->pn_entity_graph[i_part] = pn_entity_graph;
  gen_gnum->pentity_graph  [i_part] = pentity_graph;


}

/**
 *
 * \brief Compute
 *
 * \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
 *
 */

void
PDM_gnum_compute
(
 PDM_gen_gnum_t  *gen_gnum
)
{
  //Detect if geometric or topologic -- works if a procs holds no partitions
  int from_coords        = (gen_gnum->coords          != NULL);
  int from_parent        = (gen_gnum->parent          != NULL);
  int from_parent_nuplet = (gen_gnum->nuplet          != 0   );
  int from_graph         = (gen_gnum->pn_entity_graph != NULL &&
                            gen_gnum->pentity_graph   != NULL);
  int from_pcg           = (gen_gnum->pcg             != NULL);
  int from_coords_g, from_parent_g, from_parent_nuplet_g, from_graph_g, from_pcg_g;
  PDM_MPI_Allreduce(&from_coords       , &from_coords_g       , 1, PDM_MPI_INT, PDM_MPI_SUM, gen_gnum->comm);
  PDM_MPI_Allreduce(&from_parent       , &from_parent_g       , 1, PDM_MPI_INT, PDM_MPI_SUM, gen_gnum->comm);
  PDM_MPI_Allreduce(&from_parent_nuplet, &from_parent_nuplet_g, 1, PDM_MPI_INT, PDM_MPI_SUM, gen_gnum->comm);
  PDM_MPI_Allreduce(&from_graph        , &from_graph_g        , 1, PDM_MPI_INT, PDM_MPI_SUM, gen_gnum->comm);
  PDM_MPI_Allreduce(&from_pcg          , &from_pcg_g          , 1, PDM_MPI_INT, PDM_MPI_SUM, gen_gnum->comm);

  assert (from_coords_g * from_parent_g == 0);
  assert (from_pcg_g    * from_graph_g  == 0);

  if (from_coords_g != 0) {
    _gnum_from_coords_compute (gen_gnum);
  } else if (from_parent_g != 0 && from_parent_nuplet_g == 0) {
    _gnum_from_parent_compute (gen_gnum);
  } else if (from_parent_nuplet_g != 0 && // Filière "classique" mais avec optimisation ou nuplet
             from_graph_g == 0 &&
             from_pcg_g   == 0) {
    if(gen_gnum->nuplet == 1) { // Cas from parent mais optimiser
      _gnum_from_parent_compute_opt(gen_gnum);
    } else {
      _gnum_from_parent_compute_nuplet(gen_gnum);
    }
  } else if (from_graph_g != 0 || from_pcg_g != 0) {

    if (from_parent_nuplet_g != 0) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_gnum_compute, nuplet with part_graph_comm is not possible \n");
    }

    int build_pcg = 0;
    if(from_pcg_g == 0) {
      build_pcg = 1;
    }

    _gnum_from_comm_graph(gen_gnum, build_pcg);

  }
}



/**
 *
 * \brief Get global ids for a given partition
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
 * \param [in]   i_part       Current partition
 *
 * \return     Array of global ids
 *
 */

PDM_g_num_t *
PDM_gnum_get
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part
)
{
  gen_gnum->results_is_getted = PDM_TRUE;

  return gen_gnum->g_nums[i_part];
}



/**
 *
 * \brief Free
 *
 * \param [in]   gen_gnum         Pointer to \ref PDM_gen_gnum object
 *
 */

void
PDM_gnum_free
(
PDM_gen_gnum_t *gen_gnum
)
{

  if (gen_gnum->coords != NULL) {
    PDM_free(gen_gnum->coords);
  }

  if (gen_gnum->char_length != NULL) {
    PDM_free(gen_gnum->char_length);
  }

  if (gen_gnum->parent != NULL) {
    PDM_free(gen_gnum->parent);
  }

  if (gen_gnum->pn_entity_graph != NULL) {
    PDM_free(gen_gnum->pn_entity_graph);
  }
  if (gen_gnum->pentity_graph != NULL) {
    PDM_free(gen_gnum->pentity_graph);
  }

  if(( gen_gnum->owner == PDM_OWNERSHIP_KEEP ) ||
     ( gen_gnum->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !gen_gnum->results_is_getted)){
    for (int i = 0; i < gen_gnum->n_part; i++) {
      PDM_free(gen_gnum->g_nums[i]);
    }
  }

  PDM_free(gen_gnum->g_nums);
  PDM_free(gen_gnum->n_elts);

  PDM_free(gen_gnum);

}



/**
 *
 * \brief Get number of elements in a partition
 *
 * \param [in]   gen_gnum     Pointer to \ref PDM_gen_gnum object
 * \param [in]   i_part       Current partition
 *
 * \return     Number of elements
 *
 */

int
PDM_gnum_n_elt_get
(
       PDM_gen_gnum_t *gen_gnum,
 const int             i_part
)
{
  assert(gen_gnum         != NULL);
  assert(gen_gnum->n_elts != NULL);

  return gen_gnum->n_elts[i_part];
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
