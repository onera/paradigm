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
 * TODO : write module description here
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_handles.h"
#include "pdm_dmesh.h"
#include "pdm_gnum_location.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart.h"

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

typedef struct  {
  int  n_bound;
  int  n_join;
  int *face_bound_idx;
  int *face_join_idx;
  int *face_bound;
  int *face_join;
  PDM_g_num_t *face_bound_ln_to_gn;
  PDM_g_num_t *face_join_ln_to_gn;

} _boundsAndJoins_t;

/**
 * \struct _pdm_multipart_t
 * \brief  This structure describe a multipart. In addition to splitting
 *         parameters, it stores the multiples blocks and part as well as
 *         the global numbering.
 *
 */

typedef struct  {

  int               n_zone;           /*!< Number of initial zones */
  const int        *n_part;          /*!< Number of partitions per proc in each zone */
  PDM_bool_t        merge_blocks;     /*!< Merge before partitionning or not */
  PDM_part_split_t  split_method;     /*!< Partitioning method */
  PDM_MPI_Comm      comm;             /*!< MPI communicator */
  int               *part_zone_distribution;  /*!< Number of part in each zone (distribution,
                                           size = n_zone + 1)                         */
  int               *gpart_to_proc;     /*!< For each global part id, proc storing this
                                           part and localId of part in this process   */
  int               *dmeshes_ids;      /*!< Ids of distributed blocks (size = n_zone)  */
  int               *part_ids;         /*!< Ids of partitions built on each block of this
                                           process (size = n_zone)                    */
  int               *n_bounds_and_joins; /*!< Number of boundaries and joins in each zone
                                           (size = 2*n_zone, global data)             */
  _boundsAndJoins_t **pbounds_and_joins;/*!< partitionned boundary and join data in each
                                           zone/part                                  */
} _pdm_multipart_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_multiparts   = NULL;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return multipart object from it identifier
 *
 * \param [in]   multipartId    multipart identifier
 *
 */

static _pdm_multipart_t *
_get_from_id
(
 int  id
)
{

  _pdm_multipart_t *multipart = (_pdm_multipart_t *) PDM_Handles_get (_multiparts, id);

  if (multipart == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_multipart error : Bad identifier\n");
  }

  return multipart;
}

static int
 _search_rank
(
 PDM_g_num_t   elt,
 PDM_g_num_t  *array,
 int            id1,
 int            id2
)
{
  if (elt >= array[id2]) {
    PDM_printf("PPART error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
           elt, array[id1], array[id2]);
    abort();
  }

  if (elt < array[id1]) {
    PDM_printf("PPART error : Element not in initial distributed array "
           PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
           elt, array[id1], array[id2]);
    abort();
  }

  if (id2 == id1 + 1) {
    return id1;
  }

  else {

    while(array[id1] == array[id1+1]) id1++;

    int midId = (id2 + id1) / 2;

    if (elt == array[id1])
      return id1;
    else if (elt < array[midId])
      return _search_rank(elt, array, id1, midId);
    else if (elt >= array[midId])
      return _search_rank(elt, array, midId, id2);
  }
  return -1;
}

static void
_set_dface_tag_from_joins
(
 const int          dn_face,
 const int          n_join,
 const int         *djoin_gids,
 const int         *dface_join_idx,
 const PDM_g_num_t *dface_join,
 int               *dface_tag,
 const PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // 1. Construct face distribution -- will be needed to find owner of faces
  PDM_g_num_t * dface_proc = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  int *dn_face_proc = (int *) malloc((n_rank) * sizeof(int));

  PDM_MPI_Allgather((void *) &dn_face, 1, PDM_MPI_INT, (void *) dn_face_proc, 1, PDM_MPI_INT, comm);

  dface_proc[0] = 1;
  for (int i = 1; i < n_rank+1; i++)
    dface_proc[i] = (PDM_g_num_t) dn_face_proc[i-1] + dface_proc[i-1];
  free(dn_face_proc);

  int  nData = 2; //Face Id, JoinOppId
  int *face_to_send_n   = (int *) malloc(n_rank * sizeof(int));
  int *face_to_send_idx = (int *) malloc((n_rank+1) * sizeof(int));
  for (int i = 0; i < n_rank; i++)
    face_to_send_n[i] = 0;

  // 2. Prepare and send data
  //Count faces to send
  for (int ijoin = 0; ijoin < n_join; ijoin++) {
    for (int iface = dface_join_idx[ijoin]; iface < dface_join_idx[ijoin+1]; iface++) {
      int rank = _search_rank(dface_join[iface], dface_proc, 0, n_rank);
      face_to_send_n[rank] += nData;
    }
  }
  //Prepare variable stride
  face_to_send_idx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    face_to_send_idx[i] = face_to_send_idx[i-1] + face_to_send_n[i-1];
    face_to_send_n[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *face_to_send = (PDM_g_num_t *) malloc(face_to_send_idx[n_rank] * sizeof(PDM_g_num_t));
  for (int ijoin = 0; ijoin < n_join; ijoin++)
  {
    for (int iface = dface_join_idx[ijoin]; iface < dface_join_idx[ijoin+1]; iface++)
    {
      int rank = _search_rank(dface_join[iface], dface_proc, 0, n_rank);
      int idx   = face_to_send_idx[rank] + face_to_send_n[rank];
      face_to_send[idx  ]   = dface_join[iface];
      face_to_send[idx+1]   = djoin_gids[2*ijoin+1];
      face_to_send_n[rank] += nData;
    }
  }
  //Exchange sizes
  int *face_to_recv_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(face_to_send_n, 1, PDM_MPI_INT, face_to_recv_n, 1, PDM_MPI_INT, comm);
  int *face_to_recv_idx = (int *) malloc((n_rank+1) * sizeof(int));
  face_to_recv_idx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    face_to_recv_idx[i] = face_to_recv_idx[i-1] + face_to_recv_n[i-1];
  }
  //Exchange data
  PDM_g_num_t *face_to_recv = (PDM_g_num_t *) malloc(face_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(face_to_send,
                    face_to_send_n,
                    face_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    face_to_recv,
                    face_to_recv_n,
                    face_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  int nRecv = face_to_recv_idx[n_rank]/nData;

  free(face_to_send_n);
  free(face_to_send_idx);
  free(face_to_send);
  free(face_to_recv_n);
  free(face_to_recv_idx);

  // 3. Process received data : go back to local numerotation and flag received faces
  for (int iface = 0; iface < dn_face; iface++)
    dface_tag[iface] = -1;
  for (int i=0; i<nRecv; i++) {
    int lfaceId = face_to_recv[nData*i] - dface_proc[i_rank];
    dface_tag[lfaceId] = face_to_recv[nData*i + 1];
  }
  free(face_to_recv);
}


// static void
// _rebuild_boundaries2
// (
//  _pdm_multipart_t *_multipart
// )
// {
//   /*
//    * I/ For each block we build the gnum location
//    */
//   int* userdom_gnum_id = (int *) malloc(_multipart->n_zone * sizeof(int*) );
//   for (int zone_gid = 0; zone_gid < _multipart->n_zone; zone_gid++) {

//     // Gnum creation
//     userdom_gnum_id[zone_gid] = PDM_gnum_location_create( _multipart->n_part[zone_gid], _multipart->n_part[0], _multipart->comm);

//     for (int i_part = 0; i_part < _multipart->n_part[zone_gid]; i_part++) {

//       int n_cell, n_face, n_face_part_bound, n_vtx, n_proc, n_total_part, scell_face, sface_vtx, sface_group, n_face_group;
//       PDM_part_part_dim_get(_multipart->part_ids[zone_gid],
//                         i_part,
//                         &n_cell,
//                         &n_face,
//                         &n_face_part_bound,
//                         &n_vtx,
//                         &n_proc,
//                         &n_total_part,
//                         &scell_face,
//                         &sface_vtx,
//                         &sface_group,
//                         &n_face_group);

//       int n_bound = _multipart->n_bounds_and_joins[2*zone_gid];
//       int n_join  = _multipart->n_bounds_and_joins[2*zone_gid+1];
//       assert(n_face_group == n_bound + n_join);

//       int          *cell_tag;
//       int          *cell_face_idx;
//       int          *cell_face;
//       PDM_g_num_t *cell_ln_to_gn;
//       int          *face_tag;
//       int          *face_cell;
//       int          *face_vtx_idx;
//       int          *face_vtx;
//       PDM_g_num_t *face_ln_to_gn;
//       int          *face_part_bound_proc_idx;
//       int          *face_part_bound_part_idx;
//       int          *face_part_bound;
//       int          *vtx_tag;
//       double       *vtx;
//       PDM_g_num_t  *vtx_ln_to_gn;
//       int          *face_group_idx;
//       int          *face_group;
//       PDM_g_num_t *face_group_ln_to_gn;
//       PDM_part_part_val_get(_multipart->part_ids[zone_gid],
//                         i_part,
//                         &cell_tag,
//                         &cell_face_idx,
//                         &cell_face,
//                         &cell_ln_to_gn,
//                         &face_tag,
//                         &face_cell,
//                         &face_vtx_idx,
//                         &face_vtx,
//                         &face_ln_to_gn,
//                         &face_part_bound_proc_idx,
//                         &face_part_bound_part_idx,
//                         &face_part_bound,
//                         &vtx_tag,
//                         &vtx,
//                         &vtx_ln_to_gn,
//                         &face_group_idx,
//                         &face_group,
//                         &face_group_ln_to_gn
//                         );
//       //Test gnum location
//       if (zone_gid == 0) {
//         PDM_gnum_location_elements_set (idtest, i_part, n_face, face_ln_to_gn);
//         PDM_g_num_t *_numabs2 = malloc(sizeof(PDM_g_num_t) * 227);
//         for (int i=0; i < 227; i++)
//           _numabs2[i] = i+1;
//         PDM_gnum_location_requested_elements_set (idtest, i_part, 227, _numabs2);
//       }

//       //Retrieve boundaries and joins from face_group
//       int *pface_bound_idx = (int *) malloc((n_bound+1) * sizeof(int));
//       int *pface_join_idx  = (int *) malloc((n_join +1) * sizeof(int));
//       for (int i = 0; i < n_bound + 1; i++)
//         pface_bound_idx[i] = face_group_idx[i];
//       pface_join_idx[0] = 0;
//       for (int i = n_bound + 1; i < n_bound + n_join + 1; i++)
//         pface_join_idx[i-n_bound] = face_group_idx[i] - face_group_idx[n_bound];

//       int *pface_bound = (int *) malloc(pface_bound_idx[n_bound] * sizeof(int));
//       int *pface_join  = (int *) malloc(pface_join_idx[n_join]   * sizeof(int));
//       for (int i = 0; i < pface_bound_idx[n_bound]; i++)
//         pface_bound[i] = face_group[i];
//       for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
//         pface_join[i - pface_bound_idx[n_bound]] = face_group[i];

//       PDM_g_num_t *pface_bound_ln_to_gn = (PDM_g_num_t *) malloc(pface_bound_idx[n_bound] * sizeof(PDM_g_num_t));
//       PDM_g_num_t *pface_join_ln_to_gn  = (PDM_g_num_t *) malloc(pface_join_idx[n_join]   * sizeof(PDM_g_num_t));
//       for (int i = 0; i < pface_bound_idx[n_bound]; i++)
//         pface_bound_ln_to_gn[i] = face_group_ln_to_gn[i];
//       for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
//         pface_join_ln_to_gn[i - pface_bound_idx[n_bound]] = face_group_ln_to_gn[i];

//       // Store data in pbounds_and_joins
//       int idx = bounds_and_joins_idx[zone_gid] + i_part;
//       _multipart->pbounds_and_joins[idx] = malloc(sizeof(_boundsAndJoins_t));
//       _multipart->pbounds_and_joins[idx]->n_bound = n_bound;
//       _multipart->pbounds_and_joins[idx]->n_join = n_join;
//       _multipart->pbounds_and_joins[idx]->face_bound_idx = pface_bound_idx;
//       _multipart->pbounds_and_joins[idx]->face_join_idx  = pface_join_idx;
//       _multipart->pbounds_and_joins[idx]->face_bound    = pface_bound;
//       _multipart->pbounds_and_joins[idx]->face_join     = pface_join;
//       _multipart->pbounds_and_joins[idx]->face_bound_ln_to_gn    = pface_bound_ln_to_gn;
//       _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn     = pface_join_ln_to_gn;
//     }
//   }
//   free(userdom_gnum_id);
// }

static void
_rebuild_boundaries
(
 _pdm_multipart_t *_multipart
)
{

  printf("_rebuild_boundaries::\n");

  //Set structure : we need a to retrive pboundsAndJoin for a given zone/part
  int *bounds_and_joins_idx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  bounds_and_joins_idx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++)
    bounds_and_joins_idx[i + 1] = _multipart->n_part[i] + bounds_and_joins_idx[i];

  _multipart->pbounds_and_joins = (_boundsAndJoins_t **)
  malloc(bounds_and_joins_idx[_multipart->n_zone] * sizeof(_boundsAndJoins_t *));

  // Test GNUM location
  int idtest = PDM_gnum_location_create(_multipart->n_part[0], _multipart->n_part[0], _multipart->comm);
  // Loop over zones and part to get data
  for (int zone_gid = 0; zone_gid<_multipart->n_zone; zone_gid++) {
    for (int i_part = 0; i_part < _multipart->n_part[zone_gid]; i_part++) {
      int n_cell, n_face, n_face_part_bound, n_vtx, n_proc, n_total_part, scell_face, sface_vtx, sface_group, n_face_group;
      PDM_part_part_dim_get(_multipart->part_ids[zone_gid],
                        i_part,
                        &n_cell,
                        &n_face,
                        &n_face_part_bound,
                        &n_vtx,
                        &n_proc,
                        &n_total_part,
                        &scell_face,
                        &sface_vtx,
                        &sface_group,
                        &n_face_group);

      int n_bound = _multipart->n_bounds_and_joins[2*zone_gid];
      int n_join  = _multipart->n_bounds_and_joins[2*zone_gid+1];
      assert(n_face_group == n_bound + n_join);

      int          *cell_tag;
      int          *cell_face_idx;
      int          *cell_face;
      PDM_g_num_t *cell_ln_to_gn;
      int          *face_tag;
      int          *face_cell;
      int          *face_vtx_idx;
      int          *face_vtx;
      PDM_g_num_t *face_ln_to_gn;
      int          *face_part_bound_proc_idx;
      int          *face_part_bound_part_idx;
      int          *face_part_bound;
      int          *vtx_tag;
      double       *vtx;
      PDM_g_num_t  *vtx_ln_to_gn;
      int          *face_group_idx;
      int          *face_group;
      PDM_g_num_t *face_group_ln_to_gn;
      PDM_part_part_val_get(_multipart->part_ids[zone_gid],
                        i_part,
                        &cell_tag,
                        &cell_face_idx,
                        &cell_face,
                        &cell_ln_to_gn,
                        &face_tag,
                        &face_cell,
                        &face_vtx_idx,
                        &face_vtx,
                        &face_ln_to_gn,
                        &face_part_bound_proc_idx,
                        &face_part_bound_part_idx,
                        &face_part_bound,
                        &vtx_tag,
                        &vtx,
                        &vtx_ln_to_gn,
                        &face_group_idx,
                        &face_group,
                        &face_group_ln_to_gn
                        );
      //Test gnum location
      if (zone_gid == 0) {
        PDM_gnum_location_elements_set (idtest, i_part, n_face, face_ln_to_gn);
        PDM_g_num_t *_numabs2 = malloc(sizeof(PDM_g_num_t) * 227);
        for (int i=0; i < 227; i++)
          _numabs2[i] = i+1;
        PDM_gnum_location_requested_elements_set (idtest, i_part, 227, _numabs2);
      }

      //Retrieve boundaries and joins from face_group
      int *pface_bound_idx = (int *) malloc((n_bound+1) * sizeof(int));
      int *pface_join_idx  = (int *) malloc((n_join +1) * sizeof(int));
      for (int i = 0; i < n_bound + 1; i++)
        pface_bound_idx[i] = face_group_idx[i];
      pface_join_idx[0] = 0;
      for (int i = n_bound + 1; i < n_bound + n_join + 1; i++)
        pface_join_idx[i-n_bound] = face_group_idx[i] - face_group_idx[n_bound];

      int *pface_bound = (int *) malloc(pface_bound_idx[n_bound] * sizeof(int));
      int *pface_join  = (int *) malloc(pface_join_idx[n_join]   * sizeof(int));
      for (int i = 0; i < pface_bound_idx[n_bound]; i++)
        pface_bound[i] = face_group[i];
      for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
        pface_join[i - pface_bound_idx[n_bound]] = face_group[i];

      PDM_g_num_t *pface_bound_ln_to_gn = (PDM_g_num_t *) malloc(pface_bound_idx[n_bound] * sizeof(PDM_g_num_t));
      PDM_g_num_t *pface_join_ln_to_gn  = (PDM_g_num_t *) malloc(pface_join_idx[n_join]   * sizeof(PDM_g_num_t));
      for (int i = 0; i < pface_bound_idx[n_bound]; i++)
        pface_bound_ln_to_gn[i] = face_group_ln_to_gn[i];
      for (int i = pface_bound_idx[n_bound]; i < face_group_idx[n_face_group]; i++)
        pface_join_ln_to_gn[i - pface_bound_idx[n_bound]] = face_group_ln_to_gn[i];

      // Store data in pbounds_and_joins
      int idx = bounds_and_joins_idx[zone_gid] + i_part;
      _multipart->pbounds_and_joins[idx] = malloc(sizeof(_boundsAndJoins_t));
      _multipart->pbounds_and_joins[idx]->n_bound         = n_bound;
      _multipart->pbounds_and_joins[idx]->n_join          = n_join;
      _multipart->pbounds_and_joins[idx]->face_bound_idx    = pface_bound_idx;
      _multipart->pbounds_and_joins[idx]->face_join_idx     = pface_join_idx;
      _multipart->pbounds_and_joins[idx]->face_bound       = pface_bound;
      _multipart->pbounds_and_joins[idx]->face_join        = pface_join;
      _multipart->pbounds_and_joins[idx]->face_bound_ln_to_gn = pface_bound_ln_to_gn;
      _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn  = pface_join_ln_to_gn;

    }
  }
  //Test du gnum location
  PDM_gnum_location_compute (idtest);
  int *location_idx;
  int *location;

  PDM_gnum_location_get(idtest, 0, &location_idx, &location);
  PDM_gnum_location_free (idtest, 1);








  // Implementation avec le alltoall + trie (en cours mais a jeter ?)
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  //Step 0 ; construction de join_gid -> liste partitions partageant ce join en num globale (apres PT)
  // TODO -> ASSUME WE HAVE IT FOR NOW
  // Local ou pas local ? La zone est la même pour tt les parts partageant le join
  int *join_to_part_idx = (int *) malloc(3 * sizeof(int));
  join_to_part_idx[0] = 0;
  join_to_part_idx[1] = 1;
  join_to_part_idx[2] = 3;
  int *join_to_part = (int *) malloc(join_to_part_idx[2] * sizeof(int));
  join_to_part[0] = 0; //zone 0 part 0   (zone 0 had 2 parts)
  join_to_part[1] = 2; //zone 1 part 0   (zone 1 had 3 parts)
  join_to_part[2] = 4; //zone 1 part 2   (zone 1 had 3 parts)

  // ASSUME we have the array djoin_gids : for each zone, gives joinId, joinOppId
  int *djoin_gids = (int *) malloc(2*_multipart->n_zone * sizeof(int));
  djoin_gids[2*0] = 0;
  djoin_gids[2*0+1] = 1;
  djoin_gids[2*1] = 1;
  djoin_gids[2*1+1] = 0;

  // Step 1. Count data
  int *data_to_send_n = (int *) malloc(n_rank * sizeof(int));
  for (int i=0; i < n_rank; i++)
    data_to_send_n[i] = 0;

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int *face_join_idx = _multipart->pbounds_and_joins[idx]->face_join_idx;
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone+1]; ijoin++) {
        // Get destination and deduce procs that could require this data
        int join_gid = djoin_gids[2*izone];
        int opp_join_gid = djoin_gids[2*izone + 1];
        PDM_printf("[%i] Zone %i, i_part %i, ijoin %i (gid %i) : joinopp %i --> receiving parts are",
                   i_rank, izone, i_part, ijoin, join_gid, opp_join_gid);
        for (int i = join_to_part_idx[opp_join_gid]; i < join_to_part_idx[opp_join_gid+1]; i++) {
          int destPartition = join_to_part[i];
          int destProc = _multipart->gpart_to_proc[2*destPartition];

          PDM_printf(" %d (proc %d)", destPartition, _multipart->gpart_to_proc[2*destPartition]);
          //We have the destination, exchanged data is 3 times the lenght of point list
          // (pl value, LNToGN value, join_gid value)
          data_to_send_n[destProc] += 3*(face_join_idx[ijoin+1] - face_join_idx[ijoin]);
        }
        PDM_printf("\n");

      }
    }
  }
  // Step 2. Prepare data and performs alltoall
  // Prepare stride
  int *data_to_send_idx = (int *) malloc((n_rank+1) * sizeof(int));
  data_to_send_idx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    data_to_send_idx[i] = data_to_send_idx[i-1] + data_to_send_n[i-1];
    data_to_send_n[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *data_to_send = (PDM_g_num_t *) malloc(data_to_send_idx[n_rank] * sizeof(PDM_g_num_t));
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = bounds_and_joins_idx[izone] + i_part; //TO CHECK
      int *face_join_idx    = _multipart->pbounds_and_joins[idx]->face_join_idx;
      int *face_join       = _multipart->pbounds_and_joins[idx]->face_join;
      PDM_g_num_t *face_join_ln_to_gn = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;
      for (int ijoin = 0; ijoin < _multipart->n_bounds_and_joins[2*izone+1]; ijoin++) {
        int join_gid     = djoin_gids[2*izone];
        int opp_join_gid = djoin_gids[2*izone + 1];
        for (int i = join_to_part_idx[opp_join_gid]; i < join_to_part_idx[opp_join_gid+1]; i++) {
          int destPartition = join_to_part[i];
          int destProc = _multipart->gpart_to_proc[2*destPartition];
          int idx2 = data_to_send_idx[destProc] + data_to_send_n[destProc];
          int k = 0;
          for (int iface = face_join_idx[ijoin]; iface < face_join_idx[ijoin+1]; iface++) {
            data_to_send[idx2 + 3*k    ] = face_join[iface];
            data_to_send[idx2 + 3*k + 1] = face_join_ln_to_gn[iface];
            data_to_send[idx2 + 3*k + 2] = join_gid;
            k += 1;
          }
          data_to_send_n[destProc] += 3*k;
        }
      }
    }
  }
  //Exchange sizes
  int *data_to_recv_n   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(data_to_send_n, 1, PDM_MPI_INT, data_to_recv_n, 1, PDM_MPI_INT, _multipart->comm);
  int *data_to_recv_idx = (int *) malloc((n_rank+1) * sizeof(int));
  data_to_recv_idx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    data_to_recv_idx[i] = data_to_recv_idx[i-1] + data_to_recv_n[i-1];
  }
  //Exchange data
  PDM_g_num_t *data_to_recv = (PDM_g_num_t *) malloc(data_to_recv_idx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(data_to_send,
                    data_to_send_n,
                    data_to_send_idx,
                    PDM__PDM_MPI_G_NUM,
                    data_to_recv,
                    data_to_recv_n,
                    data_to_recv_idx,
                    PDM__PDM_MPI_G_NUM,
                    _multipart->comm);

  int nRecv = data_to_recv_idx[n_rank]/3;

  // Step 3. Search in received data the matching faces

  free(bounds_and_joins_idx);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure
 *
 * \param [in]   n_zone       Number of zones in the original mesh
 * \param [in]   n_part       Number of partition per proc in each zone
 * \param [in]   merge_blocks Merge or not the zones before splitting
 * \param [in]   split_method Choice of library used to split the mesh
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Identifier
 */

int
PDM_multipart_create
(
 const int              n_zone,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_part_split_t split_method,
 const PDM_MPI_Comm     comm
)
{
  printf("PDM_multipart_create::n_zone:: %d \n", n_zone);
  printf("PDM_multipart_create::n_part:: %d \n", n_part[0]);
  printf("PDM_multipart_create::split_method:: %d \n", split_method);
  /*
   * Search a ppart free id
   */

  if (_multiparts == NULL) {
    _multiparts = PDM_Handles_create (4);
  }

  _pdm_multipart_t *_multipart = (_pdm_multipart_t *) malloc(sizeof(_pdm_multipart_t));
  int id = PDM_Handles_store (_multiparts, _multipart);

  _multipart->n_zone      = n_zone;
  _multipart->n_part      = n_part;
  _multipart->merge_blocks= merge_blocks;
  _multipart->split_method= split_method;
  _multipart->comm        = comm;

  _multipart->dmeshes_ids      = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->part_ids         = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->n_bounds_and_joins = (int *) malloc(_multipart->n_zone * 2 * sizeof(int));

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    _multipart->dmeshes_ids[izone] = -1;
    _multipart->part_ids   [izone] = -1;
    _multipart->n_bounds_and_joins[2*izone]   = -1;
    _multipart->n_bounds_and_joins[2*izone+1] = -1;
  }

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  // Number of partitions in each zone (distribution)
  _multipart->part_zone_distribution = (int *) malloc((n_zone + 1) * sizeof(int));
  int *part_zone_distribution = _multipart->part_zone_distribution;
  part_zone_distribution[0] = 0;

  // For each zone (slot of n_rank + 1 in array), number of part per proc (distribution)
  int *dpart_proc = (int *) malloc(n_zone*(n_rank + 1) * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    dpart_proc[izone*(n_rank + 1)] = 0;
    PDM_MPI_Allgather((void *) &n_part[izone],
                      1,
                      PDM_MPI_INT,
                      (void *) (&dpart_proc[izone*(n_rank+1) + 1]),
                      1,
                      PDM_MPI_INT,
                      _multipart->comm);

    for (int i = 1; i < n_rank+1; i++) {
      dpart_proc[izone*(n_rank+1) + i] = dpart_proc[izone*(n_rank+1) + i] + dpart_proc[izone*(n_rank+1) + i-1];
    }

    part_zone_distribution[izone+1] = part_zone_distribution[izone] + dpart_proc[izone*(n_rank+1) + n_rank];
  }
  // For each global part number, owner proc and i_part in proc
  _multipart->gpart_to_proc = (int *) malloc(2*part_zone_distribution[n_zone] * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++){
    int zshift = part_zone_distribution[izone];
    for (int i = 0; i < n_rank; i++) {
      for (int j = dpart_proc[izone*(n_rank+1) + i]; j < dpart_proc[izone*(n_rank+1) + i+1]; j++) {
        _multipart->gpart_to_proc[2*(zshift + j)] = i;
        _multipart->gpart_to_proc[2*(zshift + j) + 1] = j - dpart_proc[izone*(n_rank+1) + i];

      }
    }
  }
  free(dpart_proc);

  PDM_printf("Created from PDM_multipart_create. You requested a multipart with %d zones \n", n_zone);
  return id;

}

/* TODO : copy doc of the function */
void PDM_multipart_register_block
(
 const int        mpart_id,
 const int        zone_gid,
 const int        block_data_id
)
{
  PDM_printf("In multipart %d, set zone n°%d using blockdata %d \n",
             mpart_id, zone_gid, block_data_id);

  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_gid < _multipart->n_zone);
  _multipart->dmeshes_ids[zone_gid] = block_data_id;
}

void
PDM_multipart_run_ppart
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(_multipart->comm, &i_rank);
  PDM_MPI_Comm_size(_multipart->comm, &n_rank);

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  }
  else
  {
    // 2. Loop over the blocks and call the partitionner
    for (int zone_gid = 0; zone_gid < _multipart->n_zone; zone_gid++) {
      PDM_printf("You requested no merge : partitionning zone %d/%d \n", zone_gid+1, _multipart->n_zone);
      int block_id = _multipart->dmeshes_ids[zone_gid];
      PDM_printf("block id for zone %d is %d\n", zone_gid, block_id);
      int dn_cell  = 0;
      int dn_face  = 0;
      int dn_vtx   = 0;
      int nBnd    = 0;
      int n_join   = 0;
      const double       *dvtx_coord;
      const int          *dface_vtx_idx;
      const PDM_g_num_t  *dface_vtx;
      const PDM_g_num_t  *dface_cell;
      const int          *dface_bound_idx;
      const PDM_g_num_t  *dface_bound;
      const int          *djoin_gids;
      const int          *dface_join_idx;
      const PDM_g_num_t  *dface_join;

      int n_face_group = 0;
      int          *dface_group_idx = NULL;
      PDM_g_num_t  *dface_group    = NULL;
      int          *dface_tag      = NULL;

      if (block_id >= 0)
      {
        PDM_dmesh_dims_get(block_id, &dn_cell, &dn_face, &dn_vtx, &nBnd, &n_join);
        PDM_dmesh_data_get(block_id, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                           &dface_bound_idx, &dface_bound, &djoin_gids, &dface_join_idx, &dface_join);
        //Merge face_bounds and face_joins into face_group
        if (dface_join_idx == NULL){
          int single_array[1] = {0};
          dface_join_idx = single_array;
        }
        n_face_group = nBnd + n_join;
        dface_group_idx = (int *) malloc((n_face_group + 1) * sizeof(int));
        dface_group = (PDM_g_num_t *) malloc((dface_bound_idx[nBnd] + dface_join_idx[n_join]) * sizeof(PDM_g_num_t));

        for (int i=0; i < nBnd + 1; i++)
          dface_group_idx[i] = dface_bound_idx[i];
        for (int i=0; i < dface_bound_idx[nBnd]; i++)
          dface_group[i] = dface_bound[i];

        for (int i=1; i < n_join + 1; i++)
          dface_group_idx[nBnd + i] = dface_bound_idx[nBnd] + dface_join_idx[i];
        for (int i=0; i < dface_join_idx[n_join]; i++)
          dface_group[dface_bound_idx[nBnd] + i] = dface_join[i];
      }
      // Fill global array n_bounds_and_joins. n_bound and n_join are supposed to be the same for
      // procs having distributed data, so we send it to procs having no data with reduce_max
      PDM_MPI_Allreduce(&nBnd, &_multipart->n_bounds_and_joins[2*zone_gid], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);
      PDM_MPI_Allreduce(&n_join, &_multipart->n_bounds_and_joins[2*zone_gid+1], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);

      // n_face_group and face_group_idx must also be know (even if filled with 0) for every proc
      if (block_id < 0)
      {
        n_face_group = _multipart->n_bounds_and_joins[2*zone_gid] + _multipart->n_bounds_and_joins[2*zone_gid+1];
        dface_group_idx = (int *) malloc((n_face_group + 1) * sizeof(int));
        for (int k=0; k < n_face_group + 1; k++)
          dface_group_idx[k] = 0;
        dface_cell = (PDM_g_num_t *) malloc(0); //Must be != NULL to enter in _dual_graph
      }

      dface_tag = (int *) malloc((dn_face) * sizeof(int));
      _set_dface_tag_from_joins(dn_face, n_join, djoin_gids, dface_join_idx, dface_join, dface_tag, _multipart->comm);

      int ppart_id = 0;
      int have_dcell_part = 0;
      int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
      // We probably should create a subcomm to imply only the procs sharing this block
      PDM_part_create(&ppart_id,
              _multipart->comm,
              _multipart->split_method,
              "PDM_PART_RENUM_CELL_NONE",
              "PDM_PART_RENUM_FACE_NONE",
              0,                          // n_property_cell
              NULL,                       // renum_properties_cell
              0,                          // n_property_face
              NULL,                       // renum_properties_face
              _multipart->n_part[zone_gid],
              dn_cell,
              dn_face,
              dn_vtx,
              n_face_group,
              NULL,                       // dcell_faceIdx
              NULL,                       // dcell_face
              NULL,                       // dcell_tag
              NULL,                       // dcell_weight
              have_dcell_part,
              dcell_part,                  // dcell_part
              dface_cell,
              dface_vtx_idx,
              dface_vtx,
              dface_tag,
              dvtx_coord,
              NULL,                       // dvtx_tag
              dface_group_idx,
              dface_group);
      PDM_printf("Partitionning done, ppardId is %d \n", ppart_id);
      //Store the partition id for future access
      _multipart->part_ids[zone_gid] = ppart_id;

      free(dcell_part);
      free(dface_group_idx);
      free(dface_group);
      free(dface_tag);
      if (block_id < 0)
      {
        free(dface_cell); //FIXME: Revoir la gestion du rien
      }
    }
    // Now separate joins and boundaries and we rebuild joins over the zones
    _rebuild_boundaries(_multipart);
  }
  // 3. rebuild the joins
  // _rebuild_boundaries2(_multipart);
}

void
PDM_multipart_part_dim_get
(
const   int  mpart_id,
const   int  zone_gid,
const   int  i_part,
 int        *n_cell,
 int        *n_face,
 int        *n_face_part_bound,
 int        *n_vtx,
 int        *n_proc,
 int        *n_total_part,
 int        *scell_face,
 int        *sface_vtx,
 int        *sface_bound,
 int        *n_face_bound,
 int        *sface_join,
 int        *n_face_join
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_dim_get(ppart_id,
                        i_part,
                        n_cell,
                        n_face,
                        n_face_part_bound,
                        n_vtx,
                        n_proc,
                        n_total_part,
                        scell_face,
                        sface_vtx,
                        sface_bound,
                        n_face_bound);
  // Get boundary and join data from pbounds_and_joins
  int idx = 0;
  for (int i = 0; i < zone_gid; i++)
    idx += _multipart->n_part[i];
  idx += i_part;
  //Attention au cas ou pas de face de bord
  *n_face_bound = _multipart->pbounds_and_joins[idx]->n_bound;
  *sface_bound = _multipart->pbounds_and_joins[idx]->face_bound_idx[*n_face_bound];
  *n_face_join  = _multipart->pbounds_and_joins[idx]->n_join;
  *sface_join  = _multipart->pbounds_and_joins[idx]->face_join_idx[*n_face_join];

}

void
PDM_multipart_part_val_get
(
const int            mpart_id,
const int            zone_gid,
const int            i_part,
      int          **cell_tag,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_tag,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      int          **vtx_tag,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_bound_idx,
      int          **face_bound,
      PDM_g_num_t  **face_bound_ln_to_gn,
      int          **face_join_idx,
      int          **face_join,
      PDM_g_num_t  **face_join_ln_to_gn
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_val_get(ppart_id,
                        i_part,
                        cell_tag,
                        cell_face_idx,
                        cell_face,
                        cell_ln_to_gn,
                        face_tag,
                        face_cell,
                        face_vtx_idx,
                        face_vtx,
                        face_ln_to_gn,
                        face_part_bound_proc_idx,
                        face_part_bound_part_idx,
                        face_part_bound,
                        vtx_tag,
                        vtx,
                        vtx_ln_to_gn,
                        face_bound_idx,
                        face_bound,
                        face_bound_ln_to_gn
                        );

  // Get boundary and join data from pbounds_and_joins
  int idx = 0;
  for (int i = 0; i < zone_gid; i++)
    idx += _multipart->n_part[i];
  idx += i_part;
  //Attention au cas ou pas de face de bord
  *face_bound_idx       = _multipart->pbounds_and_joins[idx]->face_bound_idx;
  *face_bound           = _multipart->pbounds_and_joins[idx]->face_bound;
  *face_bound_ln_to_gn  = _multipart->pbounds_and_joins[idx]->face_bound_ln_to_gn;
  *face_join_idx        = _multipart->pbounds_and_joins[idx]->face_join_idx;
  *face_join            = _multipart->pbounds_and_joins[idx]->face_join;
  *face_join_ln_to_gn   = _multipart->pbounds_and_joins[idx]->face_join_ln_to_gn;

}

void
PDM_multipart_part_color_get
(
const int            mpart_id,
const int            zone_gid,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **thread_color,
      int          **hyperplane_color
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);

  assert(zone_gid < _multipart->n_zone && i_part < _multipart->n_part[zone_gid]);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_part_color_get(ppart_id,
                          i_part,
                          cell_color,
                          face_color,
                          thread_color,
                          hyperplane_color);

}

void
PDM_multipart_time_get
(
const int       mpart_id,
const int       zone_gid,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zone_gid < _multipart->n_zone);
  int ppart_id = _multipart->part_ids[zone_gid];

  PDM_part_time_get(ppart_id,
                    elapsed,
                    cpu,
                    cpu_user,
                    cpu_sys);

}

/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */

void
PDM_multipart_free
(
 const int id
)
{
  _pdm_multipart_t *_multipart = _get_from_id (id);

  free(_multipart->dmeshes_ids);
  free(_multipart->n_bounds_and_joins);
  free(_multipart->part_zone_distribution);
  free(_multipart->gpart_to_proc);

  for (int izone = 0; izone<_multipart->n_zone; izone++)
    PDM_part_free(_multipart->part_ids[izone]);
  free(_multipart->part_ids);

  free (_multipart);

  PDM_Handles_handle_free (_multiparts, id, PDM_FALSE);

  const int n_multipart = PDM_Handles_n_get (_multiparts);

  if (n_multipart == 0) {
    _multiparts = PDM_Handles_free (_multiparts);
  }
  PDM_printf("Cleaned from PDM_multipart_free\n");
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
