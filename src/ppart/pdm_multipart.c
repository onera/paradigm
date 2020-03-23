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
  int  nBound;
  int  nJoin;
  int *faceBoundIdx;
  int *faceJoinIdx;
  int *faceBound;
  int *faceJoin;
  PDM_g_num_t *faceBoundLNToGN;
  PDM_g_num_t *faceJoinLNToGN;

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
  int               *partZoneDistri;  /*!< Number of part in each zone (distribution,
                                           size = n_zone + 1)                         */
  int               *gPartToProc;     /*!< For each global part id, proc storing this
                                           part and localId of part in this process   */
  int               *dmeshesIds;      /*!< Ids of distributed blocks (size = n_zone)  */
  int               *partIds;         /*!< Ids of partitions built on each block of this
                                           process (size = n_zone)                    */
  int               *nBoundsAndJoins; /*!< Number of boundaries and joins in each zone
                                           (size = 2*n_zone, global data)             */
  _boundsAndJoins_t **pBoundsAndJoins;/*!< partitionned boundary and join data in each
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
 const int          nJoin,
 const int         *dJoinGIds,
 const int         *dFaceJoinIdx,
 const PDM_g_num_t *dFaceJoin,
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
  int *faceToSendN   = (int *) malloc(n_rank * sizeof(int));
  int *faceToSendIdx = (int *) malloc((n_rank+1) * sizeof(int));
  for (int i = 0; i < n_rank; i++)
    faceToSendN[i] = 0;

  // 2. Prepare and send data
  //Count faces to send
  for (int ijoin = 0; ijoin < nJoin; ijoin++) {
    for (int iface = dFaceJoinIdx[ijoin]; iface < dFaceJoinIdx[ijoin+1]; iface++) {
      int rank = _search_rank(dFaceJoin[iface], dface_proc, 0, n_rank);
      faceToSendN[rank] += nData;
    }
  }
  //Prepare variable stride
  faceToSendIdx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    faceToSendIdx[i] = faceToSendIdx[i-1] + faceToSendN[i-1];
    faceToSendN[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[n_rank] * sizeof(PDM_g_num_t));
  for (int ijoin = 0; ijoin < nJoin; ijoin++)
  {
    for (int iface = dFaceJoinIdx[ijoin]; iface < dFaceJoinIdx[ijoin+1]; iface++)
    {
      int rank = _search_rank(dFaceJoin[iface], dface_proc, 0, n_rank);
      int idx   = faceToSendIdx[rank] + faceToSendN[rank];
      faceToSend[idx  ]   = dFaceJoin[iface];
      faceToSend[idx+1]   = dJoinGIds[2*ijoin+1];
      faceToSendN[rank] += nData;
    }
  }
  //Exchange sizes
  int *faceToRecvN   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(faceToSendN, 1, PDM_MPI_INT, faceToRecvN, 1, PDM_MPI_INT, comm);
  int *faceToRecvIdx = (int *) malloc((n_rank+1) * sizeof(int));
  faceToRecvIdx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    faceToRecvIdx[i] = faceToRecvIdx[i-1] + faceToRecvN[i-1];
  }
  //Exchange data
  PDM_g_num_t *faceToRecv = (PDM_g_num_t *) malloc(faceToRecvIdx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(faceToSend,
                    faceToSendN,
                    faceToSendIdx,
                    PDM__PDM_MPI_G_NUM,
                    faceToRecv,
                    faceToRecvN,
                    faceToRecvIdx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  int nRecv = faceToRecvIdx[n_rank]/nData;

  free(faceToSendN);
  free(faceToSendIdx);
  free(faceToSend);
  free(faceToRecvN);
  free(faceToRecvIdx);

  // 3. Process received data : go back to local numerotation and flag received faces
  for (int iface = 0; iface < dn_face; iface++)
    dface_tag[iface] = -1;
  for (int i=0; i<nRecv; i++) {
    int lfaceId = faceToRecv[nData*i] - dface_proc[i_rank];
    dface_tag[lfaceId] = faceToRecv[nData*i + 1];
  }
  free(faceToRecv);
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
//   for (int zoneGId = 0; zoneGId < _multipart->n_zone; zoneGId++) {

//     // Gnum creation
//     userdom_gnum_id[zoneGId] = PDM_gnum_location_create( _multipart->n_part[zoneGId], _multipart->n_part[0], _multipart->comm);

//     for (int i_part = 0; i_part < _multipart->n_part[zoneGId]; i_part++) {

//       int n_cell, n_face, n_face_part_bound, n_vtx, n_proc, n_total_part, scell_face, sface_vtx, sface_group, n_face_group;
//       PDM_part_part_dim_get(_multipart->partIds[zoneGId],
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

//       int nBound = _multipart->nBoundsAndJoins[2*zoneGId];
//       int nJoin  = _multipart->nBoundsAndJoins[2*zoneGId+1];
//       assert(n_face_group == nBound + nJoin);

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
//       PDM_part_part_val_get(_multipart->partIds[zoneGId],
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
//       if (zoneGId == 0) {
//         PDM_gnum_location_elements_set (idtest, i_part, n_face, face_ln_to_gn);
//         PDM_g_num_t *_numabs2 = malloc(sizeof(PDM_g_num_t) * 227);
//         for (int i=0; i < 227; i++)
//           _numabs2[i] = i+1;
//         PDM_gnum_location_requested_elements_set (idtest, i_part, 227, _numabs2);
//       }

//       //Retrieve boundaries and joins from face_group
//       int *pFaceBoundIdx = (int *) malloc((nBound+1) * sizeof(int));
//       int *pFaceJoinIdx  = (int *) malloc((nJoin +1) * sizeof(int));
//       for (int i = 0; i < nBound + 1; i++)
//         pFaceBoundIdx[i] = face_group_idx[i];
//       pFaceJoinIdx[0] = 0;
//       for (int i = nBound + 1; i < nBound + nJoin + 1; i++)
//         pFaceJoinIdx[i-nBound] = face_group_idx[i] - face_group_idx[nBound];

//       int *pFaceBound = (int *) malloc(pFaceBoundIdx[nBound] * sizeof(int));
//       int *pFaceJoin  = (int *) malloc(pFaceJoinIdx[nJoin]   * sizeof(int));
//       for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
//         pFaceBound[i] = face_group[i];
//       for (int i = pFaceBoundIdx[nBound]; i < face_group_idx[n_face_group]; i++)
//         pFaceJoin[i - pFaceBoundIdx[nBound]] = face_group[i];

//       PDM_g_num_t *pFaceBoundLNToGN = (PDM_g_num_t *) malloc(pFaceBoundIdx[nBound] * sizeof(PDM_g_num_t));
//       PDM_g_num_t *pFaceJoinLNToGN  = (PDM_g_num_t *) malloc(pFaceJoinIdx[nJoin]   * sizeof(PDM_g_num_t));
//       for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
//         pFaceBoundLNToGN[i] = face_group_ln_to_gn[i];
//       for (int i = pFaceBoundIdx[nBound]; i < face_group_idx[n_face_group]; i++)
//         pFaceJoinLNToGN[i - pFaceBoundIdx[nBound]] = face_group_ln_to_gn[i];

//       // Store data in pBoundsAndJoins
//       int idx = boundsAndJoinsIdx[zoneGId] + i_part;
//       _multipart->pBoundsAndJoins[idx] = malloc(sizeof(_boundsAndJoins_t));
//       _multipart->pBoundsAndJoins[idx]->nBound = nBound;
//       _multipart->pBoundsAndJoins[idx]->nJoin = nJoin;
//       _multipart->pBoundsAndJoins[idx]->faceBoundIdx = pFaceBoundIdx;
//       _multipart->pBoundsAndJoins[idx]->faceJoinIdx  = pFaceJoinIdx;
//       _multipart->pBoundsAndJoins[idx]->faceBound    = pFaceBound;
//       _multipart->pBoundsAndJoins[idx]->faceJoin     = pFaceJoin;
//       _multipart->pBoundsAndJoins[idx]->faceBoundLNToGN    = pFaceBoundLNToGN;
//       _multipart->pBoundsAndJoins[idx]->faceJoinLNToGN     = pFaceJoinLNToGN;
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
  int *boundsAndJoinsIdx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  boundsAndJoinsIdx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++)
    boundsAndJoinsIdx[i + 1] = _multipart->n_part[i] + boundsAndJoinsIdx[i];

  _multipart->pBoundsAndJoins = (_boundsAndJoins_t **)
  malloc(boundsAndJoinsIdx[_multipart->n_zone] * sizeof(_boundsAndJoins_t *));

  // Test GNUM location
  int idtest = PDM_gnum_location_create(_multipart->n_part[0], _multipart->n_part[0], _multipart->comm);
  // Loop over zones and part to get data
  for (int zoneGId = 0; zoneGId<_multipart->n_zone; zoneGId++) {
    for (int i_part = 0; i_part < _multipart->n_part[zoneGId]; i_part++) {
      int n_cell, n_face, n_face_part_bound, n_vtx, n_proc, n_total_part, scell_face, sface_vtx, sface_group, n_face_group;
      PDM_part_part_dim_get(_multipart->partIds[zoneGId],
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

      int nBound = _multipart->nBoundsAndJoins[2*zoneGId];
      int nJoin  = _multipart->nBoundsAndJoins[2*zoneGId+1];
      assert(n_face_group == nBound + nJoin);

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
      PDM_part_part_val_get(_multipart->partIds[zoneGId],
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
      if (zoneGId == 0) {
        PDM_gnum_location_elements_set (idtest, i_part, n_face, face_ln_to_gn);
        PDM_g_num_t *_numabs2 = malloc(sizeof(PDM_g_num_t) * 227);
        for (int i=0; i < 227; i++)
          _numabs2[i] = i+1;
        PDM_gnum_location_requested_elements_set (idtest, i_part, 227, _numabs2);
      }

      //Retrieve boundaries and joins from face_group
      int *pFaceBoundIdx = (int *) malloc((nBound+1) * sizeof(int));
      int *pFaceJoinIdx  = (int *) malloc((nJoin +1) * sizeof(int));
      for (int i = 0; i < nBound + 1; i++)
        pFaceBoundIdx[i] = face_group_idx[i];
      pFaceJoinIdx[0] = 0;
      for (int i = nBound + 1; i < nBound + nJoin + 1; i++)
        pFaceJoinIdx[i-nBound] = face_group_idx[i] - face_group_idx[nBound];

      int *pFaceBound = (int *) malloc(pFaceBoundIdx[nBound] * sizeof(int));
      int *pFaceJoin  = (int *) malloc(pFaceJoinIdx[nJoin]   * sizeof(int));
      for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
        pFaceBound[i] = face_group[i];
      for (int i = pFaceBoundIdx[nBound]; i < face_group_idx[n_face_group]; i++)
        pFaceJoin[i - pFaceBoundIdx[nBound]] = face_group[i];

      PDM_g_num_t *pFaceBoundLNToGN = (PDM_g_num_t *) malloc(pFaceBoundIdx[nBound] * sizeof(PDM_g_num_t));
      PDM_g_num_t *pFaceJoinLNToGN  = (PDM_g_num_t *) malloc(pFaceJoinIdx[nJoin]   * sizeof(PDM_g_num_t));
      for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
        pFaceBoundLNToGN[i] = face_group_ln_to_gn[i];
      for (int i = pFaceBoundIdx[nBound]; i < face_group_idx[n_face_group]; i++)
        pFaceJoinLNToGN[i - pFaceBoundIdx[nBound]] = face_group_ln_to_gn[i];

      // Store data in pBoundsAndJoins
      int idx = boundsAndJoinsIdx[zoneGId] + i_part;
      _multipart->pBoundsAndJoins[idx] = malloc(sizeof(_boundsAndJoins_t));
      _multipart->pBoundsAndJoins[idx]->nBound = nBound;
      _multipart->pBoundsAndJoins[idx]->nJoin = nJoin;
      _multipart->pBoundsAndJoins[idx]->faceBoundIdx = pFaceBoundIdx;
      _multipart->pBoundsAndJoins[idx]->faceJoinIdx  = pFaceJoinIdx;
      _multipart->pBoundsAndJoins[idx]->faceBound    = pFaceBound;
      _multipart->pBoundsAndJoins[idx]->faceJoin     = pFaceJoin;
      _multipart->pBoundsAndJoins[idx]->faceBoundLNToGN    = pFaceBoundLNToGN;
      _multipart->pBoundsAndJoins[idx]->faceJoinLNToGN     = pFaceJoinLNToGN;

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

  //Step 0 ; construction de joinGId -> liste partitions partageant ce join en num globale (apres PT)
  // TODO -> ASSUME WE HAVE IT FOR NOW
  // Local ou pas local ? La zone est la même pour tt les parts partageant le join
  int *JoinToPartIdx = (int *) malloc(3 * sizeof(int));
  JoinToPartIdx[0] = 0;
  JoinToPartIdx[1] = 1;
  JoinToPartIdx[2] = 3;
  int *JoinToPart = (int *) malloc(JoinToPartIdx[2] * sizeof(int));
  JoinToPart[0] = 0; //zone 0 part 0   (zone 0 had 2 parts)
  JoinToPart[1] = 2; //zone 1 part 0   (zone 1 had 3 parts)
  JoinToPart[2] = 4; //zone 1 part 2   (zone 1 had 3 parts)

  // ASSUME we have the array dJoinGIds : for each zone, gives joinId, joinOppId
  int *dJoinGIds = (int *) malloc(2*_multipart->n_zone * sizeof(int));
  dJoinGIds[2*0] = 0;
  dJoinGIds[2*0+1] = 1;
  dJoinGIds[2*1] = 1;
  dJoinGIds[2*1+1] = 0;

  // Step 1. Count data
  int *dataToSendN = (int *) malloc(n_rank * sizeof(int));
  for (int i=0; i < n_rank; i++)
    dataToSendN[i] = 0;

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = boundsAndJoinsIdx[izone] + i_part; //TO CHECK
      int *faceJoinIdx = _multipart->pBoundsAndJoins[idx]->faceJoinIdx;
      for (int ijoin = 0; ijoin < _multipart->nBoundsAndJoins[2*izone+1]; ijoin++) {
        // Get destination and deduce procs that could require this data
        int joinGId = dJoinGIds[2*izone];
        int oppJoiGId = dJoinGIds[2*izone + 1];
        PDM_printf("[%i] Zone %i, i_part %i, ijoin %i (gid %i) : joinopp %i --> receiving parts are",
                   i_rank, izone, i_part, ijoin, joinGId, oppJoiGId);
        for (int i = JoinToPartIdx[oppJoiGId]; i < JoinToPartIdx[oppJoiGId+1]; i++) {
          int destPartition = JoinToPart[i];
          int destProc = _multipart->gPartToProc[2*destPartition];

          PDM_printf(" %d (proc %d)", destPartition, _multipart->gPartToProc[2*destPartition]);
          //We have the destination, exchanged data is 3 times the lenght of point list
          // (pl value, LNToGN value, joinGId value)
          dataToSendN[destProc] += 3*(faceJoinIdx[ijoin+1] - faceJoinIdx[ijoin]);
        }
        PDM_printf("\n");

      }
    }
  }
  // Step 2. Prepare data and performs alltoall
  // Prepare stride
  int *dataToSendIdx = (int *) malloc((n_rank+1) * sizeof(int));
  dataToSendIdx[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    dataToSendIdx[i] = dataToSendIdx[i-1] + dataToSendN[i-1];
    dataToSendN[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *dataToSend = (PDM_g_num_t *) malloc(dataToSendIdx[n_rank] * sizeof(PDM_g_num_t));
  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    for (int i_part = 0; i_part < _multipart->n_part[izone]; i_part++) {
      int idx = boundsAndJoinsIdx[izone] + i_part; //TO CHECK
      int *faceJoinIdx    = _multipart->pBoundsAndJoins[idx]->faceJoinIdx;
      int *faceJoin       = _multipart->pBoundsAndJoins[idx]->faceJoin;
      PDM_g_num_t *faceJoinLNToGN = _multipart->pBoundsAndJoins[idx]->faceJoinLNToGN;
      for (int ijoin = 0; ijoin < _multipart->nBoundsAndJoins[2*izone+1]; ijoin++) {
        int joinGId = dJoinGIds[2*izone];
        int oppJoiGId = dJoinGIds[2*izone + 1];
        for (int i = JoinToPartIdx[oppJoiGId]; i < JoinToPartIdx[oppJoiGId+1]; i++) {
          int destPartition = JoinToPart[i];
          int destProc = _multipart->gPartToProc[2*destPartition];
          int idx2 = dataToSendIdx[destProc] + dataToSendN[destProc];
          int k = 0;
          for (int iface = faceJoinIdx[ijoin]; iface < faceJoinIdx[ijoin+1]; iface++) {
            dataToSend[idx2 + 3*k    ] = faceJoin[iface];
            dataToSend[idx2 + 3*k + 1] = faceJoinLNToGN[iface];
            dataToSend[idx2 + 3*k + 2] = joinGId;
            k += 1;
          }
          dataToSendN[destProc] += 3*k;
        }
      }
    }
  }
  //Exchange sizes
  int *dataToRecvN   = (int *) malloc(n_rank * sizeof(int));
  PDM_MPI_Alltoall(dataToSendN, 1, PDM_MPI_INT, dataToRecvN, 1, PDM_MPI_INT, _multipart->comm);
  int *dataToRecvIdx = (int *) malloc((n_rank+1) * sizeof(int));
  dataToRecvIdx[0] = 0;
  for(int i = 1; i < (n_rank+1); i++) {
    dataToRecvIdx[i] = dataToRecvIdx[i-1] + dataToRecvN[i-1];
  }
  //Exchange data
  PDM_g_num_t *dataToRecv = (PDM_g_num_t *) malloc(dataToRecvIdx[n_rank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(dataToSend,
                    dataToSendN,
                    dataToSendIdx,
                    PDM__PDM_MPI_G_NUM,
                    dataToRecv,
                    dataToRecvN,
                    dataToRecvIdx,
                    PDM__PDM_MPI_G_NUM,
                    _multipart->comm);

  int nRecv = dataToRecvIdx[n_rank]/3;

  // Step 3. Search in received data the matching faces

  free(boundsAndJoinsIdx);
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

  _multipart->dmeshesIds      = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->partIds         = (int *) malloc(_multipart->n_zone * sizeof(int));
  _multipart->nBoundsAndJoins = (int *) malloc(_multipart->n_zone * 2 * sizeof(int));

  for (int izone = 0; izone < _multipart->n_zone; izone++) {
    _multipart->dmeshesIds[izone] = -1;
    _multipart->partIds   [izone] = -1;
    _multipart->nBoundsAndJoins[2*izone]   = -1;
    _multipart->nBoundsAndJoins[2*izone+1] = -1;
  }

  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);

  // Number of partitions in each zone (distribution)
  _multipart->partZoneDistri = (int *) malloc((n_zone + 1) * sizeof(int));
  int *partZoneDistri = _multipart->partZoneDistri;
  partZoneDistri[0] = 0;

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

    partZoneDistri[izone+1] = partZoneDistri[izone] + dpart_proc[izone*(n_rank+1) + n_rank];
  }
  // For each global part number, owner proc and i_part in proc
  _multipart->gPartToProc = (int *) malloc(2*partZoneDistri[n_zone] * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++){
    int zshift = partZoneDistri[izone];
    for (int i = 0; i < n_rank; i++) {
      for (int j = dpart_proc[izone*(n_rank+1) + i]; j < dpart_proc[izone*(n_rank+1) + i+1]; j++) {
        _multipart->gPartToProc[2*(zshift + j)] = i;
        _multipart->gPartToProc[2*(zshift + j) + 1] = j - dpart_proc[izone*(n_rank+1) + i];

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
 const int        zoneGId,
 const int        block_data_id
)
{
  PDM_printf("In multipart %d, set zone n°%d using blockdata %d \n",
             mpart_id, zoneGId, block_data_id);

  _pdm_multipart_t *_multipart = _get_from_id (mpart_id);
  assert(zoneGId < _multipart->n_zone);
  _multipart->dmeshesIds[zoneGId] = block_data_id;
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
    for (int zoneGId = 0; zoneGId < _multipart->n_zone; zoneGId++) {
      PDM_printf("You requested no merge : partitionning zone %d/%d \n", zoneGId+1, _multipart->n_zone);
      int blockId = _multipart->dmeshesIds[zoneGId];
      PDM_printf("block id for zone %d is %d\n", zoneGId, blockId);
      int dn_cell  = 0;
      int dn_face  = 0;
      int dn_vtx   = 0;
      int nBnd    = 0;
      int nJoin   = 0;
      const double       *dvtx_coord;
      const int          *dface_vtx_idx;
      const PDM_g_num_t  *dface_vtx;
      const PDM_g_num_t  *dface_cell;
      const int          *dFaceBoundIdx;
      const PDM_g_num_t  *dFaceBound;
      const int          *dJoinGIds;
      const int          *dFaceJoinIdx;
      const PDM_g_num_t  *dFaceJoin;

      int n_face_group = 0;
      int          *dface_group_idx = NULL;
      PDM_g_num_t  *dface_group    = NULL;
      int          *dface_tag      = NULL;

      if (blockId >= 0)
      {
        PDM_dmesh_dims_get(blockId, &dn_cell, &dn_face, &dn_vtx, &nBnd, &nJoin);
        PDM_dmesh_data_get(blockId, &dvtx_coord, &dface_vtx_idx, &dface_vtx, &dface_cell,
                           &dFaceBoundIdx, &dFaceBound, &dJoinGIds, &dFaceJoinIdx, &dFaceJoin);
        //Merge FaceBounds and FaceJoins into face_group
        if (dFaceJoinIdx == NULL){
          int singleArray[1] = {0};
          dFaceJoinIdx = singleArray;
        }
        n_face_group = nBnd + nJoin;
        dface_group_idx = (int *) malloc((n_face_group + 1) * sizeof(int));
        dface_group = (PDM_g_num_t *) malloc((dFaceBoundIdx[nBnd] + dFaceJoinIdx[nJoin]) * sizeof(PDM_g_num_t));

        for (int i=0; i < nBnd + 1; i++)
          dface_group_idx[i] = dFaceBoundIdx[i];
        for (int i=0; i < dFaceBoundIdx[nBnd]; i++)
          dface_group[i] = dFaceBound[i];

        for (int i=1; i < nJoin + 1; i++)
          dface_group_idx[nBnd + i] = dFaceBoundIdx[nBnd] + dFaceJoinIdx[i];
        for (int i=0; i < dFaceJoinIdx[nJoin]; i++)
          dface_group[dFaceBoundIdx[nBnd] + i] = dFaceJoin[i];
      }
      // Fill global array nBoundsAndJoins. nBound and nJoin are supposed to be the same for
      // procs having distributed data, so we send it to procs having no data with reduce_max
      PDM_MPI_Allreduce(&nBnd, &_multipart->nBoundsAndJoins[2*zoneGId], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);
      PDM_MPI_Allreduce(&nJoin, &_multipart->nBoundsAndJoins[2*zoneGId+1], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);

      // n_face_group and face_group_idx must also be know (even if filled with 0) for every proc
      if (blockId < 0)
      {
        n_face_group = _multipart->nBoundsAndJoins[2*zoneGId] + _multipart->nBoundsAndJoins[2*zoneGId+1];
        dface_group_idx = (int *) malloc((n_face_group + 1) * sizeof(int));
        for (int k=0; k < n_face_group + 1; k++)
          dface_group_idx[k] = 0;
        dface_cell = (PDM_g_num_t *) malloc(0); //Must be != NULL to enter in _dual_graph
      }

      dface_tag = (int *) malloc((dn_face) * sizeof(int));
      _set_dface_tag_from_joins(dn_face, nJoin, dJoinGIds, dFaceJoinIdx, dFaceJoin, dface_tag, _multipart->comm);

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
              _multipart->n_part[zoneGId],
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
      _multipart->partIds[zoneGId] = ppart_id;

      free(dcell_part);
      free(dface_group_idx);
      free(dface_group);
      free(dface_tag);
      if (blockId < 0)
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
const   int  mpartId,
const   int  zoneGId,
const   int  i_part,
 int        *n_cell,
 int        *n_face,
 int        *n_face_part_bound,
 int        *n_vtx,
 int        *n_proc,
 int        *n_total_part,
 int        *scell_face,
 int        *sface_vtx,
 int        *sFaceBound,
 int        *n_faceBound,
 int        *sFaceJoin,
 int        *n_faceJoin
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && i_part < _multipart->n_part[zoneGId]);
  int ppart_id = _multipart->partIds[zoneGId];

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
                        sFaceBound,
                        n_faceBound);
  // Get boundary and join data from pBoundsAndJoins
  int idx = 0;
  for (int i = 0; i < zoneGId; i++)
    idx += _multipart->n_part[i];
  idx += i_part;
  //Attention au cas ou pas de face de bord
  *n_faceBound = _multipart->pBoundsAndJoins[idx]->nBound;
  *sFaceBound = _multipart->pBoundsAndJoins[idx]->faceBoundIdx[*n_faceBound];
  *n_faceJoin  = _multipart->pBoundsAndJoins[idx]->nJoin;
  *sFaceJoin  = _multipart->pBoundsAndJoins[idx]->faceJoinIdx[*n_faceJoin];

}

void
PDM_multipart_part_val_get
(
const int            mpartId,
const int            zoneGId,
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
      int          **faceBoundIdx,
      int          **faceBound,
      PDM_g_num_t  **faceBoundLNToGN,
      int          **faceJoinIdx,
      int          **faceJoin,
      PDM_g_num_t  **faceJoinLNToGN
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && i_part < _multipart->n_part[zoneGId]);
  int ppart_id = _multipart->partIds[zoneGId];

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
                        faceBoundIdx,
                        faceBound,
                        faceBoundLNToGN
                        );

  // Get boundary and join data from pBoundsAndJoins
  int idx = 0;
  for (int i = 0; i < zoneGId; i++)
    idx += _multipart->n_part[i];
  idx += i_part;
  //Attention au cas ou pas de face de bord
  *faceBoundIdx       = _multipart->pBoundsAndJoins[idx]->faceBoundIdx;
  *faceBound          = _multipart->pBoundsAndJoins[idx]->faceBound;
  *faceBoundLNToGN    = _multipart->pBoundsAndJoins[idx]->faceBoundLNToGN;
  *faceJoinIdx        = _multipart->pBoundsAndJoins[idx]->faceJoinIdx;
  *faceJoin           = _multipart->pBoundsAndJoins[idx]->faceJoin;
  *faceJoinLNToGN     = _multipart->pBoundsAndJoins[idx]->faceJoinLNToGN;

}

void
PDM_multipart_part_color_get
(
const int            mpartId,
const int            zoneGId,
const int            i_part,
      int          **cell_color,
      int          **face_color,
      int          **thread_color,
      int          **hyperplane_color
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && i_part < _multipart->n_part[zoneGId]);
  int ppart_id = _multipart->partIds[zoneGId];

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
const int       mpartId,
const int       zoneGId,
      double  **elapsed,
      double  **cpu,
      double  **cpu_user,
      double  **cpu_sys
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);
  assert(zoneGId < _multipart->n_zone);
  int ppart_id = _multipart->partIds[zoneGId];

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

  free(_multipart->dmeshesIds);
  free(_multipart->nBoundsAndJoins);
  free(_multipart->partZoneDistri);
  free(_multipart->gPartToProc);

  for (int izone = 0; izone<_multipart->n_zone; izone++)
    PDM_part_free(_multipart->partIds[izone]);
  free(_multipart->partIds);

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
