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
  int               *n_part;          /*!< Number of partitions per proc in each zone */
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
_set_dFaceTag_from_joins
(
 const int          dNFace,
 const int          nJoin,
 const int         *dJoinGIds,
 const int         *dFaceJoinIdx,
 const PDM_g_num_t *dFaceJoin,
 int               *dFaceTag,
 const PDM_MPI_Comm  comm
)
{
  int iRank;
  int nRank;

  PDM_MPI_Comm_rank(comm, &iRank);
  PDM_MPI_Comm_size(comm, &nRank);

  // 1. Construct face distribution -- will be needed to find owner of faces
  PDM_g_num_t * dFaceProc = (PDM_g_num_t *) malloc((nRank+1) * sizeof(PDM_g_num_t));
  int *dNFaceProc = (int *) malloc((nRank) * sizeof(int));

  PDM_MPI_Allgather((void *) &dNFace, 1, PDM_MPI_INT, (void *) dNFaceProc, 1, PDM_MPI_INT, comm);

  dFaceProc[0] = 1;
  for (int i = 1; i < nRank+1; i++)
    dFaceProc[i] = (PDM_g_num_t) dNFaceProc[i-1] + dFaceProc[i-1];
  free(dNFaceProc);

  int  nData = 2; //Face Id, JoinOppId
  int *faceToSendN   = (int *) malloc(nRank * sizeof(int));
  int *faceToSendIdx = (int *) malloc((nRank+1) * sizeof(int));
  for (int i = 0; i < nRank; i++)
    faceToSendN[i] = 0;

  // 2. Prepare and send data
  //Count faces to send
  for (int ijoin = 0; ijoin < nJoin; ijoin++) {
    for (int iface = dFaceJoinIdx[ijoin]; iface < dFaceJoinIdx[ijoin+1]; iface++) {
      int rank = _search_rank(dFaceJoin[iface], dFaceProc, 0, nRank);
      faceToSendN[rank] += nData;
    }
  }
  //Prepare variable stride
  faceToSendIdx[0] = 0;
  for (int i = 1; i < nRank + 1; i++) {
    faceToSendIdx[i] = faceToSendIdx[i-1] + faceToSendN[i-1];
    faceToSendN[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *faceToSend = (PDM_g_num_t *) malloc(faceToSendIdx[nRank] * sizeof(PDM_g_num_t));
  for (int ijoin = 0; ijoin < nJoin; ijoin++)
  {
    for (int iface = dFaceJoinIdx[ijoin]; iface < dFaceJoinIdx[ijoin+1]; iface++)
    {
      int rank = _search_rank(dFaceJoin[iface], dFaceProc, 0, nRank);
      int idx   = faceToSendIdx[rank] + faceToSendN[rank];
      faceToSend[idx  ]   = dFaceJoin[iface];
      faceToSend[idx+1]   = dJoinGIds[2*ijoin+1];
      faceToSendN[rank] += nData;
    }
  }
  //Exchange sizes
  int *faceToRecvN   = (int *) malloc(nRank * sizeof(int));
  PDM_MPI_Alltoall(faceToSendN, 1, PDM_MPI_INT, faceToRecvN, 1, PDM_MPI_INT, comm);
  int *faceToRecvIdx = (int *) malloc((nRank+1) * sizeof(int));
  faceToRecvIdx[0] = 0;
  for(int i = 1; i < (nRank+1); i++) {
    faceToRecvIdx[i] = faceToRecvIdx[i-1] + faceToRecvN[i-1];
  }
  //Exchange data
  PDM_g_num_t *faceToRecv = (PDM_g_num_t *) malloc(faceToRecvIdx[nRank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(faceToSend,
                    faceToSendN,
                    faceToSendIdx,
                    PDM__PDM_MPI_G_NUM,
                    faceToRecv,
                    faceToRecvN,
                    faceToRecvIdx,
                    PDM__PDM_MPI_G_NUM,
                    comm);
  int nRecv = faceToRecvIdx[nRank]/nData;

  free(faceToSendN);
  free(faceToSendIdx);
  free(faceToSend);
  free(faceToRecvN);
  free(faceToRecvIdx);

  // 3. Process received data : go back to local numerotation and flag received faces
  for (int iface = 0; iface < dNFace; iface++)
    dFaceTag[iface] = -1;
  for (int i=0; i<nRecv; i++) {
    int lfaceId = faceToRecv[nData*i] - dFaceProc[iRank];
    dFaceTag[lfaceId] = faceToRecv[nData*i + 1];
  }
  free(faceToRecv);
}

static void
_rebuild_boundaries
(
 _pdm_multipart_t *_multipart
)
{

  //Set structure : we need a to retrive pboundsAndJoin for a given zone/part
  int *boundsAndJoinsIdx = (int *) malloc((_multipart->n_zone + 1) * sizeof(int));
  boundsAndJoinsIdx[0] = 0;
  for (int i = 0; i < _multipart->n_zone; i++)
    boundsAndJoinsIdx[i + 1] = _multipart->n_part[i] + boundsAndJoinsIdx[i];

  _multipart->pBoundsAndJoins = (_boundsAndJoins_t **)
  malloc(boundsAndJoinsIdx[_multipart->n_zone] * sizeof(_boundsAndJoins_t *));

  // Loop over zones and part to get data
  for (int zoneGId = 0; zoneGId<_multipart->n_zone; zoneGId++)
  {
    for (int ipart = 0; ipart < _multipart->n_part[zoneGId]; ipart++)
    {
      int nCell, nFace, nFacePartBound, nVtx, nProc, nTPart, sCellFace, sFaceVtx, sFaceGroup, nFaceGroup;
      PDM_part_part_dim_get(_multipart->partIds[zoneGId],
                        ipart,
                        &nCell,
                        &nFace,
                        &nFacePartBound,
                        &nVtx,
                        &nProc,
                        &nTPart,
                        &sCellFace,
                        &sFaceVtx,
                        &sFaceGroup,
                        &nFaceGroup);

      int nBound = _multipart->nBoundsAndJoins[2*zoneGId];
      int nJoin  = _multipart->nBoundsAndJoins[2*zoneGId+1];
      assert(nFaceGroup == nBound + nJoin);

      int          *cellTag;
      int          *cellFaceIdx;
      int          *cellFace;
      PDM_g_num_t *cellLNToGN;
      int          *faceTag;
      int          *faceCell;
      int          *faceVtxIdx;
      int          *faceVtx;
      PDM_g_num_t *faceLNToGN;
      int          *facePartBoundProcIdx;
      int          *facePartBoundPartIdx;
      int          *facePartBound;
      int          *vtxTag;
      double       *vtx;
      PDM_g_num_t  *vtxLNToGN;
      int          *faceGroupIdx;
      int          *faceGroup;
      PDM_g_num_t *faceGroupLNToGN;
      PDM_part_part_val_get(_multipart->partIds[zoneGId],
                        ipart,
                        &cellTag,
                        &cellFaceIdx,
                        &cellFace,
                        &cellLNToGN,
                        &faceTag,
                        &faceCell,
                        &faceVtxIdx,
                        &faceVtx,
                        &faceLNToGN,
                        &facePartBoundProcIdx,
                        &facePartBoundPartIdx,
                        &facePartBound,
                        &vtxTag,
                        &vtx,
                        &vtxLNToGN,
                        &faceGroupIdx,
                        &faceGroup,
                        &faceGroupLNToGN
                        );

      //Retrieve boundaries and joins from faceGroup
      int *pFaceBoundIdx = (int *) malloc((nBound+1) * sizeof(int));
      int *pFaceJoinIdx  = (int *) malloc((nJoin +1) * sizeof(int));
      for (int i = 0; i < nBound + 1; i++)
        pFaceBoundIdx[i] = faceGroupIdx[i];
      pFaceJoinIdx[0] = 0;
      for (int i = nBound + 1; i < nBound + nJoin + 1; i++)
        pFaceJoinIdx[i-nBound] = faceGroupIdx[i] - faceGroupIdx[nBound];

      int *pFaceBound = (int *) malloc(pFaceBoundIdx[nBound] * sizeof(int));
      int *pFaceJoin  = (int *) malloc(pFaceJoinIdx[nJoin]   * sizeof(int));
      for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
        pFaceBound[i] = faceGroup[i];
      for (int i = pFaceBoundIdx[nBound]; i < faceGroupIdx[nFaceGroup]; i++)
        pFaceJoin[i - pFaceBoundIdx[nBound]] = faceGroup[i];

      PDM_g_num_t *pFaceBoundLNToGN = (PDM_g_num_t *) malloc(pFaceBoundIdx[nBound] * sizeof(PDM_g_num_t));
      PDM_g_num_t *pFaceJoinLNToGN  = (PDM_g_num_t *) malloc(pFaceJoinIdx[nJoin]   * sizeof(PDM_g_num_t));
      for (int i = 0; i < pFaceBoundIdx[nBound]; i++)
        pFaceBoundLNToGN[i] = faceGroupLNToGN[i];
      for (int i = pFaceBoundIdx[nBound]; i < faceGroupIdx[nFaceGroup]; i++)
        pFaceJoinLNToGN[i - pFaceBoundIdx[nBound]] = faceGroupLNToGN[i];

      // Store data in pBoundsAndJoins
      int idx = boundsAndJoinsIdx[zoneGId] + ipart;
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

  // Implementation avec le alltoall + trie (en cours mais a jeter ?)
  int myRank;
  int nRank;
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &nRank);

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
  int *dataToSendN = (int *) malloc(nRank * sizeof(int));
  for (int i=0; i < nRank; i++)
    dataToSendN[i] = 0;

  for (int izone = 0; izone < _multipart->n_zone; izone++)
  {
    for (int ipart = 0; ipart < _multipart->n_part[izone]; ipart++)
    {
      int idx = boundsAndJoinsIdx[izone] + ipart; //TO CHECK
      int *faceJoinIdx = _multipart->pBoundsAndJoins[idx]->faceJoinIdx;
      for (int ijoin = 0; ijoin < _multipart->nBoundsAndJoins[2*izone+1]; ijoin++)
      {
        // Get destination and deduce procs that could require this data
        int joinGId = dJoinGIds[2*izone];
        int oppJoiGId = dJoinGIds[2*izone + 1];
        PDM_printf("[%i] Zone %i, ipart %i, ijoin %i (gid %i) : joinopp %i --> receiving parts are",
                   myRank, izone, ipart, ijoin, joinGId, oppJoiGId);
        for (int i = JoinToPartIdx[oppJoiGId]; i < JoinToPartIdx[oppJoiGId+1]; i++)
        {
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
  int *dataToSendIdx = (int *) malloc((nRank+1) * sizeof(int));
  dataToSendIdx[0] = 0;
  for (int i = 1; i < nRank + 1; i++) {
    dataToSendIdx[i] = dataToSendIdx[i-1] + dataToSendN[i-1];
    dataToSendN[i-1] = 0;
  }
  //Prepare data
  PDM_g_num_t *dataToSend = (PDM_g_num_t *) malloc(dataToSendIdx[nRank] * sizeof(PDM_g_num_t));
  for (int izone = 0; izone < _multipart->n_zone; izone++)
  {
    for (int ipart = 0; ipart < _multipart->n_part[izone]; ipart++)
    {
      int idx = boundsAndJoinsIdx[izone] + ipart; //TO CHECK
      int *faceJoinIdx    = _multipart->pBoundsAndJoins[idx]->faceJoinIdx;
      int *faceJoin       = _multipart->pBoundsAndJoins[idx]->faceJoin;
      int *faceJoinLNToGN = _multipart->pBoundsAndJoins[idx]->faceJoinLNToGN;
      for (int ijoin = 0; ijoin < _multipart->nBoundsAndJoins[2*izone+1]; ijoin++)
      {
        int joinGId = dJoinGIds[2*izone];
        int oppJoiGId = dJoinGIds[2*izone + 1];
        for (int i = JoinToPartIdx[oppJoiGId]; i < JoinToPartIdx[oppJoiGId+1]; i++)
        {
          int destPartition = JoinToPart[i];
          int destProc = _multipart->gPartToProc[2*destPartition];
          int idx2 = dataToSendIdx[destProc] + dataToSendN[destProc];
          int k = 0;
          for (int iface = faceJoinIdx[ijoin]; iface < faceJoinIdx[ijoin+1]; iface++)
          {
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
  int *dataToRecvN   = (int *) malloc(nRank * sizeof(int));
  PDM_MPI_Alltoall(dataToSendN, 1, PDM_MPI_INT, dataToRecvN, 1, PDM_MPI_INT, _multipart->comm);
  int *dataToRecvIdx = (int *) malloc((nRank+1) * sizeof(int));
  dataToRecvIdx[0] = 0;
  for(int i = 1; i < (nRank+1); i++) {
    dataToRecvIdx[i] = dataToRecvIdx[i-1] + dataToRecvN[i-1];
  }
  //Exchange data
  PDM_g_num_t *dataToRecv = (PDM_g_num_t *) malloc(dataToRecvIdx[nRank]*sizeof(PDM_g_num_t));
  PDM_MPI_Alltoallv(dataToSend,
                    dataToSendN,
                    dataToSendIdx,
                    PDM__PDM_MPI_G_NUM,
                    dataToRecv,
                    dataToRecvN,
                    dataToRecvIdx,
                    PDM__PDM_MPI_G_NUM,
                    _multipart->comm);
  int nRecv = dataToRecvIdx[nRank]/3;

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

  int nRank;
  PDM_MPI_Comm_size(comm, &nRank);

  // Number of partitions in each zone (distribution)
  _multipart->partZoneDistri = (int *) malloc((n_zone + 1) * sizeof(int));
  int *partZoneDistri = _multipart->partZoneDistri;
  partZoneDistri[0] = 0;

  // For each zone (slot of nRank + 1 in array), number of part per proc (distribution)
  int *dPartProc = (int *) malloc(n_zone*(nRank + 1) * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++)
  {
    dPartProc[izone*(nRank + 1)] = 0;
    PDM_MPI_Allgather((void *) &n_part[izone],
                      1,
                      PDM_MPI_INT,
                      (void *) (&dPartProc[izone*(nRank+1) + 1]),
                      1,
                      PDM_MPI_INT,
                      _multipart->comm);

    for (int i = 1; i < nRank+1; i++) {
      dPartProc[izone*(nRank+1) + i] = dPartProc[izone*(nRank+1) + i] + dPartProc[izone*(nRank+1) + i-1];
    }

    partZoneDistri[izone+1] = partZoneDistri[izone] + dPartProc[izone*(nRank+1) + nRank];
  }
  // For each global part number, owner proc and ipart in proc
  _multipart->gPartToProc = (int *) malloc(2*partZoneDistri[n_zone] * sizeof(int));
  for (int izone = 0; izone < _multipart->n_zone; izone++)
  {
    int zshift = partZoneDistri[izone];
    for (int i = 0; i < nRank; i++) {
      for (int j = dPartProc[izone*(nRank+1) + i]; j < dPartProc[izone*(nRank+1) + i+1]; j++) {
        _multipart->gPartToProc[2*(zshift + j)] = i;
        _multipart->gPartToProc[2*(zshift + j) + 1] = j - dPartProc[izone*(nRank+1) + i];

      }
    }
  }
  free(dPartProc);

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

  int iRank;
  int nRank;
  PDM_MPI_Comm_rank(_multipart->comm, &iRank);
  PDM_MPI_Comm_size(_multipart->comm, &nRank);

  if (_multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  }
  else
  {
    // 2. Loop over the blocks and call the partitionner
    for (int zoneGId = 0; zoneGId<_multipart->n_zone; zoneGId++)
    {
      PDM_printf("You requested no merge : partitionning zone %d/%d \n", zoneGId+1, _multipart->n_zone);
      int blockId = _multipart->dmeshesIds[zoneGId];
      PDM_printf("block id for zone %d is %d\n", zoneGId, blockId);
      int dNCell  = 0;
      int dNFace  = 0;
      int dNVtx   = 0;
      int nBnd    = 0;
      int nJoin   = 0;
      double       *dVtxCoord     = NULL;
      int          *dFaceVtxIdx   = NULL;
      PDM_g_num_t  *dFaceVtx      = NULL;
      PDM_g_num_t  *dFaceCell     = NULL;
      int          *dFaceBoundIdx = NULL;
      PDM_g_num_t  *dFaceBound    = NULL;
      int          *dJoinGIds     = NULL;
      int          *dFaceJoinIdx  = NULL;
      PDM_g_num_t  *dFaceJoin     = NULL;

      int nFaceGroup = 0;
      int          *dFaceGroupIdx = NULL;
      PDM_g_num_t  *dFaceGroup    = NULL;
      int          *dFaceTag      = NULL;

      if (blockId >= 0)
      {
        PDM_dmesh_dims_get(blockId, &dNCell, &dNFace, &dNVtx, &nBnd, &nJoin);
        PDM_dmesh_data_get(blockId, &dVtxCoord, &dFaceVtxIdx, &dFaceVtx, &dFaceCell,
                           &dFaceBoundIdx, &dFaceBound, &dJoinGIds, &dFaceJoinIdx, &dFaceJoin);
        //Merge FaceBounds and FaceJoins into FaceGroup
        if (dFaceJoinIdx == NULL){
          int singleArray[1] = {0};
          dFaceJoinIdx = singleArray;
        }
        nFaceGroup = nBnd + nJoin;
        dFaceGroupIdx = (int *) malloc((nFaceGroup + 1) * sizeof(int));
        dFaceGroup = (PDM_g_num_t *) malloc((dFaceBoundIdx[nBnd] + dFaceJoinIdx[nJoin]) * sizeof(PDM_g_num_t));

        for (int i=0; i < nBnd + 1; i++)
          dFaceGroupIdx[i] = dFaceBoundIdx[i];
        for (int i=0; i < dFaceBoundIdx[nBnd]; i++)
          dFaceGroup[i] = dFaceBound[i];

        for (int i=1; i < nJoin + 1; i++)
          dFaceGroupIdx[nBnd + i] = dFaceBoundIdx[nBnd] + dFaceJoinIdx[i];
        for (int i=0; i < dFaceJoinIdx[nJoin]; i++)
          dFaceGroup[dFaceBoundIdx[nBnd] + i] = dFaceJoin[i];
      }
      // Fill global array nBoundsAndJoins. nBound and nJoin are supposed to be the same for
      // procs having distributed data, so we send it to procs having no data with reduce_max
      PDM_MPI_Allreduce(&nBnd, &_multipart->nBoundsAndJoins[2*zoneGId], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);
      PDM_MPI_Allreduce(&nJoin, &_multipart->nBoundsAndJoins[2*zoneGId+1], 1,
                        PDM_MPI_INT, PDM_MPI_MAX, _multipart->comm);

      // nFaceGroup and faceGroupIdx must also be know (even if filled with 0) for every proc
      if (blockId < 0)
      {
        nFaceGroup = _multipart->nBoundsAndJoins[2*zoneGId] + _multipart->nBoundsAndJoins[2*zoneGId+1];
        dFaceGroupIdx = (int *) malloc((nFaceGroup + 1) * sizeof(int));
        for (int k=0; k < nFaceGroup + 1; k++)
          dFaceGroupIdx[k] = 0;
        dFaceCell = (PDM_g_num_t *) malloc(0); //Must be != NULL to enter in _dual_graph
      }

      dFaceTag = (int *) malloc((dNFace) * sizeof(int));
      _set_dFaceTag_from_joins(dNFace, nJoin, dJoinGIds, dFaceJoinIdx, dFaceJoin, dFaceTag, _multipart->comm);

      int ppartId = 0;
      int have_dCellPart = 0;
      int *dCellPart = (int *) malloc(dNCell*sizeof(int));
      // We probably should create a subcomm to imply only the procs sharing this block
      PDM_part_create(&ppartId,
              _multipart->comm,
              _multipart->split_method,
              "PDM_PART_RENUM_CELL_NONE",
              "PDM_PART_RENUM_FACE_NONE",
              0,                          // nPropertyCell
              NULL,                       // renum_properties_cell
              0,                          // nPropertyFace
              NULL,                       // renum_properties_face
              _multipart->n_part[zoneGId],
              dNCell,
              dNFace,
              dNVtx,
              nFaceGroup,
              NULL,                       // dCellFaceIdx
              NULL,                       // dCellFace
              NULL,                       // dCellTag
              NULL,                       // dCellWeight
              have_dCellPart,
              dCellPart,                  // dCellPart
              dFaceCell,
              dFaceVtxIdx,
              dFaceVtx,
              dFaceTag,
              dVtxCoord,
              NULL,                       // dVtxTag
              dFaceGroupIdx,
              dFaceGroup);
      PDM_printf("Partitionning done, ppardId is %d \n", ppartId);
      //Store the partition id for future access
      _multipart->partIds[zoneGId] = ppartId;

      free(dCellPart);
      free(dFaceGroupIdx);
      free(dFaceGroup);
      free(dFaceTag);
      if (blockId < 0)
      {
        free(dFaceCell);
      }
    }
    // Now separate joins and boundaries and we rebuild joins over the zones
    _rebuild_boundaries(_multipart);
  }
  // 3. rebuild the joins
}

void
PDM_multipart_part_dim_get
(
const   int  mpartId,
const   int  zoneGId,
const   int  ipart,
 int        *nCell,
 int        *nFace,
 int        *nFacePartBound,
 int        *nVtx,
 int        *nProc,
 int        *nTPart,
 int        *sCellFace,
 int        *sFaceVtx,
 int        *sFaceBound,
 int        *nFaceBound,
 int        *sFaceJoin,
 int        *nFaceJoin
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && ipart < _multipart->n_part[zoneGId]);
  int ppartId = _multipart->partIds[zoneGId];

  PDM_part_part_dim_get(ppartId,
                        ipart,
                        nCell,
                        nFace,
                        nFacePartBound,
                        nVtx,
                        nProc,
                        nTPart,
                        sCellFace,
                        sFaceVtx,
                        sFaceBound,
                        nFaceBound);
  // Get boundary and join data from pBoundsAndJoins
  int idx = 0;
  for (int i = 0; i < zoneGId; i++)
    idx += _multipart->n_part[i];
  idx += ipart;
  //Attention au cas ou pas de face de bord
  *nFaceBound = _multipart->pBoundsAndJoins[idx]->nBound;
  *sFaceBound = _multipart->pBoundsAndJoins[idx]->faceBoundIdx[*nFaceBound];
  *nFaceJoin  = _multipart->pBoundsAndJoins[idx]->nJoin;
  *sFaceJoin  = _multipart->pBoundsAndJoins[idx]->faceJoinIdx[*nFaceJoin];

}

void
PDM_multipart_part_val_get
(
const int            mpartId,
const int            zoneGId,
const int            ipart,
      int          **cellTag,
      int          **cellFaceIdx,
      int          **cellFace,
      PDM_g_num_t  **cellLNToGN,
      int          **faceTag,
      int          **faceCell,
      int          **faceVtxIdx,
      int          **faceVtx,
      PDM_g_num_t  **faceLNToGN,
      int          **facePartBoundProcIdx,
      int          **facePartBoundPartIdx,
      int          **facePartBound,
      int          **vtxTag,
      double       **vtx,
      PDM_g_num_t  **vtxLNToGN,
      int          **faceBoundIdx,
      int          **faceBound,
      PDM_g_num_t  **faceBoundLNToGN,
      int          **faceJoinIdx,
      int          **faceJoin,
      PDM_g_num_t  **faceJoinLNToGN
)
{
   _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && ipart < _multipart->n_part[zoneGId]);
  int ppartId = _multipart->partIds[zoneGId];

  PDM_part_part_val_get(ppartId,
                        ipart,
                        cellTag,
                        cellFaceIdx,
                        cellFace,
                        cellLNToGN,
                        faceTag,
                        faceCell,
                        faceVtxIdx,
                        faceVtx,
                        faceLNToGN,
                        facePartBoundProcIdx,
                        facePartBoundPartIdx,
                        facePartBound,
                        vtxTag,
                        vtx,
                        vtxLNToGN,
                        faceBoundIdx,
                        faceBound,
                        faceBoundLNToGN
                        );

  // Get boundary and join data from pBoundsAndJoins
  int idx = 0;
  for (int i = 0; i < zoneGId; i++)
    idx += _multipart->n_part[i];
  idx += ipart;
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
const int            ipart,
      int          **cellColor,
      int          **faceColor,
      int          **threadColor,
      int          **hyperPlaneColor
)
{
  _pdm_multipart_t *_multipart = _get_from_id (mpartId);

  assert(zoneGId < _multipart->n_zone && ipart < _multipart->n_part[zoneGId]);
  int ppartId = _multipart->partIds[zoneGId];

  PDM_part_part_color_get(ppartId,
                          ipart,
                          cellColor,
                          faceColor,
                          threadColor,
                          hyperPlaneColor);

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
  int ppartId = _multipart->partIds[zoneGId];

  PDM_part_time_get(ppartId,
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
