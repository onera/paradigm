/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_part_bound.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_surf_part_t structure
 *
 * This function returns an initialized \ref PDM_surf_part_t structure
 *
 * \param [in]  nFace       Number of faces
 * \param [in]  faceVtxIdx  Index in the face -> vertex connectivity
 * \param [in]  faceVtxIdx  face -> vertex connectivity
 * \param [in]  faceLnToGn  Local face numbering to global face numbering
 * \param [in]  nVtx        Number of vertices
 * \param [in]  coords      Coordinates
 * \param [in]  vtxLnToGn   Local vertex numbering to global vertex numbering
 *
 * \return      A new initialized \ref PDM_surf_part_t structure
 *
 */

PDM_surf_part_t *
PDM_surf_part_create
(
const int         nFace,
const int        *faceVtxIdx,
const int        *faceVtx,
const PDM_g_num_t *faceLnToGn,
const int         nVtx,
const double     *coords,
const PDM_g_num_t *vtxLnToGn
)
{
  int vb=1;
  if (vb == 1)   PDM_printf ("==== PDM_surf_part_create ====\n");
  PDM_surf_part_t *_part = (PDM_surf_part_t *) malloc(sizeof(PDM_surf_part_t));

  _part->nFace                = nFace;
  _part->nGhostFace           = 0;
  _part->nTotalFace           = nFace;
  _part->sFaceVtx             = faceVtxIdx[nFace];
  _part->faceVtxIdx           = faceVtxIdx;
  _part->faceVtx              = faceVtx;
  _part->faceEdgeIdx          = NULL;
  _part->faceEdge             = NULL;
  _part->faceLnToGn           = faceLnToGn;
  _part->nVtx                 = nVtx;
  _part->coords               = coords;
  _part->vtxEdgeIdx           = NULL;
  _part->vtxEdge              = NULL;
  _part->vtxLnToGn            = vtxLnToGn;
  _part->nEdge                = -1;
  _part->nGhostEdge           = -1;
  _part->nTotalEdge           = -1;
  _part->edgeFace             = NULL;
  _part->edgeVtx              = NULL;
  _part->edgeLnToGn           = NULL;
  _part->edgePartBound        = NULL;
  _part->vtxPartBound         = NULL;
  _part->carLgthVtx           = NULL;
  _part->faceNormal           = NULL;
  _part->extents              = NULL;

  if (vb == 1)   PDM_printf ("==== PDM_surf_part_create ==== terminated ====\n");
  return _part;
}


/**
 * \brief Delete a \ref PDM_surf_part_t structure
 *
 * This function deletes a  PDM_surf_part_t structure
 *
 * \param [in]  part      part to delete
 *
 * \return     Null pointer
 */

PDM_surf_part_t *
PDM_surf_part_free
(
 PDM_surf_part_t * part
)
{
  int vb=1;
  if (vb == 1)   PDM_printf ("==== PDM_surf_part_free ==== \n");
  assert (part != NULL);

  if (part != NULL) {
    part->faceVtxIdx = NULL;
    part->faceVtx = NULL;
    if (part->faceEdgeIdx != NULL)
      free(part->faceEdgeIdx);
    if (part->faceEdge != NULL)
      free(part->faceEdge);

    part->faceLnToGn = NULL;
    part->coords = NULL;
    if (part->vtxEdgeIdx != NULL)
      free(part->vtxEdgeIdx);
    if (part->vtxEdge != NULL)
      free(part->vtxEdge);
    part->vtxLnToGn = NULL;

    if (part->edgeFace != NULL)
      free(part->edgeFace);

    if (part->edgeVtx != NULL)
      free(part->edgeVtx);

    if (part->edgeLnToGn != NULL)
      free(part->edgeLnToGn);

    if (part->edgePartBound != NULL)
      part->edgePartBound = PDM_part_bound_free(part->edgePartBound);
    if (part->vtxPartBound != NULL)
      part->vtxPartBound = PDM_part_bound_free(part->vtxPartBound);

    if (part->carLgthVtx != NULL)
      free (part->carLgthVtx);

    if (part->faceNormal != NULL)
      free (part->faceNormal);

    if (part->extents != NULL)
      free (part->extents);

    free(part);
  }
  if (vb == 1)   PDM_printf ("==== PDM_surf_part_free ==== terminated ====\n");

  return NULL;
}


/**
 * \brief Compute partition edge entities
 *
 * This function defines edges of an initial partitiob and
 * computes edge connectivities
 *
 * \param [in]  _part      Partition to compute
 *
 */

void
PDM_surf_part_build_edges
(
PDM_surf_part_t *part
)
{

  assert (part != NULL);
  int vb = 1;
  if (vb == 1)   PDM_printf ("==== PDM_surf_part_build_edges ==== \n");
  vb = 0;
  /*
   * build hash table with key = vertices sum and list of edges
   */

  int  lHashTableIdx = 2 * part->nVtx + 1;
  int *hashTableIdx  = (int *) malloc(sizeof(int) * lHashTableIdx);
  int *hashTable     = (int *) malloc(sizeof(int) * part->faceVtxIdx[part->nFace]);
  int *nHashTable    = (int *) malloc(sizeof(int) * 2 * part->nVtx);

  for (int i = 0; i <  part->faceVtxIdx[part->nFace]; i++) {
    hashTable[i] = -1;
  }

  for (int i = 0; i < lHashTableIdx; i++) {
    hashTableIdx[i] = 0;
  }

  for (int i = 0; i < 2 * part->nVtx; i++) {
    nHashTable[i] = 0;
  }

  int l_edges = 0;
  for (int i = 0; i < part->nFace; i++) {
    int idxFace = part->faceVtxIdx[i];
    for (int j = idxFace; j < part->faceVtxIdx[i+1]; j++) {
      int k = (j == part->faceVtxIdx[i+1] - 1) ? idxFace : (j + 1);
      int s_vtx = part->faceVtx[j] + part->faceVtx[k];
      hashTableIdx[s_vtx+1] += 1;
/*PDM_printf("i: %d  j: %d, k: %d, s_vtx:%d+%d: %d hashTableIdx[%d] : %d\n", i, j, k, part->faceVtx[j], part->faceVtx[k], s_vtx, s_vtx+1, hashTableIdx[s_vtx+1]); */
      l_edges += 2;
    }
  }

  hashTableIdx[0] = 0;
  for (int i = 1; i < lHashTableIdx; i++) {
    hashTableIdx[i] =  hashTableIdx[i] +  hashTableIdx[i-1];
  }
  if (vb == 1) {
    PDM_printf("hashTableIdx :\n");
    for (int i = 0; i < lHashTableIdx; i++) {
      PDM_printf("%d: %d  ", i, hashTableIdx[i]);
    }
    PDM_printf("\n");
    PDM_printf("l_edges : %d\n", l_edges);
  }

  int *listEdges = (int *) malloc(sizeof(int) * l_edges);
  int *edgeFaceUncompress = (int *) malloc(sizeof(int) * l_edges/2);
  l_edges = 0;

  for (int i = 0; i < part->nFace; i++) {
    int idxFace = part->faceVtxIdx[i];
    int nVtxFace = part->faceVtxIdx[i+1] - idxFace;
    for (int j = idxFace; j < idxFace + nVtxFace; j++) {
      int k = (j == part->faceVtxIdx[i+1] - 1) ? idxFace : (j + 1);
      int vtx1 = part->faceVtx[j];
      int vtx2 = part->faceVtx[k];
      int s_vtx = vtx1 + vtx2;

      hashTable[hashTableIdx[s_vtx] + (nHashTable[s_vtx]++)] = l_edges/2;
      edgeFaceUncompress[l_edges/2] = i;
/* PDM_printf("i: %d  idxFace: %d, nVtxFace: %d, j: %d, k: %d, vtx1: %d, vtx2: %d, s_vtx: %d hashTable[%d+%d] : %d, edgeFaceUncompress[%d]: %d\n", i, idxFace, nVtxFace, j, k, vtx1, vtx2, s_vtx, hashTableIdx[s_vtx] , (nHashTable[s_vtx]-1), l_edges/2, l_edges/2, i);  */
      listEdges[l_edges++] = vtx1;
      listEdges[l_edges++] = vtx2;
    }
  }
   if (vb == 1) {
    PDM_printf("nHashTable :\n");
    for (int i = 0; i < 2*part->nVtx; i++)
      PDM_printf("%d: %d   ", i, nHashTable[i]);
    PDM_printf("\n");
    PDM_printf("l_edges : %d\n", l_edges);
    PDM_printf("edgeFaceUncompress :\n");
    for (int i = 0; i < l_edges/2; i++)
      PDM_printf("%d: %d   ", i, edgeFaceUncompress[i]);
    PDM_printf("\n");
    PDM_printf("listEdges :\n");
    for (int i = 0; i < l_edges; i++) {
      PDM_printf("%d: %d   ", i, listEdges[i]);
    }
    PDM_printf("\n");
    PDM_printf("hashtable :\n");

    for (int i = 0; i < 2 * part->nVtx; i++) {
      PDM_printf("[%d] : ", i);
      for (int k = hashTableIdx[i]; k < hashTableIdx[i+1]; k++)
        PDM_printf(" %d", hashTable[k]);
      PDM_printf("\n");
    }
	}

  free(nHashTable);

  /*
   * Compress edges
   */

  int nEdgesFaces = l_edges/2;
  int *edgesToCompressEdges = (int *) malloc(sizeof(int) * nEdgesFaces);
  for (int i = 0; i < nEdgesFaces; i++) {
    edgesToCompressEdges[i] = -1;
  }

  /*
   * First allocation to l_edges
   */

  part->edgeVtx  = (int *) malloc(sizeof(int) * l_edges);
  part->nEdge    = 0;

  /*
   * Compute real edges
   */

  int nEdge = 0;

  part->vtxEdgeIdx  = (int *) malloc(sizeof(int) * (part->nVtx + 1));
  for (int i = 0; i < part->nVtx + 1; i++) {
    part->vtxEdgeIdx[i] = 0;
  }

  for (int i = 0; i < 2 * part->nVtx; i++) {
    int idx          = hashTableIdx[i];
    int nEdgeSameSum = (hashTableIdx[i+1] - idx);
    for (int j = idx; j < idx + nEdgeSameSum; j++) {
      int iEdge = hashTable[j];
/*PDM_printf("i: %d, j: %d, iEdge: %d ", i, j,iEdge); */
      if (iEdge != -1) {
/*PDM_printf("edgesToCompressEdges[%d] = %d  ", iEdge, nEdge); */
        edgesToCompressEdges[iEdge] = nEdge;
        int vtx1 = listEdges[2*iEdge];
        int vtx2 = listEdges[2*iEdge+1];
        part->edgeVtx[2*nEdge]      = vtx1;
        part->edgeVtx[2*nEdge + 1]  = vtx2;
        part->vtxEdgeIdx[vtx1] += 1;
        part->vtxEdgeIdx[vtx2] += 1;
        for (int k = j + 1; k < idx + nEdgeSameSum; k++) {
          int iEdge1 = hashTable[k];
          if (iEdge1 != -1) {
            int vtx11 = listEdges[2*iEdge1];
            int vtx12 = listEdges[2*iEdge1+1];

            if (((vtx1 == vtx11) && (vtx2 == vtx12)) ||
                ((vtx2 == vtx11) && (vtx1 == vtx12))) {
              hashTable[k] = -1;
/*PDM_printf("i: %d, j: %d, k: %d, edgesToCompressEdges[%d]: %d\n", i, j, k, iEdge1, nEdge); */

              edgesToCompressEdges[iEdge1] = nEdge;
            }
          }
        }
        nEdge += 1;
/*PDM_printf("\n"); */
      }
    }
  }

if (0 == 1) {
    PDM_printf ("nEdge: %d nEdgesFaces: %d \n", nEdge,nEdgesFaces);
    PDM_printf ("edgesToCompressEdges \n");
    for (int i = 0; i < nEdgesFaces; i++) {
		PDM_printf("edgesToCompressEdges[%d]: %d\n", i, edgesToCompressEdges[i]);
	}
    PDM_printf ("part->edgeVtx \n");
    for (int i = 0; i < nEdge; i++) {
		PDM_printf("[%d]: %d %d\n", i, part->edgeVtx[2*i], part->edgeVtx[2*i+1]);
	}
}

  part->nEdge = nEdge;

  for (int i = 0; i < part->nVtx; i++) {
    part->vtxEdgeIdx[i+1] = part->vtxEdgeIdx[i+1] + part->vtxEdgeIdx[i];
  }
  part->vtxEdge  = (int *) malloc(sizeof(int) * part->vtxEdgeIdx[part->nVtx]);
  int *nVtxEdge  = (int *) malloc(sizeof(int) * part->nVtx);
  for (int i = 0; i < part->nVtx; i++) {
    nVtxEdge[i] = 0;
  }
  if (0 == 1) {
    PDM_printf ("part->vtxEdgeIdx   --- PDM_surf_part_build_edges \n");
    for (int i = 0; i < part->nVtx; i++) {
      PDM_printf ("[%d] :  %d \n", i+1, part->vtxEdgeIdx[i]);
     }
  }

  for (int i = 0; i < part->nEdge; i++) {
    int vtx1 = part->edgeVtx[2*i    ] - 1;
    int vtx2 = part->edgeVtx[2*i + 1] - 1;
    part->vtxEdge[part->vtxEdgeIdx[vtx1] + nVtxEdge[vtx1]++] = i+1;
    part->vtxEdge[part->vtxEdgeIdx[vtx2] + nVtxEdge[vtx2]++] = i+1;
  }
  free(nVtxEdge);

  if (1 == 1) {
    PDM_printf ("part->vtxEdge   --- PDM_surf_part_build_edges \n");
    for (int i = 0; i < part->nVtx; i++) {
      PDM_printf ("[%d] :", i+1);
      for (int j = part->vtxEdgeIdx[i]; j < part->vtxEdgeIdx[i+1]; j++) {
        PDM_printf (" %d", part->vtxEdge[j]);
      }
      PDM_printf ("\n");
    }
  }

  free(hashTable);
  free(hashTableIdx);
  free(listEdges);

  /*
   * Re-allocation to real size
   */

  part->edgeVtx = (int *) realloc(part->edgeVtx, sizeof(int) * 2 * nEdge);

  /*
   * edge -> face connectivity
   */

  part->edgeFace = (int *) malloc(sizeof(int) * 2 * nEdge);
  for (int i = 0; i <  2 * nEdge; i++)
    part->edgeFace[i] = 0;

  for (int i = 0; i < nEdgesFaces; i++) {
    int iEdge = edgesToCompressEdges[i];
    int ifac1 = 2*iEdge;
    int ifac2 = 2*iEdge + 1;

    if (part->edgeFace[ifac1] == 0)
      part->edgeFace[ifac1] = edgeFaceUncompress[i] + 1;
    else if (part->edgeFace[ifac2] == 0)
      part->edgeFace[ifac2] = edgeFaceUncompress[i] + 1;
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error _build_edges_part : Error in edgeFace computing\n");
      abort();
    }
  }

  if (0 == 1) {
    PDM_printf("edgeface : PDM_surf_part_build_edges\n");
    for (int i = 0; i < part->nEdge; i++) {
      PDM_printf("%d: %d %d\n", i, part->edgeFace[2*i], part->edgeFace[2*i+1]);
    }
  }

  free(edgesToCompressEdges);

  /*
   * face -> edge connectivity
   */

  part->faceEdgeIdx = (int *) malloc(sizeof(int) * (part->nFace + 1));
  memcpy(part->faceEdgeIdx, part->faceVtxIdx, sizeof(int) * (part->nFace + 1));
  const int *faceEdgeIdx = part->faceEdgeIdx;
  part->faceEdge = (int *) malloc(sizeof(int) * faceEdgeIdx[part->nFace]);
  int *n_faceEdge = (int *) malloc(sizeof(int) * part->nFace);
  for (int i = 0; i < part->nFace; i++)
    n_faceEdge[i] = 0;

  for (int i = 0; i < nEdge; i++) {
    int face1 = PDM_ABS (part->edgeFace[2*i]);
    int face2 = PDM_ABS (part->edgeFace[2*i+1]);
    if (face1 > 0)
      part->faceEdge[faceEdgeIdx[face1-1] + n_faceEdge[face1-1]++] = i + 1;
    if (face2 > 0)
      part->faceEdge[faceEdgeIdx[face2-1] + n_faceEdge[face2-1]++] = i + 1;
  }

    if (vb == 1) {
    PDM_printf ("part->faceEdge   --- PDM_surf_part_build_edges \n");
    for (int i = 0; i < part->nFace; i++) {
      PDM_printf ("[%d] :", i+1);
      for (int j = part->faceEdgeIdx[i]; j < part->faceEdgeIdx[i+1]; j++) {
        PDM_printf (" %d", part->faceEdge[j]);
      }
      PDM_printf ("\n");
    }
  }

  free(edgeFaceUncompress);
  free(n_faceEdge);
  vb = 1;
  if (vb == 1)   PDM_printf ("==== PDM_surf_part_build_edges ==== terminated ====\n");

}


/**
 * \brief Return faceLnToGn
 *
 *
 * \param [in]  part      Partition to compute
 *
 */

const PDM_g_num_t *
PDM_surf_part_faceLnToGn_get
(
PDM_surf_part_t *part
)
{

  return part->faceLnToGn;
}

/**
 * \brief Dump a PDM_surf_part_t object
 *
 * This function dumps a surf part structure
 *
 * \param [in]  part     surf part to dump
 *
 */

void
PDM_surf_part_dump
(
 PDM_surf_part_t *part
 )
{
	int vb = 1;
	if (vb == 1)   PDM_printf ("==== PDM_surf_part_dump ====\n");
  if (part != NULL) {
	  PDM_printf ("PDM_surf_part_t : \n");
    PDM_printf ("  - edgeLnToGn :\n", part->edgeLnToGn);
    if (part->edgeLnToGn != NULL) {
      for (int j=0; j<part->nEdge; j++)
        PDM_printf ("%d -> %d \n", j+1, part->edgeLnToGn[j]);
    }

    PDM_printf ("  - nFace : %d \n", part->nFace);
    PDM_printf ("  - faceNormal : %d : ", part->faceNormal);
    if (part->faceNormal == NULL)
      PDM_printf ("\n");
    else{
      for (int j=0; j<3*part->nFace; j++)
        PDM_printf ("%d ", part->faceNormal[j]);
      PDM_printf ("\n");}
    PDM_printf ("  - nGhostFace : %d \n", part->nGhostFace);
    PDM_printf ("  - nTotalFace : %d \n", part->nTotalFace);
    PDM_printf ("  - sFaceVtx : %d \n", part->sFaceVtx);
		PDM_printf ("  - faceVtx : \n");
		for (int j = 0; j <part->nFace; j++) {
		  PDM_printf ("%d-> ", part->faceLnToGn[j]);
		  for (int k = (part->faceVtxIdx)[j]; k < (part->faceVtxIdx)[j+1]; k++)
        PDM_printf (" "PDM_FMT_G_NUM, part->vtxLnToGn[part->faceVtx[k]-1]);
		  PDM_printf ("\n");
		}
		PDM_printf ("  - faceEdge : \n");
		for (int j = 0; j <part->nFace; j++) {
		  PDM_printf ("%ld-> ", part->faceLnToGn[j]);
		  for (int k = (part->faceEdgeIdx)[j]; k < (part->faceEdgeIdx)[j+1]; k++)
        PDM_printf (" "PDM_FMT_G_NUM, part->edgeLnToGn[part->faceEdge[k]-1]);
		  PDM_printf ("\n");
		}
    PDM_printf ("  - faceLnToGn :\n");
    for (int j=0; j<part->nFace; j++)
      PDM_printf ("%d -> %d\n", j+1, part->faceLnToGn[j]);
    PDM_printf ("  - nVtx : %d \n", part->nVtx);
    PDM_printf ("  - coords : \n");
    for (int j=0; j<part->nVtx; j++){
      PDM_printf ("%d-> ", part->vtxLnToGn[j]);
      PDM_printf ("%f ", part->coords[3*j]);
      PDM_printf ("%f ", part->coords[3*j+1]);
      PDM_printf ("%f ", part->coords[3*j+2]);
      PDM_printf ("\n");}
    PDM_printf ("  - vtxEdge loc with  ghost edges: \n");
    for (int j = 0; j <part->nVtx; j++) {
      PDM_printf ("%ld-> ", part->vtxLnToGn[j]);
      for (int k = (part->vtxEdgeIdx[j]); k < (part->vtxEdgeIdx)[j+1]; k++)
        PDM_printf (" %d",part->vtxEdge[k]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
    PDM_printf ("  - vtxEdge without ghost edges: \n");
    for (int j = 0; j <part->nVtx; j++) {
      PDM_printf ("%ld-> ", part->vtxLnToGn[j]);
      for (int k = (part->vtxEdgeIdx[j]); k < (part->vtxEdgeIdx)[j+1]; k++)
        if (part->vtxEdge[k] <= part->nEdge) {
          PDM_printf (" "PDM_FMT_G_NUM, part->edgeLnToGn[part->vtxEdge[k]-1]);
        }
		  PDM_printf ("\n");
		}
/* PDM_printf ("  - vtxLnToGn : \n", part->vtxLnToGn); */
/* for (int j=0; j<part->nVtx; j++) */
/* 	PDM_printf ("%d -> %d\n", j+1, part->vtxLnToGn[j]); */
    PDM_printf ("  - nEdge : %d \n", part->nEdge);
    PDM_printf ("  - nGhostEdge : %d \n", part->nGhostEdge);
    PDM_printf ("  - nTotalEdge : %d \n", part->nTotalEdge);
    PDM_printf ("  - edgeFace loc : ");
    for (int j=0; j<part->nEdge; j++){
      PDM_printf ("\n%ld : ", part->edgeLnToGn[j]);
      PDM_printf (" %d", part->edgeFace[2*j]);
      if (part->edgeFace[2*j+1] > 0)
        PDM_printf (" %d", part->edgeFace[2*j+1]);
    }
    PDM_printf ("\n");
    PDM_printf ("  - edgeFace with ghost face : ");
    for (int j=0; j<part->nEdge; j++){
      PDM_printf ("\n%ld : ", part->edgeLnToGn[j]);
      if (part->edgeFace[2*j] <= part->nFace) {
        PDM_printf (" "PDM_FMT_G_NUM, part->faceLnToGn[part->edgeFace[2*j]-1]);
      }
      if (part->edgeFace[2*j+1] <= part->nFace) {
        if (part->edgeFace[2*j+1] > 0)
          PDM_printf (" "PDM_FMT_G_NUM, part->faceLnToGn[part->edgeFace[2*j+1]-1]);
      }
    }
    PDM_printf ("\n");
    PDM_printf ("  - edgeVtx loc without ghost faces : ");
    for (int j=0; j<part->nEdge; j++){
      PDM_printf ("\n%ld : ",  part->edgeLnToGn[j]);
      PDM_printf (" %d", part->edgeVtx[2*j]);
      PDM_printf (" %d", part->edgeVtx[2*j+1]);}
    PDM_printf ("\n");
    PDM_printf ("  - edgeVtx : ");
    for (int j=0; j<part->nEdge; j++){
      PDM_printf ("\n%ld : ",  part->edgeLnToGn[j]);
      PDM_printf (" "PDM_FMT_G_NUM, part->vtxLnToGn[part->edgeVtx[2*j]-1]);
      PDM_printf (" "PDM_FMT_G_NUM, part->vtxLnToGn[part->edgeVtx[2*j+1]-1]);}
    PDM_printf ("\n");
    /* PDM_printf ("  - edgeLnToGn :\n", part->edgeLnToGn); */
    /* if (part->edgeLnToGn != NULL) { */
    /* 	for (int j=0; j<part->nEdge; j++) */
    /* 		PDM_printf ("%d -> %d \n", j+1, part->edgeLnToGn[j]); */
    /* } */
    PDM_printf ("  - edgePartBound : %d \n", part->edgePartBound);
    PDM_printf ("  - vtxPartBound : %d \n", part->vtxPartBound);
    PDM_printf ("  - carLgthVtx : %d : ", part->carLgthVtx);
    if (part->carLgthVtx == NULL)
      PDM_printf ("\n");
    else{
      for (int j=0; j<part->nVtx; j++)
        PDM_printf ("%d ", part->carLgthVtx[j]);
      PDM_printf ("\n");}
    PDM_printf ("  - extents : %d : ", part->extents);
    if (part->extents == NULL)
      PDM_printf ("\n");
    else{
      for (int j=0; j<part->nVtx; j++)
        PDM_printf ("%d ", part->extents[j]);
      PDM_printf ("\n");}


    fflush(stdout);
  }
if (vb == 1)   PDM_printf ("==== PDM_surf_part_dump ==== terminated ====\n");
}


/*
pour afficher un PDM_surf_part_t

PDM_printf ("nFace : %d \n", _part->nFace);
PDM_printf ("nGhostFace : %d \n", _part->nGhostFace);
PDM_printf ("nTotalFace : %d \n", _part->nTotalFace);
PDM_printf ("sFaceVtx : %d \n", _part->sFaceVtx);
PDM_printf ("faceVtxIdx : %d \n", _part->faceVtxIdx);
PDM_printf ("faceVtx : %d \n", _part->faceVtx);
PDM_printf ("faceEdgeIdx : %d \n", _part->faceEdgeIdx);
PDM_printf ("faceEdge : %d \n", _part->faceEdge);
PDM_printf ("faceLnToGn : %d \n", _part->faceLnToGn);
PDM_printf ("nVtx : %d \n", _part->nVtx);
PDM_printf ("coords : %d \n", _part->coords);
PDM_printf ("vtxEdgeIdx : %d \n", _part->vtxEdgeIdx);
PDM_printf ("vtxEdge : %d \n", _part->vtxEdge);
PDM_printf ("vtxLnToGn : %d \n", _part->vtxLnToGn);
PDM_printf ("nEdge : %d \n", _part->nEdge);
PDM_printf ("nGhostEdge : %d \n", _part->nGhostEdge);
PDM_printf ("nTotalEdge : %d \n", _part->nTotalEdge);
PDM_printf ("edgeFace : %d \n", _part->edgeFace);
PDM_printf ("edgeVtx : %d \n", _part->edgeVtx);
PDM_printf ("edgeLnToGn : %d \n", _part->edgeLnToGn);
PDM_printf ("edgePartBound : %d \n", _part->edgePartBound);
PDM_printf ("vtxPartBound : %d \n", _part->vtxPartBound);
PDM_printf ("carLgthVtx : %d \n", _part->carLgthVtx);
PDM_printf ("faceNormal : %d \n", _part->faceNormal);
PDM_printf ("extents : %d \n", _part->extents);


  */

#ifdef __cplusplus
}
#endif /* __cplusplus */
