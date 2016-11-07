
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"

#ifdef PDM_HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef PDM_HAVE_PTSCOTCH
#include <ptscotch.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_part.h"
#include "pdm_part_priv.h"
#include "pdm_timer.h"

#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_hilbert.h"
#include "pdm_geom_elem.h"
#include "pdm_sort.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/


/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Random order 
 *
 * \param [in]   nElt    Number of elements
 * \param [out]  order   Random order
 *
 */

static void 
_random_order
(
const int nElt,
      int *order
)
{
  int *tmpArray = (int *) malloc (sizeof(int) * nElt);
  
  srand(time(NULL));
  for (int i = 0; i < nElt; i++) {
    tmpArray[i] = rand()%nElt;
    order[i] = i;
  }
  
  PDM_sort_int (tmpArray, order, nElt);
}

/**
 * \brief Renumber face to cell connectivity 
 * 
 * \param [in]      nCell        Number of cells
 * \param [in]      nFace        Number of faces
 * \param [in]      cellFaceIdx  Cell face connectivity Index
 * \param [in]      cellFace     Cell face connectivity
 * \param [in, out] faceCell     Face cell connectivity
 *
 */

static void 
_renum_faceCell 
(
const int  nCell,
const int  nFace,
const int *cellFaceIdx, 
const int *cellFace, 
      int *faceCell 
)
{

  for (int i = 0; i < 2 * nFace; i++) {
    faceCell[i] = -1;
  }

  for (int i = 0; i < nCell; ++i) {
    for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
      int idx = 2 * (cellFace[j]-1);
      if (faceCell[idx] == -1) { 
        faceCell[idx] = i + 1;
      }
      else { 
        faceCell[idx + 1] = i + 1;
      }
    }
  }
 
}


/**
 * \brief Order an array
 * 
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder        New order (size = \ref nElt
 * \param [in, out] Array         	Array to renumber
 *
 */

static void 
_order_array 
(
const int    sizeArray,
const size_t elt_size,        
const int   *newToOldOrder,
void        *array
)
{
  unsigned char *oldArray = (unsigned char *) malloc (sizeArray * elt_size);
  unsigned char *_array = (unsigned char *) array;
  
  for (int i = 0; i < sizeArray; ++i) {
    for (int j = 0; j < elt_size; ++j) {
      oldArray[elt_size * i + j] = _array[elt_size * i + j];
    }
  }
  
  for (int i = 0; i < sizeArray; ++i) {
    for (int j = 0; j < elt_size; ++j) {
      _array[elt_size * i + j] = oldArray[elt_size * newToOldOrder[i] +j];
    }
  }
  
  free(oldArray);
}


/**
 * \brief Order faceCell array
 * 
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder   New order (size = \ref nElt
 * \param [in, out] faceCell        Array to order
 *
 */

static void 
_order_faceCell 
(
int           nFace,
int          *newToOldOrder,
PDM_g_num_t *faceCell 
)
{
  int *oldFaceCell = (int*) malloc (nFace * 2 * sizeof(int));
  for(int i = 0; i < nFace * 2; ++i) {
    oldFaceCell [i] = faceCell [i];
  }
  
  for(int i = 0; i < nFace; ++i) {
    faceCell[i*2+0] = oldFaceCell[newToOldOrder[i]*2+0];
    faceCell[i*2+1] = oldFaceCell[newToOldOrder[i]*2+1];
  }
  
  free(oldFaceCell);
}


/**
 * \brief Order an array
 * 
 * \param [in]      sizeArray       Number of elements
 * \param [in]      newToOldOrder        New order (size = \ref nElt
 * \param [in, out] Array         	Array to renumber
 *
 */

static void 
_renum_array 
(
const int    sizeArray,
const int   *olToNewOrder,
int         *array
)
{
  int *oldArray = (int *) malloc (sizeof(int) * sizeArray);
  
  for (int i = 0; i < sizeArray; ++i) {
    oldArray[i] = array[i];
  }
  
  for (int i = 0; i < sizeArray; ++i) {
    array[i] = olToNewOrder[oldArray[i]-1] + 1;
  }
  
  free(oldArray);
}


/**
 * \brief Renumber connectivities 
 * 
 * \param [in]      nElt            Number of elements
 * \param [in]      newToOldOrder        New order (size = \ref nElt
 * \param [in, out] connectivityIdx	Connectivity index
 * \param [in, out] connectivities	Element connectivities
 *
 */
 
static void 
_renum_connectivities 
(
const int nElt,
const int *newToOldOrder,
int       *connectivityIdx,
int       *connectivities
)
{
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);
  
  int *oldConnectivities = (int *) malloc (connectivityIdx[nElt] * sizeof(int));
  
  for (int i = 0; i < connectivityIdx[nElt]; ++i) { 
    oldConnectivities[i] = connectivities[i];
  }
    
  int *oldConnectivityIdx = (int *) malloc ((nElt + 1) * sizeof(int));
  
  for (int i = 0; i < nElt+1; ++i) { 
    oldConnectivityIdx[i] = connectivityIdx[i];
  }

  for (int elem = 0; elem < nElt; ++elem) {
    int nbSSElem = oldConnectivityIdx[newToOldOrder[elem] + 1] - oldConnectivityIdx[newToOldOrder[elem]];
    connectivityIdx[elem+1] = nbSSElem;
  }
  
  connectivityIdx[0] = 0;
  for (int elem = 1; elem < nElt+1; ++elem) {
    connectivityIdx[elem] += connectivityIdx[elem-1];
  }
    
  for (int elem = 0; elem < nElt; ++elem) {
    int nbSSElem = oldConnectivityIdx[newToOldOrder[elem] + 1] - oldConnectivityIdx[newToOldOrder[elem]];
    
    for (int ssElem = 0; ssElem < nbSSElem; ++ssElem) {
      connectivities[connectivityIdx[elem] + ssElem] = oldConnectivities[oldConnectivityIdx[newToOldOrder[elem]]+ssElem];
    }
  }
  
  free(oldConnectivities);
  free(oldConnectivityIdx);
  
}


/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in]  part      Current PPART structure
 * \param [out] cellCenter Cell center
 *
 */

static void 
_compute_cellCenter 
(
_part_t *part,
double  *cellCenter
)
{
  if (part->nCell <= 0) {
    return;
  }
  
  int isPoly3D = (part->faceVtxIdx[1] > 2);

  if (isPoly3D) {
    const int isOriented = 0;
    double *volume = (double *) malloc (part->nCell * sizeof(double));
    int isDegenerated;  

    PDM_geom_elem_polyhedra_properties (isOriented,
                                        part->nCell,
                                        part->nFace,
                                        part->faceVtxIdx,
                                        part->faceVtx,
                                        part->cellFaceIdx,
                                        part->cellFace,
                                        part->nVtx,
                                        part->vtx,
                                        volume,
                                        cellCenter,
                                        NULL,
                                        &isDegenerated);
    free (volume);
  }
  
  else {
    double *surfaceVector = (double * ) malloc( sizeof(double) * 3 * part->nCell); 
    int isDegenerated;  
    
    int *connectivity = (int *) malloc (part->cellFaceIdx[part->nCell] 
                        * sizeof(int));

    
    int idx = 0;    
    for (int icell = 0; icell < part->nCell; icell++) {
      int faceCurr = part->cellFace[part->cellFaceIdx[icell]] - 1;
      int deb  = part->faceVtx[part->faceVtxIdx[faceCurr]];
      int next = part->faceVtx[part->faceVtxIdx[faceCurr]+1];
      connectivity[idx++] = deb;
      connectivity[idx++] = next;
      while (next != deb) {
        for (int j = part->cellFaceIdx[icell]; j < part->cellFaceIdx[icell+1]; ++j) {
          int face = part->cellFace[j] - 1;
          if (faceCurr != face) {
            int s1 = part->faceVtx[part->faceVtxIdx[face]  ];
            int s2 = part->faceVtx[part->faceVtxIdx[face]+1];
						
            if ((s1 == next) || (s2 == next)) { 
              if (s1 == next) {
                next = s2;
              }
              else if (s2 == next) {
                next = s1;
              }
              if (next != deb) {
                connectivity[idx++] = next;
              }
              faceCurr = face;
              break;
            }              
          }
        }
      }
    }

		assert (idx == part->cellFaceIdx[part->nCell]);
		
    PDM_geom_elem_polygon_properties (part->nCell,   
                                      part->cellFaceIdx, 
                                      connectivity,
                                      part->vtx,
                                      surfaceVector,
                                      cellCenter,
                                      NULL,
                                      &isDegenerated);

    free (surfaceVector);
    free (connectivity);

  }

}


/**
 *
 * \brief Perform cells renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */

static void 
_renum_cells
(
 _part_t *part, 
 int     *newToOldOrder
)
{
  /*
   * Cell Renumbering
   */
   
  _renum_connectivities (part->nCell,
                         newToOldOrder,
                         part->cellFaceIdx, 
                         part->cellFace); 
  
  if (part->cellTag != NULL) {
    _order_array (part->nCell,
                  sizeof(int),
                  newToOldOrder,
                  part->cellTag);
  }
   
  _order_array (part->nCell,
                sizeof(PDM_g_num_t),
                newToOldOrder,
                part->cellLNToGN); 
   
  _renum_faceCell (part->nCell,
                   part->nFace,
                   part->cellFaceIdx, 
                   part->cellFace, 
                   part->faceCell); 
   
}


/**
 *
 * \brief Perform faces renumbering from a new order
 *
 * \param [in,out]  part        Current partition
 * \param [in]      newToOldOrder    NewOrder
 *
 */

static void 
_renum_faces
(
_part_t *part, 
int     *newToOldOrder
)
{
  
  /** Renum FaceVtx / FaceVtxIdx **/

  _renum_connectivities (part->nFace, 
                         newToOldOrder, 
                         part->faceVtxIdx, 
                         part->faceVtx);

  /** CellFace **/
  
   int *oldToNewOrder = (int *) malloc (part->nFace * sizeof(int));
  
   for(int i = 0; i < part->nFace; i++) {
    oldToNewOrder[newToOldOrder[i]] = i;
   }
 
  _renum_array (part->cellFaceIdx[part->nCell], 
                oldToNewOrder,
                part->cellFace);
  
  free (oldToNewOrder);
    
    /** FaceTag **/

  if (part->faceTag != NULL) {
    _order_array (part->nFace,
                  sizeof(int),
                  newToOldOrder,
                  part->faceTag); 
  }
    
   /** FaceLNToGN **/

  _order_array (part->nFace,
                sizeof(PDM_g_num_t),
                newToOldOrder,
                part->faceLNToGN); // OK
    
    /** FaceCell Face **/
  
  _order_faceCell (part->nFace, 
                   newToOldOrder,
                   part->faceCell);  

}


/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_cells_hilbert 
(
_PDM_part_t* ppart
)
{
  
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    double *cellCenter = 
        (double *) malloc (part->nCell * 3 * sizeof(double ));
    PDM_hilbert_code_t *hilbertCodes = 
        (PDM_hilbert_code_t *) malloc (part->nCell * sizeof(PDM_hilbert_code_t));
      
    /** Barycentre computation **/
    
    _compute_cellCenter (part, cellCenter);
        
    double extents[3 * 2]; 
    
    /** Get EXTENTS LOCAL **/
    
    PDM_hilbert_get_coord_extents_seq(3, part->nCell, cellCenter, extents);
    
    /** Hilbert Coordinates Computation **/
    
    PDM_hilbert_encode_coords(3, PDM_HILBERT_CS, extents, part->nCell, cellCenter, hilbertCodes);
    
    /** CHECK H_CODES **/
    
    free(cellCenter);
    
    int *newToOldOrder = (int *) malloc (part->nCell * sizeof(int));
    for(int i = 0; i < part->nCell; ++i) {
      newToOldOrder [i] = i;
    }
      
    PDM_sort_double (hilbertCodes, newToOldOrder, part->nCell);
	  
    _renum_cells (part, newToOldOrder);
          
    free (hilbertCodes);
    free (newToOldOrder);
    
  }
}


/**
 *
 * \brief Perform a cells renumbering from a Hilbert curve
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_cells_random 
(
_PDM_part_t* ppart
)
{
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    const int nCell = part->nCell;
    
    int *order = (int *) malloc (sizeof(int) * nCell);
    
    _random_order (nCell, order);
    
    _renum_cells (part, order);
      
    free (order);
  }
}


/**
 *
 * \brief Perform a face random renumbering 
 *
 * \param [in,out]  ppart    Current PPART structure
 *
 */

static void 
_renum_faces_random 
(
_PDM_part_t* ppart
)
{
  for(int ipart = 0; ipart < ppart->nPart; ++ipart) {
    _part_t *part = ppart->meshParts[ipart];       
    const int nFace = part->nFace;
    
    int *order = (int *) malloc (sizeof(int) * nFace);
    
    _random_order (nFace, order);
    
    _renum_faces (part, order);
      
    free (order);
  }
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  ppart       ppart structure
 * \param [in]      entity      Mesh entity to renumber
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_cell
(
 _PDM_part_t              *ppart,
 PDM_part_renum_cell_t  method                 
)
{
  switch (method) {
  case PDM_PART_RENUM_CELL_HILBERT :
    _renum_cells_hilbert(ppart); 
    break;
  case PDM_PART_RENUM_CELL_RANDOM :
    _renum_cells_random(ppart); 
    break;
  case PDM_PART_RENUM_CELL_NONE :
    break;
  default:
    fprintf (stderr, "PDM_part_renum Error : unavailable face renumbering method\n");
    exit (1);    
  }
}



/**
 *
 * \brief Perform mesh entities renumbering
 *
 * \param [in,out]  ppart       ppart structure
 * \param [in]      entity      Mesh entity to renumber
 * \param [in]      method      Renumbering method
 *
 */

void 
PDM_part_renum_face
(
 _PDM_part_t           *ppart,
 PDM_part_renum_face_t  method                 
)
{
  switch (method) {
    case PDM_PART_RENUM_FACE_RANDOM :
      _renum_faces_random (ppart); 
      break;
    case PDM_PART_RENUM_FACE_NONE :
      break;
    default:
      fprintf (stderr, "PDM_part_renum Error : unavailable face renumbering method\n");
      exit (1);    
  }
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

