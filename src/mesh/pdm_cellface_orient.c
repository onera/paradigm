
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
#include "pdm_geom_elem.h"
#include "pdm_cellface_orient.h"
#include "pdm_hash_tab.h"
#include "pdm_error.h"

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


enum {false, true};

typedef enum {
  FACE_UNPROCESSED,
  FACE_IN_STACK,        
  FACE_UNCHANGED_CYCLE,
  FACE_CHANGED_CYCLE
} _face_state_t;


typedef enum {
  CELL_UNPROCESSED,
  CELL_IN_STACK,        
  CELL_COMPLETED
} _cell_state_t;


/**
 * \struct ????
 * \brief  ????
 * 
 *
 */

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/



/**
 * \brief Orient cell->face connectivity
 * 
 * At the output of th function, a face number in \ref cellFace is positive 
 * if surface normal is inside the cell, negative otherwise. \ref faceCell is 
 * oriented in the same way  
 *
 * \param [in]     nCell       Number of cells
 * \param [in]     nFace       Number of faces
 * \param [in]     nVtx        Number of vertices
 * \param [in]     coords      Vertices coordinates
 * \param [in]     cellFaceIdx Cell to face connectivity index (size = \ref nCell + 1)
 * \param [in, out]cellFace    Cell to face connectivity (size = cellFaceIdx[nCell])
 * \param [in, out]faceCell    face to cell connectivity (size = 2 * \ref nFace) or NULL
 * \param [in]     faceVtxIdx  face to vertex connectivity index (size = \ref nFace + 1)
 * \param [in]     faceVtx     face to vertex connectivity (size = faceVtxIdx[nFace])
 *
 */

void
PDM_cellface_orient
(
const int      nCell,
const int      nFace,
const int      nVtx,
const double  *coords,        
const int     *cellFaceIdx,
int           *cellFace,
int           *faceCell,
const int     *faceVtxIdx,
const int     *faceVtx
)
{

  if (nCell == 0) {
    return;
  }

  int *_faceCell = NULL;
  
  if (faceCell != NULL) {
    _faceCell = faceCell; 
  }
  else {
    _faceCell = malloc (sizeof(int)* 2 * nFace);
  }
  
  for (int i = 0; i < 2 * nFace; i++) {
    _faceCell[i] = 0;
  }
  
  for (int i = 0; i < nCell; i++) {
    for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
      int face = PDM_ABS (cellFace[j]);
      int iFace = 2 * (face - 1);
      if (faceCell[iFace] == 0) {
        faceCell[iFace] = i + 1;
      }
      else {
        faceCell[iFace+1] = i + 1;
      }
    }
  }
  
 /* 
  * Orient the first cell. 
  * ----------------------
  * 
  * As the oriented volume is positive, the face normals
  * are outside of the element 
  * 
  */
 
  int     isOriented = 0;
  int     nPolyhedra = 1;
  double  volume[3];
  double  center[3];

  PDM_geom_elem_polyhedra_properties (isOriented,
                                      nPolyhedra,
                                      nFace,
                                      faceVtxIdx,
                                      faceVtx,   
                                      cellFaceIdx,
                                      cellFace,
                                      nVtx,
                                      coords,
                                      volume,
                                      center,
                                      NULL,
                                      NULL);  
 /* 
  * Initialize an oriented face cell
  * --------------------------------
  * 
  * left cell : normal is inside the cell
  * right cell : normal is outside the cell
  * 
  * The orientation of the first cell is taking into account
  * 
  */
  
  int *orientedFaceCell = malloc (sizeof(int)* 2 * nFace);
    
  for (int i = 0; i < 2 * nFace; i++) {
    orientedFaceCell[i] = 0;
  }

  for (int i = cellFaceIdx[0]; i < cellFaceIdx[1]; i++) {
    int Face = cellFace[i];
    
    if (Face > 0) {
      orientedFaceCell[2 * (Face - 1)] = 1;
    }
    else {
      orientedFaceCell[2 * (PDM_ABS (Face) - 1) + 1] = 1;
    }
  }
  
 /* 
  * Other cells are oriented from the faces of the first cell
  * ---------------------------------------------------------
  * 
  */

  /* Build and initialize local structures : stack and tag arrays */
  
  int keyMax = 2 * nVtx;
  
  int nKeyPoly = 0;
  int sKeyPoly = 10;
 
  int maxNPolyFace = -1;

  PDM_hash_tab_t *hashOrient = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, &keyMax);

  int *keyPoly = (int *) malloc (sizeof(int) * sKeyPoly);
 
  for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {
    const int polyIdx   = cellFaceIdx[ipoly];
    const int nPolyFace = cellFaceIdx[ipoly + 1] - polyIdx;
    maxNPolyFace = PDM_MAX (maxNPolyFace, nPolyFace);
  }

  int *stackFace = (int *) malloc (sizeof(int) * maxNPolyFace);

  int nStackCell = -1;
  int *stackCell = (int *) malloc (sizeof(int) * nCell);
  
  int *tagFace = (int *) malloc (sizeof(int) * maxNPolyFace);
  int *tagCell = (int *) malloc (sizeof(int) * nCell);
  
  for (int i = 0; i < nCell; i++) {
    tagCell[i] = CELL_UNPROCESSED;
  }

  for (int iface = 0; iface < maxNPolyFace; iface++) {
    tagFace[iface] = FACE_UNPROCESSED; 
  }

  tagCell[0] = CELL_COMPLETED;

  /* Add neighbours of the first cell in the stack */
  
  for (int i = cellFaceIdx[0]; i < cellFaceIdx[1]; i++) {
    int iFace = 2 * (PDM_ABS (cellFace[i]) - 1);
    if (_faceCell[iFace] == 1) {
      if (_faceCell[iFace + 1] > 0) {
        int cell = PDM_ABS (_faceCell[iFace + 1]);
        if (tagCell[cell - 1] == CELL_UNPROCESSED) {
          stackCell[++nStackCell] = cell;
          tagCell[cell - 1] = CELL_IN_STACK;
        }
      }
    }
    else {
      int cell = PDM_ABS (_faceCell[iFace]);
      if (tagCell[cell - 1] == CELL_UNPROCESSED) {
        stackCell[++nStackCell] = cell;
      }
      tagCell[cell - 1] = CELL_IN_STACK;
    }
  }

  /* Orientation process */

  while (nStackCell >= 0) {
    int iCell = stackCell[nStackCell--] - 1;
    
    printf("--- icell : %d\n", iCell);
    
    if (tagCell[iCell] == CELL_COMPLETED) {
      continue;
    }
    
    const int polyIdx   = cellFaceIdx[iCell];
    const int nPolyFace = cellFaceIdx[iCell + 1] - polyIdx;

    /* Build pseudo edges of the current cell and store them into a hash table */
    
    int fistProcessedFace = -1;
    for (int iface = 0; iface < nPolyFace; iface++) {

      tagFace[iface] = FACE_UNPROCESSED;

      const int face          = PDM_ABS (cellFace[polyIdx + iface]) - 1;
      const int faceIdx       = faceVtxIdx[face];
      const int nFaceVertices = faceVtxIdx[face+1] - faceIdx;

      printf ("orientedFaceCell %d : %d %d\n", face, orientedFaceCell[2*face],
              orientedFaceCell[2*face+1]);
      
      if (orientedFaceCell[2*face] != 0) {
        assert (orientedFaceCell[2*face] != iCell + 1);
        assert (orientedFaceCell[2*face+1] != iCell + 1);
        assert (orientedFaceCell[2*face+1] == 0);
        tagFace[iface] = FACE_CHANGED_CYCLE;
        orientedFaceCell[2*face+1] = iCell + 1;
        fistProcessedFace = iface;        
        printf("pass 1\n");
      }
      else if (orientedFaceCell[2*face + 1] != 0) {
        assert (orientedFaceCell[2*face] != iCell + 1);
        assert (orientedFaceCell[2*face+1] != iCell + 1);
        assert (orientedFaceCell[2*face] == 0);        
        tagFace[iface] = FACE_UNCHANGED_CYCLE;        
        orientedFaceCell[2*face] = iCell + 1;        
        fistProcessedFace = iface;        
        printf("pass 2\n");
      }        

      for (int ivert = 0; ivert < nFaceVertices; ivert++) {
        const int vertex = faceVtx[faceIdx + ivert] - 1;

        const int inext = (ivert + 1) % nFaceVertices;
        const int vertexNext = faceVtx[faceIdx + inext] - 1;
        const int key = vertex + vertexNext;

        if (nKeyPoly >= sKeyPoly) {
          sKeyPoly *= 2; 
          keyPoly = (int *) realloc (keyPoly, sizeof(int) * sKeyPoly);
        }

        keyPoly[nKeyPoly++] = key;

        int *edge = malloc (sizeof(int)*3);
        edge[0] = vertex;
        edge[1] = vertexNext;
        edge[2] = iface;

        PDM_hash_tab_data_add (hashOrient, (void *) &key, edge);

      }
    }
    
    assert (fistProcessedFace != -1); // On vérifie qu'on a deja une face renseignée
    
    // TODO
    int nStackFace = -1;    
    
    /* Look a neighbour of this face */
    if (fistProcessedFace != -1) {
      const int face          = PDM_ABS (cellFace[polyIdx + fistProcessedFace]) - 1;
      const int faceIdx       = faceVtxIdx[face];
      const int nFaceVertices = faceVtxIdx[face+1] - faceIdx;

      for (int ivert = 0; ivert < nFaceVertices; ivert++) {
        const int inext = (ivert + 1) % nFaceVertices;

        const int vertex = faceVtx[faceIdx + ivert] - 1;
        const int vertexNext = faceVtx[faceIdx + inext] - 1;
        int key = vertex + vertexNext;

        int nData = PDM_hash_tab_n_data_get (hashOrient, &key);
        int **data = (int **) PDM_hash_tab_data_get (hashOrient, &key);

        for (int j = 0; j < nData; j++) {
          if (data[j] != NULL) {
            int *_edge = data[j];
            int isInverseEdge = (vertex == _edge[1]) && (vertexNext == _edge[0]);
            int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
            int isSameFace    = fistProcessedFace == _edge[2];

            if (!isSameFace) {
              if (isSameEdge || isInverseEdge) {
                stackFace[++nStackFace] = ivert;
                tagFace[nStackFace] = FACE_IN_STACK;
                break;
              }
            }
          }
        }
        if (nStackFace != -1) {
          break;
        }
      }
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Internal Error : No previous processed face \n");
    }
    
    nKeyPoly = 0;
    
    while (nStackFace >= 0) {

      int iFace = stackFace[nStackFace--];
      
      if ((tagFace[iFace] == FACE_UNCHANGED_CYCLE) || 
          (tagFace[iFace] == FACE_CHANGED_CYCLE)) {
        continue;
      }
      
      const int face          = PDM_ABS (cellFace[polyIdx + iFace]) - 1;
      const int faceIdx       = faceVtxIdx[face];
      const int nFaceVertices = faceVtxIdx[face+1] - faceIdx;

      for (int ivert = 0; ivert < nFaceVertices; ivert++) {
        const int inext = (ivert + 1) % nFaceVertices;

        const int vertex = faceVtx[faceIdx + ivert] - 1;
        const int vertexNext = faceVtx[faceIdx + inext] - 1;
        int key = vertex + vertexNext;

        int nData = PDM_hash_tab_n_data_get (hashOrient, &key);
        int **data = (int **) PDM_hash_tab_data_get (hashOrient, &key);

        int jCurrentEdge = -1;
        for (int j = 0; j < nData; j++) {
          if (data[j] != NULL) {
            int *_edge = data[j];
            int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
            int isSameFace    = iFace == _edge[2];
            if (isSameEdge && isSameFace) {
              jCurrentEdge = j;
              break;
            }
          }
        }

        assert (jCurrentEdge > -1);

        for (int j = 0; j < nData; j++) {
          if (data[j] != NULL) {
            int *_edge = data[j];
            int isInverseEdge = (vertex == _edge[1]) && (vertexNext == _edge[0]);
            int isSameEdge    = (vertex == _edge[0]) && (vertexNext == _edge[1]);
            int isSameFace    = iFace == _edge[2];
            
            int neighbour = _edge[2];
            
            if (!isSameFace) {
              if (isSameEdge || isInverseEdge) {

                if (tagFace[iFace] < FACE_UNCHANGED_CYCLE) { 

                  if (tagFace[neighbour] >= FACE_UNCHANGED_CYCLE) {
                    if (tagFace[neighbour] == FACE_UNCHANGED_CYCLE) {
                      if (isSameEdge) {
                        tagFace[iFace] = FACE_CHANGED_CYCLE;
                      }
                      else  {
                        tagFace[iFace] = FACE_UNCHANGED_CYCLE;
                      }
                    }
                    else {
                      if (isSameEdge) {
                        tagFace[iFace] = FACE_UNCHANGED_CYCLE;
                      }
                      else  {
                        tagFace[iFace] = FACE_CHANGED_CYCLE;
                      }                      
                    }                    

                    free (data[j]);
                    free (data[jCurrentEdge]);
                    data[j] = NULL;
                    data[jCurrentEdge] = NULL;

                  }
                }

                if (tagFace[neighbour] == FACE_UNPROCESSED) {
                  stackFace[++nStackFace] = neighbour;
                  tagFace[neighbour] = FACE_IN_STACK;
                }

                break;

              }
            }
          }
        }
      }
      
      if (tagFace[iFace] == FACE_IN_STACK) {
        printf ("Error reorient : no neighbour processed face found");
        abort();        
      }
    }

    /* Add neighbour cell in the stack */
    
    for (int iface = 0; iface < nPolyFace; iface++) {

      if (tagFace[iface] == FACE_CHANGED_CYCLE) {
        cellFace[polyIdx + iface] = -cellFace[polyIdx + iface];
      }
      tagFace[iface] = FACE_UNPROCESSED; 

      const int face          = PDM_ABS (cellFace[polyIdx + iface]) - 1;

      int nextCell = -1;
      if (faceCell[2 * face] == (iCell + 1)) {
        nextCell = faceCell[2 * face + 1] - 1;
      }
      else {
        nextCell = faceCell[2 * face] - 1;
      }

      if (tagCell[nextCell] == CELL_UNPROCESSED) {
        stackCell[++nStackCell] = nextCell + 1;
        tagCell[nextCell] = CELL_IN_STACK;
      }
    }
    
    tagCell[iCell] = CELL_COMPLETED;

  }
  
  /* Orient FaceCell */
  
  if (faceCell != NULL) {
    for (int i = 0; i < nCell; i++) {
      for (int j = cellFaceIdx[i]; j < cellFaceIdx[i+1]; j++) {
        int face = cellFace[j];
        int iFace = 2 * (PDM_ABS(face) - 1);
        if (PDM_ABS (faceCell[iFace]) == i+1) {
          if (face < 0) {
            faceCell[iFace] = -(i+1);
          }
          else {
            faceCell[iFace] = i+1;           
          }
        }
        else {
          if (face < 0) {
            faceCell[iFace+1] = -(i+1);
          }
          else {
            faceCell[iFace+1] = i+1;           
          }
        }
      }
    }
  }

  PDM_hash_tab_free(hashOrient);
  
  free (stackFace);
  free (tagFace);
  free (orientedFaceCell);
  if (faceCell == NULL) {
    free (_faceCell);
  }
  

}

#ifdef	__cplusplus
}
#endif

