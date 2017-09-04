
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

enum {
  NOT_DEFINE,
  IN_STACK,        
  UNCHANGED_CYCLE,
  CHANGED_CYCLE,
};

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
  
 /* 
  * Orient the first cell. As the oriented volume is positive, the face noramals
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
  * Other cells are oriented from the faces of the first cell
  */
  
  int *internFaceCell = malloc (sizeof(int)* 2 * nFace);
    
  for (int i = 0; i < 2 * nFace; i++) {
    internFaceCell[i] = 0;
  }

  for (int i = cellFaceIdx[0]; i < cellFaceIdx[1]; i++) {
    int Face = cellFace[i];
    
    if (Face > 0) {
      internFaceCell[2 * (Face - 1) + 1] = 1;
    }
    else {
      internFaceCell[2 * (PDM_ABS (Face) - 1)] = 1;
    }
  }
  
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

  int *stack = (int *) malloc (sizeof(int) * maxNPolyFace);
  int *tagFace = (int *) malloc (sizeof(int) * maxNPolyFace);

  for (int ipoly = 0; ipoly < nPolyhedra; ipoly++) {

    const int polyIdx   = cellFaceIdx[ipoly];
    const int nPolyFace = cellFaceIdx[ipoly + 1] - polyIdx;

    /*
     * Intialize cell center to the barycenter
     */

    nKeyPoly = 0;
    
    /* Search polyhedra vertices */

    int nPolyhedraVertices = 0;
    
    for (int iface = 0; iface < nPolyFace; iface++) {
      tagFace[iface] = NOT_DEFINE; 
    }

  }
  /* Continuer ici*/
  
  
  PDM_hash_tab_free(hashOrient);
  
  free (stack);
  free (tagFace);
  free (internFaceCell);
  
  

}

#ifdef	__cplusplus
}
#endif

