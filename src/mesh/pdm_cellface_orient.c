
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
}

#ifdef	__cplusplus
}
#endif

