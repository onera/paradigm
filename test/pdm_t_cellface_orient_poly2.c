#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_cellface_orient.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  /*
   *  Init
   */

  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

//nCell :  5
//nFace :  25
//nVtx  :  22
//VtxCoord  :  [ 0.  0.  0.  1.  0.  0.  0.  1.  0.  1.  1.  0.  0. 2.  0.  1.  2.  0.
//  0.  0.  1.  1.  0.  1.  0.  1.  1.  1.  1.  1.  0.  2.  1.  1. 2.  1.
//  0.  0.  2.  1.  0.  2.  0.  1.  2.  1.  1.  2.  0.  2.  2.  1. 2.  2.
//  3.  2.  2.  3.  0.  0.  3.  0.  2.  3.  2.  0.]
//cellFaceIdx  :  [ 0  6 15 21 27 33]
//cellFace  :  [ 2  9  3 18  1 13  9 19 14 23  4  6  7  5  8 11 14 10 24 13 12 15 19 18 16
// 17 25 21 23 24 20 25 22]
//faceVtxIdx  :  [  0   4   8  12  16  20  24  28  32  36  40  44 48  52  56  60  64  68
//  72  76  80  84  88  92  96 100]
//faceVtx  :  [ 3  4  2  1  7  9  3  1  1  2  8  7 20 22 19 21  6 22 20  2  2 20 21 14 18
// 19 22  6 14 21 19 18  2  4 10  8  7  8 14 13 13 15  9  7 13 14 16 15  7  8
// 10  9  8 10 16 14  9 11  5  3 11 12  6  5  5  6  4  3  3  4 10 9  4  6 12
// 10 17 18 12 11 15 17 11  9 15 16 18 17 10 12 18 16  9 10 16 15  9 10 12 11]

  const int      nCell = 5;
  const int      nFace = 25;
  const int      nVtx  = 22; 
  double   vtx[66] = { 0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  1.,  
  1.,  0.,  0., 2.,  0.,  1.,  2.,  0.,
  0.,  0.,  1.,  1.,  0.,  1.,  0.,  1.,  1.,  1.,  1.,  1.,  0.,  2.,  1.,  1., 2.,  1.,
  0.,  0.,  2.,  1.,  0.,  2.,  0.,  1.,  2.,  1.,  1.,  2.,  0.,  2.,  2.,  1., 2.,  2.,
  3.,  2.,  2.,  3.,  0.,  0.,  3.,  0.,  2.,  3.,  2.,  0.};        

  int  cellFaceIdx[6] = {0,  6, 15, 21, 27, 33};
  int  cellFace[33] = { 2,  9,  3, 18,  1, 13,  9, 19, 14, 23,  4,  6,  7,  5,
  8, 11, 14, 10, 24, 13, 12, 15, 19, 18, 16,
 17, 25, 21, 23, 24, 20, 25, 22};
  int  *faceCell = NULL;
  
  int faceVtxIdx[26] = {0,   4,   8,  12,  16,  20,  24,  28,  32,  36,  40,  44,
  48,  52,  56,  60,  64,  68,
  72,  76,  80,  84,  88,  92,  96, 100};

  int faceVtx[100] = {3,  4,  2,  1,  7,  9,  3,  1,  1,  2,  8,  7, 20, 22, 19, 21,  6, 22, 20,  2,  2, 20, 21, 14, 18,
 19, 22,  6, 14, 21, 19, 18,  2,  4, 10,  8,  7,  8, 14, 13, 13, 15,  9,  7, 13, 14, 16, 15,  7,  8,
 10,  9,  8, 10, 16, 14,  9, 11,  5,  3, 11, 12,  6,  5,  5,  6,  4,  3,  3,  4, 10, 9,  4,  6, 12,
 10, 17, 18, 12, 11, 15, 17, 11,  9, 15, 16, 18, 17, 10, 12, 18, 16,  9, 10, 16, 15,  9, 10, 12, 11};

  PDM_cellface_orient (nCell,
                       nFace,
                       nVtx,
                       vtx,        
                       cellFaceIdx,
                       cellFace,
                       faceCell,
                       faceVtxIdx,
                       faceVtx);

  PDM_MPI_Finalize();

  return 0;
}

