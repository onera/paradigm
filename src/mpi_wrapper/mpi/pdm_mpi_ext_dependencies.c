/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef PDM_HAVE_PARMETIS
#include <parmetis.h>
#include <metis.h>
#endif
#ifdef PDM_HAVE_PTSCOTCH
#include <ptscotch.h>
#include <scotch.h>
#endif

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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


#ifdef PDM_HAVE_PARMETIS

int 
PDM_ParMETIS_V3_PartKway 
(
const PDM_g_num_t *vtxdist, 
const PDM_g_num_t *xadj, 
const PDM_g_num_t *adjncy, 
const int *vwgt, 
const int *adjwgt, 
const int *wgtflag, 
const int *numflag, 
const int *ncon, 
const int *nparts, 
const double *tpwgts, 
const double *ubvec, 
const int *edgecut, 
int *part, 
const PDM_MPI_Comm comm
)
{
  MPI_Comm mpi_comm = *((MPI_Comm *) PDM_MPI_2_mpi_comm (comm));
  int iRank = 0;
  PDM_MPI_Comm_rank (comm, &iRank);
  
  int iSize = 0;
  PDM_MPI_Comm_size (comm, &iSize);
  
  PDM_g_num_t nNode = vtxdist[iRank+1] - vtxdist[iRank];
  PDM_g_num_t nEdge = xadj[nNode];

  idx_t _wgtflag = *wgtflag;
  idx_t _numflag = *numflag;
  idx_t _ncon = *ncon;; 
  idx_t _nparts = *nparts; 

  idx_t options[METIS_NOPTIONS]; /* Options */
  METIS_SetDefaultOptions(options);

  options[METIS_OPTION_NUMBERING] = 0; //C numbering = 0 (Fortran = 1)
  options[METIS_OPTION_MINCONN] = 1; //Minimize the maximum connectivity
  options[METIS_OPTION_CONTIG] = 1; //Force contiguous partitions
  //The graph should be compressed by combining together vertices that have identical adjacency lists.
  options[METIS_OPTION_COMPRESS] = 1; 


  //METIS provide the METIS SetDefaultOptions routine to set the options to their default values. 
  //After that, the application can just modify the options that is interested in modifying.
  //options[METIS_OPTION_NSEPS] = 10;
  //options[METIS_OPTION_UFACTOR] = 100;
  
  idx_t _edgecut = (idx_t) *edgecut;
  
	real_t *_tpwgts; 
  real_t *_ubvec;  
  for (int i = 0; i < *ncon; i++) {
    _tpwgts[i] = (real_t) tpwgts[i]; 
    _ubvec[i] = (real_t) ubvec[i];         
  }
  
  idx_t *__vtxdist, *_vtxdist;
  idx_t *__xadj, *_xadj;
  idx_t *__adjncy, *_adjncy;

  if (sizeof(PDM_g_num_t) == sizeof(idx_t)) {
    _vtxdist = (idx_t *) vtxdist;
    _xadj    = (idx_t *) xadj;
    _adjncy  = (idx_t *) adjncy;
    __vtxdist = NULL;
    __xadj    = NULL;
    __adjncy  = NULL;
  }
  
  else {
    __vtxdist = (idx_t *) malloc (sizeof(idx_t) * (iSize + 1));
    __xadj    = (idx_t *) malloc (sizeof(idx_t) * (nNode + 1));
    __adjncy  = (idx_t *) malloc (sizeof(idx_t) * nEdge);
    _vtxdist = __vtxdist;
    _xadj    = __xadj;
    _adjncy  = __adjncy;
    
    for (int i = 0; i < iSize + 1; i++) {
      __vtxdist[i] = vtxdist[i]; 
    }
      
    for (int i = 0; i < nNode + 1; i++) {
      __xadj[i] = xadj[i]; 
    }

    for (int i = 0; i < nEdge; i++) {
      __adjncy[i] =  adjncy[i]; 
    }
  }

  idx_t *__vwgt, *_vwgt; 
	idx_t *__adjwgt, *_adjwgt; 

  idx_t *__part, *_part; 

  if (sizeof(int) == sizeof(idx_t)) {
    _vwgt     = (idx_t *) vwgt;
    _adjwgt   = (idx_t *) adjwgt;
    _part     = (idx_t *) part;
    __vwgt    = NULL;
    __adjwgt  = NULL;
    __part    = NULL;
  }

  else {
    if (vwgt != NULL) { 
      __vwgt = (idx_t *) malloc (sizeof(idx_t) * nNode);
      for (int i = 0; i < nNode; i++) {
        __vwgt[i] = vwgt[i]; 
      }      
    }
    else {
      __vwgt = NULL;
    }
    
    if (adjwgt != NULL) { 
      __adjwgt = (idx_t *) malloc (sizeof(idx_t) * nEdge);
      for (int i = 0; i < nEdge; i++) {
        __adjwgt[i] = adjwgt[i]; 
      }      
    }
    else {
      __adjwgt = NULL;
    }

    __part = (idx_t *) malloc (sizeof(idx_t) * nNode);

    _vwgt   = __vwgt;
    _adjwgt = __adjwgt;
    _part   = __part;

  }

  int rval = (int) ParMETIS_V3_PartKway (_vtxdist,
                                         _xadj, 
                                         _adjncy, 
                                         _vwgt, 
                                         _adjwgt, 
                                         &_wgtflag, 
                                         &_numflag, 
                                         &_ncon, 
                                         &_nparts, 
	                                       _tpwgts, 
                                         _ubvec, 
                                         options, 
                                         &_edgecut,
                                        _part, 
	                                      &mpi_comm);
  
  if (sizeof(int) != sizeof(idx_t)) {
    for (int i = 0; i < nNode; i++) {
      part[i] = _part[i]; 
    }      
  }
  
  if (__vtxdist != NULL) {
    free (__vtxdist);
  }
  
  if (__xadj != NULL) {
    free (__xadj);
  }

  if (__adjncy != NULL) {
    free (__adjncy);
  }

  if (__part != NULL) {
    free (__part);
  }

  if (__vwgt != NULL) {
    free (__vwgt);
  }

  if (__adjwgt != NULL) {
    free (__adjwgt);
  }

  return rval;
  
}

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

void  
PDM_SCOTCH_dpart 
(
const PDM_g_num_t dNCell,
const PDM_g_num_t *dDualGraphIdx,
const PDM_g_num_t *dDualGraph,        
const int *cellWeight,
const int *edgeWeight,
const int check,        
const PDM_MPI_Comm comm,
const int  nPart,        
int *part
)
{
  SCOTCH_Dgraph graph;
  SCOTCH_Strat strat;
  int ierr = 0;

  MPI_Comm mpi_comm = *((MPI_Comm *) PDM_MPI_2_mpi_comm (comm));

  ierr = SCOTCH_dgraphInit (&graph, mpi_comm);
  if (ierr) {
    fprintf (stderr,"PPART error : Error in PT-Scotch graph initialization\n");
    exit(1);
  }
    
  SCOTCH_Num _baseval = 0; 
  SCOTCH_Num _vertlocnbr = (SCOTCH_Num) dNCell;
  SCOTCH_Num _vertlocmax = (SCOTCH_Num) dNCell; 
  SCOTCH_Num *_vertloctab, *__vertloctab; 
  SCOTCH_Num *_vendloctab, *__vendloctab; 
  SCOTCH_Num *_veloloctab, *__veloloctab; 
  SCOTCH_Num *_vlblloctab = NULL; 
  SCOTCH_Num _edgelocnbr = (SCOTCH_Num) dDualGraphIdx[dNCell];
  SCOTCH_Num _edgelocsiz = (SCOTCH_Num) dDualGraphIdx[dNCell];
  SCOTCH_Num *_edgeloctab, *__edgeloctab; 
  SCOTCH_Num *_edgegsttab = NULL; 
  SCOTCH_Num *_edloloctab, *__edloloctab;
  SCOTCH_Num *_part, *__part;
  
  if (sizeof(PDM_g_num_t) == sizeof(SCOTCH_Num)) {
    _vertloctab = (SCOTCH_Num *) dDualGraphIdx;
    _vendloctab = (SCOTCH_Num *) dDualGraphIdx + 1;
    _edgeloctab = (SCOTCH_Num *) dDualGraph;
    
    __vertloctab = NULL;
    __vendloctab = NULL;
    __edgeloctab = NULL;
  }
  
  else {
    __vertloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * (_vertlocnbr + 1));
    __vendloctab = __vertloctab + 1; 
    __edgeloctab  = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _edgelocsiz);

    for (int i = 0; i < _vertlocnbr + 1; i++) {
      __vertloctab[i] = dDualGraphIdx[i]; 
    }

    for (int i = 0; i < _edgelocsiz; i++) {
      __edgeloctab[i] = dDualGraph[i]; 
    }

    _vertloctab = __vertloctab;
    _vendloctab = __vendloctab;
    _edgeloctab = __edgeloctab;

  }

  if (sizeof(int) == sizeof(SCOTCH_Num)) {
    
    _veloloctab = (SCOTCH_Num *) cellWeight; 
    _edloloctab = (SCOTCH_Num *) edgeWeight;
    _part = part;
    
    __veloloctab = NULL; 
    __edloloctab = NULL;
    __part = NULL;
    
  }
  
  else {
          
    __veloloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _vertlocnbr);
    __part       = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _vertlocnbr);
    __edloloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * _edgelocsiz);
      
    for (int i = 0; i < _vertlocnbr; i++) {
      __veloloctab[i] = cellWeight[i]; 
    }

    for (int i = 0; i < _edgelocsiz; i++) {
      __edloloctab[i] = edgeWeight[i]; 
    }

    _veloloctab = __veloloctab; 
    _edloloctab = __edloloctab;
    _part = __part;
    
  }

  ierr = SCOTCH_dgraphBuild  (&graph,
                              _baseval, 
                              _vertlocnbr, 
                              _vertlocmax, 
                              _vertloctab, 
                              _vendloctab,
                              _veloloctab, 
                              _vlblloctab, 
                              _edgelocnbr, 
                              _edgelocsiz, 
                              _edgeloctab, 
                              _edgegsttab, 
                              _edloloctab);

  if (ierr) {
    fprintf(stderr, "PPART error : Error in SCOTCH_dgraphBuild\n");
    exit(1);
  }

  /* Checks graph */

  if (check) {

    ierr = SCOTCH_dgraphCheck (&graph);
  }

  if (ierr) {
    fprintf(stderr, "PPART error : Error in PT-Scotch graph check\n");
    exit(1);
  }

  /* Partitioning strategy : */

  SCOTCH_stratInit(&strat);

  const SCOTCH_Num _nPart = (SCOTCH_Num) nPart;

  ierr = SCOTCH_dgraphPart(&graph,
                               _nPart, /* Nombre de partitions demande */
                               &strat,
                               _part);    /* parts[i] donne le numero */

  SCOTCH_stratExit(&strat);
  SCOTCH_dgraphExit(&graph);

  if (__part != NULL) {
    for (int i = 0; i < _vertlocnbr; i++) {
      part[i] = __part[i]; 
    }
    free (__part);
  }

  if (__vertloctab != NULL) {
    free (__vertloctab);
  }
  
  if (__vendloctab != NULL) {
    free (__vendloctab);
  }
  
  if (__edgeloctab != NULL) {
    free (__edgeloctab);
  }

  if (__veloloctab != NULL) {
    free (__veloloctab);
  }
  
  if (__edloloctab != NULL) {
    free (__edloloctab);
  }
  
}

#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */
