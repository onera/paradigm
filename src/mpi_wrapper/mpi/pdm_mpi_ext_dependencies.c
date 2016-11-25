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
PDM_g_num_t *vtxdist, 
PDM_g_num_t *xadj, 
PDM_g_num_t *adjncy, 
int *vwgt, 
int *adjwgt, 
int *wgtflag, 
int *numflag, 
int *ncon, 
int *nparts, 
double *tpwgts, 
double *ubvec, 
int options[5], 
int *edgecut, 
int *part, 
PDM_MPI_Comm comm
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

  idx_t _options[5] = {options[0],
                       options[1], 
                       options[2], 
                       options[3], 
                       options[4]};
  
  idx_t _edgecut = edgecut;
  
	real_t *_tpwgts; 
  real_t *_ubvec;  
  for (int i = 0; i < *ncon; i++) {
    _tpwgts[i] = tpwgts[i]; 
    _ubvec[i] = ubvec[i];         
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
                                         _options, 
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

int PDM_METIS_PartGraphRecursive
(
int *nvtxs, 
int *ncon, 
int *xadj, 
int *adjncy, 
int *vwgt, 
int *adjwgt, 
int *nparts, 
double *tpwgts, 
double *ubvec, 
int options[5], 
int *edgecut, 
int *part
)
{

  idx_t _nvtxs = *nvtxs;

  idx_t _ncon = *ncon;
  
  int nEdge = xadj[_nvtxs];
  
  idx_t _nparts = *nparts; 

  idx_t _options[5] = {options[0],
                       options[1], 
                       options[2],
                       options[3],
                       options[4]};
  
  idx_t _edgecut = edgecut;
  
	real_t *_tpwgts; 
  real_t *_ubvec;  
  for (int i = 0; i < *ncon; i++) {
    _tpwgts[i] = tpwgts[i]; 
    _ubvec[i] = ubvec[i];         
  }

  int *_vsize = NULL;

  idx_t *__xadj, *_xadj;
  idx_t *__adjncy, *_adjncy;
  idx_t *__vwgt, *_vwgt; 
	idx_t *__adjwgt, *_adjwgt; 
  idx_t *__part, *_part; 

  if (sizeof(int) == sizeof(idx_t)) {
    _vwgt     = (idx_t *) vwgt;
    _adjwgt   = (idx_t *) adjwgt;
    _part     = (idx_t *) part;
    _xadj     = (idx_t *) xadj;
    _adjncy   = (idx_t *) adjncy;

    __vwgt    = NULL;
    __adjwgt  = NULL;
    __part    = NULL;
    __xadj    = NULL;
    __adjncy  = NULL;
  }

  else {
    __xadj    = (idx_t *) malloc (sizeof(idx_t) * (_nvtxs + 1));
    __adjncy  = (idx_t *) malloc (sizeof(idx_t) * nEdge);
    _xadj    = __xadj;
    _adjncy  = __adjncy;
      
    for (int i = 0; i < _nvtxs + 1; i++) {
      __xadj[i] = xadj[i]; 
    }

    for (int i = 0; i < nEdge; i++) {
      __adjncy[i] =  adjncy[i]; 
    }

    if (vwgt != NULL) { 
      __vwgt = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);
      for (int i = 0; i < _nvtxs; i++) {
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

    __part = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);

    _vwgt   = __vwgt;
    _adjwgt = __adjwgt;
    _part   = __part;

  }

  int rval = (int) METIS_PartGraphRecursive (&_nvtxs, 
                                             &_ncon, 
                                              _xadj, 
                                              _adjncy, 
                                              _vwgt, 
                                              _vsize, 
                                              _adjwgt, 
                                              &_nparts, 
                                              _tpwgts, 
                                              _ubvec, 
                                              _options, 
                                              &_edgecut, 
                                              &_part);

    if (sizeof(int) != sizeof(idx_t)) {
    for (int i = 0; i < _nvtxs; i++) {
      part[i] = _part[i]; 
    }      
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


int 
PDM_METIS_PartGraphKway
(
int *nvtxs, 
int *ncon, 
int *xadj, 
int *adjncy, 
int *vwgt, 
int *adjwgt, 
int *nparts, 
double *tpwgts, 
double *ubvec, 
int options[5], 
int *edgecut, 
int *part
)
{

  idx_t _nvtxs = *nvtxs;

  idx_t _ncon = *ncon;
  
  int nEdge = xadj[_nvtxs];
  
  idx_t _nparts = *nparts; 

  idx_t _options[5] = {options[0],
                       options[1], 
                       options[2],
                       options[3],
                       options[4]};
  
  idx_t _edgecut = edgecut;
  
	real_t *_tpwgts; 
  real_t *_ubvec;  
  for (int i = 0; i < *ncon; i++) {
    _tpwgts[i] = tpwgts[i]; 
    _ubvec[i] = ubvec[i];         
  }

  int *_vsize = NULL;

  idx_t *__xadj, *_xadj;
  idx_t *__adjncy, *_adjncy;
  idx_t *__vwgt, *_vwgt; 
	idx_t *__adjwgt, *_adjwgt; 
  idx_t *__part, *_part; 

  if (sizeof(int) == sizeof(idx_t)) {
    _vwgt     = (idx_t *) vwgt;
    _adjwgt   = (idx_t *) adjwgt;
    _part     = (idx_t *) part;
    _xadj     = (idx_t *) xadj;
    _adjncy   = (idx_t *) adjncy;

    __vwgt    = NULL;
    __adjwgt  = NULL;
    __part    = NULL;
    __xadj    = NULL;
    __adjncy  = NULL;
  }

  else {
    __xadj    = (idx_t *) malloc (sizeof(idx_t) * (_nvtxs + 1));
    __adjncy  = (idx_t *) malloc (sizeof(idx_t) * nEdge);
    _xadj    = __xadj;
    _adjncy  = __adjncy;
      
    for (int i = 0; i < _nvtxs + 1; i++) {
      __xadj[i] = xadj[i]; 
    }

    for (int i = 0; i < nEdge; i++) {
      __adjncy[i] =  adjncy[i]; 
    }

    if (vwgt != NULL) { 
      __vwgt = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);
      for (int i = 0; i < _nvtxs; i++) {
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

    __part = (idx_t *) malloc (sizeof(idx_t) * _nvtxs);

    _vwgt   = __vwgt;
    _adjwgt = __adjwgt;
    _part   = __part;

  }

  int rval = (int) METIS_PartGraphKway (&_nvtxs, 
                                        &_ncon, 
                                         _xadj, 
                                         _adjncy, 
                                         _vwgt, 
                                         _vsize, 
                                         _adjwgt, 
                                        &_nparts, 
                                         _tpwgts, 
                                         _ubvec, 
                                         _options, 
                                        &_edgecut, 
                                        &_part);

  if (__part != NULL) {
    for (int i = 0; i < _nvtxs; i++) {
      part[i] = __part[i]; 
    }
    free (__part);
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


int
PDM_METIS_SetDefaultOptions
(
int *options
)
{
  idx_t _options[METIS_NOPTIONS]; /* Options */
  int rval = (int) METIS_SetDefaultOptions(_options);
  
  for (int i = 0; i < METIS_NOPTIONS; i++) {
    options[i] = _options[i];
  }
  
  return rval;
}

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

int  
PDM_SCOTCH_dgraphInit   
(
PDM_SCOTCH_Dgraph graphptr,
PDM_MPI_Comm          proccomm)             /* Communicator to be used for all communications */
{
  MPI_Comm mpi_comm = *((MPI_Comm *) PDM_MPI_2_mpi_comm (proccomm));
  return SCOTCH_dgraphInit ((SCOTCH_Dgraph *) graphptr, mpi_comm);
}


int  
PDM_SCOTCH_graphInit   
(
PDM_SCOTCH_Graph graphptr
)
{
  return SCOTCH_dgraphInit ((SCOTCH_Graph *) graphptr);
}


int  
PDM_SCOTCH_dgraphBuild  
(
PDM_SCOTCH_Dgraph graphptr, 
const int  baseval, 
const PDM_g_num_t vertlocnbr, 
const PDM_g_num_t vertlocmax, 
PDM_g_num_t * const vertloctab, 
PDM_g_num_t * const vendloctab,
int * const veloloctab, // Poids cellules */
const PDM_g_num_t edgelocnbr, 
const PDM_g_num_t edgelocsiz, 
PDM_g_num_t * const edgeloctab, 
int  * const edloloctab // Poids faces */
)
{
  
  SCOTCH_Num _baseval = (SCOTCH_Num) baseval; 
  SCOTCH_Num _vertlocnbr = (SCOTCH_Num) vertlocnbr;
  SCOTCH_Num _vertlocmax = (SCOTCH_Num) vertlocmax; 
  SCOTCH_Num *_vertloctab, *__vertloctab; 
  SCOTCH_Num *_vendloctab, *__vendloctab; 
  SCOTCH_Num *_veloloctab, *__veloloctab; 
  SCOTCH_Num *_vlblloctab = NULL; 
  SCOTCH_Num _edgelocnbr = (SCOTCH_Num) edgelocnbr;
  SCOTCH_Num _edgelocsiz = (SCOTCH_Num) edgelocsiz;
  SCOTCH_Num *_edgeloctab, *__edgeloctab; 
  SCOTCH_Num *_edgegsttab = NULL; 
  SCOTCH_Num *_edloloctab, *__edloloctab;
  
  if (sizeof(PDM_g_num_t) == sizeof(SCOTCH_Num)) {
    _vertloctab = vertloctab;
    _vendloctab = vendloctab;
    _edgeloctab = edgeloctab;
    
    __vertloctab = NULL;
    __vendloctab = NULL;
    __edgeloctab = NULL;
  }
  
  else {
    __vertloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * vertlocnbr);
    __vendloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * vertlocnbr);
    __edgeloctab  = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * edgelocsiz);

    for (int i = 0; i < vertlocnbr; i++) {
      __vertloctab[i] = vertloctab[i]; 
    }
      
    for (int i = 0; i < vertlocnbr; i++) {
      __vendloctab[i] = vendloctab[i]; 
    }

    for (int i = 0; i < edgelocsiz; i++) {
      __edgeloctab[i] = edgeloctab[i]; 
    }

    _vertloctab = __vertloctab;
    _vendloctab = __vendloctab;
    _edgeloctab = __edgeloctab;

  }


  if (sizeof(int) == sizeof(SCOTCH_Num)) {
    
    _veloloctab = veloloctab; 
    _edloloctab = edloloctab;
    
    __veloloctab = NULL; 
    __edloloctab = NULL;
    
  }
  
  else {
          
    __veloloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * vertlocnbr);
    __edloloctab = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * edgelocsiz);
      
    for (int i = 0; i < vertlocnbr; i++) {
      __veloloctab[i] = veloloctab[i]; 
    }

    for (int i = 0; i < edgelocsiz; i++) {
      __edloloctab[i] = _edloloctab[i]; 
    }

    _veloloctab = __veloloctab; 
    _edloloctab = __edloloctab;
    
  }
    
    
  return SCOTCH_dgraphBuild  ((SCOTCH_Dgraph *) graphptr,
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
}

int 
SCOTCH_dgraphBuild  (
SCOTCH_Dgraph * const,
        const SCOTCH_Num, 
        const SCOTCH_Num, 
        const SCOTCH_Num, 
        SCOTCH_Num * const, 
        SCOTCH_Num * const, 
        SCOTCH_Num * const, 
        SCOTCH_Num * const, 
        const SCOTCH_Num, 
        const SCOTCH_Num, 
        SCOTCH_Num * const, 
        SCOTCH_Num * const, 
        SCOTCH_Num * const);


int  
PDM_SCOTCH_dgraphCheck  
(
const PDM_SCOTCH_Dgraph graphptr 
)
{
  return SCOTCH_dgraphCheck ((SCOTCH_Dgraph *) graphptr);
}


int  
PDM_SCOTCH_graphCheck  
(
const PDM_SCOTCH_Graph graphptr 
)
{
  return SCOTCH_graphCheck ((SCOTCH_Graph *) graphptr);
}


int  
PDM_SCOTCH_dgraphPart   
(
PDM_SCOTCH_Dgraph graph, 
const int nPart,
const int nVert,
PDM_SCOTCH_Strat stratptr,
int *part
)
{

  SCOTCH_Num _nPart = nPart;
  SCOTCH_Num *_part, *__part;
  
  if (sizeof(SCOTCH_Num) == sizeof(int)) {
    _part = part;
    __part = NULL;
  }
  
  else {
    __part = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * nPart);
    _part = __part;
  }

  int rval = (int) SCOTCH_dgraphPart ((SCOTCH_Dgraph *) graph,
                                     _nPart,
                                     stratptr,
                                     _part);
  
  if (__part != NULL) {
    for (int i = 0; i < nPart; i++) {
      part[i] = __part[i];
    }
    free (__part);
  }
  
  return rval;
  
}


int  
PDM_SCOTCH_graphPart   
(
PDM_SCOTCH_Graph graph, 
const int nPart,
const int nVert,
PDM_SCOTCH_Strat stratptr,
int *part
)
{

  SCOTCH_Num _nPart = nPart;
  SCOTCH_Num *_part, *__part;
  
  if (sizeof(SCOTCH_Num) == sizeof(int)) {
    _part = part;
    __part = NULL;
  }
  
  else {
    __part = (SCOTCH_Num *) malloc (sizeof(SCOTCH_Num) * nPart);
    _part = __part;
  }

  int rval = (int) SCOTCH_graphPart ((SCOTCH_Graph *) graph,
                                     _nPart,
                                     stratptr,
                                     _part);
  
  if (__part != NULL) {
    for (int i = 0; i < nPart; i++) {
      part[i] = __part[i];
    }
    free (__part);
  }
  
  return rval;
  
}


void 
PDM_SCOTCH_dgraphExit   
(
PDM_SCOTCH_Dgraph graphptr
)
{
 SCOTCH_dgraphExit ((SCOTCH_Dgraph *) graphptr);
}
    
void 
PDM_SCOTCH_graphExit   
(
PDM_SCOTCH_Graph graphptr
)
{
 SCOTCH_graphExit ((SCOTCH_Graph *) graphptr);
}
    
void 
PDM_SCOTCH_stratExit   
(
PDM_SCOTCH_Strat stratptr
)
{
 SCOTCH_graphExit ((SCOTCH_Strat *) stratptr);
}


PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphAlloc
(
void
)
{
  return malloc (sizeof(SCOTCH_Dgraph));
}


PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphFree
(
PDM_SCOTCH_Dgraph graph
)
{
  free (graph);
  return NULL;
}

PDM_SCOTCH_Graph 
PDM_SCOTCH_GraphAlloc
(
void
)
{
  return malloc (sizeof(SCOTCH_Graph));
}


PDM_SCOTCH_Graph 
PDM_SCOTCH_GraphFree
(
PDM_SCOTCH_Graph graph
)
{
  free (graph);
  return NULL;
}


PDM_SCOTCH_Strat 
PDM_SCOTCH_StratAlloc
(
void
)
{
  return malloc (sizeof(SCOTCH_Strat));
}


PDM_SCOTCH_Strat 
PDM_SCOTCH_StratFree
(
PDM_SCOTCH_Strat strat
)
{
  free (strat);
  return NULL;
}

#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */
