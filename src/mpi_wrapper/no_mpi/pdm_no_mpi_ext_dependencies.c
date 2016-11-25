/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"


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
idx_t *vtxdist, 
idx_t *xadj, 
idx_t *adjncy, 
idx_t *vwgt, 
idx_t *adjwgt, 
idx_t *wgtflag, 
idx_t *numflag, 
idx_t *ncon, 
idx_t *nparts, 
real_t *tpwgts, 
real_t *ubvec, 
idx_t *options, 
idx_t *edgecut, 
idx_t *part, 
PDM_MPI_Comm comm)
{
  vtxdist; 
  xadj; 
  adjncy; 
  vwgt; 
  adjwgt; 
  wgtflag; 
  numflag; 
  ncon; 
  nparts; 
  tpwgts; 
  ubvec; 
  options; 
  edgecut; 
  part; 
  comm;

  fprintf(stderr,"PDM_ParMETIS_V3_PartKway : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

int  
PDM_SCOTCH_dgraphInit   
(
PDM_SCOTCH_Dgraph graphptr,
PDM_MPI_Comm          proccomm)             /* Communicator to be used for all communications */
{
  graphptr;
  proccomm;
  
  fprintf(stderr,"PDM_SCOTCH_dgraphInit : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}


int  
PDM_SCOTCH_dgraphBuild  
(
PDM_SCOTCH_Dgraph graphptr, 
const SCOTCH_Num baseval, 
const SCOTCH_Num vertlocnbr, 
const SCOTCH_Num vertlocmax, 
SCOTCH_Num * const vertloctab, 
SCOTCH_Num * const vendloctab,
SCOTCH_Num * const veloloctab, 
SCOTCH_Num * const vlblloctab, 
const SCOTCH_Num edgelocnbr, 
const SCOTCH_Num edgelocsiz, 
SCOTCH_Num * const edgeloctab, 
SCOTCH_Num * const edgegsttab, 
SCOTCH_Num * const edloloctab
)
{
  graphptr; 
  baseval; 
  vertlocnbr; 
  vertlocmax; 
  vertloctab; 
  vendloctab;
  veloloctab; 
  vlblloctab; 
  edgelocnbr; 
  edgelocsiz; 
  edgeloctab; 
  edgegsttab; 
  edloloctab;

  fprintf(stderr,"PDM_SCOTCH_dgraphBuild : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}


int  
PDM_SCOTCH_dgraphCheck  
(
const PDM_SCOTCH_Dgraph graphptr 
)
{
  graphptr;
  
  fprintf(stderr,"PDM_SCOTCH_dgraphCheck : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}


int  
PDM_SCOTCH_dgraphPart   
(
PDM_SCOTCH_Dgraph graphptr, 
const SCOTCH_Num partnbr, 
SCOTCH_Strat * const stratptr, 
SCOTCH_Num * const termlocatab
)
{
  graphptr; 
  partnbr; 
  stratptr; 
  termlocatab;

  fprintf(stderr,"PDM_SCOTCH_dgraphPart : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}
    
void 
PDM_SCOTCH_dgraphExit   
(
PDM_SCOTCH_Dgraph graphptr
)
{
  graphptr;
  
  fprintf(stderr,"PDM_SCOTCH_dgraphExit : Unavailable function with pdm_no_mpi library\n" );
  abort();
}


PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphAlloc
(
void
)
{
  fprintf(stderr,"PDM_SCOTCH_DgraphAlloc : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return (PDM_SCOTCH_Dgraph) 0;
}


PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphFree
(
PDM_SCOTCH_Dgraph graph
)
{
  graph;
  
  fprintf(stderr,"PDM_SCOTCH_DgraphFree : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return (PDM_SCOTCH_Dgraph) 0;
}

#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */

