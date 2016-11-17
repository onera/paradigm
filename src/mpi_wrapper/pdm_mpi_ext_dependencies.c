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
#endif
#ifdef PDM_HAVE_PTSCOTCH
#include <ptscotch.h>
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
  MPI_Comm mpi_comm = *((MPI_Comm *) PDM_MPI_2_mpi_comm (comm));
  
  return ParMETIS_V3_PartKway (vtxdist, xadj, adjncy, vwgt, 
	                                 adjwgt, wgtflag, numflag, ncon, nparts, 
	                                 tpwgts, ubvec, options, edgecut, part, 
	                                 &mpi_comm);
}

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

int  
PDM_SCOTCH_dgraphInit   
(
PDM_SCOTCH_Dgraph graphptr,
MPI_Comm          proccomm)             /* Communicator to be used for all communications */
{
  MPI_Comm mpi_comm = *((MPI_Comm *) PDM_MPI_2_mpi_comm (proccomm));
  return SCOTCH_dgraphInit ((SCOTCH_Dgraph *) graphptr, mpi_comm);
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
  return SCOTCH_dgraphBuild  ((SCOTCH_Dgraph *) graphptr,
                              baseval, 
                              vertlocnbr, 
                              vertlocmax, 
                              vertloctab, 
                              vendloctab,
                              veloloctab, 
                              vlblloctab, 
                              edgelocnbr, 
                              edgelocsiz, 
                              edgeloctab, 
                              edgegsttab, 
                              edloloctab);
}


int  
PDM_SCOTCH_dgraphCheck  
(
const PDM_SCOTCH_Dgraph graphptr 
)
{
  return SCOTCH_dgraphCheck ((SCOTCH_Dgraph *) graphptr);
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
  return SCOTCH_dgraphPart ((SCOTCH_Dgraph *) graphptr,
                             partnbr,
                             stratptr,
                             termlocatab);
}
    
void 
PDM_SCOTCH_dgraphExit   
(
PDM_SCOTCH_Dgraph graphptr
)
{
 SCOTCH_dgraphExit ((SCOTCH_Dgraph *) graphptr);
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

#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */
