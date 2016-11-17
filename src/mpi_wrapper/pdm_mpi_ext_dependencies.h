/* 
 * File:   pdm_mpi_ext_dependencies.h
 * Author: equemera
 *
 * Created on November 16, 2016, 1:49 PM
 */

#ifndef PDM_MPI_EXT_DEPENDENCIES_H
#define	PDM_MPI_EXT_DEPENDENCIES_H

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

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef PDM_HAVE_PARMETIS
#include <metis.h>
#endif
#ifdef PDM_HAVE_PTSCOTCH
#include <scotch.h>
#endif

#ifdef	__cplusplus
extern "C" {
#endif

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
PDM_MPI_Comm comm
);

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

typedef void* PDM_SCOTCH_Dgraph;

int  
PDM_SCOTCH_dgraphInit   
(
PDM_SCOTCH_Dgraph, 
PDM_MPI_Comm comm
);

int  
PDM_SCOTCH_dgraphBuild  
(
PDM_SCOTCH_Dgraph, 
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
SCOTCH_Num * const
);

int  
PDM_SCOTCH_dgraphCheck  
(
const PDM_SCOTCH_Dgraph
);

int  
PDM_SCOTCH_dgraphPart   
(
PDM_SCOTCH_Dgraph, 
const SCOTCH_Num, 
SCOTCH_Strat * const, 
SCOTCH_Num * const
);

void 
PDM_SCOTCH_stratExit    
(
SCOTCH_Strat
);

void 
PDM_SCOTCH_dgraphExit   
(
PDM_SCOTCH_Dgraph
);

PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphAlloc
(
void
);

PDM_SCOTCH_Dgraph 
PDM_SCOTCH_DgraphFree
(
PDM_SCOTCH_Dgraph graphptr
);

#endif


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_MPI_EXT_DEPENDENCIES_H */

