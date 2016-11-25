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
);


int 
PDM_METIS_PartGraphRecursive
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
);

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
);

int
PDM_METIS_SetDefaultOptions
(
int *options
);

#endif
    
#ifdef PDM_HAVE_PTSCOTCH

typedef void* PDM_SCOTCH_Dgraph;
typedef void* PDM_SCOTCH_Graph;
typedef void* PDM_SCOTCH_Strat;

int  
PDM_SCOTCH_dgraphInit   
(
PDM_SCOTCH_Dgraph, 
PDM_MPI_Comm comm
);

int  
PDM_SCOTCH_graphInit   
(
PDM_SCOTCH_Graph 
);

int  
PDM_SCOTCH_dgraphBuild  
(
PDM_SCOTCH_Dgraph, 
const int  baseval, 
const PDM_g_num_t vertlocnbr, 
const PDM_g_num_t vertlocmax, 
PDM_g_num_t * const vertloctab, 
PDM_g_num_t * const vendloctab,
int * const veloloctab, // Poids cellules */
PDM_g_num_t * const vlblloctab, 
const PDM_g_num_t edgelocnbr, 
const PDM_g_num_t edgelocsiz, 
PDM_g_num_t * const edgeloctab, 
int  * const edloloctab // Poids faces */
);

int  
PDM_SCOTCH_dgraphCheck  
(
const PDM_SCOTCH_Dgraph graph
);

int  
PDM_SCOTCH_graphCheck  
(
const PDM_SCOTCH_Graph graph
);

int  
PDM_SCOTCH_dgraphPart   
(
PDM_SCOTCH_Dgraph graph, 
const int nPart, 
PDM_SCOTCH_Strat stratptr,
int *part
);

int  
PDM_SCOTCH_graphPart   
(
PDM_SCOTCH_Graph graph, 
const int nPart, 
PDM_SCOTCH_Strat stratptr,
int *part
);

void 
PDM_SCOTCH_stratExit    
(
PDM_SCOTCH_Strat
);

void 
PDM_SCOTCH_dgraphExit   
(
PDM_SCOTCH_Dgraph
);

void 
PDM_SCOTCH_graphExit   
(
PDM_SCOTCH_Graph graphptr
);
    
void 
PDM_SCOTCH_stratExit   
(
PDM_SCOTCH_Strat stratptr
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

PDM_SCOTCH_Graph 
PDM_SCOTCH_GraphAlloc
(
void
);

PDM_SCOTCH_Graph 
PDM_SCOTCH_GraphFree
(
PDM_SCOTCH_Graph graph
);

PDM_SCOTCH_Strat 
PDM_SCOTCH_StratAlloc
(
void
);

PDM_SCOTCH_Strat 
PDM_SCOTCH_StratFree
(
PDM_SCOTCH_Strat strat
);


#endif


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_MPI_EXT_DEPENDENCIES_H */

