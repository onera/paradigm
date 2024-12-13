/*
 * \file
 * \author equemera
 *
 * \date December 6, 2016, 9:50 AM
 */

#ifndef PDM_EXT_WRAPPER_H
#define	PDM_EXT_WRAPPER_H

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


#ifdef	__cplusplus
extern "C" {
#endif


#ifdef PDM_HAVE_PARMETIS

int
PDM_METIS_PartGraphRecursive
(
int    *n_vtxs,
int    *ncon,
int    *xadj,
int    *adjncy,
int    *vwgt,
int    *adjwgt,
int    *n_parts,
double *tpwgts,
double *ubvec,
int    *edgecut,
int    *part
);

int
PDM_METIS_PartGraphKway
(
int    *n_vtxs,
int    *ncon,
int    *xadj,
int    *adjncy,
int    *vwgt,
int    *adjwgt,
int    *n_parts,
double *tpwgts,
double *ubvec,
int    *edgecut,
int    *part
);

#endif

#ifdef PDM_HAVE_PTSCOTCH

void
PDM_SCOTCH_part
(
const int  n_cell,
      int *dual_graph_idx,
      int *dual_graph,
      int *cell_weight,
      int *edge_weight,
      int  check,
const int  n_part,
      int *part
);

#endif

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_EXT_WRAPPER_H */

