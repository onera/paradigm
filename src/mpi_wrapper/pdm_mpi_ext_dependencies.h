#ifndef __PDM_MPI_EXT_DEPENDENCIES_H__
#define __PDM_MPI_EXT_DEPENDENCIES_H__

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

#ifdef __cplusplus
extern "C" {
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
const int *n_parts,
const double *tpwgts,
const double *ubvec,
const int *edgecut,
      int *part,
const PDM_MPI_Comm comm
);

#endif

#ifdef PDM_HAVE_PTSCOTCH

void
PDM_SCOTCH_dpart
(
const PDM_g_num_t dn_cell,
const PDM_g_num_t *ddual_graph_idx,
const PDM_g_num_t *ddual_graph,
const int *cell_weight,
const int *edgeWeight,
const int check,
const PDM_MPI_Comm comm,
const int n_part,
int *part
);

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_EXT_DEPENDENCIES_H__ */

