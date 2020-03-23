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
#include "pdm_printf.h"
#include "pdm_error.h"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
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
const int *n_parts,
const double *tpwgts,
const double *ubvec,
const int *edgecut,
int *part,
const PDM_MPI_Comm comm
)
{
  vtxdist;
  xadj;
  adjncy;
  vwgt;
  adjwgt;
  wgtflag;
  numflag;
  ncon;
  n_parts;
  tpwgts;
  ubvec;
  edgecut;
  part;
  comm;

  PDM_error(__FILE__, __LINE__, 0,"PDM_ParMETIS_V3_PartKway : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 0;
}


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
)
{
  dn_cell;
  ddual_graph_idx;
  ddual_graph;
  cell_weight;
  edgeWeight;
  check;
  comm;
  n_part;
  part;

  PDM_error(__FILE__, __LINE__, 0,"PDM_SCOTCH_dpart : Unavailable function with pdm_no_mpi library\n" );
  abort();
}


#endif



#ifdef __cplusplus
}
#endif /* __cplusplus */

#if defined(__clang__)
#pragma clang diagnostic pop
#endif
