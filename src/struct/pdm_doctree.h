#ifndef PDM_DOCTREE_H
#define PDM_DOCTREE_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef  __cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_doctree_t PDM_doctree_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/



PDM_doctree_t*
PDM_doctree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
);


void
PDM_doctree_point_set
(
 PDM_doctree_t   *doct,
 const int        i_point_cloud,
 const int        n_points,
 const double    *coords
);




#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DOCTREE_H */
