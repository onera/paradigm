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
 PDM_MPI_Comm              comm,
 int                       dim,
 int                       n_part_cloud,
 double                   *global_extents,
 PDM_doctree_local_tree_t  local_tree_kind
);

void
PDM_doctree_build
(
 PDM_doctree_t     *doct
);


void
PDM_doctree_point_set
(
 PDM_doctree_t     *doct,
 const int          i_part_cloud,
 const int          n_points,
 const int         *pts_init_location,
 const PDM_g_num_t *pts_g_num,
 const double      *pts_coords
);


void
PDM_doctree_free
(
  PDM_doctree_t   *doct
);



#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DOCTREE_H */
