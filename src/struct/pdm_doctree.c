/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_doctree_priv.h"
#include "pdm_doctree.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_octree_seq.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/



/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_doctree_t*
PDM_doctree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
)
{
  PDM_doctree_t* doct = (PDM_doctree_t *) malloc(sizeof(PDM_doctree_t));

  doct->comm = comm;
  doct->dim  = dim;

  doct->global_depth_max          = 5;
  doct->global_points_in_leaf_max = 60;

  doct->global_octree = NULL;
  doct->local_octree  = NULL;
  doct->shmem_octree  = NULL;

  doct->comm_shared   = PDM_MPI_COMM_NULL;

  PDM_UNUSED(global_extents);

  return doct;
}

void
PDM_doctree_point_set
(
 PDM_doctree_t   *doct,
 const int        i_point_cloud,
 const int        n_points,
 const double    *coords
)
{
  PDM_UNUSED(doct);
  PDM_UNUSED(i_point_cloud);
  PDM_UNUSED(n_points);
  PDM_UNUSED(coords);


}



#ifdef __cplusplus
}
#endif /* __cplusplus */
