
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_handles.h"

#include "pdm_dmesh_partitioning.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/**
 * \def _PDM_part_MIN(a,b)
 * Computes the minimum of \a x and \a y.
 *
 */

#define _PDM_part_MIN(a,b) ((a) > (b) ? (b) : (a))

/**
 * \def _PDM_part_MAX(a,b)
 * Computes the maximum of \a x and \a y.
 *
 */

#define _PDM_part_MAX(a,b) ((a) < (b) ? (b) : (a))

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _dmesh_partitioning_t
 * \brief  Define a global numberring
 *
 */
typedef struct {


} _dmesh_partitioning_t;

/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *_dmps   = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
static _dmesh_partitioning_t *
_get_from_id
(
 int  dmpartitioning_id
)
{
  _dmesh_partitioning_t *dmesh_partitioning = (_dmesh_partitioning_t *) PDM_Handles_get (_dmps, dmpartitioning_id);

  if (dmesh_partitioning == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "pdm_dmesh_partitioning error : Bad _dmesh_partitioning identifier\n");
    exit(1);
  }

  return dmesh_partitioning;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
int
PDM_dmesh_partitioning_create
(
 const PDM_MPI_Comm              comm,
 const PDM_partitioning_method_t split_method
)
{
  /*
   * Search a gnum_from_hash_values free id
   */
  if (_dmps == NULL) {
    _dmps = PDM_Handles_create (4);
  }

  _dmesh_partitioning_t *_dmp = (_dmesh_partitioning_t *) malloc(sizeof(_dmesh_partitioning_t));
  int id = PDM_Handles_store (_dmps, _dmp);

  return id;
}


/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_set_from_dmesh
(
 const int dmesh_partitioning_id,
 const int dmesh_id
)
{
  // ---> Depuis l'interface du dmesh
  // Recopie
  // PDM_dmesh_dims_get(dmesh_id, &)
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
