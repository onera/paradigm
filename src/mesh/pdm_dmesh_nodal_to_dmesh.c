
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh.h"

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

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_Dist_t
 * \brief  Distance to a mesh surface structure
 *
 */

typedef struct {

  PDM_MPI_Comm      comm;                    /*!< MPI communicator */
  PDM_ownership_t   owner;                   /*!< Which have the responsabilities of results */
  PDM_bool_t        results_is_getted;       /*!< Flags to indicate if result is getted      */

  int               n_mesh;                  /*!< Number of meshes to manages                */

  int              *dmesh_nodal_ids;         /*!< Identifier list of dmesh_nodal struture
                                              *   (size = n_mesh )                           */

  int              *dmesh_ids;               /*!< Identifier list of dmesh structure
                                              *   (size = n_mesh)                            */

} _PDM_dmesh_nodal_to_dmesh;


/*============================================================================
 * Global variable
 *============================================================================*/

static PDM_Handles_t *dmn_to_dm_handles = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Return dmn_to_dm object from it identifier
 *
 * \param [in]   id        dmn_to_dm identifier
 *
 */
static _PDM_dmesh_nodal_to_dmesh *
_get_from_id
(
 int  id
)
{
  _PDM_dmesh_nodal_to_dmesh *dmn_to_dm = (_PDM_dmesh_nodal_to_dmesh *) PDM_Handles_get (dmn_to_dm_handles, id);

  if (dmn_to_dm == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_dmesh_nodal_to_dmesh error : Bad identifier\n");
  }

  return dmn_to_dm;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

int
PDM_dmesh_nodal_to_dmesh_create
(
      int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
)
{
  if (dmn_to_dm_handles == NULL) {
    dmn_to_dm_handles = PDM_Handles_create (4);
  }

  _PDM_dmesh_nodal_to_dmesh *dmn_to_dm = (_PDM_dmesh_nodal_to_dmesh *) malloc(sizeof(_PDM_dmesh_nodal_to_dmesh));

  int id = PDM_Handles_store (dmn_to_dm_handles, dmn_to_dm);

  dmn_to_dm->comm              = comm;
  dmn_to_dm->owner             = owner;
  dmn_to_dm->results_is_getted = PDM_FALSE;
  dmn_to_dm->n_mesh            = n_mesh;

  dmn_to_dm->dmesh_nodal_ids = (int *) malloc(sizeof(int));
  dmn_to_dm->dmesh_ids       = (int *) malloc(sizeof(int));

  return id;
}


void
PDM_dmesh_nodal_to_dmesh_free
(
const int hdl
)
{
  _PDM_dmesh_nodal_to_dmesh* dmn_to_dm = _get_from_id(hdl);

  if(( dmn_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( dmn_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !dmn_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < dmn_to_dm->n_mesh; ++i_mesh) {
      PDM_dmesh_free(dmn_to_dm->dmesh_ids[i_mesh]);
    }
  }

  free(dmn_to_dm->dmesh_nodal_ids);
  free(dmn_to_dm->dmesh_ids);
}

