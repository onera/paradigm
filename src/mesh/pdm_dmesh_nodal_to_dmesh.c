
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
#include "pdm_dmesh_nodal_to_dmesh_priv.h"
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


/*============================================================================
 * Global variable
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/

static
int
_generate_faces_from_dmesh_nodal
(
  int dmesh_nodal_id
)
{
  PDM_UNUSED(dmesh_nodal_id);

  return -1;
}



static
int
_generate_edges_from_dmesh_nodal
(
  int dmesh_nodal_id
)
{
  PDM_UNUSED(dmesh_nodal_id);

  return -1;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_dmesh_nodal_to_dmesh_t*
PDM_dmesh_nodal_to_dmesh_create
(
const int             n_mesh,
const PDM_MPI_Comm    comm,
const PDM_ownership_t owner
)
{
  _pdm_dmesh_nodal_to_dmesh_t *dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) malloc(sizeof(_pdm_dmesh_nodal_to_dmesh_t));

  dmn_to_dm->comm              = comm;
  dmn_to_dm->owner             = owner;
  dmn_to_dm->results_is_getted = PDM_FALSE;
  dmn_to_dm->n_mesh            = n_mesh;

  dmn_to_dm->dmesh_nodal_ids = (int                 *) malloc(sizeof(int                 ));
  dmn_to_dm->dmesh_ids       = (int                 *) malloc(sizeof(int                 ));
  dmn_to_dm->dmesh_nodal     = (_pdm_dmesh_nodal_t **) malloc(sizeof(_pdm_dmesh_nodal_t *));
  dmn_to_dm->dmesh           = (_pdm_dmesh_t       **) malloc(sizeof(_pdm_dmesh_t       *));

  return (PDM_dmesh_nodal_to_dmesh_t *) dmn_to_dm;
}

/**
 * \brief  Add dmesh_nodal
 */
void
PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal
(
        PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
  const int                         i_mesh,
        PDM_DMesh_nodal_t          *dmn
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;
  // _dmn_to_dm->dmesh_nodal_ids[i_mesh] = dmesh_nodal_id;
  _dmn_to_dm->dmesh_nodal[i_mesh] = (_pdm_dmesh_nodal_t*) dmn;
}

/**
 * \brief  Free
 */
void
PDM_dmesh_nodal_to_dmesh_free
(
  PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;

  if(( _dmn_to_dm->owner == PDM_OWNERSHIP_KEEP ) ||
     ( _dmn_to_dm->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !_dmn_to_dm->results_is_getted)){

    for(int i_mesh = 0; i_mesh < _dmn_to_dm->n_mesh; ++i_mesh) {
      PDM_dmesh_free( (PDM_dmesh_t*)_dmn_to_dm->dmesh[i_mesh]);
    }
  }

  free(_dmn_to_dm->dmesh_nodal);
  free(_dmn_to_dm->dmesh);
  free(_dmn_to_dm->dmesh_nodal_ids);
  free(_dmn_to_dm->dmesh_ids);
}



void
PDM_dmesh_nodal_to_dmesh_compute
(
  PDM_dmesh_nodal_to_dmesh_t*                dmn_to_dm,
  const PDM_dmesh_nodal_to_dmesh_transform_t transform_kind
)
{
  _pdm_dmesh_nodal_to_dmesh_t* _dmn_to_dm = (_pdm_dmesh_nodal_to_dmesh_t *) dmn_to_dm;

  for(int i_mesh = 0; i_mesh < _dmn_to_dm->n_mesh; ++i_mesh) {

    switch (transform_kind) {

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE:
        {
          int dmesh_id = _generate_faces_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal_ids[i_mesh]);
          _dmn_to_dm->dmesh_ids[i_mesh] = dmesh_id;
        }
        break;

      case PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE:
        {
          int dmesh_id = _generate_edges_from_dmesh_nodal(_dmn_to_dm->dmesh_nodal_ids[i_mesh]);
          _dmn_to_dm->dmesh_ids[i_mesh] = dmesh_id;
        }
        break;
    }
  }

  // Boundary management

  // Join management


}
