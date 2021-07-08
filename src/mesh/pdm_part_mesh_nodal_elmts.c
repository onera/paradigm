/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

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
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_geom_elem.h"
#include "pdm_array.h"

#include "pdm_writer.h"

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

/*----------------------------------------------------------------------------
 * Maximum number of blocks depending of block type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
 const int          mesh_dimension,
 const int          n_part,
 const PDM_MPI_Comm comm
)
{
  PDM_part_mesh_nodal_elmts_t *pmne = (PDM_part_mesh_nodal_elmts_t *) malloc (sizeof(PDM_part_mesh_nodal_elmts_t));

  pmne->comm           = comm;
  pmne->mesh_dimension = mesh_dimension;
  pmne->n_part         = n_part;


  pmne->n_sections      = 0;
  pmne->n_elmts         = (int * ) malloc( n_part * sizeof(int));
  pmne->sections_id     = NULL;
  pmne->sections_std    = NULL;
  pmne->sections_poly2d = NULL;
  pmne->sections_poly3d = NULL;

  return pmne;
}


/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */


void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t* pmne
)
{

  if(pmne->n_elmts != NULL) {
    free(pmne->n_elmts);
  }

  /*
   *  Free all sections
   */



  if(pmne->sections_id != NULL) {
    free(pmne->sections_id);
  }

  free(pmne);
}
