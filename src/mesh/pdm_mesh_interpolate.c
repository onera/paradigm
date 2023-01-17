
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mesh_interpolate_priv.h"
#include "pdm_mesh_interpolate.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/


/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a structure that compute a global mean
 *
 * \param [in]   n_part       Number of local partitions
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to \ref PDM_mesh_interpolate object
 */
PDM_mesh_interpolate_t*
PDM_mesh_interpolate_create
(
 const int            n_domain,
 const int           *n_part,
 const int           *n_group,
 const int            interp_kind,
 const PDM_MPI_Comm   comm
)
{
  PDM_mesh_interpolate_t* mi = malloc(sizeof(PDM_mesh_interpolate_t));

  mi->n_domain    = n_domain;
  mi->n_part      = malloc( n_domain * sizeof(int)); // Make a copy to avoid pb in cython
  for(int i = 0; i < mi->n_domain; ++i) {
    mi->n_part[i] = n_part[i];
  }
  mi->comm        = comm;

  mi->n_part_idx    = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->n_part_g_idx  = (int * ) malloc( (n_domain + 1) * sizeof(int));
  mi->parts = malloc(n_domain * sizeof(_part_t *));


  mi->n_part_idx  [0] = 0;
  mi->n_part_g_idx[0] = 0;
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    mi->parts[i_domain] = malloc( n_part[i_domain] * sizeof(_part_t));
    mi->n_part_idx[i_domain+1] = mi->n_part_idx[i_domain] + mi->n_part[i_domain];

    int n_part_l = n_part[i_domain];
    int n_part_g = -100;
    PDM_MPI_Allreduce(&n_part_l, &n_part_g, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    mi->n_part_g_idx[i_domain+1] = mi->n_part_idx[i_domain] + n_part_g;

  }



  /* Compute graph comm from gnum */
  /* Deduce graph with all graphe inside same domain and between domain */
  /* Compute weight */
  /* Create protocol only between join - Distant neighbor */

  return mi;
}


void
PDM_mesh_interpolate_part_domain_interface_shared_set
(
  PDM_mesh_interpolate_t      *mi,
  PDM_part_domain_interface_t *pdi
)
{
  mi->pdi = pdi;
}


void
PDM_mesh_interpolate_free
(
 PDM_mesh_interpolate_t *mi
)
{
  for(int i_domain = 0; i_domain < mi->n_domain; ++i_domain) {
    free(mi->parts[i_domain]);
  }
  free(mi->parts);
  free(mi->n_part_idx);
  free(mi->n_part_g_idx);
  free(mi->n_part);


  free(mi);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
