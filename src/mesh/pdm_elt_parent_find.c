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
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_elt_parent_find.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 * \brief Find parent in a set of elements 
 *
 * \param [in]     elt_distrib          Distribution of elements on processes 
 * \param [in]     elt_def_idx          Element definition index 
 * \param [in]     elt_def              Element definition
 * \param [in]     n_elt_to_find        Number of elements to find
 * \param [in]     elt_to_find_def_idx  Element to find definition index 
 * \param [in]     elt_to_find_def      Element to find definition
 * \param [in]     comm                 MPI Communicator
 * \param [inout]  parent               Parent element of found element, 0 otherwise

 */

void
PDM_elt_parent_find
(
 const int         *elt_distrib,
 const int         *elt_def_idx,
 const PDM_g_num_t *elt_def,
 const int          n_elt_to_find,
 const int         *elt_to_find_def_idx,
 const PDM_g_num_t *elt_to_find_def,
 const PDM_MPI_Comm comm,     
 const int         *parent
)
{

  /* Build distributed hash tab (key = sum of integer used for element definition */

  /* Find parent in distributed hash table */

  /* Return result in the same order of elt_to_find_def array */
  
  /* Push result into parent array */
  
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
