/*
 * \file
 */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_size_idx_from_stride.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/**
 *
 * \brief Return the integer size required to store an idx array from stride array
 *
 * \param [in]   stride    Stride array
 * \param [in]   s_sride   Size of the stride array
 * \param [in]   comm      MPI communicator
 *
 * \return the size of int required
 *
 */

int
PDM_size_idx_from_stride
(
 int          *stride,
 int          s_stride,
 PDM_MPI_Comm comm
)
{

	size_t sum = 0;		

  INT_MAX;
  int coeff = 4;
  int large = 0;
	for (int i = 0; i < s_stride; i++) {
		sum += (size_t) stride[i];
		if (sum >= (size_t) (INT_MAX/coeff)) {
			large = 1;
			break;
		}
	}

  int s_large = 0;
  PDM_MPI_Allreduce (&large, &s_large, 1, PDM_MPI_INT,  PDM_MPI_MAX, comm);

  if (s_large) {
  	return (int) sizeof(PDM_g_num_t);
  }
  else {
  	return (int) sizeof(int);
  }

}


#ifdef __cplusplus
}
#endif /* __cplusplus */

