/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_compare_operator.h"
#include "pdm_mem_tool.h"
#include "pdm_quick_sort.h"
#include "pdm_sort.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

int
PDM_compare_unsigned_ordered_nuplets_int
(
  const int size1,
  const int nuplet1[],
  const int size2,
  const int nuplet2[]
)
{
  if (size1 != size2) {
    return 0;
  }

  if (size1 == 2) {
    // Handle pairs separately
    if (nuplet2[0] == nuplet1[0] && nuplet2[1] == nuplet1[1]) {
      return 1;
    }
    else if (nuplet2[0] == nuplet1[1] && nuplet2[1] == nuplet1[0]) {
      return -1;
    }
    else {
      return 0;
    }
  }

  assert(size1 >= 3);

  for (int i = 0; i < size1; i++) {
    if (nuplet2[i] == nuplet1[0]) { // we found a common element ("pivot")
      int direct  = 1;
      int reverse = 1;
      for (int j = 1; j < size1; j++) {
        // Starting from pivot, traverse nuplet 2 in both directions
        // and compare to nuplet 1 element-wise
        if (nuplet2[(i+j)%size1] != nuplet1[j]) {
          direct = 0;
        }
        if (nuplet2[(i+size1-j)%size1] != nuplet1[j]) {
          reverse = 0;
        }
      }

      if (direct && !reverse) {
        return 1;
      } else if (reverse) {
        return -1;
      } else {
        return 0;
      }
    }
  }

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
