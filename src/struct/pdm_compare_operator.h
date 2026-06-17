#ifndef __PDM_COMPARE_OPERATOR_H__
#define __PDM_COMPARE_OPERATOR_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
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
 * \brief Compare two unsigned ordered nuplets of integers
 *
 * (Suboptimal for large nuplets)
 *
 * \param [in] size1    Size of first nuplet
 * \param [in] nuplet1  First nuplet (size = \p size1)
 * \param [in] size2    Size of second nuplet
 * \param [in] nuplet2  Second nuplet (size = \p size2)
 *
 * \return  1 if the nuplets are similar (up to a cyclic permutation)
 *         -1 if the nuplets are reversed (up to a cyclic permutation)
 *          0 otherwise
 */
int
PDM_compare_unsigned_ordered_nuplets_int
(
  const int size1,
  const int nuplet1[],
  const int size2,
  const int nuplet2[]
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /*  __PDM_SORT_H__ */
