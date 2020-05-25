/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_radix_sort.h"

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

static inline
PDM_g_num_t
_get_max_value_in_array_long
(
 PDM_g_num_t* v,
 int          size
)
{
  PDM_g_num_t max = v[0];
  for(int i = 1; i < size; ++i){
    if( v[i] > max){
      max = v[i];
    }
  }
  return max;
}


static inline
int*
_counting_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  /* First step - Count */
  int n_buckets = 10;
  int* count = malloc( (n_buckets + 1) * sizeof(int));

  /* Set to zero */
  for(int i = 0; i < n_buckets+1; ++i){
    count[i] = 0;
  }

  /* Count */
  for(int i = beg; i < end; ++i) {
    count[(v[i]/place)%10]++;
  }

  /* Rebuild array of displacement */
  for(int i = 1; i < n_buckets+1; ++i) {
    count[i] += count[i-1];
  }

  for(int i = end-1; i >= beg; i--){
    tmp[beg+count[(v[i]/place)%10]-1] = v[i];
    count[(v[i]/place)%10]--;
  }

  for(int i = beg; i < end; ++i) {
    v[i] = tmp[i];
  }

  return count;
}


static inline
void
_cc_radix_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  int place_init = place*10;
  place = 1;
  while(place < place_init){
    int* range = _counting_sort(v, tmp, beg, end, place);
    free(range);
    place *= 10;
  }

}



static inline
void
_std_radix_sort
(
 PDM_g_num_t* v,
 PDM_g_num_t* tmp,
 int beg,
 int end,
 int place
)
{
  int cache_size     = 5000;

  if(beg == end){
    return;
  }

  int* range = _counting_sort(v, tmp, beg, end, place);

  place /= 10;

  if(place == 0) {
    free(range);
    return;
  }

  /* Il faut aglomerer les ranges succesives pour gagner si elle sont trop petite */
  // int i  = 0;
  // int i+1 = 1;

  // fmt::print(" Range is : {0} \n", range);

  for(int i = 0; i < 11; ++i){
    // fmt::print("Manage {0} = {1} / {2}\n", i, range[i], range[i+1]);
    if( (range[i+1] - range[i]) == 0){
    } else if( (range[i+1] - range[i]) > cache_size){
      // fmt::print("Case 1 :: {0} - {1} \n ", range[i], range[i+1]);
      _std_radix_sort(v, tmp, beg+range[i], beg+range[i+1] , place);
    } else if( (range[i+1] - range[i]) <= cache_size) {
      // fmt::print("Case 3 :: {0} - {1} \n ", range[i], range[i+1]);
      _cc_radix_sort(v, tmp, beg+range[i], beg+range[i+1], place);
    }
  }

  free(range);

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout] array        Array to sort
 * \param [inout] order        new indice to old indice (or NULL)
 * \param [in]    lArray       Array length
 *
 */
void
PDM_radix_sort_long
(
 PDM_g_num_t *array,
 int         *order,
 int          lArray
)
{

  PDM_g_num_t max = _get_max_value_in_array_long(array, lArray);

  int n_step = -1;
  int lar = max;
  while(lar > 0){
    n_step++;
    lar /= 10;
  }

  // Il faut trouver le max dans la base pour reprensenter le min et le max
  //  --> si beaucoup d'écart on tente la moyenne
  int* tmp = (int *) malloc( lArray * sizeof(int));

  int place = (int) pow(10, n_step);
  std_radix_sort(array, tmp, 0, lArray, place);


  free(tmp);

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
