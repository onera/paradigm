/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_sort.h"
#include "pdm_logging.h"

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

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_long
(
PDM_g_num_t *a,
PDM_g_num_t *b
)
{
  PDM_g_num_t tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_int
(
int *a,
int *b
)
{
  int tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

/**
 *
 * \brief Swap element values
 *
 * \param [inout]   a       Pointer on first value
 * \param [inout]   b       Pointer on second value
 *
 */

static inline void
_swap_double
(
double *a,
double *b
)
{
  double tmp = *a;
  *a = *b;
  *b = tmp;
  return;
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
PDM_sort_long
(
 PDM_g_num_t *array,
 int         *order,
 int          lArray
)
{
  /* size of subarray sorted by straight insertion */
  const int M = 7;
  /* default size of the stack */
  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  PDM_g_num_t a;
  int        b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {

    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else{
        int k = (l+ir) / 2;
        _swap_long (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_long (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_long (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_long (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_long (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else{
        int k = (l+ir) / 2;
        _swap_long (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_long (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_long (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_long (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_long (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;

}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout]   array        Array to sort
 * \param [inout]   order        new indice to old indice (or NULL)
 * \param [in]      lArray       Array length
 *
 */

void
PDM_sort_int
(
 int         *array,
 int         *order,
 int          lArray
)
{
  /* size of subarray sorted by straight insertion */
  const int M = 7;
  /* default size of the stack */
  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  int  a;
  int  b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_int (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_int (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_int (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_int (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_int (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_int (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_int (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_int (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_int (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_int (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;

}


/**
 *
 * \brief Search element index in a sorted array
 *
 * \param [inout]   array        Array to sort
 * \param [inout]   order        new indice to old indice (or NULL)
 * \param [in]      lArray       Array length
 *
 */

void
PDM_sort_double
(
 double    *array,
 int       *order,
 int        lArray
)
{
  /* size of subarray sorted by straight insertion */

  const int M = 7;
  /* default size of the stack */

  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  double  a;
  int  b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  if (order != NULL) {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          b = order[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
            order[i+1] = order[i];
          }
          array[i+1] = a;
          order[i+1] = b;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_double (&(array[k]), &(array[l+1]));
        _swap_int (&(order[k]), &(order[l+1]));
        if (array[l] > array[ir]){
          _swap_double (&(array[l]), &(array[ir]));
          _swap_int (&(order[l]), &(order[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_double (&(array[l+1]), &(array[ir]));
          _swap_int (&(order[l+1]), &(order[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_double (&(array[l]), &(array[l+1]));
          _swap_int (&(order[l]), &(order[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        b = order[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_double (&(array[i]), &(array[j]));
          _swap_int (&(order[i]), &(order[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        order[l+1] = order[j];
        order[j] = b;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  else {
    for (;;) {
      if ( ir-l < M ) {
        for (j = l+1; j <= ir; ++j) {
          a = array[j];
          for(i = j-1; i >= l; --i) {
            if (array[i] <= a) {
              break;
            }
            array[i+1] = array[i];
          }
          array[i+1] = a;
        }
        if (jstack < 0) {
          break;
        }
        ir = istack[jstack--];
        l  = istack[jstack--];
      }
      else {
        int k = (l+ir) / 2;
        _swap_double (&(array[k]), &(array[l+1]));
        if (array[l] > array[ir]){
          _swap_double (&(array[l]), &(array[ir]));
        }
        if (array[l+1] > array[ir]) {
          _swap_double (&(array[l+1]), &(array[ir]));
        }
        if (array[l] > array[l+1]) {
          _swap_double (&(array[l]), &(array[l+1]));
        }
        i = l + 1;
        j = ir;
        a = array[l+1];
        for (;;) {
          do {
            ++i;
          } while (array[i] < a);
          do {
            --j;
          } while (array[j] > a);

          if (j < i) {
            break;
          }
          _swap_double (&(array[i]), &(array[j]));
        }

        array[l+1] = array[j];
        array[j] = a;
        jstack += 2;

        if (jstack >= sizeStack) {
          sizeStack *= 2;
          istack = (int *) realloc (istack, sizeof(int) * sizeStack);
        }

        if (ir-i+1 >= j-1) {
          istack[jstack  ] = ir;
          istack[jstack-1] = i;
          ir = j -1;
        }
        else {
          istack[jstack  ] = j -1;
          istack[jstack-1] = l;
          l = i;
        }
      }
    }
  }
  free (istack);
  return;
}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
void
PDM_quick_sort_long
(
 PDM_g_num_t a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    PDM_g_num_t t;
    PDM_g_num_t pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    PDM_quick_sort_long(a, l  , j-1);
    PDM_quick_sort_long(a, j+1,   r);
  }
}

/**
 *
 * \brief Unique
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
int
PDM_inpace_unique_long
(
 PDM_g_num_t a[],
 int l,
 int r
)
{
  // PDM_quick_sort_long(a, l, r); /* Less optimal than PDM_sort_long */
  int array_size = r - l;
  // printf("PDM_inpace_unique_long::array_size::%d\n", array_size);
  PDM_sort_long(&a[l], NULL, array_size);

  int new_size  = 1;
  int idx_write = l;
  PDM_g_num_t last_value = a[l];
  a[idx_write++] = last_value;
  for (int idx = l+1; idx <= r; idx++) {
    if(last_value != a[idx]){
      last_value = a[idx];
      a[idx_write++] = a[idx];
      new_size++;
    }
  }

  return new_size;
}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */
void
PDM_quick_sort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    PDM_quick_sort_int(a, l  , j-1);
    PDM_quick_sort_int(a, j+1,   r);
  }
}

/**
 *
 * \brief Quick sort 2
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 *
 */
void
PDM_quick_sort_int2
(
 int          a[],
 int          l,
 int          r,
 int          c[]
)
{
  if (l < r) {
    int j = r+1;
    int  t, v;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

      v    = c[i];
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    v    = c[l];
    c[l] = c[j];
    c[j] = v;

    PDM_quick_sort_int2(a, l  , j-1, c);
    PDM_quick_sort_int2(a, j+1,   r, c);
  }
}

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 * \param [inout]   c     Array sorted as a
 *
 */

void
PDM_quick_sort_long2
(
 PDM_g_num_t a[],
 int          l,
 int          r,
 int          c[]
)
{
  if (l < r) {
    int j = r+1;

    PDM_g_num_t t;
    int v;
    PDM_g_num_t pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

      v    = c[i];
      c[i] = c[j];
      c[j] = v;
    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    v    = c[l];
    c[l] = c[j];
    c[j] = v;

    PDM_quick_sort_long2(a, l  , j-1, c);
    PDM_quick_sort_long2(a, j+1,   r, c);
  }
}


/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_string
(
const void* a,
const void* b,
      void* ctxt
)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  // log_trace("PDM_operator_compare_string::  %d %d - idx -> %d %d \n",  i, j, us->idx[i], us->idx[j]);

  int ni = us->idx[i+1] - us->idx[i];
  int nj = us->idx[j+1] - us->idx[j];

  // log_trace("PDM_operator_compare_string:: %d %d - %d %d \n", i, j, ni, nj);

  char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];

  char* carr_i = malloc( sizeof(char) * (ni+1));
  char* carr_j = malloc( sizeof(char) * (nj+1));
  for(int k = 0; k < ni; ++k){
    carr_i[k] = (char)arr_i[k];
  }
  for(int k = 0; k < nj; ++k){
    carr_j[k] = (char)arr_j[k];
  }
  carr_i[ni] = '\0';
  carr_j[nj] = '\0';
  int i_comp = strcmp(carr_i, carr_j);
  free(carr_i);
  free(carr_j);
  if(i_comp >= 0){
    return 0;
  } else {
    return 1;
  }

  // Implementation different
  // if(ni < nj){
  //   return 0;
  // } else {
  //   char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  //   char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];
  //   int i_comp = strncmp(arr_i, arr_j, ni);
  //   if(i_comp > 0) {
  //     return 0;
  //   } else {
  //     return 1;
  //   }
  // }
}


/**
 *
 * \brief Equal operator for connectivities
 *
 */
int
PDM_operator_equal_string
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  int ni = us->idx[i+1] - us->idx[i];
  int nj = us->idx[j+1] - us->idx[j];

  // log_trace("PDM_operator_equal_string:: %d %d - %d %d - %d %d \n", i, j, ni, nj, us->idx[i], us->idx[j]);

  char* arr_i = (char *) &us->arr[us->idx[i]*sizeof(char)];
  char* arr_j = (char *) &us->arr[us->idx[j]*sizeof(char)];

  if(ni != nj){
    return 0;
  } else {
    int i_comp = strncmp(arr_i, arr_j, ni);
    if(i_comp == 0) {
      return 1;
    } else {
      return 0;
    }
  }
}

/**
 *
 * \brief Compare operator for connectivities
 *
 */
int
PDM_operator_compare_connectivity
(
const void* a,
const void* b,
      void* ctxt
)
{
  int i = *(const int *) a;
  int j = *(const int *) b;

  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;

  // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
  // log_trace("PDM_operator_compare_connectivity:: %d %d with key %lu %lu \n", i, j, us->key[i], us->key[j]);

  if( us->key[i] == us->key[j]){

    // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
    int ni = us->idx[i+1] - us->idx[i];
    int nj = us->idx[j+1] - us->idx[j];

    if(ni == nj){

      int* arr_i = (int *) &us->arr[us->idx[i]*sizeof(int)];
      int* arr_j = (int *) &us->arr[us->idx[j]*sizeof(int)];

      /* Dans notre cas on veut sort les entiers avant de les comparers */
      int* sort_arr_i = (int*) malloc( ni * sizeof(int));
      int* sort_arr_j = (int*) malloc( ni * sizeof(int));

      for(int k = 0; k < ni; ++k){
        sort_arr_i[k] = arr_i[k];
        sort_arr_j[k] = arr_j[k];
        // log_trace("PDM_operator_compare_connectivity:: %d %d \n", i, j);
      }
      PDM_quick_sort_int(sort_arr_i, 0, ni-1);
      PDM_quick_sort_int(sort_arr_j, 0, ni-1);

      // log_trace("Comparison of ::");
      // for(int k = 0; k < ni; ++k){
      //   log_trace(" %d %d ", sort_arr_i[k], sort_arr_j[k]);
      // }
      // log_trace("\n");

      for(int k = 0; k < ni; ++k){
        if(sort_arr_i[k] < sort_arr_j[k]) {
          free(sort_arr_i);
          free(sort_arr_j);
          return 1;
        } else if( sort_arr_i[k] > sort_arr_j[k] ) {
          free(sort_arr_i);
          free(sort_arr_j);
          return 0;
        }
      }
      free(sort_arr_i);
      free(sort_arr_j);
    } else {
      return ni < nj;
    }
  } else {
    return us->key[i] < us->key[j];
  }
  return 0;
}


/**
 *
 * \brief Equal operator for connectivities
 *
 */
int
PDM_operator_equal_connectivity
(
const void* a,
const void* b,
      void* ctxt)
{
  int i = *(const int *) a;
  int j = *(const int *) b;


  PDM_user_defined_sort* us = (PDM_user_defined_sort*) ctxt;
  // printf("PDM_operator_equal_connectivity:: %d %d - %d %d - %d %d \n", i, j, ni, nj, us->idx[i], us->idx[j]);

  if( us->key[i] == us->key[j]){
    int ni = us->idx[i+1] - us->idx[i];
    int nj = us->idx[j+1] - us->idx[j];

    if(ni == nj){

      int* arr_i = (int*) &us->arr[us->idx[i]*sizeof(int)];
      int* arr_j = (int*) &us->arr[us->idx[j]*sizeof(int)];

      /* Dans notre cas on veut sort les entiers avant de les comparers */
      int* sort_arr_i = (int*) malloc( ni * sizeof(int));
      int* sort_arr_j = (int*) malloc( ni * sizeof(int));

      for(int k = 0; k < ni; ++k){
        sort_arr_i[k] = arr_i[k];
        sort_arr_j[k] = arr_j[k];
      }
      PDM_quick_sort_int(sort_arr_i, 0, ni-1);
      PDM_quick_sort_int(sort_arr_j, 0, ni-1);

      for(int k = 0; k < ni; ++k){
        // printf(" \t sort_arr_i[%d] = %d | sort_arr_j[%d] = %d \n", k, sort_arr_i[k], k, sort_arr_j[k]);
        if(sort_arr_i[k] != sort_arr_j[k]) {
          free(sort_arr_i);
          free(sort_arr_j);
          return 0;
        }
      }

      free(sort_arr_i);
      free(sort_arr_j);
    } else {
      return 0;
    }
  } else {
    return 0;
  }

  return 1;
}

/**
 *
 * \brief Indirect sort of structure according to a user comparison function
 *
 * \param [inout] array        Array to sort
 * \param [in]    lArray       Array length
 * \param [in]    lArray       User compare function (return int )
 * \param [in]    lArray       Context anonymous pointer to a struct to perfomr user compare
 *
 */
void
PDM_sort_int_special
(
 int                  *array,
 int                   lArray,
 pdm_operator_compare  comp,
 void*                 context
)
{
  /* size of subarray sorted by straight insertion */
  const int M = 7;
  /* default size of the stack */
  int sizeStack = 64; /* default size of the stack */
  int jstack = -1;
  int l = 0;
  int i;
  int j;
  int ir = lArray - 1;
  PDM_g_num_t a;
  // int        b;
  int *istack = (int *) malloc (sizeof(int) * sizeStack);

  for (;;) {
    if ( ir-l < M ) {
      for (j = l+1; j <= ir; ++j) {
        a = array[j];
        for(i = j-1; i >= l; --i) {
          if (comp(&array[i], &a       , context)) {
            break;
          }
          array[i+1] = array[i];
        }
        array[i+1] = a;
      }
      if (jstack < 0) {
        break;
      }
      ir = istack[jstack--];
      l  = istack[jstack--];
    }
    else{
      int k = (l+ir) / 2;
      _swap_long (&(array[k]), &(array[l+1]));
      if ( comp(&array[ir], &array[l], context) ){
        _swap_long (&(array[l]), &(array[ir]));
      }
      if ( comp(&array[ir], &array[l+1], context) ) {
        _swap_long (&(array[l+1]), &(array[ir]));
      }
      if ( comp(&array[l+1], &array[l], context) ) {
        _swap_long (&(array[l]), &(array[l+1]));
      }
      i = l + 1;
      j = ir;
      a = array[l+1];
      for (;;) {
        do {
          ++i;
        } while (comp(&array[i], &a       , context));
        do {
          --j;
        } while (comp(&a       , &array[j], context));

        if (j < i) {
          break;
        }
        _swap_long (&(array[i]), &(array[j]));
      }

      array[l+1] = array[j];
      array[j] = a;
      jstack += 2;

      if (jstack >= sizeStack) {
        sizeStack *= 2;
        istack = (int *) realloc (istack, sizeof(int) * sizeStack);
      }

      if (ir-i+1 >= j-1) {
        istack[jstack  ] = ir;
        istack[jstack-1] = i;
        ir = j -1;
      }
      else {
        istack[jstack  ] = j -1;
        istack[jstack-1] = l;
        l = i;
      }
    }
  }
  free (istack);
  return;

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
