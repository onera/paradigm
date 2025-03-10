
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_quick_sort.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"


TEST_CASE("[pdm_quick_sort] - PDM_quick_sort_long") {

  int n = 8;
  PDM_g_num_t array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  PDM_quick_sort_long(array, 0, n-1);

  PDM_g_num_t expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};

  if(0 == 1) {
    PDM_log_trace_array_long(array, n, "array :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);

}


TEST_CASE("[pdm_quick_sort] - PDM_quick_sort_int") {

  int n = 8;
  int array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  PDM_quick_sort_int(array, 0, n-1);

  int expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};

  if(0 == 1) {
    PDM_log_trace_array_int(array, n, "array :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);

}


TEST_CASE("[pdm_quick_sort] - PDM_quick_sort_long2") {

  int n = 8;
  PDM_g_num_t array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  int *order = NULL;
  PDM_malloc(order, n, int);

  for(int i = 0; i < n; ++i) {
    order[i] = i;
  }

  PDM_quick_sort_long2(array, 0, n-1, order);

  PDM_g_num_t expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};
  int         expected_order[8] = {4, 2, 1, 6, 3, 0,  5,  7};

  if(0 == 1) {
    PDM_log_trace_array_long(array, n, "array :: ");
    PDM_log_trace_array_int (order, n, "order :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);
  CHECK_EQ_C_ARRAY(order, expected_order, n);

  PDM_free(order);
}


TEST_CASE("[pdm_quick_sort] - PDM_quick_sort_int2") {

  int n = 8;
  int array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  int *order = NULL;
  PDM_malloc(order, n, int);

  for(int i = 0; i < n; ++i) {
    order[i] = i;
  }

  PDM_quick_sort_int2(array, 0, n-1, order);

  int expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};
  int expected_order[8] = {4, 2, 1, 6, 3, 0,  5,  7};

  if(0 == 1) {
    PDM_log_trace_array_int(array, n, "array :: ");
    PDM_log_trace_array_int (order, n, "order :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);
  CHECK_EQ_C_ARRAY(order, expected_order, n);
  PDM_free(order);

}
