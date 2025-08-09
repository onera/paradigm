
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_sort.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"
#include "pdm_priv.h"


TEST_CASE("[pdm_sort] - PDM_sort_long") {

  int n = 8;
  PDM_g_num_t array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  int *order = NULL;
  PDM_malloc(order, n, int);

  for(int i = 0; i < n; ++i) {
    order[i] = i;
  }

  PDM_sort_long(array, order, n);

  PDM_g_num_t expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};
  int         expected_order[8] = {4, 2, 1, 6, 0, 3,  5,  7};

  if(0 == 1) {
    PDM_log_trace_array_long(array, n, "array :: ");
    PDM_log_trace_array_int (order, n, "order :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);
  CHECK_EQ_C_ARRAY(order, expected_order, n);

  PDM_free(order);
}

TEST_CASE("[pdm_sort] - PDM_sort_int") {

  int n = 8;
  int array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  int *order = NULL;
  PDM_malloc(order, n, int);

  for(int i = 0; i < n; ++i) {
    order[i] = i;
  }

  PDM_sort_int(array, order, n);

  int expected_array[8] = {1, 3, 4, 7, 8, 8, 12, 40};
  int expected_order[8] = {4, 2, 1, 6, 0, 3,  5,  7};

  if(0 == 1) {
    PDM_log_trace_array_int(array, n, "array :: ");
    PDM_log_trace_array_int(order, n, "order :: ");
  }

  CHECK_EQ_C_ARRAY(array, expected_array, n);
  CHECK_EQ_C_ARRAY(order, expected_order, n);

  PDM_free(order);
}


TEST_CASE("[pdm_sort] - PDM_sort_double") {

  int n = 8;
  double array[8] = {8, 4, 3, 8, 1, 12, 7, 40};

  int *order = NULL;
  PDM_malloc(order, n, int);

  for(int i = 0; i < n; ++i) {
    order[i] = i;
  }

  PDM_sort_double(array, order, n);

  double expected_array[8] = {1., 3., 4., 7., 8., 8., 12., 40.};
  int    expected_order[8] = {4, 2, 1, 6, 0, 3,  5,  7};

  if(0 == 1) {
    PDM_log_trace_array_double(array, n, "array :: ");
    PDM_log_trace_array_int   (order, n, "order :: ");
  }

  CHECK_EQ_C_ARRAY_FLOAT(array, expected_array, n, 1.e-6);
  CHECK_EQ_C_ARRAY(order, expected_order, n);

  PDM_free(order);
}
