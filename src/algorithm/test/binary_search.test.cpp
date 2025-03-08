
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_binary_search.h"

TEST_CASE("[pdm_binary_search] - PDM_binary_search_gap_long") {

  int         n_rank = 3;
  PDM_g_num_t distrib[4] = {0, 3, 30, 300};

  PDM_g_num_t gnum1 = 12;
  int pos1 = PDM_binary_search_gap_long(gnum1, distrib, n_rank+1);

  CHECK(pos1 == 1);

  PDM_g_num_t gnum2 = -1;
  int pos2 = PDM_binary_search_gap_long(gnum2, distrib, n_rank+1);

  CHECK(pos2 == -1);

  PDM_g_num_t gnum3 = 300;
  int pos3 = PDM_binary_search_gap_long(gnum3, distrib, n_rank+1);

  CHECK(pos3 == -1);
}


TEST_CASE("[pdm_binary_search] - PDM_binary_search_gap_size_t") {

  int    n_rank = 3;
  size_t distrib[4] = {0, 3, 30, 300};

  size_t gnum1 = 12;
  int pos1 = PDM_binary_search_gap_size_t(gnum1, distrib, n_rank+1);

  CHECK(pos1 == 1);

  size_t gnum2 = -1;
  int pos2 = PDM_binary_search_gap_size_t(gnum2, distrib, n_rank+1);

  CHECK(pos2 == -1);

  size_t gnum3 = 300;
  int pos3 = PDM_binary_search_gap_size_t(gnum3, distrib, n_rank+1);

  CHECK(pos3 == -1);
}



TEST_CASE("[pdm_binary_search] - PDM_binary_search_gap_int") {

  int n_rank = 3;
  int distrib[4] = {0, 3, 30, 300};

  int gnum1 = 12;
  int pos1 = PDM_binary_search_gap_int(gnum1, distrib, n_rank+1);

  CHECK(pos1 == 1);

  int gnum2 = -1;
  int pos2 = PDM_binary_search_gap_int(gnum2, distrib, n_rank+1);

  CHECK(pos2 == -1);

  int gnum3 = 300;
  int pos3 = PDM_binary_search_gap_int(gnum3, distrib, n_rank+1);

  CHECK(pos3 == -1);
}




TEST_CASE("[pdm_binary_search] - PDM_binary_search_long") {

  int n_elmt = 6;
  PDM_g_num_t array[6] = {-8, 0, 4, 8, 9, 10};

  PDM_g_num_t gnum1 = 4;
  int pos1 = PDM_binary_search_long(gnum1, array, n_elmt);

  CHECK(pos1 == 2);

  PDM_g_num_t gnum2 = 2;
  int pos2 = PDM_binary_search_long(gnum2, array, n_elmt);

  CHECK(pos2 == -1);

  PDM_g_num_t gnum3 = -8;
  int pos3 = PDM_binary_search_long(gnum3, array, n_elmt);

  CHECK(pos3 == 0);

}


TEST_CASE("[pdm_binary_search] - PDM_binary_search_int") {

  int n_elmt = 6;
  int array[6] = {-8, 0, 4, 8, 9, 10};

  int gnum1 = 4;
  int pos1 = PDM_binary_search_int(gnum1, array, n_elmt);

  CHECK(pos1 == 2);

  int gnum2 = 2;
  int pos2 = PDM_binary_search_int(gnum2, array, n_elmt);

  CHECK(pos2 == -1);

  int gnum3 = -8;
  int pos3 = PDM_binary_search_int(gnum3, array, n_elmt);

  CHECK(pos3 == 0);

}


TEST_CASE("[pdm_binary_search] - PDM_binary_search_gap_double") {

  int n_rank = 3;
  double distrib[4] = {0., 3., 30., 300};

  double gnum1 = 12.;
  int pos1 = PDM_binary_search_gap_double(gnum1, distrib, n_rank+1);

  CHECK(pos1 == 1);

  double gnum2 = -1.;
  int pos2 = PDM_binary_search_gap_double(gnum2, distrib, n_rank+1);

  CHECK(pos2 == -1);

  double gnum3 = 300.;
  int pos3 = PDM_binary_search_gap_double(gnum3, distrib, n_rank+1);

  CHECK(pos3 == -1);

}


TEST_CASE("pdm_binary_search] - PDM_search_rank") {

  int n_rank = 3;
  PDM_g_num_t distrib[4] = {0, 3, 30, 300};

  PDM_g_num_t gnum1 = 12;
  int pos1 = PDM_search_rank(gnum1, distrib, 0, n_rank);

  CHECK(pos1 == 1);

}
