#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_compare_operator.h"
#include "pdm_logging.h"

TEST_CASE("[pdm_compare_operator] - PDM_compare_unsigned_ordered_nuplets_int") {

  std::vector<int> a           = {7, 1, 3, 2, 4};
  std::vector<int> a_direct    = {2, 4, 7, 1, 3};
  std::vector<int> a_reverse   = {7, 4, 2, 3, 1};
  std::vector<int> not_a       = {3, 4, 2, 7, 1};
  std::vector<int> still_not_a = {7, 1, 3, 2, 4, 5};

  int res_a           = PDM_compare_unsigned_ordered_nuplets_int(a.size(), a.data(), a          .size(), a          .data());
  int res_a_direct    = PDM_compare_unsigned_ordered_nuplets_int(a.size(), a.data(), a_direct   .size(), a_direct   .data());
  int res_a_reverse   = PDM_compare_unsigned_ordered_nuplets_int(a.size(), a.data(), a_reverse  .size(), a_reverse  .data());
  int res_not_a       = PDM_compare_unsigned_ordered_nuplets_int(a.size(), a.data(), not_a      .size(), not_a      .data());
  int res_still_not_a = PDM_compare_unsigned_ordered_nuplets_int(a.size(), a.data(), still_not_a.size(), still_not_a.data());

  // printf("res_a           = %d\n", res_a          );
  // printf("res_a_direct    = %d\n", res_a_direct   );
  // printf("res_a_reverse   = %d\n", res_a_reverse  );
  // printf("res_not_a       = %d\n", res_not_a      );
  // printf("res_still_not_a = %d\n", res_still_not_a);

  CHECK(res_a           ==  1);
  CHECK(res_a_direct    ==  1);
  CHECK(res_a_reverse   == -1);
  CHECK(res_not_a       ==  0);
  CHECK(res_still_not_a ==  0);

}
