#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_distant_neighbor.h"
#include "pdm_logging.h"


//  Parenthesis is fields
//         Part 1                        Part 2
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   3(4) |   2(3) |           |   2(-3)|   3(-4)|
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   0(1) |   1(2) |           |   0(-1)|   1(-2)|
//   +--------+--------+           +--------+--------+
//
//

MPI_TEST_CASE("[1p] distant neighbor ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_part = 2;
  int n_entity[2] = {4, 4};

  std::vector<int> neighbor_idx_p1 = { 0, 0, 1, 2, 2};
  std::vector<int> neighbor_idx_p2 = { 0, 1, 1, 2, 2};

  std::vector<int> neighbor_desc_p1 = { 0, 1, 0,
                                        0, 1, 2};

  std::vector<int> neighbor_desc_p2 = { 0, 0, 1,
                                        0, 0, 2};

  std::vector<int*> neighbor_idx  = {neighbor_idx_p1 .data(), neighbor_idx_p2 .data()};
  std::vector<int*> neighbor_desc = {neighbor_desc_p1.data(), neighbor_desc_p2.data()};

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(pdm_comm,
                                                           n_part,
                                                           n_entity,
                                                           neighbor_idx.data(),
                                                           neighbor_desc.data());

  SUBCASE("constant stride ") {
    // Prepare exchange
    std::vector<int> field_p1 = { 1,  2,  3,  4};
    std::vector<int> field_p2 = {-1, -2, -3, -4};

    std::vector<int*> fields = {field_p1.data(), field_p2.data()};

    int** exch_fields = NULL;
    PDM_distant_neighbor_exch(dn,
                              sizeof(int),
                              PDM_STRIDE_CST,
                              1,
                              NULL,
                   (void **)  fields.data(),
                              NULL,
                  (void ***) &exch_fields);

    // PDM_log_trace_array_int(exch_fields[0], neighbor_idx_p1[n_entity[0]], "p1");
    // PDM_log_trace_array_int(exch_fields[1], neighbor_idx_p2[n_entity[1]], "p2");

    int exch_fields_expexted_p1[2] = {-1, -3};
    int exch_fields_expexted_p2[2] = { 2,  3};

    CHECK_EQ_C_ARRAY(exch_fields[0], exch_fields_expexted_p1, neighbor_idx_p1[n_entity[0]]); // Part 1
    CHECK_EQ_C_ARRAY(exch_fields[1], exch_fields_expexted_p2, neighbor_idx_p2[n_entity[1]]); // Part 2

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(exch_fields[i_part]);
    }
    free(exch_fields);
  }

  SUBCASE("variable strid ") {

    // Prepare exchange
    std::vector<int> field_n_p1 = { 2, 1, 4, 3};
    std::vector<int> field_n_p2 = { 3, 4, 2, 1};
    std::vector<int> field_p1   = { 1, 10,
                                    2,
                                    3, 30, 300, 3000,
                                    4, 40, 400};
    std::vector<int> field_p2   = {-1, -10, -100,
                                   -2, -20, -200, -2000,
                                   -3, -30,
                                   -4};

    std::vector<int*> fields   = {field_p1  .data(), field_p2  .data()};
    std::vector<int*> fields_n = {field_n_p1.data(), field_n_p2.data()};

    int** exch_fields_n = NULL;
    int** exch_fields   = NULL;
    PDM_distant_neighbor_exch(dn,
                              sizeof(int),
                              PDM_STRIDE_VAR,
                              -1,
                              fields_n.data(),
                   (void **)  fields.data(),
                             &exch_fields_n,
                  (void ***) &exch_fields);

    int exch_fields_n_expexted_p1[2] = {3,  2};
    int exch_fields_n_expexted_p2[2] = {1,  4};

    CHECK_EQ_C_ARRAY(exch_fields_n[0], exch_fields_n_expexted_p1, neighbor_idx_p1[n_entity[0]]); // Part 1
    CHECK_EQ_C_ARRAY(exch_fields_n[1], exch_fields_n_expexted_p2, neighbor_idx_p2[n_entity[1]]); // Part 2


    std::vector<int> size_fields(n_part);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      size_fields[i_part] = 0;
      for(int i_entity = 0; i_entity < neighbor_idx[i_part][n_entity[i_part]]; ++i_entity) {
        size_fields[i_part] += exch_fields_n[i_part][i_entity];
      }
    }

    CHECK( size_fields[0] == 5);
    CHECK( size_fields[1] == 5);

    int exch_fields_expexted_p1[5] = {-1, -10, -100,  -3, -30 };
    int exch_fields_expexted_p2[5] = { 2,   3,   30, 300, 3000};

    CHECK_EQ_C_ARRAY(exch_fields[0], exch_fields_expexted_p1, size_fields[0]); // Part 1
    CHECK_EQ_C_ARRAY(exch_fields[1], exch_fields_expexted_p2, size_fields[1]); // Part 2

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(exch_fields_n[i_part]);
      free(exch_fields[i_part]);
    }
    free(exch_fields);
    free(exch_fields_n);

  }

  PDM_distant_neighbor_free(dn);
}

// Same as before but neighbor not appear in sorted way
//  Parenthesis is fields
//         Part 1                        Part 2
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   0(1) |   1(2) |           |   2(-3)|   3(-4)|
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   3(4) |   2(3) |           |   0(-1)|   1(-2)|
//   +--------+--------+           +--------+--------+
//
//
MPI_TEST_CASE("[1p] distant neighbor (unsorted) ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_part = 2;
  int n_entity[2] = {4, 4};

  std::vector<int> neighbor_idx_p1 = { 0, 0, 1, 2, 2};
  std::vector<int> neighbor_idx_p2 = { 0, 1, 1, 2, 2};

  std::vector<int> neighbor_desc_p1 = { 0, 1, 2,
                                        0, 1, 0};

  std::vector<int> neighbor_desc_p2 = { 0, 0, 2,
                                        0, 0, 1};

  std::vector<int*> neighbor_idx  = {neighbor_idx_p1 .data(), neighbor_idx_p2 .data()};
  std::vector<int*> neighbor_desc = {neighbor_desc_p1.data(), neighbor_desc_p2.data()};

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(pdm_comm,
                                                           n_part,
                                                           n_entity,
                                                           neighbor_idx.data(),
                                                           neighbor_desc.data());

  // Prepare exchange
  std::vector<int> field_p1 = { 1,  2,  3,  4};
  std::vector<int> field_p2 = {-1, -2, -3, -4};

  std::vector<int*> fields = {field_p1.data(), field_p2.data()};

  int** exch_fields = NULL;
  PDM_distant_neighbor_exch(dn,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            1,
                            NULL,
                 (void **)  fields.data(),
                            NULL,
                (void ***) &exch_fields);

  // PDM_log_trace_array_int(exch_fields[0], neighbor_idx_p1[n_entity[0]], "p1");
  // PDM_log_trace_array_int(exch_fields[1], neighbor_idx_p2[n_entity[1]], "p2");

  int exch_fields_expexted_p1[2] = {-1, -3};
  int exch_fields_expexted_p2[2] = { 2,  3};

  CHECK_EQ_C_ARRAY(exch_fields[0], exch_fields_expexted_p1, neighbor_idx_p1[n_entity[0]]); // Part 1
  CHECK_EQ_C_ARRAY(exch_fields[1], exch_fields_expexted_p2, neighbor_idx_p2[n_entity[1]]); // Part 2

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(exch_fields[i_part]);
  }
  free(exch_fields);

  PDM_distant_neighbor_free(dn);

}


// Same as before but neighbor not appear in sorted way
//  Parenthesis is fields
//         Part 1                        Part 2
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   0(1) |   1(2) |           |   2(-3)|   3(-4)|
//   +--------+--------+           +--------+--------+
//   |        |        |           |        |        |
//   |   3(4) |   2(3) |           |   0(-1)|   1(-2)|
//   +--------+--------+           +--------+--------+
//
//
MPI_TEST_CASE("[1p] distant neighbor (multiple unsorted) ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int n_part = 2;
  int n_entity[2] = {4, 4};

  std::vector<int> neighbor_idx_p1 = { 0, 0, 2, 4, 4};
  std::vector<int> neighbor_idx_p2 = { 0, 2, 2, 4, 4};

  std::vector<int> neighbor_desc_p1 = { 0, 1, 2,
                                        0, 1, 0,
                                        0, 1, 0,
                                        0, 1, 2};

  std::vector<int> neighbor_desc_p2 = { 0, 0, 2,
                                        0, 0, 1,
                                        0, 0, 1,
                                        0, 0, 2};

  std::vector<int*> neighbor_idx  = {neighbor_idx_p1 .data(), neighbor_idx_p2 .data()};
  std::vector<int*> neighbor_desc = {neighbor_desc_p1.data(), neighbor_desc_p2.data()};

  PDM_distant_neighbor_t* dn = PDM_distant_neighbor_create(pdm_comm,
                                                           n_part,
                                                           n_entity,
                                                           neighbor_idx.data(),
                                                           neighbor_desc.data());

  // Prepare exchange
  std::vector<int> field_p1 = { 1,  2,  3,  4};
  std::vector<int> field_p2 = {-1, -2, -3, -4};

  std::vector<int*> fields = {field_p1.data(), field_p2.data()};

  int** exch_fields = NULL;
  PDM_distant_neighbor_exch_int(dn,
                            sizeof(int),
                            PDM_STRIDE_CST,
                            1,
                            NULL,
                            fields.data(),
                            NULL,
                           &exch_fields);

  PDM_log_trace_array_int(exch_fields[0], neighbor_idx_p1[n_entity[0]], "p1:: ");
  PDM_log_trace_array_int(exch_fields[1], neighbor_idx_p2[n_entity[1]], "p2:: ");

  // int exch_fields_expexted_p1[2] = {-1, -3};
  // int exch_fields_expexted_p2[2] = { 2,  3};

  // CHECK_EQ_C_ARRAY(exch_fields[0], exch_fields_expexted_p1, neighbor_idx_p1[n_entity[0]]); // Part 1
  // CHECK_EQ_C_ARRAY(exch_fields[1], exch_fields_expexted_p2, neighbor_idx_p2[n_entity[1]]); // Part 2

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(exch_fields[i_part]);
  }
  free(exch_fields);

  PDM_distant_neighbor_free(dn);

}
