#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_renum.h"


TEST_CASE("[pdm_part_renum] - PDM_part_renum_by_first_appearance") {
  /*
   *  3 triangles on 6 vertices (vertex 6 is not used by any triangle)
   *    e0 = (1, 2, 3)   e1 = (2, 4, 3)   e2 = (4, 5, 3)
   */
  int n_vtx = 6;
  int elmt_vtx_idx[4] = {0, 3, 6, 9};
  int elmt_vtx    [9] = {1, 2, 3,   2, 4, 3,   4, 5, 3};
  int old_to_new[6] = {-1, -1, -1, -1, -1, -1};
  int new_to_old[6] = {-1, -1, -1, -1, -1, -1};
  int n_new = 0;

  SUBCASE("identity traversal") {
    PDM_part_renum_by_first_appearance(3, NULL, elmt_vtx_idx, elmt_vtx, old_to_new, new_to_old, &n_new);
    int expected_old_to_new[6] = {0, 1, 2, 3, 4, -1};
    int expected_new_to_old[5] = {0, 1, 2, 3, 4};
    CHECK(n_new == 5);
    CHECK_EQ_C_ARRAY(old_to_new, expected_old_to_new, n_vtx);
    CHECK_EQ_C_ARRAY(new_to_old, expected_new_to_old, 5);
  }

  SUBCASE("given traversal order, then continuation with signed connectivity") {
    int order[3] = {2, 0, 1};
    PDM_part_renum_by_first_appearance(3, order, elmt_vtx_idx, elmt_vtx, old_to_new, new_to_old, &n_new);
    /* e2 : 4, 5, 3 -> 0, 1, 2 ; e0 : 1, 2 -> 3, 4 ; e1 : nothing new */
    int expected_old_to_new[6] = {3, 4, 2, 0, 1, -1};
    int expected_new_to_old[5] = {3, 4, 2, 0, 1};
    CHECK(n_new == 5);
    CHECK_EQ_C_ARRAY(old_to_new, expected_old_to_new, n_vtx);
    CHECK_EQ_C_ARRAY(new_to_old, expected_new_to_old, 5);

    /* Continue the numbering with a lower dimension entity (edge 6 -> 1, sign ignored) */
    int edge_vtx_idx[2] = {0, 2};
    int edge_vtx    [2] = {6, -1};
    PDM_part_renum_by_first_appearance(1, NULL, edge_vtx_idx, edge_vtx, old_to_new, new_to_old, &n_new);

    int expected_old_to_new2[6] = {3, 4, 2, 0, 1, 5};
    int expected_new_to_old2[6] = {3, 4, 2, 0, 1, 5};
    CHECK(n_new == 6);
    CHECK_EQ_C_ARRAY(old_to_new, expected_old_to_new2, n_vtx);
    CHECK_EQ_C_ARRAY(new_to_old, expected_new_to_old2, n_vtx);
  }
}


TEST_CASE("[pdm_part_renum] - PDM_part_renum_by_lexicographic_order") {
  SUBCASE("sorted lists, prefix first, empty first, ties keep initial order") {
    /*
     *  e0 = (3, 1, 2)  -> (1, 2, 3)
     *  e1 = (2, -5)    -> (2, 5)      (sign ignored)
     *  e2 = (1, 2)     -> (1, 2)      (prefix of e0 : before e0)
     *  e3 = ()         -> ()          (empty : first)
     *  e4 = (2, 1, 3)  -> (1, 2, 3)   (same list as e0 : after e0)
     *  e5 = (1, 4)     -> (1, 4)
     */
    int entity1_entity2_idx[7] = {0, 3, 5, 7, 7, 10, 12};
    int entity1_entity2   [12] = {3, 1, 2,   2, -5,   1, 2,   2, 1, 3,   1, 4};
    int new_to_old[6];
    PDM_part_renum_by_lexicographic_order(6, entity1_entity2_idx, entity1_entity2, new_to_old);
    int expected[6] = {3, 2, 0, 4, 5, 1};
    CHECK_EQ_C_ARRAY(new_to_old, expected, 6);
    /* Input connectivity is not modified */
    int expected_connectivity[12] = {3, 1, 2,   2, -5,   1, 2,   2, 1, 3,   1, 4};
    CHECK_EQ_C_ARRAY(entity1_entity2, expected_connectivity, 12);
  }

  SUBCASE("already ordered stays unchanged, empty input") {
    int entity1_entity2_idx[4] = {0, 2, 4, 6};
    int entity1_entity2    [6] = {1, 2,   1, 3,   2, 3};
    int new_to_old[3];
    PDM_part_renum_by_lexicographic_order(3, entity1_entity2_idx, entity1_entity2, new_to_old);
    int expected[3] = {0, 1, 2};
    CHECK_EQ_C_ARRAY(new_to_old, expected, 3);
    int empty_idx[1] = {0};
    PDM_part_renum_by_lexicographic_order(0, empty_idx, NULL, NULL);
  }
}
