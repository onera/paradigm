#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_unique.h"
#include "pdm_logging.h"
#include "pdm_morton.h"
#include "pdm_part_to_block.h"
#include "pdm_priv.h"


TEST_CASE("[pdm_morton] - PDM_morton_encode") {

  double coords[3] = {0.5, 0.25, 1.};
  PDM_morton_int_t level = 12;

  PDM_morton_code_t code = PDM_morton_encode(3, level, coords);

  // printf("L =  %3u [X, Y, Z] - [%5u %5u %5u] \n", code.L, code.X[0], code.X[1], code.X[2]);

  CHECK(code.L    == 12);
  CHECK(code.X[0] == 2048);
  CHECK(code.X[1] == 1024);
  CHECK(code.X[2] == 4095);

  SUBCASE("get_children") {

    PDM_morton_code_t  child[8];

    PDM_morton_get_children(3, code, child);

    // for(int i = 0; i < 8; ++i) {
    //   printf("L =  %3u [X, Y, Z] - [%5u %5u %5u] \n", child[i].L, child[i].X[0], child[i].X[1], child[i].X[2]);
    // }

    // All children should have level + 1 over her parents
    CHECK(child[0].L    == 13);
    CHECK(child[0].X[0] == 4096);
    CHECK(child[0].X[1] == 2048);
    CHECK(child[0].X[2] == 8190);

    CHECK(child[1].L    == 13);
    CHECK(child[1].X[0] == 4096);
    CHECK(child[1].X[1] == 2048);
    CHECK(child[1].X[2] == 8191);

    CHECK(child[2].L    == 13);
    CHECK(child[2].X[0] == 4096);
    CHECK(child[2].X[1] == 2049);
    CHECK(child[2].X[2] == 8190);

    CHECK(child[3].L    == 13);
    CHECK(child[3].X[0] == 4096);
    CHECK(child[3].X[1] == 2049);
    CHECK(child[3].X[2] == 8191);

    CHECK(child[4].L    == 13);
    CHECK(child[4].X[0] == 4097);
    CHECK(child[4].X[1] == 2048);
    CHECK(child[4].X[2] == 8190);

    CHECK(child[5].L    == 13);
    CHECK(child[5].X[0] == 4097);
    CHECK(child[5].X[1] == 2048);
    CHECK(child[5].X[2] == 8191);

    CHECK(child[6].L    == 13);
    CHECK(child[6].X[0] == 4097);
    CHECK(child[6].X[1] == 2049);
    CHECK(child[6].X[2] == 8190);

    CHECK(child[7].L    == 13);
    CHECK(child[7].X[0] == 4097);
    CHECK(child[7].X[1] == 2049);
    CHECK(child[7].X[2] == 8191);

  }

  SUBCASE(" a > b") {

    PDM_morton_code_t code1;
    PDM_morton_code_t code2;

    code1.L    = 12;
    code1.X[0] = 2048;
    code1.X[1] = 1024;
    code1.X[2] = 4095;

    _Bool res1 = PDM_morton_a_gt_b(code1, code1);

    printf("res1 = %i \n", (int)res1);

    CHECK(res1 == 0);

    code2.L    = 13;
    code2.X[0] = 2048;
    code2.X[1] = 1024;
    code2.X[2] = 4095;

    _Bool res2 = PDM_morton_a_gt_b(code1, code2);
    _Bool res3 = PDM_morton_a_gt_b(code2, code1);
    printf("res2 = %i \n", (int)res2);

    CHECK(res2 == 1);
    CHECK(res3 == 0);

  }

  SUBCASE(" a gtmin b") {

    PDM_morton_code_t code1;
    PDM_morton_code_t code2;

    code1.L    = 4;
    code1.X[0] = 3;
    code1.X[1] = 3;
    code1.X[2] = 2;

    _Bool res1 = PDM_morton_a_gtmin_b(code1, code1);

    // printf("res1 = %i \n", (int)res1);

    CHECK(res1 == 0);

    code2.L    = 6;
    code2.X[0] = 14;
    code2.X[1] = 15;
    code2.X[2] = 9;

    // en fait gt ça permet de tester si TOUT B est au-dessus de A

    _Bool res_gt_min = PDM_morton_a_gtmin_b(code2, code1);
    _Bool res_gt      = PDM_morton_a_gt_b   (code2, code1);

    // printf("res2_gt_min = %i \n", (int)res2_gt_min);
    // printf("res_gt      = %i \n", (int)res_gt);

    CHECK(res_gt_min == 1);
    CHECK(res_gt     == 0);

  }

}


MPI_TEST_CASE("[pdm_morton] - PDM_morton_local_sort", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  double coords[12] = {0.5, 0.25, 1.,
                       0.5, 0.25, 1.,
                       0.25, 0.5, 0.,
                       0., 0., 0.   };

  int n_pts = 4;

  int dim = 3;
  double extents[2*dim];

  PDM_morton_get_coord_extents(dim, n_pts, coords, extents, pdm_comm);
  PDM_extents_conformize(dim, extents, 1e-3);

  PDM_morton_code_t *morton_codes;
  PDM_malloc(morton_codes, n_pts, PDM_morton_code_t);

  const PDM_morton_int_t max_level = PDM_morton_max_level;
  double d[3];
  double s[3];
  PDM_morton_encode_coords(dim,
                           max_level,
                           extents,
                           n_pts,
                           coords,
                           morton_codes, d, s);


  PDM_morton_local_sort(n_pts, morton_codes);

  for(int i = 1; i < n_pts; ++i) {

    _Bool test = PDM_morton_a_gt_b(morton_codes[i - 1], morton_codes[i]);

    CHECK(test == 0);

    // printf("[%i] = L =  %3u [X, Y, Z] - [%5u %5u %5u] \n", morton_codes[i].L, morton_codes[i].X[0], morton_codes[i].X[1], morton_codes[i].X[2]);
  }


  PDM_free(morton_codes);

}
