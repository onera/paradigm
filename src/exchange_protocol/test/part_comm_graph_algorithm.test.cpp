#include <vector>
#include <numeric>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_comm_graph_algorithm.h"
#include "pdm_sort.h"
#include "pdm_vtk.h"
#include <functional>


MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p", 2) {

  // Corresponds to a QUAD of n_vtx_seg = 3
  /*

       p0                 p1
           6                 6
     3 |+++++++| 6     3 |+++++++| 6
       |       |         |       |
   3   |       |  7  2   |       |   7
       |   4   |         |   4   |
     2 |+++++++| 5     2 |+++++++| 5
       |       |         |       |
   1   |       |  5  1   |       |   5
       |       |         |       |
     1 |+++++++| 4     1 |+++++++| 4
           2                 3

  */


  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<int> vn_entity1      = {6, 6};
  std::vector<int> vn_entity2      = {7, 7};
  std::vector<std::vector<int>> ventity_bound = {{4, 1, 1, 1,
                                                  5, 1, 1, 2,
                                                  6, 1, 1, 3 },
                                                 {1, 0, 1, 4,
                                                  2, 0, 1, 5,
                                                  3, 0, 1, 6 }};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 2, 4, 6, 8, 10, 12, 14},
                                                        {0, 2, 4, 6, 8, 10, 12, 14}};

  std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   5, 4,   3, 6,   6, 5},
                                                    {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 2; // nombre de faces de bords attendu

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[8] = {5, 1, 1, 1,   7, 1, 1, 2};
  static int entity_bound_reorder_p1[8] = {1, 0, 1, 5,   2, 0, 1, 7};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 8);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 8);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}



MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - revert sens", 2) {

  // Correspond to a QUAD of n_vtx_seg = 3

/*
       p0                 p1
           6                 6
     3 |+++++++| 6     3 |+++++++| 6
       |       |         |       |
   3   |       |  7  2   |       |   7
       |   4   |         |   4   |
     2 |+++++++| 5     2 |+++++++| 5
       |       |         |       |
   1   |       |  5  1   |       |   5
       |       |         |       |
     1 |+++++++| 4     1 |+++++++| 4
           2                 3

*/
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  std::vector<int> vn_entity_bound = {3, 3};
  std::vector<int> vn_entity1      = {6, 6};
  std::vector<int> vn_entity2      = {7, 7};
  std::vector<std::vector<int>> ventity_bound = {{4, 1, 1, 1,
                                                  5, 1, 1, 2,
                                                  6, 1, 1, 3 },
                                                 {1, 0, 1, 4,
                                                  2, 0, 1, 5,
                                                  3, 0, 1, 6 }};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 2, 4, 6, 8, 10, 12, 14},
                                                        {0, 2, 4, 6, 8, 10, 12, 14}};

  // std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   5, 4,   3, 6,   6, 5},
  //                                                   {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

                                                                                 /*inverted        inverted   */
  //                                                                                |----|          |----|
  std::vector<std::vector<int>> ventity2_entity1 = {{1, 2,   4, 1,   2, 3,   2, 5,   4, 5,   3, 6,   5, 6},
                                                    {2, 1,   3, 2,   4, 1,   2, 5,   5, 4,   3, 6,   6, 5}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 2;

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[8] = {5, 1, 1, -1,   7, 1, 1, -2};
  static int entity_bound_reorder_p1[8] = {1, 0, 1, -5,   2, 0, 1, -7};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 8);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 8);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}




MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - 3D ", 2) {


  // Corresponds to a HEXA of n_vtx_seg = 3
  // here we only represent the boundary between the 2 parts
  //

  //              --------- +18              9+--------
  //                       /|                /|
  //                      / |               / |
  //                  15 /  |              /  |
  //              ----- +   |            6+---|---
  //                   /|   |            /|   |
  //                  / |17 |           / | 4 |
  //              12 /- |-- +17        /  |  8+-------
  //             -- +   |  /|        3+-- |--/|--
  //                |   | / |         |   | / |
  //                |15 |/  |         | 2 |/  |
  //             -- |-- +14 |         |  5+-- |----
  //                |  /|   |         |  /|   |
  //                | / |16 |         | / | 3 |
  //             11 |/ -| --+16       |/  |  7+-------
  //            --- +   |  /         2+-- |- /----
  //                |   | /           | 1 | /
  //                |14 |/            |   |/
  //   x          --|-- +             |  4+------
  //   ^  y         |  /13            |  /
  //   | +          | /               | /
  //   |/           |/                |/
  //   +--->z   --- +                1+------
  //              10

  //                        z=0.5 plane
  // all normals of boundary faces are z-positive

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // Keep for debug
  std::vector<std::vector<double>> vvtx_coords = {{0.0, 0.0, 0.0,    /*1*/
                                                   0.5, 0.0, 0.0,    /*2*/
                                                   1.0, 0.0, 0.0,    /*3*/
                                                   0.0, 0.5, 0.0,    /*4*/
                                                   0.5, 0.5, 0.0,    /*5*/
                                                   1.0, 0.5, 0.0,    /*6*/
                                                   0.0, 1.0, 0.0,    /*7*/
                                                   0.5, 1.0, 0.0,    /*8*/
                                                   1.0, 1.0, 0.0,    /*9*/
                                                   0.0, 0.0, 0.5,    /*10*/
                                                   0.5, 0.0, 0.5,    /*11*/
                                                   1.0, 0.0, 0.5,    /*12*/
                                                   0.0, 0.5, 0.5,    /*13*/
                                                   0.5, 0.5, 0.5,    /*14*/
                                                   1.0, 0.5, 0.5,    /*15*/
                                                   0.0, 1.0, 0.5,    /*16*/
                                                   0.5, 1.0, 0.5,    /*17*/
                                                   1.0, 1.0, 0.5 },  /*18*/

                                                  {0.0, 0.0, 0.5,    /*1*/
                                                   0.5, 0.0, 0.5,    /*2*/
                                                   1.0, 0.0, 0.5,    /*3*/
                                                   0.0, 0.5, 0.5,    /*4*/
                                                   0.5, 0.5, 0.5,    /*5*/
                                                   1.0, 0.5, 0.5,    /*6*/
                                                   0.0, 1.0, 0.5,    /*7*/
                                                   0.5, 1.0, 0.5,    /*8*/
                                                   1.0, 1.0, 0.5,    /*9*/
                                                   0.0, 0.0, 1.0,    /*10*/
                                                   0.5, 0.0, 1.0,    /*11*/
                                                   1.0, 0.0, 1.0,    /*12*/
                                                   0.0, 0.5, 1.0,    /*13*/
                                                   0.5, 0.5, 1.0,    /*14*/
                                                   1.0, 0.5, 1.0,    /*15*/
                                                   0.0, 1.0, 1.0,    /*16*/
                                                   0.5, 1.0, 1.0,    /*17*/
                                                   1.0, 1.0, 1.0}};  /*18*/

  std::vector<int> vn_entity_bound = {9 ,  9};
  std::vector<int> vn_entity1      = {18, 18};
  std::vector<int> vn_entity2      = {20, 20};
  std::vector<std::vector<int>> ventity_bound = {{10, 1, 1, 1,
                                                  11, 1, 1, 2,
                                                  12, 1, 1, 3,
                                                  13, 1, 1, 4,
                                                  14, 1, 1, 5,
                                                  15, 1, 1, 6,
                                                  16, 1, 1, 7,
                                                  17, 1, 1, 8,
                                                  18, 1, 1, 9 },
                                                 {1, 0, 1, 10,
                                                  2, 0, 1, 11,
                                                  3, 0, 1, 12,
                                                  4, 0, 1, 13,
                                                  5, 0, 1, 14,
                                                  6, 0, 1, 15,
                                                  7, 0, 1, 16,
                                                  8, 0, 1, 17,
                                                  9, 0, 1, 18}};


  /* ici les entity2 sont des quads -> donc composés de 4 noeuds à chaque fois*/
  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80},
                                                        {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80}};

  std::vector<std::vector<int>> ventity2_entity1 = {{2, 1, 4, 5,      /*1 */
                                                     2, 5, 6, 3,      /*2 */
                                                     1, 2, 11, 10,    /*3 */
                                                     7, 8, 5, 4,      /*4 */
                                                     1, 10, 13, 4,    /*5 */
                                                     11, 2, 3, 12,    /*6 */
                                                     8, 9, 6, 5,      /*7 */
                                                     11, 2, 5, 14,    /*8 */
                                                     14, 5, 4, 13,    /*9 */
                                                     12, 3, 6, 15,    /*10 */
                                                     15, 6, 5, 14,    /*11 */
                                                     4, 13, 16, 7,    /*12 */
                                                     14, 5, 8, 17,    /*13 */
                                                     14, 13, 10, 11,  /*14  bound*/
                                                     17, 8, 7, 16,    /*15 */
                                                     15, 6, 9, 18,    /*16 */
                                                     15, 14, 11, 12,  /*17  bound*/
                                                     18, 9, 8, 17,    /*18 */
                                                     17, 16, 13, 14,  /*19  bound*/
                                                     18, 17, 14, 15}, /*20  bound*/

                                                    {5, 4, 1, 2,      /*1   bound*/
                                                     6, 5, 2, 3,      /*2   bound*/
                                                     8, 7, 4, 5,      /*3   bound*/
                                                     11, 10, 1, 2,    /*4*/
                                                     9, 8, 5, 6,      /*5   bound*/
                                                     4, 1, 10, 13,    /*6*/
                                                     12, 11, 2, 3,    /*7*/
                                                     11, 2, 5, 14,    /*8*/
                                                     14, 5, 4, 13,    /*9*/
                                                     15, 12, 3, 6,    /*10*/
                                                     15, 6, 5, 14,    /*11*/
                                                     16, 7, 4, 13,    /*12*/
                                                     14, 5, 8, 17,    /*13*/
                                                     14, 13, 10, 11,  /*14*/
                                                     17, 8, 7, 16,    /*15*/
                                                     18, 15, 6, 9,    /*16*/
                                                     15, 14, 11, 12,  /*17*/
                                                     18, 9, 8, 17,    /*18*/
                                                     17, 16, 13, 14,  /*19*/
                                                     18, 17, 14, 15}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 4; // nombre de faces de bords attendu

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[16] = {14, 1, 1,  1,   17, 1, 1,  2,   19, 1, 1,  3,   20, 1, 1, 5};
  static int entity_bound_reorder_p1[16] = { 1, 0, 1, 14,    2, 0, 1, 17,    3, 0, 1, 19,    5, 0, 1,20};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 16);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");

      int *face_tag = (int *) malloc(pn_entity2 * sizeof(int));

      for(int i = 0; i < pn_entity2; ++i) {
        face_tag[i] = -1;
      }

      for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {
        int i_face = pentity2_graph[i_part][4*idx]-1;
        int t_rank = pentity2_graph[i_part][4*idx+1];
        face_tag[i_face] = t_rank;
      }

      const char* field_name[] = {"face_tag", 0 };
      const int*  field     [] = {face_tag};

      char filename[999];
      sprintf(filename, "out_face_graph_i_part=%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 pn_entity1,
                                 vvtx_coords[i_rank].data(),
                                 NULL,
                                 PDM_MESH_NODAL_QUAD4,
                                 pn_entity2,
                                 entity2_entity1,
                                 NULL,
                                 1,
                                 field_name,
                                 (const int **)   &field);


      free(face_tag);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}



MPI_TEST_CASE("[PDM_part_comm_graph_entity1_to_entity2] - 1 part - 2p - 3D - revert face_vtx", 2) {


  // Correspond to a HEXA of n_vtx_seg = 3

  // here we only represent the boundary between the 2 parts
  //

  //              --------- +18              9+--------
  //                       /|                /|
  //                      / |               / |
  //                  15 /  |              /  |
  //              ----- +   |            6+---|---
  //                   /|   |            /|   |
  //                  / |17 |           / | 4 |
  //              12 /- |-- +17        /  |  8+-------
  //             -- +   |  /|        3+-- |--/|--
  //                |   | / |         |   | / |
  //                |15 |/  |         | 2 |/  |
  //             -- |-- +14 |         |  5+-- |----
  //                |  /|   |         |  /|   |
  //                | / |16 |         | / | 3 |
  //             11 |/ -| --+16       |/  |  7+-------
  //            --- +   |  /         2+-- |- /----
  //                |   | /           | 1 | /
  //                |14 |/            |   |/
  //   x          --|-- +             |  4+------
  //   ^  y         |  /13            |  /
  //   | +          | /               | /
  //   |/           |/                |/
  //   +--->z   --- +                1+------
  //              10

  //                        z=0.5 plane
  // p1 : all normals of boundary faces are z-negative
  // p2 :  "    "     "     "      "     "  z-positive
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_part = 1;

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // Keep for debug
  std::vector<std::vector<double>> vvtx_coords = {{0.0, 0.0, 0.0,
                                                   0.5, 0.0, 0.0,
                                                   1.0, 0.0, 0.0,
                                                   0.0, 0.5, 0.0,
                                                   0.5, 0.5, 0.0,
                                                   1.0, 0.5, 0.0,
                                                   0.0, 1.0, 0.0,
                                                   0.5, 1.0, 0.0,
                                                   1.0, 1.0, 0.0,
                                                   0.0, 0.0, 0.5,
                                                   0.5, 0.0, 0.5,
                                                   1.0, 0.0, 0.5,
                                                   0.0, 0.5, 0.5,
                                                   0.5, 0.5, 0.5,
                                                   1.0, 0.5, 0.5,
                                                   0.0, 1.0, 0.5,
                                                   0.5, 1.0, 0.5,
                                                   1.0, 1.0, 0.5 },
                                                  {0.0, 0.0, 0.5,
                                                   0.5, 0.0, 0.5,
                                                   1.0, 0.0, 0.5,
                                                   0.0, 0.5, 0.5,
                                                   0.5, 0.5, 0.5,
                                                   1.0, 0.5, 0.5,
                                                   0.0, 1.0, 0.5,
                                                   0.5, 1.0, 0.5,
                                                   1.0, 1.0, 0.5,
                                                   0.0, 0.0, 1.0,
                                                   0.5, 0.0, 1.0,
                                                   1.0, 0.0, 1.0,
                                                   0.0, 0.5, 1.0,
                                                   0.5, 0.5, 1.0,
                                                   1.0, 0.5, 1.0,
                                                   0.0, 1.0, 1.0,
                                                   0.5, 1.0, 1.0,
                                                   1.0, 1.0, 1.0}};

  std::vector<int> vn_entity_bound = {9 ,  9};
  std::vector<int> vn_entity1      = {18, 18};
  std::vector<int> vn_entity2      = {20, 20};
  std::vector<std::vector<int>> ventity_bound = {{10, 1, 1, 1,
                                                  11, 1, 1, 2,
                                                  12, 1, 1, 3,
                                                  13, 1, 1, 4,
                                                  14, 1, 1, 5,
                                                  15, 1, 1, 6,
                                                  16, 1, 1, 7,
                                                  17, 1, 1, 8,
                                                  18, 1, 1, 9 },
                                                 {1, 0, 1, 10,
                                                  2, 0, 1, 11,
                                                  3, 0, 1, 12,
                                                  4, 0, 1, 13,
                                                  5, 0, 1, 14,
                                                  6, 0, 1, 15,
                                                  7, 0, 1, 16,
                                                  8, 0, 1, 17,
                                                  9, 0, 1, 18}};

  std::vector<std::vector<int>> ventity2_entity1_idx = {{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80},
                                                        {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80}};

  std::vector<std::vector<int>> ventity2_entity1 = {{2, 1, 4, 5,      // 1
                                                     2, 5, 6, 3,      // 2
                                                     1, 2, 11, 10,    // 3
                                                     7, 8, 5, 4,      // 4
                                                     1, 10, 13, 4,    // 5
                                                     11, 2, 3, 12,    // 6
                                                     8, 9, 6, 5,      // 7
                                                     11, 2, 5, 14,    // 8
                                                     14, 5, 4, 13,    // 9
                                                     12, 3, 6, 15,    // 10
                                                     15, 6, 5, 14,    // 11
                                                     4, 13, 16, 7,    // 12
                                                     14, 5, 8, 17,    // 13
                                                     11, 10, 13, 14,  // 14 original : 14, 13, 10, 11,
                                                     17, 8, 7, 16,    // 15
                                                     15, 6, 9, 18,    // 16
                                                     12, 11, 14, 15,  // 17 original : 15, 14, 11, 12
                                                     18, 9, 8, 17,    // 18
                                                     14, 13, 16, 17,  // 19 origianl : 17, 16, 13, 14
                                                     15, 14, 17, 18}, // 20 original : 18, 17, 14, 15
                                                    {5, 4, 1, 2,
                                                     6, 5, 2, 3,
                                                     8, 7, 4, 5,
                                                     11, 10, 1, 2,
                                                     9, 8, 5, 6,
                                                     4, 1, 10, 13,
                                                     12, 11, 2, 3,
                                                     11, 2, 5, 14,
                                                     14, 5, 4, 13,
                                                     15, 12, 3, 6,
                                                     15, 6, 5, 14,
                                                     16, 7, 4, 13,
                                                     14, 5, 8, 17,
                                                     14, 13, 10, 11,
                                                     17, 8, 7, 16,
                                                     18, 15, 6, 9,
                                                     15, 14, 11, 12,
                                                     18, 9, 8, 17,
                                                     17, 16, 13, 14,
                                                     18, 17, 14, 15}};

  int n_entity_bound       = vn_entity_bound     [i_rank];
  int *entity_bound        = ventity_bound       [i_rank].data();
  int *entity2_entity1_idx = ventity2_entity1_idx[i_rank].data();
  int *entity2_entity1     = ventity2_entity1    [i_rank].data();
  int pn_entity1           = vn_entity1          [i_rank];
  int pn_entity2           = vn_entity2          [i_rank];

  int  *pn_entity2_graph = NULL;
  int **pentity2_graph   = NULL;
  PDM_part_comm_graph_entity1_to_entity2(pdm_comm,
                                         n_part,
                                         &n_entity_bound,
                                         &entity_bound,
                                         0,
                                         NULL,
                                         &pn_entity1,
                                         &pn_entity2,
                                         &entity2_entity1_idx,
                                         &entity2_entity1,
                                         &pn_entity2_graph,
                                         &pentity2_graph,
                                         NULL);

  int pn_entity2_graph_expected = 4;

  CHECK(pn_entity2_graph_expected == pn_entity2_graph[0]);

  static int entity_bound_reorder_p0[16] = {14, 1, 1,  -1, 17, 1, 1,  -2, 19, 1, 1, -3, 20, 1, 1, -5};
  static int entity_bound_reorder_p1[16] = { 1, 0, 1, -14,  2, 0, 1, -17,  3, 0, 1,-19,  5, 0, 1,-20};

  MPI_CHECK_EQ_C_ARRAY(0, pentity2_graph[0], entity_bound_reorder_p0, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pentity2_graph[0], entity_bound_reorder_p1, 16);

  if(1 == 0) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(pentity2_graph[i_part], 4 * pn_entity2_graph[i_part], "pentity2_graph ::");

      int *face_tag = (int *) malloc(pn_entity2 * sizeof(int));

      for(int i = 0; i < pn_entity2; ++i) {
        face_tag[i] = -1;
      }

      for(int idx = 0; idx < pn_entity2_graph[i_part]; ++idx) {
        int i_face = pentity2_graph[i_part][4*idx]-1;
        int t_rank = pentity2_graph[i_part][4*idx+1];
        face_tag[i_face] = t_rank;
      }

      const char* field_name[] = {"face_tag", 0 };
      const int*  field     [] = {face_tag};

      char filename[999];
      sprintf(filename, "out_face_graph_i_part=%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                 pn_entity1,
                                 vvtx_coords[i_rank].data(),
                                 NULL,
                                 PDM_MESH_NODAL_QUAD4,
                                 pn_entity2,
                                 entity2_entity1,
                                 NULL,
                                 1,
                                 field_name,
                                 (const int **)   &field);


      free(face_tag);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pentity2_graph[i_part]);
  }
  free(pentity2_graph);
  free(pn_entity2_graph);

}


MPI_TEST_CASE("[PDM_part_comm_graph] - selected_entity1_to_selected_entity2 - 2p", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /**
   *         /         7 ----- 8 ----- 9 ---- 10
   *         |         |       |       |.......|
   *         |         |   2   |   3   |.. 4 ..|
   *         |         |       |       |.......|
   * Rank 1 <          3 ----- 4 ----- 5 ----- 6
   *         |        /       /       /|.......|
   *         |       /       /       / |.. 1 ..|
   *         |      /       /       /  |.......|
   *         \     /       /       /   1 ----- 2
   *         /    9 ---- 10 ---- 11   /       /
   *         |    |       |       |  /       /
   *         |    |   4   |   5   | /       /
   *         |    |       |       |/       /
   * Rank 0 <     5 ----- 6 ----- 7 ----- 8
   *         |    |       |       |.......|
   *         |    |   1   |   2   |.. 3 ..|
   *         |    |       |       |.......|
   *         \    1 ----- 2 ----- 3 ----- 4
   */


  // Element->vertex connectivity
  std::vector<std::vector<int>> velt_vtx_idx = {
    {0, 4, 8, 12, 16, 20},
    {0, 4, 8, 12, 16}
  };
  std::vector<std::vector<int>> velt_vtx = {
    {1, 2, 6, 5,
     2, 3, 7, 6,
     3, 4, 8, 7,
     5, 6, 10, 9,
     6, 7, 11, 10},
    {1, 2, 6, 5,
     3, 4, 8, 7,
     4, 5, 9, 8,
     5, 6, 10, 9}
  };

  // Selected elements
  std::vector<int> vn_selected_elt = {1, 2};
  std::vector<std::vector<int>> vselected_elt = {
    {3},
    {1, 4}
  };

  // Inter-partition vtx communication graph
  std::vector<int> vn_bound_vtx = {5, 5};
  std::vector<std::vector<int>> vbound_vtx = {
    { 7, 1, 1, 1,
      8, 1, 1, 2,
      9, 1, 1, 3,
     10, 1, 1, 4,
     11, 1, 1, 5},
    {1, 0, 1, 7,
     2, 0, 1, 8,
     3, 0, 1, 9,
     4, 0, 1, 10,
     5, 0, 1, 11}
  };

  int n_part = 1;

  // Create Part Comm Graph for vertices
  int  n_bound_vtx = vn_bound_vtx[i_rank];
  int *bound_vtx   = vbound_vtx  [i_rank].data();
  PDM_part_comm_graph_t *pcg_vtx = PDM_part_comm_graph_create(n_part,
                                                              &n_bound_vtx,
                                                              &bound_vtx,
                                                              PDM_OWNERSHIP_USER,
                                                              pdm_comm);

  // Deduce selected vertices from selected elements
  int  n_selected_elt = vn_selected_elt[i_rank];
  int *selected_elt   = vselected_elt  [i_rank].data();
  int *elt_vtx_idx    = velt_vtx_idx   [i_rank].data();
  int *elt_vtx        = velt_vtx       [i_rank].data();

  int  *n_selected_vtx = NULL;
  int **selected_vtx   = NULL;
  PDM_part_comm_graph_selected_entity1_to_selected_entity2(&n_selected_elt,
                                                           &selected_elt,
                                                           &elt_vtx_idx,
                                                           &elt_vtx,
                                                           pcg_vtx,
                                                           &n_selected_vtx,
                                                           &selected_vtx);

  // Check selected vertices
  PDM_sort_int(selected_vtx[0], NULL, n_selected_vtx[0]);

  std::vector<std::vector<int>> vexpected_selected_vtx = {
    {3, 4, 7, 8, 11},
    {1, 2, 5, 6, 9, 10}
  };

  CHECK(n_selected_vtx[0] == vexpected_selected_vtx[i_rank].size());
  MPI_CHECK_EQ_C_ARRAY(i_rank, selected_vtx[0], vexpected_selected_vtx[i_rank].data(), n_selected_vtx[0]);

  // Free memory
  PDM_free(n_selected_vtx);
  PDM_free(selected_vtx[0]);
  PDM_free(selected_vtx);

  PDM_part_comm_graph_free(pcg_vtx);
}




MPI_TEST_CASE("[PDM_part_comm_graph_concatenate] - 1 part - 1 perio - 2p", 2) {
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  int i_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  /*
   *              9 +---+---+---+ 12
   *                |           |
   *              5 +   rank 1  + 8
   *                |           |
   *                +---+---+---+
   *                1   2   3   4
   * interface -1                   interface +1
   *                9  10  11  12
   *                +---+---+---+
   *                |           |
   *              5 +   rank 0  + 8
   *                |           |
   *              1 +---+---+---+ 4
   *
   * --- rank 0 ---
   *  1 -> (0, 1,  4) through interface -1
   *
   *  4 -> (0, 1,  1) through interface  1
   *
   *  5 -> (0, 1,  8) through interface -1
   *
   *  8 -> (0, 1,  5) through interface  1
   *
   *  9 -> (1, 1,  1) through interface 0
   *    -> (0, 1, 12) through interface -1
   *    -> (1, 1,  4) through interface -1
   *
   * 10 -> (1, 1,  2) through interface 0
   *
   * 11 -> (1, 1,  3) through interface 0
   *
   * 12 -> (1, 1,  4) through interface 0
   * 12 -> (0, 1,  9) through interface 1
   * 12 -> (1, 1,  1) through interface 1
   *
   * owners : 1, 5, 9, 10, 11, 12
   *
   *
   * --- rank 1 ---
   *  1 -> (0, 1,  9) through interface 0
   *    -> (1, 1,  4) through interface -1
   *    -> (0, 1, 12) through interface -1
   *
   *  2 -> (0, 1, 10) through interface 0
   *
   *  3 -> (0, 1, 11) through interface 0
   *
   *  4 -> (0, 1, 12) through interface 0
   *    -> (1, 1,  1) through interface 1
   *    -> (0, 1,  9) through interface 1
   *
   *  5 -> (1, 1,  8) through interface -1
   *
   *  8 -> (1, 1,  5) through interface  1
   *
   *  9 -> (1, 1, 12) through interface -1
   *
   * 12 -> (1, 1,  9) through interface  1
   *
   * owners : 5, 9
   *
   */

  /* Part */
  int n_part = 1;

  /* Comm graph intra*/
  std::vector<int> vn_entity_intra = {4, 4};
  std::vector<std::vector<int>> ventity_intra = {{9,  1, 1,  1,
                                                  10, 1, 1,  2,
                                                  11, 1, 1,  3,
                                                  12, 1, 1,  4},
                                                 {1,  0, 1,  9,
                                                  2,  0, 1, 10,
                                                  3,  0, 1, 11,
                                                  4,  0, 1, 12}};

  int n_entity_intra = vn_entity_intra[i_rank];
  int *entity_intra  = ventity_intra  [i_rank].data();

  PDM_part_comm_graph_t *pcg_intra = PDM_part_comm_graph_create(n_part,
                                                               &n_entity_intra,
                                                               &entity_intra,
                                                               PDM_OWNERSHIP_USER,
                                                               pdm_comm);

  /* Comm graph perio*/
  std::vector<int> vn_entity_perio = {8, 8};
  std::vector<std::vector<int>> ventity_perio = {{1,  0, 1,  4,
                                                  4,  0, 1,  1,
                                                  5,  0, 1,  8,
                                                  8,  0, 1,  5,
                                                  9,  0, 1, 12,
                                                  9,  1, 1,  4,
                                                  12, 0, 1,  9,
                                                  12, 1, 1,  1},
                                                 {1,  1, 1,  4,
                                                  1,  0, 1, 12,
                                                  4,  1, 1,  1,
                                                  4,  0, 1,  9,
                                                  5,  1, 1,  8,
                                                  8,  1, 1,  5,
                                                  9,  1, 1, 12,
                                                  12, 1, 1,  9}};
  std::vector<std::vector<int>> ventity_perio_nplt = {{-1,  1, -1, 1, -1, -1,  1, 1},
                                                      {-1, -1,  1, 1, -1,  1, -1, 1}};

  int n_entity_perio     = vn_entity_perio   [i_rank];
  int *entity_perio      = ventity_perio     [i_rank].data();
  int *entity_perio_nplt = ventity_perio_nplt[i_rank].data();

  PDM_part_comm_graph_t *pcg_perio = PDM_part_comm_graph_with_nuplet_create(n_part,
                                                                           &n_entity_perio,
                                                                           &entity_perio,
                                                                            PDM_OWNERSHIP_USER,
                                                                            1,
                                                                           &entity_perio_nplt,
                                                                            PDM_OWNERSHIP_USER,
                                                                            PDM_TRUE,
                                                                            pdm_comm);

  PDM_part_comm_graph_t *pcgs[2] = {pcg_intra, pcg_perio};
  PDM_part_comm_graph_t *pcg_full = PDM_part_comm_graph_concatenate(pdm_comm,
                                                                    2,
                                                                    pcgs);

  /* Comm graph */
  std::vector<int> vn_entity_full = {12, 12};
  std::vector<std::vector<int>> ventity_full = {{
    9,  1, 1,  1,
    10, 1, 1,  2,
    11, 1, 1,  3,
    12, 1, 1,  4,
    1,  0, 1,  4,
    4,  0, 1,  1,
    5,  0, 1,  8,
    8,  0, 1,  5,
    9,  0, 1, 12,
    9,  1, 1,  4,
    12, 0, 1,  9,
    12, 1, 1,  1
  },
  {
    1,  0, 1,  9,
    2,  0, 1, 10,
    3,  0, 1, 11,
    4,  0, 1, 12,
    1,  1, 1,  4,
    1,  0, 1, 12,
    4,  1, 1,  1,
    4,  0, 1,  9,
    5,  1, 1,  8,
    8,  1, 1,  5,
    9,  1, 1, 12,
    12, 1, 1,  9
  }};
  std::vector<std::vector<int>> ventity_full_nplt = {{0, 0,  0, 0, -1,  1, -1, 1, -1, -1,  1, 1},
                                                     {0, 0, 0,  0, -1, -1,  1, 1, -1,  1, -1, 1}};

  int *entity_full      = NULL;
  int *entity_full_nplt = NULL;
  int n_entity_full = PDM_part_comm_graph_entity_graph_get(pcg_full,
                                                           0,
                                                          &entity_full,
                                                           PDM_OWNERSHIP_BAD_VALUE);
  PDM_part_comm_graph_entity_nuplet_get(pcg_full,
                                        0,
                                       &entity_full_nplt,
                                        PDM_OWNERSHIP_BAD_VALUE);

  int *expected_entity_full      = ventity_full     [i_rank].data();
  int *expected_entity_full_nplt = ventity_full_nplt[i_rank].data();

  MPI_CHECK_EQ_C_ARRAY(0, entity_full     , expected_entity_full     , 4*n_entity_full);
  MPI_CHECK_EQ_C_ARRAY(0, entity_full_nplt, expected_entity_full_nplt,   n_entity_full);
  MPI_CHECK_EQ_C_ARRAY(1, entity_full     , expected_entity_full     , 4*n_entity_full);
  MPI_CHECK_EQ_C_ARRAY(1, entity_full_nplt, expected_entity_full_nplt,   n_entity_full);

  PDM_part_comm_graph_free(pcg_intra);
  PDM_part_comm_graph_free(pcg_perio);
  PDM_part_comm_graph_free(pcg_full);
}