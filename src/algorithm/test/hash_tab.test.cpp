
#include <vector>
#include "doctest/doctest.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_doctest.h"
#include "pdm_hash_tab.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"


template<typename T> PDM_hash_tab_key_t get_key_type();

template<>
PDM_hash_tab_key_t get_key_type<int>() {
  return PDM_HASH_TAB_KEY_INT;
}

#ifdef PDM_LONG_G_NUM
template<>
PDM_hash_tab_key_t get_key_type<PDM_g_num_t>() {
  return PDM_HASH_TAB_KEY_LONG;
}
#endif

template<typename T>
void
test_hash_tab()
{
  auto key_type = get_key_type<T>();

  T key_max = 100;
  PDM_hash_tab_t* ht = PDM_hash_tab_create(key_type,
                                           (void *) &key_max);

  CHECK( *((int*)PDM_hash_tab_key_max_get (ht)) == key_max);
  CHECK( PDM_hash_tab_key_type_get(ht) == key_type);

  T key1  = 8;
  int data1 = 80;
  PDM_hash_tab_data_add(ht, &key1, &data1);

  T key2  = 8;
  int data2 = 800;
  PDM_hash_tab_data_add(ht, &key2, &data2);

  T key3  = 1;
  int data3 = 51;
  PDM_hash_tab_data_add(ht, &key3, &data3);

  /* Check for key 1 */
  T check_key1 = 8;
  const int n_data_check_key1 = PDM_hash_tab_n_data_get (ht, (void *) &check_key1);
  CHECK(n_data_check_key1 == 2);

  /* Check for key 2 */
  T check_key2 = 1;
  const int n_data_check_key2 = PDM_hash_tab_n_data_get (ht, (void *) &check_key2);
  CHECK(n_data_check_key2 == 1);

  int **data = (int **) PDM_hash_tab_data_get (ht, &check_key1);

  CHECK(data[0][0] == 80 );
  CHECK(data[1][0] == 800);

  int n_used_key = PDM_hash_tab_n_used_keys_get(ht);

  CHECK(n_used_key == 2);

  PDM_g_num_t* key_info = PDM_hash_tab_used_keys_get(ht);

  PDM_g_num_t expected_key_info[2] = {8, 1};

  CHECK_EQ_C_ARRAY(key_info, expected_key_info, n_used_key);

  // Ne marche pas en gnum car la fonction est pas général
  // int res_key = 0;
  // int found = PDM_hash_tab_check_collision(ht, 5, key_max, &res_key);
  // CHECK(res_key == 5);
  // CHECK(found   == 0);

  // Purge hash tab
  PDM_hash_tab_purge(ht, PDM_FALSE);  // If we set TRUE, data is free

  T check_key1_purge = 8;
  const int n_data_check_key1_purge = PDM_hash_tab_n_data_get (ht, (void *) &check_key1_purge);
  CHECK(n_data_check_key1_purge == 0);

  T check_key2_purge = 1;
  const int n_data_check_key2_purge = PDM_hash_tab_n_data_get (ht, (void *) &check_key2_purge);
  CHECK(n_data_check_key2_purge == 0);

  PDM_hash_tab_free(ht);
}

TEST_CASE("[pdm_hash_tab] - PDM_HASH_TAB_KEY_INT") {
  test_hash_tab<int>();
}

TEST_CASE("[pdm_hash_tab] - PDM_HASH_TAB_KEY_LONG") {
  test_hash_tab<PDM_g_num_t>();
}



TEST_CASE("[pdm_hash_tab] - recover_edge") {
  // Issued from QUAD4 dcube_nodal with n_vtx_seg=3

  std::vector<int> face_vtx_idx = {0, 4, 8, 12, 16};
  std::vector<int> face_vtx     = {1, 2, 5, 4, 2, 3, 6, 5, 4, 5, 8, 7, 5, 6, 9, 8};

  std::vector<int> ridge_vtx     = {1, 2, 2, 3, 8, 7, 9, 8, 4, 1, 7, 4, 3, 6, 6, 9};

  int n_quad  = 4;
  int n_ridge = 8;
  int n_vtx   = 9;
  int key_max = n_vtx / 2;

  int *work_data = NULL;

  int n_data_work_max = ( 2 * ( 2 + 1 ) ) * n_ridge + ( 4 * (2 + 1)) * n_quad; // 4 arêtes + 1 voisins
  PDM_malloc(work_data, n_data_work_max, int);

  PDM_hash_tab_t* ht = PDM_hash_tab_create(PDM_HASH_TAB_KEY_INT, (void *) &key_max);
  int n_work = 0;
  for(int i = 0; i < n_quad; ++i) {
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      int jp = (j == face_vtx_idx[i+1] - 1) ? face_vtx_idx[i] : j+1;
      int key = ( face_vtx[j] + face_vtx[jp] ) % key_max;

      int *data = &work_data[n_work];

      data[0] = face_vtx[j];
      data[1] = face_vtx[jp];
      data[2] = (i+1);

      // Add in hash_tab all information to sort it after all
      PDM_hash_tab_data_add(ht, &key, (void*) data);

      n_work += 3;

    }
  }

  // Add ridge
  for(int i = 0; i < n_ridge; ++i) {
    int i_vtx1 = ridge_vtx[2*i  ];
    int i_vtx2 = ridge_vtx[2*i+1];
    int key = ( i_vtx1 + i_vtx2 ) % key_max;

    int *data = &work_data[n_work];

    data[0] = i_vtx1;
    data[1] = i_vtx2;
    data[2] = -(i+1);

    // Add in hash_tab all information to sort it after all
    PDM_hash_tab_data_add(ht, &key, (void*) data);
    n_work += 3;
  }

  // Examine hash table to allocated memory
  int n_conflict = PDM_hash_tab_n_used_keys_get(ht);
  PDM_g_num_t *conflict_keys = PDM_hash_tab_used_keys_get(ht);

  int n_data_tot = 0;
  int n_data_max = 0;
  for(int i = 0; i < n_conflict; ++i) {
    int _key = (int) conflict_keys[i];
    int n_data = PDM_hash_tab_n_data_get(ht, &_key);
    n_data_tot += n_data;
    n_data_max = PDM_MAX(n_data_max, n_data);
  }

  int         *edge_vtx    = NULL;
  int         *edge_face   = NULL;
  int         *is_solved   = NULL;
  PDM_g_num_t *is_bnd_edge = NULL;
  PDM_malloc(edge_vtx   , 2 * n_data_tot, int        );
  PDM_malloc(edge_face  , 2 * n_data_tot, int        );
  PDM_malloc(is_solved  ,     n_data_max, int        );
  PDM_malloc(is_bnd_edge,     n_data_tot, PDM_g_num_t);

  // Loop over hash table to solve conflict
  int n_edge = 0;
  for(int i = 0; i < n_conflict; ++i) {
    int _key = (int) conflict_keys[i];
    int n_data = PDM_hash_tab_n_data_get(ht, &_key);

    int **data = (int **) PDM_hash_tab_data_get(ht, &_key);

    for(int i_data = 0; i_data < n_data; ++i_data) {
      is_solved[i_data] = 0;
    }

    // printf("conflict = %i - n_data = %i \n", i, n_data);

    for(int i_data = 0; i_data < n_data; ++i_data) {
      if(is_solved[i_data] == 1) {
        continue;
      }

      // Attention en 3d une arête appartient a plusieurs cellules !!!! --> Pas notre cas ici because 2d
      int c1_vtx1 = data[i_data][0];
      int c1_vtx2 = data[i_data][1];
      int c1_elmt = data[i_data][2];
      int c1_sgn  = PDM_SIGN(c1_elmt);

      // printf(" \t i_data = %i - (%i/%i) - %i \n", i_data, c1_vtx1, c1_vtx2, c1_elmt);

      for(int i_data_opp = i_data+1; i_data_opp < n_data; ++i_data_opp) {
        if(is_solved[i_data_opp] == 1) {
          continue;
        }
        int c2_vtx1 = data[i_data_opp][0];
        int c2_vtx2 = data[i_data_opp][1];
        int c2_elmt = data[i_data_opp][2];
        int c2_sgn  = PDM_SIGN(c2_elmt);

        // printf(" \t\t i_data_opp = %i - (%i/%i) - %i \n", i_data, c2_vtx1, c2_vtx2, c2_elmt);

        if( ((c1_vtx1 == c2_vtx1) && (c1_vtx2 == c2_vtx2)) ||
            ((c1_vtx1 == c2_vtx2) && (c1_vtx2 == c2_vtx1))) {

          is_bnd_edge[n_edge] = 0;
          if(c1_sgn == 1 && c2_sgn == 1) { // Interior faces
            edge_vtx[2*n_edge  ] = c1_vtx1;
            edge_vtx[2*n_edge+1] = c1_vtx2;

            edge_face[2*n_edge  ] = c1_elmt;
            edge_face[2*n_edge+1] = c2_elmt;
          } else {

            if(c1_sgn == 1 && c2_sgn == -1) {
              edge_vtx[2*n_edge  ] = c1_vtx1;
              edge_vtx[2*n_edge+1] = c1_vtx2;

              edge_face[2*n_edge  ] = c1_elmt;
              edge_face[2*n_edge+1] = 0;

              is_bnd_edge[n_edge] = c2_elmt;
            } else if (c1_sgn == -1 && c2_sgn == 1) {
              edge_vtx[2*n_edge  ] = c2_vtx1;
              edge_vtx[2*n_edge+1] = c2_vtx2;

              edge_face[2*n_edge  ] = c2_elmt;
              edge_face[2*n_edge+1] = 0;

              is_bnd_edge[n_edge] = c1_elmt;
            }
          }

          // printf(" \t\t\t Match !!!! %i \n", n_edge);
          n_edge++;

          is_solved[i_data_opp] = 1;
        }
      }

      is_solved[i_data] = 1;
    }

    // Check if all solved
    for(int i_data = 0; i_data < n_data; ++i_data) {
      CHECK(is_solved[i_data] == 1);
    }
  }

  CHECK(n_edge == 12);

  if(0 == 1) {

    double vtx_coords[27] = {0.,0.,0.,5.e-01,0.,0.,1.,0.,0.,0.,5.e-01,0.,5.e-01,5.e-01,0.,1.,5.e-01,0.,0.,1.,0.,5.e-01,1.,0.,1.,1.,0.};

    PDM_vtk_write_std_elements("debug_edge.vtk",
                               n_vtx,
                               vtx_coords,
                               NULL,
                               PDM_MESH_NODAL_BAR2,
                               n_edge,
                               edge_vtx,
                               is_bnd_edge,
                               0,
                               NULL,
                               NULL);
  }

  if(0 == 1) {
    PDM_log_trace_array_int (edge_vtx   , 2 * n_edge, "edge_vtx    ::");
    PDM_log_trace_array_int (edge_face  , 2 * n_edge, "edge_face   ::");
    PDM_log_trace_array_long(is_bnd_edge,     n_edge, "is_bnd_edge ::");
  }

  int         expected_edge_vtx   [24] = {1, 2, 2, 5, 6, 5, 8, 7, 7, 4, 6, 9, 5, 4, 4, 1, 2, 3, 3, 6, 5, 8, 9, 8};
  int         expected_edge_face  [24] = {1, 0, 1, 2, 2, 4, 3, 0, 3, 0, 4, 0, 1, 3, 1, 0, 2, 0, 2, 0, 3, 4, 4, 0};
  PDM_g_num_t expected_is_bnd_edge[12] = {-1, 0, 0, -3, -6, -8, 0, -5, -2, -7, 0, -4};

  CHECK_EQ_C_ARRAY(edge_vtx   , expected_edge_vtx   , 2 * n_edge);
  CHECK_EQ_C_ARRAY(edge_face  , expected_edge_face  , 2 * n_edge);
  CHECK_EQ_C_ARRAY(is_bnd_edge, expected_is_bnd_edge,     n_edge);


  PDM_free(edge_vtx);
  PDM_free(edge_face);
  PDM_free(is_solved);
  PDM_free(is_bnd_edge);
  PDM_free(work_data);

  PDM_hash_tab_free(ht);
}
