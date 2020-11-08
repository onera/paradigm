#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_multi_block_to_part.h"

MPI_TEST_CASE("[1p] multi_block_to_part",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  static int n_block = 2;
  static PDM_g_num_t multi_distrib_idx[3] = {0, 6, 9};
  static PDM_g_num_t block_distrib1_idx[2] = {0, 6};
  static PDM_g_num_t block_distrib2_idx[2] = {0, 3};

  PDM_g_num_t** block_distrib_idx = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
  block_distrib_idx[0] = block_distrib1_idx;
  block_distrib_idx[1] = block_distrib2_idx;

  static const int n_part         = 1;
  static const int n_elmts_part_0 = 6;
  PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {9, 7, 5, 2, 1, 3};

  // Convenient array
  int*          n_elmts  = (int         * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t** ln_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  ln_to_gn[0] = ln_to_gn_p0;
  n_elmts[0]  = n_elmts_part_0;

  PDM_multi_block_to_part_t* mbtp =
    PDM_multi_block_to_part_create(multi_distrib_idx,
                                   n_block,
          (const PDM_g_num_t**)    block_distrib_idx,
          (const PDM_g_num_t**)    ln_to_gn,
                                   n_elmts,
                                   n_part,
                                   pdm_comm);


  SUBCASE("constant stride = 1") {
    // Stride cst with 2 blk_data reprensenting a implicite block distribution of 6+3 elemnts
    static const int darray_size1 = 6;
    static const int darray_size2 = 3;
    static PDM_g_num_t darray1[darray_size1] = {5 , 4 , 6, 1, 3, 2};
    static PDM_g_num_t darray2[darray_size2] = {10, 11, 8};
    PDM_g_num_t** darray = (PDM_g_num_t **) malloc( n_block * sizeof(PDM_g_num_t*));
    darray[0] = darray1;
    darray[1] = darray2;

    /* Exchange */
    int** stride_one = (int ** ) malloc( n_block * sizeof(int *));
    for(int i_block = 0; i_block < n_block; ++i_block){
      stride_one[i_block] = (int * ) malloc( 1 * sizeof(int));
      stride_one[i_block][0] = 1;
    }

    PDM_g_num_t** parray = NULL;
    PDM_multi_block_to_part_exch2(mbtp, sizeof(PDM_g_num_t), PDM_STRIDE_CST,
                                  stride_one,
                       (void ** ) darray,
                                  NULL,
                       (void ***) &parray);

    static PDM_g_num_t parray_expected_p0[n_elmts_part_0] = {8, 10, 3, 4, 5, 6};

    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
        printf(" parray[%i][%i] = %i \n", i_part, ielmt, parray[i_part][ielmt]);
      }
    }
    CHECK_EQ_C_ARRAY(parray[0], parray_expected_p0, n_elmts_part_0);
    for(int i_part = 0; i_part < n_part; i_part++){
      free(parray[i_part]);
    }
    free(parray);
    for(int i_block = 0; i_block < n_block; ++i_block){
      free(stride_one[i_block]);
    }
    free(stride_one);
    free(darray);
  }





  // Free
  PDM_multi_block_to_part_free(mbtp);
  free(block_distrib_idx);
  free(ln_to_gn);
  free(n_elmts);

}


