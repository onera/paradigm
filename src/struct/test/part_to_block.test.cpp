#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_to_block.h"
#include "pdm_logging.h"



MPI_TEST_CASE("[pdm_part_to_block] - 1p - part_to_block",1) {


  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  static const int n_part         = 1;
  static const int n_elmts_part_0 = 7;
  PDM_g_num_t ln_to_gn_p0[n_elmts_part_0] = {2, 2, 5, 5, 1, 2, 1};

  // Convenient array
  int*          pn_elmt      = (int         * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t** pln_to_to_gn = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
  pln_to_to_gn[0] = ln_to_gn_p0;
  pn_elmt     [0] = n_elmts_part_0;

  int pstrid_cst_tab[ 7] = {0, 4         , 2   , 2   , 1, 3         , 2     };
  int pdata_cst_tab [14] = {   1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  int *pdata_cst = pdata_cst_tab;
  int *pstrid    = pstrid_cst_tab;

  SUBCASE("Merge mode") {

    /*
     * Create protocol
     */
    PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_MERGE,
                                                        1.,
                                                        pln_to_to_gn,
                                                        NULL,
                                                        pn_elmt,
                                                        n_part,
                                                        pdm_comm);


    int n_elmt_in_block         = PDM_part_to_block_n_elt_block_get     (ptb);
    PDM_g_num_t* distrib_elmt   = PDM_part_to_block_distrib_index_get   (ptb);
    const PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get      (ptb);
    const int* blk_count        = PDM_part_to_block_block_gnum_count_get(ptb);

    if(0 == 1) {
      log_trace("n_elmt_in_block = %i \n", n_elmt_in_block);
      PDM_log_trace_array_long(distrib_elmt, n_rank+1       , "distrib_elmt : ");
      PDM_log_trace_array_long(blk_gnum    , n_elmt_in_block, "blk_gnum     : ");
      PDM_log_trace_array_int (blk_count   , n_elmt_in_block, "blk_count    : ");
    }

    CHECK(n_elmt_in_block == 3);

    PDM_g_num_t distrib_elmt_expected[2] = {0, 5};
    PDM_g_num_t blk_gnum_expected    [3] = {1, 2, 5};
    int         blk_count_expected   [3] = {2, 3, 2};

    CHECK_EQ_C_ARRAY(distrib_elmt, distrib_elmt_expected, n_rank+1       );
    CHECK_EQ_C_ARRAY(blk_gnum    , blk_gnum_expected    , n_elmt_in_block);
    CHECK_EQ_C_ARRAY(blk_count   , blk_count_expected   , n_elmt_in_block);

    /*
     * Exchange
     */
    int *blk_strid    = NULL;
    int *blk_data_cst = NULL;
    int s_tot_size    = 0;
    int **tmp_preverse_data = NULL;
    int **tmp_preverse_stri = NULL;

    SUBCASE("Synchronous exch ") {

      s_tot_size = PDM_part_to_block_exch(ptb,
                                             sizeof(int),
                                             PDM_STRIDE_VAR_INTERLACED,
                                             -1,
                                             &pstrid,
                                 (void **)   &pdata_cst,
                                             &blk_strid,
                                 (void **)   &blk_data_cst);

      // Reverse
      PDM_part_to_block_reverse_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                                     blk_strid,
                                     blk_data_cst,
                                     &tmp_preverse_stri,
                         (void ***)  &tmp_preverse_data);


    }

    SUBCASE("Asynchronous exch") {
      int request_id = 0;
      PDM_part_to_block_iexch(ptb,
                              PDM_MPI_COMM_KIND_COLLECTIVE,
                              sizeof(int),
                              PDM_STRIDE_VAR_INTERLACED,
                              -1,
                              &pstrid,
                  (void **)   &pdata_cst,
                              &blk_strid,
                  (void **)   &blk_data_cst,
                              &request_id);
      s_tot_size = PDM_part_to_block_iexch_wait(ptb, request_id);


      // Reverse
      PDM_part_to_block_reverse_iexch(ptb,
                                      PDM_MPI_COMM_KIND_COLLECTIVE,
                                      sizeof(int),
                                      PDM_STRIDE_VAR_INTERLACED,
                                     -1,
                                     blk_strid,
                                     blk_data_cst,
                                     &tmp_preverse_stri,
                         (void ***)  &tmp_preverse_data,
                                     &request_id);
      PDM_part_to_block_reverse_iexch_wait(ptb, request_id);


    }

    int *preverse_data = tmp_preverse_data[0];
    free(tmp_preverse_data);

    int *preverse_stri = tmp_preverse_stri[0];
    free(tmp_preverse_stri);

    if(0 == 1) {
      PDM_log_trace_array_int(preverse_stri, n_elmts_part_0, "preverse_stri : ");
      PDM_log_trace_array_int(preverse_data, 35, "preverse_data : ");
    }

    // log_trace("s_tot_size = %i \n", s_tot_size);
    CHECK(s_tot_size == 14);

    // Count buffer size
    int blk_data_size = 0;
    for(int i = 0; i < n_elmt_in_block; ++i) {
      blk_data_size += blk_strid[i];
    }
    // log_trace("blk_data_size = %i \n", blk_data_size);
    CHECK(blk_data_size == 14);

    if(0 == 1) {
      PDM_log_trace_array_int(blk_strid   , n_elmt_in_block , "blk_strid    : ");
      PDM_log_trace_array_int(blk_data_cst,   blk_data_size , "blk_data_cst : ");
    }

    int blk_strid_expected[ 3] = {3, 7, 4};
    int blk_data_expected [14] = {9, 13, 14, 1, 2, 3, 4, 10, 11, 12, 5, 6, 7, 8 };

    CHECK_EQ_C_ARRAY(blk_strid   , blk_strid_expected, n_elmt_in_block);
    CHECK_EQ_C_ARRAY(blk_data_cst, blk_data_expected , blk_data_size  );

    free(blk_strid   );
    free(blk_data_cst);

    // Revser check
    int preverse_stri_expected[7 ] = {7, 7, 4, 4, 3, 7, 3};
    int preverse_data_expected[35] = {1, 2, 3, 4, 10, 11, 12, 1, 2, 3, 4, 10, 11, 12, 5, 6, 7, 8, 5, 6, 7, 8, 9, 13, 14, 1, 2, 3, 4, 10, 11, 12, 9, 13, 14};
    CHECK_EQ_C_ARRAY(preverse_stri, preverse_stri_expected , 7  );
    CHECK_EQ_C_ARRAY(preverse_data, preverse_data_expected , 35 );



    free(preverse_data);
    free(preverse_stri);


    PDM_part_to_block_free(ptb);

  }

  SUBCASE("Cleanup mode") {

    /*
     * Create protocol
     */
    PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                        PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                        1.,
                                                        pln_to_to_gn,
                                                        NULL,
                                                        pn_elmt,
                                                        n_part,
                                                        pdm_comm);

    int n_elmt_in_block         = PDM_part_to_block_n_elt_block_get     (ptb);
    PDM_g_num_t* distrib_elmt   = PDM_part_to_block_distrib_index_get   (ptb);
    const PDM_g_num_t* blk_gnum = PDM_part_to_block_block_gnum_get      (ptb);
    const int* blk_count        = PDM_part_to_block_block_gnum_count_get(ptb);

    if(0 == 1) {
      log_trace("n_elmt_in_block = %i \n", n_elmt_in_block);
      PDM_log_trace_array_long(distrib_elmt, n_rank+1       , "distrib_elmt : ");
      PDM_log_trace_array_long(blk_gnum    , n_elmt_in_block, "blk_gnum     : ");
      PDM_log_trace_array_int (blk_count   , n_elmt_in_block, "blk_count    : ");
    }

    CHECK(n_elmt_in_block == 3);

    PDM_g_num_t distrib_elmt_expected[2] = {0, 5};
    PDM_g_num_t blk_gnum_expected    [3] = {1, 2, 5};
    int         blk_count_expected   [3] = {2, 3, 2};

    CHECK_EQ_C_ARRAY(distrib_elmt, distrib_elmt_expected, n_rank+1       );
    CHECK_EQ_C_ARRAY(blk_gnum    , blk_gnum_expected    , n_elmt_in_block);
    CHECK_EQ_C_ARRAY(blk_count   , blk_count_expected   , n_elmt_in_block);

    /*
     * Exchange
     */
    int *blk_data_cst = NULL;
    int s_tot_size = 0;
    int **tmp_preverse_data = NULL;

    SUBCASE("Synchronous exch ") {
      s_tot_size = PDM_part_to_block_exch(ptb,
                                          sizeof(int),
                                          PDM_STRIDE_CST_INTERLACED,
                                          2,
                                          NULL,
                             (void **)   &pdata_cst,
                                          NULL,
                             (void **)   &blk_data_cst);

      // Reverse
      PDM_part_to_block_reverse_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_CST_INTERLACED,
                                     2,
                                     NULL,
                                     blk_data_cst,
                                     NULL,
                         (void ***)  &tmp_preverse_data);

    }

    SUBCASE("Asynchronous exch") {

      int request_id = 0;
      PDM_part_to_block_iexch(ptb,
                              PDM_MPI_COMM_KIND_COLLECTIVE,
                              sizeof(int),
                              PDM_STRIDE_CST_INTERLACED,
                              2,
                              NULL,
                  (void **)   &pdata_cst,
                              NULL,
                  (void **)   &blk_data_cst,
                              &request_id);

      s_tot_size = PDM_part_to_block_iexch_wait(ptb, request_id);

      // Reverse
      request_id = 0;
      PDM_part_to_block_reverse_iexch(ptb,
                                      PDM_MPI_COMM_KIND_COLLECTIVE,
                                      sizeof(int),
                                      PDM_STRIDE_CST_INTERLACED,
                                      2,
                                      NULL,
                                      blk_data_cst,
                                      NULL,
                          (void ***)  &tmp_preverse_data,
                                      &request_id);
      PDM_part_to_block_reverse_iexch_wait(ptb, request_id);

    }

    CHECK(s_tot_size == 6);

    if(0 == 1) {
      PDM_log_trace_array_int(blk_data_cst, s_tot_size, "blk_data_cst : ");
    }

    int blk_data_expected [6] = {9, 10, 1, 2, 5, 6};

    CHECK_EQ_C_ARRAY(blk_data_cst, blk_data_expected , s_tot_size  );

    int *preverse_data = tmp_preverse_data[0];
    free(tmp_preverse_data);

    if(0 == 1) {
      PDM_log_trace_array_int(preverse_data, 2 * n_elmts_part_0, "preverse_data : ");
    }

    int preverse_data_expected [14] = {1, 2, 1, 2, 5, 6, 5, 6, 9, 10, 1, 2, 9, 10};
    CHECK_EQ_C_ARRAY(preverse_data, preverse_data_expected, s_tot_size);

    free(preverse_data);

    free(blk_data_cst);

    PDM_part_to_block_free(ptb);

  }


  free(pn_elmt);
  free(pln_to_to_gn);

}



MPI_TEST_CASE("[pdm_part_to_block] - 2p - zero global weight", 2) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank (pdm_comm, &i_rank);
  PDM_MPI_Comm_size (pdm_comm, &n_rank);

  int n_part = 1;
  int          *pn_elmt    = (int          *) malloc(sizeof(int          ) * n_part);
  PDM_g_num_t **pln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  double      **pweight    = (double      **) malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **part_data  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t  *block_data = NULL;

  PDM_g_num_t distrib_expected[3];

  SUBCASE("Empty parts") {
    for (int i_part = 0; i_part < n_part; i_part++) {
      pn_elmt  [i_part] = 0;
      pln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * pn_elmt[i_part]);
      pweight  [i_part] = (double      *) malloc(sizeof(double     ) * pn_elmt[i_part]);
    }

    distrib_expected[0] = 0;
    distrib_expected[1] = 0;
    distrib_expected[2] = 0;
  }


  SUBCASE("Non-empty parts") {
    for (int i_part = 0; i_part < n_part; i_part++) {
      pn_elmt  [i_part] = 5;
      pln_to_gn[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * pn_elmt[i_part]);
      pweight  [i_part] = (double      *) malloc(sizeof(double     ) * pn_elmt[i_part]);
      for (int i = 0; i < pn_elmt[i_part]; i++) {
        pln_to_gn[i_part][i] = 1 + n_rank*i + i_rank;
        pweight  [i_part][i] = 0.;
      }
    }

    distrib_expected[0] = 0;
    distrib_expected[1] = 5;
    distrib_expected[2] = 10;
  }

  // Create ptb structure
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      pln_to_gn,
                                                      pweight,
                                                      pn_elmt,
                                                      n_part,
                                                      pdm_comm);

  // Check block distribution
  PDM_g_num_t *distrib = PDM_part_to_block_distrib_index_get(ptb);

  CHECK_EQ_C_ARRAY(distrib, distrib_expected, n_rank+1);

  // Exchange data from part to block
  for (int i_part = 0; i_part < n_part; i_part++) {
    part_data[i_part] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * pn_elmt[i_part]);
    for (int i = 0; i < pn_elmt[i_part]; i++) {
      part_data[i_part][i] = pln_to_gn[i_part][i];
    }
  }

  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
                         (void **) part_data,
                         NULL,
                         (void **) &block_data);

  PDM_part_to_block_free(ptb);

  // Check block data
  int dn_elmt = distrib_expected[i_rank+1] - distrib_expected[i_rank];
  PDM_g_num_t *block_data_expected = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * dn_elmt);
  for (int i = 0; i < dn_elmt; i++) {
    block_data_expected[i] = 1 + distrib_expected[i_rank] + i;
  }

  CHECK_EQ_C_ARRAY(block_data, block_data_expected, dn_elmt);

  // Free memory
  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pln_to_gn[i_part]);
    free(pweight  [i_part]);
    free(part_data[i_part]);
  }
  free(pn_elmt);
  free(pln_to_gn);
  free(pweight);
  free(part_data);
  free(block_data);
  free(block_data_expected);
}
