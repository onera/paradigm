#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_gnum_from_hash_values.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*
 *  Use case
 *
 *       3           4     4           2
 *       +-----------+     +-----------+
 *       |           |     |           |
 *       |           |     |           |
 *       |           |     |           |
 *       +-----------+     +-----------+
 *       5           1     1           6
 *
 *  A l'issu de l'algorithme on doit identifer 7 edges -->
 */


/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);



  int edge_vtx_p0[8] = {1, 4, /* Edge 1 */
                        4, 3,
                        3, 5,
                        5, 1};

  int edge_str_p0[4] = {2, 2, 2, 2};

  int edge_vtx_p1[8] = {2, 4, /* Edge 1 */
                        4, 1,
                        1, 6,
                        6, 2};
  int edge_str_p1[4] = {2, 2, 2, 2};

  /*
   * SetUp
   */
  int     n_part;
  int    *n_elmts;
  int    **part_stri;
  int    **part_data;
  size_t **part_key;
  if(n_rank == 1){
    n_part = 2;
    n_elmts   = (int *     ) malloc( n_part * sizeof(int    ));
    part_stri = (int **    ) malloc( n_part * sizeof(int*   ));
    part_data = (int **    ) malloc( n_part * sizeof(int*   ));
    part_key  = (size_t ** ) malloc( n_part * sizeof(size_t*));

    n_elmts[0] = 4;
    n_elmts[1] = 4;

    part_data[0] = edge_vtx_p0;
    part_data[1] = edge_vtx_p1;

    part_stri[0] = edge_str_p0;
    part_stri[1] = edge_str_p1;

  } else if( n_rank == 2) {
    PDM_printf ("\n Test is not implemented for 2 ranks \n");
    abort();
  }

  /*
   * Compute key
   */
  for(int i_part = 0; i_part < n_part; i_part++){
    part_key[i_part] = (size_t *) malloc(n_elmts[i_part] * sizeof(size_t));
    int idx = 0;
    for(int ielmt = 0; ielmt < n_elmts[i_part]; ++ielmt){
      size_t key = 0;
      for(int idata = 0; idata < part_stri[i_part][ielmt]; ++idata){
        key += part_data[i_part][idx++];
      }
      printf(" part_key[%d][%d] = %lu \n", i_part, ielmt, key);
      part_key[i_part][ielmt] = key;
    }
  }


  PDM_bool_t equilibrate = PDM_FALSE;


  int gnum_fhv_id = PDM_gnum_from_hash_values_create(n_part,
                                                     equilibrate,
                                                     sizeof(int),
                                                     PDM_MPI_COMM_WORLD);

  printf("gnum_fhv_id:: %d \n", gnum_fhv_id);

  for(int i_part = 0; i_part < n_part; ++i_part){
    PDM_gnum_set_hash_values(gnum_fhv_id,
                             i_part,
                             n_elmts[i_part],
                             part_key[i_part],
                             part_stri[i_part],
            (unsigned char*) part_data[i_part]);
  }

  PDM_gnum_from_hv_compute(gnum_fhv_id); /* Passage de part --> block */


  // size_t        *blk_hkeys = NULL;
  // unsigned char *blk_hdata = NULL;
  // int           *blk_hstri = NULL;
  // int            blk_size  = 0;
  // PDM_gnum_from_hv_get_block(&blk_hkeys, &blk_hdata, &blk_hstri);

  /*
   * User responsability to sort the block av
   */
  // int* order = ...


  /*
   * User responsability to sort the block av
   */
  // PDM_gnum_from_hv_set_order(order);

  /*
   *
   */
  // PDM_gnum_from_hv_generate_global_numebering(gnum_fhv_id); /* Passage de part --> block */



  /*
   * Free
   */
  PDM_gnum_from_hv_free(gnum_fhv_id, 0);
  free(part_stri);
  free(part_data);
  free(n_elmts);
  for(int i_part = 0; i_part < n_part; ++i_part){
    free(part_key[i_part]);
  }
  free(part_key);

  PDM_MPI_Finalize ();

  PDM_printf ("\nfin Test\n");

  return 0;

}
