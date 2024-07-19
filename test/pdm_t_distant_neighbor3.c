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
#include "pdm_distant_neighbor.h"
#include "pdm_points_merge.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"



/**
 *
 * \brief  Main
 *
 */

int
main
(
int   argc,
char *argv[]
)
{
  int i_rank;
  int n_rank;


  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  assert(n_rank == 4);

  int n_part = 1;

  int *n_elt;
  PDM_malloc(n_elt,n_part,int);
  if (i_rank == 0) {
    n_elt[0] = 2;
  }
  else if (i_rank == 1) {
    n_elt[0] = 2;
  }
  else if (i_rank == 2) {
    n_elt[0] = 3;
  }
  else {
    n_elt[0] = 3;
  }

  int **neighbor_idx;
  PDM_malloc(neighbor_idx,n_part,int *);
  int **neighbor_desc;
  PDM_malloc(neighbor_desc,n_part,int *);

  for (int i = 0; i < n_part; i++) {
    PDM_malloc(neighbor_idx[i],(n_elt[i] + 1),int);
    neighbor_idx[i][0] = 0;
  }

  if (i_rank == 0) {
    neighbor_idx[0][1] = 2;
    neighbor_idx[0][2] = 3;

    PDM_malloc(neighbor_desc[0],neighbor_idx[0][n_elt[0]] * 3,int);

    neighbor_desc[0][3*0+0] = 1;
    neighbor_desc[0][3*0+1] = 0;
    neighbor_desc[0][3*0+2] = 0;

    neighbor_desc[0][3*1+0] = 2;
    neighbor_desc[0][3*1+1] = 0;
    neighbor_desc[0][3*1+2] = 1;

    neighbor_desc[0][3*2+0] = 3;
    neighbor_desc[0][3*2+1] = 0;
    neighbor_desc[0][3*2+2] = 0;
  }

  else if (i_rank == 1) {
    neighbor_idx[0][1] = 2;
    neighbor_idx[0][2] = 3;

    PDM_malloc(neighbor_desc[0],neighbor_idx[0][n_elt[0]] * 3,int);

    neighbor_desc[0][3*0+0] = 0;
    neighbor_desc[0][3*0+1] = 0;
    neighbor_desc[0][3*0+2] = 0;

    neighbor_desc[0][3*1+0] = 2;
    neighbor_desc[0][3*1+1] = 0;
    neighbor_desc[0][3*1+2] = 2;

    neighbor_desc[0][3*2+0] = 3;
    neighbor_desc[0][3*2+1] = 0;
    neighbor_desc[0][3*2+2] = 1;
  }

  else if (i_rank == 2) {
    neighbor_idx[0][1] = 1;
    neighbor_idx[0][2] = 2;
    neighbor_idx[0][3] = 3;

    PDM_malloc(neighbor_desc[0],neighbor_idx[0][n_elt[0]] * 3,int);

    neighbor_desc[0][3*0+0] = 3;
    neighbor_desc[0][3*0+1] = 0;
    neighbor_desc[0][3*0+2] = 2;

    neighbor_desc[0][3*1+0] = 0;
    neighbor_desc[0][3*1+1] = 0;
    neighbor_desc[0][3*1+2] = 0;

    neighbor_desc[0][3*2+0] = 1;
    neighbor_desc[0][3*2+1] = 0;
    neighbor_desc[0][3*2+2] = 0;
  }

  else {
    neighbor_idx[0][1] = 1;
    neighbor_idx[0][2] = 2;
    neighbor_idx[0][3] = 3;

    PDM_malloc(neighbor_desc[0],neighbor_idx[0][n_elt[0]] * 3,int);

    neighbor_desc[0][3*0+0] = 0;
    neighbor_desc[0][3*0+1] = 0;
    neighbor_desc[0][3*0+2] = 1;

    neighbor_desc[0][3*1+0] = 1;
    neighbor_desc[0][3*1+1] = 0;
    neighbor_desc[0][3*1+2] = 1;

    neighbor_desc[0][3*2+0] = 2;
    neighbor_desc[0][3*2+1] = 0;
    neighbor_desc[0][3*2+2] = 0;
  }





  PDM_distant_neighbor_t *dngb = PDM_distant_neighbor_create(comm,
                                                             n_part,
                                                             n_elt,
                                                             neighbor_idx,
                                                             neighbor_desc);



  /* Exchange */
  int **send_n;
  PDM_malloc(send_n,n_part,int *);
  for (int i = 0; i < n_part; i++) {
    PDM_malloc(send_n[i],n_elt[i],int);
    for (int j = 0; j < n_elt[i]; j++) {
      send_n[i][j] = 1;//neighbor_idx[i][j+1] - neighbor_idx[i][j];
    }
  }

  // int **recv_n    = NULL;
  // int **recv_desc = NULL;
  // PDM_distant_neighbor_exch(dngb,
  //                           3*sizeof(int),
  //                           PDM_STRIDE_VAR_INTERLACED,
  //                           -1,
  //                           send_n,
  //                (void  **) neighbor_desc,
  //                           &recv_n,
  //                (void ***) &recv_desc);

  int shift[4];
  shift[0] = 1;
  shift[1] = 3;
  shift[2] = 5;
  shift[3] = 8;
  // if (i_rank == 0) {
  //   shift[0] = 1;
  //   shift[1] = 3;
  // }
  // else {
  //   shift[0] = 5;
  //   shift[1] = 8;
  // }

  int **send_i;
  PDM_malloc(send_i,n_part,int *);
  for (int i = 0; i < n_part; i++) {
    int i_part = i_rank*n_part + i;
    // PDM_malloc(send_i[i],n_elt[i],int);
    PDM_malloc(send_i[i],neighbor_idx[i][n_elt[i]],int);
    for (int j = 0; j < n_elt[i]; j++) {
      // log_trace("part %d (%d), elt %d : %d\n", i, i_part, j, shift[i_part] + j);
      send_i[i][j] = shift[i_part] + j;
      // for (int k = neighbor_idx[i][j]; k < neighbor_idx[i][j+1]; k++) {
      //   send_i[i][k] = shift[i_part] + j;
      // }
      // for (int l = neighbor_idx[i][j]; l < neighbor_idx[i][j+1]; l++) {
      //   int j_rank = neighbor_desc[i][3*l  ];
      //   int j_part = neighbor_desc[i][3*l+1];
      //   int j_num  = neighbor_desc[i][3*l+2];
      //   int g = shift[n_part*j_rank + j_part] + j_num;
      //   // log_trace("  (%d, %d, %d : %d)\n",
      //   //           neighbor_desc[i][3*l  ],
      //   //           neighbor_desc[i][3*l+1],
      //   //           neighbor_desc[i][3*l+2],
      //   //           g);
      // }
    }
  }


  int **recv_n = NULL;
  int **recv_i = NULL;
  PDM_distant_neighbor_exch(dngb,
                            sizeof(int),
                            PDM_STRIDE_VAR_INTERLACED,
                            -1,
                            send_n,
                 (void  **) send_i,
                            &recv_n,
                 (void ***) &recv_i);

  /* Check */
  // for (int i = 0; i < n_part; i++) {

  //   log_trace("\npart #%d\n", i);

  //   PDM_log_trace_array_int(recv_n[i], neighbor_idx[i][n_elt[i]], "recv_n : ");

  //   int idx = 0;
  //   for (int j = 0; j < n_elt[i]; j++) {
  //     log_trace("  elt %d : ", j);
  //     for (int l = neighbor_idx[i][j]; l < neighbor_idx[i][j+1]; l++) {
  //       for (int k = 0; k < recv_n[i][l]; k++) {
  //         int j_rank = neighbor_desc[i][3*l  ];
  //         int j_part = neighbor_desc[i][3*l+1];
  //         int j_num  = neighbor_desc[i][3*l+2];
  //         int g = shift[n_part*j_rank + j_part] + j_num;
  //         log_trace("%d (%d)   ", recv_i[i][idx], g);
  //         idx++;
  //         // log_trace("(%d, %d, %d) ",
  //         //           recv_desc[i][idx+1], recv_desc[i][idx+2], recv_desc[i][idx+3]);
  //         // idx += 3;
  //       }
  //     }
  //     log_trace("\n");
  //   }
  // }


  /* Free memory */
  for (int i = 0; i < n_part; i++) {
   PDM_free(neighbor_idx [i]);
   PDM_free(neighbor_desc[i]);

   PDM_free(send_n[i]);
   PDM_free(send_i[i]);
   PDM_free(recv_n[i]);
   PDM_free(recv_i[i]);
    //PDM_free(recv_desc[i]);
  }
 PDM_free(neighbor_idx );
 PDM_free(neighbor_desc);
 PDM_free(send_n);
 PDM_free(send_i);
 PDM_free(recv_n);
 PDM_free(recv_i);
 PDM_free(n_elt);
  //PDM_free(recv_desc);

  PDM_distant_neighbor_free(dngb);

  PDM_MPI_Finalize();

  return 0;
}
