/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <unistd.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_part_to_block.h"
#include "pdm_part_to_block_priv.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm.h"
#include "pdm_timer.h"
#include "pdm_priv.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_hilbert.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_order.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_size_idx_from_stride.h"

#include "pdm_io.h"
#include <string.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \enum _ptb_timer_step_t
 *
 */

typedef enum {

  MALLOC_ACTIVE_RANKS    = 0, // Initialisation step in Part-to-Block creation
  GENERATE_DISTRIB       = 1, // Block-distribution generation step in Part-to-Block creation
  BINARY_SEARCH          = 2, // Binary search step in Part-to-Block creation
  CREATE_EXCHANGE        = 3, // Collective communication step in Part-to-Block creation
  BLOCK_POST             = 4, // Post-processing step in Part-to-Block creation
  GLOBAL_WEIGHTS         = 5, // Global weight computation step in Part-to-Block creation
  CREATE_FROM_DISTRIB    = 6, // Part-to-Block creation from provided distribution
  CREATE_GEOM            = 7, // Geometric Part-to-Block creation
  DATA_EXCHANGE          = 8  // Collective communication step in Part-to-Block data exchange

} _ptb_timer_step_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

static const double  pdm_part_to_block_distrib_tol = 0.10;

/* Max. number of sub-iterations to get a well-balanced distribution */
static const int pdm_part_to_block_distrib_n_iter_max = 5;

static const int _sampling_factors[4] = {1, /* OD */
                                         2, /* 1D */
                                         2, /* 2D */
                                         4, /* 3D */};

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*
 * Static can cause pb if we're call function in multiple contexte.
 *   For example python and other C++ program
 * No static : truly global
 *  https://stackoverflow.com/questions/1856599/when-to-use-static-keyword-before-global-variables
 */

// Store timers
PDM_timer_t *t_timer[NTIMER_PTB] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

// Timer step by step
double t_elaps[NTIMER_PTB] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
double t_cpu[NTIMER_PTB] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

int min_exch_rank[2] = {INT_MAX, INT_MAX};
int max_exch_rank[2] = {-1, -1};

unsigned long long exch_data[2] = {0, 0};

// Number of Part-to-Block instances in a run
int n_ptb = 0;

// Number of create Part-to-Block instances
int n_ptb_open = 0;

/*=============================================================================
 * Static function definitions
 *============================================================================*/

void
PDM_extents_conformize(int    dim,
                       double extents[],
                       double eps)
{
  double max_range = 1.e-12;
  for (int i = 0; i < dim; i++) {
    max_range = PDM_MAX (max_range, extents[i+dim] - extents[i]);
  }

  // eps = 1.e-3
  const double epsilon = eps * max_range; // Add eps in cas of only one point ...
  for (int i = 0; i < dim; i++) {
    extents[i    ] -= 1.1 * epsilon; // On casse la symetrie !
    extents[i+dim] +=       epsilon;
  }
}

static
void
_counting_sort_long
(
 PDM_part_to_block_t *ptb
)
{
  // On a un pb si la distrib est géant
  int n_elt_block_tot = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
  int *block_n = NULL;
  PDM_malloc(block_n, n_elt_block_tot, int);

  PDM_g_num_t *block_gnum = NULL;
  PDM_malloc(block_gnum, ptb->tn_recv_data, PDM_g_num_t);

  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = ptb->sorted_recv_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    block_n[lid] += 1;
    block_gnum[i] = ptb->sorted_recv_gnum[i];
  }

  int *block_idx = NULL;
  PDM_malloc(block_idx, n_elt_block_tot + 1, int);
  block_idx[0] = 0;
  for(int i = 0; i < n_elt_block_tot; ++i) {
    block_idx[i+1] = block_idx[i] + block_n[i];
    block_n[i] = 0;
  }

  for(int i = 0; i < ptb->tn_recv_data; ++i) {
    int lid = block_gnum[i] - ptb->data_distrib_index[ptb->i_rank] - 1;
    int idx_write = block_idx[lid] + block_n[lid]++;
    ptb->sorted_recv_gnum[idx_write] = block_gnum[i];
    ptb->order[idx_write] = i;
  }

  PDM_free(block_n);
  PDM_free(block_idx);
  PDM_free(block_gnum);

}

/**
 *
 * \brief  Define active ranks
 *
 * \param [inout]   ptb          Part to block structure
 *
 */

static void
_active_ranks
(
 PDM_part_to_block_t  *ptb
)
{

  assert (ptb->active_ranks == NULL);

  if (ptb->s_comm == 1) {
    ptb->is_my_rank_active = 1;
    ptb->n_active_ranks = 1;
    PDM_malloc(ptb->active_ranks, ptb->n_active_ranks, int);
    ptb->active_ranks[0] = ptb->i_rank;
  }

  else {

    switch (ptb->t_distrib) {

    case PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC : {
      ptb->is_my_rank_active = 1;
      ptb->n_active_ranks = ptb->s_comm;
      PDM_malloc(ptb->active_ranks, ptb->n_active_ranks, int);
      for (int i = 0; i < ptb->n_active_ranks; i++) {
        ptb->active_ranks[i] = i;
      }
      break;
    }

    case PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE :
    case PDM_PART_TO_BLOCK_DISTRIB_PART_OF_NODE : {

      ptb->is_my_rank_active = 0;

      int rank_in_node = PDM_io_mpi_node_rank(ptb->comm);
      if (rank_in_node == 0) {
        ptb->is_my_rank_active = 1;
      }

      int *tag_active_ranks;
      PDM_malloc(tag_active_ranks, ptb->s_comm, int);

      PDM_MPI_Allgather((void *) &ptb->is_my_rank_active, 1, PDM_MPI_INT,
                        (void *) tag_active_ranks, 1, PDM_MPI_INT,
                        ptb->comm);

      ptb->n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->n_active_ranks += 1;
        }
      }

      PDM_malloc(ptb->active_ranks, ptb->n_active_ranks, int);
      int n_active_ranks = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        if (tag_active_ranks[i] == 1) {
          ptb->active_ranks[n_active_ranks++] = i;
        }
      }

      if (ptb->t_distrib == PDM_PART_TO_BLOCK_DISTRIB_ONE_PROC_PER_NODE) {
        break;
      }

      int n_node = ptb->n_active_ranks;
      double part_active_node = PDM_MIN (1, ptb->part_active_node);
      part_active_node = PDM_MAX (0, part_active_node);
      ptb->n_active_ranks = (int) floor (n_node * part_active_node);
      ptb->n_active_ranks = PDM_MAX (1, ptb->n_active_ranks);

      n_active_ranks = 0;
      int coeff = n_node / ptb->n_active_ranks;
      for (int i = 0; i < n_node; i++) {
        if (i % coeff == 0) {
          ptb->active_ranks[n_active_ranks++] = ptb->active_ranks[i];
        }
        if (n_active_ranks == ptb->n_active_ranks) {
          break;
        }
      }

      assert (n_active_ranks == ptb->n_active_ranks);

      break;
    }

    default : {
      PDM_error(__FILE__, __LINE__, 0,"Error cs_part_to_bloc : unknown distribute type\n");
      abort();
    }
    }

    /* Dump */

    if (0 == 1) {
      PDM_printf("active ranks : ");
      for(int i = 0; i < ptb->n_active_ranks; i++)
        PDM_printf("%i ", ptb->active_ranks[i]);
      PDM_printf("\n");
    }
  }
}


/**
 *
 * \brief Distrib data
 *
 * \param [inout] ptb              Part to block structure
 * \param [in]    n_totalData      Total number of data
 * \param [out]   data_distrib_index Element global number
 * \param [out]   s_block_min       Local number of elements
 * \param [out]   s_block_max       Number of partition
 *
 */

static void
_distrib_data
(
 PDM_part_to_block_t *ptb,
 int                  user_distrib
)
{
  double t1_elaps = PDM_timer_elapsed(t_timer[GENERATE_DISTRIB]);
  double t1_cpu   = PDM_timer_cpu    (t_timer[GENERATE_DISTRIB]);
  PDM_timer_resume(t_timer[GENERATE_DISTRIB]);

  PDM_g_num_t _id_max     = 0;
  PDM_g_num_t _id_max_max = 0;

  ptb->n_elt_proc= 0;
  for (int i = 0; i < ptb->n_part; i++) {
    ptb->n_elt_proc+= ptb->n_elt[i];
    if(user_distrib == 0) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        PDM_g_num_t gnum = PDM_ABS(ptb->gnum_elt[i][j]);
        _id_max = PDM_MAX (_id_max, gnum);
      }
    }
  }



  if(user_distrib == 0) {

    PDM_MPI_Allreduce (&_id_max,
                       &_id_max_max,
                       1,
                       PDM__PDM_MPI_G_NUM,
                       PDM_MPI_MAX,
                       ptb->comm);

    if (ptb->weight == NULL) {
      PDM_g_num_t _n_rank_data = _id_max_max / ptb->n_active_ranks;
      PDM_g_num_t _rest = _id_max_max % ptb->n_active_ranks;

      int n_rank_data = (int) (_n_rank_data);
      int rest       = (int) (_rest);

      ptb->s_block_max = n_rank_data;
      ptb->s_block_min = n_rank_data;

      if (rest != 0) {
        ptb->s_block_max += 1;
      }

      for (int i = 0; i < ptb->s_comm + 1; i++) {
        ptb->data_distrib_index[i] = 0;
      }

      int k = 0;
      int idx = 0;
      for (int i = 0; i < ptb->s_comm; i++) {
        ptb->data_distrib_index[i+1] +=  ptb->data_distrib_index[i];
        if (idx < ptb->n_active_ranks) {
          if (ptb->active_ranks[idx] == i) {
            ptb->data_distrib_index[i+1] += n_rank_data;
            if (k < rest)
              ptb->data_distrib_index[i+1] += 1;
            k += 1;
            idx++;
          }

        }
      }

    }

    else {
      const int dim = 2;
      const int  n_active_ranks = ptb->n_active_ranks;
      const int  sampling_factor = _sampling_factors[dim];

      double **weight = ptb->weight;
      int *_active_ranks = ptb->active_ranks;

      PDM_g_num_t *rank_index = NULL;
      PDM_distrib_weight(             sampling_factor,
                                      n_active_ranks,
                                      ptb->n_part,
                                      ptb->n_elt,
              (const PDM_g_num_t **)  ptb->gnum_elt,
              (const double      **)  weight,
                                      pdm_part_to_block_distrib_n_iter_max,
                                      pdm_part_to_block_distrib_tol,
                                      ptb->comm,
                                      &rank_index);

      if(0 == 1) {
        PDM_log_trace_array_long(rank_index, n_active_ranks, "rank_index :: ");
      }

      ptb->data_distrib_index[0] = 0;
      int k = 0;
      for (int i = 0; i < n_active_ranks; i++) {
        int i_activeRank = _active_ranks[i];
        while (k < i_activeRank) {
          ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
          k++;
        }
        ptb->data_distrib_index[k+1] = rank_index[i+1];
        k++;
      }
      while (k < ptb->s_comm) {
        ptb->data_distrib_index[k+1] = ptb->data_distrib_index[k];
        k++;
      }

      PDM_free(rank_index);

    }
  } // If User

  /* Affichage */
  if (1 == 0) {
    if (ptb->i_rank == 0) {
      PDM_printf("data_distrib_index : ");
      for(int i = 0; i < ptb->s_comm + 1; i++)
        PDM_printf(PDM_FMT_G_NUM" ", ptb->data_distrib_index[i]);
      PDM_printf("\n");
    }
  }

  PDM_timer_hang_on(t_timer[GENERATE_DISTRIB]);
  double t2_elaps = PDM_timer_elapsed(t_timer[GENERATE_DISTRIB]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[GENERATE_DISTRIB]);

  t_elaps[GENERATE_DISTRIB] += (t2_elaps - t1_elaps);
  t_cpu  [GENERATE_DISTRIB] += (t2_cpu   - t1_cpu  );

  double t3_elaps = PDM_timer_elapsed(t_timer[BINARY_SEARCH]);
  double t3_cpu   = PDM_timer_cpu    (t_timer[BINARY_SEARCH]);
  PDM_timer_resume(t_timer[BINARY_SEARCH]);

  PDM_malloc(ptb->n_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_data, ptb->s_comm, int);

  /* Pour chaque donnee le proc ou elle va etre envoyee */
  PDM_malloc(ptb->dest_proc, ptb->n_elt_proc, int);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    for (int j = 0; j < ptb->n_elt[i]; j++) {

      PDM_g_num_t _gnum_elt = PDM_ABS(ptb->gnum_elt[i][j]) - 1;

      int iproc = PDM_binary_search_gap_long (_gnum_elt,
                                              ptb->data_distrib_index,
                                              ptb->s_comm + 1);

      ptb->dest_proc[++idx] = iproc;
      // assert (ptb->dest_proc[idx] >= 0);
      ptb->n_send_data[iproc] += 1;
    }
  }

  PDM_timer_hang_on(t_timer[BINARY_SEARCH]);
  double t4_elaps = PDM_timer_elapsed(t_timer[BINARY_SEARCH]);
  double t4_cpu   = PDM_timer_cpu    (t_timer[BINARY_SEARCH]);

  t_elaps[BINARY_SEARCH] += (t4_elaps - t3_elaps);
  t_cpu  [BINARY_SEARCH] += (t4_cpu   - t3_cpu  );

  double t5_elaps = PDM_timer_elapsed(t_timer[CREATE_EXCHANGE]);
  double t5_cpu   = PDM_timer_cpu    (t_timer[CREATE_EXCHANGE]);
  PDM_timer_resume(t_timer[CREATE_EXCHANGE]);

  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);


  PDM_malloc(ptb->i_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_data, ptb->s_comm, int);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];

  PDM_g_num_t *send_gnum = NULL;
  PDM_malloc(send_gnum, ptb->tn_send_data, PDM_g_num_t);

  for (int i = 0; i < ptb->s_comm; i++)
    ptb->n_send_data[i] = 0;

  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[idx];
      send_gnum[ptb->i_send_data[iproc] +
                ptb->n_send_data[iproc]] = PDM_ABS(ptb->gnum_elt[i][j]);
      idx++;
      ptb->n_send_data[iproc] += 1;
    }
  }

  PDM_malloc(ptb->sorted_recv_gnum, ptb->tn_recv_data, PDM_g_num_t);

  PDM_MPI_Partofactiverank (ptb->n_send_data,
                            ptb->n_recv_data,
                            ptb->comm,
                            &(ptb->part_active_rank));

  if (ptb->p2p_factor < ptb->part_active_rank) {

    PDM_MPI_Alltoallv(send_gnum,
                      ptb->n_send_data,
                      ptb->i_send_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->sorted_recv_gnum,
                      ptb->n_recv_data,
                      ptb->i_recv_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->comm);
  }

  else {

    PDM_MPI_Alltoallv_p2p(send_gnum,
                          ptb->n_send_data,
                          ptb->i_send_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->sorted_recv_gnum,
                          ptb->n_recv_data,
                          ptb->i_recv_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->comm);

  }

  PDM_free(send_gnum);

  PDM_timer_hang_on(t_timer[CREATE_EXCHANGE]);
  double t6_elaps = PDM_timer_elapsed(t_timer[CREATE_EXCHANGE]);
  double t6_cpu   = PDM_timer_cpu    (t_timer[CREATE_EXCHANGE]);

  t_elaps[CREATE_EXCHANGE] += (t6_elaps - t5_elaps);
  t_cpu  [CREATE_EXCHANGE] += (t6_cpu   - t5_cpu  );

  double t7_elaps = PDM_timer_elapsed(t_timer[BLOCK_POST]);
  double t7_cpu   = PDM_timer_cpu    (t_timer[BLOCK_POST]);
  PDM_timer_resume(t_timer[BLOCK_POST]);

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    PDM_malloc(ptb->order, ptb->tn_recv_data, int);

    char *env_var = NULL;
    int counting_sort = 0;
    env_var = getenv ("PDM_PTB_COUNTING_SORT");
    if (env_var != NULL) {
      counting_sort = atoi(env_var);
    }

    if(counting_sort == 0) {
      for (int i = 0; i < ptb->tn_recv_data; i++) {
        ptb->order[i] = i;
      }

      PDM_sort_long (ptb->sorted_recv_gnum,
                     ptb->order,
                     ptb->tn_recv_data);
    } else {
      if(ptb->i_rank == 0) {
        printf("_counting_sort_long \n");
      }
      _counting_sort_long(ptb);
    }

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }

  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    ptb->n_elt_block = ptb->data_distrib_index[ptb->i_rank+1] - ptb->data_distrib_index[ptb->i_rank];
    // printf("ptb->n_elt_block::%d \n", ptb->n_elt_block);
    // printf(" ptb->data_distrib_index[ptb->i_rank+1]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank+1]);
    // printf(" ptb->data_distrib_index[ptb->i_rank]::"PDM_FMT_G_NUM" \n", ptb->data_distrib_index[ptb->i_rank]);
  }

  /*
   * Cleanup
   */

  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    PDM_malloc(ptb->block_gnum, ptb->tn_recv_data, PDM_g_num_t);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
      else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    PDM_realloc(ptb->block_gnum ,ptb->block_gnum , ptb->n_elt_block,PDM_g_num_t);

    PDM_malloc(ptb->block_gnum_count, ptb->n_elt_block, int);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  }

  else {

    ptb->block_gnum = ptb->sorted_recv_gnum;

    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  /*
   * Create idx_partial for reverse
   */
  ptb->enable_reverse = 1;
  if ( (ptb->t_post         != PDM_PART_TO_BLOCK_POST_NOTHING)       &&
       (ptb->t_post         != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) &&
       (ptb->enable_reverse == 1)) {
    PDM_malloc(ptb->idx_partial, ptb->tn_recv_data, int);
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      PDM_g_num_t gnum = ptb->sorted_recv_gnum[i];
      int idx_in_partial_block = PDM_binary_search_long(gnum, ptb->block_gnum, ptb->n_elt_block);
      assert(idx_in_partial_block != -1);
      ptb->idx_partial[i] = idx_in_partial_block;
    }
  }

  PDM_timer_hang_on(t_timer[BLOCK_POST]);
  double t8_elaps = PDM_timer_elapsed(t_timer[BLOCK_POST]);
  double t8_cpu   = PDM_timer_cpu    (t_timer[BLOCK_POST]);

  t_elaps[BLOCK_POST] += (t8_elaps - t7_elaps);
  t_cpu  [BLOCK_POST] += (t8_cpu   - t7_cpu  );
}

static
void
_distrib_data_hilbert
(
  PDM_part_to_block_t   *ptb,
  double               **pvtx_coords,
  double               **weight
)
{
  ptb->n_elt_proc = 0;
  int *part_idx;
  PDM_malloc(part_idx, ptb->n_part + 1, int);
  part_idx[0] = 0;
  for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
    ptb->n_elt_proc += ptb->n_elt[i_part];
    part_idx[i_part+1] = part_idx[i_part] + ptb->n_elt[i_part];
  }

  double *concat_vtx_coord = NULL;
  double *concat_weight    = NULL;
  if(ptb->n_part == 1 ) {
    concat_vtx_coord = pvtx_coords[0];
    concat_weight    = weight     [0];
  } else if (ptb->n_part > 1) {
    PDM_malloc(concat_vtx_coord, 3 * ptb->n_elt_proc, double);
    PDM_malloc(concat_weight   ,     ptb->n_elt_proc, double);

    int shift = 0;
    for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        int idx_write = 3*(shift + i_elt);
        concat_vtx_coord[idx_write    ] = pvtx_coords[i_part][3 * i_elt    ];
        concat_vtx_coord[idx_write + 1] = pvtx_coords[i_part][3 * i_elt + 1];
        concat_vtx_coord[idx_write + 2] = pvtx_coords[i_part][3 * i_elt + 2];
      }
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        concat_weight[shift+i_elt] = weight[i_part][i_elt];
      }
      shift += ptb->n_elt[i_part];
    }
  }

  /** Initialisation **/
  int dim = 3;
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_extents_conformize(dim, extents, 1e-3);

  /** Hilbert Coordinates Computation **/
  PDM_hilbert_code_t *hilbert_codes;
  PDM_malloc(hilbert_codes, ptb->n_elt_proc, PDM_hilbert_code_t);

  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, ptb->n_elt_proc, concat_vtx_coord, hilbert_codes);

  if(ptb->n_part > 1 ) {
    PDM_free(concat_vtx_coord);
  }

  PDM_hilbert_code_t *hilbert_codes_idx;
  PDM_malloc(hilbert_codes_idx, ptb->s_comm+1, PDM_hilbert_code_t);

  PDM_hilbert_build_rank_index(3,
                               ptb->s_comm,  // Number of chunk
                               ptb->n_elt_proc,
                               hilbert_codes,
                               concat_weight,
                               NULL,
                               hilbert_codes_idx, // Is the distrib
                               ptb->comm);


  if(0 == 1) {
    PDM_log_trace_array_double(hilbert_codes, ptb->n_elt_proc, "hilbert_codes :: ");
    PDM_log_trace_array_double(hilbert_codes_idx, ptb->s_comm+1, "hilbert_codes_idx :: ");
  }

  if(ptb->n_part > 1 ) {
    PDM_free(concat_weight);
  }

  PDM_malloc(ptb->n_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_data, ptb->s_comm, int);

  /* Pour chaque donnee le proc ou elle va etre envoyee */
  PDM_malloc(ptb->dest_proc, ptb->n_elt_proc, int);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */
  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  /*
   * Caution :  We need to send in the ordering of partition !!!!
   */
  int idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = hilbert_order[i_concat_elt]; // Donc dans la frame de depart

      size_t t_rank = PDM_hilbert_quantile_search(ptb->s_comm,
                                                  hilbert_codes[i_concat_elt],
                                                  hilbert_codes_idx);
      ptb->dest_proc  [idx++ ]  = t_rank;
      ptb->n_send_data[t_rank] += 1;
    }
  }

  PDM_free(hilbert_codes_idx);
  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  PDM_malloc(ptb->i_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_data, ptb->s_comm, int);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];

  PDM_g_num_t        *send_gnum  = NULL;
  PDM_hilbert_code_t *send_codes = NULL;
  PDM_malloc(send_gnum , ptb->tn_send_data, PDM_g_num_t       );
  PDM_malloc(send_codes, ptb->tn_send_data, PDM_hilbert_code_t);

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  // Dans le cas géométrique on doit envoyer le gnum + le codes
  idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {

      int t_rank         = ptb->dest_proc[idx++];
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = hilbert_order[i_concat_elt]; // Donc dans la frame de depart

      send_gnum [ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = PDM_ABS(ptb->gnum_elt[i_part][j]);
      send_codes[ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = hilbert_codes[i_concat_elt]; // Not sure
      ptb->n_send_data[t_rank] += 1;
    }
  }
  //PDM_free(hilbert_order);
  PDM_free(hilbert_codes);

  PDM_hilbert_code_t *sorted_recv_codes;
  PDM_malloc(ptb->sorted_recv_gnum, ptb->tn_recv_data, PDM_g_num_t       );
  PDM_malloc(sorted_recv_codes    , ptb->tn_recv_data, PDM_hilbert_code_t);

  PDM_MPI_Partofactiverank (ptb->n_send_data,
                            ptb->n_recv_data,
                            ptb->comm,
                            &(ptb->part_active_rank));

  if (ptb->p2p_factor < ptb->part_active_rank) {
    PDM_MPI_Alltoallv(send_gnum,
                      ptb->n_send_data,
                      ptb->i_send_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->sorted_recv_gnum,
                      ptb->n_recv_data,
                      ptb->i_recv_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->comm);
  
    PDM_MPI_Alltoallv(send_codes,
                      ptb->n_send_data,
                      ptb->i_send_data,
                      PDM_MPI_DOUBLE,
                      sorted_recv_codes,
                      ptb->n_recv_data,
                      ptb->i_recv_data,
                      PDM_MPI_DOUBLE,
                      ptb->comm);
  }

  else {
    PDM_MPI_Alltoallv_p2p(send_gnum,
                          ptb->n_send_data,
                          ptb->i_send_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->sorted_recv_gnum,
                          ptb->n_recv_data,
                          ptb->i_recv_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->comm);
  
    PDM_MPI_Alltoallv_p2p(send_codes,
                          ptb->n_send_data,
                          ptb->i_send_data,
                          PDM_MPI_DOUBLE,
                          sorted_recv_codes,
                          ptb->n_recv_data,
                          ptb->i_recv_data,
                          PDM_MPI_DOUBLE,
                          ptb->comm);
  }

  PDM_free(send_gnum);
  PDM_free(send_codes);
  PDM_free(part_idx);

  if(0 == 1) {
    PDM_log_trace_array_long  (ptb->sorted_recv_gnum, ptb->tn_recv_data, "ptb->sorted_recv_gnum :: ");
    PDM_log_trace_array_double(sorted_recv_codes    , ptb->tn_recv_data, "sorted_recv_codes :: ");
  }

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    PDM_malloc(ptb->order, ptb->tn_recv_data, int);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      ptb->order[i] = i;
    }

    PDM_hilbert_local_order(ptb->tn_recv_data, sorted_recv_codes, ptb->order);
    PDM_order_array(ptb->tn_recv_data, sizeof(PDM_g_num_t), ptb->order, ptb->sorted_recv_gnum);

    // PDM_log_trace_array_double(sorted_recv_codes    , ptb->tn_recv_data, "sorted_recv_codes :: ");

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }

  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) {
    PDM_error(__FILE__, __LINE__, 0,"Error PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM : not implemented \n");
  }

  /*
   * Cleanup
   */
  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    PDM_malloc(ptb->block_gnum, ptb->tn_recv_data, PDM_g_num_t);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      } else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    PDM_realloc(ptb->block_gnum ,ptb->block_gnum , ptb->n_elt_block,PDM_g_num_t);

    PDM_malloc(ptb->block_gnum_count, ptb->n_elt_block, int);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  } else {
    ptb->block_gnum       = ptb->sorted_recv_gnum;
    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  // Generate distribution
  PDM_distrib_compute(ptb->n_elt_block, ptb->data_distrib_index, -1, ptb->comm);

  PDM_free(sorted_recv_codes);

  /*
   * To do : ptb->enable_reverse = 1;
   */
  ptb->enable_reverse = 0;


}


static
void
_distrib_data_morton
(
  PDM_part_to_block_t   *ptb,
  double               **pvtx_coords,
  double               **weight
)
{
  ptb->n_elt_proc = 0;
  int *part_idx = NULL;
  PDM_malloc(part_idx, ptb->n_part + 1, int);
  part_idx[0] = 0;
  for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
    ptb->n_elt_proc += ptb->n_elt[i_part];
    part_idx[i_part+1] = part_idx[i_part] + ptb->n_elt[i_part];
  }

  double *concat_vtx_coord = NULL;
  double *concat_weight    = NULL;
  if(ptb->n_part == 1 ) {
    concat_vtx_coord = pvtx_coords[0];
    concat_weight    = weight     [0];
  } else {
    PDM_malloc(concat_vtx_coord, 3 * ptb->n_elt_proc, double);
    PDM_malloc(concat_weight   ,     ptb->n_elt_proc, double);

    int shift = 0;
    for(int i_part = 0; i_part < ptb->n_part; ++i_part) {
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        int idx_write = 3*(shift + i_elt);
        concat_vtx_coord[idx_write    ] = pvtx_coords[i_part][3 * i_elt    ];
        concat_vtx_coord[idx_write + 1] = pvtx_coords[i_part][3 * i_elt + 1];
        concat_vtx_coord[idx_write + 2] = pvtx_coords[i_part][3 * i_elt + 2];
      }
      for(int i_elt = 0; i_elt < ptb->n_elt[i_part]; ++i_elt) {
        concat_weight[shift+i_elt] = weight[i_part][i_elt];
      }
      shift += ptb->n_elt[i_part];
    }
  }


  /** Initialisation **/
  int dim = 3;
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  // PDM_morton_get_coord_extents_par(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_morton_get_coord_extents(dim, ptb->n_elt_proc, concat_vtx_coord, extents, ptb->comm);
  PDM_extents_conformize(dim, extents, 1e-3);

  /** morton Coordinates Computation **/
  PDM_morton_code_t *morton_codes = NULL;
  PDM_malloc(morton_codes, ptb->n_elt_proc, PDM_morton_code_t);

  const PDM_morton_int_t max_level = PDM_morton_max_level;
  double d[3];
  double s[3];
  PDM_morton_encode_coords(dim, max_level, extents, ptb->n_elt_proc, concat_vtx_coord, morton_codes, d, s);

  if(ptb->n_part != 1) {
    PDM_free(concat_vtx_coord);
  }

  PDM_morton_code_t *morton_codes_idx = NULL;
  PDM_malloc(morton_codes_idx, ptb->s_comm + 1, PDM_morton_code_t);
  PDM_morton_build_rank_index(3,
                              max_level, // Number of chunk
                              ptb->n_elt_proc,
                              morton_codes,
                              concat_weight,
                              NULL,
                              morton_codes_idx, // Is the distrib
                              ptb->comm);

  if(ptb->n_part != 1) {
    PDM_free(concat_weight);
  }

  PDM_malloc(ptb->n_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_data, ptb->s_comm, int);

  /* Pour chaque donnee le proc ou elle va etre envoyee */
  PDM_malloc(ptb->dest_proc, ptb->n_elt_proc, int);

  /* Calcul du nombre de donnees a envoyer a chaque procesus */
  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  /*
   * Caution :  We need to send in the ordering of partition !!!!
   */
  int idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = morton_order[i_concat_elt]; // Donc dans la frame de depart

      size_t t_rank = PDM_morton_quantile_search(ptb->s_comm,
                                                 morton_codes[i_concat_elt],
                                                 morton_codes_idx);
      ptb->dest_proc  [idx++ ]  = t_rank;
      ptb->n_send_data[t_rank] += 1;
    }
  }

  PDM_free(morton_codes_idx);

  PDM_MPI_Alltoall (ptb->n_send_data, 1, PDM_MPI_INT,
                    ptb->n_recv_data, 1, PDM_MPI_INT,
                    ptb->comm);

  PDM_malloc(ptb->i_send_data, ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_data, ptb->s_comm, int);

  ptb->i_send_data[0] = 0;
  ptb->i_recv_data[0] = 0;
  for (int i = 1; i < ptb->s_comm; i++) {
    ptb->i_send_data[i] = ptb->i_send_data[i-1] + ptb->n_send_data[i-1];
    ptb->i_recv_data[i] = ptb->i_recv_data[i-1] + ptb->n_recv_data[i-1];
  }

  ptb->tn_recv_data = ptb->i_recv_data[ptb->s_comm - 1] +
                      ptb->n_recv_data[ptb->s_comm - 1];

  ptb->tn_send_data = ptb->i_send_data[ptb->s_comm - 1] +
                      ptb->n_send_data[ptb->s_comm - 1];


  PDM_g_num_t       *send_gnum  = NULL;
  PDM_morton_code_t *send_codes = NULL;
  PDM_malloc(send_gnum , ptb->tn_send_data, PDM_g_num_t      );
  PDM_malloc(send_codes, ptb->tn_send_data, PDM_morton_code_t);

  for (int i = 0; i < ptb->s_comm; i++) {
    ptb->n_send_data[i] = 0;
  }

  // Dans le cas géométrique on doit envoyer le gnum + le codes
  idx = 0;
  for (int i_part = 0; i_part < ptb->n_part; i_part++) {
    for (int j = 0; j < ptb->n_elt[i_part]; j++) {

      int t_rank         = ptb->dest_proc[idx++];
      int i_concat_elt   = j + part_idx[i_part];
      // int old_concat_elt = morton_order[i_concat_elt]; // Donc dans la frame de depart

      send_gnum [ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]] = PDM_ABS(ptb->gnum_elt[i_part][j]);
      PDM_morton_copy(morton_codes[i_concat_elt], &send_codes[ptb->i_send_data[t_rank] + ptb->n_send_data[t_rank]]);

      ptb->n_send_data[t_rank] += 1;
    }
  }
  //PDM_free(morton_order);
  PDM_free(morton_codes);

  PDM_morton_code_t *sorted_recv_codes = NULL;
  PDM_malloc(ptb->sorted_recv_gnum, ptb->tn_recv_data, PDM_g_num_t      );
  PDM_malloc(sorted_recv_codes    , ptb->tn_recv_data, PDM_morton_code_t);

  PDM_MPI_Partofactiverank (ptb->n_send_data,
                            ptb->n_recv_data,
                            ptb->comm,
                            &(ptb->part_active_rank));

  if (ptb->p2p_factor < ptb->part_active_rank) {

    PDM_MPI_Alltoallv(send_gnum,
                      ptb->n_send_data,
                      ptb->i_send_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->sorted_recv_gnum,
                      ptb->n_recv_data,
                      ptb->i_recv_data,
                      PDM__PDM_MPI_G_NUM,
                      ptb->comm);
  
    PDM_MPI_Datatype mpi_morton_type;
    PDM_MPI_Type_create_contiguous(4, PDM_MPI_INT, &mpi_morton_type);
    PDM_MPI_Type_commit(&mpi_morton_type);
    PDM_MPI_Alltoallv(send_codes,
                      ptb->n_send_data,
                      ptb->i_send_data,
                      mpi_morton_type,
                      sorted_recv_codes,
                      ptb->n_recv_data,
                      ptb->i_recv_data,
                      mpi_morton_type,
                      ptb->comm);
    PDM_MPI_Type_free(&mpi_morton_type);
  }

  else {

    PDM_MPI_Alltoallv_p2p(send_gnum,
                          ptb->n_send_data,
                          ptb->i_send_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->sorted_recv_gnum,
                          ptb->n_recv_data,
                          ptb->i_recv_data,
                          PDM__PDM_MPI_G_NUM,
                          ptb->comm);
  
    PDM_MPI_Datatype mpi_morton_type;
    PDM_MPI_Type_create_contiguous(4, PDM_MPI_INT, &mpi_morton_type);
    PDM_MPI_Type_commit(&mpi_morton_type);
    PDM_MPI_Alltoallv_p2p(send_codes,
                          ptb->n_send_data,
                          ptb->i_send_data,
                          mpi_morton_type,
                          sorted_recv_codes,
                          ptb->n_recv_data,
                          ptb->i_recv_data,
                          mpi_morton_type,
                          ptb->comm);
    PDM_MPI_Type_free(&mpi_morton_type);
  }

  
  PDM_free(send_gnum);
  PDM_free(send_codes);
  PDM_free(part_idx);


  if(0 == 1) {
    PDM_log_trace_array_long  (ptb->sorted_recv_gnum, ptb->tn_recv_data, "ptb->sorted_recv_gnum :: ");
  }

  /*
   * Sort
   */
  if( ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM)
  {
    PDM_malloc(ptb->order, ptb->tn_recv_data, int);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      ptb->order[i] = i;
    }

    PDM_morton_local_order(ptb->tn_recv_data, sorted_recv_codes, ptb->order);
    PDM_order_array(ptb->tn_recv_data, sizeof(PDM_g_num_t      ), ptb->order, ptb->sorted_recv_gnum);
    // PDM_order_array(ptb->tn_recv_data, sizeof(PDM_morton_code_t), ptb->order, ptb->sorted_recv_gnum);

    if(0 == 1) {
      PDM_log_trace_array_long(ptb->sorted_recv_gnum, ptb->tn_recv_data, "sorted_recv_gnum");
      PDM_log_trace_array_int (ptb->order           , ptb->tn_recv_data, "order");
    }
  }

  ptb->n_elt_block = ptb->tn_recv_data;

  if( ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) {
    PDM_error(__FILE__, __LINE__, 0,"Error PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM : not implemented \n");
  }

  /*
   * Cleanup
   */
  if ( (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) &&
       (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) ) {

    int n_elt_block = 0;

    PDM_malloc(ptb->block_gnum, ptb->tn_recv_data, PDM_g_num_t);

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      if (i == 0) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      } else if (ptb->block_gnum[n_elt_block-1] != ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum[n_elt_block++] = ptb->sorted_recv_gnum[i];
      }
    }
    ptb->n_elt_block = n_elt_block;
    PDM_realloc(ptb->block_gnum ,ptb->block_gnum , ptb->n_elt_block,PDM_g_num_t);

    PDM_malloc(ptb->block_gnum_count, ptb->n_elt_block, int);
    n_elt_block = 0;
    if (ptb->tn_recv_data > 0)
      ptb->block_gnum_count[0] = 1;
    for (int i = 1; i < ptb->tn_recv_data; i++) {
      //If same than previous, juste increase stride counter
      if (ptb->sorted_recv_gnum[i-1] == ptb->sorted_recv_gnum[i]) {
        ptb->block_gnum_count[n_elt_block]++;
      }
      else { //Otherwise, count next
        ptb->block_gnum_count[++n_elt_block] = 1;
      }
    }
  } else {
    ptb->block_gnum       = ptb->sorted_recv_gnum;
    ptb->block_gnum_count = PDM_array_const_int(ptb->tn_recv_data, 1);
  }

  // Generate distribution
  PDM_distrib_compute(ptb->n_elt_block, ptb->data_distrib_index, -1, ptb->comm);

  PDM_free(sorted_recv_codes);

  /*
   * To do : ptb->enable_reverse = 1;
   */
  ptb->enable_reverse = 0;


}


static void
_compute_global_weights
(
 PDM_part_to_block_t *ptb
 )
{
  /* Send local weights */
  int *send_count = PDM_array_zeros_int (ptb->s_comm);
  double *part_weight = NULL;
  PDM_malloc(part_weight, ptb->tn_send_data, double);
  int idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      part_weight[ptb->i_send_data[rank] + send_count[rank]++] = ptb->weight[i][j];
    }
  }

  double *recv_weight = NULL;
  PDM_malloc(recv_weight, ptb->tn_recv_data, double);

  if (ptb->p2p_factor < ptb->part_active_rank) {
    PDM_MPI_Alltoallv (part_weight,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_DOUBLE,
                       recv_weight,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_DOUBLE,
                       ptb->comm);
  }
  else {
    PDM_MPI_Alltoallv_p2p (part_weight,
                           ptb->n_send_data,
                           ptb->i_send_data,
                           PDM_MPI_DOUBLE,
                           recv_weight,
                           ptb->n_recv_data,
                           ptb->i_recv_data,
                           PDM_MPI_DOUBLE,
                           ptb->comm);
  }

  /* Sum received weights */
  double *block_weight = NULL;
  PDM_malloc(block_weight, ptb->n_elt_block, double);
  for (int i = 0; i < ptb->n_elt_block; i++) {
    block_weight[i] = 0.;
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    block_weight[idx] += recv_weight[j];
  }

  idx = 0;
  for (int i = 0; i < ptb->tn_recv_data; i++) {
    int j = ptb->order[i];
    while (ptb->block_gnum[idx] < ptb->sorted_recv_gnum[i]) {
      idx++;
    }
    recv_weight[j] = block_weight[idx];
  }
  PDM_free(block_weight);

  /* Send back global weights */

  if (ptb->p2p_factor < ptb->part_active_rank) {

    PDM_MPI_Alltoallv (recv_weight,
                       ptb->n_recv_data,
                       ptb->i_recv_data,
                       PDM_MPI_DOUBLE,
                       part_weight,
                       ptb->n_send_data,
                       ptb->i_send_data,
                       PDM_MPI_DOUBLE,
                       ptb->comm);
  }

  else {

    PDM_MPI_Alltoallv_p2p (recv_weight,
                           ptb->n_recv_data,
                           ptb->i_recv_data,
                           PDM_MPI_DOUBLE,
                           part_weight,
                           ptb->n_send_data,
                           ptb->i_send_data,
                           PDM_MPI_DOUBLE,
                           ptb->comm);
  }

  PDM_free(recv_weight);

  /* Store global weights */
  PDM_malloc(ptb->weight_g, ptb->n_part, double *);
  PDM_array_reset_int (send_count, ptb->s_comm, 0);
  idx = 0;
  for (int i = 0; i < ptb->n_part; i++) {
    PDM_malloc(ptb->weight_g[i], ptb->n_elt[i], double);
    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int rank = ptb->dest_proc[idx++];
      ptb->weight_g[i][j] = part_weight[ptb->i_send_data[rank] + send_count[rank]++];
    }
  }
  PDM_free(send_count);
  PDM_free(part_weight);
}

/**
 *
 * \brief  Define active ranks
 *
 * \param [inout]   ptb          Part to block structure
 *
 */

static
PDM_part_to_block_t *
_ptb_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  if (n_ptb == 0) {
    t_timer[MALLOC_ACTIVE_RANKS] = PDM_timer_create (); // Warning : unused for now because negligable
    t_timer[GENERATE_DISTRIB   ] = PDM_timer_create ();
    t_timer[BINARY_SEARCH      ] = PDM_timer_create ();
    t_timer[CREATE_EXCHANGE    ] = PDM_timer_create ();
    t_timer[BLOCK_POST         ] = PDM_timer_create ();
    t_timer[GLOBAL_WEIGHTS     ] = PDM_timer_create ();
    t_timer[DATA_EXCHANGE      ] = PDM_timer_create ();
    t_timer[CREATE_FROM_DISTRIB] = PDM_timer_create ();
    t_timer[CREATE_GEOM        ] = PDM_timer_create ();
  }
  n_ptb++;
  n_ptb_open++;

  PDM_part_to_block_t *ptb;
  PDM_malloc(ptb, 1, PDM_part_to_block_t);

  ptb->t_distrib         = t_distrib;    /*!< Distribution type */
  ptb->t_post            = t_post;       /*!< Post processing type */
  ptb->n_active_ranks    = 0;            /*!< Number of active ranks */
  ptb->active_ranks      = NULL;         /*!< List of active ranks */
  ptb->comm              = comm;         /*!< MSG communicator */
  PDM_MPI_Comm_size (comm, &(ptb->s_comm));
  PDM_MPI_Comm_rank (comm, &(ptb->i_rank));
  ptb->is_my_rank_active   = 0;              /*!< Is active current rank */
  ptb->part_active_node    = part_active_node; /*!< Part of active nodes */

  ptb->n_part             = n_part;       /*!< Number of parts */
  ptb->n_elt              = n_elt;        /*!< Number of elements for any part */
  ptb->n_elt_proc         = 0;            /*!< Number of elements on the current rank */
  ptb->gnum_elt           = gnum_elt;     /*!< Global numbering of elements for any part */
  ptb->weight             = weight;
  ptb->dest_proc          = NULL;
  PDM_malloc(ptb->data_distrib_index, ptb->s_comm + 1, PDM_g_num_t);   /*!< Data distribution on ranks */

  ptb->s_block_min   = INT_MAX;
  ptb->s_block_max   = 0;

  ptb->i_send_data   = NULL;  /*!< Data to send to other processes index (size = s_comm) */
  ptb->i_recv_data   = NULL;  /*!< Received Data from other processes index (size = s_comm) */
  ptb->n_send_data   = NULL;  /*!< Number of data to send to other processes (size = s_comm) */
  ptb->n_recv_data   = NULL;  /*!< Number of received Data from other processes (size = s_comm) */

  ptb->tn_send_data      = 0;     /*!< Total number of sended data */
  ptb->tn_recv_data      = 0;     /*!< Total number of received data */
  ptb->sorted_recv_gnum  = NULL;  /*!< Sorted recv global num */
  ptb->order             = NULL;  /*!< Order */
  ptb->n_elt_block       = 0;
  ptb->block_gnum        = NULL;  /*!< Global number of reveived data (size = tn_recv_data) */
  ptb->block_gnum_count  = NULL;  /*!< Number of occurence of reveived data (size = tn_recv_data) */

  ptb->weight_g          = NULL; /*!< Global weights of elements for any part */
  ptb->enable_reverse    = 0;
  ptb->idx_partial       = NULL;

  /* Asynchone */
  ptb->max_exch_request = 10;
  ptb->next_request     = 0;
  PDM_malloc(ptb->s_data    , ptb->max_exch_request, size_t           );
  PDM_malloc(ptb->t_stride  , ptb->max_exch_request, PDM_stride_t     );
  PDM_malloc(ptb->cst_stride, ptb->max_exch_request, int              );
  ptb->wait_status      =  PDM_array_const_int(ptb->max_exch_request, 2);
  PDM_malloc(ptb->request_mpi, ptb->max_exch_request, PDM_MPI_Request  );

  PDM_malloc(ptb->send_buffer  , ptb->max_exch_request, unsigned char   *);
  PDM_malloc(ptb->recv_buffer  , ptb->max_exch_request, unsigned char   *);
  PDM_malloc(ptb->recv_stride  , ptb->max_exch_request, int             *);
  PDM_malloc(ptb->n_send_buffer, ptb->max_exch_request, int             *);
  PDM_malloc(ptb->i_send_buffer, ptb->max_exch_request, int             *);
  PDM_malloc(ptb->n_recv_buffer, ptb->max_exch_request, int             *);
  PDM_malloc(ptb->i_recv_buffer, ptb->max_exch_request, int             *);

  PDM_malloc(ptb->block_stride, ptb->max_exch_request, int            **);
  PDM_malloc(ptb->block_data  , ptb->max_exch_request, void           **);
  PDM_malloc(ptb->part_stride , ptb->max_exch_request, int           ***);
  PDM_malloc(ptb->part_data   , ptb->max_exch_request, void          ***);

  PDM_malloc(ptb->comm_kind, ptb->max_exch_request, PDM_mpi_comm_kind_t );
  PDM_malloc(ptb->win_send , ptb->max_exch_request, PDM_MPI_Win         );
  PDM_malloc(ptb->win_recv , ptb->max_exch_request, PDM_MPI_Win         );
  PDM_malloc(ptb->mpi_type , ptb->max_exch_request, PDM_MPI_Datatype    );

  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    ptb->send_buffer  [i_req] = NULL;
    ptb->recv_buffer  [i_req] = NULL;
    ptb->recv_stride  [i_req] = NULL;
    ptb->n_send_buffer[i_req] = NULL;
    ptb->i_send_buffer[i_req] = NULL;
    ptb->n_recv_buffer[i_req] = NULL;
    ptb->i_recv_buffer[i_req] = NULL;
    ptb->block_stride [i_req] = NULL;
    ptb->block_data   [i_req] = NULL;
    ptb->part_stride  [i_req] = NULL;
    ptb->part_data    [i_req] = NULL;
    ptb->win_send     [i_req] = PDM_MPI_WIN_NULL;
    ptb->win_recv     [i_req] = PDM_MPI_WIN_NULL;
    ptb->mpi_type     [i_req] = PDM_MPI_DATATYPE_NULL;
  }

  ptb->p2p_factor = 0.25;

  char host[1024];
  gethostname(host, 1023);

  if (!strncmp(host, "sator" , 5)) {
    ptb->p2p_factor = -0.1;
  }

  char *env_var = NULL;
  env_var = getenv ("PDM_PART_TO_BLOCK_P2P_FACTOR");
  if (env_var != NULL) {
    ptb->p2p_factor = atof (env_var);
  }

  /*
   * Active ranks definition
   */

  _active_ranks (ptb);

  return ptb;
}


static
void
_prepare_exchange
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                  **part_stride,
  size_t                *i_send_buffer,
  size_t                *i_recv_buffer,
  int                   *n_send_buffer,
  int                   *n_recv_buffer,
  int                  **recv_stride
)
{
  PDM_UNUSED(s_data);
  PDM_UNUSED(cst_stride);
  /*
   * Exchange Stride and build buffer properties
   */
  int *_recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *send_stride = NULL;
    PDM_malloc(send_stride, ptb->tn_send_data, int);

    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        send_stride[ptb->i_send_data[iproc] + n_send_buffer[iproc]] = part_stride[i][j];
        n_send_buffer[iproc] += 1;
      }
    }

    PDM_malloc(_recv_stride, ptb->tn_recv_data, int);
    assert (send_stride != NULL);
    assert (_recv_stride != NULL);

    *recv_stride = _recv_stride;

    if (ptb->p2p_factor < ptb->part_active_rank) {

      PDM_MPI_Alltoallv (send_stride,
                         ptb->n_send_data,
                         ptb->i_send_data,
                         PDM_MPI_INT,
                         _recv_stride,
                         ptb->n_recv_data,
                         ptb->i_recv_data,
                         PDM_MPI_INT,
                         ptb->comm);
    }

    else {

      PDM_MPI_Alltoallv_p2p (send_stride,
                             ptb->n_send_data,
                             ptb->i_send_data,
                             PDM_MPI_INT,
                             _recv_stride,
                             ptb->n_recv_data,
                             ptb->i_recv_data,
                             PDM_MPI_INT,
                             ptb->comm);

    }

    /*
     * Build buffers
     */

    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_send_data[i];
      int iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_send_buffer[i] += send_stride[k];
      }

      // n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_recv_data[i];
      iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++) {
        n_recv_buffer[i] += _recv_stride[k];
      }

      // n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }

    }

    PDM_free(send_stride);

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_send_data[i]; //  * cst_stride; //  * (int) s_data;
      i_recv_buffer[i] = ptb->i_recv_data[i]; //  * cst_stride; //  * (int) s_data;

      n_send_buffer[i] = ptb->n_send_data[i]; //  * cst_stride; //  * (int) s_data;
      n_recv_buffer[i] = ptb->n_recv_data[i]; //  * cst_stride; //  * (int) s_data;

    }
  }
}


static
void
_prepare_reverse_exchange
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                   *block_strid,
  size_t                *i_send_buffer,
  size_t                *i_recv_buffer,
  int                   *n_send_buffer,
  int                   *n_recv_buffer,
  int                  **send_stride,
  int                  **recv_stride
)
{
  PDM_UNUSED(s_data);
  PDM_UNUSED(cst_stride);
  /*
   * Exchange Stride and build buffer properties
   */
  int *_recv_stride = NULL;
  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {
      n_send_buffer[i] = 0;
      n_recv_buffer[i] = 0;
    }

    int *_send_stride = NULL;
    PDM_malloc(_send_stride, ptb->tn_recv_data, int);
    *send_stride = _send_stride;

    /* Because it's reverse we have tn_recn_data for send -_-' */
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      _send_stride[ptb->order[i]] = block_strid[ptb->idx_partial[i]];
    }
    // PDM_log_trace_array_int(_send_stride, ptb->tn_recv_data, "_send_stride :: ");

    int s_recv = ptb->i_send_data[ptb->s_comm-1] + ptb->n_send_data[ptb->s_comm-1];
    PDM_malloc(_recv_stride, s_recv, int);
    assert (_send_stride != NULL);
    assert (_recv_stride != NULL);

    *recv_stride = _recv_stride;

    if (ptb->p2p_factor < ptb->part_active_rank) {
      PDM_MPI_Alltoallv (_send_stride,
                         ptb->n_recv_data,
                         ptb->i_recv_data,
                         PDM_MPI_INT,
                         _recv_stride,
                         ptb->n_send_data,
                         ptb->i_send_data,
                         PDM_MPI_INT,
                         ptb->comm); 
    }
    else {
      PDM_MPI_Alltoallv_p2p (_send_stride,
                             ptb->n_recv_data,
                             ptb->i_recv_data,
                             PDM_MPI_INT,
                             _recv_stride,
                             ptb->n_send_data,
                             ptb->i_send_data,
                             PDM_MPI_INT,
                             ptb->comm);      
    }

    if(0 == 1) {
      PDM_log_trace_array_int(_recv_stride, s_recv, "recv_stride :: ");
    }

    /*
     * Build buffers
     */

    for (int i = 0; i < ptb->s_comm; i++) {
      int ibeg = ptb->i_recv_data[i];
      int iend = ptb->i_recv_data[i] + ptb->n_recv_data[i];

      n_send_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++){
        n_send_buffer[i] += _send_stride[k];
      }

      // n_send_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_send_buffer[i] = i_send_buffer[i-1] + n_send_buffer[i-1];
      } else {
        i_send_buffer[i] = 0;
      }

      ibeg = ptb->i_send_data[i];
      iend = ptb->i_send_data[i] + ptb->n_send_data[i];

      n_recv_buffer[i] = 0;
      for (int k = ibeg; k < iend; k++){
        n_recv_buffer[i] += _recv_stride[k];
      }

      // n_recv_buffer[i] *= (int) s_data;

      if (i > 0) {
        i_recv_buffer[i] = i_recv_buffer[i-1] + n_recv_buffer[i-1];
      } else {
        i_recv_buffer[i] = 0;
      }

    }

  }

  else if (t_stride == PDM_STRIDE_CST_INTERLACED) {

    for (int i = 0; i < ptb->s_comm; i++) {

      i_send_buffer[i] = ptb->i_recv_data[i]; // * cst_stride; // * (int) s_data;
      i_recv_buffer[i] = ptb->i_send_data[i]; // * cst_stride; // * (int) s_data;

      n_send_buffer[i] = ptb->n_recv_data[i]; // * cst_stride; // * (int) s_data;
      n_recv_buffer[i] = ptb->n_send_data[i]; // * cst_stride; // * (int) s_data;

    }

  }

  if(0 == 1) {
    PDM_log_trace_array_size_t(i_send_buffer, ptb->s_comm, "i_send_buffer ::");
    PDM_log_trace_array_size_t(i_recv_buffer, ptb->s_comm, "i_recv_buffer ::");
    PDM_log_trace_array_int   (n_send_buffer, ptb->s_comm, "n_send_buffer ::");
    PDM_log_trace_array_int   (n_recv_buffer, ptb->s_comm, "n_recv_buffer ::");
  }
}


static
void
_prepare_send_buffer
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                  **part_stride,
  void                 **part_data,
  size_t                *i_send_buffer,
  int                   *n_send_buffer,
  unsigned char         *send_buffer
)
{
  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  unsigned char **_part_data = (unsigned char **) part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_send_buffer[i] = 0;
  }

  int idx = -1;
  for (int i = 0; i < ptb->n_part; i++) {

    size_t *i_part = NULL;
    if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
      PDM_malloc(i_part, ptb->n_elt[i] + 1, size_t);
      assert (i_part != NULL);

      i_part[0] = 0;
      for (int j = 1; j < ptb->n_elt[i] + 1; j++) {
        i_part[j] = i_part[j-1] + ((size_t) part_stride[i][j-1] * s_data);
      }
    }

    for (int j = 0; j < ptb->n_elt[i]; j++) {
      int iproc = ptb->dest_proc[++idx];
      size_t s_octet_elt = 0;
      size_t i_part_elt = 0;

      if (t_stride == PDM_STRIDE_CST_INTERLACED) {
        s_octet_elt = cst_stride * (int) s_data;
        i_part_elt  = cst_stride * (int) s_data * j;
      }

      else if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
        s_octet_elt = i_part[j+1] - i_part[j];
        i_part_elt  = i_part[j];
      }

      int idx_write = (i_send_buffer[iproc] + n_send_buffer[iproc]) * s_data_tot;
      if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
        n_send_buffer[iproc] += part_stride[i][j];
      } else {
        n_send_buffer[iproc] += 1;
      }

      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[idx_write + k] =
          _part_data[i][i_part_elt + k];
      }
    }

    if (i_part != NULL)
      PDM_free(i_part);
  }
}



static
void
_prepare_reverse_send_buffer
(
  PDM_part_to_block_t   *ptb,
  size_t                 s_data,
  PDM_stride_t           t_stride,
  int                    cst_stride,
  int                   *send_stride,
  int                   *block_stride,
  void                  *block_data,
  size_t                *i_send_buffer,
  int                   *n_send_buffer,
  unsigned char         *send_buffer
)
{
  // unsigned char **_part_data = (unsigned char **) part_data;
  unsigned char *_block_data = (unsigned char *) block_data;

  int s_data_tot = s_data;
  if (t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  int s_send_buffer = (i_send_buffer[ptb->s_comm-1] + n_send_buffer[ptb->s_comm-1]) * s_data_tot;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    /* Because it's reverse we have tn_recn_data for send -_-' */
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      send_stride[ptb->order[i]] = block_stride[ptb->idx_partial[i]];
    }

    int *send_idx = NULL;
    PDM_malloc(send_idx, ptb->tn_recv_data + 1, int);
    send_idx[0] = 0;
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      send_idx[i+1] = send_idx[i] + send_stride[i];
    }

    // Same for block_stride
    int *block_idx = NULL;
    PDM_malloc(block_idx, ptb->n_elt_block + 1, int);
    block_idx[0] = 0;
    for(int i = 0; i < ptb->n_elt_block; ++i) {
      block_idx[i+1] = block_idx[i] + block_stride[i];
    }


    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      size_t s_octet_elt = send_stride[ptb->order[i]] * s_data;
      int idx_write = send_idx [ptb->order      [i]] * s_data;
      int idx_read  = block_idx[ptb->idx_partial[i]] * s_data;
      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[idx_write+k] = _block_data[idx_read+k];
      }
    }

    PDM_free(send_stride);
    PDM_free(send_idx);
    PDM_free(block_idx);

  }
  else  {
    size_t s_octet_elt = cst_stride * (int) s_data;
    for(int i = 0; i < ptb->tn_recv_data; ++i) {
      for (int k = 0; k < (int) s_octet_elt; k++) {
        send_buffer[s_octet_elt * ptb->order[i] + k] = _block_data[s_octet_elt*ptb->idx_partial[i]+k];
      }
    }
  }

  if(0 == 1) {
    int* tmp_send_buffer = (int *) send_buffer;
    PDM_log_trace_array_int(tmp_send_buffer, s_send_buffer/sizeof(int), "tmp_send_buffer ::");
  }

}



/**
 *
 * \brief  Post-treatment of the resulting buffer
 *
 * \param [inout]   ptb          Part to block structure
 *
 */
static
int
_post_treatment
(
  PDM_part_to_block_t  *ptb,
  size_t                s_data,
  PDM_stride_t          t_stride,
  int                   cst_stride,
  int                  *recv_stride,
  unsigned char        *recv_buffer,
  size_t                s_recv_buffer,
  int                 **block_stride,
  void                **block_data
)
{
  unsigned char *_block_data = NULL;
 PDM_malloc(_block_data, s_recv_buffer, unsigned char);
  assert(_block_data != NULL);
  *block_data = _block_data;

  int *i_recv_stride = NULL;
  int *i_block_stride = NULL;
  int s_block_data = ((int) sizeof(unsigned char) * s_recv_buffer) / (int) s_data;

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {
    *block_stride = NULL;
    int* _block_stride = NULL;
    if(ptb->tn_recv_data > 0){
      PDM_malloc(_block_stride, ptb->tn_recv_data, int);
    }
    *block_stride = _block_stride;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      _block_stride[i] = recv_stride[ptb->order[i]];
    }

    /*
     * Compute index in data
     */
    PDM_malloc(i_recv_stride , ptb->tn_recv_data + 1, int);
    PDM_malloc(i_block_stride, ptb->tn_recv_data + 1, int);

    i_recv_stride [0] = 0;
    i_block_stride[0] = 0;
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride [i+1] = i_recv_stride [i] + recv_stride[i];
      i_block_stride[i+1] = i_block_stride[i] + _block_stride[i];
    }

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      i_recv_stride [i+1] *= (int) s_data;
      i_block_stride[i+1] *= (int) s_data;
    }

    /*
     * Sort Buffer
     */

    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int old    = ptb->order[i];
      int id_old = i_recv_stride[old];
      for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
        _block_data[k] = recv_buffer[id_old++];
      }
    }

    PDM_free(recv_stride);
    PDM_free(i_recv_stride);

    /*
     * post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {

      int idx1 = 0;
      int idx2 = 0;

      if (ptb->tn_recv_data == 1) {
        idx2 = i_block_stride[1];
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        if (i == 1) {
          idx2 = i_block_stride[1];
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          _block_stride[idx1] = _block_stride[i];
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            for (int k = i_block_stride[i]; k < i_block_stride[i+1]; k++) {
              _block_data[idx2++] = _block_data[k];
            }
          }
        } else {
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
            _block_stride[idx1] += _block_stride[i];
          }
        }
      }

      /* Cleanup */

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        PDM_realloc(_block_data ,_block_data , idx2,unsigned char);
        *block_data = _block_data;

        PDM_realloc(_block_stride ,_block_stride , ptb->n_elt_block,int);

        *block_stride = _block_stride;
        s_block_data = idx2 / (int) s_data;
      }
    }
    PDM_free(i_block_stride);
  } else {

    /*
     * Sort Buffer
     */
    for (int i = 0; i < ptb->tn_recv_data; i++) {
      int n_octet = cst_stride * (int) s_data;
      int old = ptb->order[i];
      int id_old = old * n_octet;

      for (int k = i*n_octet; k < (i+1)*n_octet; k++) {
        _block_data[k] = recv_buffer[id_old++];
      }
    }

    /*
     * Post processing
     */

    if (ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) {
      int idx2 = 0;
      int idx1 = 0;

      assert (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE);

      if (ptb->tn_recv_data == 1) {
        idx2 =  cst_stride * (int) s_data;
      }

      for (int i = 1; i < ptb->tn_recv_data; i++) {
        int n_octet = cst_stride * (int) s_data;
        if (i == 1) {
          idx2 = n_octet;
        }
        if (ptb->block_gnum[idx1] != ptb->sorted_recv_gnum[i]) {
          idx1 += 1;
          if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
            int idx3 = i * cst_stride * (int) s_data;
            for (int k = 0; k < n_octet; k++) {
              _block_data[idx2++] = _block_data[idx3++];
            }
          }
        }
      }

      if (ptb->t_post == PDM_PART_TO_BLOCK_POST_CLEANUP) {
        PDM_realloc(_block_data, _block_data, idx2, unsigned char);
        *block_data = _block_data;
        s_block_data = idx2 / (int) s_data;
      }

    }
  }

  return s_block_data;
}

/**
 *
 * \brief  Post-treatment of the resulting buffer
 *
 * \param [inout]   ptb          Part to block structure
 *
 */
static
void
_post_treatment_reverse
(
  PDM_part_to_block_t  *ptb,
  size_t                s_data,
  PDM_stride_t          t_stride,
  int                   cst_stride,
  int                  *recv_stride,
  unsigned char        *recv_buffer,
  int                  *n_recv_buffer,
  size_t               *i_recv_buffer,
  int                ***part_stride,
  void               ***part_data
)
{
  PDM_malloc(*((unsigned char ***) part_data), ptb->n_part, unsigned char *);
  unsigned char **_part_data = (unsigned char **) *part_data;

  for (int i = 0; i <  ptb->s_comm; i++) {
    n_recv_buffer[i] = 0;
  }

  if (t_stride == PDM_STRIDE_VAR_INTERLACED) {

    PDM_malloc(*part_stride, ptb->n_part, int *);
    int **_part_stride = *part_stride;
    int **_part_idx;
    PDM_malloc(_part_idx, ptb->n_part, int *);

    int *n_recv_strid;
    PDM_malloc(n_recv_strid, ptb->s_comm, int);
    for(int i = 0; i < ptb->s_comm; ++i) {
      n_recv_strid[i] = 0;
    }

    /*
     * Post-treated recv_strid
     */
    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      PDM_malloc(_part_stride[i], ptb->n_elt[i]  , int);
      PDM_malloc(_part_idx   [i], ptb->n_elt[i]+1, int);

      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        int idx_read = ptb->i_send_data[iproc] + n_recv_strid[iproc]++;
        _part_stride[i][j] = recv_stride[idx_read];
      }

      _part_idx[i][0] = 0;
      for(int j = 0; j < ptb->n_elt[i]; j++) {
        _part_idx[i][j+1] = _part_idx[i][j] + _part_stride[i][j];
      }

      if(0 == 1) {
        PDM_log_trace_array_int(_part_idx[i], ptb->n_elt[i]+1, "_part_idx :: ");
      }

    }

    idx = -1;
    int n_octet = (int) s_data;
    for (int i = 0; i < ptb->n_part; i++) {
      PDM_malloc(_part_data[i], _part_idx[i][ptb->n_elt[i]] * n_octet, unsigned char);
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        int idx_write = n_octet * _part_idx[i][j];
        int idx_read  = (i_recv_buffer[iproc] + n_recv_buffer[iproc]) * s_data;
        n_recv_buffer[iproc] += _part_stride[i][j];
        for (int k = 0; k < (int) n_octet * _part_stride[i][j]; k++) {
          _part_data[i][idx_write + k] = recv_buffer[idx_read+k];
        }
      }
      PDM_free(_part_idx[i]);
    }
    PDM_free(n_recv_strid);
    PDM_free(_part_idx);
    PDM_free(recv_stride);

  } else { // PDM_STRIDE_CST_INTERLACED

    int n_octet = cst_stride * (int) s_data;
    int idx = -1;
    for (int i = 0; i < ptb->n_part; i++) {
      PDM_malloc(_part_data[i], ptb->n_elt[i] * n_octet, unsigned char);
      for (int j = 0; j < ptb->n_elt[i]; j++) {
        int iproc = ptb->dest_proc[++idx];
        int idx_read  = (i_recv_buffer[iproc] + n_recv_buffer[iproc]) * s_data * cst_stride;
        n_recv_buffer[iproc] += 1;
        for (int k = 0; k < (int) n_octet; k++) {
          _part_data[i][n_octet*j + k] = recv_buffer[idx_read+k];
        }
      }
    }
  }
}





/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Reset global statistic
 *
 */

void
PDM_part_to_block_global_statistic_reset
(
  void
)
{
  for (int i = 0; i < NTIMER_PTB; i++) {
    t_elaps[i] = 0;
    t_cpu[i] = 0;
  }

  for (int i = 0; i < 2; i++) {
    min_exch_rank[i] = INT_MAX;
    max_exch_rank[i] = -1;
    exch_data[i] = 0;
  }
}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm                 MPI communicator
 * \param [out]  min_exch_rank_send   Global min part of ranks used to send
 * \param [out]  min_exch_rank_recv   Global min part of ranks used to receive
 * \param [out]  max_exch_rank_send   Global max part of ranks used to send
 * \param [out]  max_exch_rank_recv   Global max part of ranks used to receive
 * \param [out]  min_exch_data_send   Global min sent data for a rank
 * \param [out]  min_exch_data_recv   Global min received data for a rank
 * \param [out]  max_exch_data_send   Global max sent data for a rank
 * \param [out]  max_exch_data_recv   Global max received data for a rank
 *
 */

void
PDM_part_to_block_global_statistic_get
(
 PDM_MPI_Comm comm,
 int *min_exch_rank_send,
 int *min_exch_rank_recv,
 int *max_exch_rank_send,
 int *max_exch_rank_recv,
 unsigned long long *min_exch_data_send,
 unsigned long long *min_exch_data_recv,
 unsigned long long *max_exch_data_send,
 unsigned long long *max_exch_data_recv
)
{
  unsigned long long max_exch_data[2];
  unsigned long long min_exch_data[2];

  PDM_MPI_Allreduce (exch_data, min_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (exch_data, max_exch_data, 2,
                     PDM_MPI_UNSIGNED_LONG_LONG, PDM_MPI_MAX, comm);

  *min_exch_data_send = min_exch_data[0];
  *min_exch_data_recv = min_exch_data[1];
  *max_exch_data_send = max_exch_data[0];
  *max_exch_data_recv = max_exch_data[1];


  int max_max_exch_rank[2];
  int min_min_exch_rank[2];

  PDM_MPI_Allreduce (min_exch_rank, min_min_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (max_exch_rank, max_max_exch_rank, 2,
                     PDM_MPI_INT, PDM_MPI_MAX, comm);

  *min_exch_rank_send = min_min_exch_rank[0];
  *min_exch_rank_recv = min_min_exch_rank[1];
  *max_exch_rank_send = max_max_exch_rank[0];
  *max_exch_rank_recv = max_max_exch_rank[1];

}


/**
 *
 * \brief Get global timer in part to block
 *
 * \param [in]   comm              MPI communicator
 * \param [out]  min_elaps         Min elapsed time
 * \param [out]  max_elaps         Max elapsed time
 * \param [out]  min_cpu           Min cpu time
 * \param [out]  max_cpu           Max cpu time
 * \param [out]  min_elaps_create  Global min elapsed for create function
 * \param [out]  max_elaps_create  Global max elapsed for create function
 * \param [out]  min_cpu_create    Global min cpu for create function
 * \param [out]  max_cpu_create    Global max cpu for create function
 * \param [out]  min_elaps_create2 Global min elapsed for create2 function
 * \param [out]  max_elaps_create2 Global max elapsed for create2 function
 * \param [out]  min_cpu_create2   Global min cpu for create2 function
 * \param [out]  max_cpu_create2   Global max cpu for create2 function
 * \param [out]  min_elaps_exch    Global min elapsed for exch function
 * \param [out]  max_elaps_exch    Global max elapsed for exch function
 * \param [out]  min_cpu_exch      Global min cpu for exch function
 * \param [out]  max_cpu_exch      Global max cpu for exch function
 *
 */

void
PDM_part_to_block_global_timer_get
(
 PDM_MPI_Comm comm,
 double       *min_elaps_create,
 double       *max_elaps_create,
 double       *min_cpu_create,
 double       *max_cpu_create,
 double       *min_elaps_create2,
 double       *max_elaps_create2,
 double       *min_cpu_create2,
 double       *max_cpu_create2,
 double       *min_elaps_exch,
 double       *max_elaps_exch,
 double       *min_cpu_exch,
 double       *max_cpu_exch
)
{

  double min_elaps[NTIMER_PTB];
  double max_elaps[NTIMER_PTB];
  double min_cpu[NTIMER_PTB];
  double max_cpu[NTIMER_PTB];

  PDM_MPI_Allreduce (t_elaps, min_elaps, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_elaps, max_elaps, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  PDM_MPI_Allreduce (t_cpu, min_cpu, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_cpu, max_cpu, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  // Part-to-Block creation with block-distribution generation
  *min_elaps_create  = min_elaps[MALLOC_ACTIVE_RANKS] + min_elaps[GENERATE_DISTRIB] + min_elaps[BINARY_SEARCH] + min_elaps[CREATE_EXCHANGE] + min_elaps[BLOCK_POST] + min_elaps[GLOBAL_WEIGHTS];
  *max_elaps_create  = max_elaps[MALLOC_ACTIVE_RANKS] + max_elaps[GENERATE_DISTRIB] + max_elaps[BINARY_SEARCH] + max_elaps[CREATE_EXCHANGE] + max_elaps[BLOCK_POST] + max_elaps[GLOBAL_WEIGHTS];
  *min_cpu_create    = min_cpu[MALLOC_ACTIVE_RANKS]   + min_cpu[GENERATE_DISTRIB]   + min_cpu[BINARY_SEARCH]   + min_cpu[CREATE_EXCHANGE]   + min_cpu[BLOCK_POST]   + min_cpu[GLOBAL_WEIGHTS]  ;
  *max_cpu_create    = max_cpu[MALLOC_ACTIVE_RANKS]   + max_cpu[GENERATE_DISTRIB]   + max_cpu[BINARY_SEARCH]   + max_cpu[CREATE_EXCHANGE]   + max_cpu[BLOCK_POST]   + max_cpu[GLOBAL_WEIGHTS]  ;
  // Part-to-Block creation with user provided block-distribution
  *min_elaps_create2 = min_elaps[CREATE_FROM_DISTRIB];
  *max_elaps_create2 = max_elaps[CREATE_FROM_DISTRIB];
  *min_cpu_create2   = min_cpu[CREATE_FROM_DISTRIB];
  *max_cpu_create2   = max_cpu[CREATE_FROM_DISTRIB];
  // Warning : Geometric Part-to-Block creation is not outputed while exchanges are counted
  // Data exchange
  *min_elaps_exch    = min_elaps[DATA_EXCHANGE];
  *max_elaps_exch    = max_elaps[DATA_EXCHANGE];
  *min_cpu_exch      = min_cpu[DATA_EXCHANGE];
  *max_cpu_exch      = max_cpu[DATA_EXCHANGE];

}

/**
 *
 * \brief Global write part-to-block step timer
 *
 * \param [in]  comm            MPI communicator
 * \param [in]  filename        File name
 *
 */

void
PDM_part_to_block_time_per_step_dump
(
 PDM_MPI_Comm  comm,
 const char   *filename
)
{
  // Write in parallel
  PDM_io_file_t *writer = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(filename,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_APPEND,
              PDM_IO_NATIVE,
              comm,
              -1.,
              &writer,
              &ierr);

  // MPI
  int n_rank = 0;
  PDM_MPI_Comm_size (comm, &n_rank);

  // Create timer statistics
  double min_elaps[NTIMER_PTB];
  double mean_elaps[NTIMER_PTB];
  double max_elaps[NTIMER_PTB];
  double min_cpu[NTIMER_PTB];
  double mean_cpu[NTIMER_PTB];
  double max_cpu[NTIMER_PTB];

  PDM_MPI_Allreduce (t_elaps, min_elaps, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_elaps, mean_elaps, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  PDM_MPI_Allreduce (t_elaps, max_elaps, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  PDM_MPI_Allreduce (t_cpu, min_cpu, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);

  PDM_MPI_Allreduce (t_cpu, mean_cpu, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  PDM_MPI_Allreduce (t_cpu, max_cpu, NTIMER_PTB,
                     PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  for (int i_step = 0; i_step < NTIMER_PTB; i_step++) {
    min_elaps[i_step]  /= n_ptb_open;
    mean_elaps[i_step] /= n_ptb_open;
    max_elaps[i_step]  /= n_ptb_open;

    min_cpu[i_step]  /= n_ptb_open;
    mean_cpu[i_step] /= n_ptb_open;
    max_cpu[i_step]  /= n_ptb_open;

    mean_elaps[i_step] /= n_rank;
    mean_cpu[i_step]   /= n_rank;
  } // end loop on timed steps

  // Global write times
  size_t s_buffer = 436; // buffer size for %.5f + 1
  char *buffer;
  PDM_malloc(buffer, s_buffer, char);

  for (int i = 0; i < (int) s_buffer; i++) {
    buffer[i] = '\0';
  }

  sprintf(buffer, "generate_distrib elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[GENERATE_DISTRIB], mean_elaps[GENERATE_DISTRIB], max_elaps[GENERATE_DISTRIB], min_cpu[GENERATE_DISTRIB], mean_cpu[GENERATE_DISTRIB], max_cpu[GENERATE_DISTRIB]);

  sprintf(buffer + strlen(buffer), "binary_search elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[BINARY_SEARCH], mean_elaps[BINARY_SEARCH], max_elaps[BINARY_SEARCH], min_cpu[BINARY_SEARCH], mean_cpu[BINARY_SEARCH], max_cpu[BINARY_SEARCH]);

  sprintf(buffer + strlen(buffer), "create_exchange elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[CREATE_EXCHANGE], mean_elaps[CREATE_EXCHANGE], max_elaps[CREATE_EXCHANGE], min_cpu[CREATE_EXCHANGE], mean_cpu[CREATE_EXCHANGE], max_cpu[CREATE_EXCHANGE]);

  sprintf(buffer + strlen(buffer), "block_post elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[BLOCK_POST], mean_elaps[BLOCK_POST], max_elaps[BLOCK_POST], min_cpu[BLOCK_POST], mean_cpu[BLOCK_POST], max_cpu[BLOCK_POST]);

  sprintf(buffer + strlen(buffer), "global_weights elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[GLOBAL_WEIGHTS], mean_elaps[GLOBAL_WEIGHTS], max_elaps[GLOBAL_WEIGHTS], min_cpu[GLOBAL_WEIGHTS], mean_cpu[GLOBAL_WEIGHTS], max_cpu[GLOBAL_WEIGHTS]);

  // Warning : Geometric Part-to-Block creation is not outputed while exchanges are counted
  sprintf(buffer + strlen(buffer), "data_exchange elaps %.5f %.5f %.5f cpu %.5f %.5f %.5f\n", min_elaps[DATA_EXCHANGE], mean_elaps[DATA_EXCHANGE], max_elaps[DATA_EXCHANGE], min_cpu[DATA_EXCHANGE], mean_cpu[DATA_EXCHANGE], max_cpu[DATA_EXCHANGE]);

  PDM_io_global_write(writer,
                      (PDM_l_num_t) sizeof(char),
                      (PDM_l_num_t) s_buffer,
                      buffer);

  PDM_free(buffer);

  // Finalize parallel write
  PDM_io_close(writer);
  PDM_io_free(writer);
}

/**
 *
 * \brief Write in parallel communication graph
 *
 * \param [in]  ptb             Part-to-Block structure
 * \param [in]  filename        File name
 *
 */

void
PDM_part_to_block_comm_graph_dump
(
 PDM_part_to_block_t *ptb,
 const char          *filename
)
{
  // Write in parallel
  PDM_io_file_t *writer = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(filename,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              ptb->comm,
              -1.,
              &writer,
              &ierr);

  // Create a node identifier
  PDM_MPI_Comm shared_comm = PDM_MPI_COMM_WORLD;

  PDM_MPI_Comm_split_type(ptb->comm, PDM_MPI_SPLIT_SHARED, &shared_comm);

  int i_shared_rank = 0;
  PDM_MPI_Comm_rank(shared_comm, &i_shared_rank);

  int bcast_buffer = 0;
  if (i_shared_rank == 0) {
    bcast_buffer = ptb->i_rank;
  }
  PDM_MPI_Bcast(&bcast_buffer, 1, PDM_MPI_INT32_T, 0, shared_comm);

  // Block write i_rank, node and number of send data
  int s_buffer = ptb->s_comm * 11 + 40 + 2 + 1; // (10 + 1 space) * n_rank + chaine + space + \n + 1
  char *buffer;
  PDM_malloc(buffer, s_buffer, char);

  for (int i = 0; i < (int) s_buffer; i++) {
    buffer[i] = '\0';
  }

  sprintf(buffer, "i_rank %10d\nnode %10d\nn_send", ptb->i_rank, bcast_buffer);

  for (int j_rank = 0; j_rank < ptb->s_comm; j_rank++) {
    sprintf(buffer + strlen(buffer), " %10d", ptb->n_send_data[j_rank]);
  } // end loop on n_rank
  sprintf(buffer + strlen(buffer), " \n");

  PDM_l_num_t one = 1;
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (ptb->i_rank+1);
  PDM_io_par_interlaced_write(writer,
                              PDM_STRIDE_VAR_INTERLACED,
                              (PDM_l_num_t *) &s_buffer,
                              (PDM_l_num_t) sizeof(char),
                              one,
                              &i_rank_gnum,
                              (const void *) buffer);

  PDM_free(buffer);

  // Finalize parallel write
  PDM_io_close(writer);
  PDM_io_free(writer);
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */

  // Warning : timing of _ptb_create not considered because induces issues

  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         weight,
                                         n_elt,
                                         n_part,
                                         comm);

  /*
   * Data distribution definition
   */

  _distrib_data (ptb, 0);

  /*
   * Compute global weight for each element
   */
  double t1_elaps = PDM_timer_elapsed(t_timer[GLOBAL_WEIGHTS]);
  double t1_cpu   = PDM_timer_cpu    (t_timer[GLOBAL_WEIGHTS]);
  PDM_timer_resume(t_timer[GLOBAL_WEIGHTS]);

  if (ptb->weight != NULL) {// && ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
    _compute_global_weights (ptb);
  }

  PDM_timer_hang_on(t_timer[GLOBAL_WEIGHTS]);
  double t2_elaps = PDM_timer_elapsed(t_timer[GLOBAL_WEIGHTS]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[GLOBAL_WEIGHTS]);

  t_elaps[GLOBAL_WEIGHTS] += (t2_elaps - t1_elaps);
  t_cpu  [GLOBAL_WEIGHTS] += (t2_cpu   - t1_cpu  );

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_create_from_distrib
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                         part_active_node,
 PDM_g_num_t                 **gnum_elt,
 const PDM_g_num_t            *data_distrib_index,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */

  // Warning : timing of _ptb_create not considered because induces issues

  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         NULL,
                                         n_elt,
                                         n_part,
                                         comm);

  double t1_elaps = PDM_timer_elapsed(t_timer[CREATE_FROM_DISTRIB]);
  double t1_cpu = PDM_timer_cpu(t_timer[CREATE_FROM_DISTRIB]);
  PDM_timer_resume(t_timer[CREATE_FROM_DISTRIB]);

  for(int i_rank = 0; i_rank < ptb->s_comm+1; i_rank++){
    ptb->data_distrib_index[i_rank] = data_distrib_index[i_rank];
  }

  /*
   * Data distribution definition
   */
  _distrib_data (ptb, 1);

  PDM_timer_hang_on(t_timer[CREATE_FROM_DISTRIB]);
  double t2_elaps = PDM_timer_elapsed(t_timer[CREATE_FROM_DISTRIB]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[CREATE_FROM_DISTRIB]);

  t_elaps[CREATE_FROM_DISTRIB] += (t2_elaps - t1_elaps);
  t_cpu  [CREATE_FROM_DISTRIB] += (t2_cpu   - t1_cpu);

  int n_rank_recv = 0;
  int n_rank_send = 0;


  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  return (PDM_part_to_block_t *) ptb;
}


/**
 *
 * \brief Create a partitioning to block redistribution
 *
 * \param [in]   t_distrib       Distribution type
 * \param [in]   t_post          Post processing type
 * \param [in]   part_active_node  Part of active nodes (\ref PDM_writer_BLOCK_DISTRIB_PART_OF_NODE mode)
 * \param [in]   gnum_elt        Element global number
 * \param [in]   weight          Weight of elements (or NULL)
 * \param [in]   n_elt           Local number of elements
 * \param [in]   n_part          Number of partition
 * \param [in]   comm            MPI communicator
 *
 * \return   Initialized PDM_part_to_block_t
 *
 */

PDM_part_to_block_t *
PDM_part_to_block_geom_create
(
 PDM_part_to_block_distrib_t   t_distrib,
 PDM_part_to_block_post_t      t_post,
 double                        part_active_node,
 PDM_part_geom_t               geom_kind,
 double                      **pvtx_coords,
 PDM_g_num_t                 **gnum_elt,
 double                      **weight,
 int                          *n_elt,
 int                           n_part,
 PDM_MPI_Comm                  comm
)
{
  /*
   * Common creation
   */

  // Warning : timing of _ptb_create not considered because induces issues

  PDM_part_to_block_t* ptb = _ptb_create(t_distrib,
                                         t_post,
                                         part_active_node,
                                         gnum_elt,
                                         NULL,
                                         n_elt,
                                         n_part,
                                         comm);

  double t1_elaps = PDM_timer_elapsed(t_timer[CREATE_GEOM]);
  double t1_cpu = PDM_timer_cpu(t_timer[CREATE_GEOM]);
  PDM_timer_resume(t_timer[CREATE_GEOM]);

  if(geom_kind == PDM_PART_GEOM_HILBERT ) {
    _distrib_data_hilbert(ptb, pvtx_coords, weight);
  } else if (geom_kind == PDM_PART_GEOM_MORTON ){
    _distrib_data_morton(ptb, pvtx_coords, weight);
  }

  /*
   * Compute global weight for each element
   */
  if (ptb->weight != NULL) {// && ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) {
    _compute_global_weights (ptb);
  }

  int n_rank_recv = 0;
  int n_rank_send = 0;

  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i && ptb->n_recv_data[i] > 0) {
      n_rank_recv += 1;
    }
    if (ptb->i_rank != i && ptb->n_send_data[i] > 0) {
      n_rank_send += 1;
    }
  }

  max_exch_rank[0] = PDM_MAX(max_exch_rank[0], n_rank_send);
  max_exch_rank[1] = PDM_MAX(max_exch_rank[1], n_rank_recv);
  min_exch_rank[0] = PDM_MIN(min_exch_rank[0], n_rank_send);
  min_exch_rank[1] = PDM_MIN(min_exch_rank[1], n_rank_recv);

  PDM_timer_hang_on(t_timer[CREATE_GEOM]);
  double t2_elaps = PDM_timer_elapsed(t_timer[CREATE_GEOM] );
  double t2_cpu = PDM_timer_cpu(t_timer[CREATE_GEOM]);

  t_elaps[CREATE_GEOM] += (t2_elaps - t1_elaps);
  t_cpu[CREATE_GEOM] += (t2_cpu - t1_cpu);

  return (PDM_part_to_block_t *) ptb;
}

/**
 *
 * \brief Return number of active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of active ranks
 *
 */

int
PDM_part_to_block_n_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_active_ranks;
}


/**
 *
 * \brief Return if current rank is active
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  if current rank is active
 *
 */

int
PDM_part_to_block_is_active_rank
(
 PDM_part_to_block_t *ptb
 )
{
  return ptb->is_my_rank_active;
}


/**
 *
 * \brief Return active ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  active ranks
 *
 */

int *
PDM_part_to_block_active_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->active_ranks;
}


/**
 *
 * \brief Return number of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_block_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->n_elt_block;
}


/**
 *
 * \brief Return global numbers of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers
 *
 */

PDM_g_num_t *
PDM_part_to_block_block_gnum_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum;
}

/**
 *
 * \brief Return numbers of occurence of each gnum element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global numbers counter
 *
 */

int *
PDM_part_to_block_block_gnum_count_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->block_gnum_count;
}


/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

int
PDM_part_to_block_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data,
 int                **block_stride,
 void               **block_data
)
{

  double t1_elaps = PDM_timer_elapsed(t_timer[DATA_EXCHANGE]);
  double t1_cpu = PDM_timer_cpu(t_timer[DATA_EXCHANGE]);
  PDM_timer_resume(t_timer[DATA_EXCHANGE]);

  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  size_t *i_send_buffer = NULL;
  size_t *i_recv_buffer = NULL;
  int    *n_send_buffer = NULL;
  int    *n_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, ptb->s_comm, size_t);
  PDM_malloc(i_recv_buffer, ptb->s_comm, size_t);
  PDM_malloc(n_send_buffer, ptb->s_comm, int   );
  PDM_malloc(n_recv_buffer, ptb->s_comm, int   );

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);


  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    i_send_buffer,
                    i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);

  /*
   * Prepare buffer
   */
  size_t s_send_buffer = ( i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1] ) * s_data_tot;
  size_t s_recv_buffer = ( i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1] ) * s_data_tot;

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;
  PDM_malloc(send_buffer, s_send_buffer, unsigned char);
  PDM_malloc(recv_buffer, s_recv_buffer, unsigned char);

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Data exchange
   */
  int mandatory_size = PDM_size_idx_from_stride (n_send_buffer, ptb->s_comm, ptb->comm);
  mandatory_size = PDM_MAX(PDM_size_idx_from_stride (n_send_buffer, ptb->s_comm, ptb->comm), mandatory_size);
  
  if (mandatory_size > 32) {
    PDM_MPI_Alltoallv_p2p_l(send_buffer,
                            n_send_buffer,
                            i_send_buffer,
                            mpi_type,
                            recv_buffer,
                            n_recv_buffer,
                            i_recv_buffer,
                            mpi_type,
                            ptb->comm);
  }

  else {
 
    if (ptb->p2p_factor < ptb->part_active_rank) {
     
      int *_i_send_buffer = NULL;
      int *_i_recv_buffer = NULL;
      PDM_malloc(_i_send_buffer, ptb->s_comm, int);
      PDM_malloc(_i_recv_buffer, ptb->s_comm, int);
     
      for (int i = 0; i < ptb->s_comm; i++) {
        _i_send_buffer[i] = (int) i_send_buffer[i];
        _i_recv_buffer[i] = (int) i_recv_buffer[i];
      }

      PDM_MPI_Alltoallv(send_buffer,
                        n_send_buffer,
                        _i_send_buffer,
                        mpi_type,
                        recv_buffer,
                        n_recv_buffer,
                        _i_recv_buffer,
                        mpi_type,
                        ptb->comm);

      PDM_free(_i_send_buffer); 
      PDM_free(_i_recv_buffer); 
    }

    else {

      PDM_MPI_Alltoallv_p2p_l(send_buffer,
                              n_send_buffer,
                              i_send_buffer,
                              mpi_type,
                              recv_buffer,
                              n_recv_buffer,
                              i_recv_buffer,
                              mpi_type,
                              ptb->comm);

    }
  } 

  /*
   * Statistics
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    if (ptb->i_rank != i) {
      exch_data[1] += n_recv_buffer[i];
      exch_data[0] += n_send_buffer[i];
    }
  }

  PDM_free(send_buffer);
  PDM_free(n_send_buffer);
  PDM_free(i_send_buffer);
  PDM_free(n_recv_buffer);
  PDM_free(i_recv_buffer);

  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     block_stride,
                                     block_data);
  PDM_free(recv_buffer);
  PDM_MPI_Type_free(&mpi_type);

  PDM_timer_hang_on(t_timer[DATA_EXCHANGE]);
  double t2_elaps = PDM_timer_elapsed(t_timer[DATA_EXCHANGE]);
  double t2_cpu   = PDM_timer_cpu    (t_timer[DATA_EXCHANGE]);

  t_elaps[DATA_EXCHANGE] += (t2_elaps - t1_elaps);
  t_cpu  [DATA_EXCHANGE] += (t2_cpu   - t1_cpu  );

  return s_block_data;
}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_reverse_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                 *block_stride,
 void                *block_data,
 int               ***part_stride,
 void              ***part_data
)
{
  assert(ptb->enable_reverse == 1);

  size_t *i_send_buffer = NULL;
  size_t *i_recv_buffer = NULL;
  int    *n_send_buffer = NULL;
  int    *n_recv_buffer = NULL;
  PDM_malloc(i_send_buffer, ptb->s_comm, size_t);
  PDM_malloc(i_recv_buffer, ptb->s_comm, size_t);
  PDM_malloc(n_send_buffer, ptb->s_comm, int   );
  PDM_malloc(n_recv_buffer, ptb->s_comm, int   );

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);

  /*
   * Exchange Stride and build buffer properties
   */
  int *send_stride = NULL; // Release in post-treatment
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_reverse_exchange(ptb,
                            s_data,
                            t_stride,
                            cst_stride,
                            block_stride,
                            i_send_buffer,
                            i_recv_buffer,
                            n_send_buffer,
                            n_recv_buffer,
                            &send_stride,
                            &recv_stride);

  /*
   * Prepare buffer
   */
  size_t s_send_buffer = ( i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1] ) * s_data_tot;
  size_t s_recv_buffer = ( i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1] ) * s_data_tot;

  unsigned char *send_buffer = NULL;
  unsigned char *recv_buffer = NULL;
 PDM_malloc(send_buffer, s_send_buffer, unsigned char);
 PDM_malloc(recv_buffer, s_recv_buffer, unsigned char);

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_reverse_send_buffer(ptb,
                               s_data,
                               t_stride,
                               cst_stride,
                               send_stride,
                               block_stride,
                               block_data,
                               i_send_buffer,
                               n_send_buffer,
                               send_buffer);

  /*
   * Data exchange
   */
  
  int mandatory_size = PDM_size_idx_from_stride (n_send_buffer, ptb->s_comm, ptb->comm);
  mandatory_size = PDM_MAX(PDM_size_idx_from_stride (n_send_buffer, ptb->s_comm, ptb->comm), mandatory_size);

  if (mandatory_size > 32) {
    PDM_MPI_Alltoallv_p2p_l(send_buffer,
                            n_send_buffer,
                            i_send_buffer,
                            mpi_type,
                            recv_buffer,
                            n_recv_buffer,
                            i_recv_buffer,
                            mpi_type,
                            ptb->comm);
  }

  else {
 
    if (ptb->p2p_factor < ptb->part_active_rank) {
     
      int *_i_send_buffer = NULL;
      int *_i_recv_buffer = NULL;
      PDM_malloc(_i_send_buffer, ptb->s_comm, int);
      PDM_malloc(_i_recv_buffer, ptb->s_comm, int);
     
      for (int i = 0; i < ptb->s_comm; i++) {
        _i_send_buffer[i] = (int) i_send_buffer[i];
        _i_recv_buffer[i] = (int) i_recv_buffer[i];
      }

      PDM_MPI_Alltoallv(send_buffer,
                        n_send_buffer,
                        _i_send_buffer,
                        mpi_type,
                        recv_buffer,
                        n_recv_buffer,
                        _i_recv_buffer,
                        mpi_type,
                        ptb->comm);

      PDM_free(_i_send_buffer); 
      PDM_free(_i_recv_buffer); 
    }

    else {

      PDM_MPI_Alltoallv_p2p_l(send_buffer,
                              n_send_buffer,
                              i_send_buffer,
                              mpi_type,
                              recv_buffer,
                              n_recv_buffer,
                              i_recv_buffer,
                              mpi_type,
                              ptb->comm);

    }
  } 

  _post_treatment_reverse(ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          recv_stride,
                          recv_buffer,
                          n_recv_buffer,
                          i_recv_buffer,
                          part_stride,
                          part_data);
  PDM_MPI_Type_free(&mpi_type);

  PDM_free(recv_buffer);
  PDM_free(send_buffer);

  PDM_free(i_send_buffer);
  PDM_free(i_recv_buffer);
  PDM_free(n_send_buffer);
  PDM_free(n_recv_buffer);
}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                 **part_stride,
       void                **part_data,
       int                 **block_stride,
       void                **block_data,
       int                  *request
)
{
  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;
  request_id %= ptb->max_exch_request;
  *request = request_id;

  assert(ptb->wait_status[request_id] == 2);

  PDM_malloc(ptb->i_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_buffer[request_id], ptb->s_comm, int);

  ptb->block_stride[request_id] = block_stride;
  ptb->block_data  [request_id] = block_data;

  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;
  ptb->comm_kind [request_id] = k_comm;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = NULL;
  size_t *tmp_i_recv_buffer = NULL;
  PDM_malloc(tmp_i_send_buffer, ptb->s_comm, size_t);
  PDM_malloc(tmp_i_recv_buffer, ptb->s_comm, size_t);

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);
  ptb->mpi_type[request_id] = mpi_type;

  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    tmp_i_send_buffer,
                    tmp_i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = ( tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1] ) * s_data_tot;
  int s_recv_buffer = ( tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1] ) * s_data_tot;

  if(k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    PDM_MPI_Win_allocate(s_send_buffer, sizeof(unsigned char), ptb->comm, &ptb->send_buffer[request_id], &ptb->win_send[request_id]);
    PDM_MPI_Win_allocate(s_recv_buffer, sizeof(unsigned char), ptb->comm, &ptb->recv_buffer[request_id], &ptb->win_recv[request_id]);

  } else {

    PDM_malloc(ptb->send_buffer[request_id], s_send_buffer, unsigned char);
    PDM_malloc(ptb->recv_buffer[request_id], s_recv_buffer, unsigned char);

  }

  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       tmp_i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  PDM_free(tmp_i_send_buffer);
  PDM_free(tmp_i_recv_buffer);

  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_P2P k_comm is not implemented yet\n");
    abort();
  } else if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    PDM_MPI_Ialltoallv(send_buffer,
                       n_send_buffer,
                       i_send_buffer,
                       mpi_type,
                       recv_buffer,
                       n_recv_buffer,
                       i_recv_buffer,
                       mpi_type,
                       ptb->comm,
                       &ptb->request_mpi[request_id]);
  } else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    // double t1 = PDM_MPI_Wtime();
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
    PDM_MPI_Get_ialltoallv(ptb->win_send[request_id],
                           ptb->win_recv[request_id],
                           send_buffer,
                           n_send_buffer,
                           i_send_buffer,
                           mpi_type,
                           recv_buffer,
                           n_recv_buffer,
                           i_recv_buffer,
                           mpi_type,
                           ptb->comm);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_MPI_Get_ialltoallv + fence dt = %12.5e \n", dt);

  }


  ptb->wait_status[request_id] = 0;
}



/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */

void
PDM_part_to_block_reverse_iexch
(
       PDM_part_to_block_t  *ptb,
 const PDM_mpi_comm_kind_t   k_comm,
       size_t                s_data,
       PDM_stride_t          t_stride,
       int                   cst_stride,
       int                  *block_stride,
       void                 *block_data,
       int                ***part_stride,
       void               ***part_data,
       int                  *request
)
{
  assert(ptb->enable_reverse == 1);

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;
  request_id %= ptb->max_exch_request;
  *request = request_id;

  assert(ptb->wait_status[request_id] == 2);

  PDM_malloc(ptb->i_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_buffer[request_id], ptb->s_comm, int);

  ptb->part_stride[request_id] = part_stride;
  ptb->part_data  [request_id] = part_data;


  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;
  ptb->comm_kind [request_id] = k_comm;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = NULL;
  size_t *tmp_i_recv_buffer = NULL;
  PDM_malloc(tmp_i_send_buffer, ptb->s_comm, size_t);
  PDM_malloc(tmp_i_recv_buffer, ptb->s_comm, size_t);

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);
  ptb->mpi_type[request_id] = mpi_type;

  /*
   * Exchange Stride and build buffer properties
   */
  int *send_stride = NULL; // Release in post-treatment
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_reverse_exchange(ptb,
                            s_data,
                            t_stride,
                            cst_stride,
                            block_stride,
                            tmp_i_send_buffer,
                            tmp_i_recv_buffer,
                            n_send_buffer,
                            n_recv_buffer,
                            &send_stride,
                            &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = (tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1]) * s_data_tot;
  int s_recv_buffer = (tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1]) * s_data_tot;

  if(k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    PDM_MPI_Win_allocate(s_send_buffer, sizeof(unsigned char), ptb->comm, &ptb->send_buffer[request_id], &ptb->win_send[request_id]);
    PDM_MPI_Win_allocate(s_recv_buffer, sizeof(unsigned char), ptb->comm, &ptb->recv_buffer[request_id], &ptb->win_recv[request_id]);

  } else {

    PDM_malloc(ptb->send_buffer[request_id], s_send_buffer, unsigned char);
    PDM_malloc(ptb->recv_buffer[request_id], s_recv_buffer, unsigned char);

  }

  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_reverse_send_buffer(ptb,
                               s_data,
                               t_stride,
                               cst_stride,
                               send_stride,
                               block_stride,
                               block_data,
                               tmp_i_send_buffer,
                               n_send_buffer,
                               send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  PDM_free(tmp_i_send_buffer);
  PDM_free(tmp_i_recv_buffer);

  if (k_comm == PDM_MPI_COMM_KIND_P2P) {
    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();
  } else if(k_comm == PDM_MPI_COMM_KIND_COLLECTIVE) {
    PDM_MPI_Ialltoallv(send_buffer,
                       n_send_buffer,
                       i_send_buffer,
                       mpi_type,
                       recv_buffer,
                       n_recv_buffer,
                       i_recv_buffer,
                       mpi_type,
                       ptb->comm,
                       &ptb->request_mpi[request_id]);
  } else if (k_comm == PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_P2P k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE) {

    printf ("Error PDM_part_to_block_iexch : "
            " PDM_MPI_COMM_KIND_WIN_SHARED_AND_NEIGHBOR_COLLECTIVE k_comm is not implemented yet\n");
    abort();

  } else if (k_comm == PDM_MPI_COMM_KIND_WIN_RMA) {

    // double t1 = PDM_MPI_Wtime();
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
    PDM_MPI_Get_ialltoallv(ptb->win_send[request_id],
                           ptb->win_recv[request_id],
                           send_buffer,
                           n_send_buffer,
                           i_send_buffer,
                           mpi_type,
                           recv_buffer,
                           n_recv_buffer,
                           i_recv_buffer,
                           mpi_type,
                           ptb->comm);
    // double dt = PDM_MPI_Wtime() - t1;
    // log_trace("PDM_MPI_Get_ialltoallv + fence dt = %12.5e \n", dt);

  }


  ptb->wait_status[request_id] = 0;
}




/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 * \return       Size of highest block
 */
int
PDM_part_to_block_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
  } else {
    int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
    assert(code == PDM_MPI_SUCCESS);
  }

  PDM_MPI_Type_free(&ptb->mpi_type[request_id]);

  ptb->wait_status[request_id] = 1;
  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  /*
   *  Post-treatment
   */
  size_t s_recv_buffer = ( ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1]) * s_data_tot;

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_send[request_id]);
    ptb->win_send[request_id] = PDM_MPI_WIN_NULL;
  } else {
    PDM_free(ptb->send_buffer  [request_id]);
  }
  PDM_free(ptb->n_send_buffer[request_id]);
  PDM_free(ptb->i_send_buffer[request_id]);
  PDM_free(ptb->n_recv_buffer[request_id]);
  PDM_free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  // recv_stride is free inside
  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     ptb->block_stride[request_id],
                                     ptb->block_data  [request_id]);

  /*
   * Free
   */

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_recv[request_id]);
    ptb->win_recv[request_id] = PDM_MPI_WIN_NULL;
  } else {
    PDM_free(ptb->recv_buffer  [request_id]);
  }
  ptb->recv_stride [request_id] = NULL;
  ptb->recv_buffer [request_id] = NULL;
  ptb->block_stride[request_id] = NULL;
  ptb->block_data  [request_id] = NULL;


  ptb->wait_status[request_id] = 2;

  return s_block_data;
}



/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 * \return       Size of highest block
 */
void
PDM_part_to_block_reverse_iexch_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_fence(0, ptb->win_send[request_id]);
    PDM_MPI_Win_fence(0, ptb->win_recv[request_id]);
  } else {
    int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
    assert(code == PDM_MPI_SUCCESS);
  }

  PDM_MPI_Type_free(&ptb->mpi_type[request_id]);

  ptb->wait_status[request_id] = 1;
  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];

  /*
   *  Post-treatment
   */
  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_send[request_id]);
    ptb->win_send[request_id] = PDM_MPI_WIN_NULL;
  } else {
    PDM_free(ptb->send_buffer  [request_id]);
  }
  PDM_free(ptb->n_send_buffer[request_id]);
  PDM_free(ptb->i_send_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;

  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  size_t *tmp_i_recv_buffer;
  PDM_malloc(tmp_i_recv_buffer, ptb->s_comm+1, size_t);
  tmp_i_recv_buffer[0] = 0;
  for(int i = 0; i < ptb->s_comm; ++i) {
    tmp_i_recv_buffer[i+1] = tmp_i_recv_buffer[i] + ptb->n_recv_buffer[request_id][i];
  }

  _post_treatment_reverse(ptb,
                          s_data,
                          t_stride,
                          cst_stride,
                          recv_stride,
                          recv_buffer,
                          ptb->n_recv_buffer[request_id],
                          tmp_i_recv_buffer,
                          ptb->part_stride  [request_id],
                          ptb->part_data    [request_id]);

  PDM_free(tmp_i_recv_buffer);
  PDM_free(ptb->n_recv_buffer[request_id]);
  PDM_free(ptb->i_recv_buffer[request_id]);
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  /*
   * Free
   */

  if(ptb->comm_kind[request_id] == PDM_MPI_COMM_KIND_WIN_RMA) {
    PDM_MPI_Win_free(&ptb->win_recv[request_id]);
    ptb->win_recv[request_id] = PDM_MPI_WIN_NULL;
  } else {
    PDM_free(ptb->recv_buffer  [request_id]);
  }
  ptb->recv_stride [request_id] = NULL;
  ptb->recv_buffer [request_id] = NULL;
  ptb->part_stride [request_id] = NULL;
  ptb->part_data   [request_id] = NULL;

  ptb->wait_status[request_id] = 2;

}

/**
 *
 * \brief Initialize a data exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   s_data       Data size
 * \param [in]   t_stride     Stride type
 * \param [in]   var_stride   Variable stride (size = n_part) only for \ref PDM_writer_STRIDE_VAR
 * \param [in]   cst_stride   Stride only for \ref PDM_writer_STRIDE_CST
 *
 * \return       Size of highest block
 *
 */
int
PDM_part_to_block_async_exch
(
 PDM_part_to_block_t *ptb,
 size_t               s_data,
 PDM_stride_t         t_stride,
 int                  cst_stride,
 int                **part_stride,
 void               **part_data
)
{

  if ((ptb->t_post == PDM_PART_TO_BLOCK_POST_MERGE) &&
      (t_stride ==  PDM_STRIDE_CST_INTERLACED)) {
    PDM_error(__FILE__, __LINE__, 0,"PDM_part_to_block_exch : PDM_writer_STRIDE_CST is not compatible PDM_writer_POST_MERGE post\n");
    abort ();
  }

  /*
   * Take next id for message and buffer
   */
  int request_id = ptb->next_request++;
  request_id %= ptb->max_exch_request;

  assert(ptb->wait_status[request_id] == 2);

  PDM_malloc(ptb->i_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->i_recv_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_send_buffer[request_id], ptb->s_comm, int);
  PDM_malloc(ptb->n_recv_buffer[request_id], ptb->s_comm, int);

  /* Store additionnal information necessary for post-process */
  ptb->s_data    [request_id] = s_data;
  ptb->t_stride  [request_id] = t_stride;
  ptb->cst_stride[request_id] = cst_stride;

  /* Short cut */
  int* i_send_buffer = ptb->i_send_buffer[request_id];
  int* i_recv_buffer = ptb->i_recv_buffer[request_id];
  int* n_send_buffer = ptb->n_send_buffer[request_id];
  int* n_recv_buffer = ptb->n_recv_buffer[request_id];

  size_t *tmp_i_send_buffer = NULL;
  size_t *tmp_i_recv_buffer = NULL;
  PDM_malloc(tmp_i_send_buffer, ptb->s_comm, size_t);
  PDM_malloc(tmp_i_recv_buffer, ptb->s_comm, size_t);

  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  PDM_MPI_Datatype mpi_type;
  PDM_MPI_Type_create_contiguous(s_data_tot, PDM_MPI_BYTE, &mpi_type);
  PDM_MPI_Type_commit(&mpi_type);
  ptb->mpi_type[request_id] = mpi_type;


  /*
   * Exchange Stride and build buffer properties
   */
  int *recv_stride = NULL; // Release in post-treatment
  _prepare_exchange(ptb,
                    s_data,
                    t_stride,
                    cst_stride,
                    part_stride,
                    tmp_i_send_buffer,
                    tmp_i_recv_buffer,
                    n_send_buffer,
                    n_recv_buffer,
                    &recv_stride);
  ptb->recv_stride[request_id] = recv_stride;

  /*
   * Data exchange
   */
  int s_send_buffer = ( tmp_i_send_buffer[ptb->s_comm - 1] + n_send_buffer[ptb->s_comm -1] ) * s_data_tot;
  int s_recv_buffer = ( tmp_i_recv_buffer[ptb->s_comm - 1] + n_recv_buffer[ptb->s_comm -1] ) * s_data_tot;

  PDM_malloc(ptb->send_buffer[request_id], s_send_buffer, unsigned char);
  PDM_malloc(ptb->recv_buffer[request_id], s_recv_buffer, unsigned char);

  /* Shortcut */
  unsigned char *send_buffer = ptb->send_buffer[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  assert (send_buffer != NULL);
  assert (recv_buffer != NULL);

  _prepare_send_buffer(ptb,
                       s_data,
                       t_stride,
                       cst_stride,
                       part_stride,
                       part_data,
                       tmp_i_send_buffer,
                       n_send_buffer,
                       send_buffer);

  /*
   * Copy back
   */
  for (int i = 0; i < ptb->s_comm; i++) {
    i_send_buffer[i] = (int ) tmp_i_send_buffer[i];
    i_recv_buffer[i] = (int ) tmp_i_recv_buffer[i];
  }

  PDM_free(tmp_i_send_buffer);
  PDM_free(tmp_i_recv_buffer);

  PDM_MPI_Ialltoallv(send_buffer,
                     n_send_buffer,
                     i_send_buffer,
                     mpi_type,
                     recv_buffer,
                     n_recv_buffer,
                     i_recv_buffer,
                     mpi_type,
                     ptb->comm,
                     &ptb->request_mpi[request_id]);

  ptb->wait_status[request_id] = 0;

  return request_id;
}


/**
 *
 * \brief Wait for an exchange
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 *
 */
void
PDM_part_to_block_async_wait
(
 PDM_part_to_block_t *ptb,
 int                  request_id
)
{
  // printf("PDM_part_to_block_async_wait::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 0);

  int code = PDM_MPI_Wait(&ptb->request_mpi[request_id]);
  assert(code == PDM_MPI_SUCCESS);

  ptb->wait_status[request_id] = 1;
  PDM_MPI_Type_free(&ptb->mpi_type[request_id]);

}

/**
 *
 * \brief Get the raw exchange buffer and stride and deallocate memory
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_get_raw
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  // printf("PDM_part_to_block_asyn_get_raw::request_id::%d \n", request_id);

  assert(ptb->wait_status[request_id] == 1);

  PDM_free(ptb->send_buffer  [request_id]);
  PDM_free(ptb->n_send_buffer[request_id]);
  PDM_free(ptb->i_send_buffer[request_id]);
  PDM_free(ptb->n_recv_buffer[request_id]);
  PDM_free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;

  /* Mv pointer */
  *block_stride = ptb->recv_stride[request_id];
  *block_data   = ptb->recv_buffer[request_id];

  /* Nulliffy - User is now owner of the excahnge data */
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  return ptb->tn_recv_data;

}


/**
 *
 * \brief Post-treatment of the recv buffer
 *
 * \param [in]   ptb          Part to block structure
 * \param [in]   request_id   Internal id of the current exchange
 * \param [out]  block_stride Block stride
 * \param [out]  block_data   Block data
 *
 */
int
PDM_part_to_block_asyn_post_treatment
(
 PDM_part_to_block_t *ptb,
 int                  request_id,
 int                **block_stride,
 void               **block_data
)
{
  assert(ptb->wait_status[request_id] == 1);

  size_t       s_data     = ptb->s_data    [request_id];
  PDM_stride_t t_stride   = ptb->t_stride  [request_id];
  int          cst_stride = ptb->cst_stride[request_id];
  int s_data_tot = s_data;
  if(t_stride == PDM_STRIDE_CST_INTERLACED) {
    s_data_tot = s_data * cst_stride;
  }

  // size_t s_send_buffer = ptb->i_send_buffer[request_id][ptb->s_comm - 1] + ptb->n_send_buffer[request_id][ptb->s_comm -1];
  size_t s_recv_buffer = ( ptb->i_recv_buffer[request_id][ptb->s_comm - 1] + ptb->n_recv_buffer[request_id][ptb->s_comm -1] )*s_data_tot;

  PDM_free(ptb->send_buffer  [request_id]);
  PDM_free(ptb->n_send_buffer[request_id]);
  PDM_free(ptb->i_send_buffer[request_id]);
  PDM_free(ptb->n_recv_buffer[request_id]);
  PDM_free(ptb->i_recv_buffer[request_id]);

  ptb->send_buffer  [request_id] = NULL;
  ptb->n_send_buffer[request_id] = NULL;
  ptb->i_send_buffer[request_id] = NULL;
  ptb->n_recv_buffer[request_id] = NULL;
  ptb->i_recv_buffer[request_id] = NULL;


  int           *recv_stride = ptb->recv_stride[request_id];
  unsigned char *recv_buffer = ptb->recv_buffer[request_id];

  // recv_stride is free inside
  int s_block_data = _post_treatment(ptb,
                                     s_data,
                                     t_stride,
                                     cst_stride,
                                     recv_stride,
                                     recv_buffer,
                                     s_recv_buffer,
                                     block_stride,
                                     block_data);

  /*
   * Free
   */
  PDM_free(ptb->recv_buffer  [request_id]);
  ptb->recv_stride[request_id] = NULL;
  ptb->recv_buffer[request_id] = NULL;

  ptb->wait_status[request_id] = 2;

  return s_block_data;
}


/**
 *
 * \brief Free a part to block structure
 *
 * \param [inout] ptb         Part to block structure
 *
 * \return       NULL
 */
PDM_part_to_block_t *
PDM_part_to_block_free
(
 PDM_part_to_block_t *ptb
)
{

  if (ptb->active_ranks != NULL) {
    PDM_free(ptb->active_ranks);
    ptb->active_ranks = NULL;
  }
  if (ptb->dest_proc != NULL) {
    PDM_free(ptb->dest_proc);
    ptb->dest_proc = NULL;
  }
  if (ptb->data_distrib_index != NULL) {
    PDM_free(ptb->data_distrib_index);
    ptb->data_distrib_index = NULL;
  }
  if (ptb->i_send_data != NULL) {
    PDM_free(ptb->i_send_data);
    ptb->i_send_data = NULL;
  }
  if (ptb->i_recv_data != NULL) {
    PDM_free(ptb->i_recv_data);
    ptb->i_recv_data = NULL;
  }
  if (ptb->n_send_data != NULL) {
    PDM_free(ptb->n_send_data);
    ptb->n_send_data = NULL;
  }
  if (ptb->n_recv_data != NULL) {
    PDM_free(ptb->n_recv_data);
    ptb->n_recv_data = NULL;
  }
  if (ptb->sorted_recv_gnum != NULL) {
    PDM_free(ptb->sorted_recv_gnum);
    ptb->sorted_recv_gnum = NULL;
  }

  if (ptb->enable_reverse == 1 && ptb->idx_partial != NULL) {
    PDM_free(ptb->idx_partial);
    ptb->idx_partial = NULL;
  }
  if (ptb->order != NULL) {
    PDM_free(ptb->order);
    ptb->order = NULL;
  }

  // if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING) && (ptb->block_gnum != NULL)) {
  //  PDM_free(ptb->block_gnum);
  //   ptb->block_gnum = NULL;
  // }
  if ((ptb->t_post != PDM_PART_TO_BLOCK_POST_NOTHING      ) &&
      (ptb->t_post != PDM_PART_TO_BLOCK_POST_MERGE_UNIFORM) && (ptb->block_gnum != NULL)) {
    PDM_free(ptb->block_gnum);
    ptb->block_gnum = NULL;
  }
  if (ptb->block_gnum_count != NULL) {
    PDM_free(ptb->block_gnum_count);
    ptb->block_gnum_count=NULL;
  }

  if (ptb->weight_g != NULL) {
    for (int i = 0; i < ptb->n_part; i++) {
      PDM_free(ptb->weight_g[i]);
    }
    PDM_free(ptb->weight_g);
  }

  /* This one check if all buffer has been correctly move or delete */
  for(int i_req = 0; i_req < ptb->max_exch_request; ++i_req) {
    assert(ptb->send_buffer  [i_req] == NULL);
    assert(ptb->recv_buffer  [i_req] == NULL);
    assert(ptb->recv_stride  [i_req] == NULL);
    assert(ptb->n_send_buffer[i_req] == NULL);
    assert(ptb->i_send_buffer[i_req] == NULL);
    assert(ptb->n_recv_buffer[i_req] == NULL);
    assert(ptb->i_recv_buffer[i_req] == NULL);
    assert(ptb->block_stride [i_req] == NULL);
    assert(ptb->block_data   [i_req] == NULL);
    assert(ptb->part_stride  [i_req] == NULL);
    assert(ptb->part_data    [i_req] == NULL);
    assert(ptb->win_send     [i_req] == PDM_MPI_WIN_NULL);
    assert(ptb->win_recv     [i_req] == PDM_MPI_WIN_NULL);
    assert(ptb->mpi_type     [i_req] == PDM_MPI_DATATYPE_NULL);
  }

  PDM_free(ptb->s_data       );
  PDM_free(ptb->t_stride     );
  PDM_free(ptb->cst_stride   );
  PDM_free(ptb->wait_status  );
  PDM_free(ptb->request_mpi  );
  PDM_free(ptb->send_buffer  );
  PDM_free(ptb->recv_buffer  );
  PDM_free(ptb->recv_stride  );
  PDM_free(ptb->n_send_buffer);
  PDM_free(ptb->i_send_buffer);
  PDM_free(ptb->n_recv_buffer);
  PDM_free(ptb->i_recv_buffer);
  PDM_free(ptb->block_stride );
  PDM_free(ptb->block_data   );
  PDM_free(ptb->part_stride  );
  PDM_free(ptb->part_data    );

  PDM_free(ptb->comm_kind);
  PDM_free(ptb->win_send );
  PDM_free(ptb->win_recv );
  PDM_free(ptb->mpi_type );

  PDM_free(ptb);

  n_ptb--;
  if (n_ptb == 0) {
    PDM_timer_free(t_timer[MALLOC_ACTIVE_RANKS   ]);
    PDM_timer_free(t_timer[GENERATE_DISTRIB      ]);
    PDM_timer_free(t_timer[BINARY_SEARCH         ]);
    PDM_timer_free(t_timer[CREATE_EXCHANGE       ]);
    PDM_timer_free(t_timer[BLOCK_POST            ]);
    PDM_timer_free(t_timer[GLOBAL_WEIGHTS        ]);
    PDM_timer_free(t_timer[CREATE_FROM_DISTRIB   ]);
    PDM_timer_free(t_timer[CREATE_GEOM           ]);
    PDM_timer_free(t_timer[DATA_EXCHANGE         ]);
  }

  return NULL;
}


/**
 *
 * \brief Return block distribution index
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Distribution (size = communicator size + 1)
 */

PDM_g_num_t *
PDM_part_to_block_distrib_index_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->data_distrib_index;
}


/**
 *
 * \brief Return processus destination
 *
 * \param [in] ptb         Part to block structure
 *
 * \return  Destination (size = sum of partition elements)
 */

PDM_l_num_t *
PDM_part_to_block_destination_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->dest_proc;
}

PDM_g_num_t*
PDM_part_to_block_adapt_partial_block_to_block
(
 PDM_part_to_block_t  *ptb,
 int                 **block_n,
 PDM_g_num_t           n_g_block
)
{
  PDM_g_num_t *_block_distrib_idx;
  PDM_malloc(_block_distrib_idx, ptb->s_comm + 1, PDM_g_num_t);

  for (int i = 0; i < ptb->s_comm + 1; i++) {
    _block_distrib_idx[i] = ptb->data_distrib_index[i];
  }
  int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);

  int block_n_elt_tot = _block_distrib_idx[ptb->i_rank+1] - _block_distrib_idx[ptb->i_rank];
  int* block_n_tmp = PDM_array_zeros_int(block_n_elt_tot);

  int* _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt; i1++) {
    int i = (int) (ptb->block_gnum[i1] - _block_distrib_idx[ptb->i_rank] - 1);
    // printf(" ptb->block_gnum[%i] = %i --> %i (%i)\n", i1, ptb->block_gnum[i1], _block_distrib_idx[ptb->i_rank], i);
    block_n_tmp[i] = _block_n[i1];
  }

  PDM_realloc(*block_n ,*block_n , block_n_elt_tot,int);
  _block_n = *block_n;
  for (int i1 = 0; i1 < block_n_elt_tot; i1++) {
    _block_n[i1] = block_n_tmp[i1];
  }

  PDM_free(block_n_tmp);

  PDM_g_num_t old_max = _block_distrib_idx[ptb->s_comm];
  PDM_g_num_t new_max = n_g_block;
  int diff_last = (int) (new_max - old_max);

  _block_distrib_idx[ptb->s_comm] = new_max;

  if (ptb->i_rank == (ptb->s_comm - 1)) {
    int new_size = block_n_elt_tot + diff_last;
    PDM_realloc(*block_n ,*block_n , new_size,int);
    _block_n = *block_n;
    for (int i = block_n_elt_tot; i < new_size; i++) {
      _block_n[i] = 0;
    }
  }

  return _block_distrib_idx;
}





/**
 *
 * \brief Return global weights of element in the current process
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Global weights
 *
 */

double **
PDM_part_to_block_global_weight_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->weight_g;
}



/**
 *
 * \brief Get number of MPI ranks
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return  Number of MPI ranks
 *
 */

int
PDM_part_to_block_n_ranks_get
(
 PDM_part_to_block_t *ptb
)
{
  return ptb->s_comm;
}


/**
 *
 * \brief Return total number of element in the current process (summed over all partitions)
 *
 * \param [in]   ptb          Part to block structure
 *
 * \return Total number of element in the current process
 *
 */

int
PDM_part_to_block_n_elt_proc_get
(
 PDM_part_to_block_t *ptb
 )
 {
  return ptb->n_elt_proc;
 }

#undef _MIN
#undef _MAX

#ifdef __cplusplus
}
#endif /* __cplusplus */
