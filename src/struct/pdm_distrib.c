/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#include "pdm_distrib.h"
#include "pdm_binary_search.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 * \brief   Evaluate a distribution array.
 *
 * \param [in]  n_ranges     Number of ranges in the distribution
 * \param [in]  distribution Number of elements associated to each range of the distribution
 * \param [in]  optim        Optimal count in each range
 *
 * \return  a fit associated to the distribution. If fit = 0, distribution is perfect.
 *
 */
static double
_evaluate_distribution(int          n_ranges,
                       double      *distribution,
                       double       optim)
{
  int  i;
  double  d_low = 0, d_up = 0, fit = 0;

  /*
     d_low is the max gap between the distribution count and the optimum when
     distribution is lower than optimum.
     d_up is the max gap between the distribution count and the optimum when
     distribution is greater than optimum.
  */

  for (i = 0; i < n_ranges; i++) {

    if (distribution[i] > optim)
      d_up = PDM_MAX(d_up, distribution[i] - optim);
    else
      d_low = PDM_MAX(d_low, optim - distribution[i]);

  }

  fit = (d_up + d_low) / optim;

#if 0 && defined(DEBUG) && !defined(NDEBUG)
  if (cs_glob_rank_id <= 0)
    PDM_printf( "<DISTRIBUTION EVALUATION> optim: %g, fit: %g\n",
               optim, fit);
#endif

  return  fit;
}


/**
 * \brief Define a global distribution associated to a sampling array i.e. count
 * the number of elements in each range.
 *
 *   \param [in]    dim           2D or 3D
 *   \param [in]    n_ranks       number of ranks (= number of ranges)
 *   \param [in]    gsum_weight   global sum of all weightings
 *   \param [in]    n_codes       local number of Hilbert codes
 *   \param [in]    hilbert_codes local list of Hilbert codes to distribute
 *   \param [in]    weight        weighting related to each code
 *   \param [in]    order         ordering array
 *   \param [in]    sampling      sampling array
 *   \param [inout] c_freq        pointer to the cumulative frequency array
 *   \param [inout] g_distrib     pointer to a distribution array
 *   \param [in]    comm          mpi communicator
 */

static void
_define_rank_distrib(const int             sampling_factor,
                     const int             n_part,
                     const int            *n_elt,
                     const PDM_g_num_t   **gnum_elt,
                     const double        **weight,
                     const int             n_ranks,
                     const double          gsum_weight,
                     const PDM_g_num_t     sampling[],
                           double          cfreq[],
                           double          g_distrib[],
                           PDM_MPI_Comm    comm)
{
  int  id, rank_id;

  const int  n_samples       = sampling_factor * n_ranks;

  /* Initialization */
  double   *l_distrib = (double  *) malloc (n_samples * sizeof(double));

  for (id = 0; id < n_samples; id++) {
    l_distrib[id] = 0;
    g_distrib[id] = 0;
  }

  for (int i = 0; i < n_part; i++) {
    for (int j = 0; j < n_elt[i]; j++) {
      PDM_g_num_t _gnum_elt = PDM_ABS(gnum_elt[i][j]) - 1;
      int i_sample = PDM_binary_search_gap_long(_gnum_elt, sampling, n_samples + 1);
      l_distrib[i_sample] += weight[i][j];
    }
  }

  /* Define the global distribution */
  PDM_MPI_Allreduce(l_distrib, g_distrib, n_samples,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  free(l_distrib);

  /* Define the cumulative frequency related to g_distribution */
  cfreq[0] = 0.;
  for (id = 0; id < n_samples; id++) {
    double _g_distrib  = (double)g_distrib[id];
    double _gsum_weight = (double)gsum_weight;
    cfreq[id+1] = cfreq[id] + _g_distrib/_gsum_weight;
  }
  cfreq[n_samples] = 1.0;

#if 0 && defined(DEBUG) && !defined(DEBUG) /* For debugging purpose only */

  if (cs_glob_rank_id <= 0) {

    FILE  *dbg_file = NULL;
    int  len;
    static int  loop_id1 = 0;

    len = strlen("DistribOutput_l.dat")+1+2;
    char  *rfilename = (char *) malloc (len * sizeof(char));
    sprintf(rfilename, "DistribOutput_l%02d.dat", loop_id1);

    loop_id1++;

    dbg_file = fopen(rfilename, "w");

    fprintf(dbg_file,
            "# Sample_id  |  OptCfreq  |  Cfreq  |  Sampling  |"
            "Global Distrib\n");
    for (i = 0; i < n_samples; i++)
      fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
              i, (double)i/(double)n_samples, cfreq[i],
              (double)(sampling[i]), distrib[i]);
    fprintf(dbg_file, "%8d %15.5f %15.10f %15.10f %10u\n",
            i, 1.0, 1.0, 1.0, 0);

    fclose(dbg_file);
    free(rfilename);

  }

#endif /* debugging output */

  /* Convert global distribution from n_samples to n_ranks */

  for (rank_id = 0; rank_id < n_ranks; rank_id++) {

    double   sum = 0.;
    int   shift = rank_id * sampling_factor;

    for (id = 0; id < sampling_factor; id++) {
      sum += g_distrib[shift + id];
    }
    g_distrib[rank_id] = sum;

  } /* End of loop on ranks */

#if 0 && defined(DEBUG) && !defined(NDEBUG) /* Sanity check in debug */
  {
    PDM_g_num_t   sum = 0;
    for (rank_id = 0; rank_id < n_ranks; rank_id++)
      sum += g_distrib[rank_id];

    if (sum != gsum_weight)
      PDM_error(__FILE__, __LINE__, 0,
                "Error while computing global distribution.\n"
                "sum = %u and gsum_weight = %u\n",
                sum, gsum_weight);
    exit(1);
  }
#endif /* sanity check */

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
void
PDM_distrib_compute
(
 const int           dnelt,
       PDM_g_num_t  *elt_distrib,
       int           offset,
 const PDM_MPI_Comm  comm
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Compute distribution for element */

  // PDM_g_num_t* elt_distrib = (PDM_g_num_t *) malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t  _dnelt      = (PDM_g_num_t) dnelt;

  PDM_MPI_Allgather((void *) &_dnelt,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&elt_distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  elt_distrib[0] = 1+offset;
  for (int i = 1; i < n_rank+1; i++) {
    elt_distrib[i] +=  elt_distrib[i-1];
  }

  /* Verbose */
  if (1 == 0) {
    PDM_printf("elt_distrib : "PDM_FMT_G_NUM,  elt_distrib[0]);
    for (int i = 1; i < n_rank+1; i++) {
      PDM_printf(" "PDM_FMT_G_NUM, elt_distrib[i]);
    }
    PDM_printf("\n");
  }
}

/**
 * \brief Compute distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
PDM_g_num_t*
PDM_compute_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const int              dn_entity
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t* dentity_proc = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  /*
   * Exchange
   */

  PDM_g_num_t _dn_entity = (PDM_g_num_t) dn_entity;
  PDM_MPI_Allgather((void *) &_dn_entity,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&dentity_proc[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    comm);

  dentity_proc[0] = 0;
  for (int i = 1; i < n_rank+1; i++) {
    dentity_proc[i] = dentity_proc[i] + dentity_proc[i-1];
  }

  return dentity_proc;
}


/**
 * \brief Compute uniform distribution distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
int
PDM_compute_uniform_dn_entity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t step      = n_g_entity/n_rank;
  PDM_g_num_t remainder = n_g_entity%n_rank;

  int dn_elmt = step;
  if (i_rank < remainder) { /* Distribute the remainder */
    dn_elmt += 1;
  }

  return dn_elmt;
}

/**
 * \brief Compute uniform distribution distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
PDM_g_num_t*
PDM_compute_uniform_entity_distribution
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t      n_g_entity
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t* dentity_proc = (PDM_g_num_t *) malloc( (n_rank+1) * sizeof(PDM_g_num_t));

  PDM_g_num_t step      = n_g_entity/n_rank;
  PDM_g_num_t remainder = n_g_entity%n_rank;

  dentity_proc[0] = 0;
  for (int i = 1; i < n_rank + 1; i++) {
    dentity_proc[i] = step;
    const int i1 = i - 1;
    if (i1 < remainder) { /* Distribute the remainder */
      dentity_proc[i]  += 1;
    }
  }

  for (int i = 1; i < n_rank + 1; i++) {
    dentity_proc[i] += dentity_proc[i-1];
  }

  return dentity_proc;
}


/**
 * \brief Compute uniform distribution distribution from dNelmt
 *
 * \param [in]     elt_distrib          Distribution of elements on processes
 * \param [in]     dnelt                Number of element on current process
 * \param [in]     comm                 MPI Communicator
 */
PDM_g_num_t*
PDM_compute_uniform_entity_distribution_from_partition
(
 const PDM_MPI_Comm     comm,
 const int              n_part,
 const int             *n_elmts,
 const PDM_g_num_t    **ln_to_gn
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Compute the max
   */
  PDM_g_num_t _id_max = 0;
  PDM_g_num_t n_g_entity = 0;

  for(int i_part = 0; i_part < n_part; ++i_part) {
    for(int i_elmt = 0; i_elmt < n_elmts[i_part]; ++i_elmt) {
      _id_max = PDM_MAX (_id_max, ln_to_gn[i_part][i_elmt]);
    }
  }

  // double t1 = PDM_MPI_Wtime();
  PDM_MPI_Allreduce (&_id_max,
                     &n_g_entity,
                     1,
                     PDM__PDM_MPI_G_NUM,
                     PDM_MPI_MAX,
                     comm);
  // double t2 = PDM_MPI_Wtime();
  // double dt = t2-t1;
  // printf("[%i] dt = %12.5e \n ", i_rank, dt);

  return PDM_compute_uniform_entity_distribution(comm, n_g_entity);
}


void
PDM_distrib_weight
(
  const int            sampling_factor,
  const int            n_active_ranks,
  const int            n_part,
  const int           *n_elmts,
  const PDM_g_num_t  **ln_to_gn,
  const double       **weight,
        PDM_MPI_Comm   comm
)
{

  PDM_g_num_t _id_max     = 0;
  PDM_g_num_t _id_max_max = 0;

  for (int i = 0; i < n_part; i++) {
    for (int j = 0; j < n_elmts[i]; j++) {
      PDM_g_num_t gnum = PDM_ABS(ln_to_gn[i][j]);
      _id_max = PDM_MAX (_id_max, gnum);
    }
  }

  PDM_MPI_Allreduce (&_id_max,
                     &_id_max_max,
                     1,
                     PDM__PDM_MPI_G_NUM,
                     PDM_MPI_MAX,
                     comm);


  const int n_samples = sampling_factor * n_active_ranks;

  PDM_g_num_t *sampling = malloc(sizeof(PDM_g_num_t) * (n_samples + 1));

  double  lsum_weight = 0.;
  for (int i = 0; i < n_part; i++) {
    for (int j = 0; j < n_elmts[i]; j++) {
      lsum_weight += weight[i][j];
    }
  }

  double  gsum_weight = 0.;
  PDM_MPI_Allreduce(&lsum_weight, &gsum_weight, 1,
                    PDM_MPI_DOUBLE, PDM_MPI_SUM, comm);

  double optim = gsum_weight / n_active_ranks;

  /* Define a naive sampling (uniform distribution) */
  PDM_g_num_t _n_sample_data = _id_max_max / n_samples;
  PDM_g_num_t _samplerest    = _id_max_max % n_samples;

  sampling[0] = 0;
  int k = 0;
  for (int i = 0; i < n_samples; i++) {
    sampling[i+1] = sampling[i];
    sampling[i+1] += _n_sample_data;
    if (k < _samplerest) {
      sampling[i+1] += 1;
      k += 1;
    }
  }

  /* Define the distribution associated to the current sampling array */

  double *distrib = (double *) malloc (sizeof(double) * n_samples);
  double  *cfreq = (double *) malloc (sizeof(double) * (n_samples + 1));

  _define_rank_distrib(sampling_factor,
                       n_part,
                       n_elmts,
                       ln_to_gn,
                       weight,
                       n_active_ranks,
                       gsum_weight,
                       sampling,
                       cfreq,
                       distrib,
                       comm);

  /* Initialize best choice */

  double fit = _evaluate_distribution(n_active_ranks, distrib, optim);
  double best_fit = fit;

  PDM_g_num_t  *best_sampling = (PDM_g_num_t  *) malloc (sizeof(PDM_g_num_t) * (n_samples + 1));

  for (int i = 0; i < (n_samples + 1); i++) {
    best_sampling[i] = sampling[i];
  }





}



#ifdef __cplusplus
}
#endif /* __cplusplus */
