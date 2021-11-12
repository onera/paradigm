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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_gnum.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
dconnectivity_to_extract_dconnectivity
(
 const PDM_MPI_Comm    comm,
       int             n_selected_entity1,
       PDM_g_num_t    *select_entity1,
       PDM_g_num_t    *entity1_distribution,
       int            *dentity1_entity2_idx,
       PDM_g_num_t    *dentity1_entity2
)
{

  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * Create implcite global numbering
   */
  int dn_entity1 = entity1_distribution[i_rank+1] - entity1_distribution[i_rank];
  // PDM_g_num_t* entity1_ln_to_gn = malloc(dn_entity1 * sizeof(PDM_g_num_t));
  // for(int i = 0; i < dn_entity1; ++i) {
  //   entity1_ln_to_gn[i] = entity1_distribution[i_rank] + i + 1;
  // }

  /*
   *  Create global numbering from parent
   */
  PDM_gen_gnum_t* gen_gnum_entity1 = PDM_gnum_create(3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);


  PDM_gnum_set_from_parents (gen_gnum_entity1,
                             0,
                             n_selected_entity1,
                             select_entity1);

  PDM_gnum_compute (gen_gnum_entity1);


  PDM_g_num_t* extract_entity1_ln_to_gn = PDM_gnum_get (gen_gnum_entity1, 0);

  if(1 == 1) {
    PDM_log_trace_array_long(select_entity1          , n_selected_entity1, "select_entity1:: ");
    PDM_log_trace_array_long(extract_entity1_ln_to_gn, n_selected_entity1, "extract_entity1_ln_to_gn:: ");
  }


  PDM_gnum_free(gen_gnum_entity1);
  // free(entity1_ln_to_gn);

  /*
   * Caution we need a result independant of parallelism
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(entity1_distribution,
                               (const PDM_g_num_t **) &select_entity1,
                                                      &n_selected_entity1,
                                                      1,
                                                      comm);

  /*
   * Prepare data
   */
  int* dentity1_entity2_n = (int *) malloc( sizeof(int) * dn_entity1);
  for(int i_elmt = 0; i_elmt < dn_entity1; ++i_elmt){
    dentity1_entity2_n[i_elmt] = dentity1_entity2_idx[i_elmt+1] - dentity1_entity2_idx[i_elmt];
  }

  /*
   * Exchange
   */
  int**         tmp_dextract_entity1_entity2_n;
  PDM_g_num_t** tmp_dextract_entity1_entity2;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dentity1_entity2_n,
             (void *  )   dentity1_entity2,
             (int  ***)  &tmp_dextract_entity1_entity2_n,
             (void ***)  &tmp_dextract_entity1_entity2);

  int**         dextract_entity1_entity2_n = tmp_dextract_entity1_entity2_n[0];
  PDM_g_num_t** dextract_entity1_entity2   = tmp_dextract_entity1_entity2[0];
  free(tmp_dextract_entity1_entity2_n);
  free(tmp_dextract_entity1_entity2);
  free(dentity1_entity2_n);


  PDM_block_to_part_free(btp);

  free(dextract_entity1_entity2_n);
  free(dextract_entity1_entity2  );

  free(extract_entity1_ln_to_gn);
}





/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{

  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int           dn_cell;
  int           dn_face;
  int           dn_vtx;
  int           n_face_group;
  PDM_g_num_t  *dface_cell      = NULL;
  int          *dface_vtx_idx   = NULL;
  PDM_g_num_t  *dface_vtx       = NULL;
  double       *dvtx_coord      = NULL;
  int          *dface_group_idx = NULL;
  PDM_g_num_t  *dface_group     = NULL;
  int           dface_vtx_l;
  int           dface_group_l;

  PDM_dcube_t* dcube = PDM_dcube_gen_init(comm,
                                          n_vtx_seg,
                                          length,
                                          0.,
                                          0.,
                                          0.,
                                          PDM_OWNERSHIP_KEEP);

  PDM_dcube_gen_dim_get(dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);

  PDM_dcube_gen_data_get(dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face);
  /*
   *  Choice of extraction
   */
  dconnectivity_to_extract_dconnectivity(comm,
                                         dface_group_idx[n_face_group],
                                         dface_group,
                                         face_distribution,
                                         dface_vtx_idx,
                                         dface_vtx);


  free(face_distribution);

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
