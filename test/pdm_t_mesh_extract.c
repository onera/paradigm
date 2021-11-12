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
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"

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
       PDM_g_num_t    *dentity1_entity2,
       PDM_g_num_t   **extract_entity1_distribution,
       PDM_g_num_t   **extract_entity2_distribution,
       int           **dextract_entity1_entity2_idx,
       PDM_g_num_t   **dextract_entity1_entity2,
       PDM_g_num_t   **dparent_entity1_g_num,
       PDM_g_num_t   **dparent_entity2_g_num
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

  if(0 == 1) {
    PDM_log_trace_array_long(select_entity1          , n_selected_entity1, "select_entity1:: ");
    PDM_log_trace_array_long(extract_entity1_ln_to_gn, n_selected_entity1, "extract_entity1_ln_to_gn:: ");
  }

  PDM_gnum_free(gen_gnum_entity1);

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
  int**         tmp_pextract_entity1_entity2_n;
  PDM_g_num_t** tmp_pextract_entity1_entity2;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dentity1_entity2_n,
             (void *  )   dentity1_entity2,
             (int  ***)  &tmp_pextract_entity1_entity2_n,
             (void ***)  &tmp_pextract_entity1_entity2);

  int*         pextract_entity1_entity2_n = tmp_pextract_entity1_entity2_n[0];
  PDM_g_num_t* pextract_entity1_entity2   = tmp_pextract_entity1_entity2[0];
  free(tmp_pextract_entity1_entity2_n);
  free(tmp_pextract_entity1_entity2);
  free(dentity1_entity2_n);

  PDM_block_to_part_free(btp);

  /*
   *  Remap inside true block to ensure parallelism independant
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                                      &extract_entity1_ln_to_gn,
                                                      NULL,
                                                      &n_selected_entity1,
                                                      1,
                                                      comm);


  int*         _dextract_entity1_entity2_n = NULL;
  PDM_g_num_t* _dextract_entity1_entity2   = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          1,
                          &pextract_entity1_entity2_n,
                (void **) &pextract_entity1_entity2,
                          &_dextract_entity1_entity2_n,
                (void **) &_dextract_entity1_entity2);

  PDM_g_num_t* _dparent_entity1_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                (void **) &select_entity1,
                          NULL,
                (void **) &_dparent_entity1_g_num);

  free(pextract_entity1_entity2_n);
  free(pextract_entity1_entity2  );

  /*
   * Post-Treatment
   */
  PDM_g_num_t* distrib_extract_entity1 = PDM_part_to_block_distrib_index_get(ptb);
  int dn_extract_entity1 = distrib_extract_entity1[i_rank+1] - distrib_extract_entity1[i_rank];

  int* _dextract_entity1_entity2_idx = malloc( (dn_extract_entity1 + 1) * sizeof(int));
  _dextract_entity1_entity2_idx[0] = 0;
  for(int i = 0; i < dn_extract_entity1; ++i) {
    _dextract_entity1_entity2_idx[i+1] = _dextract_entity1_entity2_idx[i] + _dextract_entity1_entity2_n[i];
  }
  free(_dextract_entity1_entity2_n);

  /*
   *  Create global numbering from parent entity2
   */
  PDM_gen_gnum_t* gen_gnum_entity2 = PDM_gnum_create(3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);


  PDM_gnum_set_from_parents (gen_gnum_entity2,
                             0,
                             _dextract_entity1_entity2_idx[dn_extract_entity1],
                             _dextract_entity1_entity2);

  PDM_gnum_compute (gen_gnum_entity2);


  PDM_g_num_t* extract_entity2_ln_to_gn = PDM_gnum_get (gen_gnum_entity2, 0);

  if(0 == 1) {
    PDM_log_trace_array_long(_dextract_entity1_entity2, _dextract_entity1_entity2_idx[dn_extract_entity1], "select_entity2:: ");
    PDM_log_trace_array_long(extract_entity2_ln_to_gn, _dextract_entity1_entity2_idx[dn_extract_entity1], "extract_entity2_ln_to_gn:: ");
  }


  PDM_gnum_free(gen_gnum_entity2);


  PDM_part_to_block_t *ptb_entity2 = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                              1.,
                                                              &extract_entity2_ln_to_gn,
                                                              NULL,
                                                              &_dextract_entity1_entity2_idx[dn_extract_entity1],
                                                              1,
                                                              comm);


  PDM_g_num_t* _dparent_entity2_g_num = NULL;
  PDM_part_to_block_exch (ptb_entity2,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
                (void **) &_dextract_entity1_entity2,
                          NULL,
                (void **) &_dparent_entity2_g_num);


  PDM_g_num_t* distrib_extract_entity2 = PDM_part_to_block_distrib_index_get(ptb_entity2);
  int dn_extract_entity2 = distrib_extract_entity2[i_rank+1] - distrib_extract_entity2[i_rank];

  if(0 == 1) {
    PDM_log_trace_array_long(_dparent_entity2_g_num, dn_extract_entity2, "_dparent_entity2_g_num:: ");
  }

  /*
   * Update connectivity with the new ones
   */
  for(int i = 0; i < _dextract_entity1_entity2_idx[dn_extract_entity1]; ++i) {
    _dextract_entity1_entity2[i] = extract_entity2_ln_to_gn[i];
  }

  free(extract_entity1_ln_to_gn);
  free(extract_entity2_ln_to_gn);
  *dextract_entity1_entity2_idx = _dextract_entity1_entity2_idx;
  *dextract_entity1_entity2     = _dextract_entity1_entity2;
  *dparent_entity1_g_num        = _dparent_entity1_g_num;
  *dparent_entity2_g_num        = _dparent_entity2_g_num;

  *extract_entity1_distribution = malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  *extract_entity2_distribution = malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t *_extract_entity1_distribution = *extract_entity1_distribution;
  PDM_g_num_t *_extract_entity2_distribution = *extract_entity2_distribution;

  for(int i = 0; i < n_rank+1; ++i) {
    _extract_entity1_distribution[i] = distrib_extract_entity1[i];
    _extract_entity2_distribution[i] = distrib_extract_entity2[i];
  }

  PDM_part_to_block_free(ptb);
  PDM_part_to_block_free(ptb_entity2);
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
  PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx );
  /*
   *  Choice of extraction
   */

  PDM_g_num_t   *extract_face_distribution = NULL;
  PDM_g_num_t   *extract_vtx_distribution  = NULL;
  int           *dextract_face_vtx_idx     = NULL;
  PDM_g_num_t   *dextract_face_vtx         = NULL;
  PDM_g_num_t   *dparent_face_g_num        = NULL;
  PDM_g_num_t   *dparent_vtx_g_num         = NULL;

  dconnectivity_to_extract_dconnectivity(comm,
                                         dface_group_idx[n_face_group],
                                         dface_group,
                                         face_distribution,
                                         dface_vtx_idx,
                                         dface_vtx,
                                         &extract_face_distribution,
                                         &extract_vtx_distribution,
                                         &dextract_face_vtx_idx,
                                         &dextract_face_vtx,
                                         &dparent_face_g_num,
                                         &dparent_vtx_g_num);


  int dn_extract_face = extract_face_distribution[i_rank+1] - extract_face_distribution[i_rank];
  int dn_extract_vtx  = extract_vtx_distribution [i_rank+1] - extract_vtx_distribution [i_rank];
  if(1 == 1) {
    PDM_log_trace_array_long(extract_face_distribution, n_rank+1, "extract_face_distribution:: ");
    PDM_log_trace_array_long(extract_vtx_distribution , n_rank+1, "extract_vtx_distribution::  ");

    PDM_log_trace_array_int(dextract_face_vtx_idx, dn_extract_face+1                     , "dextract_face_vtx_idx:: ");
    PDM_log_trace_array_long(dextract_face_vtx   , dextract_face_vtx_idx[dn_extract_face], "dextract_face_vtx:: "    );
    PDM_log_trace_array_long(dparent_face_g_num  , dn_extract_face                       , "dparent_face_g_num:: "   );
    PDM_log_trace_array_long(dparent_vtx_g_num   , dn_extract_vtx                        , "dparent_vtx_g_num:: "    );

  }

  double** tmp_dextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        vtx_distribution,
                                        dvtx_coord,
                                        &dn_extract_vtx,
                 (const PDM_g_num_t **) &dparent_vtx_g_num,
                                        &tmp_dextract_vtx_coord);

  double* dextract_vtx_coord = tmp_dextract_vtx_coord[0];
  free(tmp_dextract_vtx_coord);

  /*
   *  Visulisation
   */
  PDM_g_num_t* extract_face_ln_to_gn = malloc(dn_extract_face * sizeof(PDM_g_num_t));
  // PDM_g_num_t* extract_vtx_ln_to_gn  = malloc(dn_extract_vtx  * sizeof(PDM_g_num_t));

  for(int i = 0; i < dn_extract_face; ++i) {
    extract_face_ln_to_gn[i] = extract_face_distribution[i_rank] + i + 1;
  }

  // for(int i = 0; i < dn_extract_vtx; ++i) {
  //   extract_vtx_ln_to_gn[i] = extract_vtx_distribution[i_rank] + i + 1;
  // }

  int pn_extract_vtx = -1;
  PDM_g_num_t *pextract_vtx_ln_to_gn = NULL;
  int         *pextract_face_vtx_idx = NULL;
  int         *pextract_face_vtx     = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           extract_face_distribution,
                                                           dextract_face_vtx_idx,
                                                           dextract_face_vtx,
                                                           dn_extract_face,
                                                           extract_face_ln_to_gn,
                                                           &pn_extract_vtx,
                                                           &pextract_vtx_ln_to_gn,
                                                           &pextract_face_vtx_idx,
                                                           &pextract_face_vtx);

  double** tmp_pextract_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        extract_vtx_distribution,
                                        dextract_vtx_coord,
                                        &pn_extract_vtx,
                 (const PDM_g_num_t **) &pextract_vtx_ln_to_gn,
                                        &tmp_pextract_vtx_coord);
  double* pextract_vtx_coord = tmp_pextract_vtx_coord[0];
  free(tmp_pextract_vtx_coord);

  char filename[999];
  sprintf(filename, "export_face_%i.vtk", i_rank);
  PDM_vtk_write_polydata(filename,
                         pn_extract_vtx,
                         pextract_vtx_coord,
                         pextract_vtx_ln_to_gn,
                         dn_extract_face,
                         pextract_face_vtx_idx,
                         pextract_face_vtx,
                         extract_face_ln_to_gn,
                         NULL);

  free(dextract_vtx_coord);
  free(extract_face_ln_to_gn);
  // free(extract_vtx_ln_to_gn );

  free(extract_face_distribution);
  free(extract_vtx_distribution );
  free(dextract_face_vtx_idx    );
  free(dextract_face_vtx        );
  free(dparent_face_g_num       );
  free(dparent_vtx_g_num        );

  free(pextract_vtx_ln_to_gn);
  free(pextract_face_vtx_idx);
  free(pextract_face_vtx    );
  free(pextract_vtx_coord   );

  free(face_distribution);
  free(vtx_distribution );

  PDM_dcube_gen_free(dcube);

  PDM_MPI_Finalize();

  return 0;
}
