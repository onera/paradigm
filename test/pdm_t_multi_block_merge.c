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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multi_block_merge.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


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
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
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

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 4;
  double             length  = 1.;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);
  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_dcube_nodal_t* dcube1 = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube1, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube1);

  PDM_dmesh_nodal_t*  dmn1 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube1);
  /*
   * Define distribution of cell
   */
  PDM_dmesh_nodal_dump_vtk(dmn1, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_dcube1_");

  /*
   * Define distibution of vtx
   */
  PDM_dcube_nodal_t* dcube2 = PDM_dcube_nodal_gen_create(comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         1.,
                                                         0.,
                                                         0.,
                                                         PDM_MESH_NODAL_QUAD4,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube2, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube2);

  PDM_dmesh_nodal_t*  dmn2 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube2);

  PDM_dmesh_nodal_dump_vtk(dmn2, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_dcube2_");

  /*
   * Concatenate blocks
   */
  int n_block = 2;
  PDM_g_num_t** block_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_distrib_idx[0] = PDM_dmesh_nodal_vtx_distrib_get(dmn1);
  block_distrib_idx[1] = PDM_dmesh_nodal_vtx_distrib_get(dmn2);

  int* n_selected = malloc(n_block * sizeof(int));
  n_selected[0] = PDM_DMesh_nodal_n_vtx_get(dmn1);
  n_selected[1] = PDM_DMesh_nodal_n_vtx_get(dmn2);

  PDM_g_num_t** selected_g_num = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    selected_g_num[i_block] = malloc(n_selected[i_block] * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_selected[i_block]; ++i) {
      selected_g_num[i_block][i] = block_distrib_idx[i_block][i_rank] + i + 1;
    }
  }

  /*
   * Setup graph
   */
  int         **dmerge_idx      = malloc(n_block * sizeof(PDM_g_num_t *));
  int         **dmerge_block_id = malloc(n_block * sizeof(int         *));
  PDM_g_num_t **dmerge_g_num    = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    int dn_vtx = block_distrib_idx[i_block][i_rank+1] - block_distrib_idx[i_block][i_rank];

    dmerge_idx     [i_block] = malloc( ( dn_vtx + 1 ) * sizeof(int        ));
    dmerge_block_id[i_block] = malloc( ( dn_vtx     ) * sizeof(int        ));
    dmerge_g_num   [i_block] = malloc( ( dn_vtx     ) * sizeof(PDM_g_num_t));

    dmerge_idx     [i_block][0] = 0;

    PDM_g_num_t vtx_g_num_next = 1;
    if(i_block == 1) {
      vtx_g_num_next = 5;
    }

    int idx_write = 0;
    for(int j = 0; j < dn_vtx; ++j) {
      PDM_g_num_t vtx_g_num = block_distrib_idx[i_block][i_rank] + j + 1;
      PDM_g_num_t indi      = vtx_g_num % ( n_vtx_seg + 1 );

      dmerge_idx[i_block][j+1] = dmerge_idx[i_block][j];

      if(indi == 0 && i_block == 0) {
        dmerge_idx     [i_block][j+1] = dmerge_idx[i_block][j] + 1;
        dmerge_block_id[i_block][idx_write] = 1;
        dmerge_g_num   [i_block][idx_write] = vtx_g_num_next;
        vtx_g_num_next += n_vtx_seg+1;

        idx_write++;
      } else if(indi == 1 && i_block == 1) {
        dmerge_idx     [i_block][j+1] = dmerge_idx[i_block][j] + 1;
        dmerge_block_id[i_block][idx_write] = 0;
        dmerge_g_num   [i_block][idx_write] = vtx_g_num_next;
        vtx_g_num_next += n_vtx_seg+1;
        idx_write++;
      }

      // PDM_g_num_t g_num = n_vtx_seg * i;
      // printf("vtx_g_num = %i | indi = %i  \n", vtx_g_num, indi);
    }

    dmerge_block_id[i_block] = realloc(dmerge_block_id[i_block], idx_write * sizeof(int        ));
    dmerge_g_num   [i_block] = realloc(dmerge_g_num   [i_block], idx_write * sizeof(PDM_g_num_t));

    PDM_log_trace_array_int (dmerge_block_id[i_block], idx_write, "dmerge_block_id :: ");
    PDM_log_trace_array_long(dmerge_g_num   [i_block], idx_write, "dmerge_g_num    :: ");


  }



  // PDM_multi_block_merge_t* mbm = PDM_multi_block_merge_create(comm);

  for(int i_block = 0; i_block < n_block ; ++i_block) {
    free(selected_g_num[i_block]);
    free(dmerge_block_id[i_block]);
    free(dmerge_g_num[i_block]);
    free(dmerge_idx[i_block]);
  }
  free(selected_g_num);
  free(dmerge_idx     );
  free(dmerge_block_id);
  free(dmerge_g_num   );

  free(block_distrib_idx);
  free(n_selected);


  PDM_dcube_nodal_gen_free(dcube1);
  PDM_dcube_nodal_gen_free(dcube2);

  PDM_MPI_Finalize ();
  return 0;
}
