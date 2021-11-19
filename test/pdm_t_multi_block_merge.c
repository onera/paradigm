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

  PDM_g_num_t        n_vtx_seg = 10;
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

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t*  dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  /*
   * Define distribution of cell
   */
  PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_in_");

  /*
   * Define distibution of vtx
   */
  PDM_g_num_t n_vtx_seg_add = 4;
  PDM_dcube_nodal_t* dcube_add = PDM_dcube_nodal_gen_create(comm,
                                                            n_vtx_seg_add,
                                                            n_vtx_seg_add,
                                                            n_vtx_seg_add,
                                                            length,
                                                            0.,
                                                            1.,
                                                            0.,
                                                            PDM_MESH_NODAL_TRIA3,
                                                            1,
                                                            PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube_add, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube_add);

  PDM_dmesh_nodal_t*  dmn_add = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube_add);

  PDM_dmesh_nodal_dump_vtk(dmn_add, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_add_in_");


  PDM_dcube_nodal_gen_free(dcube);
  PDM_dcube_nodal_gen_free(dcube_add);

  PDM_MPI_Finalize ();
  return 0;
}
