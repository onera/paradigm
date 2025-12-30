#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_mem_tool.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_multipart.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -n_part <level>  Number of partition                        .\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -t               Element kind .\n\n"
     "  -h               This message.\n\n");
  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_a,
 int                   *post
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_a = atol(argv[i]);
        *n_vtx_a = (PDM_g_num_t) _n_vtx_a;
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
PDM_part_mesh_nodal_t*
_generate_mesh
(
  PDM_MPI_Comm    pdm_comm,
  int             n_vtx_seg
)
{
  int              n_part       = 1;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;

  /* Warmup */
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        2.,//lenght
                                                        -1,//-0.5*length,//-0.1,
                                                        -1,//-0.5*length,//-0.2,
                                                        -1,//-0.5*length,//0.07,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                split_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t* pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmesh_nodal, PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_complete_part_comm_graph(pmesh_nodal);

  PDM_multipart_free(mpart);
  PDM_dcube_nodal_gen_free(dcube);
  return pmesh_nodal;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

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
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_vtx_a = 10;

  int post = 0;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &post);


  PDM_part_mesh_nodal_t* pmn = _generate_mesh(comm, n_vtx_a);

  PDM_part_mesh_nodal_free(pmn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
