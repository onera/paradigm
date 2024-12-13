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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_generate_mesh.h"
#include "pdm_triangulate.h"


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
     "  -n     <n> Number of vertices per side of the square mesh.\n\n"
     "  -visu      Set meshes vtk output.\n\n"
     "  -h         This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc       Number of arguments
 * \param [in]    argv       Arguments
 * \param [inout] n_vtx_seg  Number of vertices per side of the square mesh
 * \param [inout] visu       Set meshes vtk output
 *
 */

static void
_read_args
(
  int            argc,
  char         **argv,
  PDM_g_num_t   *n_vtx_seg,
  int           *visu
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }

    i++;
  }
}


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
  // Set default values
  PDM_g_num_t n_vtx_seg = 5;
  int         visu      = 0;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &visu);

  // Initialize MPI
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  // Generate polygonal mesh
  int          *pn_vtx         = NULL;
  int          *pn_edge        = NULL;
  int          *pn_face        = NULL;
  double      **pvtx_coord     = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  int         **pface_vtx      = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_POLY_2D,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   n_vtx_seg,
                                   n_vtx_seg,
                                   1,
                                   PDM_SPLIT_DUAL_WITH_HILBERT,
                                   1.,
                                   &pn_vtx,
                                   &pn_edge,
                                   &pn_face,
                                   &pvtx_coord,
                                   &pedge_vtx,
                                   &pface_edge_idx,
                                   &pface_edge,
                                   &pface_vtx,
                                   &pvtx_ln_to_gn,
                                   &pedge_ln_to_gn,
                                   &pface_ln_to_gn);

  // Triangulate faces
  int *face_tria_idx = NULL;
  int *tria_vtx      = NULL;
  int n_tria = PDM_triangulate_faces(pn_face       [0],
                                     pface_edge_idx[0],
                                     pface_vtx     [0],
                                     pvtx_coord    [0],
                                     &face_tria_idx,
                                     &tria_vtx,
                                     NULL);

  if (visu) {
    PDM_g_num_t *tria_parent = NULL;
    PDM_malloc(tria_parent, n_tria, PDM_g_num_t);
    for (int i_face = 0; i_face < pn_face[0]; i_face++) {
      for (int i_tria = face_tria_idx[i_face]; i_tria < face_tria_idx[i_face+1]; i_tria++) {
        tria_parent[i_tria] = pface_ln_to_gn[0][i_face];
      }
    }

    char filename[999];
    sprintf(filename, "triangulate_faces_%d.vtk", i_rank);
    PDM_vtk_write_std_elements(filename,
                               pn_vtx       [0],
                               pvtx_coord   [0],
                               pvtx_ln_to_gn[0],
                               PDM_MESH_NODAL_TRIA3,
                               n_tria,
                               tria_vtx,
                               tria_parent,
                               0,
                               NULL,
                               NULL);
    free(tria_parent);
  }


  // Free memory
  free(face_tria_idx);
  free(tria_vtx     );

  free(pvtx_coord    [0]);
  free(pedge_vtx     [0]);
  free(pface_edge_idx[0]);
  free(pface_edge    [0]);
  free(pface_vtx     [0]);
  free(pvtx_ln_to_gn [0]);
  free(pedge_ln_to_gn[0]);
  free(pface_ln_to_gn[0]);
  free(pn_vtx        );
  free(pn_edge       );
  free(pn_face       );
  free(pvtx_coord    );
  free(pedge_vtx     );
  free(pface_edge_idx);
  free(pface_edge    );
  free(pface_vtx     );
  free(pvtx_ln_to_gn );
  free(pedge_ln_to_gn);
  free(pface_ln_to_gn);

  // Finalize MPI
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}