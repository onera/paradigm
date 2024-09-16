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
#include "pdm_writer_priv.h"
#include "pdm_generate_mesh.h"
#include "pdm_extract_part.h"

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
  printf
    ("\n"
     "  -h                This message.\n\n"
     "  -n_part_in  <n>   Number of partitions in initial mesh.\n\n"
     "  -n_part_out <n>   Number of partitions in extracted mesh.\n\n"
     "  -visu             Enable exports for visualization.\n\n"
     "  -n          <n>   Number of vertices per side.\n\n"
     );

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 */

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_part_in,
  int                   *n_part_out,
  int                   *visu,
  PDM_g_num_t           *n_vtx_seg
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n_part_in") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part_in = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-n_part_out") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part_out = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg = (PDM_g_num_t) atoi(argv[i]);
      }
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

int main
(
  int   argc,
  char *argv[]
)
{
  /*
   *  Init MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Read args
   */
  int         n_part_in  = 1;
  int         n_part_out = 1;
  int         visu       = 0;
  PDM_g_num_t n_vtx_seg  = 10;

  _read_args(argc,
             argv,
             &n_part_in,
             &n_part_out,
             &visu,
             &n_vtx_seg);



  int          *pn_vtx                 = NULL;
  int          *pn_edge                = NULL;
  int          *pn_face                = NULL;
  int          *pn_cell                = NULL;
  double      **pvtx_coord             = NULL;
  int         **pedge_vtx              = NULL;
  int         **pface_edge_idx         = NULL;
  int         **pface_edge             = NULL;
  int         **pface_vtx              = NULL;
  int         **pcell_face_idx         = NULL;
  int         **pcell_face             = NULL;
  PDM_g_num_t **pvtx_ln_to_gn          = NULL;
  PDM_g_num_t **pedge_ln_to_gn         = NULL;
  PDM_g_num_t **pface_ln_to_gn         = NULL;
  PDM_g_num_t **pcell_ln_to_gn         = NULL;
  int          *pn_surface             = NULL;
  int         **psurface_face_idx      = NULL;
  int         **psurface_face          = NULL;
  PDM_g_num_t **psurface_face_ln_to_gn = NULL;
  int          *pn_ridge               = NULL;
  int         **pridge_edge_idx        = NULL;
  int         **pridge_edge            = NULL;
  PDM_g_num_t **pridge_edge_ln_to_gn   = NULL;
  PDM_generate_mesh_parallelepiped_ngon(comm,
                                        PDM_MESH_NODAL_HEXA8,
                                        1,
                                        NULL,
                                        0.,
                                        0.,
                                        0.,
                                        1.,
                                        1.,
                                        1.,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_part_in,
                                        PDM_SPLIT_DUAL_WITH_HILBERT,
                                        &pn_vtx,
                                        &pn_edge,
                                        &pn_face,
                                        &pn_cell,
                                        &pvtx_coord,
                                        &pedge_vtx,
                                        &pface_edge_idx,
                                        &pface_edge,
                                        &pface_vtx,
                                        &pcell_face_idx,
                                        &pcell_face,
                                        &pvtx_ln_to_gn,
                                        &pedge_ln_to_gn,
                                        &pface_ln_to_gn,
                                        &pcell_ln_to_gn,
                                        &pn_surface,
                                        &psurface_face_idx,
                                        &psurface_face,
                                        &psurface_face_ln_to_gn,
                                        &pn_ridge,
                                        &pridge_edge_idx,
                                        &pridge_edge,
                                        &pridge_edge_ln_to_gn);


  int  *n_extract    = NULL;
  int **extract_lnum = NULL;

  PDM_malloc(n_extract,    n_part_in, int  );
  PDM_malloc(extract_lnum, n_part_in, int *);
  for (int i_part = 0; i_part < n_part_in; i_part++) {
    n_extract[i_part] = 0;
    PDM_malloc(extract_lnum[i_part], pn_cell[i_part], int);

    for (int i_cell = 0; i_cell < pn_cell[i_part]; i_cell++) {

      double _min =  HUGE_VAL;
      double _max = -HUGE_VAL;

      for (int idx_face = pcell_face_idx[i_part][i_cell]; idx_face < pcell_face_idx[i_part][i_cell+1]; idx_face++) {
        int i_face = PDM_ABS(pcell_face[i_part][idx_face]) - 1;

        for (int idx_vtx = pface_edge_idx[i_part][i_face]; idx_vtx < pface_edge_idx[i_part][i_face+1]; idx_vtx++) {
          int i_vtx = pface_vtx[i_part][idx_vtx] - 1;

          double x = pvtx_coord[i_part][3*i_vtx  ];
          double y = pvtx_coord[i_part][3*i_vtx+1];
          double z = pvtx_coord[i_part][3*i_vtx+2];

          double f = x + y + z;

          _min = PDM_MIN(_min, f);
          _max = PDM_MAX(_max, f);
        }
      }

      if (_max > 0.9 && _min < 1.1) {
        extract_lnum[i_part][n_extract[i_part]++] = i_cell + 1;
      }

    }

  }





  PDM_extract_part_t *extrp = PDM_extract_part_create(3,
                                                      n_part_in,
                                                      n_part_out,
                                                      PDM_EXTRACT_PART_KIND_REEQUILIBRATE,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_TRUE,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);

  for (int i_part = 0; i_part < n_part_in; i_part++) {
    PDM_extract_part_part_set(extrp,
                              i_part,
                              pn_cell       [i_part],
                              pn_face       [i_part],
                              0,//pn_edge       [i_part],
                              pn_vtx        [i_part],
                              pcell_face_idx[i_part],
                              pcell_face    [i_part],
                              NULL,//pface_edge_idx[i_part],
                              NULL,//pface_edge    [i_part],
                              NULL,//pedge_vtx     [i_part],
                              pface_edge_idx[i_part],
                              pface_vtx     [i_part],
                              pcell_ln_to_gn[i_part],
                              pface_ln_to_gn[i_part],
                              NULL,//pedge_ln_to_gn[i_part],
                              pvtx_ln_to_gn [i_part],
                              pvtx_coord    [i_part]);

    PDM_extract_part_selected_lnum_set(extrp,
                                       i_part,
                                       n_extract   [i_part],
                                       extract_lnum[i_part]);
  }


  PDM_extract_part_compute(extrp);

  for (int i_part = 0; i_part < n_part_in; i_part++) {
    PDM_free(extract_lnum[i_part]);
  }
  PDM_free(n_extract   );
  PDM_free(extract_lnum);



  if (visu) {
    int          *extract_n_cell        = NULL;
    int          *extract_n_face        = NULL;
    int          *extract_n_vtx         = NULL;
    int         **extract_cell_face_idx = NULL;
    int         **extract_cell_face     = NULL;
    int         **extract_face_vtx_idx  = NULL;
    int         **extract_face_vtx      = NULL;
    double      **extract_vtx_coord     = NULL;
    PDM_g_num_t **extract_cell_ln_to_gn = NULL;
    PDM_g_num_t **extract_vtx_ln_to_gn  = NULL;
    PDM_malloc(extract_n_cell       , n_part_out, int          );
    PDM_malloc(extract_n_face       , n_part_out, int          );
    PDM_malloc(extract_n_vtx        , n_part_out, int          );
    PDM_malloc(extract_cell_face_idx, n_part_out, int         *);
    PDM_malloc(extract_cell_face    , n_part_out, int         *);
    PDM_malloc(extract_face_vtx_idx , n_part_out, int         *);
    PDM_malloc(extract_face_vtx     , n_part_out, int         *);
    PDM_malloc(extract_vtx_coord    , n_part_out, double      *);
    PDM_malloc(extract_cell_ln_to_gn, n_part_out, PDM_g_num_t *);
    PDM_malloc(extract_vtx_ln_to_gn , n_part_out, PDM_g_num_t *);

    for (int i_part = 0; i_part < n_part_out; i_part++) {
      extract_n_cell[i_part] = PDM_extract_part_connectivity_get(extrp,
                                                                 i_part,
                                                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                                 &extract_cell_face    [i_part],
                                                                 &extract_cell_face_idx[i_part],
                                                                 PDM_OWNERSHIP_KEEP);

      extract_n_face[i_part] = PDM_extract_part_connectivity_get(extrp,
                                                                 i_part,
                                                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                                 &extract_face_vtx    [i_part],
                                                                 &extract_face_vtx_idx[i_part],
                                                                 PDM_OWNERSHIP_KEEP);

      extract_n_vtx[i_part] = PDM_extract_part_vtx_coord_get(extrp,
                                                             i_part,
                                                             &extract_vtx_coord[i_part],
                                                             PDM_OWNERSHIP_KEEP);

      PDM_extract_part_ln_to_gn_get(extrp,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &extract_cell_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

      PDM_extract_part_ln_to_gn_get(extrp,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    &extract_vtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);
    }

    writer_wrapper(comm,
                   "extract_part_reequilibrate",
                   "extracted_mesh",
                   n_part_out,
                   extract_n_vtx,
                   extract_vtx_coord,
                   extract_vtx_ln_to_gn,
                   extract_n_cell,
                   extract_face_vtx_idx,
                   extract_face_vtx,
                   extract_cell_ln_to_gn,
                   -1,
                   extract_n_face,
                   extract_cell_face_idx,
                   extract_cell_face,
                   "Ensight",
                   0,
                   NULL,
                   NULL,
                   0,
                   NULL,
                   NULL);


    PDM_free(extract_n_cell       );
    PDM_free(extract_n_face       );
    PDM_free(extract_n_vtx        );
    PDM_free(extract_cell_face_idx);
    PDM_free(extract_cell_face    );
    PDM_free(extract_face_vtx_idx );
    PDM_free(extract_face_vtx     );
    PDM_free(extract_vtx_coord    );
    PDM_free(extract_cell_ln_to_gn);
    PDM_free(extract_vtx_ln_to_gn );
  }


  /* Free memory */
  PDM_extract_part_free(extrp);

  for (int i_part = 0; i_part < n_part_in; i_part++) {
    PDM_free(pvtx_coord             [i_part]);
    PDM_free(pedge_vtx              [i_part]);
    PDM_free(pface_edge_idx         [i_part]);
    PDM_free(pface_edge             [i_part]);
    PDM_free(pface_vtx              [i_part]);
    PDM_free(pcell_face_idx         [i_part]);
    PDM_free(pcell_face             [i_part]);
    PDM_free(pvtx_ln_to_gn          [i_part]);
    PDM_free(pedge_ln_to_gn         [i_part]);
    PDM_free(pface_ln_to_gn         [i_part]);
    PDM_free(pcell_ln_to_gn         [i_part]);
    PDM_free(psurface_face_idx      [i_part]);
    PDM_free(psurface_face          [i_part]);
    PDM_free(psurface_face_ln_to_gn [i_part]);
    PDM_free(pridge_edge_idx        [i_part]);
    PDM_free(pridge_edge            [i_part]);
    PDM_free(pridge_edge_ln_to_gn   [i_part]);
  }
  PDM_free(pn_vtx                );
  PDM_free(pn_edge               );
  PDM_free(pn_face               );
  PDM_free(pn_cell               );
  PDM_free(pvtx_coord            );
  PDM_free(pedge_vtx             );
  PDM_free(pface_edge_idx        );
  PDM_free(pface_edge            );
  PDM_free(pface_vtx             );
  PDM_free(pcell_face_idx        );
  PDM_free(pcell_face            );
  PDM_free(pvtx_ln_to_gn         );
  PDM_free(pedge_ln_to_gn        );
  PDM_free(pface_ln_to_gn        );
  PDM_free(pcell_ln_to_gn        );
  PDM_free(pn_surface            );
  PDM_free(psurface_face_idx     );
  PDM_free(psurface_face         );
  PDM_free(psurface_face_ln_to_gn);
  PDM_free(pn_ridge              );
  PDM_free(pridge_edge_idx       );
  PDM_free(pridge_edge           );
  PDM_free(pridge_edge_ln_to_gn  );

  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
