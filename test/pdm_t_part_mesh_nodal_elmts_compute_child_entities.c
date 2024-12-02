#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"

#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_writer_priv.h"
#include "pdm_vtk.h"
#include "pdm_array.h"

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
     "  -n                    <n>  Number of vertices on the cube side.\n\n"
     "  -n_part               <n>  Number of partitions per MPI rank.\n\n"
     "  -elt_type             <t>  Element type.\n\n"
     "  -compute_parent_child      Compute only ascending link.\n\n"
     "  -no_child_pmne             Do not use child pmne.\n\n"
     "  -v                         Enable verbose mode.\n\n"
     "  -parmetis                  Call ParMETIS.\n\n"
     "  -pt-scotch                 Call PT-Scotch.\n\n"
     "  -hilbert                   Call PT-Hilbert.\n\n"
     "  -h                         This message.\n\n");

  exit(exit_code);
}


/*
 * Read arguments from the command line
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_seg,
 int                   *n_part,
 PDM_Mesh_nodal_elt_t  *elt_type,
 PDM_split_dual_t      *part_method,
 PDM_bool_t            *compute_parent_child,
 PDM_bool_t            *no_child_pmne,
 int                   *verbose
 )
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
        long n = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) n;
      }
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-elt_type") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-compute_parent_child") == 0) {
      *compute_parent_child = PDM_TRUE;
    }
    else if (strcmp(argv[i], "-no_child_pmne") == 0) {
      *no_child_pmne = PDM_TRUE;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
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

  PDM_g_num_t          n_vtx_seg            = 10;
  int                  n_part               = 1;
  PDM_Mesh_nodal_elt_t elt_type             = PDM_MESH_NODAL_HEXA8;
  PDM_split_dual_t     split_method         = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_bool_t           compute_parent_child = PDM_FALSE;
  int                  verbose              = 0;
  PDM_bool_t           no_child_pmne        = 0;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_part,
             &elt_type,
             &split_method,
             &compute_parent_child,
             &no_child_pmne,
             &verbose);


  /* Initialize MPI */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);


  /* Generate block-distributed nodal mesh */
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        elt_type,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_build(dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);


  /* Split mesh */
  int n_domain = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part,
                                                PDM_FALSE,
                                                split_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  PDM_multipart_compute(mpart);


  /* Let's go */
  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(elt_type);

  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmn, PDM_OWNERSHIP_KEEP);

  int *pn_vtx = NULL;
  PDM_malloc(pn_vtx, n_part, int);
  for (int i_part = 0; i_part < n_part; i_part++) {
    pn_vtx[i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
  }

  PDM_part_mesh_nodal_elmts_t *pmne_ridge = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, PDM_GEOMETRY_KIND_RIDGE   );
  PDM_part_mesh_nodal_elmts_t *pmne_surf  = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC);

  if (no_child_pmne) {
    pmne_ridge = NULL;
  }

  /* Compute link Surface <-> Ridge */
  int **ridge_edge_surf_face_idx = NULL;
  int **ridge_edge_surf_face     = NULL;
  int  *n_surf_edge              = NULL;
  int **surf_edge_vtx_idx        = NULL;
  int **surf_edge_vtx            = NULL;
  int **surf_face_edge_idx       = NULL;
  int **surf_face_edge           = NULL;
  log_trace("surf ridge\n");
  PDM_part_mesh_nodal_elmts_compute_child_parent(pmne_surf,
                                                 pmne_ridge,
                                                 PDM_MESH_ENTITY_EDGE,
                                                 compute_parent_child,
                                                 pn_vtx,
                                                 &ridge_edge_surf_face_idx,
                                                 &ridge_edge_surf_face,
                                                 &n_surf_edge,
                                                 &surf_edge_vtx_idx,
                                                 &surf_edge_vtx,
                                                 &surf_face_edge_idx,
                                                 &surf_face_edge);

  if (verbose) {
    log_trace("\n=== Surface <-> Ridge ====\n");
    for (int i_part = 0; i_part < n_part; i_part++) {
      log_trace("-- i_part = %d\n", i_part);

      if (!no_child_pmne) {

        // PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE, "pmn_ridge");

        int n_ridge_edge = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_ridge, i_part);
        log_trace("n_ridge_edge = %d\n", n_ridge_edge);
        PDM_log_trace_connectivity_int(ridge_edge_surf_face_idx[i_part],
                                       ridge_edge_surf_face    [i_part],
                                       n_ridge_edge,
                                       "ridge_edge_surf_face : ");

        int *ridge_edge_vtx_idx;
        int *ridge_edge_vtx;
        PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne_ridge,
                                                       i_part,
                                                       &ridge_edge_vtx_idx,
                                                       &ridge_edge_vtx);
        PDM_log_trace_connectivity_int(ridge_edge_vtx_idx,
                                       ridge_edge_vtx,
                                       n_ridge_edge,
                                       "ridge_edge_vtx : ");
        free(ridge_edge_vtx_idx);
        free(ridge_edge_vtx    );

        if (compute_parent_child) {
          log_trace("n_surf_edge  = %d\n", n_surf_edge[i_part]);
          PDM_log_trace_connectivity_int(surf_edge_vtx_idx[i_part],
                                         surf_edge_vtx    [i_part],
                                         n_surf_edge      [i_part],
                                         "surf_edge_vtx : ");
        }

      }


      if (compute_parent_child) {
        int n_surf_face = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_surf, i_part);
        PDM_log_trace_connectivity_int(surf_face_edge_idx[i_part],
                                       surf_face_edge    [i_part],
                                       n_surf_face,
                                       "surf_face_edge : ");
      }
    }
  }

  int **surf_face_cell_idx = NULL;
  int **surf_face_cell     = NULL;
  int  *n_vol_face         = NULL;
  int **vol_face_vtx_idx   = NULL;
  int **vol_face_vtx       = NULL;
  int **cell_vol_face_idx  = NULL;
  int **cell_vol_face      = NULL;

  int **ridge_edge_cell_idx = NULL;
  int **ridge_edge_cell     = NULL;
  int  *n_vol_edge          = NULL;
  int **vol_edge_vtx_idx    = NULL;
  int **vol_edge_vtx        = NULL;
  int **cell_vol_edge_idx   = NULL;
  int **cell_vol_edge       = NULL;


  if (mesh_dimension == 3) {
    PDM_part_mesh_nodal_elmts_t *pmne_vol = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, PDM_GEOMETRY_KIND_VOLUMIC);

    if (no_child_pmne) {
      pmne_surf = NULL;
    }

    /* Compute link Volume <-> Surface */
    log_trace("vol surf\n");
    PDM_part_mesh_nodal_elmts_compute_child_parent(pmne_vol,
                                                   pmne_surf,
                                                   PDM_MESH_ENTITY_FACE,
                                                   compute_parent_child,
                                                   pn_vtx,
                                                   &surf_face_cell_idx,
                                                   &surf_face_cell,
                                                   &n_vol_face,
                                                   &vol_face_vtx_idx,
                                                   &vol_face_vtx,
                                                   &cell_vol_face_idx,
                                                   &cell_vol_face);

    if (verbose) {
      log_trace("\n=== Volume <-> Surface ====\n");
      for (int i_part = 0; i_part < n_part; i_part++) {
        log_trace("-- i_part = %d\n", i_part);

        if (!no_child_pmne) {

          // PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "pmn_surf");

          int n_surf_face = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_surf, i_part);
          log_trace("n_surf_face = %d\n", n_surf_face);
          PDM_log_trace_connectivity_int(surf_face_cell_idx[i_part],
                                         surf_face_cell    [i_part],
                                         n_surf_face,
                                         "surf_face_cell : ");

          int *surf_face_vtx_idx;
          int *surf_face_vtx;
          PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne_surf,
                                                         i_part,
                                                         &surf_face_vtx_idx,
                                                         &surf_face_vtx);
          PDM_log_trace_connectivity_int(surf_face_vtx_idx,
                                         surf_face_vtx,
                                         n_surf_face,
                                         "surf_face_vtx : ");
          free(surf_face_vtx_idx);
          free(surf_face_vtx    );

          if (compute_parent_child) {
            log_trace("n_vol_face  = %d\n", n_vol_face[i_part]);
            PDM_log_trace_connectivity_int(vol_face_vtx_idx[i_part],
                                           vol_face_vtx    [i_part],
                                           n_vol_face      [i_part],
                                           "vol_face_vtx : ");
          }

        }


        if (compute_parent_child) {
          int n_cell = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_vol, i_part);
          PDM_log_trace_connectivity_int(cell_vol_face_idx[i_part],
                                         cell_vol_face    [i_part],
                                         n_cell,
                                         "cell_vol_face : ");
        }
      }
    }

    /* Compute link Volume <-> Ridge */
    log_trace("vol ridge\n");
    PDM_part_mesh_nodal_elmts_compute_child_parent(pmne_vol,
                                                   pmne_ridge,
                                                   PDM_MESH_ENTITY_EDGE,
                                                   compute_parent_child,
                                                   pn_vtx,
                                                   &ridge_edge_cell_idx,
                                                   &ridge_edge_cell,
                                                   &n_vol_edge,
                                                   &vol_edge_vtx_idx,
                                                   &vol_edge_vtx,
                                                   &cell_vol_edge_idx,
                                                   &cell_vol_edge);
    if (verbose) {
      log_trace("\n=== Volume <-> Ridge ====\n");
      for (int i_part = 0; i_part < n_part; i_part++) {
        log_trace("-- i_part = %d\n", i_part);

        if (!no_child_pmne) {

          int n_ridge_edge = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_ridge, i_part);
          log_trace("n_ridge_edge = %d\n", n_ridge_edge);
          PDM_log_trace_connectivity_int(ridge_edge_cell_idx[i_part],
                                         ridge_edge_cell    [i_part],
                                         n_ridge_edge,
                                         "ridge_edge_cell : ");

          int *ridge_edge_vtx_idx;
          int *ridge_edge_vtx;
          PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne_ridge,
                                                         i_part,
                                                         &ridge_edge_vtx_idx,
                                                         &ridge_edge_vtx);
          PDM_log_trace_connectivity_int(ridge_edge_vtx_idx,
                                         ridge_edge_vtx,
                                         n_ridge_edge,
                                         "ridge_edge_vtx : ");
          free(ridge_edge_vtx_idx);
          free(ridge_edge_vtx    );

          if (compute_parent_child) {
            log_trace("n_vol_edge  = %d\n", n_vol_edge[i_part]);
            PDM_log_trace_connectivity_int(vol_edge_vtx_idx[i_part],
                                           vol_edge_vtx    [i_part],
                                           n_vol_edge      [i_part],
                                           "vol_edge_vtx : ");
          }

        }


        if (compute_parent_child) {
          int n_cell = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne_vol, i_part);
          PDM_log_trace_connectivity_int(cell_vol_edge_idx[i_part],
                                         cell_vol_edge    [i_part],
                                         n_cell,
                                         "cell_edge : ");
        }
      }
    }
  }


  /*
   * Visu
   */
  if (verbose && compute_parent_child) {
    double      **pvtx_coord     = NULL;
    PDM_g_num_t **pvtx_ln_to_gn  = NULL;
    int          *pn_cell        = NULL;
    PDM_g_num_t **pcell_ln_to_gn = NULL;
    PDM_malloc(pvtx_coord    , n_part, double      *);
    PDM_malloc(pvtx_ln_to_gn , n_part, PDM_g_num_t *);
    PDM_malloc(pn_cell       , n_part, int          );
    PDM_malloc(pcell_ln_to_gn, n_part, PDM_g_num_t *);

    for (int i_part = 0; i_part < n_part; i_part++) {
      pvtx_coord    [i_part] = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      pvtx_ln_to_gn [i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part);

      if (mesh_dimension == 3) {
        PDM_part_mesh_nodal_elmts_t *pmne_vol = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, PDM_GEOMETRY_KIND_VOLUMIC);

        pn_cell       [i_part] = PDM_part_mesh_nodal_elmts_n_elmts_get        (pmne_vol, i_part);
        pcell_ln_to_gn[i_part] = PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne_vol, i_part, PDM_OWNERSHIP_KEEP);
      }


      char name[999];
      sprintf(name, "pmesh_nodal_to_pmesh_surf_edge_%d.vtk", i_rank*n_part + i_part);
      PDM_vtk_write_std_elements(name,
                                 pn_vtx       [i_part],
                                 pvtx_coord   [i_part],
                                 pvtx_ln_to_gn[i_part],
                                 PDM_MESH_NODAL_BAR2,
                                 n_surf_edge  [i_part],
                                 surf_edge_vtx[i_part],
                                 NULL,
                                 0,
                                 NULL,
                                 NULL);
    }

    if (mesh_dimension == 3) {
      PDM_writer_elt_geom_t cell_type = (PDM_writer_elt_geom_t) -1;
      writer_wrapper(comm,
                     "pmesh_nodal_to_pmesh_3d",
                     "3d",
                     n_part,
                     pn_vtx,
                     pvtx_coord,
                     pvtx_ln_to_gn,
                     pn_cell,
                     vol_face_vtx_idx,
                     vol_face_vtx,
                     pcell_ln_to_gn,
                     cell_type,
                     n_vol_face,
                     cell_vol_face_idx,
                     cell_vol_face,
                     "Ensight",
                     0,
                     NULL,
                     NULL,
                     0,
                     NULL,
                     NULL);

      for (int i_part = 0; i_part < n_part; i_part++) {
        char name[999];
        sprintf(name, "pmesh_nodal_to_pmesh_vol_edge_%d.vtk", i_rank*n_part + i_part);
        PDM_vtk_write_std_elements(name,
                                   pn_vtx       [i_part],
                                   pvtx_coord   [i_part],
                                   pvtx_ln_to_gn[i_part],
                                   PDM_MESH_NODAL_BAR2,
                                   n_vol_edge   [i_part],
                                   vol_edge_vtx [i_part],
                                   NULL,
                                   0,
                                   NULL,
                                   NULL);
      }
    }
    PDM_free(pvtx_coord   );
    PDM_free(pvtx_ln_to_gn);
    PDM_free(pn_cell       );
    PDM_free(pcell_ln_to_gn);
  }


  /*
   * Free memory
   */
  PDM_free(pn_vtx);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(ridge_edge_surf_face_idx[i_part]);
    PDM_free(ridge_edge_surf_face    [i_part]);
    if (compute_parent_child) {
      PDM_free(surf_edge_vtx_idx [i_part]);
      PDM_free(surf_edge_vtx     [i_part]);
      PDM_free(surf_face_edge_idx[i_part]);
      PDM_free(surf_face_edge    [i_part]);
    }
    if (mesh_dimension == 3) {
      PDM_free(surf_face_cell_idx[i_part]);
      PDM_free(surf_face_cell    [i_part]);
      if (compute_parent_child) {
        PDM_free(vol_face_vtx_idx [i_part]);
        PDM_free(vol_face_vtx     [i_part]);
        PDM_free(cell_vol_face_idx[i_part]);
        PDM_free(cell_vol_face    [i_part]);
      }

      PDM_free(ridge_edge_cell_idx[i_part]);
      PDM_free(ridge_edge_cell    [i_part]);
      if (compute_parent_child) {
        PDM_free(vol_edge_vtx_idx [i_part]);
        PDM_free(vol_edge_vtx     [i_part]);
        PDM_free(cell_vol_edge_idx[i_part]);
        PDM_free(cell_vol_edge    [i_part]);
      }
    }
  }
  PDM_free(n_surf_edge             );
  PDM_free(ridge_edge_surf_face_idx);
  PDM_free(ridge_edge_surf_face    );
  PDM_free(surf_edge_vtx_idx       );
  PDM_free(surf_edge_vtx           );
  PDM_free(surf_face_edge_idx      );
  PDM_free(surf_face_edge          );

  if (mesh_dimension == 3) {
    PDM_free(n_vol_face        );
    PDM_free(surf_face_cell_idx);
    PDM_free(surf_face_cell    );
    PDM_free(vol_face_vtx_idx  );
    PDM_free(vol_face_vtx      );
    PDM_free(cell_vol_face_idx );
    PDM_free(cell_vol_face     );

    PDM_free(n_vol_edge         );
    PDM_free(ridge_edge_cell_idx);
    PDM_free(ridge_edge_cell    );
    PDM_free(vol_edge_vtx_idx   );
    PDM_free(vol_edge_vtx       );
    PDM_free(cell_vol_edge_idx  );
    PDM_free(cell_vol_edge      );
  }


  /* Free memory */
  PDM_part_mesh_nodal_free(pmn);
  PDM_dcube_nodal_gen_free(dcube);
  PDM_multipart_free(mpart);

  /* Finalize */
  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("The End :D\n");
    fflush(stdout);
  }

  return EXIT_SUCCESS;
}
