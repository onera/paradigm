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
#include "pdm_partitioning_algorithm.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_dmesh.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_part_mesh_nodal_to_pmesh.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_distrib.h"

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
     "  -pt-scotch       Call PT-Scotch.\n\n"
     "  -t               Element type.\n\n"
     "  -old             Use old algorithm.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 */

static void
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length,
           int                   *n_part,
           int                   *post,
           PDM_Mesh_nodal_elt_t  *elt_type,
           int                   *part_method,
           int                   *old_algo,
           int                   *do_gnums)
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
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-old") == 0) {
      *old_algo = 1;
    }
    else if (strcmp(argv[i], "-do_gnums") == 0) {
      *do_gnums = 1;
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

  PDM_g_num_t        n_vtx_seg  = 10;
  double             length     = 1.;
  int                n_part     = 1;
  int                post       = 0;
  int                old_algo   = 0;
  int                do_gnums   = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
             &elt_type,
     (int *) &part_method,
             &old_algo,
             &do_gnums);

  assert(elt_type != PDM_MESH_NODAL_POLY_3D); // TODO: poly_vol_gen en dcube_nodal_gen

  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         0.,
                                                         0.,
                                                         0.,
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);

  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);

  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  // if(post) {
  //   if(dim == 3) {
  //     PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "dmn_volumic");
  //   }
  //   PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "dmn_surfacic");
  //   PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_RIDGE   , "dmn_ridge");
  // }

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_domain,
                                                   &n_part_domains,
                                                   PDM_FALSE,
                                                   part_method,
                                                   PDM_PART_SIZE_HOMOGENEOUS,
                                                   NULL,
                                                   comm,
                                                   PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE",
                                                     NULL,
                                                     "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart_id, 0, dmn);
  PDM_multipart_compute(mpart_id);

  PDM_part_mesh_nodal_t *pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart_id, 0, &pmesh_nodal, PDM_OWNERSHIP_KEEP);


  PDM_MPI_Barrier(comm);
  double t_start = PDM_MPI_Wtime();

  PDM_part_mesh_t *pm = NULL;
  if (old_algo) {
    PDM_dmesh_nodal_to_dmesh_transform_t transform_kind = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE;
    if(dim == 2) {
      transform_kind = PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE;
    }

    pm = PDM_part_mesh_nodal_to_part_mesh(pmesh_nodal,
                                          transform_kind,
                                          PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
  }
  else {

    PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmesh_nodal,
                                                                                            PDM_TRUE,
                                                                                            PDM_OWNERSHIP_USER);
    // Connectivities
    if (dim == 3) {
      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                           PDM_CONNECTIVITY_TYPE_CELL_FACE);
    }

    // PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
    //                                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE);

    // PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
    //                                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX);

    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX);


    // Global IDs
    if (do_gnums) {
      if (dim == 3) {
        PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm,
                                                       PDM_MESH_ENTITY_CELL);
      }
      PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm,
                                                     PDM_MESH_ENTITY_FACE);

      // PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm,
      //                                                PDM_MESH_ENTITY_EDGE);

      PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm,
                                                     PDM_MESH_ENTITY_VTX);
    }

    // Groups
    PDM_part_mesh_nodal_to_part_mesh_groups_enable(pmn_to_pm,
                                                   PDM_BOUND_TYPE_FACE);

    // PDM_part_mesh_nodal_to_part_mesh_groups_enable(pmn_to_pm,
    //                                                PDM_BOUND_TYPE_EDGE);


    PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

    PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                   &pm,
                                                   PDM_OWNERSHIP_USER);

    PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);
  }

  double t_end = PDM_MPI_Wtime();

  double delta_t = t_end - t_start;

  double min_delta_t, max_delta_t;
  PDM_MPI_Allreduce(&delta_t, &min_delta_t, 1, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(&delta_t, &max_delta_t, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  if (i_rank == 0) {
    printf("min/max elapsed time :  %8.2e / %8.2e\n", min_delta_t, max_delta_t)    ;
  }

  if(post) {
    // Part mesh nodal
    if(dim == 3) {
      PDM_part_mesh_nodal_dump_vtk(pmesh_nodal, PDM_GEOMETRY_KIND_VOLUMIC , "pmn_volumic");
    }
    PDM_part_mesh_nodal_dump_vtk(pmesh_nodal, PDM_GEOMETRY_KIND_SURFACIC, "pmn_surfacic");
    PDM_part_mesh_nodal_dump_vtk(pmesh_nodal, PDM_GEOMETRY_KIND_RIDGE   , "pmn_ridge");

    // Part mesh
    if(dim == 3) {

      PDM_part_mesh_dump_ensight(pm,
                                 "pmn_to_pm",
                                 "pmesh",
                                 PDM_TRUE);

      int n_face_group = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_FACE);
      for(int i_part = 0; i_part < n_part; ++i_part) {

        int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmesh_nodal, i_part);
        double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmesh_nodal, i_part);
        PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmesh_nodal, i_part);

        PDM_g_num_t *face_ln_to_gn = NULL;
        int         *face_vtx_idx  = NULL;
        int         *face_vtx      = NULL;
        int pn_face = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_FACE);
        PDM_part_mesh_entity_ln_to_gn_get(pm,
                                          i_part,
                                          PDM_MESH_ENTITY_FACE,
                                          &face_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);
        PDM_part_mesh_connectivity_get(pm,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       &face_vtx,
                                       &face_vtx_idx,
                                       PDM_OWNERSHIP_KEEP);

        int owner_face_vtx = 0;

        if (face_vtx == NULL) {
          owner_face_vtx = 1;

          int *face_edge;
          PDM_part_mesh_connectivity_get(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         &face_edge,
                                         &face_vtx_idx,
                                         PDM_OWNERSHIP_KEEP);

          int *edge_vtx_idx;
          int *edge_vtx;
          PDM_part_mesh_connectivity_get(pm,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         &edge_vtx,
                                         &edge_vtx_idx,
                                         PDM_OWNERSHIP_KEEP);

          PDM_compute_face_vtx_from_face_and_edge(pn_face,
                                                  face_vtx_idx,
                                                  face_edge,
                                                  edge_vtx,
                                                  &face_vtx);
        }

        int *face_flags;
        PDM_malloc(face_flags, pn_face, int);
        for(int i_face = 0; i_face < pn_face; ++i_face) {
          face_flags[i_face] = -1;
        }

        for(int i_group = 0; i_group < n_face_group; ++i_group ) {
          int          pn_bound        = 0;
          int         *pbound          = NULL;
          PDM_g_num_t *pbound_ln_to_gn = NULL;
          PDM_part_mesh_bound_get(pm,
                                  i_part,
                                  i_group,
                                  PDM_BOUND_TYPE_FACE,
                                  &pn_bound,
                                  &pbound,
                                  &pbound_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
          // PDM_log_trace_array_int(pbound, pn_bound, "pbound : ");
          for(int idx_face = 0; idx_face < pn_bound; ++idx_face) {
            int i_face = pbound[idx_face]-1;
            face_flags[i_face] = i_group;
          }
        }

        char filename[999];
        sprintf(filename, "out_pmesh_nodal_to_pmesh_%i_%i.vtk", i_part, i_rank);
        PDM_vtk_write_polydata(filename,
                               n_vtx,
                               vtx_coord,
                               vtx_ln_to_gn,
                               pn_face,
                               face_vtx_idx,
                               face_vtx,
                               face_ln_to_gn,
                               face_flags);

        PDM_free(face_flags);


        if (owner_face_vtx) {
          PDM_free(face_vtx);
        }
      }

    }

    else {

      int n_edge_group = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_EDGE);
      for(int i_part = 0; i_part < n_part; ++i_part) {

        int n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmesh_nodal, i_part);
        double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmesh_nodal, i_part);
        PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmesh_nodal, i_part);

        PDM_g_num_t *edge_ln_to_gn = NULL;
        int         *edge_vtx_idx  = NULL;
        int         *edge_vtx      = NULL;
        int pn_edge = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);
        PDM_part_mesh_entity_ln_to_gn_get(pm,
                                          i_part,
                                          PDM_MESH_ENTITY_EDGE,
                                          &edge_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);
        PDM_part_mesh_connectivity_get(pm,
                                       i_part,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       &edge_vtx,
                                       &edge_vtx_idx,
                                       PDM_OWNERSHIP_KEEP);
        // assert(edge_vtx_idx == NULL);
        int *edge_flags;
        PDM_malloc(edge_flags, pn_edge, int);
        for(int i_edge = 0; i_edge < pn_edge; ++i_edge) {
          edge_flags[i_edge] = -1;
        }

        for(int i_group = 0; i_group < n_edge_group; ++i_group ) {
          int          pn_bound        = 0;
          int         *pbound          = NULL;
          PDM_g_num_t *pbound_ln_to_gn = NULL;
          PDM_part_mesh_bound_get(pm,
                                  i_part,
                                  i_group,
                                  PDM_BOUND_TYPE_EDGE,
                                  &pn_bound,
                                  &pbound,
                                  &pbound_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
          // PDM_log_trace_array_int(pbound, pn_bound, "pbound : ");
          for(int idx_edge = 0; idx_edge < pn_bound; ++idx_edge) {
            int i_edge = pbound[idx_edge]-1;
            edge_flags[i_edge] = i_group;
          }
        }

        char filename[999];
        sprintf(filename, "out_pmesh_nodal_to_pmesh_%i_%i.vtk", i_part, i_rank);
        const char* field_name[] = {"edge_color", 0 };
        const int *edge_field[2] = {edge_flags};
        PDM_vtk_write_std_elements(filename,
                                   n_vtx,
                                   vtx_coord,
                                   vtx_ln_to_gn ,
                                   PDM_MESH_NODAL_BAR2,
                                   pn_edge,
                                   edge_vtx    ,
                                   edge_ln_to_gn,
                                   1,
                                   field_name,
                                   edge_field);

        PDM_free(edge_flags);
      }
    }




  }




  PDM_multipart_free(mpart_id);
  PDM_dcube_nodal_gen_free(dcube);
  PDM_part_mesh_nodal_free(pmesh_nodal);
  PDM_part_mesh_free(pm);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}


