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
#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_sort.h"
#include "pdm_distrib.h"
#include "pdm_error.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_polygon.h"
#include "pdm_array.h"
#include "pdm_triangle.h"
#include "pdm_geom_elem.h"

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
 * \param [inout]   part_method Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length,
           int                   *n_part,
           int                   *post,
           int                   *part_method,
           PDM_Mesh_nodal_elt_t  *elt_type)
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

  PDM_g_num_t          n_vtx_seg    = 10;
  double               length       = 1.;
  int                  n_part       = 1;
  int                  post         = 0;
  PDM_split_dual_t     part_method  = PDM_SPLIT_DUAL_WITH_HILBERT;
  PDM_Mesh_nodal_elt_t elt_type     = PDM_MESH_NODAL_TRIA3;

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method,
             &elt_type);
  /*
   *  Init
   */

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

  /*
   *  Create distributed cube
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

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


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if (1) {
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    int dn_vtx = vtx_distrib[i_rank+1] - vtx_distrib[i_rank];

    double noise = 0.2 / (double) (n_vtx_seg - 1);
    for (int i = 0; i < dn_vtx; i++) {
      dvtx_coord[3*i  ] += noise * (rand() / (double) RAND_MAX - 0.5);
      dvtx_coord[3*i+1] += noise * (rand() / (double) RAND_MAX - 0.5);
    }
  }

  if(post) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_RIDGE   , "out_ridge"   );
  }

  /*
   * Partitionnement
   */
  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart, -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);

  for (int i_domain = 0; i_domain < n_domain; i_domain++){
    for (int i_part = 0; i_part < n_part; i_part++){

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                       i_domain,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                       &face_edge_idx,
                                                       &face_edge,
                                                       PDM_OWNERSHIP_KEEP);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                       i_domain,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx_idx,
                                                       &edge_vtx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* face_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* edge_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t* vtx_ln_to_gn = NULL;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      i_domain,
                                      i_part,
                                      PDM_MESH_ENTITY_VTX,
                                      &vtx_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);


      double *vtx = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx,
                                                   PDM_OWNERSHIP_KEEP);

      /*
       *  Compute additionnal connectivity
       */
      PDM_malloc(edge_vtx_idx, n_edge + 1, int);
      for(int i_edge = 0; i_edge < n_edge+1; ++i_edge) {
        edge_vtx_idx[i_edge] = 2 * i_edge;
      }

      int* face_vtx_idx = NULL;
      int* face_vtx     = NULL;
      PDM_combine_connectivity(n_face,
                               face_edge_idx,
                               face_edge,
                               edge_vtx_idx,
                               edge_vtx,
                               &face_vtx_idx,
                               &face_vtx);

      int* vtx_face_idx = NULL;
      int* vtx_face     = NULL;
      PDM_connectivity_transpose(n_face,
                                 n_vtx,
                                 face_vtx_idx,
                                 face_vtx,
                                 &vtx_face_idx,
                                 &vtx_face);



      /* Flag boundary edges (for error-checking) */
      // int *_edge_face_idx = NULL;
      // int *_edge_face     = NULL;
      // PDM_connectivity_transpose(n_face,
      //                            n_edge,
      //                            face_edge_idx,
      //                            face_edge,
      //                            &_edge_face_idx,
      //                            &_edge_face);


      int    *upwind_face_out    = NULL;
      int    *downwind_face_out  = NULL;
      int    *upwind_edge_out    = NULL;
      int    *downwind_edge_out  = NULL;
      double *upwind_point_out   = NULL;
      double *downwind_point_out = NULL;

      // Use to build unit test
      if(0 == 1) {
        PDM_log_trace_array_int   (face_edge_idx, n_face+1             , "face_edge_idx ::");
        PDM_log_trace_array_int   (face_edge    , face_edge_idx[n_face], "face_edge     ::");
        PDM_log_trace_array_int   (edge_vtx    , edge_vtx_idx[n_edge]  , "edge_vtx      ::");
        PDM_log_trace_array_int   (vtx_face_idx, n_vtx+1               , "vtx_face_idx  ::");
        PDM_log_trace_array_int   (vtx_face    , vtx_face_idx[n_vtx]   , "vtx_face      ::");
        PDM_log_trace_array_double(vtx         , 3 * n_vtx             , "vtx_coords    :: ");
      }

      PDM_geom_elem_edge_upwind_and_downwind_2d(0, // XY-plane
                                                face_ln_to_gn,
                                                NULL,
                                                face_edge_idx,
                                                face_edge,
                                                n_edge,
                                                edge_vtx,
                                                vtx_face_idx,
                                                vtx_face,
                                                vtx,
                                                &upwind_face_out,
                                                &downwind_face_out,
                                                &upwind_edge_out,
                                                &downwind_edge_out,
                                                &upwind_point_out,
                                                &downwind_point_out);

      if(0 == 1) {

        PDM_log_trace_array_int   (upwind_face_out   ,     n_edge, "upwind_face_out    ::");
        PDM_log_trace_array_int   (downwind_face_out ,     n_edge, "downwind_face_out  ::");
        PDM_log_trace_array_int   (upwind_edge_out   ,     n_edge, "upwind_edge_out    ::");
        PDM_log_trace_array_int   (downwind_edge_out ,     n_edge, "downwind_edge_out  ::");
        PDM_log_trace_array_double(upwind_point_out  , 3 * n_edge, "upwind_point_out   ::");
        PDM_log_trace_array_double(downwind_point_out, 3 * n_edge, "downwind_point_out ::");

      }


      /* Visualisation */
      if (post) {

        const char *field_name[] = {"upwind_face_out", "downwind_face_out", "upwind_edge_out", "downwind_edge_out"};
        const int  *field     [] = {upwind_face_out, downwind_face_out, upwind_edge_out, downwind_edge_out};

        PDM_vtk_write_std_elements("check_edges.vtk",
                                   n_vtx,
                                   vtx,
                                   vtx_ln_to_gn,
                                   PDM_MESH_NODAL_BAR2,
                                   n_edge,
                                   edge_vtx,
                                   edge_ln_to_gn,
                                   4,
                                   field_name,
                                   field);

        // Prepare dump points
        int n_pts = 0;
        double *pts_intersect = NULL;
        PDM_malloc(pts_intersect, 3 * n_edge * 2, double);
        for(int i_edge = 0; i_edge < n_edge; ++i_edge) {
          if(upwind_edge_out[i_edge] != -1) {
            pts_intersect[3*n_pts  ] = upwind_point_out[3*i_edge  ];
            pts_intersect[3*n_pts+1] = upwind_point_out[3*i_edge+1];
            pts_intersect[3*n_pts+2] = upwind_point_out[3*i_edge+2];
            n_pts++;
          }
          if(downwind_edge_out[i_edge] != -1) {
            pts_intersect[3*n_pts  ] = downwind_point_out[3*i_edge  ];
            pts_intersect[3*n_pts+1] = downwind_point_out[3*i_edge+1];
            pts_intersect[3*n_pts+2] = downwind_point_out[3*i_edge+2];
            n_pts++;
          }
        }

        PDM_vtk_write_point_cloud("check_intersect.vtk",
                                  n_pts,
                                  pts_intersect,
                                  NULL,
                                  NULL);
        PDM_free(pts_intersect);


      }


      PDM_free(upwind_face_out);
      PDM_free(downwind_face_out);
      PDM_free(upwind_edge_out);
      PDM_free(downwind_edge_out);
      PDM_free(upwind_point_out);
      PDM_free(downwind_point_out);


      PDM_free(vtx_face_idx);
      PDM_free(vtx_face);
      PDM_free(face_vtx_idx);
      PDM_free(face_vtx);
      PDM_free(face_vtx_idx);
      PDM_free(face_vtx);
      PDM_free(edge_vtx_idx);

    }
  }


  PDM_multipart_free(mpart);


  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}
