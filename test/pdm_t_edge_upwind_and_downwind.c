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
#include "pdm_plugin.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_polygon.h"
#include "pdm_array.h"

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
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length,
           int           *n_part,
           int           *post,
           int           *part_method)
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static
void
_setup_edge_upwind_and_downwind
(
  int     n_cell,
  int     n_face,
  int     n_edge,
  int     n_vtx,
  int    *cell_face_idx,
  int    *cell_face,
  int    *face_vtx_idx,
  int    *face_vtx,
  int    *vtx_cell_idx,
  int    *vtx_cell,
  int    *edge_vtx,
  double *vtx_coord
)
{
  /*
   *  Compute normale and center_interface
   */
  double* face_center = (double *) malloc(3 * n_face * sizeof(double));
  double* face_surf   = (double *) malloc(3 * n_face * sizeof(double));
  for(int i_face = 0; i_face < n_face; ++i_face) {

    face_center[3*i_face  ] = 0.;
    face_center[3*i_face+1] = 0.;
    face_center[3*i_face+2] = 0.;

    face_surf  [3*i_face  ] = 0.;
    face_surf  [3*i_face+1] = 0.;
    face_surf  [3*i_face+2] = 0.;

    double pond = 1./(face_vtx_idx[i_face+1] - face_vtx_idx[i_face]);
    int n_vtx_on_face = face_vtx_idx[i_face+1] - face_vtx_idx[i_face];

    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = PDM_ABS(face_vtx[idx_vtx])-1;
      face_center[3*i_face  ] += vtx_coord[3*i_vtx  ];
      face_center[3*i_face+1] += vtx_coord[3*i_vtx+1];
      face_center[3*i_face+2] += vtx_coord[3*i_vtx+2];
    }
    face_center[3*i_face  ] = face_center[3*i_face  ] * pond;
    face_center[3*i_face+1] = face_center[3*i_face+1] * pond;
    face_center[3*i_face+2] = face_center[3*i_face+2] * pond;

    double xf = face_center[3*i_face  ];
    double yf = face_center[3*i_face+1];
    double zf = face_center[3*i_face+2];

    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {

      const int i_vtx1 = face_vtx[idx_vtx] - 1;
      const int i_vtx2 = face_vtx[(idx_vtx + 1) % n_vtx_on_face] - 1;

      double x2 = vtx_coord[3*i_vtx1  ];
      double y2 = vtx_coord[3*i_vtx1+1];
      double z2 = vtx_coord[3*i_vtx1+2];

      double x3 = vtx_coord[3*i_vtx2  ];
      double y3 = vtx_coord[3*i_vtx2+1];
      double z3 = vtx_coord[3*i_vtx2+2];

      double ux = x2 - xf;
      double uy = y2 - yf;
      double uz = z2 - zf;

      double vx = x3 - xf;
      double vy = y3 - yf;
      double vz = z3 - zf;

      face_surf[3*i_face  ] += 0.5 * (uy * vz - uz * vy);
      face_surf[3*i_face+1] += 0.5 * (uz * vx - ux * vz);
      face_surf[3*i_face+2] += 0.5 * (ux * vy - uy * vx);
    }
  }

  /*
   * Begin algo
   */

  int *idx_vtx_to_reset = malloc( n_vtx * sizeof(int));
  int *vtx_flags        = malloc( n_vtx * sizeof(int));

  int *idx_face_to_reset = malloc( n_face * sizeof(int));
  int *face_flags        = malloc( n_face * sizeof(int));

  int *idx_cell_to_reset = malloc( n_face * sizeof(int) );
  int *idx_cell_keep     = malloc( n_cell * sizeof(int) );
  int *idx_face_keep     = malloc( n_face * sizeof(int) );

  double *pts_in_face    = malloc( 3*n_face * sizeof(int));

  double epsilon = 1.e-12;
  for(int i_edge = 0; i_edge < n_edge; ++i_edge) {

    int i_vtx1 = edge_vtx[2*i_edge  ];
    int i_vtx2 = edge_vtx[2*i_edge+1];

    double edge_dirx = vtx_coord[3*i_vtx2  ] - vtx_coord[3*i_vtx1  ];
    double edge_diry = vtx_coord[3*i_vtx2+1] - vtx_coord[3*i_vtx1+1];
    double edge_dirz = vtx_coord[3*i_vtx2+2] - vtx_coord[3*i_vtx1+2];

    log_trace(" -------------------------------------------- \n");
    for(int idx_vtx = 0; idx_vtx < 2; ++idx_vtx) {

      int lvtx    = 0;
      int lface   = 0;
      int lselect = 0;

      int i_vtx = edge_vtx[2*i_edge+idx_vtx]-1;

      for(int idx_cell = vtx_cell_idx[i_vtx]; idx_cell < vtx_cell_idx[i_vtx+1]; ++idx_cell) {

        int i_cell = PDM_ABS(vtx_cell[idx_cell]) - 1;

        for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {
          int i_face = PDM_ABS(cell_face[idx_face])-1;
          if(face_flags[i_face] == 1) {
            continue;
          }

          // Filrage des faces pour prendre qui contienne le sommet courant
          // Here : Intersection of line with the plane
          int is_a_face_candidate = 1;
          for(int idx_v = face_vtx_idx[i_face]; idx_v < face_vtx_idx[i_face+1]; ++idx_v) {
            int i_vtxc = PDM_ABS(face_vtx[idx_v]) - 1;
            if(i_vtx == i_vtxc) {
              is_a_face_candidate = 0;
              break;
            }
          }

          if(is_a_face_candidate == 1) {
            idx_face_to_reset[lface ] = i_face;
            face_flags       [i_face] = 1;
            idx_cell_to_reset[lface ] = i_cell;
            lface++;
          }
        } /* End loop face */
      } /* End loop cell */

      PDM_log_trace_array_int(idx_face_to_reset, lface, "idx_face_to_reset ::");
      PDM_log_trace_array_int(idx_cell_to_reset, lface, "idx_cell_to_reset ::");

      /* Reloop to make up the test */
      for(int idx_face = 0; idx_face < lface; ++idx_face) {
        int i_face = idx_face_to_reset[idx_face];

        double xf = face_center[3*i_face  ];
        double yf = face_center[3*i_face+1];
        double zf = face_center[3*i_face+2];

        double sn = sqrt(  face_surf[3*i_face  ] * face_surf[3*i_face  ]
                         + face_surf[3*i_face+1] * face_surf[3*i_face+1]
                         + face_surf[3*i_face+2] * face_surf[3*i_face+2]);
        double isn = 1./sn;
        double nx = face_surf[3*i_face  ] * isn;
        double ny = face_surf[3*i_face+1] * isn;
        double nz = face_surf[3*i_face+2] * isn;

        double origin_x = vtx_coord[3*i_vtx1  ];
        double origin_y = vtx_coord[3*i_vtx1+1];
        double origin_z = vtx_coord[3*i_vtx1+2];

        double px = xf - origin_x;
        double py = yf - origin_y;
        double pz = zf - origin_z;

        double denom = edge_dirx * nx + edge_diry * ny + edge_dirz * nz;
        double numer =        px * nx +        py * ny +        pz * nz;  // distance signé entre origin et

        double intersection_x = 0.;
        double intersection_y = 0.;
        double intersection_z = 0.;
        double t              = 0.;

        int stat = -1;
        if(PDM_ABS(denom) < epsilon) {   // Epsilon is here to avoid division by 0
          if(PDM_ABS(numer) < epsilon) { // Le edge contenu dans le plan de la face
            stat = 1;
            intersection_x = origin_x;
            intersection_y = origin_y;
            intersection_z = origin_z;
          } else {   // Edge parallèle mais pas dans le même plan
            stat = 0;
          }
        } else {  // Ca intersecte nickel
          stat = 1;
          t = numer/denom;
          intersection_x = origin_x + t * edge_dirx;
          intersection_y = origin_y + t * edge_diry;
          intersection_z = origin_z + t * edge_dirz;
        }

        int keep = 0;
        if(stat == 1) {

          if((idx_vtx == 0 && t < 0. ) ||
             (idx_vtx == 1 && t > 1.)) {
            keep = 1;
            idx_cell_keep[lselect] = idx_cell_to_reset[idx_face];
            idx_face_keep[lselect] = idx_face_to_reset[idx_face];

            pts_in_face[3*lselect  ] = intersection_x;
            pts_in_face[3*lselect+1] = intersection_y;
            pts_in_face[3*lselect+2] = intersection_z;

            lselect++;
          }
        }
      } /* End loop on selected face */

      PDM_log_trace_array_int(idx_cell_keep, lselect, "idx_cell_keep ::");
      PDM_log_trace_array_int(idx_face_keep, lselect, "idx_face_keep ::");

      // Debug
      for(int idx_face = 0; idx_face < lselect; ++idx_face) {

        int i_face = idx_face_keep[idx_face];

        log_trace("face_vtx = ");
        for(int index_vtx = face_vtx_idx[i_face]; index_vtx < face_vtx_idx[i_face+1]; ++index_vtx) {
          log_trace("%i ", face_vtx[index_vtx]);
        }
        log_trace("\n");

        log_trace("face_center[%i]  = %12.5e / %12.5e / %12.5e \n", i_face, face_center[3*i_face], face_center[3*i_face+1], face_center[3*i_face+2]);
        log_trace("face_surf  [%i]  = %12.5e / %12.5e / %12.5e \n", i_face, face_surf  [3*i_face], face_surf  [3*i_face+1], face_surf  [3*i_face+2]);
        log_trace("pts_in_face[%i]  = %12.5e / %12.5e / %12.5e \n", i_face, pts_in_face[3*idx_face], pts_in_face  [3*idx_face+1], pts_in_face[3*idx_face+2]);
      }

      /* Reset */
      for(int i_reset = 0; i_reset < lface; ++i_reset) {
        face_flags[idx_face_to_reset[i_reset ]] = 0;
      }

    } /* End edge upwind / downwind */
  }



  /*
   *  Output vtk : tag elem by edge ?
   */

  free(idx_vtx_to_reset );
  free(vtx_flags        );
  free(idx_face_to_reset);
  free(face_flags       );
  free(idx_cell_to_reset);
  free(idx_cell_keep    );
  free(idx_face_keep    );
  free(pts_in_face    );
  free(face_center);
  free(face_surf  );
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
  double             length    = 1.;
  int                n_part    = 1;
  int                post      = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length,
             &n_part,
             &post,
     (int *) &part_method);
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
                                                         // PDM_MESH_NODAL_HEXA8,
                                                         PDM_MESH_NODAL_TETRA4,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  if(1 == 1) {
    // PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_VOLUMIC , "out_volumic");
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /*
   * Partitionnement
   */
  int n_zone = 1;
  int n_part_zones = n_part;
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_zone,
               &n_part_zones,
               PDM_FALSE,
               part_method,
               PDM_PART_SIZE_HOMOGENEOUS,
               NULL,
               comm,
               PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_NONE",
                                                     NULL,
                                                     "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart_id, 0, dmn);
  PDM_multipart_run_ppart(mpart_id);

  for (int i_zone = 0; i_zone < n_zone; i_zone++){
    for (int i_part = 0; i_part < n_part; i_part++){
      int  n_proc, tn_part;
      int  n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int  scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int *n_elt;

      PDM_multipart_part_dim_get(mpart_id,
                                 i_zone,
                                 i_part,
                                 &n_section,
                                 &n_elt,
                                 &n_cell,
                                 &n_face,
                                 &n_part_joins,
                                 &n_vtx,
                                 &n_proc,
                                 &tn_part,
                                 &scell_face,
                                 &sface_vtx,
                                 &sface_bound,
                                 &n_bounds,
                                 &sface_join,
                                 &n_joins);
      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, i_zone, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

      int *vtx_part_bound_proc_idx = NULL;
      int *vtx_part_bound_part_idx = NULL;
      int *vtx_part_bound          = NULL;
      PDM_multipart_part_graph_comm_vtx_data_get(mpart_id,
                                                 i_zone,
                                                 i_part,
                                                 &vtx_part_bound_proc_idx,
                                                 &vtx_part_bound_part_idx,
                                                 &vtx_part_bound);

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int n_face2 = PDM_multipart_part_connectivity_get(mpart_id,
                                                        i_zone,
                                                        i_part,
                                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                        &face_edge,
                                                        &face_edge_idx,
                                                        PDM_OWNERSHIP_KEEP);
      assert(n_face == n_face2);
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      int n_edge = PDM_multipart_part_connectivity_get(mpart_id,
                                                       i_zone,
                                                       i_part,
                                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                       &edge_vtx,
                                                       &edge_vtx_idx,
                                                       PDM_OWNERSHIP_KEEP);
      assert(edge_vtx_idx == NULL);
      PDM_g_num_t* edge_ln_to_gn = NULL;
      int n_edge2 = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                    i_zone,
                                                    i_part,
                                                    PDM_MESH_ENTITY_EDGE,
                                                    &edge_ln_to_gn,
                                                    PDM_OWNERSHIP_KEEP);
      assert(n_edge2 == n_edge);

      /*
       *  Compute additionnal connectivity
       */
      assert(face_vtx_idx == NULL);
      assert(face_vtx     == NULL);

      int* tmp_edge_vtx_idx = malloc( (n_edge+1) * sizeof(int));
      tmp_edge_vtx_idx[0] = 0;
      for(int i = 0; i < n_edge; ++i) {
        tmp_edge_vtx_idx[i+1] = tmp_edge_vtx_idx[i] + 2;
      }

      PDM_combine_connectivity(n_face, face_edge_idx, face_edge, tmp_edge_vtx_idx, edge_vtx, &face_vtx_idx, &face_vtx);
      free(tmp_edge_vtx_idx);

      int* cell_vtx_idx = NULL;
      int* cell_vtx     = NULL;
      PDM_combine_connectivity(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, &cell_vtx_idx, &cell_vtx);

      int* vtx_cell_idx = NULL;
      int* vtx_cell     = NULL;
      PDM_connectivity_transpose(n_cell, n_vtx, cell_vtx_idx, cell_vtx, &vtx_cell_idx, &vtx_cell);

      _setup_edge_upwind_and_downwind(n_cell,
                                      n_face,
                                      n_edge,
                                      n_vtx,
                                      cell_face_idx,
                                      cell_face,
                                      face_vtx_idx,
                                      face_vtx,
                                      vtx_cell_idx,
                                      vtx_cell,
                                      edge_vtx,
                                      vtx);



      free(vtx_cell_idx);
      free(vtx_cell);
      free(cell_vtx_idx);
      free(cell_vtx);
      free(face_vtx_idx);
      free(face_vtx);

    }
  }


  PDM_multipart_free(mpart_id);


  PDM_dcube_nodal_gen_free(dcube);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;

}
