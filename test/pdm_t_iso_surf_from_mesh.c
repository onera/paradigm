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
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_distrib.h"
#include "pdm_multipart.h"
#include "pdm_dist_cloud_surf.h"
#include "pdm_writer.h"
#include "pdm_iso_surface.h"

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
_read_args(int                    argc,
           char                 **argv,
           PDM_g_num_t           *n_vtx_seg,
           char                 **filename,
           PDM_Mesh_nodal_elt_t  *elt_type,
           double                *level,
           int                   *n_part,
           int                   *part_method)
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
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *level = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *elt_type = (PDM_Mesh_nodal_elt_t) atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
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


static void
_compute_distance_field_and_gradient
(
 PDM_MPI_Comm    comm,
 int             n_part,
 int            *n_vtx,
 double        **vtx_coord,
 PDM_g_num_t   **vtx_ln_to_gn,
 int             sm_n_face,
 int             sm_n_vtx,
 int            *sm_face_vtx_idx,
 int            *sm_face_vtx,
 double         *sm_vtx_coord,
 PDM_g_num_t    *sm_face_ln_to_gn,
 PDM_g_num_t    *sm_vtx_ln_to_gn,
 double       ***field,
 double       ***gradient
 )
{
  int n_point_cloud = 1;

  PDM_dist_cloud_surf_t *dist = PDM_dist_cloud_surf_create (PDM_MESH_NATURE_MESH_SETTED,
                                                            n_point_cloud,
                                                            comm,
                                                            PDM_OWNERSHIP_KEEP);

  int sm_n_part = 1;
  PDM_dist_cloud_surf_surf_mesh_global_data_set (dist,
                                                 sm_n_part);

  PDM_dist_cloud_surf_surf_mesh_part_set (dist,
                                          0,
                                          sm_n_face,
                                          sm_face_vtx_idx,
                                          sm_face_vtx,
                                          sm_face_ln_to_gn,
                                          sm_n_vtx,
                                          sm_vtx_coord,
                                          sm_vtx_ln_to_gn);


  PDM_dist_cloud_surf_n_part_cloud_set (dist,
                                        0,
                                        n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_dist_cloud_surf_cloud_set (dist,
                                   0,
                                   i_part,
                                   n_vtx[i_part],
                                   vtx_coord[i_part],
                                   vtx_ln_to_gn[i_part]);
  }

  PDM_dist_cloud_surf_compute (dist);

  *field    = (double **) malloc(sizeof(double *) * n_part);
  *gradient = (double **) malloc(sizeof(double *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    double      *distance         = NULL;
    double      *projected        = NULL;
    PDM_g_num_t *closest_elt_gnum = NULL;
    PDM_dist_cloud_surf_get (dist,
                             0,
                             i_part,
                             &distance,
                             &projected,
                             &closest_elt_gnum);

    (*field)[i_part]    = (double *) malloc(sizeof(double) * n_vtx[i_part]);
    (*gradient)[i_part] = (double *) malloc(sizeof(double) * n_vtx[i_part] * 3);

    for (int i = 0; i < n_vtx[i_part]; i++) {

      (*field)[i_part][i] = sqrt(distance[i]);

      double vec[3];
      for (int j = 0; j < 3; j++) {
        vec[j] = vtx_coord[i_part][3*i + j] - projected[3*i + j];
      }

      double mag = PDM_MODULE(vec);
      if (mag > 1.e-16) {
        mag = 1. / mag;
        for (int j = 0; j < 3; j++) {
          (*gradient)[i_part][3*i + j] = vec[j] * mag;
        }
      }

    }
  }

  PDM_dist_cloud_surf_free (dist);
}


static PDM_multipart_t *
_generate_volume_mesh
(
 PDM_MPI_Comm                  comm,
 const PDM_g_num_t             n_vtx_seg,
 const double                  g_extents[],
 const PDM_Mesh_nodal_elt_t    elt_type,
 const PDM_split_dual_t        part_method,
 const int                     n_part
 )
{
  double length = 0.;
  for (int i = 0; i < 3; i++) {
    length = PDM_MAX(length, g_extents[3+i] - g_extents[i]);
  }

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         g_extents[0],
                                                         g_extents[1],
                                                         g_extents[2],
                                                         elt_type,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);


  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);
  PDM_dmesh_nodal_generate_distribution(dmn);

  /*
   * Partitionnement
   */
  int n_zone = 1;
  // int n_part_zones = {n_part};
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;
  printf("n_part = %d\n", n_part);
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                n_part_zones,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);
  PDM_multipart_run_ppart(mpart);

  PDM_dcube_nodal_gen_free(dcube);

  return mpart;
}




static void
_read_surface_mesh
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int           *n_vtx,
 int           *n_face,
 int          **face_vtx_idx,
 int          **face_vtx,
 double       **vtx_coord,
 PDM_g_num_t  **face_ln_to_gn,
 PDM_g_num_t  **vtx_ln_to_gn
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // assert(n_rank == 1);

  if (i_rank == 0) {
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
    }

    PDM_g_num_t gn_vtx, gn_face;

    fscanf(f, PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n", &gn_vtx, &gn_face);

    *n_vtx  = (int) gn_vtx;
    *n_face = (int) gn_face;

    if (i_rank == 0) {
      printf("gsm_n_vtx = "PDM_FMT_G_NUM", gsm_n_face = "PDM_FMT_G_NUM"\n", gn_vtx, gn_face);
    }

    *vtx_coord = malloc (sizeof(double) * (*n_vtx) * 3);
    for (int i = 0; i < *n_vtx; i++) {
      fscanf(f, "%lf %lf %lf\n",
             *vtx_coord + 3*i, *vtx_coord + 3*i+1, *vtx_coord + 3*i+2);
    }

    *face_vtx_idx = malloc (sizeof(int) * (*n_face + 1));
    for (int i = 0; i <= *n_face; i++) {
      fscanf(f, "%d", *face_vtx_idx + i);
    }

    *face_vtx = malloc (sizeof(int) * (*face_vtx_idx)[*n_face]);
    for (int i = 0; i < (*face_vtx_idx)[*n_face]; i++) {
      fscanf(f, "%d", *face_vtx + i);
    }

    fclose(f);



    *vtx_ln_to_gn  = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*n_vtx));
    for (int i = 0; i < *n_vtx; i++) {
      (*vtx_ln_to_gn)[i] = i + 1; //use distrib
    }

    *face_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*n_face));
    for (int i = 0; i < *n_face; i++) {
      (*face_ln_to_gn)[i] = i + 1; //use distrib
    }
  }

  else {
    *n_vtx  = 0;
    *n_face = 0;

    *vtx_coord = malloc (sizeof(double) * (*n_vtx) * 3);
    *face_vtx_idx = malloc (sizeof(int) * (*n_face + 1));
    (*face_vtx_idx)[0] = 0;
    *face_vtx = malloc (sizeof(int) * (*face_vtx_idx)[*n_face]);

    *vtx_ln_to_gn  = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*n_vtx));
    *face_ln_to_gn = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (*n_face));
  }

}


static void
_compute_face_vtx
(
 const int   n_face,
 int        *pface_edge_idx,
 int        *pface_edge,
 int        *pedge_vtx,
 int       **pface_vtx
 )
{
  *pface_vtx = (int *) malloc(sizeof(int) * pface_edge_idx[n_face]);

  for (int i = 0; i < n_face; i++) {

    int *_pface_vtx = *pface_vtx + pface_edge_idx[i];

    int cur_vtx, next_vtx;
    int cur_edge = pface_edge[pface_edge_idx[i]];
    if (cur_edge < 0) {
      cur_edge = -cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge+1];
      next_vtx = pedge_vtx[2*cur_edge  ];
    } else {
      cur_edge = cur_edge - 1;
      cur_vtx  = pedge_vtx[2*cur_edge  ];
      next_vtx = pedge_vtx[2*cur_edge+1];
    }

    for (int ivtx = 0; ivtx < pface_edge_idx[i+1] - pface_edge_idx[i]; ivtx++) {
      _pface_vtx[ivtx] = cur_vtx;

      for (int iedg = pface_edge_idx[i]; iedg < pface_edge_idx[i+1]; iedg++) {
        cur_edge = pface_edge[iedg];
        int vtx1, vtx2;
        if (cur_edge < 0) {
          cur_edge = -cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge+1];
          vtx2 = pedge_vtx[2*cur_edge  ];
        } else {
          cur_edge = cur_edge - 1;
          vtx1 = pedge_vtx[2*cur_edge  ];
          vtx2 = pedge_vtx[2*cur_edge+1];
        }

        if (vtx1 == next_vtx) {
          cur_vtx  = next_vtx;
          next_vtx = vtx2;
          break;
        }
      }
    }

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

  PDM_g_num_t           n_vtx_seg = 10;
  char                 *filemesh  = NULL;
  PDM_Mesh_nodal_elt_t  elt_type  = PDM_MESH_NODAL_TETRA4;
  double                level     = 1.e-2;
  int                   n_part    = 1;
  #ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &filemesh,
             &elt_type,
             &level,
             &n_part,
     (int *) &part_method);

  assert(filemesh != NULL);

  /*
   *  Init
   */

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   *  Read surface mesh
   */
  if (i_rank == 0) {
    printf("-- Read surface mesh\n");
    fflush(stdout);
  }
  int          sm_n_vtx  = 0;
  int          sm_n_face = 0;
  int         *sm_face_vtx_idx  = NULL;
  int         *sm_face_vtx      = NULL;
  double      *sm_vtx_coord     = NULL;
  PDM_g_num_t *sm_face_ln_to_gn = NULL;
  PDM_g_num_t *sm_vtx_ln_to_gn  = NULL;

  _read_surface_mesh(comm,
                     filemesh,
                     &sm_n_vtx,
                     &sm_n_face,
                     &sm_face_vtx_idx,
                     &sm_face_vtx,
                     &sm_vtx_coord,
                     &sm_face_ln_to_gn,
                     &sm_vtx_ln_to_gn);

  if (0) {
    char filename[999];
    sprintf(filename, "surface_mesh_%2.2d.vtk", i_rank);
    PDM_vtk_write_polydata(filename,
                           sm_n_vtx,
                           sm_vtx_coord,
                           sm_vtx_ln_to_gn,
                           sm_n_face,
                           sm_face_vtx_idx,
                           sm_face_vtx,
                           sm_face_ln_to_gn,
                           NULL);
  }



  /*
   *  Compute global extents of surface mesh
   */
  double l_extents[6] = {
    HUGE_VAL, HUGE_VAL, HUGE_VAL,
    -HUGE_VAL, -HUGE_VAL, -HUGE_VAL
  };

  for (int i = 0; i < sm_n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      double x = sm_vtx_coord[3*i + j];
      l_extents[j]   = PDM_MIN(l_extents[j],   x);
      l_extents[3+j] = PDM_MAX(l_extents[3+j], x);
    }

  }

  double g_extents[6];
  PDM_MPI_Allreduce(l_extents,   g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(l_extents+3, g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double eps = 1e-2;
  double g_range = 0.;
  double g_center[3];
  for (int j = 0; j < 3; j++) {
    g_center[j] = 0.5*(g_extents[3+j] + g_extents[j]);
    g_range = PDM_MAX(g_range, g_extents[3+j] - g_extents[j]);
  }

  g_range = g_range * (0.5 + eps);
  for (int j = 0; j < 3; j++) {
    g_extents[j]   = g_center[j] - g_range;
    g_extents[j+3] = g_center[j] + g_range;
  }

  printf("g_extents = %f %f %f  %f %f %f\n",
         g_extents[0], g_extents[1], g_extents[2],
         g_extents[3], g_extents[4], g_extents[5]);

  /*
   *  Generate volume mesh
   */
  if (i_rank == 0) {
    printf("-- Generate volume mesh\n");
    fflush(stdout);
  }
  PDM_multipart_t *mpart = _generate_volume_mesh (comm,
                                                  n_vtx_seg,
                                                  g_extents,
                                                  elt_type,
                                                  part_method,
                                                  n_part);


  /*
   *  Compute distance function and gradient
   */
  if (i_rank == 0) {
    printf("-- Compute distance field and gradient\n");
    fflush(stdout);
  }

  double **field    = NULL;
  double **gradient = NULL;

  int          *pn_vtx        = (int *)          malloc(sizeof(int)           * n_part);
  double      **pvtx_coord    = (double **)      malloc(sizeof(double *)      * n_part);
  PDM_g_num_t **pvtx_ln_to_gn = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int i_part = 0; i_part < n_part; i_part++) {
    pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart,
                                                      0,
                                                      i_part,
                                                      pvtx_coord + i_part,
                                                      PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    pvtx_ln_to_gn + i_part,
                                    PDM_OWNERSHIP_KEEP);
  }

  _compute_distance_field_and_gradient (comm,
                                        n_part,
                                        pn_vtx,
                                        pvtx_coord,
                                        pvtx_ln_to_gn,
                                        sm_n_face,
                                        sm_n_vtx,
                                        sm_face_vtx_idx,
                                        sm_face_vtx,
                                        sm_vtx_coord,
                                        sm_face_ln_to_gn,
                                        sm_vtx_ln_to_gn,
                                        &field,
                                        &gradient);

  for (int i_part = 0; i_part < n_part; i_part++) {
    for (int i = 0; i < pn_vtx[i_part]; i++) {
      field[i_part][i] -= level;
    }
  }



  if (1) {
    PDM_writer_t *id_cs = PDM_writer_create("Ensight",
                                            PDM_WRITER_FMT_ASCII,
                                            PDM_WRITER_TOPO_CST,
                                            PDM_WRITER_OFF,
                                            "test_isosurf3d",
                                            "isosurf3d",
                                            PDM_MPI_COMM_WORLD,
                                            PDM_IO_KIND_MPI_SIMPLE,
                                            1.,
                                            NULL);

    int id_geom = PDM_writer_geom_create(id_cs,
                                         "isosurf3d_geom",
                                         n_part);

    int id_var_field = PDM_writer_var_create(id_cs,
                                             PDM_WRITER_ON,
                                             PDM_WRITER_VAR_SCALAR,
                                             PDM_WRITER_VAR_VERTICES,
                                             "field");

    int id_var_gradient = PDM_writer_var_create(id_cs,
                                                PDM_WRITER_ON,
                                                PDM_WRITER_VAR_VECTOR,
                                                PDM_WRITER_VAR_VERTICES,
                                                "gradient");

    PDM_writer_step_beg(id_cs, 0.);


    int **face_vtxNb  = (int **) malloc(sizeof(int *) * n_part);
    int **cell_faceNb = (int **) malloc(sizeof(int *) * n_part);
    int **face_vtx    = (int **) malloc(sizeof(int *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {

     int *cell_face_idx;
     int *cell_face;
     int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                      0,
                                                      i_part,
                                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                      &cell_face,
                                                      &cell_face_idx,
                                                      PDM_OWNERSHIP_KEEP);

      int *face_edge_idx;
      int *face_edge;
      int n_face = PDM_multipart_part_connectivity_get(mpart,
                                    0,
                                    i_part,
                                    PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                    &face_edge,
                                    &face_edge_idx,
                                    PDM_OWNERSHIP_KEEP);
      int *edge_vtx_idx;
      int *edge_vtx;
      PDM_multipart_part_connectivity_get(mpart,
                                          0,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &edge_vtx,
                                          &edge_vtx_idx,
                                          PDM_OWNERSHIP_KEEP);


      PDM_g_num_t *cell_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &cell_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *face_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &face_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);

      PDM_g_num_t *edge_ln_to_gn;
      PDM_multipart_part_ln_to_gn_get(mpart,
                                      0,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &edge_ln_to_gn,
                                      PDM_OWNERSHIP_KEEP);



      face_vtxNb[i_part]  = (int *) malloc(sizeof(int) * n_face);
      cell_faceNb[i_part] = (int *) malloc(sizeof(int) * n_cell);

      _compute_face_vtx(n_face,
                        face_edge_idx,
                        face_edge,
                        edge_vtx,
                        &face_vtx[i_part]);

      for (int i = 0; i < n_cell; i++) {
        cell_faceNb[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      for (int i = 0; i < n_face; i++) {
        face_vtxNb[i_part][i] = face_edge_idx[i+1] - face_edge_idx[i];
      }

      PDM_writer_geom_coord_set(id_cs,
                                id_geom,
                                i_part,
                                pn_vtx[i_part],
                                pvtx_coord[i_part],
                                pvtx_ln_to_gn[i_part]);

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                           id_geom,
                                           i_part,
                                           n_cell,
                                           n_face,
                                           face_edge_idx,
                                           face_vtxNb[i_part],
                                           face_vtx[i_part],
                                           cell_face_idx,
                                           cell_faceNb[i_part],
                                           cell_face,
                                           cell_ln_to_gn);
    }

    PDM_writer_geom_write(id_cs,
                          id_geom);

    PDM_real_t **val_field    = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_gradient = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
    for (int i_part = 0; i_part < n_part; i_part++) {

      val_field[i_part]    = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_vtx[i_part]);
      val_gradient[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * pn_vtx[i_part] * 3);

      for (int i = 0; i < pn_vtx[i_part]; i++) {
        val_field[i_part][i] = (PDM_real_t) field[i_part][i];
        for (int j = 0; j < 3; j++) {
          val_gradient[i_part][3*i+j] = (PDM_real_t) gradient[i_part][3*i+j];
        }
      }

      PDM_writer_var_set(id_cs,
                         id_var_field,
                         id_geom,
                         i_part,
                         val_field[i_part]);

      PDM_writer_var_set(id_cs,
                         id_var_gradient,
                         id_geom,
                         i_part,
                         val_gradient[i_part]);
    }

    PDM_writer_var_write(id_cs,
                         id_var_field);

    PDM_writer_var_write(id_cs,
                         id_var_gradient);


    PDM_writer_step_end(id_cs);

    PDM_writer_free(id_cs);

    for (int i = 0; i < n_part; i++) {
      free(val_field[i]);
      free(val_gradient[i]);
      free(face_vtxNb[i]);
      free(cell_faceNb[i]);
    }
    free(val_field);
    free(val_gradient);
    free(face_vtxNb);
    free(cell_faceNb);
  }


  /*
   *  Compute iso-surface from partitionned mesh
   */
  PDM_iso_surface_t* isos = PDM_iso_surface_create(3,
                                                   PDM_ISO_SURFACE_KIND_FIELD,
                                                   1,
                                                   PDM_OWNERSHIP_KEEP,
                                                   comm);


  for (int i_part = 0; i_part < n_part; i_part++) {

    int *cell_face_idx;
    int *cell_face;
    int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                     &cell_face,
                                                     &cell_face_idx,
                                                     PDM_OWNERSHIP_KEEP);

    int *face_edge_idx;
    int *face_edge;
    int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                     &face_edge,
                                                     &face_edge_idx,
                                                     PDM_OWNERSHIP_KEEP);

    int *edge_vtx_idx;
    int *edge_vtx;
    int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                     0,
                                                     i_part,
                                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                     &edge_vtx,
                                                     &edge_vtx_idx,
                                                     PDM_OWNERSHIP_KEEP);


    PDM_g_num_t *cell_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    &cell_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *face_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    &face_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *edge_ln_to_gn;
    PDM_multipart_part_ln_to_gn_get(mpart,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    &edge_ln_to_gn,
                                    PDM_OWNERSHIP_KEEP);

    PDM_iso_surface_part_set (isos,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              pn_vtx[i_part],
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              NULL,//face_vtx_idx,
                              NULL,//face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              pvtx_ln_to_gn[i_part],
                              pvtx_coord[i_part]);

    PDM_iso_surface_part_field_set (isos,
                                    i_part,
                                    field[i_part]);

    PDM_iso_surface_part_gradient_field_set (isos,
                                             i_part,
                                             gradient[i_part]);
  }


  PDM_iso_surface_compute(isos);


  char name[999];
  sprintf(name, "iso_distance_l_%3.3f", level);
  PDM_iso_surface_write(isos, name);

  PDM_iso_surface_free(isos);

  /*
   *  Free memory
   */
  for (int i = 0; i < n_part; i++) {
    // free(cell_face_idx[i]);
    // free(cell_face[i]);
    // free(face_edge_idx[i]);
    // free(face_edge[i]);
    // free(edge_vtx[i]);
    // free(vtx_coord[i]);
    // free(cell_ln_to_gn[i]);
    // free(face_ln_to_gn[i]);
    // free(edge_ln_to_gn[i]);
    // free(vtx_ln_to_gn[i]);
    // free(field[i]);
    // free(gradient[i]);
  }
  // free(n_cell);
  // free(n_face);
  // free(n_edge);
  // free(n_vtx);
  // free(cell_face_idx);
  // free(cell_face);
  // free(face_edge_idx);
  // free(face_edge);
  // free(edge_vtx);
  // free(vtx_coord);
  // free(cell_ln_to_gn);
  // free(face_ln_to_gn);
  // free(edge_ln_to_gn);
  // free(vtx_ln_to_gn);
  free(field);
  free(gradient);

  PDM_multipart_free(mpart);

  free(sm_face_vtx_idx);
  free(sm_face_vtx);
  free(sm_vtx_coord);
  free(sm_face_ln_to_gn);
  free(sm_vtx_ln_to_gn);


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();
}
