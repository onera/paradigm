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
#include "pdm_writer.h"

#include "pdm_generate_mesh.h"
#include "pdm_mesh_location.h"


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
     "  -h               This message.\n\n");

  exit(exit_code);
}


/*
 * Read arguments from the command line
 */

static void
_read_args
(
 int           argc,
 char        **argv,
 PDM_g_num_t  *src_n_vtx_seg,
 PDM_g_num_t  *tgt_n_vtx_seg,
 int          *src_n_part,
 int          *tgt_n_part,
 double       *tgt_xmin,
 double       *tgt_ymin,
 int          *nodal
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-src_n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long n = atol(argv[i]);
        *src_n_vtx_seg = (PDM_g_num_t) n;
      }
    }
    else if (strcmp(argv[i], "-tgt_n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long n = atol(argv[i]);
        *tgt_n_vtx_seg = (PDM_g_num_t) n;
      }
    }
    else if (strcmp(argv[i], "-src_n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *src_n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tgt_n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tgt_n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tgt_xmin") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tgt_xmin = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-tgt_ymin") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *tgt_ymin = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nodal") == 0) {
      *nodal = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double
_eval_field
(
 const double x,
 const double y,
 const double z
 )
{
  return 1 + 2*x + 3*y + 4*z;
}


static void
_visu_2d
(
 PDM_MPI_Comm    comm,
 const char     *directory,
 const char     *name,
 int             n_part,
 int            *n_vtx,
 double        **vtx_coord,
 PDM_g_num_t   **vtx_ln_to_gn,
 int            *n_face,
 int           **face_vtx_idx,
 int           **face_vtx,
 PDM_g_num_t   **face_ln_to_gn,
 int             n_vtx_field,
 const char    **vtx_field_names,
 double       ***vtx_field_values // [i_field][i_part][i_vtx]
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        directory,
                                        name,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       name,
                                       n_part);

  /* Create variables */
  int id_var_part = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "i_part");

  int id_var_elt_gnum = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "elt_gnum");

  int *id_var_vtx_field = NULL;
  if (n_vtx_field > 0) {
    id_var_vtx_field = malloc(sizeof(int) * n_vtx_field);

    for (int i = 0; i < n_vtx_field; i++) {
      id_var_vtx_field[i] = PDM_writer_var_create(wrt,
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAR,
                                                  PDM_WRITER_VAR_VERTICES,
                                                  vtx_field_names[i]);
    }
  }


  PDM_writer_step_beg(wrt, 0.);

  /* Write geometry */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              i_part,
                              n_vtx       [i_part],
                              vtx_coord   [i_part],
                              vtx_ln_to_gn[i_part],
                              PDM_OWNERSHIP_USER);

    PDM_writer_geom_faces_facesom_add(wrt,
                                      id_geom,
                                      i_part,
                                      n_face       [i_part],
                                      face_vtx_idx [i_part],
                                      NULL,
                                      face_vtx     [i_part],
                                      face_ln_to_gn[i_part]);
  }

  PDM_writer_geom_write(wrt, id_geom);


  /* Write "i_part" and "elt_gnum" variables */
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_real_t *val_part = malloc(sizeof(PDM_real_t) * n_face[i_part]);
    PDM_real_t *val_gnum = malloc(sizeof(PDM_real_t) * n_face[i_part]);

    for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
      val_part[i_face] = i_rank*n_part + i_part;
      val_gnum[i_face] = face_ln_to_gn[i_part][i_face];
    }

    PDM_writer_var_set(wrt,
                       id_var_part,
                       id_geom,
                       i_part,
                       val_part);
    free(val_part);

    PDM_writer_var_set(wrt,
                       id_var_elt_gnum,
                       id_geom,
                       i_part,
                       val_gnum);
    free(val_gnum);
  }


  PDM_writer_var_write(wrt, id_var_part);
  PDM_writer_var_write(wrt, id_var_elt_gnum);


   /* Write node-based variables */
  if (n_vtx_field > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_real_t *val = malloc(sizeof(PDM_real_t) * n_vtx[i_part]);

      for (int i_field = 0; i_field < n_vtx_field; i_field++) {
        for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
          val[i_vtx] = vtx_field_values[i_field][i_part][i_vtx];
        }

        PDM_writer_var_set(wrt,
                           id_var_vtx_field[i_field],
                           id_geom,
                           i_part,
                           val);
      }
      free(val);
    }

    for (int i_field = 0; i_field < n_vtx_field; i_field++) {
      PDM_writer_var_write(wrt, id_var_vtx_field[i_field]);
    }
  }

  PDM_writer_step_end(wrt);


  if (n_vtx_field > 0) {
    free(id_var_vtx_field);
  }

  PDM_writer_free(wrt);
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

  PDM_g_num_t src_n_vtx_seg = 10;
  PDM_g_num_t tgt_n_vtx_seg = 10;
  int         src_n_part    = 1;
  int         tgt_n_part    = 1;
  double      tgt_xmin      = 0.;
  double      tgt_ymin      = 0.;
  int         nodal         = 0;

  /*
   * Read args
   */

  _read_args(argc,
             argv,
             &src_n_vtx_seg,
             &tgt_n_vtx_seg,
             &src_n_part,
             &tgt_n_part,
             &tgt_xmin,
             &tgt_ymin,
             &nodal);


  /*
   * Initialize MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   * Generate and partition the source mesh
   */
  int          *src_n_vtx         = NULL;
  int          *src_n_edge        = NULL;
  int          *src_n_face        = NULL;
  double      **src_vtx_coord     = NULL;
  int         **src_edge_vtx      = NULL;
  int         **src_face_edge_idx = NULL;
  int         **src_face_edge     = NULL;
  int         **src_face_vtx      = NULL;
  PDM_g_num_t **src_vtx_ln_to_gn  = NULL;
  PDM_g_num_t **src_edge_ln_to_gn = NULL;
  PDM_g_num_t **src_face_ln_to_gn = NULL;

  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_POLY_2D,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   src_n_vtx_seg,
                                   src_n_vtx_seg,
                                   src_n_part,
                                   PDM_SPLIT_DUAL_WITH_PARMETIS,
                                   0.2,
                                   &src_n_vtx,
                                   &src_n_edge,
                                   &src_n_face,
                                   &src_vtx_coord,
                                   &src_edge_vtx,
                                   &src_face_edge_idx,
                                   &src_face_edge,
                                   &src_face_vtx,
                                   &src_vtx_ln_to_gn,
                                   &src_edge_ln_to_gn,
                                   &src_face_ln_to_gn);

  /*
   * Then we need to generate and partition a target point cloud
   * For nicer visualization we will generate a second mesh, and use
   * its nodes as the target point cloud.
   */
  int          *tgt_n_vtx         = NULL;
  int          *tgt_n_edge        = NULL;
  int          *tgt_n_face        = NULL;
  double      **tgt_vtx_coord     = NULL;
  int         **tgt_edge_vtx      = NULL;
  int         **tgt_face_edge_idx = NULL;
  int         **tgt_face_edge     = NULL;
  int         **tgt_face_vtx      = NULL;
  PDM_g_num_t **tgt_vtx_ln_to_gn  = NULL;
  PDM_g_num_t **tgt_edge_ln_to_gn = NULL;
  PDM_g_num_t **tgt_face_ln_to_gn = NULL;

  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_QUAD4,
                                   tgt_xmin,
                                   tgt_ymin,
                                   0.,
                                   1.,
                                   1.,
                                   tgt_n_vtx_seg,
                                   tgt_n_vtx_seg,
                                   tgt_n_part,
                                   PDM_SPLIT_DUAL_WITH_HILBERT,
                                   0.,
                                   &tgt_n_vtx,
                                   &tgt_n_edge,
                                   &tgt_n_face,
                                   &tgt_vtx_coord,
                                   &tgt_edge_vtx,
                                   &tgt_face_edge_idx,
                                   &tgt_face_edge,
                                   &tgt_face_vtx,
                                   &tgt_vtx_ln_to_gn,
                                   &tgt_edge_ln_to_gn,
                                   &tgt_face_ln_to_gn);





  /* Create the PDM_mesh_location_t object */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,
                                     tgt_n_part);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    PDM_mesh_location_cloud_set(mesh_loc,
                                0,
                                i_part,
                                tgt_n_vtx       [i_part],
                                tgt_vtx_coord   [i_part],
                                tgt_vtx_ln_to_gn[i_part]);
  }


  /*
   * Set the source mesh
   * Here you have essentially two options :
   *  - you can either define the mesh with *nodal* connectivity (i.e. Finite-Element style)
   *  - or with "descending" connectivity (i.e. Finite-Volume style)
   */
  PDM_mesh_location_mesh_n_part_set(mesh_loc, src_n_part);

  if (nodal) {
    for (int i_part = 0; i_part < src_n_part; i_part++) {
      PDM_mesh_location_nodal_part_set_2d(mesh_loc,
                                          i_part,
                                          src_n_face       [i_part],
                                          src_face_edge_idx[i_part],
                                          src_face_vtx     [i_part],
                                          src_face_ln_to_gn[i_part],
                                          src_n_vtx        [i_part],
                                          src_vtx_coord    [i_part],
                                          src_vtx_ln_to_gn [i_part]);
    }
  }
  else {
    for (int i_part = 0; i_part < src_n_part; i_part++) {
      PDM_mesh_location_part_set_2d(mesh_loc,
                                    i_part,
                                    src_n_face       [i_part],
                                    src_face_edge_idx[i_part],
                                    src_face_edge    [i_part],
                                    src_face_ln_to_gn[i_part],
                                    src_n_edge       [i_part],
                                    src_edge_vtx     [i_part],
                                    src_n_vtx        [i_part],
                                    src_vtx_coord    [i_part],
                                    src_vtx_ln_to_gn [i_part]);
    }
  }

  /* Set the geometric tolerance (optional) */
  double tolerance = 1e-3;
  PDM_mesh_location_tolerance_set(mesh_loc, tolerance);

  /* Set the location preconditioning method (optional) */
  PDM_mesh_location_method_set(mesh_loc,
                               PDM_MESH_LOCATION_OCTREE);

  /* Compute location */
  PDM_mesh_location_compute(mesh_loc);

  /* Dump elapsed and CPU times */
  PDM_mesh_location_dump_times(mesh_loc);



  /*
   * Get the list of (un)located target points
   */
  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                    0,
                                                    i_part);

    int *located = PDM_mesh_location_located_get(mesh_loc,
                                                 0,
                                                 i_part);

    PDM_log_trace_array_int(located, n_located, "located : ");


    int n_unlocated = PDM_mesh_location_n_unlocated_get(mesh_loc,
                                                        0,
                                                        i_part);

    int *unlocated = PDM_mesh_location_unlocated_get(mesh_loc,
                                                     0,
                                                     i_part);

    PDM_log_trace_array_int(unlocated, n_unlocated, "unlocated : ");
  }

  /*
   * Now that we have located the target points in the source mesh,
   * we can exchange data between the two.
   * To complete this exercise, we will interpolate two fields from
   * the source mesh to the target cloud.
   * The first field is cell-based : we can simply use the face global ids for such a field,
   * and check it matches the location data.
   * The second field is node-based : we can use the node coordinates.
   */

  /*
   * First, compute the spatially interpolated fields on the source side.
   * For the first field, the interpolation is straightforward : the target value is simply the same as the host source.
   * The second field interpolation is trickier as you will need the cell->vertex connectivity built during the location computation to link the interpolation weights to the appropriate source nodes.
   */

  /* Interpolate first field (cell-based) */
  double **src_send_field1 = malloc(sizeof(double *) * tgt_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {

    int         *src_to_tgt_idx          = NULL;
    PDM_g_num_t *points_gnum             = NULL;
    double      *points_coords           = NULL;
    double      *points_uvw              = NULL;
    int         *points_weights_idx      = NULL;
    double      *points_weights          = NULL;
    double      *points_dist2            = NULL;
    double      *points_projected_coords = NULL;

    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        0,
                                        i_part,
                                        &src_to_tgt_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords);

    int n_pts = src_to_tgt_idx[src_n_face[i_part]];

    src_send_field1[i_part] = malloc(sizeof(double) * n_pts);
    for (int i_elt = 0; i_elt < src_n_face[i_part]; i_elt++) {
      for (int i_pt = src_to_tgt_idx[i_elt]; i_pt < src_to_tgt_idx[i_elt+1]; i_pt++) {
        src_send_field1[i_part][i_pt] = src_face_ln_to_gn[i_part][i_elt];
      }
    }
  }


  /* Interpolate second field (node-based) */
  double **src_send_field2 = malloc(sizeof(double *) * tgt_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {

    int         *src_to_tgt_idx          = NULL;
    PDM_g_num_t *points_gnum             = NULL;
    double      *points_coords           = NULL;
    double      *points_uvw              = NULL;
    int         *points_weights_idx      = NULL;
    double      *points_weights          = NULL;
    double      *points_dist2            = NULL;
    double      *points_projected_coords = NULL;

    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        0,
                                        i_part,
                                        &src_to_tgt_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords);

    int *cell_vtx_idx = NULL;
    int *cell_vtx     = NULL;
    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      i_part,
                                      &cell_vtx_idx,
                                      &cell_vtx);

    int n_pts = src_to_tgt_idx[src_n_face[i_part]];

    src_send_field2[i_part] = malloc(sizeof(double) * n_pts);
    for (int i_elt = 0; i_elt < src_n_face[i_part]; i_elt++) {
      for (int i_pt = src_to_tgt_idx[i_elt]; i_pt < src_to_tgt_idx[i_elt+1]; i_pt++) {
        src_send_field2[i_part][i_pt] = 0;

        int elt_n_vtx = cell_vtx_idx[i_elt+1] - cell_vtx_idx[i_elt];
        assert(points_weights_idx[i_pt+1] - points_weights_idx[i_pt] == elt_n_vtx);

        for (int i_vtx = 0; i_vtx < elt_n_vtx; i_vtx++) {
          int vtx_id = cell_vtx[cell_vtx_idx[i_elt] + i_vtx] - 1;
          src_send_field2[i_part][i_pt] += src_vtx_coord[i_part][3*vtx_id] * points_weights[points_weights_idx[i_pt] + i_vtx];
        }

      }
    }
  }

  /*
   * Now, use the PartToPart object to exchange the interpolated fields from the source mesh to the target cloud.
   * This ParToPart object was built when computing the location and can be accessed from the MeshLocation object
   */

  /* Get PartToPart object (it is now owned by the user) */
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  /* Initiate exchange of first field */
  int request1 = -1;
  double **tgt_recv_field1 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) src_send_field1,
                         NULL,
              (void ***) &tgt_recv_field1,
                         &request1);

  /* Initiate exchange of second field */
  int request2 = -1;
  double **tgt_recv_field2 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) src_send_field2,
                         NULL,
              (void ***) &tgt_recv_field2,
                         &request2);


  /* Finalize both exchanges */
  PDM_part_to_part_iexch_wait(ptp, request1);
  PDM_part_to_part_iexch_wait(ptp, request2);

  /*
   * Finally, visualize the interpolated target fields.
   * (Beware of unlocated points!)
   */

  double **tgt_field[3];
  tgt_field[0] = malloc(sizeof(double *) * tgt_n_part);
  tgt_field[1] = malloc(sizeof(double *) * tgt_n_part);
  tgt_field[2] = malloc(sizeof(double *) * tgt_n_part);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    tgt_field[0][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);
    tgt_field[1][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);
    tgt_field[2][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);

    double *tgt_field1 = tgt_field[0][i_part];
    double *tgt_field2 = tgt_field[1][i_part];
    double *is_located = tgt_field[2][i_part];

    int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                    0,
                                                    i_part);

    int *located = PDM_mesh_location_located_get(mesh_loc,
                                                 0,
                                                 i_part);

    int n_unlocated = PDM_mesh_location_n_unlocated_get(mesh_loc,
                                                        0,
                                                        i_part);

    int *unlocated = PDM_mesh_location_unlocated_get(mesh_loc,
                                                     0,
                                                     i_part);

    for (int i = 0; i < n_unlocated; i++) {
      int vtx_id = unlocated[i] - 1;
      is_located[vtx_id] = 0;
      tgt_field1[vtx_id] = -1;
      tgt_field2[vtx_id] = -1;
    }

    for (int i = 0; i < n_located; i++) {
      int vtx_id = located[i] - 1;
      is_located[vtx_id] = 1;
      tgt_field1[vtx_id] = tgt_recv_field1[i_part][i];
      tgt_field2[vtx_id] = tgt_recv_field2[i_part][i];

      double error = fabs(tgt_field2[vtx_id] - tgt_vtx_coord[i_part][3*vtx_id]);
      if (error > 1e-9) {
        printf("!! error vtx "PDM_FMT_G_NUM" : %e\n",
               tgt_vtx_ln_to_gn[i_part][vtx_id],
               error);
      }
    }

  }




  const char *field_name[] = {
    "field1",
    "field2",
    "is_located"
  };

  _visu_2d(comm,
           "mesh_location_sol",
           "src_mesh",
           src_n_part,
           src_n_vtx,
           src_vtx_coord,
           src_vtx_ln_to_gn,
           src_n_face,
           src_face_edge_idx,
           src_face_vtx,
           src_face_ln_to_gn,
           0,
           NULL,
           NULL);


  _visu_2d(comm,
           "mesh_location_sol",
           "tgt_mesh",
           tgt_n_part,
           tgt_n_vtx,
           tgt_vtx_coord,
           tgt_vtx_ln_to_gn,
           tgt_n_face,
           tgt_face_edge_idx,
           tgt_face_vtx,
           tgt_face_ln_to_gn,
           3,
           field_name,
           tgt_field);


  /* Free memory */
  PDM_mesh_location_free(mesh_loc);

  PDM_part_to_part_free(ptp);

  for (int i_part = 0; i_part < src_n_part; i_part++) {
    free(src_vtx_coord    [i_part]);
    free(src_edge_vtx     [i_part]);
    free(src_face_edge_idx[i_part]);
    free(src_face_edge    [i_part]);
    free(src_face_vtx     [i_part]);
    free(src_vtx_ln_to_gn [i_part]);
    free(src_edge_ln_to_gn[i_part]);
    free(src_face_ln_to_gn[i_part]);
    free(src_send_field1  [i_part]);
    free(src_send_field2  [i_part]);
  }
  free(src_n_vtx        );
  free(src_n_edge       );
  free(src_n_face       );
  free(src_vtx_coord    );
  free(src_edge_vtx     );
  free(src_face_edge_idx);
  free(src_face_edge    );
  free(src_face_vtx     );
  free(src_vtx_ln_to_gn );
  free(src_edge_ln_to_gn);
  free(src_face_ln_to_gn);
  free(src_send_field1  ); // can be free'd right after PDM_part_to_part_iexch_wait(ptp, &request1);
  free(src_send_field2  ); // can be free'd right after PDM_part_to_part_iexch_wait(ptp, &request2);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    free(tgt_vtx_coord    [i_part]);
    free(tgt_edge_vtx     [i_part]);
    free(tgt_face_edge_idx[i_part]);
    free(tgt_face_edge    [i_part]);
    free(tgt_face_vtx     [i_part]);
    free(tgt_vtx_ln_to_gn [i_part]);
    free(tgt_edge_ln_to_gn[i_part]);
    free(tgt_face_ln_to_gn[i_part]);
    free(tgt_recv_field1  [i_part]);
    free(tgt_recv_field2  [i_part]);
    free(tgt_field[0][i_part]);
    free(tgt_field[1][i_part]);
    free(tgt_field[2][i_part]);
  }
  free(tgt_n_vtx        );
  free(tgt_n_edge       );
  free(tgt_n_face       );
  free(tgt_vtx_coord    );
  free(tgt_edge_vtx     );
  free(tgt_face_edge_idx);
  free(tgt_face_edge    );
  free(tgt_face_vtx     );
  free(tgt_vtx_ln_to_gn );
  free(tgt_edge_ln_to_gn);
  free(tgt_face_ln_to_gn);
  free(tgt_recv_field1  );
  free(tgt_recv_field2  );
  free(tgt_field[0]);
  free(tgt_field[1]);
  free(tgt_field[2]);

  PDM_MPI_Finalize();


  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return EXIT_SUCCESS;
}
