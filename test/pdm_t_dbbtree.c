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
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_inside_cloud_surf.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"

#include "pdm_vtk.h"

#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

struct _pdm_mpi_double_int_t {
  double val;
  int    rank;
};
typedef struct _pdm_mpi_double_int_t PDM_MPI_double_int_t;

static
void
end_timer_and_print(const char* msg, PDM_MPI_Comm comm, double t1){

  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  double t2 = PDM_MPI_Wtime();

  double delta_t = t2 - t1;

  PDM_MPI_double_int_t l_info;

  l_info.val  = delta_t;
  l_info.rank = i_rank;

  PDM_MPI_double_int_t g_max_info;
  PDM_MPI_double_int_t g_min_info;


  PDM_MPI_Allreduce (&l_info,
                     &g_max_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MAXLOC,
                     comm);

  PDM_MPI_Allreduce (&l_info,
                     &g_min_info,
                     1,
                     PDM_MPI_DOUBLE_INT,
                     PDM_MPI_MINLOC,
                     comm);

  if(i_rank == 0) {
    printf("[%i] %s : duration min/max -> %12.5e [on rank = %i] %12.5e [on rank = %i] \n",
           n_rank, msg, g_min_info.val, g_min_info.rank, g_max_info.val, g_max_info.rank);
  }
}

static void _rotate (const int  n_pts,
                     double    *coord)
{
  double R[3][3] = {{0.9362934, -0.2896295, 0.1986693},
                    {0.3129918,  0.9447025, -0.0978434},
                    {-0.1593451,  0.1537920,  0.9751703}};

  for (int i = 0; i < n_pts; i++) {
    double x = coord[3*i];
    double y = coord[3*i+1];
    double z = coord[3*i+2];

    for (int j = 0; j < 3; j++) {
      coord[3*i+j] = R[j][0]*x + R[j][1]*y + R[j][2]*z;
    }
  }
}

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
           PDM_g_num_t   *n_vtx_seg,
           PDM_g_num_t   *n_g_pts_clouds,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
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
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_g_pts_clouds = atol(argv[i]);
        *n_g_pts_clouds = (PDM_g_num_t) _n_g_pts_clouds;
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
      *method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *method = PDM_PART_SPLIT_HILBERT;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}





static
void
_generate_surface_mesh
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
 const PDM_split_dual_t    part_method,
 const int                 n_part,
       PDM_dmesh_nodal_t **_dmn,
       PDM_multipart_t   **_mpart
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_sphere_surf_gen_nodal(comm,
                            nu,
                            nv,
                            x_center,
                            y_center,
                            z_center,
                            radius,
                            &dmn);


  const PDM_g_num_t *distrib_vtx = PDM_DMesh_nodal_distrib_vtx_get (dmn);
  double *dvtx_coord = PDM_DMesh_nodal_vtx_get (dmn);

  int dn_vtx = distrib_vtx[i_rank+1] - distrib_vtx[i_rank];

  // _rotate (dn_vtx, dvtx_coord);


  int n_zone = 1;
  // int n_part_zones = {n_part};
  int *n_part_zones = (int *) malloc(sizeof(int) * n_zone);
  n_part_zones[0] = n_part;

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

  free(n_part_zones);

  *_mpart = mpart;
  *_dmn   = dmn;

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

  PDM_g_num_t           n_vtx_seg      = 100;
  PDM_g_num_t           n_g_pts_clouds = 10;
  double                length         = 1.;
  int                   n_part         = 1;
  int                   post           = 1;

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
             &n_g_pts_clouds,
             &length,
             &n_part,
             &post,
     (int *) &part_method);

  double radius         = length;//2*length;

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
   *  Generate cloud
   */
  if (i_rank == 0) {
    printf("-- Generate volume mesh\n");
    fflush(stdout);
  }
  int          n_pts_clouds;
  double      *pts_coord;
  PDM_g_num_t *pts_g_num;
  PDM_point_cloud_gen_random (comm,
                              n_g_pts_clouds,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_pts_clouds,
                              &pts_coord,
                              &pts_g_num);


  /*
   *  Generate surface mesh
   */
  if (i_rank == 0) {
    printf("-- Generate surface mesh\n");
    fflush(stdout);
  }
  PDM_g_num_t nu = 2*n_vtx_seg;
  PDM_g_num_t nv =   n_vtx_seg;

  PDM_dmesh_nodal_t     *dmn_surf   = NULL;
  PDM_multipart_t       *mpart_surf = NULL;
  _generate_surface_mesh (comm,
                          nu,
                          nv,
                          0.,
                          0.,
                          0.,
                          0.8*length,
                          part_method,
                          n_part,
                          &dmn_surf,
                          &mpart_surf);

  int          *surf_pn_vtx          = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_face         = (int          *) malloc(sizeof(int          ) * n_part);
  int          *surf_pn_edge         = (int          *) malloc(sizeof(int          ) * n_part);
  int         **surf_pface_edge_idx  = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_edge      = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pedge_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  int         **surf_pface_vtx       = (int         **) malloc(sizeof(int         *) * n_part);
  double      **surf_pvtx_coord      = (double      **) malloc(sizeof(double      *) * n_part);
  PDM_g_num_t **surf_pvtx_ln_to_gn   = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
  PDM_g_num_t **surf_pface_ln_to_gn  = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);

  PDM_surf_mesh_t* surf_mesh  = PDM_surf_mesh_create (n_part, comm);
  for (int i_part = 0; i_part < n_part; i_part++) {
    surf_pn_vtx[i_part] = PDM_multipart_part_vtx_coord_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           &surf_pvtx_coord[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                    0,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &surf_pvtx_ln_to_gn[i_part],
                                    PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_ln_to_gn_get(mpart_surf,
                                                           0,
                                                           i_part,
                                                           PDM_MESH_ENTITY_FACE,
                                                           &surf_pface_ln_to_gn[i_part],
                                                           PDM_OWNERSHIP_KEEP);

    surf_pn_face[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                               &surf_pface_edge[i_part],
                                                               &surf_pface_edge_idx[i_part],
                                                               PDM_OWNERSHIP_KEEP);

    int* tmp_pedge_vtx_idx = NULL;
    surf_pn_edge[i_part] = PDM_multipart_part_connectivity_get(mpart_surf,
                                                               0,
                                                               i_part,
                                                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                               &surf_pedge_vtx[i_part],
                                                               &tmp_pedge_vtx_idx,
                                                               PDM_OWNERSHIP_KEEP);
    assert(tmp_pedge_vtx_idx == NULL);


    /*
     * Compute face_vtx
     */
    PDM_compute_face_vtx_from_face_and_edge(surf_pn_face       [i_part],
                                            surf_pface_edge_idx[i_part],
                                            surf_pface_edge    [i_part],
                                            surf_pedge_vtx     [i_part],
                                            &surf_pface_vtx    [i_part]);

    if(0 == 1) {
      PDM_log_trace_array_long(surf_pface_ln_to_gn[i_part], surf_pn_face[i_part]                            , "surf_pface_ln_to_gn :: ");
      PDM_log_trace_array_int (surf_pface_edge_idx[i_part], surf_pn_face[i_part]+1                          , "surf_pface_edge_idx :: ");
      // PDM_log_trace_array_int(surf_pface_vtx[i_part]     , surf_pface_edge_idx[i_part][surf_pn_face[i_part]], "surf_pface_vtx     :: ");
    }

    PDM_surf_mesh_part_input (surf_mesh,
                              i_part,
                              surf_pn_face       [i_part],
                              surf_pface_edge_idx[i_part],
                              surf_pface_vtx     [i_part],
                              surf_pface_ln_to_gn[i_part],
                              surf_pn_vtx        [i_part],
                              surf_pvtx_coord    [i_part],
                              surf_pvtx_ln_to_gn [i_part]);
  }


  /*
   * Prepare for bbox-tree
   */
  int                *part_n_elt       = malloc (sizeof(int          ) * n_part);
  const double      **part_elt_extents = malloc (sizeof(double      *) * n_part);
  const PDM_g_num_t **part_elt_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part);
  const int         **part_elt_vtx_idx = malloc (sizeof(int         *) * n_part);
  const int         **part_elt_vtx     = malloc (sizeof(int         *) * n_part);
  int                *part_n_vtx       = malloc (sizeof(int          ) * n_part);
  const PDM_g_num_t **part_vtx_g_num   = malloc (sizeof(PDM_g_num_t *) * n_part);
  const double      **part_vtx_coord   = malloc (sizeof(double      *) * n_part);


  PDM_surf_mesh_compute_faceExtentsMesh (surf_mesh, 1e-8);
  for (int i_part = 0; i_part < n_part; i_part++) {
    part_n_elt      [i_part] = PDM_surf_mesh_part_n_face_get      (surf_mesh, i_part);
    part_elt_g_num  [i_part] = PDM_surf_mesh_part_face_g_num_get  (surf_mesh, i_part);
    part_elt_extents[i_part] = PDM_surf_mesh_part_extents_get     (surf_mesh, i_part);
    part_elt_vtx_idx[i_part] = PDM_surf_mesh_part_face_vtx_idx_get(surf_mesh, i_part);
    part_elt_vtx    [i_part] = PDM_surf_mesh_part_face_vtx_get    (surf_mesh, i_part);
    part_n_vtx      [i_part] = PDM_surf_mesh_part_n_vtx_get       (surf_mesh, i_part);
    part_vtx_g_num  [i_part] = PDM_surf_mesh_part_vtx_g_num_get   (surf_mesh, i_part);
    part_vtx_coord  [i_part] = PDM_surf_mesh_part_vtx_get         (surf_mesh, i_part);


    if (1 == 1) {
      char filename[999];
      sprintf(filename, "faces_all_%2.2d_%2.2d.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements(filename,
                                part_n_vtx      [i_part],
                                part_vtx_coord  [i_part],
                                NULL,
                                PDM_MESH_NODAL_TRIA3,
                                part_n_elt[i_part],
                                part_elt_vtx[i_part],
                                NULL,
                                0,
                                NULL,
                                NULL);

    }
  }

  /* Compute local extents */
  double my_extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL, -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < part_n_elt[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        my_extents[j]   = PDM_MIN (my_extents[j],   part_elt_extents[ipart][6*i + j]);
        my_extents[j+3] = PDM_MAX (my_extents[j+3], part_elt_extents[ipart][6*i + 3 + j]);
      }
    }
  }

  /* Compute global extents */
  double global_extents[6];
  PDM_MPI_Allreduce (my_extents,   global_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce (my_extents+3, global_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i  ] -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  /*
   * Create ray
   */
  int n_ray = n_pts_clouds;
  double* ray_coord = (double * ) malloc( 2 * 3 * n_ray * sizeof(double));
  PDM_g_num_t* ray_g_num   = pts_g_num;
  for(int i = 0; i < n_ray; ++i) {

    // ray_coord[6*i  ] = -1.;
    // ray_coord[6*i+1] = -1.;
    // ray_coord[6*i+2] = -1.;

    // ray_coord[6*i+3] =  1.;
    // ray_coord[6*i+4] =  1.;
    // ray_coord[6*i+5] =  1.;

    ray_coord[6*i  ] = global_extents[0];
    ray_coord[6*i+1] = global_extents[1];
    ray_coord[6*i+2] = global_extents[2];

    ray_coord[6*i+3] = global_extents[3];
    ray_coord[6*i+4] = global_extents[4];
    ray_coord[6*i+5] = global_extents[5];

    // ray_coord[6*i  ] = -1.;
    // ray_coord[6*i+1] =  1.;
    // ray_coord[6*i+2] = -1.;

    // ray_coord[6*i+3] =  1.;
    // ray_coord[6*i+4] = -1.;
    // ray_coord[6*i+5] =  1.;

  }

  double t1=PDM_MPI_Wtime();

  PDM_dbbtree_t *dbbt = PDM_dbbtree_create (comm, 3, global_extents);

  PDM_box_set_t  *surf_mesh_boxes = PDM_dbbtree_boxes_set_for_intersect_line(dbbt,
                                                                             n_part,
                                                                             part_n_elt,
                                                        (const double **)    part_elt_extents,
                                                    (const PDM_g_num_t **)   part_elt_g_num,
                                                                             n_ray,
                                                                             ray_coord);

  end_timer_and_print("Compute extent and build dbbtree",comm,t1);

  t1 = PDM_MPI_Wtime();

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "box_tree_%i.vtk", i_rank);
    PDM_dbbtree_box_tree_write_vtk(filename,
                                   dbbt,
                                   -1,
                                   0);
  }



  /*
   *  Intersect rays with bounding boxes of surface mesh faces
   */
  int           redistrib_n_part    = 0;
  int          *redistrib_n_box     = NULL;
  PDM_g_num_t **redistrib_box_g_num = NULL;
  int         **box_ray_idx         = NULL;
  PDM_g_num_t **box_ray_g_num       = NULL;
  PDM_dbbtree_lines_intersect_boxes2(dbbt,
                                     n_ray,
                                     ray_g_num,
                                     ray_coord,
                                     &redistrib_n_part,
                                     &redistrib_n_box,
                                     &redistrib_box_g_num,
                                     &box_ray_idx,
                                     &box_ray_g_num);

  end_timer_and_print("PDM_dbbtree_lines_intersect_boxes2 ", comm, t1);

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "ray_%i.vtk", i_rank);
    PDM_vtk_write_lines(filename,
                        n_ray,
                        ray_coord,
                        ray_g_num,
                        NULL);
  }



  free(ray_coord);

  for(int i_part = 0; i_part < redistrib_n_part; ++i_part) {
    free(redistrib_box_g_num[i_part]);
    free(box_ray_idx        [i_part]);
    free(box_ray_g_num      [i_part]);
  }
  free(redistrib_box_g_num);
  free(box_ray_idx        );
  free(box_ray_g_num      );
  free(redistrib_n_box);

  PDM_dbbtree_free    (dbbt);
  PDM_box_set_destroy (&surf_mesh_boxes);

  free(pts_coord);
  free(pts_g_num);
  free(surf_pn_edge       );
  free(surf_pn_vtx        );
  free(surf_pn_face       );
  free(surf_pface_edge    );
  free(surf_pface_edge_idx);
  for(int i_part = 0; i_part < n_part; i_part++) {
    free(surf_pface_vtx[i_part]);
  }
  free(surf_pface_vtx     );
  free(surf_pedge_vtx     );
  free(surf_pvtx_coord    );
  free(surf_pvtx_ln_to_gn );
  free(surf_pface_ln_to_gn);


  free(part_n_elt      );
  free(part_elt_extents);
  free(part_elt_g_num  );
  free(part_elt_vtx_idx);
  free(part_elt_vtx    );
  free(part_n_vtx      );
  free(part_vtx_g_num  );
  free(part_vtx_coord  );

  PDM_multipart_free(mpart_surf);
  PDM_DMesh_nodal_free(dmn_surf);
  PDM_surf_mesh_free(surf_mesh);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
  PDM_printf ("-- End\n");
  }

  PDM_MPI_Finalize ();

  return 0;
}
