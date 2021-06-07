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
#include "pdm_mesh_location.h"
#include "pdm_geom_elem.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_error.h"

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
     "  -t      <level>  Bounding boxes tolerance.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -p      <level>  Number of points to locate.\n\n"
     "  -octree          Use octree-based method.\n\n"
     "  -dbbree          Use dbbtree-based method.\n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}



/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_vtx_seg,
           double        *length,
           double        *separation,
           int           *deform,
           double        *tolerance,
           double        *marge,
           int           *n_part,
           PDM_g_num_t   *n_pts,
           int           *post,
           int           *part_method,
           PDM_mesh_location_method_t *loc_method)
{
  int i = 1;

  PDM_UNUSED (post);

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
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *marge = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-p") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_pts = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_PART_SPLIT_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_PART_SPLIT_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_PART_SPLIT_HILBERT;
    }
    else if (strcmp(argv[i], "-octree") == 0) {
      *loc_method = PDM_MESH_LOCATION_OCTREE;
    }
    else if (strcmp(argv[i], "-dbbtree") == 0) {
      *loc_method = PDM_MESH_LOCATION_DBBTREE;
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static int
_cube_mesh
(
 const int              n_part,
 const PDM_part_split_t part_method,
 const PDM_g_num_t      n_vtx_seg,
 const double           xmin,
 const double           ymin,
 const double           zmin,
 const double           length,
 const int              deform
 )
{
  PDM_dcube_t *dcube = PDM_dcube_gen_init (PDM_MPI_COMM_WORLD,
                                           n_vtx_seg,
                                           length,
                                           xmin,
                                           ymin,
                                           zmin,
                                           PDM_OWNERSHIP_KEEP);

  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  int          n_face_group;
  PDM_g_num_t *dface_cell = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx = NULL;
  double      *dvtx_coord = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group = NULL;
  int          dface_vtx_l;
  int          dface_group_l;


  PDM_dcube_gen_dim_get (dcube,
                         &n_face_group,
                         &dn_cell,
                         &dn_face,
                         &dn_vtx,
                         &dface_vtx_l,
                         &dface_group_l);

  PDM_dcube_gen_data_get (dcube,
                          &dface_cell,
                          &dface_vtx_idx,
                          &dface_vtx,
                          &dvtx_coord,
                          &dface_group_idx,
                          &dface_group);

  if (deform) {
    for (int i = 0; i < dn_vtx; i++) {
      double x = dvtx_coord[3*i];
      double z = dvtx_coord[3*i + 2];

      dvtx_coord[3*i]     += 0.1 * z * z;
      dvtx_coord[3*i + 2] += 0.2 * cos(PDM_PI * x);
    }
  }

  /*
   *  Create mesh partitiions
   */
  int ppart_id = 0;
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));

  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  PDM_part_create (&ppart_id,
                   PDM_MPI_COMM_WORLD,
                   part_method,
                   "PDM_PART_RENUM_CELL_NONE",
                   "PDM_PART_RENUM_FACE_NONE",
                   n_property_cell,
                   renum_properties_cell,
                   n_property_face,
                   renum_properties_face,
                   n_part,
                   dn_cell,
                   dn_face,
                   dn_vtx,
                   n_face_group,
                   NULL,
                   NULL,
                   NULL,
                   NULL,
                   have_dcell_part,
                   dcell_part,
                   dface_cell,
                   dface_vtx_idx,
                   dface_vtx,
                   NULL,
                   dvtx_coord,
                   NULL,
                   dface_group_idx,
                   dface_group);

  free(dcell_part);

  PDM_dcube_gen_free(dcube);

  return ppart_id;
}



static void _export_point_cloud
(
 char         *filename,
 int           n_part,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num,
 PDM_g_num_t **location
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\npoints\nASCII\nDATASET UNSTRUCTURED_GRID\n");

  int n_pts_t = 0;
  for (int ipart = 0; ipart < n_part; ipart++) {
    n_pts_t += n_pts[ipart];
  }

  fprintf(f, "POINTS %d double\n", n_pts_t);
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(f, "%f ", coord[ipart][3*i + j]);
      }
      fprintf(f, "\n");
    }
  }

  fprintf(f, "CELLS %d %d\n", n_pts_t, 2*n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts_t);
  for (int i = 0; i < n_pts_t; i++) {
    fprintf(f, "1\n");
  }

  fprintf(f, "CELL_DATA %d\n", n_pts_t);
  fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
  for (int ipart = 0; ipart < n_part; ipart++) {
    for (int i = 0; i < n_pts[ipart]; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[ipart][i]);
    }
  }

  if (location != NULL) {
    fprintf(f, "FIELD FieldData 1\n");
    fprintf(f, "location 1 %d int\n", n_pts_t);
    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int i = 0; i < n_pts[ipart]; i++) {
        fprintf(f, ""PDM_FMT_G_NUM"\n", location[ipart][i]);
      }
    }
  }

  fclose(f);
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

  PDM_g_num_t n_vtx_seg  = 10;
  double      length     = 1.;
  double      separation = 2.;
  int         deform     = 0;
  double      tolerance  = 1e-6;
  double      marge      = 0.;
  int         n_part     = 1;
  int         post       = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t part_method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t part_method  = PDM_PART_SPLIT_PTSCOTCH;
#else
  PDM_part_split_t part_method  = PDM_PART_SPLIT_HILBERT;
#endif
#endif

  PDM_g_num_t n_pts = 10;
  PDM_mesh_location_method_t loc_method = PDM_MESH_LOCATION_OCTREE;

  /*
   *  Read args
   */

  _read_args (argc,
              argv,
              &n_vtx_seg,
              &length,
              &separation,
              &deform,
              &tolerance,
              &marge,
              &n_part,
              &n_pts,
              &post,
              (int *) &part_method,
              &loc_method);


  /*
   *  Init
   */
  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);


  const double xmin = 0;
  const double ymin = 0;
  const double zmin = 0;

  /*
   *  Source cube
   */
  int ppart_src = _cube_mesh (n_part,
                              part_method,
                              n_vtx_seg,
                              xmin,
                              ymin,
                              zmin,
                              length,
                              deform);

  /*
   *  Target cube
   */
  int ppart_tgt = _cube_mesh (n_part,
                              part_method,
                              n_vtx_seg,
                              xmin + separation*length,
                              ymin,
                              zmin,
                              length,
                              deform);

  /*
   *  Mesh location structure initialization
   */
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED,
                                         1,
                                         PDM_MPI_COMM_WORLD);


  /* Set target point cloud */
  PDM_mesh_location_n_part_cloud_set (mesh_loc,
                                      0,
                                      n_part);

  double **cell_volume = malloc (sizeof(double *) * n_part);
  double **cell_center = malloc (sizeof(double *) * n_part);

  int *n_tgt = malloc (sizeof(double *) * n_part);
  PDM_g_num_t **tgt_g_num = malloc (sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get (ppart_tgt,
                           ipart,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_t_part,
                           &s_cell_face,
                           &s_face_vtx,
                           &s_face_group,
                           &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart_tgt,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    n_tgt[ipart] = n_cell;
    tgt_g_num[ipart] = malloc (sizeof(PDM_g_num_t) * n_cell);
    memcpy (tgt_g_num[ipart], cell_ln_to_gn, sizeof(PDM_g_num_t) * n_cell);

    const int is_oriented = 0;
    cell_volume[ipart] = malloc(sizeof(double) * n_cell);
    cell_center[ipart] = malloc(sizeof(double) * 3 * n_cell);

    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_cell,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        vtx,
                                        cell_volume[ipart],
                                        cell_center[ipart],
                                        NULL,
                                        NULL);

    PDM_mesh_location_cloud_set (mesh_loc,
                                 0,
                                 ipart,
                                 n_cell,
                                 cell_center[ipart],
                                 cell_ln_to_gn);
  }


  /* Set source mesh */
  PDM_mesh_location_mesh_global_data_set (mesh_loc,
                                          n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {

    int n_cell;
    int n_face;
    int n_face_part_bound;
    int n_vtx;
    int n_proc;
    int n_t_part;
    int s_cell_face;
    int s_face_vtx;
    int s_face_group;
    int n_edge_group2;

    PDM_part_part_dim_get (ppart_src,
                           ipart,
                           &n_cell,
                           &n_face,
                           &n_face_part_bound,
                           &n_vtx,
                           &n_proc,
                           &n_t_part,
                           &s_cell_face,
                           &s_face_vtx,
                           &s_face_group,
                           &n_edge_group2);

    int         *cell_tag;
    int         *cell_face_idx;
    int         *cell_face;
    PDM_g_num_t *cell_ln_to_gn;
    int         *face_tag;
    int         *face_cell;
    int         *face_vtx_idx;
    int         *face_vtx;
    PDM_g_num_t *face_ln_to_gn;
    int         *face_part_boundProcIdx;
    int         *face_part_boundPartIdx;
    int         *face_part_bound;
    int         *vtx_tag;
    double      *vtx;
    PDM_g_num_t *vtx_ln_to_gn;
    int         *face_group_idx;
    int         *face_group;
    PDM_g_num_t *face_group_ln_to_gn;

    PDM_part_part_val_get (ppart_src,
                           ipart,
                           &cell_tag,
                           &cell_face_idx,
                           &cell_face,
                           &cell_ln_to_gn,
                           &face_tag,
                           &face_cell,
                           &face_vtx_idx,
                           &face_vtx,
                           &face_ln_to_gn,
                           &face_part_boundProcIdx,
                           &face_part_boundPartIdx,
                           &face_part_bound,
                           &vtx_tag,
                           &vtx,
                           &vtx_ln_to_gn,
                           &face_group_idx,
                           &face_group,
                           &face_group_ln_to_gn);

    PDM_mesh_location_part_set (mesh_loc,
                                ipart,
                                n_cell,
                                cell_face_idx,
                                cell_face,
                                cell_ln_to_gn,
                                n_face,
                                face_vtx_idx,
                                face_vtx,
                                face_ln_to_gn,
                                n_vtx,
                                vtx,
                                vtx_ln_to_gn);
  }


  /* Set location parameters */
  PDM_mesh_location_tolerance_set (mesh_loc,
                                   tolerance);

  PDM_mesh_location_method_set (mesh_loc,
                                loc_method);



  /*
   *  Compute location
   */
  if (i_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  PDM_mesh_location_compute (mesh_loc);

  PDM_mesh_location_dump_times (mesh_loc);


  PDM_g_num_t **tgt_location = malloc (sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    int n_located = PDM_mesh_location_n_located_get (mesh_loc,
                                                     0,//i_point_cloud,
                                                     ipart);

    int *located = PDM_mesh_location_located_get (mesh_loc,
                                                  0,//i_point_cloud,
                                                  ipart);

    int n_unlocated = PDM_mesh_location_n_unlocated_get (mesh_loc,
                                                         0,//i_point_cloud,
                                                         ipart);

    int *unlocated = PDM_mesh_location_unlocated_get (mesh_loc,
                                                      0,//i_point_cloud,
                                                      ipart);

    printf("[%d] n_located = %d, n_unlocated = %d\n", i_rank, n_located, n_unlocated);

    tgt_location[ipart] = PDM_array_const_gnum (n_tgt[ipart], 0);

    PDM_g_num_t *p_location    = NULL;
    double      *p_dist2  = NULL;
    double      *p_proj_coord  = NULL;
    PDM_mesh_location_point_location_get (mesh_loc,
                                          0,//i_point_cloud,
                                          ipart,
                                          &p_location,
                                          &p_dist2,
                                          &p_proj_coord);

    for (int j = 0; j < n_located; j++) {
      int i = located[j] - 1;
      tgt_location[ipart][i] = p_location[j];
    }

    for (int j = 0; j < n_unlocated; j++) {
      int i = unlocated[j] - 1;
      tgt_location[ipart][i] = -1;
    }

    /* Check results */
    if (!deform) {
      const double location_tolerance = 1.e-6;

      const PDM_g_num_t n_cell_seg = n_vtx_seg - 1;
      const double cell_side = length / ((double) n_cell_seg);

      for (int k1 = 0; k1 < n_located; k1++) {
        int ipt = located[k1] - 1;
        double *p = cell_center[ipart] + 3*ipt;

        int i = (int) floor (p[0] / cell_side);
        int j = (int) floor (p[1] / cell_side);
        int k = (int) floor (p[2] / cell_side);

        PDM_g_num_t box_gnum = 1 + i + n_cell_seg*(j + n_cell_seg*k);

        if (p[0] < -tolerance || p[0] > length + tolerance ||
            p[1] < -tolerance || p[1] > length + tolerance ||
            p[2] < -tolerance || p[2] > length + tolerance) {
          box_gnum = -1;
        }

        if (p_location[k1] != box_gnum) {
          /*double *cp = p_proj_coord + 3*ipt;
          printf("%d ("PDM_FMT_G_NUM") (%.15lf %.15lf %.15lf): ("PDM_FMT_G_NUM") | ("PDM_FMT_G_NUM") proj : (%.15lf %.15lf %.15lf)\n",
                 ipt, pts_gnum[ipt],
                 p[0], p[1], p[2],
                 p_location[ipt], box_gnum,
                 cp[0], cp[1], cp[2]);
                 printf("\n");*/

          //-->>
          double cell_min[3] = {cell_side * i,     cell_side * j,     cell_side * k};
          double cell_max[3] = {cell_side * (i+1), cell_side * (j+1), cell_side * (k+1)};
          /*printf("cell min = (%.15lf %.15lf %.15lf)\ncell max = (%.15lf %.15lf %.15lf)\n",
            cell_min[0], cell_min[1], cell_min[2],
            cell_max[0], cell_max[1], cell_max[2]);*/

          double dist = HUGE_VAL;
          for (int idim = 0; idim < 3; idim++) {
            double _dist1 = PDM_ABS (p[idim] - cell_min[idim]);
            double _dist2 = PDM_ABS (p[idim] - cell_max[idim]);
            double _dist = PDM_MIN (_dist1, _dist2);
            dist = PDM_MIN (dist, _dist);
          }
          printf("[%d] "PDM_FMT_G_NUM" distance = %e , location = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM"\n\n",
                 i_rank, tgt_g_num[ipart][ipt], dist, p_location[k1], box_gnum);
          assert (dist < location_tolerance);
          //<<--
        }
      }



      for (int k1 = 0; k1 < n_unlocated; k1++) {
        int ipt = unlocated[k1] - 1;

        double x = cell_center[ipart][3*ipt];
        assert (x > xmin + length);
      }

    }
  }

  if (post) {
    char filename[999];
    sprintf(filename, "tgt_location_%3.3d.vtk", i_rank);

    _export_point_cloud (filename,
                       n_part,
                       n_tgt,
                       cell_center,
                       tgt_g_num,
                       tgt_location);
  }



  for (int ipart = 0; ipart < n_part; ipart++) {
    free (cell_center[ipart]);
    free (cell_volume[ipart]);
    free (tgt_g_num[ipart]);
    free (tgt_location[ipart]);
  }
  free (n_tgt);
  free (cell_center);
  free (cell_volume);
  free (tgt_g_num);
  free (tgt_location);

  PDM_mesh_location_free (mesh_loc,
                          0);

  PDM_part_free (ppart_src);
  PDM_part_free (ppart_tgt);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
