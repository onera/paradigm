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
#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "pdm_poly_vol_gen.h"
#include "pdm_dmesh.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
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
_read_args(int           argc,
           char        **argv,
           PDM_g_num_t  *nx,
           PDM_g_num_t  *ny,
           PDM_g_num_t  *nz,
           double       *lengthx,
           double       *lengthy,
           double       *lengthz,
           int          *n_part,
           int          *post,
           int          *method)
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
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
        *lengthy = atof(argv[i]);
        *lengthz = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthx = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-ly") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthy = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-lz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *lengthz = atof(argv[i]);
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
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_write_point_cloud
(
 const char        *filename,
 const char        *header,
 const int          n_pts,
 const double       coord[],
 const PDM_g_num_t  g_num[],
 const double       field[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "CELLS %d %d\n", n_pts, 2*n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1 %d\n", i);
  }

  fprintf(f, "CELL_TYPES %d\n", n_pts);
  for (int i = 0; i < n_pts; i++) {
    fprintf(f, "1\n");
  }

  if (g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS gnum int\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, ""PDM_FMT_G_NUM"\n", g_num[i]);
    }
  }

  else if (field != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_pts);
    fprintf(f, "SCALARS field double\n LOOKUP_TABLE default\n");
    for (int i = 0; i < n_pts; i++) {
      fprintf(f, "%f\n", field[i]);
    }
  }

  fclose(f);
}


static void
_write_polydata
(
 const char        *filename,
 const char        *header,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const int          face_vtx[],
 const PDM_g_num_t  face_g_num[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " %d", face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
    }
  }


  fclose(f);
}


static void
_write_polydata_gnum
(
 const char        *filename,
 const char        *header,
 const int          n_vtx,
 const double       vtx_coord[],
 const PDM_g_num_t  vtx_g_num[],
 const int          n_face,
 const int          face_vtx_idx[],
 const PDM_g_num_t  face_vtx[],
 const PDM_g_num_t  face_g_num[]
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  if (header != NULL) {
    fprintf(f, "%s\n", header);
  } else {
    fprintf(f, "point cloud\n");
  }
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET POLYDATA\n");

  fprintf(f, "POINTS %d double\n", n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", vtx_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  fprintf(f, "POLYGONS %d %d\n", n_face, n_face + face_vtx_idx[n_face]);
  for (int i = 0; i < n_face; i++) {
    fprintf(f, "%d", face_vtx_idx[i+1] - face_vtx_idx[i]);
    for (int j = face_vtx_idx[i]; j < face_vtx_idx[i+1]; j++) {
      fprintf(f, " "PDM_FMT_G_NUM, face_vtx[j] - 1);
    }
    fprintf(f, "\n");
  }


  if (vtx_g_num != NULL) {
    fprintf(f, "POINT_DATA %d\n", n_vtx);
    fprintf(f, "SCALARS vtx_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_vtx; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", vtx_g_num[i]);
    }
  }

  if (face_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_face);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_face; i++) {
      fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[i]);
    }
  }


  fclose(f);
}


static void
write_dual_edges
(
 const char        *filename,
 const int          n_cell,
 const int         *cell_face_idx,
 const int         *cell_face,
 const double      *cell_center,
 const int          n_face,
 const double      *face_center,
 const PDM_g_num_t *face_g_num
 )
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "dual edges\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", n_cell + n_face);
  for (int i = 0; i < n_cell; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", cell_center[3*i+j]);
    }
    fprintf(f, "\n");
  }

  for (int i = 0; i < n_face; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", face_center[3*i+j]);
    }
    fprintf(f, "\n");
  }

  int n_edge = cell_face_idx[n_cell];
  fprintf(f, "CELLS %d %d\n", n_edge, 3*n_edge);
  for (int i = 0; i < n_cell; i++) {
    for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
      int iface = PDM_ABS(cell_face[j]) - 1;
      fprintf(f, "2 %d %d\n", i, n_cell + iface);
    }
  }

  fprintf(f, "CELL_TYPES %d\n", n_edge);
  for (int i = 0; i < n_edge; i++) {
    fprintf(f, "3\n");
  }

  if (face_g_num != NULL) {
    fprintf(f, "CELL_DATA %d\n", n_edge);
    fprintf(f, "SCALARS face_gnum long 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int iface = PDM_ABS(cell_face[j]) - 1;
        fprintf(f, PDM_FMT_G_NUM"\n", face_g_num[iface]);
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
  PDM_g_num_t nx = 10;
  PDM_g_num_t ny = 10;
  PDM_g_num_t nz = 10;

  double lengthx = 1.;
  double lengthy = 1.;
  double lengthz = 1.;

  int           n_part   = 1;
  int           post    = 0;
#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */
  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &lengthx,
             &lengthy,
             &lengthz,
             &n_part,
             &post,
             (int *) &method);


  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);


  /*
   *  Create distributed mesh
   */
  double xmin = 0.;
  double ymin = 0.;
  double zmin = 0.;

  int randomize = 0;
  int random_seed = 0;

  PDM_g_num_t  ng_cell;
  PDM_g_num_t  ng_face;
  PDM_g_num_t  ng_vtx;
  int          n_face_group;
  int          dn_cell;
  int          dn_face;
  int          dn_vtx;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  PDM_poly_vol_gen (comm,
                    xmin,
                    ymin,
                    zmin,
                    lengthx,
                    lengthy,
                    lengthz,
                    nx,
                    ny,
                    nz,
                    randomize,
                    random_seed,
                    &ng_cell,
                    &ng_face,
                    &ng_vtx,
                    &n_face_group,
                    &dn_cell,
                    &dn_face,
                    &dn_vtx,
                    &dface_cell,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dvtx_coord,
                    &dface_group_idx,
                    &dface_group);

  if (i_rank == 0) printf("ng_cell = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM", ng_vtx = "PDM_FMT_G_NUM"\n", ng_cell, ng_face, ng_vtx);

  if (n_rank == 1) {
    char filename[999];
    sprintf(filename, "dvtx_coord_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "dvtx",
                        dn_vtx,
                        dvtx_coord,
                        NULL,
                        NULL);

    sprintf(filename, "dface_vtx_%3.3d.vtk", i_rank);
    PDM_g_num_t *distrib_face = PDM_compute_entity_distribution (comm, dn_face);
    PDM_g_num_t *face_g_num = malloc (sizeof(PDM_g_num_t) * dn_face);
    for (int i = 0; i < dn_face; i++) {
      face_g_num[i] = 1 + i + distrib_face[i_rank];
    }
    _write_polydata_gnum (filename,
                          "dface",
                          dn_vtx,
                          dvtx_coord,
                          NULL,
                          dn_face,
                          dface_vtx_idx,
                          dface_vtx,
                          face_g_num);
    free (face_g_num);

    PDM_g_num_t *face_data = PDM_array_const_gnum (dn_face, -1);
    for (int i = 0; i < n_face_group; i++) {
      for (int j = dface_group_idx[i]; j < dface_group_idx[i+1]; j++) {
        int ifac = (int) (dface_group[j] - 1);
        //assert (ifac >= 0);
        if (ifac >= 0) {
          face_data[ifac] = (PDM_g_num_t) i;
        }
      }
    }

    sprintf(filename, "dface_group_%3.3d.vtk", i_rank);
    _write_polydata_gnum (filename,
                          "dface",
                          dn_vtx,
                          dvtx_coord,
                          NULL,
                          dn_face,
                          dface_vtx_idx,
                          dface_vtx,
                          face_data);
    free (face_data);
  }


  if (1) {
    PDM_g_num_t n_octo_z_cst = nx * ny;
    PDM_g_num_t n_quadH_z_cst = (nx - 1) * (ny - 1);
    PDM_g_num_t n_tria_z_cst = 2*(nx - 1) + 2*(ny - 1) + 4;//2 * (nx + ny);
    PDM_g_num_t n_quadV_z_cst = nx*(ny+1) + ny*(nx+1) + 4*nx*ny + 2*(nx-1) + 2*(ny-1) + 8;
    PDM_g_num_t n_faceH_z_cst = n_octo_z_cst + n_quadH_z_cst + n_tria_z_cst;
    PDM_g_num_t n_face_z_cst = n_faceH_z_cst + n_quadV_z_cst;
    PDM_g_num_t *distrib_face = PDM_compute_entity_distribution (comm, dn_face);

    //printf("dface_cell:\n");
    for (int i = 0; i < dn_face; i++) {
      PDM_g_num_t g = distrib_face[i_rank] + i;
      PDM_g_num_t k = g / n_face_z_cst;
      PDM_g_num_t r = g % n_face_z_cst;
      /*printf("k = "PDM_FMT_G_NUM", r = "PDM_FMT_G_NUM", "PDM_FMT_G_NUM" : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",
        k, r, g + 1, dface_cell[2*i], dface_cell[2*i+1]);*/
      if (dface_cell[2*i] > ng_cell || dface_cell[2*i+1] > ng_cell) {
        printf("dface_cell outside bounds (face %ld, k = %ld, r = %ld)\n", g+1, k, r);
      }
    }
    PDM_MPI_Barrier (comm);
  }

  if (1) {
    //if (i_rank == 1) printf("dface_group_idx[4,5,6] = %d %d %d\n", dface_group_idx[4], dface_group_idx[5], dface_group_idx[6]);
    for (int i = 0; i < n_face_group; i++) {
      for (int j = dface_group_idx[i]; j < dface_group_idx[i+1]; j++) {
        if (dface_group[j] < 0 || dface_group[j] > ng_face) {
          printf("[%d] group %d, face %d = %ld\n", i_rank, i, j, dface_group[j]);
        }
      }
    }
    PDM_MPI_Barrier (comm);
  }

  #if 1
  /*
   *  Create mesh partitions
   */
  /* Initialize multipart */
  int *n_part_zones = (int *) malloc(sizeof(int));
  n_part_zones[0] = n_part;

  int mpart_id = PDM_multipart_create(1, n_part_zones, PDM_FALSE,
                                      method, PDM_PART_SIZE_HOMOGENEOUS,
                                      NULL, comm, PDM_OWNERSHIP_KEEP);

  /* Generate mesh */
  int n_join = 0;
  PDM_dmesh_t *dmesh = dmesh = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                                dn_cell,
                                                dn_face,
                                                -1, // dn_edge
                                                dn_vtx,
                                                n_face_group,
                                                n_join,
                                                comm);

  int *djoins_ids = malloc (sizeof(int) * n_join);
  int *dface_join_idx = malloc (sizeof(int) * (n_join + 1));
  dface_join_idx[0] = 0;
  PDM_g_num_t *dface_join = malloc (sizeof(PDM_g_num_t) * dface_join_idx[n_join]);

  PDM_dmesh_set (dmesh,
                 dvtx_coord,
                 dface_vtx_idx,
                 dface_vtx,
                 dface_cell,
                 dface_group_idx,
                 dface_group,
                 djoins_ids,
                 dface_join_idx,
                 dface_join);

  PDM_multipart_register_block (mpart_id, 0, dmesh);

  /* Connection between zones */
  int n_total_joins = 0;
  int *join_to_opposite = malloc(sizeof(int) * n_total_joins);
  PDM_multipart_register_joins (mpart_id, n_total_joins, join_to_opposite);

  /* Run */
  PDM_multipart_run_ppart (mpart_id);


  if (post) {
    /* Write faces and dual edges */
    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_proc, tn_part;
      int n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int* n_elt;

      PDM_multipart_part_dim_get(mpart_id, 0, i_part, &n_section, &n_elt,
                                 &n_cell, &n_face, &n_part_joins, &n_vtx, &n_proc, &tn_part,
                                 &scell_face, &sface_vtx, &sface_bound, &n_bounds, &sface_join, &n_joins);

      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, 0, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

      char filename[999];
      sprintf(filename, "faces_%3.3d.vtk", i_rank*n_part + i_part);
      _write_polydata (filename,
                       "faces",
                       n_vtx,
                       vtx,
                       vtx_ln_to_gn,
                       n_face,
                       face_vtx_idx,
                       face_vtx,
                       face_ln_to_gn);


      PDM_g_num_t n_octo_z_cst = nx * ny;
      PDM_g_num_t n_quadH_z_cst = (nx - 1) * (ny - 1);
      PDM_g_num_t n_tria_z_cst = 2*(nx - 1) + 2*(ny - 1) + 4;
      PDM_g_num_t n_quadV_z_cst = nx*(ny+1) + ny*(nx+1) + 4*nx*ny + 2*(nx-1) + 2*(ny-1) + 8;
      PDM_g_num_t n_faceH_z_cst = n_octo_z_cst + n_quadH_z_cst + n_tria_z_cst;
      PDM_g_num_t n_face_z_cst = n_faceH_z_cst + n_quadV_z_cst;

      PDM_g_num_t idx_quadH = n_octo_z_cst;
      PDM_g_num_t idx_tria1 = idx_quadH + n_quadH_z_cst; // bottom row
      PDM_g_num_t idx_tria2 = idx_tria1 + nx - 1;        // top row
      PDM_g_num_t idx_tria3 = idx_tria2 + nx - 1;        // left column
      PDM_g_num_t idx_tria4 = idx_tria3 + ny - 1;        // right column
      PDM_g_num_t idx_tria5 = idx_tria4 + ny - 1;        // corners

      double stepx = lengthx / (double) (3*nx);
      double stepy = lengthy / (double) (3*ny);
      double stepz = 0.;
      if (nz > 0) {
        stepz = lengthz / (double) nz;
      }

      double *cell_center = malloc (sizeof(double) * n_cell * 3);
      for (int icell = 0; icell < n_cell; icell++) {
        PDM_g_num_t g = cell_ln_to_gn[icell] - 1;
        PDM_g_num_t k = g / n_faceH_z_cst;
        PDM_g_num_t r = g % n_faceH_z_cst;
        PDM_g_num_t i, j;

        cell_center[3*icell + 2] = zmin + (k + 0.5) * stepz;
        if (r < idx_quadH) {
          // Octagon
          j = r / nx;
          i = r % nx;

          cell_center[3*icell]     = xmin + (3*i + 1.5) * stepx;
          cell_center[3*icell + 1] = ymin + (3*j + 1.5) * stepy;
        }

        else if (r < idx_tria1) {
          // Quadrangle
          PDM_g_num_t s = r - idx_quadH;

          j = s / (nx - 1);
          i = s % (nx - 1);

          cell_center[3*icell]     = xmin + 3*(i + 1) * stepx;
          cell_center[3*icell + 1] = ymin + 3*(j + 1) * stepy;
        }

        else if (r < idx_tria2) {
          // Triangle (bottom row)
          PDM_g_num_t s = r - idx_tria1;

          i = s % (nx - 1);

          cell_center[3*icell]     = xmin + 3*(i + 1) * stepx;
          cell_center[3*icell + 1] = ymin + stepy/3.;
        }

        else if (r < idx_tria3) {
          // Triangle (top row)
          r -= idx_tria2;

          i = r % (nx - 1);

          cell_center[3*icell]     = xmin + 3*(i + 1) * stepx;
          cell_center[3*icell + 1] = ymin + lengthy - stepy/3.;
        }

        else if (r < idx_tria4) {
          // Triangle (left column)
          PDM_g_num_t s = r - idx_tria3;

          j = s % (ny - 1);

          cell_center[3*icell]     = xmin + stepx/3.;
          cell_center[3*icell + 1] = ymin + 3*(j + 1) * stepy;
        }

        else if (r < idx_tria5) {
          // Triangle (right column)
          r -= idx_tria4;

          j = r % (ny - 1);

          cell_center[3*icell]     = xmin + lengthx - stepx/3.;
          cell_center[3*icell + 1] = ymin + 3*(j + 1) * stepy;
        }

        else if (r < idx_tria5 + 1) {
          // Triangle (bottom-left corner)
          cell_center[3*icell]     = xmin + stepx/3.;
          cell_center[3*icell + 1] = ymin + stepy/3.;
        }

        else if (r < idx_tria5 + 2) {
          // Triangle (bottom-right corner)
          cell_center[3*icell]     = xmin + lengthx - stepx/3.;
          cell_center[3*icell + 1] = ymin + stepy/3.;
        }

        else if (r < idx_tria5 + 3) {
          // Triangle (top-left corner)
          cell_center[3*icell]     = xmin + stepx/3.;
          cell_center[3*icell + 1] = ymin + lengthy - stepy/3.;
        }

        else if (r < idx_tria5 + 4) {
          // Triangle (top-right corner)
          cell_center[3*icell]     = xmin + lengthx - stepx/3.;
          cell_center[3*icell + 1] = ymin + lengthy - stepy/3.;
        }

        else {
          //assert(0);
          printf("[%d] part %d, cell %d ("PDM_FMT_G_NUM"), k = %ld, r = %ld\n", i_rank, i_part, icell, cell_ln_to_gn[icell], k, r);
          for (int l = 0; l < 3; l++) {
            cell_center[3*icell + l] = 0.;
          }
        }
      }

      double *face_center = malloc (sizeof(double) * n_face * 3);
      double *surf_vector = malloc (sizeof(double) * n_face * 3);
      double *char_length = malloc (sizeof(double) * n_face);
      int *is_degenerate = malloc (sizeof(int) * n_face);
      PDM_geom_elem_polygon_properties (n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        vtx,
                                        surf_vector,
                                        face_center,
                                        char_length,
                                        is_degenerate);
      free (surf_vector);
      free (char_length);
      free (is_degenerate);

      sprintf(filename, "dual_edges_%3.3d.vtk", i_rank*n_part + i_part);
      write_dual_edges (filename,
                        n_cell,
                        cell_face_idx,
                        cell_face,
                        cell_center,
                        n_face,
                        face_center,
                        face_ln_to_gn);
      free (cell_center);
      free (face_center);
    }

  }

  if (post) {
    /* Prepare writer */
    int id_cs = PDM_writer_create ("Ensight",
                                   PDM_WRITER_FMT_ASCII,
                                   PDM_WRITER_TOPO_CONSTANTE,
                                   PDM_WRITER_OFF,
                                   "test_polyvol",
                                   "polyvol",
                                   PDM_MPI_COMM_WORLD,
                                   PDM_IO_ACCES_MPI_SIMPLE,
                                   1.,
                                   NULL);

    int geom_id = PDM_writer_geom_create (id_cs,
                                          "mesh",
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_OFF,
                                          n_part);

    PDM_writer_step_beg (id_cs, 0.);

    // Cell local id
    /*int id_var_cell_id = PDM_writer_var_create (id_cs,
                                                PDM_WRITER_OFF,
                                                PDM_WRITER_VAR_SCALAIRE,
                                                PDM_WRITER_VAR_ELEMENTS,
                                                "cell_id");*/

    /* Write geometry */
    int **face_vtx_n  = malloc (sizeof(int *) * n_part);
    int **cell_face_n = malloc (sizeof(int *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      int n_proc, tn_part;
      int n_cell, n_face, n_vtx, n_bounds, n_joins, n_part_joins;
      int scell_face, sface_vtx, sface_bound, sface_join;
      int  n_section;
      int* n_elt;

      PDM_multipart_part_dim_get(mpart_id, 0, i_part, &n_section, &n_elt,
                                 &n_cell, &n_face, &n_part_joins, &n_vtx, &n_proc, &tn_part,
                                 &scell_face, &sface_vtx, &sface_bound, &n_bounds, &sface_join, &n_joins);

      double       *vtx;
      int          *cell_face_idx, *cell_face, *face_cell, *face_vtx_idx, *face_vtx;
      int          *face_bound_idx, *face_bound, *face_join_idx, *face_join;
      int          *face_part_bound_proc_idx, *face_part_bound_part_idx, *face_part_bound;
      PDM_g_num_t  *cell_ln_to_gn, *face_ln_to_gn, *vtx_ln_to_gn, *face_bound_ln_to_gn, *face_join_ln_to_gn;
      int          *cell_tag, *face_tag, *vtx_tag;
      int         **elt_vtx_idx;
      int         **elt_vtx;
      PDM_g_num_t **elt_section_ln_to_gn;

      PDM_multipart_part_val_get(mpart_id, 0, i_part, &elt_vtx_idx, &elt_vtx, &elt_section_ln_to_gn,
                                 &cell_tag, &cell_face_idx, &cell_face, &cell_ln_to_gn,
                                 &face_tag, &face_cell, &face_vtx_idx, &face_vtx, &face_ln_to_gn,
                                 &face_part_bound_proc_idx, &face_part_bound_part_idx, &face_part_bound,
                                 &vtx_tag, &vtx, &vtx_ln_to_gn, &face_bound_idx, &face_bound,
                                 &face_bound_ln_to_gn, &face_join_idx, &face_join, &face_join_ln_to_gn);

      PDM_writer_geom_coord_set (id_cs,
                                 geom_id,
                                 i_part,
                                 n_vtx,
                                 vtx,
                                 vtx_ln_to_gn);

      face_vtx_n[i_part] = malloc (sizeof(int) * n_face);
      for (int i = 0; i < n_face; i++) {
        face_vtx_n[i_part][i] = face_vtx_idx[i+1] - face_vtx_idx[i];
      }

      cell_face_n[i_part] = malloc (sizeof(int) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        cell_face_n[i_part][i] = cell_face_idx[i+1] - cell_face_idx[i];
      }

      PDM_writer_geom_cell3d_cellface_add (id_cs,
                                          geom_id,
                                          i_part,
                                          n_cell,
                                          n_face,
                                          face_vtx_idx,
                                          face_vtx_n[i_part],
                                          face_vtx,
                                          cell_face_idx,
                                          cell_face_n[i_part],
                                          cell_face,
                                          cell_ln_to_gn);

    }

    PDM_writer_geom_write (id_cs, geom_id);

    // write variables...

    PDM_writer_step_end (id_cs);

    PDM_writer_geom_data_free (id_cs, geom_id);
    PDM_writer_geom_free (id_cs, geom_id);
    PDM_writer_free (id_cs);
  }

  /*
   *  Finalize
   */
  free (dface_cell);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dvtx_coord);
  free (dface_group_idx);
  free (dface_group);
  free (dface_join_idx);
  free (dface_join);
  free (djoins_ids);
  free (join_to_opposite);

  PDM_multipart_free (mpart_id);
  free (dmesh);
  free (n_part_zones);

#endif

  PDM_MPI_Finalize();

  return 0;
}
