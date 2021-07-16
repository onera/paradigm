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
#include "pdm_part.h"
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


static void
_cell_center_from_g_num
(
 const int           nx,
 const int           ny,
 const int           nz,
 const double        xmin,
 const double        ymin,
 const double        zmin,
 const double        lengthx,
 const double        lengthy,
 const double        lengthz,
 const int           n_cell,
 const PDM_g_num_t  *cell_g_num,
 double            **cell_center
 )
{
  PDM_g_num_t n_octo_z_cst = nx * ny;
  PDM_g_num_t n_quadH_z_cst = (nx - 1) * (ny - 1);
  PDM_g_num_t n_tria_z_cst = 2*(nx - 1) + 2*(ny - 1) + 4;
  PDM_g_num_t n_faceH_z_cst = n_octo_z_cst + n_quadH_z_cst + n_tria_z_cst;

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

  *cell_center = malloc (sizeof(double) * n_cell * 3);

  for (int icell = 0; icell < n_cell; icell++) {
    PDM_g_num_t g = cell_g_num[icell] - 1;
    double *_cell_center = *cell_center + 3*icell;
    PDM_g_num_t k = g / n_faceH_z_cst;
    PDM_g_num_t r = g % n_faceH_z_cst;
    PDM_g_num_t i, j;

    _cell_center[2] = zmin + (k + 0.5) * stepz;
    if (r < idx_quadH) {
      // Octagon
      j = r / nx;
      i = r % nx;

      _cell_center[0] = xmin + (3*i + 1.5) * stepx;
      _cell_center[1] = ymin + (3*j + 1.5) * stepy;
    }

    else if (r < idx_tria1) {
      // Quadrangle
      PDM_g_num_t s = r - idx_quadH;

      j = s / (nx - 1);
      i = s % (nx - 1);

      _cell_center[0] = xmin + 3*(i + 1) * stepx;
      _cell_center[1] = ymin + 3*(j + 1) * stepy;
    }

    else if (r < idx_tria2) {
      // Triangle (bottom row)
      PDM_g_num_t s = r - idx_tria1;

      i = s % (nx - 1);

      _cell_center[0] = xmin + 3*(i + 1) * stepx;
      _cell_center[1] = ymin + stepy/3.;
    }

    else if (r < idx_tria3) {
      // Triangle (top row)
      r -= idx_tria2;

      i = r % (nx - 1);

      _cell_center[0] = xmin + 3*(i + 1) * stepx;
      _cell_center[1] = ymin + lengthy - stepy/3.;
    }

    else if (r < idx_tria4) {
      // Triangle (left column)
      PDM_g_num_t s = r - idx_tria3;

      j = s % (ny - 1);

      _cell_center[0] = xmin + stepx/3.;
      _cell_center[1] = ymin + 3*(j + 1) * stepy;
    }

    else if (r < idx_tria5) {
      // Triangle (right column)
      r -= idx_tria4;

      j = r % (ny - 1);

      _cell_center[0] = xmin + lengthx - stepx/3.;
      _cell_center[1] = ymin + 3*(j + 1) * stepy;
    }

    else if (r < idx_tria5 + 1) {
      // Triangle (bottom-left corner)
      _cell_center[0] = xmin + stepx/3.;
      _cell_center[1] = ymin + stepy/3.;
    }

    else if (r < idx_tria5 + 2) {
      // Triangle (bottom-right corner)
      _cell_center[0] = xmin + lengthx - stepx/3.;
      _cell_center[1] = ymin + stepy/3.;
    }

    else if (r < idx_tria5 + 3) {
      // Triangle (top-left corner)
      _cell_center[0] = xmin + stepx/3.;
      _cell_center[1] = ymin + lengthy - stepy/3.;
    }

    else if (r < idx_tria5 + 4) {
      // Triangle (top-right corner)
      _cell_center[0] = xmin + lengthx - stepx/3.;
      _cell_center[1] = ymin + lengthy - stepy/3.;
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
  int         *dcell_face_idx  = NULL;
  PDM_g_num_t *dcell_face      = NULL;
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
                    &dcell_face_idx,
                    &dcell_face,
                    &dface_cell,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dvtx_coord,
                    &dface_group_idx,
                    &dface_group);

  if (i_rank == 0) printf("ng_cell = "PDM_FMT_G_NUM", ng_face = "PDM_FMT_G_NUM", ng_vtx = "PDM_FMT_G_NUM"\n", ng_cell, ng_face, ng_vtx);

  /*
   *  Create mesh partitions
   */
  int have_dcell_part = 0;

  int *dcell_part = (int *) malloc(dn_cell*sizeof(int));
  int *renum_properties_cell = NULL;
  int *renum_properties_face = NULL;
  int n_property_cell = 0;
  int n_property_face = 0;

  int ppart_id = 0;

  PDM_part_create(&ppart_id,
                  PDM_MPI_COMM_WORLD,
                  method,
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



  if (post) {
    /* Write faces and dual edges */
    for (int i_part = 0; i_part < n_part; i_part++) {

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

      PDM_part_part_dim_get(ppart_id,
                            i_part,
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

      PDM_part_part_val_get (ppart_id,
                             i_part,
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


      double *cell_center = NULL;
      _cell_center_from_g_num (nx, ny, nz,
                               xmin, ymin, zmin,
                               lengthx, lengthy, lengthz,
                               n_cell,
                               cell_ln_to_gn,
                               &cell_center);

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

#if 1
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

    int id_geom = PDM_writer_geom_create (id_cs,
                                          "mesh",
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_OFF,
                                          n_part);

    PDM_writer_step_beg (id_cs, 0.);

    // Cell local id
    int id_var_cell_g_num = PDM_writer_var_create (id_cs,
                                                   PDM_WRITER_OFF,
                                                   PDM_WRITER_VAR_SCALAIRE,
                                                   PDM_WRITER_VAR_ELEMENTS,
                                                   "cell_g_num");

    int id_var_num_part = PDM_writer_var_create (id_cs,
                                                 PDM_WRITER_OFF,
                                                 PDM_WRITER_VAR_SCALAIRE,
                                                 PDM_WRITER_VAR_ELEMENTS,
                                                 "num_part");

    /* Write geometry */
    int **face_vtx_n  = malloc (sizeof(int *) * n_part);
    int **cell_face_n = malloc (sizeof(int *) * n_part);

    PDM_real_t **val_cell_g_num = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);
    PDM_real_t **val_num_part = (PDM_real_t **) malloc(sizeof(PDM_real_t *) * n_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
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

      PDM_part_part_dim_get(ppart_id,
                            i_part,
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

      PDM_part_part_val_get (ppart_id,
                             i_part,
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


      PDM_writer_geom_coord_set (id_cs,
                                 id_geom,
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
                                           id_geom,
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

      val_cell_g_num[i_part]  = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      val_num_part[i_part] = (PDM_real_t *) malloc(sizeof(PDM_real_t) * n_cell);
      for (int i = 0; i < n_cell; i++) {
        val_cell_g_num[i_part][i] = (PDM_real_t) cell_ln_to_gn[i];
        val_num_part[i_part][i] = (PDM_real_t) (i_rank*n_part + i_part);
      }
    }

    PDM_writer_geom_write (id_cs,
                           id_geom);

    // write variables
    for (int i_part = 0; i_part < n_part; i_part++) {

      PDM_writer_var_set (id_cs,
                          id_var_cell_g_num,
                          id_geom,
                          i_part,
                          val_cell_g_num[i_part]);
      PDM_writer_var_set (id_cs,
                          id_var_num_part,
                          id_geom,
                          i_part,
                          val_num_part[i_part]);
    }

    PDM_writer_var_write (id_cs,
                          id_var_cell_g_num);
    PDM_writer_var_write (id_cs,
                          id_var_num_part);

    PDM_writer_var_free (id_cs,
                         id_var_cell_g_num);
    PDM_writer_var_free (id_cs,
                         id_var_num_part);

    for (int i_part = 0; i_part < n_part; i_part++) {
      free (val_cell_g_num[i_part]);
      free (val_num_part[i_part]);
    }
    free (val_cell_g_num);
    free (val_num_part);


    PDM_writer_step_end (id_cs);

    PDM_writer_geom_data_free (id_cs, id_geom);
    PDM_writer_geom_free (id_cs, id_geom);
    PDM_writer_free (id_cs);
  }
#endif

  /*
   *  Finalize
   */
  free (dface_cell);
  free (dface_vtx_idx);
  free (dface_vtx);
  free (dvtx_coord);
  free (dface_group_idx);
  free (dface_group);

  PDM_MPI_Finalize();

  return 0;
}
