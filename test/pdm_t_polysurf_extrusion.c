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
#include "pdm_poly_vol_gen.h"

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
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
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

  int i_rank, n_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &n_rank);

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
  int         *dface_cell_idx  = NULL;
  PDM_g_num_t *dface_cell      = NULL;
  int         *dface_vtx_idx   = NULL;
  PDM_g_num_t *dface_vtx       = NULL;
  double      *dvtx_coord      = NULL;
  int         *dface_group_idx = NULL;
  PDM_g_num_t *dface_group     = NULL;

  PDM_poly_vol_gen (PDM_MPI_COMM_WORLD,
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
                    &dface_cell_idx,
                    &dface_cell,
                    &dface_vtx_idx,
                    &dface_vtx,
                    &dvtx_coord,
                    &dface_group_idx,
                    &dface_group);

  printf("[%d] dn_vtx = %d\n", i_rank, dn_vtx);

  if (1) {
    char filename[999];
    sprintf(filename, "dvtx_coord_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "dvtx",
                        dn_vtx,
                        dvtx_coord,
                        NULL,
                        NULL);

    sprintf(filename, "dface_vtx_%3.3d.vtk", i_rank);
    _write_polydata_gnum (filename,
                          "dface",
                          dn_vtx,
                          dvtx_coord,
                          NULL,
                          dn_face,
                          dface_vtx_idx,
                          dface_vtx,
                          NULL);
  }


  /*
   *  Create mesh partitions
   */

  //...



  /*
   *  Finalize
   */
  PDM_MPI_Finalize();

  return 0;
}
