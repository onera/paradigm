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
#include "pdm_mesh_nodal.h"

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
 * \param [inout]   filename   VTK mesh file
 * \param [inout]   tolerance  Bounding boxes tolerance
 * \param [inout]   n_part     Number of partitions par process
 * \param [inout]   post       Ensight outputs status
 * \param [inout]   method     Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int                          argc,
           char                       **argv,
           char                       **filename,
           double                      *tolerance,
           int                         *n_part,
           PDM_g_num_t                 *n_pts,
           int                         *post,
           int                         *part_method,
           PDM_mesh_location_method_t  *loc_method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *filename = argv[i];

    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
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
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}




static void
_random_cloud
(
 const int      n_pts,
 const double   xyz_min[3],
 const double   xyz_max[3],
 const int      n_procs,
 const int      my_rank,
 double       **coord,
 int           *n_pts_l
 )
{
  double length[3] = {xyz_max[0] - xyz_min[0],
                      xyz_max[1] - xyz_min[1],
                      xyz_max[2] - xyz_min[2]};


  *n_pts_l = (int) (n_pts/n_procs);
  *coord = malloc (sizeof(double) * 3 * (*n_pts_l));
  double *_coord = *coord;
  double x;
  int idx = 0;
  for (PDM_g_num_t i = 0; i < n_procs*(*n_pts_l); i++) {
    for (int idim = 0; idim < 3; idim++) {
      x = xyz_min[idim] + length[idim] * (double) rand() / ((double) RAND_MAX);
      if (i%n_procs == my_rank) {
        _coord[idx++] = x;
      }
    }
  }
}





static PDM_Mesh_nodal_elt_t _pdm_elt_type (const int vtk_type) {

  switch (vtk_type) {
  case 1:
    return PDM_MESH_NODAL_POINT;

  case 3:
    return PDM_MESH_NODAL_BAR2;

  case 5:
    return PDM_MESH_NODAL_TRIA3;

  case 7:
    return PDM_MESH_NODAL_POLY_2D;

  case 9:
    return PDM_MESH_NODAL_QUAD4;

  case 10:
    return PDM_MESH_NODAL_TETRA4;

  case 12:
    return PDM_MESH_NODAL_HEXA8;

  case 13:
    return PDM_MESH_NODAL_PRISM6;

  case 14:
    return PDM_MESH_NODAL_PYRAMID5;

  default:
    PDM_error (__FILE__, __LINE__, 0, "No match with VTK element type");
    return PDM_MESH_NODAL_POLY_3D;
  }

}


static void
_read_vtk_unstructured_grid_seq
(
 const char            *filename,
 int                   *n_vtx,
 PDM_g_num_t          **vtx_g_num,
 double               **vtx_coords,
 int                   *n_cells,
 PDM_g_num_t          **cell_g_num,
 int                  **cell_vtx_idx,
 int                  **cell_vtx,
 PDM_Mesh_nodal_elt_t **cell_type
 )
{
  char line[100], elt_line[100];

  double *_vtx_coords = NULL;

  int n_cell_vtx = 0;
  int *_cell_vtx_idx = NULL;
  int *_cell_vtx = NULL;
  PDM_Mesh_nodal_elt_t *_cell_type = NULL;

  int *id_cell = NULL;

  int elt_vtx[8];

  FILE *f = fopen(filename, "r");

  if (f == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Unable to open VTK file");
  }

  int count = 0;

  while (fgets(line, sizeof(line), f) != NULL && count < 1000000000) {

    if(strcmp(line,"\n") || strcmp(line,"\r\n")) {

      /* Points */
      if (strstr(line, "POINTS") != NULL) {
        int stat = sscanf(line,
                          "%*[^0123456789]%d%*[^0123456789]",
                          n_vtx);
        assert (stat);

        *vtx_g_num = malloc (sizeof(PDM_g_num_t) * (*n_vtx));
        *vtx_coords = malloc (sizeof(double) * (*n_vtx) * 3);
        _vtx_coords = *vtx_coords;

        for (int ivtx = 0; ivtx < *n_vtx; ivtx++) {
          fscanf(f, "%lf %lf %lf",
                 _vtx_coords + 3*ivtx,
                 _vtx_coords + 3*ivtx + 1,
                 _vtx_coords + 3*ivtx + 2);

          (*vtx_g_num)[ivtx] = 1 + ivtx;
        }
      }


      /* Cells */
      else if (strstr(line, "CELLS") != NULL) {
        int n_elts = 0;
        int stat = sscanf(line,
                          "%*[^0123456789]%d%*[^0123456789]%d",
                          &n_elts,
                          &n_cell_vtx);
        assert (stat);

        n_cell_vtx -= n_elts;

        *cell_g_num = malloc (sizeof(PDM_g_num_t) * n_elts);

        *cell_vtx_idx = malloc (sizeof(int) * (n_elts + 1));
        _cell_vtx_idx = *cell_vtx_idx;
        _cell_vtx_idx[0] = 0;

        *cell_vtx = malloc (sizeof(int) * n_cell_vtx);
        _cell_vtx = *cell_vtx;

        id_cell = malloc (sizeof(int) * n_elts);

        int icell = 0;
        for (int ielt = 0; ielt < n_elts; ielt++) {
          char *ptr = fgets(elt_line, sizeof(elt_line), f);
          assert (ptr != NULL);

          int n_vtx_elt;
          int stat2 = sscanf(elt_line,
                             "%d%*[^0123456789]",
                             &n_vtx_elt);
          assert (stat2);


          int unused;
          switch (n_vtx_elt) {
          case 1:
            stat2 = sscanf(elt_line, "%d %d", &unused,
                           &(elt_vtx[0]));
            break;

          case 2:
            stat2 = sscanf(elt_line, "%d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]));
            break;

          case 3:
            stat2 = sscanf(elt_line, "%d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]),
                           &(elt_vtx[2]));
            break;

          case 4:
            stat2 = sscanf(elt_line, "%d %d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]),
                           &(elt_vtx[2]),
                           &(elt_vtx[3]));
            break;

          case 5:
            stat2 = sscanf(elt_line, "%d %d %d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]),
                           &(elt_vtx[2]),
                           &(elt_vtx[3]),
                           &(elt_vtx[4]));
            break;

          case 6:
            /*stat2 = sscanf(elt_line, "%d %d %d %d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]),
                           &(elt_vtx[2]),
                           &(elt_vtx[3]),
                           &(elt_vtx[4]),
                           &(elt_vtx[5]));*/
            //for prisms, VTK and PDM vertex orderings are different
            stat2 = sscanf(elt_line, "%d %d %d %d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[2]),
                           &(elt_vtx[1]),
                           &(elt_vtx[3]),
                           &(elt_vtx[5]),
                           &(elt_vtx[4]));
            break;

          case 8:
            stat2 = sscanf(elt_line, "%d %d %d %d %d %d %d %d %d", &unused,
                           &(elt_vtx[0]),
                           &(elt_vtx[1]),
                           &(elt_vtx[2]),
                           &(elt_vtx[3]),
                           &(elt_vtx[4]),
                           &(elt_vtx[5]),
                           &(elt_vtx[6]),
                           &(elt_vtx[7]));
            break;

          default:
            printf("Unknown element with %d vertices (ignored)\n", n_vtx_elt);
            id_cell[ielt] = -1;
            continue;
          }

          assert (stat2);

          for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
            _cell_vtx[_cell_vtx_idx[icell] + ivtx] = elt_vtx[ivtx];
          }
          _cell_vtx_idx[icell+1] = _cell_vtx_idx[icell] + n_vtx_elt;

          id_cell[ielt] = icell;
          (*cell_g_num)[icell] = icell + 1;
          icell++;
        }

        *n_cells = icell;
        if (*n_cells < n_elts) {
          *cell_g_num = realloc (*cell_g_num, sizeof(PDM_g_num_t) * (*n_cells));
          *cell_vtx = realloc (*cell_vtx, sizeof(int) * _cell_vtx_idx[*n_cells]);
          *cell_vtx_idx = realloc (*cell_vtx_idx, sizeof(int) * (*n_cells + 1));
        }
      }

      /* Cell types */
      else if (strstr(line, "CELL_TYPES") != NULL) {
        int n_elts;
        int stat = sscanf(line,
                          "%*[^0123456789]%d",
                          &n_elts);
        assert (stat);

        *cell_type = malloc (sizeof(PDM_Mesh_nodal_elt_t) * (*n_cells));
        _cell_type = *cell_type;

        int icell = 0;
        for (int ielt = 0; ielt < *n_cells; ielt++) {
          int vtk_type;
          fscanf(f, "%d", &vtk_type);

          if (id_cell[ielt] < 0) {
            continue;
          }

          _cell_type[icell] =  _pdm_elt_type (vtk_type);
          icell++;
        }
      }

    }


    else {
      continue;
    }

    count++;
  }

  fclose(f);
}





static int
_mesh_nodal_from_VTK_unstructured_grid
(
 const char* filename
 )
{
  /*
   *  Read VTK file
   */
  int n_vtx;
  PDM_g_num_t *vtx_g_num = NULL;
  double *vtx_coords = NULL;

  int n_cell;
  PDM_g_num_t *cell_g_num = NULL;
  int *cell_vtx_idx = NULL;
  int *cell_vtx = NULL;
  PDM_Mesh_nodal_elt_t *cell_type = NULL;

  _read_vtk_unstructured_grid_seq (filename,
                                   &n_vtx,
                                   &vtx_g_num,
                                   &vtx_coords,
                                   &n_cell,
                                   &cell_g_num,
                                   &cell_vtx_idx,
                                   &cell_vtx,
                                   &cell_type);


  /*
   *  Build mesh nodal
   */
  const int n_parts = 1;
  /* Create mesh nodal structure */
  int mesh_nodal_id = PDM_Mesh_nodal_create (n_parts,
                                             PDM_MPI_COMM_WORLD);

  /* Set mesh nodal vertices */
  PDM_Mesh_nodal_coord_set (mesh_nodal_id,
                            0,//id_part,
                            n_vtx,
                            vtx_coords,
                            vtx_g_num);


  /* Set mesh nodal blocks */
  /*   count nb of elements for each type */
  const int n_types = 10;
  int n_elt[n_types];
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT; type <= PDM_MESH_NODAL_POLY_3D; type++) {
    n_elt[type] = 0;
  }
  for (int icell = 0; icell < n_cell; icell++) {
    PDM_Mesh_nodal_elt_t type = cell_type[icell];

    //-->>
    if (type == PDM_MESH_NODAL_POLY_2D ||
        type == PDM_MESH_NODAL_POLY_3D) {
      //not supported yet...
      continue;
    }
    //<<--

    n_elt[type]++;
  }


  /*   prepare blocks */
  const PDM_bool_t st_free_data = PDM_TRUE;
  const int order = 1;
  PDM_l_num_t **block_connec = malloc (sizeof(PDM_l_num_t *) * n_types);
  PDM_g_num_t **block_g_num  = malloc (sizeof(PDM_g_num_t *) * n_types);
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT; type <= PDM_MESH_NODAL_POLY_3D; type++) {

    if (n_elt[type] < 1) {
      continue;
    }

    int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (type, order);

    block_connec[type] = malloc (sizeof(PDM_l_num_t) * n_elt[type] * n_vtx_elt);
    block_g_num[type]  = malloc (sizeof(PDM_g_num_t) * n_elt[type]);
  }


  /*   fill blocks */
  int idx[n_types];
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT; type <= PDM_MESH_NODAL_POLY_3D; type++) {
    idx[type] = 0;
  }

  for (int icell = 0; icell < n_cell; icell++) {
    PDM_Mesh_nodal_elt_t type = cell_type[icell];

    if (n_elt[type] < 1) {
      continue;
    }

    int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (type, order);
    for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
      block_connec[type][n_vtx_elt*idx[type] + ivtx] =
        (PDM_l_num_t) (cell_vtx[cell_vtx_idx[icell] + ivtx] + 1);
    }

    block_g_num[type][idx[type]] = cell_g_num[icell];
    idx[type]++;
  }

  /*   add non-empty blocks */
  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT; type <= PDM_MESH_NODAL_POLY_3D; type++) {

    if (n_elt[type] < 1) {
      continue;
    }

    int id_block = PDM_Mesh_nodal_block_add (mesh_nodal_id,
                                             st_free_data,
                                             type);

    PDM_Mesh_nodal_block_std_set (mesh_nodal_id,
                                  id_block,
                                  0,//id_part,
                                  n_elt[type],
                                  block_connec[type],
                                  block_g_num[type],
                                  NULL);
  }


  return mesh_nodal_id;
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
  char   *filename  = NULL;
  double  tolerance = 1e-6;
  int     n_part    = 1;
  int     post      = 0;
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

  _read_args(argc,
             argv,
             &filename,
             &tolerance,
             &n_part,
             &n_pts,
             &post,
             (int *) &part_method,
             &loc_method);



  if (filename == NULL) {
    filename = malloc (sizeof(char) * 100);
    sprintf(filename, "cube_tetra_pyram.vtk");
  }

  /*
   *  Init
   */

  int my_rank;
  int n_procs;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &my_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_procs);


  /************************
   *
   * Create mesh nodal structure from VTK file
   *
   ************************/

  /* Read VTK file */
  int n_vtx;
  PDM_g_num_t *vtx_g_num = NULL;
  double *vtx_coords = NULL;

  int n_cells;
  PDM_g_num_t *cell_g_num = NULL;
  int *cell_vtx_idx = NULL;
  int *cell_vtx = NULL;
  PDM_Mesh_nodal_elt_t *cell_type = NULL;

  if (my_rank == 0) {
    printf("-- Read %s\n", filename);
    fflush(stdout);
  }
  _read_vtk_unstructured_grid_seq (filename,
                                   &n_vtx,
                                   &vtx_g_num,
                                   &vtx_coords,
                                   &n_cells,
                                   &cell_g_num,
                                   &cell_vtx_idx,
                                   &cell_vtx,
                                   &cell_type);

  double _xyz_min[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  double _xyz_max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int ivtx = 0; ivtx < n_vtx; ivtx++) {
    for (int idim = 0; idim < 3; idim++) {
      _xyz_min[idim] = PDM_MIN (_xyz_min[idim], vtx_coords[3*ivtx + idim]);
      _xyz_max[idim] = PDM_MAX (_xyz_max[idim], vtx_coords[3*ivtx + idim]);
    }
  }

  double xyz_min[3];
  double xyz_max[3];
  PDM_MPI_Allreduce (&_xyz_min, &xyz_min, 3, PDM_MPI_DOUBLE, PDM_MPI_MIN, PDM_MPI_COMM_WORLD);
  PDM_MPI_Allreduce (&_xyz_max, &xyz_max, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, PDM_MPI_COMM_WORLD);

  /*if (my_rank == 0) {
    printf("%f < x < %f\n", xyz_min[0], xyz_max[0]);
    printf("%f < y < %f\n", xyz_min[1], xyz_max[1]);
    printf("%f < z < %f\n", xyz_min[2], xyz_max[2]);
    }*/

  if (my_rank == 0) {
    printf("-- Build mesh nodal\n");
    fflush(stdout);
  }

  const int n_parts = 1;
  /* Create mesh nodal structure */
  int mesh_nodal_id = PDM_Mesh_nodal_create (n_parts,
                                             PDM_MPI_COMM_WORLD);


  /* Set mesh nodal vertices */
  PDM_Mesh_nodal_coord_set (mesh_nodal_id,
                            0,//id_part,
                            n_vtx,
                            vtx_coords,
                            vtx_g_num);



  /* Set mesh nodal blocks */
#define NO_NODE 0
#define NO_LINE 0
#define NO_SURF 0
#define NO_POLY 0

  /*   count nb of elements for each type */
  const int n_types = 10;
  int n_elt[n_types];
  for (PDM_Mesh_nodal_elt_t itype = PDM_MESH_NODAL_POINT; itype <= PDM_MESH_NODAL_POLY_3D; itype++) {
    n_elt[itype] = 0;
  }
  for (int ielt = 0; ielt < n_cells; ielt++) {
    n_elt[cell_type[ielt]]++;
  }


  /*   prepare non-empty blocks */
  const int id_offset = 1;
  const int order = 1;
  int         **block_connec = malloc (sizeof(int *)         * n_types);
  PDM_g_num_t **block_g_num  = malloc (sizeof(PDM_g_num_t *) * n_types);

  for (PDM_Mesh_nodal_elt_t itype = PDM_MESH_NODAL_POINT; itype <= PDM_MESH_NODAL_POLY_3D; itype++) {
    if (n_elt[itype] < 1) {
      continue;
    }

#if NO_NODE
    if (itype == PDM_MESH_NODAL_POINT) {
      continue;
    }
#endif

#if NO_LINE
    if (itype == PDM_MESH_NODAL_BAR2) {
      continue;
    }
#endif

#if NO_POLY
    if (itype == PDM_MESH_NODAL_POLY_2D ||
        itype == PDM_MESH_NODAL_POLY_3D) {
      //not supported yet...
      continue;
    }
#endif

#if NO_SURF
    if (itype == PDM_MESH_NODAL_TRIA3 ||
        itype == PDM_MESH_NODAL_QUAD4) {
      continue;
    }
#endif

    int n_vtx_elt = PDM_Mesh_nodal_n_vertices_element (itype, order);

    block_connec[itype] = malloc (sizeof(int) * n_elt[itype] * n_vtx_elt);
    block_g_num[itype] =  malloc (sizeof(PDM_g_num_t) * n_elt[itype]);

    int icell = 0;
    for (int ielt = 0; ielt < n_cells; ielt++) {
      if (cell_type[ielt] == itype) {
        block_g_num[itype][icell] = cell_g_num[ielt];

        assert (cell_vtx_idx[ielt+1] - cell_vtx_idx[ielt] == n_vtx_elt);

        for (int ivtx = 0; ivtx < n_vtx_elt; ivtx++) {
          block_connec[itype][n_vtx_elt*icell + ivtx] =
            cell_vtx[cell_vtx_idx[ielt] + ivtx] + id_offset;
        }

        icell++;
      }
    }
  }



  /*   add non-empty blocks */
  const PDM_bool_t st_free_data = PDM_TRUE;

  for (PDM_Mesh_nodal_elt_t itype = PDM_MESH_NODAL_POINT; itype <= PDM_MESH_NODAL_POLY_3D; itype++) {
    if (n_elt[itype] < 1) {
      continue;
    }

#if NO_NODE
    if (itype == PDM_MESH_NODAL_POINT) {
      continue;
    }
#endif

#if NO_LINE
    if (itype == PDM_MESH_NODAL_BAR2) {
      continue;
    }
#endif

#if NO_POLY
    if (itype == PDM_MESH_NODAL_POLY_2D ||
        itype == PDM_MESH_NODAL_POLY_3D) {
      //not supported yet...
      continue;
    }
#endif

#if NO_SURF
    if (itype == PDM_MESH_NODAL_TRIA3 ||
        itype == PDM_MESH_NODAL_QUAD4) {
      continue;
    }
#endif

    int id_block = PDM_Mesh_nodal_block_add (mesh_nodal_id,
                                             st_free_data,
                                             itype);


    PDM_Mesh_nodal_block_std_set (mesh_nodal_id,
                                  id_block,
                                  0,//id_part,
                                  n_elt[itype],
                                  block_connec[itype],
                                  block_g_num[itype],
                                  NULL);
  }

  if (my_rank == 0) {
    int n_blocks = PDM_Mesh_nodal_n_blocks_get (mesh_nodal_id);
    printf("n_blocks = %d\n", n_blocks);
    fflush(stdout);
  }



  /************************
   *
   * Define point cloud
   *
   ************************/
  if (my_rank == 0) {
    printf("-- Point cloud\n");
    fflush(stdout);
  }

  int n_pts_l;
  double *pts_coords = NULL;
  _random_cloud (n_pts,
                 xyz_min,
                 xyz_max,
                 n_procs,
                 my_rank,
                 &pts_coords,
                 &n_pts_l);

  int id_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *char_length = malloc(sizeof(double) * n_pts_l);
  double length = 0;
  for (int idim = 0; idim < 3; idim++) {
    length = PDM_MAX (length, xyz_max[idim] - xyz_min[idim]);
  }

  for (int i = 0; i < n_pts_l; i++) {
    char_length[i] = length * 1.e-6;
  }

  PDM_gnum_set_from_coords (id_gnum, 0, n_pts_l, pts_coords, char_length);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *pts_gnum = PDM_gnum_get(id_gnum, 0);

  PDM_gnum_free (id_gnum, 1);


#if 0
  for (int ipt = 0; ipt < n_pts_l; ipt++) {
    printf("[%d] (%ld) (%f, %f, %f)\n",
           my_rank, pts_gnum[ipt], pts_coords[3*ipt], pts_coords[3*ipt+1], pts_coords[3*ipt+2]);
  }
#endif


  /************************
   *
   * Create mesh location structure
   *
   ************************/
  if (my_rank == 0) {
    printf("-- Locate\n");
    fflush(stdout);
  }

  int location_id = PDM_mesh_location_create (PDM_MESH_NATURE_NODAL,//?
                                              1,//const int n_point_cloud,
                                              PDM_MPI_COMM_WORLD);

  PDM_mesh_location_shared_nodal_mesh_set (location_id,
                                           mesh_nodal_id);


  /* Set point cloud(s) */
  PDM_mesh_location_n_part_cloud_set (location_id,
                                      0,//i_point_cloud,
                                      1);//n_part

  PDM_mesh_location_cloud_set (location_id,
                               0,//i_point_cloud,
                               0,//i_part,
                               n_pts_l,
                               pts_coords,
                               pts_gnum);

  /* Set location parameters */
  PDM_mesh_location_tolerance_set (location_id,
                                   tolerance);

  PDM_mesh_location_method_set (location_id,
                                loc_method);


  /************************
   *
   * Compute location
   *
   ************************/
  PDM_mesh_location_compute (location_id);

  PDM_mesh_location_dump_times (location_id);

  PDM_g_num_t *location_elt_gnum = NULL;
  PDM_mesh_location_get (location_id,
                         0,//i_point_cloud,
                         0,//i_part,
                         &location_elt_gnum);


  if (my_rank == 0) {
    printf("--- Results ---\n");
  }
  for (int ipt = 0; ipt < n_pts_l; ipt++) {
    printf("Point (%ld) : location = (%ld)\n", pts_gnum[ipt], location_elt_gnum[ipt]);
  }


  /*
   *  Finalize
   */
  PDM_mesh_location_free (location_id,
                          0);

  PDM_MPI_Finalize();

  if (my_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
