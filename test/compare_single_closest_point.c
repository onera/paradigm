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
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_octree.h"
#include "pdm_para_octree.h"
#include "pdm_timer.h"

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
     "  -t      <level>  Number of Target points (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc       Number of arguments
 * \param [in]      argv       Arguments
 * \param [inout]   n_faceSeg  Number of vertices on the cube side
 * \param [inout]   length     Cube length
 * \param [inout]   nTgt       Number of Target points
 * \param [inout]   n_part     Number of partitions par process
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *n_faceSeg,
           double        *length,
           PDM_g_num_t   *nTgt,
           int           *n_max_per_leaf,
           int           *surf_source,
           int           *randomize,
           int           *repeat_last,
           int           *n_methods)
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
        long _n_faceSeg = atol(argv[i]);
        *n_faceSeg = (PDM_g_num_t) _n_faceSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
      }
    }
    else if (strcmp(argv[i], "-mpl") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_max_per_leaf = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-surf") == 0) {
      *surf_source = 1;
    }
    else if (strcmp(argv[i], "-rand") == 0) {
      *randomize = 1;
    }
    else if (strcmp(argv[i], "-repeat") == 0) {
      *repeat_last = 1;
    }
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *n_methods = atoi(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


static inline double _rand01(void) {
  return (double) rand() / (double) RAND_MAX;
}

static void
_gen_cloud_random
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_pts,
 const double        origin[3],
 const double        length,
 int                *_n_pts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  // Define distribution
  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  distrib[0] = 0;
  PDM_g_num_t step = n_pts / n_rank;
  PDM_g_num_t remainder = n_pts % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] = step;
    const int i1 = i - 1;
    if (i1 < remainder) {
      distrib[i]++;
    }
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distrib[i] += distrib[i-1];
  }

  PDM_g_num_t dn_pts = distrib[i_rank+1] - distrib[i_rank];
  *_n_pts = (int) dn_pts;

  *g_num = malloc (sizeof(PDM_g_num_t) * dn_pts);
  *coord = malloc (sizeof(double)      * dn_pts * 3);
  for (int i = 0; i < *_n_pts; i++) {
    (*g_num)[i] = 1 + i + distrib[i_rank];
    for (int j = 0; j < 3; j++) {
      (*coord)[3*i+j] = origin[j] + length * _rand01();
    }
  }
}


static void
_gen_cloud_vol
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_faceSeg,
 const double        origin[3],
 const double        length,
 const int           randomize,
 int                *npts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t n_cell = n_faceSeg * n_faceSeg * n_faceSeg;

  PDM_g_num_t *distribCell = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t n_faceFace = n_faceSeg * n_faceSeg;

  // Define distribution
  distribCell[0] = 0;
  PDM_g_num_t stepCell = n_cell / n_rank;
  PDM_g_num_t remainderCell = n_cell % n_rank;

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i]  = stepCell;
    const int i1 = i - 1;
    if (i1 < remainderCell)
      distribCell[i] += 1;
  }

  for (int i = 1; i < n_rank + 1; i++) {
    distribCell[i] += distribCell[i-1];
  }

  PDM_g_num_t _dn_cell = distribCell[i_rank+1] - distribCell[i_rank];
  *npts = (int) _dn_cell;

  const double step = length / (double) n_faceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dn_cell);
  *coord = malloc (sizeof(double)      * _dn_cell * 3);

  if (randomize) {
    for (int i = 0; i < *npts; i++) {
      (*g_num)[i] = 1 + i + distribCell[i_rank];
      for (int j = 0; j < 3; j++) {
        (*coord)[3*i+j] = origin[j] + length * _rand01();
      }
    }
  }

  else {
    int _npts = 0;
    for (PDM_g_num_t g = distribCell[i_rank]; g < distribCell[i_rank+1]; g++) {
      PDM_g_num_t i = g % n_faceSeg;
      PDM_g_num_t j = ((g - i) % n_faceFace) / n_faceSeg;
      PDM_g_num_t k = (g - i - n_faceSeg * j) / n_faceFace;

      (*coord)[3 * _npts    ] = (i + 0.5) * step + origin[0];
      (*coord)[3 * _npts + 1] = (j + 0.5) * step + origin[1];
      (*coord)[3 * _npts + 2] = (k + 0.5) * step + origin[2];
      (*g_num)[_npts++] = g + 1;//1 + i + n_faceSeg * j + n_faceFace * k;
    }
  }

  free (distribCell);
}


static void
_gen_cloud_surf
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_faceSeg,
 const double        origin[3],
 const double        length,
 const int           randomize,
 int                *npts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank, i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);


  PDM_g_num_t n_faceFace = n_faceSeg * n_faceSeg;
  PDM_g_num_t n_face     = 6 * n_faceFace;

  PDM_g_num_t step_rank = n_face / n_rank;
  PDM_g_num_t remainder = n_face % n_rank;

  PDM_g_num_t *distrib = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));
  distrib[0] = 0;
  for (int i = 0; i < n_rank; i++) {
    int n = step_rank;
    if (i < remainder) n++;
    distrib[i+1] = distrib[i] + n;
  }

  PDM_g_num_t _dn_face = distrib[i_rank+1] - distrib[i_rank];

  const double step = length / (double) n_faceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dn_face);
  *coord = malloc (sizeof(double)      * _dn_face * 3);

  int _npts = 0;
  for (PDM_g_num_t g = distrib[i_rank]; g < distrib[i_rank+1]; g++) {
    int k = g / n_faceFace;
    int dir0 = k / 3;
    int dir1 = (dir0 + 1) % 3;
    int dir2 = (dir0 + 2) % 3;
    int sgn = k % 2;

    PDM_g_num_t i = g % n_faceSeg;
    PDM_g_num_t j = ((g - i) % n_faceFace) / n_faceSeg;

    (*coord)[3 * _npts + dir0] = origin[dir0] + sgn * length;
    if (randomize) {
      (*coord)[3 * _npts + dir1] = origin[dir1] + length * _rand01();
      (*coord)[3 * _npts + dir2] = origin[dir2] + length * _rand01();
    }

    else {
      (*coord)[3 * _npts + dir1] = (i + 0.5) * step + origin[dir1];
      (*coord)[3 * _npts + dir2] = (j + 0.5) * step + origin[dir2];
    }
    (*g_num)[_npts++] = g + 1;
  }

  *npts = _npts;

  free (distrib);
}


static void
_closest_point_par
(
 PDM_MPI_Comm         comm,
 const int            n_max_per_leaf,
 const int            local_search_fun,
 const int            n_part_src,
 const int           *n_src,
 const double       **src_coord,
 const PDM_g_num_t  **src_g_num,
 const int            n_part_tgt,
 const int           *n_tgt,
 const double       **tgt_coord,
 const PDM_g_num_t  **tgt_g_num,
 PDM_g_num_t       ***closest_point_g_num,
 double            ***closest_point_dist2
 )
{
  /* Build parallel octree */
  const int depth_max = 31;

  int octree_id = PDM_para_octree_create (n_part_src,
                                          depth_max,
                                          n_max_per_leaf,
                                          0,//use_neighbours,
                                          comm);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_para_octree_point_cloud_set (octree_id,
                                     i_part,
                                     n_src[i_part],
                                     src_coord[i_part],
                                     src_g_num[i_part]);
  }

  /* Compute global extents of source and target point clouds */
  double local_min[3] = { HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double local_max[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
  for (int i_part = 0; i_part < n_part_src; i_part++) {
    const double *x = src_coord[i_part];
    for (int i = 0; i < n_src[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    const double *x = tgt_coord[i_part];
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        if (*x < local_min[j])
          local_min[j] = *x;
        if (*x > local_max[j])
          local_max[j] = *x;
        x++;
      }
    }
  }

  double global_extents[6];
  PDM_MPI_Allreduce(local_min, global_extents,     3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(local_max, global_extents + 3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  /* Break symmetry */
  double max_range = 0.;
  for (int i = 0; i < 3; i++) {
    max_range = PDM_MAX (max_range, global_extents[i+3] - global_extents[i]);
  }
  for (int i = 0; i < 3; i++) {
    global_extents[i]   -= max_range * 1.1e-3;
    global_extents[i+3] += max_range * 1.0e-3;
  }

  PDM_para_octree_build (octree_id, global_extents);
  //PDM_para_octree_dump (octree_id);
  PDM_para_octree_dump_times (octree_id);




  /* Concatenate partitions */
  int _n_tgt = 0;

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    _n_tgt += n_tgt[i_part];
  }

  double      *_tgt_coord = malloc (sizeof(double)      * _n_tgt * 3);
  PDM_g_num_t *_tgt_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  PDM_g_num_t *_closest_src_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  double      *_closest_src_dist2 = malloc (sizeof(double)      * _n_tgt);

  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        _tgt_coord[_n_tgt + 3*i + j] = tgt_coord[i_part][3*i + j];
      }
      _tgt_g_num[_n_tgt + i] = tgt_g_num[i_part][i];
    }
    _n_tgt += n_tgt[i_part];
  }


  /* Search closest source points */
  PDM_para_octree_single_closest_point (octree_id,
                                        (const _local_search_fun_t) local_search_fun,
                                        _n_tgt,
                                        _tgt_coord,
                                        _tgt_g_num,
                                        _closest_src_g_num,
                                        _closest_src_dist2);

  /* Restore partitions */
  free (_tgt_coord);
  free (_tgt_g_num);

  *closest_point_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *closest_point_dist2 = (double **) malloc (sizeof(double *) * n_part_tgt);
  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    (*closest_point_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_tgt[i_part]);
    (*closest_point_dist2)[i_part] = (double *) malloc (sizeof(double) * n_tgt[i_part]);

    for (int i = 0; i < n_tgt[i_part]; i++) {
      (*closest_point_g_num)[i_part][i] = _closest_src_g_num[_n_tgt + i];
      (*closest_point_dist2)[i_part][i] = _closest_src_dist2[_n_tgt + i];
    }
    _n_tgt += n_tgt[i_part];
  }
  free (_closest_src_g_num);
  free (_closest_src_dist2);

  /* Free parallel octree */
  PDM_para_octree_free (octree_id);
}



static void
_closest_point_seq
(
 PDM_MPI_Comm         comm,
 const int            n_part_src,
 const int           *n_src,
 const double       **src_coord,
 const PDM_g_num_t  **src_g_num,
 const int            n_part_tgt,
 const int           *n_tgt,
 const double       **tgt_coord,
 const PDM_g_num_t  **tgt_g_num,
 PDM_g_num_t       ***closest_point_g_num,
 double            ***closest_point_dist2
 )
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  const double tolerance = 1e-4;
  const int depth_max = 31;
  const int points_in_leaf_max = 4;

  int octree_id = PDM_octree_create (n_part_src,
                                     depth_max,
                                     points_in_leaf_max,
                                     tolerance,
                                     comm);

  for (int i_part = 0; i_part < n_part_src; i_part++) {
    PDM_octree_point_cloud_set (octree_id,
                                i_part,
                                n_src[i_part],
                                src_coord[i_part],
                                src_g_num[i_part]);
  }

  /* Build octree */
  PDM_timer_t *timer = PDM_timer_create ();
  double t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  PDM_octree_build (octree_id);

  PDM_timer_hang_on (timer);
  double t_end = PDM_timer_elapsed (timer) - t_begin;
  t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  double t_max;
  PDM_MPI_Reduce (&t_end, &t_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, 0, PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf ("Octree build elapsed time = %.3gs\n", t_max);
  }

  /* Concatenate partitions */
  int _n_tgt = 0;

  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    _n_tgt += n_tgt[i_part];
  }

  double      *_tgt_coord = malloc (sizeof(double)      * _n_tgt * 3);
  PDM_g_num_t *_tgt_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  PDM_g_num_t *_closest_src_g_num = malloc (sizeof(PDM_g_num_t) * _n_tgt);
  double      *_closest_src_dist2 = malloc (sizeof(double)      * _n_tgt);

  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    for (int i = 0; i < n_tgt[i_part]; i++) {
      for (int j = 0; j < 3; j++) {
        _tgt_coord[_n_tgt + 3*i + j] = tgt_coord[i_part][3*i + j];
      }
      _tgt_g_num[_n_tgt + i] = tgt_g_num[i_part][i];
    }
    _n_tgt += n_tgt[i_part];
  }

  /* Search closest source points */
  PDM_timer_hang_on (timer);
  t_begin = PDM_timer_elapsed (timer);
  PDM_timer_resume (timer);

  PDM_octree_closest_point (octree_id,
                            _n_tgt,
                            _tgt_coord,
                            _tgt_g_num,
                            _closest_src_g_num,
                            _closest_src_dist2);

  PDM_timer_hang_on (timer);
  t_end = PDM_timer_elapsed (timer) - t_begin;
  PDM_timer_resume (timer);

  PDM_MPI_Reduce (&t_end, &t_max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, 0, PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf ("Octree closest point elapsed time = %.3gs\n", t_max);
  }
  PDM_timer_free (timer);

  /* Restore partitions */
  free (_tgt_coord);
  free (_tgt_g_num);

  *closest_point_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *closest_point_dist2 = (double **) malloc (sizeof(double *) * n_part_tgt);
  _n_tgt = 0;
  for (int i_part = 0; i_part < n_part_tgt; i_part++) {
    (*closest_point_g_num)[i_part] = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * n_tgt[i_part]);
    (*closest_point_dist2)[i_part] = (double *) malloc (sizeof(double) * n_tgt[i_part]);

    for (int i = 0; i < n_tgt[i_part]; i++) {
      (*closest_point_g_num)[i_part][i] = _closest_src_g_num[_n_tgt + i];
      (*closest_point_dist2)[i_part][i] = _closest_src_dist2[_n_tgt + i];
    }
    _n_tgt += n_tgt[i_part];
  }
  free (_closest_src_g_num);
  free (_closest_src_dist2);

  /* Free octree */
  PDM_octree_free (octree_id);
}



static void
_write_point_cloud
(
 const char        *filename,
 const char        *header,
 const int          n_pts,
 const double       coord[],
 const PDM_g_num_t  g_num[]
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

  fclose(f);
}


static void
_read_point_cloud
(
 const char   *filename,
 int          *n_pts,
 double      **coord,
 PDM_g_num_t **g_num
 )
{
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Unable to open %s", filename);
  }

  char line[999];
  while (fgets(line, sizeof(line), f) != NULL) {
    if (strcmp(line,"\n") == 0 && strcmp(line,"\r\n") == 0) {
      continue;
    }

    if (strstr(line, "POINTS") != NULL) {
      int stat = sscanf(line,
                        "%*[^0123456789]%d%*[^0123456789]",
                        n_pts);
      assert (stat);

      *coord = malloc (sizeof(double) * (*n_pts) * 3);
      for (int i = 0; i < *n_pts; i++) {
        fscanf(f, "%lf %lf %lf",
               *coord + 3*i,
               *coord + 3*i + 1,
               *coord + 3*i + 2);
      }
    }

    if (strstr(line, "CELL_DATA") != NULL) {

      *g_num = malloc (sizeof(PDM_g_num_t) * (*n_pts));
      fgets(line, sizeof(line), f);
      fgets(line, sizeof(line), f);
      for (int i = 0; i < *n_pts; i++) {
        fscanf(f, PDM_FMT_G_NUM, *g_num + i);
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
  PDM_MPI_Init (&argc, &argv);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  srand (time(NULL) + i_rank);


  /*
   *  Set default values
   */
  PDM_g_num_t n_face_seg     = 10;
  double      length         = 1.;
  PDM_g_num_t n_tgt          = 10;
  int         n_max_per_leaf = 10;
  int         surf_source    = 0;
  int         randomize      = 0;
  int         repeat_last    = 0;
  int         n_methods      = 5;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &n_face_seg,
              &length,
              &n_tgt,
              &n_max_per_leaf,
              &surf_source,
              &randomize,
              &repeat_last,
              &n_methods);

  /* Random cube origin */
  double origin[3] = {0., 0., 0.};
  if (0) {
    if (i_rank == 0) {
      for (int i = 0; i < 3; i++) {
        origin[i] = 2.*_rand01() - 1.;
      }
    }

    PDM_MPI_Bcast (origin, 3, PDM_MPI_DOUBLE, 0, PDM_MPI_COMM_WORLD);
  }

  /* Random cube dimensions */
  double dimensions[3] = {1., 1., 1.};
  if (0) {
    if (i_rank == 0) {
      for (int i = 0; i < 3; i++) {
        double r = 2.*_rand01() - 1.;
        dimensions[i] = pow(2., r);
      }
      //printf("dimensions = %f %f %f\n", dimensions[0], dimensions[1], dimensions[2]);
    }

    PDM_MPI_Bcast (dimensions, 3, PDM_MPI_DOUBLE, 0, PDM_MPI_COMM_WORLD);
  }



  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  double      *tgt_coord = NULL;
  PDM_g_num_t *tgt_g_num = NULL;
  int _n_tgt;

  if (repeat_last) {
    char filename[999];
    sprintf(filename, "last/tgt_%3.3d.vtk", i_rank);

    _read_point_cloud (filename,
                       &_n_tgt,
                       &tgt_coord,
                       &tgt_g_num);
  }

  else {
    _gen_cloud_random (PDM_MPI_COMM_WORLD,
                       n_tgt,
                       origin,
                       length,
                       &_n_tgt,
                       &tgt_g_num,
                       &tgt_coord);

    for (int i = 0; i < _n_tgt; i++) {
      for (int j = 0; j < 3; j++) {
        tgt_coord[3*i+j] = origin[j] + dimensions[j] * (tgt_coord[3*i+j] - origin[j]);
      }
    }
  }

  /*
   *  Define the source point cloud
   */
  int n_part_src = 1;
  int _n_src;
  double *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;

  if (repeat_last) {
    char filename[999];
    sprintf(filename, "last/src_%3.3d.vtk", i_rank);

    _read_point_cloud (filename,
                       &_n_src,
                       &src_coord,
                       &src_g_num);
  }

  else {

    if (surf_source) {
      _gen_cloud_surf (PDM_MPI_COMM_WORLD,
                       n_face_seg,
                       origin,
                       length,
                       randomize,
                       &_n_src,
                       &src_g_num,
                       &src_coord);
    }

    else {
      _gen_cloud_vol (PDM_MPI_COMM_WORLD,
                      n_face_seg,
                      origin,
                      length,
                      randomize,
                      &_n_src,
                      &src_g_num,
                      &src_coord);
    }

    for (int i = 0; i < _n_src; i++) {
      for (int j = 0; j < 3; j++) {
        src_coord[3*i+j] = origin[j] + dimensions[j] * (src_coord[3*i+j] - origin[j]);
      }
    }

  }



  if (!repeat_last) {
    char filename[999];

    sprintf(filename, "last/tgt_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "tgt",
                        _n_tgt,
                        tgt_coord,
                        tgt_g_num);

    sprintf(filename, "last/src_%3.3d.vtk", i_rank);
    _write_point_cloud (filename,
                        "src",
                        _n_src,
                        src_coord,
                        src_g_num);
  }


  /*
   *  Compare methods
   */
  double *elapsed = malloc (sizeof(double) * n_methods);
  double t_begin, t_end;
  PDM_timer_t *timer = PDM_timer_create ();

  PDM_g_num_t **closest_point_g_num = NULL;
  double      **closest_point_dist2 = NULL;

  double elapsed_min = HUGE_VAL;
  for (int method = 0; method < n_methods; method++) {

    if (i_rank == 0) {
      printf("\n\nMethod %d\n", method);
    }

    if (method > 0) PDM_timer_hang_on (timer);
    t_begin = PDM_timer_elapsed (timer);
    PDM_timer_resume (timer);

    /* Compute closest points */
    if (method == 0) {
      _closest_point_seq (PDM_MPI_COMM_WORLD,
                          n_part_src,
                          (const int *) &_n_src,
                          (const double **) &src_coord,
                          (const PDM_g_num_t **) &src_g_num,
                          n_part_tgt,
                          (const int *) &_n_tgt,
                          (const double **) &tgt_coord,
                          (const PDM_g_num_t **) &tgt_g_num,
                          &closest_point_g_num,
                          &closest_point_dist2);
    }

    else if (method > 0) {
      _closest_point_par (PDM_MPI_COMM_WORLD,
                          n_max_per_leaf,
                          method - 1,
                          n_part_src,
                          (const int *) &_n_src,
                          (const double **) &src_coord,
                          (const PDM_g_num_t **) &src_g_num,
                          n_part_tgt,
                          (const int *) &_n_tgt,
                          (const double **) &tgt_coord,
                          (const PDM_g_num_t **) &tgt_g_num,
                          &closest_point_g_num,
                          &closest_point_dist2);
    }

    PDM_timer_hang_on (timer);
    t_end = PDM_timer_elapsed (timer) - t_begin;
    PDM_timer_resume (timer);

    PDM_MPI_Reduce (&t_end,
                    elapsed + method, 1,
                    PDM_MPI_DOUBLE,
                    PDM_MPI_MAX,
                    0,
                    PDM_MPI_COMM_WORLD);

    if (i_rank == 0) {
      printf ("\nTotal elapsed time = %.3gs\n", elapsed[method]);
      if (elapsed[method] < elapsed_min) {
        elapsed_min = elapsed[method];
      }
    }

    /* Check */
    if (!randomize) {
      if (i_rank == 0) {
        printf("-- Check\n");
        fflush(stdout);
      }

      int marge = 3;
      int ijk0;
      int ijk_lo[3];
      int ijk_hi[3];
      double cell_side = length / ((double) n_face_seg);
      double coord[3];

      PDM_g_num_t _n_wrong = 0;
      for (int itgt = 0; itgt < _n_tgt; itgt++) {
        /*printf ("[%d] pt ("PDM_FMT_G_NUM") closest src : ("PDM_FMT_G_NUM") at dist2 = %f\n",
          i_rank, tgt_g_num[itgt], closest_point_g_num[0][itgt], closest_point_dist2[0][itgt]);*/
        PDM_g_num_t true_closest_src_g_num;
        double      true_closest_src_dist2 = HUGE_VAL;

        // find i,j,k of the cell that contains the target point
        // and define search region
        if (surf_source) {
          int k;
          double max_val = -HUGE_VAL;
          for (int idim = 0; idim < 3; idim++) {
            double val1 = PDM_ABS (tgt_coord[3*itgt+idim] - origin[idim]);
            double val2 = PDM_ABS (tgt_coord[3*itgt+idim] - origin[idim] - length);

            if (val1 > max_val) {
              k = 2*idim;
              max_val = val1;
            }

            if (val2 > max_val) {
              k = 2*idim+1;
              max_val = val2;
            }

          }

          int dir0 = k / 3;
          int dir1 = (dir0 + 1) % 3;
          int dir2 = (dir0 + 2) % 3;
          int sgn = k % 2;

          coord[dir0] = origin[dir0] + sgn * length * dimensions[dir0];
          for (int j = -marge; j < marge; j++) {
            int _j = PDM_MAX (0, PDM_MIN (j, n_face_seg-1));
            coord[dir2] = (_j + 0.5) * cell_side * dimensions[dir2] + origin[dir2];
            for (int i = -marge; i < marge; i++) {
              int _i = PDM_MAX (0, PDM_MIN (i, n_face_seg-1));
              coord[dir1] = (_i + 0.5) * cell_side * dimensions[dir1] + origin[dir1];

              PDM_g_num_t g_num = 1 + i + n_face_seg*(j + n_face_seg*k);
              double dist2 = 0.;
              for (int idim = 0; idim < 3; idim++) {
                double delta = tgt_coord[3*itgt+idim] - coord[idim];
                dist2 += delta * delta;
              }

              if (dist2 < true_closest_src_dist2) {
                true_closest_src_dist2 = dist2;
                true_closest_src_g_num = g_num;
              }
            }
          }

        }

        else {
          for (int idim = 0; idim < 3; idim++) {
            ijk0 = (int) floor(tgt_coord[3*itgt+idim] / (cell_side * dimensions[idim]));
            ijk0 = PDM_MIN (PDM_MAX (ijk0, 0), n_face_seg-1);

            if (ijk0 < marge) {
              ijk_lo[idim] = 0;
              ijk_hi[idim] = 2*marge;
            } else if (ijk0 > n_face_seg - marge) {
              ijk_hi[idim] = n_face_seg;
              ijk_lo[idim] = n_face_seg - 2*marge;
            } else {
              ijk_lo[idim] = ijk0 - marge;
              ijk_hi[idim] = ijk0 + marge;
            }
          }

          // inspect search region
          for (int k = ijk_lo[2]; k < ijk_hi[2]; k++) {
            coord[2] = (k + 0.5) * cell_side * dimensions[2] + origin[2];
            for (int j = ijk_lo[1]; j < ijk_hi[1]; j++) {
              coord[1] = (j + 0.5) * cell_side * dimensions[1] + origin[1];
              for (int i = ijk_lo[0]; i < ijk_hi[0]; i++) {
                coord[0] = (i + 0.5) * cell_side * dimensions[0] + origin[0];

                PDM_g_num_t g_num = 1 + i + n_face_seg*(j + n_face_seg*k);
                double dist2 = 0.;
                for (int idim = 0; idim < 3; idim++) {
                  double delta = tgt_coord[3*itgt+idim] - coord[idim];
                  dist2 += delta * delta;
                }

                if (dist2 < true_closest_src_dist2) {
                  true_closest_src_dist2 = dist2;
                  true_closest_src_g_num = g_num;
                }

              }
            }
          }
        }

        // check
        //if (true_closest_src_g_num != closest_point_g_num[0][itgt]) {
        if (closest_point_dist2[0][itgt] >  true_closest_src_dist2 &&
            closest_point_g_num[0][itgt] != true_closest_src_g_num) {
          _n_wrong++;
          double d0 = sqrt(true_closest_src_dist2);
          double d1 = sqrt(closest_point_dist2[0][itgt]);

          printf ("[%d] ERROR pt ("PDM_FMT_G_NUM") : dist2 = %3.g | %3.g (rel. err. = %3.g), gnum = "PDM_FMT_G_NUM" | "PDM_FMT_G_NUM"\n",
                  i_rank,
                  tgt_g_num[itgt],
                  closest_point_dist2[0][itgt],
                  true_closest_src_dist2,
                  (d1 - d0)/d0,
                  closest_point_g_num[0][itgt],
                  true_closest_src_g_num);
        }
      }

      PDM_g_num_t n_wrong;
      PDM_MPI_Reduce (&_n_wrong,
                      &n_wrong,
                      1,
                      PDM__PDM_MPI_G_NUM,
                      PDM_MPI_SUM,
                      0,
                      PDM_MPI_COMM_WORLD);
      if (i_rank == 0) {
        PDM_g_num_t wrong_percentage = 100 * n_wrong;
        if (n_tgt > 0) {
          wrong_percentage /= n_tgt;
        }

        printf("\nn_wrong = "PDM_FMT_G_NUM" / "PDM_FMT_G_NUM" ("PDM_FMT_G_NUM"%%)\n",
               n_wrong,
               n_tgt,
               wrong_percentage);
        fflush(stdout);
      }
    }

    for (int i = 0; i < n_part_tgt; i++) {
      free (closest_point_g_num[i]);
      free (closest_point_dist2[i]);
    }
    free (closest_point_g_num);
    free (closest_point_dist2);
  }
  PDM_timer_free (timer);


  /*
   *  Summary
   */
  if (i_rank == 0) {
    printf("\n\n\n");
    for (int method = 0; method < n_methods; method++) {
      printf ("method %d: elapsed = %.3gs, relative to min = %.3f\n",
              method, elapsed[method], elapsed[method] / elapsed_min);
    }
  }



  /*
   *  Finalize
   */
  free (elapsed);
  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
