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
#include "pdm_closest_points.h"


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
        long _n_faceSeg = atol(argv[i]);
        *n_faceSeg = (PDM_g_num_t) _n_faceSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
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
    else if (strcmp(argv[i], "-m") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *method = atoi(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



static void
_gen_cloud
(
 const int           n_pts,
 const double        origin[3],
 const double        length,
 const int           n_rank,
 const int           i_rank,
 double            **pts_coord,
 int                *_n_pts
 )
{
  *_n_pts = (int) (n_pts/n_rank);
  if (i_rank < n_pts%n_rank) {
    (*_n_pts)++;
  }
  *pts_coord = malloc (sizeof(double) * 3 * (*_n_pts));

  int idx = 0;
  for (int i = 0; i < n_pts; i++) {
    for (int j = 0; j < 3; j++) {
      double x = origin[j] + length * (double) rand() / ((double) RAND_MAX);
      if (i%n_rank == i_rank) {
        (*pts_coord)[idx++] = x;
      }
    }
  }
}



static void
_gen_cube_cell_centers
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   n_faceSeg,
 const double        origin[3],
 const double        length,
 int                *npts,
 PDM_g_num_t       **g_num,
 double            **coord
 )
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  PDM_g_num_t *distribCell = (PDM_g_num_t *) malloc((n_rank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t n_cell     = n_faceSeg * n_faceSeg * n_faceSeg;
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

  const double step = length / (double) n_faceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dn_cell);
  *coord = malloc (sizeof(double)      * _dn_cell * 3);

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

  *npts = _npts;

  free (distribCell);
}





static void
_single_closest_point
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
  /* Build parallel octree */
  const int depth_max = 31;
  const int points_in_leaf_max = 1;
  const int build_leaf_neighbours = 1;

  int octree_id = PDM_para_octree_create (n_part_src,
                                          depth_max,
                                          points_in_leaf_max,
                                          build_leaf_neighbours,
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
  PDM_octree_build (octree_id);

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
  PDM_octree_closest_point (octree_id,
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

  /* Free octree */
  PDM_octree_free (octree_id);
}


static void
_closest_point
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
  /* Init closest points structure */
  int id = PDM_closest_points_create (comm,
                                      1,//n_closest_points
                                      PDM_OWNERSHIP_KEEP);

  /* Set point clouds */
  PDM_closest_points_n_part_cloud_set (id,
                                       n_part_src,
                                       n_part_tgt);

  for (int ipart = 0; ipart < n_part_tgt; ipart++) {
    PDM_closest_points_tgt_cloud_set (id,
                                      ipart,
                                      n_tgt[ipart],
                                      tgt_coord[ipart],
                                      tgt_g_num[ipart]);
  }

  for (int ipart = 0; ipart < n_part_src; ipart++) {
    PDM_closest_points_src_cloud_set (id,
                                      ipart,
                                      n_src[ipart],
                                      src_coord[ipart],
                                      src_g_num[ipart]);
  }

  /* Compute closest point */
  PDM_closest_points_compute (id);

  PDM_closest_points_dump_times (id);

  *closest_point_g_num = (PDM_g_num_t **) malloc (sizeof(PDM_g_num_t *) * n_part_tgt);
  *closest_point_dist2 = (double **) malloc (sizeof(double *) * n_part_tgt);
  for (int ipart = 0; ipart < n_part_tgt; ipart++) {
    PDM_closest_points_get (id,
                            ipart,
                            *closest_point_g_num + ipart,
                            *closest_point_dist2 + ipart);
  }
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


  /*
   *  Set default values
   */
  PDM_g_num_t n_face_seg = 10;
  double      length     = 1.;
  PDM_g_num_t n_tgt      = 10;
  int         method     = 0;

  /*
   *  Read args
   */
  _read_args (argc,
              argv,
              &n_face_seg,
              &length,
              &n_tgt,
              &method);

  double origin[3] = {0., 0., 0.};
  if (1) {
    for (int i = 0; i < 3; i++) {
      origin[3] = 2. * (double) rand() / ((double) RAND_MAX) - 1.;
    }
  }

  /*
   *  Define the target point cloud
   */
  int n_part_tgt = 1;
  double *tgt_coord = NULL;
  int _n_tgt;
  _gen_cloud (n_tgt,
              origin,
              length,
              n_rank,
              i_rank,
              &tgt_coord,
              &_n_tgt);
  int id_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

  double *tgt_char_length = malloc (sizeof(double) * _n_tgt);

  for (int i = 0; i < _n_tgt; i++) {
    tgt_char_length[i] = length * 1.e-6;
  }

  PDM_gnum_set_from_coords (id_gnum, 0, _n_tgt, tgt_coord, tgt_char_length);

  PDM_gnum_compute (id_gnum);

  PDM_g_num_t *tgt_g_num = PDM_gnum_get (id_gnum, 0);

  PDM_gnum_free (id_gnum);



  /*
   *  Define the source point cloud
   */
  int n_part_src = 1;
  int _n_src;
  double *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;

  _gen_cube_cell_centers (PDM_MPI_COMM_WORLD,
                          n_face_seg,
                          origin,
                          length,
                          &_n_src,
                          &src_g_num,
                          &src_coord);


  /*
   *  Compute closest point
   */
  PDM_g_num_t **closest_point_g_num = NULL;
  double      **closest_point_dist2 = NULL;
  if (method == 0) {
    if (i_rank == 0) printf ("Method : New \n");
    _single_closest_point (PDM_MPI_COMM_WORLD,
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
  } else if (method == 1) {
    if (i_rank == 0) printf ("Method : kNN with k = 1\n");
    _closest_point (PDM_MPI_COMM_WORLD,
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
  } else {
    if (i_rank == 0) printf ("Method : serial octree\n");
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



  /*
   *  Check results
   */
  if (i_rank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  int marge = 2;
  int ijk0;
  int ijk_lo[3];
  int ijk_hi[3];
  double cell_side = length / ((double) n_face_seg);
  double cell_ctr[3];

  PDM_g_num_t _n_wrong = 0;
  for (int itgt = 0; itgt < _n_tgt; itgt++) {
    /*printf ("[%d] pt ("PDM_FMT_G_NUM") closest src : ("PDM_FMT_G_NUM") at dist2 = %f\n",
      i_rank, tgt_g_num[itgt], closest_point_g_num[0][itgt], closest_point_dist2[0][itgt]);*/

    // find i,j,k of the cell that contains the target point
    // and define search region
    for (int idim = 0; idim < 3; idim++) {
      ijk0 = (int) floor(tgt_coord[3*itgt+idim] / cell_side);
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
    PDM_g_num_t true_closest_src_g_num;
    double      true_closest_src_dist2 = HUGE_VAL;
    for (int k = ijk_lo[2]; k < ijk_hi[2]; k++) {
      cell_ctr[2] = (k + 0.5) * cell_side;
      for (int j = ijk_lo[1]; j < ijk_hi[1]; j++) {
        cell_ctr[1] = (j + 0.5) * cell_side;
        for (int i = ijk_lo[0]; i < ijk_hi[0]; i++) {
          cell_ctr[0] = (i + 0.5) * cell_side;

          PDM_g_num_t g_num = 1 + i + n_face_seg*(j + n_face_seg*k);
          double dist2 = 0.;
          for (int idim = 0; idim < 3; idim++) {
            double delta = tgt_coord[3*itgt+idim] - cell_ctr[idim];
            dist2 += delta * delta;
          }

          if (dist2 < true_closest_src_dist2) {
            true_closest_src_dist2 = dist2;
            true_closest_src_g_num = g_num;
          }

        }
      }
    }

    // check
    if (true_closest_src_g_num != closest_point_g_num[0][itgt]) {
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




  /*
   *  Finalize
   */
  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);
  for (int i = 0; i < n_part_tgt; i++) {
    free (closest_point_g_num[i]);
    free (closest_point_dist2[i]);
  }
  free (closest_point_g_num);
  free (closest_point_dist2);

  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
