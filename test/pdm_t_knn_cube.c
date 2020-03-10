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
#include "pdm_part.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_closest_points.h"
#include "pdm_dcube_gen.h"
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
     "  -c      <level>  Number of closest points (default : 10).\n\n"
     "  -t      <level>  Number of Target points (default : 10).\n\n"
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
 * \param [inout]   nFaceSeg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   nClosest Number of closest points
 * \param [inout]   nTgt     Number of Target points
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nFaceSeg,
           double        *length,
           int           *nClosest,
           PDM_g_num_t   *nTgt,
           int           *n_part,
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
        long _nFaceSeg = atol(argv[i]);
        *nFaceSeg = (PDM_g_num_t) _nFaceSeg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nClosest = atoi(argv[i]);
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
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
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
_gen_clouds_random
(
 const int         nPts,
 const double      length,
 const int         numProcs,
 const int         myRank,
 const PDM_g_num_t nFaceSeg,
 double          **pts_coord,
 int              *nPts_l
 )
{
  *nPts_l = (int) (nPts/numProcs);
  *pts_coord = malloc (sizeof(double) * 3 * (*nPts_l));
  double *_pts_coord = *pts_coord;

  double offset = 0.5 * length / ((double) nFaceSeg);
  double length2 = length - 2*offset;

  int idx = 0;
  for (int i = 0; i < numProcs*(*nPts_l); i++) {
    for (int j = 0; j < 3; j++) {
      double x = offset + length2 * (double) rand() / ((double) RAND_MAX);
      if (i%numProcs == myRank) {
        _pts_coord[idx++] = x;
      }
    }
  }
}


static void
_gen_cube_cell_centers
(
 PDM_MPI_Comm       comm,
 const PDM_g_num_t  nFaceSeg,
 const double       length,
 const double       zero_x,
 const double       zero_y,
 const double       zero_z,
 int               *npts,
 PDM_g_num_t      **g_num,
 double           **coord
 )
{
  int nRank;
  int myRank;

  PDM_MPI_Comm_size(comm, &nRank);
  PDM_MPI_Comm_rank(comm, &myRank);

  PDM_g_num_t *distribCell = (PDM_g_num_t *) malloc((nRank + 1) * sizeof(PDM_g_num_t));

  PDM_g_num_t nCell     = nFaceSeg * nFaceSeg * nFaceSeg;
  PDM_g_num_t nFaceFace = nFaceSeg * nFaceSeg;

  // Define distribution
  distribCell[0] = 0;
  PDM_g_num_t stepCell = nCell / nRank;
  PDM_g_num_t remainderCell = nCell % nRank;

  for (int i = 1; i < nRank + 1; i++) {
    distribCell[i]  = stepCell;
    const int i1 = i - 1;
    if (i1 < remainderCell)
      distribCell[i] += 1;
  }

  for (int i = 1; i < nRank + 1; i++) {
    distribCell[i] += distribCell[i-1];
  }

  PDM_g_num_t _dNCell = distribCell[myRank+1] - distribCell[myRank];

  const double step = length / (double) nFaceSeg;

  *g_num = malloc (sizeof(PDM_g_num_t) * _dNCell);
  *coord = malloc (sizeof(double)      * _dNCell * 3);

  int _npts = 0;
  for (PDM_g_num_t g = distribCell[myRank]; g < distribCell[myRank+1]; g++) {
    PDM_g_num_t i = g % nFaceSeg;
    PDM_g_num_t j = ((g - i) % nFaceFace) / nFaceSeg;
    PDM_g_num_t k = (g - i - nFaceSeg * j) / nFaceFace;

    (*coord)[3 * _npts    ] = (i + 0.5) * step + zero_x;
    (*coord)[3 * _npts + 1] = (j + 0.5) * step + zero_y;
    (*coord)[3 * _npts + 2] = (k + 0.5) * step + zero_z;
    (*g_num)[_npts++] = g + 1;//1 + i + nFaceSeg * j + nFaceFace * k;
  }

  *npts = _npts;

  free (distribCell);
}





/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  int myRank;
  int numProcs;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(PDM_MPI_COMM_WORLD, &myRank);
  PDM_MPI_Comm_size(PDM_MPI_COMM_WORLD, &numProcs);

  /*
   *  Set default values
   */

  PDM_g_num_t  nFaceSeg  = 10;
  double        length  = 1.;
  int           nPart   = 1;
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#else
#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#endif
#endif

  int n_closest_points = 10;
  PDM_g_num_t nTgt = 10;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nFaceSeg,
             &length,
             &n_closest_points,
             &nTgt,
             &nPart,
             (int *) &method);


  /* Define the target point cloud */
  double *tgt_coords = NULL;
  int _nTgt_l;
  _gen_clouds_random (nTgt,
                      length,
                      numProcs,
                      myRank,
                      nFaceSeg,
                      &tgt_coords,
                      &_nTgt_l);

  int id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *tgt_char_length = malloc(sizeof(double) * _nTgt_l);

  for (int i = 0; i < _nTgt_l; i++) {
    tgt_char_length[i] = length * 1.e-6;
  }

  PDM_gnum_set_from_coords (id, 0, _nTgt_l, tgt_coords, tgt_char_length);

  PDM_gnum_compute (id);

  PDM_g_num_t *tgt_gnum = PDM_gnum_get(id, 0);

  PDM_gnum_free (id, 1);



  /*
   *  Define the source point cloud (cell centers of cube)
   */
  nPart = 1;
  int _nSrc_l;
  double *src_coords = NULL;
  PDM_g_num_t *src_gnum = NULL;
  _gen_cube_cell_centers (PDM_MPI_COMM_WORLD,
                          nFaceSeg,
                          length,
                          0.,
                          0.,
                          0.,
                          &_nSrc_l,
                          &src_gnum,
                          &src_coords);

  n_closest_points = PDM_MIN (n_closest_points, nFaceSeg*nFaceSeg*nFaceSeg);


  /* Init closest points structure */
  int id2 = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                       n_closest_points);

  PDM_closest_points_n_part_cloud_set (id2,
                                       nPart,
                                       1);

  // set tgt point cloud
  /*printf("[%d] nTgt = %d, gnum, coords =\n", myRank, _nTgt_l);
  for (int i = 0; i < _nTgt_l; i++) {
    printf(" %ld\t[%f %f %f]\n",
           tgt_gnum[i], tgt_coords[3*i], tgt_coords[3*i+1], tgt_coords[3*i+2]);
  }
  printf("\n\n\n\n");*/
  PDM_closest_points_tgt_cloud_set (id2,
                                    0,
                                    _nTgt_l,
                                    tgt_coords,
                                    tgt_gnum);

  // set src point cloud
  /*printf("[%d] nSrc = %d, gnum, coords =\n", myRank, _nSrc_l);
  for (int i = 0; i < _nSrc_l; i++) {
    printf(" %ld\t[%f %f %f]\n",
           src_gnum[i], src_coords[3*i], src_coords[3*i+1], src_coords[3*i+2]);
  }
  printf("\n");*/
  PDM_closest_points_src_cloud_set (id2,
                                    0,
                                    _nSrc_l,
                                    src_coords,
                                    src_gnum);

  /* Compute closest points */
  PDM_closest_points_compute (id2);

  PDM_closest_points_dump_times (id2);

  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (id2,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);

#if 1
  /* Check results */
  if (myRank == 0) {
    printf("-- Check\n");
    fflush(stdout);
  }

  PDM_g_num_t n_wrong = 0;
  PDM_g_num_t n_tgt = 0;

  PDM_g_num_t *true_closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_closest_points);
  double      *true_closest_src_dist = malloc (sizeof(double) * n_closest_points);
  int n_cells_radius = (int) ceil(0.5 * pow((double) n_closest_points, 1./3.));

  int ijk0;
  int ijk_lo[3];
  int ijk_hi[3];
  double cell_side = length / ((double) nFaceSeg);
  double cell_ctr[3];

  for (int itgt = 0; itgt < _nTgt_l; itgt++) {
    /*printf("[%d] %d/%d\n", myRank, itgt, _nTgt_l);*/

    n_tgt++;
    int wrong = 0;

    for (int l = 0; l < n_closest_points; l++) {
      true_closest_src_dist[l] = HUGE_VAL;
    }

    // find i,j,k of the cell that contains the target point
    // and define search region
    for (int idim = 0; idim < 3; idim++) {
      ijk0 = (int) floor(tgt_coords[3*itgt+idim] / cell_side);
      ijk0 = PDM_MIN (PDM_MAX (ijk0, 0), nFaceSeg-1);

      if (ijk0 < n_cells_radius) {
        ijk_lo[idim] = 0;
        ijk_hi[idim] = 2*n_cells_radius;
      } else if (ijk0 > nFaceSeg - n_cells_radius) {
        ijk_hi[idim] = nFaceSeg;
        ijk_lo[idim] = nFaceSeg - 2*n_cells_radius;
      } else {
        ijk_lo[idim] = ijk0 - n_cells_radius;
        ijk_hi[idim] = ijk0 + n_cells_radius;
      }
    }

    // inspect search region
    for (int k = ijk_lo[2]; k < ijk_hi[2]; k++) {
      cell_ctr[2] = (k + 0.5) * cell_side;
      for (int j = ijk_lo[1]; j < ijk_hi[1]; j++) {
        cell_ctr[1] = (j + 0.5) * cell_side;
        for (int i = ijk_lo[0]; i < ijk_hi[0]; i++) {
          cell_ctr[0] = (i + 0.5) * cell_side;

          PDM_g_num_t gnum = 1 + i + nFaceSeg*j + nFaceSeg*nFaceSeg*k;

          double dist = 0;
          for (int idim = 0; idim < 3; idim++) {
            double delta = tgt_coords[3*itgt+idim] - cell_ctr[idim];
            dist += delta * delta;
          }

          // insertion sort
          if (dist < true_closest_src_dist[n_closest_points-1]) {

            int l = n_closest_points - 1;
            while (l > 0 && dist < true_closest_src_dist[l-1]) {
              true_closest_src_gnum[l] = true_closest_src_gnum[l-1];
              true_closest_src_dist[l] = true_closest_src_dist[l-1];
              l--;
            }

            true_closest_src_gnum[l] = gnum;
            true_closest_src_dist[l] = dist;
          }
        }
      }
    }

    // check
    for (int l = 0; l < n_closest_points; l++) {
      /*printf("(%ld) %d/%d : %ld\t%f\n",
             tgt_gnum[itgt],
             l+1,
             n_closest_points,
             closest_src_gnum[n_closest_points*itgt + l],
             closest_src_dist[n_closest_points*itgt + l]);*/

      //if (closest_src_gnum[n_closest_points*itgt + l] != true_closest_src_gnum[l]) {
      /*if (closest_src_dist[n_closest_points*itgt + l] > true_closest_src_dist[l]) {
        printf("*** ERROR (%ld) [%f %f %f]: %f / %f (relative err. = %f)\t\t%ld / %ld\n",
        tgt_gnum[itgt],
        tgt_coords[3*itgt], tgt_coords[3*itgt+1], tgt_coords[3*itgt+2],
        closest_src_dist[n_closest_points*itgt + l],
        true_closest_src_dist[l],
        (sqrt(closest_src_dist[n_closest_points*itgt + l]) - sqrt(true_closest_src_dist[l]))/sqrt(true_closest_src_dist[l]),
        closest_src_gnum[n_closest_points*itgt + l],
        true_closest_src_gnum[l]);
        }*/
      //assert (closest_src_gnum[n_closest_points*itgt + l] == true_closest_src_gnum[l]);
      //assert (closest_src_dist[n_closest_points*itgt + l] <= true_closest_src_dist[l]);

      if (closest_src_dist[n_closest_points*itgt + l] > true_closest_src_dist[l]) {
        printf("[%d] (%ld) [%f %f %f]: %f / %f (relative err. = %f)\t\t%ld / %ld\n",
               myRank,
               tgt_gnum[itgt],
               tgt_coords[3*itgt], tgt_coords[3*itgt+1], tgt_coords[3*itgt+2],
               closest_src_dist[n_closest_points*itgt + l],
               true_closest_src_dist[l],
               (sqrt(closest_src_dist[n_closest_points*itgt + l]) - sqrt(true_closest_src_dist[l]))/sqrt(true_closest_src_dist[l]),
               closest_src_gnum[n_closest_points*itgt + l],
               true_closest_src_gnum[l]);

        wrong = 1;
      }
    }

    if (wrong) {
      n_wrong++;
    }
  }
  free (true_closest_src_gnum);
  free (true_closest_src_dist);

  /*PDM_g_num_t wrong_percentage = 0;
  if (n_tgt > 0) {
    wrong_percentage = 100 * n_wrong / n_tgt;
  }
  printf("[%d] n_wrong = %ld / %ld (%ld%%)\n",
         myRank,
         n_wrong,
         n_tgt,
         wrong_percentage);*/

  PDM_g_num_t n_wrong_total, n_tgt_total;

  PDM_MPI_Reduce (&n_tgt,
                  &n_tgt_total,
                  1,
                  PDM__PDM_MPI_G_NUM,
                  PDM_MPI_SUM,
                  0,
                  PDM_MPI_COMM_WORLD);

  PDM_MPI_Reduce (&n_wrong,
                  &n_wrong_total,
                  1,
                  PDM__PDM_MPI_G_NUM,
                  PDM_MPI_SUM,
                  0,
                  PDM_MPI_COMM_WORLD);

  if (myRank == 0) {
    PDM_g_num_t wrong_percentage_total = 0;
    if (n_tgt_total > 0) {
      wrong_percentage_total = 100 * n_wrong_total / n_tgt_total;
    }
    printf("n_wrong = %ld / %ld (%ld%%)\n",
           n_wrong_total,
           n_tgt_total,
           wrong_percentage_total);

    assert (n_wrong_total < 1);
  }

#endif







#if 0
  /* Export as vtk */
  if (myRank == 0) {
    printf("-- Export for visu\n");
    fflush(stdout);
  }

  char filename[9999];
  FILE *f;
  /*sprintf(filename, "/home/bandrieu/workspace/paradigma-dev/test/para_octree/src_%4.4d.vtk", myRank);

  f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "src\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", _nSrc_l);
  for (int i = 0; i < _nSrc_l; i++) {
    for (int j = 0; j < 3; j++)
      fprintf(f, "%f ", src_coords[3*i+j]);
    fprintf(f, "\n");
  }


  fprintf(f, "CELLS %d %d\n", _nSrc_l, 2*_nSrc_l);
  for (int i = 0; i < _nSrc_l; i++)
    fprintf(f, "1 %d\n", i);

  fprintf(f, "CELL_TYPES %d\n", _nSrc_l);
  for (int i = 0; i < _nSrc_l; i++)
    fprintf(f, "1\n");

    fclose(f);*/



  // send everything to rank 0
  int *recv_src_count = malloc(sizeof(int) * numProcs);
  PDM_MPI_Gather (&_nSrc_l,       1, PDM_MPI_INT,
                  recv_src_count, 1, PDM_MPI_INT,
                  0, PDM_MPI_COMM_WORLD);

  int *recv_tgt_count = malloc(sizeof(int) * numProcs);
  PDM_MPI_Gather (&_nTgt_l,       1, PDM_MPI_INT,
                  recv_tgt_count, 1, PDM_MPI_INT,
                  0, PDM_MPI_COMM_WORLD);

  int *recv_src_g_num_shift = malloc(sizeof(int) * (numProcs+1));
  int *recv_tgt_g_num_shift = malloc(sizeof(int) * (numProcs+1));
  int *recv_src_coord_shift = malloc(sizeof(int) * (numProcs+1));
  int *recv_tgt_coord_shift = malloc(sizeof(int) * (numProcs+1));
  int *recv_closest_src_shift = malloc(sizeof(int) * (numProcs+1));
  if (myRank == 0) {
    recv_src_g_num_shift[0] = 0;
    recv_tgt_g_num_shift[0] = 0;
    recv_src_coord_shift[0] = 0;
    recv_tgt_coord_shift[0] = 0;
    recv_closest_src_shift[0] = 0;
    for (int i = 0; i < numProcs; i++) {
      recv_src_g_num_shift[i+1] = recv_src_g_num_shift[i] + recv_src_count[i];
      recv_tgt_g_num_shift[i+1] = recv_tgt_g_num_shift[i] + recv_tgt_count[i];
    }
  }

  PDM_g_num_t *recv_src_g_num = NULL;
  if (myRank == 0)
    recv_src_g_num = malloc(sizeof(PDM_g_num_t) * recv_src_g_num_shift[numProcs]);

  PDM_MPI_Gatherv (src_gnum, _nSrc_l, PDM__PDM_MPI_G_NUM,
                   recv_src_g_num, recv_src_count, recv_src_g_num_shift, PDM__PDM_MPI_G_NUM,
                   0, PDM_MPI_COMM_WORLD);


  PDM_g_num_t *recv_tgt_g_num = NULL;
  if (myRank == 0)
    recv_tgt_g_num = malloc(sizeof(PDM_g_num_t) * recv_tgt_g_num_shift[numProcs]);
  PDM_MPI_Gatherv (tgt_gnum, _nTgt_l, PDM__PDM_MPI_G_NUM,
                   recv_tgt_g_num, recv_tgt_count, recv_tgt_g_num_shift, PDM__PDM_MPI_G_NUM,
                   0, PDM_MPI_COMM_WORLD);

  if (myRank == 0) {
    for (int i = 0; i < numProcs; i++) {
      recv_src_count[i] *= 3;
      recv_tgt_count[i] *= 3;
      recv_src_coord_shift[i+1] = recv_src_coord_shift[i] + recv_src_count[i];
      recv_tgt_coord_shift[i+1] = recv_tgt_coord_shift[i] + recv_tgt_count[i];
    }
  }

  double *recv_src_coord = NULL;
  if (myRank == 0)
    recv_src_coord = malloc(sizeof(double) * recv_src_coord_shift[numProcs]);
  PDM_MPI_Gatherv (src_coords, 3*_nSrc_l, PDM_MPI_DOUBLE,
                   recv_src_coord, recv_src_count, recv_src_coord_shift, PDM_MPI_DOUBLE,
                   0, PDM_MPI_COMM_WORLD);

  double *recv_tgt_coord = NULL;
  if (myRank == 0)
    recv_tgt_coord = malloc(sizeof(double) * recv_tgt_coord_shift[numProcs]);
  PDM_MPI_Gatherv (tgt_coords, 3*_nTgt_l, PDM_MPI_DOUBLE,
                   recv_tgt_coord, recv_tgt_count, recv_tgt_coord_shift, PDM_MPI_DOUBLE,
                   0, PDM_MPI_COMM_WORLD);

  if (myRank == 0) {
    for (int i = 0; i < numProcs; i++) {
      recv_tgt_count[i] = (recv_tgt_count[i] / 3) * n_closest_points;
      recv_closest_src_shift[i+1] = recv_closest_src_shift[i] + recv_tgt_count[i];
    }
  }

  PDM_g_num_t *recv_closest_src_gnum = NULL;
  double *recv_closest_src_dist = NULL;
  if (myRank == 0) {
    recv_closest_src_gnum = malloc(sizeof(PDM_g_num_t) * recv_closest_src_shift[numProcs]);
    recv_closest_src_dist = malloc(sizeof(double) * recv_closest_src_shift[numProcs]);
  }
  PDM_MPI_Gatherv (closest_src_gnum, n_closest_points*_nTgt_l, PDM__PDM_MPI_G_NUM,
                   recv_closest_src_gnum, recv_tgt_count, recv_closest_src_shift, PDM__PDM_MPI_G_NUM,
                   0, PDM_MPI_COMM_WORLD);

  PDM_MPI_Gatherv (closest_src_dist, n_closest_points*_nTgt_l, PDM_MPI_DOUBLE,
                   recv_closest_src_dist, recv_tgt_count, recv_closest_src_shift, PDM_MPI_DOUBLE,
                   0, PDM_MPI_COMM_WORLD);


  if (myRank == 0) {
    double *g_src_coord = malloc(sizeof(double) * recv_src_coord_shift[numProcs]);
    for (int i = 0; i < recv_src_g_num_shift[numProcs]; i++) {
      PDM_g_num_t ii = recv_src_g_num[i] - 1;
      for (int j = 0; j < 3; j++)
        g_src_coord[3*ii + j] = recv_src_coord[3*i + j];
    }
    free (recv_src_coord);
    free (recv_src_g_num);


    double *g_tgt_coord = malloc(sizeof(double) * recv_tgt_coord_shift[numProcs]);
    for (int i = 0; i < recv_tgt_g_num_shift[numProcs]; i++) {
      PDM_g_num_t ii = recv_tgt_g_num[i] - 1;
      for (int j = 0; j < 3; j++)
        g_tgt_coord[3*ii + j] = recv_tgt_coord[3*i + j];
    }
    free (recv_tgt_coord);

    PDM_g_num_t *g_closest_src_gnum = malloc(sizeof(PDM_g_num_t) * recv_closest_src_shift[numProcs]);
    double *g_closest_src_dist = malloc(sizeof(double) * recv_closest_src_shift[numProcs]);
    for (int i = 0; i < recv_tgt_g_num_shift[numProcs]; i++) {
      PDM_g_num_t ii = recv_tgt_g_num[i] - 1;
      for (int j = 0; j < n_closest_points; j++) {
        g_closest_src_gnum[n_closest_points*ii + j] = recv_closest_src_gnum[n_closest_points*i + j];
        g_closest_src_dist[n_closest_points*ii + j] = recv_closest_src_dist[n_closest_points*i + j];
      }
    }

    free (recv_tgt_g_num);


    int n_g_pts = recv_src_g_num_shift[numProcs] + recv_tgt_g_num_shift[numProcs];



    sprintf(filename, "/home/bandrieu/workspace/paradigma-dev/test/para_octree/valid_knn.vtk");


    f = fopen(filename, "w");

    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "validation_knn\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d double\n", n_g_pts+1);
    for (int i = 0; i < recv_src_g_num_shift[numProcs]; i++) {
      for (int j = 0; j < 3; j++)
        fprintf(f, "%f ", g_src_coord[3*i + j]);
      fprintf(f, "\n");
    }
    for (int i = 0; i < recv_tgt_g_num_shift[numProcs]; i++) {
      for (int j = 0; j < 3; j++)
        fprintf(f, "%f ", g_tgt_coord[3*i + j]);
      fprintf(f, "\n");
    }

    for (int j = 0; j < 3; j++)
      fprintf(f, "%f ", -99999.);
    fprintf(f, "\n");

    int n_edges = n_closest_points * recv_tgt_g_num_shift[numProcs];
    int n_nodes = recv_src_g_num_shift[numProcs];
    int n_cells = n_edges + n_nodes;

    fprintf(f, "CELLS %d %d\n", n_cells, 3*n_edges + 2*n_nodes);
    for (int i = 0; i < recv_src_g_num_shift[numProcs]; i++) {
      fprintf(f, "1 %d\n", i);
    }

    for (int i = 0; i < recv_tgt_g_num_shift[numProcs]; i++) {
      for (int j = 0; j < n_closest_points; j++) {
        PDM_g_num_t gnum = g_closest_src_gnum[n_closest_points*i + j];
        if (gnum < 1)
          gnum = n_g_pts - 1;

        fprintf(f, "2 %d %ld\n",
                i + recv_src_g_num_shift[numProcs],
                gnum - 1);
        /* fprintf(f, "2 %d %d\n", i + recv_src_g_num_shift[numProcs], ?); */
      }
    }

    fprintf(f, "CELL_TYPES %d\n", n_cells);
    for (int i = 0; i < n_nodes; i++)
      fprintf(f, "1\n");
    for (int i = 0; i < n_edges; i++)
      fprintf(f, "4\n");

    fprintf(f, "CELL_DATA %d\n", n_cells);
    fprintf(f, "SCALARS distance float\n LOOKUP_TABLE default\n");
    for (int i = 0; i < recv_src_g_num_shift[numProcs]; i++) {
      fprintf(f, "%f\n", 0.);
    }
    for (int i = 0; i < recv_tgt_g_num_shift[numProcs]; i++) {
      for (int j = 0; j < n_closest_points; j++) {
        PDM_g_num_t gnum = g_closest_src_gnum[n_closest_points*i + j];
        if (gnum > 0) {
          fprintf(f, "%f\n", g_closest_src_dist[n_closest_points*i + j]);
        } else {
          fprintf(f, "%f\n", 0.0);
        }
      }
    }

    fclose(f);
  }

  free (recv_src_count);
  free (recv_tgt_count);
  free (recv_src_g_num_shift);
  free (recv_tgt_g_num_shift);
  free (recv_src_coord_shift);
  free (recv_tgt_coord_shift);

#endif






  /* Free */
  PDM_closest_points_free (id2,
                           0);

  free (tgt_coords);
  free (tgt_char_length);
  free (tgt_gnum);

  free (src_coords);
  free (src_gnum);


  PDM_MPI_Finalize();

  if (myRank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  return 0;
}
