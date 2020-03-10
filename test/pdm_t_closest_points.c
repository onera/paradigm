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
#include "pdm_closest_points.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

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
_usage
(
 int exit_code
 )
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c       <level> Number of closest points (default : 10).\n\n"
     "  -s       <level> Number of Source points (default : 10).\n\n"
     "  -t       <level> Number of Target points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -clumps          Source points distributed in clumps around target points (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nClosest   Number of closest points
 * \param [inout] nSrc   Number of Source points
 * \param [inout] nTgt   Number of Target points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 int           *nClosest,
 PDM_g_num_t   *nSrc,
 PDM_g_num_t   *nTgt,
 double        *radius,
 int           *local,
 int           *rand,
 int           *clumps
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nClosest = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nSrc = atol(argv[i]);
        *nSrc = (PDM_g_num_t) _nSrc;
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

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else if (strcmp(argv[i], "-clumps") == 0) {
      *clumps = 1;
    }


    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

/**
 *
 * \brief  Random value
 *
 *  \return a random double in [-1, 1]
 */

static double
_random01
(
)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS (rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


static void
_gen_clouds_random2
(
 const int      nSrc,
 const int      nTgt_l,
 const double   radius,
 const int      numProcs,
 const int      myRank,
 double       **src_coord,
 double       **tgt_coord,
 int           *nSrc_l
 )
{
  *nSrc_l = (int) (nSrc/numProcs);
  *src_coord = malloc (sizeof(double) * 3 * (*nSrc_l));
  double *_src_coord = *src_coord;
  double x;
  int idx = 0;
  for (int i = 0; i < numProcs*(*nSrc_l); i++) {
    for (int j = 0; j < 3; j++) {
      x = _random01() * radius;
      if (i%numProcs == myRank) {
        _src_coord[idx++] = x;
      }
    }
  }


  *tgt_coord = malloc (sizeof(double) * 3 * nTgt_l);
  double *_tgt_coord = *tgt_coord;
  for (int i = 0; i < nTgt_l; i++) {
    for (int j = 0; j < 3; j++) {
      _tgt_coord[3*i+j] = _random01() * radius;
    }
  }
}


static void
_gen_clouds_random
(
 const int      nSrc_l,
 const int      nTgt_l,
 const double   radius,
 double       **src_coord,
 double       **tgt_coord
 )
{

  *src_coord = malloc (sizeof(double) * 3 * nSrc_l);
  double *_src_coord = *src_coord;
  for (int i = 0; i < nSrc_l; i++) {
    for (int j = 0; j < 3; j++) {
      _src_coord[3*i+j] = _random01() * radius;
    }
  }

  *tgt_coord = malloc (sizeof(double) * 3 * nTgt_l);
  double *_tgt_coord = *tgt_coord;
  for (int i = 0; i < nTgt_l; i++) {
    for (int j = 0; j < 3; j++) {
      _tgt_coord[3*i+j] = _random01() * radius;
    }
  }
}











static void
_gen_clouds_clumps
(
 const int      n_closest_points,
 const int      nTgt,
 const double   radius,
 const int      numProcs,
 const int      myRank,
 const double   clump_scale,
 double       **src_coord,
 double       **tgt_coord,
 int           *nTgt_l,
 int           *nSrc_l
 )
{
  *nTgt_l = (int) nTgt/numProcs;
  *nSrc_l = n_closest_points * (*nTgt_l);

  const double clump_radius = clump_scale * radius / pow(nTgt, 1./3);

  *tgt_coord = malloc (sizeof(double) * 3 * (*nTgt_l));
  *src_coord = malloc (sizeof(double) * 3 * (*nSrc_l));

  double *_tgt_coord = *tgt_coord;
  double *_src_coord = *src_coord;
  double x;

  int ii = 0;
  int idx = 0;
  for (int i = 0; i < (*nTgt_l)*numProcs; i++) {
    for (int j = 0; j < 3; j++) {
      x = _random01() * radius;
      if (i%numProcs == myRank) {
        _tgt_coord[3*ii+j] = x;
      }
    }

    for (int k = 0; k < n_closest_points; k++) {
      for (int j = 0; j < 3; j++) {
        x = _random01() * clump_radius;

        if (i%numProcs == myRank)
          _src_coord[idx++] = _tgt_coord[3*ii+j] + x;
      }
    }

    if (i%numProcs == myRank)
      ii++;
  }
}


/*============================================================================
 * Public function definitions
 *============================================================================*/
/**
 *
 * \brief  Main
 *
 */

int
main
(
 int argc,
 char *argv[]
 )
{

  PDM_MPI_Init (&argc, &argv);

  int myRank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &myRank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);

  int n_closest_points = 10;
  PDM_g_num_t nSrc = 10;
  PDM_g_num_t nTgt = 10;
  double radius = 10.;
  int local = 0;
  int rand = 0;
  int clumps = 0;

  _read_args(argc,
             argv,
             &n_closest_points,
             &nSrc,
             &nTgt,
             &radius,
             &local,
             &rand,
             &clumps);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(myRank);
  }


  /* Define the numbers of Source/Target points */

  int _nSrc_l = 0, _nTgt_l = 0;
  if (local) {
    _nSrc_l = (int) nSrc;
    _nTgt_l = (int) nTgt;
  }
  else {
    _nSrc_l = (int) (nSrc/numProcs);
    _nTgt_l = (int) (nTgt/numProcs);
    if (myRank < nSrc%numProcs) {
      _nSrc_l += 1;
    }
    if (myRank < nTgt%numProcs) {
      _nTgt_l += 1;
    }
  }

  double *src_coords = NULL;
  double *tgt_coords = NULL;


  if (clumps) {
    const double clump_scale = 0.3;

    _gen_clouds_clumps (n_closest_points,
                        nTgt,
                        radius,
                        numProcs,
                        myRank,
                        clump_scale,
                        &src_coords,
                        &tgt_coords,
                        &_nTgt_l,
                        &_nSrc_l);
  } else {
    /*_gen_clouds_random (_nSrc_l,
      _nTgt_l,
      radius,
      &src_coords,
      &tgt_coords);*/
    _gen_clouds_random2 (nSrc,
                         _nTgt_l,
                         radius,
                         numProcs,
                         myRank,
                         &src_coords,
                         &tgt_coords,
                         &_nSrc_l);
  }


  if ( n_closest_points > _nSrc_l ) n_closest_points = (int) _nSrc_l;

  /* Source points definition  */
  int id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *src_char_length = malloc(sizeof(double) * _nSrc_l);

  for (int i = 0; i < _nSrc_l; i++) {
    src_char_length[i] = radius * 1.e-6;
  }

  PDM_gnum_set_from_coords (id, 0, _nSrc_l, src_coords, src_char_length);

  PDM_gnum_compute (id);

  PDM_g_num_t *src_gnum = PDM_gnum_get(id, 0);

  PDM_gnum_free (id, 1);



  /* Target points definition */
  id = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD);

  double *tgt_char_length = malloc(sizeof(double) * _nTgt_l);

  for (int i = 0; i < _nTgt_l; i++) {
    tgt_char_length[i] = radius * 1.e-6;
  }

  PDM_gnum_set_from_coords (id, 0, _nTgt_l, tgt_coords, tgt_char_length);

  PDM_gnum_compute (id);

  PDM_g_num_t *tgt_gnum = PDM_gnum_get(id, 0);

  PDM_gnum_free (id, 1);



  int id2 = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                       n_closest_points);

  PDM_closest_points_n_part_cloud_set (id2,
                                       1,
                                       1);

  PDM_closest_points_src_cloud_set (id2,
                                    0,
                                    _nSrc_l,
                                    src_coords,
                                    src_gnum);

  PDM_closest_points_tgt_cloud_set (id2,
                                    0,
                                    _nTgt_l,
                                    tgt_coords,
                                    tgt_gnum);


  PDM_closest_points_compute (id2);


  PDM_closest_points_dump_times (id2);

  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (id2,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


#if 0
  if (myRank >= 0) {
    printf("\n\n============================\n\n");

    for (int i = 0; i < _nTgt_l; i++) {
      printf("Target point #%d (%ld) [%f, %f, %f]\n", i, tgt_gnum[i],
             tgt_coords[3*i], tgt_coords[3*i+1], tgt_coords[3*i+2]);
      for (int j = 0; j < n_closest_points; j++)
        printf("\t%d:\t%ld\t%f\n",
               j+1,
               closest_src_gnum[n_closest_points*i + j],
               closest_src_dist[n_closest_points*i + j]);
      printf("\n\n");
    }

    printf("============================\n\n");
  }
#endif


#if 0
  /* Export as vtk */
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


  PDM_closest_points_free (id2,
                           0);





  /* Free */

  free (src_coords);
  free (src_char_length);
  free (src_gnum);
  free (tgt_coords);
  free (tgt_char_length);
  free (tgt_gnum);

  if (myRank == 0) {

    PDM_printf ("\nfin Test\n");

  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;
}
