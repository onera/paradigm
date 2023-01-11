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
#include "pdm_part_to_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_closest_points.h"
#include "pdm_version.h"
#include "pdm_mesh_nodal.h"
#include "pdm_vtk.h"

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
 int           *visu
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

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}




static double _idw_interp
(
 const int     n,
       double *values,
       double *distances
 )
{
  double v = 0;
  double s = 0;
  for (int i = 0; i < n; i++) {
    double w = 1./sqrt(PDM_MAX(1e-16, distances[i]));
    v += w*values[i];
    s += w;
  }

  return v/s;
}



/*============================================================================
 * Public function definitions
 *============================================================================*/
/**
 *
 * \brief  Main
 *
 */
// @@@param[n_proc] : 1,2,3,4
// @@@param[c] : 1,2,3,10
// @@@param[s] : 10000, 20000
// @@@param[t] : 10000, 20000
int
main
(
 int argc,
 char *argv[]
 )
{

  PDM_MPI_Init (&argc, &argv);

  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  if (i_rank == 0) {
    char *version = PDM_version_get();

    printf("Version de ParaDiGM : %s\n", version);
    free(version);
  }

  int         n_closest_points = 10;
  PDM_g_num_t gn_src           = 10;
  PDM_g_num_t gn_tgt           = 10;
  double      radius           = 10.;
  int         visu             = 0;

  _read_args(argc,
             argv,
             &n_closest_points,
             &gn_src,
             &gn_tgt,
             &radius,
             &visu);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank           : %d\n", n_rank);
    PDM_printf ("  - n_closest_points : %d\n", n_closest_points);
    PDM_printf ("  - n_src            : "PDM_FMT_G_NUM"\n", gn_src);
    PDM_printf ("  - n_tgt            : "PDM_FMT_G_NUM"\n", gn_tgt);
    PDM_printf ("  - radius           : %f\n", radius);
  }


  /* Generate src and tgt point clouds */
  int          n_src;
  double      *src_coord;
  PDM_g_num_t *src_g_num;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              gn_src,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_src,
                              &src_coord,
                              &src_g_num);


  int          n_tgt;
  double      *tgt_coord;
  PDM_g_num_t *tgt_g_num;
  PDM_point_cloud_gen_random (comm,
                              123456789, // seed
                              0,         // geometric_g_num
                              gn_tgt,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_tgt,
                              &tgt_coord,
                              &tgt_g_num);


  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         n_closest_points,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set (clsp,
                                       1,
                                       1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coord,
                                    src_g_num);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num);


  PDM_closest_points_compute (clsp);


  PDM_closest_points_dump_times (clsp);

  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  if (0 == 1) {
    printf("\n\n============================\n\n");

    for (int i = 0; i < n_tgt; i++) {
      printf("Target point #%d ("PDM_FMT_G_NUM") [%f, %f, %f]\n", i, tgt_g_num[i],
             tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2]);
      for (int j = 0; j < n_closest_points; j++)
        printf("\t%d:\t"PDM_FMT_G_NUM"\t%f\n",
               j+1,
               closest_src_gnum[n_closest_points*i + j],
               closest_src_dist[n_closest_points*i + j]);
      printf("\n\n");
    }


    printf("============================\n\n");
  }


  /* Check ptp */
  if (1) {
    PDM_part_to_part_t *ptp = NULL;
    PDM_closest_points_part_to_part_get(clsp,
                                        &ptp,
                                        PDM_OWNERSHIP_KEEP);


    // define field on src cloud
    double *src_field = malloc(sizeof(double) * n_src);
    for (int i = 0; i < n_src; i++) {
      // src_field[i] = src_coord[3*i];
      src_field[i] = cos(PDM_MODULE(src_coord+3*i));
    }

    // exchange to tgt cloud
    double **recv_field = NULL;
    int request = -1;
    PDM_part_to_part_iexch(ptp,
                           PDM_MPI_COMM_KIND_P2P,
                           PDM_STRIDE_CST_INTERLACED,
                           PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                           1,
                           sizeof(double),
                           NULL,
          (const void  **) &src_field,
                           NULL,
          (      void ***) &recv_field,
                           &request);

    PDM_part_to_part_iexch_wait(ptp, request);

    // interpolate tgt field from closest src points (IDW)
    double *tgt_field = malloc(sizeof(double) * n_tgt);
    for (int i = 0; i < n_tgt; i++) {
      tgt_field[i] = _idw_interp(n_closest_points,
                                 &recv_field[0]   [n_closest_points*i],
                                 &closest_src_dist[n_closest_points*i]);
    }
    free(recv_field[0]);
    free(recv_field);

    if (visu) {
      char filename[999];

      int          n_pts[2] = {n_src, n_tgt};
      double      *coord[2] = {src_coord, tgt_coord};
      double      *field[2] = {src_field, tgt_field};
      PDM_g_num_t *g_num[2] = {src_g_num, tgt_g_num};

      const char *field_name[] = {"field"};

      for (int i = 0; i < 2; i++) {
        if (i == 0) {
          sprintf(filename, "src_cloud_%d.vtk", i_rank);
        }
        else {
          sprintf(filename, "tgt_cloud_%d.vtk", i_rank);
        }

        PDM_vtk_write_std_elements_double(filename,
                                          n_pts[i],
                                          coord[i],
                                          g_num[i],
                                          PDM_MESH_NODAL_POINT,
                                          n_pts[i],
                                          NULL,
                                          g_num[i],
                                          1,
                                          field_name,
                        (const double **) &field[i]);
      }
    }
    free(src_field);
    free(tgt_field);
  }




  PDM_closest_points_free (clsp);





  /* Free */

  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);

  if (i_rank == 0) {

    PDM_printf ("-- End\n");

  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;
}
