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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_array.h"

#include "pdm_mesh_intersection_vol_vol_atomic.h"

/*============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private functions
 *============================================================================*/

 // for main

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}

static void
_read_args
(
 int    argc,
 char **argv
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}



/*============================================================================
 * Main
 *============================================================================*/

int main(int argc, char *argv[])
{
  // Init
  PDM_MPI_Init(&argc, &argv);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

  // Triangle: B->C->D->B cyclic linked list
  // Geogebra
  // double pt0[3] = {1.5, 1, 0};
  // double pt1[3] = {0, 0, 0};
  // double pt2[3] = {0.8, 0.3, 0.4};
  // inside
  // double pt0[3] = {0.5, 0, 0};
  // double pt2[3] = {0, 0.5, 0};
  // double pt1[3] = {0, 0, 0.5};
  // XYZ
  // double pt0[3] = {1, 0, 0};
  // double pt1[3] = {0, 1, 0};
  // double pt2[3] = {0, 0, 1};
  // dbg
  // double pt0[3] = {0, 0, 0};
  // double pt1[3] = {-1, 1, 0};
  // double pt2[3] = {1, -1, 1};
  // dbg 2
  // double pt0[3] = {0.5, 0.5,0.5};
  // double pt1[3] = {0.5,-0.5,0.5};
  // double pt2[3] = {1.5,-0.5,0.5};
  // double pt0[3] = {0.000000,-1.000000,1.500000};
  // double pt1[3] = {-1.000000,-1.000000,1.500000};
  // double pt2[3] = {-1.000000,0.000000,0.500000};
  // double pt0[3] = {0.400000,-0.400000,0.900000};
  // double pt1[3] = {0.400000,0.600000,-0.100000};
  // double pt2[3] = {1.400000,-0.400000,-0.100000};
  // double pt0[3] = {1.100000,-0.000000,-0.000000};
  // double pt1[3] = {0.100000,1.000000,-0.000000};
  // double pt2[3] = {0.100000,-0.000000,1.000000};
  // double pt0[3] = {0.9, -0.1, 0.2};
  // double pt1[3] = {0.9, 0.5, -0.4};
  // double pt2[3] = {0.3, 0.5, 0.2};
  double pt0[3] = {0.000000,1.000000,0.000000};
  double pt1[3] = {1.000000,0.000000,0.000000};
  double pt2[3] = {0.160000,0.160000,0.760000};

  double triaB_coord[9];
  memcpy(triaB_coord + 3*0, pt0, sizeof(double)*3);
  memcpy(triaB_coord + 3*1, pt1, sizeof(double)*3);
  memcpy(triaB_coord + 3*2, pt2, sizeof(double)*3);


  double  *vtx_coordA = NULL;
  int      n_vtxA     = 0;
  int     *face_vtxA  = NULL;
  int      n_faceA    = 0;
  double  *vtx_coordB = NULL;
  int      n_vtxB     = 0;
  int     *face_vtxB  = NULL;
  int      n_faceB    = 0;
  double vol1 = PDM_mesh_intersection_vol_vol_atomic_compute(triaB_coord,
                                                             &vtx_coordA,
                                                             &n_vtxA,
                                                             &face_vtxA,
                                                             &n_faceA,
                                                             &vtx_coordB,
                                                             &n_vtxB,
                                                             &face_vtxB,
                                                             &n_faceB);

  log_trace("vol1 = %20.16f\n", vol1);

  double vol2 = PDM_mesh_intersection_vol_vol_atomic_compute2(triaB_coord);

  log_trace("vol2 = %20.16f\n", vol2);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}

