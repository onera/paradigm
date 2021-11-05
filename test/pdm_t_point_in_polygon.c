#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_plane.h"
#include "pdm_polygon.h"
#include "pdm_predicate.h"



static void
_read_args (int            argc,
            char         **argv,
            int           *n_vtx,
            int           *n_test,
            double        *offset)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {
    if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        exit(EXIT_FAILURE);
      else
        *n_vtx = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        exit(EXIT_FAILURE);
      else
        *n_test = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-o") == 0) {
      i++;
      if (i >= argc)
        exit(EXIT_FAILURE);
      else
        *offset = atof(argv[i]);
    }
    i++;
  }
}



static void
_compute_bounds_and_normal 
(
const int     n_vtx,
const double *vtx_coord,
double        bounds[6],
double        normal[3]
)
{
  for (int j = 0; j < 3; j++) {
    bounds[2*j  ] =  DBL_MAX;
    bounds[2*j+1] = -DBL_MAX;
  }

  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      bounds[2*j]   = PDM_MIN (bounds[2*j],   vtx_coord[3*i + j]);
      bounds[2*j+1] = PDM_MAX (bounds[2*j+1], vtx_coord[3*i + j]);
    }
  }

  PDM_plane_normal (n_vtx, vtx_coord, normal);
}





#if 1
int main (int argc, char *argv[])
{
  PDM_predicate_exactinit();

  #define n_vtxA 4
  double vtx_coordA[3*n_vtxA] = {
8.97127790e+00,  9.03787915e+00,  0.00000000e+00,
8.96783146e+00,  9.04104976e+00,  0.00000000e+00,
8.97045526e+00,  9.04344279e+00,  0.00000000e+00,
8.97487487e+00,  9.04019175e+00,  0.00000000e+00};

  #define n_vtxB 4
  double vtx_coordB[3*n_vtxB] = {
8.96972e+00,  9.03665051e+00,  0.00000000e+00,
8.96456e+00,  9.04059081e+00,  0.00000000e+00,
8.96974e+00,  9.04397358e+00,  0.00000000e+00,
8.97600e+00,  9.03935956e+00,  0.00000000e+00};

  PDM_polygon_status_t statusA[n_vtxA] = {
PDM_POLYGON_INSIDE,
PDM_POLYGON_INSIDE,
PDM_POLYGON_INSIDE,
PDM_POLYGON_OUTSIDE};

  PDM_polygon_status_t statusB[n_vtxB] = {
PDM_POLYGON_OUTSIDE,
PDM_POLYGON_OUTSIDE,
PDM_POLYGON_OUTSIDE,
PDM_POLYGON_OUTSIDE};

  double boundsA[6], normalA[6];
  _compute_bounds_and_normal (n_vtxA, vtx_coordA, boundsA, normalA);

  double boundsB[6], normalB[6];
  _compute_bounds_and_normal (n_vtxB, vtx_coordB, boundsB, normalB);

  for (int i = 0; i < n_vtxA; i++) {
    PDM_polygon_status_t status1 = PDM_polygon_point_in (vtx_coordA + 3*i,
                                                         n_vtxB,
                                                         vtx_coordB,
                                                         boundsB,
                                                         normalB);

    PDM_polygon_status_t status2 = PDM_polygon_point_in_new (vtx_coordA + 3*i,
                                                             n_vtxB,
                                                             vtx_coordB,
                                                             boundsB,
                                                             normalB);
    printf("A %d : %d  -  %d   /  %d\n", i, (int) status1, (int) status2, (int) statusA[i]);
  }


  for (int i = 0; i < n_vtxB; i++) {
    PDM_polygon_status_t status1 = PDM_polygon_point_in (vtx_coordB + 3*i,
                                                         n_vtxA,
                                                         vtx_coordA,
                                                         boundsA,
                                                         normalA);

    PDM_polygon_status_t status2 = PDM_polygon_point_in_new (vtx_coordB + 3*i,
                                                             n_vtxA,
                                                             vtx_coordA,
                                                             boundsA,
                                                             normalA);
    printf("B %d : %d  -  %d   /  %d\n", i, (int) status1, (int) status2, (int) statusB[i]);
  }

  return 0;
}

#else

int main (int argc, char *argv[])
{
  PDM_predicate_exactinit();

  int    n_vtx  = 6;
  int    n_test = 1000;
  double offset = 0.;

  _read_args (argc,
              argv,
              &n_vtx,
              &n_test,
              &offset);


  /* Create polygon */
  double *vtx_coord = malloc (sizeof(double) * n_vtx * 3);
  for (int i = 0; i < n_vtx; i++) {
    double t = 2.* PDM_PI * i / (double) n_vtx;
    vtx_coord[3*i    ] = cos(t);
    vtx_coord[3*i + 1] = sin(t);
    vtx_coord[3*i + 2] = 0.;
  }


  printf("vtx_coord =\n");
  for (int i = 0; i < n_vtx; i++) {
    printf("%f %f %f\n", vtx_coord[3*i], vtx_coord[3*i+1], vtx_coord[3*i+2]);
  }


  double bounds[6] = {DBL_MAX, -DBL_MAX, 
                      DBL_MAX, -DBL_MAX, 
                      DBL_MAX, -DBL_MAX};
  for (int i = 0; i < n_vtx; i++) {
    for (int j = 0; j < 3; j++) {
      bounds[2*j]   = PDM_MIN (bounds[2*j],   vtx_coord[3*i + j]);
      bounds[2*j+1] = PDM_MAX (bounds[2*j+1], vtx_coord[3*i + j]);
    }
  }
  /*printf("bounds = %f %f %f %f %f %f\n", 
         bounds[0], bounds[1], bounds[2],
         bounds[3], bounds[4], bounds[5]);*/

  double normal[3];
  PDM_plane_normal (n_vtx, vtx_coord, normal);


  /* Check vertices */
  for (int i = 0; i < n_vtx; i++) {
    PDM_polygon_status_t status1 = PDM_polygon_point_in (vtx_coord + 3*i,
                                                         n_vtx,
                                                         vtx_coord,
                                                         bounds,
                                                         normal);

    PDM_polygon_status_t status2 = PDM_polygon_point_in_new (vtx_coord + 3*i,
                                                             n_vtx,
                                                             vtx_coord,
                                                             bounds,
                                                             normal);

    printf("vtx %d : %d  -  %d\n", 
           i, status1 == PDM_POLYGON_INSIDE, status2 == PDM_POLYGON_INSIDE);
  }

  /* Check edges */
  for (int i = 0; i < n_vtx; i++) {
    int ip = (i+1) % n_vtx;

    double vec[3] = {vtx_coord[3*ip  ] - vtx_coord[3*i  ],
                     vtx_coord[3*ip+1] - vtx_coord[3*i+1],
                     vtx_coord[3*ip+2] - vtx_coord[3*i+2]};
    double vec_in[3];
    PDM_CROSS_PRODUCT(vec_in, normal, vec);

    int n1 = 0, n2 = 0;
    for (int j = 0; j < n_test; j++) {
      double t = (double) rand() / (double) RAND_MAX;
      double pt[3] = {(1-t)*vtx_coord[3*i  ] + t*vtx_coord[3*ip  ],
                      (1-t)*vtx_coord[3*i+1] + t*vtx_coord[3*ip+1],
                      (1-t)*vtx_coord[3*i+2] + t*vtx_coord[3*ip+2]};
      for (int k = 0; k < 3; k++) {
        pt[k] += offset * vec_in[k];
      }


      PDM_polygon_status_t status1 = PDM_polygon_point_in (pt,
                                                           n_vtx,
                                                           vtx_coord,
                                                           bounds,
                                                           normal);

      PDM_polygon_status_t status2 = PDM_polygon_point_in_new (pt,
                                                               n_vtx,
                                                               vtx_coord,
                                                               bounds,
                                                               normal);
      //printf("%f %f %f  : %d  - %d\n", pt[0], pt[1], pt[2], (int) status1, (int) status2);

      n1 += (status1 == PDM_POLYGON_INSIDE);
      n2 += (status2 == PDM_POLYGON_INSIDE);
    }

    printf("edge %d : %d/%d  -  %d/%d\n", i, n1, n_test, n2, n_test);
  }



  free (vtx_coord);

  return 0;
}
#endif
