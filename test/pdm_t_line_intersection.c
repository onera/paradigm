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

#include "pdm_line.h"
#include "pdm_predicate.h"


static PDM_line_intersect_t
_PDM_line_intersection_2drobust
(
 const double  a[3],
 const double  b[3],
 const double  c[3],
 const double  d[3],
 double       *u,
 double       *v
 )
{
  double ha = PDM_predicate_orient2d (c, d, a);
  double hb = PDM_predicate_orient2d (c, d, b);

  /*if ((ha < 0 && hb < 0) || (ha > 0 && hb > 0)) {
    return PDM_LINE_INTERSECT_NO;
  }*/


  double hc = PDM_predicate_orient2d (a, b, c);
  double hd = PDM_predicate_orient2d (a, b, d);

  printf("ha = %g, hb = %g, hc = %g, hf = %g\n", ha, hb, hc, hd);

  if ((ha < 0 && hb < 0) || (ha > 0 && hb > 0)) {
    return PDM_LINE_INTERSECT_NO;
  }
  if ((hc < 0 && hd < 0) || (hc > 0 && hd > 0)) {
    return PDM_LINE_INTERSECT_NO;
  }

  if (ha == hb) {
    // ===> ha = hb = 0
    // =?=> hc = hd = 0
    printf("ha = hb = %g, hc = %g, hd = %g\n", ha, hc, hd);
    return PDM_LINE_INTERSECT_ON_LINE;
  }


  *u = ha / (ha - hb);
  *v = hc / (hc - hd);
  printf("  u = %g, v = %g\n", *u, *v);

  /*
   * Check parametric coordinates for intersection.
   */

  if ( (0.0 <= *u) && (*u <= 1.0) && (0.0 <= *v) && (*v <= 1.0) ) {
    return PDM_LINE_INTERSECT_YES;
  }
  else {
    return PDM_LINE_INTERSECT_NO;
  }
}







static void
_read_args (int            argc,
            char         **argv,
            int           *n_test,
            double        *scale,
            double        *offset)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {
    if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        exit(EXIT_FAILURE);
      else
        *n_test = atoi(argv[i]);
    }
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc)
        exit(EXIT_FAILURE);
      else
        *scale = atof(argv[i]);
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



static const char *intersect_status[3] = {"NO", "YES", "ON_LINE"};


int main (int argc, char *argv[])
{
  PDM_predicate_exactinit();

  int    n_test = 1000;
  double scale  = 1.;
  double offset = 0.;
  _read_args (argc, 
              argv,
              &n_test,
              &scale,
              &offset);


  double u1, v1, u2, v2;
  double pts[12];
  double *a1 = pts;
  double *a2 = pts + 3;
  double *b1 = pts + 6;
  double *b2 = pts + 9;

  int n_diff = 0;

  for (int i = 0; i < n_test; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 2; k++) {
        pts[3*j + k] = scale * (2.*(double) rand() / (double) RAND_MAX - 1.);
      }
      pts[3*j + 1] = 0.;
      pts[3*j + 2] = 0.;
    }

    for (int j = 0; j < 2; j++) {
      double r = (double) rand() / (double) RAND_MAX;
      for (int k = 0; k < 2; k++) {
        pts[3*(j+2) + k] = (1-r)*a1[k] + r*a2[k];
      }
    }

    PDM_line_intersect_t status1 = PDM_line_intersection (a1,
                                                          a2,
                                                          b1,
                                                          b2,
                                                          &u1,
                                                          &v1);

    PDM_line_intersect_t status2 = PDM_line_intersection_2drobust (a1,
                                                                   a2,
                                                                   b1,
                                                                   b2,
                                                                   &u2,
                                                                   &v2);

    printf("%s - %s\n", intersect_status[status1], intersect_status[status2]);

    n_diff += (status1 != status2);
  }

  printf("n_diff = %d / %d\n", n_diff, n_test);

  return 0;
}

