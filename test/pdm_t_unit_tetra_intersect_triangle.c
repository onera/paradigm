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

/*============================================================================
 * Global variables
 *============================================================================*/

static int test_debug = 1;

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct Element Element;
struct Element
{
  double coord[3];
  Element *next;
};

typedef struct List List;
struct List
{
    Element *head;
};

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
 int            argc,
 char         **argv
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

// for cyclic linked list

static void
_print_cll
(
 List *cll
)
{
  if (cll == NULL)
  {
    exit(EXIT_FAILURE);
  }

  Element *current = cll->head;
  printf("head->");
  while (1) {
    printf("[(%f,%f,%f),%p]->", current->coord[0], current->coord[1], current->coord[2], (void *) current->next);
    if (current->next == cll->head) break;
    current = current->next;
  }
  printf("head\n");
}

// to factorize

static double
_plane_equation_function
(
 double coord[3],
 int    i
)
{
  switch (i) {
  // OYZ
  case 0:
    return coord[0];
  // X+Y=1
  case 1:
    return 1 - coord[0] - coord[1];
  // OZX
  case 2:
    return coord[1];
  // OXY
  case 3:
    return coord[2];
  // X+Y+Z=1
  case 4:
    return 1 - coord[0] - coord[1] - coord[2];
  default:
    PDM_error(__FILE__, __LINE__, 0, "Only 5 planes are considered\n");
    return 0;
  }
}

// for intersection

// --> if A or outside is NULL, it implies Vcolumn = 0
static void
_determine_A_outside
(
 Element  **cll_storage,
 int        idx,
 List      *cll,
 List     **A,
 List     **outside
)
{

  // check if the triangle is inside the tetrahedron
  Element *current = cll->head;
  while (1) {
    // 0<= x, y, z <= 1 and 1-x-y-z >= 0
    int cond0 = current->coord[0] >= 0 && current->coord[0] <= 1;
    int cond1 = current->coord[1] >= 0 && current->coord[1] <= 1;
    int cond2 = current->coord[2] >= 0 && current->coord[2] <= 1;
    int cond3 = 1-current->coord[0]-current->coord[1]-current->coord[2] >= 0;
    // not totally in unit tetrahedron
    if (!(cond0 && cond1 && cond2 && cond3)) {
      break;
    }
    if (current->next == cll->head) break;
    current = current->next;
  }
  if (current->next == cll->head) {

    if (test_debug) {
      printf("triangle in tetrahedra\n");
    }

    // A = cll;
    *A       = cll;
    *outside = NULL;
    return;
  }

  // intersect with all 5 planes
  for (int i = 0; i < 5; i++) {
    // maximum 2 intersections
    int count_intersect = 0;
    current = cll->head;

    Element *I10 = cll_storage[idx++];
    Element *I11 = cll_storage[idx++];
    Element *I20 = cll_storage[idx++];
    Element *I21 = cll_storage[idx++];

    while (count_intersect < 2) {
      double *coord1 = current->coord;
      double *coord2 = current->next->coord;

      // segment in plane
      if (_plane_equation_function(coord1, i) == 0 && _plane_equation_function(coord2, i) == 0) {

        if (test_debug) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) is on plane %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

        // polygon in plane
        if (_plane_equation_function(current->next->next->coord, i) == 0) {

          if (test_debug) {
            printf("polygon is in plane %d\n", i);
          }


          free(A);
          free(outside);
          *A       = NULL;
          *outside = NULL;
          return;
        }

        // only segment in plane
        else {

          // outside
          if (_plane_equation_function(current->next->next->coord, i) <= 0) {

            if (test_debug) {
              printf("polygon is outside %d\n", i);
            }

            free(A);
            free(outside);
            *A       = NULL;
            *outside = NULL;
            return;

          }

          // inside
          else {

            if (test_debug) {
              printf("polygon is inside %d, cll unchanged\n", i);
            }

            free(outside);
            *A       = cll;
            *outside = NULL;

            // count_intersect += 2;
            // idx -= 4;
            break; // skip to next plane

          }
        }
      } // end if segment in plane

      // intersection
      else if ((_plane_equation_function(coord1, i) >= 0 && _plane_equation_function(coord2, i) <= 0) ||
               (_plane_equation_function(coord1, i) <= 0 && _plane_equation_function(coord2, i) >= 0)) {

        if (test_debug) {
          printf("intersection with %d\n", i);
        }

        double t = _plane_equation_function(coord2, i)/(_plane_equation_function(coord2, i)-_plane_equation_function(coord1, i));

        double pt[3] = {t*coord1[0]+(1-t)*coord2[0],
                        t*coord1[1]+(1-t)*coord2[1],
                        t*coord1[2]+(1-t)*coord2[2]};

        // current-> (outside) intersection (inside) -> current->next
        if (_plane_equation_function(current->coord, i) <= 0) {
          if (count_intersect == 0) {
            *outside   = cll;
            (*A)->head = I11;

            I11->next = current->next;
            current->next = I10->next;
            I10->next = NULL;

            I10->coord[0] = pt[0];
            I10->coord[1] = pt[1];
            I10->coord[2] = pt[2];

            I11->coord[0] = pt[0];
            I11->coord[1] = pt[1];
            I11->coord[2] = pt[2];

          } else {

            I10->next = I20;
            I20->next = current->next;
            current->next = I21->next;
            I21->next = (*A)->head;

            I20->coord[0] = pt[0];
            I20->coord[1] = pt[1];
            I20->coord[2] = pt[2];

            I21->coord[0] = pt[0];
            I21->coord[1] = pt[1];
            I21->coord[2] = pt[2];

          }
        }

        // current-> (inside) intersection (outside) -> current->next
        else {
          if (count_intersect == 0) {
            (*outside)->head = I11;
            *A               = cll;

            I11->next = current->next;
            current->next = I10->next;
            I10->next = NULL;

            I10->coord[0] = pt[0];
            I10->coord[1] = pt[1];
            I10->coord[2] = pt[2];

            I11->coord[0] = pt[0];
            I11->coord[1] = pt[1];
            I11->coord[2] = pt[2];

          } else {

            I10->next = I20;
            I20->next = current->next;
            current->next = I11->next;
            I11->next = (*outside)->head;

            I20->coord[0] = pt[0];
            I20->coord[1] = pt[1];
            I20->coord[2] = pt[2];

            I21->coord[0] = pt[0];
            I21->coord[1] = pt[1];
            I21->coord[2] = pt[2];

          }

        }
      }

      // nothing
      else {

        if (test_debug) {
          printf("no intersection with %d\n", i);
        }

      }

      // to correcly loop over cll
      cll = *A;
      if (current->next == cll->head) break;
      current = current->next;

    } // end while loop
  } // end for loop

}

/*============================================================================
 * Main
 *============================================================================*/

int main(int argc, char *argv[])
{
  // Init
  PDM_MPI_Init(&argc, &argv);

  // Malloc for cll
  int max_size = 10*2 + 3;
  int idx      = 0;
  Element **cll_storage = malloc(sizeof(Element *) * max_size);
  for (int i = 0; i < max_size; i++) {
    cll_storage[i] = malloc(sizeof(Element));
  }

  // Triangle: B->C->D->B cyclic linked list
  // Geogebra
  double pt0[3] = {1.5, 1, 0};
  double pt1[3] = {0, 0, 0};
  double pt2[3] = {0.8, 0.3, 0.4};
  // inside
  // double pt0[3] = {0.5, 0, 0};
  // double pt1[3] = {0, 0.5, 0};
  // double pt2[3] = {0, 0, 0.5};

  Element *B = cll_storage[idx++];
  memcpy(B->coord, pt0, sizeof(double)*3);

  Element *C = cll_storage[idx++];
  memcpy(C->coord, pt1, sizeof(double)*3);
  B->next = C;

  Element *D = cll_storage[idx++];
  memcpy(D->coord, pt2, sizeof(double)*3);
  D->next = B;
  C->next = D;

  List *cll = malloc(sizeof(List));
  cll->head = B;

  // debug
  if (test_debug) {
    printf("Triangle:\n");
    _print_cll(cll);
  }

  // Determine A and outside (B before projection)
  List *A       = NULL;//malloc(sizeof(List
  List *outside = malloc(sizeof(List));
  _determine_A_outside(cll_storage, idx, cll, &A, &outside);

  // debug
  if (test_debug) {
    printf("A: ");
    if (A == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(A);
    }

    printf("outside: ");
    if (outside == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(outside);
    }
  }

  // free
  for (int i = 0; i < max_size; i++) {
    free(cll_storage[i]);
  }
  free(cll_storage);
  free(cll);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}

