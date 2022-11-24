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

typedef enum
{
  TWO_EXTREMA_INTERSECTION,
  ONE_EXTREMA_INTERSECTION,
  ONE_EXTREMA_ONE_SHARP_INTERSECTION,
  ONE_SHARP_INTERSECTION,
  TWO_SHARP_INTERSECTION
} intersection_t;

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

// --> if outside is NULL, it implies Vcolumn = 0
static void
_determine_A_outside
(
 Element  **cll_storage,
 int        idx,
 List      *cll,
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

    // cll remains as is
    free(*outside);
    *outside = NULL;
    return;
  }

  // intersect with all 5 planes
  for (int i = 0; i < 5; i++) {
    // maximum 2 intersections
    int intersect_idx                = 0;
    int extrema_prev_is_in           = -1;
    intersection_t intersection_type = -1;

    current = cll->head;

    Element *in[2]  = {cll_storage[idx++], cll_storage[idx++]};
    Element *out[2] = {cll_storage[idx++], cll_storage[idx++]};

    Element *prev_in  = NULL;
    Element *prev_out = NULL;
    int idx_prev_in  = -1;
    int idx_prev_out = -1;

    while (intersect_idx < 2) {

      // current segment coordinates
      double *coord1 = current->coord;
      double *coord2 = current->next->coord;

      // current function values
      double f1 = _plane_equation_function(coord1, i);
      double f2 = _plane_equation_function(coord2, i);
      double f3 = _plane_equation_function(current->next->next->coord, i);

      // segment in plane
      if ((f1 == 0) && (f2 == 0)) {

        if (test_debug) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) is on plane %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

        // polygon in plane
        if (f3 == 0) {

          if (test_debug) {
            printf("polygon is in plane %d\n", i);
          }

          free(cll);
          free(*outside);
          cll      = NULL;
          *outside = NULL;
          return;
        }

        // only segment in plane
        else {

          // outside
          if (f3 <= 0) {

            if (test_debug) {
              printf("polygon is outside %d\n", i);
            }

            free(cll);
            free(*outside);
            cll      = NULL;
            *outside = NULL;
            return;

          }

          // inside
          else {

            if (test_debug) {
              printf("polygon is inside %d, cll unchanged\n", i);
            }

            intersection_type = TWO_EXTREMA_INTERSECTION;

            if (i == 4) {
              // cll remains as is
              free(*outside);
              *outside = NULL;
            }
            break; // skip to next plane

          }
        }
      } // end if segment in plane

      // intersection
      else if ((f1 >= 0 && f2 <= 0) ||
               (f1 <= 0 && f2 >= 0)) {

        if (test_debug) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) intersection with %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

        double t = f2/(f2-f1);

        if (t == 0 || t == 1) {

          if (intersection_type == -1) {
            printf("ONE_EXTREMA_INTERSECTION \n");
            intersection_type = ONE_EXTREMA_INTERSECTION;
          }

          if (intersect_idx == 1) {
            if (intersection_type == ONE_SHARP_INTERSECTION) {
              intersection_type = ONE_EXTREMA_ONE_SHARP_INTERSECTION;
              printf("ONE_EXTREMA_ONE_SHARP_INTERSECTION \n");
            }
          } // end if one intersection

          // now first extrema intersection
          memcpy(in[intersect_idx]->coord,  coord1, sizeof(double) * 3);
          memcpy(out[intersect_idx]->coord, coord1, sizeof(double) * 3);

          if (t == 1) { // at begin of segment t == 1
            // from inside to outside
            // can not be one because segment on plane already dealt with
            if (f2 < 0) {
              extrema_prev_is_in = 1;
              printf("from in to out (intersect_idx = %d)\n", intersect_idx);
              idx_prev_in              = intersect_idx;
              prev_in                  = current;
              out[intersect_idx]->next = current->next;
              in[intersect_idx]->next  = in[(intersect_idx+1)%2];
            }
            // from outside to inside
            else { // >= 0
              extrema_prev_is_in = 0;
              printf("from out to in (intersect_idx = %d)\n", intersect_idx);
              idx_prev_out             = intersect_idx;
              prev_out                 = current;
              in[intersect_idx]->next  = current->next;
              out[intersect_idx]->next = out[(intersect_idx+1)%2];
            }

            // increase only for t == 1, easy to handle case begin segment
            intersect_idx += 1;
          } // end if t == 1
        } // end if one extrema intersection

        else {

          if (intersect_idx == 1) {
            if (intersection_type == ONE_EXTREMA_INTERSECTION) {
              intersection_type = ONE_EXTREMA_ONE_SHARP_INTERSECTION;
              printf("ONE_EXTREMA_ONE_SHARP_INTERSECTION \n");
            }
            else if (intersection_type == ONE_SHARP_INTERSECTION) {
              intersection_type = TWO_SHARP_INTERSECTION;
              printf("TWO_SHARP_INTERSECTION \n");
            }
          } // end if one intersection

          else if (intersect_idx == 0) {
            intersection_type = ONE_SHARP_INTERSECTION;
            printf("ONE_SHARP_INTERSECTION \n");
          } // end if no intersection

          double pt[3] = {t*coord1[0]+(1-t)*coord2[0],
                          t*coord1[1]+(1-t)*coord2[1],
                          t*coord1[2]+(1-t)*coord2[2]};

          memcpy(in[intersect_idx]->coord,  pt, sizeof(double) * 3);
          memcpy(out[intersect_idx]->coord, pt, sizeof(double) * 3);

          // current-> (outside) intersection (inside) -> current->next
          if (f1 < 0) {
            printf("from out to in (intersect_idx = %d)\n", intersect_idx);
            idx_prev_out             = intersect_idx;
            prev_out                 = current;
            in[intersect_idx]->next  = current->next;
            out[intersect_idx]->next = out[(intersect_idx+1)%2];
          }
          // current-> (inside) intersection (outside) -> current->next
          else if (f1 > 0) {
            printf("from in to out (intersect_idx = %d)\n", intersect_idx);
            idx_prev_in              = intersect_idx;
            prev_in                  = current;
            out[intersect_idx]->next = current->next;
            in[intersect_idx]->next  = in[(intersect_idx+1)%2];
          }

          intersect_idx += 1;
        } // end if sharp intersection

      } // end if intersection

      // nothing
      else {

        if (test_debug) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) has no intersection with %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

      }

      // to correcly loop over cll
      if (current->next == cll->head) break;
      current = current->next;

    } // end while loop

    // connect intersection points to the linked list

    printf("in[0]: %p/ in[1]: %p\n", (void *) in[0], (void *) in[1]);
    printf("out[0]: %p/ out[1]: %p\n", (void *) out[0], (void *) out[1]);
    if (intersection_type == TWO_SHARP_INTERSECTION || intersection_type == ONE_EXTREMA_ONE_SHARP_INTERSECTION) {
      if (intersection_type == TWO_SHARP_INTERSECTION) {
        prev_in->next  = in[idx_prev_in];
        prev_out->next = out[idx_prev_out];
      }

      else if (intersection_type == ONE_EXTREMA_ONE_SHARP_INTERSECTION) {

        if (extrema_prev_is_in) {
          prev_in->next  = in[idx_prev_in]->next;
          in[idx_prev_in] = prev_in;
          prev_out->next = out[idx_prev_out];
        }
        else {
          prev_in->next = in[idx_prev_in];
          prev_out->next = out[idx_prev_out]->next;
          out[idx_prev_out] = prev_out;

        }

      }

      // update A and outside
      cll->head        = in[0];
      (*outside)->head = out[0];
    }

    if (test_debug) {
       printf("cll at plane %d: ", i);
      _print_cll(cll);
    }

  } // end for loop

}

// for volume

static double
_prism_volume
(
 double coord1[3],
 double coord2[3],
 double coord3[3]
)
{
  double parallelepiped_area = 0.5 * ((coord2[0] - coord1[0])*(coord3[1] - coord1[1]) - (coord3[0] - coord1[0])*(coord2[1] - coord1[1]));
  return ((1./3.) * (coord1[2] + coord2[2] + coord3[2]) * parallelepiped_area);
}

static double
_column_volume
(
 List *cll
)
{
  double volume = 0;

  Element *current = cll->head->next;

  while (current->next != cll->head) {

    double prism_volume = _prism_volume(cll->head->coord, current->coord, current->next->coord);
    volume += prism_volume;

    current = current->next;
  }

  return volume;
}

// for vtk

static void
_cll_to_polydata
(
  List *cll,
  char filename[999]
)
{
  int size_min = 3;

  int  n_vtx        = 0;
  double *vtx_coord = malloc(sizeof(double) * size_min * 3);

  Element *current = cll->head;

  while (1) {

    if (n_vtx > size_min -1) {
      size_min *= 2;
      vtx_coord = realloc(vtx_coord, sizeof(double) * size_min * 3);
    }

    memcpy(vtx_coord + n_vtx * 3, current->coord, sizeof(double) * 3);
    n_vtx++;

    if (current->next == cll->head) break;
    current = current->next;
  }

  int face_vtx_idx[2] = {0, n_vtx};

  int *face_vtx = malloc(sizeof(int) * n_vtx);
  for (int i = 0; i < n_vtx; i++) {
    face_vtx[i] = i + 1;
  }


  PDM_vtk_write_polydata(filename,
                         n_vtx,
                         vtx_coord,
                         NULL,
                         1,
                         face_vtx_idx,
                         face_vtx,
                         NULL,
                         NULL);

  free(vtx_coord);
  free(face_vtx);
}

/*============================================================================
 * Main
 *============================================================================*/

int main(int argc, char *argv[])
{
  // Init
  PDM_MPI_Init(&argc, &argv);

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")

  // Malloc for cll
  int max_size = 10*2 + 3;
  int idx      = 0;
  Element **cll_storage = malloc(sizeof(Element *) * max_size);
  for (int i = 0; i < max_size; i++) {
    cll_storage[i] = malloc(sizeof(Element));
  }

  // Triangle: B->C->D->B cyclic linked list
  // Geogebra
  // double pt0[3] = {1.5, 1, 0};
  // double pt1[3] = {0, 0, 0};
  // double pt2[3] = {0.8, 0.3, 0.4};
  // inside
  double pt0[3] = {0.5, 0, 0};
  double pt2[3] = {0, 0.5, 0};
  double pt1[3] = {0, 0, 0.5};
  // XYZ
  // double pt0[3] = {1, 0, 0};
  // double pt1[3] = {0, 1, 0};
  // double pt2[3] = {0, 0, 1};


  Element *ptA = cll_storage[idx++];
  memcpy(ptA->coord, pt0, sizeof(double)*3);

  Element *ptB = cll_storage[idx++];
  memcpy(ptB->coord, pt1, sizeof(double)*3);
  ptA->next = ptB;

  Element *ptC = cll_storage[idx++];
  memcpy(ptC->coord, pt2, sizeof(double)*3);
  ptC->next = ptA;
  ptB->next = ptC;

  List *cll = malloc(sizeof(List));
  cll->head = ptA;

  // debug
  if (test_debug) {
    printf("Triangle:\n");
    _print_cll(cll);
  }

  // Determine A and outside (B before projection)
  List *outside = malloc(sizeof(List));
  _determine_A_outside(cll_storage, idx, cll, &outside);

  // debug
  if (test_debug) {
    printf("A: ");
    if (cll == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(cll);
    }

    printf("outside: ");
    if (outside == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(outside);
    }
  }

  // projection to go from outside to B
  List *B = outside;

  if (B != NULL) {
    Element *current = B->head;

    while (1) {

      double zB = 1 - current->coord[0] - current->coord[1];
      current->coord[2] = zB;

      if (current->next == B->head) break;
      current = current->next;
    }
  }

  // debug
  if (test_debug) {
    printf("B: ");
    if (B == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(B);
    }
  }

  // vtk
  if (test_debug) {
    if (cll != NULL) {
      char filename[999] = "A.vtk";
      _cll_to_polydata(cll, filename);
    }

    if (B != NULL) {
      char filename[999] = "B.vtk";
      _cll_to_polydata(B, filename);
    }
  }

  // compute columns
  double volumeA = -1;
  double volumeB = -1;

  if (cll == NULL) {
    volumeA = 0;
  }
  else {
    volumeA = _column_volume(cll);
  }

  if (B == NULL) {
    volumeB = 0;
  }
  else {
    volumeB = _column_volume(B);
  }

  // debug
  if (test_debug) {
    printf("volumeA = %f\n", volumeA);
    printf("volumeB = %f\n", volumeB);
  }

  // free
  for (int i = 0; i < max_size; i++) {
    free(cll_storage[i]);
  }
  free(cll_storage);
  if (cll != NULL) free(cll);
  if (outside != NULL) free(outside);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}

