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

static int dbg = 1;

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

typedef enum
{
  EXTREMUM,
  SHARP
} intersect_t;

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
  fflush(stdout);
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

// v2 of _determine_A_outside
static void
_determine_A_outside2
(
 Element  **cll_storage,
 int        idx,
 List     **cll,
 List     **outside
)
{

  // check if the triangle is inside the tetrahedron
  Element *current = (*cll)->head;
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
    if (current->next == (*cll)->head) break;
    current = current->next;
  }
  if (current->next == (*cll)->head) {

    if (dbg) {
      printf("triangle in tetrahedra\n");
    }

    // cll remains as is
    free(*outside);
    *outside = NULL;
    return;
  }

  // intersect with all 5 planes
  for (int i = 0; i < 5; i++) {

    if (dbg) {
      printf("Plane %d\n", i);
    }

    int intersect_idx = 0;

    // for linked lists
    current = (*cll)->head->next;

    intersect_t intersect_type[2] = {-1, -1};
    Element *in[2]  = {cll_storage[idx++], cll_storage[idx++]};
    Element *out[2] = {cll_storage[idx++], cll_storage[idx++]};

    Element *prev_in  = NULL;
    Element *prev_out = NULL;
    int idx_prev_in  = -1;
    int idx_prev_out = -1;

    // function value
    double fp = _plane_equation_function((*cll)->head->coord, i);
    double fc = _plane_equation_function(current->coord, i);
    double fn = 0;

    while (intersect_idx < 2) {

      if (dbg) {
        printf("--> SEGMENT (%f,%f,%f)-(%f,%f,%f)\n", current->coord[0], current->coord[1], current->coord[2], current->next->coord[0], current->next->coord[1], current->next->coord[2]);
      }

      // function value for next
      fn = _plane_equation_function(current->next->coord, i);

      // **********BEGIN********** //

      if (fc == 0) {

        if (fn == 0) {

          if (fp < 0) {

            if (dbg) {
              printf("polygon outside plane %d\n", i);
            }

            if (i == 4) {
              (*outside)->head = (*cll)->head;
              free(*cll);
              *cll = NULL;
              return;
            } // plane X+Y+Z=1

            else {
              free(*cll);
              free(*outside);
              *cll     = NULL;
              *outside = NULL;
              return;
            } // other planes

          } // polygon outside

          else if (fp > 0) {

            if (dbg) {
              printf("polygon inside plane %d\n", i);
            }

            if (i == 4) {
              // cll remains as is
              free(*outside);
              *outside = NULL;
              return;
            } // plane X+Y+Z=1

            else {
              break;
            } // other planes

          } // polygon inside

          else {

            if (dbg) {
              printf("polygon on plane %d\n", i);
            }

            if (i == 4) {
              // cll remains as is
              free(*outside);
              *outside = NULL;
              return;
            } // plane X+Y+Z=1

            else {
              free(*cll);
              free(*outside);
              *cll     = NULL;
              *outside = NULL;
              return;
            } // other planes
          } // polygon on (fp == 0)

        } // current and next on plane

        else {

          if (fp == 0) {

            if (fn < 0) {

              if (dbg) {
                printf("polygon outside plane %d\n", i);
              }

              if (i == 4) {
                (*outside)->head = (*cll)->head;
                free(*cll);
                *cll = NULL;
                return;
              } // plane X+Y+Z=1

              else {
               free(*cll);
               free(*outside);
               *cll     = NULL;
               *outside = NULL;
               return;
              } // other planes

            } // polygon outside

            else {

              if (dbg) {
                printf("polygon inside plane %d\n", i);
              }

              if (i == 4) {
                // cll remains as is
                free(*outside);
                *outside = NULL;
                return;
              } // plane X+Y+Z=1

              else {
                break;
              } // other planes

            } // polygon inside

          } // current and previous on plane

          else {

            if (fp*fn > 0) {

              if (fp > 0) {

                if (dbg) {
                  printf("polygon inside plane %d\n", i);
                }

                if (i == 4) {
                  // cll remains as is
                  free(*outside);
                  *outside = NULL;
                  return;
                } // plane X+Y+Z=1

                else {
                  break;
                } // other planes

              } // polygon inside

              else {

                if (dbg) {
                  printf("polygon outside plane %d\n", i);
                }

                if (i == 4) {
                  (*outside)->head = (*cll)->head;
                  free(*cll);
                  *cll = NULL;
                  return;
                } // plane X+Y+Z=1

                else {
                  free(*cll);
                  free(*outside);
                  *cll     = NULL;
                  *outside = NULL;
                  return;
                } // other planes

              } // polygon outside (fp < 0)

            } // next and previous on same side

            else {

              if (dbg) {
                printf("one extremum intersection with plane %d\n", i);
              }

              // coordinates of intersection
              memcpy(in[intersect_idx]->coord,  current->coord, sizeof(double) * 3);
              memcpy(out[intersect_idx]->coord, current->coord, sizeof(double) * 3);

              if (fp < 0 && fn > 0) {
                idx_prev_out             = intersect_idx;
                prev_out                 = current;
                in[intersect_idx]->next  = current->next;
                out[intersect_idx]->next = out[(intersect_idx+1)%2];
              } // from outside to inside

              else if (fp > 0 && fn < 0) {
                idx_prev_in              = intersect_idx;
                prev_in                  = current;
                out[intersect_idx]->next = current->next;
                in[intersect_idx]->next  = in[(intersect_idx+1)%2];
              } // from inside to outside

              intersect_type[intersect_idx++] = EXTREMUM;

            } // next and previous on different side (fn*fp < 0)

          } // only current on plane

        } // next not on plane

      } // current on plane

      else if (fc*fn < 0) {

        if (dbg) {
          printf("one sharp intersection with plane %d\n", i);
        }

        // coordinates of intersection
        double t = fc/(fc-fn);

        double pt[3] = {t*current->next->coord[0]+(1-t)*current->coord[0],
                        t*current->next->coord[1]+(1-t)*current->coord[1],
                        t*current->next->coord[2]+(1-t)*current->coord[2]};

        memcpy(in[intersect_idx]->coord,  pt, sizeof(double) * 3);
        memcpy(out[intersect_idx]->coord, pt, sizeof(double) * 3);

        if (fc < 0 && fn > 0) {
          idx_prev_out             = intersect_idx;
          prev_out                 = current;
          in[intersect_idx]->next  = current->next;
          out[intersect_idx]->next = out[(intersect_idx+1)%2];
        } // from outside to inside

        else if (fc > 0 && fn < 0) {
          idx_prev_in              = intersect_idx;
          prev_in                  = current;
          out[intersect_idx]->next = current->next;
          in[intersect_idx]->next  = in[(intersect_idx+1)%2];
        } // from inside to outside

        intersect_type[intersect_idx++] = SHARP;

      } // sharp intersection


      // ***********END*********** //

      // while loop end condition
      if (current->next == (*cll)->head->next) break;

      // update
      fp = fc;
      fc = fn;
      current = current->next;

    } // end while loop on cll

    // not intersection at all
    if (intersect_idx == 0) {
      printf("not intersection at all : fc = %f\n", fc);
      if (fc < 0) {
        if (i == 4) {
          (*outside)->head = (*cll)->head;
          free(*cll);
          *cll = NULL;
          return;
        } // plane X+Y+Z=1

        else {
          free(*cll);
          free(*outside);
          *cll     = NULL;
          *outside = NULL;
          return;
        } // other planes
      } // current outside
    } // 0 intersection

    // reconnect linked lists
    if (i == 4) {
      free(*outside);
      *outside = NULL;
    }
    if (intersect_idx == 2) {
      for (int j = 0; j < 2; j++) {
        if (j == idx_prev_in) {
          if (intersect_type[j] == SHARP) {
            prev_in->next  = in[idx_prev_in];
          }

          if (intersect_type[j] == EXTREMUM) {
            prev_in->next  = in[idx_prev_in]->next;
            in[idx_prev_in] = prev_in;
          }
        } // from in to out

        if (j == idx_prev_out) {
          if (intersect_type[j] == SHARP) {
            prev_out->next = out[idx_prev_out];
          }

          if (intersect_type[j] == EXTREMUM) {
            prev_out->next = out[idx_prev_out]->next;
            out[idx_prev_out] = prev_out;
          }
        } // from out to in
      } // loop on intersection points

      // update A and B
      (*cll)->head     = in[0];
      *outside = malloc(sizeof(List));
      (*outside)->head = out[0];
    } // 2 intersections

    if (dbg) {
       printf("cll at plane %d: ", i);
      _print_cll(*cll);
    }

  } // end for loop on planes

}

// --> if outside is NULL, it implies Vcolumn = 0
static void
_determine_A_outside
(
 Element  **cll_storage,
 int        idx,
 List     **cll,
 List     **outside
)
{

  // check if the triangle is inside the tetrahedron
  Element *current = (*cll)->head;
  int is_inside = 1;
  while (1) {
    // 0<= x, y, z <= 1 and 1-x-y-z >= 0
    int cond0 = current->coord[0] >= 0 && current->coord[0] <= 1;
    int cond1 = current->coord[1] >= 0 && current->coord[1] <= 1;
    int cond2 = current->coord[2] >= 0 && current->coord[2] <= 1;
    int cond3 = 1-current->coord[0]-current->coord[1]-current->coord[2] >= 0;
    // not totally in unit tetrahedron
    if (!(cond0 && cond1 && cond2 && cond3)) {
      is_inside = 0;
      break;
    }
    if (current->next == (*cll)->head) break;
    current = current->next;
  }
  if (is_inside) {

    if (dbg) {
      log_trace("triangle in tetrahedra\n");
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

    current = (*cll)->head;

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

        if (dbg) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) is on plane %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

        // polygon in plane
        if (f3 == 0) {

          if (dbg) {
            printf("polygon is in plane %d\n", i);
          }

          free(*cll);
          free(*outside);
          *cll      = NULL;
          *outside = NULL;
          return;
        }

        // only segment in plane
        else {

          // outside
          if (f3 <= 0) {

            if (dbg) {
              printf("polygon is outside %d\n", i);
            }

            free(*cll);
            free(*outside);
            *cll      = NULL;
            *outside = NULL;
            return;

          }

          // inside
          else {

            if (dbg) {
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

        if (dbg) {
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

        if (dbg) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) has no intersection with %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
        }

        if (f1 < 0) {

          if (dbg) {
            printf("is outside \n");
          }

        }

        if (f1 > 0) {

          if (dbg) {
            printf("is inside \n");
          }

        }

      }

      // to correcly loop over cll
      if (current->next == (*cll)->head) break;
      current = current->next;

    } // end while loop

    // no intersection at all
    if (intersect_idx == 0) {
       if (_plane_equation_function((*cll)->head->coord, i) < 0) {

          if (dbg) {
            printf("is totally outside \n");
          }

          free(*cll);
          free(*outside);
          *cll      = NULL;
          *outside  = NULL;
          return;
        }
    }

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
      (*cll)->head     = in[0];
      (*outside)->head = out[0];
    }

    if (dbg) {
       printf("cll at plane %d: ", i);
      _print_cll(*cll);
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








static void
_print_cll2
(
 int     head,
 double *coord,
 int    *next
 )
{
  if (head < 0) {
    log_trace("NULL\n");
    return;
  }

  int current = head;
  log_trace("head->");
  while (1) {
    log_trace("[(%f,%f,%f),%d]->",
              coord[3*current+0], coord[3*current+1], coord[3*current+2], next[current]);
    current = next[current];
    if (current == head) {
      break;
    }
  }
  log_trace("head\n");
}


static void
_determine_A_outside3
(
 int    *cll_head_in,
 int    *cll_head_out,
 double *cll_coord,
 int    *cll_next
 )
{
  // int dbg = 1;

  /* Check if the initial polygon (triangle) is inside the unit tetrahedron */
  int current = *cll_head_in;
  int inside = 1;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (cll_coord[3*current+j] < 0 || cll_coord[3*current+j] > 1) {
        log_trace("i = %d, outside %d\n", i, j);
        inside = 0;
        break;
      }
    }

    if (!inside) {
      break;
    }

    log_trace("i = %d : %f\n", i, 1 - cll_coord[3*current] - cll_coord[3*current+1] - cll_coord[3*current+2]);
    if (1 - cll_coord[3*current] - cll_coord[3*current+1] - cll_coord[3*current+2] < 0) {
      inside = 0;
      break;
    }

    current = cll_next[current];
  }

  if (inside) {
    if (dbg) {
      log_trace("triangle in tetrahedra\n");
    }

    // cll remains as is
    return;
  }



  int prev[2] = {-1, -1};

  intersect_t intersect_type[2] = {-1, -1};

  /* Intersect/clip with all 5 planes */
  int idx = 3;
  for (int iplane = 0; iplane < 5; iplane++) {

    if (*cll_head_in < 0) {
      return;
    }

    if (dbg) {
      log_trace(">>>\n");
      log_trace("Plane %d\n", iplane);
    }

    int idx0 = idx;

    int intersect_idx = 0;

    current = cll_next[*cll_head_in];

    int start = current;


    // function value
    double *cc = cll_coord + 3*current;
    double fp = _plane_equation_function(&cll_coord[3*(*cll_head_in)], iplane);
    double fc = _plane_equation_function(&cll_coord[3*current],        iplane);
    double fn = 0;

    while (intersect_idx < 2) {

      int next = cll_next[current];
      double *cn = cll_coord + 3*next;
      // function value for next
      fn = _plane_equation_function(cn, iplane);
      if (dbg) {
        log_trace("--> SEGMENT (%f,%f,%f)-(%f,%f,%f)\n",
                  cc[0], cc[1], cc[2],
                  cn[0], cn[1], cn[2]);
        log_trace("fp/fc/fn = %f / %f / %f\n", fp, fc, fn);
      }

      // **********BEGIN********** //
      if (fc == 0) {

        if (fn == 0) {

          if (fp <= 0) {
            // Segment cn on plane, rest outside or on plane
            if (dbg) {
              log_trace("Segment cn on plane, rest outside or on plane\n");
            }
            *cll_head_out = *cll_head_in;
            *cll_head_in  = -1;
            break;
          }
          else {
            // Segment cn on plane, rest inside
            if (dbg) {
              log_trace("Segment cn on plane, rest inside\n");
            }
            *cll_head_out = -1;
            break;
          }

        } // End if fc == 0, fn == 0

        else {

          if (fp == 0) {

            if (fn < 0) {
              // Segment pc on plane, rest outside
              if (dbg) {
                log_trace("Segment pc on plane, rest outside\n");
              }
              *cll_head_out = *cll_head_in;
              *cll_head_in  = -1;
              break;
            }
            else {
              // Segment pc on plane, rest inside
              if (dbg) {
                log_trace("Segment pc on plane, rest inside\n");
              }
              *cll_head_out = -1;
              break;
            }

          } // End if fc == 0, fp == 0

          else {

            if (fp*fn > 0) {
              if (fp > 0) {
                // c on plane, rest inside
                if (dbg) {
                  log_trace("c on plane, rest inside\n");
                }
                *cll_head_out = -1;
                break;
              }
              else {
                // c on plane, rest outside
                if (dbg) {
                  log_trace("c on plane, rest outside\n");
                }
                *cll_head_out = *cll_head_in;
                *cll_head_in  = -1;
                break;
              }
            }

            else {
              // c on plane, p and n in opposite sides
              if (dbg) {
                log_trace("c on plane, p and n in opposite sides\n");
              }
              prev[intersect_idx] = current;
              if (intersect_idx == 0) {
                if (fn > 0) {
                  // entering
                  *cll_head_out = current;
                  *cll_head_in  = next;
                }
                else {
                  // exiting
                  *cll_head_in  = current;
                  *cll_head_out = next;
                }
                cll_next[idx] = idx0 + 3;
              }
              else {
                cll_next[idx] = idx0 + 1;
              }
              cll_next[idx+1] = next;

              for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 3; j++) {
                  cll_coord[3*idx+j] = cc[j];
                }
                idx++;
              }
              intersect_type[intersect_idx++] = EXTREMUM;
            }

          } // End if fc == 0, fp != 0

        } // End if fc == 0, fn != 0

      } // End if fc == 0

      else if (fc * fn < 0) {
        // Regular intersection between c and n
        if (dbg) {
          log_trace("Regular intersection between c and n\n");
        }

        prev[intersect_idx] = current;
        if (intersect_idx == 0) {
          if (fn > 0) {
            // entering
            *cll_head_out = current;
            *cll_head_in  = next;
          }
          else {
            // exiting
            *cll_head_in  = current;
            *cll_head_out = next;
          }
          cll_next[idx] = idx0 + 3;
        }
        else {
          cll_next[idx] = idx0 + 1;
        }
        cll_next[idx+1] = next;

        double t = fc/(fc - fn);
        for (int i = 0; i < 2; i++) {
          for (int j = 0; j < 3; j++) {
            cll_coord[3*idx+j] = (1-t)*cc[j] + t*cn[j];
          }
          idx++;
        }

        intersect_type[intersect_idx++] = SHARP;
      }

      // ***********END*********** //

      // while loop end condition
      current = cll_next[current];
      if (current == start) break;

      // update
      fp = fc;
      fc = fn;
      cc = cll_coord + 3*current;

    } // end while loop on cll

    if (dbg) {
      log_trace("  intersect_idx = %d\n", intersect_idx);
    }


    // no intersection at all
    if (intersect_idx == 0) {
      if (fc < 0) {
        // Polygon outside
        if (dbg) {
          log_trace("Polygon outside (no intersection at all)\n");
        }
        *cll_head_out = *cll_head_in;
        *cll_head_in  = -1;
        return;
      }
      else {
        // Polygon inside
        if (dbg) {
          log_trace("Polygon inside (no intersection at all)\n");
        }
        if (iplane < 4) {
          *cll_head_out = -1;
        }
      }
    }

    else if (intersect_idx == 2) {

      for (int i = 0; i < 2; i++) {
        if (intersect_type[i] == EXTREMUM) {
          cll_next[prev[i]] = cll_next[idx0 + 2*i];
        }
        else {
          cll_next[prev[i]] = idx0 + 2*i;
        }
      }

    }

    else {
      // error?!
    }

    if (dbg) {
       log_trace("IN at plane %d: ", iplane);
      _print_cll2(*cll_head_in, cll_coord, cll_next);
      log_trace("OUT at plane %d: ", iplane);
      _print_cll2(*cll_head_out, cll_coord, cll_next);
      log_trace("<<<\n");
    }


  } // End of loop on planes


}


static double
_column_volume2
(
 int     head,
 double *coord,
 int    *next
)
{
  // int dbg = 1;

  double volume = 0;
  if (head < 0) {
    return volume;
  }

  int current = next[head];

  while (next[current] != head) {

    double prism_volume = _prism_volume(&coord[3*head],
                                        &coord[3*current],
                                        &coord[3*next[current]]);
    if (dbg) {
      log_trace("            volumeK += %f\n", prism_volume);
    }
    volume += prism_volume;

    current = next[current];
  }

  return volume;
}

static void
_cll_to_polydata2
(
 int      head,
 double  *coord,
 int     *next,
 const char *filename
 // double **vtx_coord,
 // int    **face_vtx
 )
{
  int size_min = 3;

  int n_vtx = 0;
  double *vtx_coord = malloc(sizeof(double) * size_min * 3);

  int current = head;

  while (1) {

    if (n_vtx > size_min -1) {
      size_min *= 2;
      vtx_coord = realloc(vtx_coord, sizeof(double) * size_min * 3);
    }

    memcpy(vtx_coord + n_vtx * 3, coord + 3*current, sizeof(double) * 3);
    n_vtx++;

    if (next[current] == head) break;
    current = next[current];
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

  // return n_vtx;
}


static double
_reference_volume_compute2
(
 double triaB_coord[9]
)
{
  // int dbg = 1;

  double coord[23*3];
  int    next[23];

  /* Initial polygon == triangle */
  for (int i = 0; i < 3; i++) {
    memcpy(&coord[3*i], &triaB_coord[3*i], sizeof(double)*3);
    next[i] = (i+1)%3;
  }


  int head_in  = 0;
  int head_out = -1;

  // Determine A and B before projection
  _determine_A_outside3(&head_in,
                        &head_out,
                        coord,
                        next);

  // debug
  if (dbg && 1) {
    log_trace("A: ");
    _print_cll2(head_in,
                coord,
                next);

    log_trace("B: ");
    _print_cll2(head_out,
                coord,
                next);
  }

  // Projection of B
  if (head_out >= 0) {
    int current = head_out;

    while (1) {
      double zB = 1 - coord[3*current+0] - coord[3*current+1];
      coord[3*current+2] = zB;

      current = next[current];
      if (current == head_out) break;
    }
  }

  if (dbg && 1) {
    if (head_in >= 0) {
      char filename[999] = "A2.vtk";
      _cll_to_polydata2(head_in,  coord, next, filename);
    }

    if (head_out >= 0) {
      char filename[999] = "B2.vtk";
      _cll_to_polydata2(head_out, coord, next, filename);
    }
  }

  // Determine volume
  double volumeA = _column_volume2(head_in,  coord, next);
  double volumeB = _column_volume2(head_out, coord, next);

  return volumeA + volumeB;
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

  double vol = _reference_volume_compute2(triaB_coord);



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
  if (dbg) {
    printf("Triangle:\n");
    _print_cll(cll);

    char filename[999] = "T.vtk";
    _cll_to_polydata(cll, filename);
  }

  // Determine A and outside (B before projection)
  List *outside = malloc(sizeof(List));
  _determine_A_outside2(cll_storage, idx, &cll, &outside);

  // debug
  if (dbg) {
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
  if (dbg) {
    printf("B: ");
    if (B == NULL) {
      printf("NULL\n");
    } else {
      _print_cll(B);
    }
  }

  // vtk
  if (dbg) {
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
  if (dbg) {
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

