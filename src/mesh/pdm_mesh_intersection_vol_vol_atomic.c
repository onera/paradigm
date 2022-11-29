#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm_mesh_intersection_vol_vol_atomic.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm.h"
#include "pdm_priv.h"

/*============================================================================
 * Global variables
 *============================================================================*/

static int dbg = 1;

/*============================================================================
 * Type definitions
 *============================================================================*/

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


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
/*=============================================================================
 * Static function definitions
 *============================================================================*/
const double epsilon = 1.e-15;

static inline int _is_zero
(
 const double x
 )
{
  return PDM_ABS(x) < epsilon;
}


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
    log_error("Only 5 planes are considered\n");
    return 0;
  }
}




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
  log_trace("head->");
  while (1) {
    log_trace("[(%f,%f,%f),%p]->", current->coord[0], current->coord[1], current->coord[2], (void *) current->next);
    if (current->next == cll->head) break;
    current = current->next;
  }
  log_trace("head\n");
}



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
      log_trace("triangle in tetrahedra\n");
    }

    // cll remains as is
    free(*outside);
    *outside = NULL;
    return;
  }

  // intersect with all 5 planes
  for (int i = 0; i < 5; i++) {

    if (dbg) {
      log_trace("Plane %d\n", i);
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
        log_trace("--> SEGMENT (%f,%f,%f)-(%f,%f,%f)\n", current->coord[0], current->coord[1], current->coord[2], current->next->coord[0], current->next->coord[1], current->next->coord[2]);
      }

      // function value for next
      fn = _plane_equation_function(current->next->coord, i);

      // **********BEGIN********** //

      if (fc == 0) {

        if (fn == 0) {

          if (fp < 0) {

            if (dbg) {
              log_trace("polygon outside plane %d\n", i);
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
              log_trace("polygon inside plane %d\n", i);
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
              log_trace("polygon on plane %d\n", i);
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
                log_trace("polygon outside plane %d\n", i);
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
                log_trace("polygon inside plane %d\n", i);
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
                  log_trace("polygon inside plane %d\n", i);
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
                  log_trace("polygon outside plane %d\n", i);
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
                log_trace("one extremum intersection with plane %d\n", i);
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
          log_trace("one sharp intersection with plane %d\n", i);
        }

        // coordinates of intersection
        double t = fc/(fc-fn);
        t = PDM_MIN(1-epsilon, t);

        double pt[3] = {t*current->next->coord[0]+(1-t)*current->coord[0],
                        t*current->next->coord[1]+(1-t)*current->coord[1],
                        t*current->next->coord[2]+(1-t)*current->coord[2]};

        if (dbg) {
          log_trace("t = %20.16f\n", t);
          log_trace("pt = %20.16f %20.16f %20.16f \n", pt[0], pt[1], pt[2]);
        }

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
      if (dbg) {
        log_trace("not intersection at all : fc = %f\n", fc);
      }
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
      (*cll)->head = in[0];
      free(*outside);
      *outside = malloc(sizeof(List));
      (*outside)->head = out[0];
    } // 2 intersections

    if (dbg) {
       log_trace("cll at plane %d: ", i);
      _print_cll(*cll);
    }

  } // end for loop on planes

}



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
  const char *filename
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

  if (head < 0) {
    return;
  }

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




static void
_determine_A_outside2
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

    *cll_head_out = -1;
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
        log_trace("fp/fc/fn = %20.16f / %20.16f / %20.16f\n", fp, fc, fn);
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
        t = PDM_MIN(1-epsilon, t);
        if (dbg) {
          log_trace("t = %20.16f\n", t);
          log_trace("fc - fn = %20.16f\n", fc - fn);
          log_trace("pt = %20.16f %20.16f %20.16f \n",
                    (1-t)*cc[0] + t*cn[0],
                    (1-t)*cc[1] + t*cn[1],
                    (1-t)*cc[2] + t*cn[2]);
        }
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
        if (iplane == 4) {
          *cll_head_out = *cll_head_in;
        }
        else {
          *cll_head_out = -1;
        }
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

      char filename[999];

      sprintf(filename, "in_%d.vtk", iplane);
      _cll_to_polydata2(*cll_head_in, cll_coord, cll_next, filename);

      sprintf(filename, "out_%d.vtk", iplane);
      _cll_to_polydata2(*cll_head_out, cll_coord, cll_next, filename);
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


PDM_GCC_SUPPRESS_WARNING_POP

/*=============================================================================
 * Public function definitions
 *============================================================================*/


double
PDM_mesh_intersection_vol_vol_atomic_compute
(
 double    triaB_coord[9],
 double  **vtx_coordA,
 int      *n_vtxA,
 int     **face_vtxA,
 int      *n_faceA,
 double  **vtx_coordB,
 int      *n_vtxB,
 int     **face_vtxB,
 int      *n_faceB
)
{
  PDM_UNUSED(vtx_coordA);
  PDM_UNUSED(n_vtxA);
  PDM_UNUSED(face_vtxA);
  PDM_UNUSED(n_faceA);
  PDM_UNUSED(vtx_coordB);
  PDM_UNUSED(n_vtxB);
  PDM_UNUSED(face_vtxB);
  PDM_UNUSED(n_faceB);

  // malloc List structure
  List *A = malloc(sizeof(List));
  List *B = malloc(sizeof(List));

  // int dbg = 0;
  int idx = 0;

  // Triangle cll
  int max_size = 10*2 + 3;
  Element **cll_storage = malloc(sizeof(Element *) * max_size);
  for (int i = 0; i < max_size; i++) {
    cll_storage[i] = malloc(sizeof(Element));
  }
  Element *pt0 = cll_storage[idx++];
  memcpy(pt0->coord, triaB_coord,     sizeof(double)*3);

  Element *pt1 = cll_storage[idx++];
  memcpy(pt1->coord, triaB_coord + 3, sizeof(double)*3);
  pt0->next = pt1;

  Element *pt2 = cll_storage[idx++];
  memcpy(pt2->coord, triaB_coord + 6, sizeof(double)*3);
  pt2->next = pt0;
  pt1->next = pt2;

  A->head = pt0;

  // debug
  if (dbg && 1) {
    log_trace("Triangle:\n");
    _print_cll(A);
  }

  // Determine A and B before projection
  _determine_A_outside(cll_storage, idx, &A, &B);

  // debug
  if (dbg && 1) {
    log_trace("A: ");
    if (A == NULL) {
      log_trace("NULL\n");
    } else {
      _print_cll(A);
    }

    log_trace("B: ");
    if (B == NULL) {
      log_trace("NULL\n");
    } else {
      _print_cll(B);
    }

    if (B != NULL) {
      _cll_to_polydata(B, "B_noproj.vtk");
    }
  }

  // Projection of B
  if (B != NULL) {
    Element *current = B->head;

    while (1) {

      double zB = 1 - current->coord[0] - current->coord[1];
      current->coord[2] = zB;

      if (current->next == B->head) break;
      current = current->next;
    }
  }

  // vtk
  if (dbg && 1) {

    if (A != NULL) {
      _cll_to_polydata(A, "A.vtk");
    }

    if (B != NULL) {
      _cll_to_polydata(B, "B.vtk");
    }
  }

  // Determine volume
  double volumeA = -1;
  double volumeB = -1;

  if (A == NULL) {
    volumeA = 0;
  }
  else {
    volumeA = _column_volume(A);
  }

  if (B == NULL) {
    volumeB = 0;
  }
  else {
    volumeB = _column_volume(B);
  }

  // free List structure
  if (A != NULL) free(A);
  if (B != NULL) free(B);
  for (int i = 0; i < max_size; i++) {
    free(cll_storage[i]);
  }
  free(cll_storage);

  return volumeA + volumeB;
}


double
PDM_mesh_intersection_vol_vol_atomic_compute2
(
 double triaB_coord[9]
 )
{
  // int dbg = 1;
  if (dbg) {
    log_trace("initial triangle :\n");
    for (int i = 0; i < 3; i++) {
      log_trace("%20.16f, %20.16f, %20.16f\n",
                triaB_coord[3*i], triaB_coord[3*i+1], triaB_coord[3*i+2]);
    }
  }

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
  _determine_A_outside2(&head_in,
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

    if (head_out >= 0) {
      _cll_to_polydata2(head_out, coord, next, "B_noproj2.vtk");
    }
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
      _cll_to_polydata2(head_in,  coord, next, "A2.vtk");
    }

    if (head_out >= 0) {
      _cll_to_polydata2(head_out, coord, next, "B2.vtk");
    }
  }

  // Determine volume
  double volumeA = _column_volume2(head_in,  coord, next);
  double volumeB = _column_volume2(head_out, coord, next);

  return volumeA + volumeB;
}
