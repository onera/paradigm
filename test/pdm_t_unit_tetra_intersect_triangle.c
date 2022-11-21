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
 * Type definitions
 *============================================================================*/

typedef struct Element Element;
struct Element
{
  double coord[3];
  Element *next;

  int intersect;
};

typedef struct List List;
struct List
{
    Element *head;
};

/*============================================================================
 * Private functions
 *============================================================================*/

// for linked list

static void
_print_linked_list
(
 List *ll
)
{
  if (ll == NULL)
  {
    exit(EXIT_FAILURE);
  }

  Element *current = ll->head;
  printf("head->");
  while (current != NULL) {
    printf("[(%f,%f,%f),%d,%p]->", current->coord[0], current->coord[1], current->coord[2], current->intersect, (void *) current->next);
    current = current->next;
  }
  printf("NULL\n");
}

static void
_free
(
 List *ll
)
{
  if (ll->head != NULL) {
    Element *current = ll->head;
    Element *prev = NULL;
    while (current != NULL) {
      prev = current;
      current = current->next;
      free(prev);
    }
  }
  free(ll);
}

static void
_split_list
(
 List *out,
 List *in
)
{
  Element *current_out = out->head;
  while (current_out->intersect == 0) {
    current_out = current_out->next;
  }

  // out intersection found, tag as no more intersection
  current_out->intersect = 0;

  // create copy for in list
  Element *intersect_copy_in = malloc(sizeof(Element));
  memcpy(intersect_copy_in->coord, current_out->coord, sizeof(double)*3);
  intersect_copy_in->next = current_out->next;
  intersect_copy_in->intersect = 0;
  in->head = intersect_copy_in;

  Element *current_in = in->head;
  while (current_in->intersect == 0) {
    current_in = current_in->next;
  }

  // in intersection found, tag as no more intersection
  current_in->intersect = 0;

  // create tail for in
  Element *in_tail = malloc(sizeof(Element));
  memcpy(in_tail->coord, in->head->coord, sizeof(double)*3);
  in_tail->next = NULL;
  in_tail->intersect = 0;

  // create copy for out list
  Element *intersect_copy_out = malloc(sizeof(Element));
  memcpy(intersect_copy_out->coord, current_in->coord, sizeof(double)*3);
  intersect_copy_out->next = current_in->next;
  intersect_copy_out->intersect = 0;

  // set next to tail
  current_in->next = in_tail;

  // set out next
  current_out->next = intersect_copy_out;
}

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

/*============================================================================
 * Main
 *============================================================================*/

int main(int argc, char *argv[])
{
  // Init
  PDM_MPI_Init(&argc, &argv);

  // // create linked list
  // double tab[3] = {0, 0, 0};

  // Element *head = malloc(sizeof(Element));
  // memcpy(head->coord, tab, sizeof(double)*3);
  // head->intersect = 0;

  // tab[0] = 1;

  // Element *i1 = malloc(sizeof(Element));
  // memcpy(i1->coord, tab, sizeof(double)*3);
  // i1->intersect = 1;
  // head->next = i1;

  // tab[0] = 2;

  // Element *A = malloc(sizeof(Element));
  // memcpy(A->coord, tab, sizeof(double)*3);
  // A->intersect = 0;
  // i1->next = A;

  // tab[0] = 3;

  // Element *i2 = malloc(sizeof(Element));
  // memcpy(i2->coord, tab, sizeof(double)*3);
  // i2->intersect = 1;
  // A->next = i2;

  // tab[0] = 0;

  // Element *head1 = malloc(sizeof(Element));
  // memcpy(head1->coord, tab, sizeof(double)*3);
  // head1->intersect = 0;
  // head1->next = NULL;
  // i2->next = head1;

  // List *ll = malloc(sizeof(List));
  // ll->head = head;

  // // print it
  // _print_linked_list(ll);

  // // split linked list
  // List *in = malloc(sizeof(List));
  // _split_list(ll, in);

  // printf("after split:\n");
  // _print_linked_list(ll);
  // _print_linked_list(in);

  // // free
  // _free(ll);
  // _free(in);

  // Triangle
  // Geogebra
  double pt0[3] = {1.5, 1, 0};
  double pt1[3] = {0, 0, 0};
  double pt2[3] = {0.8, 0.3, 0.4};
  // inside
  // double pt0[3] = {0.5, 0, 0};
  // double pt1[3] = {0, 0.5, 0};
  // double pt2[3] = {0, 0, 0.5};

  Element *A0 = malloc(sizeof(Element));
  memcpy(A0->coord, pt0, sizeof(double)*3);
  A0->intersect = 0;

  Element *B = malloc(sizeof(Element));
  memcpy(B->coord, pt1, sizeof(double)*3);
  B->intersect = 0;
  A0->next = B;

  Element *C = malloc(sizeof(Element));
  memcpy(C->coord, pt2, sizeof(double)*3);
  C->intersect = 0;
  B->next = C;

  Element *A1 = malloc(sizeof(Element));
  memcpy(A1->coord, pt0, sizeof(double)*3);
  A1->intersect = 0;
  A1->next = NULL;
  C->next = A1;

  List *ll = malloc(sizeof(List));
  ll->head = A0;

  printf("Triangle:\n");
  _print_linked_list(ll);

  // Determine A and B
  List *A = malloc(sizeof(List));
  List *outside = malloc(sizeof(List));

  for (int i = 0; i < 1; i++) {

    // check if inside tetrahedron, then ll == A and break
    Element *current = ll->head;
    while (current != NULL) {
      int cond0 = current->coord[0] >= 0 && current->coord[0] <= 1;
      int cond1 = current->coord[1] >= 0 && current->coord[1] <= 1;
      int cond2 = current->coord[2] >= 0 && current->coord[2] <= 1;
      int cond3 = 1-current->coord[0]-current->coord[1]-current->coord[2] >= 0;
      // not totally in unit tetrahedron
      if (!(cond0 && cond1 && cond2 && cond3)) {
        break;
      }
      current = current->next;
    }
    // is totaly inside
    if (current == NULL) {
      printf("at step %d is totally inside the unit tetrahedron\n", i);
      A = ll;
      break;
    }

    // check weather head is in or out
    int cond_head_out = -1;
    switch (i) {
    // OYZ
    case 0:
      cond_head_out = ll->head->coord[0] < 0;
      break;
    // X+Y=1
    case 1:
      cond_head_out = (ll->head->coord[0] + ll->head->coord[1]) > 1;
      break;
    // OZX
    case 2:
      cond_head_out = ll->head->coord[2] < 0;
      break;
    // OXY
    case 3:
      cond_head_out = ll->head->coord[3] < 0;
      break;
    // X+Y+Z=1
    case 4:
      cond_head_out = (ll->head->coord[0] + ll->head->coord[1] + ll->head->coord[2]) > 1;
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "Only 5 planes\n");
      break;
    }

    // find the maximum two intersections for each plane
    int count_intersect = 0;
    current = ll->head;
    // --> loop over the polyhedra segments
    while (current->next != NULL && count_intersect < 2) {
      double *coord1 = current->coord;
      double *coord2 = current->next->coord;

      switch (i) {
      // OYZ
      case 0:
        // on plane ?
        // TO DO: make robust
        if (coord1[0] == 0 && coord2[0] == 0) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) is on plane %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);

          // polygon in plane ?
          if (current->next->next->coord[0] == 0) {
            printf("polygon is in plane %d\n", i);
            break; // TO DO: output that Vcolumn is 0
          } else {
            Element *i1 = malloc(sizeof(Element));
            memcpy(i1->coord, current->coord, sizeof(double) * 3);
            i1->intersect = 1;

            Element *i2 = malloc(sizeof(Element));
            memcpy(i2->coord, current->next->coord, sizeof(double) * 3);
            i2->intersect = 1;

            i2->next      = current->next;
            i1->next      = i2;
            current->next = i1;

            count_intersect += 2;

            current = i2->next;
          }
        }

        // has intersection ?
        else if ((coord1[0] >= 0 && coord2[0] <= 0) || (coord1[0] <= 0 && coord2[0] >= 0)) {
          printf("segment (%f,%f,%f)-(%f,%f,%f) intersects plane %d\n", coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2], i);
          double t = coord2[0]/(coord2[0]-coord1[0]);
          Element *i1 = malloc(sizeof(Element));
          i1->coord[0] = t*coord1[0]+(1-t)*coord2[0];
          i1->coord[1] = t*coord1[1]+(1-t)*coord2[1];
          i1->coord[2] = t*coord1[2]+(1-t)*coord2[2];
          i1->next  = current->next;
          i1->intersect = 1;
          current->next = i1;

          count_intersect += 1;

          current = i1->next;
        }

        else {
          current = current->next;
        }
        break;
      // X+Y=1
      case 1:
        break;
      // OZX
      case 2:
        break;
      // OXY
      case 3:
        break;
      // X+Y+Z=1
      case 4:
        break;
      default:
        PDM_error(__FILE__, __LINE__, 0, "Only 5 planes\n");
        break;
      }
    }

    // --> split and keep depending on head in/out
    // TO DO: gestion des fuites mémoires ici !!!
    List *in = NULL;
    List *out = NULL;

    if (cond_head_out) {
      in = malloc(sizeof(List));
      out = ll;
      _split_list(out, in);
    } else {
      out = malloc(sizeof(List));
      in = ll;
      _split_list(in, out);
    }

    printf("After split:\n");
    _print_linked_list(in);
    _print_linked_list(out);

    if (i < 4) {
      _free(out);
      ll = in;
    } else { // if X+Y+Z=1
      A = in;
      outside = out;
    }
  }

  // free
  _free(ll);
  free(A);
  free(outside);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}
