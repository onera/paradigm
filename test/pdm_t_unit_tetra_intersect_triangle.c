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

  // Planes
  double **normals = malloc(sizeof(double *) * 5);
  // OYZ
  double n0[3] = {1, 0, 0};
  normals[0] = n0;
  // X+Y = 1
  double n1[3] = {-1, -1, 0};
  normals[1] = n1;
  // OZX
  double n2[3] = {0, 1, 0};
  normals[2] = n2;
  // OXY
  double n3[3] = {0, 0, 1};
  normals[3] = n3;
  // X+Y+Z = 1
  double n4[4] = {-1, -1, -1};
  normals[4] = n4;
  double **points  = malloc(sizeof(double *) * 5);
  double z[3] = {0, 0, 1};
  double x[3] = {1, 0, 0};
  // OYZ
  points[0] = z; // Z
  // X+Y = 1
  points[1] = x; // X
  // OZX
  points[2] = z; // Z
  // OXY
  points[3] = x; // X
  // X+Y+Z = 1
  points[4] = z; // Z

  // Determine A and B
  List *A = malloc(sizeof(List));
  List *outside = malloc(sizeof(List));

  for (int i = 0; i < 5; i++) {
    double *_n = normals[i];
    double *_pt = points[i];

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

    // find the maximum two intersections for each plane
  }

  // free
  _free(ll);
  free(normals);
  free(points);
  free(A);
  free(outside);

  // Finalize
  PDM_MPI_Finalize();

  return 0;
}
