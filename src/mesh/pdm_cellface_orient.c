
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_array.h"
#include "pdm_geom_elem.h"
#include "pdm_cellface_orient.h"
#include "pdm_hash_tab.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/


enum {false, true};

typedef enum {
  FACE_UNPROCESSED,
  FACE_IN_STACK,
  FACE_UNCHANGED_CYCLE,
  FACE_CHANGED_CYCLE
} _face_state_t;


typedef enum {
  CELL_UNPROCESSED,
  CELL_IN_STACK,
  CELL_COMPLETED
} _cell_state_t;

/*============================================================================
 * Static global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/*=============================================================================
 * Public function definitions
 *============================================================================*/



/**
 * \brief Orient cell->face connectivity
 *
 * At the output of th function, a face number in \ref cell_face is positive
 * if surface normal is inside the cell, negative otherwise. \ref face_cell is
 * oriented in the same way
 *
 * \param [in]      n_cell         Number of cells
 * \param [in]      n_face         Number of faces
 * \param [in]      n_vtx          Number of vertices
 * \param [in]      coords         Vertices coordinates
 * \param [in]      cell_face_idx  Cell to face connectivity index (size = \ref n_cell + 1)
 * \param [in, out] cell_face      Cell to face connectivity (size = cell_face_idx[n_cell])
 * \param [in, out] face_cell      Face to cell connectivity (size = 2 * \ref n_face) or NULL
 * \param [in]      face_vtx_idx   Face to vertex connectivity index (size = \ref n_face + 1)
 * \param [in]      face_vtx       Face to vertex connectivity (size = face_vtx_idx[n_face])
 *
 */

void
PDM_cellface_orient
(
const int      n_cell,
const int      n_face,
const int      n_vtx,
const double  *coords,
const int     *cell_face_idx,
int           *cell_face,
int           *face_cell,
const int     *face_vtx_idx,
const int     *face_vtx
)
{
  if (n_cell == 0) {
    return;
  }

  PDM_timer_t *t1 = PDM_timer_create();
  PDM_timer_resume(t1);

  int *_face_cell = NULL;

  if (face_cell != NULL) {
    _face_cell = face_cell;
  }
  else {
    _face_cell = PDM_array_zeros_int(2*n_face);

    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = PDM_ABS (cell_face[j]);
        int i_face = 2 * (face - 1);
        if (_face_cell[i_face] == 0) {
          _face_cell[i_face] = i + 1;
        }
        else {
          _face_cell[i_face+1] = i + 1;
        }
      }
    }
  }

  if (1 == 0) {
    printf("_face_cell : ");
    for (int i = 0; i < n_face; i++) {
      printf("%d : %d %d\n", i+1, _face_cell[2*i], _face_cell[2*i+1]);
    }
    printf("\n");
  }

  int *orientedface_cell = PDM_array_zeros_int(2*n_face);


  int key_max = 2 * n_vtx;
  PDM_hash_tab_t *hash_orient = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, &key_max);

  int max_n_poly_face = -1;
  int max_edges = -1;

  for (int ipoly = 0; ipoly < n_cell; ipoly++) {
    const int poly_idx   = cell_face_idx[ipoly];
    const int npoly_face = cell_face_idx[ipoly + 1] - poly_idx;
    max_n_poly_face = PDM_MAX (max_n_poly_face, npoly_face);
    int n_edge_cell = 0;
    for (int i = poly_idx; i < poly_idx + npoly_face; i++) {
      int face = PDM_ABS(cell_face[i]) - 1;
      n_edge_cell += face_vtx_idx[face+1] - face_vtx_idx[face];

    }
    max_edges = PDM_MAX (max_edges, n_edge_cell);
  }

  int n_processed_face = 0;
  int n_stack_cell     = -1;
  int *stack_face      = NULL;
  int *stack_cell      = NULL;
  int *processed_face  = NULL;
  PDM_malloc(stack_face    , max_n_poly_face, int);
  PDM_malloc(stack_cell    , n_cell         , int);
  PDM_malloc(processed_face, max_n_poly_face, int);

  int *tag_cell = PDM_array_const_int(n_cell         , CELL_UNPROCESSED);
  int *tag_face = PDM_array_const_int(max_n_poly_face, FACE_UNPROCESSED);

  tag_cell[0] = CELL_COMPLETED;

  int n_edges = 0;
  const int n_data_edge = 3;
  int *edges = NULL;
  PDM_malloc(edges, max_edges * n_data_edge, int);

 /*
  * Orient the first cell of the first component
  * --------------------------------------------
  *
  * As the oriented volume is positive, the face normals
  * are outside of the element
  *
  */

  int first_cell_comp = 0;

  while (first_cell_comp != -1) {

    int     is_oriented = 0;
    int     n_polyhedra = 1;
    double  volume[3];
    double  center[3];

    int *_cell_face_idx = (int *) cell_face_idx + first_cell_comp;

    PDM_geom_elem_polyhedra_properties (is_oriented,
                                        n_polyhedra,
                                        n_face,
                                        face_vtx_idx,
                                        face_vtx,
                                        _cell_face_idx,
                                        cell_face,
                                        n_vtx,
                                        coords,
                                        volume,
                                        center,
                                        NULL,
                                        NULL);

   /*
    * Initialize an oriented face cell
    * --------------------------------
    *
    * left cell : normal is inside the cell
    * right cell : normal is outside the cell
    *
    * The orientation of the first cell is taking into account
    *
    */

    tag_cell[first_cell_comp] = CELL_COMPLETED;

    for (int i = cell_face_idx[first_cell_comp]; i < cell_face_idx[first_cell_comp+1]; i++) {
      int i_face = cell_face[i];

      if (i_face > 0) {
        orientedface_cell[2 * (i_face - 1)] = first_cell_comp + 1;
      }
      else {
        orientedface_cell[2 * (PDM_ABS (i_face) - 1) + 1] = first_cell_comp + 1;
      }
    }


   /*
    * Other cells are oriented from the faces of the first cell
    * ---------------------------------------------------------
    *
    */

    /* Add neighbours of the first cell in the stack */

    for (int i = cell_face_idx[first_cell_comp]; i < cell_face_idx[first_cell_comp+1]; i++) {
      int i_face = 2 * (PDM_ABS (cell_face[i]) - 1);

      // if (_face_cell[i_face] == 1) {
      if (_face_cell[i_face] == first_cell_comp+1) {
        if (_face_cell[i_face + 1] > 0) {
          int cell = PDM_ABS (_face_cell[i_face + 1]);
          if (tag_cell[cell - 1] == CELL_UNPROCESSED) {
            stack_cell[++n_stack_cell] = cell;
            tag_cell[cell - 1] = CELL_IN_STACK;
          }
        }
      }
      else {
        int cell = PDM_ABS (_face_cell[i_face]);
        if (tag_cell[cell - 1] == CELL_UNPROCESSED) {
          stack_cell[++n_stack_cell] = cell;
        }
        tag_cell[cell - 1] = CELL_IN_STACK;
      }
    }

    /* Orientation process */

    while (n_stack_cell >= 0) {
      if (1 == 0) {
        printf("orientedface_cell : ");
        for (int i = 0; i < n_face; i++) {
          printf("%d : %d %d\n", i+1,orientedface_cell[2*i], orientedface_cell[2*i+1]);
        }
        printf("\n");

        printf("\n_stack_cell : ");
        for (int i = 0; i < n_stack_cell; i++) {
          printf(" %d", stack_cell[i]);
        }
        printf("\n");
      }

      n_edges = 0;

      int i_cell = stack_cell[n_stack_cell--] - 1;

      if (1 == 0) {
        printf("i_cell : %d\n", i_cell);
      }

      if (tag_cell[i_cell] == CELL_COMPLETED) {
        continue;
      }

      const int poly_idx   = cell_face_idx[i_cell];
      const int npoly_face = cell_face_idx[i_cell + 1] - poly_idx;

      /* Build pseudo edges of the current cell and store them into a hash table */

      n_processed_face = 0;

      for (int iface = 0; iface < npoly_face; iface++) {

        tag_face[iface] = FACE_UNPROCESSED;

        const int face          = PDM_ABS (cell_face[poly_idx + iface]) - 1;

        if (orientedface_cell[2*face] != 0) {
          assert (orientedface_cell[2*face  ] != i_cell + 1);
          assert (orientedface_cell[2*face+1] != i_cell + 1);
          assert (orientedface_cell[2*face+1] == 0        );
          tag_face[iface] = FACE_CHANGED_CYCLE;
          orientedface_cell[2*face+1] = i_cell + 1;
          processed_face[n_processed_face++] = iface;
        }
        else if (orientedface_cell[2*face + 1] != 0) {
          assert (orientedface_cell[2*face  ] != i_cell + 1);
          assert (orientedface_cell[2*face+1] != i_cell + 1);
          assert (orientedface_cell[2*face  ] == 0        );
          tag_face[iface] = FACE_UNCHANGED_CYCLE;
          orientedface_cell[2*face] = i_cell + 1;
          processed_face[n_processed_face++] = iface;
        }
      }

      if (n_processed_face == 0) {
        PDM_error (__FILE__, __LINE__, 0, "Error reorient : no processed face found\n");
      }

      for (int iface = 0; iface < npoly_face; iface++) {

        const int face          = PDM_ABS (cell_face[poly_idx + iface]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int vertex = face_vtx[faceIdx + ivert] - 1;

          const int inext = (ivert + 1) % n_faceVertices;
          const int vertex_next = face_vtx[faceIdx + inext] - 1;
          const int key = vertex + vertex_next;

          int *edge = edges + n_data_edge * n_edges;
          edge[0] = vertex;
          edge[1] = vertex_next;
          edge[2] = iface;

          n_edges += 1;

          PDM_hash_tab_data_add (hash_orient, (void *) &key, edge);

        }
      }

      /* FIXME Check edges connexion (ignore edges connected with more 2 faces) */

      int nStack_face = -1;

      /* Look for a neighbour of this face */

      for (int i = 0; i < n_processed_face; i++) {
        int _current_processed_face = processed_face[i];
        const int face          = PDM_ABS (cell_face[poly_idx + _current_processed_face]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int inext = (ivert + 1) % n_faceVertices;

          const int vertex = face_vtx[faceIdx + ivert] - 1;
          const int vertex_next = face_vtx[faceIdx + inext] - 1;
          int key = vertex + vertex_next;

          int nData = PDM_hash_tab_n_data_get (hash_orient, &key);

          int **data = (int **) PDM_hash_tab_data_get (hash_orient, &key);

          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isInverseEdge = (vertex == _edge[1]) && (vertex_next == _edge[0]);
              int isSameEdge    = (vertex == _edge[0]) && (vertex_next == _edge[1]);
              int isSameFace    = _current_processed_face == _edge[2];
              if (!isSameFace) {
                if (isSameEdge || isInverseEdge) {
                  if (tag_face[ _edge[2]] == FACE_UNPROCESSED) {
                    stack_face[++nStack_face] = _edge[2];
                    tag_face[ _edge[2]] = FACE_IN_STACK;
                  }
                  break;
                }
              }
            }
          }
        }
      }
      if (1 == 0) {
        printf("\nstack_face : ");
        for (int i = 0; i < nStack_face+1; i++) {
          printf(" %d", stack_face[i]);
        }
        printf("\n");
      }

      while (nStack_face >= 0) {

        int iFace = stack_face[nStack_face--];

        if ((tag_face[iFace] == FACE_UNCHANGED_CYCLE) ||
            (tag_face[iFace] == FACE_CHANGED_CYCLE)) {
          continue;
        }

        const int face          = PDM_ABS (cell_face[poly_idx + iFace]) - 1;
        const int faceIdx       = face_vtx_idx[face];
        const int n_faceVertices = face_vtx_idx[face+1] - faceIdx;

        for (int ivert = 0; ivert < n_faceVertices; ivert++) {
          const int inext = (ivert + 1) % n_faceVertices;

          const int vertex = face_vtx[faceIdx + ivert] - 1;
          const int vertex_next = face_vtx[faceIdx + inext] - 1;
          int key = vertex + vertex_next;

          int nData = PDM_hash_tab_n_data_get (hash_orient, &key);
          int **data = (int **) PDM_hash_tab_data_get (hash_orient, &key);

          int jCurrentEdge = -1;
          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isSameEdge    = (vertex == _edge[0]) && (vertex_next == _edge[1]);
              int isSameFace    = iFace == _edge[2];
              if (isSameEdge && isSameFace) {
                jCurrentEdge = j;
                break;
              }
            }
          }

          assert (jCurrentEdge > -1);

          for (int j = 0; j < nData; j++) {
            if (data[j] != NULL) {
              int *_edge = data[j];
              int isInverseEdge = (vertex == _edge[1]) && (vertex_next == _edge[0]);
              int isSameEdge    = (vertex == _edge[0]) && (vertex_next == _edge[1]);
              int isSameFace    = iFace == _edge[2];

              int neighbour = _edge[2];

              if (!isSameFace) {
                if (isSameEdge || isInverseEdge) {

                  if (tag_face[iFace] < FACE_UNCHANGED_CYCLE) {

                    if (tag_face[neighbour] >= FACE_UNCHANGED_CYCLE) {
                      if (tag_face[neighbour] == FACE_UNCHANGED_CYCLE) {
                        if (isSameEdge) {
                          tag_face[iFace] = FACE_CHANGED_CYCLE;
                          orientedface_cell[2*face+1] = i_cell+1;
                        }
                        else  {
                          tag_face[iFace] = FACE_UNCHANGED_CYCLE;
                          orientedface_cell[2*face] = i_cell+1;
                        }
                      }
                      else {
                        if (isSameEdge) {
                          tag_face[iFace] = FACE_UNCHANGED_CYCLE;
                          orientedface_cell[2*face] = i_cell+1;
                        }
                        else  {
                          tag_face[iFace] = FACE_CHANGED_CYCLE;
                          orientedface_cell[2*face+1] = i_cell+1;
                        }
                      }
                    }
                  }

                  if (tag_face[neighbour] == FACE_UNPROCESSED) {
                    stack_face[++nStack_face] = neighbour;
                    tag_face[neighbour] = FACE_IN_STACK;
                  }

                  break;

                }
              }
            }
          }
        }

        if (tag_face[iFace] == FACE_IN_STACK) {
          // printf(" oooo %i \n", iFace);
          // printf(" oooo %i \n", face);
          PDM_error (__FILE__, __LINE__, 0, "Error reorient : no neighbour processed face found\n");
        }
      }

      /* Add cell neighbours in the stack */

      for (int iface = 0; iface < npoly_face; iface++) {

        if (!((tag_face[iface] == FACE_UNCHANGED_CYCLE) ||
              (tag_face[iface] == FACE_CHANGED_CYCLE))) {
          printf(" oooo %i \n", iface);
          printf(" Link to : %i %i \n", face_cell[2*iface], face_cell[2*iface+1]);
          PDM_error (__FILE__, __LINE__, 0, "Error reorient : a face of polyhedron is not processed\n");
        }

        if (tag_face[iface] == FACE_CHANGED_CYCLE) {
          cell_face[poly_idx + iface] = -cell_face[poly_idx + iface];
        }

        tag_face[iface] = FACE_UNPROCESSED;

        const int face          = PDM_ABS (cell_face[poly_idx + iface]) - 1;

        int next_cell = -1;
        if (_face_cell[2 * face] == (i_cell + 1)) {
          if (_face_cell[2 * face + 1] != 0) {
            next_cell = _face_cell[2 * face + 1] - 1;
          }
        }
        else {
          next_cell = _face_cell[2 * face] - 1;
        }

        if (next_cell != -1) {
          if (tag_cell[next_cell] == CELL_UNPROCESSED) {
            stack_cell[++n_stack_cell] = next_cell + 1;
            tag_cell[next_cell] = CELL_IN_STACK;
          }
        }
      }

      PDM_hash_tab_purge(hash_orient, PDM_FALSE);

      tag_cell[i_cell] = CELL_COMPLETED;

    }

    int icheck = first_cell_comp;
    first_cell_comp = -1;
    for (int k = icheck; k < n_cell; k++) {
      if (tag_cell[k] == CELL_UNPROCESSED) {
        first_cell_comp = k;
        break;
      }
    }
  }

  PDM_hash_tab_free(hash_orient);

  /* Orient face_cell */
  if (face_cell != NULL) {
    for (int i = 0; i < n_cell; i++) {
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        int face = cell_face[j];
        int i_face = 2 * (PDM_ABS(face) - 1);
        if (PDM_ABS (face_cell[i_face]) == i+1) {
          if (face < 0) {
            face_cell[i_face] = -(i+1);
          }
          else {
            face_cell[i_face] = i+1;
          }
        }
        else {
          if (face < 0) {
            face_cell[i_face+1] = -(i+1);
          }
          else {
            face_cell[i_face+1] = i+1;
          }
        }
      }
    }
  }

  PDM_free(edges            );
  PDM_free(stack_face       );
  PDM_free(tag_face         );
  PDM_free(stack_cell       );
  PDM_free(tag_cell         );
  PDM_free(orientedface_cell);
  PDM_free(processed_face   );

  if (face_cell == NULL) {
    PDM_free(_face_cell);
  }

  PDM_timer_hang_on (t1);
  double et1 = PDM_timer_elapsed (t1);
  PDM_timer_free (t1);

  if (0 == 1) {
    printf("elapsed time cell_face_orient : %12.5e\n", et1);
  }
}

#ifdef	__cplusplus
}
#endif
