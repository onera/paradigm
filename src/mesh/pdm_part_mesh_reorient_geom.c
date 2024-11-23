/*----------------------------------------------------------------------------
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_priv.h"
#include "pdm_plane.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_array.h"
#include "pdm_part_mesh_reorient_geom.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Fix mesh orientation geometrically
 *
 * \param [in]    mesh_dim               Mesh dimension (2 or 3)
 * \param [in]    mesh_inv               Mesh invariance direction (2D)
 * \param [in]    n_part                 Number of partition
 * \param [in]    n_cell                 Number of cells
 * \param [in]    n_face                 Number of faces
 * \param [in]    n_face_group           Number of face groups
 * \param [inout] cell_face_idx          Cell face connectivity index or NULL
 *                                       (size : n_cell + 1, numbering : 0 to n-1)
 * \param [inout] cell_face              Cell face connectivity or NULL
 *                                       (size : cell_face_idx[n_cell], numbering : 1 to n)
 * \param [inout] face_cell              Face cell connectivity or NULL
 *                                       (size : 2 * n_face, numbering : 1 to n)
 * \param [in]    face_vtx_idx           Face to vertex connectivity index
 *                                       (size : n_face + 1, numbering : 0 to n-1)
 * \param [inout] face_vtx               Face to vertex connectivity
 *                                       (size : face_vtx_idx[n_face], numbering : 1 to n)
 * \param [in]    vtx_coord              Vertex coordinates
 *                                       (size : 3*n_vtx)
 * \param [in]    face_group_idx         Index of faces list of each group
 *                                       (size = n_face_group + 1) or NULL
 * \param [in]    face_group             Faces list of each group
 *                                       (size = face_group[face_group_idx[n_face_group]], numbering : 1 to n)
 *                                       or NULL
 * \param [in]    comm                   MPI communicator
 *
 */

int
PDM_part_mesh_reorient_geom
(
 const int                    mesh_dim,
 const double                *mesh_inv,
 const int                    n_part,
 const int                   *n_cell,
 const int                   *n_face,
 const int                   *n_face_group,
       int                  **cell_face_idx,
       int                  **cell_face,
       int                  **face_cell,
       int                  **face_vtx_idx,
       int                  **face_vtx,
       double               **vtx_coord,
 const int                  **face_group_idx,
 const int                  **face_group,
 const int                    comm
)
{

  // init
  int  fixed  = 0;
  int  orient = 0;
  int _orient = 0;

  int **face_cell_signed_idx = NULL;
  int **face_cell_signed     = NULL;

  int **face_is_group;
  PDM_malloc (face_is_group, n_part, int*);
  for (int i_part = 0; i_part < n_part; i_part++) {
    face_is_group[i_part] = PDM_array_const_int(n_face[i_part], 0);
  }

  int cell_face_is_input = 0;
  if (cell_face != NULL) {
    cell_face_is_input = 1;
  }
  int face_cell_is_input = 0;
  if (face_cell != NULL) {
    face_cell_is_input = 1;
  }
  assert(cell_face_is_input == 1 || face_cell_is_input == 1);

  // convert cell_face/face_cell to signed connectivities
  if (face_cell_is_input == 0) {

    // transpose cell_face if not available
    PDM_part_connectivity_transpose(n_part,
                                    n_cell,
                                    n_face,
                                    cell_face_idx,
                                    cell_face,
                                   &face_cell_signed_idx,
                                   &face_cell_signed);

  } else {

    // convert face_cell only
    PDM_malloc (face_cell_signed_idx, n_part, int*);

    PDM_malloc (face_cell_signed, n_part, int*);

    for (int i_part = 0; i_part < n_part; i_part++) {

      int* _face_cell_signed_idx;
      PDM_malloc (_face_cell_signed_idx, (n_face[i_part]+1),int);

      int* _face_cell_signed;
      PDM_malloc (_face_cell_signed, 2 * n_face[i_part],int);

      _face_cell_signed_idx[0] = 0;
      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
        _face_cell_signed_idx[i_face+1] = _face_cell_signed_idx[i_face];
        if (face_cell[i_part][2*i_face  ] != 0) {
          _face_cell_signed[_face_cell_signed_idx[i_face+1]++] =  PDM_ABS(face_cell[i_part][2*i_face]);
        }
        if (face_cell[i_part][2*i_face+1] != 0) {
          _face_cell_signed[_face_cell_signed_idx[i_face+1]++] = -PDM_ABS(face_cell[i_part][2*i_face+1]);
        }
      }

      PDM_realloc(_face_cell_signed, _face_cell_signed, _face_cell_signed_idx[n_face[i_part]],int);

      face_cell_signed_idx[i_part] = _face_cell_signed_idx;
      face_cell_signed    [i_part] = _face_cell_signed;

    }

  }

  // check groups if available
  if (face_group != NULL) {

    for (int i_part = 0; i_part < n_part; i_part++) {

      int *_face_cell_idx  = face_cell_signed_idx[i_part];
      int *_face_cell      = face_cell_signed    [i_part];
      int *_face_is_group  = face_is_group       [i_part];

      for (int i_face_group = 0; i_face_group < n_face_group[i_part]; i_face_group++) {
        for (int i_face = face_group_idx[i_part][i_face_group]; i_face < face_group_idx[i_part][i_face_group+1]; i_face++) {
          int j_face = face_group[i_part][i_face]-1;
          assert(_face_cell_idx[j_face+1]-_face_cell_idx[j_face] == 1);
          _face_is_group[j_face] = 1;
          int idx_face = _face_cell_idx[j_face];
          if (_face_cell[idx_face] < 0) {
            orient++;
            _face_cell[idx_face] *= -1; // positive face_cell for inner cell
          }
        }
      }

    }

    PDM_MPI_Allreduce(&orient, &_orient, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
    if (_orient > 0) fixed = 1;

  }

  // update cell_face connectivity
  int **lcell_face_idx = NULL;
  int **lcell_face     = NULL;
  if (cell_face_is_input == 0) {
    PDM_malloc (lcell_face_idx, n_part, int*);
    PDM_malloc (lcell_face, n_part, int*);

  } else {
    lcell_face_idx = cell_face_idx;
    lcell_face     = cell_face;
  }
  if (cell_face_is_input == 0 || orient > 0) {

    for (int i_part = 0; i_part < n_part; i_part++) {

      int* _cell_face_idx = NULL;
      int* _cell_face     = NULL;

      PDM_connectivity_transpose(n_face[i_part],
                                 n_cell[i_part],
                                 face_cell_signed_idx[i_part],
                                 face_cell_signed    [i_part],
                                &_cell_face_idx,
                                &_cell_face);

      if (cell_face_is_input == 0) {
        lcell_face_idx[i_part] = _cell_face_idx;
        lcell_face    [i_part] = _cell_face;
      } else {
        for (int i_cell = 0; i_cell < n_cell[i_part]; i_cell++) {
          cell_face_idx[i_part][i_cell+1] = _cell_face_idx[i_cell+1];
          for (int idx_face = _cell_face_idx[i_cell]; idx_face < _cell_face_idx[i_cell+1]; idx_face++) {
            cell_face[i_part][idx_face] = _cell_face[idx_face];
          }
        }
        PDM_free(_cell_face_idx);
        PDM_free(_cell_face    );
      }

    }

  }

   orient = 0;
  _orient = 0;

  // check face_cell connectivity
  for (int i_part = 0; i_part < n_part; i_part++) {

    double *face_center;
    PDM_malloc (face_center, 3 * n_face[i_part], double);

    double *cell_center;
    PDM_malloc (cell_center, 3 * n_cell[i_part], double);

    int *_face_cell_idx = face_cell_signed_idx[i_part];
    int *_face_cell     = face_cell_signed    [i_part];
    int *_face_is_group = face_is_group       [i_part];

    if (mesh_dim == 3) {

      // face center
      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
        face_center[3*i_face  ] = 0.;
        face_center[3*i_face+1] = 0.;
        face_center[3*i_face+2] = 0.;
        int n_vtx_on_face = face_vtx_idx[i_part][i_face+1] - face_vtx_idx[i_part][i_face];
        for (int idx_vtx = face_vtx_idx[i_part][i_face]; idx_vtx < face_vtx_idx[i_part][i_face+1]; idx_vtx++) {
          int i_vtx = face_vtx[i_part][idx_vtx]-1;
          face_center[3*i_face  ] += vtx_coord[i_part][3*i_vtx  ];
          face_center[3*i_face+1] += vtx_coord[i_part][3*i_vtx+1];
          face_center[3*i_face+2] += vtx_coord[i_part][3*i_vtx+2];
        }
        face_center[3*i_face  ] = face_center[3*i_face  ] / n_vtx_on_face;
        face_center[3*i_face+1] = face_center[3*i_face+1] / n_vtx_on_face;
        face_center[3*i_face+2] = face_center[3*i_face+2] / n_vtx_on_face;
      }

    } else {

      // edge center
      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
        face_center[3*i_face  ] = 0.;
        face_center[3*i_face+1] = 0.;
        face_center[3*i_face+2] = 0.;
        for (int idx_vtx = 2*i_face; idx_vtx < 2*(i_face+1); idx_vtx++) {
          int i_vtx = face_vtx[i_part][idx_vtx]-1;
          face_center[3*i_face  ] += vtx_coord[i_part][3*i_vtx  ];
          face_center[3*i_face+1] += vtx_coord[i_part][3*i_vtx+1];
          face_center[3*i_face+2] += vtx_coord[i_part][3*i_vtx+2];
        }
        face_center[3*i_face  ] = face_center[3*i_face  ] / 2;
        face_center[3*i_face+1] = face_center[3*i_face+1] / 2;
        face_center[3*i_face+2] = face_center[3*i_face+2] / 2;
      }

    }

    // cell center
    for (int i_cell = 0; i_cell < n_cell[i_part]; i_cell++) {
      cell_center[3*i_cell  ] = 0.;
      cell_center[3*i_cell+1] = 0.;
      cell_center[3*i_cell+2] = 0.;
      int n_face_on_cell = lcell_face_idx[i_part][i_cell+1] - lcell_face_idx[i_part][i_cell];
      for (int idx_face = lcell_face_idx[i_part][i_cell]; idx_face < lcell_face_idx[i_part][i_cell+1]; idx_face++) {
        int i_face = PDM_ABS(lcell_face[i_part][idx_face])-1;
        cell_center[3*i_cell  ] += face_center[3*i_face  ];
        cell_center[3*i_cell+1] += face_center[3*i_face+1];
        cell_center[3*i_cell+2] += face_center[3*i_face+2];
      }
      cell_center[3*i_cell  ] = cell_center[3*i_cell  ] / n_face_on_cell;
      cell_center[3*i_cell+1] = cell_center[3*i_cell+1] / n_face_on_cell;
      cell_center[3*i_cell+2] = cell_center[3*i_cell+2] / n_face_on_cell;
    }

    if (mesh_dim == 3) {

      int n_vtx_on_face_max = 0;

      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
        n_vtx_on_face_max = PDM_MAX(n_vtx_on_face_max, face_vtx_idx[i_part][i_face+1]-face_vtx_idx[i_part][i_face]);
      }

      double *_vtx_on_face;
      PDM_malloc (_vtx_on_face, 3 * n_vtx_on_face_max, double);

      // orient geometrically with simplified normal definition
      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {

        double normal[3];
        double tmp[3];
        int idx_write = 0;
        for (int idx_vtx = face_vtx_idx[i_part][i_face]; idx_vtx < face_vtx_idx[i_part][i_face+1]; idx_vtx++) {
          int i_vtx = face_vtx[i_part][idx_vtx]-1;
          _vtx_on_face[3*idx_write  ] = vtx_coord[i_part][3*i_vtx  ];
          _vtx_on_face[3*idx_write+1] = vtx_coord[i_part][3*i_vtx+1];
          _vtx_on_face[3*idx_write+2] = vtx_coord[i_part][3*i_vtx+2];
          idx_write++;
        }
        PDM_plane_normal(face_vtx_idx[i_part][i_face+1]-face_vtx_idx[i_part][i_face],
                         _vtx_on_face,
                         normal);
        int idx_face = _face_cell_idx[i_face];
        int i_cell   = PDM_ABS(_face_cell[idx_face])-1;
        tmp[0] = face_center[3*i_face  ] - cell_center[3*i_cell  ];
        tmp[1] = face_center[3*i_face+1] - cell_center[3*i_cell+1];
        tmp[2] = face_center[3*i_face+2] - cell_center[3*i_cell+2];
        double scal = normal[0]*tmp[0] + normal[1]*tmp[1] + normal[2]*tmp[2];
        // internal face
        if (_face_cell_idx[i_face+1]-_face_cell_idx[i_face] == 2) {
          assert(_face_cell[idx_face]*_face_cell[idx_face+1] < 0);
          if (_face_cell[idx_face] > 0 && scal < 0.) {
            orient++;
            _face_cell[idx_face  ] = -PDM_ABS(_face_cell[idx_face  ]); // negative face_cell right
            _face_cell[idx_face+1] =  PDM_ABS(_face_cell[idx_face+1]); // positive face_cell left
          } else if (_face_cell[idx_face] < 0 && scal > 0.) {
            orient++;
            _face_cell[idx_face  ] =  PDM_ABS(_face_cell[idx_face  ]); // positive face_cell left
            _face_cell[idx_face+1] = -PDM_ABS(_face_cell[idx_face+1]); // negative face_cell right
          }
        // boundary face
        } else {
          // groupe face
          if (_face_is_group[i_face] == 1) {
            if (scal < 0.) {
              orient++;
              int idx_vtx = face_vtx_idx[i_part][i_face];
              int n_vtx_on_face = face_vtx_idx[i_part][i_face+1] - face_vtx_idx[i_part][i_face];
              for (int i_vtx = 0; i_vtx < n_vtx_on_face/2; i_vtx++) { // revert face_vtx
                int j_vtx = face_vtx[i_part][idx_vtx+i_vtx];
                face_vtx[i_part][idx_vtx+i_vtx                ] = face_vtx[i_part][idx_vtx+n_vtx_on_face-1-i_vtx];
                face_vtx[i_part][idx_vtx+n_vtx_on_face-1-i_vtx] = j_vtx;
              }
            }
          // partitioning face
          } else {
            if (_face_cell[idx_face] > 0 && scal < 0.) {
              orient++;
              _face_cell[idx_face] = -PDM_ABS(_face_cell[idx_face]); // negative face_cell right
            } else if (_face_cell[idx_face] < 0 && scal > 0.) {
              orient++;
              _face_cell[idx_face] =  PDM_ABS(_face_cell[idx_face]); // positive face_cell left
            }
          }
        }

      }

      PDM_free(_vtx_on_face);

    } else {

      // orient geometrically with exact normal definition
      for (int i_edge = 0; i_edge < n_face[i_part]; i_edge++) {

        double edge[3];
        double normal[3];
        double tmp[3];
        int i_vtx = face_vtx[i_part][2*i_edge  ]-1;
        int j_vtx = face_vtx[i_part][2*i_edge+1]-1;
        edge[0] = vtx_coord[i_part][3*j_vtx  ] - vtx_coord[i_part][3*i_vtx  ];
        edge[1] = vtx_coord[i_part][3*j_vtx+1] - vtx_coord[i_part][3*i_vtx+1];
        edge[2] = vtx_coord[i_part][3*j_vtx+2] - vtx_coord[i_part][3*i_vtx+2];
        normal[0] = mesh_inv[2]*edge[1] - mesh_inv[1]*edge[2];
        normal[1] = mesh_inv[0]*edge[2] - mesh_inv[2]*edge[0];
        normal[2] = mesh_inv[1]*edge[0] - mesh_inv[0]*edge[1];
        int idx_face = _face_cell_idx[i_edge];
        int i_face   = PDM_ABS(_face_cell[idx_face])-1;
        tmp[0] = face_center[3*i_edge  ] - cell_center[3*i_face  ];
        tmp[1] = face_center[3*i_edge+1] - cell_center[3*i_face+1];
        tmp[2] = face_center[3*i_edge+2] - cell_center[3*i_face+2];
        double scal = normal[0]*tmp[0] + normal[1]*tmp[1] + normal[2]*tmp[2];
        // internal edge
        if (_face_cell_idx[i_edge+1]-_face_cell_idx[i_edge] == 2) {
          assert(_face_cell[idx_face]*_face_cell[idx_face+1] < 0);
          if (_face_cell[idx_face] > 0 && scal < 0.) {
            orient++;
            _face_cell[idx_face  ] = -PDM_ABS(_face_cell[idx_face  ]); // negative face_cell right
            _face_cell[idx_face+1] =  PDM_ABS(_face_cell[idx_face+1]); // positive face_cell left
          } else if (_face_cell[idx_face] < 0 && scal > 0.) {
            orient++;
            _face_cell[idx_face  ] =  PDM_ABS(_face_cell[idx_face  ]); // positive face_cell left
            _face_cell[idx_face+1] = -PDM_ABS(_face_cell[idx_face+1]); // negative face_cell right
          }
        // boundary edge
        } else {
          // groupe edge
          if (_face_is_group[i_edge] == 1) {
            if (scal < 0.) {
              orient++;
              face_vtx[i_part][2*i_edge  ] = j_vtx+1; // revert face_vtx
              face_vtx[i_part][2*i_edge+1] = i_vtx+1;
            }
          // partitioning edge
          } else {
            if (_face_cell[idx_face] > 0 && scal < 0.) {
              orient++;
              _face_cell[idx_face] = -PDM_ABS(_face_cell[idx_face]); // negative face_cell right
            } else if (_face_cell[idx_face] < 0 && scal > 0.) {
              orient++;
              _face_cell[idx_face] =  PDM_ABS(_face_cell[idx_face]); // positive face_cell left
            }
          }
        }

      }

    }

    PDM_free(face_center);
    PDM_free(cell_center);

  }

  PDM_MPI_Allreduce(&orient, &_orient, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  if (_orient > 0) fixed = 1;

  // convert face_cell to unsigned connectivities if available in input
  if (face_cell_is_input == 1) {

    for (int i_part = 0; i_part < n_part; i_part++) {

      int* _face_cell_signed_idx = face_cell_signed_idx[i_part];
      int* _face_cell_signed     = face_cell_signed    [i_part];

      for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
        face_cell[i_part][2*i_face  ] = 0;
        face_cell[i_part][2*i_face+1] = 0;
        for (int idx_face = _face_cell_signed_idx[i_face]; idx_face < _face_cell_signed_idx[i_face+1]; idx_face++) {
          if (_face_cell_signed[idx_face] > 0) {
            face_cell[i_part][2*i_face  ] =  _face_cell_signed[idx_face]; // positive face_cell left
          } else {
            face_cell[i_part][2*i_face+1] = -_face_cell_signed[idx_face]; // negative face_cell right
          }
        }
      }

    }

  }

  // update cell_face connectivity if available in input
  if (cell_face_is_input == 1 && orient > 0) {

    for (int i_part = 0; i_part < n_part; i_part++) {

      int* _cell_face_idx = NULL;
      int* _cell_face     = NULL;

      PDM_connectivity_transpose(n_face[i_part],
                                 n_cell[i_part],
                                 face_cell_signed_idx[i_part],
                                 face_cell_signed    [i_part],
                                &_cell_face_idx,
                                &_cell_face);

      for (int i_cell = 0; i_cell < n_cell[i_part]; i_cell++) {
        cell_face_idx[i_part][i_cell+1] = _cell_face_idx[i_cell+1];
        for (int idx_face = _cell_face_idx[i_cell]; idx_face < _cell_face_idx[i_cell+1]; idx_face++) {
          cell_face[i_part][idx_face] = _cell_face[idx_face];
        }
      }

      PDM_free(_cell_face_idx);
      PDM_free(_cell_face    );

    }

  }

  // free memory
  if (cell_face_is_input == 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(lcell_face_idx[i_part]);
      PDM_free(lcell_face    [i_part]);
    }
    PDM_free(lcell_face_idx);
    PDM_free(lcell_face    );
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(face_cell_signed_idx[i_part]);
    PDM_free(face_cell_signed    [i_part]);
    PDM_free(face_is_group       [i_part]);
  }
  PDM_free(face_cell_signed_idx);
  PDM_free(face_cell_signed    );
  PDM_free(face_is_group       );

  return fixed;
}
#ifdef	__cplusplus
}
#endif
