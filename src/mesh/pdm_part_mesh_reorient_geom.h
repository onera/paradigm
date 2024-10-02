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
);

#ifdef	__cplusplus
}
#endif
