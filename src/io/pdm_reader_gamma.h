#ifndef __PDM_READER_GAMMA_H__
#define __PDM_READER_GAMMA_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a dmesh nodal from a file in ASCII GAMMA mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 * \param[in]  fix_orientation_2d  Ensure positive area for 2d faces (\warning mesh should be in (x,y) plane with pointing upward)
 * \param[in]  fix_orientation_3d  Ensure positive volume for 3d cells
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */

PDM_dmesh_nodal_t *
PDM_reader_gamma_dmesh_nodal
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int            fix_orientation_2d,
 int            fix_orientation_3d
);


void
PDM_write_meshb
(
  const char         *filename,
  const int          *n_elt_table,
        int         **tag_table,
        PDM_g_num_t **vtx_connect_table,
  const double       *vtx_coords
);


void
PDM_write_gamma_sol
(
  const char   *filename,
  const int     n_vtx,
  const int     n_field,
  const double *fields
);

void
PDM_read_gamma_sol
(
  const char   *filename,
  const int     n_vtx,
  const int     n_field,
        double *fields
);


void
PDM_write_gamma_matsym
(
  const char   *filename,
  const int     n_vtx,
  const double *fields
);


/**
 * \brief Read solution file in Gamma Mesh Format
 *
 * \param [in]  filename       Solution file name
 * \param [out] n_field        Number of fields
 * \param [out] field_stride   Field strides (size = \p n_field)
 * \param [out] field_values   Field values (size = \p n_field, for each field \p i, size = \p n_vtx * \p field_stride[i])
 *
 * \return Number of vertices
 */

int
PDM_read_gamma_sol_at_vertices
(
 const char   *filename,
 int          *n_field,
 int         **field_stride,
 double     ***field_values
 );


#ifdef __cplusplus
}
#endif
#endif
