/*
 * \file
 */

#ifndef __PDM_READER_GAMMA_H__
#define __PDM_READER_GAMMA_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */



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
PDM_write_meshb(
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
