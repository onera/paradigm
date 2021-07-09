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

#include <stdio.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "pdm_poly_vol_gen.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"


void
PDM_poly_vol_gen
(
 PDM_MPI_Comm  comm,
 double        xmin,
 double        ymin,
 double        zmin,
 double        lengthx,
 double        lengthy,
 double        lengthz,
 PDM_g_num_t   nx,
 PDM_g_num_t   ny,
 PDM_g_num_t   nz,
 int           randomize,
 int           random_seed,
 PDM_g_num_t  *ng_cell,
 PDM_g_num_t  *ng_face,
 PDM_g_num_t  *ng_vtx,
 int          *n_face_group,
 int          *dn_cell,
 int          *dn_face,
 int          *dn_vtx,
 int         **dcell_face_idx,
 PDM_g_num_t **dcell_face,
 int         **dface_cell_idx,
 PDM_g_num_t **dface_cell,
 int         **dface_vtx_idx,
 PDM_g_num_t **dface_vtx,
 double      **dvtx_coord,
 int         **dface_group_idx,
 PDM_g_num_t **dface_group
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t n_vtx1 = 2*nx;
  PDM_g_num_t n_vtx2 = nx + 1;
  PDM_g_num_t n_vtx3 = n_vtx1 + 2*n_vtx2;
  PDM_g_num_t n_vtx_z_cst = (ny+1)*2*nx + (nx+1)*2*ny + 4;


  PDM_g_num_t n_octo_z_cst = nx * ny;
  PDM_g_num_t n_quadH_z_cst = (nx - 1) * (ny - 1);
  PDM_g_num_t n_tria_z_cst = 2*(nx - 1) + 2*(ny - 1) + 4;//2 * (nx + ny);
  PDM_g_num_t n_quadV_z_cst = nx*(ny+1) + ny*(nx+1) + 4*nx*ny + 2*(nx-1) + 2*(ny-1) + 8;
  PDM_g_num_t n_faceH_z_cst = n_octo_z_cst + n_quadH_z_cst + n_tria_z_cst;
  PDM_g_num_t n_face_z_cst = n_faceH_z_cst + n_quadV_z_cst;
  PDM_g_num_t ng_face_lim = 2*(n_faceH_z_cst + 2*(nx+1) + 2*(ny+1));

  *ng_vtx = n_vtx_z_cst * (nz + 1);
  *ng_face = (nz + 1)*(n_face_z_cst) - n_quadV_z_cst;//nz*n_face_z_cst + n_faceH_z_cst;
  *ng_cell = n_faceH_z_cst * nz;
  *n_face_group = 6;

  /* Define distributions */
  PDM_g_num_t *distrib_vtx  = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_cell = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  PDM_g_num_t *distrib_face_lim = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_rank + 1));
  distrib_vtx[0]      = 0;
  distrib_face[0]     = 0;
  distrib_cell[0]     = 0;
  distrib_face_lim[0] = 0;

  PDM_g_num_t step_vtx       = *ng_vtx / n_rank;
  PDM_g_num_t remainder_vtx  = *ng_vtx % n_rank;

  PDM_g_num_t step_face      = *ng_face / n_rank;
  PDM_g_num_t remainder_face = *ng_face % n_rank;

  PDM_g_num_t step_cell      = *ng_cell / n_rank;
  PDM_g_num_t remainder_cell = *ng_cell % n_rank;

  PDM_g_num_t step_face_lim      = ng_face_lim / n_rank;
  PDM_g_num_t remainder_face_lim = ng_face_lim % n_rank;

  for (int i = 0; i < n_rank; i++) {
    distrib_vtx[i+1] = distrib_vtx[i] + step_vtx;
    if (i < remainder_vtx) {
      distrib_vtx[i+1]++;
    }

    distrib_face[i+1] = distrib_face[i] + step_face;
    if (i < remainder_face) {
      distrib_face[i+1]++;
    }

    distrib_cell[i+1] = distrib_cell[i] + step_cell;
    if (i < remainder_cell) {
      distrib_cell[i+1]++;
    }

    distrib_face_lim[i+1] = distrib_face_lim[i] + step_face_lim;
    if (i < remainder_face_lim) {
      distrib_face_lim[i+1]++;
    }
  }
  *dn_vtx  = (int) (distrib_vtx[i_rank+1]  - distrib_vtx[i_rank]);
  *dn_face = (int) (distrib_face[i_rank+1] - distrib_face[i_rank]);
  *dn_cell = (int) (distrib_cell[i_rank+1] - distrib_cell[i_rank]);
  int dn_face_lim = (int) (distrib_face_lim[i_rank+1] - distrib_face_lim[i_rank]);

  /*
   *  Vertices
   */
  *dvtx_coord = malloc (sizeof(double) * (*dn_vtx) * 3);
  double *_dvtx_coord = *dvtx_coord;

  double stepx = lengthx / (double) (3*nx);
  double stepy = lengthy / (double) (3*ny);
  double stepz = 0.;
  if (nz > 0) {
    stepz = lengthz / (double) nz;
  }

  for (int ivtx = 0; ivtx < *dn_vtx; ivtx++) {
    PDM_g_num_t g = distrib_vtx[i_rank] + ivtx;
    PDM_g_num_t k = g / n_vtx_z_cst;
    PDM_g_num_t r = g % n_vtx_z_cst;
    PDM_g_num_t i, j;
    double x, y;

    if (r > n_vtx_z_cst - 5) {
      // Corner
      r -= (n_vtx_z_cst - 4);
      j = r / 2;
      i = r % 2;

      x = xmin + i*lengthx;
      y = ymin + j*lengthy;
    }
    else {
      PDM_g_num_t jj = r / n_vtx3;
      PDM_g_num_t s = r % n_vtx3;

      if (s < n_vtx1) {
        // row n_vtx1
        i = 3*(s/2) + 1 + s%2;
        j = 3*jj;
      }
      else {
        // row n_vtx2
        s -= n_vtx1;
        j = 3*jj + s / n_vtx2 + 1;
        i = 3*(s % n_vtx2);
      }
      x = xmin + i*stepx;
      y = ymin + j*stepy;
    }

    _dvtx_coord[3*ivtx    ] = x;
    _dvtx_coord[3*ivtx + 1] = y;
    _dvtx_coord[3*ivtx + 2] = zmin + k*stepz;
  }


  /*
   *  Faces
   */
  PDM_g_num_t idx_octo  = 0;
  PDM_g_num_t idx_quadH = idx_octo + n_octo_z_cst;
  PDM_g_num_t idx_tria1 = idx_quadH + n_quadH_z_cst;
  PDM_g_num_t idx_tria2 = idx_tria1 + nx - 1;
  PDM_g_num_t idx_tria3 = idx_tria2 + nx - 1;
  PDM_g_num_t idx_tria4 = idx_tria3 + ny - 1;
  PDM_g_num_t idx_tria5 = idx_tria4 + ny - 1;
  PDM_g_num_t idx_quadV = n_faceH_z_cst;

  *dface_vtx_idx = malloc (sizeof(int) * (*dn_face + 1));
  int *_dface_vtx_idx = *dface_vtx_idx;
  _dface_vtx_idx[0] = 0;

  int s_face_vtx = 8*(*dn_face);
  *dface_vtx = malloc (sizeof(PDM_g_num_t) * s_face_vtx);

  for (int ifac = 0; ifac < (*dn_face); ifac++) {

    PDM_g_num_t *_dface_vtx = *dface_vtx + _dface_vtx_idx[ifac];

    PDM_g_num_t g = distrib_face[i_rank] + ifac;
    PDM_g_num_t k = g / n_face_z_cst;
    PDM_g_num_t r = g % n_vtx_z_cst;
    PDM_g_num_t i, j;

    if (r < idx_quadH) {
      // Octagon
      j = r / nx;
      i = r % nx;

      PDM_g_num_t idx = k*n_vtx_z_cst + j*n_vtx3 + 1;
      _dface_vtx[0] = idx + 2*i;
      _dface_vtx[1] = idx + 2*i + 1;
      _dface_vtx[2] = idx + n_vtx1 + i + 1;
      _dface_vtx[3] = idx + n_vtx1 + n_vtx2 + i + 1;
      _dface_vtx[4] = idx + n_vtx3 + 2*i + 1;
      _dface_vtx[5] = idx + n_vtx3 + 2*i;
      _dface_vtx[6] = idx + n_vtx1 + n_vtx2 + i;
      _dface_vtx[7] = idx + n_vtx1 + i;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 8;
    }

    else if (r < idx_tria1) {
      // Quadrangle
      r -= idx_quadH;

      j = r / (nx - 1);
      i = r % (nx - 1);

      PDM_g_num_t idx = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idx + i + 1;
      _dface_vtx[1] = idx + n_vtx2 + 2*(i+1);
      _dface_vtx[2] = idx + n_vtx1 + n_vtx2 + i + 1;
      _dface_vtx[3] = idx + n_vtx2 + 2*i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;
    }

    else if (r < idx_tria2) {
      // Triangle (row ymin)
      r -= idx_tria1;

      i = r % (nx - 1);

      PDM_g_num_t idx = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idx + 2*i + 1;
      _dface_vtx[1] = idx + 2*(i+1);
      _dface_vtx[2] = idx + n_vtx1 + i + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria3) {
      // Triangle (row ymax)
      r -= idx_tria2;

      i = r % (nx - 1);

      PDM_g_num_t idx = k*n_vtx_z_cst + ny*n_vtx3 - n_vtx2 + 1;
      _dface_vtx[0] = idx + n_vtx2 + 2*i + 1;
      _dface_vtx[1] = idx + i + 1;
      _dface_vtx[2] = idx + n_vtx2 + 2*(i+1);
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria4) {
      // Triangle (column xmin)
      r -= idx_tria3;

      j = r % (ny - 1);

      PDM_g_num_t idx = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idx;
      _dface_vtx[1] = idx + n_vtx2;
      _dface_vtx[2] = idx + n_vtx2 + n_vtx1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria5) {
      // Triangle (column xmax)
      r -= idx_tria4;

      j = r % (ny - 1);

      PDM_g_num_t idx = k*n_vtx_z_cst + n_vtx1 + n_vtx2 + j*n_vtx3 + 1;
      _dface_vtx[0] = idx + nx;
      _dface_vtx[1] = idx + n_vtx2 + n_vtx1 + nx;
      _dface_vtx[2] = idx + n_vtx2 + 2*(nx - 1) + 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria5 + 1) {
      // Triangle (bottom-left corner)
      PDM_g_num_t idx = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idx + n_vtx_z_cst - 4;
      _dface_vtx[1] = idx;
      _dface_vtx[2] = idx + n_vtx1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria5 + 2) {
      // Triangle (bottom-right corner)
      PDM_g_num_t idx = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idx + n_vtx1 - 1;
      _dface_vtx[1] = idx + n_vtx_z_cst - 3;
      _dface_vtx[2] = idx + n_vtx1 + n_vtx2 - 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria5 + 3) {
      // Triangle (top-left corner)
      PDM_g_num_t idx = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idx + (ny-1)*n_vtx3 + n_vtx1 + n_vtx2;
      _dface_vtx[1] = idx + ny*n_vtx3;
      _dface_vtx[2] = idx + n_vtx_z_cst - 2;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else if (r < idx_tria5 + 4) {
      // Triangle (top-right corner)
      PDM_g_num_t idx = k*n_vtx_z_cst + 1;
      _dface_vtx[0] = idx + (ny-1)*n_vtx3 + n_vtx1 + 2*n_vtx2 - 1;
      _dface_vtx[1] = idx + n_vtx_z_cst - 1;
      _dface_vtx[2] = idx + n_vtx_z_cst - 5;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 3;
    }

    else {
      _dface_vtx[0] = 1;
      _dface_vtx[1] = 1;
      _dface_vtx[2] = 1;
      _dface_vtx[3] = 1;
      _dface_vtx_idx[ifac+1] = _dface_vtx_idx[ifac] + 4;
    }

    if (r < idx_quadV && k == nz) {
      // flip face
      int n_vtx = _dface_vtx_idx[ifac+1] - _dface_vtx_idx[ifac];
      for (int l = 0; l < n_vtx/2; l++) {
        PDM_g_num_t tmp = _dface_vtx[l];
        _dface_vtx[l] = _dface_vtx[n_vtx - l - 1];
        _dface_vtx[n_vtx - l - 1] = tmp;
      }
    }
  }
}
