/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2023       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
#include <stdlib.h>

#include "pdm.h"
#include "pdm_mpi.h"

#include "pdm_error.h"

#include "pdm_mesh_nodal.h"
#include "pdm_isosurface.h"
#include "pdm_isosurface_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


void
_check_entry_mesh_coherence
(
 PDM_isosurface_t *isos,
 int               entry_mesh_type
)
{
  if (isos->entry_mesh_type==0) {
    isos->entry_mesh_type=entry_mesh_type;
  }
  else if (isos->entry_mesh_type!=entry_mesh_type) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t:entry_mesh_type already set to %d.\n", isos->entry_mesh_type);
  }
}


PDM_isosurface_t *
PDM_isosurface_create
(
 PDM_MPI_Comm             comm,
 int                      mesh_dimension,
 PDM_Mesh_nodal_elt_t     elt_type
)
{
  PDM_isosurface_t *isos = (PDM_isosurface_t *) malloc(sizeof(PDM_isosurface_t));

  // > Save communicator
  isos->comm = comm;

  // > Entry mesh information
  isos->is_dist_or_part = -1; 
  isos->entry_mesh_type =  0; 
  isos->entry_mesh_dim  =  mesh_dimension;

  // > Isosurface mesh information
  isos->iso_elt_type = elt_type; 

  // > Isosurfaces informations
  isos->n_isosurface = 0;
  isos->n_isovalues       = NULL;
  isos->  isovalues       = NULL;
  isos->kind            = NULL;
  isos->field_function  = NULL;
  isos->eq_coeffs       = NULL;
  isos->use_gradient    = NULL;
  isos->iso_func        = NULL;

  // > Distributed entry data
  isos->dmesh        = NULL;
  isos->dmesh_nodal  = NULL;

  isos->distrib_cell = NULL;
  isos->distrib_face = NULL;
  isos->distrib_edge = NULL;
  isos->distrib_vtx  = NULL;

  isos->dvtx_coord     = NULL;
  isos->dcell_face_idx = NULL;
  isos->dcell_face     = NULL;
  isos->dface_edge_idx = NULL;
  isos->dface_edge     = NULL;
  isos->dface_vtx_idx  = NULL;
  isos->dface_vtx      = NULL;
  isos->dedge_vtx      = NULL;

  isos->dgroup_face_idx = NULL;
  isos->dgroup_face     = NULL;
  isos->dgroup_edge_idx = NULL;
  isos->dgroup_edge     = NULL;
  isos->dgroup_vtx_idx  = NULL;
  isos->dgroup_vtx      = NULL;

  isos->dfield    = NULL;
  isos->dgradient = NULL;

  // > Partitioned entry data
  isos->pmesh       = NULL;
  isos->pmesh_nodal = NULL;

  isos->n_cell = NULL;
  isos->n_face = NULL;
  isos->n_edge = NULL;
  isos->n_vtx  = NULL;

  isos->cell_gnum = NULL;
  isos->face_gnum = NULL;
  isos->edge_gnum = NULL;
  isos-> vtx_gnum = NULL;

  isos->vtx_coord     = NULL;
  isos->cell_face_idx = NULL;
  isos->cell_face     = NULL;
  isos->face_edge_idx = NULL;
  isos->face_edge     = NULL;
  isos->face_vtx_idx  = NULL;
  isos->face_vtx      = NULL;
  isos->edge_vtx      = NULL;

  isos->n_group_face   = NULL;
  isos->group_face_idx = NULL;
  isos->group_face     = NULL;
  isos->n_group_edge   = NULL;
  isos->group_edge_idx = NULL;
  isos->group_edge     = NULL;
  isos->n_group_vtx    = NULL;
  isos->group_vtx_idx  = NULL;
  isos->group_vtx      = NULL;

  isos->field    = NULL;
  isos->gradient = NULL;


  return isos;
}


int
PDM_isosurface_add
(
 PDM_isosurface_t       *isos,
 PDM_iso_surface_kind_t  kind,
 int                     n_isovalues,
 double                 *isovalues
)
{

  int id_isosurface = isos->n_isosurface;
  isos->n_isosurface++;
  
  printf("isos->n_isosurface = %d\n", isos->n_isosurface);
  isos->kind           = realloc(isos->kind          , sizeof(PDM_iso_surface_kind_t          ) * isos->n_isosurface);
  isos->eq_coeffs      = realloc(isos->eq_coeffs     , sizeof(double                         *) * isos->n_isosurface);
  isos->field_function = realloc(isos->field_function, sizeof(PDM_isosurface_field_function_t ) * isos->n_isosurface);
  
  if (isos->is_dist_or_part==0) { // is distributed
    isos->dfield         = realloc(isos->dfield        , sizeof(double                         *) * isos->n_isosurface);
  } else if (isos->is_dist_or_part==1) { // is partitioned
    isos->field          = realloc(isos->field         , sizeof(double                        **) * isos->n_isosurface);
  }
  isos->n_isovalues    = realloc(isos->n_isovalues   , sizeof(int                             ) * isos->n_isosurface);
  isos->isovalues      = realloc(isos->isovalues     , sizeof(double                         *) * isos->n_isosurface);
  isos->use_gradient   = realloc(isos->use_gradient  , sizeof(int                            *) * isos->n_isosurface);

  isos->kind       [id_isosurface] = kind;
  isos->n_isovalues[id_isosurface] = n_isovalues;
  isos->  isovalues[id_isosurface] = malloc(sizeof(double) * n_isovalues);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[id_isosurface][i] = isovalues[i];
  }

  if        (kind==PDM_ISO_SURFACE_KIND_FIELD) {
    isos->eq_coeffs     [id_isosurface] = NULL;
  }
  else if (kind==PDM_ISO_SURFACE_KIND_PLANE) {
    int n_coeff = 4;
    isos->eq_coeffs[id_isosurface] = malloc(sizeof(double) * n_coeff);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_SPHERE) {
    int n_coeff = 4;
    isos->eq_coeffs[id_isosurface] = malloc(sizeof(double) * n_coeff);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_ELLIPSE) {
    int n_coeff = 6;
    isos->eq_coeffs[id_isosurface] = malloc(sizeof(double) * n_coeff);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_QUADRIC) {
    int n_coeff = 10;
    isos->eq_coeffs[id_isosurface] = malloc(sizeof(double) * n_coeff);
  }
  else if (kind==PDM_ISO_SURFACE_KIND_HEART) {
    isos->eq_coeffs[id_isosurface] = NULL;
  }
  else if (kind==PDM_ISO_SURFACE_KIND_FUNCTION) {
    isos->eq_coeffs[id_isosurface] = NULL;
  }

  return id_isosurface;
}


void
PDM_isosurface_equation_set
(
 PDM_isosurface_t *isos,
 int               id_isosurface,
 double           *coeff,
 int               use_gradient
)
{
  if        (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_PLANE) {
    int n_coeff = 4;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_SPHERE) {
    int n_coeff = 4;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_ELLIPSE) {
    int n_coeff = 6;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_QUADRIC) {
    int n_coeff = 10;
    for (int i_coeff=0; i_coeff<n_coeff; ++i_coeff) {
      isos->eq_coeffs[id_isosurface][i_coeff] = coeff[i_coeff];
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_equation_set method cause its kind is %d.\n", id_isosurface, isos->kind[id_isosurface]);
  }

  isos->use_gradient[id_isosurface] = use_gradient;

}


void
PDM_isosurface_field_function_set
(
 PDM_isosurface_t                *isos,
 int                              id_isosurface,
 PDM_isosurface_field_function_t  func
)
{
  if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FUNCTION) {
    isos->field_function[id_isosurface] = func;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface n°%d doesn't support PDM_isosurface_field_function_set method cause its kind is %d.\n", id_isosurface, isos->kind[id_isosurface]);
  }

}


void
PDM_isosurface_reset
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{

}


void
PDM_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{

}


void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
)
{

}


void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
)
{
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    printf("id_iso = %d\n", id_iso);
    if (isos->n_isovalues[id_iso]>0) {
      free(isos->isovalues[id_iso]);
      if (isos->kind[id_iso]==PDM_ISO_SURFACE_KIND_FIELD) {
        if (isos->is_dist_or_part==0) { // is distributed
          free(isos->dfield[id_iso]);
        } else if (isos->is_dist_or_part==1) { // is partitioned
          printf("\t free isos->field[id_iso]\n");
          free(isos->field[id_iso]); // TODO: check if field is desallocated by iso or not
        }
      } else  {
        free(isos->eq_coeffs[id_iso]);
      }
    }
  }

  free(isos->kind);
  free(isos->n_isovalues);
  free(isos->isovalues);
  free(isos->eq_coeffs);
  free(isos->use_gradient);
  free(isos->field_function);
  if (isos->is_dist_or_part==0) { // is distributed
    free(isos->dfield);
  } else if (isos->is_dist_or_part==1) { // is partitioned
          printf("\t free isos->field[id_iso]\n");
    free(isos->field);
  }

  free(isos);
}


#ifdef  __cplusplus
}
#endif