#ifndef __PDM_ISO_SURFACE_H__
#define __PDM_ISO_SURFACE_H__

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

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

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_iso_surface_t PDM_iso_surface_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Build an iso surface struct
 *
 *
 * \return     Identifier
 */
PDM_iso_surface_t*
PDM_iso_surface_create
(
 const int                    dim,
       PDM_iso_surface_kind_t iso_kind,
 const int                    n_part,
       PDM_ownership_t        ownership,
       PDM_MPI_Comm           comm
);


void
PDM_iso_surface_compute
(
  PDM_iso_surface_t        *isos
);

// See with Eric et Bastien : par type ou une fonction avec 1000 arguments ?
void
PDM_iso_surface_part_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);

void
PDM_iso_surface_part_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *field
);

void
PDM_iso_surface_part_gradient_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *gradient_field
);

// See with Eric et Bastien : par type ou une fonction avec 1000 arguments ?
void
PDM_iso_surface_dconnectivity_set
(
  PDM_iso_surface_t        *isos,
  PDM_connectivity_type_t   connectivity_type,
  PDM_g_num_t              *dconnect,
  int                      *dconnect_idx
);

void
PDM_iso_surface_distrib_set
(
  PDM_iso_surface_t        *isos,
  PDM_mesh_entities_t       entity_type,
  PDM_g_num_t              *distrib_entity
);

void
PDM_iso_surface_dvtx_coord_set
(
  PDM_iso_surface_t *isos,
  double            *dvtx_coord
);

void
PDM_iso_surface_dfield_set
(
  PDM_iso_surface_t *isos,
  double            *dfield
);

void
PDM_iso_surface_dgrad_field_set
(
  PDM_iso_surface_t *isos,
  double            *dgrad_field
);


void
PDM_iso_surface_get
(
  PDM_iso_surface_t  *isos,
  PDM_g_num_t       **diso_face_vtx,
  int               **diso_face_vtx_idx,
  double            **diso_vtx_coord
);


// iso_vtx_to_cell -> Pour chaque iso_vtx, la connectivité avec la cellule qui contient ce vtx
// iso_vtx_to_vtx  -> Pour chaque iso_vtx, la connectivité (a trié) avec tous les vtx qui on contribué à l'iso
// @Bastien : Est-ce que tu aurai un poids a associé a cette connectivité ?
void
PDM_iso_surface_graph_comm_get
(
  PDM_iso_surface_t   *isos,
  PDM_g_num_t        **iso_vtx_to_cell,
  PDM_g_num_t        **iso_vtx_to_vtx,
  int                **iso_vtx_to_vtx_idx
);

void
PDM_iso_surface_plane_equation_set
(
  PDM_iso_surface_t        *isos,
  double                    a,
  double                    b,
  double                    c,
  double                    d
);

void
PDM_iso_surface_free
(
  PDM_iso_surface_t  *isos
);


void
PDM_iso_surface_write
(
 PDM_iso_surface_t  *isos,
 const char         *name
 );

void
PDM_iso_surface_eval_field_and_gradient_set
(
 PDM_iso_surface_t *isos,
 void (*eval_field_and_gradient) (const double, const double, const double,
                                   double *,
                                   double *, double *, double *)
 );


void
PDM_iso_surface_surface_get
(
 PDM_iso_surface_t  *isos,
 int                *n_vtx,
 int                *n_elt,
 int               **elt_vtx_idx,
 int               **elt_vtx,
 double            **vtx_coord,
 PDM_g_num_t       **elt_ln_to_gn,
 PDM_g_num_t       **vtx_ln_to_gn,
 PDM_g_num_t       **elt_parent_g_num
 );

void
PDM_iso_surface_dump_times
(
 PDM_iso_surface_t  *isos
 );

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_ISO_SURFACE_H__ */
