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
#include <math.h>

#include "pdm.h"
#include "pdm_mpi.h"

#include "pdm_error.h"
#include "pdm_logging.h"

#include "pdm_array.h"
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

static
inline
double
_plane_field
(
  double x,
  double y,
  double z,
  double *plane_equation
)
{
  return  plane_equation[0] * x
        + plane_equation[1] * y
        + plane_equation[2] * z
        - plane_equation[3];
}



static
inline
double
_sphere_field
(
  double x,
  double y,
  double z,
  double *sphere_equation
)
{
  return   pow(x-sphere_equation[0], 2.)
         + pow(y-sphere_equation[1], 2.)
         + pow(z-sphere_equation[2], 2.)
         - pow(  sphere_equation[3], 2.) ;
}



static
inline
double
_ellipse_field
(
  double x,
  double y,
  double z,
  double *ellipse_equation
)
{
  return   pow( (x-ellipse_equation[0]) / ellipse_equation[3], 2.)
         + pow( (y-ellipse_equation[1]) / ellipse_equation[4], 2.)
         + pow( (z-ellipse_equation[2]) / ellipse_equation[5], 2.)
         -         ellipse_equation[6];
}



static
inline
double
_quadric_field
(
  double x,
  double y,
  double z,
  double *quadric_equation
)
{
  return   quadric_equation[6] * pow( (x-quadric_equation[0]) / quadric_equation[3], 2.)
         + quadric_equation[7] * pow( (y-quadric_equation[1]) / quadric_equation[4], 2.)
         + quadric_equation[8] * pow( (z-quadric_equation[2]) / quadric_equation[5], 2.)
         - quadric_equation[9];
}



static
inline
double
_heart_field
(
  double x,
  double y,
  double z,
  double *heart_equation
)
{
  PDM_UNUSED(heart_equation);
  double a=1.;
  double b=2.;
  // return pow(x*x + y*y -1,3.) - x*x*y*y*y; // 2D
  return pow( pow(       x,2.)
            + pow((1.+b)*y,2.)
            + pow(       z,2.)
            -            1.   ,3.)
        -   x*x*z*z*z
        - a*y*y*z*z*z;
        // +           1.;
}



static void
_compute_iso_field
(
  PDM_isosurface_t *isos,
  int               id_isosurface
)
{
  int debug = 1;

  PDM_part_mesh_nodal_t *pmn = isos->pmesh_nodal;
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);

  // > Allocate if 
  if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_PLANE    ||
      isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_SPHERE   ||
      isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_ELLIPSE  ||
      isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_QUADRIC  ||
      isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_HEART    ||
      isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FUNCTION ) {
    isos->field[id_isosurface] = malloc(sizeof(double *) * n_part); 
  }
  else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FIELD) {
    if (isos->field[id_isosurface]==NULL) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: iso n°%d is PDM_ISO_SURFACE_KIND_FIELD, but field isn't defined.\n", id_isosurface);
    }
  }

  double **vtx_field = isos->field[id_isosurface];


  for (int i_part=0; i_part<n_part; ++i_part) {
    int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);


    if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FIELD) {
      if (isos->field[id_isosurface][i_part]==NULL) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: iso n°%d is PDM_ISO_SURFACE_KIND_FIELD, but field of part n°%d isn't defined.\n", id_isosurface, i_part);
      }
    } else {
      
      log_trace("id_isosurface = %d ; i_part = %d ; n_vtx = %d ; \n", id_isosurface, i_part, n_vtx);
      vtx_field[i_part] = malloc(sizeof(double) * n_vtx); 

      for (int i_vtx=0; i_vtx<n_vtx; ++i_vtx) {
        
        double x = vtx_coord[3*i_vtx  ];
        double y = vtx_coord[3*i_vtx+1];
        double z = vtx_coord[3*i_vtx+2];

        if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_PLANE   ) {
          vtx_field[i_part][i_vtx] = _plane_field(x, y, z, isos->eq_coeffs[id_isosurface]);
        }

        else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_SPHERE  ) {
          vtx_field[i_part][i_vtx] = _sphere_field(x, y, z, isos->eq_coeffs[id_isosurface]);
        }

        else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_ELLIPSE ) {
          vtx_field[i_part][i_vtx] = _ellipse_field(x, y, z, isos->eq_coeffs[id_isosurface]);
        }

        else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_QUADRIC ) {
          vtx_field[i_part][i_vtx] = _quadric_field(x, y, z, isos->eq_coeffs[id_isosurface]);
        }

        else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_HEART   ) {
          vtx_field[i_part][i_vtx] = _heart_field(x, y, z, isos->eq_coeffs[id_isosurface]);
        }

        else if (isos->kind[id_isosurface]==PDM_ISO_SURFACE_KIND_FUNCTION) {
          isos->iso_func[id_isosurface](x, y, z, &vtx_field[i_part][i_vtx]);
        }
      }
    }

    if (debug==1) {
      PDM_log_trace_array_double(vtx_field[i_part], n_vtx, "_compute_iso_field:: vtx_field [i_part]::");
    }
  }
}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/



static void
_free_iso_vtx
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_vtx_coord[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_vtx_coord[id_iso][i_part]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_coord[id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_gnum[id_iso][i_part]);
      }
      free(isos->iso_vtx_lparent_idx   [id_iso][i_part]);
      free(isos->iso_vtx_lparent       [id_iso][i_part]);
      if (isos->iso_owner_vtx_parent_weight[id_iso][i_part]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_vtx_parent_weight[id_iso][i_part]);
      }
      free(isos->isovalue_vtx_idx[id_iso][i_part]);
    }
  }
  if (isos->iso_n_vtx[id_iso]!=NULL) {
    free(isos->iso_n_vtx            [id_iso]);
  }
  if (isos->iso_vtx_coord[id_iso]!=NULL) {
    free(isos->iso_vtx_coord[id_iso]);
  }
  if (isos->iso_vtx_gnum[id_iso]!=NULL) {
    free(isos->iso_vtx_gnum[id_iso]);
  }
  if (isos->iso_vtx_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_vtx_lparent_idx[id_iso]);
  }
  if (isos->iso_vtx_lparent[id_iso]!=NULL) {
    free(isos->iso_vtx_lparent[id_iso]);
  }
  if (isos->iso_vtx_parent_weight[id_iso]!=NULL) {
    free(isos->iso_vtx_parent_weight[id_iso]);
  }
  if (isos->isovalue_vtx_idx[id_iso]!=NULL) {
    free(isos->isovalue_vtx_idx[id_iso]);
  }

  isos->iso_n_vtx            [id_iso] = NULL;
  isos->iso_vtx_coord        [id_iso] = NULL;
  isos->iso_vtx_gnum         [id_iso] = NULL;
  isos->iso_vtx_lparent_idx  [id_iso] = NULL;
  isos->iso_vtx_lparent      [id_iso] = NULL;
  isos->iso_vtx_parent_weight[id_iso] = NULL;
  isos->isovalue_vtx_idx     [id_iso] = NULL;
}


static void
_free_iso_edge
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_connec[id_iso][i_part][PDM_CONNECTIVITY_TYPE_EDGE_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_edge_vtx[id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_EDGE]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_edge_gnum[id_iso][i_part]);
      }
      free(isos->iso_edge_lparent_idx[id_iso][i_part]);
      free(isos->iso_edge_lparent    [id_iso][i_part]);
      free(isos->iso_edge_group_idx  [id_iso][i_part]);
      free(isos->iso_edge_group_lnum [id_iso][i_part]);
      free(isos->iso_edge_group_gnum [id_iso][i_part]);
      free(isos->isovalue_edge_idx   [id_iso][i_part]);
    }
  }
  if (isos->iso_n_edge[id_iso]!=NULL) {
    free(isos->iso_n_edge[id_iso]);
  }
  if (isos->iso_edge_vtx[id_iso]!=NULL) {
    free(isos->iso_edge_vtx[id_iso]);
  }
  if (isos->iso_edge_gnum[id_iso]!=NULL) {
    free(isos->iso_edge_gnum[id_iso]);
  }
  if (isos->iso_edge_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_edge_lparent_idx[id_iso]);
  }
  if (isos->iso_edge_lparent[id_iso]!=NULL) {
    free(isos->iso_edge_lparent[id_iso]);
  }
  if (isos->iso_edge_group_idx[id_iso]!=NULL) {
    free(isos->iso_edge_group_idx[id_iso]);
  }
  if (isos->iso_edge_group_lnum[id_iso]!=NULL) {
    free(isos->iso_edge_group_lnum[id_iso]);
  }
  if (isos->iso_edge_group_gnum[id_iso]!=NULL) {
    free(isos->iso_edge_group_gnum[id_iso]);
  }
  if (isos->isovalue_edge_idx[id_iso]!=NULL) {
    free(isos->isovalue_edge_idx[id_iso]);
  }

  isos->iso_n_edge          [id_iso] = NULL;
  isos->iso_edge_vtx        [id_iso] = NULL;
  isos->iso_edge_gnum       [id_iso] = NULL;
  isos->iso_edge_lparent_idx[id_iso] = NULL;
  isos->iso_edge_lparent    [id_iso] = NULL;
  isos->iso_n_edge_group    [id_iso] = 0;
  isos->iso_edge_group_idx  [id_iso] = NULL;
  isos->iso_edge_group_lnum [id_iso] = NULL;
  isos->iso_edge_group_gnum [id_iso] = NULL;
  isos->isovalue_edge_idx   [id_iso] = NULL;
}


static void
_free_iso_face
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      if (isos->iso_owner_connec[id_iso][i_part][PDM_CONNECTIVITY_TYPE_FACE_VTX]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_face_vtx_idx   [id_iso][i_part]);
        free(isos->iso_face_vtx       [id_iso][i_part]);
      }
      if (isos->iso_owner_gnum[id_iso][i_part][PDM_MESH_ENTITY_FACE]==PDM_OWNERSHIP_KEEP) {
        free(isos->iso_face_gnum      [id_iso][i_part]);
      }
      free(isos->iso_face_lparent_idx[id_iso][i_part]);
      free(isos->iso_face_lparent    [id_iso][i_part]);
      free(isos->isovalue_face_idx   [id_iso][i_part]);
    }
  }

  if (isos->iso_n_face[id_iso]!=NULL) {
    free(isos->iso_n_face[id_iso]);
  }
  if (isos->iso_face_gnum[id_iso]!=NULL) {
    free(isos->iso_face_gnum[id_iso]);
  }
  if (isos->iso_face_vtx_idx[id_iso]!=NULL) {
    free(isos->iso_face_vtx_idx[id_iso]);
  }
  if (isos->iso_face_vtx[id_iso]!=NULL) {
    free(isos->iso_face_vtx[id_iso]);
  }
  if (isos->iso_face_lparent_idx[id_iso]!=NULL) {
    free(isos->iso_face_lparent_idx[id_iso]);
  }
  if (isos->iso_face_lparent[id_iso]!=NULL) {
    free(isos->iso_face_lparent[id_iso]);
  }
  if (isos->isovalue_face_idx[id_iso]!=NULL) {
    free(isos->isovalue_face_idx[id_iso]);
  }

  isos->iso_n_face          [id_iso] = NULL;
  isos->iso_face_vtx_idx    [id_iso] = NULL;
  isos->iso_face_vtx        [id_iso] = NULL;
  isos->iso_face_gnum       [id_iso] = NULL;
  isos->iso_face_lparent_idx[id_iso] = NULL;
  isos->iso_face_lparent    [id_iso] = NULL;
  isos->isovalue_face_idx   [id_iso] = NULL;
}


static void
_free_field
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  if (isos->field[id_iso]!=NULL) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) { // TODO: fix this n_part (volumic)
      if (isos->kind[id_iso]!=PDM_ISO_SURFACE_KIND_FIELD && isos->field[id_iso][i_part]!=NULL) {
        free(isos->field[id_iso][i_part]);
      }
    }
    free(isos->field[id_iso]);
  }
  isos->field[id_iso] = NULL;
}


static void
_free_owner
(
  PDM_isosurface_t *isos,
  int               id_iso
)
{
  free(isos->iso_owner_vtx_coord        [id_iso]);
  free(isos->iso_owner_vtx_parent_weight[id_iso]);
  for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
    if (isos->iso_owner_gnum[id_iso]!=NULL) {
      free(isos->iso_owner_gnum   [id_iso][i_part]);
    }
    if (isos->iso_owner_connec[id_iso]!=NULL) {
      free(isos->iso_owner_connec [id_iso][i_part]);
    }
    if (isos->iso_owner_lparent[id_iso]!=NULL) {
      free(isos->iso_owner_lparent[id_iso][i_part]);
    }
  }
  if (isos->iso_owner_gnum[id_iso]!=NULL) {
    free(isos->iso_owner_gnum   [id_iso]);
  }
  if (isos->iso_owner_connec[id_iso]!=NULL) {
    free(isos->iso_owner_connec [id_iso]);
  }
  if (isos->iso_owner_lparent[id_iso]!=NULL) {
    free(isos->iso_owner_lparent[id_iso]);
  }
  if (isos->iso_owner_edge_bnd[id_iso]!=NULL) {
    free(isos->iso_owner_edge_bnd[id_iso]);
  }
  isos->iso_owner_vtx_coord        [id_iso] = NULL;
  isos->iso_owner_vtx_parent_weight[id_iso] = NULL;
  isos->iso_owner_gnum             [id_iso] = NULL;
  isos->iso_owner_connec           [id_iso] = NULL;
  isos->iso_owner_lparent          [id_iso] = NULL;
  isos->iso_owner_edge_bnd         [id_iso] = NULL;
  isos->iso_owner_ptp              [id_iso] = NULL;
}


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
 int                      mesh_dimension
 // PDM_Mesh_nodal_elt_t     elt_type
)
{
  PDM_isosurface_t *isos = (PDM_isosurface_t *) malloc(sizeof(PDM_isosurface_t));

  // > Save communicator
  isos->comm = comm;

  // > Entry mesh information
  isos->is_dist_or_part = -1; 
  isos->entry_mesh_type =  0; 
  isos->entry_mesh_dim  =  mesh_dimension;

  // // > Isosurface mesh information
  // isos->iso_elt_type = elt_type; 

  // > Link with entry mesh
  isos->compute_ptp = NULL;

  // > Isosurfaces informations
  isos->n_isosurface = 0;
  isos->n_isovalues     = NULL;
  isos->  isovalues     = NULL;
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

  isos->n_part = -1;
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


  // > Partitioned output data
  isos->iso_n_vtx             = NULL;
  isos->iso_vtx_coord         = NULL;
  isos->iso_vtx_gnum          = NULL;
  isos->iso_vtx_lparent_idx   = NULL;
  isos->iso_vtx_lparent       = NULL;
  isos->iso_vtx_parent_weight = NULL;
  isos->isovalue_vtx_idx      = NULL;

  isos->iso_n_edge            = NULL;
  isos->iso_edge_vtx          = NULL;
  isos->iso_edge_gnum         = NULL;
  isos->iso_edge_lparent_idx  = NULL;
  isos->iso_edge_lparent      = NULL;
  isos->iso_n_edge_group      = NULL;
  isos->iso_edge_group_idx    = NULL;
  isos->iso_edge_group_lnum   = NULL;
  isos->iso_edge_group_gnum   = NULL;
  isos->isovalue_edge_idx     = NULL;

  isos->iso_n_face            = NULL;
  isos->iso_face_vtx_idx      = NULL;
  isos->iso_face_vtx          = NULL;
  isos->iso_face_gnum         = NULL;
  isos->iso_face_lparent_idx  = NULL;
  isos->iso_face_lparent      = NULL;
  isos->isovalue_face_idx     = NULL;

  // > Part_to_part between iso entities and entry mesh entites
  isos->iso_ptp_vtx  = NULL;
  isos->iso_ptp_edge = NULL;
  isos->iso_ptp_face = NULL;

  // > Owners
  isos->iso_owner_vtx_coord         = NULL;
  isos->iso_owner_vtx_parent_weight = NULL;
  isos->iso_owner_gnum              = NULL;
  isos->iso_owner_connec            = NULL;
  isos->iso_owner_lparent           = NULL;
  isos->iso_owner_edge_bnd          = NULL;
  isos->iso_owner_ptp               = NULL;


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
  // TODO: check that difference between isovalues>ISOSURFACE_EPS
  
  int id_isosurface = isos->n_isosurface;
  isos->n_isosurface++;
  
  printf("isos->n_isosurface = %d\n", isos->n_isosurface);
  isos->kind           = realloc(isos->kind          , sizeof(PDM_iso_surface_kind_t          ) * isos->n_isosurface);
  isos->eq_coeffs      = realloc(isos->eq_coeffs     , sizeof(double                         *) * isos->n_isosurface);
  isos->field_function = realloc(isos->field_function, sizeof(PDM_isosurface_field_function_t ) * isos->n_isosurface);
  
  isos->field_function = realloc(isos->field_function, sizeof(PDM_isosurface_field_function_t ) * isos->n_isosurface);
  
  isos->compute_ptp = realloc(isos->compute_ptp, sizeof(int *) * isos->n_isosurface);
  isos->compute_ptp[id_isosurface] = PDM_array_zeros_int(PDM_MESH_ENTITY_MAX);

  // TODO: status if its ok to allocate d and p variable (because set functions order not imposed ?)
  isos->dfield         = realloc(isos->dfield        , sizeof(double                         *) * isos->n_isosurface);
  isos->field          = realloc(isos->field         , sizeof(double                        **) * isos->n_isosurface);
  isos->field[id_isosurface] = NULL;
  if (kind==PDM_ISO_SURFACE_KIND_FIELD) {
    int n_part = 0;
    if (isos->entry_mesh_type==-1) {
      n_part = isos->n_part;
      if (n_part==-1) {
        PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: entry_mesh_type = %d but n_part isn't defined.\n", isos->entry_mesh_type);
      }
    }
    else if (isos->entry_mesh_type==-2) {
      n_part = PDM_part_mesh_n_part_get(isos->pmesh);
    }
    else if (isos->entry_mesh_type==-3) {
      n_part = PDM_part_mesh_nodal_n_part_get(isos->pmesh_nodal);
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_t: Impossible to defined manually field without setting mesh first.\n", isos->entry_mesh_type);
    }
    isos->field[id_isosurface] = malloc(sizeof(double *) * n_part);
  }

  isos->n_isovalues    = realloc(isos->n_isovalues   , sizeof(int                             ) * isos->n_isosurface);
  isos->isovalues      = realloc(isos->isovalues     , sizeof(double                         *) * isos->n_isosurface);
  isos->use_gradient   = realloc(isos->use_gradient  , sizeof(int                            *) * isos->n_isosurface);

  // > Partitioned iso vertices
  isos->iso_n_vtx             = realloc(isos->iso_n_vtx            , sizeof(int           *) * isos->n_isosurface);
  isos->iso_vtx_coord         = realloc(isos->iso_vtx_coord        , sizeof(double       **) * isos->n_isosurface);
  isos->iso_vtx_gnum          = realloc(isos->iso_vtx_gnum         , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_vtx_lparent_idx   = realloc(isos->iso_vtx_lparent_idx  , sizeof(int          **) * isos->n_isosurface);
  isos->iso_vtx_lparent       = realloc(isos->iso_vtx_lparent      , sizeof(int          **) * isos->n_isosurface);
  isos->iso_vtx_parent_weight = realloc(isos->iso_vtx_parent_weight, sizeof(double       **) * isos->n_isosurface);
  isos->isovalue_vtx_idx      = realloc(isos->isovalue_vtx_idx     , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_vtx            [id_isosurface] = NULL;
  isos->iso_vtx_coord        [id_isosurface] = NULL;
  isos->iso_vtx_gnum         [id_isosurface] = NULL;
  isos->iso_vtx_lparent_idx  [id_isosurface] = NULL;
  isos->iso_vtx_lparent      [id_isosurface] = NULL;
  isos->iso_vtx_parent_weight[id_isosurface] = NULL;
  isos->isovalue_vtx_idx     [id_isosurface] = NULL;

  // > Partitioned iso edges
  isos->iso_n_edge           = realloc(isos->iso_n_edge          , sizeof(int           *) * isos->n_isosurface);
  isos->iso_edge_vtx         = realloc(isos->iso_edge_vtx        , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_gnum        = realloc(isos->iso_edge_gnum       , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_edge_lparent_idx = realloc(isos->iso_edge_lparent_idx, sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_lparent     = realloc(isos->iso_edge_lparent    , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_edge_group     = realloc(isos->iso_n_edge_group    , sizeof(int            ) * isos->n_isosurface);
  isos->iso_edge_group_idx   = realloc(isos->iso_edge_group_idx  , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_group_lnum  = realloc(isos->iso_edge_group_lnum , sizeof(int          **) * isos->n_isosurface);
  isos->iso_edge_group_gnum  = realloc(isos->iso_edge_group_gnum , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->isovalue_edge_idx    = realloc(isos->isovalue_edge_idx   , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_edge          [id_isosurface] = NULL;
  isos->iso_edge_vtx        [id_isosurface] = NULL;
  isos->iso_edge_gnum       [id_isosurface] = NULL;
  isos->iso_edge_lparent_idx[id_isosurface] = NULL;
  isos->iso_edge_lparent    [id_isosurface] = NULL;
  isos->iso_n_edge_group    [id_isosurface] = 0;
  isos->iso_edge_group_idx  [id_isosurface] = NULL;
  isos->iso_edge_group_lnum [id_isosurface] = NULL;
  isos->iso_edge_group_gnum [id_isosurface] = NULL;
  isos->isovalue_edge_idx   [id_isosurface] = NULL;

  // > Partitioned iso faces
  isos->iso_n_face           = realloc(isos->iso_n_face          , sizeof(int           *) * isos->n_isosurface);
  isos->iso_face_vtx_idx     = realloc(isos->iso_face_vtx_idx    , sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_vtx         = realloc(isos->iso_face_vtx        , sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_gnum        = realloc(isos->iso_face_gnum       , sizeof(PDM_g_num_t  **) * isos->n_isosurface);
  isos->iso_face_lparent_idx = realloc(isos->iso_face_lparent_idx, sizeof(int          **) * isos->n_isosurface);
  isos->iso_face_lparent     = realloc(isos->iso_face_lparent    , sizeof(int          **) * isos->n_isosurface);
  isos->isovalue_face_idx    = realloc(isos->isovalue_face_idx   , sizeof(int          **) * isos->n_isosurface);
  isos->iso_n_face          [id_isosurface] = NULL;
  isos->iso_face_vtx_idx    [id_isosurface] = NULL;
  isos->iso_face_vtx        [id_isosurface] = NULL;
  isos->iso_face_gnum       [id_isosurface] = NULL;
  isos->iso_face_lparent_idx[id_isosurface] = NULL;
  isos->iso_face_lparent    [id_isosurface] = NULL;
  isos->isovalue_face_idx   [id_isosurface] = NULL;

  // > Part_to_part between iso entities and entry mesh entites
  isos->iso_ptp_vtx  = realloc(isos->iso_ptp_vtx , sizeof(PDM_part_to_part_t *) * isos->n_isosurface);
  isos->iso_ptp_edge = realloc(isos->iso_ptp_edge, sizeof(PDM_part_to_part_t *) * isos->n_isosurface);
  isos->iso_ptp_face = realloc(isos->iso_ptp_face, sizeof(PDM_part_to_part_t *) * isos->n_isosurface);

  // > Partitioned owners
  isos->iso_owner_vtx_coord         = realloc(isos->iso_owner_vtx_coord         , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_vtx_parent_weight = realloc(isos->iso_owner_vtx_parent_weight , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_gnum              = realloc(isos->iso_owner_gnum              , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_connec            = realloc(isos->iso_owner_connec            , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_lparent           = realloc(isos->iso_owner_lparent           , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_edge_bnd          = realloc(isos->iso_owner_edge_bnd          , sizeof(PDM_ownership_t **) * isos->n_isosurface);
  isos->iso_owner_ptp               = realloc(isos->iso_owner_ptp               , sizeof(PDM_ownership_t  *) * isos->n_isosurface);
  isos->iso_owner_vtx_coord        [id_isosurface] = NULL;
  isos->iso_owner_vtx_parent_weight[id_isosurface] = NULL;
  isos->iso_owner_gnum             [id_isosurface] = NULL;
  isos->iso_owner_connec           [id_isosurface] = NULL;
  isos->iso_owner_lparent          [id_isosurface] = NULL;
  isos->iso_owner_edge_bnd         [id_isosurface] = NULL;
  isos->iso_owner_ptp              [id_isosurface] = NULL;

  isos->kind       [id_isosurface] = kind;
  isos->n_isovalues[id_isosurface] = n_isovalues;
  isos->  isovalues[id_isosurface] = malloc(sizeof(double) * n_isovalues);
  for (int i=0; i<n_isovalues; ++i) {
    isos->isovalues[id_isosurface][i] = isovalues[i];
  }

  if (kind==PDM_ISO_SURFACE_KIND_FIELD) {
    isos->eq_coeffs[id_isosurface] = NULL;
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
  // > Distributed
  if (isos->is_dist_or_part==0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_isosurface_reset not implemented for dist_entry\n");
  }
  // > Partitioned
  else if (isos->is_dist_or_part==1) {
    for (int i_part=0; i_part<isos->iso_n_part; ++i_part) {
      _free_iso_vtx (isos, id_isosurface);
      _free_iso_edge(isos, id_isosurface);
      _free_iso_face(isos, id_isosurface);
      _free_owner(isos, id_isosurface);
      _free_field(isos, id_isosurface);
    }
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d is invalid.\n", isos->is_dist_or_part);
  }
}


void
PDM_isosurface_compute
(
 PDM_isosurface_t *isos,
 int               id_isosurface
)
{
  int debug = 1;

  if (debug==1) {
    log_trace("PDM_isosurface:: compute isosurface n°%d\n", id_isosurface);
  }

  if (isos->is_dist_or_part==0) { // Distributed entry
    if (isos->entry_mesh_type<0) {
      PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d incoherent with isos->entry_mesh_type = %d < 0.\n", isos->is_dist_or_part, isos->entry_mesh_type);
    } else if (isos->entry_mesh_type==1) { // Dist mesh alamano

    } else if (isos->entry_mesh_type==2) { // Dist mesh

    } else if (isos->entry_mesh_type==3) { // Dist mesh nodal

    } else {
      PDM_error(__FILE__, __LINE__, 0, "Isosurface isos->entry_mesh_type = %d is invalid for distributed entry.\n", isos->entry_mesh_type);
    }
  } else if (isos->is_dist_or_part==1) { // Partitioned entry
    if (isos->entry_mesh_type>0) {
      PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d incoherent with isos->entry_mesh_type = %d > 0.\n", isos->is_dist_or_part, isos->entry_mesh_type);
    } else if (isos->entry_mesh_type==-1) { // Part mesh alamano

    } else if (isos->entry_mesh_type==-2) { // Part mesh

    } else if (isos->entry_mesh_type==-3) { // Part mesh nodal

      if (debug == 1) {
        int mesh_dim = PDM_part_mesh_nodal_mesh_dimension_get(isos->pmesh_nodal);
        if (mesh_dim>1) {
          PDM_part_mesh_nodal_dump_vtk(isos->pmesh_nodal,
                                       PDM_GEOMETRY_KIND_SURFACIC,
                                       "pmn_surfacic_entry_mesh");
        }
        if (mesh_dim>2) {
          PDM_part_mesh_nodal_dump_vtk(isos->pmesh_nodal,
                                       PDM_GEOMETRY_KIND_VOLUMIC,
                                       "pmn_volumic_entry_mesh");
        }
      }

      /*
       * TODO: Extract part sur les cellules
       */

      
      /*
       * Compute field once for all to avoid multiple if in code
       * and if user function is costful, would be costful once
       * TODO: compute field on extract part mesh
       */

      _compute_iso_field(isos, id_isosurface);

      PDM_isosurface_marching_algo(isos,
                                   id_isosurface);

    } else {
      PDM_error(__FILE__, __LINE__, 0, "Isosurface isos->entry_mesh_type = %d is invalid for partitioned entry.\n", isos->entry_mesh_type);
    }
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Isosurface is_dist_or_part = %d is invalid.\n", isos->is_dist_or_part);
  }

}


void
PDM_isosurface_dump_times
(
 PDM_isosurface_t *isos
)
{
  PDM_UNUSED(isos);
}


void
PDM_isosurface_free
(
  PDM_isosurface_t  *isos
)
{
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->n_isovalues[id_iso]>0) {
      free(isos->isovalues[id_iso]);
      if (isos->kind[id_iso]==PDM_ISO_SURFACE_KIND_FIELD) {
        if (isos->is_dist_or_part==0) { // is distributed
          free(isos->dfield[id_iso]);
        }
      } else  {
        free(isos->eq_coeffs[id_iso]);
      }
    }
  }

  free(isos->n_isovalues);
  free(isos->isovalues);
  free(isos->eq_coeffs);
  free(isos->use_gradient);
  free(isos->field_function);
  free(isos->dfield);
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_field(isos, id_iso);
  }
  free(isos->field);
  free(isos->kind);

  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_iso_vtx (isos, id_iso);
    _free_iso_edge(isos, id_iso);
    _free_iso_face(isos, id_iso);
  }

  free(isos->iso_n_vtx);
  free(isos->iso_vtx_coord);
  free(isos->iso_vtx_gnum);
  free(isos->iso_vtx_lparent_idx);
  free(isos->iso_vtx_lparent);
  free(isos->iso_vtx_parent_weight);
  free(isos->isovalue_vtx_idx);
  
  free(isos->iso_n_edge);
  free(isos->iso_edge_vtx);
  free(isos->iso_edge_gnum);
  free(isos->iso_edge_lparent_idx);
  free(isos->iso_edge_lparent);
  free(isos->iso_n_edge_group);
  free(isos->iso_edge_group_idx);
  free(isos->iso_edge_group_lnum);
  free(isos->iso_edge_group_gnum);
  free(isos->isovalue_edge_idx);

  free(isos->iso_n_face);
  free(isos->iso_face_vtx_idx);
  free(isos->iso_face_vtx);
  free(isos->iso_face_gnum);
  free(isos->iso_face_lparent_idx);
  free(isos->iso_face_lparent);
  free(isos->isovalue_face_idx);


  // > Part_to_part between iso entities and entry mesh entites
  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_VTX]==1 &&
        isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_VTX]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_vtx[id_iso]);
    }
    else if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_EDGE]==1 &&
             isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_EDGE]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_edge[id_iso]);
    }
    else if (isos->  compute_ptp[id_iso][PDM_MESH_ENTITY_FACE]==1 &&
             isos->iso_owner_ptp[id_iso][PDM_MESH_ENTITY_FACE]==PDM_OWNERSHIP_KEEP) {
      PDM_part_to_part_free(isos->iso_ptp_face[id_iso]);
    }
  }

  free(isos->iso_ptp_vtx);
  free(isos->iso_ptp_edge);
  free(isos->iso_ptp_face);



  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    _free_owner(isos, id_iso);
  }
  free(isos->iso_owner_vtx_coord);
  free(isos->iso_owner_vtx_parent_weight);
  free(isos->iso_owner_gnum);
  free(isos->iso_owner_connec);
  free(isos->iso_owner_lparent);
  free(isos->iso_owner_edge_bnd);
  free(isos->iso_owner_ptp);



  for (int id_iso=0; id_iso<isos->n_isosurface; ++id_iso) {
    free(isos->compute_ptp[id_iso]);
  }
  free(isos->compute_ptp);


  free(isos);
}


#ifdef  __cplusplus
}
#endif