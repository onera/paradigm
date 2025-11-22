/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_surf_mesh.h"
#include "pdm_array.h"
#include "pdm_binary_search.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_part_to_block.h"
#include "pdm_priv.h"
#include "pdm_surf_mesh_priv.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief Compute global number of entities
 *
 * \param [in]  mesh         Mesh 
 *
 */

static 
void
_n_g_enttities_compute
(
PDM_surf_mesh_t *mesh
)
{

  if (mesh->nGVtx < 0) {

    PDM_g_num_t n_g_num_vtx = 0;

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_surf_part_t *part = mesh->part[i];
      for (int j = 0; j < part->n_vtx; j++) {
        n_g_num_vtx = PDM_MAX (n_g_num_vtx, part->vtx_ln_to_gn[j]);
      }
    }

    PDM_MPI_Allreduce(&n_g_num_vtx, &(mesh->nGVtx), 1,
                      PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, mesh->comm);
  }

  if (mesh->nGFace < 0) {

    PDM_g_num_t n_g_num_face = 0;

    for (int i = 0; i < mesh->n_part; i++) {
      PDM_surf_part_t *part = mesh->part[i];
      for (int j = 0; j < part->n_face; j++) {
        n_g_num_face = PDM_MAX (n_g_num_face, part->face_ln_to_gn[j]);
      }
    }

    PDM_MPI_Allreduce(&n_g_num_face, &(mesh->nGFace), 1,
                      PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, mesh->comm);
  }
 
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


PDM_surf_mesh_t *
PDM_surf_mesh_create
(
const int    n_part,
PDM_MPI_Comm comm
)
{

  PDM_surf_mesh_t *mesh;
  PDM_malloc(mesh,1,PDM_surf_mesh_t);

  mesh->nGFace  = -1;
  mesh->nGVtx   = -1;

  mesh->comm    = comm;
  mesh->n_part   = n_part;
  PDM_malloc(mesh->part,n_part ,PDM_surf_part_t *);

  PDM_MPI_Allreduce ((void *)&n_part, (void *)&(mesh->nGPart), 1,
                     PDM_MPI_INT, PDM_MPI_SUM, mesh->comm);

  mesh->gMinCarLgthVtx =  DBL_MAX;
  mesh->gMaxCarLgthVtx = -DBL_MAX;

  return (PDM_surf_mesh_t *) mesh;
}


PDM_surf_mesh_t *
PDM_surf_mesh_free
(
PDM_surf_mesh_t *mesh
)
{
  if (mesh != NULL) {

    if (mesh->part != NULL) {
      for (int i = 0; i < mesh->n_part; i++) {
        mesh->part[i] = PDM_surf_part_free(mesh->part[i]);
      }
      PDM_free(mesh->part);
      mesh->part = NULL;
    }

    PDM_free(mesh);
  }

  return NULL;
}

void
PDM_surf_mesh_part_input
(
 PDM_surf_mesh_t      *mesh,
 const int            i_part,
 const int            n_face,
 const int           *face_vtx_idx,
 const int           *face_vtx,
 const PDM_g_num_t   *face_ln_to_gn,
 const int            n_vtx,
 const double        *coords,
 const PDM_g_num_t   *vtx_ln_to_gn
)
{
  assert (mesh != NULL);

  mesh->part[i_part] = PDM_surf_part_create(n_face,
                                            face_vtx_idx,
                                            face_vtx,
                                            face_ln_to_gn,
                                            n_vtx,
                                            coords,
                                            vtx_ln_to_gn);
}


void
PDM_surf_mesh_compute_faceExtentsMesh
(
 PDM_surf_mesh_t *mesh,
 double           tolerance
)
{
  assert (mesh != NULL);

  for (int i = 0; i < mesh->n_part; i++) {

    PDM_surf_part_t *part = mesh->part[i];
    const int n_face = part->n_face;
    const int *face_vtx = part->face_vtx;
    const int *face_vtx_idx = part->face_vtx_idx;
    const double *coords = part->coords;

    PDM_malloc(part->extents,6 * n_face,double);

    /*
     * TODO : Optimization : Split this boucle by blocks (vectorization)
     */

    double *_extents = part->extents;
    for (int j = 0; j < n_face; j++) {

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[k1]   = DBL_MAX;
        _extents[3+k1] = -DBL_MAX;
      }

      for (int k = face_vtx_idx[j]; k < face_vtx_idx[j+1]; k++) {
        int iVtx = face_vtx[k] - 1;
        double *_coords = (double *) coords + 3 * iVtx;

        for (int k1 = 0; k1 < 3; k1++) {
          _extents[k1]   = PDM_MIN (_coords[k1], _extents[k1]);
          _extents[3+k1] = PDM_MAX (_coords[k1], _extents[3+k1]);
        }

      }

      double delta = -DBL_MAX;

      for (int k1 = 0; k1 < 3; k1++) {
        delta = PDM_MAX (delta, fabs (_extents[k1+3] - _extents[k1]));
      }

      delta *= tolerance;

      for (int k1 = 0; k1 < 3; k1++) {
        _extents[k1]   +=  - delta;
        _extents[3+k1] +=    delta;
      }

      _extents += 6;
    }
  }
}

int
PDM_surf_mesh_n_part_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);

  return mesh->n_part;
}


int
PDM_surf_mesh_part_n_face_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->n_face;
}


int
PDM_surf_mesh_part_n_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->n_vtx;
}


const double *
PDM_surf_mesh_part_extents_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);


  PDM_surf_part_t *part  =  mesh->part[i_part];

  return part->extents;
}


const PDM_g_num_t *
PDM_surf_mesh_part_face_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->face_ln_to_gn;
}


const PDM_g_num_t *
PDM_surf_mesh_part_vtx_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->vtx_ln_to_gn;
}


const PDM_g_num_t *
PDM_surf_mesh_part_edge_g_num_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->edgeLnToGn;
}


const int *
PDM_surf_mesh_part_face_edge_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->faceEdge;
}



const int *
PDM_surf_mesh_part_face_vtx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->face_vtx;
}


const int *
PDM_surf_mesh_part_face_vtx_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  =  mesh->part[i_part];

  return part->face_vtx_idx;
}

const int *
PDM_surf_mesh_part_face_edge_idx_get
(
 PDM_surf_mesh_t      *mesh,
 int                   i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->faceEdgeIdx;
}



const double *
PDM_surf_mesh_part_vtx_get
(
 PDM_surf_mesh_t *mesh,
 int              i_part
)
{
  assert (mesh != NULL);

  PDM_surf_part_t *part  = mesh->part[i_part];

  return part->coords;
}


PDM_g_num_t
PDM_surf_mesh_n_g_vtx_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);
  
  _n_g_enttities_compute(mesh);

  return mesh->nGVtx;
}


PDM_g_num_t
PDM_surf_mesh_n_g_face_get
(
 PDM_surf_mesh_t      *mesh
)
{
  assert (mesh != NULL);

  _n_g_enttities_compute(mesh);

  return mesh->nGFace;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
