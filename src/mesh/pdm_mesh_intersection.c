/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mesh_intersection_priv.h"
#include "pdm_mesh_intersection.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"
#include "pdm_extract_part.h"
#include "pdm_extract_part_priv.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"
#include "pdm_writer.h"
#include "pdm_gnum.h"
#include "pdm_box_priv.h"
#include "pdm_unique.h"
#include "pdm_triangulate.h"
#include "pdm_geom_elem.h"

#include "pdm_mesh_intersection_vol_vol_atomic.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_export_vtk_1d
(
 const char               *pattern,
       PDM_extract_part_t *extrp_mesh
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp_mesh->comm, &i_rank);


  for(int i_part = 0; i_part < extrp_mesh->n_part_out; ++i_part) {

    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;
    int n_edge = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_EDGE  , &edge_ln_to_gn, PDM_OWNERSHIP_KEEP);
    int n_vtx  = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_VERTEX, &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp_mesh, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

    int  *edge_vtx      = NULL;
    int  *edge_vtx_idx  = NULL;
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_KEEP);

    char filename[999];
    sprintf(filename, "%s_%i_%i.vtk", pattern, i_part, i_rank);
    PDM_vtk_write_std_elements(filename,
                               n_vtx,
                               vtx_coord,
                               vtx_ln_to_gn,
                               PDM_MESH_NODAL_BAR2,
                               n_edge,
                               edge_vtx,
                               edge_ln_to_gn,
                               0,
                               NULL,
                               NULL);
  }
}

static
void
_export_vtk_2d
(
 const char               *pattern,
       PDM_extract_part_t *extrp_mesh
)
{
  int i_rank;
  PDM_MPI_Comm_rank(extrp_mesh->comm, &i_rank);

  for(int i_part = 0; i_part < extrp_mesh->n_part_out; ++i_part) {

    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;
    int n_face = PDM_extract_part_parent_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_FACE  , &face_ln_to_gn, PDM_OWNERSHIP_KEEP);
    int n_vtx  = PDM_extract_part_ln_to_gn_get(extrp_mesh, i_part, PDM_MESH_ENTITY_VERTEX, &vtx_ln_to_gn , PDM_OWNERSHIP_KEEP);

    double *vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp_mesh, i_part, &vtx_coord, PDM_OWNERSHIP_KEEP);

    int  *face_edge     = NULL;
    int  *face_edge_idx = NULL;
    int  *edge_vtx      = NULL;
    int  *edge_vtx_idx  = NULL;
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_KEEP);
    PDM_extract_part_connectivity_get(extrp_mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_KEEP);

    int *face_vtx = NULL;
    PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);

    char filename[999];
    sprintf(filename, "%s_%i_%i.vtk", pattern, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           vtx_coord,
                           vtx_ln_to_gn,
                           n_face,
                           face_edge_idx,
                           face_vtx,
                           face_ln_to_gn,
                           NULL);


    free(face_vtx);
  }
}

static void
_export_vtk_3d
(
 const char               *name_chr,
       PDM_extract_part_t *extrp
)
{
  int i_rank;

  PDM_MPI_Comm_rank(extrp->comm, &i_rank);

  int          *pn_extract_cell        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_face        = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int          *pn_extract_vtx         = (int          *) malloc(extrp->n_part_out * sizeof(int          ));
  int         **pextract_cell_face     = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_cell_face_idx = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx      = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  int         **pextract_face_vtx_idx  = (int         **) malloc(extrp->n_part_out * sizeof(int         *));
  double      **pextract_vtx           = (double      **) malloc(extrp->n_part_out * sizeof(double      *));
  PDM_g_num_t **pextract_cell_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_face_ln_to_gn = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));
  PDM_g_num_t **pextract_vtx_ln_to_gn  = (PDM_g_num_t **) malloc(extrp->n_part_out * sizeof(PDM_g_num_t *));


  for(int i_part = 0; i_part < extrp->n_part_out; ++i_part) {

    pn_extract_cell[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_CELL);

    pn_extract_face[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                            i_part,
                                                            PDM_MESH_ENTITY_FACE);

    pn_extract_vtx[i_part] = PDM_extract_part_n_entity_get(extrp,
                                                           i_part,
                                                           PDM_MESH_ENTITY_VERTEX);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                      &pextract_cell_face[i_part],
                                      &pextract_cell_face_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                      &pextract_face_vtx[i_part],
                                      &pextract_face_vtx_idx[i_part],
                                      PDM_OWNERSHIP_KEEP);

    PDM_ownership_t owner = PDM_OWNERSHIP_KEEP;
    if(pextract_face_vtx[i_part] == NULL) {
      owner = PDM_OWNERSHIP_KEEP;

      int *face_edge     = NULL;
      int *face_edge_idx = NULL;
      int *edge_vtx     = NULL;
      int *edge_vtx_idx = NULL;
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &face_edge,
                                        &face_edge_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_extract_part_connectivity_get(extrp,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx,
                                        &edge_vtx_idx,
                                        PDM_OWNERSHIP_KEEP);
      PDM_compute_face_vtx_from_face_and_edge(pn_extract_face[i_part], face_edge_idx, face_edge, edge_vtx, &pextract_face_vtx[i_part]);

      pextract_face_vtx_idx[i_part] = malloc((pn_extract_face[i_part]+1) * sizeof(int));
      for(int i = 0; i < pn_extract_face[i_part]+1; ++i) {
        pextract_face_vtx_idx[i_part][i] = face_edge_idx[i];
      }
    }


    PDM_extract_part_vtx_coord_get(extrp,
                                   i_part,
                                   &pextract_vtx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_CELL,
                                  &pextract_cell_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &pextract_face_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_VERTEX,
                                  &pextract_vtx_ln_to_gn[i_part],
                                  PDM_OWNERSHIP_KEEP);

    // log_trace(" %s --> %i \n", name_chr, pn_extract_cell[i_part]);
    /* Vtk en légende */
    char filename[999];
    sprintf(filename, "%s_%3.3d_%3.3d.vtk", name_chr, i_part, i_rank);
    PDM_vtk_write_polydata(filename,
                           pn_extract_vtx[i_part],
                           pextract_vtx[i_part],
                           pextract_vtx_ln_to_gn[i_part],
                           pn_extract_face[i_part],
                           pextract_face_vtx_idx[i_part],
                           pextract_face_vtx[i_part],
                           pextract_face_ln_to_gn[i_part],
                           NULL);

    if (owner == PDM_OWNERSHIP_KEEP) {
      free(pextract_face_vtx    [i_part]);
      free(pextract_face_vtx_idx[i_part]);
    }

  }

  // La visu concatene merge les valeurs donc on voit pas grand choses
  // _visu (name_chr,
  //        extrp->n_part_out,
  //        pn_extract_cell,
  //        pn_extract_face,
  //        pn_extract_vtx,
  //        pextract_cell_face_idx,
  //        pextract_cell_face,
  //        pextract_cell_ln_to_gn,
  //        pextract_face_vtx_idx,
  //        pextract_face_vtx,
  //        pextract_face_ln_to_gn,
  //        pextract_vtx,
  //        pextract_vtx_ln_to_gn);



  free(pn_extract_cell       );
  free(pn_extract_face       );
  free(pn_extract_vtx        );
  free(pextract_cell_face    );
  free(pextract_cell_face_idx);
  free(pextract_face_vtx     );
  free(pextract_face_vtx_idx );
  free(pextract_vtx          );
  free(pextract_cell_ln_to_gn);
  free(pextract_face_ln_to_gn);
  free(pextract_vtx_ln_to_gn );
}

static
void
_compute_extents_3d
(
 int      n_cell,
 int     *cell_face_idx,
 int     *cell_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over cell */
  for(int i_cell = 0; i_cell < n_cell; ++i_cell ) {

    double *_extents = box_extents + 6 * i_cell;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_face = cell_face_idx[i_cell]; idx_face < cell_face_idx[i_cell+1]; ++idx_face) {

      int i_face = PDM_ABS(cell_face[idx_face])-1;

      for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
        int i_vtx = face_vtx[idx_vtx]-1;

        for (int i_dim = 0; i_dim < 3; i_dim++) {
          double x = vtx_coord[3*i_vtx + i_dim];

          if (x < _extents[i_dim]) {
            _extents[i_dim] = x;
          }
          if (x > _extents[3+i_dim]) {
            _extents[3+i_dim] = x;
          }
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}


static
void
_compute_extents_2d_from_face_vtx
(
 int      n_face,
 int     *face_vtx_idx,
 int     *face_vtx,
 double  *vtx_coord,
 double  *box_extents,
 double  *global_extents
)
{
  const double tolerance   = 1.e-12;
  const double eps_extents = 1.e-7;
  const int dim = 3;

  /* Loop over face */
  for(int i_face = 0; i_face < n_face; ++i_face ) {

    double *_extents = box_extents + 6 * i_face;

    /* Init */
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[  i_dim] =  HUGE_VAL;
      _extents[3+i_dim] = -HUGE_VAL;
    }

    /* Loop over face and vtx to compute bbox */
    for(int idx_vtx = face_vtx_idx[i_face]; idx_vtx < face_vtx_idx[i_face+1]; ++idx_vtx) {
      int i_vtx = face_vtx[idx_vtx]-1;

      for (int i_dim = 0; i_dim < 3; i_dim++) {
        double x = vtx_coord[3*i_vtx + i_dim];

        if (x < _extents[i_dim]) {
          _extents[i_dim] = x;
        }
        if (x > _extents[3+i_dim]) {
          _extents[3+i_dim] = x;
        }
      }
    }

    double delta = 0.;
    for (int i_dim = 0; i_dim < 3; i_dim++) {
      double x = _extents[3+i_dim] - _extents[i_dim];

      if (delta < x) {
        delta = x;
      }
    }

    if (delta > eps_extents) {
      delta *= tolerance;
    } else {
      delta = eps_extents;
    }

    for (int i_dim = 0; i_dim < 3; i_dim++) {
      _extents[i_dim]   -= delta;
      _extents[3+i_dim] += delta;
    }

    for (int k = 0; k < dim; k++) {
      global_extents[k]       = PDM_MIN(_extents[k    ], global_extents[k    ]);
      global_extents[dim + k] = PDM_MAX(_extents[dim+k], global_extents[dim+k]);
    }

  } /* End loop cell */
}

static
void
_compute_part_mesh_extents
(
  PDM_part_mesh_t   *mesh,
  int                dim_mesh,
  double            *global_extents,
  double          ***extents_out
)
{
  int n_part = mesh->n_part;
  double **extents = malloc(n_part * sizeof(double *));
  if(dim_mesh == 3) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *cell_face     = NULL;
      int    *cell_face_idx = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_vtx      = NULL;
      double *vtx_coord     = NULL;

      int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_USER);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      extents[i_part] = malloc(6 * n_cell * sizeof(double));
      if(face_vtx == NULL) {
        int    *face_edge_idx  = NULL;
        int    *face_edge      = NULL;
        int    *edge_vtx_idx   = NULL;
        int    *edge_vtx       = NULL;
        int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge , &face_edge_idx, PDM_OWNERSHIP_USER);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx  , &edge_vtx_idx , PDM_OWNERSHIP_USER);
        PDM_compute_face_vtx_from_face_and_edge(n_face, face_edge_idx, face_edge, edge_vtx, &face_vtx);
        _compute_extents_3d(n_cell, cell_face_idx, cell_face, face_edge_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
        free(face_vtx);
      } else {
        _compute_extents_3d(n_cell, cell_face_idx, cell_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
      }
    }
  } else if(dim_mesh == 2) {
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int    *face_vtx      = NULL;
      int    *face_vtx_idx  = NULL;
      int    *face_edge_idx = NULL;
      int    *face_edge     = NULL;
      int    *edge_vtx_idx  = NULL;
      int    *edge_vtx      = NULL;
      double *vtx_coord     = NULL;

      // A gerer le cas mixte face_vtx ou face_edge + edge_vtx

      int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx, PDM_OWNERSHIP_USER);

      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged
      extents[i_part] = malloc(6 * n_face * sizeof(double));

      if(face_vtx != NULL) {
        _compute_extents_2d_from_face_vtx(n_face, face_vtx_idx, face_vtx, vtx_coord, extents[i_part], global_extents);
      } else {
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE , &face_edge, &face_edge_idx, PDM_OWNERSHIP_USER);
        assert(face_edge != NULL);
        PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
        assert(edge_vtx_idx == NULL);
        int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

        edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
        for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
          edge_vtx_idx[i_edge] = 2 * i_edge;
        }
        _compute_extents_3d(n_face, face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);
        // _compute_extents_2d_from_face_edge(n_face, face_edge_idx, face_edge, edge_vtx, vtx_coord, extents[i_part], global_extents);
        free(edge_vtx_idx);
      }
    }

  } else {
    int    *edge_vtx_idx  = NULL;
    int    *edge_vtx      = NULL;
    double *vtx_coord     = NULL;
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);

      extents[i_part] = malloc(6 * n_edge * sizeof(double));
      PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER); // Il faudrait un unchanged

      PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX  , &edge_vtx, &edge_vtx_idx, PDM_OWNERSHIP_USER);
      assert(edge_vtx_idx == NULL);
      edge_vtx_idx = malloc((n_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < n_edge+1; ++i_edge){
        edge_vtx_idx[i_edge] = 2 * i_edge;
      }

      _compute_extents_2d_from_face_vtx(n_edge, edge_vtx_idx, edge_vtx, vtx_coord, extents[i_part], global_extents);

      free(edge_vtx_idx);
    }
  }
  *extents_out = extents;
}

static
void
_select_elements_by_global_bbox
(
  PDM_part_mesh_t *mesh,
  int              dim_mesh,
  double         **box_extents,
  double          *g_mesh_global_extents,
  int            **n_extract_elmt_out,
  double        ***extract_box_extents_out,
  int           ***extract_elmt_init_location_out,
  PDM_g_num_t   ***extract_elmt_ln_to_gn_out
)
{
  int n_part = mesh->n_part;
  int i_rank;
  PDM_MPI_Comm_rank(mesh->comm, &i_rank);

  int          *n_extract_elmt             = malloc(n_part * sizeof(int         *));
  double      **extract_box_extents        = malloc(n_part * sizeof(double      *));
  int         **extract_elmt_init_location = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **extract_elmt_ln_to_gn      = malloc(n_part * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < n_part; ++i_part) {

    int n_entity = 0;
    PDM_g_num_t* entity_ln_to_gn = NULL;
    if(dim_mesh == 3) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    } else if(dim_mesh == 2) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    } else if(dim_mesh == 1) {
      n_entity = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE, &entity_ln_to_gn, PDM_OWNERSHIP_USER);
    }

    // char filename[999];
    // sprintf(filename, "titi_%i.vtk", n_entity);
    // PDM_vtk_write_boxes(filename,
    //                     n_entity,
    //                     box_extents[i_part],
    //                     NULL);

    n_extract_elmt[i_part] = 0;
    extract_box_extents       [i_part] = malloc(6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = malloc(3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = malloc(    n_entity * sizeof(PDM_g_num_t));

    double *_box_extents = box_extents[i_part];

    for(int i = 0; i < n_entity; ++i) {

      double *box_min = _box_extents + 6*i;
      double *box_max = box_min + 3;

      int intersect = 1;
      for (int j = 0; j < 3; j++) {
        if (box_min[j] > g_mesh_global_extents[j+3] ||
            box_max[j] < g_mesh_global_extents[j  ]) {
          intersect = 0;
          break;
        }
      }

      if (intersect) {
        for (int j = 0; j < 6; j++) {
          extract_box_extents  [i_part][6*n_extract_elmt[i_part]+j] = _box_extents[6*i+j];
        }
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]  ] = i_rank;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+1] = i_part;
        extract_elmt_init_location[i_part][3*n_extract_elmt[i_part]+2] = i;

        extract_elmt_ln_to_gn[i_part][n_extract_elmt[i_part]] = entity_ln_to_gn[i];

        n_extract_elmt[i_part]++;
      }
    }
    extract_box_extents       [i_part] = realloc(extract_box_extents       [i_part], 6 * n_entity * sizeof(double     ));
    extract_elmt_init_location[i_part] = realloc(extract_elmt_init_location[i_part], 3 * n_entity * sizeof(int        ));
    extract_elmt_ln_to_gn     [i_part] = realloc(extract_elmt_ln_to_gn     [i_part],     n_entity * sizeof(PDM_g_num_t));

  }

  *n_extract_elmt_out             = n_extract_elmt;
  *extract_box_extents_out        = extract_box_extents;
  *extract_elmt_init_location_out = extract_elmt_init_location;
  *extract_elmt_ln_to_gn_out      = extract_elmt_ln_to_gn;
}


static
void
_redistrib_boxes
(
 PDM_MPI_Comm    comm,
 PDM_box_set_t  *boxes_mesh_a,
 PDM_box_set_t  *boxes_mesh_b,
 int            *box_a_to_box_b_idx,
 int            *box_a_to_box_b,
 int           **redistribute_box_a_to_box_b_idx,
 int           **redistribute_box_a_to_box_b
)
{

  int n_elt_mesh_a = PDM_box_set_get_size (boxes_mesh_a);
  int n_elt_mesh_b = PDM_box_set_get_size (boxes_mesh_b);

  PDM_g_num_t *gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_a);
  PDM_g_num_t *gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_mesh_b);

  /*****************************************************************************
   *                                                                           *
   *  Transfer intersection information from partitions to blocks              *
   * with PDM_part_to_block_exch function                                      *
   *                                                                           *
   *  Results :                                                                *
   *      - block_a_boxes_b_idx                                                *
   *      - block_a_boxes_b_gnum_data                                          *
   *                                                                           *
   ****************************************************************************/

  /*
   * Tentative Bruno :
   *   - Ponderate work
   *   - Extract only cell with job
   */
  double* weight = (double *) malloc( n_elt_mesh_a * sizeof(double));
  // PDM_g_num_t* extract_mesh_a_g_num = (PDM_g_num_t *) malloc( n_elt_mesh_a * sizeof(PDM_g_num_t));
  for (int i = 0; i < n_elt_mesh_a; i++) {
    weight[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  // TODO : Geometric to better locality
  PDM_part_to_block_t *ptb_boxes_a = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                              PDM_PART_TO_BLOCK_POST_MERGE,
                                                              1.,
                                            (PDM_g_num_t **) &gnum_elt_mesh_a,
                                                             &weight,
                                                              &n_elt_mesh_a,
                                                              1,
                                                              comm);

  int n_elt_block_a = PDM_part_to_block_n_elt_block_get (ptb_boxes_a);
  free(weight);

  PDM_g_num_t *block_gnum_a = PDM_part_to_block_block_gnum_get (ptb_boxes_a);

  int *part_stride_a = (int *) malloc (sizeof(int) * n_elt_mesh_a);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    part_stride_a[i] = box_a_to_box_b_idx[i+1] - box_a_to_box_b_idx[i];
  }

  /*
   * Exchange connectivity box_a_to_box_b
   */
  PDM_g_num_t *box_a_to_box_b_g_num = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) *
                                                              box_a_to_box_b_idx[n_elt_mesh_a]);

  for (int k = 0; k < box_a_to_box_b_idx[n_elt_mesh_a]; k++) {
    box_a_to_box_b_g_num[k] = gnum_elt_mesh_b[box_a_to_box_b[k]];
  }

  int         *block_a_boxes_b_stride;
  PDM_g_num_t *block_a_boxes_b_gnum_data;

  PDM_part_to_block_exch (ptb_boxes_a,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         0,
                         &part_stride_a,
               (void **) &box_a_to_box_b_g_num,
                         &block_a_boxes_b_stride,
               (void **) &block_a_boxes_b_gnum_data);
  free(box_a_to_box_b_g_num);


  // if (1) {
  //   int idx = 0;
  //   log_trace("--- Avant ---\n");
  //   for (int i = 0; i < n_elt_block_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", block_gnum_a[i]);
  //     for (int j = 0; j < block_a_boxes_b_stride[i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, block_a_boxes_b_gnum_data[idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }

  /* Remove duplicates */
  if (1) {
    int idx_read  = 0;
    int idx_write = 0;
    for (int i = 0; i < n_elt_block_a; i++) {

      int n = block_a_boxes_b_stride[i];

      int m = PDM_inplace_unique_long(block_a_boxes_b_gnum_data + idx_read,
                                      NULL,
                                      0,
                                      n-1);

      for (int j = 0; j < m; j++) {
        block_a_boxes_b_gnum_data[idx_write++] = block_a_boxes_b_gnum_data[idx_read+j];
      }
      block_a_boxes_b_stride[i] = m;

      idx_read += n;
    }
  }

  // if (1) {
  //   int idx = 0;
  //   log_trace("--- Après ---\n");
  //   for (int i = 0; i < n_elt_block_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", block_gnum_a[i]);
  //     for (int j = 0; j < block_a_boxes_b_stride[i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, block_a_boxes_b_gnum_data[idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }




  /*****************************************************************************
   *                                                                           *
   * Redistribute boxes_mesh_a intersections to ensure a good load balacing          *
   * in comm MPI communicator                                              *
   * This step removes intersections found many times on different ranks       *
   *                                                                           *
   * After this step, data are stored in a block with n_elt_block_a, block_gnum_a,
   * part_stride_a                                                              *
   *                                                                           *
   ****************************************************************************/
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  /*
   * Redistribute data boxes_mesh_a from blockB distribution with a PDM_box_distrib_t
   * structure
   * TODO: - Build a new PDM_box_distrib_t structure more simple
   *       - Hide PDM_box_distrib_t attributes
   */

  PDM_l_num_t *destination = PDM_part_to_block_destination_get (ptb_boxes_a);

  PDM_g_num_t n_g_elmt_mesh_a = PDM_box_set_get_global_size(boxes_mesh_a);

  PDM_box_distrib_t *distrib_a = PDM_box_distrib_create(n_elt_mesh_a,
                                                        n_g_elmt_mesh_a,
                                                        1, // Don't use in this case
                                                        comm);

  PDM_g_num_t n_g_elmt_mesh_b = PDM_box_set_get_global_size(boxes_mesh_b);

  PDM_box_distrib_t *distrib_b = PDM_box_distrib_create(n_elt_mesh_b,
                                                        n_g_elmt_mesh_b,
                                                        1, // Don't use in this case
                                                        comm);


  int *count_elts_a = (int *) malloc (sizeof(int) * n_rank);
  int *count_elts_b = (int *) malloc (sizeof(int) * n_rank);

  for (int i = 0; i < n_rank + 1; i++) {
    distrib_a->index[i] = 0;
    distrib_b->index[i] = 0;
  }

  for (int i = 0; i < n_rank; i++) {
    count_elts_a[i] = 0;
    count_elts_b[i] = 0;
  }

  for (int i = 0; i < n_elt_mesh_a; i++) {
    int t_rank = destination[i] + 1;

    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      distrib_a->index[t_rank]++;
      distrib_b->index[t_rank] += part_stride_a[i];
    }
  }

  for (int i = 0; i < n_rank; i++) {
    distrib_a->index[i+1] += distrib_a->index[i];
    distrib_b->index[i+1] += distrib_b->index[i];
  }

  distrib_a->list = (int *) malloc (sizeof(int) * distrib_a->index[n_rank]);
  distrib_b->list = (int *) malloc (sizeof(int) * distrib_b->index[n_rank]);

  for (int i = 0; i < n_elt_mesh_a; i++) {
    if(part_stride_a[i] > 0 ) { // To see with Eric and Bastien --> I use it to extract only the intersect part
      int t_rank = destination[i]; // EQU + 1; mais ce n est pas necessaire
      int idx_a = distrib_a->index[t_rank] + (count_elts_a[t_rank]++);
      distrib_a->list[idx_a] = i;
      int idx_b = distrib_b->index[t_rank] + count_elts_b[t_rank];
      count_elts_b[t_rank] += part_stride_a[i];
      int k=0;
      for (int j = box_a_to_box_b_idx[i]; j < box_a_to_box_b_idx[i+1]; j++) {
        distrib_b->list[idx_b+k++] = box_a_to_box_b[j];
      }
    }
  }


  free (part_stride_a);
  free (count_elts_a);
  free (count_elts_b);

  PDM_box_distrib_clean (distrib_a);
  PDM_box_distrib_clean (distrib_b);

  PDM_box_set_redistribute (distrib_a, boxes_mesh_a);
  PDM_box_set_redistribute (distrib_b, boxes_mesh_b);

  PDM_box_distrib_destroy (&distrib_a);
  PDM_box_distrib_destroy (&distrib_b);

  PDM_box_set_remove_duplicate (boxes_mesh_a);
  PDM_box_set_remove_duplicate (boxes_mesh_b);

  /*
   * All boxes are redistribute we need to update box_a_to_box_b array
   *    - Caution if morton / hilbert the array block_gnum_a IS NOT order
   */
  n_elt_mesh_a    = PDM_box_set_get_size (boxes_mesh_a); // Caution not the same of the fist call because redistibute
  n_elt_mesh_b    = PDM_box_set_get_size (boxes_mesh_b); // Caution not the same of the fist call because redistibute

  gnum_elt_mesh_a = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_a);
  gnum_elt_mesh_b = (PDM_g_num_t *) PDM_box_set_get_g_num(boxes_mesh_b);

  PDM_block_to_part_t *btp = PDM_block_to_part_create_from_sparse_block(block_gnum_a,
                                                                        n_elt_block_a,
                                                (const PDM_g_num_t **)  &gnum_elt_mesh_a,
                                                                        &n_elt_mesh_a,
                                                                        1,
                                                                        comm);

  PDM_part_to_block_free (ptb_boxes_a);
  int         **tmp_redistribute_box_a_to_box_b_n     = NULL;
  PDM_g_num_t **tmp_redistribute_box_a_to_box_b_g_num = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         block_a_boxes_b_stride,
                         block_a_boxes_b_gnum_data,
                         &tmp_redistribute_box_a_to_box_b_n,
           (void ***)    &tmp_redistribute_box_a_to_box_b_g_num);
  free (block_a_boxes_b_stride);
  free (block_a_boxes_b_gnum_data);

  // if (1) {
  //   int idx = 0;
  //   log_trace("--- 2 ---\n");
  //   for (int i = 0; i < n_elt_mesh_a; i++) {
  //     log_trace("faceA "PDM_FMT_G_NUM": faceB ", gnum_elt_mesh_a[i]);
  //     for (int j = 0; j < tmp_redistribute_box_a_to_box_b_n[0][i]; j++) {
  //       log_trace(" "PDM_FMT_G_NUM, tmp_redistribute_box_a_to_box_b_g_num[0][idx++]);
  //     }
  //     log_trace("\n");
  //   }
  // }

  PDM_block_to_part_free(btp);

  int         *redistribute_box_a_to_box_b_n     = tmp_redistribute_box_a_to_box_b_n    [0];
  PDM_g_num_t *redistribute_box_a_to_box_b_g_num = tmp_redistribute_box_a_to_box_b_g_num[0];
  free(tmp_redistribute_box_a_to_box_b_n    );
  free(tmp_redistribute_box_a_to_box_b_g_num);


  /*
   * Translate in frame of B
   */

  int         *order              = (int         *) malloc(n_elt_mesh_b * sizeof(int        ));
  PDM_g_num_t *gnum_elt_mesh_b_cp = (PDM_g_num_t *) malloc(n_elt_mesh_b * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_elt_mesh_b; ++i ) {
    order             [i] = i;
    gnum_elt_mesh_b_cp[i] = gnum_elt_mesh_b[i];
  }

  PDM_sort_long(gnum_elt_mesh_b_cp, order, n_elt_mesh_b);


  int *_redistribute_box_a_to_box_b_idx = (int *) malloc((n_elt_mesh_a+1) * sizeof(int));
  _redistribute_box_a_to_box_b_idx[0] = 0;
  int n_tot_connect = 0;
  for(int i = 0; i < n_elt_mesh_a; ++i) {
    n_tot_connect += redistribute_box_a_to_box_b_n[i];
    _redistribute_box_a_to_box_b_idx[i+1] = _redistribute_box_a_to_box_b_idx[i] + redistribute_box_a_to_box_b_n[i];
  }

  int *_redistribute_box_a_to_box_b = (int *) malloc( n_tot_connect * sizeof(int));

  for(int i = 0; i < n_tot_connect; ++i) {
    int pos = PDM_binary_search_long(redistribute_box_a_to_box_b_g_num[i], gnum_elt_mesh_b_cp, n_elt_mesh_b);
    _redistribute_box_a_to_box_b[i] = order[pos];
  }

  *redistribute_box_a_to_box_b_idx = _redistribute_box_a_to_box_b_idx;
  *redistribute_box_a_to_box_b     = _redistribute_box_a_to_box_b;

  free(order);
  free(gnum_elt_mesh_b_cp);
  free(redistribute_box_a_to_box_b_n    );
  free(redistribute_box_a_to_box_b_g_num);
}


static
PDM_extract_part_t*
_create_extract_part
(
 PDM_part_mesh_t *mesh,
 int              dim_mesh,
 PDM_box_set_t   *boxes_meshes
)
{
  int n_part_out = 1;
  PDM_extract_part_t* extrp_mesh = PDM_extract_part_create(dim_mesh,
                                                           mesh->n_part,
                                                           n_part_out,
                                                           PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT, // Not used
                                                           PDM_FALSE,                   // compute_child_gnum
                                                           PDM_OWNERSHIP_KEEP,
                                                           mesh->comm);

  int n_elt_mesh = PDM_box_set_get_size (boxes_meshes);

  // printf("n_elt_mesh = %i  \n", n_elt_mesh);

  PDM_g_num_t *gnum_elt_mesh = (PDM_g_num_t *) PDM_box_set_get_g_num (boxes_meshes);

  int *init_location_elt_mesh = (int  *) PDM_box_set_origin_get(boxes_meshes);


  for(int i_part = 0; i_part < mesh->n_part; ++i_part) {

    int n_cell = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_CELL  );
    int n_face = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_FACE  );
    int n_edge = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_EDGE  );
    int n_vtx  = PDM_part_mesh_n_entity_get(mesh, i_part, PDM_MESH_ENTITY_VERTEX);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    PDM_g_num_t *face_ln_to_gn = NULL;
    PDM_g_num_t *edge_ln_to_gn = NULL;
    PDM_g_num_t *vtx_ln_to_gn  = NULL;

    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_CELL  , &cell_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_FACE  , &face_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_EDGE  , &edge_ln_to_gn, PDM_OWNERSHIP_USER);
    PDM_part_mesh_entity_ln_to_gn_get(mesh, i_part, PDM_MESH_ENTITY_VERTEX, &vtx_ln_to_gn , PDM_OWNERSHIP_USER);

    int *cell_face     = NULL;
    int *cell_face_idx = NULL;
    int *face_vtx      = NULL;
    int *face_vtx_idx  = NULL;
    int *face_edge     = NULL;
    int *face_edge_idx = NULL;
    int *edge_vtx      = NULL;
    int *edge_vtx_idx  = NULL;

    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, &cell_face, &cell_face_idx, PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , &face_vtx , &face_vtx_idx , PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, &face_edge, &face_edge_idx, PDM_OWNERSHIP_USER);
    PDM_part_mesh_connectivity_get(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , &edge_vtx , &edge_vtx_idx , PDM_OWNERSHIP_USER);

    double *vtx_coord = NULL;
    PDM_part_mesh_vtx_coord_get(mesh, i_part, &vtx_coord, PDM_OWNERSHIP_USER);

    PDM_extract_part_part_set(extrp_mesh,
                              i_part,
                              n_cell,
                              n_face,
                              n_edge,
                              n_vtx ,
                              cell_face_idx,
                              cell_face,
                              face_edge_idx,
                              face_edge,
                              edge_vtx,
                              face_vtx_idx,
                              face_vtx,
                              cell_ln_to_gn,
                              face_ln_to_gn,
                              edge_ln_to_gn,
                              vtx_ln_to_gn,
                              vtx_coord);
  }




  /*  Setup target frame */
  PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, gnum_elt_mesh, init_location_elt_mesh);
  // PDM_g_num_t *target_g_num = malloc(sizeof(PDM_g_num_t) * n_elt_mesh);
  // memcpy(target_g_num, gnum_elt_mesh, sizeof(PDM_g_num_t) * n_elt_mesh);
  // PDM_extract_part_target_set(extrp_mesh, 0, n_elt_mesh, target_g_num, init_location_elt_mesh);

  PDM_extract_part_compute(extrp_mesh);

  return extrp_mesh;
}


 static PDM_ownership_t
_get_extracted_mesh_vol
(
 PDM_extract_part_t  *extrp,
 int                 *n_cell,
 int                 *n_face,
 int                 *n_vtx,
 int                **cell_face_idx,
 int                **cell_face,
 int                **face_vtx_idx,
 int                **face_vtx,
 double             **vtx_coord,
 PDM_g_num_t        **cell_ln_to_gn,
 PDM_g_num_t        **vtx_ln_to_gn
 )
 {
  *n_cell = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                              cell_face,
                                              cell_face_idx,
                                              PDM_OWNERSHIP_KEEP);

  *n_face = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                              face_vtx,
                                              face_vtx_idx,
                                              PDM_OWNERSHIP_KEEP);

  PDM_ownership_t owner = PDM_OWNERSHIP_KEEP;
  if (*face_vtx == NULL) {
    owner = PDM_OWNERSHIP_USER;

    int *face_edge = NULL;
    PDM_extract_part_connectivity_get(extrp,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                      &face_edge,
                                      face_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

    int *edge_vtx     = NULL;
    int *edge_vtx_idx = NULL;
    PDM_extract_part_connectivity_get(extrp,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

    PDM_compute_face_vtx_from_face_and_edge(*n_face,
                                            *face_vtx_idx,
                                            face_edge,
                                            edge_vtx,
                                            face_vtx);
  }

  *n_vtx = PDM_extract_part_vtx_coord_get(extrp,
                                          0,
                                          vtx_coord,
                                          PDM_OWNERSHIP_KEEP);

  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       0,
                                       PDM_MESH_ENTITY_CELL,
                                       cell_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);

  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       0,
                                       PDM_MESH_ENTITY_VERTEX,
                                       vtx_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);

  return owner;
 }


static void
_export_ensight3d
(
 PDM_MPI_Comm  comm,
 const char   *name,
 int           n_cell,
 int           n_face,
 int           n_vtx,
 int          *cell_face_idx,
 int          *cell_face,
 int          *face_vtx_idx,
 int          *face_vtx,
 double       *vtx_coord,
 PDM_g_num_t  *cell_ln_to_gn,
 PDM_g_num_t  *vtx_ln_to_gn
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);


  PDM_gen_gnum_t *gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 1., comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gnum_vtx,
                            0,
                            n_vtx,
                            vtx_ln_to_gn);

  PDM_gnum_compute(gnum_vtx);

  PDM_g_num_t *extract_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx,
                                                   0);


  PDM_gen_gnum_t *gnum_cell = PDM_gnum_create(3, 1, PDM_FALSE, 1., comm, PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_from_parents(gnum_cell,
                            0,
                            n_cell,
                            cell_ln_to_gn);

  PDM_gnum_compute(gnum_cell);

  PDM_g_num_t *extract_cell_ln_to_gn = PDM_gnum_get(gnum_cell,
                                                    0);


  PDM_writer_t *wrt = PDM_writer_create("Ensight",
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_CST,
                                        PDM_WRITER_OFF,
                                        name,
                                        name,
                                        PDM_MPI_COMM_WORLD,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       "iso_surface_volume_mesh",
                                       1);

  int id_var_rank = PDM_writer_var_create(wrt,
                                          PDM_WRITER_ON,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "i_rank");

  int id_var_gnum = PDM_writer_var_create(wrt,
                                          PDM_WRITER_ON,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "cell_g_num");

  PDM_writer_step_beg(wrt, 0.);

  int *cell_face_n = malloc(sizeof(int) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    cell_face_n[i] = cell_face_idx[i+1] - cell_face_idx[i];
  }

  int *face_vtx_n  = malloc(sizeof(int) * n_face);
  for (int i = 0; i < n_face; i++) {
    face_vtx_n[i] = face_vtx_idx[i+1] - face_vtx_idx[i];
  }


  PDM_writer_geom_coord_set(wrt,
                            id_geom,
                            0,
                            n_vtx,
                            vtx_coord,
                            extract_vtx_ln_to_gn,
                            PDM_OWNERSHIP_USER);

  PDM_writer_geom_cell3d_cellface_add(wrt,
                                      id_geom,
                                      0,
                                      n_cell,
                                      n_face,
                                      face_vtx_idx,
                                      face_vtx_n,
                                      face_vtx,
                                      cell_face_idx,
                                      cell_face_n,
                                      cell_face,
                                      extract_cell_ln_to_gn);

  PDM_writer_geom_write(wrt,
                        id_geom);

  PDM_real_t *val_rank = malloc(sizeof(PDM_real_t) * n_cell);
  PDM_real_t *val_gnum = malloc(sizeof(PDM_real_t) * n_cell);
  for (int i = 0; i < n_cell; i++) {
    val_rank[i] = i_rank;
    val_gnum[i] = cell_ln_to_gn[i];
  }

  PDM_writer_var_set(wrt,
                     id_var_rank,
                     id_geom,
                     0,
                     val_rank);

  PDM_writer_var_write(wrt,
                       id_var_rank);
  PDM_writer_var_free(wrt,
                      id_var_rank);


  PDM_writer_var_set(wrt,
                     id_var_rank,
                     id_geom,
                     0,
                     val_gnum);

  PDM_writer_var_write(wrt,
                       id_var_gnum);
  PDM_writer_var_free(wrt,
                      id_var_gnum);

  PDM_writer_step_end(wrt);

  PDM_writer_free(wrt);

  free(val_rank   );
  free(val_gnum   );
  free(cell_face_n);
  free(face_vtx_n );

  PDM_gnum_free(gnum_vtx );
  PDM_gnum_free(gnum_cell);
}

static void
_dump_elementary_vol_vol
(
 PDM_MPI_Comm  comm,
 int           cellA_id,
 int           faceA_id,
 int           cellB_id,
 int           triaA_id,
 int           faceB_id,
 int           triaB_id,
 double       *tetraA_coord,
 double       *triaB_coord
 )
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  char filename[999];
  sprintf(filename, "vol_vol_A%d_%d_%d_B%d_%d_%d_rank_%d.vtk",
          cellA_id,
          faceA_id,
          triaA_id,
          cellB_id,
          faceB_id,
          triaB_id,
          i_rank);


  FILE *f = fopen(filename, "w");

  fprintf(f, "# vtk DataFile Version 2.0\n");
  fprintf(f, "mesh\n");
  fprintf(f, "ASCII\n");
  fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

  fprintf(f, "POINTS %d double\n", 4 + 3);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", tetraA_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fprintf(f, "%.20lf ", triaB_coord[3*i+j]);
    }
    fprintf(f, "\n");
  }

  // fprintf(f, "CELLS %d %d\n", 2, 5 + 4);
  // fprintf(f, "4 0 1 2 3\n3 4 5 6\n");

  // fprintf(f, "CELL_TYPES 2\n10\n5\n");

  // fprintf(f, "CELL_DATA 2\n");
  // fprintf(f, "FIELD field 4\n");
  // fprintf(f, "mesh_id 1 2 int\n0 1\n");
  // fprintf(f, "cell_id 1 2 int\n%d %d\n", cellA_id, cellB_id);
  // fprintf(f, "face_id 1 2 int\n%d %d\n", faceA_id, faceB_id);
  // fprintf(f, "tria_id 1 2 int\n%d %d\n", triaA_id, triaB_id);
  fprintf(f, "CELLS %d %d\n", 1 + 3 + 2 + 1, 2 + 3*3 + 2*4 + 4);
  fprintf(f, "1 0\n2 0 1\n2 0 2\n2 0 3\n3 1 2 3\n3 0 1 2\n3 4 5 6\n");

  fprintf(f, "CELL_TYPES 7\n1\n3\n3\n3\n5\n5\n5\n");

  fprintf(f, "CELL_DATA 7\n");
  fprintf(f, "FIELD field 4\n");
  fprintf(f, "mesh_id 1 7 int\n0 0 0 0 0 0 1\n");
  fprintf(f, "cell_id 1 7 int\n%d %d %d %d %d %d %d\n", cellA_id, cellA_id, cellA_id, cellA_id, cellA_id, cellA_id, cellB_id);
  fprintf(f, "face_id 1 7 int\n%d %d %d %d %d %d %d\n", faceA_id, faceA_id, faceA_id, faceA_id, faceA_id, faceA_id, faceB_id);
  fprintf(f, "tria_id 1 7 int\n%d %d %d %d %d %d %d\n", triaA_id, triaA_id, triaA_id, triaA_id, triaA_id, triaA_id, triaB_id);

  fclose(f);
}
static
void
_mesh_intersection_vol_vol
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  /* method : 0 -> with pointers, 1 -> without */
  int method = 1;

  int dbg = 1;
  int dbg_Karmijn = 0;

  int *cellA_cellB_idx = redistribute_box_a_to_box_b_idx;
  int *cellA_cellB     = redistribute_box_a_to_box_b;

  /* Get connectivities and coordinates */
  int          n_cellA = 0;
  int          n_faceA = 0;
  int          n_vtxA  = 0;
  int         *cellA_faceA_idx = NULL;
  int         *cellA_faceA     = NULL;
  int         *faceA_vtxA_idx  = NULL;
  int         *faceA_vtxA      = NULL;
  double      *vtxA_coord      = NULL;
  PDM_g_num_t *cellA_ln_to_gn  = NULL;
  PDM_g_num_t *vtxA_ln_to_gn   = NULL;
  PDM_ownership_t owner_face_vtxA = _get_extracted_mesh_vol(extrp_mesh_a,
                                                            &n_cellA,
                                                            &n_faceA,
                                                            &n_vtxA,
                                                            &cellA_faceA_idx,
                                                            &cellA_faceA,
                                                            &faceA_vtxA_idx,
                                                            &faceA_vtxA,
                                                            &vtxA_coord,
                                                            &cellA_ln_to_gn,
                                                            &vtxA_ln_to_gn);

  int          n_cellB = 0;
  int          n_faceB = 0;
  int          n_vtxB  = 0;
  int         *cellB_faceB_idx = NULL;
  int         *cellB_faceB     = NULL;
  int         *faceB_vtxB_idx  = NULL;
  int         *faceB_vtxB      = NULL;
  double      *vtxB_coord      = NULL;
  PDM_g_num_t *cellB_ln_to_gn  = NULL;
  PDM_g_num_t *vtxB_ln_to_gn   = NULL;
  PDM_ownership_t owner_face_vtxB = _get_extracted_mesh_vol(extrp_mesh_b,
                                                            &n_cellB,
                                                            &n_faceB,
                                                            &n_vtxB,
                                                            &cellB_faceB_idx,
                                                            &cellB_faceB,
                                                            &faceB_vtxB_idx,
                                                            &faceB_vtxB,
                                                            &vtxB_coord,
                                                            &cellB_ln_to_gn,
                                                            &vtxB_ln_to_gn);

  if (dbg && 0) {
    PDM_log_trace_connectivity_int(cellA_cellB_idx,
                                   cellA_cellB,
                                   n_cellA,
                                   "cellA_cellB : ");
  }



  /*
   * Panic vtk
   */
  if (dbg) {
    _export_vtk_3d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_3d("extrp_mesh_b", extrp_mesh_b);

    _export_ensight3d(mi->comm,
                      "vol_vol_meshA",
                      n_cellA,
                      n_faceA,
                      n_vtxA,
                      cellA_faceA_idx,
                      cellA_faceA,
                      faceA_vtxA_idx,
                      faceA_vtxA,
                      vtxA_coord,
                      cellA_ln_to_gn,
                      vtxA_ln_to_gn);

    _export_ensight3d(mi->comm,
                      "vol_vol_meshB",
                      n_cellB,
                      n_faceB,
                      n_vtxB,
                      cellB_faceB_idx,
                      cellB_faceB,
                      faceB_vtxB_idx,
                      faceB_vtxB,
                      vtxB_coord,
                      cellB_ln_to_gn,
                      vtxB_ln_to_gn);
  }



  /*
   *  Triangulate all faces of mesh B
   */
  int max_face_vtx_n = 0;
  int s_triaA_vtxA = 0;
  for (int i = 0; i < n_faceA; i++) {
    int face_vtx_n = faceA_vtxA_idx[i+1] - faceA_vtxA_idx[i];
    max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    s_triaA_vtxA   = PDM_MAX(s_triaA_vtxA,   face_vtx_n);
  }
  s_triaA_vtxA = 3*(s_triaA_vtxA - 2);

  int *faceB_triaB_idx = malloc(sizeof(int) * (n_faceB + 1));
  faceB_triaB_idx[0] = 0;
  for (int i = 0; i < n_faceB; i++) {
    int face_vtx_n = faceB_vtxB_idx[i+1] - faceB_vtxB_idx[i];
    max_face_vtx_n = PDM_MAX(max_face_vtx_n, face_vtx_n);
    int face_tria_n = face_vtx_n - 2;
    faceB_triaB_idx[i+1] = faceB_triaB_idx[i] + face_tria_n;
  }

  PDM_triangulate_state_t *tri_state = PDM_triangulate_state_create(max_face_vtx_n);

  int *triaB_vtxB = malloc(sizeof(int) * faceB_triaB_idx[n_faceB] * 3);
  int *triaA_vtxA = malloc(sizeof(int) * s_triaA_vtxA);

  for (int faceB_id = 0; faceB_id < n_faceB; faceB_id++) {

    int *_face_vtx = faceB_vtxB + faceB_vtxB_idx[faceB_id];
    int face_vtx_n = faceB_vtxB_idx[faceB_id+1] - faceB_vtxB_idx[faceB_id];

    int *_tria_vtx  = triaB_vtxB + 3*faceB_triaB_idx[faceB_id];

    int n_tria;
    if (face_vtx_n == 3) {
      /* Triangular face */
      n_tria = 1;
      memcpy(_tria_vtx, _face_vtx, sizeof(int) * 3);
    }
    else if (face_vtx_n == 4) {
      /* Quadrilateral face */
      n_tria = PDM_triangulate_quadrangle(3,
                                          vtxB_coord,
                                          NULL,
                                          _face_vtx,
                                          _tria_vtx);
    }
    else {
      /* Polygonal face */
      n_tria = PDM_triangulate_polygon(3,
                                       face_vtx_n,
                                       vtxB_coord,
                                       NULL,
                                       _face_vtx,
                                       PDM_TRIANGULATE_MESH_DEF,
                                       _tria_vtx,
                                       tri_state);
    }

    assert(n_tria == faceB_triaB_idx[faceB_id+1] - faceB_triaB_idx[faceB_id]);

  } // End of loop on faces B



  if (dbg) {

    int *triaB_faceB_id = malloc(sizeof(int) * faceB_triaB_idx[n_faceB]);
    for (int faceB_id = 0; faceB_id < n_faceB; faceB_id++) {
      for (int triaB_id = faceB_triaB_idx[faceB_id]; triaB_id < faceB_triaB_idx[faceB_id+1]; triaB_id++) {
        triaB_faceB_id[triaB_id] = faceB_id;
      }
    }

    int i_rank;
    PDM_MPI_Comm_rank(mi->comm, &i_rank);
    char filename[999];
    sprintf(filename, "triaB_%d.vtk", i_rank);

    const char *field_name[]   = {"faceB_id"};
    const int  *field_value[1] = {triaB_faceB_id};

    PDM_vtk_write_std_elements(filename,
                               n_vtxB,
                               vtxB_coord,
                               vtxB_ln_to_gn,
                               PDM_MESH_NODAL_TRIA3,
                               faceB_triaB_idx[n_faceB],
                               triaB_vtxB,
                               NULL,
                               1,
                               field_name,
                               field_value);
    free(triaB_faceB_id);
  }



  double *cellA_volume = malloc(sizeof(double) * n_cellA);
  double *cellA_center = malloc(sizeof(double) * n_cellA * 3);
  PDM_geom_elem_polyhedra_properties(1,
                                     n_cellA,
                                     n_faceA,
                                     faceA_vtxA_idx,
                                     faceA_vtxA,
                                     cellA_faceA_idx,
                                     cellA_faceA,
                                     n_vtxA,
                                     vtxA_coord,
                                     cellA_volume,
                                     cellA_center,
                                     NULL,
                                     NULL);



  double *cellA_cellB_volume = malloc(sizeof(double) * cellA_cellB_idx[n_cellA]);

  //-->>
  for (int i = 0; i < cellA_cellB_idx[n_cellA]; i++) {
    cellA_cellB_volume[i] = 0;
  }
  //<<--

  double tetraA_coord[12];
  double triaB_coord[9];

  // malloc for atomic computation
  // int max_size = 10*2 + 3;
  // Element **cll_storage = malloc(sizeof(Element *) * max_size);
  // for (int i = 0; i < max_size; i++) {
  //   cll_storage[i] = malloc(sizeof(Element));
  // }

  // for vtk
  // TO DO: add reallocs
  int min_sizeA = 2000;
  int min_sizeB = 2000;

  double *vtk_vtx_coordA    = NULL;
  int     vtk_n_vtxA        = 0;
  int    *vtk_face_vtx_idxA = NULL;
  int    *vtk_face_vtxA     = NULL;
  int     vtk_n_faceA       = 0;
  double *vtk_vtx_coordB    = NULL;
  int     vtk_n_vtxB        = 0;
  int    *vtk_face_vtx_idxB = NULL;
  int    *vtk_face_vtxB     = NULL;
  int     vtk_n_faceB       = 0;
  if (dbg_Karmijn) {
    vtk_vtx_coordA    = malloc(sizeof(double) * min_sizeA * 3);
    vtk_face_vtx_idxA = malloc(sizeof(int) * (min_sizeA + 1));
    vtk_face_vtx_idxA[0] = 0;
    vtk_face_vtxA     = malloc(sizeof(int) * min_sizeA);
    vtk_vtx_coordB    = malloc(sizeof(double) * min_sizeB * 3);
    vtk_face_vtx_idxB = malloc(sizeof(int) * (min_sizeB + 1));
    vtk_face_vtx_idxB[0] = 0;
    vtk_face_vtxB     = malloc(sizeof(int) * min_sizeB);
  }

  for (int cellA_id = 0; cellA_id < n_cellA; cellA_id++) {

    if (cellA_id != 2) continue;

    /* Compute a 'center' point to tetrahedrize current cell A */
    // Reference vertex = 1st vertex of first face
    // Use cell center instead????
    int ref_vtxA_id = -1;
    int first_face_id = PDM_ABS(cellA_faceA[cellA_faceA_idx[cellA_id]]) - 1;
    ref_vtxA_id = faceA_vtxA[faceA_vtxA_idx[first_face_id]];


    if (dbg) {
      log_trace("cellA_id %d ("PDM_FMT_G_NUM") : ref_vtxA_id = %d\n",
                cellA_id, cellA_ln_to_gn[cellA_id], ref_vtxA_id);
    }

    memcpy(&tetraA_coord[0], &vtxA_coord[3*(ref_vtxA_id-1)], sizeof(double)*3);
    // memcpy(&tetraA_coord[0], &cellA_center[3*cellA_id], sizeof(double)*3);
    // tetraA_coord[0] = 1.234569;
    // tetraA_coord[1] = -10.3147;
    // tetraA_coord[2] = 3.4242424213;

    /* Initialize intersection volumes to zero */
    for (int icellB = cellA_cellB_idx[cellA_id]; icellB < cellA_cellB_idx[cellA_id+1]; icellB++) {
      cellA_cellB_volume[icellB] = 0.;
    }

    /* Loop on faces A */
    for (int ifaceA = cellA_faceA_idx[cellA_id]; ifaceA < cellA_faceA_idx[cellA_id+1]; ifaceA++) {

      int faceA_id   = PDM_ABS (cellA_faceA[ifaceA]) - 1;
      int faceA_sign = PDM_SIGN(cellA_faceA[ifaceA]);

      if (dbg) {
        log_trace("  faceA_id %d, sign = %d\n", faceA_id, faceA_sign);
      }

      /* Triangulate current face A */
      int *_face_vtx = faceA_vtxA + faceA_vtxA_idx[faceA_id];
      int face_vtx_n = faceA_vtxA_idx[faceA_id+1] - faceA_vtxA_idx[faceA_id];

      int n_tria;
      if (face_vtx_n == 3) {
        /* Triangular face */
        n_tria = 1;
        memcpy(triaA_vtxA, _face_vtx, sizeof(int) * 3);
      }
      else if (face_vtx_n == 4) {
        /* Quadrilateral face */
        n_tria = PDM_triangulate_quadrangle(3,
                                            vtxA_coord,
                                            NULL,
                                            _face_vtx,
                                            triaA_vtxA);
      }
      else {
        /* Polygonal face */
        n_tria = PDM_triangulate_polygon(3,
                                         face_vtx_n,
                                         vtxA_coord,
                                         NULL,
                                         _face_vtx,
                                         PDM_TRIANGULATE_MESH_DEF,
                                         triaA_vtxA,
                                         tri_state);
      }


      /* Loop on triangles A */
      for (int triaA_id = 0; triaA_id < n_tria; triaA_id++) {

        int *_triaA_vtxA = triaA_vtxA + 3*triaA_id;
        if (dbg) {
          PDM_log_trace_array_int(_triaA_vtxA, 3, "    _triaA_vtxA : ");
        }

        /* Ignore if current triangle contains reference vertex */
        if (_triaA_vtxA[0] == ref_vtxA_id ||
            _triaA_vtxA[1] == ref_vtxA_id ||
            _triaA_vtxA[2] == ref_vtxA_id) {
            continue;
        }

        if (faceA_sign > 0) {
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            int vtxA_id = _triaA_vtxA[ivtx] - 1;
            memcpy(&tetraA_coord[3*(ivtx+1)], &vtxA_coord[3*vtxA_id], sizeof(double)*3);
          }
        }
        else {
          for (int ivtx = 0; ivtx < 3; ivtx++) {
            int vtxA_id = _triaA_vtxA[2-ivtx] - 1;
            memcpy(&tetraA_coord[3*(ivtx+1)], &vtxA_coord[3*vtxA_id], sizeof(double)*3);
          }
        }

        // log_trace("tetraA_coord :\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n",
        //           tetraA_coord[3*0], tetraA_coord[3*0+1], tetraA_coord[3*0+2],
        //           tetraA_coord[3*1], tetraA_coord[3*1+1], tetraA_coord[3*1+2],
        //           tetraA_coord[3*2], tetraA_coord[3*2+1], tetraA_coord[3*2+2],
        //           tetraA_coord[3*3], tetraA_coord[3*3+1], tetraA_coord[3*3+2]);


        double mat[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            mat[i][j] = tetraA_coord[3*(j+1)+i] - tetraA_coord[i];
          }
        }

        double det = mat[0][0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
        -            mat[1][0]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
        +            mat[2][0]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]);

        log_trace("det = %f\n", det);

        PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
        if (det == 0.) {
          continue;
        }
        PDM_GCC_SUPPRESS_WARNING_POP

        double idet = 1./det;


        /* Loop on candidate cells B */
        for (int icellB = cellA_cellB_idx[cellA_id]; icellB < cellA_cellB_idx[cellA_id+1]; icellB++) {
          int cellB_id = cellA_cellB[icellB];

          if (cellB_id != 7) continue;

          if (dbg) {
            log_trace("      icellB %d, cellB_id %d\n", icellB, cellB_id);
          }

          /* Loop on faces B */
          for (int ifaceB = cellB_faceB_idx[cellB_id]; ifaceB < cellB_faceB_idx[cellB_id+1]; ifaceB++) {

            int faceB_id   = PDM_ABS (cellB_faceB[ifaceB]) - 1;
            int faceB_sign = PDM_SIGN(cellB_faceB[ifaceB]);

            if (dbg) {
              log_trace("        faceB_id %d, sign = %d\n", faceB_id, faceB_sign);
            }

            /* Loop on triangles B */
            for (int triaB_id = faceB_triaB_idx[faceB_id]; triaB_id < faceB_triaB_idx[faceB_id+1]; triaB_id++) {

              int *_triaB_vtxB = triaB_vtxB + 3*triaB_id;

              if (dbg) {
                log_trace("          triaB_id %d : ", triaB_id);
                PDM_log_trace_array_int(_triaB_vtxB, 3, "");
              }

              if (faceB_sign > 0) {
                for (int ivtx = 0; ivtx < 3; ivtx++) {
                  int vtxB_id = _triaB_vtxB[ivtx] - 1;
                  memcpy(&triaB_coord[3*ivtx], &vtxB_coord[3*vtxB_id], sizeof(double)*3);
                }
              }
              else {
                for (int ivtx = 0; ivtx < 3; ivtx++) {
                  int vtxB_id = _triaB_vtxB[2-ivtx] - 1;
                  memcpy(&triaB_coord[3*ivtx], &vtxB_coord[3*vtxB_id], sizeof(double)*3);
                }
              }
              // log_trace("triaB_coord :\n%f %f %f\n%f %f %f\n%f %f %f\n",
              //           triaB_coord[3*0], triaB_coord[3*0+1], triaB_coord[3*0+2],
              //           triaB_coord[3*1], triaB_coord[3*1+1], triaB_coord[3*1+2],
              //           triaB_coord[3*2], triaB_coord[3*2+1], triaB_coord[3*2+2]);

              if (dbg && 0) {
                _dump_elementary_vol_vol(mi->comm,
                                         cellA_id,
                                         faceA_id,
                                         cellB_id,
                                         triaA_id,
                                         faceB_id,
                                         triaB_id,
                                         tetraA_coord,
                                         triaB_coord);
              }


              /* Transform triangle to current tetra's local frame */
              for (int ivtx = 0; ivtx < 3; ivtx++) {
                double rhs[3];
                for (int i = 0; i < 3; i++) {
                  rhs[i] = triaB_coord[3*ivtx + i] - tetraA_coord[i];
                }

                triaB_coord[3*ivtx    ] = idet *
                (  rhs[0]*(mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])
                 - rhs[1]*(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])
                 + rhs[2]*(mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2]));

                triaB_coord[3*ivtx + 1] = idet *
                (  mat[0][0]*(rhs[1]*mat[2][2] - rhs[2]*mat[1][2])
                 - mat[1][0]*(rhs[0]*mat[2][2] - rhs[2]*mat[0][2])
                 + mat[2][0]*(rhs[0]*mat[1][2] - rhs[1]*mat[0][2]));

                triaB_coord[3*ivtx + 2] = idet *
                (  mat[0][0]*(mat[1][1]*rhs[2] - mat[2][1]*rhs[1])
                 - mat[1][0]*(mat[0][1]*rhs[2] - mat[2][1]*rhs[0])
                 + mat[2][0]*(mat[0][1]*rhs[1] - mat[1][1]*rhs[0]));
              }

              if (dbg && 1) {
                // check inverse transform
                for (int i = 0; i < 3; i++) {
                  int vtxB_id = _triaB_vtxB[i] - 1;
                  if (faceB_sign < 0) {
                    vtxB_id = _triaB_vtxB[2-i] - 1;
                  }
                  double err = 0;
                  for (int j = 0; j < 3; j++) {
                    double u = triaB_coord[3*i  ];
                    double v = triaB_coord[3*i+1];
                    double w = triaB_coord[3*i+2];
                    double delta = vtxB_coord[3*vtxB_id+j] -
                    ((1-u-v-w) * tetraA_coord[  j] +
                     u         * tetraA_coord[3+j] +
                     v         * tetraA_coord[6+j] +
                     w         * tetraA_coord[9+j]);
                    err += delta*delta;
                  }

                  if (err > 1e-12) {
                    log_trace("!!! error = %e\n", sqrt(err));
                  }
                }
              }


              /* Perform elementary computation */
              double *local_vtx_coordA    = NULL;
              int     local_n_vtxA        = 0;
              int    *local_face_vtxA     = NULL;
              double *local_vtx_coordB    = NULL;
              int     local_n_vtxB        = 0;
              int    *local_face_vtxB     = NULL;

              // double volume = _reference_volume_compute(cll_storage,
              //                                           triaB_coord,
              //                                           &local_vtx_coordA,
              //                                           &local_n_vtxA,
              //                                           &local_face_vtxA,
              //                                           &vtk_n_faceA,
              //                                           &local_vtx_coordB,
              //                                           &local_n_vtxB,
              //                                           &local_face_vtxB,
              //                                           &vtk_n_faceB);
              double volume = 0;

              if (method == 0) {
                volume = PDM_mesh_intersection_vol_vol_atomic_compute(triaB_coord,
                                                                      &local_vtx_coordA,
                                                                      &local_n_vtxA,
                                                                      &local_face_vtxA,
                                                                      &vtk_n_faceA,
                                                                      &local_vtx_coordB,
                                                                      &local_n_vtxB,
                                                                      &local_face_vtxB,
                                                                      &vtk_n_faceB);
              }
              else {
                volume = PDM_mesh_intersection_vol_vol_atomic_compute2(triaB_coord);
              }


              if (dbg && 1) {
                log_trace("vtk_n_faceA = %d, vtk_n_faceB = %d\n", vtk_n_faceA, vtk_n_faceB);
              }

              // for vtk
              if (dbg_Karmijn) {
                if (local_n_vtxA > 0) {
                  printf("local_n_vtxA = %d\n", local_n_vtxA);
                  memcpy(vtk_vtx_coordA + 3* vtk_n_vtxA, local_vtx_coordA, sizeof(double) * 3 * local_n_vtxA);
                  vtk_face_vtx_idxA[vtk_n_faceA+1] = vtk_face_vtx_idxA[vtk_n_faceA] + local_n_vtxA;
                  printf("vtk_face_vtx_idxA[vtk_n_vtxA] = %d\n", vtk_face_vtx_idxA[vtk_n_vtxA]);
                  memcpy(vtk_face_vtxA + vtk_face_vtx_idxA[vtk_n_faceA], local_face_vtxA, sizeof(int) * local_n_vtxA);
                  vtk_n_vtxA += local_n_vtxA;
                }

                if (local_n_vtxB > 0) {
                  memcpy(vtk_vtx_coordB + 3* vtk_n_vtxB, local_vtx_coordB, sizeof(double) * 3 * local_n_vtxB);
                  vtk_face_vtx_idxB[vtk_n_faceB+1] = vtk_face_vtx_idxB[vtk_n_faceB] + local_n_vtxB;
                  memcpy(vtk_face_vtxB + vtk_face_vtx_idxB[vtk_n_faceB], local_face_vtxB, sizeof(int) * local_n_vtxB);
                  vtk_n_vtxB += local_n_vtxB;
                }
              }

              // free
              free(local_vtx_coordA);
              free(local_face_vtxA);
              free(local_vtx_coordB);
              free(local_face_vtxB);

              if (dbg) {
                log_trace("            volume Karmijn = %20.16f\n", volume);
                if (1) {//volume != 0) {
                  log_trace("********\n");
                  for (int i = 0; i < 3; i++) {
                    log_trace("double pt%d[3] = {%20.16f,%20.16f,%20.16f};\n",
                              i, triaB_coord[3*i], triaB_coord[3*i+1], triaB_coord[3*i+2]);
                  }
                  log_trace("********\n");
                }
              }

              /* Add elementray volume contribution */
              log_trace("sign(det) = %d / faceA_sign = %d\n", PDM_SIGN(det), faceA_sign);
              cellA_cellB_volume[icellB] += volume * PDM_ABS(det);

            } // End of loop on triangles of current face B

          } // End of loop on faces of current cell B

        } // End of loop on candidate cells B for current cell A

      } // End of loop on triangles of current face A

    } // End of loop on faces of current cell A


    if (dbg) {
      for (int icellB = cellA_cellB_idx[cellA_id]; icellB < cellA_cellB_idx[cellA_id+1]; icellB++) {
        int cellB_id = cellA_cellB[icellB];
        if (cellA_cellB_volume[icellB] > 0.00001) {log_trace("---> BIG");}
        log_trace("cellA %d ("PDM_FMT_G_NUM") cellB %d ("PDM_FMT_G_NUM"), volume = %20.16f\n",
                  cellA_id, cellA_ln_to_gn[cellA_id],
                  cellB_id, cellB_ln_to_gn[cellB_id],
                  cellA_cellB_volume[icellB]);
      }
    }

  } // End of loop on cells A

  // for vtk
  if (dbg_Karmijn) {
    if (vtk_n_vtxA > 0) {
      char filename[999] = "A.vtk";
      PDM_vtk_write_polydata(filename,
                             vtk_n_vtxA,
                             vtk_vtx_coordA,
                             NULL,
                             vtk_n_faceA,
                             vtk_face_vtx_idxA,
                             vtk_face_vtxA,
                             NULL,
                             NULL);

      free(vtk_vtx_coordA);
      free(vtk_face_vtx_idxA);
      free(vtk_face_vtxA);
    }

    if (vtk_n_vtxB > 0) {
      char filename[999] = "B.vtk";
      PDM_vtk_write_polydata(filename,
                             vtk_n_vtxB,
                             vtk_vtx_coordB,
                             NULL,
                             vtk_n_faceB,
                             vtk_face_vtx_idxB,
                             vtk_face_vtxB,
                             NULL,
                             NULL);

      free(vtk_vtx_coordB);
      free(vtk_face_vtx_idxB);
      free(vtk_face_vtxB);
    }
  }

  // free
  // for (int i = 0; i < max_size; i++) {
  //   free(cll_storage[i]);
  // }
  // free(cll_storage);

  free(faceB_triaB_idx);
  free(triaB_vtxB);
  free(triaA_vtxA);
  PDM_triangulate_state_destroy(tri_state);



  if (dbg) {
    // Crude check
    double l_total_volume_AB = 0;
    for (int i = 0; i < cellA_cellB_idx[n_cellA]; i++) {
      l_total_volume_AB += cellA_cellB_volume[i];
    }

    double l_total_volume_A  = 0;
    // double *cellA_volume = malloc(sizeof(double) * n_cellA);
    // double *cellA_center = malloc(sizeof(double) * n_cellA * 3);
    PDM_geom_elem_polyhedra_properties(1,
                                       n_cellA,
                                       n_faceA,
                                       faceA_vtxA_idx,
                                       faceA_vtxA,
                                       cellA_faceA_idx,
                                       cellA_faceA,
                                       n_vtxA,
                                       vtxA_coord,
                                       cellA_volume,
                                       cellA_center,
                                       NULL,
                                       NULL);
    // free(cellA_center);

    for (int cellA_id = 0; cellA_id < n_cellA; cellA_id++) {
      l_total_volume_A += cellA_volume[cellA_id];
    }


    double g_total_volume_AB;
    PDM_MPI_Allreduce(&l_total_volume_AB, &g_total_volume_AB, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    double g_total_volume_A;
    PDM_MPI_Allreduce(&l_total_volume_A, &g_total_volume_A, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    log_trace("total volume of A inter B : local = %20.16f, global = %20.16f (%3.3f%%)\n",
              l_total_volume_AB, g_total_volume_AB,
              100*g_total_volume_AB / g_total_volume_A);

    // cas cube, translation (0.5,0.5,0.5)
    double exact = 0.5;
    log_trace("error : absolute = %e, relative = %e\n",
              PDM_ABS(g_total_volume_AB - exact),
              PDM_ABS(g_total_volume_AB - exact)/exact);
  }

  free(cellA_center);
  free(cellA_volume);


  if (owner_face_vtxA == PDM_OWNERSHIP_USER) {
    free(faceA_vtxA);
  }

  if (owner_face_vtxB == PDM_OWNERSHIP_USER) {
    free(faceB_vtxB);
  }

  // Do not free...
  free(cellA_cellB_volume);
}

static
void
_mesh_intersection_vol_surf
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);

  /*
   * Panic vtk
   */
  if(0 == 1) {
    _export_vtk_3d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_2d("extrp_mesh_b", extrp_mesh_b);
  }

}



static void
_get_extracted_mesh_surf
(
 PDM_extract_part_t  *extrp,
 int                 *n_face,
 int                 *n_edge,
 int                 *n_vtx,
 int                **face_edge_idx,
 int                **face_edge,
 int                **edge_vtx,
 double             **vtx_coord,
 PDM_g_num_t        **face_ln_to_gn
 )
 {
  *n_face = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                              face_edge,
                                              face_edge_idx,
                                              PDM_OWNERSHIP_KEEP);
  int *edge_vtx_idx = NULL;
  *n_edge = PDM_extract_part_connectivity_get(extrp,
                                              0,
                                              PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                              edge_vtx,
                                              &edge_vtx_idx,
                                              PDM_OWNERSHIP_KEEP);

  *n_vtx = PDM_extract_part_vtx_coord_get(extrp,
                                          0,
                                          vtx_coord,
                                          PDM_OWNERSHIP_KEEP);

  PDM_extract_part_parent_ln_to_gn_get(extrp,
                                       0,
                                       PDM_MESH_ENTITY_FACE,
                                       face_ln_to_gn,
                                       PDM_OWNERSHIP_KEEP);
 }


static inline void
_vector_ab
(
      double ab[3],
const double a[3],
const double b[3]
)
{
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];
}

static void
_polygon_geom_properties
(
 int     n_edge,
 int    *face_edge,
 int    *edge_vtx,
 double *vtx_coord,
 double *normal,
 double *barycenter
 )
{
  for (int i = 0; i < 3; i++) {
    normal    [i] = 0;
    barycenter[i] = 0;
  }

  for (int iedge = 0; iedge < n_edge; iedge++) {
    int edge_id   = PDM_ABS (face_edge[iedge]) - 1;

    int vtx_id0 = edge_vtx[2*edge_id  ] - 1;
    int vtx_id1 = edge_vtx[2*edge_id+1] - 1;

    for (int i = 0; i < 3; i++) {
      barycenter[i] += vtx_coord[3*vtx_id0+i] + vtx_coord[3*vtx_id1+i];
    }
  }

  double normalization = 1./(2. * n_edge);

  for (int i = 0; i < 3; i++) {
    barycenter[i] *= normalization;
  }

  for (int iedge = 0; iedge < n_edge; iedge++) {
    int edge_id   = PDM_ABS (face_edge[iedge]) - 1;
    int edge_sign = PDM_SIGN(face_edge[iedge]);

    double vec[2][3];
    for (int j = 0; j < 2; j++) {
      int vtx_id = edge_vtx[2*edge_id+j] - 1;
      _vector_ab(vec[j], barycenter, &vtx_coord[3*vtx_id]);
    }

    double cross[3];
    PDM_CROSS_PRODUCT(cross, vec[0], vec[1]);
    for (int i = 0; i < 3; i++) {
      normal[i] += edge_sign * cross[i];
    }
  }
}


PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
static inline void
_clip1
(
       double *uc,
       double *vc,
 const double  ud,
 const double  vd
 )
{
  if (*uc < 0) {
    if (*uc == ud) {
      *vc = 0;
    }
    else {
      *vc -= (*uc)*(vd - (*vc))/(ud - (*uc));
    }
    *uc = 0;
  }

  if (*vc < 0) {
    if (*vc == vd) {
      *uc = 0;
    }
    else {
      *uc -= (*vc)*(ud - (*uc))/(vd - (*vc));
    }
    *vc = 0;
  }
}

static inline void
_clip2
(
 double *u,
 double *v
 )
{
  double w = (*u) + (*v);
  if (w > 1) {
    double iw = 1./w;
    *u *= iw;
    *v *= iw;
  }
}

static
void
_mesh_intersection_surf_surf
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  int dbg = 1;

  /*
   * Panic vtk
   */
  if (dbg) {
    _export_vtk_2d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_2d("extrp_mesh_b", extrp_mesh_b);
  }

  /*
   * TODO : Faire la géometrie :D
   */
  int *faceA_faceB_idx = redistribute_box_a_to_box_b_idx;
  int *faceA_faceB     = redistribute_box_a_to_box_b;

  /* Get connectivities and coordinates */
  int          n_faceA         = 0;
  int          n_edgeA         = 0;
  int          n_vtxA          = 0;
  int         *faceA_edgeA_idx = NULL;
  int         *faceA_edgeA     = NULL;
  int         *edgeA_vtxA      = NULL;
  double      *vtxA_coord      = NULL;
  PDM_g_num_t *faceA_ln_to_gn  = NULL;
  _get_extracted_mesh_surf(extrp_mesh_a,
                           &n_faceA,
                           &n_edgeA,
                           &n_vtxA,
                           &faceA_edgeA_idx,
                           &faceA_edgeA,
                           &edgeA_vtxA,
                           &vtxA_coord,
                           &faceA_ln_to_gn);

  int          n_faceB         = 0;
  int          n_edgeB         = 0;
  int          n_vtxB          = 0;
  int         *faceB_edgeB_idx = NULL;
  int         *faceB_edgeB     = NULL;
  int         *edgeB_vtxB      = NULL;
  double      *vtxB_coord      = NULL;
  PDM_g_num_t *faceB_ln_to_gn  = NULL;
  _get_extracted_mesh_surf(extrp_mesh_b,
                           &n_faceB,
                           &n_edgeB,
                           &n_vtxB,
                           &faceB_edgeB_idx,
                           &faceB_edgeB,
                           &edgeB_vtxB,
                           &vtxB_coord,
                           &faceB_ln_to_gn);


  double *faceA_faceB_volume = malloc(sizeof(double) * faceA_faceB_idx[n_faceA]);

  for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {

    double faceA_normal[3];
    double faceA_center[3];
    _polygon_geom_properties(faceA_edgeA_idx[faceA_id+1] - faceA_edgeA_idx[faceA_id],
                             faceA_edgeA + faceA_edgeA_idx[faceA_id],
                             edgeA_vtxA,
                             vtxA_coord,
                             faceA_normal,
                             faceA_center);

    /* Unit normal */
    double mag = PDM_DOT_PRODUCT(faceA_normal, faceA_normal);
    if (mag <= 0) {
      // degenerate polygon
      continue;
    }

    double imag = 1./sqrt(mag);
    for (int i = 0; i < 3; i++) {
      faceA_normal[i] *= imag;
    }

    for (int ifaceB = faceA_faceB_idx[faceA_id]; ifaceB < faceA_faceB_idx[faceA_id+1]; ifaceB++) {
      int faceB_id = faceA_faceB[ifaceB];

      int dbg_pair = 0;// dbg && (faceA_id == 5 && faceB_id == 0);

      double area = 0.;

      double faceB_normal[3];// Compute only once and store?
      double faceB_center[3];// Not used
      _polygon_geom_properties(faceB_edgeB_idx[faceB_id+1] - faceB_edgeB_idx[faceB_id],
                               faceB_edgeB + faceB_edgeB_idx[faceB_id],
                               edgeB_vtxB,
                               vtxB_coord,
                               faceB_normal,
                               faceB_center);

      int signAB = (int) PDM_SIGN(PDM_DOT_PRODUCT(faceA_normal, faceB_normal));
      if (dbg_pair) {
        log_trace("faceA %d ("PDM_FMT_G_NUM") faceB %d ("PDM_FMT_G_NUM"), signAB = %d\n",
                  faceA_id, faceA_ln_to_gn[faceA_id], faceB_id, faceB_ln_to_gn[faceB_id], signAB);
        log_trace("faceA_center = %f %f %f\n", faceA_center[0], faceA_center[1], faceA_center[2]);
      }

      for (int iedgeA = faceA_edgeA_idx[faceA_id]; iedgeA < faceA_edgeA_idx[faceA_id+1]; iedgeA++) {
        int edgeA_id   = PDM_ABS (faceA_edgeA[iedgeA]) - 1;
        int edgeA_sign = PDM_SIGN(faceA_edgeA[iedgeA]);

        int vtxA_id0 = edgeA_vtxA[2*edgeA_id  ] - 1;
        int vtxA_id1 = edgeA_vtxA[2*edgeA_id+1] - 1;
        if (edgeA_sign < 0) {
            int tmp = vtxA_id0;
            vtxA_id0 = vtxA_id1;
            vtxA_id1 = tmp;
          }

        if (dbg_pair) {
          log_trace("  edgeA (%d)%d: %d %d\n", edgeA_sign, edgeA_id, vtxA_id0, vtxA_id1);
        }

        double *a = vtxA_coord + 3*vtxA_id0;
        double *b = vtxA_coord + 3*vtxA_id1;

        double ka[3], kb[3];
        _vector_ab(ka, faceA_center, a);
        _vector_ab(kb, faceA_center, b);

        double kaka = PDM_DOT_PRODUCT(ka, ka);
        double kakb = PDM_DOT_PRODUCT(ka, kb);
        double kbkb = PDM_DOT_PRODUCT(kb, kb);

        double det = kaka*kbkb - kakb*kakb;

        if (det <= 0.) {
          // points k, a and b are collinear, skip edge ab
          continue;
        }

        double idet = 1./det;

        double normal_kab[3];
        PDM_CROSS_PRODUCT(normal_kab, ka, kb);

        double area_kab = 0.5*PDM_DOT_PRODUCT(normal_kab, faceA_normal);


        for (int iedgeB = faceB_edgeB_idx[faceB_id]; iedgeB < faceB_edgeB_idx[faceB_id+1]; iedgeB++) {
          int edgeB_id   = PDM_ABS (faceB_edgeB[iedgeB]) - 1;
          int edgeB_sign = PDM_SIGN(faceB_edgeB[iedgeB]);

          int vtxB_id0 = edgeB_vtxB[2*edgeB_id  ] - 1;
          int vtxB_id1 = edgeB_vtxB[2*edgeB_id+1] - 1;
          if (edgeB_sign < 0) {
            int tmp = vtxB_id0;
            vtxB_id0 = vtxB_id1;
            vtxB_id1 = tmp;
          }

          if (dbg_pair) {
            log_trace("    edgeB (%d)%d: %d %d\n", edgeB_sign, edgeB_id, vtxB_id0, vtxB_id1);
          }

          double *c = vtxB_coord + 3*vtxB_id0;
          double *d = vtxB_coord + 3*vtxB_id1;

          /* Compute the barycentric coordinates of c and d in kab */
          double kc[3], kd[3];
          _vector_ab(kc, faceA_center, c);
          _vector_ab(kd, faceA_center, d);

          double kcka = PDM_DOT_PRODUCT(kc, ka);
          double kckb = PDM_DOT_PRODUCT(kc, kb);
          double kdka = PDM_DOT_PRODUCT(kd, ka);
          double kdkb = PDM_DOT_PRODUCT(kd, kb);

          double uc = kcka*kbkb - kckb*kakb;
          double ud = kdka*kbkb - kdkb*kakb;

          if (dbg_pair) {
            log_trace("      uc = %f, ud = %f\n", uc, ud);
          }

          if (uc <= 0 && ud <= 0) {
            continue;
          }

          double vc = kaka*kckb - kakb*kcka;
          double vd = kaka*kdkb - kakb*kdka;

          if (dbg_pair) {
            log_trace("      vc = %f, vd = %f\n", vc, vd);
          }

          if (vc <= 0 && vd <= 0) {
            continue;
          }

          if (uc*vd - vc*ud == 0) {
            // points k, c and d are collinear, skip edge cd
            continue;
          }

          uc *= idet;
          vc *= idet;
          ud *= idet;
          vd *= idet;

          /* Compute intersection between ab and cd */
          int intersect = 0;
          double m00 = -1;
          double m01 = uc - ud;
          double m10 = 1;
          double m11 = vc - vd;
          double s;

          double det2 = m00*m11 - m01*m10;
          if (det2 != 0) {
            double idet2 = 1./det2;
            double r0 = uc - 1;
            double r1 = vc;
            s        = (r0*m11 - r1*m01) * idet2;
            double t = (m00*r1 - m10*r0) * idet2;
            if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
              intersect = 1;
            }
          }

          /* Clip kcd by 'quarter space' {u,v >= 0} */
          _clip1(&uc, &vc, ud, vd);
          _clip1(&ud, &vd, uc, vc);

          /* Clip kcd by triangle kab {u+v <= 1} */
          _clip2(&uc, &vc);
          _clip2(&ud, &vd);

          if (dbg_pair) {
            log_trace("      final uc = %f, vc = %f\n", uc, vc);
            log_trace("      final ud = %f, vd = %f\n", ud, vd);
          }

          /* Add contribution */
          double f = 0;
          if (intersect) {
            f = (1-s)*(vd - vc) + s*(uc - ud);
            if (dbg_pair) {
              log_trace("    c2 = %f %f %f\n",
                        faceA_center[0] + uc*ka[0] + vc*kb[0],
                        faceA_center[1] + uc*ka[1] + vc*kb[1],
                        faceA_center[2] + uc*ka[2] + vc*kb[2]);
              log_trace("    x  = %f %f %f\n",
                        faceA_center[0] + (1-s)*ka[0] + s*kb[0],
                        faceA_center[1] + (1-s)*ka[1] + s*kb[1],
                        faceA_center[2] + (1-s)*ka[2] + s*kb[2]);
              log_trace("    d2 = %f %f %f\n",
                        faceA_center[0] + ud*ka[0] + vd*kb[0],
                        faceA_center[1] + ud*ka[1] + vd*kb[1],
                        faceA_center[2] + ud*ka[2] + vd*kb[2]);
            }
          }
          else {
            f = uc*vd - vc*ud;
            if (dbg_pair) {
              log_trace("    c2 = %f %f %f\n",
                        faceA_center[0] + uc*ka[0] + vc*kb[0],
                        faceA_center[1] + uc*ka[1] + vc*kb[1],
                        faceA_center[2] + uc*ka[2] + vc*kb[2]);
              log_trace("    d2 = %f %f %f\n",
                        faceA_center[0] + ud*ka[0] + vd*kb[0],
                        faceA_center[1] + ud*ka[1] + vd*kb[1],
                        faceA_center[2] + ud*ka[2] + vd*kb[2]);
            }
          }

          if (dbg_pair) {
            log_trace("    f = %f\n", f);
            log_trace("    area += %f\n", f*area_kab);
          }
          area += f*area_kab;

        } // End of loop on current faceB's edges

      } // End of loop on current faceA's edges

      faceA_faceB_volume[ifaceB] = signAB * area;
      if (0) {//dbg) {
        log_trace("faceA %d ("PDM_FMT_G_NUM") faceB %d ("PDM_FMT_G_NUM"), volume = %20.16f (%3.3f%)\n",
                  faceA_id, faceA_ln_to_gn[faceA_id],
                  faceB_id, faceB_ln_to_gn[faceB_id],
                  faceA_faceB_volume[ifaceB],
                  100*faceA_faceB_volume[ifaceB]*imag*2);
      }

    } // End of loop on faces B

  } // End of loop on faces A



  if (dbg) {
    // Crude check
    double l_total_area_AB = 0;
    for (int i = 0; i < faceA_faceB_idx[n_faceA]; i++) {
      l_total_area_AB += faceA_faceB_volume[i];
    }

    double l_total_area_A  = 0;
    for (int faceA_id = 0; faceA_id < n_faceA; faceA_id++) {
      double faceA_normal[3];
      double faceA_center[3];
      _polygon_geom_properties(faceA_edgeA_idx[faceA_id+1] - faceA_edgeA_idx[faceA_id],
                               faceA_edgeA + faceA_edgeA_idx[faceA_id],
                               edgeA_vtxA,
                               vtxA_coord,
                               faceA_normal,
                               faceA_center);

      l_total_area_A += 0.5*PDM_MODULE(faceA_normal);
    }

    double g_total_area_AB;
    PDM_MPI_Allreduce(&l_total_area_AB, &g_total_area_AB, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    double g_total_area_A;
    PDM_MPI_Allreduce(&l_total_area_A, &g_total_area_A, 1,
                      PDM_MPI_DOUBLE, PDM_MPI_SUM, mi->comm);

    log_trace("total area of A inter B : local = %20.16f, global = %20.16f (%3.3f%%)\n",
              l_total_area_AB, g_total_area_AB,
              100*g_total_area_AB / g_total_area_A);

    // cas plan, translation (0.5,0.5,0) + rotation PI/5
    double exact = 0.0875401518835469;
    log_trace("error : absolute = %e, relative = %e\n",
              PDM_ABS(g_total_area_AB - exact),
              PDM_ABS(g_total_area_AB - exact)/exact);
  }

  // Do not free...
  free(faceA_faceB_volume);

}
PDM_GCC_SUPPRESS_WARNING_POP


static
void
_mesh_intersection_surf_line
(
 PDM_mesh_intersection_t *mi,
 PDM_extract_part_t      *extrp_mesh_a,
 PDM_extract_part_t      *extrp_mesh_b,
 int                     *redistribute_box_a_to_box_b_idx,
 int                     *redistribute_box_a_to_box_b
)
{
  PDM_UNUSED(mi);
  PDM_UNUSED(extrp_mesh_a);
  PDM_UNUSED(extrp_mesh_b);
  PDM_UNUSED(redistribute_box_a_to_box_b_idx);
  PDM_UNUSED(redistribute_box_a_to_box_b);
  if(0 == 1) {
    _export_vtk_2d("extrp_mesh_a", extrp_mesh_a);
    _export_vtk_1d("extrp_mesh_b", extrp_mesh_b);
  }

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_mesh_intersection_t*
PDM_mesh_intersection_create
(
 const PDM_mesh_intersection_kind_t intersection_kind,
 const int                          dim_mesh_a,
 const int                          dim_mesh_b,
 const int                          n_part_mesh_a,
 const int                          n_part_mesh_b,
 const double                       project_coeff,
       PDM_MPI_Comm                 comm
)
{
  PDM_mesh_intersection_t *mi = (PDM_mesh_intersection_t *) malloc(sizeof(PDM_mesh_intersection_t));

  mi->comm = comm;
  mi->intersect_kind = intersection_kind;
  mi->n_part_mesh_a  = n_part_mesh_a;
  mi->n_part_mesh_b  = n_part_mesh_b;
  mi->dim_mesh_a     = dim_mesh_a;
  mi->dim_mesh_b     = dim_mesh_b;
  mi->project_coef   = project_coeff;

  mi->mesh_a = PDM_part_mesh_create(n_part_mesh_a, comm);
  mi->mesh_b = PDM_part_mesh_create(n_part_mesh_b, comm);

  return mi;
}

void
PDM_mesh_intersection_compute
(
  PDM_mesh_intersection_t  *mi
)
{
  /*
   * Compute extents of mesh_a and mesh_b
   */
  double **extents_mesh_a = NULL;
  double **extents_mesh_b = NULL;
  double global_extents[6] = { -HUGE_VAL, -HUGE_VAL, -HUGE_VAL, HUGE_VAL,  HUGE_VAL,  HUGE_VAL};
  double mesh_global_extents[2][6] = {{ HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL},
                                      {HUGE_VAL,  HUGE_VAL,  HUGE_VAL,
                                       -HUGE_VAL, -HUGE_VAL, -HUGE_VAL}};
  double g_mesh_global_extents[2][6];
  double g_global_extents        [6];
  _compute_part_mesh_extents(mi->mesh_a, mi->dim_mesh_a, mesh_global_extents[0], &extents_mesh_a);
  _compute_part_mesh_extents(mi->mesh_b, mi->dim_mesh_b, mesh_global_extents[1], &extents_mesh_b);

  /*
   * Global extents exchange
   */
  const int dim = 3;
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh], g_mesh_global_extents[i_mesh], dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MIN, mi->comm);
    PDM_MPI_Allreduce(mesh_global_extents[i_mesh]+dim, g_mesh_global_extents[i_mesh]+dim, dim,
                      PDM_MPI_DOUBLE, PDM_MPI_MAX, mi->comm);
  }

  /* Union or intersection of global extents */
  for(int i_mesh = 0; i_mesh < 2; ++i_mesh) {
    for (int k = 0; k < 3; k++) {
      // Union
      // global_extents[k]     = PDM_MIN(mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      // global_extents[3 + k] = PDM_MAX(mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
      // Intersection
      global_extents[k]     = PDM_MAX(g_mesh_global_extents[i_mesh][k  ], global_extents[k  ]);
      global_extents[3 + k] = PDM_MIN(g_mesh_global_extents[i_mesh][3+k], global_extents[3+k]);
    }
  }
  for(int i = 0; i < 6; ++i) {
    g_global_extents[i] = global_extents[i];
  }
  double max_range = -HUGE_VAL;
  double min_range =  HUGE_VAL;

  for (int k = 0; k < dim; k++) {
    max_range = PDM_MAX(max_range, (g_global_extents[dim+k] - g_global_extents[k]));
    min_range = PDM_MIN(min_range, (g_global_extents[dim+k] - g_global_extents[k]));
  }

  for (int k = 0; k < dim; k++) {
    g_global_extents[k]     += -max_range * 1.1e-3; // On casse la symetrie !
    g_global_extents[dim+k] +=  max_range * 1e-3;
  }


  /*
   * Extraction - In option ?
   */
  int           n_mesh = 2;
  int           n_part                    [n_mesh];
  int          *n_extract_elmt            [n_mesh];
  double      **extract_box_extents       [n_mesh];
  int         **extract_elmt_init_location[n_mesh];
  PDM_g_num_t **extract_elmt_ln_to_gn     [n_mesh];

  n_part[0] = mi->n_part_mesh_a;
  n_part[1] = mi->n_part_mesh_b;

  _select_elements_by_global_bbox(mi->mesh_a, mi->dim_mesh_a,
                                  extents_mesh_a, g_mesh_global_extents[1], // On enleve tout ce qui est en dehors de B
                                  &n_extract_elmt[0],
                                  &extract_box_extents[0],
                                  &extract_elmt_init_location[0],
                                  &extract_elmt_ln_to_gn[0]);
  _select_elements_by_global_bbox(mi->mesh_b, mi->dim_mesh_b,
                                  extents_mesh_b, g_mesh_global_extents[0], // On enleve tout ce qui est en dehors de A
                                  &n_extract_elmt[1],
                                  &extract_box_extents[1],
                                  &extract_elmt_init_location[1],
                                  &extract_elmt_ln_to_gn[1]);

  for(int i_part = 0; i_part < mi->n_part_mesh_a; ++i_part) {
    free(extents_mesh_a[i_part]);
  }
  for(int i_part = 0; i_part < mi->n_part_mesh_b; ++i_part) {
    free(extents_mesh_b[i_part]);
  }
  free(extents_mesh_a);
  free(extents_mesh_b);

  // Attention le dbtree fait  le init_location sauf que la il faut le forcer !!!!!!
  PDM_dbbtree_t *dbbtree_mesh_a = PDM_dbbtree_create (mi->comm, dim, g_global_extents);

  PDM_box_set_t  *boxes_mesh_a = PDM_dbbtree_boxes_set_with_init_location(dbbtree_mesh_a,
                                                                          mi->mesh_a->n_part,
                                                                          n_extract_elmt            [0],
                                                  (const int         **)  extract_elmt_init_location[0],
                                                  (const double      **)  extract_box_extents       [0],
                                                  (const PDM_g_num_t **)  extract_elmt_ln_to_gn     [0]);

  /*
   * Intersect with B
   */
  int *box_a_to_box_b_idx = NULL;
  int *box_a_to_box_b     = NULL;
  PDM_box_set_t  *boxes_mesh_b =  PDM_dbbtree_intersect_boxes_with_init_location_set(dbbtree_mesh_a,
                                                                                     mi->mesh_b->n_part,
                                                                                     n_extract_elmt            [1],
                                                              (const int         **) extract_elmt_init_location[1],
                                                              (const double      **) extract_box_extents       [1],
                                                              (const PDM_g_num_t **) extract_elmt_ln_to_gn     [1],
                                                                                     &box_a_to_box_b_idx,
                                                                                     &box_a_to_box_b);

  /* Free extraction */
  for(int i_mesh = 0; i_mesh < n_mesh; ++i_mesh) {
    for(int i_part = 0; i_part < n_part[i_mesh]; ++i_part) {
      free(extract_elmt_init_location[i_mesh][i_part]);
      free(extract_box_extents       [i_mesh][i_part]);
      free(extract_elmt_ln_to_gn     [i_mesh][i_part]);
    }
    free(n_extract_elmt            [i_mesh]);
    free(extract_elmt_init_location[i_mesh]);
    free(extract_box_extents       [i_mesh]);
    free(extract_elmt_ln_to_gn     [i_mesh]);
  }


  /*
   *  Redistrib all boxes (inplace) like overlay before extracting mesh
   */
  int *redistribute_box_a_to_box_b_idx = NULL;
  int *redistribute_box_a_to_box_b     = NULL;
  _redistrib_boxes(mi->comm,
                   boxes_mesh_a,
                   boxes_mesh_b,
                   box_a_to_box_b_idx,
                   box_a_to_box_b,
                   &redistribute_box_a_to_box_b_idx,
                   &redistribute_box_a_to_box_b);
  free(box_a_to_box_b_idx);
  free(box_a_to_box_b);


  /*
   * Extract part
   */
  PDM_extract_part_t* extrp_mesh_a = _create_extract_part(mi->mesh_a,
                                                          mi->dim_mesh_a,
                                                          boxes_mesh_a);
  PDM_extract_part_t* extrp_mesh_b = _create_extract_part(mi->mesh_b,
                                                          mi->dim_mesh_b,
                                                          boxes_mesh_b);

  PDM_dbbtree_free (dbbtree_mesh_a);
  // PDM_box_set_destroy (&boxes_mesh_a);
  // PDM_box_set_destroy (&boxes_mesh_b);

  /*
   * Geometry begin here ...
   */
  if(mi->dim_mesh_a == 3 && mi->dim_mesh_b == 3) {
    _mesh_intersection_vol_vol(mi,
                               extrp_mesh_a,
                               extrp_mesh_b,
                               redistribute_box_a_to_box_b_idx,
                               redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 3 && mi->dim_mesh_b == 2) {
    // On suppose que l'utilisateur met A = Vol et B = Surf
    _mesh_intersection_vol_surf(mi,
                                extrp_mesh_a,
                                extrp_mesh_b,
                                redistribute_box_a_to_box_b_idx,
                                redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 2 && mi->dim_mesh_b == 2) {
    _mesh_intersection_surf_surf(mi,
                                 extrp_mesh_a,
                                 extrp_mesh_b,
                                 redistribute_box_a_to_box_b_idx,
                                 redistribute_box_a_to_box_b);
  } else if(mi->dim_mesh_a == 2 && mi->dim_mesh_b == 1) {
    // On suppose que l'utilisateur met A = Vol et B = Surf
    _mesh_intersection_surf_line(mi,
                                 extrp_mesh_a,
                                 extrp_mesh_b,
                                 redistribute_box_a_to_box_b_idx,
                                 redistribute_box_a_to_box_b);
  } else {
    PDM_error(__FILE__, __LINE__, 0,
              "PDM_mesh_intersection_compute error : Cannot handle meshA with dim = %i and meshB = %i \n", mi->dim_mesh_a, mi->dim_mesh_b);
  }
  PDM_box_set_destroy (&boxes_mesh_a);
  PDM_box_set_destroy (&boxes_mesh_b);

  free(redistribute_box_a_to_box_b_idx);
  free(redistribute_box_a_to_box_b    );

  PDM_extract_part_free(extrp_mesh_a);
  PDM_extract_part_free(extrp_mesh_b);

}


void
PDM_mesh_intersection_part_set
(
  PDM_mesh_intersection_t  *mi,
  PDM_ol_mesh_t             i_mesh,
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
)
{
  PDM_part_mesh_t* mesh = mi->mesh_a;
  if(i_mesh == PDM_OL_MESH_B) {
    mesh = mi->mesh_b;
  }

  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_CELL  , n_cell);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_FACE  , n_face);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_EDGE  , n_edge);
  PDM_part_mesh_n_entity_set(mesh, i_part, PDM_MESH_ENTITY_VERTEX, n_vtx );

  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_CELL_FACE, cell_face, cell_face_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_EDGE, face_edge, face_edge_idx, PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_FACE_VTX , face_vtx , face_vtx_idx , PDM_OWNERSHIP_USER);
  PDM_part_mesh_connectivity_set(mesh, i_part, PDM_CONNECTIVITY_TYPE_EDGE_VTX , edge_vtx , NULL         , PDM_OWNERSHIP_USER);

  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_CELL  , cell_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_FACE  , face_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_EDGE  , edge_ln_to_gn, PDM_OWNERSHIP_USER);
  PDM_part_mesh_entity_ln_to_gn_set(mesh, i_part, PDM_MESH_ENTITY_VERTEX, vtx_ln_to_gn , PDM_OWNERSHIP_USER);

  PDM_part_mesh_vtx_coord_set(mesh, i_part, vtx_coord, PDM_OWNERSHIP_USER);
}



void
PDM_mesh_intersection_free
(
 PDM_mesh_intersection_t* mi
)
{

  PDM_part_mesh_free(mi->mesh_a);
  PDM_part_mesh_free(mi->mesh_b);

  free(mi);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
