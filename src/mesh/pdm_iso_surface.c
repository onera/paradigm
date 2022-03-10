/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_timer.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_gnum.h"
#include "pdm_iso_surface.h"
#include "pdm_iso_surface_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*=============================================================================
 * Private function definitions
 *============================================================================*/


static
void
_iso_line_dist
(
  PDM_iso_surface_t        *isos
)
{
  PDM_UNUSED(isos);

}

static
void
_iso_surf_dist
(
  PDM_iso_surface_t        *isos
)
{
  PDM_UNUSED(isos);

}

static
void
_iso_surface_dist
(
  PDM_iso_surface_t        *isos
)
{
  PDM_UNUSED(isos);

  /*
   * Select gnum that contains iso-surface
   */
  int i_rank;
  PDM_MPI_Comm_rank(isos->comm, &i_rank);
  assert(isos->distrib_edge != NULL);
  int dn_edge = isos->distrib_edge[i_rank+1] - isos->distrib_edge[i_rank];

  PDM_g_num_t* edge_ln_to_gn = (PDM_g_num_t * ) malloc( dn_edge * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_edge; ++i) {
    edge_ln_to_gn[i] = isos->distrib_edge[i_rank] + i + 1;
  }


  int          pn_vtx           = 0;
  PDM_g_num_t *pvtx_ln_to_gn    = NULL;
  int         *pedge_vtx_idx    = NULL;
  int         *pedge_vtx        = NULL;

  int* dedge_vtx_idx = malloc( (dn_edge + 1) * sizeof(int));
  dedge_vtx_idx[0] = 0;
  for(int i = 0; i < dn_edge; ++i) {
    dedge_vtx_idx[i+1] = dedge_vtx_idx[i] + 1;
  }

  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           isos->distrib_edge,
                                                           dedge_vtx_idx,
                                                           isos->dedge_vtx,
                                                           dn_edge,
                                     (const PDM_g_num_t *) edge_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pedge_vtx_idx,
                                                           &pedge_vtx);
  free(pedge_vtx_idx);
  free(edge_ln_to_gn);

  PDM_block_to_part_t* btp_vtx = PDM_block_to_part_create(isos->distrib_vtx,
                                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                                          &pn_vtx,
                                                          1,
                                                          isos->comm);

  int cst_stride = 1;
  double **tmp_pvtx_coord = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
            (void *  )   isos->dvtx_coord,
            (int  ***)   NULL,
            (void ***)   tmp_pvtx_coord);
  double* pvtx_coord = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);

  double **tmp_pfield = NULL;
  PDM_block_to_part_exch(btp_vtx,
                         sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &cst_stride,
            (void *  )   isos->dfield,
            (int  ***)   NULL,
            (void ***)   tmp_pfield);
  double* pfield = tmp_pfield[0];
  free(tmp_pfield);


  PDM_block_to_part_free(btp_vtx);

  /*
   *  Loop on edge to tag all edge
   */
  int    *dedge_tag    = (int    * ) malloc(    dn_edge * sizeof(int   ));
  double *dedge_center = (double * ) malloc(3 * dn_edge * sizeof(double));
  for(int i = 0; i < dn_edge; ++i) {

    int i_vtx1 = pedge_vtx[2*i  ]-1;
    int i_vtx2 = pedge_vtx[2*i+1]-1;

    dedge_tag[i] = 0;

    // Besoin des coordonnés si call back
    double x1 = pvtx_coord[3*i_vtx1  ];
    double y1 = pvtx_coord[3*i_vtx1+1];
    double z1 = pvtx_coord[3*i_vtx1+2];

    double x2 = pvtx_coord[3*i_vtx2  ];
    double y2 = pvtx_coord[3*i_vtx2+1];
    double z2 = pvtx_coord[3*i_vtx1+2];

    double val1 = pfield[i_vtx1];
    double val2 = pfield[i_vtx2];

    int sgn1 = PDM_SIGN(val1);
    int sgn2 = PDM_SIGN(val2);

    if(sgn1 * sgn2 < 0) {
      dedge_tag[i] = 1;
    }

    dedge_center[3*i  ] = 0.5 * (x1 + x2);
    dedge_center[3*i+1] = 0.5 * (y1 + y2);
    dedge_center[3*i+2] = 0.5 * (z1 + z2);
  }

  free(edge_ln_to_gn);
  free(pvtx_ln_to_gn);
  free(pedge_vtx_idx);
  free(pedge_vtx);
  free(pvtx_coord);
  free(pfield);

  PDM_g_num_t *dentity_edge     = NULL;
  int         *dentity_edge_idx = NULL;
  PDM_g_num_t *distrib_entity   = NULL;
  int          dn_entity        = -1;

  if(isos->dim == 3) {
    PDM_deduce_combine_connectivity(isos->comm,
                                    isos->distrib_cell,
                                    isos->distrib_face,
                                    isos->dcell_face_idx,
                                    isos->dcell_face,
                                    isos->dface_edge_idx,
                                    isos->dface_edge,
                                    1,
                                    &dentity_edge_idx,
                                    &dentity_edge);
    distrib_entity = isos->distrib_cell;
  } else {
    dentity_edge     = isos->dface_edge;
    dentity_edge_idx = isos->dface_edge_idx;
    distrib_entity   = isos->distrib_face;
  }

  dn_entity = distrib_entity[i_rank+1] - distrib_entity[i_rank];

  /*
   *  Deduce all entity concerns by the iso surface (could be optimize)
   */
  PDM_block_to_part_t* btp = PDM_block_to_part_create(isos->distrib_edge,
                               (const PDM_g_num_t **) &dentity_edge,
                                                      &dentity_edge_idx[dn_entity],
                                                      1,
                                                      isos->comm);

  int strid_one = 1;
  int **tmp_dentity_edge_tag = NULL;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_tag,
            (int  ***)   NULL,
            (void ***)  &tmp_dentity_edge_tag);
  int *dentity_edge_tag = tmp_dentity_edge_tag[0];
  free(tmp_dentity_edge_tag);
  free(dedge_tag);

  double **tmp_dentity_edge_center = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         &strid_one,
            (void *  )   dedge_center,
            (int  ***)   NULL,
            (void ***)  &tmp_dentity_edge_center);
  double *dentity_edge_center = tmp_dentity_edge_center[0];
  free(tmp_dentity_edge_center);
  free(dedge_center);

  int         *dentity_tag            = malloc(     dn_entity * sizeof(int        ));
  PDM_g_num_t *entity_to_extract_gnum = malloc(     dn_entity * sizeof(PDM_g_num_t));
  double      *dentity_center         = malloc( 3 * dn_entity * sizeof(double     ));
  int  n_entity_tag = 0;
  int idx_write   = 0;
  for(int i = 0; i < dn_entity; ++i) {
    dentity_tag[i] = 0;

    for(int idx_entity = dentity_edge_idx[i]; idx_entity < dentity_edge_idx[i+1]; ++idx_entity) {
      if(dentity_edge_tag[idx_entity] == 1) {
        dentity_tag[i] = 1;
        entity_to_extract_gnum[n_entity_tag++] = distrib_entity[i_rank] + i + 1;
        break;
      }
    }

    if(dentity_tag[i] == 1) {
      dentity_center[3*idx_write  ] = 0.;
      dentity_center[3*idx_write+1] = 0.;
      dentity_center[3*idx_write+2] = 0.;

      double inv = 1./((double) (dentity_edge_idx[i+1] - dentity_edge_idx[i]));
      for(int idx_entity = dentity_edge_idx[i]; idx_entity < dentity_edge_idx[i+1]; ++idx_entity) {
        dentity_center[3*idx_write  ] += dentity_edge_center[3*idx_entity  ];
        dentity_center[3*idx_write+1] += dentity_edge_center[3*idx_entity+1];
        dentity_center[3*idx_write+2] += dentity_edge_center[3*idx_entity+2];
      }
      dentity_center[3*idx_write  ] = dentity_center[3*idx_write  ] * inv;
      dentity_center[3*idx_write+1] = dentity_center[3*idx_write+1] * inv;
      dentity_center[3*idx_write+2] = dentity_center[3*idx_write+2] * inv;

      idx_write++;
    }
  }
  free(dentity_edge_center);
  PDM_block_to_part_free(btp);

  free(dentity_edge_tag);
  free(dentity_edge_center);

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "out_iso_surf_equi_entity_coord_%2.2d.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_entity_tag,
                              dentity_center,
                              NULL,
                              NULL);
  }

  /*
   * Rebuild partition that contains entity and reequilibrate
   */
  PDM_gen_gnum_t* gnum_equi = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_coords(gnum_equi, 0, n_entity_tag, dentity_center, NULL);
  PDM_gnum_compute(gnum_equi);
  PDM_g_num_t* child_equi_entity_gnum = PDM_gnum_get(gnum_equi, 0);
  PDM_gnum_free(gnum_equi);
  free(dentity_center);

  /*
   * Equilibrage avec le part_to_block
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                       &child_equi_entity_gnum,
                                                       NULL,
                                                       &n_entity_tag,
                                                       1,
                                                       isos->comm);

  int n_entity_equi = PDM_part_to_block_n_elt_block_get (ptb);
  PDM_g_num_t *block_entity_equi_shild_g_num = PDM_part_to_block_block_gnum_get (ptb);

  PDM_g_num_t *block_entity_equi_parent_g_num = NULL;
  PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST_INTERLACED,
                          1,
                          NULL,
               (void **) &entity_to_extract_gnum,
                          NULL,
               (void **) &block_entity_equi_parent_g_num);

  free(entity_to_extract_gnum);

  /*
   * A refléchir pour le 3D
   */
  int          pn_edge_equi               = 0;
  PDM_g_num_t *pequi_parent_edge_ln_to_gn = NULL;
  int         *pequi_entity_edge_idx      = NULL;
  int         *pequi_entity_edge          = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           distrib_entity,
                                                           dentity_edge_idx,
                                                           dentity_edge,
                                                           n_entity_equi,
                                     (const PDM_g_num_t *) block_entity_equi_parent_g_num,
                                                           &pn_edge_equi,
                                                           &pequi_parent_edge_ln_to_gn,
                                                           &pequi_entity_edge_idx,
                                                           &pequi_entity_edge);

  PDM_gen_gnum_t* gnum_edge = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_edge, 0, pn_edge_equi, pequi_parent_edge_ln_to_gn);
  PDM_gnum_compute(gnum_edge);
  PDM_g_num_t* pequi_edge_ln_to_gn = PDM_gnum_get(gnum_edge, 0);
  PDM_gnum_free(gnum_edge);

  int          pn_vtx_equi               = 0;
  PDM_g_num_t *pequi_parent_vtx_ln_to_gn = NULL;
  int         *pequi_edge_vtx_idx        = NULL;
  int         *pequi_edge_vtx            = NULL;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(isos->comm,
                                                           isos->distrib_edge,
                                                           dedge_vtx_idx,
                                                           isos->dedge_vtx,
                                                           pn_edge_equi,
                                     (const PDM_g_num_t *) pequi_parent_edge_ln_to_gn,
                                                           &pn_vtx_equi,
                                                           &pequi_parent_vtx_ln_to_gn,
                                                           &pequi_edge_vtx_idx,
                                                           &pequi_edge_vtx);
  free(dedge_vtx_idx);

  PDM_gen_gnum_t* gnum_vtx = PDM_gnum_create(3, 1, PDM_FALSE, 0., isos->comm, PDM_OWNERSHIP_USER);
  PDM_gnum_set_from_parents(gnum_vtx, 0, pn_vtx_equi, pequi_parent_vtx_ln_to_gn);
  PDM_gnum_compute(gnum_vtx);
  PDM_g_num_t* pequi_vtx_ln_to_gn = PDM_gnum_get(gnum_vtx, 0);
  PDM_gnum_free(gnum_vtx);

  double **tmp_pequi_vtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(isos->comm,
                                        1,
                                        isos->distrib_vtx,
                                        isos->dvtx_coord,
                                        &pn_vtx_equi,
                 (const PDM_g_num_t **) &pequi_parent_vtx_ln_to_gn,
                                        &tmp_pequi_vtx_coord);
  double* pequi_vtx_coord = tmp_pequi_vtx_coord[0];
  free(tmp_pequi_vtx_coord);

  if(isos->dim == 2) {
    _iso_line_dist(isos);
  } else {
    _iso_surf_dist(isos);
  }


  free(pequi_edge_ln_to_gn);
  free(pequi_vtx_ln_to_gn);
  free(pequi_vtx_coord);
  free(pequi_parent_vtx_ln_to_gn);
  free(pequi_edge_vtx_idx);
  free(pequi_edge_vtx);
  free(pequi_parent_edge_ln_to_gn);
  free(pequi_entity_edge_idx);
  free(pequi_entity_edge);
  free(block_entity_equi_parent_g_num);
  PDM_part_to_block_free(ptb);


  if(isos->dim == 3) {
    free(dentity_edge);
    free(dentity_edge_idx);
  }


  free(dedge_center);
  free(dedge_tag);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_iso_surface_t*
PDM_iso_surface_create
(
 const int             dim,
 const int             n_part,
       PDM_ownership_t ownership,
       PDM_MPI_Comm    comm
)
{
  PDM_iso_surface_t *isos = (PDM_iso_surface_t *) malloc(sizeof(PDM_iso_surface_t));

  isos->dim       = dim;
  isos->n_part    = n_part;
  isos->ownership = ownership;
  isos->comm      = comm;
  isos->is_dist   = -1;

  isos->n_cell         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_face         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_edge         = (int          *) malloc(n_part * sizeof(int          ));
  isos->n_vtx          = (int          *) malloc(n_part * sizeof(int          ));

  isos->pcell_face     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pcell_face_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge     = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_edge_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pedge_vtx      = (int         **) malloc(n_part * sizeof(int         *));
  isos->cell_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->face_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->edge_ln_to_gn  = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));
  isos->vtx_ln_to_gn   = (PDM_g_num_t **) malloc(n_part * sizeof(PDM_g_num_t *));

  isos->pface_vtx_idx = (int         **) malloc(n_part * sizeof(int         *));
  isos->pface_vtx     = (int         **) malloc(n_part * sizeof(int         *));

  isos->pvtx_coord      = (double **) malloc(n_part * sizeof(double *));
  isos->pfield          = (double **) malloc(n_part * sizeof(double *));
  isos->pgradient_field = (double **) malloc(n_part * sizeof(double *));

  return isos;
}

void
PDM_iso_surface_compute
(
  PDM_iso_surface_t        *isos
)
{
  if(isos->is_dist == 0) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_iso_surface_compute Not implemented with partition layout\n");
  } else {
    _iso_surface_dist(isos);
  }



}

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
)
{
  isos->n_cell        [i_part] = n_cell;
  isos->n_face        [i_part] = n_face;
  isos->n_edge        [i_part] = n_edge;
  isos->n_vtx         [i_part] = n_vtx;
  isos->pcell_face    [i_part] = cell_face;
  isos->pcell_face_idx[i_part] = cell_face_idx;
  isos->pface_edge    [i_part] = face_edge;
  isos->pface_edge_idx[i_part] = face_edge_idx;
  isos->pedge_vtx     [i_part] = edge_vtx;
  isos->cell_ln_to_gn [i_part] = cell_ln_to_gn;
  isos->face_ln_to_gn [i_part] = face_ln_to_gn;
  isos->edge_ln_to_gn [i_part] = edge_ln_to_gn;
  isos->vtx_ln_to_gn  [i_part] = vtx_ln_to_gn;
  isos->pface_vtx_idx [i_part] = face_vtx_idx;
  isos->pface_vtx     [i_part] = face_vtx;
  isos->pvtx_coord    [i_part] = vtx_coord;
}


void
PDM_iso_surface_part_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *field
)
{
  isos->is_dist = 0;
  isos->pfield[i_part] = field;
}

void
PDM_iso_surface_part_gradient_field_set
(
  PDM_iso_surface_t        *isos,
  int                       i_part,
  double                   *gradient_field
)
{
  isos->is_dist = 0;
  isos->pgradient_field[i_part] = gradient_field;
}

void
PDM_iso_surface_dconnectivity_set
(
  PDM_iso_surface_t        *isos,
  PDM_connectivity_type_t   connectivity_type,
  PDM_g_num_t              *dconnect,
  int                      *dconnect_idx
)
{
  isos->is_dist = 1;
  switch (connectivity_type) {
   case PDM_CONNECTIVITY_TYPE_CELL_FACE:
     isos->dcell_face     = dconnect;
     isos->dcell_face_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_EDGE:
     isos->dface_edge     = dconnect;
     isos->dface_edge_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_FACE_VTX:
     isos->dface_vtx     = dconnect;
     isos->dface_vtx_idx = dconnect_idx;
     break;
   case PDM_CONNECTIVITY_TYPE_EDGE_VTX:
     isos->dedge_vtx     = dconnect;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid connectivity_type for iso_surface %d\n", connectivity_type);
    break;
   }

}

void
PDM_iso_surface_distrib_set
(
  PDM_iso_surface_t        *isos,
  PDM_mesh_entities_t       entity_type,
  PDM_g_num_t              *distrib_entity
)
{
  isos->is_dist = 1;
  switch (entity_type) {
   case PDM_MESH_ENTITY_CELL:
     isos->distrib_cell = distrib_entity;
     break;
   case PDM_MESH_ENTITY_FACE:
     isos->distrib_face = distrib_entity;
     break;
   case PDM_MESH_ENTITY_EDGE:
     isos->distrib_edge = distrib_entity;
     break;
   case PDM_MESH_ENTITY_VERTEX:
     isos->distrib_vtx = distrib_entity;
     break;
   default:
    PDM_error(__FILE__, __LINE__, 0, "invalid entity_type for iso_surface %d\n", entity_type);
    break;
   }
}


void
PDM_iso_surface_dvtx_coord_set
(
  PDM_iso_surface_t *isos,
  double            *dvtx_coord
)
{
  isos->is_dist     = 1;
  isos->dvtx_coord  = dvtx_coord;
}

void
PDM_iso_surface_dfield_set
(
  PDM_iso_surface_t *isos,
  double            *dfield
)
{
  isos->is_dist = 1;
  isos->dfield  = dfield;
}

void
PDM_iso_surface_dgrad_field_set
(
  PDM_iso_surface_t *isos,
  double            *dgrad_field
)
{
  isos->is_dist         = 1;
  isos->dgradient_field = dgrad_field;
}

void
PDM_iso_surface_free
(
  PDM_iso_surface_t        *isos
)
{
  free(isos->n_cell        );
  free(isos->n_face        );
  free(isos->n_edge        );
  free(isos->n_vtx         );
  free(isos->pcell_face    );
  free(isos->pcell_face_idx);
  // Si pface_edge a été calculé il faut le free
  free(isos->pface_edge    );
  free(isos->pface_edge_idx);

  free(isos->pface_vtx     );
  free(isos->pface_vtx_idx );
  free(isos->pedge_vtx     );
  free(isos->cell_ln_to_gn );
  free(isos->face_ln_to_gn );
  free(isos->edge_ln_to_gn );
  free(isos->vtx_ln_to_gn  );

  free(isos->pvtx_coord     );
  free(isos->pfield         );
  free(isos->pgradient_field);

  free(isos);
}
