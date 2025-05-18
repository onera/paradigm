/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_geom.h"
#include "pdm_part_mesh_geom.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_error.h"
#include "pdm_predicate.h"
#include "pdm_mem_tool.h"
#include "pdm_logging.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
int
_check_if_all_simplices
(
  PDM_part_mesh_nodal_t  *pmn
)
{
  int all_simplices = 1;

  // Nodal
  PDM_geometry_kind_t  geom_kind    = (pmn->mesh_dimension == 2) ? PDM_GEOMETRY_KIND_SURFACIC : PDM_GEOMETRY_KIND_VOLUMIC;
  PDM_Mesh_nodal_elt_t simplex_type = (pmn->mesh_dimension == 2) ? PDM_MESH_NODAL_TRIA3       : PDM_MESH_NODAL_TETRA4;

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for (int i_section = 0; i_section < n_section; i_section++) {
    PDM_Mesh_nodal_elt_t elt_type = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                               sections_id[i_section]);

    if (elt_type != simplex_type) {
      all_simplices = 0;
      break;
    }
  }

  return all_simplices;
}

/**
 * \brief Compute *local* measure of dual volumes.
 */
static void
_compute_dual_volume_simplex
(
  int                 dim,
  int                 n_part,
  int                *n_elt,
  int                *n_vtx,
  int               **elt_vtx,
  double            **vtx_coord,
  double           ***out_vtx_volume
)
{
  double **vtx_volume = NULL;
  PDM_malloc(vtx_volume, n_part, double *);

  int elt_size = dim + 1;
  double factor = (dim == 2) ? 1./6. : 1./24.; // vol / (det * n_vtx_per_elt)

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_calloc(vtx_volume[i_part], n_vtx[i_part], double);

    for (int i_elt = 0; i_elt < n_elt[i_part]; i_elt++) {

      int *_elt_vtx = elt_vtx[i_part] + elt_size*i_elt;

      double elt_vol = 0;
      if (dim == 2) {
        elt_vol = factor * PDM_predicate_orient2d(vtx_coord[i_part] + 3*(_elt_vtx[0]-1),
                                                  vtx_coord[i_part] + 3*(_elt_vtx[1]-1),
                                                  vtx_coord[i_part] + 3*(_elt_vtx[2]-1));
      }
      else {
        elt_vol = factor * PDM_predicate_orient3d(vtx_coord[i_part] + 3*(_elt_vtx[0]-1),
                                                  vtx_coord[i_part] + 3*(_elt_vtx[1]-1),
                                                  vtx_coord[i_part] + 3*(_elt_vtx[2]-1),
                                                  vtx_coord[i_part] + 3*(_elt_vtx[3]-1));
      }

      for (int i = 0; i < elt_size; i++) {
        vtx_volume[i_part][_elt_vtx[i]-1] += elt_vol;
      }

    } // End loop on elements
  } // End loop on parts

  *out_vtx_volume = vtx_volume;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

void
PDM_part_mesh_nodal_dual_volume_compute
(
  PDM_part_mesh_nodal_t   *pmn,
  double                ***dual_vol
)
{
  /**
   *  Prevoir une syncrho en option ?
   *  OU
   *  Adpater la gradation pour utiliser habilement le pcg pour ne pas sommer 2 fois les contributions de complexité
   */

  int all_simplices = _check_if_all_simplices(pmn);

  if (all_simplices) {
    PDM_geometry_kind_t geom_kind = (pmn->mesh_dimension == 2) ? PDM_GEOMETRY_KIND_SURFACIC : PDM_GEOMETRY_KIND_VOLUMIC;

    int     *n_vtx     = NULL;
    int     *n_elt     = NULL;
    int    **elt_vtx   = NULL;
    double **vtx_coord = NULL;
    PDM_malloc(n_vtx    , pmn->n_part, int     );
    PDM_malloc(n_elt    , pmn->n_part, int     );
    PDM_malloc(elt_vtx  , pmn->n_part, int    *);
    PDM_malloc(vtx_coord, pmn->n_part, double *);

    for (int i_part = 0; i_part <pmn->n_part; i_part++) {
      int *elt_vtx_idx = NULL;
      n_elt[i_part] = PDM_part_mesh_nodal_cell_vtx_connect_get(pmn,
                                                               geom_kind,
                                                               i_part,
                                                               &elt_vtx_idx,
                                                               &elt_vtx[i_part]);
      n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get    (pmn, i_part);
      vtx_coord[i_part] = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

      PDM_free(elt_vtx_idx);
    }

    _compute_dual_volume_simplex(pmn->mesh_dimension,
                                 pmn->n_part,
                                 n_elt,
                                 n_vtx,
                                 elt_vtx,
                                 vtx_coord,
                                 dual_vol);

    for (int i_part = 0; i_part < pmn->n_part; i_part++) {
      PDM_free(elt_vtx[i_part]);
    }
    PDM_free(n_vtx);
    PDM_free(n_elt);
    PDM_free(elt_vtx);
    PDM_free(vtx_coord);

    if(pmn->pcg[PDM_MESH_ENTITY_VTX] == NULL) {
      PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_MESH_ENTITY_VTX);
    }

    // Synchro volume :
    PDM_part_comm_graph_all_reduce(pmn->pcg[PDM_MESH_ENTITY_VTX],
                                   PDM_MPI_DOUBLE,
                                   PDM_MPI_SUM,
             (unsigned char **)    *dual_vol);


  } else {

    /*
     * Dual volume computation need all downing connectivity
     *   - We use part_mesh_nodal_to_part_mesh to express all connectity then compute the dual volume
     */
    PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                            PDM_FALSE,
                                                                                            PDM_OWNERSHIP_USER);

    if (pmn->mesh_dimension == 3) {
      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_CELL_FACE);
    }
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_FACE_EDGE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm, PDM_CONNECTIVITY_TYPE_EDGE_VTX);
    PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_VTX);

    PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

    PDM_part_mesh_t *pmesh = NULL;
    PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm, &pmesh, PDM_OWNERSHIP_KEEP);

    // Compute dual volume via part_mesh
    PDM_part_mesh_dual_volume_compute(pmesh, dual_vol);

    PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);
  }

}


#ifdef __cplusplus
}
#endif /* __cplusplus */
