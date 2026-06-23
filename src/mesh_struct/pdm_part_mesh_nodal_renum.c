/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_array.h"
#include "pdm_cuthill.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_order.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts_priv.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_renum.h"
#include "pdm_part_renum.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Macro definition
 *============================================================================*/



void
PDM_part_mesh_nodal_elmts_renum_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old_order
)
{
  // CHECK_PMNE(pmne)

  int  *n_entity2       = NULL;
  int **entity2_vtx     = NULL;
  int **entity2_vtx_idx = NULL;
  PDM_malloc(n_entity2      , pmne->n_part, int  );
  PDM_malloc(entity2_vtx    , pmne->n_part, int *);
  PDM_malloc(entity2_vtx_idx, pmne->n_part, int *);

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
    n_entity2[i_part] = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne, i_part, &entity2_vtx_idx[i_part], &entity2_vtx[i_part]);
  }



  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
    PDM_free(entity2_vtx    [i_part]);
    PDM_free(entity2_vtx_idx[i_part]);
  }
  PDM_free(entity2_vtx    );
  PDM_free(entity2_vtx_idx);
  PDM_free(n_entity2);


}


void
PDM_part_mesh_nodal_elmts_reorder_sections_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                           i_part,
  int                         **new_to_old_order_elmt_by_section
)
{
  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];

    PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                            id_section);

    int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                            id_section,
                                                            i_part);


    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      PDM_error(__FILE__, __LINE__, 0,
              "Error in PDM_part_mesh_nodal_elmts_reorder_sections_elmts : POLY_2D is not implemented \n");
    } else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
      PDM_error(__FILE__, __LINE__, 0,
              "Error in PDM_part_mesh_nodal_elmts_reorder_sections_elmts : POLY_3D is not implemented \n");
    } else {

      int         *connec;
      PDM_g_num_t *numabs;
      int         *_parent_num;
      PDM_g_num_t *parent_entity_g_num;
      int          order;
      const char  *ho_ordering;

      PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                   id_section,
                                                   i_part,
                                                   &connec,
                                                   &numabs,
                                                   &_parent_num,
                                                   &parent_entity_g_num,
                                                   &order,
                                                   &ho_ordering,
                                                   PDM_OWNERSHIP_BAD_VALUE);
      int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

      // PDM_log_trace_array_int(new_to_old_order_elmt_by_section[i_section], n_elt, "new_to_old_order_elmt_by_section[i_section]");
      // PDM_log_trace_array_int(connec, n_elt * n_vtx_elt, "connect avant ::");

      PDM_order_array(n_elt, sizeof(int) * n_vtx_elt, new_to_old_order_elmt_by_section[i_section], connec);
      // PDM_log_trace_array_int(connec, n_elt * n_vtx_elt, "connect après ::");
      if(numabs != NULL) {
        PDM_order_array(n_elt, sizeof(PDM_g_num_t), new_to_old_order_elmt_by_section[i_section], numabs);
      }
      if(_parent_num != NULL) {
        PDM_order_array(n_elt, sizeof(int), new_to_old_order_elmt_by_section[i_section], _parent_num);
      }
      if(parent_entity_g_num != NULL) {
        PDM_order_array(n_elt, sizeof(PDM_g_num_t), new_to_old_order_elmt_by_section[i_section], parent_entity_g_num);
      }

    }
  }

}

void
PDM_part_mesh_nodal_elmts_reorder_elmts
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  int                         **new_to_old_order
)
{
  int n_section = pmne->n_section;

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

    int *section_id = NULL;
    PDM_malloc(section_id, pmne->n_elmts[i_part], int);
    int idx_elmt = 0;
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, i_section, i_part);
      for(int i_elt = 0; i_elt < n_elt_in_section; ++i_elt) {
        section_id[idx_elmt++] = i_section;
      }
    }

    PDM_order_array(pmne->n_elmts[i_part],
                    sizeof(int),
                    new_to_old_order[i_part],
                    section_id);

    int **new_to_old_order_elmt_by_section = NULL;
    PDM_malloc(new_to_old_order_elmt_by_section, n_section, int *);

    int *idx_write_section = PDM_array_zeros_int(n_section);
    int *section_elmt_idx  = PDM_array_zeros_int(n_section+1);
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt_in_section = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne, i_section, i_part);
      section_elmt_idx[i_section+1] = section_elmt_idx[i_section] + n_elt_in_section;
      PDM_malloc(new_to_old_order_elmt_by_section[i_section], n_elt_in_section, int);
    }

    for(int i = 0; i < pmne->n_elmts[i_part]; ++i) {
      int i_section = section_id[i]; // Dans le nouvelle ordre
      new_to_old_order_elmt_by_section[i_section][idx_write_section[i_section]++] = new_to_old_order[i_part][i] - section_elmt_idx[i_section];
    }

    PDM_part_mesh_nodal_elmts_reorder_sections_elmts(pmne, i_part, new_to_old_order_elmt_by_section);

    for (int i_section = 0; i_section < n_section; i_section++) {
      PDM_free(new_to_old_order_elmt_by_section[i_section]);
    }
    PDM_free(new_to_old_order_elmt_by_section);
    PDM_free(idx_write_section);
    PDM_free(section_elmt_idx);
    PDM_free(section_id);

  }


}



void
PDM_part_mesh_nodal_elmts_reorder_vtx
(
  PDM_part_mesh_nodal_elmts_t  *pmne,
  const int                   **old_to_new_order_vtx
)
{

  int  n_section   = PDM_part_mesh_nodal_elmts_n_section_get  (pmne);
  int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
    for(int i_section = 0; i_section < n_section; ++i_section) {

      int id_section = sections_id[i_section];

      PDM_Mesh_nodal_elt_t t_elt = PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                              id_section);

      int n_elt = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                              id_section,
                                                              i_part);


      if (t_elt == PDM_MESH_NODAL_POLY_2D) {
        PDM_error(__FILE__, __LINE__, 0,
                "Error in PDM_part_mesh_nodal_elmts_reorder_vtx : POLY_2D is not implemented \n");
      } else if (t_elt == PDM_MESH_NODAL_POLY_3D) {
        PDM_error(__FILE__, __LINE__, 0,
                "Error in PDM_part_mesh_nodal_elmts_reorder_vtx : POLY_3D is not implemented \n");
      } else {

        int         *connec;
        PDM_g_num_t *numabs;
        int         *_parent_num;
        PDM_g_num_t *parent_entity_g_num;
        int          order;
        const char  *ho_ordering;

        PDM_part_mesh_nodal_elmts_section_std_ho_get(pmne,
                                                     id_section,
                                                     i_part,
                                                     &connec,
                                                     &numabs,
                                                     &_parent_num,
                                                     &parent_entity_g_num,
                                                     &order,
                                                     &ho_ordering,
                                                     PDM_OWNERSHIP_BAD_VALUE);
        int n_vtx_elt = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);

        PDM_part_renum_array(n_elt * n_vtx_elt, old_to_new_order_vtx[i_part], connec);

      }
    }
  }
}



void
PDM_part_mesh_nodal_reorder_vtx
(
  PDM_part_mesh_nodal_t  *pmn,
  int                   **old_to_new_order_vtx,
  int                   **new_to_old_order_vtx
)
{

  int dimension = PDM_part_mesh_nodal_mesh_dimension_get(pmn);

  for (int i_dim = 1; i_dim <= dimension; i_dim++) {
    PDM_mesh_entities_t entity_type = PDM_dimension_to_entity_type(i_dim);
    PDM_geometry_kind_t geom_kind   = PDM_entity_type_to_geometry_kind(entity_type);

    int n_group = PDM_part_mesh_nodal_n_group_get(pmn, geom_kind);

    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

    PDM_part_mesh_nodal_elmts_reorder_vtx(pmne, (const int **) old_to_new_order_vtx);

  }

  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    PDM_Mesh_nodal_vtx_t *vtx = pmn->vtx[i_part];

    PDM_order_array (vtx->n_vtx,
                     3 * sizeof(double),
                     new_to_old_order_vtx[i_part],
                     vtx->_coords);

    if(vtx->_numabs) {
      PDM_order_array (vtx->n_vtx,
                       sizeof(PDM_g_num_t),
                       new_to_old_order_vtx[i_part],
                       vtx->_numabs);
    }

    if(vtx->_numparent) {
      PDM_order_array (vtx->n_vtx,
                       sizeof(int),
                       new_to_old_order_vtx[i_part],
                       vtx->_numparent);
    }
  }

}


// 1 - Part_mesh
// 2 - Part_mesh_nodal (plusieurs possible -> pas géré - RCM / HYPERPLANE / Est-ce quu'on gère sommet + elemts en même temps ou decolerer)
// 3 - Part_mesh & Part_mesh_nodal :
//       -

void
PDM_part_mesh_nodal_renum
(
  PDM_part_mesh_nodal_t *pmn, // A voir si on donne le part_mesh en entré directement
  const char            *renum_entity_method    [4],
  const int             *renum_entity_properties[4]
)
{
  PDM_UNUSED(renum_entity_method);
  PDM_UNUSED(renum_entity_properties);

  int n_part = pmn->n_part;

  // int renum_method_cell = PDM_part_renum_method_cell_idx_get(renum_entity_method[3]);
  // int renum_method_face = PDM_part_renum_method_face_idx_get(renum_entity_method[2]);
  // int renum_method_edge = PDM_part_renum_method_edge_idx_get(renum_entity_method[1]);
  // int renum_method_vtx  = PDM_part_renum_method_vtx_idx_get (renum_entity_method[0]);

  // PDM_PART_RENUM_VTX_NONE
  // PDM_PART_RENUM_NODAL_VTX_NONE

  // Si dédié aux nodal => pas besoin du part_mesh
  if(0 == 1) {
    // Create part_mesh & il faudra le garder
    // Renum part_mesh
    // PDM_part_mesh_t* pm = NULL;
    // PDM_part_mesh_renum(pm,
    //                     renum_entity_method,
    //                     renum_entity_properties);
    abort();

    // Apply renum into part_mesh_nodal

  } else {

    int dimension = PDM_part_mesh_nodal_mesh_dimension_get(pmn);

    int  *n_vtx        = NULL;
    int  *n_elmt       = NULL;
    int **elmt_vtx     = NULL;
    int **elmt_vtx_idx = NULL;
    PDM_malloc(n_vtx       , pmn->n_part, int  );
    PDM_malloc(n_elmt      , pmn->n_part, int  );
    PDM_malloc(elmt_vtx    , pmn->n_part, int *);
    PDM_malloc(elmt_vtx_idx, pmn->n_part, int *);

    PDM_part_mesh_nodal_elmts_t* pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_principal_dim_get(pmn);

    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
      n_vtx [i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      n_elmt[i_part] = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne, i_part, &elmt_vtx_idx[i_part], &elmt_vtx[i_part]);
    }

    int **vtx_elmt_idx = NULL;
    int **vtx_elmt     = NULL;
    PDM_part_connectivity_transpose(n_part,
                                    n_elmt,
                                    n_vtx,
                                    elmt_vtx_idx,
                                    elmt_vtx,
                                    &vtx_elmt_idx,
                                    &vtx_elmt);

    int **elmt_elmt_idx = NULL;
    int **elmt_elmt     = NULL;
    PDM_part_combine_connectivity(n_part,
                                  n_elmt,
                                  elmt_vtx_idx,
                                  elmt_vtx,
                                  vtx_elmt_idx,
                                  vtx_elmt,
                                  &elmt_elmt_idx,
                                  &elmt_elmt);

    int n_section = PDM_part_mesh_nodal_n_section_get(pmn);

    // BFS
    int **new_to_old_order_elmt = NULL;
    PDM_malloc(new_to_old_order_elmt, pmne->n_part, int *);
    for(int i_part = 0; i_part < pmne->n_part; ++i_part) {

      for(int i = 0; i < elmt_elmt_idx[i_part][n_elmt[i_part]]; ++i) {
        elmt_elmt[i_part][i] -= 1;
      }
      // PDM_log_trace_connectivity_int(elmt_elmt_idx[i_part], elmt_elmt[i_part], n_elmt[i_part], "elmt_elmt :");

      PDM_malloc(new_to_old_order_elmt[i_part], n_elmt[i_part], int);
      // PDM_malloc(section_id           [i_part], n_elmt[i_part], int);
      PDM_cuthill_generate(n_elmt[i_part], elmt_elmt_idx[i_part], elmt_elmt[i_part], new_to_old_order_elmt[i_part]);

      // PDM_log_trace_array_int(new_to_old_order_elmt[i_part], n_elmt[i_part], "new_to_old_order_elmt ::");


    }

    // Deduction du new_to_old_order par section
    PDM_part_mesh_nodal_elmts_reorder_elmts(pmne, new_to_old_order_elmt);


    // PDM_part_mesh_nodal_elmts_reorder_elmts(pmne, elmts_new_to_old_order);

    // Step 2 :
    //   - Renum vertices dans l'ordre d'apparition des élements
    int **new_to_old_order_vtx = NULL;
    int **old_to_new_order_vtx = NULL;
    PDM_malloc(new_to_old_order_vtx, pmne->n_part, int *);
    PDM_malloc(old_to_new_order_vtx, pmne->n_part, int *);
    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {

      new_to_old_order_vtx[i_part] = PDM_array_const_int(n_vtx[i_part], -1);
      old_to_new_order_vtx[i_part] = PDM_array_const_int(n_vtx[i_part], -1);

      int n_vtx_new = 0;
      for(int i_elt = 0; i_elt < n_elmt[i_part]; ++i_elt) {
        int old_order = new_to_old_order_elmt[i_part][i_elt]; // Par section plutot ?
        for(int idx_vtx = elmt_vtx_idx[i_part][old_order]; idx_vtx < elmt_vtx_idx[i_part][old_order+1]; ++idx_vtx) {
          int old_vtx_id = elmt_vtx[i_part][idx_vtx]-1;

          if (old_to_new_order_vtx[i_part][old_vtx_id] == -1) {
            old_to_new_order_vtx[i_part][old_vtx_id] = n_vtx_new;
            new_to_old_order_vtx[i_part][n_vtx_new] = old_vtx_id;
            n_vtx_new++;
          }
        }
      }
    }

    PDM_part_mesh_nodal_reorder_vtx(pmn, old_to_new_order_vtx, new_to_old_order_vtx);

    for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
      PDM_free(new_to_old_order_vtx[i_part]);
      PDM_free(old_to_new_order_vtx[i_part]);
    }
    PDM_free(new_to_old_order_vtx);
    PDM_free(old_to_new_order_vtx);

    // Step 3 :
    //   - Renum all lower dimension (On prends le min des sommets connecté)

    for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
      PDM_free(elmt_vtx     [i_part]);
      PDM_free(elmt_vtx_idx [i_part]);
      PDM_free(vtx_elmt_idx [i_part]);
      PDM_free(vtx_elmt     [i_part]);
      PDM_free(elmt_elmt_idx[i_part]);
      PDM_free(elmt_elmt    [i_part]);
    }
    PDM_free(elmt_vtx    );
    PDM_free(elmt_vtx_idx);
    PDM_free(vtx_elmt_idx);
    PDM_free(vtx_elmt    );
    PDM_free(elmt_elmt_idx);
    PDM_free(elmt_elmt);
    PDM_free(n_elmt);
    PDM_free(n_vtx);


    for(int i_part = 0; i_part < pmne->n_part; ++i_part) {
      PDM_free(new_to_old_order_elmt[i_part]);
    }
    PDM_free(new_to_old_order_elmt);

    // for (int i_dim = 1; i_dim <= dimension; i_dim++) {
    //   PDM_mesh_entities_t entity_type = PDM_dimension_to_entity_type(i_dim);
    //   PDM_geometry_kind_t geom_kind   = PDM_entity_type_to_geometry_kind(entity_type);

    //   int n_group = PDM_part_mesh_nodal_n_group_get(pmn, geom_kind);

    //   PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);




    //   int **elmts_new_to_old_order = NULL;
    //   PDM_part_mesh_nodal_elmts_renum_elmts(pmne, &elmts_new_to_old_order);

    //   PDM_part_mesh_nodal_elmts_reorder_elmts(pmne, elmts_new_to_old_order);
    // }



  }


}

#ifdef __cplusplus
}
#endif /* __cplusplus */
