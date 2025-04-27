/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_array.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_unique.h"
#include "pdm_part_comm_graph.h"
#include "pdm_priv.h"
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


static void
_find_local_parent_group
(
  PDM_part_mesh_nodal_elmts_t *pmne,
  int                          n_part,
  int                        **entity1_entity2_idx,
  int                        **entity1_entity2,
  int                         *n_entity2,
  int                       ***entity2_tag_n,
  int                       ***entity2_tag_idx,
  int                       ***entity2_tag
)
{
  int **_entity2_tag_n   = NULL;
  int **_entity2_tag_idx = NULL;
  int **_entity2_tag     = NULL;
  PDM_malloc(_entity2_tag_n  , n_part, int *);
  PDM_malloc(_entity2_tag_idx, n_part, int *);
  PDM_malloc(_entity2_tag    , n_part, int *);

  int n_group_entity1 = PDM_part_mesh_nodal_elmts_n_group_get(pmne);

  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_calloc(_entity2_tag_n[i_part], n_entity2[i_part], int);

    for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
      int          n_entity1_group        = 0;
      int         *entity1_group          = NULL;
      PDM_g_num_t *entity1_group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                         &n_entity1_group,
                                         &entity1_group,
                                         &entity1_group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);

      for(int i_entity1 = 0; i_entity1 < n_entity1_group; ++i_entity1) {
        int entity1 = entity1_group[i_entity1]-1;
        for (int i_entity2=entity1_entity2_idx[i_part][entity1  ];
                 i_entity2<entity1_entity2_idx[i_part][entity1+1]; ++i_entity2) {
          int entity2 = PDM_ABS(entity1_entity2[i_part][i_entity2])-1;
          _entity2_tag_n[i_part][entity2]++;
        }
      }
    }

    PDM_malloc(_entity2_tag_idx[i_part], n_entity2[i_part]+1, int);
    _entity2_tag_idx[i_part][0] = 0;
    for(int i_vtx = 0; i_vtx < n_entity2[i_part]; ++i_vtx) {
      _entity2_tag_idx[i_part][i_vtx+1] = _entity2_tag_idx[i_part][i_vtx] + _entity2_tag_n[i_part][i_vtx];
      _entity2_tag_n  [i_part][i_vtx  ] = 0;
    }

    PDM_malloc(_entity2_tag[i_part], _entity2_tag_idx[i_part][n_entity2[i_part]], int);
    for(int i_group = 0; i_group < n_group_entity1; ++i_group) {
      int          n_entity1_group        = 0;
      int         *entity1_group          = NULL;
      PDM_g_num_t *entity1_group_ln_to_gn = NULL;
      PDM_part_mesh_nodal_elmts_group_get(pmne,
                                          i_part,
                                          i_group,
                                         &n_entity1_group,
                                         &entity1_group,
                                         &entity1_group_ln_to_gn,
                                          PDM_OWNERSHIP_BAD_VALUE);

      for(int i_entity1 = 0; i_entity1 < n_entity1_group; ++i_entity1) {
        int entity1 = entity1_group[i_entity1]-1;
        for (int i_entity2=entity1_entity2_idx[i_part][entity1  ];
                 i_entity2<entity1_entity2_idx[i_part][entity1+1]; ++i_entity2) {
          int entity2 = PDM_ABS(entity1_entity2[i_part][i_entity2])-1;
          int i_write = _entity2_tag_idx[i_part][entity2] + _entity2_tag_n[i_part][entity2]++;
          _entity2_tag[i_part][i_write] = i_group+1;
        }
      }
    }

  }

  *entity2_tag_n   = _entity2_tag_n;
  *entity2_tag_idx = _entity2_tag_idx;
  *entity2_tag     = _entity2_tag;

}


static void
_transpose_group_information
(
  PDM_MPI_Comm    comm,
  int             n_part,
  int            *n_entity,
  int           **entity_group_idx,
  PDM_g_num_t   **entity_group,
  int            *out_n_group,
  int          ***out_group_entity_idx,
  int          ***out_group_entity
)
{
  int debug_verbose = 0;

  int **_entity_group     = NULL;
  int **_group_entity_idx = NULL;
  int **_group_entity     = NULL;
  PDM_malloc(_entity_group    , n_part, int *);
  PDM_malloc(_group_entity_idx, n_part, int *);
  PDM_malloc(_group_entity    , n_part, int *);

  int ln_group = 0;
  int gn_group = 0;
  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_malloc(_entity_group[i_part], entity_group_idx[i_part][n_entity[i_part]], int);

    for (int i_group=0; i_group<entity_group_idx[i_part][n_entity[i_part]]; ++i_group) {
      ln_group = PDM_MAX(entity_group[i_part][i_group], ln_group);
      _entity_group[i_part][i_group] = (int) entity_group[i_part][i_group];
    }

    if (debug_verbose==1) {
      int entity_group_size = entity_group_idx[i_part][n_entity[i_part]];
      PDM_log_trace_array_int(entity_group_idx[i_part], n_entity[i_part]+1, " entity_group_idx :: ");
      PDM_log_trace_array_int(_entity_group   [i_part], entity_group_size , "_entity_group     :: ");
    }
  }

  PDM_MPI_Allreduce(&ln_group, &gn_group, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  for (int i_part=0; i_part<n_part; ++i_part) {
    PDM_connectivity_transpose(n_entity[i_part],
                               gn_group,
                               entity_group_idx [i_part],
                               _entity_group    [i_part],
                              &_group_entity_idx[i_part],
                              &_group_entity    [i_part]);

    if (debug_verbose==1) {
      int group_entity_size = _group_entity_idx[i_part][gn_group];
      PDM_log_trace_array_int(_group_entity_idx[i_part], gn_group+1       , "group_entity_idx :: ");
      PDM_log_trace_array_int(_group_entity    [i_part], group_entity_size, "group_entity     :: ");
    }

    PDM_free(_entity_group[i_part]);
  }
  PDM_free(_entity_group);

  *out_n_group = gn_group;
  *out_group_entity_idx = _group_entity_idx;
  *out_group_entity     = _group_entity;
}


static void
_generate_group_entity_gnum
(
  PDM_MPI_Comm    comm,
  int             n_part,
  PDM_g_num_t   **entity_gnum,
  int             n_group,
  int           **group_entity_idx,
  int           **group_entity,
  PDM_g_num_t  ***out_group_entity_gnum
)
{
  PDM_g_num_t **group_entity_gnum = NULL;
  PDM_malloc(group_entity_gnum, n_part, PDM_g_num_t *);

  for (int i_part=0; i_part<n_part; ++i_part) {

    PDM_malloc(group_entity_gnum[i_part], group_entity_idx[i_part][n_group], PDM_g_num_t);

    for (int i_entity=0; i_entity<group_entity_idx[i_part][n_group]; ++i_entity) {
      group_entity_gnum[i_part][i_entity] = entity_gnum[i_part][group_entity[i_part][i_entity]-1];
    }
  }

  for (int i_group=0; i_group<n_group; ++i_group) {

    PDM_gen_gnum_t* gen_group_gnum = PDM_gnum_create(3,
                                                     n_part,
                                                     PDM_TRUE,
                                                     1e-4,
                                                     comm,
                                                     PDM_OWNERSHIP_USER);
    PDM_gnum_set_parents_nuplet(gen_group_gnum, 1);
    for (int i_part=0; i_part<n_part; ++i_part) {
      int n_group_entity = group_entity_idx[i_part][i_group+1]-group_entity_idx[i_part][i_group];
      PDM_gnum_set_from_parents(gen_group_gnum, i_part, n_group_entity, &group_entity_gnum[i_part][group_entity_idx[i_part][i_group]]);
    }

    PDM_gnum_compute(gen_group_gnum);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int n_group_entity = group_entity_idx[i_part][i_group+1]-group_entity_idx[i_part][i_group];
      PDM_g_num_t *group_gnum = PDM_gnum_get(gen_group_gnum, i_part);
      for (int i_entity=0; i_entity<n_group_entity; ++i_entity) {
        group_entity_gnum[i_part][group_entity_idx[i_part][i_group]+i_entity] = group_gnum[i_entity];
      }
      PDM_free(group_gnum);
    }

    PDM_gnum_free(gen_group_gnum);
  }

  *out_group_entity_gnum = group_entity_gnum;
}


static void
_generate_group_gnum
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int           *n_entity,
  int          **entity_tag_idx,
  int          **entity_tag,
  int         ***out_entity_group_idx,
  PDM_g_num_t ***out_entity_group
)
{

  int           l_n_parent_max          = 0;
  int         **group_parent_n          = NULL;
  int         **entity_group_idx        = NULL;
  // int         **entity_parent_group_idx = NULL;
  // int         **entity_parent_group     = NULL;
  PDM_g_num_t **group_parent_nplt       = NULL;
  PDM_malloc(group_parent_n         , n_part, int         *);
  PDM_malloc(entity_group_idx       , n_part, int         *);
  // PDM_malloc(entity_parent_group_idx, n_part, int         *);
  // PDM_malloc(entity_parent_group    , n_part, int         *);
  PDM_malloc(group_parent_nplt      , n_part, PDM_g_num_t *);

  for (int i_part=0; i_part<n_part; ++i_part) {
    /**
     * All parent group are known by entities and coherent over all procs, unique parents and :
     *  - count l_n_parent_max for gnum gnum
     *  - count n_parent for each entity
     *  - count straddling entity in entity_group_idx
     */
    // int entity_parent_group_size = 0;

    PDM_malloc(group_parent_n  [i_part], n_entity[i_part]  , int);
    PDM_malloc(entity_group_idx[i_part], n_entity[i_part]+1, int); entity_group_idx[i_part][0] = 0;
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      group_parent_n[i_part][i_entity  ] = 0;
      entity_group_idx [i_part][i_entity+1] = entity_group_idx[i_part][i_entity];
      int beg = entity_tag_idx[i_part][i_entity  ];
      int end = entity_tag_idx[i_part][i_entity+1];
      if(end - beg > 0){
        int n_unique   = PDM_inplace_unique(entity_tag[i_part], beg, end-1);
        l_n_parent_max = PDM_MAX(l_n_parent_max, n_unique);
        if(n_unique > 1) {
          // entity_parent_group_size  += n_unique;
          group_parent_n  [i_part][i_entity  ] = n_unique;
          entity_group_idx[i_part][i_entity+1] = entity_group_idx[i_part][i_entity]+1;
        }
      }
    }
  }
  int n_parent_max = -1;
  PDM_MPI_Allreduce(&l_n_parent_max, &n_parent_max, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);

  /**
   * Prepare nuplet from group gnum computation
   *
   * Here we could store parent group at same time
   */
  for (int i_part=0; i_part<n_part; ++i_part) {
    // int i_write_parent = 0;
    int i_write_nuplet = 0;
    // PDM_malloc(entity_parent_group_idx[i_part], n_entity[i_part]+1                                     , int        ); entity_parent_group_idx[i_part][0] = 0;
    // PDM_malloc(entity_parent_group    [i_part], entity_parent_group_size                               , int        );
    PDM_malloc(group_parent_nplt      [i_part], n_parent_max*entity_group_idx[i_part][n_entity[i_part]], PDM_g_num_t);
    for(int i_entity = 0; i_entity < n_entity[i_part]; ++i_entity) {
      // entity_parent_group_idx[i_part][i_entity+1] = entity_parent_group_idx[i_part][i_entity];
      int beg = entity_tag_idx[i_part][i_entity  ];
      if (group_parent_n[i_part][i_entity]>0) {
        // entity_parent_group_idx[i_part][i_entity+1] += group_parent_n[i_part][i_entity];
        for (int i_read=beg; i_read<beg+group_parent_n[i_part][i_entity]; ++i_read) {
          group_parent_nplt[i_part][i_write_nuplet++] = (PDM_g_num_t) entity_tag[i_part][i_read];
          // entity_parent_group [i_part][i_write_parent++] = entity_tag[i_part][i_read];
        }
        for (int i_read=beg+group_parent_n[i_part][i_entity]; i_read<beg+n_parent_max; ++i_read) {
          group_parent_nplt[i_part][i_write_nuplet++] = 0;
        }
      }
    }
    PDM_free(group_parent_n[i_part]);
  }
  PDM_free(group_parent_n);


  /**
   * Generate ids for groups
   */
  PDM_gen_gnum_t* gen_group_id = PDM_gnum_create(3,
                                                 n_part,
                                                 PDM_TRUE,
                                                 1e-4,
                                                 comm,
                                                 PDM_OWNERSHIP_USER);
  PDM_gnum_set_parents_nuplet(gen_group_id, n_parent_max);
  for (int i_part=0; i_part<n_part; ++i_part) {
    PDM_gnum_set_from_parents(gen_group_id, i_part, entity_group_idx[i_part][n_entity[i_part]], group_parent_nplt[i_part]);
  }
  PDM_gnum_compute(gen_group_id);

  PDM_g_num_t **group_id  = NULL;
  // int         **entity_group = NULL;
  PDM_malloc(group_id , n_part, PDM_g_num_t *);
  // PDM_malloc(entity_group, n_part, int         *);
  for (int i_part=0; i_part<n_part; ++i_part) {
    group_id[i_part] = PDM_gnum_get(gen_group_id, i_part);
    // PDM_malloc(entity_group[i_part], entity_group_idx[i_part][n_entity[i_part]], int);
    PDM_free(group_parent_nplt[i_part]);

  }

  PDM_gnum_free(gen_group_id);
  PDM_free(group_parent_nplt);

  *out_entity_group_idx = entity_group_idx;
  *out_entity_group     = group_id;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


void
PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_mesh_entities_t    entity_type
)
{
  if (pmn->pcg[entity_type] != NULL) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum: pmn->pcg[entity_type=%d]!=NULL\n", entity_type);
  }

  int          *n_entity    = NULL;
  PDM_g_num_t **entity_gnum = NULL;
  PDM_malloc(n_entity   , pmn->n_part, int         );
  PDM_malloc(entity_gnum, pmn->n_part, PDM_g_num_t*);

  if (entity_type==PDM_MESH_ENTITY_VTX) {
    for (int i_part = 0; i_part < pmn->n_part; ++i_part) {
      n_entity   [i_part] = PDM_part_mesh_nodal_n_vtx_get    (pmn, i_part);
      entity_gnum[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    }
  } else if (entity_type==PDM_MESH_ENTITY_EDGE || entity_type==PDM_MESH_ENTITY_FACE) {

    PDM_geometry_kind_t geom_kind = PDM_entity_type_to_geometry_kind(entity_type);
    PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn, geom_kind);

    for (int i_part = 0; i_part < pmn->n_part; ++i_part) {
      n_entity   [i_part] = PDM_part_mesh_nodal_elmts_n_elmts_get(pmne, i_part);
      entity_gnum[i_part] = PDM_part_mesh_nodal_elmts_g_num_get_from_part(pmne, i_part, PDM_OWNERSHIP_KEEP);
    }
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum: invalid entity_type (=%d)\n", entity_type);
  }

  // Retrieve partition boundary entities using global IDs
  int n_rank;
  PDM_MPI_Comm_size(pmn->comm, &n_rank);

  PDM_g_num_t *part_distribution = PDM_compute_entity_distribution(pmn->comm, pmn->n_part);

  int  *n_entity_part_bound        = NULL;
  int **entity_proc_bound_idx      = NULL;
  int **entity_part_bound_idx      = NULL;
  int **entity_part_bound          = NULL;
  int **entity_part_bound_priority = NULL;
  PDM_part_generate_entity_graph_comm(pmn->comm,
                                      part_distribution,
                                      NULL,
                                      pmn->n_part,
                                      n_entity,
               (const PDM_g_num_t **) entity_gnum,
                                      NULL,
                                      &entity_proc_bound_idx,
                                      &entity_part_bound_idx,
                                      &entity_part_bound,
                                      &entity_part_bound_priority);
  PDM_malloc(n_entity_part_bound, pmn->n_part, int);
  for (int i_part = 0; i_part < pmn->n_part; i_part++) {
    n_entity_part_bound[i_part] = entity_part_bound_idx[i_part][part_distribution[n_rank]];

    PDM_free(entity_proc_bound_idx     [i_part]);
    PDM_free(entity_part_bound_idx     [i_part]);
    PDM_free(entity_part_bound_priority[i_part]);
  }
  PDM_free(entity_proc_bound_idx);
  PDM_free(entity_part_bound_idx);
  PDM_free(entity_part_bound_priority);
  PDM_free(part_distribution);

  PDM_free(n_entity);
  PDM_free(entity_gnum);

  // Build part comm graph
  pmn->pcg[entity_type] = PDM_part_comm_graph_create(pmn->n_part,
                                                     n_entity_part_bound,
                                                     entity_part_bound,
                                                     PDM_OWNERSHIP_KEEP,
                                                     pmn->comm);
  pmn->pcg_ownership[entity_type] = PDM_OWNERSHIP_KEEP;

  PDM_free(n_entity_part_bound);
  PDM_free(entity_part_bound);
}

void
PDM_part_mesh_nodal_part_comm_graph_deduce_from_vtx
(
  PDM_part_mesh_nodal_t *pmn,
  PDM_geometry_kind_t    geom_kind
)
{

  if (pmn->pcg[PDM_MESH_ENTITY_VTX] != NULL) {
    PDM_error (__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_graph_comm_deduce_from_vtx: pmn->pcg[PDM_MESH_ENTITY_VTX]!=NULL is mandatory in order to deduce other \n");
  }

  PDM_part_mesh_nodal_elmts_t* pmne = NULL;
  if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
    pmne = pmn->surfacic;
  } else if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    pmne = pmn->ridge;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners not implemented for geom_kind %d\n", geom_kind);
  }

  PDM_mesh_entities_t mesh_entity = PDM_geometry_kind_to_entity_type(geom_kind);

  if(pmn->pcg[mesh_entity] != NULL) {
    return; // Already compute
  }

  int  *n_vtx           = NULL;
  int  *n_entity2       = NULL;
  int **entity2_vtx     = NULL;
  int **entity2_vtx_idx = NULL;
  PDM_malloc(n_vtx          , pmn->n_part, int  );
  PDM_malloc(n_entity2      , pmn->n_part, int  );
  PDM_malloc(entity2_vtx    , pmn->n_part, int *);
  PDM_malloc(entity2_vtx_idx, pmn->n_part, int *);

  for(int i_part = 0; i_part < pmn->n_part; ++i_part) {
    n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    n_entity2[i_part] = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne, i_part, &entity2_vtx_idx[i_part], &entity2_vtx[i_part]);
  }

  PDM_part_comm_graph_entity1_to_part_comm_graph_entity2(pmn->pcg[PDM_MESH_ENTITY_VTX],
                                                         n_vtx,
                                                         n_entity2,
                                                         entity2_vtx_idx,
                                                         entity2_vtx,
                                                         &pmn->pcg[mesh_entity]);

  PDM_free(n_entity2      );
  PDM_free(entity2_vtx    );
  PDM_free(entity2_vtx_idx);
}



void
PDM_part_mesh_nodal_compute_straddling_entities
(
  PDM_part_mesh_nodal_t  *pmn,
  PDM_geometry_kind_t     geom_kind,
  PDM_geometry_kind_t     geom_kind_tgt
)
{
  int debug_verbose = 0;
  int debug_visu    = 0;

  if (debug_visu==1) {
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_VOLUMIC , "volume");
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "surface");
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE   , "ridge");
  }

  int i_rank = -1;
  PDM_MPI_Comm_rank(pmn->comm, &i_rank);


  if (!((geom_kind==PDM_GEOMETRY_KIND_SURFACIC && (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE ||
                                                   geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER  ) ) ||
        (geom_kind==PDM_GEOMETRY_KIND_RIDGE    &&  geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER    )  )) {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners cannot build entities of geom_kind %d from entities of geom_kind %d\n", geom_kind_tgt, geom_kind);
  }


  PDM_part_mesh_nodal_elmts_t* pmne = NULL;
  if (geom_kind==PDM_GEOMETRY_KIND_SURFACIC) {
    pmne = pmn->surfacic;
  }
  else if (geom_kind==PDM_GEOMETRY_KIND_RIDGE) {
    pmne = pmn->ridge;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners not implemented for geom_kind %d\n", geom_kind);
  }

  if (pmne==NULL) {
    if (i_rank==0) {
      printf("Warning: no entity of geom_kind %d found\n", geom_kind);
    }
    return;
  }

  /**
   * Detect locally which vertices are on multiple entity1 groups
   */
  int  *n_vtx           = NULL;
  int  *n_entity1       = NULL;
  int **entity1_vtx_idx = NULL;
  int **entity1_vtx     = NULL;
  PDM_malloc(n_vtx          , pmn->n_part, int  );
  PDM_malloc(n_entity1      , pmn->n_part, int  );
  PDM_malloc(entity1_vtx_idx, pmn->n_part, int *);
  PDM_malloc(entity1_vtx    , pmn->n_part, int *);

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    n_vtx    [i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    n_entity1[i_part] = PDM_part_mesh_nodal_elmts_cell_vtx_connect_get(pmne, i_part,
                                                                      &entity1_vtx_idx[i_part],
                                                                      &entity1_vtx[i_part]);
  }

  int **tag_vtx_n       = NULL;
  int **tag_vtx_idx     = NULL;
  int **tag_vtx         = NULL;
  _find_local_parent_group(pmne,
                           pmn->n_part,
                           entity1_vtx_idx,
                           entity1_vtx,
                           n_vtx,
                          &tag_vtx_n,
                          &tag_vtx_idx,
                          &tag_vtx);

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(entity1_vtx_idx[i_part]);
    PDM_free(entity1_vtx    [i_part]);
  }
  PDM_free(entity1_vtx_idx);
  PDM_free(entity1_vtx);


  /**
   * Exchange local information to reduce it globally
   */
  if (pmn->pcg[PDM_MESH_ENTITY_VTX]==NULL) {
    PDM_part_mesh_nodal_part_comm_graph_compute_from_gnum(pmn, PDM_MESH_ENTITY_VTX);
  }


  if (debug_verbose==1) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_int(tag_vtx_idx[i_part], tag_vtx[i_part], n_vtx[i_part], "local tag_vtx :: ");
    }
  }

  int **_tag_vtx_n = NULL;
  int **_tag_vtx   = NULL;
  PDM_part_comm_graph_gather_strided_data(pmn->pcg[PDM_MESH_ENTITY_VTX],
                                          1*sizeof(int),
                                          PDM_STRIDE_CST_INTERLACED,
                                          n_vtx,
                                          tag_vtx_n,
                             (void  **)   tag_vtx,
                                        &_tag_vtx_n,
                             (void ***) &_tag_vtx);

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_n  [i_part]);
    PDM_free(tag_vtx_idx[i_part]);
    PDM_free(tag_vtx    [i_part]);
    tag_vtx_idx[i_part] = PDM_array_new_idx_from_sizes_int(_tag_vtx_n[i_part], n_vtx[i_part]);
  }
  PDM_free(tag_vtx_n);
  PDM_free(tag_vtx);
  tag_vtx_n = _tag_vtx_n;
  tag_vtx   = _tag_vtx;


  if (debug_verbose==1) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_log_trace_connectivity_int(tag_vtx_idx[i_part], tag_vtx[i_part], n_vtx[i_part], "global tag_vtx :: ");
    }
  }



  int         **vtx_group_idx  = NULL;
  PDM_g_num_t **vtx_group_gnum = NULL;
  PDM_l_num_t **ridge_vtx_candidate = NULL;
  if (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE) {
    /**
     * If geom_kind_tgt is ridge, generate vtx group gnum is not necessary,
     * we only need to get straddling vertices to generate edges containing them
     */
    PDM_malloc(ridge_vtx_candidate, pmn->n_part, PDM_l_num_t *);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {

      PDM_malloc(ridge_vtx_candidate[i_part], n_vtx[i_part], PDM_l_num_t);

      for(int i_vtx = 0; i_vtx < n_vtx[i_part]; ++i_vtx) {
        int beg = tag_vtx_idx[i_part][i_vtx  ];
        int end = tag_vtx_idx[i_part][i_vtx+1];
        if(end - beg > 0){
          int n_unique = PDM_inplace_unique(tag_vtx[i_part], beg, end-1);
          ridge_vtx_candidate[i_part][i_vtx] = n_unique > 1;
        }
      }
    }
  } else {
    _generate_group_gnum(pmn->comm,
                         pmn->n_part,
                         n_vtx,
                         tag_vtx_idx,
                         tag_vtx,
                        &vtx_group_idx,
                        &vtx_group_gnum);
  }
  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_idx[i_part]);
    PDM_free(tag_vtx    [i_part]);
  }
  PDM_free(tag_vtx_idx);
  PDM_free(tag_vtx);


  if (geom_kind_tgt==PDM_GEOMETRY_KIND_CORNER) {

    if (pmn->corner!=NULL) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners : part_mesh_nodal already has corner section\n");
    }

    /**
     * Count global n_group while casting group id into integer
     * and transpose vtx_group for pmne storage
     */
    int **group_vtx_idx = NULL;
    int **group_vtx     = NULL;
    int  g_n_group = 0;
    _transpose_group_information(pmn->comm,
                                 pmn->n_part,
                                 n_vtx,
                                 vtx_group_idx,
                                 vtx_group_gnum,
                                &g_n_group,
                                &group_vtx_idx,
                                &group_vtx);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(vtx_group_idx [i_part]);
      PDM_free(vtx_group_gnum[i_part]);
    }
    PDM_free(vtx_group_idx);
    PDM_free(vtx_group_gnum);



    /**
     * Generate global id for group entities, then store corners in pmne
     */
    PDM_g_num_t **vtx_gnum       = NULL;
    PDM_g_num_t **group_vtx_gnum = NULL;
    PDM_malloc(vtx_gnum, pmn->n_part, PDM_g_num_t*);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      vtx_gnum[i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    }

    _generate_group_entity_gnum(pmn->comm,
                                pmn->n_part,
                                vtx_gnum,
                                g_n_group,
                                group_vtx_idx,
                                group_vtx,
                               &group_vtx_gnum);

    PDM_free(vtx_gnum);


    pmn->corner = PDM_part_mesh_nodal_elmts_create(0, pmn->n_part, pmn->comm);
    int section_corner = PDM_part_mesh_nodal_elmts_add(pmn->corner, PDM_MESH_NODAL_POINT);
    PDM_part_mesh_nodal_elmts_n_group_set(pmn->corner, g_n_group);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_part_mesh_nodal_elmts_std_set(pmn->corner,
                                        section_corner,
                                        i_part,
                                        group_vtx_idx[i_part][g_n_group],
                                        group_vtx[i_part],
                                        NULL,
                                        NULL,
                                        NULL,
                                        PDM_OWNERSHIP_KEEP);
      for (int i_group=0; i_group<g_n_group; ++i_group) {
        int n_group_vtx = group_vtx_idx[i_part][i_group+1]-group_vtx_idx[i_part][i_group];
        int         *_group_vtx      = NULL;
        PDM_g_num_t *_group_vtx_gnum = NULL;
        PDM_malloc(_group_vtx     , n_group_vtx, int);
        PDM_malloc(_group_vtx_gnum, n_group_vtx, PDM_g_num_t);
        memcpy(_group_vtx     , &group_vtx     [i_part][group_vtx_idx[i_part][i_group]], n_group_vtx*sizeof(int        ));
        memcpy(_group_vtx_gnum, &group_vtx_gnum[i_part][group_vtx_idx[i_part][i_group]], n_group_vtx*sizeof(PDM_g_num_t));
        PDM_part_mesh_nodal_elmts_group_set(pmn->corner, i_part,
                                            i_group, n_group_vtx,
                                            _group_vtx, _group_vtx_gnum,
                                            PDM_OWNERSHIP_KEEP);
      }
      PDM_free(group_vtx_idx [i_part]);
      PDM_free(group_vtx_gnum[i_part]);
    }
    PDM_free(group_vtx_gnum);
    PDM_free(group_vtx_idx);
    PDM_free(group_vtx);

  }
  else { //edge

    if (pmn->ridge!=NULL) {
      PDM_error(__FILE__, __LINE__, 0, "PDM_part_mesh_nodal_compute_topo_corners : part_mesh_nodal already has ridge section\n");
    }

    /**
     * Decompose face edges then unique them before searching parent groups
     */
    int  *n_edge = NULL;
    int **entity1_edge_idx     = NULL;
    int **entity1_edge         = NULL;
    int **entity1_edge_vtx_idx = NULL;
    int **entity1_edge_vtx     = NULL;
    int **edge_parent          = NULL;
    int **edge_parent_pos      = NULL;
    PDM_part_mesh_nodal_elmts_sections_local_decompose_edges(pmne,
                                                             ridge_vtx_candidate,
                                                            &n_edge,
                                                            &entity1_edge_idx,
                                                            &entity1_edge_vtx_idx,
                                                            &entity1_edge_vtx,
                                                            &edge_parent,
                                                            &edge_parent_pos);

    if (debug_verbose==1) {
      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        log_trace("n_edge = %d\n", n_edge[i_part]);
        PDM_log_trace_array_int(entity1_edge_idx    [i_part], n_entity1[i_part], "entity1_edge_idx::");
        PDM_log_trace_array_int(entity1_edge_vtx_idx[i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_edge_vtx_idx::");
        PDM_log_trace_array_int(entity1_edge_vtx    [i_part], entity1_edge_vtx_idx[i_part][entity1_edge_idx[i_part][n_entity1[i_part]]], "entity1_edge_vtx::");
        PDM_log_trace_array_int(edge_parent         [i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_parent::");
        PDM_log_trace_array_int(edge_parent_pos     [i_part], entity1_edge_idx    [i_part][                         n_entity1[i_part]] , "entity1_parent_pos::");
      }
    }


    int  *n_unique_edge           = NULL;
    int **unique_edge_vtx_idx     = NULL;
    int **unique_edge_vtx         = NULL;
    int **unique_edge_to_edge_idx = NULL;
    int **unique_edge_to_edge     = NULL;
    PDM_malloc(n_unique_edge          , pmn->n_part, int  );
    PDM_malloc(unique_edge_vtx_idx    , pmn->n_part, int *);
    PDM_malloc(unique_edge_vtx        , pmn->n_part, int *);
    PDM_malloc(unique_edge_to_edge_idx, pmn->n_part, int *);
    PDM_malloc(unique_edge_to_edge    , pmn->n_part, int *);
    PDM_malloc(entity1_edge           , pmn->n_part, int *);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {

      entity1_edge[i_part] = PDM_array_new_idx_from_const_stride_int(1, n_edge[i_part]-1);
      PDM_malloc(unique_edge_vtx_idx[i_part],                              n_edge[i_part]+ 1, int);
      PDM_malloc(unique_edge_vtx    [i_part], entity1_edge_vtx_idx[i_part][n_edge[i_part]]  , int);
      unique_edge_vtx_idx[i_part][0] = 0;

      PDM_part_mesh_nodal_elmts_generate_entity_connectivity(n_edge              [i_part],
                                                             entity1_edge_vtx_idx[i_part],
                                                             entity1_edge_vtx    [i_part],
                                                             edge_parent         [i_part],
                                                             0,
                                                             NULL,
                                                             NULL,
                                                            &unique_edge_to_edge_idx[i_part],
                                                            &unique_edge_to_edge    [i_part],
                                                             n_entity1              [i_part],
                                                             entity1_edge_idx       [i_part],
                                                             entity1_edge           [i_part],
                                                            &n_unique_edge          [i_part],
                                                             unique_edge_vtx_idx    [i_part],
                                                             unique_edge_vtx        [i_part],
                                                             PDM_TRUE);
      PDM_free(entity1_edge_vtx_idx   [i_part]);
      PDM_free(entity1_edge_vtx       [i_part]);
      PDM_free(edge_parent            [i_part]);
      PDM_free(edge_parent_pos        [i_part]);
      PDM_free(unique_edge_to_edge_idx[i_part]);
      PDM_free(unique_edge_to_edge    [i_part]);
    }
    PDM_free(n_edge);
    PDM_free(entity1_edge_vtx_idx);
    PDM_free(entity1_edge_vtx);
    PDM_free(edge_parent);
    PDM_free(edge_parent_pos);
    PDM_free(unique_edge_to_edge_idx);
    PDM_free(unique_edge_to_edge);


    PDM_g_num_t **unique_edge_gnum = NULL;
    if (pmn->pcg[PDM_MESH_ENTITY_EDGE]==NULL) {

      pmn->ridge = PDM_part_mesh_nodal_elmts_create(1, pmn->n_part, pmn->comm);
      int section_ridge = PDM_part_mesh_nodal_elmts_add(pmn->ridge, PDM_MESH_NODAL_BAR2);

      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        PDM_part_mesh_nodal_elmts_std_set(pmn->ridge,
                                          section_ridge,
                                          i_part,
                                          n_unique_edge   [i_part],
                                          unique_edge_vtx [i_part],
                                          NULL,
                                          NULL,
                                          NULL,
                                          PDM_OWNERSHIP_USER);
      }


      // > Get vertices part_comm_graph
      int  *pn_vtx_graph = NULL;
      int **pvtx_graph   = NULL;
      PDM_malloc(pn_vtx_graph, pmn->n_part, int  );
      PDM_malloc(pvtx_graph  , pmn->n_part, int *);
      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        pn_vtx_graph[i_part] = PDM_part_comm_graph_entity_graph_get(pmn->pcg[PDM_MESH_ENTITY_VTX],
                                                                    i_part,
                                                                   &pvtx_graph[i_part],
                                                                    PDM_OWNERSHIP_BAD_VALUE);
      }


      // > Compute edge part_comm_graph from vertices pcg + edge->vtx connectivity
      int  *pn_edge_graph = NULL;
      int **pedge_graph   = NULL;
      PDM_part_comm_graph_entity1_to_entity2(pmn->comm,
                                             pmn->n_part,
                                             pn_vtx_graph,
                                             pvtx_graph,
                                             0,
                                             NULL,
                                             n_vtx,
                                             n_unique_edge,
                                             unique_edge_vtx_idx,
                                             unique_edge_vtx,
                                            &pn_edge_graph,
                                            &pedge_graph,
                                             NULL);

      for (int i_part=0; i_part<pmn->n_part; ++i_part) {
        PDM_free(unique_edge_vtx_idx[i_part]);
      }
      PDM_free(unique_edge_vtx_idx);
      PDM_free(pn_vtx_graph);
      PDM_free(pvtx_graph);

      PDM_part_comm_graph_t *pcg_edge = PDM_part_comm_graph_create(pmn->n_part,
                                                                   pn_edge_graph,
                                                                   pedge_graph,
                                                                   PDM_OWNERSHIP_KEEP,
                                                                   pmn->comm);
      PDM_free(pn_edge_graph);
      PDM_free(pedge_graph);

      PDM_part_mesh_nodal_part_comm_graph_set(pmn, pcg_edge, PDM_MESH_ENTITY_EDGE, PDM_OWNERSHIP_KEEP);
    }



    if (debug_visu==1) {
      PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE, "edge_detected");
    }


    int **tag_edge_n       = NULL;
    int **tag_edge_idx     = NULL;
    int **tag_edge         = NULL;
    _find_local_parent_group(pmne,
                             pmn->n_part,
                             entity1_edge_idx,
                             entity1_edge,
                             n_unique_edge,
                            &tag_edge_n,
                            &tag_edge_idx,
                            &tag_edge);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(tag_edge_idx    [i_part]);
      PDM_free(entity1_edge_idx[i_part]);
      PDM_free(entity1_edge    [i_part]);
    }
    PDM_free(entity1_edge_idx);
    PDM_free(entity1_edge);


    int **_tag_edge_n = NULL;
    int **_tag_edge   = NULL;
    PDM_part_comm_graph_gather_strided_data(pmn->pcg[PDM_MESH_ENTITY_EDGE],
                                            1*sizeof(int),
                                            PDM_STRIDE_CST_INTERLACED,
                                            n_unique_edge,
                                            tag_edge_n,
                               (void  **)   tag_edge,
                                          &_tag_edge_n,
                               (void ***) &_tag_edge);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(tag_edge_n [i_part]);
      PDM_free(tag_edge   [i_part]);
      tag_edge_idx[i_part] = PDM_array_new_idx_from_sizes_int(_tag_edge_n[i_part], n_unique_edge[i_part]);
      PDM_free(_tag_edge_n[i_part]);
    }
    PDM_free(tag_edge_n);
    PDM_free(tag_edge);
    PDM_free(_tag_edge_n);
    tag_edge = _tag_edge;

    int         **edge_group_idx  = NULL;
    PDM_g_num_t **edge_group_gnum = NULL;
    _generate_group_gnum(pmn->comm,
                         pmn->n_part,
                         n_unique_edge,
                         tag_edge_idx,
                         tag_edge,
                        &edge_group_idx,
                        &edge_group_gnum);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(tag_edge_idx[i_part]);
      PDM_free(tag_edge    [i_part]);
    }
    PDM_free(tag_edge_idx);
    PDM_free(tag_edge);


    /**
     * All edge created arent ridges, so we need to extract them and regenerate a gnum
     */
    int          *n_ridge          = NULL;
    int         **ridge_vtx        = NULL;
    PDM_g_num_t **ridge_gnum       = NULL;
    int         **ridge_group_idx  = NULL;
    PDM_g_num_t **ridge_group_gnum = NULL;
    PDM_malloc(n_ridge        , pmn->n_part, int          );
    PDM_malloc(ridge_vtx      , pmn->n_part, int         *);
    PDM_malloc(ridge_gnum     , pmn->n_part, PDM_g_num_t *);
    PDM_malloc(ridge_group_idx, pmn->n_part, int         *);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      n_ridge[i_part] = 0;
      for (int i_edge=0; i_edge<n_unique_edge[i_part]; ++i_edge) {
        if (edge_group_idx[i_part][i_edge+1]-edge_group_idx[i_part][i_edge]>0) {
          n_ridge[i_part]++;
        }
      }

      PDM_g_num_t *vtx_gnum = PDM_part_mesh_nodal_vtx_g_num_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

      PDM_malloc(ridge_vtx      [i_part], 2*n_ridge[i_part]  , int        );
      PDM_malloc(ridge_gnum     [i_part], 2*n_ridge[i_part]  , PDM_g_num_t);
      PDM_malloc(ridge_group_idx[i_part],   n_ridge[i_part]+1, int        );
      ridge_group_idx[i_part][0] = 0;
      int i_write = 0;
      for (int i_edge=0; i_edge<n_unique_edge[i_part]; ++i_edge) {
        if (edge_group_idx[i_part][i_edge+1]-edge_group_idx[i_part][i_edge]>0) {
          int vtx1 = unique_edge_vtx[i_part][2*i_edge  ];
          int vtx2 = unique_edge_vtx[i_part][2*i_edge+1];
          ridge_vtx      [i_part][2*i_write  ] = vtx1;
          ridge_vtx      [i_part][2*i_write+1] = vtx2;
          ridge_gnum     [i_part][2*i_write  ] = vtx_gnum[vtx1-1]; // sure about -1 ?
          ridge_gnum     [i_part][2*i_write+1] = vtx_gnum[vtx2-1]; // sure about -1 ?
          ridge_group_idx[i_part][  i_write+1] = ridge_group_idx[i_part][i_write]+1;
          i_write++;
        }
      }
      PDM_free(unique_edge_vtx[i_part]);
      PDM_free(edge_group_idx [i_part]);
    }
    PDM_free(n_unique_edge);
    PDM_free(unique_edge_vtx);
    PDM_free(unique_edge_gnum);
    PDM_free(edge_group_idx);

    PDM_gen_gnum_t* gen_ridge_gnum = PDM_gnum_create(3,
                                                     pmn->n_part,
                                                     PDM_TRUE,
                                                     1e-4,
                                                     pmn->comm,
                                                     PDM_OWNERSHIP_USER);
    PDM_gnum_set_parents_nuplet(gen_ridge_gnum, 2);
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_gnum_set_from_parents(gen_ridge_gnum, i_part, n_ridge[i_part], ridge_gnum[i_part]);
    }

    PDM_gnum_compute(gen_ridge_gnum);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_gnum[i_part]);
      ridge_gnum[i_part] = PDM_gnum_get(gen_ridge_gnum, i_part);
    }

    PDM_gnum_free(gen_ridge_gnum);
    ridge_group_gnum = edge_group_gnum;


    /**
     * Count global n_group while casting group id into integer
     * and transpose vtx_group for pmne storage
     */
    int **group_ridge_idx = NULL;
    int **group_ridge     = NULL;
    int  g_n_group = 0;
    _transpose_group_information(pmn->comm,
                                 pmn->n_part,
                                 n_ridge,
                                 ridge_group_idx,
                                 ridge_group_gnum,
                                &g_n_group,
                                &group_ridge_idx,
                                &group_ridge);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_group_idx [i_part]);
      PDM_free(ridge_group_gnum[i_part]);
    }
    PDM_free(ridge_group_idx);
    PDM_free(ridge_group_gnum);


    /**
     * Generate global id for group entities, then store corners in pmne
     */

    PDM_g_num_t **group_ridge_gnum = NULL;
    _generate_group_entity_gnum(pmn->comm,
                                pmn->n_part,
                                ridge_gnum,
                                g_n_group,
                                group_ridge_idx,
                                group_ridge,
                               &group_ridge_gnum);


    PDM_part_mesh_nodal_elmts_free(pmn->ridge);

    pmn->ridge = PDM_part_mesh_nodal_elmts_create(1, pmn->n_part, pmn->comm);
    int section_ridge = PDM_part_mesh_nodal_elmts_add(pmn->ridge, PDM_MESH_NODAL_BAR2);
    PDM_part_mesh_nodal_elmts_n_group_set(pmn->ridge, g_n_group);

    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_part_mesh_nodal_elmts_std_set(pmn->ridge,
                                        section_ridge,
                                        i_part,
                                        n_ridge   [i_part],
                                        ridge_vtx [i_part],
                                        ridge_gnum[i_part],
                                        NULL,
                                        NULL,
                                        PDM_OWNERSHIP_KEEP);
      for (int i_group=0; i_group<g_n_group; ++i_group) {
        int n_group_ridge = group_ridge_idx[i_part][i_group+1]-group_ridge_idx[i_part][i_group];
        int         *_group_ridge      = NULL;
        PDM_g_num_t *_group_ridge_gnum = NULL;
        PDM_malloc(_group_ridge     , n_group_ridge, int);
        PDM_malloc(_group_ridge_gnum, n_group_ridge, PDM_g_num_t);
        memcpy(_group_ridge     , &group_ridge     [i_part][group_ridge_idx[i_part][i_group]], n_group_ridge*sizeof(int        ));
        memcpy(_group_ridge_gnum, &group_ridge_gnum[i_part][group_ridge_idx[i_part][i_group]], n_group_ridge*sizeof(PDM_g_num_t));
        PDM_part_mesh_nodal_elmts_group_set(pmn->ridge,
                                            i_part,
                                            i_group,
                                            n_group_ridge,
                                            _group_ridge,
                                            _group_ridge_gnum,
                                            PDM_OWNERSHIP_KEEP);
      }
      PDM_free(group_ridge_idx [i_part]);
      PDM_free(group_ridge     [i_part]);
      PDM_free(group_ridge_gnum[i_part]);

      if (debug_visu==1) {
        PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_RIDGE, "ridge_final");
      }
    }
    PDM_free(n_ridge);
    PDM_free(ridge_vtx);
    PDM_free(ridge_gnum);
    PDM_free(group_ridge_gnum);
    PDM_free(group_ridge);
    PDM_free(group_ridge_idx);
  }

  for (int i_part=0; i_part<pmn->n_part; ++i_part) {
    PDM_free(tag_vtx_n[i_part]);
  }
  PDM_free(tag_vtx_n);

  if (geom_kind_tgt==PDM_GEOMETRY_KIND_RIDGE) {
    for (int i_part=0; i_part<pmn->n_part; ++i_part) {
      PDM_free(ridge_vtx_candidate[i_part]);
    }
    PDM_free(ridge_vtx_candidate);
  }

  PDM_free(n_entity1);
  PDM_free(n_vtx);

}

#ifdef __cplusplus
}
#endif /* __cplusplus */
