#include <stddef.h>
#include <vector>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_graph_dual.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_vtk.h"


static
PDM_part_mesh_nodal_t*
_generate_mesh
(
  PDM_MPI_Comm    pdm_comm,
  int             n_vtx_seg
)
{
  int              n_part       = 1;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;

  /* Warmup */
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        2.,
                                                        -1,
                                                        -1,
                                                        -1,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                split_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t* pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmesh_nodal, PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_complete_part_comm_graph(pmesh_nodal);

  PDM_multipart_free(mpart);
  PDM_dcube_nodal_gen_free(dcube);
  return pmesh_nodal;
}

static
PDM_part_mesh_t*
_part_mesh_nodal_to_part_mesh
(
  int                    dim,
  PDM_part_mesh_nodal_t *pmn
)
{
  PDM_part_mesh_nodal_to_part_mesh_t* pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                          PDM_FALSE,
                                                                                          PDM_OWNERSHIP_USER);

  PDM_part_mesh_t *pm = NULL;
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_EDGE);
  if(dim == 3) {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX);
  } else {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX);
  }

  PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

  PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                 &pm,
                                                 PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);

  // Finalize part_mesh
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);

  return pm;
}

static
void
_warmup_node_centered
(
  PDM_part_mesh_t   *pm,
  int              **out_pn_node,
  int              **out_pn_arc,
  int             ***out_parc_node_idx,
  int             ***out_parc_node,
  int             ***out_pnode_weight,
  int             ***out_parc_weight
)
{
  int n_part = PDM_part_mesh_n_part_get(pm);

  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;
  PDM_malloc(pn_node      , n_part, int  );
  PDM_malloc(pn_arc       , n_part, int  );
  PDM_malloc(parc_node_idx, n_part, int *);
  PDM_malloc(parc_node    , n_part, int *);
  PDM_malloc(pnode_weight , n_part, int *);
  PDM_malloc(parc_weight  , n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_node[i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VTX );
    pn_arc [i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &parc_node    [i_part],
                                   &parc_node_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    if(parc_node_idx[i_part] == NULL) {
      PDM_malloc(parc_node_idx[i_part], pn_arc[i_part] + 1, int);
      for(int i = 0; i < pn_arc[i_part]+1; ++i) {
        parc_node_idx[i_part][i] = 2*i;
      }
    }

    PDM_malloc(pnode_weight[i_part], pn_node[i_part], int);
    PDM_malloc(parc_weight [i_part], pn_arc [i_part], int);

    for(int i = 0; i < pn_node[i_part]; ++i) {
      pnode_weight[i_part][i] = 1;
    }
    for(int i = 0; i < pn_arc[i_part]; ++i) {
      parc_weight[i_part][i] = 1;
    }
  }

  *out_pn_node       = pn_node;
  *out_pn_arc        = pn_arc;
  *out_parc_node_idx = parc_node_idx;
  *out_parc_node     = parc_node;
  *out_pnode_weight  = pnode_weight;
  *out_parc_weight   = parc_weight;
}



static
void
_warmup_cell_centered
(
  PDM_part_mesh_t   *pm,
  int              **out_pn_node,
  int              **out_pn_arc,
  int             ***out_pnode_arc_idx,
  int             ***out_pnode_arc,
  int             ***out_pnode_weight,
  int             ***out_parc_weight
)
{
  int n_part = PDM_part_mesh_n_part_get(pm);

  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;
  PDM_malloc(pn_node      , n_part, int  );
  PDM_malloc(pn_arc       , n_part, int  );
  PDM_malloc(pnode_arc_idx, n_part, int *);
  PDM_malloc(pnode_arc    , n_part, int *);
  PDM_malloc(pnode_weight , n_part, int *);
  PDM_malloc(parc_weight  , n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_node[i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_FACE);
    pn_arc [i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &pnode_arc    [i_part],
                                   &pnode_arc_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    PDM_malloc(pnode_weight[i_part], pn_node[i_part], int);
    PDM_malloc(parc_weight [i_part], pn_arc [i_part], int);

    for(int i = 0; i < pn_node[i_part]; ++i) {
      pnode_weight[i_part][i] = 1;
    }
    for(int i = 0; i < pn_arc[i_part]; ++i) {
      parc_weight[i_part][i] = 1;
    }
  }

  *out_pn_node       = pn_node;
  *out_pn_arc        = pn_arc;
  *out_pnode_arc_idx = pnode_arc_idx;
  *out_pnode_arc     = pnode_arc;
  *out_pnode_weight  = pnode_weight;
  *out_parc_weight   = parc_weight;
}

MPI_TEST_CASE("[PDM_part_assembly_dual_graph] Full vtx_centered ", 2) {

  int n_rank;
  PDM_MPI_Comm_size(test_comm, &n_rank);

  int n_vtx_a = 4;
  PDM_part_mesh_nodal_t* pmn = _generate_mesh(test_comm, n_vtx_a);

  PDM_part_mesh_t* pm = _part_mesh_nodal_to_part_mesh(2, pmn);

  /*
   * Node = vertices
   * Arc  = edges
   */
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pselect_node  = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;
  _warmup_node_centered(pm,
                        &pn_node,
                        &pn_arc,
                        &parc_node_idx,
                        &parc_node,
                        &pnode_weight,
                        &parc_weight);

  // Transpose
  PDM_part_connectivity_transpose(n_part,
                                  pn_arc,
                                  pn_node,
                                  parc_node_idx,
                                  parc_node,
                                  &pnode_arc_idx,
                                  &pnode_arc);

  PDM_part_comm_graph_t *pcg_node = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_VTX,
                                    &pcg_node,
                                    PDM_OWNERSHIP_KEEP);

  PDM_part_comm_graph_t *pcg_arc = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_EDGE,
                                    &pcg_arc,
                                    PDM_OWNERSHIP_KEEP);

  int           n_tot_node     = 0;
  PDM_g_num_t  *gnode_node_idx = NULL;
  PDM_g_num_t  *gnode_node     = NULL;
  int          *gnode_weight   = NULL;
  int          *garc_weight    = NULL;
  PDM_g_num_t  *distrib_node   = NULL;
  int         **part_to_graph  = NULL;
  PDM_part_assembly_dual_graph(test_comm,
                               n_part,
                               pn_node,
                               pn_arc,
                               pselect_node,
                               pnode_arc_idx,
                               pnode_arc,
                               parc_node_idx,
                               parc_node,
                               pnode_weight,
                               parc_weight,
                               pcg_node,
                               pcg_arc,
                               &n_tot_node,
                               &gnode_node_idx,
                               &gnode_node,
                               &gnode_weight,
                               &garc_weight,
                               &distrib_node,
                               &part_to_graph);

  if(0 == 1) {
    PDM_log_trace_array_long(gnode_node_idx,       n_tot_node+1              , "gnode_node_idx ::");
    PDM_log_trace_array_long(gnode_node    , (int) gnode_node_idx[n_tot_node], "gnode_node     ::");
    PDM_log_trace_array_int (gnode_weight  ,       n_tot_node                , "gnode_weight   ::");
    PDM_log_trace_array_int (garc_weight   , (int) gnode_node_idx[n_tot_node], "garc_weight    ::");
    PDM_log_trace_array_long(distrib_node  ,       n_rank+1                  , "distrib_node   ::");
  }

  MPI_CHECK(0, n_tot_node == 10);
  MPI_CHECK(1, n_tot_node == 6 );

  PDM_g_num_t expected_p0_gnode_node_idx[11] = {0, 6, 10, 14, 20, 26, 30, 33, 37, 41, 43};
  PDM_g_num_t expected_p0_gnode_node    [43] = {1, 3, 4, 12, 13, 15, 0, 4, 5, 13, 3, 6, 14, 15, 0, 2, 4, 6, 7, 15, 0, 1, 3, 5, 7, 8, 1, 4, 8, 9, 2, 3, 7, 3, 4, 6, 8, 4, 5, 7, 9, 5, 8};
  PDM_g_num_t expected_p0_gnode_weight  [10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p0_garc_weight   [43] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p0_distrib_node  [3 ] = {0, 10, 16};

  PDM_g_num_t expected_p1_gnode_node_idx[7 ] = {0, 2, 6, 10, 13, 17, 23};
  PDM_g_num_t expected_p1_gnode_node    [23] = {11, 14, 10, 12, 14, 15, 0, 11, 13, 15, 0, 1, 12, 2, 10, 11, 15, 0, 2, 3, 11, 12, 14};
  PDM_g_num_t expected_p1_gnode_weight  [6 ] = {1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p1_garc_weight   [23] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p1_distrib_node  [3 ] = {0, 10, 16};

  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node_idx, gnode_node_idx, 11);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node    , gnode_node    , 43);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_weight  , gnode_weight  , 10);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_garc_weight   , garc_weight   , 43);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_distrib_node  , distrib_node  , 3 );

  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node_idx, gnode_node_idx, 7 );
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node    , gnode_node    , 23);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_weight  , gnode_weight  , 6 );
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_garc_weight   , garc_weight   , 23);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_distrib_node  , distrib_node  , 3 );


  PDM_free(gnode_node_idx);
  PDM_free(gnode_node    );
  PDM_free(gnode_weight  );
  PDM_free(garc_weight   );
  PDM_free(distrib_node  );

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(part_to_graph[i_part]);
  }
  PDM_free(part_to_graph);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_arc_idx[i_part]);
    PDM_free(pnode_arc    [i_part]);
    PDM_free(pnode_weight [i_part]);
    PDM_free(parc_weight  [i_part]);
  }
  PDM_free(pn_node      );
  PDM_free(pn_arc       );
  PDM_free(pselect_node );
  PDM_free(pnode_arc_idx);
  PDM_free(pnode_arc    );
  PDM_free(parc_node_idx);
  PDM_free(parc_node    );
  PDM_free(pnode_weight );
  PDM_free(parc_weight  );

  PDM_part_mesh_nodal_free(pmn);
  PDM_part_mesh_free(pm);
}


MPI_TEST_CASE("[PDM_part_assembly_dual_graph] select vtx_centered ", 2) {

  int n_rank;
  PDM_MPI_Comm_size(test_comm, &n_rank);

  int n_vtx_a = 10;
  PDM_part_mesh_nodal_t* pmn = _generate_mesh(test_comm, n_vtx_a);

  PDM_part_mesh_t* pm = _part_mesh_nodal_to_part_mesh(2, pmn);

  /*
   * Node = vertices
   * Arc  = edges
   */
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pselect_node  = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;
  _warmup_node_centered(pm,
                        &pn_node,
                        &pn_arc,
                        &parc_node_idx,
                        &parc_node,
                        &pnode_weight,
                        &parc_weight);

  // Transpose
  PDM_part_connectivity_transpose(n_part,
                                  pn_arc,
                                  pn_node,
                                  parc_node_idx,
                                  parc_node,
                                  &pnode_arc_idx,
                                  &pnode_arc);

  PDM_part_comm_graph_t *pcg_node = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_VTX,
                                    &pcg_node,
                                    PDM_OWNERSHIP_KEEP);

  PDM_part_comm_graph_t *pcg_arc = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_EDGE,
                                    &pcg_arc,
                                    PDM_OWNERSHIP_KEEP);

  PDM_malloc(pselect_node, n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_node[i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VTX );
    PDM_malloc(pselect_node[i_part], pn_node[i_part], int);
    double *vtx_coords = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);
    for(int i = 0; i < pn_node[i_part]; ++i) {
      if( (vtx_coords[3*i  ] > -0.25 && vtx_coords[3*i  ] < 0.25) &&
          (vtx_coords[3*i+1] > -0.25 && vtx_coords[3*i+1] < 0.25)){
        pselect_node[i_part][i] = 1;
      } else {
        pselect_node[i_part][i] = 0;
      }
    }
  }

  int           n_tot_node     = 0;
  PDM_g_num_t  *gnode_node_idx = NULL;
  PDM_g_num_t  *gnode_node     = NULL;
  int          *gnode_weight   = NULL;
  int          *garc_weight    = NULL;
  PDM_g_num_t  *distrib_node   = NULL;
  int         **part_to_graph  = NULL;
  PDM_part_assembly_dual_graph(test_comm,
                               n_part,
                               pn_node,
                               pn_arc,
                               pselect_node,
                               pnode_arc_idx,
                               pnode_arc,
                               parc_node_idx,
                               parc_node,
                               pnode_weight,
                               parc_weight,
                               pcg_node,
                               pcg_arc,
                               &n_tot_node,
                               &gnode_node_idx,
                               &gnode_node,
                               &gnode_weight,
                               &garc_weight,
                               &distrib_node,
                               &part_to_graph);

  if(0 == 1) {
    PDM_log_trace_array_long(gnode_node_idx,       n_tot_node+1              , "gnode_node_idx ::");
    PDM_log_trace_array_long(gnode_node    , (int) gnode_node_idx[n_tot_node], "gnode_node     ::");
    PDM_log_trace_array_int (gnode_weight  ,       n_tot_node                , "gnode_weight   ::");
    PDM_log_trace_array_int (garc_weight   , (int) gnode_node_idx[n_tot_node], "garc_weight    ::");
    PDM_log_trace_array_long(distrib_node  ,       n_rank+1                  , "distrib_node   ::");
  }

  MPI_CHECK(0, n_tot_node == 3);
  MPI_CHECK(1, n_tot_node == 1);

  PDM_g_num_t expected_p0_gnode_node_idx[4] = {0, 2, 5, 8};
  PDM_g_num_t expected_p0_gnode_node    [8] = {1, 2, 0, 2, 3, 0, 1, 3 };
  PDM_g_num_t expected_p0_gnode_weight  [3] = {1, 1, 1};
  PDM_g_num_t expected_p0_garc_weight   [8] = {1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p0_distrib_node  [3] = {0, 3, 4 };

  PDM_g_num_t expected_p1_gnode_node_idx[2] = {0, 2};
  PDM_g_num_t expected_p1_gnode_node    [2] = {1, 2};
  PDM_g_num_t expected_p1_gnode_weight  [1] = {1};
  PDM_g_num_t expected_p1_garc_weight   [2] = {1, 1};
  PDM_g_num_t expected_p1_distrib_node  [3] = {0, 3, 4};

  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node_idx, gnode_node_idx, 4);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node    , gnode_node    , 8);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_weight  , gnode_weight  , 3);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_garc_weight   , garc_weight   , 8);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_distrib_node  , distrib_node  , 3);

  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node_idx, gnode_node_idx, 2);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node    , gnode_node    , 2);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_weight  , gnode_weight  , 1);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_garc_weight   , garc_weight   , 2);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_distrib_node  , distrib_node  , 3);


  PDM_free(gnode_node_idx);
  PDM_free(gnode_node    );
  PDM_free(gnode_weight  );
  PDM_free(garc_weight   );
  PDM_free(distrib_node  );

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(part_to_graph[i_part]);
  }
  PDM_free(part_to_graph);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_arc_idx[i_part]);
    PDM_free(pnode_arc    [i_part]);
    PDM_free(pnode_weight [i_part]);
    PDM_free(parc_weight  [i_part]);
    PDM_free(pselect_node [i_part]);
  }
  PDM_free(pn_node      );
  PDM_free(pn_arc       );
  PDM_free(pselect_node );
  PDM_free(pnode_arc_idx);
  PDM_free(pnode_arc    );
  PDM_free(parc_node_idx);
  PDM_free(parc_node    );
  PDM_free(pnode_weight );
  PDM_free(parc_weight  );

  PDM_part_mesh_nodal_free(pmn);
  PDM_part_mesh_free(pm);


}



MPI_TEST_CASE("[PDM_part_assembly_dual_graph] Full cell_centered ", 2) {

  int n_rank;
  PDM_MPI_Comm_size(test_comm, &n_rank);

  int n_vtx_a = 4;
  PDM_part_mesh_nodal_t* pmn = _generate_mesh(test_comm, n_vtx_a);

  PDM_part_mesh_t* pm = _part_mesh_nodal_to_part_mesh(2, pmn);

  /*
   * Node = Faces (because 2d)
   * Arc  = edges
   */
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pselect_node  = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;
  _warmup_cell_centered(pm,
                        &pn_node,
                        &pn_arc,
                        &pnode_arc_idx,
                        &pnode_arc,
                        &pnode_weight,
                        &parc_weight);

  // Transpose
  PDM_part_connectivity_transpose(n_part,
                                  pn_node,
                                  pn_arc,
                                  pnode_arc_idx,
                                  pnode_arc,
                                  &parc_node_idx,
                                  &parc_node);

  PDM_part_comm_graph_t *pcg_node = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_FACE,
                                    &pcg_node,
                                    PDM_OWNERSHIP_KEEP);

  PDM_part_comm_graph_t *pcg_arc = NULL;
  PDM_part_mesh_part_comm_graph_get(pm,
                                    PDM_MESH_ENTITY_EDGE,
                                    &pcg_arc,
                                    PDM_OWNERSHIP_KEEP);

  int           n_tot_node     = 0;
  PDM_g_num_t  *gnode_node_idx = NULL;
  PDM_g_num_t  *gnode_node     = NULL;
  int          *gnode_weight   = NULL;
  int          *garc_weight    = NULL;
  PDM_g_num_t  *distrib_node   = NULL;
  int         **part_to_graph  = NULL;
  PDM_part_assembly_dual_graph(test_comm,
                               n_part,
                               pn_node,
                               pn_arc,
                               pselect_node,
                               pnode_arc_idx,
                               pnode_arc,
                               parc_node_idx,
                               parc_node,
                               pnode_weight,
                               parc_weight,
                               pcg_node,
                               pcg_arc,
                               &n_tot_node,
                               &gnode_node_idx,
                               &gnode_node,
                               &gnode_weight,
                               &garc_weight,
                               &distrib_node,
                               &part_to_graph);

  if(0 == 1) {
    PDM_log_trace_array_long(gnode_node_idx,       n_tot_node+1              , "gnode_node_idx ::");
    PDM_log_trace_array_long(gnode_node    , (int) gnode_node_idx[n_tot_node], "gnode_node     ::");
    PDM_log_trace_array_int (gnode_weight  ,       n_tot_node                , "gnode_weight   ::");
    PDM_log_trace_array_int (garc_weight   , (int) gnode_node_idx[n_tot_node], "garc_weight    ::");
    PDM_log_trace_array_long(distrib_node  ,       n_rank+1                  , "distrib_node   ::");
  }

  MPI_CHECK(0, n_tot_node == 9);
  MPI_CHECK(1, n_tot_node == 9);

  PDM_g_num_t expected_p0_gnode_node_idx[10] = {0, 3, 6, 8, 10, 12, 15, 17, 20, 21};
  PDM_g_num_t expected_p0_gnode_node    [21] = {1, 5, 17, 0, 2, 14, 1, 7, 4, 16, 3, 5, 0, 4, 6, 5, 7, 2, 6, 8, 7};
  PDM_g_num_t expected_p0_gnode_weight  [9 ] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p0_garc_weight   [21] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p0_distrib_node  [3 ] = {0, 9, 18};

  PDM_g_num_t expected_p1_gnode_node_idx[10] = {0, 1, 4, 6, 9, 11, 13, 15, 18, 21};
  PDM_g_num_t expected_p1_gnode_node    [21] = {10, 9, 11, 15, 10, 12, 11, 13, 17, 12, 14, 1, 13, 10, 16, 3, 15, 17, 0, 12, 16};
  PDM_g_num_t expected_p1_gnode_weight  [9 ] = {1, 1, 1, 1, 1, 1, 1, 1, 1 };
  PDM_g_num_t expected_p1_garc_weight   [21] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  PDM_g_num_t expected_p1_distrib_node  [3 ] = {0, 9, 18};

  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node_idx, gnode_node_idx, 10 );
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_node    , gnode_node    , 21);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_gnode_weight  , gnode_weight  , 9 );
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_garc_weight   , garc_weight   , 21);
  MPI_CHECK_EQ_C_ARRAY(0, expected_p0_distrib_node  , distrib_node  , 3 );

  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node_idx, gnode_node_idx, 10 );
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_node    , gnode_node    , 21);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_gnode_weight  , gnode_weight  , 9 );
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_garc_weight   , garc_weight   , 21);
  MPI_CHECK_EQ_C_ARRAY(1, expected_p1_distrib_node  , distrib_node  , 3 );


  PDM_free(gnode_node_idx);
  PDM_free(gnode_node    );
  PDM_free(gnode_weight  );
  PDM_free(garc_weight   );
  PDM_free(distrib_node  );

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(part_to_graph[i_part]);
  }
  PDM_free(part_to_graph);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(parc_node_idx[i_part]);
    PDM_free(parc_node    [i_part]);
    PDM_free(pnode_weight [i_part]);
    PDM_free(parc_weight  [i_part]);
  }
  PDM_free(pn_node      );
  PDM_free(pn_arc       );
  PDM_free(pselect_node );
  PDM_free(pnode_arc_idx);
  PDM_free(pnode_arc    );
  PDM_free(parc_node_idx);
  PDM_free(parc_node    );
  PDM_free(pnode_weight );
  PDM_free(parc_weight  );

  PDM_part_mesh_nodal_free(pmn);
  PDM_part_mesh_free(pm);
}


// Same but with pcg_node == NULL (cell centered)

MPI_TEST_CASE("[PDM_transfer_entity1_part_id_to_entity2_part_id] ", 2) {


}

