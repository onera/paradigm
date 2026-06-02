#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_graph_dual.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_sort.h"
#include "pdm_unique.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number vtx in side of mesh A (default : 10).\n\n"
     "  -n_part <level>  Number of partition                        .\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -t               Element kind .\n\n"
     "  -h               This message.\n\n");
  exit (exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int                    argc,
 char                 **argv,
 PDM_g_num_t           *n_vtx_a,
 int                   *post
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_vtx_a = atol(argv[i]);
        *n_vtx_a = (PDM_g_num_t) _n_vtx_a;
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
PDM_part_mesh_nodal_t*
_generate_mesh
(
  PDM_MPI_Comm    pdm_comm,
  int             n_vtx_seg
)
{
  int              n_part       = 1;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_HILBERT;

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

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
  int   argc,
  char *argv[]
)
{
  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_g_num_t n_vtx_a = 10;

  int post = 0;

  _read_args(argc,
             argv,
             &n_vtx_a,
             &post);

  PDM_part_mesh_nodal_t* pmn = _generate_mesh(comm, n_vtx_a);

  if(post) {
    PDM_part_mesh_nodal_dump_vtk(pmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }

  /*
   * Differents modes de splitting
   *   - Split les sommets (idéal pour l'adaptation de maillage )
   *   - Split les cellules (le plus classiques )
   *   - On peut partir d'un pmn ou d'un pm
   */

  /* Warm-up */
  int n_part = PDM_part_mesh_nodal_n_part_get(pmn);
  int  *pn_elt       = NULL;
  int  *pn_vtx       = NULL;
  int **pelt_vtx_idx = NULL;
  int **pelt_vtx     = NULL;
  int **pvtx_vtx_idx = NULL;
  int **pvtx_vtx     = NULL;
  PDM_malloc(pn_elt      , n_part, int  );
  PDM_malloc(pn_vtx      , n_part, int  );
  PDM_malloc(pelt_vtx_idx, n_part, int *);
  PDM_malloc(pelt_vtx    , n_part, int *);
  PDM_malloc(pvtx_vtx_idx, n_part, int *);
  PDM_malloc(pvtx_vtx    , n_part, int *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    pn_vtx[i_part] = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
    PDM_geometry_kind_t principal_kind = PDM_part_mesh_nodal_principal_geom_kind_get(pmn);
    pn_elt[i_part] = PDM_part_mesh_nodal_cell_vtx_connect_get(pmn,
                                                              principal_kind,
                                                              i_part,
                                                              &pelt_vtx_idx[i_part],
                                                              &pelt_vtx    [i_part]);

    // Vtx-vtx connectivity
    int *pvtx_elt_idx = NULL;
    int *pvtx_elt     = NULL;
    PDM_connectivity_transpose(pn_elt       [i_part],
                               pn_vtx       [i_part],
                               pelt_vtx_idx [i_part],
                               pelt_vtx     [i_part],
                               &pvtx_elt_idx,
                               &pvtx_elt);

    PDM_combine_connectivity(pn_vtx       [i_part],
                             pvtx_elt_idx,
                             pvtx_elt,
                             pelt_vtx_idx [i_part],
                             pelt_vtx     [i_part],
                             &pvtx_vtx_idx[i_part],
                             &pvtx_vtx    [i_part]);

    PDM_free(pvtx_elt_idx);
    PDM_free(pvtx_elt    );
  }

  /*
   * Graph assembly
   */
  PDM_part_mesh_nodal_to_part_mesh_t* pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                          PDM_FALSE,
                                                                                          PDM_OWNERSHIP_USER);

  // Connectivities
  int dim = 2;
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

  // Hook
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);


  int  *pn_node       = NULL;
  int  *pn_arc        = NULL;
  int **pselect_node  = NULL;
  int **pnode_arc_idx = NULL;
  int **pnode_arc     = NULL;
  int **parc_node_idx = NULL;
  int **parc_node     = NULL;
  int **pnode_weight  = NULL;
  int **parc_weight   = NULL;

  PDM_malloc(pn_node      , n_part, int  );
  PDM_malloc(pn_arc       , n_part, int  );
  // PDM_malloc(pnode_arc_idx, n_part, int *);
  // PDM_malloc(pnode_arc    , n_part, int *);
  PDM_malloc(parc_node_idx, n_part, int *);
  PDM_malloc(parc_node    , n_part, int *);
  PDM_malloc(pnode_weight , n_part, int *);
  PDM_malloc(parc_weight  , n_part, int *);

  int with_select = 0;
  if(with_select == 1) {
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

      // PDM_log_trace_array_int(pselect_node[i_part], pn_node[i_part], "pselect_node ::");
    }
  }


  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_node[i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VTX );
    pn_arc [i_part] = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &parc_node    [i_part],
                                   &parc_node_idx[i_part],
                                   PDM_OWNERSHIP_KEEP);

    double *vtx_coords = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

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

    // This is stupid but why not
    double vdir[3] = {1., 0., 0.};
    for(int i = 0; i < pn_arc[i_part]; ++i) {
      int i_vtx1 = parc_node[i_part][2*i  ]-1;
      int i_vtx2 = parc_node[i_part][2*i+1]-1;

      double v[3] = {vtx_coords[3*i_vtx2  ] - vtx_coords[3*i_vtx1  ],
                     vtx_coords[3*i_vtx2+1] - vtx_coords[3*i_vtx1+1],
                     vtx_coords[3*i_vtx2+2] - vtx_coords[3*i_vtx1+2]};
      double mod = PDM_MODULE(v);
      v[0] = v[0] / mod;
      v[1] = v[1] / mod;
      v[2] = v[2] / mod;

      double vdot = PDM_DOT_PRODUCT(vdir, v);
      int idot = (int) ( PDM_ABS(vdot) * 10) ;

      if(PDM_ABS(vdot) > 0.9999) {
        idot = 100;
      } else {
        idot = 0;
      }

      // log_trace("i_arc = %i / i_vtx1 = %i / i_vtx2 = %i --> idot = %i (%12.5e) \n ", i, i_vtx1, i_vtx2, idot, vdot);

      parc_weight[i_part][i] += idot;
    }

    // PDM_log_trace_array_int(parc_weight[i_part], pn_arc[i_part], "parc_weight ::");
  }

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

  /*
   * Tester avec cell_vtx + vtx_cell aussi -> Shortcut for mesh adaptation + quality
   */
  int           n_tot_node     = 0;
  PDM_g_num_t  *gnode_node_idx = NULL;
  PDM_g_num_t  *gnode_node     = NULL;
  int          *gnode_weight   = NULL;
  int          *garc_weight    = NULL;
  PDM_g_num_t  *distrib_node   = NULL;
  int         **part_to_graph  = NULL;
  PDM_part_assembly_dual_graph(comm,
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

  // for(int i = 0; i < gnode_node_idx[n_tot_node]; ++i) {
  //   gnode_node[i] -= 1;
  // }

  int *node_part_id = NULL;
  PDM_malloc(node_part_id, n_tot_node, int);
  PDM_para_graph_split(PDM_SPLIT_DUAL_WITH_PTSCOTCH,
                       distrib_node,
                       gnode_node_idx,
                       gnode_node,
                       gnode_weight,
                       garc_weight,
                       n_rank,
                       NULL,
                       node_part_id,
                       comm);

  // PDM_log_trace_array_int(node_part_id, n_tot_node, "node_part_id ::");

  PDM_free(gnode_node_idx);
  PDM_free(gnode_node    );
  PDM_free(gnode_weight  );
  PDM_free(garc_weight   );
  PDM_free(distrib_node  );

  // L'ordre du graphe suit le "select"
  // Attention car on compte pas les is_owner=0 (et ils sont pas dans le graphe)
  // Bon dans tous les cas on essayera de ce rammener au entités principale

  // Synchro color and deconcatenate
  int **vtx_id  = NULL;
  PDM_malloc(vtx_id , n_part, int    *);
  for(int i_part = 0; i_part < n_part; ++i_part) {

    PDM_malloc( vtx_id[i_part], pn_node[i_part], int   );

    for(int i = 0; i < pn_node[i_part]; ++i) {
      int l_node = part_to_graph[i_part][i];
      if(l_node != -1) {
        vtx_id [i_part][i] = node_part_id[l_node];
      } else {
        vtx_id [i_part][i] = -1;
      }
    }
  }

  PDM_part_comm_graph_all_reduce(pcg_node,
                                 PDM_MPI_INT,
                                 1,
                                 PDM_MPI_MAX,
             (unsigned char **)  vtx_id);

  if(0 == 1) { // Usefull to setup unit test easily
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_int(vtx_id[i_part], pn_node[i_part], "vtx_id ::");
      PDM_log_trace_array_int(pelt_vtx_idx[i_part], pn_elt[i_part]+1, "pelt_vtx_idx ::");
      PDM_log_trace_array_int(pelt_vtx[i_part], pelt_vtx_idx[i_part][pn_elt[i_part]], "pelt_vtx ::");
    }
  }

  /*
   * Update delt_id
   */
  int **elt_id = NULL;
  PDM_transfer_entity1_part_id_to_entity2_part_id(n_part,
                                                  vtx_id,
                                                  pn_elt,
                                                  pelt_vtx_idx,
                                                  pelt_vtx,
                                                  &elt_id);



  if(post) {

    double **delt_id = NULL;
    double **dvtx_id = NULL;

    PDM_malloc(delt_id, n_part, double *);
    PDM_malloc(dvtx_id, n_part, double *);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_malloc(delt_id[i_part], pn_elt [i_part], double);
      PDM_malloc(dvtx_id[i_part], pn_node[i_part], double);

      for(int i = 0; i < pn_elt[i_part]; ++i) {
        delt_id[i_part][i] = elt_id[i_part][i];
      }
    }

    for(int i_part = 0; i_part < n_part; ++i_part) {
      for(int i = 0; i < pn_node[i_part]; ++i) {
        dvtx_id[i_part][i] = vtx_id[i_part][i];
      }
    }

    const char    *elt_field_name [] = {"delt_id"};
    double       **elt_field      [] = {delt_id};
    const char    *field_vtx_name [] = {"dvtx_id"};
    double       **field_vtx      [] = {dvtx_id};
    char filename[999];
    sprintf(filename, "repart_pmn");
    PDM_part_mesh_nodal_dump_vtk_with_fields(pmn,
                                             PDM_GEOMETRY_KIND_SURFACIC,
                                             filename,
                                             1,
                                             elt_field_name,
                   (const double ***)        elt_field,
                                             1,
                                             field_vtx_name,
                   (const double ***)        field_vtx);

    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_free(delt_id[i_part]);
      PDM_free(dvtx_id[i_part]);
    }
    PDM_free(delt_id);
    PDM_free(dvtx_id);


  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(vtx_id [i_part]);
    PDM_free(elt_id [i_part]);
  }
  PDM_free(vtx_id);
  PDM_free(elt_id);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(part_to_graph[i_part]);
  }
  PDM_free(part_to_graph);
  PDM_free(node_part_id);


  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnode_arc_idx[i_part]);
    PDM_free(pnode_arc    [i_part]);
    PDM_free(pnode_weight [i_part]);
    PDM_free(parc_weight  [i_part]);
  }

  /* Free all */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pelt_vtx_idx[i_part]);
    PDM_free(pelt_vtx    [i_part]);
    PDM_free(pvtx_vtx_idx[i_part]);
    PDM_free(pvtx_vtx    [i_part]);
  }
  PDM_free(pelt_vtx_idx);
  PDM_free(pelt_vtx    );
  PDM_free(pvtx_vtx_idx);
  PDM_free(pvtx_vtx    );
  PDM_free(pn_elt);
  PDM_free(pn_vtx);

  if(with_select == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_free(pselect_node[i_part]);
    }
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

  PDM_part_mesh_free(pm);
  PDM_part_mesh_nodal_free(pmn);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
