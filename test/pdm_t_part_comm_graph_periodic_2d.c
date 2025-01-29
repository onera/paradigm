#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_comm_graph.h"
#include "pdm_part_extension_algorithm.h"
#include "pdm_printf.h"
#include "pdm_priv.h"
#include "pdm_vtk.h"

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
  PDM_printf(
    "\n"
    "  Usage: \n\n"
    "  -n           <int>  Number vtx along each side (default: 10).\n\n"
    "  -n_part      <int>  Number of partitions per MPI rank (default: 1).\n\n"
    "  -part_method <int>  Partitioning method (default: 3 (Hilbert)).\n\n"
    "  -pj                 Perodic in y-direction.\n\n"
    "  -rotation           Perodic by rotation.\n\n"
    "  -visu               Enable output for visualization.\n\n"
    "  -h                  This message.\n\n");
  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc         Number of arguments
 * \param [in]    argv         Arguments
 * \param [inout] n_vtx_seg    Number of vtx along each side of the domain
 * \param [inout] n_part       Number of partitions per MPI rank
 * \param [inout] part_method  Partitioning method
 * \param [inout] periodic_j   Perodic in y-direction
 * \param [inout] rotation     Perodic by rotation
 * \param [inout] visu         Enable output for visualization
 *
 */
static void
_read_args
(
  int                    argc,
  char                 **argv,
  PDM_g_num_t           *n_vtx_seg,
  int                   *n_part,
  PDM_split_dual_t      *part_method,
  int                   *periodic_j,
  PDM_bool_t            *rotation,
  PDM_bool_t            *visu
)
{
  int i = 1;

  /* Parse and check command line */
  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long n = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) n;
      }
    }

    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-part_method") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *part_method = (PDM_split_dual_t) atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-pj") == 0) {
      *periodic_j = 1;
    }

    else if (strcmp(argv[i], "-rotation") == 0) {
      *rotation = PDM_TRUE;
    }

    else if (strcmp(argv[i], "-visu") == 0) {
      *visu = PDM_TRUE;
    }

    else {
      _usage(EXIT_FAILURE);
    }

    i++;
  }
}



/**
 *
 * \brief  Main
 *
 */
int main
(
  int   argc,
  char *argv[]
)
{
  /* Command line arguments */
  PDM_g_num_t      n_vtx_seg   = 10;
  int              n_part      = 1;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  int              periodic_j  = 0;
  PDM_bool_t       rotation    = PDM_FALSE;
  PDM_bool_t       visu        = PDM_FALSE;
  _read_args(argc,
             argv,
             &n_vtx_seg,
             &n_part,
             &part_method,
             &periodic_j,
             &rotation,
             &visu);
             
  /* Initialize MPI */
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Generate block-distributed mesh with periodic joints */
  int n_domain_i = 1;
  int n_domain_j = 1;
  int n_domain_k = 1;

  int n_domain   = n_domain_i*n_domain_j*n_domain_k;

  int periodic_i = 1;
  int periodic_k = 0;

  PDM_dcube_nodal_t      **dcube    = NULL;
  PDM_domain_interface_t  *dom_itrf = NULL;

  PDM_dcube_nodal_cart_topo(comm,
                            n_domain_i,
                            n_domain_j,
                            n_domain_k,
                            periodic_i,
                            periodic_j,
                            periodic_k,
                            n_vtx_seg,
                            n_vtx_seg,
                            1,
                            2.,
                            -1.,
                            -1.,
                            0.,
                            PDM_MESH_NODAL_TRIA3,
                            1,
                            &dcube,
                            &dom_itrf,
                            PDM_OWNERSHIP_KEEP);

  /* Partitioning */
  int *n_part_per_domain = PDM_array_const_int(n_domain, n_part);

  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                n_part_per_domain,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i_dom]);
    PDM_dmesh_nodal_generate_distribution(dmn);
    PDM_multipart_dmesh_nodal_set(mpart, i_dom, dmn);
  }

  PDM_multipart_compute(mpart);


  /* Get partitioned mesh */
  PDM_part_mesh_nodal_t **pmn = NULL;
  PDM_malloc(pmn, n_domain, PDM_part_mesh_nodal_t *);
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn[i_dom],
                                      PDM_OWNERSHIP_USER);
  }

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_dcube_nodal_gen_free(dcube[i_dom]);
  }
  PDM_free(dcube);


  /* Get partitioned domain interface */
  int n_interface = PDM_domain_interface_n_interface_get(dom_itrf);

  PDM_part_domain_interface_t *pdom_itrf = PDM_part_domain_interface_create(n_interface,
                                                                            n_domain,
                                                                            n_part_per_domain,
                                                                            (n_domain > 1),
                                                                            PDM_OWNERSHIP_KEEP,
                                                                            comm);
  
  int          *itrf_dn  = NULL;
  PDM_g_num_t **itrf_ids = NULL;
  int         **itrf_dom = NULL;  
  PDM_domain_interface_get(dom_itrf,
                           PDM_BOUND_TYPE_VTX,
                           &itrf_dn,
                           &itrf_ids,
                           &itrf_dom);

  int          **n_vtx        = NULL;
  PDM_g_num_t ***vtx_ln_to_gn = NULL;
  PDM_malloc(n_vtx,        n_domain, int          *);
  PDM_malloc(vtx_ln_to_gn, n_domain, PDM_g_num_t **);
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_malloc(n_vtx       [i_dom], n_part_per_domain[i_dom], int          );
    PDM_malloc(vtx_ln_to_gn[i_dom], n_part_per_domain[i_dom], PDM_g_num_t *);
  
    for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {
      n_vtx       [i_dom][i_part] = PDM_part_mesh_nodal_n_vtx_get    (pmn[i_dom], i_part);
      vtx_ln_to_gn[i_dom][i_part] = PDM_part_mesh_nodal_vtx_g_num_get(pmn[i_dom], i_part, PDM_OWNERSHIP_KEEP);
    }
  }

  PDM_ddomain_interface_to_pdomain_interface(comm,
                                             n_interface,
                                             n_domain,
                                             (n_domain > 1),
                                             PDM_BOUND_TYPE_VTX,
                                             itrf_dn,
                                             itrf_ids,
                                             itrf_dom,
                                             n_part_per_domain,
                                             n_vtx,
                                             vtx_ln_to_gn,
                                             pdom_itrf);
  
  PDM_domain_interface_free(dom_itrf);

  /* Get vtx_part_bound */
  int ln_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    ln_part += n_part_per_domain[i_dom];
  }

  int gn_part = 0;
  PDM_MPI_Allreduce(&ln_part, &gn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int ***vtx_part_bound_proc_idx = NULL;
  int ***vtx_part_bound_part_idx = NULL;
  int ***vtx_part_bound          = NULL;
  PDM_malloc(vtx_part_bound_proc_idx, n_domain, int **);
  PDM_malloc(vtx_part_bound_part_idx, n_domain, int **);
  PDM_malloc(vtx_part_bound,          n_domain, int **);
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_malloc(vtx_part_bound_proc_idx[i_dom], n_part_per_domain[i_dom], int *);
    PDM_malloc(vtx_part_bound_part_idx[i_dom], n_part_per_domain[i_dom], int *);
    PDM_malloc(vtx_part_bound         [i_dom], n_part_per_domain[i_dom], int *);
    for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &vtx_part_bound_proc_idx[i_dom][i_part],
                                        &vtx_part_bound_part_idx[i_dom][i_part],
                                        &vtx_part_bound         [i_dom][i_part],
                                        PDM_OWNERSHIP_KEEP);
    }
  }

  /* Make unified graph */
  int **unified_vtx_graph_idx   = NULL;
  int **unified_vtx_graph_trplt = NULL;
  int **unified_vtx_graph_itrf  = NULL;
  PDM_part_extension_build_entity1_graph(pdom_itrf,
                                         PDM_BOUND_TYPE_VTX,
                                         n_domain,
                                         n_part_per_domain,
                                         n_vtx,
                                         vtx_ln_to_gn,
                                         vtx_part_bound_proc_idx,
                                         vtx_part_bound_part_idx,
                                         vtx_part_bound,
                                         NULL,
                                         1,
                                         &unified_vtx_graph_idx,
                                         &unified_vtx_graph_trplt,
                                         &unified_vtx_graph_itrf,
                                         comm);

  /* Create Part Comm Graph */
  int  *n_unified_vtx_graph = NULL;
  int **unified_vtx_graph   = NULL;
  PDM_malloc(n_unified_vtx_graph, ln_part, int  );
  PDM_malloc(unified_vtx_graph,   ln_part, int *);
  int j_part = 0; 
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {

      for (int i_vtx = 0; i_vtx < n_vtx[i_dom][i_part]; i_vtx++) {
        unified_vtx_graph_idx[j_part][i_vtx+1] /= 3;
      }
      n_unified_vtx_graph[j_part] = unified_vtx_graph_idx[j_part][n_vtx[i_dom][i_part]];

      // Triplet to quadruplet
      PDM_malloc(unified_vtx_graph[j_part], n_unified_vtx_graph[j_part] * 4, int);
      int i = 0;
      for (int i_vtx = 0; i_vtx < n_vtx[i_dom][i_part]; i_vtx++) {
        for (int idx = unified_vtx_graph_idx[j_part][i_vtx]; idx < unified_vtx_graph_idx[j_part][i_vtx+1]; idx++) {
          unified_vtx_graph[j_part][i++] = i_vtx + 1;
          unified_vtx_graph[j_part][i++] = unified_vtx_graph_trplt[j_part][3*idx  ];
          unified_vtx_graph[j_part][i++] = unified_vtx_graph_trplt[j_part][3*idx+1] + 1;
          unified_vtx_graph[j_part][i++] = unified_vtx_graph_trplt[j_part][3*idx+2] + 1;
        }
      }

      if (visu) {
        log_trace("j_part %d, unified graph in quadruplets :\n", j_part);
        for (int i = 0; i < n_unified_vtx_graph[j_part]; i++) {
          log_trace("%4d %2d %2d %4d  itrf %d\n", 
                    unified_vtx_graph[j_part][4*i  ],
                    unified_vtx_graph[j_part][4*i+1],
                    unified_vtx_graph[j_part][4*i+2],
                    unified_vtx_graph[j_part][4*i+3],
                    unified_vtx_graph_itrf[j_part][i]);
        }
      }

      j_part++;
    }
  }

  PDM_part_comm_graph_t *pcg_vtx = PDM_part_comm_graph_with_nuplet_create(ln_part,
                                                                          n_unified_vtx_graph,
                                                                          unified_vtx_graph,
                                                                          PDM_OWNERSHIP_USER,
                                                                          1,
                                                                          unified_vtx_graph_itrf,
                                                                          PDM_OWNERSHIP_USER,
                                                                          PDM_TRUE,
                                                                          comm);

  if (visu) {
    /* Visualize mesh */
    for (int i_dom = 0; i_dom < n_domain; i_dom++) {
      char prefix[999];
      sprintf(prefix, "periodic_2d_pmn_domain_%d", i_dom);
      PDM_part_mesh_nodal_dump_vtk(pmn[i_dom],
                                   PDM_GEOMETRY_KIND_SURFACIC,
                                   prefix);

      sprintf(prefix, "periodic_2d_pmn_ridge_domain_%d", i_dom);
      PDM_part_mesh_nodal_dump_vtk(pmn[i_dom],
                                   PDM_GEOMETRY_KIND_RIDGE,
                                   prefix);
    }
  }


  /* Deduce part comm graph of edges */
  int  *pn_vtx        = NULL;
  int  *pn_edge       = NULL;
  int **pedge_vtx_idx = NULL;
  int **pedge_vtx     = NULL;
  PDM_malloc(pn_vtx,        ln_part, int  );
  PDM_malloc(pn_edge,       ln_part, int  );
  PDM_malloc(pedge_vtx_idx, ln_part, int *);
  PDM_malloc(pedge_vtx,     ln_part, int *);

  j_part = 0;
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {
      pn_vtx [j_part] = n_vtx[i_dom][i_part];
      pn_edge[j_part] = PDM_part_mesh_nodal_cell_vtx_connect_get(pmn[i_dom],
                                                                 PDM_GEOMETRY_KIND_RIDGE,
                                                                 i_part,
                                                                 &pedge_vtx_idx[j_part],
                                                                 &pedge_vtx    [j_part]);
      
      j_part++;
    }
  }

  int  *n_unified_edge_graph = NULL;
  int **unified_edge_graph   = NULL;
  int **unified_edge_itrf    = NULL;
  PDM_part_comm_graph_entity1_to_entity2(comm,
                                         ln_part,
                                         n_unified_vtx_graph,
                                         unified_vtx_graph,
                                         1,
                                         unified_vtx_graph_itrf,
                                         pn_vtx,
                                         pn_edge,
                                         pedge_vtx_idx,
                                         pedge_vtx,
                                         &n_unified_edge_graph,
                                         &unified_edge_graph,
                                         &unified_edge_itrf);

  if (visu) {
    j_part = 0; 
    for (int i_dom = 0; i_dom < n_domain; i_dom++) {
      for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {
        log_trace("j_part %d, unified edge graph in quadruplets :\n", j_part);
        for (int i = 0; i < n_unified_edge_graph[j_part]; i++) {
          log_trace("%4d %2d %2d %4d through interface %d\n", 
                    unified_edge_graph[j_part][4*i  ],
                    unified_edge_graph[j_part][4*i+1],
                    unified_edge_graph[j_part][4*i+2],
                    unified_edge_graph[j_part][4*i+3],
                    unified_edge_itrf [j_part][  i  ]);
        }

        j_part++;
      }
    }
  }

  PDM_part_comm_graph_t *pcg_edge = PDM_part_comm_graph_with_nuplet_create(ln_part,
                                                                           n_unified_edge_graph,
                                                                           unified_edge_graph,
                                                                           PDM_OWNERSHIP_USER,
                                                                           1,
                                                                           unified_edge_itrf,
                                                                           PDM_OWNERSHIP_USER,
                                                                           PDM_TRUE,
                                                                           comm);

  if (visu) {
    double **send_data = NULL;
    double **recv_data = NULL;
    PDM_malloc(send_data, ln_part, double *);

    j_part = 0; 
    for (int i_dom = 0; i_dom < n_domain; i_dom++) {
      for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {
        
        double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn[i_dom], i_part, PDM_OWNERSHIP_KEEP);

        PDM_malloc(send_data[j_part], n_unified_edge_graph[j_part] * 3 * 2, double);
        for (int i_bnd = 0; i_bnd < n_unified_edge_graph[j_part]; i_bnd++) {
          int i_edge = unified_edge_graph[j_part][4*i_bnd] - 1;
          int i_itrf = PDM_ABS(unified_edge_itrf[j_part][i_bnd]);

          for (int idx_vtx = 0; idx_vtx < 2; idx_vtx++) {
            double *p = &send_data[j_part][6*i_bnd + 3*idx_vtx];
            int i_vtx = pedge_vtx[j_part][2*i_edge+idx_vtx] - 1;
            memcpy(p, &vtx_coord[3*i_vtx], sizeof(double) * 3);

            if (i_itrf == 1) {
              p[0] -= 2*PDM_SIGN(unified_edge_itrf[j_part][i_bnd]);
            }
            if (i_itrf == 2) {
              p[1] -= 2*PDM_SIGN(unified_edge_itrf[j_part][i_bnd]);
            }
          }
        }

        j_part++;
      }
    }

    PDM_part_comm_graph_exch(pcg_edge, 
                             sizeof(double) * 3 * 2,
                             PDM_STRIDE_CST_INTERLACED,
                             1,
                             NULL,
                  (void  **) send_data,
                             NULL,
                  (void ***) &recv_data);

    j_part = 0; 
    for (int i_dom = 0; i_dom < n_domain; i_dom++) {
      for (int i_part = 0; i_part < n_part_per_domain[i_dom]; i_part++) {

        char name[999];
        sprintf(name, "check_perio_edge_%d_%d.vtk", j_part, i_rank);

        const int *is_owner = PDM_part_comm_graph_owner_get(pcg_edge, j_part);

        int *connec   = NULL;
        int *opp_rank = NULL;
        int *opp_part = NULL;
        int *opp_lnum = NULL;
        PDM_malloc(connec,   2*n_unified_edge_graph[j_part], int);
        PDM_malloc(opp_rank,   n_unified_edge_graph[j_part], int);
        PDM_malloc(opp_part,   n_unified_edge_graph[j_part], int);
        PDM_malloc(opp_lnum,   n_unified_edge_graph[j_part], int);
        for (int i = 0; i < n_unified_edge_graph[j_part]; i++) {
          connec[2*i  ] = 1 + 2*i;
          connec[2*i+1] = 1 + 2*i+1;

          opp_rank[i] = unified_edge_graph[j_part][4*i+1];
          opp_part[i] = unified_edge_graph[j_part][4*i+2];
          opp_lnum[i] = unified_edge_graph[j_part][4*i+3];
        }

        const char *elt_field_name [] = {"itrf", "opp_rank", "opp_part", "opp_lnum", "is_owner"};
        const int  *elt_field_value[] = {unified_edge_itrf[j_part], opp_rank, opp_part, opp_lnum, is_owner};

        PDM_vtk_write_std_elements(name, 
                                   2*n_unified_edge_graph[j_part],
                                   recv_data[j_part],
                                   NULL,
                                   PDM_MESH_NODAL_BAR2,
                                   n_unified_edge_graph[j_part],
                                   connec,
                                   NULL,
                                   5,
                                   elt_field_name,
                                   elt_field_value);
        PDM_free(connec);
        PDM_free(opp_rank);
        PDM_free(opp_part);
        PDM_free(opp_lnum);

        PDM_free(send_data[j_part]);
        PDM_free(recv_data[j_part]);
        j_part++;
      }
    }
    PDM_free(send_data);
    PDM_free(recv_data);
  }


  /* Free memory */
  PDM_part_comm_graph_free(pcg_vtx);
  PDM_part_comm_graph_free(pcg_edge);
  PDM_multipart_free(mpart);
  PDM_part_domain_interface_free(pdom_itrf);
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_part_mesh_nodal_free(pmn[i_dom]);
    PDM_free(n_vtx       [i_dom]);
    PDM_free(vtx_ln_to_gn[i_dom]);

  }
  for (int i_part = 0; i_part < ln_part; i_part++) {
    PDM_free(unified_vtx_graph_idx  [i_part]);
    PDM_free(unified_vtx_graph_trplt[i_part]);
    PDM_free(unified_vtx_graph_itrf [i_part]);

    PDM_free(unified_vtx_graph[i_part]);

    PDM_free(pedge_vtx_idx[i_part]);
    PDM_free(pedge_vtx    [i_part]);

    PDM_free(unified_edge_graph[i_part]);
    PDM_free(unified_edge_itrf [i_part]);
  }
  PDM_free(n_unified_vtx_graph);
  PDM_free(unified_vtx_graph);
  PDM_free(pn_vtx);
  PDM_free(pn_edge);
  PDM_free(pedge_vtx_idx);
  PDM_free(pedge_vtx);
  PDM_free(n_unified_edge_graph);
  PDM_free(unified_edge_graph);
  PDM_free(unified_edge_itrf);

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    PDM_free(vtx_part_bound_proc_idx[i_dom]);
    PDM_free(vtx_part_bound_part_idx[i_dom]);
    PDM_free(vtx_part_bound         [i_dom]);
  }
  PDM_free(vtx_part_bound_proc_idx);
  PDM_free(vtx_part_bound_part_idx);
  PDM_free(vtx_part_bound);

  PDM_free(pmn);
  PDM_free(n_vtx);
  PDM_free(vtx_ln_to_gn);
  PDM_free(unified_vtx_graph_idx);
  PDM_free(unified_vtx_graph_trplt);
  PDM_free(unified_vtx_graph_itrf);
  PDM_free(n_part_per_domain);

  /* Finalize MPI */
  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("The End\n");
  }

  return EXIT_SUCCESS;
}