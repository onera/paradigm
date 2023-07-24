#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_multipart.h"
#include "pdm_array.h"
#include "pdm_part_to_part.h"
#include "pdm_extract_part.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_part_connectivity_transform.h"

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
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
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
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_g_pts,
 int           *n_dom_i,
 int           *periodic_i,
 int           *n_depth
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
        long _n_g_pts = atol(argv[i]);
        *n_g_pts = (PDM_g_num_t) _n_g_pts;
      }
    }
    else if (strcmp(argv[i], "-ni") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_i = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-depth") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_depth = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
PDM_g_num_t*
_compute_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_MPI_Comm     comm
)
{

  PDM_g_num_t *shift_by_domain_loc = PDM_array_const_gnum(n_domain, 0);
  PDM_g_num_t *shift_by_domain     = (PDM_g_num_t *) malloc((n_domain+1) * sizeof(PDM_g_num_t));

  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {

      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        shift_by_domain_loc[i_domain] = PDM_MAX(shift_by_domain_loc[i_domain], _pentity_ln_to_gn[i]);
      }
    }
  }

  shift_by_domain[0] = 0;
  PDM_MPI_Allreduce(shift_by_domain_loc, &shift_by_domain[1], n_domain, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_array_accumulate_gnum(shift_by_domain, n_domain+1);

  free(shift_by_domain_loc);

  return shift_by_domain;
}

// static
// void
// _shift_ln_to_gn
// (
//   int          n_entity,
//   PDM_g_num_t *entity_ln_to_gn,
//   PDM_g_num_t  shift,
//   int          sens
// )
// {
//   for(int i = 0; i < n_entity; ++i) {
//     entity_ln_to_gn[i] = entity_ln_to_gn[i] + sens * shift;
//   }
// }

static
void
_offset_ln_to_gn_by_domain
(
  int              n_domain,
  int             *n_part,
  int            **pn_entity,
  PDM_g_num_t   ***pentity_ln_to_gn,
  PDM_g_num_t     *shift_by_domain,
  int              sens
)
{
  for(int i_domain = 0; i_domain < n_domain; ++i_domain) {
    for(int i_part = 0; i_part < n_part[i_domain]; ++i_part) {
      int          _pn_entity        = pn_entity       [i_domain][i_part];
      PDM_g_num_t *_pentity_ln_to_gn = pentity_ln_to_gn[i_domain][i_part];
      for(int i = 0; i < _pn_entity; ++i) {
        _pentity_ln_to_gn[i] = _pentity_ln_to_gn[i] + sens * shift_by_domain[i_domain];
      }
    }
  }
}


static
void
_part_extension
(
  int                          n_depth,
  int                          n_domain,
  int*                         n_part,
  PDM_MPI_Comm                 comm,
  PDM_multipart_t             *mpart,
  PDM_part_domain_interface_t *pdi
)
{
  PDM_UNUSED(n_depth);
  PDM_UNUSED(mpart);
  PDM_UNUSED(pdi);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int* n_part_g = malloc(n_domain * sizeof(int));
  PDM_MPI_Allreduce(n_part, n_part_g, n_domain, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int          **pn_vtx             = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int          *pflat_n_vtx         = (int           *) malloc( n_domain * sizeof(int           ));
  int          *pflat_n_edge        = (int           *) malloc( n_domain * sizeof(int           ));
  int         **pflat_edge_vtx_idx  = (int          **) malloc( n_domain * sizeof(int          *));
  int         **pflat_edge_vtx      = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      pflat_n_vtx       [ln_part_tot+i_part] = pn_vtx       [i_dom][i_part];
      pflat_vtx_ln_to_gn[ln_part_tot+i_part] = pvtx_ln_to_gn[i_dom][i_part];

      pflat_n_edge      [ln_part_tot+i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                                               i_dom,
                                                                               i_part,
                                                                               PDM_MESH_ENTITY_EDGE,
                                                                               &pflat_edge_ln_to_gn[ln_part_tot+i_part],
                                                                               PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &pflat_edge_vtx    [ln_part_tot+i_part],
                                          &pflat_edge_vtx_idx[ln_part_tot+i_part],
                                          PDM_OWNERSHIP_KEEP);

      assert(pflat_edge_vtx_idx[ln_part_tot+i_part] == NULL);
      pflat_edge_vtx_idx[ln_part_tot+i_part] = malloc((pflat_n_edge[ln_part_tot+i_part] + 1) * sizeof(int));
      for(int i_edge = 0; i_edge < pflat_n_edge[ln_part_tot+i_part]+1; ++i_edge) {
        pflat_edge_vtx_idx[ln_part_tot+i_part][i_edge] = 2 * i_edge;
      }


    }
    ln_part_tot += n_part[i_dom];
  }

  int         **pflat_vtx_edge_idx = NULL;
  int         **pflat_vtx_edge     = NULL;
  PDM_part_connectivity_transpose(ln_part_tot,
                                  pflat_n_edge,
                                  pflat_n_vtx,
                                  pflat_edge_vtx_idx,
                                  pflat_edge_vtx,
                                  &pflat_vtx_edge_idx,
                                  &pflat_vtx_edge);

  PDM_g_num_t* shift_by_domain_vtx = _compute_offset_ln_to_gn_by_domain(n_domain,
                                                                        n_part,
                                                                        pn_vtx,
                                                                        pvtx_ln_to_gn,
                                                                        comm);

  /* Shift ln_to_gn */
  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                             shift_by_domain_vtx,
                             1);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      int *ppart_vtx_proc_idx = NULL;
      int *ppart_vtx_part_idx = NULL;
      int *ppart_vtx          = NULL; // (i_vtx, i_proc, i_part, i_vtx_opp)
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_BOUND_TYPE_VTX,
                                        &ppart_vtx_proc_idx,
                                        &ppart_vtx_part_idx,
                                        &ppart_vtx,
                                        PDM_OWNERSHIP_KEEP);

      // PDM_log_trace_array_int(ppart_vtx, 4 * ppart_vtx_part_idx[n_part_g[i_dom]], "ppart_vtx");

    }
  }

  //
  int         **part1_to_part2_idx         = malloc(ln_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet_idx = NULL; //malloc(ln_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet     = malloc(ln_part_tot * sizeof(int *));

  // Count
  int li_part = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {
      int *ppart_vtx_proc_idx = NULL;
      int *ppart_vtx_part_idx = NULL;
      int *ppart_vtx          = NULL; // (i_vtx, i_proc, i_part, i_vtx_opp)
      PDM_multipart_part_graph_comm_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_BOUND_TYPE_VTX,
                                        &ppart_vtx_proc_idx,
                                        &ppart_vtx_part_idx,
                                        &ppart_vtx,
                                        PDM_OWNERSHIP_KEEP);

      PDM_log_trace_array_int(ppart_vtx, 4 * ppart_vtx_part_idx[n_part_g[i_dom]], "ppart_vtx");

      part1_to_part2_idx[li_part] = malloc((pflat_n_vtx[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pflat_n_vtx[li_part]);

      int n_entity_bound = ppart_vtx_part_idx[n_part_g[i_dom]];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_vtx[4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      PDM_log_trace_array_int(part1_to_part2_n, pflat_n_vtx[li_part], "part1_to_part2_n ::");

      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]];
      part1_to_part2_triplet[li_part] = malloc(n_connect_tot * sizeof(int));

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_vtx[4*idx_entity]-1;
        int i_proc_opp   = ppart_vtx[4*idx_entity+1];
        int i_part_opp   = ppart_vtx[4*idx_entity+2]-1; // A gere -> Le shift
        int i_entity_opp = ppart_vtx[4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + part1_to_part2_n[i_entity]++;
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp;
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      PDM_log_trace_array_int(part1_to_part2_idx    [li_part],     pflat_n_vtx[li_part], "part1_to_part2_idx ::");
      PDM_log_trace_array_int(part1_to_part2_triplet[li_part],     n_connect_tot       , "part1_to_part2_triplet ::");

      free(part1_to_part2_n);
      li_part += 1;
    }
  }


  // Attention à l'histoire du x3 dans l'idx
  PDM_part_to_part_t* ptp = PDM_part_to_part_create_from_num2_triplet((const PDM_g_num_t **) pflat_vtx_ln_to_gn,
                                                                      (const int          *) pflat_n_vtx,
                                                                      ln_part_tot,
                                                                      (const int          *) pflat_n_vtx,
                                                                      ln_part_tot,
                                                                      (const int         **) part1_to_part2_idx,
                                                                      (const int         **) part1_to_part2_triplet_idx,
                                                                      (const int         **) part1_to_part2_triplet,
                                                                      comm);


  int  *n_ref_lnum2 = NULL;
  int **ref_lnum2   = NULL;
  PDM_part_to_part_ref_lnum2_get(ptp, &n_ref_lnum2, &ref_lnum2);

  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");
    }
  }

  PDM_part_to_part_free(ptp);

  int          *pextract_n_edge                = NULL;
  int         **pextract_vtx_edge_idx          = NULL;
  int         **pextract_vtx_edge              = NULL;
  PDM_g_num_t **pextract_edge_parent_ln_to_gn  = NULL;
  int         **pextract_edge_to_edge_location = NULL;
  PDM_part_to_part_t* ptp_vtx_edge = NULL;


  PDM_pconnectivity_to_pconnectivity_from_location_keep(comm,
                                                        ln_part_tot,
                               (const int            *) pflat_n_vtx,
                               (const int           **) pflat_vtx_edge_idx,
                               (const int           **) pflat_vtx_edge,
                               (const PDM_g_num_t   **) pflat_edge_ln_to_gn,
                                                        ln_part_tot,
                               (const int            *) pflat_n_vtx,
                               (const PDM_g_num_t   **) pflat_vtx_ln_to_gn,
                               (const int           **) part1_to_part2_idx,
                               (const int           **) part1_to_part2_triplet,
                                                        &pextract_n_edge,
                                                        &pextract_vtx_edge_idx,
                                                        &pextract_vtx_edge,
                                                        &pextract_edge_parent_ln_to_gn,
                                                        &pextract_edge_to_edge_location,
                                                        &ptp_vtx_edge);

  if(1 == 1) {
    PDM_log_trace_array_int(pextract_n_edge, ln_part_tot, "pextract_n_edge");
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_long(pextract_edge_parent_ln_to_gn [i_part],     pextract_n_edge[i_part], "pextract_edge_parent_ln_to_gn");
      PDM_log_trace_array_int (pextract_edge_to_edge_location[i_part], 3 * pextract_n_edge[i_part], "pextract_edge_to_edge_location");
    }
  }


  PDM_part_to_part_free(ptp_vtx_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_vtx_edge_idx         [i_part]);
    free(pextract_vtx_edge             [i_part]);
    free(pextract_edge_parent_ln_to_gn [i_part]);
    free(pextract_edge_to_edge_location[i_part]);
  }
  free(pextract_n_edge);
  free(pextract_vtx_edge_idx);
  free(pextract_vtx_edge);
  free(pextract_edge_parent_ln_to_gn);
  free(pextract_edge_to_edge_location);

  // Objectif :
  //   - Faire un part_to_part qui englobe les raccords entre domain et entre partition
  //   - ATTENTION au piege du shift

  // A faire :
  //  - Shifter tout les ln_to_gn

  // On échange le numero de domain pour retrouver le lien partition -> domaine opposé !
  // On échange le kind aussi pour dissocié périodicité/interface entre domaine/interface entre partition

  /* Unshift ln_to_gn */
  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                             shift_by_domain_vtx,
                             -1);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(part1_to_part2_idx        [i_part]);
    // free(part1_to_part2_triplet_idx[i_part]);
    free(part1_to_part2_triplet    [i_part]);

    free(pflat_vtx_edge_idx[i_part] );
    free(pflat_vtx_edge    [i_part] );
    free(pflat_edge_vtx_idx[i_part] );
  }
  free(part1_to_part2_idx        );
  // free(part1_to_part2_triplet_idx);
  free(part1_to_part2_triplet    );

  free(shift_by_domain_vtx);
  free(n_part_g);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    free(pn_vtx       [i_dom]);
    free(pvtx_ln_to_gn[i_dom]);
  }
  free(pn_vtx);
  free(pvtx_ln_to_gn);
  free(pflat_n_vtx       );
  free(pflat_n_edge      );
  free(pflat_vtx_ln_to_gn);
  free(pflat_edge_ln_to_gn);
  free(pflat_edge_vtx_idx );
  free(pflat_edge_vtx     );
  free(pflat_vtx_edge_idx );
  free(pflat_vtx_edge     );
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
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t n_g_pts   = 10;
  int         n_dom_i    = 1;
  int         periodic_i = 0;
  int         n_depth    = 1;
  _read_args(argc,
             argv,
             &n_g_pts,
             &n_dom_i,
             &periodic_i,
             &n_depth);

  double      **dvtx_coord   = NULL;
  PDM_g_num_t **distrib_edge = NULL;
  PDM_g_num_t **distrib_vtx  = NULL;
  PDM_g_num_t **dedge_vtx    = NULL;

  PDM_domain_interface_t *dom_itrf = NULL;
  PDM_generate_cart_topo_lines(comm,
                               n_dom_i,
                               periodic_i,
                               0.,
                               0.,
                               0.,
                               1.,
                               n_g_pts,
                               &distrib_edge,
                               &distrib_vtx,
                               &dedge_vtx,
                               &dvtx_coord,
                               &dom_itrf);

  PDM_dmesh_t **dm = malloc(n_dom_i * sizeof(PDM_dmesh_t *));

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {

    int dn_edge = distrib_edge[i_dom][i_rank+1] - distrib_edge[i_dom][i_rank];
    int dn_vtx  = distrib_vtx [i_dom][i_rank+1] - distrib_vtx [i_dom][i_rank];

    /*
     * Create dmesh
     */
    dm[i_dom] = PDM_dmesh_create(PDM_OWNERSHIP_KEEP,
                                 0,
                                 0,
                                 dn_edge,
                                 dn_vtx,
                                 comm);

    PDM_dmesh_connectivity_set(dm[i_dom],
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               dedge_vtx[i_dom],
                               NULL,
                               PDM_OWNERSHIP_USER);

    PDM_dmesh_vtx_coord_set(dm[i_dom],
                            dvtx_coord[i_dom],
                            PDM_OWNERSHIP_USER);
  }


  /*
   * Mulitpart
   */
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
  // PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
  int* n_part = malloc(n_dom_i * sizeof(int));
  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
    n_part[i_dom] = 1;
  }
  PDM_multipart_t* mpart = PDM_multipart_create(n_dom_i,
                                                n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
    PDM_multipart_register_block(mpart, i_dom, dm[i_dom]);
  }

  PDM_multipart_run_ppart(mpart);

  int n_domain = n_dom_i;

  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

    for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);



    }
  }

  if(0 == 1) {
    for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {
      for(int i_part = 0; i_part < n_part[i_dom]; ++i_part) {

        PDM_g_num_t* lpvtx_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        i_dom,
                                        i_part,
                                        PDM_MESH_ENTITY_VERTEX,
                                        &lpvtx_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t* pedge_ln_to_gn = NULL;
        int pn_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                      i_dom,
                                                      i_part,
                                                      PDM_MESH_ENTITY_EDGE,
                                                      &pedge_ln_to_gn,
                                                      PDM_OWNERSHIP_KEEP);

        int *pedge_vtx     = NULL;
        int *pedge_vtx_idx = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            i_dom,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                            &pedge_vtx,
                                            &pedge_vtx_idx,
                                            PDM_OWNERSHIP_KEEP);

        double *pvtx_coord = NULL;
        int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     i_dom,
                                                     i_part,
                                                     &pvtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
        char filename[999];
        sprintf(filename, "out_part_vtx_i_dom=%i_%i.vtk", i_dom, i_rank);
        PDM_vtk_write_std_elements (filename,
                                    n_vtx,
                                    pvtx_coord,
                                    lpvtx_ln_to_gn,
                                    PDM_MESH_NODAL_BAR2,
                                    pn_edge,
                                    pedge_vtx,
                                    pedge_ln_to_gn,
                                    0,
                                    NULL,
                                    NULL);

      }
    }
  }

  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_itrf,
                                                                                   n_part,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pn_vtx,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pvtx_ln_to_gn);

  PDM_part_domain_interface_free(pdi);
  PDM_domain_interface_free(dom_itrf);


  /*
   * Extention de partition
   */

  /*
   *  - Step 1 : n_dom = 1, n_proc = 2 + rank1
   *  - Step 2 : n_dom = 1, n_proc = 2 + rank2
   *  - Step 3 : n_dom = 1, n_proc = 1 + rank1 + perio
   *  - Step 4 : n_dom = 1, n_proc = 1 + rank2 + perio
   *  - Step 5 : n_dom = 1, n_proc = 1 + rank1 + perio mais qui se recouvre plusieurs fois (genre 2 cellules)
   *  - Step 6 : n_dom = 2, n_proc = 1 + rank1
   *  - Step 7 : n_dom = 2, n_proc = 2 + rank1
   *  - Step 8 : n_dom = 2, n_proc = 2 + rank1 + perio
   *  - Step 9 : n_dom = 2, n_proc = 2 + rank1 + perio mais qui se recouvre plusieurs fois (genre 2 cellules)
   *
   */
  _part_extension(n_depth,
                  n_dom_i,
                  n_part,
                  comm,
                  mpart,
                  pdi);



  PDM_multipart_free(mpart);
  free(n_part);

  for(int i_dom = 0; i_dom < n_dom_i; ++i_dom) {

    free (dvtx_coord  [i_dom]);
    free (distrib_vtx [i_dom]);
    free (distrib_edge[i_dom]);
    free (dedge_vtx   [i_dom]);
    PDM_dmesh_free(dm[i_dom]);
  }
  free (dvtx_coord);
  free (distrib_vtx);
  free (distrib_edge);
  free (dedge_vtx);
  free (dm);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {
    free(pn_vtx       [i_dom]);
    free(pvtx_ln_to_gn[i_dom]);
  }
  free(pn_vtx);
  free(pvtx_ln_to_gn);


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
