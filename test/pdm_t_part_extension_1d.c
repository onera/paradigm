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
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_unique.h"
#include "pdm_order.h"

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

  int          **pn_vtx              = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn       = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int           *pflat_n_vtx         = (int           *) malloc( n_domain * sizeof(int           ));
  int           *pflat_n_edge        = (int           *) malloc( n_domain * sizeof(int           ));
  int          **pflat_edge_vtx_idx  = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pflat_edge_vtx      = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t  **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  double       **pflat_vtx_coords    = (double       **) malloc( n_domain * sizeof(double       *));

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

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

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

  int          *pn_vtx_num             = NULL;
  int         **pvtx_num               = NULL;
  int         **pvtx_opp_location_idx  = NULL;
  int         **pvtx_opp_location      = NULL;
  int         **pvtx_opp_interface_idx = NULL;
  int         **pvtx_opp_interface     = NULL;
  int         **pvtx_opp_sens          = NULL;
  PDM_g_num_t **pvtx_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         PDM_BOUND_TYPE_VTX,
                                         pflat_n_vtx,
                                         pflat_vtx_ln_to_gn,
                                         &pn_vtx_num,
                                         &pvtx_num,
                                         &pvtx_opp_location_idx,
                                         &pvtx_opp_location,
                                         &pvtx_opp_interface_idx,
                                         &pvtx_opp_interface,
                                         &pvtx_opp_sens,
                                         &pvtx_opp_gnum);

  /*
   * On rajoute l'échange pour le gnum courant
   * Et avec le part_to_part on recupère le gnum opposé
   * Pour chaque entité d'interface
   * Puis on utilise le ptp du extract_part pour tout transferer le gnum et les triplet et le sens
   */



  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pvtx_num[i_part], pn_vtx_num[i_part], "pvtx_num ::");
      PDM_log_trace_graph_nuplet_int(pvtx_opp_location_idx[i_part],
                                     pvtx_opp_location    [i_part],
                                     3,
                                     pn_vtx_num[i_part], "pvtx_opp_location ::");
      PDM_log_trace_array_int(pvtx_opp_interface[i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_interface ::");
      PDM_log_trace_array_int(pvtx_opp_sens     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_sens ::");
    }
  }

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
  int         **part1_to_part2_interface   = malloc(ln_part_tot * sizeof(int *));

  /*
   * Pour la recursion, il faut idéalement calculer le graphe avec les partition shifter avec
   * PDM_part_generate_entity_graph_comm avec pentity_hint pour accelerer ou choisir
   * Typiquement pour la deuxième passe on prends que tous les sommets étendu au step d'avant non unifié
   */

  /*
   * Dans tout les cas il ne faut pas oublié de faire le graph de comm proprement
   * avec les connexions entres partitions et entre domaine
   * avec PDM_part_domain_interface_as_graph (à reimplementer)
   */

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

      // PDM_log_trace_array_int(ppart_vtx, 4 * ppart_vtx_part_idx[n_part_g[i_dom]], "ppart_vtx");

      part1_to_part2_idx[li_part] = malloc((pflat_n_vtx[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pflat_n_vtx[li_part]);

      int n_entity_bound = ppart_vtx_part_idx[n_part_g[i_dom]];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_vtx[4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        int n_opp = pvtx_opp_location_idx[li_part][idx_entity+1] - pvtx_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }

      // PDM_log_trace_array_int(part1_to_part2_n, pflat_n_vtx[li_part], "part1_to_part2_n ::");

      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_vtx[4*idx_entity]-1;
        int i_proc_opp   = ppart_vtx[4*idx_entity+1];
        int i_part_opp   = ppart_vtx[4*idx_entity+2]-1; // A gere -> Le shift
        int i_entity_opp = ppart_vtx[4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp;
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }


      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        for(int idx_opp = pvtx_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pvtx_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pvtx_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pvtx_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pvtx_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pvtx_opp_interface[li_part][  idx_opp  ];

        }
      }

      // PDM_log_trace_array_int(part1_to_part2_idx      [li_part], pflat_n_vtx[li_part], "part1_to_part2_idx       ::");
      // PDM_log_trace_array_int(part1_to_part2_triplet  [li_part], n_connect_tot       , "part1_to_part2_triplet   ::");
      // PDM_log_trace_array_int(part1_to_part2_interface[li_part], n_connect_tot/3     , "part1_to_part2_interface ::");

      free(part1_to_part2_n);
      li_part += 1;
    }
  }

  /*
   * Pour les periodicités il faut le geré côté part1
   * car côté part2 on a un pb d'extraction de plusieurs gnum mais de periodicité différentes
   * Une fois qu'on a fait le job côté part2 on génére une numérotation asbolu avec des doublés
   * On doit également gardé le doublé existant pour faire l'unification au step d'après
   * Exemple :
   *    - 1 cellules periodiques par noeuds
   *       -> On recupère (1, -1) et (1, 1)
   *       -> On fait un gnum des extensions de partitions : donc 1, 2
   *       -> On shift par rapport à l'interne donc 2, 3
   *    Si on réaplique l'algo pour du rang 2 :
   *       -> On recupère ()
   * ATTENTION : on peut retomber sur la cellule de base
   */


  // ONLY to debug the ref_lnum2
  if(1 == 1) {
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

    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");
    }

    PDM_part_to_part_free(ptp);
  }

  int          *pextract_n_edge                = NULL;
  int         **pextract_vtx_edge_idx          = NULL;
  int         **pextract_vtx_edge              = NULL;
  PDM_g_num_t **pextract_edge_parent_ln_to_gn  = NULL;
  int         **pextract_edge_to_edge_location = NULL;
  PDM_part_to_part_t* ptp_vtx_edge = NULL;

  /*
   *  Il faudrait prendre part1 = uniquement les noeuds frontière (N. Delinger Idea)
   *  Reduction de part1 (n_vtx) -> part (n_vtx concerné par l'extension)
   */

  /*
   * Ici on a 1 seul edge alors qu'on aimerai en avoir 2 !!! Because same gnum for edge inside connectivty vtx_edge
   *   On peut hack la fonction (sinon tout refaire)
   *   Il faut se créer un  vtx_edge + edge_ln_to_gn en dedoublant
   */


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

  /*
   * Echange des numero d'interfaces via le ptp_vtx_edge
   */
  int exch_request = -1;
  int **pextract_edge_interface = NULL;
  PDM_part_to_part_iexch(ptp_vtx_edge,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(int),
                         NULL,
        (const void **)  part1_to_part2_interface,
                         NULL,
            (void ***)   &pextract_edge_interface,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp_vtx_edge, exch_request);

  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(pextract_edge_interface[i_part], pextract_n_edge[i_part], "pextract_edge_interface");
    }
  }


  // PDM_log_trace_array_int(pflat_n_edge, ln_part_tot, "pflat_n_edge");

  PDM_part_to_part_free(ptp_vtx_edge);

  /*
   * Call extract_part
   */
  PDM_extract_part_t* extrp = PDM_extract_part_create(1,
                                                      ln_part_tot,
                                                      ln_part_tot,
                                                      PDM_EXTRACT_PART_KIND_FROM_TARGET,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT,
                                                      PDM_TRUE,
                                                      PDM_OWNERSHIP_KEEP,
                                                      comm);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_extract_part_part_set(extrp,
                              i_part,
                              0,
                              0,
                              pflat_n_edge[i_part],
                              pflat_n_vtx [i_part],
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              pflat_edge_vtx[i_part],
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              pflat_edge_ln_to_gn[i_part],
                              pflat_vtx_ln_to_gn [i_part],
                              pflat_vtx_coords   [i_part]);

    PDM_extract_part_target_set(extrp,
                                i_part,
                                pextract_n_edge               [i_part],
                                pextract_edge_parent_ln_to_gn [i_part],
                                pextract_edge_to_edge_location[i_part]);
  }

  PDM_extract_part_compute(extrp);

  /*
   * On doit créer une nouvelle numérotation absolu qui prend en compte les interfaces
   * MAIS en gardant l'ancienne
   * Donc il faut faire un pdm_gnum avec 2 x le nombre de partition pour tout unifier
   *  A faire pour les entités principales puis descendre, mais attention il faut merger en descendant
   */
  PDM_gen_gnum_t* gnum_edge = PDM_gnum_create(3,
                                              2 * ln_part_tot,
                                              PDM_TRUE,
                                              1.e-6,
                                              comm,
                                              PDM_OWNERSHIP_KEEP);
  PDM_gnum_set_parents_nuplet(gnum_edge, 2); // Si transformation multiple on est ken

  PDM_g_num_t **pflat_edge_ln_to_gn_and_interface = malloc((2 * ln_part_tot) * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    pflat_edge_ln_to_gn_and_interface[i_part] = malloc(2 * pflat_n_edge[i_part] * sizeof(PDM_g_num_t));

    for(int i = 0; i < pflat_n_edge[i_part]; ++i) {
      pflat_edge_ln_to_gn_and_interface[i_part][2*i  ] = pflat_edge_ln_to_gn[i_part][i];
      pflat_edge_ln_to_gn_and_interface[i_part][2*i+1] = 0;
    }

    PDM_gnum_set_from_parents(gnum_edge,
                              i_part,
                              pflat_n_edge[i_part],
                              pflat_edge_ln_to_gn_and_interface[i_part]);

    pflat_edge_ln_to_gn_and_interface[ln_part_tot+i_part] = malloc(2 * pextract_n_edge[i_part] * sizeof(PDM_g_num_t));
    for(int i = 0; i < pextract_n_edge[i_part]; ++i) {
      pflat_edge_ln_to_gn_and_interface[ln_part_tot+i_part][2*i  ] = pextract_edge_parent_ln_to_gn[i_part][i];
      pflat_edge_ln_to_gn_and_interface[ln_part_tot+i_part][2*i+1] = pextract_edge_interface      [i_part][i];
    }

    PDM_gnum_set_from_parents(gnum_edge,
                              ln_part_tot+i_part,
                              pextract_n_edge[i_part],
                              pflat_edge_ln_to_gn_and_interface[ln_part_tot+i_part]);

  }

  PDM_gnum_compute(gnum_edge);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_g_num_t* ln_to_gn1 = PDM_gnum_get(gnum_edge, i_part);
    PDM_log_trace_array_long(ln_to_gn1, pflat_n_edge[i_part], "ln_to_gn1 ::");

    PDM_g_num_t* ln_to_gn2 = PDM_gnum_get(gnum_edge, i_part+ln_part_tot);
    PDM_log_trace_array_long(ln_to_gn2, pextract_n_edge[i_part], "ln_to_gn2 ::");

  }


  for(int i_part = 0; i_part < 2 * ln_part_tot; ++i_part) {
    free(pflat_edge_ln_to_gn_and_interface[i_part]);
  }
  free(pflat_edge_ln_to_gn_and_interface);


  /*
   * Hook le ptp des vtx du extrp_part
   */
  PDM_part_to_part_t* ptp_vtx = NULL;
  PDM_extract_part_part_to_part_get(extrp,
                                    PDM_MESH_ENTITY_VERTEX,
                                    &ptp_vtx,
                                    PDM_OWNERSHIP_KEEP);

  // En fait on veut le parent_ln_to_gn



  /*
   * Unified all points
   */
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int pn_extract_vtx = PDM_extract_part_n_entity_get(extrp,
                                                       i_part,
                                                       PDM_MESH_ENTITY_VERTEX);
    int pn_extract_edge = PDM_extract_part_n_entity_get(extrp,
                                                       i_part,
                                                       PDM_MESH_ENTITY_EDGE);
    int* pextract_edge_vtx     = NULL;
    int* pextract_edge_vtx_idx = NULL;
    PDM_extract_part_connectivity_get(extrp,
                                      i_part,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &pextract_edge_vtx,
                                      &pextract_edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *pextract_vtx_parent_ln_to_gn  = NULL;
    PDM_extract_part_parent_ln_to_gn_get(extrp,
                                         i_part,
                                         PDM_MESH_ENTITY_VERTEX,
                                         &pextract_vtx_parent_ln_to_gn,
                                         PDM_OWNERSHIP_KEEP);

    PDM_g_num_t *pextract_edge_ln_to_gn = NULL;
    PDM_extract_part_ln_to_gn_get(extrp,
                                  i_part,
                                  PDM_MESH_ENTITY_EDGE,
                                  &pextract_edge_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

    double *pextract_vtx_coord = NULL;
    PDM_extract_part_vtx_coord_get(extrp,
                                  i_part,
                                  &pextract_vtx_coord,
                                  PDM_OWNERSHIP_KEEP);

    if(1 == 1) {
      printf("pn_extract_vtx  = %i\n", pn_extract_vtx );
      printf("pn_extract_edge = %i\n", pn_extract_edge);
      PDM_log_trace_array_int (pextract_edge_vtx, 2 * pn_extract_edge, "pextract_edge_vtx ::");
      PDM_log_trace_array_long(pextract_vtx_parent_ln_to_gn , pn_extract_vtx , "pextract_vtx_parent_ln_to_gn  ::");
      PDM_log_trace_array_long(pextract_edge_parent_ln_to_gn[i_part], pn_extract_edge, "pextract_edge_parent_ln_to_gn ::");
    }

    /* Check edge */
    PDM_g_num_t* sorted_edge_ln_to_gn = malloc( pflat_n_edge[i_part] * sizeof(PDM_g_num_t));
    PDM_g_num_t* sorted_vtx_ln_to_gn  = malloc( pflat_n_vtx [i_part] * sizeof(PDM_g_num_t));

    PDM_g_num_t* full_edge_ln_to_gn = PDM_gnum_get(gnum_edge, i_part);
    for(int i_edge = 0; i_edge < pflat_n_edge[i_part]; ++i_edge) {
      // sorted_edge_ln_to_gn[i_edge] = pflat_edge_ln_to_gn[i_part][i_edge];
      sorted_edge_ln_to_gn[i_edge] = full_edge_ln_to_gn[i_edge];
    }
    PDM_sort_long(sorted_edge_ln_to_gn, NULL, pflat_n_edge[i_part]);

    int* vtx_order = malloc( pflat_n_vtx [i_part] * sizeof(int));
    for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
      sorted_vtx_ln_to_gn[i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      vtx_order          [i_vtx] = i_vtx;
    }
    PDM_sort_long(sorted_vtx_ln_to_gn, vtx_order, pflat_n_vtx[i_part]);

    // Cascade
    int pn_extract_edge_keep = 0;
    int pn_extract_vtx_keep  = 0;
    int* pextract_edge_vtx_keep_idx   = malloc((pn_extract_edge + 1) * sizeof(int));
    int* to_keep                      = malloc((pn_extract_edge    ) * sizeof(int));
    int* to_keep_vtx                  = malloc((pn_extract_vtx     ) * sizeof(int));

    for(int i_vtx = 0; i_vtx < pn_extract_vtx; ++i_vtx) {
      to_keep_vtx[i_vtx] = -1;
    }

    pextract_edge_vtx_keep_idx[0] = 0;

    PDM_g_num_t* full_extract_edge_ln_to_gn = PDM_gnum_get(gnum_edge, ln_part_tot+i_part);
    for(int i_edge = 0; i_edge < pn_extract_edge; ++i_edge) {
      // int pos = PDM_binary_search_long(pextract_edge_parent_ln_to_gn[i_part][i_edge],
      int pos = PDM_binary_search_long(full_extract_edge_ln_to_gn[i_edge],
                                       sorted_edge_ln_to_gn,
                                       pflat_n_edge[i_part]);
      to_keep[i_edge] = 0;
      if(pos == -1) {
        to_keep[i_edge] = 1;
        pextract_edge_vtx_keep_idx[pn_extract_edge_keep+1] = pextract_edge_vtx_keep_idx[pn_extract_edge_keep];
        pextract_edge_vtx_keep_idx[pn_extract_edge_keep+1] += 2;
        pn_extract_edge_keep++;

        int i_vtx1 = PDM_ABS (pextract_edge_vtx[2*i_edge  ]);
        int i_vtx2 = PDM_ABS (pextract_edge_vtx[2*i_edge+1]);

        if(to_keep_vtx[i_vtx1-1] == -1) {
          to_keep_vtx[pn_extract_vtx_keep++] = i_vtx1-1;
        }
        if(to_keep_vtx[i_vtx2-1] == -1) {
          to_keep_vtx[pn_extract_vtx_keep++] = i_vtx2-1;
        }


      }
    }

    // PDM_log_trace_array_int(pextract_edge_vtx_keep_idx, pn_extract_edge_keep+1, "pextract_edge_vtx_keep_idx ::");

    PDM_g_num_t *pn_extract_edge_vtx_gnum_keep = malloc(pextract_edge_vtx_keep_idx[pn_extract_edge_keep] * sizeof(PDM_g_num_t));

    pn_extract_edge_keep = 0;
    for(int i_edge = 0; i_edge < pn_extract_edge; ++i_edge) {
      if(to_keep[i_edge] == 1) {
        // Copy
        int idx_write = pextract_edge_vtx_keep_idx[pn_extract_edge_keep];
        int i_vtx1 = PDM_ABS (pextract_edge_vtx[2*i_edge  ]);
        int i_vtx2 = PDM_ABS (pextract_edge_vtx[2*i_edge+1]);
        int sgn1   = PDM_SIGN(pextract_edge_vtx[2*i_edge  ]);
        int sgn2   = PDM_SIGN(pextract_edge_vtx[2*i_edge+1]);
        pn_extract_edge_vtx_gnum_keep[idx_write  ] = sgn1 * pextract_vtx_parent_ln_to_gn[i_vtx1-1];
        pn_extract_edge_vtx_gnum_keep[idx_write+1] = sgn2 * pextract_vtx_parent_ln_to_gn[i_vtx2-1];

        pn_extract_edge_keep++;
      }
    }

    // PDM_log_trace_array_long(pn_extract_edge_vtx_gnum_keep , 2 * pn_extract_edge , "pn_extract_edge_vtx_gnum_keep  ::");

    PDM_g_num_t *extract_vtx_sorted_ln_to_gn = malloc(pn_extract_vtx_keep * sizeof(PDM_g_num_t));
    for(int idx_vtx = 0; idx_vtx < pn_extract_vtx_keep; ++idx_vtx) {
      int i_vtx = to_keep_vtx[idx_vtx];
      extract_vtx_sorted_ln_to_gn[i_vtx] = pextract_vtx_parent_ln_to_gn[i_vtx];
      // to_keep_vtx[idx_vtx] = -1;
    }
    PDM_sort_long(extract_vtx_sorted_ln_to_gn, to_keep_vtx, pn_extract_vtx_keep);



    // int n_extract_vtx_sorted = PDM_inplace_unique_long(extract_vtx_sorted_ln_to_gn, NULL, 0, pextract_edge_vtx_keep_idx[pn_extract_edge_keep]-1);
    // extract_vtx_sorted_ln_to_gn = realloc(extract_vtx_sorted_ln_to_gn, n_extract_vtx_sorted * sizeof(PDM_g_num_t));

    PDM_log_trace_array_long(extract_vtx_sorted_ln_to_gn, pn_extract_vtx_keep , "extract_vtx_sorted_ln_to_gn ::");


    /*
     * Tri indirect des informations recu par le ptb du extract_part
     *  - On cherche après les couples (gnum_opp, -interface) dans l'information transferer du part_domain_interface
     *  - A faire si interface multiple -> Hardcore
     * Cet algo gère l'unification avec le local il faut faire pareil avec la partie étendu aussi (coin
     */


    // Reconstruct vtx
    int *pextract_edge_vtx_keep = malloc(pextract_edge_vtx_keep_idx[pn_extract_edge_keep] * sizeof(int));
    int pn_extract_vtx_keep2 = 0;
    int *pextract_vtx_num = PDM_array_const_int(pn_extract_vtx_keep, -1);
    int *to_keep_vtx2     = PDM_array_const_int(pn_extract_vtx_keep, -1);
    for(int i_vtx = 0; i_vtx < pextract_edge_vtx_keep_idx[pn_extract_edge_keep]; ++i_vtx) {
      PDM_g_num_t gnum = PDM_ABS (pn_extract_edge_vtx_gnum_keep[i_vtx]);
      int sgn          = PDM_SIGN(pn_extract_edge_vtx_gnum_keep[i_vtx]);
      int pos          = PDM_binary_search_long(gnum, sorted_vtx_ln_to_gn, pflat_n_vtx[i_part]);
      if(pos != -1) {
        int old_pos = vtx_order[pos];
        pextract_edge_vtx_keep[i_vtx] = sgn * ( old_pos + 1 );
      } else {

        int pos2 = PDM_binary_search_long(gnum, extract_vtx_sorted_ln_to_gn, pn_extract_vtx_keep);
        assert(pos2 != -1);
        if(pextract_vtx_num[pos2] == -1) {
          pextract_vtx_num[pos2] = pflat_n_vtx[i_part] + pn_extract_vtx_keep2;
          to_keep_vtx2[pn_extract_vtx_keep2++] = to_keep_vtx[pos2];
        }

        pextract_edge_vtx_keep[i_vtx] = sgn * ( pextract_vtx_num[pos2] + 1 );
      }
    }

    if(1 == 1) {
      PDM_log_trace_array_long(pextract_edge_vtx_keep , 2 * pn_extract_edge , "pextract_edge_vtx_keep  ::");
      PDM_log_trace_array_long(pflat_edge_vtx[i_part] , 2 * pflat_n_edge[i_part] , "pflat_edge_vtx  ::");
      PDM_log_trace_array_int(pextract_vtx_num , pn_extract_vtx_keep  , "pextract_vtx_num  ::");
      PDM_log_trace_array_int(to_keep_vtx2     , pn_extract_vtx_keep2 , "to_keep_vtx2      ::");
      printf("pn_extract_vtx_keep2 = %i \n", pn_extract_vtx_keep2);
    }


    /*
     * Hook coordinates and ln_to_gn
     */
    int n_extended_vtx  = pn_extract_vtx_keep2;
    int n_extended_edge = pn_extract_edge_keep;

    PDM_g_num_t *extended_vtx_ln_to_gn = malloc(    n_extended_vtx * sizeof(PDM_g_num_t ));
    double      *extended_vtx_coord    = malloc(3 * n_extended_vtx * sizeof(double      ));

    for(int i_vtx = 0; i_vtx < n_extended_vtx; ++i_vtx) {
      int i_keep_vtx = to_keep_vtx2[i_vtx];
      extended_vtx_ln_to_gn[i_vtx] = pextract_vtx_parent_ln_to_gn[i_keep_vtx];
      extended_vtx_coord[3*i_vtx  ] = pextract_vtx_coord[3*i_keep_vtx  ];
      extended_vtx_coord[3*i_vtx+1] = pextract_vtx_coord[3*i_keep_vtx+1];
      extended_vtx_coord[3*i_vtx+2] = pextract_vtx_coord[3*i_keep_vtx+2];
    }

    /*
     * Concat with current partition and visu
     */
    if(1 == 1) {

      int pn_concat_vtx  = pflat_n_vtx [i_part] + n_extended_vtx;
      int pn_concat_edge = pflat_n_edge[i_part] + n_extended_edge;

      int         *concat_edge_vtx      = malloc(2 * pn_concat_edge * sizeof(int         ));
      PDM_g_num_t *concat_edge_ln_to_gn = malloc(    pn_concat_edge * sizeof(PDM_g_num_t ));
      double      *concat_vtx_coord     = malloc(3 * pn_concat_vtx  * sizeof(double      ));
      PDM_g_num_t *concat_vtx_ln_to_gn  = malloc(    pn_concat_vtx  * sizeof(PDM_g_num_t ));

      for(int i_edge = 0; i_edge < pflat_n_edge[i_part]; ++i_edge) {
        concat_edge_vtx[2*i_edge  ] = pflat_edge_vtx[i_part][2*i_edge  ];
        concat_edge_vtx[2*i_edge+1] = pflat_edge_vtx[i_part][2*i_edge+1];
        concat_edge_ln_to_gn[i_edge] = pflat_edge_ln_to_gn[i_part][i_edge];
      }

      for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
        concat_vtx_coord[3*i_vtx  ] = pflat_vtx_coords[i_part][3*i_vtx  ];
        concat_vtx_coord[3*i_vtx+1] = pflat_vtx_coords[i_part][3*i_vtx+1];
        concat_vtx_coord[3*i_vtx+2] = pflat_vtx_coords[i_part][3*i_vtx+2];
        concat_vtx_ln_to_gn[i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      }


      for(int i_edge = 0; i_edge < n_extended_edge; ++i_edge) {
        concat_edge_vtx[2*(pflat_n_edge[i_part]+i_edge)  ] = pextract_edge_vtx_keep[2*i_edge  ];
        concat_edge_vtx[2*(pflat_n_edge[i_part]+i_edge)+1] = pextract_edge_vtx_keep[2*i_edge+1];
        concat_edge_ln_to_gn[i_edge] = -1; // pflat_edge_ln_to_gn[i_part][i_edge];
      }

      for(int i_vtx = 0; i_vtx < n_extended_vtx; ++i_vtx) {
        concat_vtx_coord[3*(pflat_n_vtx[i_part]+i_vtx)  ] = extended_vtx_coord[3*i_vtx  ];
        concat_vtx_coord[3*(pflat_n_vtx[i_part]+i_vtx)+1] = extended_vtx_coord[3*i_vtx+1];
        concat_vtx_coord[3*(pflat_n_vtx[i_part]+i_vtx)+2] = extended_vtx_coord[3*i_vtx+2];
        concat_vtx_ln_to_gn[pflat_n_vtx[i_part]+i_vtx] = extended_vtx_ln_to_gn[i_vtx];
      }

      char filename[999];
      sprintf(filename, "out_part_concate_vtx_i_part=%i_%i.vtk", i_part, i_rank);
      PDM_vtk_write_std_elements (filename,
                                  pn_concat_vtx,
                                  concat_vtx_coord,
                                  concat_vtx_ln_to_gn,
                                  PDM_MESH_NODAL_BAR2,
                                  pn_concat_edge,
                                  concat_edge_vtx,
                                  concat_edge_ln_to_gn,
                                  0,
                                  NULL,
                                  NULL);

      free(concat_edge_vtx);
      free(concat_edge_ln_to_gn);
      free(concat_vtx_coord);
      free(concat_vtx_ln_to_gn);

    }

    free(extended_vtx_ln_to_gn);
    free(extended_vtx_coord   );


    free(to_keep_vtx2);
    free(vtx_order);
    free(to_keep);
    free(to_keep_vtx);
    free(extract_vtx_sorted_ln_to_gn);
    free(pextract_vtx_num);
    free(sorted_edge_ln_to_gn);
    free(sorted_vtx_ln_to_gn);
    free(pextract_edge_vtx_keep);
    free(pextract_edge_vtx_keep_idx);
    free(pn_extract_edge_vtx_gnum_keep);

  }

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_edge_interface[i_part]);
  }
  free(pextract_edge_interface);

  PDM_extract_part_free(extrp);
  PDM_gnum_free(gnum_edge);

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
    free(part1_to_part2_interface  [i_part]);

    free(pflat_vtx_edge_idx[i_part] );
    free(pflat_vtx_edge    [i_part] );
    free(pflat_edge_vtx_idx[i_part] );
  }
  free(part1_to_part2_idx        );
  // free(part1_to_part2_triplet_idx);
  free(part1_to_part2_triplet    );
  free(part1_to_part2_interface  );

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
  free(pflat_vtx_coords);
  free(pflat_edge_ln_to_gn);
  free(pflat_edge_vtx_idx );
  free(pflat_edge_vtx     );
  free(pflat_vtx_edge_idx );
  free(pflat_vtx_edge     );

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_num             [i_part]);
    free(pvtx_opp_location_idx[i_part]);
    free(pvtx_opp_location    [i_part]);
    free(pvtx_opp_interface   [i_part]);
    free(pvtx_opp_sens        [i_part]);
    free(pvtx_opp_gnum        [i_part]);
  }
  free(pn_vtx_num            );
  free(pvtx_num              );
  free(pvtx_opp_location_idx );
  free(pvtx_opp_location     );
  free(pvtx_opp_interface_idx);
  free(pvtx_opp_interface    );
  free(pvtx_opp_sens         );
  free(pvtx_opp_gnum         );

}


static
void
_part_extension2
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

  int          **pn_vtx              = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge             = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn       = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn      = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int           *pflat_n_vtx         = (int           *) malloc( n_domain * sizeof(int           ));
  int           *pflat_n_edge        = (int           *) malloc( n_domain * sizeof(int           ));
  int          **pflat_edge_vtx_idx  = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pflat_edge_vtx      = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t  **pflat_vtx_ln_to_gn  = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  PDM_g_num_t  **pflat_edge_ln_to_gn = (PDM_g_num_t  **) malloc( n_domain * sizeof(PDM_g_num_t  *));
  double       **pflat_vtx_coords    = (double       **) malloc( n_domain * sizeof(double       *));

  int ln_part_tot = 0;
  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    pn_vtx        [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pvtx_ln_to_gn [i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));
    pn_edge       [i_dom]  = (int          *) malloc( n_part[i_dom] * sizeof(int          ));
    pedge_ln_to_gn[i_dom]  = (PDM_g_num_t **) malloc( n_part[i_dom] * sizeof(PDM_g_num_t *));

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

      pn_edge       [i_dom][i_part] = pflat_n_edge       [ln_part_tot+i_part];
      pedge_ln_to_gn[i_dom][i_part] = pflat_edge_ln_to_gn[ln_part_tot+i_part];

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

      PDM_multipart_part_vtx_coord_get(mpart,
                                       i_dom,
                                       i_part,
                                       &pflat_vtx_coords[ln_part_tot+i_part],
                                       PDM_OWNERSHIP_KEEP);

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

  PDM_g_num_t* shift_by_domain_edge = _compute_offset_ln_to_gn_by_domain(n_domain,
                                                                         n_part,
                                                                         pn_edge,
                                                                         pedge_ln_to_gn,
                                                                         comm);

  /* Shift ln_to_gn */
  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                             shift_by_domain_vtx,
                             1);

  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_edge,
                             pedge_ln_to_gn,
                             shift_by_domain_edge,
                             1);

  /*
   * Changer les noms ?
   */
  int          *pn_vtx_num             = NULL;
  int         **pvtx_num               = NULL;
  int         **pvtx_opp_location_idx  = NULL;
  int         **pvtx_opp_location      = NULL;
  int         **pvtx_opp_interface_idx = NULL;
  int         **pvtx_opp_interface     = NULL;
  int         **pvtx_opp_sens          = NULL;
  PDM_g_num_t **pvtx_opp_gnum          = NULL;

  PDM_part_domain_interface_view_by_part(pdi,
                                         PDM_BOUND_TYPE_VTX,
                                         pflat_n_vtx,
                                         pflat_vtx_ln_to_gn,
                                         &pn_vtx_num,
                                         &pvtx_num,
                                         &pvtx_opp_location_idx,
                                         &pvtx_opp_location,
                                         &pvtx_opp_interface_idx,
                                         &pvtx_opp_interface,
                                         &pvtx_opp_sens,
                                         &pvtx_opp_gnum);

  int          *pn_vtx_opp_gnum_and_itrf = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pvtx_opp_gnum_and_itrf   = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    // int *_pvtx_opp_location = pvtx_opp_location_idx[i_part];
    int *_pvtx_opp_location_idx = pvtx_opp_location_idx[i_part];
    int n_connect_tot = _pvtx_opp_location_idx[pn_vtx_num[i_part]];

    // On garde le lien courant (indirect sort of pvtx_num) + également keep le sens_opp
    pvtx_opp_gnum_and_itrf[i_part] = malloc(2 * n_connect_tot * sizeof(PDM_g_num_t));
    PDM_g_num_t* _pvtx_opp_gnum_and_itrf = pvtx_opp_gnum_and_itrf[i_part];

    for(int idx_entity = 0; idx_entity < pn_vtx_num[i_part]; ++idx_entity) {
      int i_entity = pvtx_num[i_part][idx_entity];
      for(int idx_opp = _pvtx_opp_location_idx[idx_entity]; idx_opp < _pvtx_opp_location_idx[idx_entity+1]; ++idx_opp) {

        _pvtx_opp_gnum_and_itrf[2*idx_opp  ] = pvtx_opp_gnum     [i_part][idx_opp];
        _pvtx_opp_gnum_and_itrf[2*idx_opp+1] = pvtx_opp_interface[i_part][idx_opp];
      }
    }

    PDM_log_trace_array_long(pvtx_opp_gnum_and_itrf[i_part], 2 * n_connect_tot, "pvtx_opp_gnum_and_itrf ::");

    int *vtx_opp_position = malloc(n_connect_tot * sizeof(int));
    pn_vtx_opp_gnum_and_itrf[i_part] = PDM_order_inplace_unique_and_order_long(n_connect_tot, 2, pvtx_opp_gnum_and_itrf[i_part], vtx_opp_position);


    PDM_log_trace_array_int(vtx_opp_position, pn_vtx_opp_gnum_and_itrf[i_part], "vtx_opp_position = ");
    free(vtx_opp_position);


  }



  /*
   * On rajoute l'échange pour le gnum courant
   * Et avec le part_to_part on recupère le gnum opposé
   * Pour chaque entité d'interface
   * Puis on utilise le ptp du extract_part pour tout transferer le gnum et les triplet et le sens
   */



  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

      PDM_log_trace_array_int(pvtx_num[i_part], pn_vtx_num[i_part], "pvtx_num ::");
      PDM_log_trace_graph_nuplet_int(pvtx_opp_location_idx[i_part],
                                     pvtx_opp_location    [i_part],
                                     3,
                                     pn_vtx_num[i_part], "pvtx_opp_location ::");
      PDM_log_trace_array_int (pvtx_opp_interface[i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_interface ::");
      PDM_log_trace_array_int (pvtx_opp_sens     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_sens ::");
      PDM_log_trace_array_long(pvtx_opp_gnum     [i_part], pvtx_opp_location_idx[i_part][pn_vtx_num[i_part]], "pvtx_opp_gnum ::");
    }
  }

  //
  int         **part1_to_part2_idx           = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_triplet_idx   = NULL; //malloc(ln_part_tot * sizeof(int *));
  int         **part1_to_part2_triplet       = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_interface     = malloc(ln_part_tot * sizeof(int         *));

  int         **part1_to_part2_edge_n        = malloc(ln_part_tot * sizeof(int         *));
  int         **part1_to_part2_edge_vtx_n    = malloc(ln_part_tot * sizeof(int         *));
  PDM_g_num_t **part1_to_part2_edge_vtx_gnum = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  PDM_g_num_t **part1_to_part2_edge_gnum     = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  /*
   * Pour la recursion, il faut idéalement calculer le graphe avec les partition shifter avec
   * PDM_part_generate_entity_graph_comm avec pentity_hint pour accelerer ou choisir
   * Typiquement pour la deuxième passe on prends que tous les sommets étendu au step d'avant non unifié
   */

  /*
   * Dans tout les cas il ne faut pas oublié de faire le graph de comm proprement
   * avec les connexions entres partitions et entre domaine
   * avec PDM_part_domain_interface_as_graph (à reimplementer)
   */

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

      // PDM_log_trace_array_int(ppart_vtx, 4 * ppart_vtx_part_idx[n_part_g[i_dom]], "ppart_vtx");

      part1_to_part2_idx[li_part] = malloc((pflat_n_vtx[li_part] + 1) *sizeof(int));
      part1_to_part2_idx[li_part][0] = 0;

      int *part1_to_part2_n = PDM_array_zeros_int(pflat_n_vtx[li_part]);

      int n_entity_bound = ppart_vtx_part_idx[n_part_g[i_dom]];

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity = ppart_vtx[4*idx_entity]-1;
        part1_to_part2_n[i_entity] += 1;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        int n_opp = pvtx_opp_location_idx[li_part][idx_entity+1] - pvtx_opp_location_idx[li_part][idx_entity];
        part1_to_part2_n[i_entity] += n_opp;
      }


      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {
        part1_to_part2_idx[li_part][i_entity+1] = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
        part1_to_part2_n[i_entity] = 0;
      }

      int n_connect_tot = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]];
      part1_to_part2_triplet  [li_part] = malloc(n_connect_tot   * sizeof(int));
      part1_to_part2_interface[li_part] = malloc(n_connect_tot/3 * sizeof(int));

      printf("n_connect_tot = %i \n", n_connect_tot);

      for(int idx_entity = 0; idx_entity < n_entity_bound; ++idx_entity) {
        int i_entity     = ppart_vtx[4*idx_entity]-1;
        int i_proc_opp   = ppart_vtx[4*idx_entity+1];
        int i_part_opp   = ppart_vtx[4*idx_entity+2]-1; // A gere -> Le shift
        int i_entity_opp = ppart_vtx[4*idx_entity+3]-1;

        int idx_write = part1_to_part2_idx[li_part][i_entity] + part1_to_part2_n[i_entity];
        part1_to_part2_triplet[li_part][idx_write  ] = i_proc_opp;
        part1_to_part2_triplet[li_part][idx_write+1] = i_part_opp;
        part1_to_part2_triplet[li_part][idx_write+2] = i_entity_opp;

        idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
        part1_to_part2_interface[li_part][idx_write] = 0;

        // part1_to_part2_triplet[li_part][idx_entity+1] = part1_to_part2_triplet[li_part][idx_entity] + 3;
      }

      /* From interface */
      for(int idx_entity = 0; idx_entity < pn_vtx_num[li_part]; ++idx_entity) {
        int i_entity = pvtx_num[li_part][idx_entity];
        for(int idx_opp = pvtx_opp_location_idx[li_part][idx_entity  ];
                idx_opp < pvtx_opp_location_idx[li_part][idx_entity+1]; ++idx_opp) {

          int idx_write = part1_to_part2_idx[li_part][i_entity] + 3 * part1_to_part2_n[i_entity];
          part1_to_part2_triplet  [li_part][idx_write  ] = pvtx_opp_location [li_part][3*idx_opp  ];
          part1_to_part2_triplet  [li_part][idx_write+1] = pvtx_opp_location [li_part][3*idx_opp+1];
          part1_to_part2_triplet  [li_part][idx_write+2] = pvtx_opp_location [li_part][3*idx_opp+2];

          // Il faudra le faire en stride variable si periodicité composé
          idx_write = part1_to_part2_idx[li_part][i_entity]/3 + part1_to_part2_n[i_entity]++;
          part1_to_part2_interface[li_part][idx_write  ] = pvtx_opp_interface[li_part][  idx_opp  ];

        }
      }

      PDM_log_trace_array_int(part1_to_part2_idx      [li_part], pflat_n_vtx[li_part]+1, "part1_to_part2_idx       ::");
      PDM_log_trace_array_int(part1_to_part2_triplet  [li_part], n_connect_tot         , "part1_to_part2_triplet   ::");
      PDM_log_trace_array_int(part1_to_part2_interface[li_part], n_connect_tot/3       , "part1_to_part2_interface ::");

      /*
       * Creation des buffers d'envoi des connectivités
       */
      int n_send = part1_to_part2_idx[li_part][pflat_n_vtx[li_part]]/3;
      part1_to_part2_edge_vtx_n[li_part] = malloc(n_send * sizeof(int));
      part1_to_part2_edge_n    [li_part] = malloc(n_send * sizeof(int));
      int n_send_edge     = 0;
      int n_send_edge_vtx = 0;

      int *_pflat_vtx_edge_idx        = pflat_vtx_edge_idx       [li_part];
      int *_part1_to_part2_edge_vtx_n = part1_to_part2_edge_vtx_n[li_part];
      int *_part1_to_part2_edge_n     = part1_to_part2_edge_n    [li_part];

      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          _part1_to_part2_edge_vtx_n[idx] = 0;
          _part1_to_part2_edge_n    [idx] = 0;

          /* Copy */
          for(int idx_edge = _pflat_vtx_edge_idx[i_entity]; idx_edge < _pflat_vtx_edge_idx[i_entity+1]; ++idx_edge ) {
            int i_edge = PDM_ABS(pflat_vtx_edge[li_part][idx_edge])-1;
            _part1_to_part2_edge_n    [idx] += 1;
            _part1_to_part2_edge_vtx_n[idx] += pflat_edge_vtx_idx[li_part][i_edge+1] - pflat_edge_vtx_idx[li_part][i_edge];
          }
          n_send_edge     += _part1_to_part2_edge_n    [idx];
          n_send_edge_vtx += _part1_to_part2_edge_vtx_n[idx];
        }
      }

      printf("n_send_edge     = %i \n", n_send_edge);
      printf("n_send_edge_vtx = %i \n", n_send_edge_vtx);

      part1_to_part2_edge_gnum    [li_part] = malloc(n_send_edge     * sizeof(PDM_g_num_t));
      part1_to_part2_edge_vtx_gnum[li_part] = malloc(n_send_edge_vtx * sizeof(PDM_g_num_t));
      PDM_g_num_t *_part1_to_part2_edge_gnum     = part1_to_part2_edge_gnum    [li_part];
      PDM_g_num_t *_part1_to_part2_edge_vtx_gnum = part1_to_part2_edge_vtx_gnum[li_part];

      /* Fill loop */
      n_send_edge     = 0;
      n_send_edge_vtx = 0;
      for(int i_entity = 0; i_entity < pflat_n_vtx[li_part]; ++i_entity) {

        for(int idx = part1_to_part2_idx  [li_part][i_entity]/3; idx < part1_to_part2_idx  [li_part][i_entity+1]/3; ++idx) {

          /* Copy */
          for(int idx_edge = _pflat_vtx_edge_idx[i_entity]; idx_edge < _pflat_vtx_edge_idx[i_entity+1]; ++idx_edge ) {
            int i_edge = PDM_ABS(pflat_vtx_edge[li_part][idx_edge])-1;
            _part1_to_part2_edge_gnum    [n_send_edge++    ] = pflat_edge_ln_to_gn[li_part][i_edge];

            for(int idx_vtx = pflat_edge_vtx_idx[li_part][i_edge]; idx_vtx < pflat_edge_vtx_idx[li_part][i_edge+1]; ++idx_vtx) {
              int i_vtx = PDM_ABS (pflat_edge_vtx[li_part][idx_vtx])-1;
              int sgn   = PDM_SIGN(pflat_edge_vtx[li_part][idx_vtx]);
              _part1_to_part2_edge_vtx_gnum[n_send_edge_vtx++] = sgn * pflat_vtx_ln_to_gn[li_part][i_vtx];
            }
          }
        }
      }


      free(part1_to_part2_n);
      li_part += 1;
    }
  }


  /*
   * Create part_to_part to exchange all data in opposit part
   */
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

  /*
   * Exchange edge gnum
   */
  int exch_request = -1;
  int         **pextract_edge_n    = NULL;
  PDM_g_num_t **pextract_edge_gnum = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(PDM_g_num_t),
        (const int **)   part1_to_part2_edge_n,
        (const void **)  part1_to_part2_edge_gnum,
                         &pextract_edge_n,
            (void ***)   &pextract_edge_gnum,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);


  /*
   * Exchange connectivity
   */
  int         **pextract_edge_vtx_n    = NULL;
  PDM_g_num_t **pextract_edge_vtx_gnum = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_VAR_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(PDM_g_num_t),
        (const int **)   part1_to_part2_edge_vtx_n,
        (const void **)  part1_to_part2_edge_vtx_gnum,
                         &pextract_edge_vtx_n,
            (void ***)   &pextract_edge_vtx_gnum,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(part1_to_part2_edge_n   [i_part]);
    free(part1_to_part2_edge_gnum[i_part]);
    free(part1_to_part2_edge_vtx_n   [i_part]);
    free(part1_to_part2_edge_vtx_gnum[i_part]);
  }
  free(part1_to_part2_edge_vtx_n   );
  free(part1_to_part2_edge_vtx_gnum);
  free(part1_to_part2_edge_n   );
  free(part1_to_part2_edge_gnum);

  int **pextract_edge_interface = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(int),
                         NULL,
        (const void **)  part1_to_part2_interface,
                         NULL,
            (void ***)   &pextract_edge_interface,
                         &exch_request);
  PDM_part_to_part_iexch_wait(ptp, exch_request);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_log_trace_array_int(ref_lnum2[i_part], n_ref_lnum2[i_part], "ref_lnum2 :");
  }

  int **pextract_edge_vtx_idx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    pextract_edge_vtx_idx[i_part] = malloc((n_ref_lnum2[i_part] + 1) * sizeof(int));

    pextract_edge_vtx_idx[i_part][0] = 0;
    for(int i_edge = 0; i_edge < n_ref_lnum2[i_part]; ++i_edge) {
      pextract_edge_vtx_idx[i_part][i_edge+1] = pextract_edge_vtx_idx[i_part][i_edge] + pextract_edge_vtx_n[i_part][i_edge];
    }

  }

  /* Verbose all recv data */
  if(1 == 1) {
    for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
      PDM_log_trace_array_int(pextract_edge_interface[i_part], n_ref_lnum2[i_part], "pextract_edge_interface : ");
      PDM_log_trace_array_int(pextract_edge_vtx_n    [i_part], n_ref_lnum2[i_part], "pextract_edge_vtx_n     : ");
      PDM_log_trace_array_int(pextract_edge_n        [i_part], n_ref_lnum2[i_part], "pextract_edge_n         : ");

      int n_recv_edge     = 0;
      int n_recv_edge_vtx = 0;
      for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
        n_recv_edge     += pextract_edge_n    [i_part][i_ref];
        n_recv_edge_vtx += pextract_edge_vtx_n[i_part][i_ref];
      }
      PDM_log_trace_array_long(pextract_edge_vtx_gnum    [i_part], n_recv_edge_vtx, "pextract_edge_vtx_gnum : ");
      PDM_log_trace_array_long(pextract_edge_gnum        [i_part], n_recv_edge    , "pextract_edge_gnum     : ");
    }
  }

  /*
   * Post-treatment
   *   - Create new global numbering with new entity
   *   - Merge connectivity
   *   - Create new global entity for descending connectivity
   */
  PDM_gen_gnum_t* gen_gnum_edge = PDM_gnum_create(3,
                                                  ln_part_tot,
                                                  PDM_TRUE,
                                                  1.e-6,
                                                  comm,
                                                  PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_edge, 2);


  int          *pn_edge_only_by_interface        = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pedge_ln_to_gn_only_by_interface = malloc(ln_part_tot * sizeof(PDM_g_num_t *));

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    pn_edge_only_by_interface[i_part] = 0;

    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        pn_edge_only_by_interface[i_part]++;
      }
    }

    pedge_ln_to_gn_only_by_interface[i_part] = malloc(2 * pn_edge_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    PDM_g_num_t *_pedge_ln_to_gn_only_by_interface = pedge_ln_to_gn_only_by_interface[i_part];
    pn_edge_only_by_interface[i_part] = 0;

    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        int idx_write = pn_edge_only_by_interface[i_part]++;
        _pedge_ln_to_gn_only_by_interface[2*idx_write  ] = _pextract_edge_gnum     [i_ref];
        _pedge_ln_to_gn_only_by_interface[2*idx_write+1] = _pextract_edge_interface[i_ref];
      }
    }

    PDM_gnum_set_from_parents(gen_gnum_edge,
                              i_part,
                              pn_edge_only_by_interface[i_part],
                              pedge_ln_to_gn_only_by_interface[i_part]);

  }

  PDM_gnum_compute(gen_gnum_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    PDM_g_num_t* extented_edge_ln_to_gn = PDM_gnum_get(gen_gnum_edge, i_part);
    PDM_log_trace_array_long(extented_edge_ln_to_gn, pn_edge_only_by_interface[i_part], "extented_edge_ln_to_gn ::");

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    /* Update */
    int idx_read = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {
        _pextract_edge_gnum[i_ref] = extented_edge_ln_to_gn[idx_read++] + shift_by_domain_edge[n_domain];
      }
    }

    PDM_log_trace_array_long(_pextract_edge_gnum, n_ref_lnum2[i_part], "_pextract_edge_gnum (Update) : ");

  }


  PDM_gnum_free(gen_gnum_edge);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pedge_ln_to_gn_only_by_interface[i_part]);
  }

  free(pn_edge_only_by_interface);
  free(pedge_ln_to_gn_only_by_interface);


  /*
   * ATTENTION : On peut avoir 2 fois le même edge ou face (surtout au deuxieme passage )
   *   Faire un unique puis une indirection (on garde les buffers tel quel c'est plus simple je pense )
   *   Pour chaque sommet connecté par interrface on envoie le gnum pour faire la pair d'arrivé
   *   On doit reconstruire les infos d'interface pour le step d'après (exemple edge --> C'est ultra MANDATORY pour le périodique infini )
   *   Si tout les adjacences sont bonnes on s'arrête sinon on reprendre avec la nouvelle numérotation absolu et les nouvelles interfaces
   *   Attention au multidomaine !!!
   */


  /*
   * Echange vtx_gnum + tri indirect des doublet
   */

  /*
   * Create relation d'interface par edge
   */

  /*
   * Test d'arrêt par adjacence
   */

  /*
   * Creation numero absolu vtx
   */
  PDM_gen_gnum_t* gen_gnum_vtx = PDM_gnum_create(3,
                                                 ln_part_tot,
                                                 PDM_TRUE,
                                                 1.e-6,
                                                 comm,
                                                 PDM_OWNERSHIP_KEEP);

  PDM_gnum_set_parents_nuplet(gen_gnum_vtx, 2);

  /*
   *  Il faut mettre QUE les nouveu vertex qui ne sont pas locaux OU une translation explicit via les raccords existe
   *    - On prends chaque couple (vtx, interface)
   *    - Si edge interne -> Rien a faire
   *    - Si le sommet est dans une liste d'interface tradusible -> On le traduit
   *    - Sinon on le rajoute dans les gnum a calculer
   */

  int          *pn_vtx_only_by_interface        = malloc(ln_part_tot * sizeof(int          ));
  PDM_g_num_t **pvtx_ln_to_gn_only_by_interface = malloc(ln_part_tot * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* _pvtx_opp_gnum_and_itrf = pvtx_opp_gnum_and_itrf[i_part];

    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    int         *_pextract_edge_vtx_idx  = pextract_edge_vtx_idx [i_part];
    PDM_g_num_t *_pextract_edge_vtx_gnum = pextract_edge_vtx_gnum[i_part];

    pn_vtx_only_by_interface[i_part] = 0;

    /* Count */
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, pn_vtx_opp_gnum_and_itrf[i_part]);
          if(pos != -1) {
            pn_vtx_only_by_interface[i_part]++;
          }
          printf("Search (%i/%i) -> pos = %i \n", vtx_g_num, _pextract_edge_interface[i_ref], pos);
        }
      }
    }

    /* Fill */
    pvtx_ln_to_gn_only_by_interface[i_part] = malloc(2 * pn_vtx_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    pn_vtx_only_by_interface[i_part] = 0;
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {
      if(_pextract_edge_interface[i_ref] != 0) {

        /* Blinder d'assert */
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, pn_vtx_opp_gnum_and_itrf[i_part]);
          if(pos != -1) {
            int idx_write = pn_vtx_only_by_interface[i_part]++;
            pvtx_ln_to_gn_only_by_interface[i_part][2*idx_write  ] = vtx_g_num;
            pvtx_ln_to_gn_only_by_interface[i_part][2*idx_write+1] = -_pextract_edge_interface[i_ref];
          }
        }
      }
    }

    PDM_gnum_set_from_parents(gen_gnum_vtx,
                              i_part,
                              pn_vtx_only_by_interface[i_part],
                              pvtx_ln_to_gn_only_by_interface[i_part]);
  }
  PDM_gnum_compute(gen_gnum_vtx);

  /*
   * Unify edge_vtx
   */
  int **pextract_edge_vtx = malloc(ln_part_tot * sizeof(int *));
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {

    PDM_g_num_t* extented_vtx_ln_to_gn = PDM_gnum_get(gen_gnum_vtx, i_part);
    PDM_log_trace_array_long(extented_vtx_ln_to_gn, pn_vtx_only_by_interface[i_part], "extented_vtx_ln_to_gn ::");

    /* Sort unique the new gnum to unify */
    int         *extended_vtx_order           = malloc( pn_vtx_only_by_interface[i_part] * sizeof(int        ));
    PDM_g_num_t *extented_vtx_ln_to_gn_sorted = malloc( pn_vtx_only_by_interface[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pn_vtx_only_by_interface[i_part]; ++i_vtx) {
      extented_vtx_ln_to_gn_sorted[i_vtx] = extented_vtx_ln_to_gn[i_vtx];
      extended_vtx_order          [i_vtx] = i_vtx;
    }

    PDM_sort_long(extented_vtx_ln_to_gn_sorted, extended_vtx_order, pn_vtx_only_by_interface[i_part]);


    /*
     * Sort current part vtx_ln_to_gn
     */
    int         *vtx_order            = malloc( pflat_n_vtx[i_part] * sizeof(int        ));
    PDM_g_num_t *pvtx_ln_to_gn_sorted = malloc( pflat_n_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < pflat_n_vtx[i_part]; ++i_vtx) {
      pvtx_ln_to_gn_sorted[i_vtx] = pflat_vtx_ln_to_gn[i_part][i_vtx];
      vtx_order           [i_vtx] = i_vtx;
    }
    PDM_sort_long(pvtx_ln_to_gn_sorted, vtx_order, pflat_n_vtx[i_part]);


    PDM_g_num_t *_pvtx_opp_gnum_and_itrf  = pvtx_opp_gnum_and_itrf [i_part];
    int         *_pextract_edge_interface = pextract_edge_interface[i_part];
    PDM_g_num_t *_pextract_edge_gnum      = pextract_edge_gnum     [i_part];

    int         *_pextract_edge_vtx_idx  = pextract_edge_vtx_idx [i_part];
    PDM_g_num_t *_pextract_edge_vtx_gnum = pextract_edge_vtx_gnum[i_part];

    pextract_edge_vtx[i_part] = malloc( _pextract_edge_vtx_idx[n_ref_lnum2[i_part]] * sizeof(int));
    int         *_pextract_edge_vtx  = pextract_edge_vtx [i_part];

    /* */
    for(int i_ref = 0; i_ref < n_ref_lnum2[i_part]; ++i_ref) {

      /*
       * First case :
       *   - Edge is move by interface
       */
      if(_pextract_edge_interface[i_ref] != 0) {
        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          PDM_g_num_t gnum_to_find[2] = {vtx_g_num, -_pextract_edge_interface[i_ref]};
          int pos = PDM_order_binary_search_long(gnum_to_find, _pvtx_opp_gnum_and_itrf, 2, pn_vtx_opp_gnum_and_itrf[i_part]);
          printf(" --- Search (%i/%i) -> pos = %i \n", vtx_g_num, _pextract_edge_interface[i_ref], pos);
          if(pos != -1) {
            _pextract_edge_vtx[idx_edge] = ( pflat_n_vtx[i_part] + pos + 1); // ATTENTION SIGN

          } else {
            int pos_int = PDM_binary_search_long(vtx_g_num, pvtx_ln_to_gn_sorted, pflat_n_vtx[i_part]);
            assert(pos_int != -1);
            int orig_pos = vtx_order[pos_int];
            _pextract_edge_vtx[idx_edge] = ( orig_pos + 1);

            printf("\t --- Search interior (%i) -> pos_int = %i \n", vtx_g_num, pos_int);

          }
        }
      } else { // Second case : Interior edge

        for(int idx_edge = _pextract_edge_vtx_idx[i_ref]; idx_edge < _pextract_edge_vtx_idx[i_ref+1]; ++idx_edge) {

          PDM_g_num_t vtx_g_num = _pextract_edge_vtx_gnum[idx_edge];
          int pos_int = PDM_binary_search_long(vtx_g_num, pvtx_ln_to_gn_sorted, pflat_n_vtx[i_part]);
          assert(pos_int != -1);
          int orig_pos = vtx_order[pos_int];
          _pextract_edge_vtx[idx_edge] = ( orig_pos + 1);

        }
      }
    }

    PDM_log_trace_connectivity_int(_pextract_edge_vtx_idx, _pextract_edge_vtx, n_ref_lnum2[i_part], "pextract_edge_vtx ::");




    free(extended_vtx_order          );
    free(extented_vtx_ln_to_gn_sorted);

    free(vtx_order);
    free(pvtx_ln_to_gn_sorted);

  }



  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_edge_vtx[i_part]);
  }
  free(pextract_edge_vtx);


  PDM_gnum_free(gen_gnum_vtx);

  /*
   * Pour faire la cascade IL FAUT garder les opp_location !!!!
   */

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_opp_gnum_and_itrf[i_part]);
  }
  free(pvtx_opp_gnum_and_itrf);
  free(pn_vtx_opp_gnum_and_itrf);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_ln_to_gn_only_by_interface[i_part]);
  }
  free(pn_vtx_only_by_interface       );
  free(pvtx_ln_to_gn_only_by_interface);


  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pextract_edge_interface[i_part]);
    free(pextract_edge_vtx_n    [i_part]);
    free(pextract_edge_vtx_idx  [i_part]);
    free(pextract_edge_vtx_gnum [i_part]);
    free(pextract_edge_n        [i_part]);
    free(pextract_edge_gnum     [i_part]);
  }
  free(pextract_edge_interface);
  free(pextract_edge_vtx_n    );
  free(pextract_edge_vtx_idx  );
  free(pextract_edge_vtx_gnum );
  free(pextract_edge_n        );
  free(pextract_edge_gnum     );



  PDM_part_to_part_free(ptp);

  /* Unshift ln_to_gn */
  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                             shift_by_domain_vtx,
                             -1);

  _offset_ln_to_gn_by_domain(n_domain,
                             n_part,
                             pn_edge,
                             pedge_ln_to_gn,
                             shift_by_domain_edge,
                             -1);

  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(part1_to_part2_idx        [i_part]);
    // free(part1_to_part2_triplet_idx[i_part]);
    free(part1_to_part2_triplet    [i_part]);
    free(part1_to_part2_interface  [i_part]);

    free(pflat_vtx_edge_idx[i_part] );
    free(pflat_vtx_edge    [i_part] );
    free(pflat_edge_vtx_idx[i_part] );
  }
  free(part1_to_part2_idx        );
  // free(part1_to_part2_triplet_idx);
  free(part1_to_part2_triplet    );
  free(part1_to_part2_interface  );

  free(pflat_vtx_edge_idx );
  free(pflat_vtx_edge     );

  free(shift_by_domain_vtx);
  free(shift_by_domain_edge);
  free(n_part_g);

  for(int i_dom = 0; i_dom < n_domain; ++i_dom) {

    free(pn_edge       [i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
  }

  free(pn_edge             );
  free(pedge_ln_to_gn      );
  free(pn_vtx             );
  free(pvtx_ln_to_gn      );
  free(pflat_n_vtx        );
  free(pflat_n_edge       );
  free(pflat_edge_vtx_idx );
  free(pflat_edge_vtx     );
  free(pflat_vtx_ln_to_gn );
  free(pflat_edge_ln_to_gn);
  free(pflat_vtx_coords   );
  for(int i_part = 0; i_part < ln_part_tot; ++i_part) {
    free(pvtx_num             [i_part]);
    free(pvtx_opp_location_idx[i_part]);
    free(pvtx_opp_location    [i_part]);
    free(pvtx_opp_interface   [i_part]);
    free(pvtx_opp_sens        [i_part]);
    free(pvtx_opp_gnum        [i_part]);
  }
  free(pn_vtx_num            );
  free(pvtx_num              );
  free(pvtx_opp_location_idx );
  free(pvtx_opp_location     );
  free(pvtx_opp_interface_idx);
  free(pvtx_opp_interface    );
  free(pvtx_opp_sens         );
  free(pvtx_opp_gnum         );



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
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_IMPLICIT;
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

  PDM_domain_interface_free(dom_itrf);


  /*
   * Extention de partition
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
  if(0 == 1) {
    _part_extension(n_depth,
                    n_dom_i,
                    n_part,
                    comm,
                    mpart,
                    pdi);
  }
  _part_extension2(n_depth,
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

  PDM_part_domain_interface_free(pdi);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
