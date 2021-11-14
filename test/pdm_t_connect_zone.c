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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_dcube_gen.h"
#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

static
void
_exchange_point_list
(
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t  **dface_join_opp,
 PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_UNUSED(group_join_to_zone_cur);
  PDM_UNUSED(group_join_to_zone_opp);

  /*
   * We have for all extraction zone, we need to exchange id
   */
  PDM_g_num_t **distrib_join   = malloc(n_group_join * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **dface_join_opp = malloc(n_group_join * sizeof(PDM_g_num_t *));
  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    distrib_join[i_group_join] = PDM_compute_entity_distribution(comm, dn_face_join);
  }

  *dface_join_opp = malloc(dface_join_idx[n_group_join] * sizeof(PDM_g_num_t));
  PDM_g_num_t* _dface_join_opp = *dface_join_opp;

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    /*
     * Each id in dface_join can be seen as a glabal numbering, in the same order of the opposite window
     * The current part is the implicit numbering
     * The block is a ptr of the opposite part
     */
    // int i_group_zone_opp = group_join_to_zone_opp[i_group_join];
    int i_group_join_opp = group_join_to_join_opp[i_group_join];
    // printf(" i_zone           = %i \n", i_zone);
    // printf(" i_group_join     = %i \n", i_group_join);
    // printf(" i_group_zone_opp = %i \n", i_group_zone_opp);
    // printf(" i_group_join_opp = %i \n", i_group_join_opp);
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    PDM_g_num_t* distrib_join_cur = distrib_join[i_group_join    ];
    PDM_g_num_t* distrib_join_opp = distrib_join[i_group_join_opp];

    PDM_g_num_t *join_ln_to_gn = malloc(dn_face_join * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_face_join; ++i) {
      join_ln_to_gn[i] = distrib_join_cur[i_rank] + i + 1;
    }

    /*
     * Exchange
     */
    PDM_g_num_t *blk_dface_join_cur = &dface_join[dface_join_idx[i_group_join    ]];
    PDM_g_num_t *blk_dface_join_opp = &dface_join[dface_join_idx[i_group_join_opp]];

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_join_opp,
                                 (const PDM_g_num_t **) &join_ln_to_gn,
                                                        &dn_face_join,
                                                        1,
                                                        comm);

    int cst_stride = 1;
    PDM_g_num_t* sub_dface_join_opp = dface_join_opp[dface_join_idx[i_group_join]];
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) blk_dface_join_opp,
                           NULL,
               (void ** ) &sub_dface_join_opp);
    // dface_join_opp[i_group_join] = tmp_dface_join_opp[0];
    // free(tmp_dface_join_opp);

    if(1 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur, dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(sub_dface_join_opp, dn_face_join, "dface_join_opp :: ");
    }

    PDM_block_to_part_free(btp);
    free(join_ln_to_gn);
  }
}






static
void
_deduce_descending_join
(
 int            n_zone,
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 PDM_g_num_t  **extract_face_distribution,
 PDM_g_num_t  **extract_vtx_distribution,
 int          **dextract_face_vtx_idx,
 PDM_g_num_t  **dextract_face_vtx,
 PDM_g_num_t  **dparent_face_g_num,
 PDM_g_num_t  **dparent_vtx_g_num,
 PDM_g_num_t  **pextract_old_to_new,
 double       **dextract_vtx_coord,
 PDM_MPI_Comm   comm
)
{
  PDM_UNUSED(n_zone);
  PDM_UNUSED(dface_join_idx);
  PDM_UNUSED(dface_join);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(extract_face_distribution);
  PDM_UNUSED(extract_vtx_distribution);
  PDM_UNUSED(dextract_face_vtx_idx);
  PDM_UNUSED(dextract_face_vtx);
  PDM_UNUSED(dparent_face_g_num);
  PDM_UNUSED(dparent_vtx_g_num);
  PDM_UNUSED(pextract_old_to_new);
  PDM_UNUSED(dextract_vtx_coord);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   * We have for all extraction zone, we need to exchange id
   */
  PDM_g_num_t **distrib_join   = malloc(n_group_join * sizeof(PDM_g_num_t *));
  PDM_g_num_t **dface_join_opp = malloc(n_group_join * sizeof(PDM_g_num_t *));
  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    distrib_join[i_group_join] = PDM_compute_entity_distribution(comm, dn_face_join);
  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    /*
     * Each id in dface_join can be seen as a glabal numbering, in the same order of the opposite window
     * The current part is the implicit numbering
     * The block is a ptr of the opposite part
     */
    int i_group_zone_opp = group_join_to_zone_opp[i_group_join];
    int i_group_join_opp = group_join_to_join_opp[i_group_join];
    // printf(" i_zone           = %i \n", i_zone);
    // printf(" i_group_join     = %i \n", i_group_join);
    // printf(" i_group_zone_opp = %i \n", i_group_zone_opp);
    // printf(" i_group_join_opp = %i \n", i_group_join_opp);
    int dn_face_join = dface_join_idx[i_group_join+1] - dface_join_idx[i_group_join];
    PDM_g_num_t* distrib_join_cur = distrib_join[i_group_join    ];
    PDM_g_num_t* distrib_join_opp = distrib_join[i_group_join_opp];

    PDM_g_num_t *join_ln_to_gn = malloc(dn_face_join * sizeof(PDM_g_num_t));
    for(int i = 0; i < dn_face_join; ++i) {
      join_ln_to_gn[i] = distrib_join_cur[i_rank] + i + 1;
    }

    /*
     * Exchange
     */
    PDM_g_num_t *blk_dface_join_cur = &dface_join[dface_join_idx[i_group_join    ]];
    PDM_g_num_t *blk_dface_join_opp = &dface_join[dface_join_idx[i_group_join_opp]];

    PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_join_opp,
                                 (const PDM_g_num_t **) &join_ln_to_gn,
                                                        &dn_face_join,
                                                        1,
                                                        comm);

    int cst_stride = 1;
    PDM_g_num_t** tmp_dface_join_opp = NULL;
    PDM_block_to_part_exch2(btp,
                            sizeof(PDM_g_num_t),
                            PDM_STRIDE_CST,
                            &cst_stride,
                   (void *) blk_dface_join_opp,
                            NULL,
                 (void ***) &tmp_dface_join_opp);
    dface_join_opp[i_group_join] = tmp_dface_join_opp[0];
    free(tmp_dface_join_opp);

    if(1 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur          , dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(dface_join_opp[i_group_join], dn_face_join, "dface_join_opp :: ");
    }


    PDM_block_to_part_free(btp);

    free(join_ln_to_gn);

  }

  PDM_g_num_t **dedge_distrib  = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_vtx_idx  = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_vtx      = malloc(n_zone * sizeof(PDM_g_num_t *));
  int         **dedge_face_idx = malloc(n_zone * sizeof(int         *));
  PDM_g_num_t **dedge_face     = malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t **key_ln_to_gn = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));
  PDM_g_num_t **data_send    = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t*));

  // PDM_g_num_t key_mod = ; // Comment definir un bon key_mod ?
  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    /*
     * Generate edge numbering
     */
    int dn_face = extract_face_distribution[i_zone][i_rank+1] - extract_face_distribution[i_zone][i_rank];
    int n_edge_elt_tot = dextract_face_vtx_idx[i_zone][dn_face];

    PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(     n_edge_elt_tot    * sizeof(PDM_g_num_t) );
    int*         tmp_parent_elmt_pos    = (int         *) malloc(     n_edge_elt_tot    * sizeof(int        ) );
    int*         tmp_dface_edge_vtx_idx = (int         *) malloc( ( n_edge_elt_tot + 1) * sizeof(int        ) );
    PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc( 2 * n_edge_elt_tot    * sizeof(PDM_g_num_t) );

    int n_elmt_current = 0;
    int n_edge_current = 0;
    tmp_dface_edge_vtx_idx[0] = 0;
    PDM_poly2d_decomposes_edges(dn_face,
                                &n_elmt_current,
                                &n_edge_current,
                                extract_face_distribution[i_zone][i_rank],
                                -1,
                                dextract_face_vtx[i_zone],
                                dextract_face_vtx_idx[i_zone],
                                tmp_dface_edge_vtx_idx,
                                tmp_dface_edge_vtx,
                                tmp_dface_edge,
                                NULL,
                                NULL,
                                tmp_parent_elmt_pos);
    assert(n_edge_current == n_edge_elt_tot);

    /*
     *  Compute edges connectivity
     */
    int  dn_edge = -1;
    PDM_generate_entitiy_connectivity_raw(comm,
                                          extract_vtx_distribution[i_zone][n_rank],
                                          n_edge_elt_tot,
                                          tmp_dface_edge,
                                          tmp_dface_edge_vtx_idx,
                                          tmp_dface_edge_vtx,
                                          &dn_edge,
                                          &dedge_distrib [i_zone],
                                          &dedge_vtx_idx [i_zone],
                                          &dedge_vtx     [i_zone],
                                          &dedge_face_idx[i_zone],
                                          &dedge_face    [i_zone]);
    free(tmp_parent_elmt_pos    );

    if(1 == 1) {
      PDM_log_trace_array_long(dedge_vtx   [i_zone], dedge_vtx_idx [i_zone][dn_edge], "dedge_vtx :: ");
      PDM_log_trace_array_long(dedge_face  [i_zone], dedge_face_idx[i_zone][dn_edge], "dedge_face :: ");
      // PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_face, "key_ln_to_gn :: ");
    }


    /*
     * Remplissage des clés + préparation buffer d'envoi
     *    - Une clé par edge pour unifier !!
     *    - Maybe besoin de signé le dgroup_face pour identifié les vtx aprés !!!
     *  Attnetion les l'ordre des faces dans le dfaces_group n'est pas le même que dans l'extraction !!!!!
     *  Il faut appliqué le pextract_old_to_new
     *
     *   -> Le mieux est peut-être de prétraité pour ne plus avoir le pb aprés ...
     */
    key_ln_to_gn[i_zone] = malloc(dn_edge * sizeof(PDM_g_num_t)); // Toutes les edges ont une clé car tout vient de l'extraction
    int idx_write = 0;


    // abort();

    // PDM_g_num_t *_dface_join_cur = dface_join    [i_zone];
    // for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

    //   PDM_g_num_t *_dface_join_opp = dface_join_opp[i_zone][i_group_join];
    //   int k = 0;
    //   for(int j = dface_join_idx[i_group_join]; j < dface_join_idx[i_group_join+1]; ++j) {
    //     // key_ln_to_gn[i_zone][idx_write++] = (_dface_join_cur[j] + _dface_join_opp[k++]) % key_mod + 1;
    //     // key_ln_to_gn[i_zone][idx_write++] = shift_group + (_dface_join_cur[j] + _dface_join_opp[k++]) % key_mod + 1;
    //     key_ln_to_gn[i_zone][idx_write++] = (_dface_join_cur[j] + _dface_join_opp[k++]);
    //   }
    // }

    if(1 == 1) {
      PDM_log_trace_array_long(key_ln_to_gn[i_zone], dn_edge, "key_ln_to_gn :: ");
    }


  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    free(distrib_join  [i_group_join]);
    free(dface_join_opp[i_group_join]);
  }
  free(distrib_join);
  free(dface_join_opp);
  free(data_send);

  // for(int i_zone = 0; i_zone < n_zone; ++i_zone) {

  // }

  // for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {

  //   int i_group_zone_opp = group_join_to_zone_opp[i_group_join];
  //   int i_group_join_opp = group_join_to_join_opp[i_group_join];

  // }


  /*
   * Pour les coins coins il faut probablement un echange en plus entre les graph 2 à 2 puis unifier
   * Maybe utiliser le PDM_part_dentity_group_to_pentity_group
   * A la fin de l'algo un vtx est potentiellement lié à plusieurs groupes.
   * On peut merger avec un part_to_block tout les tags d'un vertex en appliquant sur les PL et PLD
   * --> Ca gère le coin
   */


  for(int i_zone = 0; i_zone < n_zone; ++i_zone) {
    free(dedge_distrib [i_zone]);
    free(dedge_vtx_idx [i_zone]);
    free(dedge_vtx     [i_zone]);
    free(dedge_face_idx[i_zone]);
    free(dedge_face    [i_zone]);
    free(key_ln_to_gn  [i_zone]);
  }

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);
  free(key_ln_to_gn);
}

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -h               This message.\n\n");

  exit(exit_code);
}

/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t  *n_vtx_seg,
           double        *length)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n_vtx_seg = atol(argv[i]);
        *n_vtx_seg = (PDM_g_num_t) _n_vtx_seg;
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}


/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_g_num_t        n_vtx_seg = 10;
  double             length    = 1.;
  int                n_zone    = 3;


  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm     comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /* Alloc for distributed mesh */
  int *dn_cell       = (int *) malloc(n_zone * sizeof(int));
  int *dn_face       = (int *) malloc(n_zone * sizeof(int));
  int *dn_vtx        = (int *) malloc(n_zone * sizeof(int));
  int *n_face_group  = (int *) malloc(n_zone * sizeof(int));
  int *dface_vtx_s   = (int *) malloc(n_zone * sizeof(int));
  int *dface_group_s = (int *) malloc(n_zone * sizeof(int));

  PDM_g_num_t  **dface_cell      = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_vtx_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_vtx       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  double       **dvtx_coord      = (double      **) malloc(n_zone * sizeof(double      *));
  int          **dface_group_idx = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_group     = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));
  int          **dface_bnd_idx   = (int         **) malloc(n_zone * sizeof(int         *));
  PDM_g_num_t  **dface_bnd       = (PDM_g_num_t **) malloc(n_zone * sizeof(PDM_g_num_t *));

  PDM_g_num_t  **extract_face_distribution = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **extract_vtx_distribution  = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  int          **dextract_face_vtx_idx     = (int         **) malloc( n_zone * sizeof(int         *));
  PDM_g_num_t  **dextract_face_vtx         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_face_g_num        = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **dparent_vtx_g_num         = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pextract_old_to_new       = (PDM_g_num_t **) malloc( n_zone * sizeof(PDM_g_num_t *));
  double       **dextract_vtx_coord        = (double      **) malloc( n_zone * sizeof(double      *));

  int n_group_join = 2*(n_zone-1);
  int *group_join_to_zone_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_zone_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));

  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  int tmp_i_zone = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join % 2 == 0) {
      group_join_to_join_opp[i_join] = i_join + 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone + 1;
    } else {
      group_join_to_join_opp[i_join] = i_join - 1;
      group_join_to_zone_opp[i_join] = tmp_i_zone++;
    }
  }

  PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");
  PDM_log_trace_array_int(group_join_to_zone_opp, n_group_join, "group_join_to_zone_opp :: ");

  for(int i_group_join = 0; i_group_join < n_group_join+1; ++i_group_join) {
    dface_join_idx[i_group_join] = 0;
  }

  PDM_dcube_t **dcube = (PDM_dcube_t **) malloc(n_zone * sizeof(PDM_dcube_t *));
  int tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {

    dcube[i_zone] = PDM_dcube_gen_init(comm, n_vtx_seg, length, i_zone, 0., 0., PDM_OWNERSHIP_KEEP);
    PDM_dcube_gen_dim_get(dcube         [i_zone],
                          &n_face_group [i_zone],
                          &dn_cell      [i_zone],
                          &dn_face      [i_zone],
                          &dn_vtx       [i_zone],
                          &dface_vtx_s  [i_zone],
                          &dface_group_s[i_zone]);

    PDM_dcube_gen_data_get(dcube          [i_zone],
                          &dface_cell     [i_zone],
                          &dface_vtx_idx  [i_zone],
                          &dface_vtx      [i_zone],
                          &dvtx_coord     [i_zone],
                          &dface_group_idx[i_zone],
                          &dface_group    [i_zone]);

    /*
     * Les faces groups du dcube sont : zmin, zmax, xmin, xmax, ymin, ymax
     * Il faut les séparer en faces de bords et faces raccord, sachant que
     * les zones sont alignées selon X
     */
    int n_bnd = 4;
    int n_jn  = 2;
    if (i_zone == 0){
      n_bnd++;
      n_jn-- ;
    }
    if (i_zone == n_zone-1){
      n_bnd++;
      n_jn-- ;
    }

    // Join numbering (left to right, increasing i_zone)
    printf("n_jn = %i \n", n_jn);
    printf("tmp_i_group_join = %i \n", tmp_i_group_join);

    dface_bnd_idx [i_zone] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_zone][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      int group_size = dface_group_idx[i_zone][igroup+1] - dface_group_idx[i_zone][igroup];
      if (copy_to_bnd) { //Its a boundary
        dface_bnd_idx[i_zone][i_bnd++] = group_size;
      } else { //Its a join
        group_join_to_zone_cur[i_jn-1] = i_zone;
        dface_join_idx[i_jn++] = group_size;
      }
    }
    for (int i = 0; i < n_bnd; i++) {
      dface_bnd_idx[i_zone][i+1] = dface_bnd_idx[i_zone][i+1] + dface_bnd_idx[i_zone][i];
    }

    printf("i_jn = %i \n", i_jn);

    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      if (copy_to_bnd){ //Its a boundary
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_bnd[i_zone][i_bnd++] = dface_group[i_zone][i];
        }
      } else { //Its a join
        int k = 0;
        for (int i = dface_group_idx[i_zone][igroup]; i < dface_group_idx[i_zone][igroup+1]; i++) {
          dface_join[dface_join_idx[i_jn]+k++] = dface_group[i_zone][i];
        }
        i_jn++;
      }
    }

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;
  }

  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx (END) :: ");
  PDM_log_trace_array_int(group_join_to_zone_cur, n_group_join, "group_join_to_zone_cur :: ");

  /*
   * Setup dface_join_opp + group_id in current layout
   */



  /*
   * Algorithm begin - Extract faces
   */
  tmp_i_group_join = 0;
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    int n_jn  = 2;
    if (i_zone == 0       ){n_jn--;}
    if (i_zone == n_zone-1){n_jn--;}

    /*
     *  Now we have all joins create we need to extract them
     */
    PDM_g_num_t* face_distribution = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
    PDM_g_num_t* vtx_distribution  = PDM_compute_entity_distribution(comm, dn_vtx [i_zone]);

    /*
     *  By construction of actual join, each join of same zone are concatenate
     */
    int dn_l_face_join = dface_join_idx[tmp_i_group_join+n_jn] - dface_join_idx[tmp_i_group_join];
    PDM_g_num_t* l_dface_join = &dface_join[dface_join_idx[tmp_i_group_join]];

    PDM_dconnectivity_to_extract_dconnectivity(comm,
                                               dn_l_face_join,
                                               l_dface_join,
                                               face_distribution,
                                               dface_vtx_idx[i_zone],
                                               dface_vtx[i_zone],
                                               &extract_face_distribution[i_zone],
                                               &extract_vtx_distribution[i_zone],
                                               &dextract_face_vtx_idx[i_zone],
                                               &dextract_face_vtx[i_zone],
                                               &dparent_face_g_num[i_zone],
                                               &dparent_vtx_g_num[i_zone],
                                               &pextract_old_to_new[i_zone]);

    int dn_extract_vtx  = extract_vtx_distribution[i_zone][i_rank+1] - extract_vtx_distribution[i_zone][i_rank];

    double** tmp_dextract_vtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          1,
                                          vtx_distribution,
                                          dvtx_coord[i_zone],
                                          &dn_extract_vtx,
                   (const PDM_g_num_t **) &dparent_vtx_g_num[i_zone],
                                          &tmp_dextract_vtx_coord);

    dextract_vtx_coord[i_zone] = tmp_dextract_vtx_coord[0];
    free(tmp_dextract_vtx_coord);

    free(face_distribution);
    free(vtx_distribution);

    /*
     *  Go to nexts join
     */
    tmp_i_group_join += n_jn;

  }




  _deduce_descending_join(n_zone,
                          n_group_join,
                          group_join_to_zone_cur,
                          group_join_to_zone_opp,
                          group_join_to_join_opp,
                          dface_join_idx,
                          dface_join,
                          dface_vtx_idx,
                          dface_vtx,
                          extract_face_distribution,
                          extract_vtx_distribution,
                          dextract_face_vtx_idx,
                          dextract_face_vtx,
                          dparent_face_g_num,
                          dparent_vtx_g_num,
                          pextract_old_to_new,
                          dextract_vtx_coord,
                          comm);

  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx [i_zone]);
    free(dface_bnd     [i_zone]);
    free(extract_face_distribution[i_zone]);
    free(extract_vtx_distribution [i_zone]);
    free(dextract_face_vtx_idx    [i_zone]);
    free(dextract_face_vtx        [i_zone]);
    free(dparent_face_g_num       [i_zone]);
    free(dparent_vtx_g_num        [i_zone]);
    free(pextract_old_to_new      [i_zone]);
    free(dextract_vtx_coord       [i_zone]);
    PDM_dcube_gen_free(dcube[i_zone]);
  }
  free(dcube);
  free(dn_cell);
  free(dn_face);
  free(dn_vtx);
  free(n_face_group);
  free(dface_group_s);
  free(dface_vtx_s);
  free(dface_cell);
  free(dface_vtx_idx);
  free(dface_vtx);
  free(dvtx_coord);
  free(dface_group_idx);
  free(dface_group);
  free(dface_bnd_idx);
  free(dface_bnd);
  free(dface_join_idx);
  free(dface_join);
  free(group_join_to_zone_cur);
  free(group_join_to_zone_opp);
  free(group_join_to_join_opp);

  free(extract_face_distribution);
  free(extract_vtx_distribution );
  free(dextract_face_vtx_idx    );
  free(dextract_face_vtx        );
  free(dparent_face_g_num       );
  free(dparent_vtx_g_num        );
  free(pextract_old_to_new      );
  free(dextract_vtx_coord       );

  PDM_MPI_Finalize();

  return 0;
}
