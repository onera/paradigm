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
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_sort.h"
#include "pdm_array.h"
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

    int i_group_join_opp = group_join_to_join_opp[i_group_join];
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
    PDM_g_num_t* sub_dface_join_opp = &_dface_join_opp[dface_join_idx[i_group_join]];
    PDM_block_to_part_exch(btp,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_CST,
                           &cst_stride,
                  (void *) blk_dface_join_opp,
                           NULL,
               (void ** ) &sub_dface_join_opp);

    if(0 == 1) {
      PDM_log_trace_array_long(blk_dface_join_cur, dn_face_join, "dface_join_cur :: ");
      PDM_log_trace_array_long(sub_dface_join_opp, dn_face_join, "dface_join_opp :: ");
    }

    PDM_block_to_part_free(btp);
    free(join_ln_to_gn);
  }

  for(int i_group_join = 0; i_group_join < n_group_join; ++i_group_join) {
    free(distrib_join[i_group_join]);
  }
  free(distrib_join);
}


static void
merge_zone
(
 int            n_zone,
 int            n_group_join,
 int           *group_join_to_zone_cur,
 int           *group_join_to_zone_opp,
 int           *group_join_to_join_opp,
 int           *dface_join_idx,
 PDM_g_num_t   *dface_join,
 PDM_g_num_t   *dface_join_opp,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 PDM_g_num_t  **dface_cell,
 PDM_MPI_Comm   comm
)
{

  /*
   * On a plusieurs connectivités distribué  :
   *   -> On veut unifier ce qui est commun
   *   -> Et actualiser les connectivités
   *
   *  PDM_multi_block_to_part
   *    -> Avec face_ln_to_gn = implicite face_ln_to_gn - face à supprimer
   *          --> PDM_redistribute
   *
   */


  /*
   * Par connectivité ascendante ou descendane ?
   *    --> face_cell
   *    --> Par ascendance --> Trop dure je pense
   *   Dans pdm_mesh_adpation c'est en ascendant mais on part des vertex .
   *   Subtile on reconstruit dans ce sens --> vtx -> edge -> face -> cell
   *   Mais on construit les connectivités descendantes --> edge_vtx -> face_edge...
   */


  /*
   * Par descendance :
   *   --> Calcul du dcell_face --> Attention au signe !!!!!!!
   *   --> Quand on transpose la connectvity il faut d'aborder transformé le dface_cell en orienté
   *   --> This one : PDM_setup_connectivity_idx
   *
   *   PMD_multiblock_to_part avec dcell_face + cell_ln_to_gn concateante
   *     --> Donc un dcell_face re-repartie mais faux en face
   *     --> Si on est des ouf --> Reorder_block_hilbert
   *   --> Pour update le dcell_face
   *         --> mbtp avec ln_to_gn = face_to_keep
   *   --> On utilise le dface_jon + dface_join_opp
   *   Attention au desequilibre de charge possibleeeeeeeee
   *   Pour l'update des numero de face --> block_to_part avec un block variable de numero de faces old_to_new
   *     --> Si la face est supprimé strid = 0
   *
   */

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
  int                n_zone    = 2;


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

  int n_group_join = 2*(n_zone-1);
  int          *dface_join_idx  = (int         *) malloc((n_group_join + 1) * sizeof(int        ));
  PDM_g_num_t  *dface_join      = NULL; // A allouer propremet

  //int n_group_join = 2*2*(n_zone-1); // Try 2*2 jns
  int *group_join_to_zone_cur = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_zone_opp = (int *) malloc( n_group_join * sizeof(int));
  int *group_join_to_join_opp = (int *) malloc( n_group_join * sizeof(int));


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
  // Test 2*2 jns
  /*group_join_to_join_opp[0] = 2;*/
  /*group_join_to_join_opp[1] = 3;*/
  /*group_join_to_join_opp[2] = 0;*/
  /*group_join_to_join_opp[3] = 1;*/
  /*group_join_to_zone_opp[0] = 1;*/
  /*group_join_to_zone_opp[1] = 1;*/
  /*group_join_to_zone_opp[2] = 0;*/
  /*group_join_to_zone_opp[3] = 0;*/


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
    //Try different setup with 2*2 jns
    /*
    n_jn = 2;
    n_face_group[i_zone] += 1;
    int *dface_group_idx_new = malloc((n_face_group[i_zone]+2) * sizeof(int));
    //For zone0, jn is 3 ; for zone1, jn is 2
    int jn_group = -1;
    if (i_zone == 0) jn_group = 3;
    else if (i_zone == 1) jn_group = 2;
    for (int i = 0; i < n_face_group[i_zone]+1; i++) {
      if (i <= jn_group) {
        dface_group_idx_new[i] = dface_group_idx[i_zone][i];
      }
      else if (i==jn_group+1) {
        int nface_this_group = dface_group_idx[i_zone][i] - dface_group_idx[i_zone][i-1];
        dface_group_idx_new[i] = dface_group_idx[i_zone][i-1] + (nface_this_group / 2);
        dface_group_idx_new[i+1] = dface_group_idx[i_zone][i-1] + nface_this_group;
      }
      else {
        dface_group_idx_new[i] = dface_group_idx[i_zone][i-1];
      }
    }
    PDM_log_trace_array_int(dface_group_idx[i_zone], 6+1, "dfacegroupidx :");
    PDM_log_trace_array_int(dface_group_idx_new, 7+1, "dfacegroupidxnew :");
    PDM_log_trace_array_long(dface_group[i_zone], dface_group_idx[i_zone][6], "dfacegroup :");
    free(dface_group_idx[i_zone]);
    dface_group_idx[i_zone] = dface_group_idx_new;
    */
    //

    // Join numbering (left to right, increasing i_zone)

    dface_bnd_idx [i_zone] = (int *) malloc((n_bnd        + 1) * sizeof(int));

    // First pass to count and allocate
    int i_bnd = 1;
    int i_jn  = tmp_i_group_join+1;
    dface_bnd_idx[i_zone][0]  = 0;
    // dface_join_idx[0] = 0;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      /*int copy_to_bnd = -1; //Try 2*2 jns*/
      /*if (i_zone == 0) copy_to_bnd = (igroup != 3) && (igroup !=4);*/
      /*if (i_zone == 1) copy_to_bnd = (igroup != 2) && (igroup !=3);*/
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


    for (int i = tmp_i_group_join; i < i_jn-1; i++) {
      dface_join_idx[i+1] = dface_join_idx[i+1] + dface_join_idx[i];
    }

    /* A bit stupid but complicated to made it in ohter way for a clear test */
    dface_join = realloc(dface_join, dface_join_idx[i_jn-1] * sizeof(PDM_g_num_t));

    // Second pass to copy
    dface_bnd [i_zone] = (PDM_g_num_t *) malloc(dface_bnd_idx [i_zone][n_bnd        ] * sizeof(PDM_g_num_t));
    i_bnd = 0;
    i_jn  = tmp_i_group_join;
    for (int igroup = 0; igroup < n_face_group[i_zone]; igroup++) {
      int copy_to_bnd = (igroup != 2 || (igroup == 2 && i_zone == 0)) && (igroup != 3 || (igroup == 3 && i_zone == n_zone-1));
      /*int copy_to_bnd = -1; //Try 2*2 jns*/
      /*if (i_zone == 0) copy_to_bnd = (igroup != 3) && (igroup !=4);*/
      /*if (i_zone == 1) copy_to_bnd = (igroup != 2) && (igroup !=3);*/
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

  /*
   * Setup dface_join_opp + group_id in current layout
   */
  PDM_g_num_t *dface_join_opp = NULL;
  _exchange_point_list(n_group_join,
                       group_join_to_join_opp,
                       dface_join_idx,
                       dface_join,
                       &dface_join_opp,
                       comm);

  log_trace("Global join data (%d)\n", n_group_join);
  PDM_log_trace_array_int(group_join_to_zone_cur, n_group_join, "group_join_to_zone_cur :: ");
  PDM_log_trace_array_int(group_join_to_join_opp, n_group_join, "group_join_to_join_opp :: ");
  PDM_log_trace_array_int(group_join_to_zone_opp, n_group_join, "group_join_to_zone_opp :: ");

  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "dface_join_idx :: ");
  PDM_log_trace_array_long(dface_join    , dface_join_idx[n_group_join], "dface_join     :: ");
  PDM_log_trace_array_long(dface_join_opp, dface_join_idx[n_group_join], "dface_join_opp :: ");

  // New version begins

  // Extract all the jn faces
  
  PDM_g_num_t *face_per_block_offset = (PDM_g_num_t *) malloc((n_zone+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *vtx_per_block_offset  = (PDM_g_num_t *) malloc((n_zone+1) * sizeof(PDM_g_num_t));
  face_per_block_offset[0] = 0;
  vtx_per_block_offset[0] = 0;
  PDM_MPI_Allreduce(dn_face, &face_per_block_offset[1], n_zone, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_MPI_Allreduce(dn_vtx , &vtx_per_block_offset[1] , n_zone, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);
  PDM_array_accumulate_gnum(face_per_block_offset, n_zone+1);
  PDM_array_accumulate_gnum(vtx_per_block_offset , n_zone+1);

  //int n_face_join = dface_join_idx[n_group_join]; //
  int n_face_join = 0; // Put faces only once
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    if (i_join <= group_join_to_join_opp[i_join])
      n_face_join += 2*(dface_join_idx[i_join+1] - dface_join_idx[i_join]);
  }

  PDM_g_num_t *multi_gnum        = malloc(n_face_join * sizeof(PDM_g_num_t));

  int idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    log_trace("Treat jn %i\n", i_join);
    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        multi_gnum[idx++] = dface_join[i_face_jn] + face_per_block_offset[i_zone_cur];
        multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];
      }
      /*for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {*/
        /*multi_gnum[idx++] = dface_join_opp[i_face_jn] + face_per_block_offset[i_zone_opp];*/
      /*}*/
    }
  }
  assert (idx == n_face_join);

  PDM_log_trace_array_long(face_per_block_offset, n_zone+1, "face_per_block_offset :: ");
  PDM_log_trace_array_long(vtx_per_block_offset,  n_zone+1, "vtx_per_block_offset :: ");
  PDM_log_trace_array_long(multi_gnum, n_face_join, "multi_gnum :: ");

  PDM_g_num_t **all_face_distribution = malloc(n_zone * sizeof(PDM_g_num_t*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    all_face_distribution[i_zone] = PDM_compute_entity_distribution(comm, dn_face[i_zone]);
  }
  PDM_multi_block_to_part_t *mptb = PDM_multi_block_to_part_create(face_per_block_offset,
                                                                   n_zone,
                                            (const PDM_g_num_t **) all_face_distribution,
                                            (const PDM_g_num_t **)&multi_gnum,
                                                                  &n_face_join,
                                                                   1,
                                                                   comm);
  //Prepare data to send : face -> vtx connectivity 
  int **face_vtx_n       = malloc(n_zone * sizeof(int*));
  PDM_g_num_t **face_vtx_shifted = malloc(n_zone * sizeof(int*));
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    face_vtx_n[i_zone]       = malloc(dn_face[i_zone] * sizeof(int));
    face_vtx_shifted[i_zone] = malloc(dface_vtx_idx[i_zone][dn_face[i_zone]] * sizeof(PDM_g_num_t));

    for (int i_face = 0; i_face < dn_face[i_zone]; i_face++) {
      face_vtx_n[i_zone][i_face] = dface_vtx_idx[i_zone][i_face+1] - dface_vtx_idx[i_zone][i_face];
      for (int j = dface_vtx_idx[i_zone][i_face]; j < dface_vtx_idx[i_zone][i_face+1]; j++) {
        face_vtx_shifted[i_zone][j] = dface_vtx[i_zone][j] + vtx_per_block_offset[i_zone];
      }
    }
    //PDM_log_trace_array_long(face_vtx_shifted[i_zone], dface_vtx_idx[i_zone][dn_face[i_zone]], "face_vtx_shifted :: ");
  }
  //int **face_vtx = malloc(n_zone * sizeof(int*));



  int         **part_stride = NULL;
  PDM_g_num_t **part_data   = NULL;
  PDM_multi_block_to_part_exch2(mptb,
                                sizeof(PDM_g_num_t),
                                PDM_STRIDE_VAR,
                                face_vtx_n,
                 (void **)      face_vtx_shifted,
                               &part_stride,
                 (void ***)    &part_data);

  int         *face_vtx_both_n   = part_stride[0];
  int         *face_vtx_both_idx = PDM_array_new_idx_from_sizes_int(face_vtx_both_n, n_face_join);
  PDM_g_num_t *face_vtx_both     = part_data[0];

  free(part_data);
  free(part_stride);
  PDM_multi_block_to_part_free(mptb);

  int n_recv = 0;
  for (int i = 0; i < n_face_join; i++)
    n_recv += face_vtx_both_n[i];

  log_trace("Face vtx received after MBTP\n");
  PDM_log_trace_array_int (face_vtx_both_idx, n_face_join+1, "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx_both, n_recv, "face_vtx :: ");


  int *face_vtx_idx = malloc(((n_face_join / 2 )+1) * sizeof(int));
  PDM_g_num_t *face_vtx = malloc((n_recv / 2) * sizeof(PDM_g_num_t));
  PDM_g_num_t *face_vtx_opp = malloc((n_recv / 2 ) * sizeof(PDM_g_num_t));
  face_vtx_idx[0] = 0;

  idx = 0;
  int recv_idx = 0;
  int data_idx = 0;
  int recv_data_idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int n_face_this_jn = dface_join_idx[i_join+1] - dface_join_idx[i_join];
    if (i_join <= i_join_opp) {
      for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
        int n_edge = face_vtx_both_n[2*recv_idx]; //Jump opposite face stride
        // Update idx array
        face_vtx_idx[idx+1] = face_vtx_idx[idx] + n_edge;
        recv_idx++;
        idx++;
        // Updata data
        memcpy(&face_vtx    [data_idx], &face_vtx_both[recv_data_idx], n_edge * sizeof(PDM_g_num_t));
        memcpy(&face_vtx_opp[data_idx], &face_vtx_both[recv_data_idx+n_edge], n_edge * sizeof(PDM_g_num_t));
        recv_data_idx += 2*n_edge;
        data_idx += n_edge;
      }
    }
  }

  //Old version, when jn were not interlaced
  
  /*
  //We have jn_idx, prepare jn_data_idx to allow easier access
  int *jn_data_idx = malloc(n_group_join * sizeof(int));
  jn_data_idx[0] = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    jn_data_idx[i_join+1] = 0;
    for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
      jn_data_idx[i_join+1] += part_stride[0][k];
    }
    jn_data_idx[i_join+1] += jn_data_idx[i_join];
  }
  PDM_log_trace_array_int(dface_join_idx, n_group_join+1, "jn_idx :: ");
  PDM_log_trace_array_int(jn_data_idx, n_group_join+1, "jn_data_idx :: ");

  //Reorder all received face->edge connectivity to have pairs of faces
  int *face_edge_idx = malloc(((n_face_join / 2 )+1) * sizeof(int));
  PDM_g_num_t *face_edge = malloc((n_recv / 2) * sizeof(PDM_g_num_t));
  PDM_g_num_t *face_edge_opp = malloc((n_recv / 2 ) * sizeof(PDM_g_num_t));
  face_edge_idx[0] = 0;
  idx = 0;
  int data_idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int n_face_this_jn = dface_join_idx[i_join+1] - dface_join_idx[i_join];

    log_trace("Post treat jn %i\n", i_join);
    if (i_join <= i_join_opp) {
      //Copy idx
      for (int k = dface_join_idx[i_join]; k < dface_join_idx[i_join+1]; k++) {
        face_edge_idx[idx+1] = face_edge_idx[idx] + part_stride[0][k];
        idx++;
      }
      //Copy data
      int data_size = jn_data_idx[i_join+1] - jn_data_idx[i_join];
      log_trace("data size is %d\n", data_size);
      assert (jn_data_idx[i_join_opp+1] - jn_data_idx[i_join_opp] == data_size);
      memcpy(&face_edge[data_idx], &part_data[0][jn_data_idx[i_join]], data_size * sizeof(PDM_g_num_t));
      memcpy(&face_edge_opp[data_idx], &part_data[0][jn_data_idx[i_join_opp]], data_size * sizeof(PDM_g_num_t));
      data_idx += data_size;
    }
  }
  */
  assert (idx == n_face_join / 2 );
  assert (data_idx == n_recv / 2);

  log_trace("Split face_vtx & face_vtx donor \n");
  PDM_log_trace_array_int(face_vtx_idx, idx+1, "face_vtx_idx :: ");
  PDM_log_trace_array_long(face_vtx, n_recv/2, "face_vtx :: ");
  PDM_log_trace_array_long(face_vtx_opp, n_recv/2, "face_vtx_opp :: ");

  //Now we have some pairs of faces (each pair appears only one) + face_vtx for this pairs
  
  //To get the number of unique vertex, we can do a part to block using face_vtx_both as gnum
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &face_vtx_both,
                                                      NULL,
                                                     &face_vtx_both_idx[n_face_join],
                                                      1,
                                                      comm);
  int n_vtx;
  int n_vtx_loc = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_MPI_Allreduce(&n_vtx_loc, &n_vtx, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  PDM_part_to_block_free(ptb);

  //Decompose into edges
  int ext_dn_face = n_face_join;
  int n_edge_elt = face_vtx_both_idx[ext_dn_face];
  int n_elmt_current = 0;
  int n_edge_current = 0;
  PDM_g_num_t *face_distri = PDM_compute_entity_distribution(comm, ext_dn_face);

  PDM_g_num_t* tmp_dface_edge         = (PDM_g_num_t *) malloc(n_edge_elt       * sizeof(PDM_g_num_t) );
  int*         tmp_parent_elmt_pos    = (int         *) malloc(n_edge_elt       * sizeof(int        ) );
  int*         tmp_dface_edge_vtx_idx = (int         *) malloc((n_edge_elt + 1) * sizeof(int        ) );
  PDM_g_num_t* tmp_dface_edge_vtx     = (PDM_g_num_t *) malloc(2*n_edge_elt     * sizeof(PDM_g_num_t) );

  tmp_dface_edge_vtx_idx[0] = 0;
  PDM_poly2d_decomposes_edges(ext_dn_face,
                              &n_elmt_current,
                              &n_edge_current,
                              face_distri[i_rank],
                              -1,
                              face_vtx_both,
                              face_vtx_both_idx,
                              tmp_dface_edge_vtx_idx, //Numéro de sommet des edges
                              tmp_dface_edge_vtx,
                              tmp_dface_edge,  // Numéro des edges pour une face
                              NULL,
                              NULL,
                              tmp_parent_elmt_pos);
  assert(n_edge_current == n_edge_elt);

  /*PDM_log_trace_array_long(tmp_dface_edge, n_edge_elt, "dface_edge :: ");*/
  /*PDM_log_trace_connectivity_long(tmp_dface_edge_vtx_idx, tmp_dface_edge_vtx, n_edge_elt, "dface_edge :: ");*/

  PDM_g_num_t dn_edge;
  PDM_g_num_t *dedge_distrib  = NULL;
  int         *dedge_vtx_idx  = NULL;
  PDM_g_num_t *dedge_vtx      = NULL;
  int         *dedge_face_idx = NULL;
  PDM_g_num_t *dedge_face     = NULL;
  PDM_generate_entitiy_connectivity_raw(comm,
                                        n_vtx, // n vtx tot
                                        n_edge_elt,
                                        tmp_dface_edge,
                                        tmp_dface_edge_vtx_idx,
                                        tmp_dface_edge_vtx,
                                        &dn_edge,
                                        &dedge_distrib ,
                                        &dedge_vtx_idx ,
                                        &dedge_vtx     ,
                                        &dedge_face_idx,
                                        &dedge_face    );
  free(tmp_parent_elmt_pos    );

  log_trace("Edges rebuild\n");
  PDM_log_trace_array_long(dedge_distrib, n_rank+1, "dedge_distri ::");
  //PDM_log_trace_array_long(dedge_vtx, 2*dn_edge, "dedge_vtx ::");
  //PDM_log_trace_connectivity_long(dedge_face_idx, dedge_face, dn_edge, "dedge_face ::");

  // Count the number of internal & external edges
  int dn_internal_edge = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {
    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    dn_internal_edge += (int) (n_face_this_edge > 1);
  }
  int dn_external_edge = dn_edge - dn_internal_edge;
  log_trace("dn internal edges is %i \n", dn_internal_edge);
  log_trace("dn external edges is %i \n", dn_external_edge);



  // Transport data to edges from face to build the key


  //Prepare numbering
  PDM_g_num_t *dedge_face_abs = (PDM_g_num_t *) malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int         *dedge_face_sgn = (int         *) malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  for(int i = 0; i < dedge_face_idx[dn_edge]; ++i) {
    dedge_face_abs[i] = PDM_ABS (dedge_face[i]);
    dedge_face_sgn[i] = PDM_SIGN(dedge_face[i]);
  }

  //Prepare data on multi_gnum (use same ordering than before) to transfert to edges
  int *dextract_face_group_idNEW = malloc(ext_dn_face*sizeof(int));
  PDM_g_num_t *dextract_face_joinNEW     = malloc(ext_dn_face*sizeof(PDM_g_num_t));
  PDM_g_num_t *dextract_face_join_oppNEW = malloc(ext_dn_face*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int i_join_opp = group_join_to_join_opp[i_join];
    int i_zone_cur = group_join_to_zone_cur[i_join];
    int i_zone_opp = group_join_to_zone_opp[i_join];

    if (i_join <= i_join_opp) {
      for (int i_face_jn = dface_join_idx[i_join]; i_face_jn < dface_join_idx[i_join+1]; i_face_jn++) {
        //Face data
        dextract_face_joinNEW    [idx] = dface_join[i_face_jn]; //Here unshifted
        dextract_face_join_oppNEW[idx] = dface_join_opp[i_face_jn];
        dextract_face_group_idNEW[idx++] = i_join;
        //Opp face data
        dextract_face_joinNEW    [idx] = dface_join_opp[i_face_jn];
        dextract_face_join_oppNEW[idx] = dface_join[i_face_jn];
        dextract_face_group_idNEW[idx++] = i_join_opp;
      }
    }
  }
  assert (idx == n_face_join);


  PDM_g_num_t* dedge_face_join       = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  PDM_g_num_t* dedge_face_join_opp   = malloc(dedge_face_idx[dn_edge] * sizeof(PDM_g_num_t));
  int        * dedge_face_group_id   = malloc(dedge_face_idx[dn_edge] * sizeof(int        ));
  PDM_block_to_part_t *btp = PDM_block_to_part_create(face_distri,
                               (const PDM_g_num_t **) &dedge_face_abs,
                                                      &dedge_face_idx[dn_edge],
                                                      1,
                                                      comm);
  int cst_stride = 1;
  PDM_block_to_part_exch(btp,
                         sizeof(int),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_group_idNEW,
                         NULL,
             (void ** ) &dedge_face_group_id);
  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_joinNEW,
                         NULL,
             (void ** ) &dedge_face_join);

  PDM_block_to_part_exch(btp,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST,
                         &cst_stride,
                (void *) dextract_face_join_oppNEW,
                         NULL,
             (void ** ) &dedge_face_join_opp);

  PDM_block_to_part_free(btp);

  log_trace("Transport data on edges\n");
  PDM_log_trace_array_int (dedge_face_idx,      dn_edge+1,               "dedge_face_idx      ::");
  PDM_log_trace_array_int (dedge_face_group_id, dedge_face_idx[dn_edge], "dedge_face_group_id ::");
  PDM_log_trace_array_long(dedge_face_join    , dedge_face_idx[dn_edge], "dedge_face_join     ::");
  PDM_log_trace_array_long(dedge_face_join_opp, dedge_face_idx[dn_edge], "dedge_face_join_opp ::");

  // Begin hash table
  PDM_g_num_t *key_ln_to_gn = malloc(dn_internal_edge * sizeof(PDM_g_num_t)); 
  int         *stride_one   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_two   = malloc(dn_internal_edge * sizeof(int        ));
  int         *stride_four  = malloc(dn_internal_edge * sizeof(int        ));

  int         *zone_id              = malloc(  dn_internal_edge * sizeof(int        ));
  PDM_g_num_t *data_send_connect    = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_edge_g_num = malloc(  dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_group      = malloc(4*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_sens       = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));
  PDM_g_num_t *data_send_face_g_num = malloc(2*dn_internal_edge * sizeof(PDM_g_num_t));

  int i_int_edge = 0;
  int idx_write2 = 0;
  int idx_write4 = 0;
  for(int i_edge = 0; i_edge < dn_edge; ++i_edge) {

    int n_face_this_edge = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    if (n_face_this_edge == 1) {
      continue;
    }
    assert (n_face_this_edge == 2);

    stride_one [i_int_edge] = 1;
    stride_two [i_int_edge] = 2;
    stride_four[i_int_edge] = 4;
    //Retrive zone id using group id of any of two faces
    data_send_edge_g_num[i_int_edge] = dedge_distrib[i_rank] + i_edge + 1;
    int group_id = dedge_face_group_id[dedge_face_idx[i_edge]];
    zone_id[i_int_edge] = group_join_to_zone_cur[group_id];

    int key = 0;
    for(int j = dedge_face_idx[i_edge]; j < dedge_face_idx[i_edge+1]; ++j) { //Do it for the two faces data
      key += (dedge_face_join[j] + dedge_face_join_opp[j]);

      data_send_face_g_num[idx_write2]   = dedge_face[j];
      idx_write2++;
      //data_send_sens      [idx_write2++] = dedge_face_group_sens[j];

      data_send_connect[idx_write4] = dedge_face_join    [j];
      data_send_group[idx_write4++] = dedge_face_group_id[j];

      data_send_connect[idx_write4] = dedge_face_join_opp[j];
      data_send_group[idx_write4++] = group_join_to_join_opp[dedge_face_group_id[j]]; // This is group join opp
    }
    key_ln_to_gn[i_int_edge] = key;

    i_int_edge++;
  }
  assert(idx_write2 == 2*dn_internal_edge);
  assert(idx_write4 == 4*dn_internal_edge);

  PDM_log_trace_array_long(key_ln_to_gn, dn_internal_edge, "key_ln_to_gn :: ");

  //Attention, pb d'équilibrage car les clés sont réparties vers la fin ... Un proc risque
  // de se retrouver avec tt les clés
                       ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                    (PDM_g_num_t **) &key_ln_to_gn,
                                                      NULL,
                                                     &dn_internal_edge,
                                                      1,
                                                      comm);
  // Get protocol data
  int blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  // This one will contain 2 for non conflicting hash key, more otherwise
  int *gnum_n_occurences = PDM_part_to_block_block_gnum_count_get(ptb);

  int gnum_n_occurences_tot = 0;
  for (int k = 0; k < blk_size; k++) {
    gnum_n_occurences_tot += gnum_n_occurences[k];
  }
  int *blk_entity_idx    = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  int *blk_data_face_idx = (int *) malloc((gnum_n_occurences_tot + 1)*sizeof(int));
  for (int k = 0; k < gnum_n_occurences_tot + 1; k++) {
    blk_entity_idx[k]    = 4*k;
    blk_data_face_idx[k] = 2*k;
  }
  
  // Exch data : zone_id, edge gnum, face & facedonor id, face & facedonor group id, face_sens and faces_gnum

  int *unused_recv_stride = NULL;
  int *blk_zone_id   = NULL;
  int exch_size = PDM_part_to_block_exch(ptb,
                                         sizeof(int),
                                         PDM_STRIDE_VAR,
                                         -1,
                               (int **)  &stride_one,
                               (void **) &zone_id,
                                         &unused_recv_stride,
                               (void **) &blk_zone_id);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  PDM_g_num_t* blk_edge_g_num           = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_one,
                           (void **) &data_send_edge_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_edge_g_num);
  free(unused_recv_stride); // Same as gnum_n_occurences 
  assert (exch_size == gnum_n_occurences_tot);

  /*
  PDM_g_num_t* blk_data_sens   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_sens,
                                     &unused_recv_stride,
                           (void **) &blk_data_sens);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences 
  assert (exch_size == 2*gnum_n_occurences_tot);
  */

  PDM_g_num_t* blk_data_face_g_num   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_two,
                           (void **) &data_send_face_g_num,
                                     &unused_recv_stride,
                           (void **) &blk_data_face_g_num);
  free(unused_recv_stride); // Same as 2*gnum_n_occurences
  assert (exch_size == 2*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_connect   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                          (int **)  &stride_four,
                          (void **) &data_send_connect,
                                    &unused_recv_stride,
                          (void **) &blk_data_connect);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);

  PDM_g_num_t* blk_data_group   = NULL;
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(int),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int **)  &stride_four,
                           (void **) &data_send_group,
                                     &unused_recv_stride,
                           (void **) &blk_data_group);
  free(unused_recv_stride); // Same as 4*gnum_n_occurences 
  assert (exch_size == 4*gnum_n_occurences_tot);


  PDM_log_trace_array_int(gnum_n_occurences   , blk_size               , "gnum_n_occurences   :: ");
  PDM_log_trace_array_int(blk_zone_id         , gnum_n_occurences_tot  , "blk_zone_id         :: ");
  PDM_log_trace_array_long(blk_edge_g_num     , gnum_n_occurences_tot  , "blk_edge_g_num      :: ");
  PDM_log_trace_array_long(blk_data_face_g_num, 2*gnum_n_occurences_tot, "blk_data_face_g_num :: ");
  //PDM_log_trace_array_long(blk_data_sens      , 2*gnum_n_occurences_tot, "blk_data_sens       :: ");
  PDM_log_trace_array_long(blk_data_connect   , 4*gnum_n_occurences_tot, "blk_data_connect    :: ");
  PDM_log_trace_array_long(blk_data_group     , 4*gnum_n_occurences_tot, "blk_data_group      :: ");


  free(key_ln_to_gn        );
  free(data_send_connect   );
  free(data_send_group     );
  free(data_send_edge_g_num);
  free(data_send_sens      );
  free(stride_one          );
  free(stride_two          );
  free(stride_four         );
  free(data_send_face_g_num);
  free(zone_id             );

  // Post treatemement : resolve conflicting keys
  int n_max_entity_per_key = 0;
  for(int i = 0; i < blk_size; ++i) {
    n_max_entity_per_key = PDM_MAX(gnum_n_occurences   [i], n_max_entity_per_key);
  }
  int n_max_connec = 4*n_max_entity_per_key;


  //Number of edge to treat on each zone (Usefull ?)
  int* zone_id_n = PDM_array_zeros_int(n_zone);
  for (int i_edge = 0; i_edge < gnum_n_occurences_tot; i_edge++) {
    zone_id_n[blk_zone_id[i_edge]]++;
  }

  PDM_log_trace_array_int(zone_id_n        , n_zone                , "zone_id_n         :: ");
  log_trace("n_max_entity_per_key = %i \n", n_max_entity_per_key);
  log_trace("n_max_connec         = %i \n", n_max_connec);
  int*          already_treat    = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  int*          same_entity_idx  = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  //int*          sens_entity      = (int         *) malloc(  n_max_entity_per_key * sizeof(int          ) );
  PDM_g_num_t *results_edge     = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));
  PDM_g_num_t *results_edge_opp = malloc(gnum_n_occurences_tot * sizeof(PDM_g_num_t));

  /* Reset to fill */
  PDM_array_reset_int(zone_id_n, n_zone, 0);

  idx  = 0;
  int idx_w = 0;
  for(int i_key = 0; i_key < blk_size; ++i_key) {

    int n_matching_edge = gnum_n_occurences[i_key];

    log_trace(" i_key = %i | n_matching_edge = %i \n", i_key, n_matching_edge);

    /* Reset */
    PDM_array_reset_int(already_treat, n_matching_edge, -1);

    /* Loop over all entitys in conflict and sort all */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {

      //Each internal edge comes with 2 faces -> 4 data
      int beg    = 4*(idx+i_entity);

      // Caution inplace sort !!!!
      PDM_sort_long(&blk_data_connect[beg], NULL, 4);
      PDM_sort_int (&blk_data_group  [beg], NULL, 4);

      if(0 == 1) {
        PDM_log_trace_array_long(&blk_data_connect[beg], 4, "blk_data_connect (sort) :: ");
        PDM_log_trace_array_int (&blk_data_group  [beg], 4, "blk_data_group   (sort) :: ");
      }
    }

    /*
     *  Identify pair or invalid other
     */
    for(int i_entity = 0; i_entity < n_matching_edge; ++i_entity) {
      int beg1    = 4*(idx+i_entity);

      int n_same         = 0;
      int i_entity2_same = -1;

      if(already_treat[i_entity] == 1) {
        continue;
      }

      for(int i_entity2 = i_entity+1; i_entity2 < n_matching_edge; ++i_entity2) {

        if(already_treat[i_entity2] == -1) {
          int beg2    = 4*(idx+i_entity2);

          if(!PDM_array_are_equal_int(&blk_data_group[beg1], &blk_data_group[beg2], 4)) {
            continue;
          }

          if(!PDM_array_are_equal_int(&blk_data_connect[beg1], &blk_data_connect[beg2], 4)) {
            continue;
          }

          already_treat[i_entity2] = 1;
          same_entity_idx[n_same++] = i_entity2;
          i_entity2_same = i_entity2;

        }
      } /* End for i_entity2 */
      assert(n_same == 1);

      //log_trace("i_entity = %i | i_entity2_same = %i | n_same = %i \n", i_entity, i_entity2_same, n_same);

      /*
       * Renvoie des resultats :
       *    - Par edge --> Il faut donc également le numero de zones
       */
      int edge_idx     = (idx+i_entity);
      int edge_idx_opp = (idx+i_entity2_same);

      int i_zone_cur = blk_zone_id[edge_idx];
      int i_zone_opp = blk_zone_id[edge_idx_opp];

      // Set data for edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx_opp];
      //results_edge    [i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx];
      //results_edge_opp[i_zone_cur][zone_id_n[i_zone_cur]  ] = blk_edge_g_num[edge_idx_opp];
      //results_zone_opp[i_zone_cur][zone_id_n[i_zone_cur]++] = blk_zone_id[edge_idx_opp];

      // Set data for opposite edge
      results_edge    [idx_w] = blk_edge_g_num[edge_idx_opp];
      results_edge_opp[idx_w++] = blk_edge_g_num[edge_idx];
      //results_edge    [i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx_opp];
      //results_edge_opp[i_zone_opp][zone_id_n[i_zone_opp]  ] = blk_edge_g_num[edge_idx];
      //results_zone_opp[i_zone_opp][zone_id_n[i_zone_opp]++] = blk_zone_id[edge_idx];

      // Renvoi de tout les edges candidats à travers les faces ???
      already_treat[i_entity] = 1;
    }

    idx  += n_matching_edge;
  }
  assert (idx_w == gnum_n_occurences_tot);

  PDM_part_to_block_free(ptb); // Needed for gnum_n_occurences

  log_trace("Conflict resolved, gnum are\n");
  PDM_log_trace_array_long(results_edge, gnum_n_occurences_tot, "edge gnum ::");
  PDM_log_trace_array_long(results_edge_opp, gnum_n_occurences_tot, "edge gnum opp::");

  free(zone_id_n);


  // Send back result on edge distribution (we were on key distribution)
                       ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_NOTHING,
                                                       1.,
                                                      &results_edge,
                                                       dedge_distrib,
                                                      &gnum_n_occurences_tot,
                                                       1,
                                                       comm);

  assert (PDM_part_to_block_n_elt_block_get(ptb) == dn_internal_edge);
  PDM_g_num_t *dedge_gnum     = malloc(PDM_part_to_block_n_elt_block_get(ptb) * sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_gnum_opp = NULL;

  PDM_g_num_t *dedge_gnum_tmp = PDM_part_to_block_block_gnum_get(ptb);
  memcpy(dedge_gnum, dedge_gnum_tmp, PDM_part_to_block_n_elt_block_get(ptb) * sizeof(PDM_g_num_t));

  PDM_part_to_block_exch(ptb,
                        sizeof(PDM_g_num_t),
                        PDM_STRIDE_CST,
                        1,
                        NULL,
             (void **) &(results_edge_opp),
                        NULL,
              (void **) &dedge_gnum_opp);
  PDM_part_to_block_free(ptb);
  
  // Rebuild idx array (no data on external edges)
  /*int *dedge_gnum_idx = (int *) malloc ((dn_edge + 1) * sizeof(int)); //This one will be useless probably*/
  /*int count = 0;*/
  /*dedge_gnum_idx[0] = 0;*/
  /*for (int i = 0; i < dn_edge[i_zone]; i++) {*/
    /*if (i + dedge_distrib[i_zone][i_rank] + 1 == dedge_gnum[count]) {*/
      /*count++;*/
    /*}*/
    /*dedge_gnum_idx[i+1] = count;*/
  /*}*/
  log_trace("Internal edge matches after conflict resolution \n");
  PDM_log_trace_array_long(dedge_gnum, dn_internal_edge, "dedge gnum :: ");
  PDM_log_trace_array_long(dedge_gnum_opp, dn_internal_edge, "dedge gnum_opp :: ");

  log_trace("Generate dface->edge from dedge->face\n");
  int         *dface_edge_idx = NULL;
  PDM_g_num_t *dface_edge     = NULL;
  PDM_dconnectivity_transpose(comm,
                              dedge_distrib,
                              face_distri,
                              dedge_face_idx,
                              dedge_face,
                              1,
                             &dface_edge_idx,
                             &dface_edge);

  PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, ext_dn_face, "dface_edge :: ");



  //BlockToPart avec lngn == dface_edge pour aller chercher des infos des edges (dedge opp)
  PDM_g_num_t *dface_edge_abs = (PDM_g_num_t *) malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  for (int i = 0; i < dface_edge_idx[ext_dn_face]; i++)
    dface_edge_abs[i] = PDM_ABS(dface_edge[i]);
                       btp = PDM_block_to_part_create(dedge_distrib,
                               (const PDM_g_num_t **) &dface_edge_abs,
                                                      &dface_edge_idx[ext_dn_face],
                                                      1,
                                                      comm);

  //Prepare data to send on edge block
  idx = 0;
  int *dedge_gnum_n = malloc(dn_edge*sizeof(int));
  for (int i = 0; i < dn_edge; i++) {
    if (i + dedge_distrib[i_rank] + 1 == dedge_gnum[idx]) {
      dedge_gnum_n[i] = 1;
      idx++;
    }
    else {
      dedge_gnum_n[i] = 0;
    }
  }
  //PDM_log_trace_array_int(dedge_gnum_n, dn_edge, "dedge_gnum_n");


  int         **recv_stride_tmp = NULL;
  PDM_g_num_t **recv_data_tmp   = NULL;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          dedge_gnum_n,
                          dedge_gnum_opp,
                         &recv_stride_tmp,
              (void ***) &recv_data_tmp);
  int         *pedge_gnum_n   = recv_stride_tmp[0];
  PDM_g_num_t *pedge_gnum_opp = recv_data_tmp[0];
  free(recv_stride_tmp);
  free(recv_data_tmp);

  int stride2 = 2;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                         &stride2,
                          dedge_vtx,
                          NULL,
              (void ***) &recv_data_tmp);
  PDM_g_num_t *pedge_vtx = recv_data_tmp[0];
  free(recv_data_tmp);
  PDM_log_trace_array_long(pedge_vtx, 2*dface_edge_idx[ext_dn_face], "pedge_vtx");



  PDM_g_num_t *face_edge_wopp = malloc(dface_edge_idx[ext_dn_face]*sizeof(PDM_g_num_t));
  idx = 0;
  for (int i = 0; i < dface_edge_idx[ext_dn_face]; i++) {
    if (pedge_gnum_n[i] == 1)
      face_edge_wopp[i] = pedge_gnum_opp[idx++];
    else
      face_edge_wopp[i] = 0;
  }

  PDM_log_trace_connectivity_long(dface_edge_idx, face_edge_wopp, ext_dn_face, "dface_edge :: ");

  free(dface_edge_abs);
  free(dedge_gnum_n);
  free(pedge_gnum_n);
  free(pedge_gnum_opp);
  PDM_block_to_part_free(btp);

  PDM_g_num_t *p_all_edge_gnum     = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_edge_gnum_opp = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx      = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_opp  = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  PDM_g_num_t *p_all_vtx_group  = malloc(dface_edge_idx[ext_dn_face] * sizeof(PDM_g_num_t));
  // Avec la construction des faces de bord, on a des paires faces / face opp
  int glob_idx = 0;
  for (int i_face = 0; i_face < ext_dn_face/2; i_face++) {
    log_trace("\niface %d :\n", i_face);
    int face_len = dface_edge_idx[2*i_face+1] - dface_edge_idx[2*i_face];
    log_trace("face data\n");
    PDM_log_trace_array_long(&face_vtx[face_vtx_idx[i_face]], face_len, "face_vtx");
    PDM_log_trace_array_long(&dface_edge[dface_edge_idx[2*i_face]], face_len, "face_edge");
    PDM_log_trace_array_long(&pedge_vtx[2*dface_edge_idx[2*i_face]], 2*face_len, "edge vertices");
    PDM_log_trace_array_long(&face_edge_wopp[dface_edge_idx[2*i_face]], face_len, "face_edge_wopp");
    log_trace("face opp data\n");
    PDM_log_trace_array_long(&face_vtx_opp[face_vtx_idx[i_face]], face_len, "face_vtx");
    PDM_log_trace_array_long(&dface_edge[dface_edge_idx[2*i_face+1]], face_len, "face_opp_edge_opp");
    PDM_log_trace_array_long(&pedge_vtx[2*dface_edge_idx[2*i_face+1]], 2*face_len, "edge vertices opp");
    // Search any received edge (we should have at least one)
    PDM_g_num_t opp_edge_key = 0;
    int idx;
    int sign;
    for (idx = 0; idx < face_len; idx++) {
      opp_edge_key = face_edge_wopp[dface_edge_idx[2*i_face] + idx];
      sign = PDM_SIGN(dface_edge[dface_edge_idx[2*i_face]+idx]);
      if (opp_edge_key != 0)
        break;
    }
    PDM_g_num_t edge_key = PDM_ABS(dface_edge[dface_edge_idx[2*i_face]+idx]);
    //Search idx of opposite edge in opposite 
    int opp_idx;
    int opp_sign;
    for (opp_idx = 0; opp_idx < face_len; opp_idx++) {
      int candidate = dface_edge[dface_edge_idx[2*i_face+1]+opp_idx];
      opp_sign = PDM_SIGN(candidate);
      candidate = PDM_ABS(candidate);
      if (candidate == opp_edge_key)
        break;
    }

    log_trace("Reference edge will be %d (at position %d) with sign is %d -- Opp edge is %d (at position %d) with sign %d\n", edge_key, idx, sign, opp_edge_key, opp_idx, opp_sign);
    
    // Reorder edges & opp_edges
    PDM_g_num_t *ordered_edge     = malloc(face_len * sizeof(PDM_g_num_t));
    PDM_g_num_t *ordered_edge_opp = malloc(face_len * sizeof(PDM_g_num_t));
    PDM_g_num_t *ordered_vtx      = malloc(face_len * sizeof(PDM_g_num_t));
    PDM_g_num_t *ordered_vtx_opp  = malloc(face_len * sizeof(PDM_g_num_t));

    int next_vtx, next_vtx_opp, cur_vtx, cur_vtx_opp;
    if (sign == 1) {
      cur_vtx = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*idx];
      next_vtx = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*idx+1];
    }
    else {
      cur_vtx = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*idx+1];
      next_vtx = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*idx];
    }
    if (sign == 1) {
      cur_vtx_opp  = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*opp_idx+1];
      next_vtx_opp = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*opp_idx];
    }
    else {
      cur_vtx_opp  = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*opp_idx];
      next_vtx_opp = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*opp_idx+1]; //Invert looping order for opposite face
    }
    for (int i = 0; i < face_len; i++) {
      //Fill
      log_trace("Cur vtx %d and opp %d \n", cur_vtx, cur_vtx_opp);
      ordered_edge[i] = PDM_ABS(dface_edge[dface_edge_idx[2*i_face]+idx]);
      ordered_edge_opp[i] = PDM_ABS(dface_edge[dface_edge_idx[2*i_face+1]+opp_idx]);
      ordered_vtx[i]     = cur_vtx;
      ordered_vtx_opp[i] = cur_vtx_opp;
      //This is for face
      for (int j = 0; j < face_len; j++) {
        if (j != idx) {
          int vtx1 = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*j];
          int vtx2 = pedge_vtx[2*dface_edge_idx[2*i_face] + 2*j+1];
          if (vtx1 == next_vtx) {
            idx = j;
            cur_vtx  = vtx1;
            next_vtx = vtx2;
            break;
          }
          else if (vtx2 == next_vtx) {
            idx = j;
            cur_vtx  = vtx2;
            next_vtx = vtx1;
            break;
          }
        }
      }
      //This is for opposite face
      for (int j = 0; j < face_len; j++) {
        if (j != opp_idx) {
          int vtx1 = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*j];
          int vtx2 = pedge_vtx[2*dface_edge_idx[2*i_face+1] + 2*j+1];
          if (vtx1 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx1;
            next_vtx_opp = vtx2;
            break;
          }
          else if (vtx2 == next_vtx_opp) {
            opp_idx = j;
            cur_vtx_opp  = vtx2;
            next_vtx_opp = vtx1;
            break;
          }
        }
      }
      //ordered_edge_opp[i] = PDM_ABS(dface_edge[dface_edge_idx[2*i_face+1]+opp_idx]);
      //Search next edge indices
    }
    PDM_log_trace_array_long(ordered_edge,     face_len, "ordered edges");
    PDM_log_trace_array_long(ordered_edge_opp, face_len, "ordered edges_opp");
    PDM_log_trace_array_long(ordered_vtx,     face_len, "ordered vtx");
    PDM_log_trace_array_long(ordered_vtx_opp, face_len, "ordered vtx_opp");
    for (int i = 0; i < face_len; i++) {
      p_all_edge_gnum[glob_idx]     = ordered_edge[i];
      p_all_edge_gnum_opp[glob_idx] = ordered_edge_opp[i];
      p_all_vtx[glob_idx]     = ordered_vtx[i];
      p_all_vtx_opp[glob_idx] = ordered_vtx_opp[i];
      p_all_vtx_group[glob_idx] = dextract_face_group_idNEW[2*i_face];


      p_all_edge_gnum[glob_idx + (dface_edge_idx[ext_dn_face]/2)]     = ordered_edge_opp[i];
      p_all_edge_gnum_opp[glob_idx + (dface_edge_idx[ext_dn_face]/2)] = ordered_edge[i];
      p_all_vtx[glob_idx + (dface_edge_idx[ext_dn_face]/2)]     = ordered_vtx_opp[i];
      p_all_vtx_opp[glob_idx + (dface_edge_idx[ext_dn_face]/2)] = ordered_vtx[i];
      p_all_vtx_group[glob_idx + (dface_edge_idx[ext_dn_face]/2)] = dextract_face_group_idNEW[2*i_face+1];
      glob_idx++;
    }
    free(ordered_edge);
    free(ordered_edge_opp);
    free(ordered_vtx);
    free(ordered_vtx_opp);
  }
  log_trace("edge matching on face distribution\n");
  PDM_log_trace_array_long(p_all_edge_gnum,     dface_edge_idx[ext_dn_face], "p_all_edge_gnum     ::");
  PDM_log_trace_array_long(p_all_edge_gnum_opp, dface_edge_idx[ext_dn_face], "p_all_edge_gnum_opp ::");
  PDM_log_trace_array_long(p_all_vtx,     dface_edge_idx[ext_dn_face], "p_all_vtx     ::");
  PDM_log_trace_array_long(p_all_vtx_opp, dface_edge_idx[ext_dn_face], "p_all_vtx_opp ::");

  //Finally, send this back to edge distribution
  stride_one = PDM_array_const_int(dface_edge_idx[ext_dn_face], 1);
                       ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                       PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                       1.,
                                                      &p_all_edge_gnum,
                                                       dedge_distrib,
                                                      &(dface_edge_idx[ext_dn_face]),
                                                       1,
                                                       comm);
  unused_recv_stride = NULL;
  PDM_g_num_t *dall_edge_gnum_opp = NULL;

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_edge_gnum_opp,
                                     &unused_recv_stride,
                           (void **) &dall_edge_gnum_opp);
  assert (exch_size == dn_edge);
  free(unused_recv_stride);


  PDM_log_trace_array_long(dall_edge_gnum_opp, exch_size, "recv dall_edge_gnum_opp");
  PDM_part_to_block_free(ptb);
  free(stride_one);


  //We have everything on edges. We can build edge - edgeopp array
  int *all_order    = (int *) malloc(dedge_face_idx[dn_edge]*sizeof(int));
  int *group_id_tmp = (int *) malloc(dedge_face_idx[dn_edge]*sizeof(int));
  memcpy(group_id_tmp, dedge_face_group_id, dedge_face_idx[dn_edge]*sizeof(int));

  for (int i_edge =  0; i_edge < dn_edge; i_edge++) {
    int n_connected_faces = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    int start_idx = dedge_face_idx[i_edge];
    //order array must be initialized
    for (int j = 0; j < n_connected_faces; j++)
      all_order[start_idx + j] = j;
    PDM_sort_int(&group_id_tmp[start_idx], &all_order[start_idx], n_connected_faces);
  }
  free(group_id_tmp);

  //First pass to count
  int *dedge_group_n   = PDM_array_zeros_int(n_group_join);
  for (int i_edge =  0; i_edge < dn_edge; i_edge++) {
    int n_connected_faces = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    int start_idx = dedge_face_idx[i_edge];
    int *_dedge_face_group_id = &dedge_face_group_id[start_idx];
    int *loc_order = &all_order[start_idx];

    dedge_group_n[_dedge_face_group_id[loc_order[0]]]++; //Edge should be related to at least one face, so its ok
    for (int i = 1; i < n_connected_faces; i++) {
      if (_dedge_face_group_id[loc_order[i]] != _dedge_face_group_id[loc_order[i-1]])
        dedge_group_n[_dedge_face_group_id[loc_order[i]]]++;
    }
  }

  //Second pass to fill
  int *dedge_group_idx =  PDM_array_new_idx_from_sizes_int(dedge_group_n, n_group_join);
  PDM_g_num_t *dedge_group     = malloc(dedge_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_g_num_t *dedge_group_opp = malloc(dedge_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_array_reset_int(dedge_group_n, n_group_join, 0); // reset to count

  for (int i_edge =  0; i_edge < dn_edge; i_edge++) {
    int n_connected_faces = dedge_face_idx[i_edge+1] - dedge_face_idx[i_edge];
    int start_idx = dedge_face_idx[i_edge];
    int *_dedge_face_group_id = &dedge_face_group_id[start_idx];
    int *loc_order = &all_order[start_idx];

    int i_group = _dedge_face_group_id[loc_order[0]]; //Edge should be related to at least one face, so its ok
    dedge_group[dedge_group_idx[i_group] + dedge_group_n[i_group]] = dedge_distrib[i_rank] + i_edge + 1;
    dedge_group_opp[dedge_group_idx[i_group] + dedge_group_n[i_group]] = dall_edge_gnum_opp[i_edge];
    dedge_group_n[i_group]++;
    for (int i = 1; i < n_connected_faces; i++) {
      if (_dedge_face_group_id[loc_order[i]] != _dedge_face_group_id[loc_order[i-1]]) {
        i_group = _dedge_face_group_id[loc_order[i]];
        dedge_group[dedge_group_idx[i_group] + dedge_group_n[i_group]] = dedge_distrib[i_rank] + i_edge + 1;
        dedge_group_opp[dedge_group_idx[i_group] + dedge_group_n[i_group]] = dall_edge_gnum_opp[i_edge];//
        dedge_group_n[i_group]++;
      }
    }
  }
  free(dedge_group_n);
  free(all_order);

  log_trace("Edge & Edge opp per group in edge global numbering\n");
  PDM_log_trace_array_int(dedge_group_idx, n_group_join+1, "dedge_group_idx ::");
  PDM_log_trace_array_long(dedge_group, dedge_group_idx[n_group_join], "dedge_group     ::");
  PDM_log_trace_array_long(dedge_group_opp, dedge_group_idx[n_group_join], "dedge_group_opp ::");


  //Todo : we could shift back to position of vtx in extraction to have a better
  //balance of edge distribution
  stride_one = PDM_array_const_int(dface_edge_idx[ext_dn_face], 1);
                       ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                     &p_all_vtx,
                                                      NULL,
                                                     &(dface_edge_idx[ext_dn_face]),
                                                      1,
                                                      comm);

  PDM_g_num_t *vtx_distriF = PDM_part_to_block_distrib_index_get(ptb);
  blk_size = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *dall_vtx   = PDM_part_to_block_block_gnum_get(ptb);

  int *recv_stride = NULL;
  PDM_g_num_t *dall_vtx_opp   = NULL;
  PDM_g_num_t *dall_vtx_group = NULL;
  //For vertex we can do the send directly to a (new) vertex distribution
  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_opp,
                                     &recv_stride,
                           (void **) &dall_vtx_opp);

  exch_size = PDM_part_to_block_exch(ptb,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR,
                                     -1,
                           (int  **) &stride_one,
                           (void **) &p_all_vtx_group,
                                     &unused_recv_stride, //Same 
                           (void **) &dall_vtx_group);
  free(unused_recv_stride);

  PDM_log_trace_array_long(vtx_distriF     , n_rank+1, "vtx_distri F     ");
  PDM_log_trace_array_long(dall_vtx     , blk_size, "dall_vtx     ");
  PDM_log_trace_array_long(recv_stride     , blk_size, "recv stride     ");
  PDM_log_trace_array_long(dall_vtx_opp , exch_size, "recv dall_vtx_opp ");
  PDM_log_trace_array_long(dall_vtx_group , exch_size, "recv dall_vtx_group ");
  free(stride_one);

  //Post treat vertex data
  all_order    = malloc(exch_size*sizeof(int));
  group_id_tmp = malloc(exch_size * sizeof(int));
  memcpy(group_id_tmp, dall_vtx_group, exch_size*sizeof(int));

  int start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    //order array must be initialized
    for (int j = 0; j < n_recv; j++)
      all_order[start_vtx + j] = j;
    PDM_sort_int(&group_id_tmp[start_vtx], &all_order[start_vtx], n_recv);
    start_vtx += n_recv;
  }
  free(group_id_tmp);
  PDM_log_trace_array_int(all_order, exch_size, "all order ");
  
  //First pass to count
  int *dvtx_group_n   = PDM_array_zeros_int(n_group_join);
  start_vtx = 0;
  for (int i_vtx = 0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    int *_dall_vtx_group = &dall_vtx_group[start_vtx];
    int *loc_order = &all_order[start_vtx];

    dvtx_group_n[_dall_vtx_group[loc_order[0]]]++; //Vtx should be related to at least one face, so its ok
    for (int i = 1; i < n_recv; i++) {
      if (_dall_vtx_group[loc_order[i]] != _dall_vtx_group[loc_order[i-1]])
        dvtx_group_n[_dall_vtx_group[loc_order[i]]]++;
    }
    start_vtx += n_recv;
  }
  //Second pass to fill
  int *dvtx_group_idx =  PDM_array_new_idx_from_sizes_int(dvtx_group_n, n_group_join);
  PDM_g_num_t *dvtx_group     = malloc(dvtx_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_g_num_t *dvtx_group_opp = malloc(dvtx_group_idx[n_group_join]*sizeof(PDM_g_num_t));
  PDM_array_reset_int(dvtx_group_n, n_group_join, 0); // reset to count

  start_vtx = 0;
  for (int i_vtx =  0; i_vtx < blk_size; i_vtx++) {
    int n_recv = recv_stride[i_vtx];
    int *_dall_vtx_group = &dall_vtx_group[start_vtx];
    PDM_g_num_t opp_vtx_gnum = dall_vtx_opp[start_vtx]; //Will be the same even if vertex appears more than once
    int *loc_order = &all_order[start_vtx];

    int i_group = _dall_vtx_group[loc_order[0]];
    dvtx_group[dvtx_group_idx[i_group] + dvtx_group_n[i_group]]     = dall_vtx[i_vtx];
    dvtx_group_opp[dvtx_group_idx[i_group] + dvtx_group_n[i_group]] = opp_vtx_gnum;
    dvtx_group_n[i_group]++;
    for (int i = 1; i < n_recv; i++) {
      if (_dall_vtx_group[loc_order[i]] != _dall_vtx_group[loc_order[i-1]]) {
        i_group = _dall_vtx_group[loc_order[i]];
        dvtx_group[dvtx_group_idx[i_group] + dvtx_group_n[i_group]] = dall_vtx[i_vtx];
        dvtx_group_opp[dvtx_group_idx[i_group] + dvtx_group_n[i_group]] = opp_vtx_gnum;
        dvtx_group_n[i_group]++;
      }
    }
    start_vtx += n_recv;
  }
  free(dvtx_group_n);
  free(all_order);
  PDM_part_to_block_free(ptb);

  //Ultimate step : go back to original vtx numbering. All we have to do is retrieve zone
  // and substract zone offset  for (int i_vtx = 0; i_vtx < dvtx_group_idx[n_group_join]; i_vtx++) {
  for (int i_join = 0; i_join < n_group_join; i_join++) {
    int zone_cur_offset = vtx_per_block_offset[group_join_to_zone_cur[i_join]];
    int zone_opp_offset = vtx_per_block_offset[group_join_to_zone_opp[i_join]];
    for (int i_vtx = dvtx_group_idx[i_join]; i_vtx < dvtx_group_idx[i_join+1]; i_vtx++) {
      dvtx_group[i_vtx]     -= zone_cur_offset;
      dvtx_group_opp[i_vtx] -= zone_opp_offset;
    }
  }
  log_trace("Vtx & vtx opp per group in vtx global numbering\n");
  PDM_log_trace_array_int(dvtx_group_idx, n_group_join+1, "dvtx_group_idx ::");
  PDM_log_trace_array_long(dvtx_group, dvtx_group_idx[n_group_join], "dvtx_group     ::");
  PDM_log_trace_array_long(dvtx_group_opp, dvtx_group_idx[n_group_join], "dvtx_group_opp ::");



  free(dface_edge_idx);
  free(dface_edge);

  free(dextract_face_joinNEW);
  free(dextract_face_join_oppNEW);
  free(dextract_face_group_idNEW);

  free(dedge_face_abs);
  free(dedge_face_sgn);

  free(dedge_face_join);
  free(dedge_face_join_opp);
  free(dedge_face_group_id);


  free(face_distri);

  free(face_vtx_both_n);
  free(face_vtx_both_idx);
  free(face_vtx_both);

  free(face_per_block_offset);
  free(vtx_per_block_offset);
  free(multi_gnum);

  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(all_face_distribution[i_zone]);
    free(face_vtx_n[i_zone]);
    free(face_vtx_shifted[i_zone]);
  }
  free(all_face_distribution);
  free(face_vtx_n);
  free(face_vtx_shifted);

  free(face_vtx_idx);
  free(face_vtx);
  free(face_vtx_opp);

  free(dedge_distrib);
  free(dedge_vtx_idx);
  free(dedge_vtx);
  free(dedge_face_idx);
  free(dedge_face);


  /* Free memory */
  for (int i_zone = 0; i_zone < n_zone; i_zone++) {
    free(dface_bnd_idx [i_zone]);
    free(dface_bnd     [i_zone]);
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
  free(dface_join_opp);

  free(group_join_to_zone_cur);
  free(group_join_to_zone_opp);
  free(group_join_to_join_opp);




  PDM_MPI_Finalize();

  return 0;
}
