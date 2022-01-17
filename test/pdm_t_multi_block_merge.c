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
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_priv.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_distrib.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multi_block_merge.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"

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
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of vertices on the cube side.\n\n"
     "  -l      <level>  Cube length.\n\n"
     "  -n_part <level>  Number of partitions par process.\n\n"
     "  -post            Ensight outputs (only if n_part == 1). \n\n"
     "  -parmetis        Call ParMETIS.\n\n"
     "  -pt-scocth       Call PT-Scotch.\n\n"
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

int interface_to_graph
(
  const int      n_interface,
  int           *interface_dn,
  PDM_g_num_t  **interface_ids,
  int          **interface_dom,
  int          **graph_idx,
  PDM_g_num_t  **graph_ids,
  int          **graph_dom,
  PDM_MPI_Comm   comm
)
{
  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);

  // Step 0 : retrieve some data. We need (to offset gnums)
  //   - the number of involved blocks
  //   - the max id occuring in each block
  int n_zone = -1;
  int max_zone_loc = 0;
  for (int itrf = 0; itrf < n_interface; itrf++) {
    for (int k = 0; k < 2*interface_dn[itrf]; k++) {
      max_zone_loc = PDM_MAX(max_zone_loc, interface_dom[itrf][k]);
    }
  }
  PDM_MPI_Allreduce(&max_zone_loc, &n_zone, 1, PDM_MPI_INT, PDM_MPI_MAX, comm);
  n_zone++; //Because zone numbering start at 0
  log_trace("max zone loc %d\n", n_zone);

  PDM_g_num_t *max_per_zone_loc = PDM_array_const_gnum(n_zone, 0);
  PDM_g_num_t *max_per_zone     = malloc((n_zone+1) * sizeof(PDM_g_num_t));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    for (int k = 0; k < interface_dn[itrf]; k++) {
      int dom    = interface_dom[itrf][2*k];
      int domopp = interface_dom[itrf][2*k+1];
      max_per_zone_loc[dom] = PDM_MAX(max_per_zone_loc[dom], interface_ids[itrf][2*k]);
      max_per_zone_loc[domopp] = PDM_MAX(max_per_zone_loc[domopp], interface_ids[itrf][2*k+1]);
    }
  }
  PDM_log_trace_array_long(max_per_zone_loc, n_zone, "max per zone loc");
  max_per_zone[0] = 0;
  PDM_MPI_Allreduce(max_per_zone_loc, &max_per_zone[1], n_zone, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, comm);
  PDM_log_trace_array_long(max_per_zone, n_zone+1, "max per itrf");
  PDM_array_accumulate_gnum(max_per_zone, n_zone+1);
  PDM_log_trace_array_long(max_per_zone, n_zone+1, "max per zone");
  free(max_per_zone_loc);
  
  // Prepare first PtB with multiple partitions.
  // Use (shifted) ids as gnum and send tuple (shited) id, opp_id
  PDM_g_num_t **interface_ids_shifted = malloc(n_interface * sizeof(PDM_g_num_t*));
  PDM_g_num_t **send_data             = malloc(n_interface * sizeof(PDM_g_num_t*));
  double      **weight                = malloc(n_interface * sizeof(double*     ));
  int         **stride_one            = malloc(n_interface * sizeof(int*        ));
  PDM_g_num_t  *interface_dn_twice    = malloc(n_interface * sizeof(PDM_g_num_t));
  for (int itrf = 0; itrf < n_interface; itrf++) {
    stride_one[itrf]            = malloc(2*interface_dn[itrf]*sizeof(int        ));
    interface_ids_shifted[itrf] = malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    send_data[itrf]             = malloc(2*interface_dn[itrf]*sizeof(PDM_g_num_t));
    weight[itrf]                = malloc(2*interface_dn[itrf]*sizeof(double     ));
    interface_dn_twice[itrf]    = 2*interface_dn[itrf];
    for (int k = 0; k < interface_dn[itrf]; k++) {
      interface_ids_shifted[itrf][2*k]   = interface_ids[itrf][2*k] + max_per_zone[interface_dom[itrf][2*k]];
      interface_ids_shifted[itrf][2*k+1] = interface_ids[itrf][2*k+1] + max_per_zone[interface_dom[itrf][2*k+1]];
      send_data[itrf][2*k]   = interface_ids[itrf][2*k+1] + max_per_zone[interface_dom[itrf][2*k+1]];
      send_data[itrf][2*k+1] = interface_ids[itrf][2*k] + max_per_zone[interface_dom[itrf][2*k]];
      weight[itrf][2*k]   = 1.;
      weight[itrf][2*k+1] = 1.;
      stride_one[itrf][2*k]   = 1;
      stride_one[itrf][2*k+1] = 1;
    }
    PDM_log_trace_array_long(interface_ids_shifted[itrf], 2*interface_dn[itrf], "shifted");
    PDM_log_trace_array_long(send_data[itrf], 2*interface_dn[itrf], "send");
  }
  
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      interface_ids_shifted,
                                                      weight,
                                                      interface_dn_twice,
                                                      n_interface,
                                                      comm);
  // Save distribution & gnum from first PtB. We will use it for following PtBs
  int n_gnum = PDM_part_to_block_n_elt_block_get(ptb);
  PDM_g_num_t *distri = malloc((n_rank+1) * sizeof(PDM_g_num_t));
  PDM_g_num_t *gnum   = malloc(n_gnum    * sizeof(PDM_g_num_t));
  memcpy(gnum,   PDM_part_to_block_block_gnum_get(ptb),        n_gnum*sizeof(PDM_g_num_t));
  memcpy(distri, PDM_part_to_block_distrib_index_get(ptb), (n_rank+1)*sizeof(PDM_g_num_t));

  int         *recv_stride = NULL;
  PDM_g_num_t *recv_data   = NULL;
  int n_connected_l = PDM_part_to_block_exch(ptb,
                                             sizeof(PDM_g_num_t),
                                             PDM_STRIDE_VAR,
                                             -1,
                                             stride_one,
                                   (void **) send_data,
                                             &recv_stride,
                                   (void **) &recv_data);
  log_trace("n recv %d\n", n_connected_l);
  PDM_log_trace_array_int(recv_stride, n_gnum, "recv stride");
  PDM_log_trace_array_long(PDM_part_to_block_block_gnum_get(ptb), n_gnum, "gnum");
  PDM_log_trace_array_long(recv_data, n_connected_l, "recv data");

  PDM_part_to_block_free(ptb);
  for (int itrf = 0; itrf < n_interface; itrf++) {
    free(stride_one           [itrf]);
    free(send_data            [itrf]);
    free(weight               [itrf]);
    free(interface_ids_shifted[itrf]);
  }
  free(stride_one);
  free(weight);
  free(send_data);
  free(interface_dn_twice);
  free(interface_ids_shifted);

  int n_connected;
  PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  log_trace("Total size of graph : %d \n", n_connected);


  /* After first exchange, we received for each gnum a list (usually of size one) of connected
   * gnum. The idea is to consider this received list as a partition, and to send to each
   * entity of this partition : the neighbors in the received list (if stride was > 1) + the original gnum.
   *
   * After some iterations this will group the related ids together
  */

  int n_connected_prev = 0;
  while(n_connected_prev != n_connected) {
    log_trace("\nStart it\n");
    n_connected_prev = n_connected;

    ptb = PDM_part_to_block_create2(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                    PDM_PART_TO_BLOCK_POST_MERGE,
                                    1.,
                                   &recv_data,
                                    distri,
                                   &n_connected_l,
                                    1,
                                    comm);

    
    int *send_stride = malloc(n_connected_l*sizeof(int));
    int w_idx = 0;
    int n_data = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        send_stride[w_idx++] = recv_stride[k];
        n_data += recv_stride[k];
      }
    }
    assert (w_idx == n_connected_l);
    PDM_g_num_t *send_data = malloc(n_data*sizeof(PDM_g_num_t));
    w_idx = 0;
    int r_idx = 0;
    for (int k = 0; k < n_gnum; k++) {
      for (int j = 0; j < recv_stride[k]; j++) {
        //Add gnum
        send_data[w_idx++] = gnum[k];
        //Add others
        for (int i = 0; i < j; i++)
          send_data[w_idx++] = recv_data[r_idx + i];
        for (int i = j+1; i < recv_stride[k]; i++)
          send_data[w_idx++] = recv_data[r_idx + i];
      }
      r_idx += recv_stride[k];
    }
    assert (r_idx == n_connected_l);
    assert (w_idx == n_data);
    PDM_log_trace_array_int(send_stride, n_connected_l, "send stride");
    PDM_log_trace_array_long(send_data, n_data, "send_data");

    int         *recv_stride_next = NULL;
    PDM_g_num_t *recv_data_next   = NULL;
    PDM_part_to_block_exch(ptb,
                           sizeof(PDM_g_num_t),
                           PDM_STRIDE_VAR,
                           -1,
                           &send_stride,
                 (void **) &send_data,
                           &recv_stride_next,
                 (void **) &recv_data_next);

    free(send_stride);
    free(send_data);
    PDM_part_to_block_free(ptb);
    
    // Post treat recv data to remove duplicated per gnum and count size of graph
    int start = 0;
    n_connected_l = 0;
    for (int i = 0; i < n_gnum; i ++) {
      int n_unique = PDM_inplace_unique_long(recv_data_next, NULL, start, start+recv_stride_next[i]-1);
      //Compress array at the same time (let meanless data at the end of array)
      memcpy(&recv_data_next[n_connected_l], &recv_data_next[start], n_unique*sizeof(PDM_g_num_t));
      start += recv_stride_next[i];
      recv_stride_next[i] = n_unique;
      n_connected_l += n_unique;
    }

    PDM_log_trace_array_long(gnum, n_gnum, "gnum");
    PDM_log_trace_array_int(recv_stride_next, n_gnum, "recv stride");
    PDM_log_trace_array_long(recv_data_next, n_connected_l, "recv data");
    PDM_MPI_Allreduce(&n_connected_l, &n_connected, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
    log_trace("Total size of graph : %d \n", n_connected);

    free(recv_stride);
    free(recv_data); // To free after PTB because it was used as lngn
    recv_data = recv_data_next;
    recv_stride = recv_stride_next;
  }
  
  // When iteration are completed, all the connections are known by every id.
  // Last step is to compress the graph and to redistribute it
  // To do that we take for each group of related id the min of it as lngn
  int n_keys = 0;
  n_connected_l = 0;
  int r_idx = 0;
  int *is_key_gr = PDM_array_const_int(n_gnum, 1);
  for (int k = 0; k < n_gnum; k++) {
    for (int j = 0; j < recv_stride[k]; j++) {
      if (recv_data[r_idx+j] < gnum[k]) {
        is_key_gr[k] = 0;
        break;
      }
    }
    if (is_key_gr[k])
      n_connected_l += recv_stride[k] + 1;
    r_idx  += recv_stride[k];
    n_keys += is_key_gr[k];
  }

  PDM_g_num_t *lngn_gr        = malloc(n_keys*sizeof(PDM_g_num_t));
  int         *send_stride_gr = malloc(n_keys*sizeof(int        ));
  double      *weight_gr      = malloc(n_keys*sizeof(double     ));
  PDM_g_num_t *send_data_gr   = malloc(n_connected_l *sizeof(PDM_g_num_t));
  int w_idx = 0;
  int w_idx2 = 0;
  r_idx = 0;
  for (int k = 0; k < n_gnum; k++) {
    if (is_key_gr[k]) {
      lngn_gr[w_idx]        = gnum[k];
      send_stride_gr[w_idx] = recv_stride[k] + 1; //Include gnum in send data so we have directly graph
      weight_gr[w_idx]      = (double) (recv_stride[k] + 1);
      w_idx++;
      send_data_gr[w_idx2++] = gnum[k];
      memcpy(&send_data_gr[w_idx2], &recv_data[r_idx], recv_stride[k]*sizeof(PDM_g_num_t));
      w_idx2 += recv_stride[k];
    }
    r_idx += recv_stride[k];
  }
  log_trace("\n");
  PDM_log_trace_array_long(lngn_gr, n_keys, "keys");
  PDM_log_trace_array_int(send_stride_gr, n_keys, "send stride last");
  PDM_log_trace_array_long(send_data_gr, n_connected_l, "send data last");
  
  // Data of previous iteration is not usefull anymore
  free(gnum);
  free(distri);
  free(recv_stride);
  free(recv_data);
  
  //TODO In fact we just want to do a block to block, but PTB + weights compute distribution for us
  //and we are lazy
  ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                 PDM_PART_TO_BLOCK_POST_MERGE,
                                 1.,
                                &lngn_gr,
                                &weight_gr,
                                &n_keys,
                                 1,
                                 comm);
  int graph_dn = PDM_part_to_block_n_elt_block_get(ptb);

  int         *graph_size = NULL;
  PDM_g_num_t *graph_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR,
                         -1,
                         &send_stride_gr,
               (void **) &send_data_gr,
                         &graph_size,
               (void **) &graph_gnum);

  PDM_part_to_block_free(ptb);
  free(send_stride_gr);
  free(send_data_gr);
  free(weight_gr);
  free(lngn_gr);
  free(is_key_gr);

  int* _graph_idx = PDM_array_new_idx_from_sizes_int(graph_size, graph_dn);
  int *_graph_dom = malloc(_graph_idx[graph_dn]*sizeof(int));

  // Retrieve domain and local gnum in domain
  for (int i = 0; i < _graph_idx[graph_dn]; i++) {
    _graph_dom[i] = PDM_binary_search_gap_long(graph_gnum[i]-1, max_per_zone, n_zone+1);
    graph_gnum[i] -= max_per_zone[_graph_dom[i]];
  }
  free(graph_size);
  free(max_per_zone);

  *graph_idx = _graph_idx;
  *graph_ids =  graph_gnum;
  *graph_dom = _graph_dom;
  return graph_dn;
}


/**
 *
 * \brief  Main
 *
 */
int main(int argc, char *argv[])
{

  /*
   *  Set default values
   */

  PDM_g_num_t        n_vtx_seg = 4;
  double             length  = 1.;

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);
  /*
   *  Init
   */
  int i_rank;
  int n_rank;

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_dcube_nodal_t* dcube1 = PDM_dcube_nodal_gen_create(comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        length,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube1, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube1);

  PDM_dmesh_nodal_t*  dmn1 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube1);
  /*
   * Define distribution of cell
   */
  PDM_dmesh_nodal_dump_vtk(dmn1, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_dcube1_");

  /*
   * Define distibution of vtx
   */
  PDM_dcube_nodal_t* dcube2 = PDM_dcube_nodal_gen_create(comm,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         n_vtx_seg,
                                                         length,
                                                         1.,
                                                         0.,
                                                         0.,
                                                         PDM_MESH_NODAL_QUAD4,
                                                         1,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_dcube_nodal_gen_ordering_set (dcube2, "PDM_HO_ORDERING_CGNS");
  PDM_dcube_nodal_gen_build (dcube2);

  PDM_dmesh_nodal_t*  dmn2 = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube2);

  PDM_dmesh_nodal_dump_vtk(dmn2, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic_dcube2_");

  /*
   * Concatenate blocks
   */
  int n_block = 2;
  PDM_g_num_t** block_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_distrib_idx[0] = PDM_dmesh_nodal_vtx_distrib_get(dmn1);
  block_distrib_idx[1] = PDM_dmesh_nodal_vtx_distrib_get(dmn2);

  int* n_selected = malloc(n_block * sizeof(int));
  n_selected[0] = PDM_DMesh_nodal_n_vtx_get(dmn1);
  n_selected[1] = PDM_DMesh_nodal_n_vtx_get(dmn2);

  PDM_g_num_t** selected_g_num = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    selected_g_num[i_block] = malloc(n_selected[i_block] * sizeof(PDM_g_num_t));
    PDM_log_trace_array_long(block_distrib_idx[i_block], n_rank+1, "block_distrib_idx ::");
    for(int i = 0; i < n_selected[i_block]; ++i) {
      selected_g_num[i_block][i] = block_distrib_idx[i_block][i_rank] + i + 1;
    }
  }

  /*
   * Setup graph
   */
  int         **dmerge_idx      = malloc(n_block * sizeof(int         *));
  int         **dmerge_block_id = malloc(n_block * sizeof(int         *));
  PDM_g_num_t **dmerge_g_num    = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    int dn_vtx = block_distrib_idx[i_block][i_rank+1] - block_distrib_idx[i_block][i_rank];

    dmerge_idx     [i_block] = malloc( ( dn_vtx + 1 ) * sizeof(int        ));
    dmerge_block_id[i_block] = malloc( ( dn_vtx     ) * sizeof(int        ));
    dmerge_g_num   [i_block] = malloc( ( dn_vtx     ) * sizeof(PDM_g_num_t));

    dmerge_idx     [i_block][0] = 0;

    // PDM_g_num_t vtx_g_num_next = 1;
    PDM_g_num_t step = (block_distrib_idx[i_block][i_rank] + 1) / (n_vtx_seg + 1);
    PDM_g_num_t rem  = (block_distrib_idx[i_block][i_rank] + 1) % (n_vtx_seg + 1);
    // PDM_g_num_t vtx_g_num_next = step + 1;
    // PDM_g_num_t vtx_g_num_next = step*(n_vtx_seg+1) + rem;
    PDM_g_num_t vtx_g_num_next = step * ( n_vtx_seg + 1 ) + 1;

    // PDM_g_num_t vtx_g_num_next = (block_distrib_idx[i_block][i_rank] + 1) % ;
    if(i_block == 1) {
      vtx_g_num_next = (step+1) * (n_vtx_seg + 1) ; //+rem;
    }
    log_trace("dist = %i | step = %i | rem = %i -> %i \n", block_distrib_idx[i_block][i_rank], step, rem, vtx_g_num_next);

    int idx_write = 0;
    for(int j = 0; j < dn_vtx; ++j) {
      PDM_g_num_t vtx_g_num = block_distrib_idx[i_block][i_rank] + j + 1;
      PDM_g_num_t indi      = vtx_g_num % ( n_vtx_seg + 1 );

      dmerge_idx[i_block][j+1] = dmerge_idx[i_block][j];

      if(indi == 0 && i_block == 0) {
        dmerge_idx     [i_block][j+1] = dmerge_idx[i_block][j] + 1;
        dmerge_block_id[i_block][idx_write] = 1;
        dmerge_g_num   [i_block][idx_write] = vtx_g_num_next;
        vtx_g_num_next += n_vtx_seg+1;
        idx_write++;
      } else if(indi == 1 && i_block == 1) {
        dmerge_idx     [i_block][j+1] = dmerge_idx[i_block][j] + 1;
        dmerge_block_id[i_block][idx_write] = 0;
        dmerge_g_num   [i_block][idx_write] = vtx_g_num_next;
        vtx_g_num_next += n_vtx_seg+1;
        idx_write++;
      }

      // PDM_g_num_t g_num = n_vtx_seg * i;
      // printf("vtx_g_num = %i | indi = %i  \n", vtx_g_num, indi);
    }

    dmerge_block_id[i_block] = realloc(dmerge_block_id[i_block], idx_write * sizeof(int        ));
    dmerge_g_num   [i_block] = realloc(dmerge_g_num   [i_block], idx_write * sizeof(PDM_g_num_t));

    PDM_log_trace_array_int (dmerge_block_id[i_block], idx_write, "dmerge_block_id :: ");
    PDM_log_trace_array_long(dmerge_g_num   [i_block], idx_write, "dmerge_g_num    :: ");

  }

  PDM_multi_block_merge_t* mbm = PDM_multi_block_merge_create(block_distrib_idx,
                                                              n_block,
                                                              n_selected,
                                                              selected_g_num,
                                                              dmerge_idx,
                                                              dmerge_block_id,
                                                              dmerge_g_num,
                                                              comm);


  // Exchange

  double** dvtx_coord = malloc(n_block * sizeof(double *));
  dvtx_coord[0] = PDM_DMesh_nodal_vtx_get(dmn1);
  dvtx_coord[1] = PDM_DMesh_nodal_vtx_get(dmn2);

  int** stride_one = malloc(n_block * sizeof(int *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    stride_one[i_block] = PDM_array_const_int(1, 1);
  }

  double *dmerge_vtx_coord = NULL;
  PDM_multi_block_merge_exch(mbm,
                             3 * sizeof(double),
                             PDM_STRIDE_CST,
                             stride_one,
                 (void * )   dvtx_coord,
                             NULL,
                 (void **)   &dmerge_vtx_coord);


  //
  //
  //
  // PDM_multi_block_merge_exch(mbm_elt,
  //                            3 * sizeof(double),
  //                            PDM_STRIDE_CST,
  //                            stride_one,
  //                (void * )   dcell_vtx,
  //                            NULL,
  //                (void **)   &dmerge_dcell_vtx);


  // origin_cell get_orgin_block (size= n_dmerge_cell)

  // orgin_vtx = 4 * s_orgini_cell

  // origin = 

  // Creer dans PDM_multi_block_merge une fonction qui applique la nouvelle numerotation
  // à un tableau contenant des références à l'ancienne numerotation  
  //
  // Transformer indication numerotation en doublon / numabs origin
  //
  // PDM_multi_block_merge_apply_array(mbm,
  //                            size_dmerge_dcell_vtx,
  //                            dmerge_vtx_origi_block,
  //                            dmerge_dcell_vtx,
  //                            dmerge_dcell_new_vtx);



  free(dvtx_coord);
  exit(1);

  /*
   * Same protocol for cells
   */

  PDM_g_num_t** block_elmt_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_elmt_distrib_idx[0] = (PDM_g_num_t *) PDM_DMesh_nodal_distrib_section_get(dmn1, PDM_GEOMETRY_KIND_SURFACIC, 0);
  block_elmt_distrib_idx[1] = (PDM_g_num_t *) PDM_DMesh_nodal_distrib_section_get(dmn2, PDM_GEOMETRY_KIND_SURFACIC, 0);

  int* n_elmt_selected = malloc(n_block * sizeof(int));
  n_elmt_selected[0] = PDM_DMesh_nodal_section_n_elt_get(dmn1, PDM_GEOMETRY_KIND_SURFACIC, 0);
  n_elmt_selected[1] = PDM_DMesh_nodal_section_n_elt_get(dmn2, PDM_GEOMETRY_KIND_SURFACIC, 0);

  PDM_g_num_t** selected_elmt_g_num = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    selected_elmt_g_num[i_block] = malloc(n_elmt_selected[i_block] * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_elmt_selected[i_block]; ++i) {
      selected_elmt_g_num[i_block][i] = block_elmt_distrib_idx[i_block][i_rank] + i + 1;
    }
  }

  /*
   * Setup graph
   */
  int         **dmerge_elmt_idx      = malloc(n_block * sizeof(int         *));
  int         **dmerge_elmt_block_id = malloc(n_block * sizeof(int         *));
  PDM_g_num_t **dmerge_elmt_g_num    = malloc(n_block * sizeof(PDM_g_num_t *));
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    dmerge_elmt_idx     [i_block] = malloc( ( n_elmt_selected[i_block] + 1 ) * sizeof(int        ));
    dmerge_elmt_block_id[i_block] = malloc( ( 0 ) * sizeof(int        ));
    dmerge_elmt_g_num   [i_block] = malloc( ( 0 ) * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_elmt_selected[i_block]+1; ++i) {
      dmerge_elmt_idx[i_block][i] = 0;
    }
  }

  PDM_multi_block_merge_t* mbm_elmt = PDM_multi_block_merge_create(block_elmt_distrib_idx,
                                                                   n_block,
                                                                   n_elmt_selected,
                                                                   selected_elmt_g_num,
                                                                   dmerge_elmt_idx,
                                                                   dmerge_elmt_block_id,
                                                                   dmerge_elmt_g_num,
                                                                   comm);

  /*
   * Exchange + Update
   */

  PDM_g_num_t *old_vtx_distrib;
  int         *dold_vtx_to_new_idx;
  PDM_g_num_t *dold_vtx_to_new;
  PDM_multi_block_merge_get_old_to_new(mbm, &old_vtx_distrib, &dold_vtx_to_new_idx, &dold_vtx_to_new);

  PDM_g_num_t* multi_block_vtx_distrib = PDM_multi_block_merge_get_multi_distrib(mbm);


  PDM_g_num_t** block_elmt_vtx = malloc(n_block * sizeof(PDM_g_num_t *));
  block_elmt_vtx[0] = PDM_DMesh_nodal_section_std_get(dmn1, PDM_GEOMETRY_KIND_SURFACIC, 0);
  block_elmt_vtx[1] = PDM_DMesh_nodal_section_std_get(dmn2, PDM_GEOMETRY_KIND_SURFACIC, 0);
  PDM_g_num_t** block_elmt_shift_distrib_idx = malloc(n_block * sizeof(PDM_g_num_t *));
  int strid_cst = 4;
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    stride_one[i_block][0] = strid_cst; // Because QUAD
    block_elmt_shift_distrib_idx[i_block] = malloc( 4 * n_elmt_selected[i_block] * sizeof(PDM_g_num_t));


    for(int i = 0; i < 4 * n_elmt_selected[i_block]; ++i) {
      block_elmt_shift_distrib_idx[i_block][i] = multi_block_vtx_distrib[i_block] + block_elmt_vtx[i_block][i];
    }

  }

  PDM_g_num_t *dmerge_elmt_vtx = NULL;
  PDM_multi_block_merge_exch_and_update_child_g_num(mbm_elmt,
                                                    old_vtx_distrib,
                                                    dold_vtx_to_new_idx,
                                                    dold_vtx_to_new,
                                                    PDM_STRIDE_CST,
                                                    strid_cst,
                                                    stride_one,
                                     (void *)       block_elmt_shift_distrib_idx,
                                                    NULL,
                                     (void **)      &dmerge_elmt_vtx);

  /*
   * Visualisation
   */
  int dn_merge_vtx               = PDM_multi_block_merge_get_n_block(mbm);
  PDM_g_num_t* distrib_merge_vtx = PDM_multi_block_merge_get_distrib(mbm);

  int dn_merge_elmt               = PDM_multi_block_merge_get_n_block(mbm_elmt);
  PDM_g_num_t* distrib_merge_elmt = PDM_multi_block_merge_get_distrib(mbm_elmt);

  assert(dn_merge_vtx  == distrib_merge_vtx [i_rank+1] - distrib_merge_vtx [i_rank]);
  assert(dn_merge_elmt == distrib_merge_elmt[i_rank+1] - distrib_merge_elmt[i_rank]);

  printf("dn_merge_elmt = %i \n", dn_merge_elmt);

  PDM_g_num_t *merge_elmt_ln_to_gn = malloc( dn_merge_elmt      * sizeof(PDM_g_num_t));
  int         *dconnec_idx         = malloc((dn_merge_elmt + 1) * sizeof(int        ));

  dconnec_idx[0] = 0;
  for(int i = 0; i < dn_merge_elmt; ++i) {
    merge_elmt_ln_to_gn[i] = distrib_merge_elmt[i_rank] + i + 1;
    dconnec_idx[i+1] = dconnec_idx[i] + 4; // Because QUAD
  }

  PDM_log_trace_connectivity_long(dconnec_idx, dmerge_elmt_vtx, dn_merge_elmt, "dmerge_elmt_vtx :: ");


  PDM_g_num_t *pvtx_ln_to_gn;
  int         *pcell_vtx_idx;
  int         *pcell_vtx;
  int          pn_vtx;
  PDM_part_dconnectivity_to_pconnectivity_sort_single_part(comm,
                                                           distrib_merge_elmt,
                                                           dconnec_idx,
                                                           dmerge_elmt_vtx,
                                                           dn_merge_elmt,
                                  (const PDM_g_num_t *)    merge_elmt_ln_to_gn,
                                                           &pn_vtx,
                                                           &pvtx_ln_to_gn,
                                                           &pcell_vtx_idx,
                                                           &pcell_vtx);

  double** tmp_pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        1,
                                        distrib_merge_vtx,
                                        dmerge_vtx_coord,
                                        &pn_vtx,
                 (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                        &tmp_pvtx_coord);
  double* pvtx_coord_out = tmp_pvtx_coord[0];
  free(tmp_pvtx_coord);


  char filename_elmt[999];
  sprintf(filename_elmt, "merge_mesh_%2.2d.vtk", i_rank);
  PDM_vtk_write_std_elements(filename_elmt,
                             pn_vtx,
                             pvtx_coord_out,
                             pvtx_ln_to_gn,
                             PDM_MESH_NODAL_QUAD4,
                             dn_merge_elmt,
                             pcell_vtx,
                             merge_elmt_ln_to_gn,
                             0,
                             NULL,
                             NULL);

  PDM_g_num_t* merge_vtx_ln_to_gn = malloc(dn_merge_vtx * sizeof(PDM_g_num_t));
  for(int i = 0; i < dn_merge_vtx; ++i) {
    merge_vtx_ln_to_gn[i] = distrib_merge_vtx[i_rank] + i + 1;
  }
  char filename[999];
  sprintf(filename, "debug_dvtx_coord_merge_%2.2d.vtk", i_rank);
  PDM_vtk_write_point_cloud(filename,
                            dn_merge_vtx,
                            dmerge_vtx_coord,
                            merge_vtx_ln_to_gn,
                            NULL);

  free(dmerge_vtx_coord);
  free(merge_vtx_ln_to_gn);
  free(merge_elmt_ln_to_gn);
  free(dconnec_idx);
  free(pvtx_ln_to_gn);
  free(pcell_vtx_idx);
  free(pcell_vtx);
  free(pvtx_coord_out);

  PDM_multi_block_merge_free(mbm);
  PDM_multi_block_merge_free(mbm_elmt);

  free(dmerge_elmt_vtx);
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    free(selected_g_num[i_block]);
    free(dmerge_block_id[i_block]);
    free(dmerge_g_num[i_block]);
    free(dmerge_idx[i_block]);
    free(selected_elmt_g_num[i_block]);
    free(dmerge_elmt_block_id[i_block]);
    free(dmerge_elmt_g_num[i_block]);
    free(dmerge_elmt_idx[i_block]);
  }
  for(int i_block = 0; i_block < n_block ; ++i_block) {
    free(stride_one[i_block]);
    free(block_elmt_shift_distrib_idx[i_block]);
  }
  free(stride_one);
  free(block_elmt_vtx);
  free(block_elmt_shift_distrib_idx);
  free(selected_g_num);
  free(dmerge_idx     );
  free(dmerge_block_id);
  free(dmerge_g_num   );
  free(selected_elmt_g_num);
  free(dmerge_elmt_idx     );
  free(dmerge_elmt_block_id);
  free(dmerge_elmt_g_num   );

  free(block_distrib_idx);
  free(block_elmt_distrib_idx);
  free(n_selected);
  free(n_elmt_selected);


  PDM_dcube_nodal_gen_free(dcube1);
  PDM_dcube_nodal_gen_free(dcube2);

  PDM_MPI_Finalize ();
  return 0;
}
