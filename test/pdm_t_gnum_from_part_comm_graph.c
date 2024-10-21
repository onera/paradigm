/*============================================================================
 * Test for mesh deformation with parallel clustering
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_printf.h"
#include "pdm_binary_search.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_timer.h"
#include "pdm_distrib.h"
#include "pdm_gnum.h"
#include "pdm_part_to_part.h"
#include "pdm_part_to_block.h"
#include "pdm_para_octree.h"
#include "pdm_closest_points.h"
#include "pdm_dist_cloud_surf.h"

#include "pdm_extract_part.h"
#include "pdm_vtk.h"
#include "pdm_generate_mesh.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_logging.h"
#include "pdm_distant_neighbor.h"
#include "pdm_part_comm_graph.h"

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
_read_args(int                     argc,
           char                  **argv,
           PDM_g_num_t           *n_vtx_seg,
           double                *length)
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

/*=============================================================================
 * Public function definitions
 *============================================================================*/

static
void
_generate_gnum_old
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int           *pn_entity,
  PDM_g_num_t  **pentity_ln_to_gn,
  PDM_g_num_t ***out_pentity_ln_to_gn
)
{
  PDM_gen_gnum_t *gnum = PDM_gnum_create(3,
                                         n_part,
                                         PDM_FALSE,
                                         1.,
                                         comm,
                                         PDM_OWNERSHIP_USER);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_gnum_set_from_parents(gnum,
                              i_part,
                              pn_entity[i_part],
                              pentity_ln_to_gn[i_part]);
  }

  PDM_gnum_compute(gnum);

  *out_pentity_ln_to_gn = malloc(sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part; i_part++) {
    (*out_pentity_ln_to_gn)[i_part] = PDM_gnum_get(gnum, i_part);
  }

  PDM_gnum_free(gnum);
}


static
void
_generate_gnum
(
  PDM_MPI_Comm   comm,
  int            n_part,
  int           *pn_entity,
  int          **pentity_bound_part_idx,
  int          **pentity_bound,
  PDM_g_num_t ***pentity_ln_to_gn
)
{
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_g_part = 0;
  PDM_MPI_Allreduce(&n_part, &n_g_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  int *pn_entity_bound = NULL;
  PDM_malloc(pn_entity_bound, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_entity_bound[i_part] = pentity_bound_part_idx[i_part][n_g_part];
  }

  /* Prepare protocol */
  PDM_part_comm_graph_t *ptpgc = PDM_part_comm_graph_create(n_part,
                                                            pn_entity_bound,
                                                            pentity_bound,
                                                            comm);
  PDM_free(pn_entity_bound);


  PDM_g_num_t **_pentity_ln_to_gn = malloc(n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    _pentity_ln_to_gn[i_part] = malloc(pn_entity[i_part] * sizeof(PDM_g_num_t));
    for(int i = 0; i < pn_entity[i_part]; ++i) {
      _pentity_ln_to_gn[i_part][i] = -1;
    }
  }

  /* Count interface entity */
  int *n_l_entity_owner_bound = malloc(n_part * sizeof(int));
  int *n_l_entity_ghost_bound = malloc(n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    n_l_entity_owner_bound[i_part] = 0;
    n_l_entity_ghost_bound[i_part] = 0;
    int n_entity_bound = pentity_bound_part_idx[i_part][n_g_part];
    for(int i = 0; i < n_entity_bound; ++i) {
      int i_entity = pentity_bound[i_part][4*i  ]-1;
      int t_rank   = pentity_bound[i_part][4*i+1];
      int t_part   = pentity_bound[i_part][4*i+2]-1;
      if(_pentity_ln_to_gn[i_part][i_entity] == -1) {
        if(i_rank < t_rank || (i_rank == t_rank && i_part < t_part)) {
          n_l_entity_owner_bound[i_part]++;
          _pentity_ln_to_gn[i_part][i_entity] = -2;
        } else {
          n_l_entity_ghost_bound[i_part]++;
          _pentity_ln_to_gn[i_part][i_entity] = -3;
        }
      }
    }
  }

  int *n_l_entity_interior = malloc(n_part * sizeof(int));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    n_l_entity_interior[i_part] = pn_entity[i_part] - n_l_entity_owner_bound[i_part] - n_l_entity_ghost_bound[i_part];
  }

  /* Syncho rank */
  PDM_g_num_t l_shift = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    l_shift += n_l_entity_interior[i_part] + n_l_entity_owner_bound[i_part];
  }

  PDM_g_num_t g_shift = 0;
  PDM_MPI_Exscan(&l_shift, &g_shift, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, comm);

  if(0 == 1) {
    PDM_log_trace_array_int(pn_entity             , n_part, "pn_entity              ::");
    PDM_log_trace_array_int(n_l_entity_interior   , n_part, "n_l_entity_interior    ::");
    PDM_log_trace_array_int(n_l_entity_owner_bound, n_part, "n_l_entity_owner_bound ::");
    PDM_log_trace_array_int(n_l_entity_ghost_bound, n_part, "n_l_entity_ghost_bound ::");

    log_trace("l_shift = %i \n", l_shift);
    log_trace("g_shift = %i \n", g_shift);
  }

  /* Generation ln_to_gn */
  PDM_g_num_t l_part_shift = 0;
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int idx_write_interior = 0;
    int idx_write_bound    = 0;
    for(int i = 0; i < pn_entity[i_part]; ++i) {
      if(_pentity_ln_to_gn[i_part][i] == -1) { // Interior
        _pentity_ln_to_gn[i_part][i] = g_shift + idx_write_interior + 1;
        idx_write_interior++;
      } else if (_pentity_ln_to_gn[i_part][i] == -2) {
        _pentity_ln_to_gn[i_part][i] = g_shift + n_l_entity_interior[i_part] + idx_write_bound + 1;
        idx_write_bound++;
      }
    }
    l_part_shift += n_l_entity_owner_bound[i_part];
  }

  /* Prepare send */
  int         **send_gnum_n = malloc(n_part * sizeof(int         *));
  PDM_g_num_t **send_gnum   = malloc(n_part * sizeof(PDM_g_num_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_entity_bound = pentity_bound_part_idx[i_part][n_g_part];
    send_gnum  [i_part] = malloc(n_entity_bound * sizeof(PDM_g_num_t));
    send_gnum_n[i_part] = malloc(n_entity_bound * sizeof(int        ));

    int idx_write = 0;
    for(int i = 0; i < n_entity_bound; ++i) {
      int i_entity = pentity_bound[i_part][4*i  ]-1;
      int t_rank   = pentity_bound[i_part][4*i+1];
      int t_part   = pentity_bound[i_part][4*i+2]-1;
      if(_pentity_ln_to_gn[i_part][i_entity] > 0) {
        if(i_rank < t_rank || (i_rank == t_rank && i_part < t_part)) {
          send_gnum_n[i_part][i] = 1;
          send_gnum  [i_part][idx_write] = _pentity_ln_to_gn[i_part][i_entity];
          idx_write++;
        } else {
          send_gnum_n[i_part][i] = 0;
        }
      } else {
        send_gnum_n[i_part][i] = 0;
      }
    }

    // PDM_log_trace_array_int (send_gnum_n[i_part], n_entity_bound, "send_gnum_n ::" );
    // PDM_log_trace_array_long(send_gnum  [i_part], idx_write, "send_gnum ::" );
  }

  PDM_g_num_t **recv_gnum   = NULL;
  int         **recv_gnum_n = NULL;
  PDM_part_comm_graph_exch(ptpgc,
                                   sizeof(PDM_g_num_t),
                                   PDM_STRIDE_VAR_INTERLACED,
                                   1,
                                   send_gnum_n,
                      (void **)    send_gnum,
                                   &recv_gnum_n,
                      (void ***)   &recv_gnum);


  /*
   * Copy
   */
  for(int i_part = 0; i_part < n_part; ++i_part) {
    int n_entity_bound = pentity_bound_part_idx[i_part][n_g_part];

    if(0 == 1) {
      int n_recv = 0;
      for(int i = 0; i < n_entity_bound; ++i) {
        n_recv += recv_gnum_n[i_part][i];
      }
      PDM_log_trace_array_long(_pentity_ln_to_gn[i_part], pn_entity[i_part], "_pentity_ln_to_gn ::");
      PDM_log_trace_array_long(recv_gnum[i_part], n_recv, "recv_gnum ::");
      PDM_log_trace_array_int (recv_gnum_n[i_part], n_entity_bound, "recv_gnum_n ::");
    }

    int idx_read = 0;
    for(int i = 0; i < n_entity_bound; ++i) {
      int i_entity = pentity_bound[i_part][4*i  ]-1;
      if(recv_gnum_n[i_part][i] == 1) {
        // log_trace("i_entity = %i replace by recv_gnum[%i] = %i (before = %i) \n", i_entity, idx_read, recv_gnum[i_part][idx_read], _pentity_ln_to_gn[i_part][i_entity] );
        _pentity_ln_to_gn   [i_part][i_entity] = recv_gnum [i_part][idx_read];
        idx_read++;
      }
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(send_gnum_n[i_part]);
    free(send_gnum  [i_part]);
    free(recv_gnum_n[i_part]);
    free(recv_gnum  [i_part]);
  }
  free(send_gnum_n);
  free(recv_gnum_n);
  free(send_gnum);
  free(recv_gnum);
  free(n_l_entity_interior);
  free(n_l_entity_owner_bound);
  free(n_l_entity_ghost_bound);

  *pentity_ln_to_gn = _pentity_ln_to_gn;

  PDM_part_comm_graph_free(ptpgc);
}



/*=============================================================================
 * Test-specific function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int main(int argc, char *argv[])
{
  PDM_MPI_Init(&argc, &argv);

  PDM_g_num_t n_vtx_seg               = 30;
  double      length                  = 1.;
  double      zero_x                  = -length/2;
  double      zero_y                  = -length/2;
  double      zero_z                  = -length/2;
  int         n_part                  = 1;

  _read_args(argc,
             argv,
             &n_vtx_seg,
             &length);

  /*
   * Generate a cube
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int n_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int          *pn_vtx                 = NULL;
  int          *pn_edge                = NULL;
  int          *pn_face                = NULL;
  int          *pn_cell                = NULL;
  int          *pn_surface             = NULL;
  int          *pn_ridge               = NULL;
  double      **pvtx_coord             = NULL;
  int         **pedge_vtx              = NULL;
  int         **pface_edge_idx         = NULL;
  int         **pface_edge             = NULL;
  int         **pface_vtx              = NULL;
  int         **pcell_face_idx         = NULL;
  int         **pcell_face             = NULL;
  int         **psurface_face_idx      = NULL;
  int         **psurface_face          = NULL;
  int         **pridge_edge_idx        = NULL;
  int         **pridge_edge            = NULL;
  PDM_g_num_t **pvtx_ln_to_gn          = NULL;
  PDM_g_num_t **pedge_ln_to_gn         = NULL;
  PDM_g_num_t **pface_ln_to_gn         = NULL;
  PDM_g_num_t **pcell_ln_to_gn         = NULL;
  PDM_g_num_t **psurface_face_ln_to_gn = NULL;
  PDM_g_num_t **pridge_edge_ln_to_gn   = NULL;

  PDM_generate_mesh_parallelepiped_ngon(comm,
                                        PDM_MESH_NODAL_HEXA8,
                                        1,
                                        NULL,
                                        zero_x,
                                        zero_y,
                                        zero_z,
                                        length,
                                        length,
                                        length,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        2,
                                        n_part,
                                        PDM_SPLIT_DUAL_WITH_PTSCOTCH,
                                       &pn_vtx,
                                       &pn_edge,
                                       &pn_face,
                                       &pn_cell,
                                       &pvtx_coord,
                                       &pedge_vtx,
                                       &pface_edge_idx,
                                       &pface_edge,
                                       &pface_vtx,
                                       &pcell_face_idx,
                                       &pcell_face,
                                       &pvtx_ln_to_gn,
                                       &pedge_ln_to_gn,
                                       &pface_ln_to_gn,
                                       &pcell_ln_to_gn,
                                       &pn_surface,
                                       &psurface_face_idx,
                                       &psurface_face,
                                       &psurface_face_ln_to_gn,
                                       &pn_ridge,
                                       &pridge_edge_idx,
                                       &pridge_edge,
                                       &pridge_edge_ln_to_gn);


  // PDM_log_trace_array_int(n_owned_vtx, n_part, "n_owned_vtx ::");
  // free(n_owned_vtx);
  // for (int i_part = 0; i_part < n_part; i_part++) {
  //   free(lnum_owned_vtx[i_part]);
  // }
  // free(lnum_owned_vtx);


  int         **pvtx_bound_proc_idx  = NULL;
  int         **pvtx_bound_part_idx  = NULL;
  int         **pvtx_bound           = NULL;
  int         **pinternal_vtx_priority        = NULL;
  PDM_g_num_t *distrib_partition = PDM_compute_entity_distribution(comm, n_part);
  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      NULL,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pvtx_bound_proc_idx,
                                     &pvtx_bound_part_idx,
                                     &pvtx_bound,
                                     &pinternal_vtx_priority);


  int         **pedge_bound_proc_idx  = NULL;
  int         **pedge_bound_part_idx  = NULL;
  int         **pedge_bound           = NULL;
  int         **pinternal_edge_priority        = NULL;
  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      NULL,
                                      n_part,
                                      pn_edge,
               (const PDM_g_num_t **) pedge_ln_to_gn,
                                      NULL,
                                     &pedge_bound_proc_idx,
                                     &pedge_bound_part_idx,
                                     &pedge_bound,
                                     &pinternal_edge_priority);


  int         **pface_bound_proc_idx  = NULL;
  int         **pface_bound_part_idx  = NULL;
  int         **pface_bound           = NULL;
  int         **pinternal_face_priority        = NULL;
  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      NULL,
                                      n_part,
                                      pn_face,
               (const PDM_g_num_t **) pface_ln_to_gn,
                                      NULL,
                                     &pface_bound_proc_idx,
                                     &pface_bound_part_idx,
                                     &pface_bound,
                                     &pinternal_face_priority);

  free(distrib_partition);


  /* Start test with vtx */
  PDM_MPI_Barrier(comm);
  double t1 = PDM_MPI_Wtime();

  PDM_g_num_t **pold_ln_to_gn = NULL;
  _generate_gnum_old(comm,
                     n_part,
                     pn_vtx,
                     pvtx_ln_to_gn,
                     &pold_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pold_ln_to_gn[i_part]);
  }
  free(pold_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_old = PDM_MPI_Wtime() - t1;


  PDM_MPI_Barrier(comm);
  double t2 = PDM_MPI_Wtime();

  int n_g_part = 0;
  PDM_MPI_Allreduce(&n_part, &n_g_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

  PDM_g_num_t **pnew_vtx_ln_to_gn = NULL;
  _generate_gnum(comm,
                 n_part,
                 pn_vtx,
                 pvtx_bound_part_idx,
                 pvtx_bound,
                 &pnew_vtx_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_new = PDM_MPI_Wtime() - t2;

  if(i_rank == 0) {
    printf("vtx : dt_old = %12.5e / dt_new = %12.5e -> %12.5e\n", dt_old, dt_new, dt_old/dt_new);
  }

  if(1 == 1) {
    char filename[999];
    for(int i_part = 0; i_part < n_part; ++i_part) {

      int *pvtx_owner = malloc(pn_vtx[i_part] * sizeof(int));
      for(int i = 0; i < pn_vtx[i_part]; ++i) {
        pvtx_owner[i] = -1;
      }

      int n_entity_bound = pvtx_bound_part_idx[i_part][n_g_part];
      for(int i = 0; i < n_entity_bound; ++i) {
        int i_entity = pvtx_bound[i_part][4*i]-1;
        int t_rank   = pvtx_bound[i_part][4*i+1];
        int t_part   = pvtx_bound[i_part][4*i+2]-1;
        if(pvtx_owner[i_entity] != 0) {
          if(i_rank < t_rank || (i_rank == t_rank && i_part < t_part)) {
            pvtx_owner[i_entity] = 1;
          } else {
            pvtx_owner[i_entity] = 0;
          }
        }
      }


      sprintf(filename, "point_cloud_before_%04d_%04d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_vtx[i_part],
                                pvtx_coord[i_part],
                                NULL,
                                pvtx_owner);
      sprintf(filename, "point_cloud_after_%04d_%04d.vtk", i_part, i_rank);
      PDM_vtk_write_point_cloud(filename,
                                pn_vtx[i_part],
                                pvtx_coord[i_part],
                                pnew_vtx_ln_to_gn[i_part],
                                pvtx_owner);

      free(pvtx_owner);
    }
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    // PDM_log_trace_array_long(pnew_vtx_ln_to_gn[i_part], pn_vtx[i_part], "pnew_vtx_ln_to_gn : ");
    free(pnew_vtx_ln_to_gn[i_part]);
  }
  free(pnew_vtx_ln_to_gn);

  /* Edges */
  PDM_MPI_Barrier(comm);
  double t1_edge = PDM_MPI_Wtime();

  PDM_g_num_t **pold_edge_ln_to_gn = NULL;
  _generate_gnum_old(comm,
                     n_part,
                     pn_edge,
                     pedge_ln_to_gn,
                     &pold_edge_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pold_edge_ln_to_gn[i_part]);
  }
  free(pold_edge_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_old_edge = PDM_MPI_Wtime() - t1_edge;


  PDM_MPI_Barrier(comm);
  double t2_edge = PDM_MPI_Wtime();

  PDM_g_num_t **pnew_edge_ln_to_gn = NULL;
  _generate_gnum(comm,
                 n_part,
                 pn_edge,
                 pedge_bound_part_idx,
                 pedge_bound,
                 &pnew_edge_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_new_edge = PDM_MPI_Wtime() - t2_edge;

  if(i_rank == 0) {
    printf("edge : dt_old = %12.5e / dt_new = %12.5e -> %12.5e\n", dt_old_edge, dt_new_edge, dt_old_edge/dt_new_edge);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pnew_edge_ln_to_gn[i_part]);
  }
  free(pnew_edge_ln_to_gn);


  /* Face */
  PDM_MPI_Barrier(comm);
  double t1_face = PDM_MPI_Wtime();

  PDM_g_num_t **pold_face_ln_to_gn = NULL;
  _generate_gnum_old(comm,
                     n_part,
                     pn_face,
                     pface_ln_to_gn,
                     &pold_face_ln_to_gn);

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pold_face_ln_to_gn[i_part]);
  }
  free(pold_face_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_old_face = PDM_MPI_Wtime() - t1_face;


  PDM_MPI_Barrier(comm);
  double t2_face = PDM_MPI_Wtime();

  PDM_g_num_t **pnew_face_ln_to_gn = NULL;
  _generate_gnum(comm,
                 n_part,
                 pn_face,
                 pface_bound_part_idx,
                 pface_bound,
                 &pnew_face_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_new_face = PDM_MPI_Wtime() - t2_face;

  if(i_rank == 0) {
    printf("face : dt_old = %12.5e / dt_new = %12.5e -> %12.5e\n", dt_old_face, dt_new_face, dt_old_face/dt_new_face);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    free(pnew_face_ln_to_gn[i_part]);
  }
  free(pnew_face_ln_to_gn);



  if(1 == 0) {
    /*
     *
     */
    int *pn_vtx_bound = NULL;
    PDM_malloc(pn_vtx_bound, n_part, int);
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_vtx_bound[i_part] = pvtx_bound_part_idx[i_part][n_g_part];
    }
    PDM_part_comm_graph_t *ptpgc = PDM_part_comm_graph_create(n_part,
                                                                              pn_vtx_bound,
                                                                              pvtx_bound,
                                                                              comm);
    PDM_free(pn_vtx_bound);


    /* Prepare send */
    int         **send_gnum_n   = malloc(n_part * sizeof(int         *));
    PDM_g_num_t **send_gnum_var = malloc(n_part * sizeof(PDM_g_num_t *));
    PDM_g_num_t **send_gnum     = malloc(n_part * sizeof(PDM_g_num_t *));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      int n_entity_bound = pvtx_bound_part_idx[i_part][n_g_part];
      send_gnum    [i_part] = malloc(n_entity_bound * sizeof(PDM_g_num_t));
      send_gnum_var[i_part] = malloc(n_entity_bound * sizeof(PDM_g_num_t));
      send_gnum_n  [i_part] = malloc(n_entity_bound * sizeof(int        ));

      int idx_write = 0;
      for(int i = 0; i < n_entity_bound; ++i) {
        int i_entity = pvtx_bound[i_part][4*i]-1;
        send_gnum[i_part][i] = pvtx_ln_to_gn[i_part][i_entity];

        int t_rank = pvtx_bound[i_part][4*i+1];
        if(i_rank < t_rank) {
          send_gnum_n  [i_part][i] = 1;
          send_gnum_var[i_part][idx_write++] = pvtx_ln_to_gn[i_part][i_entity];
        } else {
          send_gnum_n  [i_part][i] = 0;
        }

      }
    }

    PDM_g_num_t **recv_gnum = NULL;
    PDM_part_comm_graph_exch(ptpgc,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_CST_INTERLACED,
                                     1,
                                     NULL,
                         (void **)   send_gnum,
                                     NULL,
                         (void ***)  &recv_gnum);

    if(0 == 1) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        int n_entity_bound = pvtx_bound_part_idx[i_part][n_g_part];
        PDM_log_trace_array_long(recv_gnum[i_part], n_entity_bound, "recv_gnum ::");
      }
    }


    PDM_g_num_t **recv_gnum_var = NULL;
    int         **recv_gnum_n   = NULL;
    PDM_part_comm_graph_exch(ptpgc,
                                     sizeof(PDM_g_num_t),
                                     PDM_STRIDE_VAR_INTERLACED,
                                     1,
                                     send_gnum_n,
                        (void **)    send_gnum,
                                     &recv_gnum_n,
                        (void ***)   &recv_gnum_var);

    if(0 == 1) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        int n_entity_bound = pvtx_bound_part_idx[i_part][n_g_part];

        int n_recv = 0;
        for(int i = 0; i < n_entity_bound; ++i) {
          n_recv += recv_gnum_n[i_part][i];
        }
        PDM_log_trace_array_int (recv_gnum_n  [i_part], n_entity_bound, "recv_gnum_n ::");
        PDM_log_trace_array_long(recv_gnum_var[i_part], n_recv, "recv_gnum_var ::");
      }
    }

    for(int i_part = 0; i_part < n_part; ++i_part) {
      free(send_gnum_n  [i_part]);
      free(send_gnum_var[i_part]);
      free(recv_gnum_n  [i_part]);
      free(recv_gnum_var[i_part]);
      free(send_gnum    [i_part]);
      free(recv_gnum    [i_part]);
    }
    free(send_gnum_n);
    free(send_gnum_var);
    free(recv_gnum_n);
    free(recv_gnum_var);
    free(send_gnum);
    free(recv_gnum);

    PDM_part_comm_graph_free(ptpgc);
  }



  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pvtx_bound_proc_idx[i_part]);
    free(pvtx_bound_part_idx[i_part]);
    free(pvtx_bound         [i_part]);
    free(pinternal_vtx_priority      [i_part]);
  }
  free(pvtx_bound_proc_idx);
  free(pvtx_bound_part_idx);
  free(pvtx_bound         );
  free(pinternal_vtx_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pedge_bound_proc_idx[i_part]);
    free(pedge_bound_part_idx[i_part]);
    free(pedge_bound         [i_part]);
    free(pinternal_edge_priority      [i_part]);
  }
  free(pedge_bound_proc_idx);
  free(pedge_bound_part_idx);
  free(pedge_bound         );
  free(pinternal_edge_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pface_bound_proc_idx[i_part]);
    free(pface_bound_part_idx[i_part]);
    free(pface_bound         [i_part]);
    free(pinternal_face_priority      [i_part]);
  }
  free(pface_bound_proc_idx);
  free(pface_bound_part_idx);
  free(pface_bound         );
  free(pinternal_face_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    free(pvtx_coord            [i_part]);
    free(pedge_vtx             [i_part]);
    free(pface_edge_idx        [i_part]);
    free(pface_edge            [i_part]);
    free(pface_vtx             [i_part]);
    free(pcell_face_idx        [i_part]);
    free(pcell_face            [i_part]);
    free(psurface_face_idx     [i_part]);
    free(psurface_face         [i_part]);
    free(pridge_edge_idx       [i_part]);
    free(pridge_edge           [i_part]);
    free(pvtx_ln_to_gn         [i_part]);
    free(pedge_ln_to_gn        [i_part]);
    free(pface_ln_to_gn        [i_part]);
    free(pcell_ln_to_gn        [i_part]);
    free(psurface_face_ln_to_gn[i_part]);
    free(pridge_edge_ln_to_gn  [i_part]);
  }
  free(pn_vtx                );
  free(pn_edge               );
  free(pn_face               );
  free(pn_cell               );
  free(pn_surface            );
  free(pn_ridge              );
  free(pvtx_coord            );
  free(pedge_vtx             );
  free(pface_edge_idx        );
  free(pface_edge            );
  free(pface_vtx             );
  free(pcell_face_idx        );
  free(pcell_face            );
  free(psurface_face_idx     );
  free(psurface_face         );
  free(pridge_edge_idx       );
  free(pridge_edge           );
  free(pvtx_ln_to_gn         );
  free(pedge_ln_to_gn        );
  free(pface_ln_to_gn        );
  free(pcell_ln_to_gn        );
  free(psurface_face_ln_to_gn);
  free(pridge_edge_ln_to_gn  );

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
