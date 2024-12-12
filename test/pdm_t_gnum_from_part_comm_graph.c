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
  PDM_MPI_Comm             comm,
  int                      n_part,
  int                     *pn_entity,
  PDM_part_comm_graph_t   *pcg,
  PDM_g_num_t           ***out_pentity_ln_to_gn
)
{
  PDM_gen_gnum_t *gnum = PDM_gnum_create(3,
                                         n_part,
                                         PDM_FALSE,
                                         1.,
                                         comm,
                                         PDM_OWNERSHIP_USER);

  PDM_gnum_set_from_part_comm_graph(gnum,
                                    pn_entity,
                                    pcg);
  PDM_gnum_compute(gnum);

  *out_pentity_ln_to_gn = malloc(sizeof(PDM_g_num_t *));
  for (int i_part = 0; i_part < n_part; i_part++) {
    (*out_pentity_ln_to_gn)[i_part] = PDM_gnum_get(gnum, i_part);
  }

  PDM_gnum_free(gnum);
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
  int i_rank;
  PDM_MPI_Comm_size(comm, &n_rank);
  PDM_MPI_Comm_rank(comm, &i_rank);

  int n_g_part = 0;
  PDM_MPI_Allreduce(&n_part, &n_g_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);

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
                                        n_vtx_seg,
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
  /*
   * Vertices
   */
  int         **pvtx_bound_proc_idx    = NULL;
  int         **pvtx_bound_part_idx    = NULL;
  int         **pvtx_bound             = NULL;
  int         **pinternal_vtx_priority = NULL;
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

  int *pn_vtx_bound = NULL;
  PDM_malloc(pn_vtx_bound, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_vtx_bound[i_part] = pvtx_bound_part_idx[i_part][n_g_part];
  }
  PDM_part_comm_graph_t *pgc_vtx = PDM_part_comm_graph_create(n_part,
                                                              pn_vtx_bound,
                                                              pvtx_bound,
                                                              comm);
  PDM_free(pn_vtx_bound);

  /*
   * Edges
   */
  int         **pedge_bound_proc_idx    = NULL;
  int         **pedge_bound_part_idx    = NULL;
  int         **pedge_bound             = NULL;
  int         **pinternal_edge_priority = NULL;
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
  int *pn_edge_bound = NULL;
  PDM_malloc(pn_edge_bound, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_edge_bound[i_part] = pedge_bound_part_idx[i_part][n_g_part];
  }
  PDM_part_comm_graph_t *pgc_edge = PDM_part_comm_graph_create(n_part,
                                                               pn_edge_bound,
                                                               pedge_bound,
                                                               comm);
  PDM_free(pn_edge_bound);

  /*
   * Faces
   */
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


  int *pn_face_bound = NULL;
  PDM_malloc(pn_face_bound, n_part, int);
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_face_bound[i_part] = pface_bound_part_idx[i_part][n_g_part];
  }
  PDM_part_comm_graph_t *pgc_face = PDM_part_comm_graph_create(n_part,
                                                               pn_face_bound,
                                                               pface_bound,
                                                               comm);
  PDM_free(pn_face_bound);

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
    PDM_free(pold_ln_to_gn[i_part]);
  }
  PDM_free(pold_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_old = PDM_MPI_Wtime() - t1;


  PDM_MPI_Barrier(comm);
  double t2 = PDM_MPI_Wtime();

  PDM_g_num_t **pnew_vtx_ln_to_gn = NULL;
  _generate_gnum(comm,
                 n_part,
                 pn_vtx,
                 pgc_vtx,
                 &pnew_vtx_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_new = PDM_MPI_Wtime() - t2;

  if(i_rank == 0) {
    printf("vtx : dt_old = %12.5e / dt_new = %12.5e -> %12.5e\n", dt_old, dt_new, dt_old/dt_new);
  }

  if(0 == 1) {
    char filename[999];
    for(int i_part = 0; i_part < n_part; ++i_part) {

      const int *pvtx_graph_owner = PDM_part_comm_graph_owner_get(pgc_vtx, i_part);
      int *pvtx_owner = NULL;
      PDM_malloc(pvtx_owner, pn_vtx[i_part], int);
      for(int i = 0; i < pn_vtx[i_part]; ++i) {
        pvtx_owner[i] = -1;
      }

      int n_entity_bound = pvtx_bound_part_idx[i_part][n_g_part];
      for(int i = 0; i < n_entity_bound; ++i) {
        int i_entity = pvtx_bound[i_part][4*i]-1;
        if(pvtx_owner[i_entity] != 0) {
          pvtx_owner[i_entity] = pvtx_graph_owner[i];
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

      PDM_free(pvtx_owner);
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
                 pgc_edge,
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
                 pgc_face,
                 &pnew_face_ln_to_gn);

  PDM_MPI_Barrier(comm);
  double dt_new_face = PDM_MPI_Wtime() - t2_face;

  if(i_rank == 0) {
    printf("face : dt_old = %12.5e / dt_new = %12.5e -> %12.5e\n", dt_old_face, dt_new_face, dt_old_face/dt_new_face);
  }

  for(int i_part = 0; i_part < n_part; ++i_part) {
    PDM_free(pnew_face_ln_to_gn[i_part]);
  }
  PDM_free(pnew_face_ln_to_gn);


  PDM_part_comm_graph_free(pgc_vtx);
  PDM_part_comm_graph_free(pgc_edge);
  PDM_part_comm_graph_free(pgc_face);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pvtx_bound_proc_idx   [i_part]);
    PDM_free(pvtx_bound_part_idx   [i_part]);
    PDM_free(pvtx_bound            [i_part]);
    PDM_free(pinternal_vtx_priority[i_part]);
  }
  PDM_free(pvtx_bound_proc_idx   );
  PDM_free(pvtx_bound_part_idx   );
  PDM_free(pvtx_bound            );
  PDM_free(pinternal_vtx_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pedge_bound_proc_idx   [i_part]);
    PDM_free(pedge_bound_part_idx   [i_part]);
    PDM_free(pedge_bound            [i_part]);
    PDM_free(pinternal_edge_priority[i_part]);
  }
  PDM_free(pedge_bound_proc_idx   );
  PDM_free(pedge_bound_part_idx   );
  PDM_free(pedge_bound            );
  PDM_free(pinternal_edge_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pface_bound_proc_idx   [i_part]);
    PDM_free(pface_bound_part_idx   [i_part]);
    PDM_free(pface_bound            [i_part]);
    PDM_free(pinternal_face_priority[i_part]);
  }
  PDM_free(pface_bound_proc_idx);
  PDM_free(pface_bound_part_idx);
  PDM_free(pface_bound         );
  PDM_free(pinternal_face_priority);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_free(pvtx_coord            [i_part]);
    PDM_free(pedge_vtx             [i_part]);
    PDM_free(pface_edge_idx        [i_part]);
    PDM_free(pface_edge            [i_part]);
    PDM_free(pface_vtx             [i_part]);
    PDM_free(pcell_face_idx        [i_part]);
    PDM_free(pcell_face            [i_part]);
    PDM_free(psurface_face_idx     [i_part]);
    PDM_free(psurface_face         [i_part]);
    PDM_free(pridge_edge_idx       [i_part]);
    PDM_free(pridge_edge           [i_part]);
    PDM_free(pvtx_ln_to_gn         [i_part]);
    PDM_free(pedge_ln_to_gn        [i_part]);
    PDM_free(pface_ln_to_gn        [i_part]);
    PDM_free(pcell_ln_to_gn        [i_part]);
    PDM_free(psurface_face_ln_to_gn[i_part]);
    PDM_free(pridge_edge_ln_to_gn  [i_part]);
  }
  PDM_free(pn_vtx                );
  PDM_free(pn_edge               );
  PDM_free(pn_face               );
  PDM_free(pn_cell               );
  PDM_free(pn_surface            );
  PDM_free(pn_ridge              );
  PDM_free(pvtx_coord            );
  PDM_free(pedge_vtx             );
  PDM_free(pface_edge_idx        );
  PDM_free(pface_edge            );
  PDM_free(pface_vtx             );
  PDM_free(pcell_face_idx        );
  PDM_free(pcell_face            );
  PDM_free(psurface_face_idx     );
  PDM_free(psurface_face         );
  PDM_free(pridge_edge_idx       );
  PDM_free(pridge_edge           );
  PDM_free(pvtx_ln_to_gn         );
  PDM_free(pedge_ln_to_gn        );
  PDM_free(pface_ln_to_gn        );
  PDM_free(pcell_ln_to_gn        );
  PDM_free(psurface_face_ln_to_gn);
  PDM_free(pridge_edge_ln_to_gn  );

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }
  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}
