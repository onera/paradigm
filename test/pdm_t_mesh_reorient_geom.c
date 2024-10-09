#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <pdm_mpi.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_generate_mesh.h"
#include "pdm_array.h"
#include "pdm_part_mesh_reorient_geom.h"

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

  PDM_g_num_t n_vtx_seg = 10;
  int         n_part    = 2;

  // Initialize MPI
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);

  /* 1) 2D */

  // Generate partitioned mesh

  double *inv = malloc(3 * sizeof(double));
  inv[0] = 0.;
  inv[1] = 0.;
  inv[2] = 1.;

  int          *n_vtx         = NULL;
  int          *n_edge        = NULL;
  int          *n_face        = NULL;
  double      **vtx_coord     = NULL;
  int         **edge_vtx      = NULL;
  int         **face_edge_idx = NULL;
  int         **face_edge     = NULL;
  int         **face_vtx      = NULL;
  PDM_g_num_t **vtx_ln_to_gn  = NULL;
  PDM_g_num_t **edge_ln_to_gn = NULL;
  PDM_g_num_t **face_ln_to_gn = NULL;
  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_POLY_2D,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   n_vtx_seg,
                                   n_vtx_seg,
                                   n_part,
                                   PDM_SPLIT_DUAL_WITH_HILBERT,
                                   0.,
                                   &n_vtx,
                                   &n_edge,
                                   &n_face,
                                   &vtx_coord,
                                   &edge_vtx,
                                   &face_edge_idx,
                                   &face_edge,
                                   &face_vtx,
                                   &vtx_ln_to_gn,
                                   &edge_ln_to_gn,
                                   &face_ln_to_gn);

  int reorient = PDM_part_mesh_reorient_geom(2,
                                             inv,
                                             n_part,
                                             n_face,
                                             n_edge,
                                             0,
                                             face_edge_idx,
                                             face_edge,
                                             NULL,
                                             NULL,
                                             edge_vtx,
                                             vtx_coord,
                                             NULL,
                                             NULL,
                                             comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 2D 1 :D\n");
    fflush(stdout);
  }

  int **edge_face = (int**) malloc(n_part * sizeof(int*));
  for (int i_part = 0; i_part < n_part; i_part++) {
    int* _edge_face = PDM_array_const_int(2 * n_edge[i_part], 0);
    for (int i_face = 0; i_face < n_face[i_part]; i_face++) {
      for (int idx_edge = face_edge_idx[i_part][i_face]; idx_edge < face_edge_idx[i_part][i_face+1]; idx_edge++) {
        int i_edge = PDM_ABS(face_edge[i_part][idx_edge])-1;
        if (face_edge[i_part][idx_edge] > 0) {
          _edge_face[2*i_edge+1] = i_face+1;
        }
        if (face_edge[i_part][idx_edge] < 0) {
          _edge_face[2*i_edge  ] = i_face+1;
        }
      }
    }
    edge_face[i_part] = _edge_face;
  }

  reorient = PDM_part_mesh_reorient_geom(2,
                                         inv,
                                         n_part,
                                         n_face,
                                         n_edge,
                                         0,
                                         face_edge_idx,
                                         face_edge,
                                         edge_face,
                                         NULL,
                                         edge_vtx,
                                         vtx_coord,
                                         NULL,
                                         NULL,
                                         comm);
  assert(reorient == 1);

  if (i_rank == 0) {
    printf("Ok 2D 2 :D\n");
    fflush(stdout);
  }

  reorient = PDM_part_mesh_reorient_geom(2,
                                         inv,
                                         n_part,
                                         n_face,
                                         n_edge,
                                         0,
                                         NULL,
                                         NULL,
                                         edge_face,
                                         NULL,
                                         edge_vtx,
                                         vtx_coord,
                                         NULL,
                                         NULL,
                                         comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 2D 3 :D\n");
    fflush(stdout);
  }

  reorient = PDM_part_mesh_reorient_geom(2,
                                         inv,
                                         n_part,
                                         n_face,
                                         n_edge,
                                         0,
                                         face_edge_idx,
                                         face_edge,
                                         NULL,
                                         NULL,
                                         edge_vtx,
                                         vtx_coord,
                                         NULL,
                                         NULL,
                                         comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 2D 4 :D\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(vtx_coord    [i_part]);
    free(edge_vtx     [i_part]);
    free(face_edge_idx[i_part]);
    free(face_edge    [i_part]);
    free(face_vtx     [i_part]);
    free(edge_face    [i_part]);
    free(vtx_ln_to_gn [i_part]);
    free(edge_ln_to_gn[i_part]);
    free(face_ln_to_gn[i_part]);
  }
  free(n_vtx        );
  free(n_edge       );
  free(n_face       );
  free(vtx_coord    );
  free(edge_vtx     );
  free(face_edge_idx);
  free(face_edge    );
  free(edge_face    );
  free(face_vtx     );
  free(vtx_ln_to_gn );
  free(edge_ln_to_gn);
  free(face_ln_to_gn);
  free(inv);

  /* 2) 3D */

  // Generate partitioned mesh

  int          *n_cell                = NULL;
  int          *n_surface             = NULL;
  int          *n_ridge               = NULL;
  int         **cell_face_idx         = NULL;
  int         **cell_face             = NULL;
  int         **surface_face_idx      = NULL;
  int         **surface_face          = NULL;
  int         **ridge_edge_idx        = NULL;
  int         **ridge_edge            = NULL;
  PDM_g_num_t **cell_ln_to_gn         = NULL;
  PDM_g_num_t **surface_face_ln_to_gn = NULL;
  PDM_g_num_t **ridge_edge_ln_to_gn   = NULL;
  PDM_generate_mesh_parallelepiped_ngon(comm,
                                        PDM_MESH_NODAL_TETRA4,
                                        1,
                                        NULL,
                                        0.,
                                        0.,
                                        0.,
                                        10.,
                                        10.,
                                        10.,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_vtx_seg,
                                        n_part,
                                        PDM_SPLIT_DUAL_WITH_HILBERT,
                                        &n_vtx,
                                        &n_edge,
                                        &n_face,
                                        &n_cell,
                                        &vtx_coord,
                                        &edge_vtx,
                                        &face_edge_idx,
                                        &face_edge,
                                        &face_vtx,
                                        &cell_face_idx,
                                        &cell_face,
                                        &vtx_ln_to_gn,
                                        &edge_ln_to_gn,
                                        &face_ln_to_gn,
                                        &cell_ln_to_gn,
                                        &n_surface,
                                        &surface_face_idx,
                                        &surface_face,
                                        &surface_face_ln_to_gn,
                                        &n_ridge,
                                        &ridge_edge_idx,
                                        &ridge_edge,
                                        &ridge_edge_ln_to_gn);

  reorient = PDM_part_mesh_reorient_geom(3,
                                         NULL,
                                         n_part,
                                         n_cell,
                                         n_face,
                                         n_surface,
                                         cell_face_idx,
                                         cell_face,
                                         NULL,
                                         face_edge_idx,
                                         face_vtx,
                                         vtx_coord,
                           (const int**) surface_face_idx,
                           (const int**) surface_face,
                                         comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 3D 1 :D\n");
    fflush(stdout);
  }

  int **face_cell = (int**) malloc(n_part * sizeof(int*));
  for (int i_part = 0; i_part < n_part; i_part++) {
    int* _face_cell = PDM_array_const_int(2 * n_face[i_part], 0);
    for (int i_cell = 0; i_cell < n_cell[i_part]; i_cell++) {
      for (int idx_face = cell_face_idx[i_part][i_cell]; idx_face < cell_face_idx[i_part][i_cell+1]; idx_face++) {
        int i_face = PDM_ABS(cell_face[i_part][idx_face])-1;
        if (cell_face[i_part][idx_face] > 0) {
          _face_cell[2*i_face+1] = i_cell+1;
        }
        if (cell_face[i_part][idx_face] < 0) {
          _face_cell[2*i_face  ] = i_cell+1;
        }
      }
    }
    face_cell[i_part] = _face_cell;
  }

  reorient = PDM_part_mesh_reorient_geom(3,
                                         NULL,
                                         n_part,
                                         n_cell,
                                         n_face,
                                         n_surface,
                                         cell_face_idx,
                                         cell_face,
                                         face_cell,
                                         face_edge_idx,
                                         face_vtx,
                                         vtx_coord,
                           (const int**) surface_face_idx,
                           (const int**) surface_face,
                                         comm);
  assert(reorient == 1);

  if (i_rank == 0) {
    printf("Ok 3D 2 :D\n");
    fflush(stdout);
  }

  reorient = PDM_part_mesh_reorient_geom(3,
                                         NULL,
                                         n_part,
                                         n_cell,
                                         n_face,
                                         n_surface,
                                         NULL,
                                         NULL,
                                         face_cell,
                                         face_edge_idx,
                                         face_vtx,
                                         vtx_coord,
                           (const int**) surface_face_idx,
                           (const int**) surface_face,
                                         comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 3D 3 :D\n");
    fflush(stdout);
  }

  reorient = PDM_part_mesh_reorient_geom(3,
                                         NULL,
                                         n_part,
                                         n_cell,
                                         n_face,
                                         n_surface,
                                         cell_face_idx,
                                         cell_face,
                                         NULL,
                                         face_edge_idx,
                                         face_vtx,
                                         vtx_coord,
                           (const int**) surface_face_idx,
                           (const int**) surface_face,
                                         comm);
  assert(reorient == 0);

  if (i_rank == 0) {
    printf("Ok 3D 4 :D\n");
    fflush(stdout);
  }

  for (int i_part = 0; i_part < n_part; i_part++) {
    free(vtx_coord            [i_part]);
    free(edge_vtx             [i_part]);
    free(face_edge_idx        [i_part]);
    free(face_edge            [i_part]);
    free(face_vtx             [i_part]);
    free(cell_face_idx        [i_part]);
    free(cell_face            [i_part]);
    free(surface_face_idx     [i_part]);
    free(surface_face         [i_part]);
    free(ridge_edge_idx       [i_part]);
    free(ridge_edge           [i_part]);
    free(vtx_ln_to_gn         [i_part]);
    free(edge_ln_to_gn        [i_part]);
    free(face_ln_to_gn        [i_part]);
    free(cell_ln_to_gn        [i_part]);
    free(surface_face_ln_to_gn[i_part]);
    free(ridge_edge_ln_to_gn  [i_part]);
    free(face_cell            [i_part]);
  }
  free(n_vtx                );
  free(n_edge               );
  free(n_face               );
  free(n_cell               );
  free(n_surface            );
  free(n_ridge              );
  free(vtx_coord            );
  free(edge_vtx             );
  free(face_edge_idx        );
  free(face_edge            );
  free(face_vtx             );
  free(cell_face_idx        );
  free(cell_face            );
  free(surface_face_idx     );
  free(surface_face         );
  free(ridge_edge_idx       );
  free(ridge_edge           );
  free(vtx_ln_to_gn         );
  free(edge_ln_to_gn        );
  free(face_ln_to_gn        );
  free(cell_ln_to_gn        );
  free(surface_face_ln_to_gn);
  free(ridge_edge_ln_to_gn  );
  free(face_cell            );

  // Finalize
  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("The End :D\n");
    fflush(stdout);
  }

  return EXIT_SUCCESS;
}
