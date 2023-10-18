---
jupytext:
  text_representation:
    extension: '.md'
    format_name: myst
    format_version: '0.7'
    jupytext_version: 1.4.0+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Exercise 2 : Localization of a point cloud inside a mesh

+++

*(Load custom magics)*

```{code-cell} ipython3
import os, sys
module_path = os.path.abspath(os.path.join('../../utils'))
if module_path not in sys.path:
    sys.path.append(module_path)
```

```{code-cell}
%reload_ext visu_magics
%reload_ext code_magics
```

The goal of this exercise is to get used
Your job is to fill the code blocks left blank using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/prepro_algo/index.html#python-api).

+++

## Include the required headers

```{code-cell}
%%code_block -p exercise2 -i 1
// Required headers
#include "pdm_generate_mesh.h"
#include "pdm_mesh_location.h"
#include "pdm_writer_priv.h"

```

## Generate a partitioned "source" mesh

```{code-cell}
%%code_block -p exercise2 -i 2
int main(int argc, char *argv[])
{
  PDM_g_num_t src_n_vtx_seg = 10;
  PDM_g_num_t tgt_n_vtx_seg = 10;
  int         src_n_part    = 1;
  int         tgt_n_part    = 1;
  double      tgt_xmin      = 0.3;
  double      tgt_ymin      = 0.3;

  /*
   * Initialize MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank, n_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Generate and partition the source mesh
  int          *src_n_vtx         = NULL;
  int          *src_n_edge        = NULL;
  int          *src_n_face        = NULL;
  double      **src_vtx_coord     = NULL;
  int         **src_edge_vtx      = NULL;
  int         **src_face_edge_idx = NULL;
  int         **src_face_edge     = NULL;
  int         **src_face_vtx      = NULL;
  PDM_g_num_t **src_vtx_ln_to_gn  = NULL;
  PDM_g_num_t **src_edge_ln_to_gn = NULL;
  PDM_g_num_t **src_face_ln_to_gn = NULL;
  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_POLY_2D,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   src_n_vtx_seg,
                                   src_n_vtx_seg,
                                   src_n_part,
                                   PDM_SPLIT_DUAL_WITH_PARMETIS,
                                   &src_n_vtx,
                                   &src_n_edge,
                                   &src_n_face,
                                   &src_vtx_coord,
                                   &src_edge_vtx,
                                   &src_face_edge_idx,
                                   &src_face_edge,
                                   &src_face_vtx,
                                   &src_vtx_ln_to_gn,
                                   &src_edge_ln_to_gn,
                                   &src_face_ln_to_gn);

```

## Generate a partitioned "target" mesh
We will use its vertices as a point cloud.

```{code-cell}
%%code_block -p exercise2 -i 3
  int          *tgt_n_vtx         = NULL;
  int          *tgt_n_edge        = NULL;
  int          *tgt_n_face        = NULL;
  double      **tgt_vtx_coord     = NULL;
  int         **tgt_edge_vtx      = NULL;
  int         **tgt_face_edge_idx = NULL;
  int         **tgt_face_edge     = NULL;
  int         **tgt_face_vtx      = NULL;
  PDM_g_num_t **tgt_vtx_ln_to_gn  = NULL;
  PDM_g_num_t **tgt_edge_ln_to_gn = NULL;
  PDM_g_num_t **tgt_face_ln_to_gn = NULL;

  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_QUAD4,
                                   tgt_xmin,
                                   tgt_ymin,
                                   0.,
                                   1.,
                                   1.,
                                   tgt_n_vtx_seg,
                                   tgt_n_vtx_seg,
                                   tgt_n_part,
                                   PDM_SPLIT_DUAL_WITH_HILBERT,
                                   &tgt_n_vtx,
                                   &tgt_n_edge,
                                   &tgt_n_face,
                                   &tgt_vtx_coord,
                                   &tgt_edge_vtx,
                                   &tgt_face_edge_idx,
                                   &tgt_face_edge,
                                   &tgt_face_vtx,
                                   &tgt_vtx_ln_to_gn,
                                   &tgt_edge_ln_to_gn,
                                   &tgt_face_ln_to_gn);

```

## Create the `MeshLocation` object

```{code-cell}
%%code_block -p exercise2 -i 4
  // Create the PDM_mesh_location_t object
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

```

## Set the target point cloud

```{code-cell}
%%code_block -p exercise2 -i 5
  // Set target point cloud
  PDM_mesh_location_n_part_cloud_set(mesh_loc,
                                     0,
                                     tgt_n_part);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    PDM_mesh_location_cloud_set(mesh_loc,
                                0,
                                i_part,
                                tgt_n_vtx       [i_part],
                                tgt_vtx_coord   [i_part],
                                tgt_vtx_ln_to_gn[i_part]);
  }
```

## Set the source mesh
Here you have essentially two options :
- you can either define the mesh with "nodal" connectivity (i.e. Finite-Element style)
- or with "descending" connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -p exercise2 -i 6
  // Set source mesh
  PDM_mesh_location_mesh_n_part_set(mesh_loc, src_n_part);

  int nodal = 1;

  if (nodal) {
    for (int i_part = 0; i_part < src_n_part; i_part++) {
      PDM_mesh_location_nodal_part_set_2d(mesh_loc,
                                          i_part,
                                          src_n_face       [i_part],
                                          src_face_edge_idx[i_part],
                                          src_face_vtx     [i_part],
                                          src_face_ln_to_gn[i_part],
                                          src_n_vtx        [i_part],
                                          src_vtx_coord    [i_part],
                                          src_vtx_ln_to_gn [i_part]);
    }
  }
  else {
    for (int i_part = 0; i_part < src_n_part; i_part++) {
      PDM_mesh_location_part_set_2d(mesh_loc,
                                    i_part,
                                    src_n_face       [i_part],
                                    src_face_edge_idx[i_part],
                                    src_face_edge    [i_part],
                                    src_face_ln_to_gn[i_part],
                                    src_n_edge       [i_part],
                                    src_edge_vtx     [i_part],
                                    src_n_vtx        [i_part],
                                    src_vtx_coord    [i_part],
                                    src_vtx_ln_to_gn [i_part]);
    }
  }
```

## Set some optional parameters

```{code-cell}
%%code_block -p exercise2 -i 7
  // Set the geometric tolerance (optional)
  double tolerance = 1e-3;
  PDM_mesh_location_tolerance_set(mesh_loc, tolerance);

  // Set the location preconditioning method (optional)
  PDM_mesh_location_method_set(mesh_loc,
                               PDM_MESH_LOCATION_OCTREE);

```

## Compute the localization

```{code-cell}
%%code_block -p exercise2 -i 8
  // Compute location
  PDM_mesh_location_compute(mesh_loc);

  // Dump elapsed and CPU times
  PDM_mesh_location_dump_times(mesh_loc);

```

## Results

Now that we have located the target points in the source mesh, we can exchange data between the two.
To complete this exercise, we will interpolate two fields from the source mesh to the target cloud:
1. a cell-based field: we can simply use the face global ids for such a field, and check it matches the location data.
2. a node-based field: we can use the node coordinates.

First, compute the spatially interpolated fields on the source side.
For the first field, the interpolation is straightforward : the target value is simply the same as the host source.
The second field interpolation is trickier as you will need the cell->vertex connectivity built during the location computation to link the interpolation weights to the appropriate source nodes.

### Interpolate the first field (cell-based)

```{code-cell}
%%code_block -p exercise2 -i 10
  // Interpolate first field (cell-based)
  double **src_send_field1 = malloc(sizeof(double *) * tgt_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {

    int         *src_to_tgt_idx          = NULL;
    PDM_g_num_t *points_gnum             = NULL;
    double      *points_coords           = NULL;
    double      *points_uvw              = NULL;
    int         *points_weights_idx      = NULL;
    double      *points_weights          = NULL;
    double      *points_dist2            = NULL;
    double      *points_projected_coords = NULL;

    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        0,
                                        i_part,
                                        &src_to_tgt_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords);

    int n_pts = src_to_tgt_idx[src_n_face[i_part]];

    src_send_field1[i_part] = malloc(sizeof(double) * n_pts);
    for (int i_elt = 0; i_elt < src_n_face[i_part]; i_elt++) {
      for (int i_pt = src_to_tgt_idx[i_elt]; i_pt < src_to_tgt_idx[i_elt+1]; i_pt++) {
        src_send_field1[i_part][i_pt] = src_face_ln_to_gn[i_part][i_elt];
      }
    }
  }

```

### Interpolate the second field (node-based)

```{code-cell}
%%code_block -p exercise2 -i 11
  // Interpolate second field (node-based)
  double **src_send_field2 = malloc(sizeof(double *) * tgt_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {

    int         *src_to_tgt_idx          = NULL;
    PDM_g_num_t *points_gnum             = NULL;
    double      *points_coords           = NULL;
    double      *points_uvw              = NULL;
    int         *points_weights_idx      = NULL;
    double      *points_weights          = NULL;
    double      *points_dist2            = NULL;
    double      *points_projected_coords = NULL;

    PDM_mesh_location_points_in_elt_get(mesh_loc,
                                        0,
                                        i_part,
                                        &src_to_tgt_idx,
                                        &points_gnum,
                                        &points_coords,
                                        &points_uvw,
                                        &points_weights_idx,
                                        &points_weights,
                                        &points_dist2,
                                        &points_projected_coords);

    int *cell_vtx_idx = NULL;
    int *cell_vtx     = NULL;
    PDM_mesh_location_cell_vertex_get(mesh_loc,
                                      i_part,
                                      &cell_vtx_idx,
                                      &cell_vtx);

    int n_pts = src_to_tgt_idx[src_n_face[i_part]];

    src_send_field2[i_part] = malloc(sizeof(double) * n_pts);
    for (int i_elt = 0; i_elt < src_n_face[i_part]; i_elt++) {
      for (int i_pt = src_to_tgt_idx[i_elt]; i_pt < src_to_tgt_idx[i_elt+1]; i_pt++) {
        src_send_field2[i_part][i_pt] = 0;

        int elt_n_vtx = cell_vtx_idx[i_elt+1] - cell_vtx_idx[i_elt];
        assert(points_weights_idx[i_pt+1] - points_weights_idx[i_pt] == elt_n_vtx);

        for (int i_vtx = 0; i_vtx < elt_n_vtx; i_vtx++) {
          int vtx_id = cell_vtx[cell_vtx_idx[i_elt] + i_vtx] - 1;
          src_send_field2[i_part][i_pt] += src_vtx_coord[i_part][3*vtx_id] * points_weights[points_weights_idx[i_pt] + i_vtx];
        }

      }
    }
  }

```


### Exchange the interpolated fields from source to target

Now, use the PartToPart object to exchange the interpolated fields from the source mesh to the target cloud.
This ParToPart object was built when computing the location and can be accessed from the MeshLocation object

```{code-cell}
%%code_block -p exercise2 -i 12
  // Get PartToPart object (it is now owned by the user)
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);

  // Initiate exchange of first field
  int request1 = -1;
  double **tgt_recv_field1 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) src_send_field1,
                         NULL,
              (void ***) &tgt_recv_field1,
                         &request1);

  // Initiate exchange of second field
  int request2 = -1;
  double **tgt_recv_field2 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                         1,
                         sizeof(double),
                         NULL,
        (const void  **) src_send_field2,
                         NULL,
              (void ***) &tgt_recv_field2,
                         &request2);


  // Finalize both exchanges
  PDM_part_to_part_iexch_wait(ptp, request1);
  PDM_part_to_part_iexch_wait(ptp, request2);

```

### Check the interpolated received on the target side

Finally, visualize the interpolated target fields.
(Beware of unlocated points!)

```{code-cell}
%%code_block -p exercise2 -i 13
  double **tgt_field[3];
  tgt_field[0] = malloc(sizeof(double *) * tgt_n_part);
  tgt_field[1] = malloc(sizeof(double *) * tgt_n_part);
  tgt_field[2] = malloc(sizeof(double *) * tgt_n_part);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    tgt_field[0][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);
    tgt_field[1][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);
    tgt_field[2][i_part] = malloc(sizeof(double) * tgt_n_vtx[i_part]);

    double *tgt_field1 = tgt_field[0][i_part];
    double *tgt_field2 = tgt_field[1][i_part];
    double *is_located = tgt_field[2][i_part];

    int n_located = PDM_mesh_location_n_located_get(mesh_loc,
                                                    0,
                                                    i_part);

    int *located = PDM_mesh_location_located_get(mesh_loc,
                                                 0,
                                                 i_part);

    int n_unlocated = PDM_mesh_location_n_unlocated_get(mesh_loc,
                                                        0,
                                                        i_part);

    int *unlocated = PDM_mesh_location_unlocated_get(mesh_loc,
                                                     0,
                                                     i_part);

    for (int i = 0; i < n_unlocated; i++) {
      int vtx_id = unlocated[i] - 1;
      is_located[vtx_id] = 0;
      tgt_field1[vtx_id] = -1;
      tgt_field2[vtx_id] = -1;
    }

    for (int i = 0; i < n_located; i++) {
      int vtx_id = located[i] - 1;
      is_located[vtx_id] = 1;
      tgt_field1[vtx_id] = tgt_recv_field1[i_part][i];
      tgt_field2[vtx_id] = tgt_recv_field2[i_part][i];

      double error = fabs(tgt_field2[vtx_id] - tgt_vtx_coord[i_part][3*vtx_id]);
      if (error > 1e-9) {
        printf("!! error vtx "PDM_FMT_G_NUM" : %e\n",
               tgt_vtx_ln_to_gn[i_part][vtx_id],
               error);
      }
    }

  }



  const char *field_name[] = {
    "field1",
    "field2",
    "is_located"
  };

  writer_wrapper(comm,
                 "visu",
                 "src_mesh",
                 src_n_part,
                 src_n_vtx,
                 src_vtx_coord,
                 src_vtx_ln_to_gn,
                 src_n_face,
                 src_face_edge_idx,
                 src_face_vtx,
                 src_face_ln_to_gn,
                 -1,
                 0,
                 NULL,
                 NULL,
                 "Ensight",
                 0,
                 NULL,
                 NULL,
                 0,//3,
                 field_name,
                 NULL);//src_field);


  writer_wrapper(comm,
                 "visu",
                 "tgt_mesh",
                 tgt_n_part,
                 tgt_n_vtx,
                 tgt_vtx_coord,
                 tgt_vtx_ln_to_gn,
                 tgt_n_face,
                 tgt_face_edge_idx,
                 tgt_face_vtx,
                 tgt_face_ln_to_gn,
                 -1,
                 0,
                 NULL,
                 NULL,
                 "Ensight",
                 0,
                 NULL,
                 NULL,
                 3,
                 field_name,
                 tgt_field);

```

### Free memory

```{code-cell}
%%code_block -p exercise2 -i 14
  // Free memory
  PDM_mesh_location_free(mesh_loc);

  PDM_part_to_part_free(ptp);

  for (int i_part = 0; i_part < src_n_part; i_part++) {
    free(src_vtx_coord    [i_part]);
    free(src_edge_vtx     [i_part]);
    free(src_face_edge_idx[i_part]);
    free(src_face_edge    [i_part]);
    free(src_face_vtx     [i_part]);
    free(src_vtx_ln_to_gn [i_part]);
    free(src_edge_ln_to_gn[i_part]);
    free(src_face_ln_to_gn[i_part]);
    free(src_send_field1  [i_part]);
    free(src_send_field2  [i_part]);
  }
  free(src_n_vtx        );
  free(src_n_edge       );
  free(src_n_face       );
  free(src_vtx_coord    );
  free(src_edge_vtx     );
  free(src_face_edge_idx);
  free(src_face_edge    );
  free(src_face_vtx     );
  free(src_vtx_ln_to_gn );
  free(src_edge_ln_to_gn);
  free(src_face_ln_to_gn);
  free(src_send_field1  ); // can be free'd right after PDM_part_to_part_iexch_wait(ptp, &request1);
  free(src_send_field2  ); // can be free'd right after PDM_part_to_part_iexch_wait(ptp, &request2);

  for (int i_part = 0; i_part < tgt_n_part; i_part++) {
    free(tgt_vtx_coord    [i_part]);
    free(tgt_edge_vtx     [i_part]);
    free(tgt_face_edge_idx[i_part]);
    free(tgt_face_edge    [i_part]);
    free(tgt_face_vtx     [i_part]);
    free(tgt_vtx_ln_to_gn [i_part]);
    free(tgt_edge_ln_to_gn[i_part]);
    free(tgt_face_ln_to_gn[i_part]);
    free(tgt_recv_field1  [i_part]);
    free(tgt_recv_field2  [i_part]);
    free(tgt_field[0][i_part]);
    free(tgt_field[1][i_part]);
    free(tgt_field[2][i_part]);
  }
  free(tgt_n_vtx        );
  free(tgt_n_edge       );
  free(tgt_n_face       );
  free(tgt_vtx_coord    );
  free(tgt_edge_vtx     );
  free(tgt_face_edge_idx);
  free(tgt_face_edge    );
  free(tgt_face_vtx     );
  free(tgt_vtx_ln_to_gn );
  free(tgt_edge_ln_to_gn);
  free(tgt_face_ln_to_gn);
  free(tgt_recv_field1  );
  free(tgt_recv_field2  );
  free(tgt_field[0]);
  free(tgt_field[1]);
  free(tgt_field[2]);
```


### Finalize

```{code-cell}
%%code_block -p exercise2 -i 15
  PDM_MPI_Finalize();

  if (i_rank == 0) {
    printf("-- End :)\n");
    fflush(stdout);
  }

  return EXIT_SUCCESS;
}

```

## Run the code
Moment of truth!

```{code-cell}
%merge_code_blocks -l c -p exercise2 -n 2 -c
```


## Visualize the results

```{code-cell}
%%visualize -nl -sv
visu/SRC_MESH.case
visu/TGT_MESH.case : is_located
```
