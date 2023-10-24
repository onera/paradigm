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

In this second exercise we will focus on the **Mesh Location** feature.
It consists in computing the location of one or more partitioned point clouds (referred to as the *targets*) inside a partitioned mesh (referred to as the *source*).
<!-- Some target points may be unlocated if they lie outside the source mesh. -->
A mapping between the source mesh elements and the target points they contain is computed, which consists in
  - geometric data (distances, barycentric and parametric coordinates, ...) ;
  - an MPI communication graph as the associated entities are, in general, distributed on different processes.

This mapping is typically used for interpolating data from the source mesh to the point clouds in applications such as coupling between non-matching grids, <span style="color:red">**(autres exemples...)**</span>.


The aim of this exercise is to perform such an interpolation.

Your task is to fill in the empty code cells using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/prepro_algo/mesh_location.html#mesh-location).

*Note: For easier visualization, we will study a two-dimensional case but the feature is also available in three dimensions.*

+++

## Load magic commands
As usual we start by loading the custom magic commands.

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



+++

## Include the required headers and initialize MPI

To start, we include the required C headers and initialize MPI.

```{code-cell}
%%code_block -p exercise_2 -i 1

// Required headers
#include "pdm_mesh_location.h"
#include "pdm_generate_mesh.h"

int main(int argc, char *argv[])
{
  // Initialize MPI
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
```

## Generate a partitioned "source" mesh

We then generate the partitioned source mesh.

By now you should be capable of partitioning a mesh using **ParaDiGM** (if not, you should definitely take a look at [**Exercise 1**](../02_Exercise_1/exercise_1.ipynb)).
To gain some time, let's use the [*PDM_generate_mesh*](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_formation/user_manual/simple_mesh_gen/generate_mesh.html) service to generate a partitioned mesh in a single function call.

Here we generate a square mesh composed of polygonal elements.

```{code-cell}
%%code_block -p exercise_2 -i 2

  // Generate partitioned source mesh
  PDM_g_num_t src_n_vtx_seg = 10; // number of vertices along each side of the square
  int         src_n_part    = 1;  // number of partitions per MPI rank

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
                                   0.8,
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

We then generate a second partitioned mesh.
We will use its vertices as a target point cloud.
This second mesh is deliberately offset so that some target points lie outside the source mesh.
These points may not be located.
We will see later how to deal with these *unlocated* points.

```{code-cell}
%%code_block -p exercise_2 -i 3

  // Generate target source mesh
  PDM_g_num_t tgt_n_vtx_seg = 8;    // number of vertices along each side of the square
  int         tgt_n_part    = 1;    // number of partitions per MPI rank
  double      tgt_xmin      = 0.25; // x-offset
  double      tgt_ymin      = 0.25; // y-offset

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
                                   0.,
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

## Create the `PDM_mesh_location_t` object

Now that we have all the required inputs, let's create an instance of the `PDM_mesh_location_t` structure.

```{code-cell}
%%code_block -p exercise_2 -i 4

  // Create the PDM_mesh_location_t object
  // EXO
  PDM_mesh_location_t *mesh_loc = PDM_mesh_location_create(1,
                                                           comm,
                                                           PDM_OWNERSHIP_KEEP);

```

## Set the target point cloud

Now let's provide the target point cloud to the structure.

```{code-cell}
%%code_block -p exercise_2 -i 5

  // Set target point cloud
  // EXO
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

Now let's provide the source mesh to the structure.

Here you have essentially two options:
- you can either define the mesh with "nodal" connectivity (i.e. Finite-Element style)
- or with "descending" connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -p exercise_2 -i 6

  // Set source mesh
  // EXO
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

The location algorithm uses a preconditioning stage which consists in associating candidate elements and points before computed the exact location.
Three preconditioning methods are available:
- `PDM_MESH_LOCATION_OCTREE`
- `PDM_MESH_LOCATION_DBBTREE`
- `PDM_MESH_LOCATION_LOCATE_ALL_TGT`

The first two methods use bounding boxes and distributed tree data structures to find the candidate pairs efficiently.
Theses boxes can be expanded using a **relative geometric tolerance** allowing for safer candidate detection.
The third method uses a combination of both trees to ensure all target points are "located", i.e. associated to the nearest source element.

(By default, the `PDM_MESH_LOCATION_OCTREE` method is used, with a relative tolerance equal to zero.)

```{code-cell}
%%code_block -p exercise_2 -i 7

  // Set the location preconditioning method (optional)
  // EXO
  PDM_mesh_location_method_set(mesh_loc,
                               PDM_MESH_LOCATION_OCTREE);

  // Set the geometric tolerance (optional)
  // EXO
  double tolerance = 1e-6;
  PDM_mesh_location_tolerance_set(mesh_loc, tolerance);

```

## Compute the localization

Now that everything is ready, we can compute the localization.
Once the calculation is complete, we can display the elapsed time and CPU time.

```{code-cell}
%%code_block -p exercise_2 -i 8

  // Compute location
  // EXO
  PDM_mesh_location_compute(mesh_loc);

  // Dump elapsed and CPU times
  // EXO
  PDM_mesh_location_dump_times(mesh_loc);

```

## Results




<span style="color:red">
Si vous vous rappelez, pour construire un objet Part-to-Part, on doit spécifier:
- les partitions côté "1"
- les partitions côté "2"
- le lien/graphe 1 -> 2 (en g_num)

Ici, 1 = source et 2 = cible.
Le lien 1->2 est donc "la liste des points cibles contenus dans chaque élement du maillage source" (pt_in_elt aka src_to_tgt)

Rappel data_def_order

recouvrir échange field1 par calcul src_field2
</span>

Now that we have located the target points in the source mesh, we can exchange data between the two.
To complete this exercise, we will interpolate two fields from the source mesh to the target cloud:
  1. a cell-based field: we can simply use the face global ids for such a field, and check it matches the location data.
  2. a node-based field: we can use the node coordinates.

First, compute the spatially interpolated fields on the source side.
For the first field, the interpolation is straightforward : the target value is simply the same as the host source.
The second field interpolation is trickier as you will need the cell->vertex connectivity built during the location computation to link the interpolation weights to the appropriate source nodes.

### Retrieve the `PDM_part_to_part_t` instance

<span style="color:red">EXPLIQUER</span>

```{code-cell}
%%code_block -p exercise_2 -i 9

  // Get PDM_part_to_part_t object (it is now owned by the user)
  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_location_part_to_part_get(mesh_loc,
                                     0,
                                     &ptp,
                                     PDM_OWNERSHIP_USER);
```

### Exchange the first field from source to target

<span style="color:red">EXPLIQUER</span>

```{code-cell}
%%code_block -p exercise_2 -i 10

  // Initiate exchange of first field (source elements global ids)
  int request1 = -1;
  PDM_g_num_t **tgt_recv_field1 = NULL;
  PDM_part_to_part_iexch(ptp,
                         PDM_MPI_COMM_KIND_P2P,
                         PDM_STRIDE_CST_INTERLACED,
                         PDM_PART_TO_PART_DATA_DEF_ORDER_PART1,
                         1,
                         sizeof(PDM_g_num_t),
                         NULL,
        (const void  **) src_face_ln_to_gn,
                         NULL,
              (void ***) &tgt_recv_field1,
                         &request1);

```

### Interpolate the second field (node-based)

<span style="color:red">EXPLIQUER</span>

```{code-cell}
%%code_block -p exercise_2 -i 11

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


### Exchange the second interpolated fields

Now, use the [`PDM_part_to_p<!-- art_t`](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/comm_graph/ptp.html) object to exchange the interpolated fields from the source mesh to the target cloud.
This `PDM_part_to_part_t` object was built when computing the location and can be accessed from the `PDM_Mesh_location_t` object. -->

<span style="color:red">EXPLIQUER</span>

```{code-cell}
%%code_block -p exercise_2 -i 12

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
```


### Check the interpolated received on the target side

Finally, check and visualize the received interpolated target fields.

<span style="color:red">
EXPLIQUER
- wait
- comment adresser tableaux reçus

EXO
- wait
- remplir tableaux pré-alloués
</span>

*(Watch out for unlocated points!)*

```{code-cell}
%%code_block -p exercise_2 -i 13

  // Finalize both exchanges
  PDM_part_to_part_iexch_wait(ptp, request1);
  PDM_part_to_part_iexch_wait(ptp, request2);

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


  double **src_elt_field_values = malloc(sizeof(double *) * src_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {
    src_elt_field_values[i_part] = malloc(sizeof(double) * src_n_face[i_part]);
    for (int i_elt = 0; i_elt < src_n_face[i_part]; i_elt++) {
      src_elt_field_values[i_part][i_elt] = (double) src_face_ln_to_gn[i_part][i_elt];
    }
  }

  double **src_vtx_field_values = malloc(sizeof(double *) * src_n_part);
  for (int i_part = 0; i_part < src_n_part; i_part++) {
    src_vtx_field_values[i_part] = malloc(sizeof(double) * src_n_vtx[i_part]);
    for (int i_vtx = 0; i_vtx < src_n_vtx[i_part]; i_vtx++) {
      src_vtx_field_values[i_part][i_vtx] = src_vtx_coord[i_part][3*i_vtx];
    }
  }

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
                 1,
                 &field_name[0],
                 &src_elt_field_values,
                 1,
                 &field_name[1],
                 &src_vtx_field_values);

  for (int i_part = 0; i_part < src_n_part; i_part++) {
    free(src_elt_field_values[i_part]);
    free(src_vtx_field_values[i_part]);
  }
  free(src_elt_field_values);
  free(src_vtx_field_values);


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

## Free memory

```{code-cell}
%%code_block -p exercise_2 -i 14
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


## Finalize

```{code-cell}
%%code_block -p exercise_2 -i 15
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
%merge_code_blocks -l c -p exercise_2 -n 2 -c
```


## Visualize the results

```{code-cell}
%%visualize -nl -sv
visu/SRC_MESH.case
visu/TGT_MESH.case : is_located : points
```
