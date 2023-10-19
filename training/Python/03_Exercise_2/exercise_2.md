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

Your job is to fill the code cells left blank using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/prepro_algo/index.html#python-api).

+++

## Load the required Python modules

```{code-cell}
%%code_block -p exercise_2 -i 1
# Load modules
import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

```

## Generate a partitioned "source" mesh


By now, you know how to partition a mesh.

```{code-cell}
%%code_block -p exercise_2 -i 2
# Generate source mesh
src_n_vtx_seg = 10
src_n_part    = 1

src_mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                            elt_type    = PDM._PDM_MESH_NODAL_POLY_2D,
                                            xmin        = 0.,
                                            ymin        = 0.,
                                            zmin        = 0.,
                                            lengthx     = 1.,
                                            lengthy     = 1.,
                                            n_x         = src_n_vtx_seg,
                                            n_y         = src_n_vtx_seg,
                                            n_part      = src_n_part,
                                            part_method = PDM._PDM_SPLIT_DUAL_WITH_PARMETIS)

src_n_vtx         = src_mesh["pn_vtx"]         # Number of vertices in each partition
src_n_face        = src_mesh["pn_face"]        # Number of faces in each partition
src_vtx_coord     = src_mesh["pvtx_coord"]     #
src_face_vtx_idx  = src_mesh["pface_edge_idx"] #
src_face_vtx      = src_mesh["pface_vtx"]      #
src_face_edge     = src_mesh["pface_edge"]     #
src_edge_vtx      = src_mesh["pedge_vtx"]      #
src_vtx_ln_to_gn  = src_mesh["pvtx_ln_to_gn"]  #
src_face_ln_to_gn = src_mesh["pface_ln_to_gn"] #

```

## Generate a partitioned "target" mesh
We will use its vertices as a point cloud.

```{code-cell}
%%code_block -p exercise_2 -i 3
# Generate target mesh
tgt_n_vtx_seg = 10
tgt_n_part    = 1

tgt_mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                            elt_type    = PDM._PDM_MESH_NODAL_QUAD4,
                                            xmin        = 0.3,
                                            ymin        = 0.3,
                                            zmin        = 0,
                                            lengthx     = 1.,
                                            lengthy     = 1.,
                                            n_x         = tgt_n_vtx_seg,
                                            n_y         = tgt_n_vtx_seg,
                                            n_part      = tgt_n_part,
                                            part_method = PDM._PDM_SPLIT_DUAL_WITH_PARMETIS)
tgt_n_vtx         = tgt_mesh["pn_vtx"]
tgt_n_face        = tgt_mesh["pn_face"]
tgt_vtx_coord     = tgt_mesh["pvtx_coord"]
tgt_face_vtx_idx  = tgt_mesh["pface_edge_idx"]
tgt_face_vtx      = tgt_mesh["pface_vtx"]
tgt_vtx_ln_to_gn  = tgt_mesh["pvtx_ln_to_gn"]
tgt_face_ln_to_gn = tgt_mesh["pface_ln_to_gn"]

```

## Create the `MeshLocation` object

```{code-cell}
%%code_block -p exercise_2 -i 4
# Create MeshLocation instance
mesh_loc = PDM.MeshLocation(1,
                            comm)

```

## Set the target point cloud

```{code-cell}
%%code_block -p exercise_2 -i 5
# Set target point cloud
mesh_loc.n_part_cloud_set(0,
                          tgt_n_part)

for i_part in range(tgt_n_part):
  mesh_loc.cloud_set(0,
                     i_part,
                     tgt_vtx_coord   [i_part],
                     tgt_vtx_ln_to_gn[i_part])
```

## Set the source mesh
Here you have essentially two options :
- you can either define the mesh with "nodal" connectivity (i.e. Finite-Element style)
- or with "descending" connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -p exercise_2 -i 6
# Set source mesh
nodal = True
mesh_loc.mesh_n_part_set(src_n_part)
for i_part in range(src_n_part):
  if nodal:
    mesh_loc.nodal_part_set_2d(i_part,
                               src_face_vtx_idx [i_part],
                               src_face_vtx     [i_part],
                               src_face_ln_to_gn[i_part],
                               src_vtx_coord    [i_part],
                               src_vtx_ln_to_gn [i_part])
  else:
    mesh_loc.part_set_2d(i_part,
                         src_face_vtx_idx [i_part],
                         src_face_edge    [i_part],
                         src_face_ln_to_gn[i_part],
                         src_edge_vtx     [i_part],
                         src_vtx_coord    [i_part],
                         src_vtx_ln_to_gn [i_part])

```

## Set some optional parameters

```{code-cell}
%%code_block -p exercise_2 -i 7
# Geometric tolerance
mesh_loc.tolerance_set(1e-6)

# Preconditioning method
mesh_loc.method_set(0)

```

## Compute the localization

```{code-cell}
%%code_block -p exercise_2 -i 8
# Compute localization
mesh_loc.compute()

# Dump elapsed and CPU times
mesh_loc.dump_times()

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
%%code_block -p exercise_2 -i 10
# Interpolate first field
src_send_field1 = []

for i_part in range(src_n_part):
  src_result = mesh_loc.points_in_elt_get(0, i_part)
  src_to_tgt_idx = src_result["elt_pts_inside_idx"]
  n_pts = src_to_tgt_idx[src_n_face[i_part]]

  field1 = np.zeros(n_pts, dtype=np.double)
  for i_elt in range(src_n_face[i_part]):
    for i_pt in range(src_to_tgt_idx[i_elt], src_to_tgt_idx[i_elt+1]):
      field1[i_pt] = src_face_ln_to_gn[i_part][i_elt]
  src_send_field1.append(field1)
```

### Interpolate the second field (node-based)

```{code-cell}
%%code_block -p exercise_2 -i 11
# Interpolate second field
src_send_field2 = []

for i_part in range(src_n_part):
  src_result = mesh_loc.points_in_elt_get(0, i_part)
  src_to_tgt_idx = src_result["elt_pts_inside_idx"]
  n_pts = src_to_tgt_idx[src_n_face[i_part]]

  # Interpolate second field (node-based)
  src_connect = mesh_loc.cell_vertex_get(i_part)
  src_cell_vtx_idx = src_connect["cell_vtx_idx"]
  src_cell_vtx     = src_connect["cell_vtx"]

  weights_idx = src_result["points_weights_idx"]
  weights     = src_result["points_weights"]

  field2 = np.zeros(n_pts, dtype=np.double)
  for i_elt in range(src_n_face[i_part]):
    for i_pt in range(src_to_tgt_idx[i_elt], src_to_tgt_idx[i_elt+1]):
      field2[i_pt] = 0

      elt_n_vtx = src_cell_vtx_idx[i_elt+1] - src_cell_vtx_idx[i_elt]
      assert(weights_idx[i_pt+1] - weights_idx[i_pt] == elt_n_vtx)
      for i_vtx in range(elt_n_vtx):
        vtx_id = src_cell_vtx[src_cell_vtx_idx[i_elt] + i_vtx] - 1
        field2[i_pt] += src_vtx_coord[i_part][3*vtx_id] * weights[weights_idx[i_pt] + i_vtx]

  src_send_field2.append(field2)

```


### Exchange the interpolated fields from source to target

Now, use the PartToPart object to exchange the interpolated fields from the source mesh to the target cloud.
This ParToPart object was built when computing the location and can be accessed from the MeshLocation object

```{code-cell}
%%code_block -p exercise_2 -i 12
# Get PartToPart object
ptp = mesh_loc.part_to_part_get(0)

# Initiate exchange of first field
src_stride = 1
request1 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                     PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                     src_send_field1,
                     part1_stride=src_stride,
                     interlaced_str=True)

# Initiate exchange of second field
request2 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                     PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                     src_send_field2,
                     part1_stride=src_stride,
                     interlaced_str=True)

# Finalize both exchanges
tgt_stride, tgt_recv_field1 = ptp.wait(request1)
tgt_stride, tgt_recv_field2 = ptp.wait(request2)

```

### Check the interpolated received on the target side

Finally, visualize the interpolated target fields.
(Beware of unlocated points!)

```{code-cell}
%%code_block -p exercise_2 -i 13
# Check interpolated fields
pis_located = []
ptgt_field1 = []
ptgt_field2 = []
for i_part in range(tgt_n_part):
  located_tgt = mesh_loc.located_get(0, i_part)

  is_located =  np.zeros(tgt_n_vtx[i_part], dtype=bool)
  tgt_field1 = -np.ones(tgt_n_vtx[i_part], dtype=np.double)
  tgt_field2 = -np.ones(tgt_n_vtx[i_part], dtype=np.double)

  for i, i_vtx in enumerate(located_tgt):
    is_located[i_vtx-1] = True
    tgt_field1[i_vtx-1] = tgt_recv_field1[i_part][i]
    tgt_field2[i_vtx-1] = tgt_recv_field2[i_part][i]
    error = abs(tgt_recv_field2[i_part][i] - tgt_vtx_coord[i_part][3*(i_vtx-1)])
    if error > 1e-9:
      print(f"!! error vtx {tgt_vtx_ln_to_gn[i_part][i_vtx]} : {error}")

  pis_located.append(is_located)
  ptgt_field1.append(tgt_field1)
  ptgt_field2.append(tgt_field2)

# Export for visualization
PDM.writer_wrapper(comm,
                   "visu",
                   "src_mesh",
                   src_vtx_coord,
                   src_vtx_ln_to_gn,
                   src_face_vtx_idx,
                   src_face_vtx,
                   src_face_ln_to_gn)

PDM.writer_wrapper(comm,
                   "visu",
                   "tgt_mesh",
                   tgt_vtx_coord,
                   tgt_vtx_ln_to_gn,
                   tgt_face_vtx_idx,
                   tgt_face_vtx,
                   tgt_face_ln_to_gn,
                   vtx_fields={
                   "is_located" : pis_located,
                   "field1"     : ptgt_field1,
                   "field2"     : ptgt_field2})

```

Finalize?



## Run the code
Moment of truth!

```{code-cell}
%merge_code_blocks -l python -p exercise_2 -n 2 -c
```


## Visualize the results

```{code-cell}
%%visualize -nl -sv
visu/SRC_MESH.case
visu/TGT_MESH.case : is_located
```
