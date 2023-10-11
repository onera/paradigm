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
```

The goal of this exercise is to get used
Your job is to fill the code blocks left blank using the API referenced [here](https://numerics.gitlab-pages.onera.net/mesh/paradigm/dev_doc_pretty/user_manual/prepro_algo/index.html#python-api).

+++

## Begin the exercise

```{code-cell}
%%code_block -l python -p exercise2 -i 1
# Load modules
import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM
from util_visu import visu_2d
```

## Generate a partitioned "source" mesh

```{code-cell}
%%code_block -l python -p exercise2 -i 2
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
%%code_block -l python -p exercise2 -i 3
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
%%code_block -l python -p exercise2 -i 4
# (...)
```

## Set the target point cloud

```{code-cell}
%%code_block -l python -p exercise2 -i 5
# (...)
```

## Set the source mesh
Here you have essentially two options :
- you can either define the mesh with "nodal" connectivity (i.e. Finite-Element style)
- or with "descending" connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -l python -p exercise2 -i 6
# (...)
```

## Set some optional parameters

```{code-cell}
%%code_block -l python -p exercise2 -i 7
# Geometric tolerance
# (...)

# Preconditioning method
# (...)
```

## Compute the localization

```{code-cell}
%%code_block -l python -p exercise2 -i 8
# (...)
```

## Dump the elapsed and CPU times

```{code-cell}
%%code_block -l python -p exercise2 -i 9
# (...)
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
%%code_block -l python -p exercise2 -i 10
# (...)
```

### Interpolate the second field (node-based)

```{code-cell}
%%code_block -l python -p exercise2 -i 11
# (...)
```


### Exchange the interpolated fields from source to target

Now, use the PartToPart object to exchange the interpolated fields from the source mesh to the target cloud.
This ParToPart object was built when computing the location and can be accessed from the MeshLocation object

```{code-cell}
%%code_block -l python -p exercise2 -i 10
# Get PartToPart object
# (...)

# Initiate exchange of first field
# (...)

# Initiate exchange of second field
# (...)

# Finalize both exchanges
# (...)
```

### Check the interpolated received on the target side

Finally, visualize the interpolated target fields.
(Beware of unlocated points!)

```{code-cell}
%%code_block -l python -p exercise2 -i 13
# (...)

# Export for visualization
# (...)
```

Finalize?



## Run the code
Moment of truth!

```{code-cell}
%merge_code_blocks -l python -p exercise2 -n 2 -v -c
```


## Visualize the results

```{code-cell}
%%visualize
/stck/bandrieu/workspace/formations/trainings/cwipi/cwipi_writer/exercise1_code1_code2/CHR.case : s_super~fancy~field1
/stck/bandrieu/workspace/formations/trainings/cwipi/cwipi_writer/exercise1_code2_code1/CHR.case : r_super~fancy~field1
```
