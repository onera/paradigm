# Exercise 1 : Mesh partitioning

Let's get started using `ParaDiGM` by focusing on its mesh partitioning solutions.
It is a mandatory step for parallel numerical simulations. First a mesh is generated.
Then it is partitionned in subdomains. Those will be mapped onto the processors of a parallel machine.

## Generate the initial mesh

In this section, `ParaDiGM` tools are used to generate the simple mesh for exercise 1 : a cube.

```{code-cell}
%%code_block -l python -p exercise_1 -i 1

#!/usr/bin/env python

import numpy as np
import mpi4py.MPI as MPI
import Pypdm.Pypdm as PDM

 # Initialize MPI environment
 comm   = MPI.COMM_WORLD
 n_rank = MPI.COMM_WORLD.size
 i_rank = MPI.COMM_WORLD.rank

# Generate block-distributed parallelepided mesh
  n_x      = 10
  n_y      = 10
  n_z      = 10
  lengthx  = 1.
  xmin     = 0.
  ymin     = 0.
  zmin     = 0.
  elt_type = PDM._PDM_MESH_NODAL_TETRA4
  order    = 1 # call PDM_dcube_nodal_gen_ordering_set if order > 1
  dcube = PDM.DCubeNodalGenerator(n_x,
                                  n_y,
                                  n_z,
                                  lengthx,
                                  xmin,
                                  ymin,
                                  zmin,
                                  elt_type,
                                  order,
                                  comm)

  dcube.compute()

  dmn = dcube.get_dmesh_nodal()

  PDM.generate_distribution(dmn)
```

Now, we have a block-distributed cube mesh.

## Mesh partitioning

For mesh partitioning, as for all other `ParaDiGM` features, there are 5 main steps:
- **create** the feature object
- **set** the data necessary to operate with that feature
- **compute**, operate the algorithm of the feature
- **get**, retreive the ouput of the algorithm
- **free** the memory allocated to operate the feature

Following this logic, let's start **creating** the mesh partitioning object.

*Note : since this is a simple example, we will not focus on the concepts of zone and part. To understand those concepts, see Annex 1*

```{code-cell}
%%code_block -l python -p exercise_1 -i 2

# Create partitioning object
n_zone = 1 # fixed
n_part = 1 # fixed
i_part = 0 # fixed
i_zone = 0 # fixed
part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT;
mpart = PDM.MultiPart(n_zone,                             # Number of zones
                      np.array([n_part]).astype(np.intc), # Number of partitions per zone
                      0,                                  # PDM_FALSE (do not fuse zones)
                      part_method,                        # Partitioning method
                      1,                                  # PDM_PART_SIZE_HOMOGENEOUS (subdomains are equaly balanced)
                      None,                               # Weight (in %) of each partition in heterogeneous case
                      comm)                               # MPI communicator

```

Here, we chose the partition the cube with the Hilbert method. This method implemented in `ParaDiGM` does not ensure the subdomain to be connected.
This method is favored within the `ParaDiGM` algorithms since it provides quickly a good load balance. To ensure the partitions are connected use
`PDM_SPLIT_DUAL_WITH_PARMETIS` or `PDM_SPLIT_DUAL_WITH_PTSCOTCH` which call the external libraries ParMETIS and PT-Scotch.

After the partitioning, it is possible the reorder the mesh entities. In this simple example, we won't do that.

```{code-cell}
%%code_block -l python -p exercise_1 -i 3

renum_cell = bytes("PDM_PART_RENUM_CELL_NONE", 'ascii')
renum_face = bytes("PDM_PART_RENUM_FACE_NONE", 'ascii')
mpart.multipart_set_reordering(-1,         # All zones
                               renum_cell,
                               renum_face,
                               None)

```

We go on to the **set** step, in which we will provide the block-distributed cube mesh to the mesh partitioning object `mpart`.

```{code-cell}
%%code_block -l python -p exercise_1 -i 4

mpart.multipart_register_dmesh_nodal(i_zone, dmn)
```

Now we can run the partitioning algorithm. This is the so calle **compute** setp.

```{code-cell}
%%code_block -l python -p exercise_1 -i 5

mpart.multipart_run_ppart()
```

## Get the partitionned mesh

The aim of this whole process is to retreive a partitionned mesh. Let's move on to the **get** step.
Depending on the numerical method, the mesh is not described and stored in the same style. For this next part of the exercise,
you can choose if you want to retreive the partitionned cube in nodal or decending connectivity. Now move on to the chosen sub-section.

### Nodal connectivity (i.e. Finite-Element style)

```{code-cell}
%%code_block -l python -p exercise_1 -i 3

```

### Descending connectivity (i.e. Finite-Volume style)

```{code-cell}
%%code_block -l python -p exercise_1 -i 4

```

## Annex 1
