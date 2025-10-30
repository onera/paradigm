#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
# import argparse
import Pypdm.Pypdm as PDM


comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size


# Generate block-distributed parallelepided mesh
n_x      = 6
n_y      = 3
n_z      = 2
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

# Create partitioning object
n_domain_before = 1 # fixed
n_part_per_rank = 1 # fixed
i_part   = 0 # fixed
i_domain = 0 # fixed
part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT;
mpart = PDM.MultiPart(n_domain_before,
                      np.array([n_part_per_rank]).astype(np.intc),
                      0,
                      part_method,
                      1,
                      np.ones(1).astype(np.double),
                      comm)

renum_cell = bytes("PDM_PART_RENUM_CELL_NONE", 'ascii')
renum_face = bytes("PDM_PART_RENUM_FACE_NONE", 'ascii')
mpart.reordering_set(-1, # i_domain
                      renum_cell,
                      None,
                      renum_face)

mpart.dmesh_nodal_set(i_domain, dmn)

mpart.compute()

pmn = mpart.part_mesh_nodal_get(i_domain)

# print(pmn.n_part_get())

pcg_vtx = pmn.part_comm_graph_get(PDM._PDM_MESH_ENTITY_VTX)

entity_graph = pcg_vtx.get_entity_graph(i_part)

print(entity_graph)
# modulo 4 les infos contenues sont :
# 4k + 0. numero d'entité (num loc) <- pointList
# 4k + 1. target rank <-P0N0
# 4k + 2. target part (1 based)
# 4k + 3. target entity number <- pointlistdonor
print("iloc i_copy_rank  i_copy_rank_part   i_copy_rank_loc")
for i in range(int(entity_graph.size//4)):
  print(entity_graph[4*i],"            ", entity_graph[4*i+1],"              ",entity_graph[4*i+2],"               ", entity_graph[4*i+3])

if i_rank == 0:
  print("End :)")
