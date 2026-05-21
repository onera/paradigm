#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

# Generate block-distributed parallelepided mesh
n_x      = 6
n_y      = 3
n_z      = 1
lengthx  = 1.
xmin     = 0.
ymin     = 0.
zmin     = 0.
elt_type = PDM._PDM_MESH_NODAL_TRIA3
order    = 1
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

dmn.generate_distribution()

# Create partitioning object
n_domain_before = 1
n_part_per_rank = 1
i_domain = 0
i_part   = 0
part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT;
mpart = PDM.MultiPart(n_domain_before,
                      np.array([n_part_per_rank]).astype(np.intc),
                      0,
                      part_method,
                      1,
                      np.ones(1).astype(np.double),
                      comm)
mpart.dmesh_nodal_set(i_domain, dmn)
mpart.compute()

pmn = mpart.part_mesh_nodal_get(i_domain)

p_vtx_coord_part = mpart.vtx_coord_get(i_domain,
                                       i_part)
p_vtx_coord = [p_vtx_coord_part]
p_nvtx = [p_vtx_coord[0].shape[0]//3]
_, face_edge = mpart.connectivity_get(i_domain,
                                      i_part,
                                      PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
p_face_edge = [face_edge]
_, edge_vtx = mpart.connectivity_get(i_domain,
                                     i_part,
                                     PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)
p_edge_vtx_idx = [2*np.arange(edge_vtx.shape[0]//2+1, dtype=np.int32)]
p_edge_vtx     = [edge_vtx]

pcg_vtx = pmn.part_comm_graph_vtx_get()

# Generate tmp fields array
p_vtx_coord_prev = [np.copy(p_vtx_coord_part)]

# Generate vtx group
p_vtx_frozen = [np.array([i+1 for i in range(p_vtx_coord[0].shape[0]//3) if p_vtx_coord[0][3*i] < 0.4], dtype=np.int32)]

# Generate edge pcg
pcg_edge = PDM.pcg_entity1_to_entity2(pcg_vtx,
                                      p_nvtx,
                                      p_edge_vtx_idx,
                                      p_edge_vtx)

# Laplacian parameters
damping = 0.9
n_iter  = 30
tol     = 0.1
iter    = 0
eps     = 1e15;

while (iter < n_iter and eps > tol):

  # Generate edge weights
  p_edge_weight = PDM.compute_cotangent_weights(p_vtx_coord,
                                                p_face_edge,
                                                p_edge_vtx,
                                                pcg_edge)

  # Swap current and previous fields
  tmp_swap         = p_vtx_coord
  p_vtx_coord_prev = p_vtx_coord
  p_vtx_coord      = tmp_swap

  # Laplacian smoothing
  eps = PDM.laplacian_smoothing_fields_one_iteration(comm,
                                 p_vtx_frozen,
                                 pcg_vtx,
                                 p_edge_vtx,
                                 p_edge_weight,
                                 pcg_edge,
                                 damping,
                                 n_iter,
                                 tol,
                                 3,
                                 p_vtx_coord_prev,
                                 p_vtx_coord)

  iter += 1

if i_rank == 0:
  print("End :)")
