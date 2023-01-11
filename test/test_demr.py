#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import sys
import Pypdm.Pypdm as PDM
from Pypdm.Pypdm import npy_pdm_gnum_dtype

import sys
sys.path.append("/stck/bandrieu/Public/adaptation/cavity_operator3_v3/")
from mod_vtk import *

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size


# n_part = 1

# sphere_mesh = PDM.sphere_surf_icosphere_gen(comm,
#                                             4,
#                                             0.,
#                                             0.,
#                                             0.,
#                                             1.)
# print(sphere_mesh)

# distrib_vtx  = sphere_mesh["distrib_vtx"]
# distrib_face = sphere_mesh["distrib_face"]

# dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank]
# dn_face = distrib_face[i_rank+1] - distrib_face[i_rank]


"""
sphere_mesh = PDM.sphere_surf_icosphere_gen_part(comm,
                                                 4,
                                                 0.,
                                                 0.,
                                                 0.,
                                                 1.,
                                                 n_part,
                                                 PDM._PDM_SPLIT_DUAL_WITH_HILBERT)
# print(sphere_mesh)

for i in range(n_part):
  coord = []
  for j in range(sphere_mesh["pn_vtx"][i]):
    coord.append([x for x in sphere_mesh["pvtx_coord"][i][3*j:3*(j+1)]])

  connec = []
  for j in range(sphere_mesh["pn_face"][i]):
    start = sphere_mesh["pface_vtx_idx"][i][j]
    end   = sphere_mesh["pface_vtx_idx"][i][j+1]
    connec.append([k for k in sphere_mesh["pface_vtx"][i][start:end] - 1])

  vtk_write_polydata("sphere_mesh_rank%d_part%d.vtk" % (i_rank, i),
                     coord,
                     connec)
"""

# Generate a distributed surface mesh
n = 3
length = 1.
dcube = PDM.DCubeNodalGenerator(n, n, n,
                                length,
                                0., 0., 0.,
                                PDM._PDM_MESH_NODAL_TRIA3,
                                1,
                                comm)

dcube.compute()

dmn_capsule = dcube.get_dmesh_nodal()
print("[{}] dmn OK".format(i_rank))

print(dmn_capsule.dmesh_nodal_get_sections(PDM._PDM_GEOMETRY_KIND_SURFACIC,
                                           comm))

# Split the mesh
mpart = PDM.MultiPart(1,
                      np.ones(1).astype(np.intc),
                      0,
                      PDM._PDM_SPLIT_DUAL_WITH_HILBERT,
                      1,#PDM_PART_SIZE_HOMOGENEOUS,
                      np.ones(1).astype(np.double),
                      comm)

mpart.multipart_register_dmesh_nodal(0, dmn_capsule)

mpart.multipart_run_ppart()
print("[{}] mpart OK".format(i_rank))

# Extract a random set of edges and redistribute
extract_fraction = 0.25

dims = mpart.multipart_dim_get(0, 0)
n_edge = dims["n_face"]

vals = mpart.multipart_val_get(0, 0)
# edge_vtx = vals["np_face_vtx"]
edges = mpart.multipart_connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)
print("edges : {}".format(edges))
edge_vtx = edges["np_entity1_entity2"]

edge_ln_to_gn = mpart.multipart_ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_EDGE)["np_entity_ln_to_gn"]
n_edge = len(edge_ln_to_gn)

print("n_edge = {}".format(n_edge))
print("edge_ln_to_gn = {}".format(edge_ln_to_gn))

faces = mpart.multipart_connectivity_get(0, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)
face_edge_idx = faces["np_entity1_entity2_idx"]
face_edge     = faces["np_entity1_entity2"]

n_face = len(face_edge_idx)-1

face_vtx = PDM.compute_face_vtx_from_face_and_edge(n_face,
                                                   face_edge_idx,
                                                   face_edge,
                                                   edge_vtx)

# n_edge = len(edges["np_entity1_entity2_idx"]) - 1

print("n_edge = {}".format(n_edge))
print("edge_vtx = {}".format(edge_vtx))

# get global number of edges
lmax_edge_ln_to_gn = np.amax(edge_ln_to_gn)

gmax_edge_ln_to_gn = comm.allreduce(lmax_edge_ln_to_gn, op=MPI.MAX)

# if i_rank == 0:
#   print("gmax_edge_ln_to_gn = {}".format(gmax_edge_ln_to_gn))

step_edge = int(extract_fraction * gmax_edge_ln_to_gn)

select_edge = []
for i in range(n_edge):
  if edge_ln_to_gn[i] % step_edge == 0:
    select_edge.append(i)
    # print("[{}] select edge {}".format(i_rank, edge_ln_to_gn[i]))


# Generate segments orthogonal to extracted edges
# use extract_part!
vtx_coord = vals["np_vtx_coord"]
vtx_ln_to_gn = mpart.multipart_ln_to_gn_get(0, 0, PDM._PDM_MESH_ENTITY_VERTEX)["np_entity_ln_to_gn"]
n_vtx = len(vtx_ln_to_gn)

edge_mid_coord = np.zeros(3*n_edge)
for i in range(n_edge):
  for j in range(2):
    vtx_id = edge_vtx[2*i+j] - 1
    edge_mid_coord[3*i:3*(i+1)] = edge_mid_coord[3*i:3*(i+1)] + 0.5*vtx_coord[3*vtx_id:3*(vtx_id+1)]


extrp = PDM.ExtractPart(1,
                        1,
                        1,
                        PDM._PDM_EXTRACT_PART_KIND_LOCAL,#REEQUILIBRATE,
                        PDM._PDM_SPLIT_DUAL_WITH_HILBERT,
                        1,
                        comm)

extrp.selected_lnum_set(0, np.array(select_edge).astype(np.intc))

extrp.part_set(0,
               0, # n_cell
               0, # n_face
               n_edge,
               n_vtx,
               None, # cell_face_idx
               None, # cell_face
               None, # face_edge_idx
               None, # face_edge
               edge_vtx,
               None, # face_vtx_idx
               None, # face_vtx
               None, # cell_ln_to_gn
               None, # face_ln_to_gn
               edge_ln_to_gn,
               vtx_ln_to_gn,
               vtx_coord)

extrp.compute()


n_segment = extrp.n_entity_get(0, PDM._PDM_MESH_ENTITY_EDGE)

segment_vtx_idx, segment_vtx = extrp.connectivity_get(0, PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

segment_ln_to_gn = extrp.ln_to_gn_get(0, PDM._PDM_MESH_ENTITY_EDGE)

extrp_vtx_coord = extrp.vtx_coord_get(0)

segment_base_coord = np.zeros(3*n_segment)
segment_coord = np.zeros((2*n_segment, 3))
segment_connec = []
for i, edge_id in enumerate(select_edge):
  for j in range(2):
    for k in range(3):
      vtx_id = segment_vtx[2*i+j] - 1
      segment_base_coord[3*i+k] += 0.5*extrp_vtx_coord[3*vtx_id+k]
  segment_coord[2*i  ] = segment_base_coord[3*i:3*(i+1)]
  segment_coord[2*i+1] = 1.1*segment_coord[2*i]
  segment_connec.append([2*i, 2*i+1])

"""

# n_segment = len(select_edge)
# segment_base_coord = np.zeros(3*n_segment)
# segment_coord = np.zeros((2*n_segment, 3))
# segment_connec = []
# for i, edge_id in enumerate(select_edge):
#   segment_base_coord[3*i:3*(i+1)] = edge_mid_coord[3*edge_id:3*(edge_id+1)]
#   segment_coord[2*i  ] = segment_base_coord[3*i:3*(i+1)]
#   segment_coord[2*i+1] = 1.1*segment_coord[2*i]
#   segment_connec.append([2*i, 2*i+1])


# gen_gnum = PDM.GlobalNumbering(3, 1, 0, 1., comm)

# gen_gnum.gnum_set_from_parent(0, n_segment, edge_ln_to_gn[select_edge])

# gen_gnum.gnum_compute()

# segment_ln_to_gn = gen_gnum.gnum_get(0)["gnum"]

"""

# Randomly redistribute the segments?
# to do

# Locate the segments on the surface edges
clsp = PDM.ClosestPoints(comm, 1)

clsp.n_part_cloud_set(1, 1)

clsp.tgt_cloud_set(0, n_segment, segment_base_coord, segment_ln_to_gn)
clsp.src_cloud_set(0, n_edge,    edge_mid_coord,     edge_ln_to_gn)

clsp.compute()

# get ptp

# Exchange data between the segments and the edges


# Export for visu (ENSIGHT/VTK?)
_edge_vtx = np.reshape(edge_vtx-1, (n_edge, 2), order="C")
_vtx_coord = np.reshape(vtx_coord, (dims["n_vtx"], 3), order="C")
vtk_write_std_elements("sphere_edges_rank%d.vtk" % i_rank,
                       _vtx_coord,
                       ELT_TYPE_EDGE,
                       _edge_vtx,
                       {"gnum" : edge_ln_to_gn})

vtk_write_std_elements("sphere_segments_rank%d.vtk" % i_rank,
                       segment_coord,
                       ELT_TYPE_EDGE,
                       segment_connec)


connec = []
for j in range(n_face):
  start = face_edge_idx[j]
  end   = face_edge_idx[j+1]
  connec.append([k for k in face_vtx[start:end] - 1])

vtk_write_polydata("sphere_mesh_rank%d.vtk" % i_rank,
                   _vtx_coord,
                   connec)

# for i in range(n_part):
#   coord = []
#   for j in range(sphere_mesh["pn_vtx"][i]):
#     coord.append([x for x in sphere_mesh["pvtx_coord"][i][3*j:3*(j+1)]])

#   connec = []
#   for j in range(sphere_mesh["pn_face"][i]):
#     start = sphere_mesh["pface_vtx_idx"][i][j]
#     end   = sphere_mesh["pface_vtx_idx"][i][j+1]
#     connec.append([k for k in sphere_mesh["pface_vtx"][i][start:end] - 1])

#   vtk_write_polydata("sphere_mesh_rank%d_part%d.vtk" % (i_rank, i),
#                      coord,
#                      connec)


print("[{}] End".format(i_rank))

