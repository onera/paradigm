#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

assert n_rank == 2, "This test is supposed to run on exactly 2 MPI ranks"

if i_rank == 0:
  entity_graph = [np.array([3, 1, 1, 1,
                            6, 1, 1, 5,
                            9, 1, 1, 9], dtype=np.int32)]
  entity_nuplet = [np.array([1, 1, 1], dtype=np.int32)]
  send_data   = [np.array([1, 2, 3], dtype=np.int32)]

  send_stri_v = [np.array([0, 1, 0], dtype=np.int32)]
  send_data_v = [np.array([1], dtype=np.int32)]
else:
  entity_graph = [np.array([1, 0, 1, 3,
                            5, 0, 1, 6,
                            9, 0, 1, 9], dtype=np.int32)]
  entity_nuplet = [np.array([1, 1, 1], dtype=np.int32)]

  send_data   = [np.array([10, 20, 30], dtype=np.int32)]

  send_stri_v = [np.array([2, 0, 0], dtype=np.int32)]
  send_data_v = [np.array([20, 10], dtype=np.int32)]

def test_constructor():
  pcg = PDM.PartCommGraph(comm, entity_graph)
  owner = pcg.owner_get(0)

  if i_rank == 0:
    assert (owner == np.array([1, 1, 1])).all()
  else:
    assert (owner == np.array([0, 0, 0])).all()

  # > Constant stride
  _, recv_data = pcg.exch(send_data)

  if i_rank == 0:
    assert (recv_data[0] == np.array([10, 20, 30])).all()
  else:
    assert (recv_data[0] == np.array([1, 2, 3])).all()

  # > Variable stride
  recv_stri, recv_data = pcg.exch(send_data_v, send_stri_v)

  if i_rank == 0:
    assert (recv_stri[0] == np.array([2, 0, 0])).all()
    assert (recv_data[0] == np.array([20, 10])).all()
  else:
    assert (recv_stri[0] == np.array([0, 1, 0])).all()
    assert (recv_data[0] == np.array([1])).all()

  # > Allreduce
  if i_rank == 0:
    pdata = [np.ones(10, dtype=np.int32)]
  else:
    pdata = [np.ones(10, dtype=np.int32)]

  pcg.allreduce(1, MPI.SUM, False, pdata)
  if i_rank == 0:
    assert( pdata == np.array([1, 1, 2, 1, 1, 2, 1, 1, 2, 1])).all()
  else:
    assert( pdata == np.array([2, 1, 1, 1, 2, 1, 1, 1, 2, 1])).all()

def test_constructor_with_nuplet():
  pcg = PDM.PartCommGraph(comm, entity_graph, entity_nuplet, False)
  owner = pcg.owner_get(0)

  if i_rank == 0:
    assert (owner == np.array([1, 1, 1])).all()
  else:
    assert (owner == np.array([0, 0, 0])).all()

  print(owner)

def test_with_part_mesh_nodal():
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

  dmn.generate_distribution()

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
  mpart.dmesh_nodal_set(i_domain, dmn)
  mpart.compute()

  pmn = mpart.part_mesh_nodal_get(i_domain)

  pcg_vtx = pmn.part_comm_graph_get(PDM._PDM_MESH_ENTITY_VTX)

  entity_graph = pcg_vtx.entity_graph_get(i_part)

  print(entity_graph)
  print("iloc i_copy_rank  i_copy_rank_part   i_copy_rank_loc")
  for i in range(int(entity_graph.size//4)):
    print(entity_graph[4*i],"            ", entity_graph[4*i+1],"              ",entity_graph[4*i+2],"               ", entity_graph[4*i+3])


def test_entity1_to_entity2():

  if i_rank == 0:
    pvtx_graph = [
      np.array([
        3, 1, 1, 1,
        4, 0, 2, 1,
        5, 0, 2, 2,
        6, 1, 1, 4,
        7, 0, 2, 4,
        8, 1, 1, 7,
        8, 1, 2, 1,
        9, 0, 2, 6,
        10, 0, 2, 7,
        10, 1, 2, 4
      ], dtype=np.int32),
      np.array([
        1, 0, 1, 4,
        2, 0, 1, 5,
        4, 0, 1, 7,
        6, 0, 1, 9,
        7, 0, 1, 10,
        7, 1, 2, 4,
        10, 1, 2, 7
      ], dtype=np.int32),
    ]

    pn_vtx = [10, 10]

    pedge_vtx = [
      np.array([
        1, 2,
        2, 3,
        1, 4,
        2, 5,
        3, 6,
        4, 5,
        5, 6,
        5, 7,
        6, 8,
        7, 8,
        7, 9,
        8, 10,
        9, 10
      ], dtype=np.int32),
      np.array([
        1, 2,
        1, 3,
        2, 4,
        3, 4,
        3, 5,
        4, 6,
        5, 6,
        6, 7,
        5, 8,
        6, 9,
        7, 10,
        8, 9,
        9, 10
      ], dtype=np.int32)
    ]
  else:
    pvtx_graph = [
      np.array([
        1, 0, 1, 3,
        4, 0, 1, 6,
        7, 0, 1, 8,
        7, 1, 2, 1,
        8, 1, 2, 2,
        9, 1, 2, 3
      ], dtype=np.int32),
      np.array([
        1, 0, 1, 8,
        1, 1, 1, 7,
        2, 1, 1, 8,
        3, 1, 1, 9,
        4, 0, 2, 7,
        4, 0, 1, 10,
        7, 0, 2, 10
      ], dtype=np.int32),
    ]

    pn_vtx = [9, 9]

    pedge_vtx = [
      np.array([
        1, 2,
        2, 3,
        1, 4,
        2, 5,
        3, 6,
        4, 5,
        5, 6,
        4, 7,
        5, 8,
        6, 9,
        7, 8,
        8, 9
      ], dtype=np.int32),
      np.array([
        1, 2,
        2, 3,
        1, 4,
        2, 5,
        3, 6,
        4, 5,
        5, 6,
        4, 7,
        5, 8,
        6, 9,
        7, 8,
        8, 9
      ], dtype=np.int32)
    ]


  n_part = 2
  pedge_vtx_idx = [
    2*np.arange(ev.size//2 + 1, dtype=np.int32) for ev in pedge_vtx
  ]

  pcg_vtx = PDM.PartCommGraph(comm,
                              pvtx_graph)

  pcg_edge = pcg_vtx.entity1_to_entity2(pn_vtx,
                                        pedge_vtx_idx,
                                        pedge_vtx)
  # pcg_edge = PDM.pcg_entity1_to_entity2(pcg_vtx,
  #                                       pn_vtx,
  #                                       pedge_vtx_idx,
  #                                       pedge_vtx)

  pedge_graph = [pcg_edge.entity_graph_get(i_part) for i_part in range(n_part)]

  if i_rank == 0:
    assert (pedge_graph[0] == np.array([ 6, 0, 2, 1,
                                         8, 0, 2, 3,
                                        11, 0, 2, 6,
                                        13, 0, 2, 8,
                                         5, 1, 1, 3,
                                         9, 1, 1, 8,
                                        12, 1, 2, 3])).all()
    assert (pedge_graph[1] == np.array([ 1, 0, 1,  6,
                                         3, 0, 1,  8,
                                         6, 0, 1, 11,
                                         8, 0, 1, 13,
                                        11, 1, 2,  8])).all()
  else:
    assert (pedge_graph[0] == np.array([ 3, 0, 1, 5,
                                         8, 0, 1, 9,
                                        11, 1, 2, 1,
                                        12, 1, 2, 2])).all()
    assert (pedge_graph[1] == np.array([3, 0, 1, 12,
                                        8, 0, 2, 11,
                                        1, 1, 1, 11,
                                        2, 1, 1, 12])).all()


test_constructor()
test_constructor_with_nuplet()
test_with_part_mesh_nodal()
test_entity1_to_entity2()

if i_rank == 0:
  print("End :)")