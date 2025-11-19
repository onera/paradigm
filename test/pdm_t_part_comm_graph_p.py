#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

comm = MPI.COMM_WORLD

i_rank = comm.rank
n_rank = comm.size

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


def test_constructor_with_nuplet():
  pcg = PDM.PartCommGraph.nuplet_create(comm, entity_graph, entity_nuplet, 0)
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

test_constructor()
test_constructor_with_nuplet()
test_with_part_mesh_nodal()



