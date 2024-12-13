#!/usr/bin/env python

# Imports
import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM
from math import *

def run_test():

  # MPI
  comm = MPI.COMM_WORLD
  i_rank = comm.rank
  n_rank = comm.size

  # Generate mesh
  n_part = 2
  mesh = PDM.generate_mesh_parallelepiped_ngon(comm,
                                               PDM._PDM_MESH_NODAL_TETRA4,
                                               1,
                                               "",
                                               0.,
                                               0.,
                                               0.,
                                               1.,
                                               1.,
                                               1.,
                                               10,
                                               10,
                                               10,
                                               n_part,
                                               PDM._PDM_SPLIT_DUAL_WITH_HILBERT)

  reorient = PDM.part_mesh_reorient_geom(3,
                                         np.array([0., 0., 0.]),
                                         mesh["pn_cell"],
                                         mesh["pn_face"],
                                         mesh["pn_surface"],
                                         mesh["pcell_face_idx"],
                                         mesh["pcell_face"],
                                         None,
                                         mesh["pface_edge_idx"],
                                         mesh["pface_vtx"],
                                         mesh["pvtx_coord"],
                                         mesh["psurface_face_idx"],
                                         mesh["psurface_face"],
                                         comm);
  assert reorient == 0

  if i_rank == 0:
    print("Ok 3D 1 :D", flush=True)

  mesh["pface_cell"] = []
  for i_part in range(n_part):
    mesh["pface_cell"].append(np.zeros(2*mesh["pn_face"][i_part], dtype=np.int32))
    for i_cell in range(mesh["pn_cell"][i_part]):
      for idx_face in range(mesh["pcell_face_idx"][i_part][i_cell], mesh["pcell_face_idx"][i_part][i_cell+1]):
        i_face = np.abs(mesh["pcell_face"][i_part][idx_face])-1
        if (mesh["pcell_face"][i_part][idx_face] > 0):
          mesh["pface_cell"][i_part][2*i_face+1] = i_cell+1;
        if (mesh["pcell_face"][i_part][idx_face] < 0):
          mesh["pface_cell"][i_part][2*i_face  ] = i_cell+1;

  reorient = PDM.part_mesh_reorient_geom(3,
                                         np.array([0., 0., 0.]),
                                         mesh["pn_cell"],
                                         mesh["pn_face"],
                                         mesh["pn_surface"],
                                         mesh["pcell_face_idx"],
                                         mesh["pcell_face"],
                                         mesh["pface_cell"],
                                         mesh["pface_edge_idx"],
                                         mesh["pface_vtx"],
                                         mesh["pvtx_coord"],
                                         mesh["psurface_face_idx"],
                                         mesh["psurface_face"],
                                         comm)
  assert reorient == 1

  if i_rank == 0:
    print("Ok 3D 2 :D", flush=True)

  reorient = PDM.part_mesh_reorient_geom(3,
                                         np.array([0., 0., 0.]),
                                         mesh["pn_cell"],
                                         mesh["pn_face"],
                                         mesh["pn_surface"],
                                         None,
                                         None,
                                         mesh["pface_cell"],
                                         mesh["pface_edge_idx"],
                                         mesh["pface_vtx"],
                                         mesh["pvtx_coord"],
                                         mesh["psurface_face_idx"],
                                         mesh["psurface_face"],
                                         comm);
  assert reorient == 0

  if i_rank == 0:
    print("Ok 3D 3 :D", flush=True)

  # Check
  if i_rank == 0:
    print("End :)", flush=True)


if __name__ == '__main__':
  run_test()
