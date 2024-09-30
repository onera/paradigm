#!/usr/bin/env python
import numpy as np


def my_field_function(x, y, z):
  """
  User-defined field function
  """
  import numpy as np
  return x + y + z + 0.2*np.cos(3*y)


def run(visu):
  import mpi4py.MPI as MPI
  import Pypdm.Pypdm as PDM

  comm = MPI.COMM_WORLD

  # Generate volume mesh
  mesh = PDM.generate_mesh_parallelepiped_ngon(comm,
                                               PDM._PDM_MESH_NODAL_HEXA8,
                                               1,
                                               "",
                                               -0.5,
                                               -0.5,
                                               -0.5,
                                               1.,
                                               1.,
                                               1.,
                                               10,
                                               10,
                                               10,
                                               1,
                                               PDM._PDM_SPLIT_DUAL_WITH_HILBERT)

  # Create isosurface object
  isos = PDM.Isosurface(3, comm)

  # Add isosurface of user-defined field function
  id_iso = isos.add(PDM.Isosurface.FUNCTION, [-0.3, 0.0, 0.3])

  # Set field function
  isos.field_function_set(id_iso, my_field_function)

  isos.redistribution_set(PDM.Isosurface.REEQUILIBRATE,
                          PDM._PDM_SPLIT_DUAL_WITH_HILBERT)

  # Set mesh
  isos.n_part_set(1)

  # Connectivities
  isos.connectivity_set(0,
                        PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE,
                        mesh["pcell_face_idx"][0],
                        mesh["pcell_face"    ][0])

  isos.connectivity_set(0,
                        PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX,
                        mesh["pface_edge_idx"][0],
                        mesh["pface_vtx"     ][0])

  # Coordinates
  isos.coordinates_set(0,
                       mesh["pvtx_coord"][0])

  # Global IDs
  isos.ln_to_gn_set(0,
                    PDM._PDM_MESH_ENTITY_CELL,
                    mesh["pcell_ln_to_gn"][0])

  isos.ln_to_gn_set(0,
                    PDM._PDM_MESH_ENTITY_FACE,
                    mesh["pface_ln_to_gn"][0])

  isos.ln_to_gn_set(0,
                    PDM._PDM_MESH_ENTITY_VTX,
                    mesh["pvtx_ln_to_gn"][0])

  # Groups
  isos.group_set(0,
                 PDM._PDM_MESH_ENTITY_FACE,
                 mesh["psurface_face_idx"     ][0],
                 mesh["psurface_face"         ][0],
                 mesh["psurface_face_ln_to_gn"][0])

  # Compute isosurface
  isos.compute(id_iso)
  isos.dump_times()

  # Export isosurface
  if visu:
    iso_vtx_coord    = isos.coordinates_get(id_iso, 0)
    iso_vtx_ln_to_gn = isos.ln_to_gn_get(id_iso, 0, PDM._PDM_MESH_ENTITY_VTX)

    iso_face_vtx_idx, iso_face_vtx = isos.connectivity_get(id_iso, 0, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)
    iso_face_ln_to_gn = isos.ln_to_gn_get(id_iso, 0, PDM._PDM_MESH_ENTITY_FACE)

    PDM.writer_wrapper(comm,
                       "isosurface_user_function_p",
                       "isosurface",
                       [iso_vtx_coord],
                       [iso_vtx_ln_to_gn],
                       [iso_face_vtx_idx],
                       [iso_face_vtx],
                       [iso_face_ln_to_gn])

  if comm.rank == 0:
    print("The End :D")


if __name__ == '__main__':
  """
  Parse command line arguments
  """
  import argparse
  parser = argparse.ArgumentParser()

  parser.add_argument("-v",
                      "--visu",
                      help="Export for visualization",
                      action="store_true")

  args = parser.parse_args()

  run(visu = args.visu)
