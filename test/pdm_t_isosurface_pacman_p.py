#!/usr/bin/env python
import numpy as np
import time

def sdf_sphere(x, y, z, cx, cy, cz, r):
  """
  Sphere signed distance function
  """
  return np.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2) - r

def sdf_pacman(coord, time):
  """
  Pacman signed distance function
  """
  x = coord[0::3]
  y = coord[1::3]
  z = coord[2::3]

  t = np.mod(time, 5.)

  # head and eyes
  head = sdf_sphere(x, y-t, z,  0   , 0   , 0   , 0.48)
  leye = sdf_sphere(x, y-t, z, -0.18, 0.32, 0.18, 0.11)
  reye = sdf_sphere(x, y-t, z,  0.18, 0.32, 0.18, 0.11)

  # mouth
  yi =  0.258
  zi = -0.128

  gamma = 2.03173778224
  delta = 0.5

  T = 0.1 + 0.45*(np.cos(np.pi*(2.*t+1.)) + 1.)
  A1 = np.cos(gamma + delta*T)
  B1 = -np.sin(gamma + delta*T)
  C1 = -(A1*yi + B1*zi)

  A2 = -np.cos(gamma - T)
  B2 =  np.sin(gamma - T)
  C2 = -(A2*yi + B2*zi)

  fm1 = A1*(y-t) + B1*z + C1
  fm2 = A2*(y-t) + B2*z + C2
  mouth = np.maximum(fm1, fm2)

  # subtract mouth wedge from head
  pacman = np.maximum(head, -mouth)

  # subtract eyes from head
  eye    = np.minimum(leye, reye)
  pacman = np.maximum(pacman, -eye)

  # gums
  rgum = 0.06
  zgum = zi
  iy   = np.array(y + 0.5).astype(int)

  # rgum1 = np.maximum(0, np.minimum(0.125*(t - iy - 2.), rgum)) # gums reappear after a delay
  rgum1 = np.zeros(iy.shape)
  rgum1[np.where(iy > t - 1.1*rgum)] = rgum

  gum = sdf_sphere(x, y-iy, z, 0, 0, zgum, rgum1)

  # merge gums
  pacman = np.minimum(pacman, gum)

  return pacman


def run(n_vtx_seg, elt_type, n_step, visu, local, part_method):
  import mpi4py.MPI as MPI
  import Pypdm.Pypdm as PDM

  comm = MPI.COMM_WORLD

  # Generate volume mesh
  length_y = 5
  n_part   = 1

  if comm.rank == 0:
    print("Generate volume mesh")

  comm.Barrier()
  t_start = time.time()

  mesh = PDM.generate_mesh_parallelepiped_ngon(comm,
                                               elt_type,
                                               1,
                                               "",
                                               -0.5,
                                               -0.5,
                                               -0.5,
                                               1,
                                               length_y,
                                               1,
                                               n_vtx_seg,
                                               int(length_y*n_vtx_seg),
                                               n_vtx_seg,
                                               n_part,
                                               part_method)

  comm.Barrier()
  t_end = time.time()

  if comm.rank == 0:
    delta_t = t_end - t_start
    print("Took %.3fs" % (delta_t))


  # Create isosurface object
  isos = PDM.Isosurface(3, comm)

  id_iso_field = isos.add(PDM.Isosurface.FIELD, [0.])

  if not local:
    isos.redistribution_set(PDM.Isosurface.REEQUILIBRATE,
                            part_method)

  # Set mesh
  isos.n_part_set(n_part)

  for i_part in range(n_part):
    # Connectivities
    isos.connectivity_set(i_part,
                          PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE,
                          mesh["pcell_face_idx"][i_part],
                          mesh["pcell_face"    ][i_part])

    isos.connectivity_set(i_part,
                          PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX,
                          mesh["pface_edge_idx"][i_part],
                          mesh["pface_vtx"     ][i_part])

    # Coordinates
    isos.coordinates_set(i_part,
                         mesh["pvtx_coord"][i_part])

    # Global IDs
    isos.ln_to_gn_set(i_part,
                      PDM._PDM_MESH_ENTITY_CELL,
                      mesh["pcell_ln_to_gn"][i_part])

    isos.ln_to_gn_set(i_part,
                      PDM._PDM_MESH_ENTITY_FACE,
                      mesh["pface_ln_to_gn"][i_part])

    isos.ln_to_gn_set(i_part,
                      PDM._PDM_MESH_ENTITY_VTX,
                      mesh["pvtx_ln_to_gn"][i_part])

    # Groups
    isos.group_set(i_part,
                   PDM._PDM_MESH_ENTITY_FACE,
                   mesh["pn_surface"            ][i_part],
                   mesh["psurface_face_idx"     ][i_part],
                   mesh["psurface_face"         ][i_part],
                   mesh["psurface_face_ln_to_gn"][i_part])


  # Prepare writer
  if visu:
    writer = PDM.Writer("Ensight",
                        PDM._PDM_WRITER_FMT_BIN,
                        PDM._PDM_WRITER_TOPO_VARIABLE,
                        PDM._PDM_WRITER_OFF,
                        "isosurface_pacman_p",
                        "pacman",
                        comm,
                        PDM._PDM_IO_KIND_MPI_SIMPLE,
                        1.,
                        "")

    id_geom = writer.geom_create("pacman_surf", n_part)

  # Time loop
  comm.Barrier()
  t_start = time.time()

  t = 0.

  for i_step in range(n_step):
    if comm.rank == 0:
      print(f"Step {i_step}/{n_step}")

    # Update field values
    pfield = []
    for i_part in range(n_part):
      field = sdf_pacman(mesh["pvtx_coord"][i_part], t)
      pfield.append(field)
      isos.field_set(id_iso_field, i_part, field)

    # Compute isosurface
    isos.compute(id_iso_field)

    # Export isosurface
    if visu:
      writer.step_beg(t)
      for i_part in range(n_part):
        iso_vtx_coord    = isos.coordinates_get(id_iso_field, i_part)
        iso_vtx_ln_to_gn = isos.ln_to_gn_get(id_iso_field, i_part, PDM._PDM_MESH_ENTITY_VTX)

        iso_face_vtx_idx, iso_face_vtx = isos.connectivity_get(id_iso_field, i_part, PDM._PDM_CONNECTIVITY_TYPE_FACE_VTX)
        iso_face_ln_to_gn = isos.ln_to_gn_get(id_iso_field, i_part, PDM._PDM_MESH_ENTITY_FACE)

        writer.geom_coord_set(id_geom,
                              i_part,
                              iso_vtx_coord,
                              iso_vtx_ln_to_gn)

        writer.geom_faces_facevtx_add(id_geom,
                                      i_part,
                                      iso_face_vtx_idx,
                                      iso_face_vtx,
                                      iso_face_ln_to_gn)

      writer.geom_write(id_geom)
      writer.geom_data_reset(id_geom)
      writer.step_end()

    isos.reset(id_iso_field)

    t += 0.05

  isos.dump_times()

  comm.Barrier()
  t_end = time.time()

  if comm.rank == 0:
    delta_t = t_end - t_start
    print("Total elapsed time [s] = %.3f (%.3f / step)" % (delta_t, delta_t/n_step))

    print("The End :D")



if __name__ == '__main__':
  """
  Parse command line arguments
  """
  import argparse
  parser = argparse.ArgumentParser()

  parser.add_argument("-n",
                      "--n_vtx_seg",
                      help="Number of vertices along side of domain",
                      type=int,
                      default=10)

  parser.add_argument("-elt_type",
                      "--elt_type",
                      help="Volume mesh element type",
                      type=int,
                      default=8)

  parser.add_argument("-n_step",
                      "--n_step",
                      help="Number of time steps",
                      type=int,
                      default=10)

  parser.add_argument("-v",
                      "--visu",
                      help="Export for visualization",
                      action="store_true")

  parser.add_argument("-local",
                      "--local",
                      help="Do not redistribute the isosurface mesh",
                      action="store_true")

  parser.add_argument("-part_method",
                      "--part_method",
                      help="Partitioning method (1: ParMetis, 2: PT-Scotch, 3: Hilbert, 4: Implicit)",
                      type=int,
                      default=3)

  # TODO: wrap nodal mesh generation
  # parser.add_argument("-nodal",
  #                     "--nodal",
  #                     help="Use nodal mesh",
  #                     action="store_true")

  args = parser.parse_args()

  run(n_vtx_seg   = args.n_vtx_seg,
      elt_type    = args.elt_type,
      n_step      = args.n_step,
      visu        = args.visu,
      local       = args.local,
      part_method = args.part_method)
