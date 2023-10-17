#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM

import argparse

def visu_2d(comm,
            directory,
            name,
            vtx_coord,
            vtx_ln_to_gn,
            face_vtx_idx,
            face_vtx,
            face_ln_to_gn,
            vtx_fields=None):

  wrt = PDM.Writer("Ensight",
                   PDM._PDM_WRITER_FMT_ASCII,
                   PDM._PDM_WRITER_TOPO_VARIABLE,
                   PDM._PDM_WRITER_OFF,
                   directory,
                   name,
                   comm,
                   PDM._PDM_IO_KIND_MPI_SIMPLE,
                   1.,
                   "")

  n_part = len(vtx_coord)

  id_geom = wrt.geom_create(name,
                            n_part)

  id_var_rank = wrt.var_create(PDM._PDM_WRITER_OFF,
                               PDM._PDM_WRITER_VAR_SCALAR,
                               PDM._PDM_WRITER_VAR_ELEMENTS,
                               "i_part")

  id_var_elt_gnum = wrt.var_create(PDM._PDM_WRITER_OFF,
                                   PDM._PDM_WRITER_VAR_SCALAR,
                                   PDM._PDM_WRITER_VAR_ELEMENTS,
                                   "elt_gnum")

  id_var_vtx_fields = dict()
  if vtx_fields is not None:
    for name in vtx_fields:
      id_var_vtx_fields[name] = wrt.var_create(PDM._PDM_WRITER_OFF,
                                               PDM._PDM_WRITER_VAR_SCALAR,
                                               PDM._PDM_WRITER_VAR_VERTICES,
                                               name)

  wrt.step_beg(0.)

  for i_part in range(n_part):
    wrt.geom_coord_set(id_geom,
                       i_part,
                       vtx_coord   [i_part],
                       vtx_ln_to_gn[i_part])

    wrt.geom_faces_facevtx_add(id_geom,
                               i_part,
                               face_vtx_idx [i_part],
                               face_vtx     [i_part],
                               face_ln_to_gn[i_part])

  wrt.geom_write(id_geom)

  for i_part in range(n_part):
    wrt.var_set(id_var_rank,
                id_geom,
                i_part,
                (comm.rank*n_part + i_part)*np.ones(len(face_ln_to_gn[i_part]), dtype=np.double))

    wrt.var_set(id_var_elt_gnum,
                id_geom,
                i_part,
                face_ln_to_gn[i_part].astype(np.double))

    if vtx_fields is not None:
      for name in vtx_fields:
        wrt.var_set(id_var_vtx_fields[name],
                    id_geom,
                    i_part,
                    vtx_fields[name][i_part].astype(np.double))


  wrt.var_write(id_var_rank)
  wrt.var_write(id_var_elt_gnum)

  if vtx_fields is not None:
    for name in vtx_fields:
      wrt.var_write(id_var_vtx_fields[name])

  wrt.step_end()







def run_tuto():
  comm = MPI.COMM_WORLD



  parser = argparse.ArgumentParser()

  parser.add_argument("-n",          "--n_vtx_seg",  type=int, default=10)
  parser.add_argument("-src_n_part", "--src_n_part", type=int, default=1)
  parser.add_argument("-tgt_n_part", "--tgt_n_part", type=int, default=1)
  parser.add_argument("-nodal",      "--nodal",      action="store_true")

  args = parser.parse_args()

  n          = args.n_vtx_seg
  src_n_part = args.src_n_part
  tgt_n_part = args.tgt_n_part
  nodal      = args.nodal



  """
  Generate and partition the source mesh
  """
  src_mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                              elt_type    = PDM._PDM_MESH_NODAL_POLY_2D,
                                              xmin        = 0.,
                                              ymin        = 0.,
                                              zmin        = 0.,
                                              lengthx     = 1.,
                                              lengthy     = 1.,
                                              n_x         = n,
                                              n_y         = n,
                                              n_part      = src_n_part,
                                              part_method = PDM._PDM_SPLIT_DUAL_WITH_PARMETIS)

  src_n_vtx        = src_mesh["pn_vtx"]
  src_n_face       = src_mesh["pn_face"]
  src_vtx_coord    = src_mesh["pvtx_coord"]
  src_face_vtx_idx = src_mesh["pface_edge_idx"]
  src_face_vtx     = src_mesh["pface_vtx"]
  src_face_edge    = src_mesh["pface_edge"]
  src_edge_vtx     = src_mesh["pedge_vtx"]

  src_vtx_ln_to_gn  = src_mesh["pvtx_ln_to_gn"]
  src_face_ln_to_gn = src_mesh["pface_ln_to_gn"]

  """
  Randomize the coordinates just for fun
  """
  for i_part in range(src_n_part):
    for i_vtx in range(src_n_vtx[i_part]):
      # make sure the deformation is consistent for vertices on partition boundaries
      np.random.seed(src_vtx_ln_to_gn[i_part][i_vtx])
      for i in range(2):
        x = src_vtx_coord[i_part][3*i_vtx + i]
        if x > 1e-6 and x < 1-1e-6:
          src_vtx_coord[i_part][3*i_vtx + i] = x + 0.07*(2*np.random.rand() - 1)/(n - 1)


  visu_2d(comm,
          "mesh_location_sol_p",
          "src_mesh",
          src_vtx_coord,
          src_vtx_ln_to_gn,
          src_face_vtx_idx,
          src_face_vtx,
          src_face_ln_to_gn)



  """
  Then we need to generate and partition a target point cloud
  For nicer visualization we will generate a second mesh, and use
  its nodes as the target point cloud.
  """
  tgt_mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                              elt_type    = PDM._PDM_MESH_NODAL_QUAD4,
                                              xmin        = 0.3,
                                              ymin        = 0.3,
                                              zmin        = 0,
                                              lengthx     = 1.,
                                              lengthy     = 1.,
                                              n_x         = n,
                                              n_y         = n,
                                              n_part      = tgt_n_part,
                                              part_method = PDM._PDM_SPLIT_DUAL_WITH_PARMETIS)

  tgt_n_vtx        = tgt_mesh["pn_vtx"]
  tgt_n_face       = tgt_mesh["pn_face"]
  tgt_vtx_coord    = tgt_mesh["pvtx_coord"]
  tgt_face_vtx_idx = tgt_mesh["pface_edge_idx"]
  tgt_face_vtx     = tgt_mesh["pface_vtx"]

  tgt_vtx_ln_to_gn  = tgt_mesh["pvtx_ln_to_gn"]
  tgt_face_ln_to_gn = tgt_mesh["pface_ln_to_gn"]



  """
  Create the MeshLocation object
  """
  mesh_loc = PDM.MeshLocation(1,
                              comm)


  """
  Set the target point cloud
  """
  mesh_loc.n_part_cloud_set(0,
                            tgt_n_part)

  for i_part in range(tgt_n_part):
    mesh_loc.cloud_set(0,
                       i_part,
                       tgt_vtx_coord   [i_part],
                       tgt_vtx_ln_to_gn[i_part])

  """
  Set the source mesh
  Here you have essentially two options :
   - you can either define the mesh with *nodal* connectivity (i.e. Finite-Element style)
   - or with "descending" connectivity (i.e. Finite-Volume style)
  """
  mesh_loc.mesh_n_part_set(src_n_part)

  for i_part in range(src_n_part):
    if nodal:
      mesh_loc.nodal_part_set_2d(i_part,
                                 src_face_vtx_idx [i_part],
                                 src_face_vtx     [i_part],
                                 src_face_ln_to_gn[i_part],
                                 src_vtx_coord    [i_part],
                                 src_vtx_ln_to_gn [i_part])
    else:
      mesh_loc.part_set_2d(i_part,
                           src_face_vtx_idx [i_part],
                           src_face_edge    [i_part],
                           src_face_ln_to_gn[i_part],
                           src_edge_vtx     [i_part],
                           src_vtx_coord    [i_part],
                           src_vtx_ln_to_gn [i_part])

  """
  Set the geometric tolerance (optional)
  """
  mesh_loc.tolerance_set(1e-6)

  """
  Set the location preconditioning method (optional)
  """
  mesh_loc.method_set(0)

  """
  Compute location
  """
  mesh_loc.compute()

  """
  Dump elapsed and CPU times
  """
  mesh_loc.dump_times()


  """
  Get the list of (un)located target points
  """
  for i_part in range(tgt_n_part):
    located_tgt   = mesh_loc.located_get(0, i_part)
    unlocated_tgt = mesh_loc.located_get(0, i_part)

    print(f"[{comm.rank}] part {i_part}, located_tgt : {located_tgt}", flush=True)

  """
  We can also access the other location information from the target point-of-view...
  (in which source mesh element is located each target point, if any)
  """


  """
  ...as well as from the source point-of-view (which targets points are located inside each source mesh element)
  """
  # src_result = mesh_loc.points_in_elt_get(0, 0)

  # src_to_tgt_idx = src_result["elt_pts_inside_idx"]
  # n_pts = src_to_tgt_idx[src_n_face]

  """
  Now that we have located the target points in the source mesh,
  we can exchange data between the two.
  To complete this exercise, we will interpolate two fields from
  the source mesh to the target cloud.
  The first field is cell-based : we can simply use the face global ids for such a field,
  and check it matches the location data.
  The second field is node-based : we can use the node coordinates.
  """

  """
  First, compute the spatially interpolated fields on the source side.
  For the first field, the interpolation is straightforward : the target value is simply the same as the host source.
  The second field interpolation is trickier as you will need the cell->vertex connectivity built during the location computation to link the interpolation weights to the appropriate source nodes.
  """
  src_send_field1 = []
  src_send_field2 = []

  for i_part in range(src_n_part):
    src_result = mesh_loc.points_in_elt_get(0, i_part)
    src_to_tgt_idx = src_result["elt_pts_inside_idx"]
    n_pts = src_to_tgt_idx[src_n_face[i_part]]


    # Interpolate first field (cell-based)
    field1 = np.zeros(n_pts, dtype=np.double)
    for i_elt in range(src_n_face[i_part]):
      for i_pt in range(src_to_tgt_idx[i_elt], src_to_tgt_idx[i_elt+1]):
        field1[i_pt] = src_face_ln_to_gn[i_part][i_elt]
    src_send_field1.append(field1)



    # Interpolate second field (node-based)
    src_connect = mesh_loc.cell_vertex_get(i_part)
    src_cell_vtx_idx = src_connect["cell_vtx_idx"]
    src_cell_vtx     = src_connect["cell_vtx"]

    weights_idx = src_result["points_weights_idx"]
    weights     = src_result["points_weights"]

    field2 = np.zeros(n_pts, dtype=np.double)
    for i_elt in range(src_n_face[i_part]):
      for i_pt in range(src_to_tgt_idx[i_elt], src_to_tgt_idx[i_elt+1]):
        field2[i_pt] = 0

        elt_n_vtx = src_cell_vtx_idx[i_elt+1] - src_cell_vtx_idx[i_elt]
        assert(weights_idx[i_pt+1] - weights_idx[i_pt] == elt_n_vtx)
        for i_vtx in range(elt_n_vtx):
          vtx_id = src_cell_vtx[src_cell_vtx_idx[i_elt] + i_vtx] - 1
          field2[i_pt] += src_vtx_coord[i_part][3*vtx_id] * weights[weights_idx[i_pt] + i_vtx]

    src_send_field2.append(field2)



  """
  Second, use the PartToPart object to exchange the interpolated fields from the source mesh
  to the target cloud.
  This ParToPart object was built when computing the location and can be accessed from the MeshLocation object
  """

  # Get PartToPart object
  ptp = mesh_loc.part_to_part_get(0)

  # Initiate exchange of first field
  src_stride = 1
  request1 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                       PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                       src_send_field1,
                       part1_stride=src_stride,
                       interlaced_str=True)

  # Initiate exchange of second field
  request2 = ptp.iexch(PDM._PDM_MPI_COMM_KIND_P2P,
                       PDM._PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2,
                       src_send_field2,
                       part1_stride=src_stride,
                       interlaced_str=True)

  # Finalize both exchanges
  tgt_stride, tgt_recv_field1 = ptp.wait(request1)
  tgt_stride, tgt_recv_field2 = ptp.wait(request2)


  """
  Finally, visualize the interpolated target fields.
  (Beware of unlocated points!)
  """
  pis_located = []
  ptgt_field1 = []
  ptgt_field2 = []
  for i_part in range(tgt_n_part):
    located_tgt = mesh_loc.located_get(0, i_part)

    is_located = np.zeros(tgt_n_vtx[i_part], dtype=bool)
    tgt_field1 = np.zeros(tgt_n_vtx[i_part], dtype=np.double)
    tgt_field2 = np.zeros(tgt_n_vtx[i_part], dtype=np.double)
    for i, i_vtx in enumerate(located_tgt):
      is_located[i_vtx-1] = True
      tgt_field1[i_vtx-1] = tgt_recv_field1[i_part][i]
      tgt_field2[i_vtx-1] = tgt_recv_field2[i_part][i]

      error = abs(tgt_recv_field2[i_part][i] - tgt_vtx_coord[i_part][3*(i_vtx-1)])
      if error > 1e-9:
        print(f"!! error vtx {tgt_vtx_ln_to_gn[i_part][i_vtx]} : {error}")

    pis_located.append(is_located)
    ptgt_field1.append(tgt_field1)
    ptgt_field2.append(tgt_field2)

  visu_2d(comm,
       "mesh_location_sol_p",
       "tgt_mesh",
       tgt_vtx_coord,
       tgt_vtx_ln_to_gn,
       tgt_face_vtx_idx,
       tgt_face_vtx,
       tgt_face_ln_to_gn,
       {"is_located" : pis_located,
       "field1"      : ptgt_field1,
       "field2"      : ptgt_field2})



  """
  Final bonus : interpolate and exchange from target to source (if it makes any sense?)
  """


if __name__ == '__main__':
  run_tuto()
