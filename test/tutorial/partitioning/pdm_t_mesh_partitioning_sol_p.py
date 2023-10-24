#!/usr/bin/env python

def run_test():
  import numpy as np
  import argparse
  import mpi4py.MPI as MPI
  import Pypdm.Pypdm as PDM

  vtk = 0
  fe  = 0

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

  # Create partitioning object
  n_zone = 1 # fixed
  n_part = 1 # fixed
  i_part = 0 # fixed
  i_zone = 0 # fixed
  part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT;
  mpart = PDM.MultiPart(n_zone,
                        np.array([n_part]).astype(np.intc),
                        0,
                        part_method,
                        1,
                        np.ones(1).astype(np.double),
                        comm)

  renum_cell = bytes("PDM_PART_RENUM_CELL_NONE", 'ascii')
  renum_face = bytes("PDM_PART_RENUM_FACE_NONE", 'ascii')
  mpart.set_reordering(-1, # i_zone
                       renum_cell,
                       None,
                       renum_face)

  mpart.register_dmesh_nodal(i_zone, dmn)

  mpart.compute()

  # Get mesh arrrays in FE structure
  if fe:
    output = mpart.vtx_coord_get(i_zone,
                                 i_part)
    coords = output["np_vtx_coord"]

    pmn = mpart.part_mesh_nodal_get(i_zone)
    i_section = 0 # fixed
    output = PDM.part_mesh_nodal_get_sections(pmn,
                                              PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                                              i_part)
    elt_vtx      = output[i_section]["np_connec"]
    elt_ln_to_gn = output[i_section]["np_numabs"]
    n_elt        = len(elt_ln_to_gn)
    vtx_ln_to_gn = PDM.part_mesh_nodal_vtx_g_num_get(pmn, i_part)
    n_vtx        = len(vtx_ln_to_gn)

  # Get mesh arrrays in FV structure
  else :
    output = mpart.ln_to_gn_get(i_zone,
                                i_part,
                                PDM._PDM_MESH_ENTITY_VERTEX)

    vtx_ln_to_gn = output["np_entity_ln_to_gn"]
    n_vtx = len(vtx_ln_to_gn)

    output = mpart.vtx_coord_get(i_zone,
                                 i_part)

    coords = output["np_vtx_coord"]

    output = mpart.ln_to_gn_get(i_zone,
                                i_part,
                                PDM._PDM_MESH_ENTITY_EDGE)

    edge_ln_to_gn = output["np_entity_ln_to_gn"]
    n_edge = len(edge_ln_to_gn)

    output = mpart.connectivity_get(i_zone,
                                    i_part,
                                    PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    edge_vtx = output["np_entity1_entity2"]

    output = mpart.ln_to_gn_get(i_zone,
                                i_part,
                                PDM._PDM_MESH_ENTITY_FACE)

    face_ln_to_gn = output["np_entity_ln_to_gn"]
    n_face = len(face_ln_to_gn)

    output = mpart.connectivity_get(i_zone,
                                    i_part,
                                    PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)

    face_edge_idx = output["np_entity1_entity2_idx"]
    face_edge     = output["np_entity1_entity2"]

    output = mpart.ln_to_gn_get(i_zone,
                                i_part,
                                PDM._PDM_MESH_ENTITY_CELL)

    cell_ln_to_gn = output["np_entity_ln_to_gn"]
    n_cell = len(cell_ln_to_gn)

    output = mpart.connectivity_get(i_zone,
                                    i_part,
                                    PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    cell_face_idx = output["np_entity1_entity2_idx"]
    cell_face     = output["np_entity1_entity2"]

    # BONUS

    # step 1 : create
    extend_type = PDM._PDM_EXTEND_FROM_VTX
    depth       = 1
    part_ext = PDM.PartExtension(n_zone,
                                 np.array([n_part]).astype(np.intc),
                                 extend_type,
                                 depth,
                                 comm)

    # step 2 : set
    output = mpart.graph_comm_get(i_zone,
                                  i_part,
                                  PDM._PDM_BOUND_TYPE_VTX)

    vtx_part_bound_proc_idx = output["np_entity_part_bound_proc_idx"]
    vtx_part_bound_part_idx = output["np_entity_part_bound_part_idx"]
    vtx_part_bound          = output["np_entity_part_bound"]

    part_ext.set_part(i_zone,
                      i_part,
                      n_cell,
                      n_face,
                      0, # n_face_part_bound
                      0, # n_face_group
                      n_edge,
                      n_vtx,
                      cell_face_idx,
                      cell_face,
                      None, # face_cell
                      face_edge_idx,
                      face_edge,
                      None, # face_vtx_idx
                      None, # face_vtx
                      edge_vtx,
                      None, # face_group_idx
                      None, # face_group
                      None, # face_join_idx
                      None, # face_join
                      None, # face_part_bound_proc_idx
                      None, # face_part_bound_part_idx
                      None, # face_part_bound
                      vtx_part_bound_proc_idx,
                      vtx_part_bound_part_idx,
                      vtx_part_bound,
                      cell_ln_to_gn,
                      face_ln_to_gn,
                      edge_ln_to_gn,
                      vtx_ln_to_gn,
                      None, # face_group_ln_to_gn
                      coords)

    # step 3 : compute
    part_ext.compute()

    # step 4 : get
    # Cell
    cell_ext_ln_to_gn = part_ext.get_ln_to_gn(i_zone,
                                              i_part,
                                              PDM._PDM_MESH_ENTITY_CELL)

    cell_face_ext, cell_face_ext_idx = part_ext.get_connectivity(i_zone,
                                                                 i_part,
                                                                 PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    # Face
    face_ext_ln_to_gn = part_ext.get_ln_to_gn(i_zone,
                                              i_part,
                                              PDM._PDM_MESH_ENTITY_FACE)

    face_edge_ext, face_edge_ext_idx = part_ext.get_connectivity(i_zone,
                                                                 i_part,
                                                                 PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)

    # Edge
    edge_ext_ln_to_gn = part_ext.get_ln_to_gn(i_zone,
                                              i_part,
                                              PDM._PDM_MESH_ENTITY_EDGE)

    edge_vtx_ext, edge_vtx_ext_idx = part_ext.get_connectivity(i_zone,
                                                               i_part,
                                                               PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    # Vertices
    vtx_ext_ln_to_gn = part_ext.get_ln_to_gn(i_zone,
                                             i_part,
                                             PDM._PDM_MESH_ENTITY_VERTEX)

    vtx_coord_ext = part_ext.get_coord(i_zone,
                                       i_part)

    # step 5 : free (implicit in Python)

if __name__ == '__main__':
  run_test()
