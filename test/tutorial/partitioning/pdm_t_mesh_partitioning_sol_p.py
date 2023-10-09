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

  del dcube # facultatif

  # Create partitioning object
  n_zone      = 1; # fixed
  n_part      = 1; # fixed
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
  mpart.multipart_set_reordering(-1, # i_zone
                                 renum_cell,
                                 renum_face,
                                 None)

  mpart.multipart_register_dmesh_nodal(0, dmn)

  mpart.multipart_run_ppart()

  i_part = 0 # fixed
  i_zone = 0 # fixed

  # Get mesh arrrays in FE structure
  if fe:
    pmn = mpart.multipart_part_mesh_nodal_get(i_zone)

    output = PDM.part_mesh_nodal_get_sections(pmn,
                                              PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                                              i_part)

    elt_vtx      = output["np_connec"]
    elt_ln_to_gn = output["np_numabs"]
    # TO DO : get n_elt, n_vtx, coords, vtx_ln_to_gn adding it to pdm_part_mesh_nodal.pxi

  # Get mesh arrrays in FV structure
  else :
    output = mpart.multipart_ln_to_gn_get(i_part,
                                          i_zone,
                                          PDM._PDM_MESH_ENTITY_VERTEX)

    vtx_ln_to_gn = output["np_entity_ln_to_gn"]
    n_vtx = len(vtx_ln_to_gn)

    output = mpart.multipart_vtx_coord_get(i_part,
                                           i_zone)

    coords = output["np_vtx_coord"]

    output = mpart.multipart_ln_to_gn_get(i_part,
                                          i_zone,
                                          PDM._PDM_MESH_ENTITY_EDGE)

    edge_ln_to_gn = output["np_entity_ln_to_gn"]
    n_edge = len(edge_ln_to_gn)

    output = mpart.multipart_connectivity_get(i_part,
                                              i_zone,
                                              PDM._PDM_CONNECTIVITY_TYPE_EDGE_VTX)

    edge_vtx = output["np_entity1_entity2"]

    output = mpart.multipart_ln_to_gn_get(i_part,
                                          i_zone,
                                          PDM._PDM_MESH_ENTITY_FACE)

    face_ln_to_gn = output["np_entity_ln_to_gn"]
    n_face = len(face_ln_to_gn)

    output = mpart.multipart_connectivity_get(i_part,
                                              i_zone,
                                              PDM._PDM_CONNECTIVITY_TYPE_FACE_EDGE)

    face_edge_idx = output["np_entity1_entity2_idx"]
    face_edge     = output["np_entity1_entity2"]

    output = mpart.multipart_ln_to_gn_get(i_part,
                                          i_zone,
                                          PDM._PDM_MESH_ENTITY_CELL)

    cell_ln_to_gn = output["np_entity_ln_to_gn"]
    n_cell = len(cell_ln_to_gn)

    output = mpart.multipart_connectivity_get(i_part,
                                              i_zone,
                                              PDM._PDM_CONNECTIVITY_TYPE_CELL_FACE)

    cell_face_idx = output["np_entity1_entity2_idx"]
    cell_face     = output["np_entity1_entity2"]

if __name__ == '__main__':
  run_test()
