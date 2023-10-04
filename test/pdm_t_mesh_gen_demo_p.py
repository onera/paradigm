#!/usr/bin/env python

import mpi4py.MPI as MPI
import numpy as np
import Pypdm.Pypdm as PDM


import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-n",      "--n_vtx_seg", default=3)
parser.add_argument("-n_part", "--n_part",    default=1)
parser.add_argument("-v",      "--visu",      action="store_true")

args = parser.parse_args()

n      = int(args.n_vtx_seg)
n_part = int(args.n_part)
visu   = args.visu

comm = MPI.COMM_WORLD

i_rank = comm.rank


"""
UV sphere
"""
dmn_capsule = PDM.sphere_surf_gen_nodal(comm,
                                        nu       = 2*n,
                                        nv       = n,
                                        x_center = 0,
                                        y_center = 0,
                                        z_center = 0,
                                        radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_SURFACIC,
                       "visu_sphere_")

"""
Icosphere
"""
dmn_capsule = PDM.sphere_surf_icosphere_gen_nodal(comm,
                                                  n,
                                                  x_center = 0,
                                                  y_center = 0,
                                                  z_center = 0,
                                                  radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_SURFACIC,
                       "visu_icosphere_")


"""
Icoball
"""
dmn_capsule = PDM.sphere_vol_icosphere_gen_nodal(comm,
                                                 n,
                                                 x_center = 0,
                                                 y_center = 0,
                                                 z_center = 0,
                                                 radius   = 1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                       "visu_icoball_")

"""
Hollow ball
"""
dmn_capsule = PDM.sphere_vol_hollow_gen_nodal(comm,
                                              n,
                                              n_layer         = 2*n,
                                              x_center        = 0,
                                              y_center        = 0,
                                              z_center        = 0,
                                              radius_int      = 1,
                                              radius_ext      = 2,
                                              geometric_ratio = 1.1)

if visu:
  dmn_capsule.dump_vtk(PDM._PDM_GEOMETRY_KIND_VOLUMIC,
                       "visu_hollow_ball_")


"""
Dcube nodal
"""
dico = {
  # "_PDM_MESH_NODAL_POINT"    : PDM._PDM_MESH_NODAL_POINT   ,
  # "_PDM_MESH_NODAL_BAR2"     : PDM._PDM_MESH_NODAL_BAR2    ,
  "_PDM_MESH_NODAL_TRIA3"    : PDM._PDM_MESH_NODAL_TRIA3   ,
  "_PDM_MESH_NODAL_QUAD4"    : PDM._PDM_MESH_NODAL_QUAD4   ,
  "_PDM_MESH_NODAL_POLY_2D"  : PDM._PDM_MESH_NODAL_POLY_2D ,
  "_PDM_MESH_NODAL_TETRA4"   : PDM._PDM_MESH_NODAL_TETRA4  ,
  "_PDM_MESH_NODAL_PYRAMID5" : PDM._PDM_MESH_NODAL_PYRAMID5,
  "_PDM_MESH_NODAL_PRISM6"   : PDM._PDM_MESH_NODAL_PRISM6  ,
  "_PDM_MESH_NODAL_HEXA8"    : PDM._PDM_MESH_NODAL_HEXA8   ,
  # "_PDM_MESH_NODAL_POLY_3D"  : PDM._PDM_MESH_NODAL_POLY_3D ,
  }


for t in dico:
  dcube = PDM.DCubeNodalGenerator(nx     = n,
                                  ny     = n,
                                  nz     = n,
                                  length = 1,
                                  zero_x = 0,
                                  zero_y = 0,
                                  zero_z = 0,
                                  t_elmt = dico[t],
                                  order  = 1,
                                  comm   = comm)

  dcube.compute()

  dmn_capsule = dcube.get_dmesh_nodal()

  if t == "_PDM_MESH_NODAL_TRIA3" \
  or t == "_PDM_MESH_NODAL_QUAD4" \
  or t == "_PDM_MESH_NODAL_POLY_2D":
    geom_kind = PDM._PDM_GEOMETRY_KIND_SURFACIC
  else:
    geom_kind = PDM._PDM_GEOMETRY_KIND_VOLUMIC

  if visu:
    dmn_capsule.dump_vtk(geom_kind,
                         f"visu_dcube_{t[16:].lower()}")


"""
Ngon
"""
mesh = PDM.generate_mesh_rectangle_ngon(comm        = comm,
                                        elt_type    = PDM._PDM_MESH_NODAL_POLY_2D,
                                        xmin        = 0.,
                                        ymin        = 0.,
                                        zmin        = 0.,
                                        lengthx     = 1.,
                                        lengthy     = 1.,
                                        n_x         = n,
                                        n_y         = n,
                                        n_part      = n_part,
                                        part_method = PDM._PDM_SPLIT_DUAL_WITH_HILBERT)

# Randomize
for i_part in range(n_part):
  for i_vtx in range(mesh["pn_vtx"][i_part]):
    np.random.seed(mesh["pvtx_ln_to_gn"][i_part][i_vtx])
    for i in range(2):
      x = mesh["pvtx_coord"][i_part][3*i_vtx + i]
      if x > 1e-6 and x < 1-1e-6:
        mesh["pvtx_coord"][i_part][3*i_vtx + i] = x + 0.1*(2*np.random.rand() - 1)/(n - 1)

wrt = PDM.Writer("Ensight",
                 PDM._PDM_WRITER_FMT_ASCII,
                 PDM._PDM_WRITER_TOPO_VARIABLE,
                 PDM._PDM_WRITER_OFF,
                 "mesh_gen_demo",
                 "rectange_ngon",
                 comm,
                 PDM._PDM_IO_KIND_MPI_SIMPLE,
                 1.,
                 "")

id_geom = wrt.geom_create("rectange_ngon",
                          n_part)

id_var_part = wrt.var_create(PDM._PDM_WRITER_OFF,
                             PDM._PDM_WRITER_VAR_SCALAR,
                             PDM._PDM_WRITER_VAR_ELEMENTS,
                             "i_part")

wrt.step_beg(0.)

for i_part in range(n_part):
  wrt.geom_coord_set(id_geom,
                     i_part,
                     # mesh["pn_vtx"][i_part],
                     mesh["pvtx_coord"][i_part],
                     mesh["pvtx_ln_to_gn"][i_part])

  wrt.geom_cell2d_cellface_add(id_geom,
                               i_part,
                               # mesh["pn_face"][i_part],
                               # mesh["pn_edge"][i_part],
                               # None,
                               # None,
                               mesh["pedge_vtx"][i_part],
                               mesh["pface_edge_idx"][i_part],
                               # None,
                               mesh["pface_edge"][i_part],
                               mesh["pface_ln_to_gn"][i_part])

  wrt.var_set(id_var_part,
              id_geom,
              i_part,
              (n_part*i_rank + i_part)*np.ones(mesh["pn_face"][i_part], dtype=np.double))

wrt.geom_write(id_geom)

wrt.var_write(id_var_part)

wrt.step_end()

"""
Polyvol gen
"""
# ...
