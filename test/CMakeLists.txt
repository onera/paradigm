# --------------------------
# Construction des Tests
# --------------------------
# > A gerer, attention a Robin
# if (PDM_ENABLE_CUDA)
#   list(APPEND LINK_LIBRARIES pdm_gpu)
# endif()

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
if (PDM_ENABLE_EXTENSION_PDMA)
  file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/meshes DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/../extensions/paradigma/test/")
endif()

# Test pdm_t_part_cube
# --------------

set (MPIEXEC_MAX_NUMPROCS 2)
set (MPIEXEC_EXECUTABLE mpirun)

include(TestCreate)

# NB: Faire le passage dans le -genv de MPI n'est pas équivalent à faire l'export avant ...
execute_process(COMMAND sh -c "${MPIEXEC} --version | grep -c \"Open MPI\"" OUTPUT_VARIABLE IS_OPENMPI OUTPUT_STRIP_TRAILING_WHITESPACE)
set (MPIEXEC_GENV_COMMAND      "")
set (MPIEXEC_GENV_PRELOAD      "")
set (MPIEXEC_GENV_PRELOAD_PATH "")
if (CMAKE_BUILD_TYPE STREQUAL "Sanitize")
  execute_process(COMMAND gcc -print-file-name=libasan.so OUTPUT_VARIABLE PRELOAD_ASAN OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(IS_OPENMPI)
    set(MPIEXEC_GENV_COMMAND      "-x")
    set(MPIEXEC_GENV_PRELOAD      "")
    set(MPIEXEC_GENV_PRELOAD_PATH "LD_PRELOAD=${PRELOAD_ASAN}:${PDM_BINARY_DLCLOSE_DIR}/libdlclose.so")
  else()
    set(MPIEXEC_GENV_COMMAND      "-genv")
    set(MPIEXEC_GENV_PRELOAD      "LD_PRELOAD")
    set(MPIEXEC_GENV_PRELOAD_PATH "${PRELOAD_ASAN}:${PDM_BINARY_DLCLOSE_DIR}/libdlclose.so")
  endif()

endif()

if(PDM_ENABLE_UNIT_TEST)
  file(GLOB_RECURSE __test_files
       CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../src/*.test.cpp")
  test_cpp_unit_create(pdm_t_unit_test ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK SOURCES ${__test_files} )
endif()

if (PARMETIS_FOUND OR PTSCOTCH_FOUND)
  test_c_create (pdm_t_plane_agglo 2  PDM_LIST_TEST PDM_LIST_NRANK)
  if (PARMETIS_FOUND)
    test_c_create (pdm_t_partitioning_dcube ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
    test_c_create (pdm_t_part_dcube_agglo ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
  endif()
endif()

# gnum_utils
test_c_create (pdm_t_gnum_nuplet               ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_gnum_from_part_comm_graph ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_gnum                                            4 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_gen_gnum                                        2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_gnum_location                                   4 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_gnum_with_merge                                 4 PDM_LIST_TEST PDM_LIST_NRANK)

# meshgen
test_c_create (pdm_t_generate_mesh      2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_polysurf_extrusion 2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dcube_nodal_gen    2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_cube_nodal    1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_icoball            1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_sphere_vol_gen     1 PDM_LIST_TEST PDM_LIST_NRANK)

# exchange_protocol
test_c_create (pdm_t_part1_to_selected_part2_from_mesh_location 4 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_distant_neighbor3                          4 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_distant_neighbor                           1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_distant_neighbor_nomatch                   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_distant_neighbor2                          2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_to_block_geom                         1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_to_part                               2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_comm_graph_periodic_2d                2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_pcg_group_vtx                              3 PDM_LIST_TEST PDM_LIST_NRANK)

# ppart
test_c_create (pdm_t_part_dcube               ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_dcube_nodal         ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_extension                                 2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_partitioning_dcube_nodal                       2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_multipart                ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_multipart_0d                                   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_multipart_1d                                   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_multipart_1d_nodal                             1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_extension_2d_draft                        2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_extension_3d_draft                        2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_join_domains                                   2 PDM_LIST_TEST PDM_LIST_NRANK)


# mesh_topology
test_c_create (pdm_t_dmesh_nodal_to_dmesh   2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_to_dmesh_nodal   2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_find_ridges_3d   2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_dmesh_extract_nodal    2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_extract_nodal_2d 2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_extract_3d       2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_extract_2d       2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dmesh_extract_1d       2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_extract_part                     1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_vtx                 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_nodal               1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_from_target         1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_reequilibrate       1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_nodal_from_target   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_nodal_local         1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_part_nodal_reequilibrate 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_extract_surf_part                1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_cellface_orient                   ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_domain_interface_vtx_to_edge                       2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_domain_interface_vtx_translate                          2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_compute_part_edges                           1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_mesh_nodal_to_part_mesh                 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_part_mesh_nodal_elmts_compute_child_entities 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_mesh_check              1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_multi_block_merge       1 PDM_LIST_TEST PDM_LIST_NRANK)

# utils
test_c_create (pdm_t_global_mean   ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_global_reduce ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)

# geometry
test_c_create (pdm_t_dist                     2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_nodal               2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_ho                  2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_strip               2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_rocket              2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_axitransbump        2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_bl                  1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_cloud_surf_from_stl 2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_point_merge              2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_knn_cube                     ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_mesh_location_dcube          ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_location_disjoint_cubes ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_location_multisection   ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_location_poly_surf      ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_location_from_file                            2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_location_ho                                   2 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_plane_meshes_intersect                             2 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_field_cell_to_vtx       1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_field_cell_to_vtx_nodal 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_laplacian_smoothing 2 PDM_LIST_TEST PDM_LIST_NRANK)

# io
test_c_create (pdm_t_writer_ensight_3d               ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_writer_ensight_moving_mesh_3d   ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_writer_ensight_variable_mesh_3d ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_reader_gamma                                          1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_reader_vtk                                            1 PDM_LIST_TEST PDM_LIST_NRANK)

#
test_c_create (pdm_t_mesh_intersection_surf_surf 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_intersection_vol_surf  1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_intersection_vol_line  1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_intersection_vol_vol   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_intersection_surf_line 1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_edge_upwind_and_downwind    1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_edge_upwind_and_downwind_2d 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_overlay_from_vol_mesh       1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_plane_intersect_unit        2 PDM_LIST_TEST PDM_LIST_NRANK)

# geometrees
test_c_create (pdm_t_par_octree                ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_para_octree                                     1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_octree                    ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_doctree                   ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_kdtree                    ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_point_tree                ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_kdtree_points_inside_ball ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_par_octree_shared         ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)

# algo
test_c_create (pdm_t_closest_points                 ${MPIEXEC_MAX_NUMPROCS} PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_plane_meshes_intersect_unit                          1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_intersect_line_box                                   1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_closest_points_1_point 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_dist_ridge             1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_point_location_nodal        1  PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_projection_point_cloud_surf 1 PDM_LIST_TEST PDM_LIST_NRANK)


test_c_create (pdm_t_intersect_tetrahedra     1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_triangles_intersect      1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_unit_tetra_intersect_triangle      1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_intersection_surf_surf_atomic 1 PDM_LIST_TEST PDM_LIST_NRANK)


test_c_create (pdm_t_isosurface_2d_nodal 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_isosurface_3d_nodal 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_isosurface_2d_ngon  1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_isosurface_3d_ngon  1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_triangulate        1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mesh_reorient_geom 3 PDM_LIST_TEST PDM_LIST_NRANK)

# misc
test_c_create (pdm_t_mesh_adapt   1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_newton_ho    1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_svd          1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_mixing_plane 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_ho_ordering  1 PDM_LIST_TEST PDM_LIST_NRANK)

test_c_create (pdm_t_point_in_polygon            1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_line_intersection           1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_ho_location                 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_subdivision_bezier_triangle 1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_projection_back_mesh        1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_subdivision_bezier_curve    1 PDM_LIST_TEST PDM_LIST_NRANK)
test_c_create (pdm_t_ho_bezier_intersect         1 PDM_LIST_TEST PDM_LIST_NRANK)

if (PDM_ENABLE_PYTHON_BINDINGS)
   test_python_create (pdm_t_part_to_part_p             2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (test_demr                        2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_mesh_gen_demo_p            2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_gnum_p                     2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_block_to_part_p            3 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_part_comm_graph_p          2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_part_to_block_p            3 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_isosurface_pacman_p        2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_isosurface_user_function_p 2 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_mesh_reorient_geom_p       3 PDM_LIST_TEST PDM_LIST_NRANK)
   test_python_create (pdm_t_generate_graph_comm_p      3 PDM_LIST_TEST PDM_LIST_NRANK)
   if (PDM_ENABLE_BLASLAPACK)
     test_python_create (pdm_t_rotation_p 1 PDM_LIST_TEST PDM_LIST_NRANK)
   endif()
endif()

if(PDM_ENABLE_Fortran)
  test_fortran_create (pdm_t_closest_points_f             2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_mesh_location_f              1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_to_block_f              1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_to_block_from_distrib_f 2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_block_to_part_f              1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_dcube_gen_f                  1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_dcube_nodal_gen_f            2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_plane_intersect_unit_f       1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_global_reduce_f              2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_dist_f                       2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_gnum_f                       2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_f                       2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_writer_f                     2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_mesh_check_f                 1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_to_part_f               2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_gnum_location_f              2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_extension_f             2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_extension_draft_f       2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_sphere_surf_gen_f            1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_multipart_f                  2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_part_mesh_nodal_cellvtx_f    2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_generate_mesh_f              2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_intersection_prepro          2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_extract_part_f               2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_extract_part_nodal_f         2 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_point_in_convex_hull_f       1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_isosurface_3d_ngon_f         1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_isosurface_3d_nodal_f        1 PDM_LIST_TEST PDM_LIST_NRANK)
  test_fortran_create (pdm_t_isosurface_nodal_f           1 PDM_LIST_TEST PDM_LIST_NRANK)
endif()


add_subdirectory(tutorial)

set (PDM_LIST_TEST  ${PDM_LIST_TEST}  PARENT_SCOPE )
set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )
