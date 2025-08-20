.. _change_log:

Release notes
#############


Version 2.6.0 (December 2024)
-----------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_dmesh(_nodal)` : Add function to find topological ridges
- :pdmkw:`PDM_dmesh_extract` : Extract groups (in nodal form)
- :pdmkw:`PDM_mem_tool` : Overload malloc, calloc, realloc and free with macros for better memory management
- :pdmkw:`PDM_gnum` : Complete Fortran API
- :pdmkw:`PDM_part_mesh_nodal` : Complete Fortran API
- :pdmkw:`PDM_multipart` : Complete Python API
- :pdmkw:`PDM_extract_part` : Local and reequilibrate modes for nodal meshes + Automatic extraction of lower-dimension elements
- :pdmkw:`PDM_geom_elemt` : Add up/downwind edges in 2D (in any axis-aligned plane)
- :pdmkw:`PDM_part_mesh_nodal_to_part_mesh` : Customizable conversion of part_mesh_nodal to part_mesh (connectivities, gnums, groups, link elmt->entity)
- :pdmkw:`PDM_part_comm_graph` : Symmetric communication graph for inter-partition exchanges (without gnums)
- :pdmkw:`PDM_isosurface` : Generation of iso-surface/line meshes in 3d/2d partitioned or distributed meshes (nodal or ngon)
- :pdmkw:`PDM_domain_interface` : Complete Fortran API
- :pdmkw:`PDM_part_extension` : Support infinite periodicity, with multiple periodic boundaries
- :pdmkw:`PDM_part_to_part` : add collective mode for exchanges and reverse exchanges
- :pdmkw:`PDM_part_mesh_reorient_geom` : reorient a mesh geometrically

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_extract_part` : Better management of ownerships
- Transfer Gitlab CI from Spiro to Juno

❌ Removed
^^^^^^^^^^
- :pdmkw:`PDM_mesh_nodal` : Replace struct with :pdmkw:`PDM_part_mesh_nodal` + Remove redundant functions
- :pdmkw:`PDM_part_mesh_nodal_to_pmesh` : Replace free function by struct :pdmkw:`PDM_part_mesh_nodal_to_part_mesh`

🔧 Fixed
^^^^^^^^
- :pdmkw:`PDM_reader_gamma` : Fix read of meshes with multiple elements
- :pdmkw:`PDM_mesh_intersection` : Manage cases with zero candidates for intersection
- :pdmkw:`PDM_distrib` : Fix division by zero in empty distributions
- :pdmkw:`PDM_mesh_location` : Fix when all target points are outside of the global source bounding box
- :pdmkw:`PDM_dcube_gen` : Fix orientation of boundary faces
- :pdmkw:`PDM_extract_part` : Allow n_part_out != n_part_in in reequilibrate mode
- :pdmkw:`PDM_domain_interface` : Fix face to vtx conversion
- :pdmkw:`PDM_io` : Support relative paths
- :pdmkw:`PDM_part_to_part` : Fix Fortran wrapping of part_to_part_reverse_iexch

|

Version 2.5.0 (May 2024)
------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_reader_gamma` : Add prisms, hexahedra, pyramids, quadrilaterals (C and Python API)
- :pdmkw:`PDM_dmesh_nodal` : Mesh reorientation + Find ridges in surface mesh + Fortran API
- Documentation : :pdmkw:`PDM_io`, :pdmkw:`PDM_writer`
- :pdmkw:`PDM_extract_part` : Entities renumbering
- :pdmkw:`PDM_multipart` : Entities renumbering
- :pdmkw:`PDM_part_to_block` : Timers + export of the topology of the communication graph + automatic switchover to point-to-point communication if the communication graph is sparse + improved exchange capacity (> 2GB)
- :pdmkw:`PDM_block_to_part` : Automatic switchover to point-to-point communication if the communication graph is sparse + improved exchange capacity (> 2GB)
- :pdmkw:`PDM_domain_interface` : :pdmkw:`PDM_domain_interface_translate_entity1_entity2` in Python API
- :pdmkw:`PDM_part_to_part` : improved exchange capacity (> 2GB)

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_block_to_block` : Improvement (Binary search replaced by analytical formula)
- :pdmkw:`PDM_multipart` : Add some checks and conditions for entities reordering
- :pdmkw:`PDM_multipart` : For a single MPI rank, local and global numbering are identical
- :pdmkw:`PDM_extract_part` : 1-based Entities numbering

🔧 Fixed
^^^^^^^^
- :pdmkw:`PDM_part_to_block` : handle global sum of weights equal to zero
- :pdmkw:`PDM_part_to_part` : Memory leak (local copy of the MPI communicator) + request storage
- :pdmkw:`PDM_mesh_intersection` : No candidate for the intersection
- :pdmkw:`PDM_part`/:pdmkw:`pdm_multipart` : Hilbert with number of ranks > number of mesh entities
- :pdmkw:`PDM_extract_part` : parent local numbering
- :pdmkw:`PDM_dmesh_nodal_to_dmesh` : Python API

|

Version 2.4.0 (January 2024)
----------------------------

🌱 Added
^^^^^^^^
- Training : Jupyter notebooks in 3 languages (C, Fortran, Python)
- Documentation : Partial Sphinx documentation
- Fortran API : :pdmkw:`PDM_dmesh_nodal`, :pdmkw:`PDM_dist_cloud_surf`, :pdmkw:`PDM_mesh_intersection`, :pdmkw:`PDM_part_connectivity_transform`, :pdmkw:`PDM_dcube_nodal_gen`
- Python API : :pdmkw:`PDM_distrib`, :pdmkw:`PDM_multi_block_to_part`
- 7 Tests + 2 Tutorials

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_mesh_location` : API
- :pdmkw:`PDM_multipart` : API + Add 1D and 0D
- :pdmkw:`PDM_part_extension` : API

|

Version 2.3.0 (June 2023)
-------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_generate_mesh` : Unique API for simple mesh generation
- :pdmkw:`PDM_dmesh_extract` : Extract mesh in distributed blocks

⚠️ Changed
^^^^^^^^^^
- CI : Improvment
- :pdmkw:`PDM_mesh_location` : Add :pdmkw:`PDM_MESh_LOCATION_LOCATE_ALL_TGT` mode (All points are located)
- :pdmkw:`PDM_multipart` : Add 0D, 1D and 2D meshes
- :pdmkw:`PDM_extract_part` : Add 0D and 1D and meshes + Fix npart=0 case + groups preservation
- :pdmkw:`PDM_part_to_part` : Add API to setup a user buffer + API Fortran
- :pdmkw:`PDM_iso_surface` : Groups preservation + Fix polygon mode
- Cython : Improvment (ownership)

|

Version 2.2.0 (March 2023)
--------------------------

🌱 Added
^^^^^^^^
- migration to the internal gitlab server
- CI with LeakSanitizer
- :pdmkw:`pdm_run` : script to simplify the launch of tests
- :pdmkw:`PDM_field_cell_to_vtx` : Python interface
- :pdmkw:`PDM_global_mean` : Python interface
- :pdmkw:`PDM_inria_mesh_fmt` : Python interface
- :pdmkw:`PDM_mesh_intersection` : Python interface
- :pdmkw:`PDM_sphere_surf_gen` : Python interface
- :pdmkw:`PDM_iso_surface.pxi` (ParaDiGMA) : Python interfaceq
- :pdmkw:`PDM_field_cell_to_vtx` : Transfer a cell field to a vertex field
- :pdmkw:`PDM_mesh_nodal` : Fortran interface
- :pdmkw:`PDM_part_mesh_nodal` : Fortran interface
- :pdmkw:`PDM_multipart` : Fortran interface
- 30 Tests (including 3 for ParaDiGMA)

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_dist_cloud_surf` : :pdmkw:`PDM_part_mesh_nodal` is accepted in input
- :pdmkw:`PDM_ho_location` : Implement a Newton method for all elements (for all orders)
- :pdmkw:`PDM_mesh_intersection` : :pdmkw:`PDM_part_mesh_nodal` is accepted in input
- :pdmkw:`PDM_part_mesh_nodal` : Replace :pdmkw:`PDM_part_mesh` (more general and generic structure)
- :pdmkw:`PDM_part_mesh_nodal_to_pmesh` : Compute faces and (or) edges from a :pdmkw:`PDM_part_mesh_nodal`
- :pdmkw:`PDM_lagrange_to_bezier` : Implement all standard high-order elements
- :pdmkw:`PDM_mesh_location` : Add high order elements + switch to :pdmkw:`PDM_part_mesh_nodal`
- :pdmkw:`PDM_point_location` : Use the Newton method primarily
- :pdmkw:`PDM_writer` : Switch to :pdmkw:`PDM_part_mesh_nodal` structure
- :pdmkw:`PDM_part` : Call :pdmkw:`PDM_multi_part` from :pdmkw:`PDM_part` if the environment variable :pdmkw:`PDM_USE_MULTIPART` is defined
- :pdmkw:`PDM_part_renum` : Renumbering of groups of entities (face, edge, vertex)
- :pdmkw:`PDM_multi_part` : :pdmkw:`PDM_part_mesh_nodal` as output + HO
- :pdmkw:`PDM_partitioning_algorithm` : Support HO elements
- :pdmkw:`PDM_partitioning_nodal_algorithm` : Support HO elements
- :pdmkw:`PDM_box_tree` : An user init location is accepted in input
- :pdmkw:`PDM_extract_part` : :pdmkw:`PDM_part_mesh_nodal` is accepted in input
- :pdmkw:`PDM_gnum` : Build global numbering from tuple + Optimization
- :pdmkw:`PDM_part_mesh` : Manage group (Face, edge, vertex) + boundaries
- :pdmkw:`PDM_iso_surface` (ParaDiGMA) : Transfer boundary data to iso_surface (beta)

🧊 Deprecated
^^^^^^^^^^^^^
- :pdmkw:`PDM_mesh_nodal` : Replace with :pdmkw:`PDM_part_mesh_nodal`

❌ Removed
^^^^^^^^^^
- Nuga sources (ParaDiGMA)

🔧 Fixed
^^^^^^^^
- :pdmkw:`PDM_para_octree` (fix multiple points)
- :pdmkw:`PDM_part_extension` (fix n_part > 1)
- :pdmkw:`PDM_part_to_block` (fix n_part > 1)
- :pdmkw:`PDM_point_tree_seq`

|

Version 2.1.0 (December 2022)
-----------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_ho_bezier` : Elementary functions for high-order Bezier elements
- :pdmkw:`PDM_ho_bezier_basis` : Bezier basis functions
- :pdmkw:`PDM_lagrange_to_bezier` : Conversion for Lagarange basis to Bezier basis
- :pdmkw:`PDM_box_gen` : Generate sets of boxes
- :pdmkw:`PDM_sphere_vol_gen` : Generate 3D mesh of a sphere
- :pdmkw:`PDM_reader_gamma` : Read GAMMA mesh format (partial implementation)
- :pdmkw:`PDM_reader_stl` : Read STL mesh format

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_dmesh_nodal_reorder` : remove 'order' argument (minor change of API)
- :pdmkw:`PDM_mesh_location` : Add the communication graph to the results
- :pdmkw:`PDM_mesh_location` : New optimized algorithm
- :pdmkw:`PDM_Mesh_nodal` : cell3d_cellface_add (minor change of API)
- :pdmkw:`PDM_multipart`  : get_part_mesh_nodal (minor change of API (add ownership))
- :pdmkw:`PDM_part_to_block` : Add reverse exchange
- :pdmkw:`PDM_part_to_part` : Add :pdmkw:`PDM_part_to_part_create_from_num2_triplet`
- :pdmkw:`PDM_extract_part` : Improved and optimized
- :pdmkw:`PDM_dist_cloud_surf` : Improved memory efficient and performance (optional, set :pdmkw:`PDM_DIST_CLOUD_SURF_OPTIM` to 1 to use it)
- :pdmkw:`PDM_writer` : Concatenate all time step into a single file per field

🔧 Fixed
^^^^^^^^
- bugfix : :pdmkw:`PDM_part_to_part`
- bugfix : :pdmkw:`PDM_extract_part`

|

Version 2.0.0 (June 2022)
-------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_vtk` : Write ASCII VTK files for visulization with ParaView
- :pdmkw:`PDM_dgeom_elem` : Multiple method to compute center of entity with distributed conectivities and coordinates (usefull for hilbert ordering)
- :pdmkw:`PDM_dmesh_nodal_elmts` : Sub-structure of dmesh_nodal that contains for one geom kind (volumic/surfacic/ridge/corner) the associate element connectivtities
- :pdmkw:`PDM_domain_interface` : Structure and algorithme to manage interface between domain (for exemple between 2 mesh with 2 separate global numbering). Algorithm is useful to deduce other connectivities betwenn domain (for exemple face -> vtx)
- :pdmkw:`PDM_ho_ordering` : Generic definition of high order elements
- :pdmkw:`PDM_part_domain_interface` : Same as :pdmkw:`PDM_domain_interface` but in partitioned view (usefull for part_extension)
- :pdmkw:`PDM_predicate` : Numerically robust geometric predicates
- :pdmkw:`PDM_point_cloud_gen` : Generate a random point cloud
- :pdmkw:`PDM_sphere_surf_gen` : Generate a distributed suface mesh of a sphere
- :pdmkw:`PDM_extract_part` : Build a child mesh from the extraction of selected elements in the parent mesh
- :pdmkw:`PDM_global_reduce` : Global reduce applied to a partitionned field
- :pdmkw:`PDM_part_to_part` : Browse an MPI communication graph defined by two global numberings
- :pdmkw:`PDM_pointer_array` : Fortran derived type for interfacing C pointers of arrays
- :pdmkw:`PDM_hkey` : Compute an hash key from a list of integer
- :pdmkw:`PDM_memory_stats` : Structure to monitor the memory uses by snapshot

⚠️ Changed
^^^^^^^^^^
- New Fortran API based on iso-c-binding
- Adds functionalities in Fortran API
- Adds functionalities in Python API

❌ Removed
^^^^^^^^^^
- Removes handles in some functionalities

|

Version 1.13.0 (July 2021)
--------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_part_extension`: Generate extension of a existing partition for graph, usefull to setup ghost cell
- :pdmkw:`PDM_distant_neighbor`: Exchange protocol in order to communicate by triplet (i_proc, i_part, i_entity)
- :pdmkw:`PDM_mesh_adapt`: First API of mesh adaptation
- :pdmkw:`PDM_dcube_nodal_gen`: Cube generator with element connectivity (3D : Hexa/Pyra/Tetra and 2D : Quad/Tri )
- :pdmkw:`PDM_part1_to_selected_part2`: Exchange protocol between two partitioned numbering
- :pdmkw:`PDM_interpolate_from_mesh_location`: Simple inplementation of exhange after a localisation (can be replace by :pdmkw:`PDM_part1_to_selected_part2`)
- :pdmkw:`PDM_dmesh_nodal`: Low level structure to describe a mesh by elemt in distributed maner
- :pdmkw:`PDM_dmesh_nodal_to_dmesh`: Algorithm to deduce from mesh_nodal connectivity a descending connectivity (ex: cell_vtx ---> face_cell)
- :pdmkw:`PDM_poly_vol_gen`: Polyhedral mesh generation (:pdmkw:`PDM_poly_surf` extruded)
- :pdmkw:`PDM_array`: Generic functions for arrays
- :pdmkw:`PDM_partitioning_nodal_algorithm`: A collection of algorithm to help partitionning of dmesh_nodal (Beta)

⚠️ Changed
^^^^^^^^^^
- ParaDiGMA is an extension of ParaDiGM and becomes a git submodule of ParaDiGM
- :pdmkw:`PDM_closest_points`: Greatly improves robustness and efficiency
- :pdmkw:`PDM_mesh_location`: Greatly improves robustness and efficiency
- :pdmkw:`PDM_para_octree`: Greatly improves robustness and efficiency
- :pdmkw:`PDM_dbbtree`: Greatly improves robustness and efficiency
- :pdmkw:`PDM_overlay`: Greatly improves robustness and efficiency
- :pdmkw:`PDM_dist_cloud_surf`: Greatly improves robustness and efficiency
- CMake 3.0 -> CMake 3.19
- Adds functionalities in Fortran API
- Adds functionalities in Python API
- Removes handles in some functionalities

|

Version 1.12.5 (July 2021)
--------------------------

🔧 Fixed
^^^^^^^^
bugfix : :pdmkw:`pdm_writer_ensight`

|

Version 1.12.4 (March 2021)
---------------------------

🔧 Fixed
^^^^^^^^
bugfix : :pdmkw:`pdm_hilbert` 2D

|

Version 1.12.3 (March 2021)
---------------------------

🔧 Fixed
^^^^^^^^
bugfix : renum cacheblocking

|

Version 1.12.2 (February 2021)
------------------------------

🔧 Fixed
^^^^^^^^
bugfix : :pdmkw:`pdm_box_tree` (normalization)
bugfix : :pdmkw:`pdm_dbbtree` (normalization)

|

Version 1.12.1 (January 2021)
-----------------------------

|

Version 1.12.0 (November 2020)
------------------------------

🌱 Added
^^^^^^^^
- API :

  - :pdmkw:`PDM_overlay` : Intersection of two plane meshes
  - :pdmkw:`PDM_closest_points`     : Parallel knn algorithm
  - :pdmkw:`PDM_dmesh_partitioning` : New implementation of parallel meshe partitioning (old implementation : :pdmkw:`PDM_part`)
  - :pdmkw:`PDM_multipart` : Parallel meshes partitioning with preservation of links between meshes
  - :pdmkw:`PDM_mesh_location` : A partitioned point cloud location in a partitionned mesh
  - :pdmkw:`PDM_distant_neighbor` : Exchange data between distant neighbors (only collective communication for the moment)

- Internal functions :

  - :pdmkw:`PDM_mean_values` : Polygon and Polyedra mean values (generalization of barycentric coordinates)
  - :pdmkw:`PDM_ho_basis` : High order basis for tetrahedron, hexahedron, prism, pyramid, triangle and quadrangle
  - :pdmkw:`PDM_ho_location` : Point location into a high order curved element
  - :pdmkw:`PDM_logging` : Managing a logging file
  - :pdmkw:`PDM_dconnectivity_transform` : Set of function to perform distributed operation on connectivity (combine, transform, reverse)
  - :pdmkw:`PDM_dmesh_nodal_elements_utils` : Operations about a dmesh_nodal mesh
  - :pdmkw:`PDM_point_location` : Point location in a :pdmkw:`PDM_mesh_nodal`
  - :pdmkw:`PDM_tetrahedron` : Elementary functions about tetrahedron
  - :pdmkw:`PDM_triangulate` : Polygon triangulation
  - :pdmkw:`PDM_para_graph_dual` : Parallel building of a dual graph mesh
  - :pdmkw:`PDM_partitioning_algorithm` : Steps of the partitioning algorithm called from :pdmkw:`PDM_dmesh_partitioning`
  - :pdmkw:`PDM_compare_operator` : Implementation of comparaison operator for string and connectivity (useful for :pdmkw:`PDM_gnum_from_hash_values`)
  - :pdmkw:`PDM_dmesh` : Distributed mesh defined from a top-down connectivity
  - :pdmkw:`PDM_equal_operator` : Implementation of equal operator for string and connectivity (useful for :pdmkw:`PDM_gnum_from_hash_values`)
  - :pdmkw:`PDM_gnum_from_hash_values` : Structure to build a global numbering from keys, useful to building edge global numebring or face global numbering. The idea is to use a distribution of key for each rank, sort locally and resend to original distribution

- Fortran interface : :pdmkw:`PDM_closest_points`, :pdmkw:`PDM_mesh_location`, :pdmkw:`PDM_block_to_part`, :pdmkw:`PDM_part_to_block` :pdmkw:`PDM_overlay`
- Cython interface : :pdmkw:`PDM_closest_points`, :pdmkw:`PDM_dcube_gen`, :pdmkw:`PDM_distant_neighbor`, :pdmkw:`PDM_dmesh`, :pdmkw:`PDM_gnum`, :pdmkw:`PDM_gnum_location`, :pdmkw:`PDM_mesh_location`, :pdmkw:`PDM_multi_part` :pdmkw:`PDM_points_merge` :pdmkw:`PDM_overlay`
- C tests : :pdmkw:`pdm_t_closest_points`, :pdmkw:`pdm_t_create_edge_gnum`, :pdmkw:`pdm_t_distant_neighbor`, :pdmkw:`pdm_t_distant_neighbor_nomatch`, :pdmkw:`pdm_t_gnum_from_hash_values`, :pdmkw:`pdm_t_knn_cube`, :pdmkw:`pdm_t_mesh_location_dcube`, :pdmkw:`pdm_t_mesh_location_poly_surf`, :pdmkw:`pdm_t_multipart`, :pdmkw:`pdm_t_order`, :pdmkw:`pdm_t_partitioning_dcube` :pdmkw:`pdm_t_plane_intersect_unit`, :pdmkw:`pdm_t_plane_meshes_intersect_unit`, :pdmkw:`pdm_t_plane_meshes_intersect`
- Fortran tests : :pdmkw:`pdm_t_block_to_part_f`, :pdmkw:`pdm_t_mesh_location_f`, :pdmkw:`pdm_t_part_to_block_f` :pdmkw:`pdm_t_plane_intersect_unit_f`
- Unitests from DocTest : :pdmkw:`pdm_dconnectivity_transform.test`, :pdmkw:`pdm_dmesh_nodal.test`, :pdmkw:`pdm_dmesh_nodal_elements_utils`, block_to_part.test

|

Version 1.11.0 (December 2019)
------------------------------

🌱 Added
^^^^^^^^
- :pdmkw:`PDM_para_octree` : Parallel linear octree (partial implementation)

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`PDM_block_to_part` : Point to point communication if buffers are too large
- :pdmkw:`PDM_dbbtree` : Change default parameters
- :pdmkw:`PDM_mesh_nodal` : Add cell centers
- :pdmkw:`PDM_MPI` : Add large allToall (sdispls and rdispls in int64)

❌ Removed
^^^^^^^^^^
- :pdmkw:`PDM_box_tree` :  Remove extents attribute

🔧 Fixed
^^^^^^^^
- Windows portage
- Corrections for :pdmkw:`PDM_block_to_part`
- Corrections for :pdmkw:`PDM_dist_cloud_surf`
- :pdmkw:`PDM_MPI_Comme_split` : memory leaks

|

Version 1.10.0 (June 2019)
--------------------------

🌱 Added
^^^^^^^^
- ChangeLog file
- :pdmkw:`PDM_gnum_location` : Give the location of an entity (process, number of partition, local number) from its global numbering

⚠️ Changed
^^^^^^^^^^
- :pdmkw:`pdm_mesh_dist` -> :pdmkw:`pdm_dist_cloud_surf`

🔧 Fixed
^^^^^^^^
- Correction for multigrid
- Correction for part_to_block