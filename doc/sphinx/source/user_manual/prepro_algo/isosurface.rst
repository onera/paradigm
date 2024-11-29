.. _isosurface:

Iso-surfaces
============

**pdm_isosurface** provides a distributed and partitioned API for isosurface mesh generation
from 2D or 3D meshes. Entry mesh can be defined through a **nodal mesh** (volumic and boundaries elements)
or a **mesh** (cell, face, (edge) and vertices connectivities description). Depending of entry, algorithm won't be the same.

If entry mesh is a nodal mesh fully composed of tetrahedras and triangles, a edge based variation of
`marching tetrahedras <https://fr.wikipedia.org/wiki/Marching_tetrahedra>`_ algorithm will be used.

For multi-elements nodal meshes and for meshes, the `ngon algorithm <https://www.sciencedirect.com/science/article/pii/S0021999121004745>`_ 
will be used.

In any case, if entry is a 2D mesh, elements
will be triangulated to use the edge based variation of `marching tetrahedras <https://fr.wikipedia.org/wiki/Marching_tetrahedra>`_ algorithm.


C API
-----


Initialization
""""""""""""""

.. doxygenfunction:: PDM_isosurface_create

Input mesh definition
"""""""""""""""""""""

Partitioned
~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_n_part_set
.. doxygenfunction:: PDM_isosurface_connectivity_set
.. doxygenfunction:: PDM_isosurface_vtx_coord_set
.. doxygenfunction:: PDM_isosurface_ln_to_gn_set
.. doxygenfunction:: PDM_isosurface_group_set

.. doxygenfunction:: PDM_isosurface_part_mesh_set

.. doxygenfunction:: PDM_isosurface_part_mesh_nodal_set

Block-distributed
~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_dconnectivity_set
.. doxygenfunction:: PDM_isosurface_dvtx_coord_set
.. doxygenfunction:: PDM_isosurface_distrib_set
.. doxygenfunction:: PDM_isosurface_dgroup_set

.. doxygenfunction:: PDM_isosurface_dmesh_set

.. doxygenfunction:: PDM_isosurface_dmesh_nodal_set

Iso-surface settings
""""""""""""""""""""

.. doxygenfunction:: PDM_isosurface_redistribution_set

.. doxygenfunction:: PDM_isosurface_n_part_out_set

.. doxygenfunction:: PDM_isosurface_add

.. doxygenfunction:: PDM_isosurface_isovalues_set

.. doxygenenum:: PDM_iso_surface_kind_t

.. doxygenfunction:: PDM_isosurface_equation_set
.. doxygenfunction:: PDM_isosurface_field_function_set

.. doxygentypedef:: PDM_isosurface_field_function_t

.. doxygenfunction:: PDM_isosurface_set_tolerance


Partitioned discrete field
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_field_set

Block-distributed discrete field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_dfield_set

Iso-surface computation
"""""""""""""""""""""""

.. doxygenfunction:: PDM_isosurface_reset
.. doxygenfunction:: PDM_isosurface_compute
.. doxygenfunction:: PDM_isosurface_dump_times

Outputs
"""""""

.. .. todo::

..    - sortie en part_mesh_nodal/dmesh_nodal?


Partitioned
~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_connectivity_get
.. doxygenfunction:: PDM_isosurface_vtx_coord_get
.. doxygenfunction:: PDM_isosurface_ln_to_gn_get
.. doxygenfunction:: PDM_isosurface_group_get
.. doxygenfunction:: PDM_isosurface_local_parent_get
.. doxygenfunction:: PDM_isosurface_parent_weight_get
.. doxygenfunction:: PDM_isosurface_isovalue_entity_idx_get

Block-distributed
~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_distrib_get
.. doxygenfunction:: PDM_isosurface_dconnectivity_get
.. doxygenfunction:: PDM_isosurface_dvtx_coord_get
.. doxygenfunction:: PDM_isosurface_dgroup_get
.. doxygenfunction:: PDM_isosurface_dparent_weight_get

Communication graphs
~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_part_to_part_enable
.. doxygenfunction:: PDM_isosurface_part_to_part_get


Finalization
""""""""""""

.. doxygenfunction:: PDM_isosurface_free




Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. todo:: TO DO

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: Isosurface

    Python structure to perform isosurface and slice construction. Once initialized, all the following
    methods apply to a :class:`Isosurface` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.Isosurface.__init__

    .. rubric:: Isosurface general inputs

    .. automethod:: Pypdm.Pypdm.Isosurface.tolerance_set
    .. automethod:: Pypdm.Pypdm.Isosurface.add
    .. automethod:: Pypdm.Pypdm.Isosurface.isovalues_set
    .. automethod:: Pypdm.Pypdm.Isosurface.equation_set
    .. automethod:: Pypdm.Pypdm.Isosurface.field_function_set
    .. automethod:: Pypdm.Pypdm.Isosurface.compute
    .. automethod:: Pypdm.Pypdm.Isosurface.reset
    .. automethod:: Pypdm.Pypdm.Isosurface.part_to_part_enable

    .. rubric:: Partitioned source mesh definition

    .. automethod:: Pypdm.Pypdm.Isosurface.mesh_n_part_set
    .. automethod:: Pypdm.Pypdm.Isosurface.n_part_out_set
    .. automethod:: Pypdm.Pypdm.Isosurface.connectivity_set
    .. automethod:: Pypdm.Pypdm.Isosurface.coordinates_set
    .. automethod:: Pypdm.Pypdm.Isosurface.ln_to_gn_set
    .. automethod:: Pypdm.Pypdm.Isosurface.group_set
    .. automethod:: Pypdm.Pypdm.Isosurface.part_mesh_set
    .. automethod:: Pypdm.Pypdm.Isosurface.part_mesh_nodal_set
    .. automethod:: Pypdm.Pypdm.Isosurface.redistribution_set
    .. automethod:: Pypdm.Pypdm.Isosurface.field_set

    .. rubric:: Distributed source mesh definition

    .. automethod:: Pypdm.Pypdm.Isosurface.dconnectivity_set
    .. automethod:: Pypdm.Pypdm.Isosurface.dcoordinates_set
    .. automethod:: Pypdm.Pypdm.Isosurface.distribution_set
    .. automethod:: Pypdm.Pypdm.Isosurface.dgroup_set
    .. automethod:: Pypdm.Pypdm.Isosurface.dmesh_set
    .. automethod:: Pypdm.Pypdm.Isosurface.dmesh_nodal_set
    .. automethod:: Pypdm.Pypdm.Isosurface.dfield_set

    .. rubric:: Partitioned output mesh get

    .. automethod:: Pypdm.Pypdm.Isosurface.connectivity_get
    .. automethod:: Pypdm.Pypdm.Isosurface.coordinates_get
    .. automethod:: Pypdm.Pypdm.Isosurface.ln_to_gn_get
    .. automethod:: Pypdm.Pypdm.Isosurface.group_get
    .. automethod:: Pypdm.Pypdm.Isosurface.parent_lnum_get
    .. automethod:: Pypdm.Pypdm.Isosurface.parent_weight_get
    .. automethod:: Pypdm.Pypdm.Isosurface.isovalue_idx_get
    .. automethod:: Pypdm.Pypdm.Isosurface.parent_weight_get

    .. rubric:: Distributed output mesh get

    .. automethod:: Pypdm.Pypdm.Isosurface.dconnectivity_get
    .. automethod:: Pypdm.Pypdm.Isosurface.dcoordinates_get
    .. automethod:: Pypdm.Pypdm.Isosurface.distribution_get
    .. automethod:: Pypdm.Pypdm.Isosurface.dgroup_get
    .. automethod:: Pypdm.Pypdm.Isosurface.dparent_weight_get
    .. automethod:: Pypdm.Pypdm.Isosurface.disovalue_entity_get

    .. rubric:: General output

    .. automethod:: Pypdm.Pypdm.Isosurface.part_to_part_get
    .. automethod:: Pypdm.Pypdm.Isosurface.dump_times




.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


Annexes
-------

Algorithm description
"""""""""""""""""""""

General wrapping
~~~~~~~~~~~~~~~~

Two core algorithm are implemented here, but at least six type of entries are managed. Implementation tries to simplify all combinations by reducing them to common parts as quickly as possible. The following figure summarizes the process.

.. image:: ../../../../images/isosurface_algo_graph.svg
  :width: 100%
  :align: center

Marching algo
~~~~~~~~~~~~~

This algorithm is a variation of `marching tetrahedras <https://fr.wikipedia.org/wiki/Marching_tetrahedra>`_ based on tetrahedra edges to manage degenerated cases such as when isosurface is passing exactly on an edge or a tri of tetrahedras.

First, for each element (triangle or tetrahedra), its edges are decomposed to see if any of it is crossed by isosurface or on it. Then, iso-vertices are generated for these edges while storing parent vertices information and a weight for a future interpolation. Finally, all volumic elements are traversed to build faces fully on isosurface level.

Once these informations are built, isosurface mesh can be fully generated by going though all elements and, depending of how many edges are crossed by isosurface, we can generate isosurface element according to a configuration table.
During this step, entry parent are preserved for future interpolation and the global ids from entry mesh is used to prepare
the isosurface entities global ids. In 3D, for triangles, boundary information is transfered on generated isosurface edges.

Finally, global ids are computed for groups and isosurface entities.

Ngon algo
~~~~~~~~~
This algorithm is an implementation of `Lopez et al algorithm <https://www.sciencedirect.com/science/article/pii/S0021999121004745>`_.

For each element crossed by isosurface and for each of its face, we go through edges and if it is crossed by isosurface an isosurface vertices is added on edge, linked with previous one to build an edge. At the end of the loop on element face,
the edge built must be closed, so we can define a 2D element. During the loop it is easy to preserve parent, boundary and global ids information for isosurface global ids generation and a future interpolation.