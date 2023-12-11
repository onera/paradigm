.. _isosurface:

Iso-surfaces
============

C API
-----

Enumerators
"""""""""""

.. doxygenenum:: PDM_iso_surface_kind_t

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

.. doxygenfunction:: PDM_isosurface_mesh_nodal_set

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

.. doxygenfunction:: PDM_isosurface_add
.. doxygenfunction:: PDM_isosurface_equation_set
.. doxygenfunction:: PDM_isosurface_field_function_set

.. todo::

  - enable/disable ptps?
  - elt_type, extract_kind, part_method? (create ou add?)
  - PDM_isosurface_field_gradient_function_set?

Partitioned discrete field
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_field_set
.. doxygenfunction:: PDM_isosurface_gradient_set

Block-distributed discrete field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_field_set
.. doxygenfunction:: PDM_isosurface_gradient_set

Iso-surface computation
"""""""""""""""""""""""

.. doxygenfunction:: PDM_isosurface_reset
.. doxygenfunction:: PDM_isosurface_compute
.. doxygenfunction:: PDM_isosurface_dump_times

Outputs
"""""""

Partitioned
~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_connectivity_get
.. doxygenfunction:: PDM_isosurface_vtx_coord_get
.. doxygenfunction:: PDM_isosurface_ln_to_gn_get

Block-distributed
~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_dconnectivity_get
.. doxygenfunction:: PDM_isosurface_dvtx_coord_get

Communication graphs
~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_isosurface_part_to_part_get


Finalization
""""""""""""

.. doxygenfunction:: PDM_isosurface_free
