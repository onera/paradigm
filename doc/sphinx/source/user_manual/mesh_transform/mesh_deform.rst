.. _mesh_deform:

Mesh Deformation
================

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_clustering_kind_t

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_deform_create


Set inputs
~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_deform_surf_part_set

.. doxygenfunction:: PDM_mesh_deform_cloud_part_set

Perform deformation
~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_deform_compute

.. doxygenfunction:: PDM_mesh_deform_cloud_displ_finalize

Get outputs
~~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_deform_n_aux_geom_get

.. doxygenfunction:: PDM_mesh_deform_cloud_block_n_points_get

.. doxygenfunction:: PDM_mesh_deform_cloud_block_buffer_size_get

.. doxygenfunction:: PDM_mesh_deform_cloud_block_buffer_from_surf_get

.. doxygenfunction:: PDM_mesh_deform_cloud_block_get

.. doxygenfunction:: PDM_mesh_deform_cloud_n_points_part_get

.. doxygenfunction:: PDM_mesh_deform_cloud_displ_part_get

Finalize
~~~~~~~~

.. doxygenfunction:: PDM_mesh_deform_partial_free

.. doxygenfunction:: PDM_mesh_deform_free

Timers
~~~~~~

.. doxygenfunction:: PDM_mesh_deform_dump_times

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_create

  Set inputs
  ~~~~~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_surf_part_set

  .. f:autosubroutine:: PDM_mesh_deform_cloud_part_set

  Perform deformation
  ~~~~~~~~~~~~~~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_compute

  .. f:autosubroutine:: PDM_mesh_deform_cloud_displ_finalize

  Get outputs
  ~~~~~~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_cloud_block_buffer_from_surf_get

  .. f:autosubroutine:: PDM_mesh_deform_cloud_block_get

  .. f:autosubroutine:: PDM_mesh_deform_cloud_n_points_part_get

  .. f:autosubroutine:: PDM_mesh_deform_cloud_displ_part_get

  .. f:autosubroutine:: PDM_mesh_deform_n_aux_geom_get

  .. f:autosubroutine:: PDM_mesh_deform_cloud_block_n_points_get

  .. f:autosubroutine:: PDM_mesh_deform_cloud_block_buffer_size_get

  Finalize
  ~~~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_partial_free

  .. f:autosubroutine:: PDM_mesh_deform_free

  Timers
  ~~~~~~

  .. f:autosubroutine:: PDM_mesh_deform_dump_times

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. py:class:: MeshDeform

    Python object to perform mesh deformation.
    Once initialized, all the following
    methods apply to a :class:`MeshDeform` instance.

    .. rubric:: Initialization

    .. autofunction:: Pypdm.Pypdm.MeshDeform.__cinit__

    .. rubric:: Methods summary

    .. autosummary::
      :nosignatures:

      ~Pypdm.Pypdm.MeshDeform.surf_part_set
      ~Pypdm.Pypdm.MeshDeform.cloud_part_set
      ~Pypdm.Pypdm.MeshDeform.compute
      ~Pypdm.Pypdm.MeshDeform.cloud_displ_finalize
      ~Pypdm.Pypdm.MeshDeform.cloud_n_points_part_get
      ~Pypdm.Pypdm.MeshDeform.cloud_displ_part_get
      ~Pypdm.Pypdm.MeshDeform.n_aux_geom
      ~Pypdm.Pypdm.MeshDeform.cloud_block_n_points
      ~Pypdm.Pypdm.MeshDeform.cloud_block_buffer_from_surf_get
      ~Pypdm.Pypdm.MeshDeform.cloud_block_get
      ~Pypdm.Pypdm.MeshDeform.dump_times

    .. rubric:: Set inputs

    .. automethod:: Pypdm.Pypdm.MeshDeform.surf_part_set
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_part_set

    .. rubric:: Perform deformation

    .. automethod:: Pypdm.Pypdm.MeshDeform.compute
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_displ_finalize

    .. rubric:: Get outputs

    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_n_points_part_get
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_displ_part_get
    .. automethod:: Pypdm.Pypdm.MeshDeform.n_aux_geom
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_block_n_points
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_block_buffer_from_surf_get
    .. automethod:: Pypdm.Pypdm.MeshDeform.cloud_block_get

    .. rubric:: Timers

    .. automethod:: Pypdm.Pypdm.MeshDeform.dump_times
