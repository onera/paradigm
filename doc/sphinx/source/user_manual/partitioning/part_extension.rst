.. _part_extension:

Part-Extension
==============

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_extend_type_t

Initialization
~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_create

Set inputs
~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_connectivity_set

.. doxygenfunction:: PDM_part_extension_vtx_coord_set

.. doxygenfunction:: PDM_part_extension_ln_to_gn_set

.. doxygenfunction:: PDM_part_extension_part_bound_graph_set

.. doxygenfunction:: PDM_part_extension_group_set

.. .. doxygenfunction:: PDM_part_extension_set_part

.. .. doxygenfunction:: PDM_part_extension_part_domain_interface_shared_set

Perform exchange of extended partition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_compute

Get outputs
~~~~~~~~~~~

.. doxygenfunction:: PDM_part_extension_connectivity_get

.. doxygenfunction:: PDM_part_extension_ln_to_gn_get

.. doxygenfunction:: PDM_part_extension_vtx_coord_get

.. doxygenfunction:: PDM_part_extension_group_get

.. .. doxygenfunction:: PDM_part_extension_interface_get

.. .. doxygenfunction:: PDM_part_extension_composed_interface_get


Finalize
~~~~~~~~

.. doxygenfunction:: PDM_part_extension_free

Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_create

  Set inputs
  ~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_set_part

  Perform exchange of extended partition
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_compute

  Get outputs
  ~~~~~~~~~~~

  .. f:autosubroutine PDM_part_extension_connectivity_get

  .. f:autosubroutine PDM_part_extension_ln_to_gn_get

  .. f:autosubroutine PDM_part_extension_vtx_coord_get

  .. f:autosubroutine PDM_part_extension_group_get

  Finalize
  ~~~~~~~~

  .. f:autosubroutine PDM_part_extension_free

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  ~~~~~~~~~~~~~~

  .. autoclass:: Pypdm.Pypdm.PartExtension

  Set inputs
  ~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.PartExtension.connectivity_set

  .. autofunction:: Pypdm.Pypdm.PartExtension.vtx_coord_set

  .. autofunction:: Pypdm.Pypdm.PartExtension.ln_to_gn_set

  .. autofunction:: Pypdm.Pypdm.PartExtension.part_bound_graph_set

  .. autofunction:: Pypdm.Pypdm.PartExtension.group_set

  .. .. autofunction:: Pypdm.Pypdm.PartExtension.set_part

  .. .. autofunction:: Pypdm.Pypdm.PartExtension.part_domain_interface_shared_set

  Perform exchange of extended partition
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.PartExtension.compute

  Get outputs
  ~~~~~~~~~~~

  .. autofunction:: Pypdm.Pypdm.PartExtension.connectivity_get

  .. autofunction:: Pypdm.Pypdm.PartExtension.vtx_coord_get

  .. autofunction:: Pypdm.Pypdm.PartExtension.ln_to_gn_get

  .. autofunction:: Pypdm.Pypdm.PartExtension.group_get

  .. .. autofunction:: Pypdm.Pypdm.PartExtension.get_interface

  .. .. autofunction:: Pypdm.Pypdm.PartExtension.get_composed_interface

.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
