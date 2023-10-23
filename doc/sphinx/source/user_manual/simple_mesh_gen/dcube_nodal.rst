.. _dcube_nodal:

Distributed nodal square/cube mesh
==================================


C API
-----

Initialization
""""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_create

Options
"""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_random_factor_set

.. doxygenfunction:: PDM_dcube_nodal_gen_ordering_set


Mesh generation
"""""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_build

.. doxygenfunction:: PDM_dcube_nodal_gen_dmesh_nodal_get


Multi-zone mesh generation
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: PDM_dcube_nodal_cart_topo


Finalization
""""""""""""

.. doxygenfunction:: PDM_dcube_nodal_gen_free




Fortran API
-----------

.. ifconfig:: enable_fortran_doc == 'ON'

  .. todo::
    ...

.. ifconfig:: enable_fortran_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)




Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  Initialization
  """"""""""""""

  .. autoclass:: Pypdm.Pypdm.DCubeNodalGenerator

  Options
  """""""

  .. autofunction:: Pypdm.Pypdm.DCubeNodalGenerator.set_random_factor

  .. autofunction:: Pypdm.Pypdm.DCubeNodalGenerator.set_ordering


  Mesh generation
  """""""""""""""

  .. autofunction:: Pypdm.Pypdm.DCubeNodalGenerator.compute

  .. autofunction:: Pypdm.Pypdm.DCubeNodalGenerator.get_dmesh_nodal


.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)

