.. _connec_transform:

Connectivity transformation
===========================


API
"""


.. dropdown:: Partitioned connectivity transformations

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_combine_connectivity
      .. doxygenfunction:: PDM_connectivity_transpose
      .. doxygenfunction:: PDM_compute_face_vtx_from_face_and_edge
      .. doxygenfunction:: PDM_compute_face_vtx_from_face_and_edge_unsigned
      .. doxygenfunction:: PDM_part_connectivity_to_connectivity_idx



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_combine_connectivity
        .. f:autosubroutine:: PDM_connectivity_transpose
        .. f:autosubroutine:: PDM_compute_face_vtx_from_face_and_edge

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.combine_connectivity
        .. autofunction:: Pypdm.Pypdm.connectivity_transpose
        .. autofunction:: Pypdm.Pypdm.compute_face_vtx_from_face_and_edge
        .. autofunction:: Pypdm.Pypdm.part_connectivity_to_connectivity_idx

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Block-distributed connectivity transformations

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_deduce_combine_connectivity
      .. doxygenfunction:: PDM_dconnectivity_transpose
      .. doxygenfunction:: PDM_dconnectivity_dface_vtx_from_face_and_edge
      .. doxygenfunction:: PDM_dfacecell_to_dcellface
      .. doxygenfunction:: PDM_dcellface_to_dfacecell


    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.dconnectivity_combine
        .. autofunction:: Pypdm.Pypdm.dconnectivity_transpose
        .. autofunction:: Pypdm.Pypdm.dfacecell_to_dcellface
        .. autofunction:: Pypdm.Pypdm.dcellface_to_dfacecell

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
