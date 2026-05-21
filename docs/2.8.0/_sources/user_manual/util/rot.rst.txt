.. _rotation:

Rotations in 3D
===============

This module handles rotations (and translations) of rigid bodies in 3D (using `quaternions <https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation>`_).

API
"""

.. dropdown:: Format conversion

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_rotation_axis_angle_to_euler_angles
      .. doxygenfunction:: PDM_rotation_axis_angle_to_rotation_matrix
      .. doxygenfunction:: PDM_rotation_axis_angle_to_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix

      .. doxygenfunction:: PDM_rotation_euler_angles_to_axis_angle
      .. doxygenfunction:: PDM_rotation_euler_angles_to_euler_angles
      .. doxygenfunction:: PDM_rotation_euler_angles_to_rotation_matrix
      .. doxygenfunction:: PDM_rotation_euler_angles_to_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_periodic_t_info_to_homogeneous_matrix

      .. doxygenfunction:: PDM_rotation_rotation_matrix_to_axis_angle
      .. doxygenfunction:: PDM_rotation_rotation_matrix_to_euler_angles
      .. doxygenfunction:: PDM_rotation_rotation_matrix_to_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix

      .. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_axis_angle
      .. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_euler_angles
      .. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_rotation_matrix

      .. doxygenfunction:: PDM_rotation_two_vectors_to_axis_angle
      .. doxygenfunction:: PDM_rotation_two_vectors_to_euler_angles
      .. doxygenfunction:: PDM_rotation_two_vectors_to_rotation_matrix
      .. doxygenfunction:: PDM_rotation_two_vectors_to_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_two_vectors_and_rotation_center_to_homogeneous_matrix

      .. doxygenfunction:: PDM_rotation_axes_and_origin_to_homogeneous_matrix


    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.axis_angle_to_euler_angles
        .. autofunction:: Pypdm.Pypdm.axis_angle_to_rotation_matrix
        .. autofunction:: Pypdm.Pypdm.axis_angle_to_homogeneous_matrix
        .. autofunction:: Pypdm.Pypdm.axis_angle_and_rotation_center_to_homogeneous_matrix

        .. autofunction:: Pypdm.Pypdm.euler_angles_to_axis_angle
        .. autofunction:: Pypdm.Pypdm.euler_angles_to_euler_angles
        .. autofunction:: Pypdm.Pypdm.euler_angles_to_rotation_matrix
        .. autofunction:: Pypdm.Pypdm.euler_angles_to_homogeneous_matrix
        .. autofunction:: Pypdm.Pypdm.euler_angles_and_rotation_center_to_homogeneous_matrix

        .. autofunction:: Pypdm.Pypdm.rotation_matrix_to_axis_angle
        .. autofunction:: Pypdm.Pypdm.rotation_matrix_to_euler_angles
        .. autofunction:: Pypdm.Pypdm.rotation_matrix_to_homogeneous_matrix
        .. autofunction:: Pypdm.Pypdm.rotation_matrix_and_rotation_center_to_homogeneous_matrix

        .. autofunction:: Pypdm.Pypdm.two_vectors_to_axis_angle
        .. autofunction:: Pypdm.Pypdm.two_vectors_to_euler_angles
        .. autofunction:: Pypdm.Pypdm.two_vectors_to_rotation_matrix
        .. autofunction:: Pypdm.Pypdm.two_vectors_to_homogeneous_matrix
        .. autofunction:: Pypdm.Pypdm.two_vectors_and_rotation_center_to_homogeneous_matrix

        .. autofunction:: Pypdm.Pypdm.axes_and_origin_to_homogeneous_matrix


      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



.. dropdown:: Application to coordinates & vector fields

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_rotation_apply_homogeneous_matrix
      .. doxygenfunction:: PDM_rotation_compose_homogeneous_matrices
      .. doxygenfunction:: PDM_rotation_apply_euler_angles_and_rotation_center
      .. doxygenfunction:: PDM_rotation_apply_axis_angle_and_rotation_center
      .. doxygenfunction:: PDM_rotation_apply_rotation_matrix_and_rotation_center
      .. doxygenfunction:: PDM_rotation_axis_angle_to_euler_angles


    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_coords
        .. autofunction:: Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_vector_field
        .. autofunction:: Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_coords
        .. autofunction:: Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_vector_field
        .. autofunction:: Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_coords
        .. autofunction:: Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_vector_field
        .. autofunction:: Pypdm.Pypdm.apply_homogeneous_matrix_to_coords

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)
