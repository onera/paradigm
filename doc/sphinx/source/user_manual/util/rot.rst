.. _rot:

Rotations in 3D
===============

This module handles rotations (and translations) of rigid bodies in 3D (using
 `quaternions <https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation>`_ ).
C API
-----

Format conversion
"""""""""""""""""

.. doxygenfunction:: PDM_rotation_axis_angle_to_euler_angles
.. doxygenfunction:: PDM_rotation_axis_angle_to_rotation_matrix
.. doxygenfunction:: PDM_rotation_axis_angle_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_euler_angles_to_axis_angle
.. doxygenfunction:: PDM_rotation_euler_angles_to_euler_angles
.. doxygenfunction:: PDM_rotation_euler_angles_to_rotation_matrix
.. doxygenfunction:: PDM_rotation_euler_angles_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_rotation_matrix_to_axis_angle
.. doxygenfunction:: PDM_rotation_rotation_matrix_to_euler_angles
.. doxygenfunction:: PDM_rotation_rotation_matrix_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_axis_angle
.. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_euler_angles
.. doxygenfunction:: PDM_rotation_homogeneous_matrix_to_rotation_matrix
.. doxygenfunction:: PDM_rotation_two_vectors_to_axis_angle
.. doxygenfunction:: PDM_rotation_two_vectors_to_euler_angles
.. doxygenfunction:: PDM_rotation_two_vectors_to_rotation_matrix
.. doxygenfunction:: PDM_rotation_two_vectors_to_homogeneous_matrix

.. doxygenfunction:: PDM_rotation_axis_angle_and_rotation_center_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_euler_angles_and_rotation_center_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_rotation_matrix_and_rotation_center_to_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_periodic_t_info_to_homogeneous_matrix

Application to coordinates
""""""""""""""""""""""""""

.. doxygenfunction:: PDM_rotation_apply_homogeneous_matrix
.. doxygenfunction:: PDM_rotation_compose_homogeneous_matrices
.. doxygenfunction:: PDM_rotation_apply_euler_angles_and_rotation_center
.. doxygenfunction:: PDM_rotation_apply_axis_angle_and_rotation_center
.. doxygenfunction:: PDM_rotation_apply_rotation_matrix_and_rotation_center
.. doxygenfunction:: PDM_rotation_axis_angle_to_euler_angles

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'

  .. autosummary::
    :nosignatures:

    ~Pypdm.Pypdm.axis_angle_to_euler_angles
    ~Pypdm.Pypdm.axis_angle_to_rotation_matrix
    ~Pypdm.Pypdm.euler_angles_to_axis_angle
    ~Pypdm.Pypdm.euler_angles_to_euler_angles
    ~Pypdm.Pypdm.euler_angles_to_rotation_matrix
    ~Pypdm.Pypdm.rotation_matrix_to_axis_angle
    ~Pypdm.Pypdm.rotation_matrix_to_euler_angles
    ~Pypdm.Pypdm.two_vectors_to_axis_angle
    ~Pypdm.Pypdm.two_vectors_to_euler_angles
    ~Pypdm.Pypdm.two_vectors_to_rotation_matrix
    ~Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_coords
    ~Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_vector_field
    ~Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_coords
    ~Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_vector_field
    ~Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_coords
    ~Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_vector_field
    ~Pypdm.Pypdm.axis_angle_and_rotation_center_to_homogeneous_matrix
    ~Pypdm.Pypdm.euler_angles_and_rotation_center_to_homogeneous_matrix
    ~Pypdm.Pypdm.rotation_matrix_and_rotation_center_to_homogeneous_matrix
    ~Pypdm.Pypdm.periodic_t_info_to_homogeneous_matrix

  .. rubric:: Format conversion

  .. autofunction:: Pypdm.Pypdm.axis_angle_to_euler_angles
  .. autofunction:: Pypdm.Pypdm.axis_angle_to_rotation_matrix
  .. autofunction:: Pypdm.Pypdm.euler_angles_to_axis_angle
  .. autofunction:: Pypdm.Pypdm.euler_angles_to_euler_angles
  .. autofunction:: Pypdm.Pypdm.euler_angles_to_rotation_matrix
  .. autofunction:: Pypdm.Pypdm.rotation_matrix_to_axis_angle
  .. autofunction:: Pypdm.Pypdm.rotation_matrix_to_euler_angles
  .. autofunction:: Pypdm.Pypdm.two_vectors_to_axis_angle
  .. autofunction:: Pypdm.Pypdm.two_vectors_to_euler_angles
  .. autofunction:: Pypdm.Pypdm.two_vectors_to_rotation_matrix

  .. autofunction:: Pypdm.Pypdm.axis_angle_and_rotation_center_to_homogeneous_matrix
  .. autofunction:: Pypdm.Pypdm.euler_angles_and_rotation_center_to_homogeneous_matrix
  .. autofunction:: Pypdm.Pypdm.rotation_matrix_and_rotation_center_to_homogeneous_matrix
  .. autofunction:: Pypdm.Pypdm.periodic_t_info_to_homogeneous_matrix

  .. rubric:: Application to coordinates & vector fields

  .. autofunction:: Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_coords
  .. autofunction:: Pypdm.Pypdm.apply_euler_angles_and_rotation_center_to_vector_field
  .. autofunction:: Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_coords
  .. autofunction:: Pypdm.Pypdm.apply_axis_angle_and_rotation_center_to_vector_field
  .. autofunction:: Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_coords
  .. autofunction:: Pypdm.Pypdm.apply_rotation_matrix_and_rotation_center_to_vector_field
  
.. ifconfig:: enable_python_doc == 'OFF'

  .. warning::
    Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


