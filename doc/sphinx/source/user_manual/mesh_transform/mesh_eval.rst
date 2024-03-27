.. _mesh_eval:

Mesh Evaluation
===============

C API
-----

Enumerators
~~~~~~~~~~~

.. doxygenenum:: PDM_mesh_eval_quality_t

Metric-Based
~~~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_eval_edge_length

.. doxygenfunction:: PDM_mesh_eval_elt_edge_ratio

.. doxygenfunction:: PDM_mesh_eval_triangle_quality

.. doxygenfunction:: PDM_mesh_eval_tetrahedron_quality

.. doxygenfunction:: PDM_mesh_eval_elt_anisotropy_ratio

.. doxygenfunction:: PDM_mesh_eval_triangle_skewness

.. doxygenfunction:: PDM_mesh_eval_tetrahedron_skewness

.. doxygenfunction:: PDM_mesh_eval_triangle_radii_ratio

.. doxygenfunction:: PDM_mesh_eval_dmesh_complexity

.. doxygenfunction:: PDM_mesh_eval_pmesh_complexity

.. doxygenfunction:: PDM_mesh_eval_elt_natural_metric

.. doxygenfunction:: PDM_mesh_eval_triangle_load

.. doxygenfunction:: PDM_mesh_eval_tetrahedron_load

Topological
~~~~~~~~~~~

.. doxygenfunction:: PDM_mesh_eval_pvtx_valency

Python API
----------

.. ifconfig:: enable_python_doc == 'ON'
