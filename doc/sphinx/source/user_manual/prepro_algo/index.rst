.. _prepro_algo:

###################################
Pre-/co-/post-processing algorithms
###################################



Point cloud localization inside a mesh
======================================

C API
-----

.. doxygenfile:: pdm_mesh_location.h

Fortran API
-----------
.. f:automodule:: pdm_mesh_location

Python API
----------

.. autoclass:: Pypdm.Pypdm.MeshLocation
  :members:


Other
=====
* nearest neighbors search
* distance from surface
* inside cloud surf
* mesh intersection
* extraction of mesh partitions
* iso-surfaces & slices
* overlay
* mesh adaptation/remeshing
