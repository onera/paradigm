.. _mesh_location:

Mesh location
=============

**Mesh location** is a service for locating points inside a partitioned, unstructured mesh.
All types of mesh elements are supported: standard elements (bars, triangles,
quadrangles, tetrahedra, pyramids, prisms, hexahedra), general polygons and polyhedra,
as well as high-order, curved elements.

A mapping between the source mesh elements and the target points they contain is computed, which consists in
  - geometric data (distances, barycentric and parametric coordinates, ...) ;
  - a :ref:`PDM_part_to_part <ptp>` instance to transfer data in parallel.

This mapping can be used to performed spatial interpolation, as in `CWIPI <https://github.com/onera/cwipi>`_ and `Maia <https://github.com/onera/maia>`_.

API
"""

.. dropdown:: Initialization


  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_create



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.__init__
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set target point clouds

  Multiple target point clouds can be processed with a single Mesh Location instance.
  Each of them can have a distinct number of parts per MPI rank, and has its own global numbering.

  A point cloud is defined by coordinates and global IDs (if you don't have a global numbering, check out :ref:`this page <gnum>` to see how to generate one).


  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_n_part_cloud_set
      .. doxygenfunction:: PDM_mesh_location_cloud_set



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_n_part_cloud_set
        .. f:autosubroutine:: PDM_mesh_location_cloud_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.n_part_cloud_set
          :noindex:
        .. autofunction:: Pypdm.Pypdm.MeshLocation.cloud_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set source mesh

  The source mesh can be defined in several ways.

  If you have a :ref:`Part Mesh Nodal <pmn>` instance, you can use it directly
  (note that this is the only way to define high-order, curved elements):

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_shared_nodal_mesh_set

    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_shared_nodal_mesh_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)


    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        Not yet available

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  Alternatively, you can define the source mesh by providing each part separately.
  If you do so, you should always start by setting the number of part:

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_mesh_n_part_set



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_mesh_n_part_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.mesh_n_part_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  Then you can use the appropriate setter depending on the mesh dimension and type:

  .. dropdown:: Volume mesh

    .. dropdown:: Nodal connectivity

      .. tab-set::

          .. tab-item:: C

            .. doxygenfunction:: PDM_mesh_location_nodal_part_set



          .. tab-item:: Fortran

            .. ifconfig:: enable_fortran_doc == 'ON'

              .. f:autosubroutine:: PDM_mesh_location_nodal_part_set

            .. ifconfig:: enable_fortran_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



          .. tab-item:: Python

            .. ifconfig:: enable_python_doc == 'ON'

              .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set
                :noindex:

            .. ifconfig:: enable_python_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



    .. dropdown:: Ngon connectivity

      .. tab-set::

          .. tab-item:: C

            .. doxygenfunction:: PDM_mesh_location_part_set



          .. tab-item:: Fortran

            .. ifconfig:: enable_fortran_doc == 'ON'

              .. f:autosubroutine:: PDM_mesh_location_part_set

            .. ifconfig:: enable_fortran_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



          .. tab-item:: Python

            .. ifconfig:: enable_python_doc == 'ON'

              .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set
                :noindex:

            .. ifconfig:: enable_python_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



  .. dropdown:: Surface mesh

    .. dropdown:: Nodal connectivity

      .. tab-set::

          .. tab-item:: C

            .. doxygenfunction:: PDM_mesh_location_nodal_part_set_2d



          .. tab-item:: Fortran

            .. ifconfig:: enable_fortran_doc == 'ON'

              .. f:autosubroutine:: PDM_mesh_location_nodal_part_set_2d

            .. ifconfig:: enable_fortran_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



          .. tab-item:: Python

            .. ifconfig:: enable_python_doc == 'ON'

              .. autofunction:: Pypdm.Pypdm.MeshLocation.nodal_part_set_2d
                :noindex:

            .. ifconfig:: enable_python_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



    .. dropdown:: Ngon connectivity

      .. tab-set::

          .. tab-item:: C

            .. doxygenfunction:: PDM_mesh_location_part_set_2d



          .. tab-item:: Fortran

            .. ifconfig:: enable_fortran_doc == 'ON'

              .. f:autosubroutine:: PDM_mesh_location_part_set_2d

            .. ifconfig:: enable_fortran_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



          .. tab-item:: Python

            .. ifconfig:: enable_python_doc == 'ON'

              .. autofunction:: Pypdm.Pypdm.MeshLocation.part_set_2d
                :noindex:

            .. ifconfig:: enable_python_doc == 'OFF'

              .. warning::
                Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: *Optional* parameters

  .. dropdown:: Geometric tolerance

    The location algorithms rely on bounding-box tests to quickly find candidate pairs of points and elements, before computing the exact location.
    A target point will be considered as *located* if it lies in the bounding box of at least one source mesh element (even if it does not actually lie inside the element).
    These bounding boxes can be expanded using a relative tolerance.

    By default, the tolerance is set to zero, but it can be adjusted.
    This is especially useful for non-planar, surface meshes where alignment with cartesian axes might cause some detection misses if the tolerance is set too low.
    However, keep in mind that setting a very large tolerance will have a significant impact on performance.

    We recommend keeping the tolerance between 0 and 0.1.

    .. tab-set::

      .. tab-item:: C

        .. doxygenfunction:: PDM_mesh_location_tolerance_set



      .. tab-item:: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_mesh_location_tolerance_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autoattribute:: Pypdm.Pypdm.MeshLocation.tolerance
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  .. dropdown:: Preconditioning method

    Experienced users can also choose the preconditioning method used in the first step of the location algorithm:

    .. tab-set::

      .. tab-item:: C

        .. doxygenfunction:: PDM_mesh_location_method_set



      .. tab-item:: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_mesh_location_method_set

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autoattribute:: Pypdm.Pypdm.MeshLocation.method
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Compute location

  Once the target point clouds and source mesh have been set, the location can be computed:

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_compute



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_compute

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.compute
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  Once the calculation is complete, you can optionally display the elapsed and CPU times:

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_dump_times



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_dump_times

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.dump_times
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get location data

  Once computed, the location data can retrieved either from the perspective of the target point cloud or from the perspective of the source mesh.

  .. dropdown:: Target perspective

    First, the number and IDs of located/unlocated points in each part:

    .. tab-set::

      .. tab-item:: C

        .. doxygenfunction:: PDM_mesh_location_n_located_get
        .. doxygenfunction:: PDM_mesh_location_located_get

        .. doxygenfunction:: PDM_mesh_location_n_unlocated_get
        .. doxygenfunction:: PDM_mesh_location_unlocated_get



      .. tab-item:: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_mesh_location_n_located_get
          .. f:autosubroutine:: PDM_mesh_location_located_get

          .. f:autosubroutine:: PDM_mesh_location_n_unlocated_get
          .. f:autosubroutine:: PDM_mesh_location_unlocated_get

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.MeshLocation.located_get
            :noindex:
          .. autofunction:: Pypdm.Pypdm.MeshLocation.unlocated_get
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


    Second, the location data of located points:

    .. tab-set::

      .. tab-item:: C

        .. doxygenfunction:: PDM_mesh_location_point_location_get



      .. tab-item:: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_mesh_location_point_location_get

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.MeshLocation.location_get
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)



  .. dropdown:: Source perspective

    .. tab-set::

      .. tab-item:: C

        .. doxygenfunction:: PDM_mesh_location_points_in_elt_get

        .. doxygenfunction:: PDM_mesh_location_cell_vertex_get



      .. tab-item:: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_mesh_location_points_in_elt_get

          .. f:subroutine:: pdm_mesh_location_cell_vertex_get(mloc, i_part, cell_vtx_idx, cell_vtx)

            Get the cell→vertex connectivity used for internal computations

            .. note::
              For non-standard elements, this connectivity is built by ParaDiGM and is necessary to associate
              the ``points_weights`` array (returned by **pdm_mesh_location_points_in_elt_get**)
              to the appropriate mesh vertices.

            :p c_ptr mesh_loc[in]:           Mesh location instance
            :p integer i_part[in]:           Partition identifier
            :p integer(:) cell_vtx_idx[out]: Index for cell → vertex connectivity
            :p integer(:) cell_vtx[out]:     Cell → vertex connectivity

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.MeshLocation.points_in_elt_get
            :noindex:
          .. autofunction:: Pypdm.Pypdm.MeshLocation.cell_vertex_get
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Data transfer

  Data can be transferred between the source mesh and each point cloud using the :ref:`PDM_part_to_part <ptp>` instance created by **Mesh location**.

  .. note:: *Direct* exchanges go from source to target and *reverse* exchanges go from target to source.

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_part_to_part_get



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_part_to_part_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.MeshLocation.part_to_part_get
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_mesh_location_free



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_mesh_location_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      |python_gc|




