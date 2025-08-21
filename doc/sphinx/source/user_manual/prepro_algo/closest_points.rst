.. _closest_points:

Closest points
==============

Description
"""""""""""

**Closest points** is a service for finding nearest neighbors between a source and a target point cloud.

A mapping between the associated points is computed, which consists in
  - the global IDs of closest points ;
  - the distance from each target point to its closest source points ;
  - a :ref:`PDM_part_to_part <ptp>` instance to transfer data in parallel.

This mapping can be used to performed spatial interpolation, as in `CWIPI <https://github.com/onera/cwipi>`_ and `Maia <https://github.com/onera/maia>`_.



API
"""

.. dropdown:: Initialization


  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.__init__
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)





.. dropdown:: Set point clouds

  Both the target and source point clouds can have a distinct number of parts per MPI rank, and have their own global numbering (which *must* range from 1 to the total number of points in each cloud).

  A point cloud is defined by coordinates and global IDs (if you don't have a global numbering, check out :ref:`this page <gnum>` to see how to generate one).

  First, the number of partitions of each point cloud must be defined:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_n_part_cloud_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_n_part_cloud_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.n_part_cloud_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  Then, each point cloud is defined one partition at a time:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_src_cloud_set
      .. doxygenfunction:: PDM_closest_points_tgt_cloud_set



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_src_cloud_set
        .. f:autosubroutine:: PDM_closest_points_tgt_cloud_set

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.src_cloud_set
          :noindex:
        .. autofunction:: Pypdm.Pypdm.ClosestPoints.tgt_cloud_set
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Find closest points

  Once the source and target point clouds have been set, the correspondence between them can be computed:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_compute



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_compute

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.compute
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)


  Once the calculation is complete, you can optionally display the elapsed and CPU times:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_dump_times



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_dump_times

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.dump_times
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get results

  Once computed, the correspondence data can retrieved either from the perspective of the target point cloud or from the perspective of the source point cloud.


  .. dropdown:: Target perspective

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_closest_points_get



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_closest_points_get

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.ClosestPoints.points_get
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




  .. dropdown:: Source perspective

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_closest_points_tgt_in_src_get
        .. doxygenfunction:: PDM_closest_points_tgt_in_src_dist_get



      .. tab-item:: Fortran
        :sync: Fortran

        .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_closest_points_tgt_in_src_get
          .. f:autosubroutine:: PDM_closest_points_tgt_in_src_dist_get

        .. ifconfig:: enable_fortran_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



      .. tab-item:: Python
        :sync: Python

        .. ifconfig:: enable_python_doc == 'ON'

          .. autofunction:: Pypdm.Pypdm.ClosestPoints.tgt_in_src_get
            :noindex:

        .. ifconfig:: enable_python_doc == 'OFF'

          .. warning::
            Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Data transfer

  Data can be transferred between the source mesh and each point cloud using the :ref:`PDM_part_to_part <ptp>` instance created by **Closest points**.

  .. note:: *Direct* exchanges go from source to target and *reverse* exchanges go from target to source.

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_part_to_part_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_part_to_part_get

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.ClosestPoints.part_to_part_get
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_closest_points_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_closest_points_free

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|

