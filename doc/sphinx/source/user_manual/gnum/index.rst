.. _gnum:


################
Global numbering
################


.. _gen_gnum:

Global numbering generation
===========================

**ParaDiGM** relies heavily on the notion of :ref:`global numbering <concept_global_id>`.
If your code does not use global IDs, these can be generated from geometric data.
Such global numbering is achieved by encoding Cartesian coordinates along the `Morton space-filling curve <https://en.wikipedia.org/wiki/Z-order_curve>`_.


API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_gnum_create



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_gnum_create

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.__init__
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set inputs

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_gnum_set_from_coords
      .. doxygenfunction:: PDM_gnum_set_from_parents
      .. doxygenfunction:: PDM_gnum_set_parents_nuplet



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: PDM_gnum_set_from_coords
          .. f:autosubroutine:: PDM_gnum_set_from_parents
          .. f:autosubroutine:: PDM_gnum_set_parents_nuplet

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.set_from_coords
          :noindex:
        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.set_from_parent
          :noindex:
        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.set_parents_nuplet
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Build global numbering

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_gnum_compute



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:subroutine:: pdm_gnum_compute(gen_gnum)

          Build global numbering

          :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.compute
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get outputs

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_gnum_get



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: PDM_gnum_get


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.get
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_gnum_free



    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:subroutine:: pdm_gnum_free(gen_gnum)

          Free the Global Numbering Generation object

          :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      |python_gc|





Examples
""""""""

.. dropdown:: From coordinates

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      (Extract from ``test/pdm_t_gen_gnum.c``)

      .. literalinclude:: ../../../../../test/pdm_t_gen_gnum.c
        :name: c_gen_gnum_ex
        :language: c
        :lines: 6-7, 158-186


    .. tab-item:: Fortran
      :sync: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        (Extract from ``test/pdm_t_gnum_f.F90``)

        .. literalinclude:: ../../../../../test/pdm_t_gnum_f.F90
          :name: fortran_gen_gnum_ex
          :language: fortran
          :dedent: 2
          :lines: 25, 29, 42-46, 70-99


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python
      :sync: Python

      .. ifconfig:: enable_python_doc == 'ON'

        (Extract from ``test/pdm_t_gnum_p.py``)

        .. literalinclude:: ../../../../../test/pdm_t_gnum_p.py
          :name: python_gen_gnum_ex
          :language: python
          :dedent: 2
          :lines: 63-83

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)