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

    .. tab-item:: C

      .. doxygenfunction:: PDM_gnum_create



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_gnum_create_

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.__init__
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Set inputs

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_gnum_set_from_coords
      .. doxygenfunction:: PDM_gnum_set_from_parents
      .. doxygenfunction:: PDM_gnum_set_parents_nuplet



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

          .. f:autosubroutine:: pdm_gnum_set_from_coords_
          .. f:autosubroutine:: pdm_gnum_set_from_parents_
          .. f:autosubroutine:: pdm_gnum_set_parents_nuplet_

      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

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

    .. tab-item:: C

      .. doxygenfunction:: PDM_gnum_compute



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:subroutine:: pdm_gnum_compute(gen_gnum)

          Build global numbering

          :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.compute
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Get outputs

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_gnum_get



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:autosubroutine:: pdm_gnum_get_


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      .. ifconfig:: enable_python_doc == 'ON'

        .. autofunction:: Pypdm.Pypdm.GlobalNumbering.get
          :noindex:

      .. ifconfig:: enable_python_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_python_interface>` to enable the Python API)




.. dropdown:: Finalization

  .. tab-set::

    .. tab-item:: C

      .. doxygenfunction:: PDM_gnum_free



    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        .. f:subroutine:: pdm_gnum_free(gen_gnum)

          Free the Global Numbering Generation object

          :p c_ptr gen_gnum [in]:  C pointer to PDM_gen_gnum_t object


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

      The instance is automatically freed by Python's garbage collector when it is no longer referenced





Examples
""""""""

.. dropdown:: From coordinates

  .. tab-set::

    .. tab-item:: C

      (Extract from ``test/pdm_t_gen_gnum.c``)

      .. literalinclude:: ../../../../../test/pdm_t_gen_gnum.c
        :name: c_gen_gnum_ex
        :language: c
        :dedent: 2
        :lines: 159-186


    .. tab-item:: Fortran

      .. ifconfig:: enable_fortran_doc == 'ON'

        (Extract from ``test/pdm_t_gnum_f.F90``)

        .. literalinclude:: ../../../../../test/pdm_t_gnum_f.F90
          :name: fortran_gen_gnum_ex
          :language: fortran
          :dedent: 2
          :lines: 43-46, 70-99


      .. ifconfig:: enable_fortran_doc == 'OFF'

        .. warning::
          Unavailable (refer to the :ref:`installation guide <enable_fortran_interface>` to enable the Fortran API)



    .. tab-item:: Python

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