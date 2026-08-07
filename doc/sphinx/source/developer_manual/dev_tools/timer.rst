.. _timer:

Timer
=====


Description
"""""""""""


Usage
"""""


API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_timer_create


.. dropdown:: Usage

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_timer_start
      .. doxygenfunction:: PDM_timer_end

.. dropdown:: Post-treatment

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_timer_print
      .. doxygenfunction:: PDM_timer_log
      .. doxygenfunction:: PDM_timer_dump_json
      .. doxygenfunction:: PDM_timer_gather_dump
      .. doxygenfunction:: PDM_timer_gather_dump_json
      .. doxygenfunction:: PDM_timer_get_report_string

.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_timer_free
