.. _exch_helper:

Exchange helper
===============

Description
"""""""""""

**Exchange helper** is an intermediate structure that helps developers managing MPI pointwise data exchanges.
This structure can be used to build high-level :ref:`exchange protocols <comm_graph>` or to setup raw exchanges, *eg*
in CFD codes.

Schematicaly, the following layers are used in ParaDiGM to perform pointwise data exchanges (from high-level to low-level):

1. The exchange protocol structure is in charge of establishing the *communication graph*, which is the list of connected
   processes and the number of data to send (and receive) to (and from) each connected process.
   This graph is created from input data such as global ids or (rank, part, local_id) adresses.
   The exchange protocol structure also reorders input data into sending and receiving buffers, to sort it according to destination process,
   and do the opposite operation once data is received.
2. The **exchange helper** structure provides unified APIs that allow to switch between several MPI exchange modes:
   blocking or non-blocking, point-to-point or collective or one sided, persistant or oneshot, etc.
   It operates once data have been prepared (reordered) by the exchange protocol, using the computed *communication graph*.
3. Wrapping functions provided by :file:`pdm_mpi_extended.h` file, for exemple ``PDM_MPI_Isends``, which begin N non-blocking
   sends from the current process to its connected N target processes. Theses functions are used by the exchange helper to
   shorten the implementation.
4. Lastly, raw MPI primitives, for exemple ``MPI_Isend``, which begins a non-blocking send from the current process
   to a specified target process, are called from :file:`pdm_mpi_extended.h` wrapping functions.


Exchange modes
""""""""""""""

Since the exchange helper allows to switch between several MPI primitives, this section briefly recalls the caracteristics
of the different modes. Please refer to `MPI documentation <https://www.mpi-forum.org/docs/>`_ for more details.

- **Non-blocking communications** : while blocking communications stop the program until data is fully received,
  non-blocking communications immediatly return a ``request`` object, allowing the program to continue meanwhile
  exchange is performed in the background. The developer is responsible to check if the exchange is completed,
  and eventually to wait for it when the data is actually needed.
- **Persistent communications** : this mode allows to reuse the same "commmunication channel" more than once, when the
  exchanged metadata (target process, datatype, datasize, etc.) remains the same. It typically involves a call to a ``MPI_*_init``
  function, that creates the "communication channel", and then several calls to the ``MPI_Start`` function which triggers
  the exchange with current buffer content.

  .. note:: Persistent communications are always non-blocking

- **Collective communications** : these communications involve all the processes in the current communicator
  and must be called collectively by all of them. A useful exemple is ``MPI_Alltoall`` where each rank sends and receives
  data from all the other ranks. On the contrary, point to point communication only needs to be called
  by the sending and receiving processes.

  - **Neighborhood collective communications** : this is a subset of collective communications where the actual
    *communication graph* is provided to MPI, in order to reduce the number of communications (especially useful if
    the communication graph is sparse).

  .. warning:: When the exchange helper is used in the point to point exchange mode, *at least* all the
     non-orphan processes in the communication graph must call the exchange function.
     When the exchange helper is used in a collective exchange mode, all the processes must call the exchange function.

- **One sided communications** : unlike the classical two-sided model (MPI_Send / MPI_Recv) where both
  sender and receiver must explicitly participate in the data transfer, one-sided (or RMA) communications
  allow a process to directly read from or write to the memory of another process
  without its immediate involvement. This can happen in specific memory regions exposed as *windows*.

When using the exchange helper, the underlying communication method depends on:

- the exchange function called to switch between blocking, non-blocking or persistent modes;
- the value of the ``PDM_mpi_comm_kind_t`` argument to switch between peer-to-peer, collective or
  one sided modes:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenenum:: PDM_mpi_comm_kind_t

  Note than when ``PDM_MPI_COMM_KIND_COLLECTIVE`` is selected, 
  neighborhood collective are automatically used if the communicator is associated
  to a distributed graph topology (``MPI_Topo_test`` behind the hood).
  Otherwise, standard dense collective communications are used.

Note that some combinations make no sense, see each function for details.

API
"""

.. dropdown:: Initialization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_exchange_helper_create

.. dropdown:: Blocking exchanges

  .. note:: Blocking exchanges are not available in ``PDM_MPI_COMM_KIND_WIN_RMA`` mode

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_exchange_helper_exch

.. dropdown:: Non-blocking exchanges

  Non-blocking exchanges can be either oneshot or persistent (see above).
  In both cases, started exchange must be waited with the ``PDM_exchange_helper_exch_wait`` function.

  .. dropdown:: Oneshot

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_exchange_helper_iexch

  .. dropdown:: Persistent

    Persistent communication allows to reuse the same "communication chanel" (see above).
    A persistent communication is initialized once with a ``init`` method, then uses several
    times with the ``start`` and ``wait`` methods. When the communication chanel
    is no longer needed, it must be finalized with a ``free`` method.

    .. note:: Persistent exchanges are not yet implemented with neighborhood collective mode

    .. tab-set::
      :sync-group: language

      .. tab-item:: C
        :sync: C

        .. doxygenfunction:: PDM_exchange_helper_exch_init
        .. doxygenfunction:: PDM_exchange_helper_exch_start
        .. doxygenfunction:: PDM_exchange_helper_exch_free

  The ``wait`` method is common to oneshot and persistent exchanges:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_exchange_helper_exch_wait



.. dropdown:: "Oneway" exchanges

  Oneway exchanges is an advanced mode of the exchange helper that decorrelates the send and receive part
  of the exchange. It allows developpers to finely manage the flow of the exchanges, for example to
  initiate several sends before starting to receive data.
  This can be useful for specific applications such as code coupling.

  .. note:: The oneway exchanges are non-blocking, and always use a point-to-point mode.

  Oneway exchanges involve the ``PDM_exchange_direction_t`` enum:

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenenum:: PDM_exchange_direction_t

  This enum kind indicates if data is send from provided ``buffer`` or received into provided ``buffer``.

  .. important:: The developer must ensure that each send is matched by a corresponding receive.


  The function exists in oneshot and persistent form. In both cases, the output request can be
  waited with ``PDM_exchange_helper_exch_wait``.
  For the persistent form, functions ``PDM_exchange_helper_exch_start`` and
  ``PDM_exchange_helper_exch_free`` must be used in addition (see above).

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_exchange_helper_iexch_one_way
      .. doxygenfunction:: PDM_exchange_helper_exch_one_way_init


.. dropdown:: Finalization

  .. tab-set::
    :sync-group: language

    .. tab-item:: C
      :sync: C

      .. doxygenfunction:: PDM_exchange_helper_free

