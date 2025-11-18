cdef extern from "pdm_part_comm_graph.h":

  ctypedef struct PDM_part_comm_graph_t:
    pass

  PDM_part_comm_graph_t* PDM_part_comm_graph_create(int               n_part,
                                                    int              *pn_entity_graph,
                                                    int             **pentity_graph,
                                                    PDM_ownership_t   ownership,
                                                    PDM_MPI_Comm      comm);

  PDM_part_comm_graph_t* PDM_part_comm_graph_with_nuplet_create(int               n_part,
                                                                int              *pn_entity_graph,
                                                                int             **pentity_graph,
                                                                PDM_ownership_t   ownership_graph,
                                                                int               nuplet_size,
                                                                int             **pentity_nuplet,
                                                                PDM_ownership_t   ownership_nuplet,
                                                                PDM_bool_t        is_signed,
                                                                PDM_MPI_Comm      comm);

  void PDM_part_comm_graph_exch(PDM_part_comm_graph_t   *pcg,
                                size_t                   s_data,
                                PDM_stride_t             t_stride,
                                int                      cst_stride,
                                int                    **send_entity_stride,
                                void                   **send_entity_data,
                                int                   ***recv_entity_stride,
                                void                  ***recv_entity_data);

  int* PDM_part_comm_graph_owner_get(PDM_part_comm_graph_t *pcg,
                                     int                    i_part);

  int PDM_part_comm_graph_entity_graph_get(PDM_part_comm_graph_t  *pcg,
                                           int                     i_part,
                                           int                   **entity_graph,
                                           PDM_ownership_t         ownership);

  void PDM_part_comm_graph_all_reduce(PDM_part_comm_graph_t   *pcg,
                                      PDM_MPI_Datatype         datatype,
                                      PDM_MPI_Op               op,
                                      unsigned char          **pdata);

  int PDM_part_comm_graph_entity_nuplet_get(PDM_part_comm_graph_t  *pcg,
                                            int                     i_part,
                                            int                   **entity_nuplet,
                                            PDM_ownership_t         ownership);

  PDM_part_comm_graph_t* PDM_part_comm_graph_free(PDM_part_comm_graph_t* pcg);


# ========================================================================
# ------------------------------------------------------------------------
cdef class PartCommGraphCapsule:
  """

  """
  cdef PDM_part_comm_graph_t *pcg

  def __cinit__(self, object caps):
    self.pcg = <PDM_part_comm_graph_t *> PyCapsule_GetPointer(caps, NULL);

  def get_entity_graph(self, int i_part):
    return entity_graph_get(self, i_part)

# ========================================================================
# ------------------------------------------------------------------------
cdef class PartCommGraph:

  cdef PDM_part_comm_graph_t *pcg
  cdef MPI.Comm               py_comm

  cdef int                   *_pn_entity_graph
  cdef int                  **_pentity_graph
  cdef int                  **_pentity_nuplet

  def __init__(self,
               MPI.Comm    comm,
               int         n_part,
               list        pn_entity_graph,
               list        pentity_graph):
    """
    __init__(comm, n_part)

    Create a new :py:class`PartCommGraph` instance

    Parameters:
      comm (MPI.Comm) : MPI communicator
      n_part (int) : Number of partitions
      pn_entity_grpah (int*) : Number of bound (size = \p n_part)
      pentity_graph   (int**) : Graph comm identifier (size = 4 * \p pn_entity_graph[i_part]):
        For each entity :
          - entity local number (1-based)
          - Connected process   (0-based)
          - Connected partition on the connected process (1-based)
          - Connected entity local number in the connected partition (1-based)
    """

    self.py_comm  = comm

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

    self._pn_entity_graph = list_to_int_pointer(pn_entity_graph)
    self._pentity_graph   = np_list_to_int_pointers(pentity_graph)

    self._n_part = n_part

    self.pcg =  PDM_part_comm_graph_create(self._n_part,
                                  <int *>  self._pn_entity_graph,
                                  <int **> self._pentity_graph,
                                  PDM_OWNERSHIP_USER,
                                  PDMC)

  def with_nuplet_create(self,
                         MPI.Comm comm,
                         int      n_part,
                         list     pn_entity_graph,
                         list     pentity_graph,
                         int      nuplet_size,
                         list     pentity_nuplet,
                         bint     is_signed):
    """
    Add doc
    """
    self.py_comm  = comm

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI ... FIXME: what is happening here ?
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

    self._pn_entity_graph = list_to_int_pointer(pn_entity_graph)
    self._pentity_graph   = np_list_to_int_pointers(pentity_graph)

    self._pentity_nuplet  = np_list_to_int_pointers(pentity_nuplet)
    self._nuplet_size     = nuplet_size

    self._n_part = n_part

    self.pcg = PDM_part_comm_graph_with_nuplet_create(self._n_part,
                                                      self._pn_entity_graph,
                                                      self._pentity_graph,
                                                      PDM_OWNERSHIP_USER,
                                                      self._nuplet_size,
                                                      self._pentity_nuplet,
                                                      PDM_OWNERSHIP_USER,
                                         <PDM_bool_t> is_signed,
                                         PDMC)


  def exch(self,
           list send_entity_data,
           int  send_entity_stride=1,
           bint interlaced_str=True):
    """
    Add doc
    """
    return exch(self, send_entity_data, send_entity_stride, interlaced_str)

  def owner_get(self, int i_part):
    """
    Add doc
    """
    cdef int* owner = PDM_part_comm_graph_owner_get(self.pcg,i_part) #returns an int*

    return create_numpy_i(owner, self._pn_entity_graph[i_part])

  def entity_graph_get(self, int i_part):
    """
    Add doc
    """
    return entity_graph_get(self, i_part)

  def all_reduce(self,
      MPI.Datatype    datatype,
      MPI.Op          op,
      list pdata
      ):
    """
    Add doc
    """
    all_reduce(self, datatype, op, pdata)

  def entity_nuplet_get(self):
    """
    Add doc
    """

  def _dealloc__(self):
    """
    Add doc
    """
    PDM_part_comm_graph_free(self.pcg)

# ------------------------------------------------------------------------
ctypedef fused PyPartCommGraph:
  PartCommGraph
  PartCommGraphCapsule

# ------------------------------------------------------------------------
def exch(PyPartCommGraph pypcg,
         list            send_entity_data,
         int             send_entity_stride=1,
         bint            interlaced_str=True):
  """
  Add doc
  """

  cdef int request_exch
  cdef PDM_stride_t _stride_t

  cdef int _stride_cst = 0
  cdef int** _send_entity_stride = NULL
  if isinstance(send_entity_stride, int):
    _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
    _stride_cst = send_entity_stride
  elif isinstance(send_entity_stride, list):
    _stride_t = PDM_STRIDE_VAR_INTERLACED

    #FIXME: recuperer la size de la data à envoyer ?
    # assert len(send_entity_stride) == pypcg.

    _send_entity_stride = np_list_to_int_pointers(send_entity_stride)
  else:
    raise ValueError("Invalid stride in pcg ech")

  cdef void** _send_entity_data = np_list_to_void_pointers(send_entity_data)

  ref_dtype = recover_dtype(send_entity_data, pypcg.py_comm)
  cdef size_t s_data   = ref_dtype.itemsize
  cdef size_t npy_type = ref_dtype.num

  cdef int**   _recv_entity_stride = NULL;
  cdef void** _recv_entity_data   = NULL;

  PDM_part_comm_graph_exch(pypcg.pcg,
                           s_data,
                           _stride_t,
                           _stride_cst,
               <int  ** >  _send_entity_stride,
               <void ** >  _send_entity_data,
               <int  ***> &_recv_entity_stride,
               <void ***> &_recv_entity_data)



# ------------------------------------------------------------------------
def entity_graph_get(PyPartCommGraph pypcg, int i_part):
  """
  entity_graph_get(i_part)
    Parameters:
      i_part (int) : partition identifier

    Returns:
      Graph comm identifier (`int *`)

  """
  cdef int *entity_graph = NULL

  cdef n_entity = PDM_part_comm_graph_entity_graph_get(
     pypcg.pcg,
     i_part,
    &entity_graph,
     PDM_OWNERSHIP_USER)

  return create_numpy_i(entity_graph, n_entity)

# ------------------------------------------------------------------------
def all_reduce(PyPartCommGraph pypcg,
               MPI.Datatype    datatype,
               MPI.Op          op,
               list            pdata):
  """
  Add doc
  """
  # FIXME: modification in-place <= il faut que je recupère ce tableau....
  cdef void **_pdata = np_list_to_void_pointers(pdata)
  cdef PDM_MPI_Datatype c_datatype = <MPI_Datatype> datatype.ob_mpi
  cdef PDM_MPI_Op       c_op       = <MPI_Op      > op.ob_mpi

  PDM_part_comm_graph_all_reduce(pypcg.pcg,
                                 c_datatype,
                                 c_op,
              <unsigned char **> &_pdata)

# ------------------------------------------------------------------------
def entity_nuplet_get(PyPartCommGraph pypcg, int i_part):
  """
  Add doc
  """
  cdef int *entity_nuplet = NULL
  cdef int  n_entity
  cdef int  size_of_nuplet

  n_entity = PDM_part_comm_graph_entity_nuplet_get(pypcg.pcg,
                                                   i_part,
                                                   &entity_nuplet,
                                                   PDM_OWNERSHIP_USER)

  return create_numpy_i(entity_nuplet, pypcg._pn_entity_graph[i_part] * pypcg._nuplet_size)
