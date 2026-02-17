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

  int PDM_part_comm_graph_n_part_get(PDM_part_comm_graph_t *pcg);
  int* PDM_part_comm_graph_owner_get(PDM_part_comm_graph_t *pcg,
                                     int                    i_part);

  int PDM_part_comm_graph_entity_graph_get(PDM_part_comm_graph_t  *pcg,
                                           int                     i_part,
                                           int                   **entity_graph,
                                           PDM_ownership_t         ownership);

  void PDM_part_comm_graph_all_reduce(PDM_part_comm_graph_t   *pcg,
                                      PDM_MPI_Datatype         datatype,
                                      int                      stride,
                                      PDM_MPI_Op               op,
                                      unsigned char          **pdata);

  int PDM_part_comm_graph_entity_nuplet_get(PDM_part_comm_graph_t  *pcg,
                                            int                     i_part,
                                            int                   **entity_nuplet,
                                            PDM_ownership_t         ownership);

  PDM_part_comm_graph_t* PDM_part_comm_graph_free(PDM_part_comm_graph_t* pcg);

cdef extern from "pdm_part_comm_graph_algorithm.h":

  void PDM_part_comm_graph_entity1_to_part_comm_graph_entity2(PDM_part_comm_graph_t   *ptpgc_entity1,
                                                              int                     *pn_entity1,
                                                              int                     *pn_entity2,
                                                              int                    **entity2_entity1_idx,
                                                              int                    **entity2_entity1,
                                                              PDM_part_comm_graph_t  **out_ptpgc_entity2)

# ========================================================================
# ------------------------------------------------------------------------
cdef class PartCommGraph:

  cdef PDM_part_comm_graph_t *pcg
  cdef MPI.Comm               py_comm

  def __init__(self,
               MPI.Comm    comm,
               list        pentity_graph,
               list        pentity_nuplet=None,
               bint        is_signed=True):
    """
    __init__(comm, pentity_graph, pentity_nuplet=None, is_signed=True)

    Create a new :py:class`PartCommGraph` instance

    Parameters:
      comm            (MPI.Comm) : MPI communicator
      pentity_graph   (int**)    : Graph comm identifier (size = 4 * \p pn_entity_graph[i_part]):
        For each entity :
          - entity local number (1-based)
          - Connected process   (0-based)
          - Connected partition on the connected process (1-based)
          - Connected entity local number in the connected partition (1-based)
      pentity_nuplet (list of np.ndarray[in]) : Additional nuplets (size = \p nuplet_size * \p pn_entity_graph[i_part])
      is_signed      (int**)                  : Use signed nuplets
    """

    self.py_comm  = comm

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

    cdef int **_pentity_graph  = NULL
    cdef int **_pentity_nuplet = NULL

    _n_part          = len(pentity_graph)
    cdef int* _pn_entity_graph = list_to_int_pointer([g.size // 4 for g in pentity_graph])
    _pentity_graph        = np_list_to_int_pointers(pentity_graph)

    _nuplet_size = 0
    if pentity_nuplet is None:
      self.pcg =  PDM_part_comm_graph_create(_n_part,
                                             _pn_entity_graph,
                                    <int **> _pentity_graph,
                                             PDM_OWNERSHIP_USER,
                                             PDMC)
    else:
      _pentity_nuplet = np_list_to_int_pointers(pentity_nuplet)
      for i_part in range(_n_part):
        _nuplet_size = pentity_nuplet[i_part].size // _pn_entity_graph[i_part]

      self.pcg = PDM_part_comm_graph_with_nuplet_create(_n_part,
                                                        _pn_entity_graph,
                                                        _pentity_graph,
                                                        PDM_OWNERSHIP_USER,
                                                        _nuplet_size,
                                                        _pentity_nuplet,
                                                        PDM_OWNERSHIP_USER,
                                           <PDM_bool_t> is_signed,
                                                        PDMC)
      free(_pentity_nuplet)
    free(_pentity_graph)

  @staticmethod
  cdef from_ptr(PDM_part_comm_graph_t* ptr, py_comm):
    cdef PartCommGraph obj = PartCommGraph.__new__(PartCommGraph)
    obj.py_comm  = py_comm
    obj.pcg      = ptr

    return obj

  def exch(self,
           list send_entity_data,
           send_entity_stride=1,
           bint interlaced_str=True):
    """
      exch(send_entity_data, send_entity_stride=1, interlaced_str=True)

      Exchange data between graph comm with synchronous blocking exchange

      Parameters:
        send_entity_data   (list)                      : Graph data for each part
        send_entity_stride (`int` or `list, optional`) : Stride of Part1 data
        interlaced_str    (bool, optional)             : Is the data interlaced? (default = **True**)

      Returns :
        - Recv stride (same dtype as ``send_entity_stride`` )
        - Recv data   (`list` of same dtype as ``send_entity_data``)
    """

    cdef int request_exch
    cdef PDM_stride_t _stride_t

    cdef int   _stride_cst = 0
    cdef int** _send_entity_stride = NULL
    if isinstance(send_entity_stride, int):
      _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
      _stride_cst = send_entity_stride
    elif isinstance(send_entity_stride, list):
      _stride_t = PDM_STRIDE_VAR_INTERLACED
      _send_entity_stride = np_list_to_int_pointers(send_entity_stride)
    else:
      raise ValueError("Invalid stride in pcg ech")

    cdef void** _send_entity_data = np_list_to_void_pointers(send_entity_data)

    ref_dtype = recover_dtype(send_entity_data, self.py_comm)
    cdef size_t s_data   = ref_dtype.itemsize
    cdef size_t npy_type = ref_dtype.num

    cdef int  **_recv_entity_stride = NULL
    cdef void **_recv_entity_data   = NULL

    PDM_part_comm_graph_exch(self.pcg,
                            s_data,
                            _stride_t,
                            _stride_cst,
                <int  ** >  _send_entity_stride,
                <void ** >  _send_entity_data,
                <int  ***> &_recv_entity_stride,
                <void ***> &_recv_entity_data)

    lnp_part_strid = list()
    lnp_part_data  = list()
    cdef int pn_entity
    cdef int *dummy = NULL
    for i_part in range(PDM_part_comm_graph_n_part_get(self.pcg)):
      pn_entity = PDM_part_comm_graph_entity_graph_get(self.pcg,
                                                        i_part,
                                                        &dummy,
                                                        PDM_OWNERSHIP_BAD_VALUE)
      if _stride_t == PDM_STRIDE_VAR_INTERLACED:

        strid_size = pn_entity

        np_part2_stride = create_numpy_i(_recv_entity_stride[i_part], strid_size)
        dim_np = np_part2_stride.sum()

        np_part2_data = create_numpy(_recv_entity_data[i_part], npy_type, dim_np)

        lnp_part_strid.append(np_part2_stride)
        lnp_part_data .append(np_part2_data)

      elif _stride_t == PDM_STRIDE_CST_INTERLACED:
        dim_np  = pn_entity * _stride_cst
        np_part2_data = create_numpy(_recv_entity_data[i_part], npy_type, dim_np)

        lnp_part_data .append(np_part2_data)

    free(_send_entity_data)
    if _stride_t == PDM_STRIDE_VAR_INTERLACED:
      free(_send_entity_stride)
      free(_recv_entity_stride)
    free(_recv_entity_data)

    return lnp_part_strid, lnp_part_data


  def owner_get(self, int i_part):
    """
    owner_get(i_part)

    Get the owner array computed inside the structure, useful to manage reduction of array for example

    Parameters:
      i_part (int) : Partition identifier

    Returns:
      Owner array, 0 is not owner, 1 is owner  (`np.array[np.int]`)
    """
    cdef int *dummy = NULL
    cdef pn_entity = PDM_part_comm_graph_entity_graph_get(self.pcg,
                                                          i_part,
                                                          &dummy,
                                                          PDM_OWNERSHIP_BAD_VALUE)
    
    cdef int* owner = PDM_part_comm_graph_owner_get(self.pcg,i_part) #returns an int*

    np_owner = create_numpy_i(owner, pn_entity, flag_owndata=False)
    return NPY.copy(np_owner)

  def entity_graph_get(self, int i_part):
    """
    entity_graph_get(i_part)

    Get entity graph

    Parameters:
      i_part (int) : Partition identifier

    Returns:
      Graph comm identifier (`np.array[np.int]`):
        For each entity :
          - entity local number (1-based)
          - Connected process   (0-based)
          - Connected partition on the connected process (1-based)
          - Connected entity local number in the connected partition (1-based)

    """
    cdef int *entity_graph = NULL

    cdef n_entity = PDM_part_comm_graph_entity_graph_get(self.pcg,
                                                        i_part,
                                                        &entity_graph,
                                                        PDM_OWNERSHIP_USER)

    return create_numpy_i(entity_graph, 4 * n_entity)

  def all_reduce(self,
                 int             stride,
                 MPI.Op          op,
                 list            pdata):
    """
    all_reduce(stride, op, pdata)

      Parameters:
        stride   (int)                              : Constant data stride
        op       (MPI.Op)                           : Reduction operation kind (SUM/MIN/MAX)
        pdata    (`list` of `np.ndarray[datatype]`) : Data buffer, value is modified inplace

    """
    cdef void **_pdata = np_list_to_void_pointers(pdata)
    cdef PDM_MPI_Op       c_op       = <MPI_Op      > op.ob_mpi
    cdef MPI.Datatype mpi_dtype
    ref_dtype = recover_dtype(pdata, self.py_comm)
    mpi_dtype = MPI._typedict.get(ref_dtype.char)
    # > This one not working, it seems it give a special type
    # from mpi4py.util.dtlib import from_numpy_dtype
    # mpi_dtype = from_numpy_dtype(ref_dtype)
    cdef PDM_MPI_Datatype c_datatype = <MPI_Datatype> mpi_dtype.ob_mpi

    PDM_part_comm_graph_all_reduce(self.pcg,
                                  c_datatype,
                                  stride,
                                  c_op,
                <unsigned char **> _pdata)
    free(_pdata)

  def entity_nuplet_get(self, i_part):
    """
    entity_nuplet_get(i_part)

    Get entity nuplet

    Parameters:
      i_part (int) : Partition identifier
    """
    cdef int *entity_graph  = NULL
    cdef int *entity_nuplet = NULL
    cdef int  n_entity
    cdef int  size_of_nuplet

    n_entity = PDM_part_comm_graph_entity_graph_get(self.pcg,
                                                    i_part,
                                                    &entity_graph,
                                                    PDM_OWNERSHIP_USER)
    nuplet_size = PDM_part_comm_graph_entity_nuplet_get(self.pcg,
                                                        i_part,
                                                        &entity_nuplet,
                                                        PDM_OWNERSHIP_USER)

    np_entity_nuplet = None
    if (entity_nuplet != NULL):
      np_entity_nuplet = create_numpy_i(entity_nuplet, n_entity * nuplet_size)

    return np_entity_nuplet

  def __dealloc__(self):
    """
    """
    PDM_part_comm_graph_free(self.pcg)

# ------------------------------------------------------------------------



# ------------------------------------------------------------------------
def pcg_entity1_to_entity2(PartCommGraph pypcg_entity1,
                           list          pn_entity1,
                           list          entity2_entity1_idx,
                           list          entity2_entity1):
  """
  Returns a \ref PDM_part_comm_graph_t python object

  Parameters:
    pypcg_entity1       (PartCommGraph) : \ref PDM_part_comm_graph_t structure for entity1
    pn_entity1          (list         ) : Number of entity1 (size = n_part)
    entity2_entity1_idx (list         ) : Connectivity index (size = \p pn_entity2 + 1)
    entity2_entity1     (list         ) : Connectivity array (size = \p entity2_entity1_idx[\p pn_entity2] )
  """

  cdef int  *_pn_entity1 = NULL
  cdef int  *_pn_entity2 = NULL
  cdef int **_entity2_entity1_idx = NULL
  cdef int **_entity2_entity1     = NULL

  pn_entity2 = [part_entity2_entity1_idx.size -1 for part_entity2_entity1_idx in entity2_entity1_idx]

  _pn_entity1 = list_to_int_pointer(pn_entity1)
  _pn_entity2 = list_to_int_pointer(pn_entity2)
  _entity2_entity1_idx = np_list_to_int_pointers(entity2_entity1_idx)
  _entity2_entity1     = np_list_to_int_pointers(entity2_entity1)

  cdef PDM_part_comm_graph_t *_out_ptpgc_entity2 = NULL;
  PDM_part_comm_graph_entity1_to_part_comm_graph_entity2(pypcg_entity1.pcg,
                                                         _pn_entity1,
                                                         _pn_entity2,
                                                         _entity2_entity1_idx,
                                                         _entity2_entity1,
                                                         &_out_ptpgc_entity2)

  free(_pn_entity1)
  free(_pn_entity2)
  free(_entity2_entity1_idx)
  free(_entity2_entity1)

  return PartCommGraph.from_ptr(_out_ptpgc_entity2, None) #TODO retrive comm
