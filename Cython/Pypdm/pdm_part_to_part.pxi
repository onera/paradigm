import warnings

cdef extern from "pdm_part_to_part.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_to_part_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_part_to_part_t *PDM_part_to_part_create(PDM_g_num_t   **gnum_elt1,
                                                int            *n_elt1,
                                                int             n_part1,
                                                PDM_g_num_t   **gnum_elt2,
                                                int            *n_elt2,
                                                int             n_part2,
                                                int           **part1_to_part2_idx,
                                                PDM_g_num_t   **part1_to_part2,
                                                PDM_MPI_Comm    comm);
    PDM_part_to_part_t *PDM_part_to_part_create_from_num2_triplet(const PDM_g_num_t   **gnum_elt1,
                                                                  const int            *n_elt1,
                                                                  const int             n_part1,
                                                                  const int            *n_elt2,
                                                                  const int             n_part2,
                                                                  const int           **part1_to_part2_idx,
                                                                  const int           **part1_to_part2_triplet_idx,
                                                                  const int           **part1_to_part2_triplet,
                                                                  const PDM_MPI_Comm    comm);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_iexch(PDM_part_to_part_t                *ptp,
                                const PDM_mpi_comm_kind_t          k_comm,
                                const PDM_stride_t                 t_stride,
                                const PDM_part_to_part_data_def_t  t_part1_data_def,
                                const int                          cst_stride,
                                const size_t                       s_data,
                                const int                        **part1_stride,
                                const void                       **part1_data,
                                int                             ***part2_stride,
                                void                            ***part2_data,
                                int                               *request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_iexch_wait(PDM_part_to_part_t                *ptp,
                                     const int                          request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_gnum1_come_from_get(PDM_part_to_part_t *ptp,
                                              int              ***gnum1_come_from_idx,
                                              PDM_g_num_t      ***gnum1_come_from);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_ref_lnum2_get(PDM_part_to_part_t   *ptp,
                                        int                 **n_ref_lnum2,
                                        int                ***ref_lnum2);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_unref_lnum2_get(PDM_part_to_part_t   *ptp,
                                          int                 **n_unref_lnum2,
                                          int                ***unref_lnum2);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_reverse_iexch(PDM_part_to_part_t             *ptp,
                                        PDM_mpi_comm_kind_t             k_comm,
                                        PDM_stride_t                    t_stride,
                                        PDM_part_to_part_data_def_t     t_part2_data_def,
                                        int                             cst_stride,
                                        size_t                          s_data,
                                        int                           **part2_stride,
                                        void                          **part2_data,
                                        int                          ***part1_stride,
                                        void                         ***part1_data,
                                        int                            *request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_reverse_iexch_wait(PDM_part_to_part_t                *ptp,
                                             const int                          request);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_to_part_n_part_get( PDM_part_to_part_t *ptp    ,
                                      int                *n_part1,
                                      int                *n_part2);
    void PDM_part_to_part_n_part_and_n_elt_get( PDM_part_to_part_t *ptp    ,
                                                int                *n_part1,
                                                int                *n_part2,
                                                int               **n_elt1 ,
                                                int               **n_elt2 );

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # WARNING : do not wrap PDM_part_to_part_part1_to_part2_(single_part_)get or
    # because part1_to_part2 can be NULL (if PtP is created from triplets)
    # or can reference invalid memory (if pointer is free after creation)
    # PDM_part_to_part_part1_to_part2_idx_get is safe because PtP structure owns this memory
    void PDM_part_to_part_part1_to_part2_idx_get(PDM_part_to_part_t *ptp               ,
                                                 int               **n_elt1            ,
                                                 int              ***part1_to_part2_idx);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    PDM_MPI_Comm PDM_part_to_part_comm_get(PDM_part_to_part_t *ptp);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    PDM_part_to_part_t *PDM_part_to_part_free(PDM_part_to_part_t *ptp);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



# ========================================================================
# ------------------------------------------------------------------------
cdef class PartToPart:
  # ************************************************************************
  # > Class attributes
  cdef PDM_part_to_part_t         *ptp
  
  cdef dict                        request_data

  DATA_DEF_ORDER_PART1           = _PDM_PART_TO_PART_DATA_DEF_ORDER_PART1
  DATA_DEF_ORDER_PART1_TO_PART2  = _PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2
  DATA_DEF_ORDER_PART2           = _PDM_PART_TO_PART_DATA_DEF_ORDER_PART2
  DATA_DEF_ORDER_GNUM1_COME_FROM = _PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM

  # ************************************************************************
  def __init__(self, MPI.Comm comm,
                      list part1_ln_to_gn,
                      list part2_ln_to_gn,
                      list part1_to_part2_idx,
                      list part1_to_part2):
    """
    __init__(comm, part1_ln_to_gn, part2_ln_to_gn, part1_to_part2_idx, part1_to_part2)

    Create a Part-to-Part redistribution from global ids

    Parameters:
      comm               (MPI.Comm)                               : MPI communicator
      part1_ln_to_gn     (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Element global ids in Part1
      part2_ln_to_gn     (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Element global ids in Part2
      part1_to_part2_idx (`list` of `np.ndarray[np.int32_t]`)     : Index for Part1→Part2 mapping
      part1_to_part2     (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Part1→Part2 mapping (global ids)
    """
    n_part1 = len(part1_ln_to_gn)
    n_part2 = len(part2_ln_to_gn)
    self.request_data = dict()

    assert(len(part1_to_part2_idx) == n_part1)

    # > Convert input data
    cdef PDM_MPI_Comm pdm_comm = py_comm_to_pdm_comm(comm)

    cdef int* n_elt1 = list_to_int_pointer([array.size for array in part1_ln_to_gn])
    cdef int* n_elt2 = list_to_int_pointer([array.size for array in part2_ln_to_gn])

    cdef PDM_g_num_t** _part1_ln_to_gn     = np_list_to_gnum_pointers(part1_ln_to_gn)
    cdef PDM_g_num_t** _part2_ln_to_gn     = np_list_to_gnum_pointers(part2_ln_to_gn)
    cdef int**         _part1_to_part2_idx = np_list_to_int_pointers (part1_to_part2_idx)
    cdef PDM_g_num_t** _part1_to_part2     = np_list_to_gnum_pointers(part1_to_part2)

    self.ptp = PDM_part_to_part_create(<const PDM_g_num_t **> _part1_ln_to_gn,
                                       <const int *>          n_elt1,
                                                              n_part1,
                                       <const PDM_g_num_t **> _part2_ln_to_gn,
                                       <const int *>          n_elt2,
                                                              n_part2,
                                       <const int **>         _part1_to_part2_idx,
                                       <const PDM_g_num_t **> _part1_to_part2,
                                                              pdm_comm)
    free(n_elt1) # Copied in part to part structure
    free(n_elt2)
    free(_part1_ln_to_gn)
    free(_part2_ln_to_gn)
    free(_part1_to_part2_idx)
    free(_part1_to_part2)

  # ------------------------------------------------------------------------
  @staticmethod
  cdef from_ptr(PDM_part_to_part_t* ptr):

    cdef PartToPart obj = PartToPart.__new__(PartToPart)

    obj.ptp      = ptr
    obj.request_data = dict()

    return obj


  @classmethod
  def from_triplet(cls, MPI.Comm comm,
                        list part1_ln_to_gn,
                        list n_elt2,
                        list part1_to_part2_idx,
                        list part1_to_part2_triplet_idx,
                        list part1_to_part2_triplet):
    """
    from_triplet(comm, part1_ln_to_gn, n_elt2, part1_to_part2_idx, part1_to_part2_triplet_idx, part1_to_part2_triplet)

    An alternative constructor (classmethod) that create a PDM_part_to_part object from part2 triplet.

    Parameters:
      comm                       (MPI.Comm)                               : MPI communicator
      part1_ln_to_gn             (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Element global ids in Part1
      n_elt2                     (`list` of `int`)                        : Number of elements in Part2
      part1_to_part2_idx         (`list` of `np.ndarray[np.int32_t]`)     : Index for Part1→Part2 mapping
      part1_to_part2_triplet_idx (`list` of `np.ndarray[np.int32_t]`)     : Index for Part2 triplets
      part1_to_part2             (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Part1→Part2 mapping (i_rank, i_part, lnum) triplet
    """

    cdef PartToPart obj = PartToPart.__new__(PartToPart)
    cdef PDM_MPI_Comm pdm_comm = py_comm_to_pdm_comm(comm)

    n_part1 = len(part1_ln_to_gn)
    n_part2 = len(n_elt2)

    assert(len(part1_to_part2_idx) == n_part1)

    cdef int* _n_elt1 = list_to_int_pointer([array.size for array in part1_ln_to_gn])
    cdef int* _n_elt2 = list_to_int_pointer(n_elt2)

    cdef PDM_g_num_t** _part1_ln_to_gn             = np_list_to_gnum_pointers(part1_ln_to_gn)
    cdef int**         _part1_to_part2_idx         = np_list_to_int_pointers (part1_to_part2_idx)
    cdef int**         _part1_to_part2_triplet_idx = NULL
    cdef int** _part1_to_part2_triplet             = np_list_to_int_pointers (part1_to_part2_triplet)

    if part1_to_part2_triplet_idx is not None:
      _part1_to_part2_triplet_idx = np_list_to_int_pointers (part1_to_part2_triplet_idx)

    cdef PDM_part_to_part_t *ptp = NULL
    ptp = PDM_part_to_part_create_from_num2_triplet(<const PDM_g_num_t **> _part1_ln_to_gn,
                                                    <const int          *> _n_elt1,
                                                                           n_part1,
                                                    <const int          *> _n_elt2,
                                                                           n_part2,
                                                    <const int         **> _part1_to_part2_idx,
                                                    <const int         **> _part1_to_part2_triplet_idx,
                                                    <const int         **> _part1_to_part2_triplet,
                                                                            pdm_comm)

    free(_n_elt1)
    free(_n_elt2)
    free(_part1_ln_to_gn)
    free(_part1_to_part2_idx)
    free(_part1_to_part2_triplet_idx)
    free(_part1_to_part2_triplet)

    obj.ptp = ptp
    obj.request_data = dict()

    return obj


  # --------------------------------------------------------------------
  def get_referenced_lnum2(self):
    """
    Get referenced Part2 elements

    Returns:
      Referenced Part2 elements (1-based local ids) (`list` of `np.ndarray[np_int32_t]`)
    """
    cdef int  *n_ref_lnum2 = NULL;
    cdef int **ref_lnum2   = NULL;
    cdef int n_part1, n_part2
    PDM_part_to_part_n_part_get(self.ptp, &n_part1, &n_part2)
    PDM_part_to_part_ref_lnum2_get(self.ptp,
                                  &n_ref_lnum2,
                                  &ref_lnum2);
    lnp_ref_lnum2 = list()
    for i_part in range(n_part2):
      np_ref_lnum2_id = create_numpy_i(ref_lnum2[i_part], n_ref_lnum2[i_part], flag_owndata=False)
      lnp_ref_lnum2.append(NPY.copy(np_ref_lnum2_id))
    return lnp_ref_lnum2

  # --------------------------------------------------------------------
  def get_unreferenced_lnum2(self):
    """
    Get unreferenced Part2 elements

    Returns:
      Unreferenced Part2 elements (1-based local ids) (`list` of `np.ndarray[np_int32_t]`)
    """
    cdef int  *n_unref_lnum2 = NULL;
    cdef int **unref_lnum2   = NULL;
    cdef int n_part1, n_part2
    PDM_part_to_part_n_part_get(self.ptp, &n_part1, &n_part2)
    PDM_part_to_part_unref_lnum2_get(self.ptp,
                                  &n_unref_lnum2,
                                  &unref_lnum2);
    lnp_unref_lnum2 = list()
    for i_part in range(n_part2):
      np_unref_lnum2_id = create_numpy_i(unref_lnum2[i_part], n_unref_lnum2[i_part], flag_owndata=False)
      lnp_unref_lnum2.append(NPY.copy(np_unref_lnum2_id))
    return lnp_unref_lnum2

  # --------------------------------------------------------------------
  def get_gnum1_come_from(self):
    """
    Get Part2→Part1 mapping for referenced Part2 elements

    Returns:
      Dictionary
        - ``"come_from_idx"`` (`list` of `np.ndarray[np.int32_t]`)     : Index for Part2→Part1 mapping
        - ``"come_from"``     (`list` of `np.ndarray[npy_pdm_gnum_t]`) : Part2→Part1 mapping (global ids)
    """
    # Needed to get the sizes n_ref_lnum2
    cdef int  *n_ref_lnum2 = NULL;
    cdef int **ref_lnum2   = NULL;
    cdef int n_part1, n_part2
    PDM_part_to_part_n_part_get(self.ptp, &n_part1, &n_part2)
    PDM_part_to_part_ref_lnum2_get(self.ptp,
                                  &n_ref_lnum2,
                                  &ref_lnum2);

    cdef int         **gnum1_come_from_idx = NULL
    cdef PDM_g_num_t **gnum1_come_from     = NULL
    PDM_part_to_part_gnum1_come_from_get(self.ptp,
                                        &gnum1_come_from_idx,
                                        &gnum1_come_from)
    come_from_l = list()
    for i_part in range(n_part2):
      n_ref = n_ref_lnum2[i_part]
      # Idx array is shaped as n_ref_lnum2
      sending_idx   = create_numpy_i(gnum1_come_from_idx[i_part], n_ref+1, flag_owndata=False)
      sending_array = create_numpy_g(gnum1_come_from[i_part], gnum1_come_from_idx[i_part][n_ref], flag_owndata=False)
      come_from_l.append({'come_from_idx' : NPY.copy(sending_idx), 'come_from' : NPY.copy(sending_array)})
    return come_from_l

  # --------------------------------------------------------------------
  def iexch(                            self,
            PDM_mpi_comm_kind_t         k_comm,
            PDM_part_to_part_data_def_t t_part1_data_def,
            list                        part1_data,
            part1_stride=1,
            bint interlaced_str=True):
    """
    iexch(k_comm, t_part1_data_def, part1_data, part1_stride=1, interlaced_str=True)

    Initiate a non-blocking exchange (Part1→Part2)

    Parameters:
      k_comm           (int)                       : Kind of MPI communication
      t_part1_data_def (int)                       : Kind of Part1 data definition
      part1_data       (list)                      : Part1 data
      part1_stride     (`int` or `list, optional`) : Stride of Part1 data
      interlaced_str   (bool, optional)            : Is the data interlaced? (default = **True**)

    Returns:
      Request ID (`int`)

    Admissible values for ``k_comm``:
      - 0 : Peer-to-peer (``MPI_issend``/``MPI_irecv``)
      - 1 : Collective (``MPI_Ialltoall`` and alike)

    .. note::
      Additional communication kinds will be available in the future

    Admissible values for ``t_part1_data_def``:
      - :py:attr:`PartToPart.DATA_DEF_ORDER_PART1`          : Data defined according to the `part1` arrays order
      - :py:attr:`PartToPart.DATA_DEF_ORDER_PART1_TO_PART2` : Data defined according to the `part1_to_part2` arrays order
    """
    cdef int request_exch
    cdef PDM_stride_t _stride_t

    cdef int n_part1
    cdef int* n_elt1
    PDM_part_to_part_n_part_and_n_elt_get(self.ptp, &n_part1, NULL, &n_elt1, NULL)

    cdef int   _part1_stride_cst = 0
    cdef int** _part1_stride     = NULL
    if isinstance(part1_stride, int):
      _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
      _part1_stride_cst = part1_stride
    elif isinstance(part1_stride, list):
      _stride_t = PDM_STRIDE_VAR_INTERLACED
      assert len(part1_stride) == n_part1
      for i_part in range(n_part1):
        assert_single_dim_np(part1_stride[i_part], NPY.int32, n_elt1[i_part])
      _part1_stride = np_list_to_int_pointers(part1_stride)
    else:
      raise ValueError("Invalid stride in PtB exchange")

    cdef void** _part1_data   = np_list_to_void_pointers(part1_data)

    cdef PDM_MPI_Comm pdm_comm = PDM_part_to_part_comm_get(self.ptp)
    py_comm = pdm_comm_to_py_comm(pdm_comm)
    
    ref_dtype = recover_dtype(part1_data, py_comm)
    cdef size_t s_data   = ref_dtype.itemsize
    cdef size_t npy_type = ref_dtype.num

    cdef int**  _part2_stride = NULL;
    cdef void** _part2_data   = NULL;
    PDM_part_to_part_iexch(self.ptp,
                          k_comm,
                          _stride_t,
                          t_part1_data_def,
                          _part1_stride_cst,
                          s_data,
          <const int ** >  _part1_stride,
          <const void **>  _part1_data,
                          &_part2_stride,
                <void ***> &_part2_data,
                          &request_exch)
    self.request_data[request_exch] = [<uintptr_t> _part2_stride, <uintptr_t> _part2_data, _part1_stride_cst, npy_type]



    if _stride_t == PDM_STRIDE_VAR_INTERLACED:
      free(_part1_stride)
    free(_part1_data)

    return request_exch

  # --------------------------------------------------------------------
  def wait(self, int request_id):
    """
    wait(request_id)

    Finalize a non-blocking exchange (Part1→Part2)

    Parameters:
      request_id (int) : Request ID

    Returns:
      - Part2 stride (same dtype as ``part1_stride`` in :py:func:`iexch`)
      - Part2 data   (`list` of same dtype as ``part1_data`` in :py:func:`iexch`)
    """
    PDM_part_to_part_iexch_wait(self.ptp, request_id)

    cdef int         **gnum1_come_from_idx = NULL
    cdef PDM_g_num_t **gnum1_come_from     = NULL
    PDM_part_to_part_gnum1_come_from_get(self.ptp,
                                        &gnum1_come_from_idx,
                                        &gnum1_come_from)

    cdef int  *n_ref_lnum2 = NULL;
    cdef int **ref_lnum2   = NULL;
    PDM_part_to_part_ref_lnum2_get(self.ptp,
                                  &n_ref_lnum2,
                                  &ref_lnum2);

    requested = self.request_data.pop(request_id)
    cdef uintptr_t _part2_stride_id = requested[0]
    cdef uintptr_t _part2_data_id   = requested[1]
    cdef int**  _part2_stride = <int**> _part2_stride_id
    cdef void** _part2_data = <void **> _part2_data_id
    cdef NPY.npy_intp dim_np

    # cdef int *__part2_stride = NULL
    cdef int n_part1, n_part2
    PDM_part_to_part_n_part_get(self.ptp, &n_part1, &n_part2)

    lnp_part_strid = list()
    lnp_part_data  = list()
    cdef size_t npy_type = requested[3]
    for i_part in range(n_part2):

      if(_part2_stride != NULL):

        strid_size = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]]

        np_part2_stride = create_numpy_i(_part2_stride[i_part], strid_size)
        dim_np = np_part2_stride.sum()

        np_part2_data = create_numpy(_part2_data[i_part], npy_type, dim_np)

        lnp_part_strid.append(np_part2_stride)
        lnp_part_data .append(np_part2_data)

      else:
        cst_stride = requested[2]
        dim_np  = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]] * cst_stride
        np_part2_data = create_numpy(_part2_data[i_part], npy_type, dim_np)

        lnp_part_data .append(np_part2_data)

    free(_part2_stride)
    free(_part2_data)


    return lnp_part_strid, lnp_part_data


  # --------------------------------------------------------------------
  def reverse_iexch(self,
                    PDM_mpi_comm_kind_t         k_comm,
                    PDM_part_to_part_data_def_t t_part2_data_def,
                    list                        part2_data,
                    part2_stride=1,
                    bint interlaced_str=True):
    """
    reverse_iexch(k_comm, t_part2_data_def, part2_data, part2_stride=1, interlaced_str=True)

    Initiate a non-blocking exchange (Part2→Part1)

    Parameters:
      k_comm           (int)                       : Kind of MPI communication
      t_part2_data_def (int)                       : Kind of Part2 data definition
      part2_data       (list)                      : Part2 data
      part2_stride     (`int` or `list, optional`) : Stride of Part2 data
      interlaced_str   (bool, optional)            : Is the data interlaced? (default = **True**)

    Returns:
      Request ID (`int`)

    Admissible values for ``k_comm``:
      - 0 : Peer-to-peer (``MPI_issend``/``MPI_irecv``)

    .. note::
      Additional communication kinds will be available in the future

    Admissible values for ``t_part2_data_def``:
      - :py:attr:`PartToPart.DATA_DEF_ORDER_PART2`           : Data defined according to the `part2` arrays order
      - :py:attr:`PartToPart.DATA_DEF_ORDER_GNUM1_COME_FROM` : Data defined according to the `gnum1_come_from` arrays order
    """
    cdef int          request_exch
    cdef PDM_stride_t _stride_t

    cdef int   _part2_stride_cst = 0
    cdef int** _part2_stride = NULL

    # To check stride size
    cdef int   n_elt       = 0;
    cdef int  *n_ref_lnum2 = NULL;
    cdef int **ref_lnum2   = NULL;

    # To check stride size (gnum_come_from case)
    cdef int         **gnum1_come_from_idx = NULL;
    cdef PDM_g_num_t **gnum1_come_from     = NULL;

    cdef int  n_part2
    cdef int* n_elt2
    PDM_part_to_part_n_part_and_n_elt_get(self.ptp, NULL, &n_part2, NULL, &n_elt2)

    if isinstance(part2_stride, int):
      _stride_t = PDM_STRIDE_CST_INTERLACED if interlaced_str else PDM_STRIDE_CST_INTERLEAVED
      _part2_stride_cst = part2_stride

    elif isinstance(part2_stride, list):
      _stride_t = PDM_STRIDE_VAR_INTERLACED
      assert len(part2_stride) == n_part2

      PDM_part_to_part_ref_lnum2_get(self.ptp,
                                    &n_ref_lnum2,
                                    &ref_lnum2);

      for i_part in range(n_part2):
        if (t_part2_data_def==PDM_PART_TO_PART_DATA_DEF_ORDER_PART1_TO_PART2):
          n_elt = n_ref_lnum2[i_part]
        elif(t_part2_data_def==PDM_PART_TO_PART_DATA_DEF_ORDER_GNUM1_COME_FROM):
          PDM_part_to_part_gnum1_come_from_get(self.ptp,
                                              &gnum1_come_from_idx,
                                              &gnum1_come_from);
          n_elt = gnum1_come_from_idx[i_part][n_ref_lnum2[i_part]]
        else :
          n_elt = n_elt2[i_part]
        assert_single_dim_np(part2_stride[i_part], NPY.int32, n_elt)
      _part2_stride = np_list_to_int_pointers(part2_stride)
    else:
      raise ValueError("Invalid stride in PtB exchange")

    cdef void** _part2_data   = np_list_to_void_pointers(part2_data)

    cdef PDM_MPI_Comm pdm_comm = PDM_part_to_part_comm_get(self.ptp)
    py_comm = pdm_comm_to_py_comm(pdm_comm)

    ref_dtype = recover_dtype(part2_data, py_comm)
    cdef size_t s_data   = ref_dtype.itemsize
    cdef size_t npy_type = ref_dtype.num

    cdef int**  _part1_stride = NULL;
    cdef void** _part1_data   = NULL;
    PDM_part_to_part_reverse_iexch(self.ptp,
                                  k_comm,
                                  _stride_t,
                                  t_part2_data_def,
                                  _part2_stride_cst,
                                  s_data,
                  <const int ** >  _part2_stride,
                  <const void **>  _part2_data,
                                  &_part1_stride,
                        <void ***> &_part1_data,
                                  &request_exch)

    self.request_data[request_exch] = [<uintptr_t> _part1_stride, <uintptr_t> _part1_data, _part2_stride_cst, npy_type]

    if _stride_t == PDM_STRIDE_VAR_INTERLACED:
      free(_part2_stride)
    free(_part2_data)

    return request_exch

  # --------------------------------------------------------------------
  def reverse_wait(self, int request_id):
    """
    reverse_wait(request_id)

    Finalize a non-blocking exchange (Part2→Part1)

    Parameters:
      request_id (int) : Request ID

    Returns:
      - Part1 stride (same dtype as ``part2_stride`` in :py:func:`reverse_iexch`)
      - Part1 data   (`list` of same dtype as ``part2_data`` in :py:func:`reverse_iexch`)
    """
    PDM_part_to_part_reverse_iexch_wait(self.ptp, request_id)

    requested = self.request_data.pop(request_id)
    cdef uintptr_t _part1_stride_id = requested[0]
    cdef uintptr_t _part1_data_id   = requested[1]
    cdef int**  _part1_stride = <int**> _part1_stride_id
    cdef void** _part1_data = <void **> _part1_data_id
    cdef NPY.npy_intp dim_np

    cdef int   n_part1, n_part2
    cdef int*  n_elt1
    cdef int** part1_to_part2_idx
    PDM_part_to_part_n_part_get(self.ptp, &n_part1, &n_part2)
    PDM_part_to_part_part1_to_part2_idx_get(self.ptp, &n_elt1, &part1_to_part2_idx)

    lnp_part_strid = list()
    lnp_part_data  = list()
    cdef size_t npy_type = requested[3]


    for i_part in range(n_part1):

      if(_part1_stride != NULL):

        strid_size = n_elt1[i_part]

        np_part1_stride = create_numpy_i(_part1_stride[i_part], strid_size)
        dim_np = np_part1_stride.sum()
        # print("dim_np : ", dim_np)
        # print("np_part1_stride : ", np_part1_stride)

        np_part1_data = create_numpy(_part1_data[i_part], npy_type, dim_np)

        lnp_part_strid.append(np_part1_stride)
        lnp_part_data .append(np_part1_data)
        free(_part1_data)
        free(_part1_stride)

      else:
        cst_stride = requested[2]
        dim_np  = part1_to_part2_idx[i_part][n_elt1[i_part]] * cst_stride


        np_part1_data = create_numpy(_part1_data[i_part], npy_type, dim_np)

        lnp_part_data .append(np_part1_data)

        free(_part1_data)

        if(_part1_stride != NULL):
          free(_part1_stride)


    return lnp_part_strid, lnp_part_data

  # --------------------------------------------------------------------
  def __dealloc__(self):
    PDM_part_to_part_free(self.ptp)

  # --------------------------------------------------------------------

# ------------------------------------------------------------------------
# ========================================================================

