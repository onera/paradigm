
cdef extern from "pdm_part_to_block.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    int PDM_distant_neighbor_create(PDM_MPI_Comm   comm,
                                    int            n_part,
                                    int           *n_entity,
                                    int          **neighbor_idx,
                                    int          **neighbor_desc)

    void PDM_distant_neighbor_free(int id)

    void PDM_distant_neighbor_exch(int      id,
                                   size_t         s_data,
                                   PDM_stride_t   t_stride,
                                   int            cst_stride,
                                   int          **send_entity_stride,
                                   void         **send_entity_data,
                                   int         ***recv_entity_stride,
                                   void        ***recv_entity_data);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class DistantNeighbor:
    """
       DistantNeighbor: Interface for pdm_distant_neighbor.c
    """
    # ************************************************************************
    # > Class attributes
    cdef int   _dnid
    cdef int   _n_part
    cdef int  *_n_entity
    cdef int **_neighbor_idx
    cdef int **_neighbor_desc
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm,
                        int      n_part,
                        list     neighbor_idx,
                        list     neighbor_desc):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] pneighbor_idx
        cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='fortran'] pneighbor_desc
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Some verification
        assert(len(neighbor_idx) == len(neighbor_desc) == n_part)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Store partN and parameter
        self._n_part        = n_part
        self._neighbor_idx  = <int **> malloc(self._n_part * sizeof(int **))
        self._neighbor_desc = <int **> malloc(self._n_part * sizeof(int **))
        self._n_entity      = <int *>  malloc(self._n_part * sizeof(int * ))
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Prepare
        for i_part in range(self._n_part):

          # ------------------------------------------------
          # > Get shape of array
          n_entity       = neighbor_idx[i_part].shape[0]-1
          pneighbor_idx  = neighbor_idx[i_part]
          pneighbor_desc = neighbor_desc[i_part]

          self._n_entity[i_part]      = <int  > n_entity
          self._neighbor_idx[i_part]  = <int *> pneighbor_idx.data
          self._neighbor_desc[i_part] = <int *> pneighbor_desc.data
          # ------------------------------------------------

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Create
        self._dnid = PDM_distant_neighbor_create(PDMC,
                                                 n_part,
                                                 self._n_entity,
                                                 self._neighbor_idx,
                                                 self._neighbor_desc)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    # def DistantNeighbor_Exchange(self, dict dField, dict pField, list pStrid = None):
    #     """
    #        TODOUX : 1) Exchange of variables types array
    #                 2) Assertion of type and accross MPI of the same field
    #     """
    #     # ************************************************************************
    #     # > Declaration
    #     cdef NPY.ndarray   dArray
    #     cdef NPY.ndarray   pArray
    #     cdef NPY.ndarray   partStrid
    #     cdef int           idx
    #     cdef int           blk_size
    #     cdef NPY.npy_intp *ArrayDim

    #     # > For PDM
    #     cdef size_t   s_data
    #     cdef int      strideOne
    #     cdef int    **part_stride
    #     cdef int     *block_stride
    #     cdef void    *block_data
    #     cdef void   **part_data
    #     cdef int      ndim
    #     cdef int      npyflags=-1;
    #     cdef NPY.ndarray tmpData
    #     # ************************************************************************

    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #     # > Allocate
    #     part_data = <void **> malloc(self.partN * sizeof(void **))
    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #     # > Prepare stride
    #     if(self.t_stride == 0): # Cst Stride
    #        strideOne     = 1
    #        part_stride   = NULL
    #     else:
    #        strideOne     = 0
    #        part_stride   = <int **> malloc(self.partN * sizeof(int *))
    #        assert(pStrid is not None)
    #        assert(len(pStrid) == self.partN)
    #        # for idx in xrange(self.partN):
    #        for idx, partStrid in enumerate(pStrid):
    #           part_stride[idx] = <int *> partStrid.data
    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #     ndim = -1

    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #     # > Loop on all field to build
    #     for field, partList in pField.iteritems():

    #       # print field, partList

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       # > Prepare part_data
    #       for idx, pArray in enumerate(partList):
    #         # ------------------------------------------------
    #         # > Get flow solution
    #         if(pArray.ndim == 2):
    #           assert(pArray.shape[1] == self.NbElmts[idx])
    #           ndim = 2
    #         else:
    #           if(self.t_stride == <PDM_stride_t   >(0)):
    #             assert(pArray.shape[0] == self.NbElmts[idx])
    #           ndim = 1
    #         part_data[idx] = <void *> pArray.data
    #         # ------------------------------------------------

    #         # ------------------------------------------------
    #         # > Fill s_data - How to check if array is different ?
    #         s_data     = pArray.dtype.itemsize
    #         dtype_data = pArray.dtype.num
    #         dtypep     = pArray.dtype
    #         # ------------------------------------------------

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       # > Prepare block_data
    #       block_stride  = NULL
    #       block_data    = NULL
    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       # > Exchange
    #       # print "PDM_part_to_block_exch "
    #       c_size = PDM_part_to_block_exch(self.PTB,
    #                                       s_data,
    #                                       self.t_stride,
    #                                       strideOne,
    #                                       part_stride,
    #                                       part_data,
    #                                       &block_stride,
    #                                       <void **> &block_data)
    #       # print "PDM_part_to_block_exch end "
    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       blk_size = PDM_part_to_block_n_elt_block_get(self.PTB);
    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       # > Put in dict
    #       dim           = <NPY.npy_intp> c_size
    #       if(c_size == 0):
    #         # dField[field] = None # Attention faire une caspule vide serait mieux non ?
    #         dField[field] = NPY.empty((0), dtype=dtypep)
    #         # print 'Attention in PDM_part_to_block'
    #       else:
    #         if(ndim == 2):
    #           ArrayDim    = <NPY.npy_intp *> malloc(2 * sizeof(NPY.npy_intp *))
    #           ArrayDim[0] = <NPY.npy_intp> 1
    #           ArrayDim[1] = <NPY.npy_intp> c_size

    #           # > Put in dField
    #           tmpData = NPY.PyArray_SimpleNewFromData(ndim, ArrayDim, dtype_data, <void *> block_data)
    #           PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
    #           dField[field] = tmpData
    #           # > Free
    #           free(ArrayDim)
    #         else:
    #           # dField[field] = NPY.PyArray_SimpleNewFromData(1, &dim, dtype_data,
    #           #                                               <void *> block_data)
    #           tmpData = NPY.PyArray_SimpleNewFromData(1, &dim, dtype_data,
    #                                                      <void *> block_data)
    #           PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
    #           dField[field] = tmpData
    #           # print(dField[field].flags)
    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #       # > Stride management
    #       if(self.t_stride == 1 ):
    #         dimStri = <NPY.npy_intp> blk_size
    #         # dField[field+'#Stride'] = NPY.PyArray_SimpleNewFromData(1, &dimStri, NPY.NPY_INT32, <void *> block_stride)
    #         tmpData = NPY.PyArray_SimpleNewFromData(1, &dimStri, NPY.NPY_INT32, <void *> block_stride)
    #         PyArray_ENABLEFLAGS(tmpData, NPY.NPY_OWNDATA);
    #         dField[field+'#Stride'] = tmpData
    #       # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #     # > Deallocate
    #     free(part_data)
    #     if(self.t_stride != 0): # Var Stride
    #       free(part_stride)
    #     # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # # ------------------------------------------------------------------------
    # def getDistributionCopy(self):
    #   """
    #      Return a copy of the distrisbution array compute in library
    #      Copy because remove of PTB object can made a core ...
    #   """
    #   # ************************************************************************
    #   # > Declaration
    #   cdef PDM_g_num_t* Distrib
    #   # ************************************************************************

    #   # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #   # > Get
    #   Distrib = PDM_part_to_block_distrib_index_get(self.PTB)
    #   # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #   # ::::::::::::::::::::::::::::::::::::::::::::::::::
    #   dim        = <NPY.npy_intp> (self.Size+1)
    #   DistribNPY = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, <void *> Distrib)
    #   # ::::::::::::::::::::::::::::::::::::::::::::::::::

    #   return NPY.copy(DistribNPY)


    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # ************************************************************************

      # > Free Ppart Structure
      PDM_distant_neighbor_free(self._dnid)

      free(self._n_entity)
      free(self._neighbor_idx)
      free(self._neighbor_desc)

