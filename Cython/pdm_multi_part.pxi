
cdef extern from "pdm_multipart.h":

    # -> PPART bases functions
    # ------------------------------------------------------------------
    # MPI_Comm      comm,
    int PDM_multipart_create(int              n_zone,
                             int              n_part,
                             PDM_bool_t       merge_blocks,
                             PDM_part_split_t split_method,
                             PDM_MPI_Comm     comm)

    # ------------------------------------------------------------------
    void PDM_multipart_register_block(int        mpart_id,
                                      int        zoneGId,
                                      int        block_data_id)

    # ------------------------------------------------------------------
    void PDM_multipart_run_ppart(int id);

    # ------------------------------------------------------------------
    void PDM_multipart_part_dim_get(int  mpartId,
                                    int  zoneGId,
                                    int  ipart,
                                    int *nCell,
                                    int *nFace,
                                    int *nFacePartBound,
                                    int *nVtx,
                                    int *nProc,
                                    int *nTPart,
                                    int *sCellFace,
                                    int *sFaceVtx,
                                    int *sFaceGroup,
                                    int *nFaceGroup);

    # ------------------------------------------------------------------
    void PDM_multipart_part_val_get(int            mpartId,
                                    int            zoneGId,
                                    int            ipart,
                                    int          **cellTag,
                                    int          **cellFaceIdx,
                                    int          **cellFace,
                                    PDM_g_num_t  **cellLNToGN,
                                    int          **faceTag,
                                    int          **faceCell,
                                    int          **faceVtxIdx,
                                    int          **faceVtx,
                                    PDM_g_num_t  **faceLNToGN,
                                    int          **facePartBoundProcIdx,
                                    int          **facePartBoundPartIdx,
                                    int          **facePartBound,
                                    int          **vtxTag,
                                    double       **vtx,
                                    PDM_g_num_t  **vtxLNToGN,
                                    int          **faceGroupIdx,
                                    int          **faceGroup,
                                    PDM_g_num_t  **faceGroupLNToGN)
    # ------------------------------------------------------------------
    void PDM_multipart_part_color_get(int            mpartId,
                                      int            zoneGId,
                                      int            ipart,
                                      int          **cellColor,
                                      int          **faceColor,
                                      int          **threadColor,
                                      int          **hyperPlaneColor)

    # ------------------------------------------------------------------
    void PDM_multipart_time_get(int       mpartId,
                                int       zoneGId,
                                double  **elapsed,
                                double  **cpu,
                                double  **cpu_user,
                                double  **cpu_sys);
    # ------------------------------------------------------------------
    void PDM_multipart_free(int id);

# ------------------------------------------------------------------
cdef class MultiPart:
    """
       MultiPpart
    """
    # > For Ppart
    cdef int _mpart_id
    # ------------------------------------------------------------------
    def __cinit__(self,
                  int              n_zone,
                  int              n_part,
                  PDM_bool_t       merge_blocks,
                  PDM_part_split_t split_method,
                  MPI.Comm         comm):

        """
        """
        # ~> Communicator Mpi
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi

        # -> Create PPART
        # DM_multipart_create(&_mpart_id,
        self._mpart_id = PDM_multipart_create(n_zone,
                                              n_part,
                                              merge_blocks,
                                              split_method,
                                              PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))

    # ------------------------------------------------------------------
    # def multi_part_partial_free(self):
    #     print 'multi_part_partial_free MultiPart'
    #     PDM_multipart_partial_free(self._mpart_id)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        print '__dealloc__ MultiPart '
        PDM_multipart_free(self._mpart_id)

    # ------------------------------------------------------------------
    def multipart_register_block(self, int zoneGId,
                                       int block_data_id):
        """
        """
        PDM_multipart_register_block(self._mpart_id, zoneGId, block_data_id)

    # ------------------------------------------------------------------
    def multipart_run_ppart(self):
        """
        """
        PDM_multipart_run_ppart(self._mpart_id)

    # ------------------------------------------------------------------
    def multipart_dim_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int nCell
        cdef int nProc
        cdef int nTPart
        cdef int nFace
        cdef int nFacePartBound
        cdef int nVertex
        cdef int sCellFace
        cdef int sFaceVertex
        cdef int sFaceGroup
        cdef int nFaceGroup
        # ************************************************************************

        PDM_multipart_part_dim_get(self._mpart_id,
                                   zoneGId,
                                   ipart,
                                   &nCell,
                                   &nFace,
                                   &nFacePartBound,
                                   &nVertex,
                                   &nProc,
                                   &nTPart,
                                   &sCellFace,
                                   &sFaceVertex,
                                   &sFaceGroup,
                                   &nFaceGroup)

        return {'nCell'          :nCell,
                'ipart'          :ipart,
                'nFace'          :nFace,
                'nTPart'         :nTPart,
                'nProc'          :nProc,
                'nFacePartBound' :nFacePartBound,
                'nVertex'        :nVertex,
                'sCellFace'      :sCellFace,
                'sFaceVertex'    :sFaceVertex,
                'sFaceGroup'     :sFaceGroup,
                'nFaceGroup'     :nFaceGroup}

    # ------------------------------------------------------------------
    def mutlipart_val_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellTag,
        cdef int          *cellFaceIdx,
        cdef int          *cellFace,
        cdef PDM_g_num_t  *cellLNToGN,
        cdef int          *faceTag,
        cdef int          *faceCell,
        cdef int          *faceVertexIdx,
        cdef int          *faceVertex,
        cdef PDM_g_num_t  *faceLNToGN,
        cdef int          *facePartBound,
        cdef int          *facePartBoundProcIdx,
        cdef int          *facePartBoundPartIdx,
        cdef int          *vertexTag,
        cdef double       *vertex,
        cdef PDM_g_num_t  *vertexLNToGN,
        cdef int          *faceGroupIdx,
        cdef int          *faceGroup,
        cdef PDM_g_num_t  *faceGroupLNToGN
        # ************************************************************************

        # dims = self.part_dim_get(self._mpart_id, ipart)
        dims = self.multipart_dim_get(ipart, zoneGId)

        # -> Call PPART to get info
        PDM_multipart_part_val_get(self._mpart_id,
                                   zoneGId,
                                   ipart,
                                   &cellTag,
                                   &cellFaceIdx,
                                   &cellFace,
                                   &cellLNToGN,
                                   &faceTag,
                                   &faceCell,
                                   &faceVertexIdx,
                                   &faceVertex,
                                   &faceLNToGN,
                                   &facePartBoundProcIdx,
                                   &facePartBoundPartIdx,
                                   &facePartBound,
                                   &vertexTag,
                                   &vertex,
                                   &vertexLNToGN,
                                   &faceGroupIdx,
                                   &faceGroup,
                                   &faceGroupLNToGN)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cellTag            Cell tag (size = nCell)
        if (cellTag == NULL) :
            npCellTag = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellTag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> cellTag)

        # \param [out]  cellFaceIdx        Cell to face connectivity index (size = nCell + 1)
        if (cellFaceIdx == NULL) :
            npCellFaceIdx = None
        else :
            dim = <NPY.npy_intp> (dims['nCell'] + 1)
            npCellFaceIdx = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> cellFaceIdx)

        # \param [out]  cellFace           Cell to face connectivity (size = cellFaceIdx[nCell] = lCellFace)
        if (cellFace == NULL) :
            npCellFace = None
        else :
            dim = <NPY.npy_intp> dims['sCellFace']
            npCellFace = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> cellFace)

        # \param [out]  cellLNToGN         Cell local numbering to global numbering (size = nCell)
        # dim = <NPY.npy_intp> dims['nCell']
        if (cellLNToGN == NULL) :
            npCellLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        PDM_G_NUM_NPY_INT,
                                                        <void *> cellLNToGN)

        # \param [out]  faceTag            Face tag (size = nFace)
        if (faceTag == NULL) :
            npFaceTag = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceTag = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_INT32,
                                                     <void *> faceTag)

        # \param [out]  faceCell           Face to cell connectivity  (size = 2 * nFace)
        if (faceCell == NULL) :
            npFaceCell = None
        else :
            dim = <NPY.npy_intp> (2 * dims['nFace'])
            npFaceCell = NPY.PyArray_SimpleNewFromData(1,
                                                      &dim,
                                                      NPY.NPY_INT32,
                                                      <void *> faceCell)

        # \param [out]  faceVtxIdx         Face to Vertex connectivity index (size = nFace + 1)
        if (faceVertexIdx == NULL) :
            npFaceVertexIdx = None
        else :
            dim = <NPY.npy_intp> (dims['nFace'] + 1)
            npFaceVertexIdx = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> faceVertexIdx)

        # \param [out]  faceVtx            Face to Vertex connectivity (size = faceVtxIdx[nFace])
        cdef NPY.ndarray[NPY.int32_t, ndim=1] npFaceVertex
        if (faceVertex == NULL) :
            npFaceVertex = None
        else :
            dim = <NPY.npy_intp> dims['sFaceVertex']
            npFaceVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> faceVertex)
            # PyArray_ENABLEFLAGS(npFaceVertex, NPY.NPY_OWNDATA)
            # print '*'*1000
            # print 'Take ownership'
            # print '*'*1000

        # \param [out]  faceLNToGN         Face local numbering to global numbering (size = nFace)
        if (faceLNToGN == NULL) :
            npFaceLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceLNToGN   = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          PDM_G_NUM_NPY_INT,
                                                          <void *> faceLNToGN)

        # \param [out]  facePartBound      Partitioning boundary faces
        if (facePartBound == NULL) :
            npFacePartBound = None
        else :
            dim = <NPY.npy_intp> (4 * dims['nFacePartBound'])
            npFacePartBound   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBound)

        # \param [out]  facePartBoundProcIdx  Partitioning boundary faces block distribution from processus (size = nProc + 1)
        if (facePartBoundProcIdx == NULL) :
            npfacePartBoundProcIdx = None
        else :
            dim = <NPY.npy_intp> ( dims['nProc'] + 1)
            npfacePartBoundProcIdx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBoundProcIdx)

        # \param [out]  facePartBoundPartIdx  Partitioning boundary faces block distribution from partition (size = nTPart + 1)
        if (facePartBoundPartIdx == NULL) :
            npfacePartBoundPartIdx = None
        else :
            dim = <NPY.npy_intp> ( dims['nTPart'] + 1)
            npfacePartBoundPartIdx   = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             NPY.NPY_INT32,
                                                             <void *> facePartBoundPartIdx)

        # \param [out]  vtxTag             Vertex tag (size = nVtx)
        if (vertexTag == NULL) :
            npVertexTag = None
        else :
            dim = <NPY.npy_intp> dims['nVertex']
            npVertexTag   = NPY.PyArray_SimpleNewFromData(1,
                                                         &dim,
                                                         NPY.NPY_INT32,
                                                         <void *> vertexTag)

        # \param [out]  vtx                Vertex coordinates (size = 3 * nVtx)
        if (vertex == NULL) :
            npVertex = None
        else :
            dim = <NPY.npy_intp> (3 * dims['nVertex'])
            npVertex  = NPY.PyArray_SimpleNewFromData(1,
                                                     &dim,
                                                     NPY.NPY_DOUBLE,
                                                     <void *> vertex)

        # \param [out]  vtxLNToGN          Vertex local numbering to global numbering (size = nVtx)
        if (vertexLNToGN == NULL) :
            npVertexLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['nVertex']
            npVertexLNToGN  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           PDM_G_NUM_NPY_INT,
                                                           <void *> vertexLNToGN)

        # \param [out]  faceGroupIdx       face group index (size = nFaceGroup + 1)
        if (faceGroupIdx == NULL) :
            npFaceGroupIdx = None
        else :
            dim = <NPY.npy_intp> (self._nFaceGroup + 1)
            npFaceGroupIdx  = NPY.PyArray_SimpleNewFromData(1,
                                                           &dim,
                                                           NPY.NPY_INT32,
                                                           <void *> faceGroupIdx)

        # \param [out]  faceGroup          faces for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
        if (faceGroup == NULL) :
            npFaceGroup = None
        else :
            dim = <NPY.npy_intp> dims['sFaceGroup']
            npFaceGroup = NPY.PyArray_SimpleNewFromData(1,
                                                       &dim,
                                                       NPY.NPY_INT32,
                                                       <void *> faceGroup)

        # \param [out]  faceGroupLNToGN    faces global numbering for each group (size = faceGroupIdx[nFaceGroup] = lFaceGroup)
        if (faceGroupLNToGN == NULL) :
            npFaceGroupLNToGN = None
        else :
            dim = <NPY.npy_intp> dims['sFaceGroup']
            npFaceGroupLNToGN = NPY.PyArray_SimpleNewFromData(1,
                                                             &dim,
                                                             PDM_G_NUM_NPY_INT,
                                                             <void *> faceGroupLNToGN)

        return {'npCellTag'                  : npCellTag,
                'npCellFaceIdx'              : npCellFaceIdx,
                'npCellFace'                 : npCellFace,
                'npCellLNToGN'               : npCellLNToGN,
                'npFaceTag'                  : npFaceTag,
                'npFaceCell'                 : npFaceCell,
                'npFaceVertexIdx'            : npFaceVertexIdx,
                'npFaceVertex'               : npFaceVertex,
                'npFaceLNToGN'               : npFaceLNToGN,
                'npfacePartBoundProcIdx'     : npfacePartBoundProcIdx,
                'npfacePartBoundPartIdx'     : npfacePartBoundPartIdx,
                'npFacePartBound'            : npFacePartBound,
                'npVertexTag'                : npVertexTag,
                'npVertex'                   : npVertex,
                'npVertexLNToGN'             : npVertexLNToGN,
                'npFaceGroupIdx'             : npFaceGroupIdx,
                'npFaceGroup'                : npFaceGroup,
                'npFaceGroupLNToGN'          : npFaceGroupLNToGN}

    # ------------------------------------------------------------------
    def multipart_color_get(self, int ipart, int zoneGId):
        """
           Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cellColor,
        cdef int          *faceColor
        cdef int          *threadColor
        cdef int          *hyperPlaneColor
        # ************************************************************************

        # dims = self.part_dim_get(self._mpart_id, ipart)
        dims = self.multipart_dim_get(ipart)

        # -> Call PPART to get info
        PDM_multipart_part_color_get(self._mpart_id,
                                     zoneGId,
                                     ipart,
                                     &cellColor,
                                     &faceColor,
                                     &threadColor,
                                     &hyperPlaneColor)
        # -> Begin
        cdef NPY.npy_intp dim

        # \param [out]  cellColor            Cell tag (size = nCell)
        if (cellColor == NULL):
            npCellColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npCellColor = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> cellColor)
        # \param [out]  faceColor            Cell tag (size = nFace)
        if (faceColor == NULL):
            npFaceColor = None
        else :
            dim = <NPY.npy_intp> dims['nFace']
            npFaceColor = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_INT32,
                                                        <void *> faceColor)

        # \param [out]  threadColor            Cell tag (size = nCell)
        if (threadColor == NULL):
            npThreadColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npThreadColor = NPY.PyArray_SimpleNewFromData(1,
                                                          &dim,
                                                          NPY.NPY_INT32,
                                                          <void *> threadColor)

        # \param [out]  hyperPlaneColor            Cell tag (size = nCell)
        if (hyperPlaneColor == NULL):
            npHyperPlaneColor = None
        else :
            dim = <NPY.npy_intp> dims['nCell']
            npHyperPlaneColor = NPY.PyArray_SimpleNewFromData(1,
                                                              &dim,
                                                              NPY.NPY_INT32,
                                                              <void *> hyperPlaneColor)
        return {'npCellColor'       : npCellColor,
                'npFaceColor'       : npFaceColor,
                'npThreadColor'     : npThreadColor,
                'npHyperPlaneColor' : npHyperPlaneColor}

    # ------------------------------------------------------------------
    def multipart_time_get(self, int zoneGId):
        """
        Get times
        """
        # ************************************************************************
        # > Declaration
        cdef double *elapsed
        cdef double *cpu
        cdef double *cpu_user
        cdef double *cpu_sys
        # ************************************************************************

        PDM_multipart_time_get(self._mpart_id, zoneGId, &elapsed, &cpu, &cpu_user, &cpu_sys)

        d_elapsed = {'total'              : elapsed[0],
                     'building graph'     : elapsed[1],
                     'splitting graph'    : elapsed[2],
                     'building partitions': elapsed[3]}

        d_cpu     = {'total'              : cpu[0],
                     'building graph'     : cpu[1],
                     'splitting graph'    : cpu[2],
                     'building partitions': cpu[3]}

        d_cpu_user = {'total'              : cpu_user[0],
                      'building graph'     : cpu_user[1],
                      'splitting graph'    : cpu_user[2],
                      'building partitions': cpu_user[3]}

        d_cpu_sys = {'total'              : cpu_sys[0],
                     'building graph'     : cpu_sys[1],
                     'splitting graph'    : cpu_sys[2],
                     'building partitions': cpu_sys[3]}

        return {'elapsed' : d_elapsed, 'cpu' : d_cpu, 'cpu_user' : d_cpu_user,  'cpu_sys' : d_cpu_sys}


    # ------------------------------------------------------------------
    # def multipart_stat_get(self, int zoneGId):
    #     """
    #     Get statistics
    #     """
    #     # ************************************************************************
    #     # > Declaration
    #     cdef int      cells_average,
    #     cdef int      cells_median,
    #     cdef double   cells_std_deviation,
    #     cdef int      cells_min,
    #     cdef int      cells_max,
    #     cdef int      bound_part_faces_average,
    #     cdef int      bound_part_faces_median,
    #     cdef double   bound_part_faces_std_deviation,
    #     cdef int      bound_part_faces_min,
    #     cdef int      bound_part_faces_max,
    #     cdef int      bound_part_faces_sum
    #     # ************************************************************************

    #     PDM_multipart_stat_get(self._mpart_id,
    #                            zoneGId,
    #                            &cells_average,
    #                            &cells_median,
    #                            &cells_std_deviation,
    #                            &cells_min,
    #                            &cells_max,
    #                            &bound_part_faces_average,
    #                            &bound_part_faces_median,
    #                            &bound_part_faces_std_deviation,
    #                            &bound_part_faces_min,
    #                            &bound_part_faces_max,
    #                            &bound_part_faces_sum)

    #     return {'cells_average'                  : cells_average,
    #             'cells_median'                   : cells_median,
    #             'cells_std_deviation'            : cells_std_deviation,
    #             'cells_min'                      : cells_min,
    #             'cells_max'                      : cells_max,
    #             'bound_part_faces_average'       : bound_part_faces_average,
    #             'bound_part_faces_median'        : bound_part_faces_median,
    #             'bound_part_faces_std_deviation' : bound_part_faces_std_deviation,
    #             'bound_part_faces_min'           : bound_part_faces_min,
    #             'bound_part_faces_max'           : bound_part_faces_max,
    #             'bound_part_faces_sum'           : bound_part_faces_sum}
