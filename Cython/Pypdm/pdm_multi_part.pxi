
cdef extern from "pdm_multipart.h":
    ctypedef enum PDM_part_size_t:
        PDM_PART_SIZE_HOMONEGEOUS   = 1
        PDM_PART_SIZE_HETEROGENEOUS = 2
    ctypedef struct PDM_multipart_t:
        pass

    # -> PPART bases functions
    # ------------------------------------------------------------------
    # MPI_Comm      comm,
    PDM_multipart_t* PDM_multipart_create(int              n_zone,
                                          int*             n_part,
                                          PDM_bool_t       merge_blocks,
                                          PDM_split_dual_t split_method,
                                          PDM_part_size_t  part_size_method,
                                          double*          part_fraction,
                                          PDM_MPI_Comm     comm,
                                          PDM_ownership_t  owner)

    # ------------------------------------------------------------------
    void PDM_multipart_register_block(PDM_multipart_t *mtp,
                                      int              i_zone,
                                      PDM_dmesh_t     *dmesh)

    # ------------------------------------------------------------------
    void PDM_multipart_register_dmesh_nodal(PDM_multipart_t   *mtp,
                                            int                i_zone,
                                            PDM_dmesh_nodal_t *dmesh)

    # ------------------------------------------------------------------
    void PDM_multipart_register_joins(PDM_multipart_t *mtp,
                                      int              n_total_joins,
                                      int*             matching_join_array)

    # ------------------------------------------------------------------
    void PDM_multipart_set_reordering_options(      PDM_multipart_t *mtp,
                                              const int              i_zone,
                                              const char            *renum_cell_method,
                                              const int             *renum_cell_properties,
                                              const char            *renum_face_method)
    void PDM_multipart_set_reordering_options_vtx(      PDM_multipart_t *mtp,
                                                  const int              i_zone,
                                                  const char            *renum_vtx_method)

    # ------------------------------------------------------------------
    void PDM_multipart_run_ppart(PDM_multipart_t *mtp)

    # ------------------------------------------------------------------
    void PDM_multipart_part_dim_get(PDM_multipart_t  *mtp,
                                    int               i_zone,
                                    int               ipart,
                                    int              *n_cell,
                                    int              *n_face,
                                    int              *n_face_part_bound,
                                    int              *nVtx,
                                    int              *n_proc,
                                    int              *nt_part,
                                    int              *scell_face,
                                    int              *sface_vtx,
                                    int              *s_face_bound,
                                    int              *n_face_bound)

    # ------------------------------------------------------------------
    void PDM_multipart_part_val_get(PDM_multipart_t   *mtp,
                                    int                i_zone,
                                    int                ipart,
                                    int              **cell_face_idx,
                                    int              **cell_face,
                                    PDM_g_num_t      **cell_ln_to_gn,
                                    int              **face_cell,
                                    int              **face_vtx_idx,
                                    int              **face_vtx,
                                    PDM_g_num_t      **face_ln_to_gn,
                                    int              **face_part_bound_proc_idx,
                                    int              **face_part_bound_part_idx,
                                    int              **face_part_bound,
                                    double           **vtx,
                                    PDM_g_num_t      **vtx_ln_to_gn,
                                    int              **face_bound_idx,
                                    int              **face_bound,
                                    PDM_g_num_t      **face_bound_ln_to_gn)

    # ------------------------------------------------------------------
    int PDM_multipart_part_vtx_coord_get(PDM_multipart_t  *multipart,
                                         const int         i_zone,
                                         const int         i_part,
                                         double          **vtx_coord,
                                         PDM_ownership_t   ownership)

    # ------------------------------------------------------------------
    void PDM_multipart_get_part_mesh_nodal(PDM_multipart_t   *mtp,
                                           const int   i_zone,
                                           PDM_part_mesh_nodal_t **pmesh_nodal,
                                           PDM_ownership_t         ownership)

    # ------------------------------------------------------------------
    void PDM_multipart_part_graph_comm_get(PDM_multipart_t    *multipart,
                                           const int           i_zone,
                                           const int           i_part,
                                           PDM_bound_type_t    bound_type,
                                           int               **ppart_bound_proc_idx,
                                           int               **ppart_bound_part_idx,
                                           int               **ppart_bound,
                                           PDM_ownership_t     ownership);
    # ------------------------------------------------------------------
    void PDM_multipart_part_hyperplane_color_get(PDM_multipart_t  *multipart,
                                                 int               i_zone,
                                                 int               i_part,
                                                 int             **hyperplane_color,
                                                 PDM_ownership_t   ownership);
    # ------------------------------------------------------------------
    void PDM_multipart_part_thread_color_get(PDM_multipart_t  *multipart,
                                                 int               i_zone,
                                                 int               i_part,
                                                 int             **thread_color,
                                                 PDM_ownership_t   ownership);

    # ------------------------------------------------------------------
    void PDM_multipart_part_ghost_infomation_get(PDM_multipart_t *mtp,
                                                 int              i_zone,
                                                 int              ipart,
                                                 int             **vtx_ghost_information,
                                                 PDM_ownership_t   ownership)

    # ------------------------------------------------------------------
    int PDM_multipart_part_connectivity_get(PDM_multipart_t          *mtp,
                                            int                       i_zone,
                                            int                       i_part,
                                            PDM_connectivity_type_t   connectivity_type,
                                            int                     **connect,
                                            int                     **connect_idx,
                                            PDM_ownership_t           ownership);
    # ------------------------------------------------------------------
    int PDM_multipart_part_ln_to_gn_get(PDM_multipart_t      *mtp,
                                        int                   i_zone,
                                        int                   i_part,
                                        PDM_mesh_entities_t   entity_type,
                                        PDM_g_num_t         **entity_ln_to_gn,
                                        PDM_ownership_t       ownership);
    # ------------------------------------------------------------------
    int PDM_multipart_partition_color_get(PDM_multipart_t      *mtp,
                                          int                   i_zone,
                                          int                   i_part,
                                          PDM_mesh_entities_t   entity_type,
                                          int                 **entity_ln_to_gn,
                                          PDM_ownership_t       ownership);

    # ------------------------------------------------------------------
    int PDM_multipart_part_tn_part_get(PDM_multipart_t *multipart,
                                       const int        i_zone);

    # ------------------------------------------------------------------
    void PDM_multipart_time_get(PDM_multipart_t      *mtp,
                                int                   i_zone,
                                double              **elapsed,
                                double              **cpu,
                                double              **cpu_user,
                                double              **cpu_sys)

    # ------------------------------------------------------------------
    int PDM_multipart_part_n_entity_get(PDM_multipart_t      *multipart,
                                        int                   i_zone,
                                        int                   i_part,
                                        PDM_mesh_entities_t   entity_type)
    # ------------------------------------------------------------------
    void PDM_multipart_free(PDM_multipart_t      *mtp)

# ------------------------------------------------------------------
cdef class MultiPart:
    # > For Ppart
    cdef PDM_multipart_t* _mtp
    cdef int n_rank
    # ------------------------------------------------------------------
    # Fake init (Use only for docstring)
    def __init__(self,
                  int                                           n_zone,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                  int                                           merge_blocks,
                  PDM_split_dual_t                              split_method,
                  PDM_part_size_t                               part_size_method,
                  NPY.ndarray[NPY.double_t  , mode='c', ndim=1] part_fraction,
                  MPI.Comm                                      comm):

      """
      __init__(n_zone, n_part, merge_blocks, split_method, part_size_method, part_fraction, comm)

      Build a MultiPart structure. This method allows to split multiple zones.

      Parameters:
        n_zone           (int)                     : Number of zones in the original mesh
        n_part           (np.ndarray[np.int32_t])  : Number of partition per rank in each zone
        merge_blocks     (int)                     : Merge or not the zones before splitting
        split_method     (PDM_split_dual_t)        : Choice of library used to split the mesh
        part_size_method (PDM_part_size_t)         : Choice of homogeneous or heterogeneous partitions
        part_fraction    (np.ndarray[np.double_t]) : Weight (in %) of each partition in heterogeneous case if ``part_size_method`` is set to ``PDM_PART_SIZE_HETEROGENEOUS``
        comm             (MPI.Comm)                : MPI communicator
      """

    def __cinit__(self,
                  int                                           n_zone,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_part,
                  int                                           merge_blocks,
                  PDM_split_dual_t                              split_method,
                  PDM_part_size_t                               part_size_method,
                  NPY.ndarray[NPY.double_t  , mode='c', ndim=1] part_fraction,
                  MPI.Comm                                      comm):

      """
      __cinit__(n_zone, n_part, merge_blocks, split_method, part_size_method, part_fraction, comm)

      Build a MultiPart structure. This method allows to split multiple zones.

      Parameters:
        n_zone           (int)                     : Number of zones in the original mesh
        n_part           (np.ndarray[np.int32_t])  : Number of partition per rank in each zone
        merge_blocks     (int)                     : Merge or not the zones before splitting
        split_method     (PDM_split_dual_t)        : Choice of library used to split the mesh
        part_size_method (PDM_part_size_t)         : Choice of homogeneous or heterogeneous partitions
        part_fraction    (np.ndarray[np.double_t]) : Weight (in %) of each partition in heterogeneous case if ``part_size_method`` is set to ``PDM_PART_SIZE_HETEROGENEOUS``
        comm             (MPI.Comm)                : MPI communicator
      """
      # ~> MPI communicator
      cdef MPI.MPI_Comm c_comm = comm.ob_mpi
      cdef double* part_fraction_data

      # print("MultiPart::n_zone -->", n_zone)
      # print("MultiPart::n_part -->", n_part)
      # print("MultiPart::merge_blocks -->", merge_blocks)
      # print("MultiPart::split_method -->", split_method)
      self.n_rank = comm.Get_size()

      if part_fraction is None:
        part_fraction_data = NULL
      else:
        part_fraction_data = <double *> part_fraction.data

      # -> Create PPART
      self._mtp = PDM_multipart_create(n_zone,
                                            <int*> n_part.data,
                                            <PDM_bool_t> merge_blocks,
                                            split_method,
                                            part_size_method,
                                            part_fraction_data,
                                            PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                            PDM_OWNERSHIP_USER) # Python take ownership

    # ------------------------------------------------------------------
    def __dealloc__(self):
        PDM_multipart_free(self._mtp)

    # ------------------------------------------------------------------
    def multipart_register_block(self, int i_zone,
                                       DMesh dmesh): # DMesh = DistributedMeshCaspule or DistributedMesh
      """
      multipart_register_block(i_zone, dmesh)

      Set distributed mesh data for the input zone

      Parameters:
        i_zone (int)               : Zone identifier
        dmesh  (:py:class:`DMesh`) : Distributed mesh structure
      """
      PDM_multipart_register_block(self._mtp,
                                   i_zone,
                                   dmesh._dm)

    # ------------------------------------------------------------------
    def multipart_register_dmesh_nodal(self, int i_zone,
                                       DMeshNodal dmn):
      """
      multipart_register_dmesh_nodal(i_zone, dmn)

      Set distributed mesh data for the input zone.
      The mesh is described by nodal connectivity.

      Parameters:
        i_zone (int)                    : Zone identifier
        dmn    (:py:class:`DMeshNodal`) : Distributed nodal mesh structure
      """
      PDM_multipart_register_dmesh_nodal(self._mtp,
                                         i_zone,
                                         dmn.dmn)

    # ------------------------------------------------------------------
    def multipart_register_joins(self, int n_total_joins,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] matching_join):
      """
      multipart_register_joins(n_total_joins, matching_join)

      Set connecting data between all the zones

      Parameters:
        n_total_joins (int)                    : Total number of interfaces
        matching_join (np.ndarray[np.int32_t]) :  For each global join id, give the global id of the opposite
      """
      PDM_multipart_register_joins(       self._mtp,
                                          n_total_joins,
                                   <int*> matching_join.data)
    # ------------------------------------------------------------------
    def multipart_set_reordering(self, int i_zone,
                                 char *renum_cell_method,
                                 char *renum_face_method,
                                 NPY.ndarray[NPY.int32_t, mode='c', ndim=1] renum_properties_cell):
      """
      multipart_set_reordering(self, i_zone, renum_cell_method, renum_face_method, renum_properties_cell)

      Set the reordering methods to be used after partitioning

      Parameters:
        i_zone                (int)                    : Id of zone which parameters apply (or -1 for all zones)
        renum_cell_method     (str)                    : Choice of renumbering method for cells
        renum_cell_properties (np.ndarray[np.int32_t]) : Parameters used by cache-blocking method : ``(n_cell_per_cache_wanted, is_asynchronous, is_vectorisation, n_vect_face, split_method)``
        renum_face_method     (str)                    : Choice of renumbering method for faces
      """
      cdef int *renum_properties_cell_data
      if (renum_properties_cell is None):
        renum_properties_cell_data = NULL
      else:
        renum_properties_cell_data = <int *> renum_properties_cell.data
      PDM_multipart_set_reordering_options(self._mtp,
                                           i_zone,
                                           renum_cell_method,
                                           renum_properties_cell_data,
                                           renum_face_method)
    # ------------------------------------------------------------------
    def multipart_set_reordering_vtx(self, int i_zone,
                                     char *renum_vtx_method):
      """
      multipart_set_reordering_vtx(i_zone, renum_vtx_method)

      Set the vertex reordering methods to be used after partitioning

      Parameters:
        i_zone           (int) : Id of zone which parameters apply (or -1 for all zones)
        renum_vtx_method (str) : Choice of renumbering method for vertices
      """
      PDM_multipart_set_reordering_options_vtx(self._mtp,
                                               i_zone,
                                               renum_vtx_method)
    # ------------------------------------------------------------------
    def multipart_run_ppart(self):
      """
      Construct the partitioned meshes on all zones
      """
      PDM_multipart_run_ppart(self._mtp)

    # ------------------------------------------------------------------
    def multipart_n_entity_get(self, int i_part, int i_zone, PDM_mesh_entities_t entity_type):
      """
      multipart_n_entity_get(i_part, i_zone, entity_type)

      Get the number of entities with given type

      Parameters:
        i_part      (int)                 : Partition identifier
        i_zone      (int)                 : Zone identifier
        entity_type (PDM_mesh_entities_t) : Entity type
      """
      return PDM_multipart_part_n_entity_get(self._mtp, i_zone, i_part, entity_type)

    def multipart_dim_get(self, int ipart, int i_zone):
        """
        Get partition dimensions
        """
        # ************************************************************************
        # > Declaration
        cdef int n_cell
        cdef int n_proc
        cdef int nt_part
        cdef int n_face
        cdef int n_face_part_bound
        cdef int n_vtx
        cdef int scell_face
        cdef int s_face_vtx
        cdef int s_face_bound
        cdef int n_face_bound
        # ************************************************************************

        PDM_multipart_part_dim_get(self._mtp,
                                   i_zone,
                                   ipart,
                                   &n_cell,
                                   &n_face,
                                   &n_face_part_bound,
                                   &n_vtx,
                                   &n_proc,
                                   &nt_part,
                                   &scell_face,
                                   &s_face_vtx,
                                   &s_face_bound,
                                   &n_face_bound)

        return {'n_cell'            : n_cell,
                'ipart'             : ipart,
                'n_face'            : n_face,
                'nt_part'           : nt_part,
                'n_proc'            : n_proc,
                'n_face_part_bound' : n_face_part_bound,
                'n_vtx'             : n_vtx,
                'scell_face'        : scell_face,
                's_face_vtx'        : s_face_vtx,
                's_face_bound'      : s_face_bound,
                'n_face_bound'      : n_face_bound}

    # ------------------------------------------------------------------
    def multipart_val_get(self, int ipart, int i_zone):
        """
           Get partition dimensions (deprecated)
        """
        # ************************************************************************
        # > Declaration
        cdef int          *cell_face_idx,
        cdef int          *cell_face,
        cdef PDM_g_num_t  *cell_ln_to_gn,
        cdef int          *face_cell,
        cdef int          *face_vtx_idx,
        cdef int          *face_vtx,
        cdef PDM_g_num_t  *face_ln_to_gn,
        cdef int          *face_part_bound,
        cdef int          *face_part_bound_proc_idx,
        cdef int          *face_part_bound_part_idx,
        cdef double       *vtx_coord,
        cdef PDM_g_num_t  *vtx_ln_to_gn,
        cdef int          *face_bound_idx,
        cdef int          *face_bound,
        cdef PDM_g_num_t  *face_bound_ln_to_gn
        # ************************************************************************

        # dims = self.part_dim_get(self._mtp, ipart)
        dims = self.multipart_dim_get(ipart, i_zone)

        # -> Call PPART to get info
        PDM_multipart_part_val_get(self._mtp,
                                   i_zone,
                                   ipart,
                                   &cell_face_idx,
                                   &cell_face,
                                   &cell_ln_to_gn,
                                   &face_cell,
                                   &face_vtx_idx,
                                   &face_vtx,
                                   &face_ln_to_gn,
                                   &face_part_bound_proc_idx,
                                   &face_part_bound_part_idx,
                                   &face_part_bound,
                                   &vtx_coord,
                                   &vtx_ln_to_gn,
                                   &face_bound_idx,
                                   &face_bound,
                                   &face_bound_ln_to_gn)

        return {
            'np_cell_face_idx'            : create_numpy_or_none_i(cell_face_idx,            dims['n_cell']+1),
            'np_cell_face'                : create_numpy_or_none_i(cell_face,                dims['scell_face']),
            'np_cell_ln_to_gn'            : create_numpy_or_none_g(cell_ln_to_gn,            dims['n_cell']),
            'np_face_cell'                : create_numpy_or_none_i(face_cell,                2*dims['n_face']),
            'np_face_vtx_idx'             : create_numpy_or_none_i(face_vtx_idx,             dims['n_face']+1),
            'np_face_vtx'                 : create_numpy_or_none_i(face_vtx,                 dims['s_face_vtx']),
            'np_face_ln_to_gn'            : create_numpy_or_none_g(face_ln_to_gn,            dims['n_face']),
            'np_face_part_bound_proc_idx' : create_numpy_or_none_i(face_part_bound_proc_idx, dims['n_proc']+1),
            'np_face_part_bound_part_idx' : create_numpy_or_none_i(face_part_bound_part_idx, dims['nt_part']+1),
            'np_face_part_bound'          : create_numpy_or_none_i(face_part_bound,          4*dims['n_face_part_bound']),
            'np_vtx_coord'                : create_numpy_or_none_d(vtx_coord,                3*dims['n_vtx']),
            'np_vtx_ln_to_gn'             : create_numpy_or_none_g(vtx_ln_to_gn,             dims['n_vtx']),
            'np_face_bound_idx'           : create_numpy_or_none_i(face_bound_idx,           dims['n_face_bound']+1),
            'np_face_bound'               : create_numpy_or_none_i(face_bound,               dims['s_face_bound']),
            'np_face_bound_ln_to_gn'      : create_numpy_or_none_g(face_bound_ln_to_gn,      dims['s_face_bound'])
        }


    # ------------------------------------------------------------------
    def multipart_part_mesh_nodal_get(self, int i_zone):
      """
      multipart_part_mesh_nodal_get(i_zone)

      Retrieve the partitioned nodal mesh

      Parameters:
        i_zone (int) : Zone identifier

      Returns:
        Partitioned nodal mesh object (:py:class:`PartMeshNodalCaspule`)
      """
      cdef PDM_part_mesh_nodal_t *pmesh_nodal
      PDM_multipart_get_part_mesh_nodal(self._mtp, i_zone, &pmesh_nodal, PDM_OWNERSHIP_USER)
      if pmesh_nodal == NULL:
        return None
      else:
        #See pdm_part_mesh_nodal.pxi
        py_caps = PyCapsule_New(pmesh_nodal, NULL, NULL);
        return PartMeshNodalCaspule(py_caps)

    # ------------------------------------------------------------------
    def multipart_hyper_plane_color_get(self, int i_part, int i_zone):
      """
      multipart_hyper_plane_color_get(i_part, i_zone)

      Get array containing hyperplane color

      Parameters:
        i_part (int) : Partition identifier
        i_zone (int) : Zone identifier

      Returns:
        Hyperplane color (np.ndarray[np.int32_t])
      """
      # ************************************************************************
      # > Declaration
      cdef int           n_cell
      cdef int          *hyper_plane_color
      # ************************************************************************

      n_cell = PDM_multipart_part_n_entity_get(self._mtp, i_zone, i_part, PDM_MESH_ENTITY_CELL)
      PDM_multipart_part_hyperplane_color_get(self._mtp, i_zone, i_part, &hyper_plane_color, PDM_OWNERSHIP_USER);

      return {'np_hyper_plane_color' : create_numpy_or_none_i(hyper_plane_color, n_cell)}

    # ------------------------------------------------------------------
    def multipart_thread_color_get(self, int i_part, int i_zone):
      """
      multipart_thread_color_get(i_part, i_zone)

      Get partition dimensionsGet array containing thread color - Only if specific reordering (in ParaDiGMA plugins)

      Parameters:
        i_part (int) : Partition identifier
        i_zone (int) : Zone identifier

      Returns:
        Thread color (np.ndarray[np.int32_t])
      """
      # ************************************************************************
      # > Declaration
      cdef int           n_cell
      cdef int          *thread_color
      # ************************************************************************

      n_cell = PDM_multipart_part_n_entity_get(self._mtp, i_zone, i_part, PDM_MESH_ENTITY_CELL)
      PDM_multipart_part_thread_color_get(self._mtp, i_zone, i_part, &thread_color, PDM_OWNERSHIP_USER);

      return {'np_thread_color' : create_numpy_or_none_i(thread_color, n_cell)}


    # ------------------------------------------------------------------
    def multipart_ghost_information_get(self, int i_part, int i_zone):
      """
      multipart_ghost_information_get(self, int i_part, int i_zone)

      Get array containing vtx_ghost_information, useful to have a priority on vertex between multiple partitions

      Parameters:
        i_part (int) : Partition identifier
        i_zone (int) : Zone identifier

      Returns:
        Integer that gives the current priority of vertices on current partitions (np.ndarray[np.int32_t])
      """
      # ************************************************************************
      # > Declaration
      cdef int          *vtx_ghost_information
      # ************************************************************************

      # dims = self.part_dim_get(self._mtp, i_part)
      dims = self.multipart_dim_get(i_part, i_zone)

      # -> Call PPART to get info
      PDM_multipart_part_ghost_infomation_get(self._mtp,
                                              i_zone,
                                              i_part,
                                              &vtx_ghost_information,
                                              PDM_OWNERSHIP_USER)
      # -> Begin
      cdef NPY.npy_intp dim

      return {'np_vtx_ghost_information' : create_numpy_or_none_i(vtx_ghost_information, dims['n_vtx'])}

    # ------------------------------------------------------------------
    def multipart_connectivity_get(self, int i_part, int i_zone, PDM_connectivity_type_t connectivity_type):
      """
      multipart_connectivity_get(i_part, i_zone, connectivity_type)

      Get a partitioned connectivity

      Parameters:
        i_part            (int)                     : Partition identifier
        i_zone            (int)                     : Zone identifier
        connectivity_type (PDM_connectivity_type_t) : Connectivity type

      Returns:
        Dictionary:
          - ``"np_entity1_entity2"``     (np.ndarray[np.int32_t]) : Connectivity array
          - ``"np_entity1_entity2_idx"`` (np.ndarray[np.int32_t]) : Connectivity index
      """
      # ************************************************************************
      # > Declaration
      cdef int          *entity1_entity2
      cdef int          *entity1_entity2_idx
      # ************************************************************************

      # -> Call PPART to get info
      n_entity1 = PDM_multipart_part_connectivity_get(self._mtp,
                                                      i_zone,
                                                      i_part,
                                                      connectivity_type,
                                                      &entity1_entity2,
                                                      &entity1_entity2_idx,
                                                      PDM_OWNERSHIP_USER)

      if (entity1_entity2_idx == NULL and entity1_entity2 != NULL):
          np_entity1_entity2_idx = None
          np_entity1_entity2 = create_numpy_i(entity1_entity2, 2*n_entity1)
      elif(entity1_entity2 != NULL):
          np_entity1_entity2_idx = create_numpy_i(entity1_entity2_idx, n_entity1+1)
          np_entity1_entity2 = create_numpy_i(entity1_entity2, np_entity1_entity2_idx[n_entity1])
      else:
        np_entity1_entity2_idx = None
        np_entity1_entity2     = None


      return {'np_entity1_entity2'     : np_entity1_entity2,
              'np_entity1_entity2_idx' : np_entity1_entity2_idx}


    # ------------------------------------------------------------------
    def multipart_ln_to_gn_get(self, int i_part, int i_zone, PDM_mesh_entities_t entity_type):
      """
      multipart_ln_to_gn_get(i_part, i_zone, entity_type)

      Get the global ids of entities with given type

      Parameters:
        i_part      (int)               : Partition identifier
        i_zone      (int)               : Zone identifier
        entity_type (PDM_entity_type_t) : Entity type

      Returns:
        Global ids (np.ndarray[npy_pdm_gnum_t])
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_g_num_t  *entity_ln_to_gn
      # ************************************************************************

      # -> Call PPART to get info
      n_entity1 = PDM_multipart_part_ln_to_gn_get(self._mtp,
                                                  i_zone,
                                                  i_part,
                                                  entity_type,
                                                  &entity_ln_to_gn,
                                                  PDM_OWNERSHIP_USER)

      return {'np_entity_ln_to_gn' : create_numpy_or_none_g(entity_ln_to_gn, n_entity1)}

    # ------------------------------------------------------------------
    def multipart_vtx_coord_get(self, int i_part, int i_zone):
      """
      multipart_vtx_coord_get(i_part, i_zone)

      Get vertex coordinates

      Parameters:
          i_part (int) : Partition identifier
          i_zone (int) : Zone identifier

      Returns:
        Vertex coordinates (np.ndarray[np.double_t])
      """
      # ************************************************************************
      # > Declaration
      cdef double *vtx_coord
      # ************************************************************************

      n_vtx = PDM_multipart_part_vtx_coord_get(self._mtp,
                                               i_zone,
                                               i_part,
                                               &vtx_coord,
                                               PDM_OWNERSHIP_USER)

      return {'np_vtx_coord' : create_numpy_or_none_d(vtx_coord, 3*n_vtx)}

    # ------------------------------------------------------------------
    def multipart_part_color_get(self, int i_part, int i_zone, PDM_mesh_entities_t entity_type):
      """
      multipart_part_color_get(i_part, i_zone, entity_type)

      Get the color of entities with given type

      Parameters:
        i_part      (int)               : Partition identifier
        i_zone      (int)               : Zone identifier
        entity_type (PDM_entity_type_t) : Entity type

      Returns:
        Entity color (np.ndarray[np.int32_t])
      """
      # ************************************************************************
      # > Declaration
      cdef int  *entity_color
      # ************************************************************************

      # -> Call PPART to get info
      n_entity1 = PDM_multipart_partition_color_get(self._mtp,
                                                    i_zone,
                                                    i_part,
                                                    entity_type,
                                                    &entity_color,
                                                    PDM_OWNERSHIP_USER)

      return {'np_entity_color' : create_numpy_or_none_i(entity_color, n_entity1)}

    # ------------------------------------------------------------------
    def multipart_graph_comm_get(self,
                                 int i_part,
                                 int i_zone,
                                 PDM_bound_type_t bound_type):
        """
        multipart_graph_comm_get(i_part, i_zone, bound_type)

        Returns the connection graph between partition for the requested bound type

        Parameters:
          i_part     (int)              : Partition identifier
          i_zone     (int)              : Zone identifier
          bound_type (PDM_bound_type_t) : Bound type

        Returns:
          Dictionary:
            - ``"np_entity_part_bound_proc_idx"`` (np.ndarray[np.int32_t]) : Partitioning boundary entities index from process
            - ``"np_entity_part_bound_part_idx"`` (np.ndarray[np.int32_t]) : Partitioning boundary entities index from partition
            - ``"np_entity_part_bound"``          (np.ndarray[np.int32_t]) : Partitioning boundary entities
        """
        # ************************************************************************
        # > Declaration
        cdef int          *entity_part_bound
        cdef int          *entity_part_bound_proc_idx
        cdef int          *entity_part_bound_part_idx
        # ************************************************************************
        # -> Call PPART to get info
        PDM_multipart_part_graph_comm_get(self._mtp,
                                          i_zone,
                                          i_part,
                                          bound_type,
                                          &entity_part_bound_proc_idx,
                                          &entity_part_bound_part_idx,
                                          &entity_part_bound,
                                          PDM_OWNERSHIP_USER)

        tn_part = PDM_multipart_part_tn_part_get(self._mtp, i_zone)

        if (entity_part_bound == NULL) :
            np_entity_part_bound = None
        else :
            np_entity_part_bound   = create_numpy_i(entity_part_bound, 4*entity_part_bound_part_idx[tn_part])


        return {'np_entity_part_bound_proc_idx'  : create_numpy_or_none_i(entity_part_bound_proc_idx, self.n_rank+1),
                'np_entity_part_bound_part_idx'  : create_numpy_or_none_i(entity_part_bound_part_idx, tn_part+1),
                'np_entity_part_bound'           : np_entity_part_bound}

    # ------------------------------------------------------------------
    def multipart_time_get(self, int i_zone):
        """
        Get times (not implemented)
        """
        # ************************************************************************
        # > Declaration
        cdef double *elapsed
        cdef double *cpu
        cdef double *cpu_user
        cdef double *cpu_sys
        # ************************************************************************

        PDM_multipart_time_get(self._mtp, i_zone, &elapsed, &cpu, &cpu_user, &cpu_sys)

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

