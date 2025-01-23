
cdef extern from "pdm_part_mesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_part_mesh_nodal_t:
      pass
    ctypedef struct PDM_part_mesh_nodal_elmts_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    PDM_part_mesh_nodal_t* PDM_part_mesh_nodal_create(int          mesh_dimension,
                                                      int          n_part,
                                                      PDM_MPI_Comm comm)

    void PDM_part_mesh_nodal_coord_set(PDM_part_mesh_nodal_t *pmn,
                                       int                    id_part,
                                       int                    n_vtx,
                                       double                *coords,
                                       PDM_ownership_t        owner)

    void PDM_part_mesh_nodal_vtx_gnum_set(PDM_part_mesh_nodal_t *pmn,
                                          int                    id_part,
                                          PDM_g_num_t           *numabs,
                                          PDM_ownership_t        owner)

    int PDM_part_mesh_nodal_n_part_get(PDM_part_mesh_nodal_t *pmn)
    int PDM_part_mesh_nodal_mesh_dimension_get( PDM_part_mesh_nodal_t *pmn)

    int PDM_part_mesh_nodal_n_vtx_get(PDM_part_mesh_nodal_t *pmn,
                                      int                    id_part)

    double* PDM_part_mesh_nodal_vtx_coord_get(PDM_part_mesh_nodal_t *pmn,
                                              int                    id_part,
                                              PDM_ownership_t        owner);


    PDM_g_num_t* PDM_part_mesh_nodal_vtx_g_num_get(PDM_part_mesh_nodal_t *pmn,
                                                   int                    id_part,
                                                   PDM_ownership_t        owner);

    int PDM_part_mesh_nodal_n_section_in_geom_kind_get(PDM_part_mesh_nodal_t *pmn,
                                                       PDM_geometry_kind_t    geom_kind)

    int * PDM_part_mesh_nodal_sections_id_in_geom_kind_get(PDM_part_mesh_nodal_t *pmn,
                                                           PDM_geometry_kind_t    geom_kind)

    int PDM_part_mesh_nodal_n_section_get(PDM_part_mesh_nodal_t *pmn)

    int * PDM_part_mesh_nodal_sections_id_get(PDM_part_mesh_nodal_t *pmn)


    PDM_Mesh_nodal_elt_t PDM_part_mesh_nodal_section_elt_type_get(PDM_part_mesh_nodal_t *pmn,
                                                                  int                    id_section)

    int PDM_part_mesh_nodal_section_add(PDM_part_mesh_nodal_t *pmn,
                                        PDM_Mesh_nodal_elt_t   t_elt)

    void PDM_part_mesh_nodal_section_std_set(PDM_part_mesh_nodal_t *pmn,
                                             int                    id_block,
                                             int                    id_part,
                                             int                    n_elt,
                                             int                   *connec,
                                             PDM_g_num_t           *numabs,
                                             int                   *parent_num,
                                             PDM_g_num_t           *parent_entity_g_num,
                                             PDM_ownership_t        owner)

    int PDM_part_mesh_nodal_section_n_elt_get(PDM_part_mesh_nodal_t  *pmn,
                                              int                     id_block,
                                              int                     id_part)


    void PDM_part_mesh_nodal_section_std_get(PDM_part_mesh_nodal_t  *pmn,
                                             int                     id_block,
                                             int                     id_part,
                                             int                   **connec,
                                             PDM_g_num_t           **numabs,
                                             int                   **parent_num,
                                             PDM_g_num_t           **parent_entity_g_num,
                                             PDM_ownership_t         ownership)

    int* PDM_part_mesh_nodal_section_elmt_to_entity_get(PDM_part_mesh_nodal_t *pmn,
                                                        const int              i_section,
                                                        const int              id_part,
                                                        PDM_ownership_t        ownership)

    int PDM_part_mesh_nodal_section_id_from_geom_kind_get(PDM_part_mesh_nodal_t  *pmn,
                                                          const PDM_geometry_kind_t     geom_kind,
                                                          const int                     id_section_in_geom_kind)

    void PDM_part_mesh_nodal_n_group_set(PDM_part_mesh_nodal_t  *pmn,
                                         PDM_geometry_kind_t     geom_kind,
                                         const int               n_group)

    void PDM_part_mesh_nodal_group_set(PDM_part_mesh_nodal_t  *pmn,
                                       PDM_geometry_kind_t     geom_kind,
                                       const int               i_part,
                                       const int               i_group,
                                       int                     n_group_elmt,
                                       int                    *group_elmt,
                                       PDM_g_num_t            *group_ln_to_gn,
                                       PDM_ownership_t         ownership)

    void PDM_part_mesh_nodal_group_get(PDM_part_mesh_nodal_t  *pmn,
                                       PDM_geometry_kind_t     geom_kind,
                                       const int               i_part,
                                       const int               i_group,
                                       int                    *n_group_elmt,
                                       int                   **group_elmt,
                                       PDM_g_num_t           **group_ln_to_gn,
                                       PDM_ownership_t         ownership);

    int PDM_part_mesh_nodal_n_group_get(PDM_part_mesh_nodal_t  *pmn,
                                        PDM_geometry_kind_t     geom_kind);

    void PDM_part_mesh_nodal_free( PDM_part_mesh_nodal_t* pmn)

# ------------------------------------------------------------------
cdef class PartMeshNodal:
    """
      PartMeshNodal: Mesh structure for multi-elmts description
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_part_mesh_nodal_t *pmn
    cdef object keep_alive
    cdef int n_rank
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm    comm,
                        int         n_part,
                        int         mesh_dimension = 3):
        """
        """
        self.keep_alive = []
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.n_rank = comm.Get_size()
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # > Convert mpi4py -> PDM_MPI
        cdef MPI.MPI_Comm c_comm = comm.ob_mpi
        cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.pmn = PDM_part_mesh_nodal_create(mesh_dimension, n_part, PDMC)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    def set_coordinates(self,
                        id_part,
                        NPY.ndarray[NPY.double_t  , mode='c', ndim=1] pvtx_coord,
                        NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] pvtx_ln_to_gn):
        """
        set_coordinates(id_part, pvtx_coord, pvtx_ln_to_gn)

        Define partition vertices

        Parameters:
          id_part       (int)                          : Partition identifier
          pvtx_coord    (`np.ndarray[np.double_t]`   ) : Vertex coordinates
          pvtx_ln_to_gn (`np.ndarray[npy_pdm_gnum_t]`) : Vertex global IDs
        """
        # ************************************************************************
        # > Declaration
        cdef int n_vtx
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.keep_alive.append(pvtx_coord)
        self.keep_alive.append(pvtx_ln_to_gn)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        n_vtx = pvtx_coord.shape[0]//3
        PDM_part_mesh_nodal_coord_set(self.pmn,
                                      id_part,
                                      n_vtx,
                                      np_to_double_pointer(pvtx_coord),
                                      PDM_OWNERSHIP_USER)
        PDM_part_mesh_nodal_vtx_gnum_set(self.pmn,
                                         id_part,
                                         np_to_gnum_pointer(pvtx_ln_to_gn),
                                         PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def add_section(self,
                    PDM_Mesh_nodal_elt_t elmt_type):
      """
      add_section(elmt_type)

        Add for each part, a section of elmt_type and return the internal associated id in order to call set_section

        Parameters:
          elmt_type (PDM_Mesh_nodal_elt_t) : Kind of element for current sections
      """
      id_section = PDM_part_mesh_nodal_section_add(self.pmn,
                                                   elmt_type)
      return id_section

    # ------------------------------------------------------------------------
    def set_section(self,
                    int                                           id_section,
                    int                                           id_part,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] elmt_vtx,
                    NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] numabs,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] parent_num,
                    NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] parent_entity_g_num,
                    int                                           n_elemts):
        """
        set_section(id_section, id_part, elmt_vtx, numabs, parent_num, parent_entity_g_num, n_elemts)

        For id_section and id_part, set the element connectivity and the associated parent_num

        Parameters:
          id_section (int)                          : id of the section (return by add_section)
          id_part    (int)                          : id of the part (max = n_part)
          elmt_vtx   (`np.ndarray[np.int32_t]`)     : Element connectivity
          numabs     (`np.ndarray[npy_pdm_gnum_t]`) : Global numbering
          parent_num (`np.ndarray[np.int32_t]`)     : Correspondence table with a PartMesh (if any, else None)
          numabs     (`np.ndarray[npy_pdm_gnum_t]`) : Global numbering
          n_elemts   (int)                          : Number of elements in section
        """
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        self.keep_alive.append(elmt_vtx)
        self.keep_alive.append(numabs)
        self.keep_alive.append(parent_num)
        self.keep_alive.append(parent_entity_g_num)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        PDM_part_mesh_nodal_section_std_set(self.pmn,
                                            id_section,
                                            id_part,
                                            n_elemts,
                                            np_to_int_pointer(elmt_vtx),
                                            np_to_gnum_pointer(numabs),
                                            np_to_int_pointer(parent_num),
                                            np_to_gnum_pointer(parent_entity_g_num),
                                            PDM_OWNERSHIP_USER)
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def n_group_set(self,
                    PDM_geometry_kind_t geom_kind,
                    int                 n_group):
        """
        n_group_set(geom_kind, n_group)

        Set number of group for the geom_kind

        Parameters:
          geom_kind (PDM_geometry_kind_t) : Geometry kind (0D, 1D, 2D, 3D)
          n_group   (int)                 : Number of group on current geom_kind
        """
        PDM_part_mesh_nodal_n_group_set(self.pmn,
                                        geom_kind,
                                        n_group)

    # ------------------------------------------------------------------------
    def group_set(self,
                  PDM_geometry_kind_t                           geom_kind,
                  int                                           i_part,
                  int                                           i_group,
                  NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_elmt,
                  NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] group_ln_to_gn):
        """
        group_set(geom_kind, i_part, i_group, group_elmt, group_ln_to_gn)

        Set group for a partition

        Parameters:
          geom_kind      (PDM_geometry_kind_t)          : Geometry kind (VOLUMIC, SURFACIC, RIDGE, CORNER)
          i_part         (int)                          : id of the section (return by add_section)
          i_group        (int)                          : id of the group
          group_elmt     (`np.ndarray[np.int32_t]`)     : Group->entity connectivity (1-based local ids)
          group_ln_to_gn (`np.ndarray[npy_pdm_gnum_t]`) : Group->entity connectivity (group-specific global ids)
        """
        cdef int n_group_elmt
        n_group_elmt = group_elmt.shape[0]
        self.keep_alive.append(group_elmt)
        self.keep_alive.append(group_ln_to_gn)
        PDM_part_mesh_nodal_group_set(self.pmn,
                                      geom_kind,
                                      i_part,
                                      i_group,
                                      n_group_elmt,
                                      np_to_int_pointer(group_elmt),
                                      np_to_gnum_pointer(group_ln_to_gn),
                                      PDM_OWNERSHIP_USER)

    def to_view_capsule(self):
      """
      """
      return PyCapsule_New(self.pmn, NULL, NULL);

    def dim_get(self):
      """
      Get PartMeshNodal dimension

      Returns: Mesh dimension (0, 1, 2 or 3)
      """
      return part_mesh_nodal_dim_get(self)

    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
      """
      PDM_part_mesh_nodal_free(self.pmn)


# ------------------------------------------------------------------
cdef class PartMeshNodalCapsule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_part_mesh_nodal_t* pmn
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshNodalCapsule", PyCapsule_GetName(caps))
    cdef PDM_part_mesh_nodal_t* caps_pmn = <PDM_part_mesh_nodal_t *> PyCapsule_GetPointer(caps, NULL)
    self.pmn = caps_pmn;

  def dim_get(self):
    return part_mesh_nodal_dim_get(self)

  def n_part_get(self):
    return part_mesh_nodal_n_part_get(self)

  def coord_get(self, i_part):
    """
    coord_get(i_part)

    Get coordinates of mesh
    """
    return part_mesh_nodal_vtx_coord_get(self, i_part)

  def vtx_g_num_get(self, i_part):
    """
    vtx_g_num_get(i_part)

    Get global ids of mesh vertices
    """
    return part_mesh_nodal_vtx_g_num_get(self, i_part)

  def get_sections(self, PDM_geometry_kind_t geom_kind, int i_part):
    """
    get_sections(geom_kind, i_part)

    Get all standard sections (one partition)

    Parameters:
      geom_kind (PDM_geometry_kind_t) : Geometry kind (volume, surface, ridge or corner)
      i_part    (int)                 : Partition identifier

    Returns:
      List of sections. Each section is represented as a dictionary

        - ``"pdm_type"``               (`int`)                        : Element type
        - ``"np_connec"``              (`np.ndarray[np.int32_t]`)     : Connectivity
        - ``"np_numabs"``              (`np.ndarray[npy_pdm_gnum_t]`) : Element global ids
        - ``"np_parent_num"``          (`np.ndarray[np.int32_t]`)     : Element parent local ids
        - ``"np_parent_entity_g_num"`` (`np.ndarray[npy_pdm_gnum_t]`) : Element parent global ids
    """
    return part_mesh_nodal_get_sections(self, geom_kind, i_part)

  def get_n_group(self, PDM_geometry_kind_t geom_kind):
    """
    get_n_group(geom_kind)

    Get group number for given ``geom_kind``

    Parameters:
      geom_kind (PDM_geometry_kind_t) : Geometry kind (volume, surface, ridge or corner)

    Returns:
      Group number
    """
    return part_mesh_nodal_n_group_get(self, geom_kind)

  def get_group(self, PDM_geometry_kind_t geom_kind, int i_part, int i_group):
    """
    get_group(geom_kind, i_part, i_group)

    Get partition group

    Parameters:
      geom_kind (PDM_geometry_kind_t) : Geometry kind (volume, surface, ridge or corner)
      i_part    (int)                 : Partition identifier
      i_group   (int)                 : Group identifier

    Returns:
      Tuple

        - ``"group_elmt"``             (`np.ndarray[np.int32_t]`)     : Connectivity group elements
        - ``"group_ln_to_gn"``         (`np.ndarray[npy_pdm_gnum_t]`) : Element global ids
    """
    return part_mesh_nodal_get_group(self, geom_kind, i_part, i_group)

  def __dealloc__(self):
    """
    """
    PDM_part_mesh_nodal_free(self.pmn)


ctypedef fused PMeshNodal:
  PartMeshNodal
  PartMeshNodalCapsule

def part_mesh_nodal_vtx_g_num_get(PMeshNodal pypmn, int i_part):
  """
  Get vertex global numbering for partition i_part
  """
  cdef PDM_g_num_t *vtx_ln_to_gn = NULL
  vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pypmn.pmn, i_part, PDM_OWNERSHIP_USER)

  n_vtx = PDM_part_mesh_nodal_n_vtx_get(pypmn.pmn, i_part)

  return create_numpy_g(vtx_ln_to_gn, n_vtx, True)


def part_mesh_nodal_vtx_coord_get(PMeshNodal pypmn, int i_part):
  """
  Get vertex coordinates for partition i_part
  """

  cdef double *vtx_coords = NULL
  vtx_coords = PDM_part_mesh_nodal_vtx_coord_get(pypmn.pmn, i_part, PDM_OWNERSHIP_USER)

  n_vtx = PDM_part_mesh_nodal_n_vtx_get(pypmn.pmn, i_part)

  return create_numpy_d(vtx_coords, 3 * n_vtx, True)


def part_mesh_nodal_dim_get(PMeshNodal pypmn):
  return PDM_part_mesh_nodal_mesh_dimension_get(pypmn.pmn)


def part_mesh_nodal_n_part_get(PMeshNodal pypmn):
  return PDM_part_mesh_nodal_n_part_get(pypmn.pmn)


def part_mesh_nodal_get_sections(PMeshNodal pypmn, PDM_geometry_kind_t geom_kind, int i_part):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_section
  cdef int                   n_vtx_per_elmt
  cdef int                   n_elmt_in_section
  cdef int                  *section_id
  cdef int                  *parent_num
  cdef int                  *elt2entity
  cdef int                  *connec
  cdef PDM_g_num_t          *numabs
  cdef PDM_g_num_t          *parent_entity_g_num
  cdef double               *vtx_coord
  cdef PDM_Mesh_nodal_elt_t  t_elmt
  cdef NPY.npy_intp          dim
  # ************************************************************************

  n_section  = PDM_part_mesh_nodal_n_section_in_geom_kind_get  (pypmn.pmn, geom_kind)
  section_id = PDM_part_mesh_nodal_sections_id_in_geom_kind_get(pypmn.pmn, geom_kind)

  sections = []
  for i_section in range(n_section):
    id_section_in_geom_kind = section_id[i_section]
    id_section = PDM_part_mesh_nodal_section_id_from_geom_kind_get(pypmn.pmn,
                                                                   geom_kind,
                                                                   id_section_in_geom_kind)

    t_elmt = PDM_part_mesh_nodal_section_elt_type_get(pypmn.pmn, id_section)
    assert(t_elmt != PDM_MESH_NODAL_POLY_2D)
    assert(t_elmt != PDM_MESH_NODAL_POLY_3D)

    n_elmt_in_section = PDM_part_mesh_nodal_section_n_elt_get(pypmn.pmn, id_section, i_part)

    PDM_part_mesh_nodal_section_std_get(pypmn.pmn, id_section, i_part, &connec, &numabs, &parent_num, &parent_entity_g_num, PDM_OWNERSHIP_USER)
    elt2entity = PDM_part_mesh_nodal_section_elmt_to_entity_get(pypmn.pmn, id_section, i_part, PDM_OWNERSHIP_USER)

    n_vtx_per_elmt = PDM_Mesh_nodal_n_vtx_elt_get(t_elmt, 1)

    np_connec     = create_numpy_i(connec,     n_elmt_in_section*n_vtx_per_elmt)
    np_parent_num = None
    if(parent_num != NULL):
      np_parent_num = create_numpy_i(parent_num, n_elmt_in_section)

    np_elt_entity = None
    if elt2entity != NULL:
      np_elt_entity = create_numpy_i(elt2entity, n_elmt_in_section)
    np_numabs = create_numpy_g(numabs, n_elmt_in_section)

    np_parent_entity_g_num = None
    if(parent_entity_g_num != NULL):
      np_parent_entity_g_num = create_numpy_g(parent_entity_g_num, n_elmt_in_section)

    sections.append({"n_elmt"                 : n_elmt_in_section,
                     "pdm_type"               : t_elmt,
                     "np_connec"              : np_connec,
                     "np_parent_num"          : np_parent_num,
                     "np_element_to_entity"   : np_elt_entity,
                     "np_numabs"              : np_numabs,
                     "np_parent_entity_g_num" : np_parent_entity_g_num})

  return sections

def part_mesh_nodal_n_group_get(PMeshNodal pypmn, PDM_geometry_kind_t geom_kind):
  """
  """
  return PDM_part_mesh_nodal_n_group_get(pypmn.pmn,
                                         geom_kind)

def part_mesh_nodal_get_group(PMeshNodal pypmn, PDM_geometry_kind_t geom_kind, int i_part, int i_group):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int                   n_group_elmt
  cdef int                  *group_elmt
  cdef PDM_g_num_t          *group_ln_to_gn
  # ************************************************************************


  PDM_part_mesh_nodal_group_get(pypmn.pmn,
                                geom_kind,
                                i_part,
                                i_group,
                                &n_group_elmt,
                                &group_elmt,
                                &group_ln_to_gn,
                                PDM_OWNERSHIP_USER);

  np_group_elmt = create_numpy_i(group_elmt, n_group_elmt)
  if(group_ln_to_gn != NULL):
    np_group_ln_to_gn = create_numpy_g(group_ln_to_gn, n_group_elmt)

  return np_group_elmt, np_group_ln_to_gn
