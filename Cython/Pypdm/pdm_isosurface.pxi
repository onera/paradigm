
cdef extern from "pdm_isosurface.h":
  ctypedef struct PDM_isosurface_t:
      pass

  ctypedef double (*PDM_isosurface_python_field_function_t)(void   *python_object,
                                                            int     id_iso,
                                                            double  x,
                                                            double  y,
                                                            double  z);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of function
  PDM_isosurface_t *PDM_isosurface_create(PDM_MPI_Comm comm,
                                          int          mesh_dimension);

  void PDM_isosurface_set_tolerance(PDM_isosurface_t *isos,
                                    double            tolerance);

  void PDM_isosurface_n_part_set(PDM_isosurface_t *isos,
                                 int               n_part);

  void PDM_isosurface_connectivity_set(PDM_isosurface_t        *isos,
                                       int                      i_part,
                                       PDM_connectivity_type_t  connectivity_type,
                                       int                      n_entity,
                                       int                     *connect_idx,
                                       int                     *connect);

  void PDM_isosurface_vtx_coord_set(PDM_isosurface_t *isos,
                                    int               i_part,
                                    int               n_vtx,
                                    double           *vtx_coord);

  void PDM_isosurface_ln_to_gn_set(PDM_isosurface_t    *isos,
                                   int                  i_part,
                                   PDM_mesh_entities_t  entity_type,
                                   PDM_g_num_t         *ln_to_gn);

  void PDM_isosurface_n_group_set(PDM_isosurface_t    *isos,
                                  PDM_mesh_entities_t  entity_type,
                                  int                  n_group);

  void PDM_isosurface_group_set(PDM_isosurface_t    *isos,
                                int                  i_part,
                                PDM_mesh_entities_t  entity_type,
                                int                 *group_entity_idx,
                                int                 *group_entity,
                                PDM_g_num_t         *group_entity_ln_to_gn);

  void PDM_isosurface_part_mesh_set(PDM_isosurface_t *isos,
                                    PDM_part_mesh_t  *pmesh);

  void PDM_isosurface_part_mesh_nodal_set(PDM_isosurface_t      *isos,
                                          PDM_part_mesh_nodal_t *pmn);

  void PDM_isosurface_dconnectivity_set(PDM_isosurface_t        *isos,
                                        PDM_connectivity_type_t  connectivity_type,
                                        int                     *dconnect_idx,
                                        PDM_g_num_t             *dconnect);

  void PDM_isosurface_dvtx_coord_set(PDM_isosurface_t *isos,
                                     double           *dvtx_coord);

  void PDM_isosurface_distrib_set(PDM_isosurface_t    *isos,
                                  PDM_mesh_entities_t  entity_type,
                                  PDM_g_num_t         *distrib);

  void PDM_isosurface_dgroup_set(PDM_isosurface_t    *isos,
                                 PDM_mesh_entities_t  entity_type,
                                 int                 *dgroup_entity_idx,
                                 PDM_g_num_t         *dgroup_entity);

  void PDM_isosurface_dmesh_set(PDM_isosurface_t *isos,
                                PDM_dmesh_t      *dmesh);

  void PDM_isosurface_dmesh_nodal_set(PDM_isosurface_t  *isos,
                                      PDM_dmesh_nodal_t *dmn);

  int PDM_isosurface_add(PDM_isosurface_t       *isos,
                         PDM_iso_surface_kind_t  kind,
                         int                     n_isovalues,
                         double                 *isovalues);

  void PDM_isosurface_isovalues_set(PDM_isosurface_t *isos,
                                    int               id_isosurface,
                                    int               n_isovalues,
                                    double           *isovalues);

  void PDM_isosurface_equation_set(PDM_isosurface_t *isos,
                                   int               id_isosurface,
                                   double           *coeff);

  void isosurface_python_field_function_set(PDM_isosurface_t                       *isos,
                                            int                                     id_isosurface,
                                            PDM_isosurface_python_field_function_t  func);

  void isosurface_python_object_set(PDM_isosurface_t *isos,
                                    void             *python_object);

  void PDM_isosurface_field_set(PDM_isosurface_t *isos,
                                int               id_isosurface,
                                int               i_part,
                                double           *field);

  void PDM_isosurface_dfield_set(PDM_isosurface_t *isos,
                                 int               id_isosurface,
                                 double           *dfield);


  void PDM_isosurface_redistribution_set(PDM_isosurface_t        *isos,
                                         PDM_extract_part_kind_t  extract_kind,
                                         PDM_split_dual_t         part_method);


  void PDM_isosurface_reset(PDM_isosurface_t *isos,
                            int               id_isosurface);

  void PDM_isosurface_n_part_out_set(PDM_isosurface_t *isos,
                                     int               n_part_out);


  void PDM_isosurface_compute(PDM_isosurface_t *isos,
                              int               id_isosurface);

  int PDM_isosurface_connectivity_get(PDM_isosurface_t         *isos,
                                      int                       id_isosurface,
                                      int                       i_part,
                                      PDM_connectivity_type_t   connectivity_type,
                                      int                     **connect_idx,
                                      int                     **connect,
                                      PDM_ownership_t           ownership);

  int PDM_isosurface_vtx_coord_get(PDM_isosurface_t  *isos,
                                   int                id_isosurface,
                                   int                i_part,
                                   double           **vtx_coord,
                                   PDM_ownership_t    ownership);

  int PDM_isosurface_ln_to_gn_get(PDM_isosurface_t     *isos,
                                  int                   id_isosurface,
                                  int                   i_part,
                                  PDM_mesh_entities_t   entity_type,
                                  PDM_g_num_t         **ln_to_gn,
                                  PDM_ownership_t       ownership);

  int PDM_isosurface_group_get(PDM_isosurface_t     *isos,
                               int                   id_isosurface,
                               int                   i_part,
                               PDM_mesh_entities_t   entity_type,
                               int                 **group_entity_idx,
                               int                 **group_entity,
                               PDM_g_num_t         **group_entity_ln_to_gn,
                               PDM_ownership_t       ownership);

  int PDM_isosurface_isovalue_entity_idx_get(PDM_isosurface_t     *isos,
                                             int                   id_isosurface,
                                             int                   i_part,
                                             PDM_mesh_entities_t   entity_type,
                                             int                 **isovalue_entity_idx,
                                             PDM_ownership_t       ownership);

  int PDM_isosurface_dconnectivity_get(PDM_isosurface_t         *isos,
                                       int                       id_isosurface,
                                       PDM_connectivity_type_t   connectivity_type,
                                       int                     **dconnect_idx,
                                       PDM_g_num_t             **dconnect,
                                       PDM_ownership_t           ownership);

  int PDM_isosurface_dparent_weight_get(PDM_isosurface_t     *isos,
                                        int                   id_iso,
                                        PDM_mesh_entities_t   entity_type,
                                        int                 **dvtx_parent_idx,
                                        double              **dvtx_parent_weight,
                                        PDM_ownership_t       ownership);

  int PDM_isosurface_dvtx_coord_get(PDM_isosurface_t  *isos,
                                    int                id_isosurface,
                                    double           **dvtx_coord,
                                    PDM_ownership_t    ownership);

  void PDM_isosurface_distrib_get(PDM_isosurface_t     *isos,
                                  int                   id_isosurface,
                                  PDM_mesh_entities_t   entity_type,
                                  PDM_g_num_t         **distribution);

  int PDM_isosurface_dgroup_get(PDM_isosurface_t     *isos,
                                int                   id_isosurface,
                                PDM_mesh_entities_t   entity_type,
                                int                 **dgroup_entity_idx,
                                PDM_g_num_t         **dgroup_entity,
                                PDM_ownership_t       ownership);

  int PDM_isosurface_disovalue_entity_get(PDM_isosurface_t     *isos,
                                          int                   id_isosurface,
                                          PDM_mesh_entities_t   entity_type,
                                          int                 **disovalue_entity_idx,
                                          PDM_g_num_t         **disovalue_entity,
                                          PDM_ownership_t       ownership);

  int PDM_isosurface_local_parent_get(PDM_isosurface_t     *isos,
                                      int                   id_isosurface,
                                      int                   i_part,
                                      PDM_mesh_entities_t   entity_type,
                                      int                 **entity_parent_idx,
                                      int                 **entity_parent,
                                      PDM_ownership_t       ownership);

  int PDM_isosurface_parent_weight_get(PDM_isosurface_t     *isos,
                                       int                   id_isosurface,
                                       int                   i_part,
                                       PDM_mesh_entities_t   entity_type,
                                       int                 **parent_idx,
                                       double              **parent_weight,
                                       PDM_ownership_t       ownership);

  void PDM_isosurface_part_to_part_enable(PDM_isosurface_t     *isos,
                                          int                   id_isosurface,
                                          PDM_mesh_entities_t   entity_type,
                                          int                   unify_parent_info);

  void PDM_isosurface_part_to_part_get(PDM_isosurface_t     *isos,
                                       int                   id_isosurface,
                                       PDM_mesh_entities_t   entity_type,
                                       PDM_part_to_part_t  **ptp,
                                       PDM_ownership_t       ownership);

  void PDM_isosurface_dump_times(PDM_isosurface_t *isos);

  void PDM_isosurface_free(PDM_isosurface_t  *isos);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#-------------------------------------------------------------------
# CALLBACK
cdef double callback(void   *_isos,
                     int     id_iso,
                     double  x,
                     double  y,
                     double  z):

  cdef Isosurface isos = <Isosurface> _isos

  assert(id_iso in isos.user_defined_field_function)

  cdef double value = isos.user_defined_field_function[id_iso](x, y, z)

  return value

# ------------------------------------------------------------------
cdef class Isosurface:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_isosurface_t *_isos

  # cdef list keep_alive
  cdef MPI.Comm py_comm
  cdef int      i_rank
  cdef int      n_rank
  cdef dict ptp_entity
  cdef dict user_defined_field_function
  cdef list got_parent_idx

  FIELD    = _PDM_ISO_SURFACE_KIND_FIELD
  PLANE    = _PDM_ISO_SURFACE_KIND_PLANE
  SPHERE   = _PDM_ISO_SURFACE_KIND_SPHERE
  ELLIPSE  = _PDM_ISO_SURFACE_KIND_ELLIPSE
  QUADRIC  = _PDM_ISO_SURFACE_KIND_QUADRIC
  HEART    = _PDM_ISO_SURFACE_KIND_HEART
  FUNCTION = _PDM_ISO_SURFACE_KIND_FUNCTION

  LOCAL         = _PDM_EXTRACT_PART_KIND_LOCAL
  REEQUILIBRATE = _PDM_EXTRACT_PART_KIND_REEQUILIBRATE
  FROM_TARGET   = _PDM_EXTRACT_PART_KIND_FROM_TARGET

  # ************************************************************************

  # > Common setter API
  def __init__(self, int mesh_dim, MPI.Comm comm):
    """
    Isosurface(mesh_dim, comm)

    Create the isosurface structure.

    Parameters:
      mesh_dim (int)      : Entry mesh dimension
      comm     (MPI.Comm) : MPI communicator
    """
    # self.keep_alive  = list()

    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
    cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    self.py_comm = comm
    self.i_rank  = comm.Get_rank()
    self.n_rank  = comm.Get_size()

    self._isos = PDM_isosurface_create(pdm_comm,
                                       mesh_dim)

    self.ptp_entity                  = dict()
    self.user_defined_field_function = dict()
    self.got_parent_idx              = [dict(), dict(), dict(), dict()]

    isosurface_python_object_set(self._isos,
                        <void *> self)

  def tolerance_set(self, tolerance):
    """
    tolerance_set(tolerance)

    Set isosurface tolerance. May improve resulting mesh quality.

    Parameters:
      tolerance (double) : Field tolerance for isosurface computing (Default at 0.)
    """
    PDM_isosurface_set_tolerance(self._isos, tolerance)

  def add(self, kind,
          list  isovalues):
    """
    add(kind, isovalues)

    Set isosurface kind and isovalues.

    Parameters:
      kind      (PDM_iso_surface_kind_t) : Isosurface kind
      isovalues (list of double)         : Isosurface values
    Returns:
      id_iso (int) Isosurface id
    """
    cdef int     id_iso         = -1
    cdef int     n_isovalues    = len(isovalues)
    cdef double *isovalues_data = list_to_double_pointer(isovalues)

    id_iso = PDM_isosurface_add(self._isos, kind, n_isovalues, isovalues_data)

    self.ptp_entity[id_iso] = dict()

    free(isovalues_data) # deep-copied in isos

    for entity_type in range(_PDM_MESH_ENTITY_MAX):
      self.got_parent_idx[entity_type][id_iso] = False

    return id_iso

  def isovalues_set(self, id_iso,
                    list  isovalues):
    """
    isovalues_set(id_iso, isovalues)

    Set (or reset) isovalues for specified isosurface id.

    Parameters:
      id_iso    (int)            : Isosurface id
      isovalues (list of double) : Isosurface values
    """
    cdef int     n_isovalues    = len(isovalues)
    cdef double *isovalues_data = list_to_double_pointer(isovalues)

    PDM_isosurface_isovalues_set(self._isos, id_iso, n_isovalues, isovalues_data)

    free(isovalues_data) # deep-copied in isos

  def equation_set(self, id_iso,
                   list coefficients):
    """
    equation_set(id_iso, coefficients)

    Set (or reset) equation coefficients for specified isosurface id.

    Parameters:
      id_iso       (int)            : Isosurface id
      coefficients (list of double) : Equation coefficients
    """
    cdef double *coeff_data = list_to_double_pointer(coefficients)

    PDM_isosurface_equation_set(self._isos, id_iso, coeff_data)

    free(coeff_data) # deep-copied in isos

  def field_function_set(self, id_iso, fct):
    """
    field_function_set(id_iso, fct)

    Set source field function

    Parameters:
      id_iso (int) : Isosurface id
      fct          : Field function
    """
    self.user_defined_field_function[id_iso] = fct

    isosurface_python_field_function_set(self._isos,
                                         id_iso,
                                         callback)

  def compute(self, id_iso):
    """
    compute(id_iso)

    Computed asked isosurface.

    Parameters:
      id_iso (int) : Isosurface id
    """
    PDM_isosurface_compute(self._isos, id_iso)

  def reset(self, id_iso):
    """
    reset(id_iso)

    Reset computed isosurface to allow future recomputation
    with updated isovalues or field.

    Parameters:
      id_iso (int) : Isosurface id
    """
    PDM_isosurface_reset(self._isos, id_iso)


  # > Partitionned setter API
  def n_part_set(self, n_part):
    """
    n_part(n_part)

    Set entry mesh number of partitions.

    Parameters:
      n_part (int) : Number of partitions
    """
    PDM_isosurface_n_part_set(self._isos, n_part)

  def n_part_out_set(self, n_part_out):
    """
    n_part_out_set(n_part)

    Set the number of partitions in the isosurface mesh (Optional).

    .. warning:: This function must be called prior to :py:func:`compute`

    .. note::
      By default, the number of partitions in the isosurface mesh is set to:

        - 1 in `REEQUILIBRATE` mode
        - the number of partitions in the source mesh in `LOCAL` mode (mandatory)

    Parameters:
      n_part_out (int) : Number of partitions
    """
    PDM_isosurface_n_part_out_set(self._isos, n_part_out)

  def connectivity_set(self,                     i_part,
                                                 connectivity_type,
      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connectivity_idx,
      NPY.ndarray[NPY.int32_t, mode='c', ndim=1] connectivity):
    """
    connectivity_set(i_part, connectivity_type, connectivity_idx, connectivity)

    Set partition connectivity.

    Parameters:
      i_part            (int)                     : Partition id
      connectivity_type (PDM_connectivity_type_t) : Connectivity type
      connectivity_idx  (np.ndarray[np.int32_t])  : Connectivity index
      connectivity      (np.ndarray[np.int32_t])  : Connectivity
    """
    # self.keep_alive.append(connectivity_idx)
    # self.keep_alive.append(connectivity)

    cdef int n_entity = 0
    if connectivity_type == _PDM_CONNECTIVITY_TYPE_EDGE_VTX:
      n_entity = connectivity.size // 2
    else:
      n_entity = connectivity_idx.size - 1

    cdef int *connect_idx_data = np_to_int_pointer(connectivity_idx)
    cdef int *connect_data     = np_to_int_pointer(connectivity)

    PDM_isosurface_connectivity_set(self._isos, i_part,
                                    connectivity_type, 
                                    n_entity,
                                    connect_idx_data,
                                    connect_data)

  def coordinates_set(self,                       i_part,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] coordinates):
    """
    coordinates_set(i_part, coordinates)

    Set partition coordinates.

    Parameters:
      i_part      (int)                     : Partition id
      coordinates (np.ndarray[np.double_t]) : Coordinates
    """
    # self.keep_alive.append(coordinates)

    cdef int     n_vtx      = coordinates.size // 3
    cdef double *coord_data = np_to_double_pointer(coordinates)

    PDM_isosurface_vtx_coord_set(self._isos, i_part, n_vtx, coord_data)

  def ln_to_gn_set(self,                            i_part,
                                                    entity_type,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] ln_to_gn):
    """
    ln_to_gn_set(i_part, entity_type, ln_to_gn)

    Set global numbering.

    Parameters:
      i_part      (int)                        : Partition id
      entity_type (PDM_mesh_entities_t)        : Entity type
      ln_to_gn    (np.ndarray[npy_pdm_gnum_t]) : Global ids
    """
    # self.keep_alive.append(ln_to_gn)
    cdef PDM_g_num_t *ln_to_gn_data = np_to_gnum_pointer(ln_to_gn)

    PDM_isosurface_ln_to_gn_set(self._isos, i_part,
                                entity_type, 
                                ln_to_gn_data)

  def group_set(self,                               i_part,
                                                    entity_type,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity_idx,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] group_ln_to_gn):
    """
    group_set(i_part, entity_type, group_entity_idx, group_entity, group_ln_to_gn)

    Set partition groups.

    Parameters:
      i_part           (int)                        : Partition id
      entity_type      (PDM_mesh_entities_t)        : Entity type
      group_entity_idx (np.ndarray[np.int_t])       : Groups index
      group_entity     (np.ndarray[np.int_t])       : Group entities
      group_ln_to_gn   (np.ndarray[npy_pdm_gnum_t]) : Group entities global ids
    """
    # self.keep_alive.append(group_entity_idx)
    # self.keep_alive.append(group_entity)
    # self.keep_alive.append(group_ln_to_gn)

    cdef int         *group_entity_idx_data = np_to_int_pointer (group_entity_idx)
    cdef int         *group_entity_data     = np_to_int_pointer (group_entity)
    cdef PDM_g_num_t *group_ln_to_gn_data   = np_to_gnum_pointer(group_ln_to_gn)

    cdef int n_group = group_entity_idx.size - 1

    PDM_isosurface_n_group_set(self._isos,
                               entity_type,
                               n_group)

    PDM_isosurface_group_set(self._isos, i_part,
                             entity_type,
                             group_entity_idx_data,
                             group_entity_data,
                             group_ln_to_gn_data)

  def part_mesh_set(self, part_mesh):
    """
    part_mesh_set(part_mesh)

    Set PDM_part_mesh.

    .. warning:: Not implemented yet.

    Parameters:
      part_mesh (PDM_part_mesh) : PDM_part_mesh
    """
    raise NotImplementedError() #PartMesh not defined yet
    # PDM_isosurface_part_mesh_set(self._isos, part_mesh._pm)

  def part_mesh_nodal_set(self, PartMeshNodal part_mesh_nodal):
    """
    part_mesh_nodal_set(part_mesh_nodal)

    Set PDM_part_mesh_nodal.

    Parameters:
      part_mesh_nodal (PDM_part_mesh_nodal) : PDM_part_mesh_nodal
    """
    PDM_isosurface_part_mesh_nodal_set(self._isos, part_mesh_nodal.pmn)

  def redistribution_set(self,
                         PDM_extract_part_kind_t extract_kind,
                         PDM_split_dual_t        part_method):
    """
    redistribution_set(extract_kind, part_method)

    Set reequilibrate strategy and repartionning tool.

    Parameters:
      extract_kind (PDM_extract_part_kind_t) : PDM_extract_part_kind_t
      part_method  (PDM_split_dual_t       ) : PDM_split_dual_t
    """
    PDM_isosurface_redistribution_set(self._isos, extract_kind, part_method)

  def field_set(self,                             id_iso,
                                                  i_part,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
    """
    field_set(id_iso, i_part, field)

    Set partition field.

    Parameters:
      id_iso      (int)                     : Isosurface id
      i_part      (int)                     : Partition id
      coordinates (np.ndarray[np.double_t]) : Field
    """
    cdef double *field_data = np_to_double_pointer(field)

    PDM_isosurface_field_set(self._isos, id_iso, i_part, field_data)

  # > Distributed setter API
  def dconnectivity_set(self,                       connectivity_type,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connectivity_idx,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connectivity):
    """
    dconnectivity_set(connectivity_type, connectivity_idx, connectivity)

    Set distributed connectivity.

    Parameters:
      connectivity_type (PDM_connectivity_type_t) : Connectivity type
      connectivity_idx  (np.ndarray[np.int32_t])  : Distributed connectivity index
      connectivity      (np.ndarray[np.int32_t])  : Distributed connectivity
    """
    cdef int         *connec_idx_data = np_to_int_pointer (connectivity_idx)
    cdef PDM_g_num_t *connec_data     = np_to_gnum_pointer(connectivity)

    PDM_isosurface_dconnectivity_set(self._isos, connectivity_type,
                                     connec_idx_data, connec_data)

  def dcoordinates_set(self,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] coordinates):
    """
    dcoordinates_set(coordinates)

    Set distributed coordinates.

    Parameters:
      coordinates (np.ndarray[np.double_t]) : Distributed coordinates
    """
    cdef double *coord_data = np_to_double_pointer(coordinates)

    PDM_isosurface_dvtx_coord_set(self._isos, coord_data)

  def distribution_set(self,                        entity_type,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] distribution):
    """
    distribution_set(entity_type, distribution)

    Set entity distribution.

    Parameters:
      distribution (np.ndarray[npy_pdm_gnum_t]) : Entity distribution
    """
    cdef PDM_g_num_t *distrib_data = np_to_gnum_pointer(distribution)

    PDM_isosurface_distrib_set(self._isos, entity_type, distrib_data)

  def dgroup_set(self,                              entity_type,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity_idx,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] group_entity):
    """
    dgroup_set(entity_type, group_entity_idx, group_entity)

    Set distributed groups.

    Parameters:
      entity_type      (PDM_mesh_entities_t)        : Entity type
      group_entity_idx (np.ndarray[np.int_t])       : Groups index
      group_entity     (np.ndarray[npy_pdm_gnum_t]) : Distributed group entities global ids
    """
    cdef int         *group_entity_idx_data = np_to_int_pointer (group_entity_idx)
    cdef PDM_g_num_t *group_entity_data     = np_to_gnum_pointer(group_entity)

    cdef int n_group = group_entity_idx.size - 1

    PDM_isosurface_n_group_set(self._isos,
                               entity_type,
                               n_group)

    PDM_isosurface_dgroup_set(self._isos, entity_type,
                              group_entity_idx_data, group_entity_data)

  def dmesh_set(self, DMesh dmesh):
    """
    dmesh_set(dmesh)

    Set PDM_dmesh.

    Parameters:
      dmesh (PDM_dmesh) : PDM_dmesh
    """
    PDM_isosurface_dmesh_set(self._isos, dmesh._dm)

  def dmesh_nodal_set(self, DMeshNodal dmesh_nodal):
    """
    dmesh_nodal_set(dmesh_nodal)

    Set PDM_dmesh_nodal.

    Parameters:
      dmesh_nodal (PDM_dmesh_nodal) : PDM_dmesh_nodal
    """
    PDM_isosurface_dmesh_nodal_set(self._isos, dmesh_nodal.dmn)

  def dfield_set(self,                            id_iso,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
    """
    dfield_set(id_iso, field)

    Set distributed field.

    Parameters:
      id_iso      (int)                     : Isosurface id
      coordinates (np.ndarray[np.double_t]) : Field
    """
    cdef double *field_data = np_to_double_pointer(field)

    PDM_isosurface_dfield_set(self._isos, id_iso, field_data)


  # > Partitioned getter API
  def connectivity_get(self, id_iso, i_part, connectivity_type):
    """
    connectivity_get(id_iso, i_part, connectivity_type)

    Get isosurface connectivity.

    Parameters:
      id_iso            (int)                     : Isosurface id
      i_part            (int)                     : Partition id
      connectivity_type (PDM_connectivity_type_t) : Connectivity type

    Returns:
      connectivity_idx (`np.ndarray[np.int32_t]`) : Connectivity index
      connectivity     (`np.ndarray[np.int32_t]`) : Connectivity
    """
    cdef int  n_entity         = 0
    cdef int *connectivity_idx = NULL
    cdef int *connectivity     = NULL
    n_entity = PDM_isosurface_connectivity_get(self._isos, id_iso, i_part,
                                               connectivity_type,
                                              &connectivity_idx,
                                              &connectivity,
                                               PDM_OWNERSHIP_USER)
    if connectivity_type==PDM_CONNECTIVITY_TYPE_EDGE_VTX:
      connectivity_size = 2*n_entity
    elif connectivity_type==PDM_CONNECTIVITY_TYPE_FACE_VTX:
      connectivity_size = connectivity_idx[n_entity]
    else:
      raise ValueError(f"PDM_isosurface_t: has no connectivity of type {connectivity_type}")

    np_connectivity_idx = create_numpy_i(connectivity_idx, n_entity+1       , flag_owndata=True)
    np_connectivity     = create_numpy_i(connectivity    , connectivity_size, flag_owndata=True)

    return np_connectivity_idx, np_connectivity

  def coordinates_get(self, id_iso, i_part):
    """
    coordinates_get(id_iso, i_part)

    Get isosurface coordinates.

    Parameters:
      id_iso (int) : Isosurface id
      i_part (int) : Partition id

    Returns:
      coordinates (`np.ndarray[np.double_t]`) : Coordinates
    """
    cdef int     n_vtx       = 0
    cdef double *coordinates = NULL
    n_vtx = PDM_isosurface_vtx_coord_get(self._isos, id_iso, i_part,
                                        &coordinates, PDM_OWNERSHIP_USER)

    np_coordinates = create_numpy_d(coordinates, 3*n_vtx, flag_owndata=True)

    return np_coordinates

  def ln_to_gn_get(self, id_iso, i_part, entity_type):
    """
    ln_to_gn_get(id_iso, i_part, entity_type)

    Get entity_type isosurface global ids.

    Parameters:
      id_iso      (int)                 : Isosurface id
      i_part      (int)                 : Partition id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      ln_to_gn (`np.ndarray[npy_pdm_gnum_t]`) : Global ids
    """
    cdef int          n_entity = 0
    cdef PDM_g_num_t *ln_to_gn = NULL
    n_entity = PDM_isosurface_ln_to_gn_get(self._isos, id_iso, i_part,
                                           entity_type,
                                          &ln_to_gn,
                                           PDM_OWNERSHIP_USER)

    np_ln_to_gn = create_numpy_g(ln_to_gn, n_entity, flag_owndata=True)

    return np_ln_to_gn

  def group_get(self, id_iso, i_part, entity_type):
    """
    group_get(id_iso, i_part, entity_type)

    Get entity_type isosurface groups.

    Parameters:
      id_iso      (int)                 : Isosurface id
      i_part      (int)                 : Partition id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      group_entity_idx (`np.ndarray[np.int32_t]`)     : Group index
      group_entity     (`np.ndarray[np.int32_t]`)     : Group entities
      group_ln_to_gn   (`np.ndarray[npy_pdm_gnum_t]`) : Group entities global ids
    """
    cdef int          n_group          = 0
    cdef int         *group_entity_idx = NULL
    cdef int         *group_entity     = NULL
    cdef PDM_g_num_t *group_ln_to_gn   = NULL
    n_group = PDM_isosurface_group_get(self._isos, id_iso, i_part,
                                           entity_type,
                                          &group_entity_idx,
                                          &group_entity,
                                          &group_ln_to_gn,
                                           PDM_OWNERSHIP_USER)
    group_entity_size = group_entity_idx[n_group]
    np_group_entity_idx = create_numpy_i(group_entity_idx, n_group+1        , flag_owndata=True)
    np_group_entity     = create_numpy_i(group_entity    , group_entity_size, flag_owndata=True)
    np_group_ln_to_gn   = create_numpy_g(group_ln_to_gn  , group_entity_size, flag_owndata=True)

    return np_group_entity_idx, np_group_entity, np_group_ln_to_gn

  def part_to_part_enable(self, id_iso, entity_type,
                          bint unify_parent_info=False):
    """
    Enable construction of a communication graph between source mesh entities and iso-surface entities.

    .. warning:: This function must be called prior to :py:func:`compute`

    Parameters:
      id_iso            (int)               : Isosurface id
      entity_type       (PDM_entity_type_t) : Entity type
      unify_parent_info (bool)              : Get all parent over all procs (not implemented)
    """
    PDM_isosurface_part_to_part_enable(self._isos, id_iso, entity_type, unify_parent_info)

  def part_to_part_get(self, id_iso, entity_type):
    """
    part_to_part_get(id_iso, entity_type)

    Get entity_type PDM_part_to_part between entry partitioned mesh and isosurface.

    .. note:: This function must be called prior to :py:func:`compute`

    Parameters:
      id_iso      (int)                 : Isosurface id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      PDM_part_to_part (PDM_part_to_part) : PDM_part_to_part
    """
    cdef PDM_part_to_part_t *ptp
    try:
      return self.ptp_entity[id_iso][entity_type]
    except KeyError:
      PDM_isosurface_part_to_part_get(self._isos, id_iso, entity_type,
                                     &ptp,
                                      PDM_OWNERSHIP_USER)
      py_caps_ptp = PyCapsule_New(ptp, NULL, NULL)
      self.ptp_entity[id_iso][entity_type] = PartToPartCapsule(py_caps_ptp, self.py_comm)
      return self.ptp_entity[id_iso][entity_type]

  def parent_lnum_get(self, id_iso, i_part, entity_type):
    """
    parent_lnum_get(id_iso, i_part, entity_type)

    Get entity_type isosurface entities parent local ids.

    Parameters:
      id_iso      (int)                 : Isosurface id
      i_part      (int)                 : Partition id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      parent_idx  (`np.ndarray[np.int32_t]`)     : Parent index
      parent_lnum (`np.ndarray[npy_pdm_gnum_t]`) : Parent entities local ids
    """
    cdef int  n_entity    = 0
    cdef int *parent_idx  = NULL
    cdef int *parent_lnum = NULL
    n_entity = PDM_isosurface_local_parent_get(self._isos, id_iso, i_part,
                                               entity_type,
                                              &parent_idx,
                                              &parent_lnum,
                                               PDM_OWNERSHIP_USER)

    own_idx = not self.got_parent_idx[entity_type][id_iso]
    self.got_parent_idx[entity_type][id_iso] = True

    parent_lnum_size = parent_idx[n_entity]
    np_parent_idx    = create_numpy_i(parent_idx , n_entity+1      , flag_owndata=own_idx)
    np_parent_lnum   = create_numpy_i(parent_lnum, parent_lnum_size, flag_owndata=True)

    return np_parent_idx, np_parent_lnum

  def parent_weight_get(self, id_iso, i_part, entity_type):
    """
    parent_weight_get(id_iso, i_part, entity_type)

    Get isosurface entity parent interpolation weight.

    Parameters:
      id_iso      (int)                 : Isosurface id
      i_part      (int)                 : Partition id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      parent_idx    (`np.ndarray[np.int32_t]`)  : Entity parent index
      parent_weight (`np.ndarray[np.double_t]`) : Entity parent weight
    """
    cdef int     n_entity      = 0
    cdef int    *parent_idx    = NULL
    cdef double *parent_weight = NULL

    n_entity = PDM_isosurface_parent_weight_get(self._isos, id_iso, i_part, entity_type,
                                               &parent_idx,
                                               &parent_weight,
                                                PDM_OWNERSHIP_USER)

    own_idx = not self.got_parent_idx[entity_type][id_iso]
    self.got_parent_idx[entity_type][id_iso] = True

    parent_weight_size = parent_idx[n_entity]
    np_parent_idx    = create_numpy_i(parent_idx,    n_entity+1        , flag_owndata=own_idx)
    np_parent_weight = create_numpy_d(parent_weight, parent_weight_size, flag_owndata=True)

    return np_parent_idx, np_parent_weight

  def isovalue_idx_get(self, id_iso, i_part, entity_type):
    """
    isovalue_idx_get(id_iso, i_part, entity_type)

    Get isosurface's isovalues index which indicate
    which entities belongs to which isovalues.

    Parameters:
      id_iso      (int)                 : Isosurface id
      i_part      (int)                 : Partition id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      isovalue_idx (`np.ndarray[np.int32_t]`) : Isovalue index
    """
    cdef int  n_isovalues  = 0
    cdef int *isovalue_idx = NULL
    n_isovalues = PDM_isosurface_isovalue_entity_idx_get(self._isos, id_iso, i_part,
                                                         entity_type,
                                                        &isovalue_idx,
                                                         PDM_OWNERSHIP_USER)

    np_isovalue_idx = create_numpy_i(isovalue_idx, n_isovalues+1, flag_owndata=True)

    return np_isovalue_idx


  # > Distributed getter API
  def dconnectivity_get(self, id_iso, connectivity_type):
    """
    dconnectivity_get(id_iso, connectivity_type)

    Get isosurface distributed connectivity.

    Parameters:
      id_iso            (int)                     : Isosurface id
      connectivity_type (PDM_connectivity_type_t) : Connectivity type

    Returns:
      connectivity_idx (`np.ndarray[np.int32_t]`)        : Distributed connectivity index
      connectivity     (`np.ndarray[np.npy_pdm_gnum_t]`) : Distributed connectivity
    """
    cdef int          dn_entity         = 0
    cdef int         *dconnectivity_idx = NULL
    cdef PDM_g_num_t *dconnectivity     = NULL
    dn_entity = PDM_isosurface_dconnectivity_get(self._isos, id_iso,
                                                 connectivity_type, 
                                                &dconnectivity_idx,
                                                &dconnectivity,
                                                 PDM_OWNERSHIP_USER)

    dconnectivity_size   = dconnectivity_idx[dn_entity]
    np_dconnectivity_idx = create_numpy_i(dconnectivity_idx, dn_entity+1       , flag_owndata=True)
    np_dconnectivity     = create_numpy_g(dconnectivity    , dconnectivity_size, flag_owndata=True)

    return np_dconnectivity_idx, np_dconnectivity

  def dcoordinates_get(self, id_iso):
    """
    dcoordinates_get(id_iso)

    Get isosurface distributed coordinates.

    Parameters:
      id_iso (int) : Isosurface id

    Returns:
      coordinates (`np.ndarray[np.double_t]`) : Distributed coordinates
    """
    cdef int     dn_vtx       = 0
    cdef double *dcoordinates = NULL
    dn_vtx = PDM_isosurface_dvtx_coord_get(self._isos, id_iso,
                                          &dcoordinates, PDM_OWNERSHIP_USER)

    np_dcoordinates = create_numpy_d(dcoordinates, 3*dn_vtx, flag_owndata=True)

    return np_dcoordinates

  def distribution_get(self, id_iso, entity_type):
    """
    distribution_get(id_iso, entity_type)

    Get isosurface entity distribution.

    Parameters:
      id_iso      (int)                 : Isosurface id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      distribution (`np.ndarray[np.npy_pdm_gnum_t]`) : Entity distribution
    """
    cdef PDM_g_num_t *distrib = NULL
    PDM_isosurface_distrib_get(self._isos, id_iso, entity_type,
                              &distrib)

    np_distrib = create_numpy_g(distrib, self.n_rank, flag_owndata=True)

    return np_distrib

  def dgroup_get(self, id_iso, entity_type):
    """
    dgroup_get(id_iso, entity_type)

    Get entity_type isosurface distributed groups.

    Parameters:
      id_iso      (int)                 : Isosurface id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      group_entity_idx (`np.ndarray[np.int32_t]`)     : Group index
      group_entity     (`np.ndarray[npy_pdm_gnum_t]`) : Group entities global ids
    """
    cdef int          n_group           = 0
    cdef int         *dgroup_entity_idx = NULL
    cdef PDM_g_num_t *dgroup_entity     = NULL
    n_group = PDM_isosurface_dgroup_get(self._isos, id_iso,
                                        entity_type, 
                                       &dgroup_entity_idx,
                                       &dgroup_entity,
                                        PDM_OWNERSHIP_USER)
    dgroup_entity_size = dgroup_entity_idx[n_group]
    np_dgroup_entity_idx = create_numpy_i(dgroup_entity_idx, n_group+1         , flag_owndata=True)
    np_dgroup_entity     = create_numpy_g(dgroup_entity    , dgroup_entity_size, flag_owndata=True)

    return np_dgroup_entity_idx, np_dgroup_entity

  def dparent_weight_get(self, id_iso, entity_type):
    """
    dparent_weight_get(id_iso)

    Get isosurface entity parent interpolation weight.

    Parameters:
      id_iso      (int)                 : Isosurface id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      parent_idx    (`np.ndarray[np.int32_t]`)  : Entity parent index
      parent_weight (`np.ndarray[np.double_t]`) : Entity parent weights
    """
    cdef int     n_entity      = 0
    cdef int    *parent_idx    = NULL
    cdef double *parent_weight = NULL

    n_entity = PDM_isosurface_dparent_weight_get(self._isos, id_iso, entity_type,
                                                &parent_idx,
                                                &parent_weight,
                                                 PDM_OWNERSHIP_USER)

    parent_weight_size = parent_idx[n_entity]
    np_parent_weight = create_numpy_d(parent_weight, parent_weight_size, flag_owndata=True)

    return np_parent_weight

  def disovalue_entity_get(self, id_iso, entity_type):
    """
    disovalue_entity_get(id_iso, entity_type)

    Get isosurface's isovalue id of entities which indicate
    which entities belongs to which isovalues.

    Parameters:
      id_iso      (int)                 : Isosurface id
      entity_type (PDM_mesh_entities_t) : Entity type

    Returns:
      disovalue_entity_idx (`np.ndarray[np.int32_t]`)     : Isovalue index
      disovalue_entity     (`np.ndarray[npy_pdm_gnum_t]`) : Isovalue entities global ids
    """
    cdef int          n_isovalues          = 0
    cdef int         *disovalue_entity_idx = NULL
    cdef PDM_g_num_t *disovalue_entity     = NULL

    n_isovalues = PDM_isosurface_disovalue_entity_get(self._isos, id_iso, entity_type,
                                                     &disovalue_entity_idx,
                                                     &disovalue_entity,
                                                      PDM_OWNERSHIP_USER)

    disovalue_entity_size   = disovalue_entity_idx[n_isovalues]
    np_disovalue_entity_idx = create_numpy_i(disovalue_entity_idx, n_isovalues+1        , flag_owndata=True)
    np_disovalue_entity     = create_numpy_g(disovalue_entity    , disovalue_entity_size, flag_owndata=True)

    return np_disovalue_entity_idx, np_disovalue_entity


  def dump_times(self):
    """
    Dump elapsed and CPU times
    """
    PDM_isosurface_dump_times(self._isos)

  def __dealloc__(self):
    """
    """
    PDM_isosurface_free(self._isos)
