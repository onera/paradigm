
cdef extern from "pdm_isosurface.h":
  ctypedef struct PDM_isosurface_t:
      pass

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
                                       int                     *connect_idx,
                                       int                     *connect);

  void PDM_isosurface_vtx_coord_set(PDM_isosurface_t *isos,
                                    int               i_part,
                                    double           *vtx_coord);

  void PDM_isosurface_ln_to_gn_set(PDM_isosurface_t    *isos,
                                   int                  i_part,
                                   PDM_mesh_entities_t  entity_type,
                                   int                  n_entity,
                                   PDM_g_num_t         *ln_to_gn);

  void PDM_isosurface_group_set(PDM_isosurface_t    *isos,
                                int                  i_part,
                                PDM_mesh_entities_t  entity_type,
                                int                  n_group,
                                int                 *group_entity_idx,
                                int                 *group_entity,
                                PDM_g_num_t         *group_entity_ln_to_gn);

  void PDM_isosurface_part_mesh_set(PDM_isosurface_t *isos,
                                    PDM_part_mesh_t  *pmesh);

  void PDM_isosurface_mesh_nodal_set(PDM_isosurface_t      *isos,
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
                                 int                  n_group,
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

  void PDM_isosurface_set_isovalues(PDM_isosurface_t *isos,
                                    int               id_isosurface,
                                    int               n_isovalues,
                                    double           *isovalues);

  void PDM_isosurface_equation_set(PDM_isosurface_t *isos,
                                   int               id_isosurface,
                                   double           *coeff,
                                   int               use_gradient);

  void PDM_isosurface_field_function_set(PDM_isosurface_t                *isos,
                                         int                              id_isosurface,
                                         PDM_isosurface_field_function_t  func);

  void PDM_isosurface_field_set(PDM_isosurface_t *isos,
                                int               id_isosurface,
                                int               i_part,
                                double           *field);

  void PDM_isosurface_gradient_set(PDM_isosurface_t *isos,
                                   int               id_isosurface,
                                   int               i_part,
                                   double           *gradient);

  void PDM_isosurface_dfield_set(PDM_isosurface_t *isos,
                                 int               id_isosurface,
                                 double           *dfield);

  void PDM_isosurface_dgradient_set(PDM_isosurface_t *isos,
                                    int               id_isosurface,
                                    double           *dgradient);


  void PDM_isosurface_redistribution_set(PDM_isosurface_t        *isos,
                                         PDM_extract_part_kind_t  extract_kind,
                                         PDM_split_dual_t         part_method);


  void PDM_isosurface_reset(PDM_isosurface_t *isos,
                            int               id_isosurface);


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

  int PDM_isosurface_dconnectivity_get(PDM_isosurface_t         *isos,
                                       int                       id_isosurface,
                                       PDM_connectivity_type_t   connectivity_type,
                                       int                     **dconnect_idx,
                                       PDM_g_num_t             **dconnect,
                                       PDM_ownership_t           ownership);

  int PDM_isosurface_parent_gnum_get(PDM_isosurface_t     *isos,
                                     int                   id_iso,
                                     PDM_mesh_entities_t   entity_type,
                                     int                 **parent_idx,
                                     PDM_g_num_t         **parent_gnum,
                                     PDM_ownership_t       ownership);

  void PDM_isosurface_dvtx_parent_weight_get(PDM_isosurface_t     *isos,
                                             int                   id_iso,
                                             double              **dvtx_parent_weight,
                                             PDM_ownership_t       ownership);

  int PDM_isosurface_dvtx_protocol_get(PDM_isosurface_t *isos);

  int PDM_isosurface_dvtx_coord_get(PDM_isosurface_t  *isos,
                                    int                id_isosurface,
                                    double           **dvtx_coord,
                                    PDM_ownership_t    ownership);

  void PDM_isosurface_distrib_get(PDM_isosurface_t     *isos,
                                  int                   id_isosurface,
                                  PDM_mesh_entities_t   entity_type,
                                  PDM_g_num_t         **distribution,
                                  PDM_ownership_t       ownership);

  int PDM_isosurface_dgroup_get(PDM_isosurface_t     *isos,
                                int                   id_isosurface,
                                PDM_mesh_entities_t   entity_type,
                                int                 **dgroup_entity_idx,
                                PDM_g_num_t         **dgroup_entity,
                                PDM_ownership_t       ownership);

  int PDM_isosurface_local_parent_get(PDM_isosurface_t     *isos,
                                      int                   id_isosurface,
                                      int                   i_part,
                                      PDM_mesh_entities_t   entity_type,
                                      int                 **entity_parent_idx,
                                      int                 **entity_parent,
                                      PDM_ownership_t       ownership);

  int PDM_isosurface_vtx_parent_weight_get(PDM_isosurface_t  *isos,
                                           int                id_isosurface,
                                           int                i_part,
                                           int              **vtx_parent_idx,
                                           double           **vtx_parent_weight,
                                           PDM_ownership_t    ownership);

  void PDM_isosurface_enable_part_to_part(PDM_isosurface_t     *isos,
                                          int                   id_isosurface,
                                          PDM_mesh_entities_t   entity_type,
                                          int                   unify_parent_info);

  void PDM_isosurface_part_to_part_get(PDM_isosurface_t     *isos,
                                       int                   id_isosurface,
                                       PDM_mesh_entities_t   entity_type,
                                       PDM_part_to_part_t  **ptp,
                                       PDM_ownership_t       ownership);

  void PDM_isosurface_free(PDM_isosurface_t  *isos);

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class Isosurface:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_isosurface_t *_isos
  
  cdef dict ptp_entity

  # ************************************************************************

  # > Common setter API
  def __cinit__(self, int      mesh_dim, 
                      MPI.Comm comm):
    """
    __init__(mesh_dim, comm)

    Create the structure.

    Parameters:
      mesh_dim (int)      : Entry mesh dimension
      comm     (MPI.Comm) : MPI communicator
    """

    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
    cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    self._isos = PDM_isosurface_create(pdm_comm,
                                        mesh_dim)

  def tolerance_set(self, tolerance):
    """
    set_tolerance(tolerance)

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
    cdef int  id_iso         = -1
    cdef int  n_isovalues    = len(isovalues)
    cdef int *isovalues_data = list_to_double_pointer(isovalues)

    id_iso = PDM_isosurface_add(self._isos, kind, n_isovalues, isovalues_data)

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
    cdef int  n_isovalues    = len(isovalues)
    cdef int *isovalues_data = list_to_double_pointer(isovalues)

    PDM_isosurface_set_isovalues(self._isos, id_iso, n_isovalues, isovalues_data)

  def field_function_set(self):
    """
    Not implemented.
    """
    raise NotImplementedError()

  def compute(self, id_iso):
    """
    compute(id_iso)

    Computed aksed isosurface.

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
  def n_part(self, n_part):
    """
    n_part(n_part)

    Set entry mesh number or partitions.

    Parameters:
      n_part (int) : Number of partitions
    """
    PDM_isosurface_n_part_set(self._isos, n_part)

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
    PDM_isosurface_connectivity_set(self._isos, i_part,
                                    connectivity_type, 
                                    connectivity_idx, connectivity)

  def coordinates_set(self,                       i_part,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] coordinates):
    """
    coordinates_set(i_part, coordinates)

    Set partition coordinates.

    Parameters:
      i_part      (int)                     : Partition id
      coordinates (np.ndarray[np.double_t]) : Coordinates
    """
    PDM_isosurface_vtx_coord_set(self._isos, i_part, coordinates)

  def ln_to_gn_set(self,                            i_part,
                                                    entity_type,
                                                    n_entity,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] ln_to_gn):
    """
    ln_to_gn_set(i_part, entity_type, n_entity, ln_to_gn)

    Set global numbering.

    Parameters:
      i_part      (int)                        : Partition id
      entity_type (PDM_entity_type_t)          : Entity type
      n_entity    (int)                        : Number of entities
      ln_to_gn    (np.ndarray[npy_pdm_gnum_t]) : Global ids
    """
    PDM_isosurface_ln_to_gn_set(self._isos, i_part,
                                entity_type, 
                                n_entity, ln_to_gn)

  def group_set(self,                               i_part,
                                                    entity_type,
                                                    n_group,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity_idx
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] group_ln_to_gn):
    """
    group_set(i_part, entity_type,
              n_group, group_entity_idx, group_entity, group_ln_to_gn)

    Set partition groups.

    Parameters:
      i_part           (int)                        : Partition id
      entity_type      (PDM_entity_type_t)          : Entity type
      n_entity         (int)                        : Number of entities
      group_entity_idx (np.ndarray[np.int_t])       : Groups index
      group_entity     (np.ndarray[np.int_t])       : Group entities
      group_ln_to_gn   (np.ndarray[npy_pdm_gnum_t]) : Group entities global ids
    """
    PDM_isosurface_group_set(self._isos, i_part,
                             entity_type,
                             n_group,
                             group_entity_idx,
                             group_entity,
                             group_ln_to_gn)

  def part_mesh_set(self, part_mesh):
    """
    part_mesh_set(part_mesh)

    Set PDM_part_mesh.

    Parameters:
      part_mesh (PDM_part_mesh) : PDM_part_mesh
    """
    PDM_isosurface_part_mesh_set(self._isos, part_mesh)

  def part_mesh_nodal_set(self, part_mesh_nodal):
    """
    part_mesh_nodal_set(part_mesh_nodal)

    Set PDM_part_mesh_nodal.

    Parameters:
      part_mesh_nodal (PDM_part_mesh_nodal) : PDM_part_mesh_nodal
    """
    PDM_isosurface_mesh_nodal_set(self._isos, part_mesh_nodal)

  def redistribution_set(self, extract_kind, part_method):
    """
    redistribution_set(extract_kind, part_method)

    Set reequilibrate strategy and repartionning tool.

    Parameters:
      extract_kind (PDM_extract_part_kind_t) : PDM_extract_part_kind_t
      part_method  (PDM_split_dual_t       ) : PDM_split_dual_t       
    """

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
    PDM_isosurface_field_set(self._isos, i_part, coordinates)

  def gradient_set(self,                          id_iso,
                                                  i_part,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] gradient):
    """
    Not implemented.
    """
    raise NotImplementedError()


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
    PDM_isosurface_connectivity_set(self._isos, connectivity_type, 
                                    connectivity_idx, connectivity)

  def dcoordinates_set(self,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] coordinates):
    """
    dcoordinates_set(coordinates)

    Set distributed coordinates.

    Parameters:
      coordinates (np.ndarray[np.double_t]) : Distributed coordinates
    """
    PDM_isosurface_dvtx_coord_set(self._isos, coordinates)

  def distribution_set(self,                        entity_type,
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] distribution):
    """
    distribution_set(entity_type, distribution)

    Set entity distribution.

    Parameters:
      distribution (np.ndarray[np.double_t]) : Entity distribution
    """
    PDM_isosurface_distrib_set(self._isos, entity_type, distribution)

  def dgroup_set(self,                              entity_type,
                                                    n_group,
      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] group_entity_idx
      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] group_entity):
    """
    dgroup_set(entity_type, n_group, group_entity_idx, group_entity)

    Set distributed groups.

    Parameters:
      entity_type      (PDM_entity_type_t)          : Entity type
      n_entity         (int)                        : Number of groups
      group_entity_idx (np.ndarray[np.int_t])       : Groups index
      group_entity     (np.ndarray[npy_pdm_gnum_t]) : Distributed group entities global ids
    """
    PDM_isosurface_dgroup_set(self._isos, entity_type,
                              n_group, group_entity_idx, group_entity)

  def dmesh_set(self, dmesh):
    """
    dmesh_set(dmesh)

    Set PDM_dmesh.

    Parameters:
      dmesh (PDM_dmesh) : PDM_dmesh
    """
    PDM_isosurface_part_mesh_set(self._isos, part_mesh)

  def dmesh_nodal_set(self, dmesh_nodal):
    """
    dmesh_nodal_set(dmesh_nodal)

    Set PDM_dmesh_nodal.

    Parameters:
      dmesh_nodal (PDM_dmesh_nodal) : PDM_dmesh_nodal
    """
    PDM_isosurface_mesh_nodal_set(self._isos, dmesh_nodal)

  def dfield_set(self,                            id_iso,
                                                  i_part,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
    """
    dfield_set(id_iso, i_part, field)

    Set distributed field.

    Parameters:
      id_iso      (int)                     : Isosurface id
      i_part      (int)                     : Partition id
      coordinates (np.ndarray[np.double_t]) : Field
    """
    PDM_isosurface_field_set(self._isos, i_part, coordinates)

  def dgradient_set(self,                          id_iso,
      NPY.ndarray[NPY.double_t, mode='c', ndim=1] gradient):
    """
    Not implemented.
    """
    raise NotImplementedError()


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
    cdef int  n_entity         = 0;
    cdef int *connectivity_idx = NULL;
    cdef int *connectivity     = NULL;
    n_entity = PDM_isosurface_connectivity_get(self._isos, id_iso, i_part,
                                               connectivity_type, 
                                              &connect_idx,
                                              &connect,
                                               PDM_OWNERSHIP_USER)

    connectivity_size   = connect_idx[n_entity]
    np_connectivity_idx = create_numpy_i(connect_idx, n_entity         , flag_owndata=True)
    np_connectivity     = create_numpy_i(connect    , connectivity_size, flag_owndata=True)

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
    cdef int     n_vtx       = 0;
    cdef double *coordinates = NULL;
    n_vtx = PDM_isosurface_vtx_coord_get(self._isos, id_iso, i_part,
                                        &coordinates, PDM_OWNERSHIP_USER)

    np_coordinates = create_numpy_i(coordinates, n_vtx, flag_owndata=True)

    return np_coordinates

  def ln_to_gn_get(self, id_iso, i_part, entity_type):
    """
    ln_to_gn_get(id_iso, i_part, entity_type)

    Get entity_type isosurface global ids.

    Parameters:
      id_iso      (int)               : Isosurface id
      i_part      (int)               : Partition id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      ln_to_gn (`np.ndarray[npy_pdm_gnum_t]`) : Global ids
    """
    cdef int          n_entity = 0;
    cdef PDM_g_num_t *ln_to_gn = NULL;
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
      id_iso      (int)               : Isosurface id
      i_part      (int)               : Partition id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      n_group          (int)                          : Number of groups
      group_entity_idx (`np.ndarray[np.int32_t]`)     : Group index
      group_entity     (`np.ndarray[np.int32_t]`)     : Group entities
      group_ln_to_gn   (`np.ndarray[npy_pdm_gnum_t]`) : Group entities global ids
    """
    cdef int          n_group          = 0;
    cdef int         *group_entity_idx = NULL;
    cdef int         *group_entity     = NULL;
    cdef PDM_g_num_t *group_ln_to_gn   = NULL;
    n_group = PDM_isosurface_group_get(self._isos, id_iso, i_part,
                                           entity_type, 
                                          &group_entity_idx,
                                          &group_entity,
                                          &group_ln_to_gn,
                                           PDM_OWNERSHIP_USER)
    group_entity_size = group_entity_idx[n_group]
    np_group_entity_idx = create_numpy_i(group_entity_idx, n_group          , flag_owndata=True)
    np_group_entity     = create_numpy_i(group_entity    , group_entity_size, flag_owndata=True)
    np_group_ln_to_gn   = create_numpy_i(group_ln_to_gn  , group_entity_size, flag_owndata=True)

    return n_group, np_group_entity_idx, np_group_entity, np_group_ln_to_gn

  def part_to_part_get(self, id_iso, entity_type):
    """
    part_to_part_get(id_iso, entity_type)

    Get entity_type PDM_part_to_part between entry partitioned mesh and isosurface.

    Parameters:
      id_iso      (int)               : Isosurface id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      PDM_part_to_part (PDM_part_to_part) : PDM_part_to_part
    """
    cdef PDM_part_to_part_t *ptp
    try:
      return self.ptp_entity[entity_type]
    except KeyError:
      PDM_isosurface_part_to_part_get(self._isos, id_iso, entity_type,
                                     &ptp,
                                      PDM_OWNERSHIP_USER)
      py_caps_ptp = PyCapsule_New(ptp, NULL, NULL);
      self.ptp_entity[entity_type] = PartToPartCapsule(py_caps_ptp, self.py_comm)
      return self.ptp_entity[entity_type]

  def parent_lnum_get(self, id_iso, i_part, entity_type):
    """
    parent_lnum_get(id_iso, i_part, entity_type)

    Get entity_type isosurface entities parent local ids.

    Parameters:
      id_iso      (int)               : Isosurface id
      i_part      (int)               : Partition id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      parent_idx  (`np.ndarray[np.int32_t]`)     : Parent index
      parent_lnum (`np.ndarray[npy_pdm_gnum_t]`) : Parent entities local ids
    """
    cdef int  n_entity    = 0;
    cdef int *parent_idx  = NULL;
    cdef int *parent_lnum = NULL;
    n_entity = PDM_isosurface_local_parent_get(self._isos, id_iso,
                                               entity_type, 
                                              &parent_idx,
                                              &parent_lnum,
                                               PDM_OWNERSHIP_USER)
    parent_lnum_size = parent_idx[n_entity]
    np_parent_idx    = create_numpy_i(parent_idx , n_entity        , flag_owndata=True)
    np_parent_lnum   = create_numpy_i(parent_lnum, parent_lnum_size, flag_owndata=True)

    return np_parent_idx, np_parent_lnum


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
    cdef int  dn_entity         = 0;
    cdef int *dconnectivity_idx = NULL;
    cdef int *dconnectivity     = NULL;
    dn_entity = PDM_isosurface_dconnectivity_get(self._isos, id_iso, i_part,
                                                 connectivity_type, 
                                                &dconnectivity_idx,
                                                &dconnectivity,
                                                 PDM_OWNERSHIP_USER)

    dconnectivity_size   = donnectivity_idx[n_entity]
    np_dconnectivity_idx = create_numpy_i(dconnectivity_idx, n_entity          , flag_owndata=True)
    np_dconnectivity     = create_numpy_i(dconnectivity    , dconnectivity_size, flag_owndata=True)

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
    cdef int     dn_vtx       = 0;
    cdef double *dcoordinates = NULL;
    dn_vtx = PDM_isosurface_dvtx_coord_get(self._isos, id_iso, i_part,
                                          &dcoordinates, PDM_OWNERSHIP_USER)

    np_dcoordinates = create_numpy_i(dcoordinates, n_vtx, flag_owndata=True)

    return np_dcoordinates

  def distrib_get(self, id_iso, entity_type):
    """
    Not implemented
    """
    raise NotImplementedError()

  def dgroup_get(self, id_iso, entity_type):
    """
    dgroup_get(id_iso, entity_type)

    Get entity_type isosurface distributed groups.

    Parameters:
      id_iso      (int)               : Isosurface id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      n_group          (int)                          : Number of groups
      group_entity_idx (`np.ndarray[np.int32_t]`)     : Group index
      group_entity     (`np.ndarray[npy_pdm_gnum_t]`) : Group entities global ids
    """
    cdef int          n_group          = 0;
    cdef int         *dgroup_entity_idx = NULL;
    cdef PDM_g_num_t *dgroup_entity     = NULL;
    n_group = PDM_isosurface_group_get(self._isos, id_iso, i_part,
                                       entity_type, 
                                      &dgroup_entity_idx,
                                      &dgroup_entity,
                                       PDM_OWNERSHIP_USER)
    dgroup_entity_size = dgroup_entity_idx[n_group]
    np_dgroup_entity_idx = create_numpy_i(dgroup_entity_idx, n_group           , flag_owndata=True)
    np_dgroup_entity     = create_numpy_g(dgroup_entity    , dgroup_entity_size, flag_owndata=True)

    return n_group, np_dgroup_entity_idx, np_dgroup_entity

  def parent_gnum_get(self, id_iso, entity_type):
    """
    parent_gnum_get(id_iso, entity_type)

    Get entity_type isosurface entities parent global ids.

    Parameters:
      id_iso      (int)               : Isosurface id
      entity_type (PDM_entity_type_t) : Entity type

    Returns:
      parent_idx  (`np.ndarray[np.int32_t]`)     : Parent index
      parent_gnum (`np.ndarray[npy_pdm_gnum_t]`) : Parent entities global ids
    """
    cdef int          n_entity    = 0;
    cdef int         *parent_idx  = NULL;
    cdef PDM_g_num_t *parent_gnum = NULL;
    n_entity = PDM_isosurface_parent_gnum_get(self._isos, id_iso,
                                              entity_type, 
                                             &parent_idx,
                                             &parent_gnum,
                                              PDM_OWNERSHIP_USER)
    parent_gnum_size = parent_idx[n_entity]
    np_parent_idx    = create_numpy_i(parent_idx , n_entity        , flag_owndata=True)
    np_parent_gnum   = create_numpy_g(parent_gnum, parent_gnum_size, flag_owndata=True)

    return np_parent_idx, np_parent_gnum

  def dvtx_parent_weight_get(self, id_iso):
    """
    vtx_parent_weight_get(id_iso)

    Get isosurface vertices parent weight.

    Parameters:
      id_iso (int) : Isosurface id

    Returns:
      vtx_parent_weight_get (`np.ndarray[np.double_t]`) : Vertives parent weight
    """
    cdef int          n_vtx                 = 0;
    cdef int         *vtx_parent_idx        = NULL;
    cdef PDM_g_num_t *vtx_parent_gnum       = NULL;
    cdef double      *vtx_parent_weight_get = NULL;

    n_vtx = PDM_isosurface_parent_gnum_get(self._isos, id_iso,
                                           entity_type, 
                                          &vtx_parent_gnum,
                                          &vtx_parent_weight_get,
                                           PDM_OWNERSHIP_USER)
    PDM_isosurface_dvtx_parent_weight_get(self._isos, id_iso,
                                         &vtx_parent_weight_get,
                                          PDM_OWNERSHIP_USER)

    vtx_parent_gnum_size = vtx_idx[n_vtx]
    np_vtx_parent_weight_get = create_numpy_d(vtx_parent_weight_get, vtx_parent_gnum_size, flag_owndata=True)

    return np_vtx_parent_weight_get


  def __dealloc__(self):
    """
    """
    PDM_isosurface_free(self._isos)
