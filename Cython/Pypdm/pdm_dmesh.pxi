
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_t* PDM_dmesh_create(PDM_ownership_t owner,
                                  int             dn_cell,
                                  int             dn_face,
                                  int             dn_edge,
                                  int             dn_vtx,
                                  PDM_MPI_Comm    comm)

    # void PDM_dmesh_set(PDM_dmesh_t  *dm,
    #                    double       *dvtx_coord,
    #                    int          *dface_vtx_idx,
    #                    PDM_g_num_t  *dface_vtx,
    #                    PDM_g_num_t  *dface_cell,
    #                    int          *dface_bound_idx,
    #                    PDM_g_num_t  *dface_bound,
    #                    int          *join_g_dms,
    #                    int          *dface_join_idx,
    #                    PDM_g_num_t  *dface_join)

    # void PDM_dmesh_dims_get(PDM_dmesh_t *dm,
    #                         int         *dn_cell,
    #                         int         *dn_face,
    #                         int         *dn_edge,
    #                         int         *dn_vtx,
    #                         int         *n_bnd,
    #                         int         *n_joins)

    void PDM_dmesh_data_get(PDM_dmesh_t   *dm,
                            double       **dvtx_coord,
                            int          **dface_vtx_idx,
                            PDM_g_num_t  **dface_vtx,
                            PDM_g_num_t  **dface_cell,
                            int          **dface_bound_idx,
                            PDM_g_num_t  **dface_bound,
                            int          **join_g_dms,
                            int          **dface_join_idx,
                            PDM_g_num_t  **dface_join)

    int PDM_dmesh_dn_entity_get(PDM_dmesh_t         *dmesh,
                                PDM_mesh_entities_t  entity_type)
    int PDM_dmesh_connectivity_get(PDM_dmesh_t              *dmesh,
                                   PDM_connectivity_type_t   connectivity_type,
                                   PDM_g_num_t             **connect,
                                   int                     **connect_idx,
                                   PDM_ownership_t           ownership)
    int PDM_dmesh_distrib_get(PDM_dmesh_t              *dmesh,
                               PDM_mesh_entities_t       entity,
                               PDM_g_num_t             **distrib)
    int  PDM_dmesh_bound_get(PDM_dmesh_t       *dmesh,
                             PDM_bound_type_t   bound_type,
                             PDM_g_num_t      **connect,
                             int              **connect_idx,
                             PDM_ownership_t    ownership)
    void PDM_dmesh_vtx_coord_set(PDM_dmesh_t      *dmesh,
                                 double           *dvtx_coord,
                                 PDM_ownership_t   ownership)
    void PDM_dmesh_vtx_coord_get(PDM_dmesh_t      *dmesh,
                                 double          **dvtx_coord,
                                 PDM_ownership_t   ownership)
    void  PDM_dmesh_bound_set(PDM_dmesh_t       *dmesh,
                              PDM_bound_type_t   bound_type,
                              int               n_bound,
                              PDM_g_num_t       *connect,
                              int               *connect_idx,
                              PDM_ownership_t    ownership)
    void PDM_dmesh_connectivity_set(PDM_dmesh_t              *dmesh,
                                    PDM_connectivity_type_t   connectivity_type,
                                    PDM_g_num_t              *connect,
                                    int                      *connect_idx,
                                    PDM_ownership_t           ownership);
    void PDM_dmesh_free(PDM_dmesh_t   *dm);

    void PDM_dmesh_find_topological_ridges(PDM_MPI_Comm   comm,
                                           PDM_g_num_t   *distrib_face,
                                           int           *dface_vtx_idx,
                                           PDM_g_num_t   *dface_vtx,
                                           int            n_group_face,
                                           int           *dgroup_face_idx,
                                           PDM_g_num_t   *dgroup_face,
                                           PDM_g_num_t  **out_distrib_ridge,
                                           PDM_g_num_t  **out_dridge_vtx,
                                           int           *out_n_group_ridge,
                                           int          **out_dgroup_edge_idx,
                                           PDM_g_num_t  **out_dgroup_edge,
                                           int          **out_dridge_face_group_idx,
                                           int          **out_dridge_face_group);
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
# cdef void PDM_pydmesh_free(object caps):
#   # print("PDM_pydmesh_free", PyCapsule_GetName(caps))
#   cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, <const char*> PyCapsule_GetName(caps))
#   PDM_dmesh_free(dm);

# ------------------------------------------------------------------
cdef class DistributedMeshCapsule:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, object caps):
    """
    """
    # print("DistributedMeshCapsule", PyCapsule_GetName(caps))
    cdef PDM_dmesh_t* dm = <PDM_dmesh_t *> PyCapsule_GetPointer(caps, NULL)
    self._dm = dm;

  # ------------------------------------------------------------------------
  def dmesh_vtx_coord_set(self,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
      dmesh_vtx_coord_set(self, dvtx_coord)
  def dmesh_vtx_coord_get(self):
      return dmesh_vtx_coord_get(self)
  # ------------------------------------------------------------------------
  def dmesh_connectivity_set(self, PDM_connectivity_type_t connectivity_type,
                             NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_connectivity_set(self, connectivity_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_distrib_get(self, entity_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_set(self, PDM_bound_type_t bound_type,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_bound_set(self, bound_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    PDM_dmesh_free(self._dm)

# ------------------------------------------------------------------
cdef class DistributedMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  cdef int          n_rank
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, MPI.Comm comm,
                      dn_cell,
                      dn_face,
                      dn_edge,
                      dn_vtx):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    self.n_rank = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._dm = PDM_dmesh_create(PDM_OWNERSHIP_UNGET_RESULT_IS_FREE,
                                dn_cell,
                                dn_face,
                                dn_edge,
                                dn_vtx,
                                PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm))
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  # def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx    not None,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_cell,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_bound_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_bound,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] join_g_dms,
  #                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_join_idx,
  #                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_join):
  #   """
  #   """
  #   # print("dvtx_coord", dvtx_coord)
  #   # print("dface_vtx_idx", dface_vtx_idx)
  #   # print("dface_vtx", dface_vtx)
  #   # print("dface_cell", dface_cell)
  #   # print("dface_bound_idx", dface_bound_idx)
  #   # print("dface_bound", dface_bound)

  #   PDM_dmesh_set(self._dm,
  #                 <double*>      dvtx_coord.data,
  #                 <int*>         dface_vtx_idx.data,
  #                 <PDM_g_num_t*> dface_vtx.data,
  #                 <PDM_g_num_t*> dface_cell.data,
  #                 <int*>         dface_bound_idx.data,
  #                 <PDM_g_num_t*> dface_bound.data,
  #                 <int*>         join_g_dms.data,
  #                 <int*>         dface_join_idx.data,
  #                 <PDM_g_num_t*> dface_join.data)
  # ------------------------------------------------------------------------
  def dmesh_vtx_coord_set(self,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
      dmesh_vtx_coord_set(self, dvtx_coord)
  def dmesh_vtx_coord_get(self):
      return dmesh_vtx_coord_get(self)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_get(self, PDM_connectivity_type_t connectivity_type):
    """
    """
    return dmesh_connectivity_get(self, connectivity_type)

  # ------------------------------------------------------------------------
  def dmesh_connectivity_set(self, PDM_connectivity_type_t connectivity_type,
                             NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                             NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_connectivity_set(self, connectivity_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def dmesh_distrib_get(self, PDM_mesh_entities_t entity_type):
    """
    """
    return dmesh_connectivity_get(self, entity_type)
  # ------------------------------------------------------------------------
  def dmesh_bound_get(self, PDM_bound_type_t bound_type):
    """
    """
    return dmesh_bound_get(self, bound_type)

  # ------------------------------------------------------------------------
  def dmesh_bound_set(self, PDM_bound_type_t bound_type,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None):
    """
    """
    dmesh_bound_set(self, bound_type, connect_idx, connect)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    # print("DistributedMesh::__dealloc__")
    PDM_dmesh_free(self._dm)
    # print("DistributedMesh::__dealloc__")

ctypedef fused DMesh:
  DistributedMesh
  DistributedMeshCapsule

# ------------------------------------------------------------------------
def dmesh_connectivity_get(DMesh pydm, PDM_connectivity_type_t connectivity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  # ************************************************************************
  dn_entity = PDM_dmesh_connectivity_get(pydm._dm,
                                         connectivity_type,
                                         &connect,
                                         &connect_idx,
                                         PDM_OWNERSHIP_USER)

  np_connect_idx = create_numpy_or_none_i(connect_idx, dn_entity+1)

  if (connect == NULL) :
      np_connect = None
  else :
      if(np_connect_idx is not None):
        dim = <NPY.npy_intp> connect_idx[dn_entity]
      else:
        dim = <NPY.npy_intp> 2 * dn_entity # Face cell
      np_connect = create_numpy_g(connect, dim)
  return (np_connect_idx, np_connect)

# ------------------------------------------------------------------------
def dmesh_distrib_get(DMesh pydm, PDM_mesh_entities_t entity_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef PDM_g_num_t  *distrib
  cdef NPY.npy_intp  dim
  cdef int size
  # ************************************************************************
  size = PDM_dmesh_distrib_get(pydm._dm, entity_type, &distrib)
  return create_numpy_or_none_g(distrib, size+1, flag_owndata=False)

# ------------------------------------------------------------------------
def dmesh_bound_get(DMesh pydm, PDM_bound_type_t bound_type):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int          *connect_idx
  cdef PDM_g_num_t  *connect
  cdef int           dn_entity
  cdef NPY.npy_intp  dim
  cdef int           n_bnd
  # ************************************************************************
  n_bnd = PDM_dmesh_bound_get(pydm._dm,
                              bound_type,
                              &connect,
                              &connect_idx,
                              PDM_OWNERSHIP_USER)

  np_connect_idx = create_numpy_or_none_i(connect_idx, n_bnd+1)

  if (connect == NULL) :
      np_connect = None
  else :
      np_connect = create_numpy_g(connect, connect_idx[n_bnd])
  return (np_connect_idx, np_connect)

# ------------------------------------------------------------------------
def dmesh_bound_set(DMesh pydm,
                    PDM_bound_type_t bound_type,
                    NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                    NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None
                    ):
  """
  """
  # ************************************************************************
  # > Declaration
  cdef int           n_bnd
  # ************************************************************************

  n_bnd = connect_idx.shape[0]-1

  PDM_dmesh_bound_set(pydm._dm,
                      bound_type,
                      n_bnd,
      <PDM_g_num_t *> connect.data,
              <int *> connect_idx.data,
                      PDM_OWNERSHIP_USER)

# ------------------------------------------------------------------------
def dmesh_connectivity_set(DMesh pydm,
                           PDM_connectivity_type_t entity_type,
                           NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] connect_idx,
                           NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] connect    not None
                           ):
  """
  """
  # ************************************************************************
  # > Declaration
  # ************************************************************************

  cdef PDM_g_num_t* _connect     = np_to_gnum_pointer(connect)
  cdef int*         _connect_idx = np_to_int_pointer(connect_idx)
  
  PDM_dmesh_connectivity_set(pydm._dm,
                             entity_type,
                             _connect,
                             _connect_idx,
                             PDM_OWNERSHIP_USER)

# ------------------------------------------------------------------------
def dmesh_vtx_coord_set(DMesh pydm,
                        NPY.ndarray[NPY.double_t, mode='c', ndim=1] dvtx_coord not None):
    PDM_dmesh_vtx_coord_set(pydm._dm,
                 <double *> dvtx_coord.data,
                            PDM_OWNERSHIP_USER)

def dmesh_vtx_coord_get(DMesh pydm):
    dn_vtx = PDM_dmesh_dn_entity_get(pydm._dm, PDM_MESH_ENTITY_VTX)
    cdef double* dvtx_coord = NULL
    PDM_dmesh_vtx_coord_get(pydm._dm, 
                           &dvtx_coord,
                            PDM_OWNERSHIP_USER)
    return create_numpy_d(dvtx_coord, 3*dn_vtx)


def dfind_topological_ridge(MPI.Comm comm,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] distrib_face,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx,
                            NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dgroup_face_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dgroup_face):
  """
  dfind_topological_ridge(comm, distrib_face, dface_vtx_idx, dface_vtx, dgroup_face_idx, dgroup_face)

  Retrieve ridges from block-distributed faces with groups, and build associated edges

  Parameters:
    comm            (MPI.Comm)                   : MPI communicator
    distrib_face    (np.ndarray[npy_pdm_gnum_t]) : Distribution of faces
    dface_vtx_idx   (np.ndarray[np.int32_t])     : Index for block-distributed face->vertex connectivity
    dface_vtx       (np.ndarray[npy_pdm_gnum_t]) : Block-distributed face->vertex connectivity
    dgroup_face_idx (np.ndarray[np.int32_t])     : Index for block-distributed group->face connectivity
    dgroup_face     (np.ndarray[npy_pdm_gnum_t]) : Block-distributed group->face connectivity

  Returns:
    Distribution of edges (only those on ridges)                       (np.ndarray[npy_pdm_gnum_t])
    Block-distributed edge->vertex connectivity (only those on ridges) (np.ndarray[npy_pdm_gnum_t])
    Index for block-distributed group->edge connectivity               (np.ndarray[np.int32_t])
    Block-distributed group->edge connectivity                         (np.ndarray[npy_pdm_gnum_t])
    Index for ridge->surface connectivity                              (np.ndarray[np.int32_t])
    Ridge->surface connectivity                                        (np.ndarray[np.int32_t])
  """
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_g_num_t  *distrib_ridge         = NULL
  cdef PDM_g_num_t  *dridge_vtx            = NULL
  cdef int           n_group_ridge         = 0
  cdef int          *dgroup_edge_idx       = NULL
  cdef PDM_g_num_t  *dgroup_edge           = NULL
  cdef int          *dridge_face_group_idx = NULL
  cdef int          *dridge_face_group     = NULL

  cdef int n_group_face = dgroup_face_idx.shape[0]-1

  PDM_dmesh_find_topological_ridges(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                    <PDM_g_num_t *> distrib_face.data,
                    <int         *> dface_vtx_idx.data,
                    <PDM_g_num_t *> dface_vtx.data,
                                    n_group_face,
                    <int         *> dgroup_face_idx.data,
                    <PDM_g_num_t *> dgroup_face.data,
                                    &distrib_ridge,
                                    &dridge_vtx,
                                    &n_group_ridge,
                                    &dgroup_edge_idx,
                                    &dgroup_edge,
                                    &dridge_face_group_idx,
                                    &dridge_face_group);
  n_rank = comm.Get_size()
  i_rank = comm.Get_rank()
  cdef int dn_edge = distrib_ridge[i_rank+1] - distrib_ridge[i_rank]
  np_distrib_ridge   = create_numpy_g(distrib_ridge, n_rank+1)
  np_dridge_vtx      = create_numpy_g(dridge_vtx, 2 * dn_edge)
  np_dgroup_edge_idx = create_numpy_i(dgroup_edge_idx, n_group_ridge+1)
  np_dgroup_edge     = create_numpy_g(dgroup_edge, dgroup_edge_idx[n_group_ridge])

  np_dridge_face_group_idx = create_numpy_i(dridge_face_group_idx, dn_edge+1)
  np_dridge_face_group     = create_numpy_i(dridge_face_group    , dridge_face_group_idx[dn_edge])

  return np_distrib_ridge, np_dridge_vtx, np_dgroup_edge_idx, np_dgroup_edge, np_dridge_face_group_idx, np_dridge_face_group
