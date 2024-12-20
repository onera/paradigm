
cdef extern from "pdm_partitioning_algorithm.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_part_distgroup_to_partgroup(PDM_MPI_Comm      comm,
                                         PDM_g_num_t      *entity_distribution,
                                         int               n_group,
                                         int              *dgroup_idx,
                                         PDM_g_num_t      *dgroup,
                                         int               n_part,
                                         int              *pn_entity,
                                         PDM_g_num_t     **pentity_ln_to_gn,
                                         int            ***pgroup_idx,
                                         int            ***pgroup,
                                         PDM_g_num_t    ***pgroup_ln_to_gn)
    void PDM_part_dcoordinates_to_pcoordinates(PDM_MPI_Comm    comm,
                                               int             n_part,
                                               PDM_g_num_t    *vertex_distribution,
                                               double         *dvtx_coord,
                                               int            *pn_vtx,
                                               PDM_g_num_t   **pvtx_ln_to_gn,
                                               double       ***pvtx_coord)

    void PDM_compute_dface_vtx_from_edges(PDM_MPI_Comm   comm,
                                          int            dn_face,
                                          int            dn_edge,
                                          int           *dface_edge_idx,
                                          PDM_g_num_t   *dface_edge,
                                          PDM_g_num_t   *dedge_vtx,
                                          PDM_g_num_t  **dface_vtx);

    void PDM_compute_face_edge_from_face_vtx(PDM_MPI_Comm    comm,
                                             int             n_part,
                                             int            *pn_face,
                                             int            *pn_vtx,
                                             int           **pface_vtx_idx,
                                             int           **pface_vtx,
                                             PDM_g_num_t   **pface_ln_to_gn,
                                             PDM_g_num_t   **pvtx_ln_to_gn,
                                             int          ***pface_edge_idx,
                                             int          ***pface_edge,
                                             int           **pn_edge,
                                             int          ***pedge_vtx,
                                             PDM_g_num_t  ***pedge_ln_to_gn);

    void PDM_compute_graph_comm_entity_ownership(int             n_part,
                                                 int            *n_entity,
                                                 PDM_g_num_t   **entity_ln_to_gn,
                                                 int           **n_owned_entity,
                                                 int          ***lnum_owned_entity,
                                                 PDM_MPI_Comm    comm)

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ===================================================================================
def part_distgroup_to_partgroup(MPI.Comm                                      comm,
                                NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] entity_distribution,
                                int                                           n_group,
                                NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dgroup_idx,
                                NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dgroup,
                                int                                           n_part,
                                list                                          pn_entity,
                                list                                          pentity_ln_to_gn):

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef PDM_g_num_t * entity_distribution_data = np_to_gnum_pointer(entity_distribution)

    cdef int _n_group = n_group
    cdef int _n_part = n_part

    cdef int *dgroup_idx_data = np_to_int_pointer(dgroup_idx)
    cdef PDM_g_num_t * dgroup_data = np_to_gnum_pointer(dgroup)

    cdef int * _pn_entity = list_to_int_pointer(pn_entity)
    cdef PDM_g_num_t** LNToGN = np_list_to_gnum_pointers(pentity_ln_to_gn)


    cdef int **_pgroup_idx
    cdef int **_pgroup
    cdef PDM_g_num_t **_pgroup_ln_to_gn

    PDM_part_distgroup_to_partgroup(
                                    PDMC,
                                    entity_distribution_data,
                                    _n_group,
                                    dgroup_idx_data,
                                    dgroup_data,
                                    _n_part,
                                    _pn_entity,
             <const PDM_g_num_t **> LNToGN,
                                    &_pgroup_idx,
                                    &_pgroup,
                                    &_pgroup_ln_to_gn
                                    )
    list_group_part = list()

    for i_part in range(_n_part):
        if (_pgroup_idx == NULL) :
            npgroup_idx = None
        else :
            npgroup_idx = create_numpy_i(_pgroup_idx[i_part], _n_group+1)
        if (_pgroup == NULL) :
            npgroup = None
        else :
            npgroup = create_numpy_i(_pgroup[i_part], _pgroup_idx[i_part][_n_group])

        if (_pgroup_ln_to_gn == NULL) :
            npgroup_ln_to_gn = None
        else :
            npgroup_ln_to_gn = create_numpy_g(_pgroup_ln_to_gn[i_part], _pgroup_idx[i_part][_n_group])

        dict_group_part = {'npZSRGroupIdx'    : npgroup_idx,
                           'npZSRGroup'       : npgroup,
                           'npZSRGroupLNToGN' : npgroup_ln_to_gn}
        list_group_part.append(dict_group_part)

    free(_pn_entity)
    free(_pgroup_idx)
    free(_pgroup)
    free(_pgroup_ln_to_gn)
    free(LNToGN)

    return list_group_part

cdef extern from "pdm_closest_points.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    void PDM_transform_to_parent_gnum(int           n_part_initial,
                                      int          *n_elmt_initial,
                                      PDM_g_num_t **child_ln_to_gn,
                                      PDM_g_num_t **parent_ln_to_gn,
                                      int           n_part_to_transform,
                                      int          *n_elmt_to_transform,
                                      PDM_g_num_t **gnum_to_transform,
                                      PDM_MPI_Comm  comm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ===================================================================================
#Why here ? Could be moved in closest_point
def transform_to_parent_gnum(list     gnum_to_transform,
                             list     child_ln_to_gn,
                             list     parent_ln_to_gn,
                             MPI.Comm comm):
  """
  """
  assert len(parent_ln_to_gn) == len(child_ln_to_gn)
  for child_gn, parent_gn in zip(child_ln_to_gn, parent_ln_to_gn):
    assert child_gn.shape == parent_gn.shape

  cdef int n_part_out = len(gnum_to_transform)
  cdef int n_part_ini   = len(child_ln_to_gn)

  cdef PDM_g_num_t **_gnum_to_transform = np_list_to_gnum_pointers(gnum_to_transform)
  cdef int          *_pn_elt_out        = list_to_int_pointer([t.shape[0] for t in gnum_to_transform])

  cdef PDM_g_num_t **_child_ln_to_gn  = np_list_to_gnum_pointers(child_ln_to_gn)
  cdef PDM_g_num_t **_parent_ln_to_gn = np_list_to_gnum_pointers(parent_ln_to_gn)
  cdef int          *_pn_elt_in       = list_to_int_pointer([t.shape[0] for t in parent_ln_to_gn])

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)
  PDM_transform_to_parent_gnum(                 n_part_ini,
                               <const int          *> _pn_elt_in,
                               <const PDM_g_num_t **> _child_ln_to_gn,
                               <const PDM_g_num_t **> _parent_ln_to_gn,
                                                n_part_out,
                               <const int          *> _pn_elt_out,
                               <PDM_g_num_t **> _gnum_to_transform,
                                                PDMC)

  free(_pn_elt_out)
  free(_pn_elt_in)
  free(_child_ln_to_gn)
  free(_parent_ln_to_gn)
  free(_gnum_to_transform)
# ===================================================================================


# ===================================================================================
def part_dcoordinates_to_pcoordinates(MPI.Comm                                      comm,
                                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_distribution,
                                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                                      list                                          l_pvtx_ln_to_gn):
    """
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef int n_part = len(l_pvtx_ln_to_gn)
    cdef int          *pn_vtx        = list_to_int_pointer([t.shape[0] for t in l_pvtx_ln_to_gn])
    cdef PDM_g_num_t **pvtx_ln_to_gn = np_list_to_gnum_pointers(l_pvtx_ln_to_gn)
    cdef double      **pvtx_coord


    PDM_part_dcoordinates_to_pcoordinates(PDMC,
                                          n_part,
                     <const PDM_g_num_t*> vtx_distribution.data,
                     <const double*>      dvtx_coord.data,
                     <const int *>        pn_vtx,
                     <const PDM_g_num_t**>pvtx_ln_to_gn,
                                          &pvtx_coord);

    l_pvtx_coord = list()
    for i in range(n_part):
        l_pvtx_coord.append(create_numpy_d(pvtx_coord[i], 3*pn_vtx[i]))

    free(pvtx_coord)
    free(pn_vtx)
    free(pvtx_ln_to_gn)

    return l_pvtx_coord


# ===================================================================================
def compute_dface_vtx_from_edges(MPI.Comm                                      comm,
                                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_edge_idx,
                                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_edge,
                                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dedge_vtx):
    """
    """
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef int dn_face = dface_edge_idx.shape[0]-1
    cdef int dn_edge = dedge_vtx.shape[0]//2
    cdef PDM_g_num_t *dface_vtx

    PDM_compute_dface_vtx_from_edges(PDMC,
                                     dn_face,
                                     dn_edge,
                    <int *>          dface_edge_idx.data,
                    <PDM_g_num_t *>  dface_edge.data,
                    <PDM_g_num_t *>  dedge_vtx.data,
                                    &dface_vtx)

    return create_numpy_g(dface_vtx, dface_edge.shape[0])


def compute_face_edge_from_face_vtx(MPI.Comm comm,
                                    list pn_face,
                                    list pn_vtx,
                                    list pface_vtx_idx,
                                    list pface_vtx,
                                    list pface_ln_to_gn,
                                    list pvtx_ln_to_gn):

  cdef int n_part = len(pn_face)
  assert len(pn_face) == len(pn_vtx) == len(pface_vtx_idx) == \
         len(pface_vtx) == len(pface_ln_to_gn) == len(pvtx_ln_to_gn)

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

  cdef int* _pn_face = list_to_int_pointer(pn_face)
  cdef int* _pn_vtx  = list_to_int_pointer(pn_vtx)

  cdef int** _pface_vtx_idx = np_list_to_int_pointers(pface_vtx_idx)
  cdef int** _pface_vtx     = np_list_to_int_pointers(pface_vtx)

  cdef PDM_g_num_t** _pface_ln_to_gn = np_list_to_gnum_pointers(pface_ln_to_gn)
  cdef PDM_g_num_t** _pvtx_ln_to_gn = np_list_to_gnum_pointers(pvtx_ln_to_gn)

  cdef int**         _pface_edge_idx = NULL
  cdef int**         _pface_edge     = NULL
  cdef int*          _pn_edge        = NULL
  cdef int**         _pedge_vtx      = NULL
  cdef PDM_g_num_t** _pedge_ln_to_gn = NULL

  PDM_compute_face_edge_from_face_vtx(PDMC,
                                      n_part,
                                      _pn_face,
                                      _pn_vtx,
                                      _pface_vtx_idx,
                                      _pface_vtx,
                                      _pface_ln_to_gn,
                                      _pvtx_ln_to_gn,
                                     &_pface_edge_idx,
                                     &_pface_edge,
                                     &_pn_edge,
                                     &_pedge_vtx,
                                     &_pedge_ln_to_gn)


  all_part_data = []
  for i_part in range(n_part):
    _n_face = _pn_face[i_part]
    _n_edge = _pn_edge[i_part]
    _s_face_edge = _pface_edge_idx[i_part][_n_face]

    part_data = {'np_face_edge_idx' : create_numpy_i(_pface_edge_idx[i_part], _n_face+1),
                 'np_face_edge'     : create_numpy_i(_pface_edge    [i_part], _s_face_edge),
                 'np_edge_vtx'      : create_numpy_i(_pedge_vtx     [i_part], 2*_n_edge),
                 'np_edge_ln_to_gn' : create_numpy_g(_pedge_ln_to_gn[i_part], _n_edge)}
    all_part_data.append(part_data)

  free(_pn_face)
  free(_pn_vtx)
  free(_pface_vtx_idx)
  free(_pface_vtx)
  free(_pface_ln_to_gn)
  free(_pvtx_ln_to_gn)

  free(_pface_edge_idx)
  free(_pface_edge)
  free(_pn_edge)
  free(_pedge_vtx)
  free(_pedge_ln_to_gn)

  return all_part_data

# ===================================================================================
def compute_graph_comm_entity_ownership(int                                           n_part,
                                        NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] n_entity not None,
                                        list                                          entity_ln_to_gn not None,
                                        MPI.Comm                                      comm):
  """
  compute_graph_comm_entity_ownership(n_part, n_entity, entity_ln_to_gn, comm)
  """

  # Convert mpi4py -> PDM_MPI
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm PDM_comm = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  cdef int  *_n_owned_entity    = NULL
  cdef int **_lnum_owned_entity = NULL

  cdef PDM_g_num_t **_entity_ln_to_gn = np_list_to_gnum_pointers(entity_ln_to_gn)

  PDM_compute_graph_comm_entity_ownership(n_part,
                                  <int *> n_entity.data,
                         <PDM_g_num_t **> _entity_ln_to_gn,
                                          &_n_owned_entity,
                                          &_lnum_owned_entity,
                                          PDM_comm)

  n_owned_entity = create_numpy_i(_n_owned_entity, n_part)

  lnum_owned_entity = []
  for i_part in range(n_part):
    lnum_owned_entity.append(create_numpy_i(_lnum_owned_entity[i_part], n_owned_entity[i_part]))

  free(_entity_ln_to_gn)
  free(_lnum_owned_entity)

  return n_owned_entity, lnum_owned_entity
