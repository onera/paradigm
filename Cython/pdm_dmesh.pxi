
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    int PDM_dmesh_create(int          dn_cell,
                         int          dn_face,
                         int          dn_vtx,
                         int          n_bnd,
                         int          n_join);

    void PDM_dmesh_set(int           id,
                       double       *dvtx_coord,
                       int          *dface_vtx_idx,
                       PDM_g_num_t  *dface_vtx,
                       PDM_g_num_t  *dface_cell,
                       int          *dface_bound_idx,
                       PDM_g_num_t  *dface_bound,
                       int          *join_g_ids,
                       int          *dface_join_idx,
                       PDM_g_num_t  *dface_join);

    void PDM_dmesh_dims_get(int   id,
                            int        *dn_cell,
                            int        *dn_face,
                            int        *dn_vtx,
                            int        *n_bnd,
                            int        *n_joins);

    void PDM_dmesh_data_get(int            id,
                            double       **dvtx_coord,
                            int          **dface_vtx_idx,
                            PDM_g_num_t  **dface_vtx,
                            PDM_g_num_t  **dface_cell,
                            int          **dface_bound_idx,
                            PDM_g_num_t  **dface_bound,
                            int          **join_g_ids,
                            int          **dface_join_idx,
                            PDM_g_num_t  **dface_join);

    void PDM_dmesh_free(int id);

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistributedMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef public int _id
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __init__(self, dn_cell,
                     dn_face,
                     dn_vtx,
                     n_bnd,
                     n_join):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._id = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, n_bnd, n_join)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx    not None,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_cell,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_bound_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_bound,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] join_g_ids,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_join_idx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_join):
    """
    """

    PDM_dmesh_set(self._id,
                  <double*>      dvtx_coord.data,
                  <int*>         dface_vtx_idx.data,
                  <PDM_g_num_t*> dface_vtx.data,
                  <PDM_g_num_t*> dface_cell.data,
                  <int*>         dface_bound_idx.data,
                  <PDM_g_num_t*> dface_bound.data,
                  <int*>         join_g_ids.data,
                  <int*>         dface_join_idx.data,
                  <PDM_g_num_t*> dface_join.data)


  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_dmesh_free(self._id)

