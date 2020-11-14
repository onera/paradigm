
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_t* PDM_dmesh_create(int          dn_cell,
                                  int          dn_face,
                                  int          dn_vtx,
                                  int          n_bnd,
                                  int          n_join);

    void PDM_dmesh_set(PDM_dmesh_t  *dm,
                       double       *dvtx_coord,
                       int          *dface_vtx_dmx,
                       PDM_g_num_t  *dface_vtx,
                       PDM_g_num_t  *dface_cell,
                       int          *dface_bound_dmx,
                       PDM_g_num_t  *dface_bound,
                       int          *join_g_dms,
                       int          *dface_join_dmx,
                       PDM_g_num_t  *dface_join);

    void PDM_dmesh_dims_get(PDM_dmesh_t *dm,
                            int         *dn_cell,
                            int         *dn_face,
                            int         *dn_vtx,
                            int         *n_bnd,
                            int         *n_joins);

    void PDM_dmesh_data_get(PDM_dmesh_t   *dm,
                            double       **dvtx_coord,
                            int          **dface_vtx_dmx,
                            PDM_g_num_t  **dface_vtx,
                            PDM_g_num_t  **dface_cell,
                            int          **dface_bound_dmx,
                            PDM_g_num_t  **dface_bound,
                            int          **join_g_dms,
                            int          **dface_join_dmx,
                            PDM_g_num_t  **dface_join);

    void PDM_dmesh_free(PDM_dmesh_t   *dm);

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DistributedMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_t* _dm
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, dn_cell,
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
    self._dm = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, n_bnd, n_join)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dvtx_coord   not None,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_vtx_dmx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_vtx    not None,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_cell,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_bound_dmx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_bound,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] join_g_dms,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dface_join_dmx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dface_join):
    """
    """

    PDM_dmesh_set(self._dm,
                  <double*>      dvtx_coord.data,
                  <int*>         dface_vtx_dmx.data,
                  <PDM_g_num_t*> dface_vtx.data,
                  <PDM_g_num_t*> dface_cell.data,
                  <int*>         dface_bound_dmx.data,
                  <PDM_g_num_t*> dface_bound.data,
                  <int*>         join_g_dms.data,
                  <int*>         dface_join_dmx.data,
                  <PDM_g_num_t*> dface_join.data)


  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_dmesh_free(self._dm)

