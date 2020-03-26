
cdef extern from "pdm_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    int PDM_dmesh_create(int          dNCell,
                         int          dNFace,
                         int          dNVtx,
                         int          dNBnd,
                         int          dNJoin);

    void PDM_dmesh_set(int           id,
                       double       *dVtxCoord,
                       int          *dFaceVtxIdx,
                       PDM_g_num_t  *dFaceVtx,
                       PDM_g_num_t  *dFaceCell,
                       int          *dFaceBoundIdx,
                       PDM_g_num_t  *dFaceBound,
                       int          *dJoinGIds,
                       int          *dFaceJoinIdx,
                       PDM_g_num_t  *dFaceJoin);

    void PDM_dmesh_dims_get(int   id,
                            int        *dNCell,
                            int        *dNFace,
                            int        *dNVtx,
                            int        *dNBnd,
                            int        *dNJoins);

    void PDM_dmesh_data_get(int            id,
                            double       **dVtxCoord,
                            int          **dFaceVtxIdx,
                            PDM_g_num_t  **dFaceVtx,
                            PDM_g_num_t  **dFaceCell,
                            int          **dFaceBoundIdx,
                            PDM_g_num_t  **dFaceBound,
                            int          **dJoinGIds,
                            int          **dFaceJoinIdx,
                            PDM_g_num_t  **dFaceJoin);

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
  def __init__(self, dNCell,
                     dNFace,
                     dNVtx,
                     dNBnd,
                     dNJoin):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._id = PDM_dmesh_create(dNCell, dNFace, dNVtx, dNBnd, dNJoin)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def dmesh_set(self, NPY.ndarray[NPY.double_t  , mode='c', ndim=1] dVtxCoord   not None,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceVtxIdx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceVtx    not None,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceCell,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceBoundIdx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceBound,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dJoinGIds,
                      NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] dFaceJoinIdx,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] dFaceJoin):
    """
    """

    PDM_dmesh_set(self._id,
                  <double*>      dVtxCoord.data,
                  <int*>         dFaceVtxIdx.data,
                  <PDM_g_num_t*> dFaceVtx.data,
                  <PDM_g_num_t*> dFaceCell.data,
                  <int*>         dFaceBoundIdx.data,
                  <PDM_g_num_t*> dFaceBound.data,
                  <int*>         dJoinGIds.data,
                  <int*>         dFaceJoinIdx.data,
                  <PDM_g_num_t*> dFaceJoin.data)


  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    print('PDM_dmesh_free')
    PDM_dmesh_free(self._id)

