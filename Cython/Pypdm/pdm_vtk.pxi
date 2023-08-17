cdef extern from "pdm_vtk.h":

  void PDM_vtk_write_polydata(char*        filename,
                              int          n_vtx,
                              double*      vtx_coord,
                              PDM_g_num_t* vtx_g_num,
                              int          n_face,
                              int*         face_vtx_idx,
                              int*         face_vtx,
                              PDM_g_num_t* face_g_num,
                              int*         face_color)

  void PDM_vtk_write_point_cloud(char*        filename,
                                 int          n_vtx,
                                 double*      vtx_coord,
                                 PDM_g_num_t* vtx_g_num,
                                 int*         color)

# ------------------------------------------------------------------

def vtk_write_polydata(char *                                           filename,
                       int                                              n_vtx,
                       NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord,
                       NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
                       int                                              n_face,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_vtx_idx,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_vtx,
                       NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    face_g_num,
                       NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     face_color):
    """
    vtk_write_polydata(filename, n_vtx, vtx_coord, vtx_g_num, n_face, face_vtx_idx, face_vtx, face_g_num, face_color)
    Write polygons in a file of vtk format.

    Parameters:
      filename        (str)                           :
      n_vtx           (int)                           :
      vtx_coord       (Numpy array of NPY.double_t)   :
      vtx_g_num       (Numpy array of npy_pdm_gnum_t) :
      n_face          (int)                           :
      face_vtx_idx    (Numpy array of NPY.int32_t)    :
      face_vtx        (Numpy array of  NPY.int32_t)   :
      face_g_num      (Numpy array of npy_pdm_gnum_t) :
      face_color      (Numpy array of  NPY.int32_t)   :
    """

    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           <PDM_real_t *> vtx_coord.data,
                           <PDM_g_num_t*> vtx_g_num.data,
                           n_face,
                           <int *> face_vtx_idx.data,
                           <int *> face_vtx.data,
                           <PDM_g_num_t*> face_g_num.data,
                           <int *> face_color.data)

def vtk_write_point_cloud(char *                                           filename,
                          int                                              n_vtx,
                          NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord,
                          NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
                          NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     color):

  """
  vtk_write_point_cloud(filename, n_vtx, vtx_coord, vtx_g_num, color)
  Write point cloud in a file of vtk format.

  Parameters:
    filename        (str)                           :
    n_vtx           (int)                           :
    vtx_coord       (Numpy array of NPY.double_t)   :
    vtx_g_num       (Numpy array of npy_pdm_gnum_t) :
    color           (Numpy array of  NPY.int32_t)   :
  """

  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            <PDM_real_t *> vtx_coord.data,
                            <PDM_g_num_t*> vtx_g_num.data,
                            <int *> color.data)
