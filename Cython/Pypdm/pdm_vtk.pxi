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

  void PDM_vtk_write_std_elements_double(const char                 *filename,
                                         const int                   n_vtx,
                                         const double                vtx_coord[],
                                         const PDM_g_num_t           vtx_g_num[],
                                         const PDM_Mesh_nodal_elt_t  elt_type,
                                         const int                   n_elt,
                                         const int                   elt_vtx[],
                                         const PDM_g_num_t           elt_g_num[],
                                         const int                   n_elt_ifield,
                                         const char                 *elt_ifield_name[],
                                         const int                  *elt_ifield[])

# ------------------------------------------------------------------

def vtk_write_polydata(char *                                           filename,
                       NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord,
                       NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
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

    cdef int n_vtx  = len(vtx_coord) // 3
    cdef int n_face = len(face_vtx_idx) - 1

    cdef PDM_g_num_t *c_vtx_g_num  = np_to_gnum_pointer(vtx_g_num)
    cdef PDM_g_num_t *c_face_g_num = np_to_gnum_pointer(face_g_num)
    cdef int         *c_face_color = np_to_int_pointer(face_color)

    PDM_vtk_write_polydata(filename,
                           n_vtx,
                           <PDM_real_t *> vtx_coord.data,
                           c_vtx_g_num,
                           n_face,
                           <int *> face_vtx_idx.data,
                           <int *> face_vtx.data,
                           c_face_g_num,
                           c_face_color)

def vtk_write_point_cloud(char *                                           filename,
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

  cdef int n_vtx = len(vtx_coord) // 3

  cdef PDM_g_num_t *c_vtx_g_num = np_to_gnum_pointer(vtx_g_num)
  cdef int         *c_color     = np_to_int_pointer(color)

  PDM_vtk_write_point_cloud(filename,
                            n_vtx,
                            <PDM_real_t *> vtx_coord.data,
                            c_vtx_g_num,
                            c_color)


def vtk_write_std_elements_double(char                                            *filename,
                                  NPY.ndarray[NPY.double_t, mode='c', ndim=1]      vtx_coord not None,
                                  NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  vtx_g_num,
                                  PDM_Mesh_nodal_elt_t                             elt_type,
                                  NPY.ndarray[NPY.int32_t  , mode='c', ndim=1]     elt_vtx,
                                  NPY.ndarray[npy_pdm_gnum_t  , mode='c', ndim=1]  elt_g_num,
                                  dict                                             elt_fields):
  """
  vtk_write_std_elements_double(filename, vtx_coord, vtx_g_num, elt_type, elt_vtx, elt_g_num, elt_fields)
  Write an unstructured mesh made of standard elements of the same type.

  Parameters:
    filename   (str)                        : Output file name
    vtx_coord  (np.ndarray[np.double_t])    : Vertex coordinates
    vtx_g_num  (np.ndarray[npy_pdm_gnum_t]) : Vertex global IDs (or None)
    elt_type   (int)                        : Element type
    elt_vtx    (np.ndarray[np.int32_t])     : Element->vertex connectivity
    elt_g_num  (np.ndarray[npy_pdm_gnum_t]) : Element global IDs (or None)
    elt_fields (dict)                       : Element-based fields
  """
  cdef int stride
  if   elt_type == _PDM_MESH_NODAL_POINT:
    stride = 1
  elif elt_type == _PDM_MESH_NODAL_BAR2:
    stride = 2
  elif elt_type == _PDM_MESH_NODAL_TRIA3:
    stride = 3
  elif elt_type == _PDM_MESH_NODAL_QUAD4:
    stride = 4
  elif elt_type == _PDM_MESH_NODAL_TETRA4:
    stride = 4
  elif elt_type == _PDM_MESH_NODAL_PYRAMID5:
    stride = 5
  elif elt_type == _PDM_MESH_NODAL_PRISM6:
    stride = 6
  elif elt_type == _PDM_MESH_NODAL_HEXA8:
    stride = 8
  else:
    raise ValueError(f"Invalid elt_type {elt_type}")

  if len(elt_vtx)%stride != 0:
    raise ValueError(f"Size of elt_vtx must be a multiple of {stride} for elt_type {elt_type}")

  cdef int n_vtx = len(vtx_coord) // 3
  cdef int n_elt = len(elt_vtx)   // stride

  cdef PDM_g_num_t *c_vtx_g_num = np_to_gnum_pointer(vtx_g_num)
  cdef PDM_g_num_t *c_elt_g_num = np_to_gnum_pointer(elt_g_num)

  # cdef n_elt_field = len(elt_fields)

  PDM_vtk_write_std_elements_double(filename,
                                    n_vtx,
                         <double *> vtx_coord.data,
                                    c_vtx_g_num,
                                    elt_type,
                                    n_elt,
                            <int *> elt_vtx.data,
                                    c_elt_g_num,
                                    0,#n_elt_field,
                                    NULL,#elt_field_name,
                                    NULL)#elt_field)
