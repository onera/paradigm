
cdef extern from "pdm_reader_gamma.h":
    void PDM_write_meshb(char           *filename,
                         int            *n_elt_table,
                         int           **tag_table,
                         PDM_g_num_t   **vtx_connect_table,
                         double          *vtx_coords
    )

    PDM_dmesh_nodal_t *PDM_reader_gamma_dmesh_nodal(PDM_MPI_Comm   comm,
                                                    char    *filename,
                                                    int            fix_orientation_2d,
                                                    int            fix_orientation_3d)

    void PDM_write_gamma_sol(char   *filename,
                             int     n_vtx,
                             int     n_field,
                             double    *fields)
    void PDM_read_gamma_sol(char   *filename,
                             int     n_vtx,
                             int     n_field,
                             double    *fields)

    void PDM_write_gamma_matsym(char   *filename,
                                int     n_vtx,
                                double *fields)


    int PDM_read_gamma_sol_at_vertices(char     *filename,
                                       int      *n_field,
                                       int     **field_stride,
                                       double ***field_values)

def write_meshb(char                                        *filename,
                NPY.ndarray[NPY.int32_t, mode ='c', ndim=1]  n_elt_table,
                list                                         tag_table,
                list                                         vtx_connect_table,
                NPY.ndarray[NPY.double_t, mode ='c', ndim=1] vtx_coords):
    """
    """
    # ---- Variable initialization
    # TODO : Find a way to use _PDM_MESH_NODAL_N_ELEMENT_TYPES at compile time
    #        instead of hardcoding it
    DEF size = 19

    cdef int*                                          tag_buffer[size]
    cdef NPY.ndarray[NPY.int32_t, ndim=1, mode='c']    tag_temp

    cdef PDM_g_num_t*                                  vtx_connect_buffer[size]
    cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='c'] vtx_connect_temp

    # ---- Conversion to raw pointers
    for idx in range(_PDM_MESH_NODAL_N_ELEMENT_TYPES):
        tag_temp         = tag_table[idx]
        vtx_connect_temp = vtx_connect_table[idx]

        tag_buffer[idx]         = <int *> tag_temp.data
        vtx_connect_buffer[idx] = <PDM_g_num_t *> vtx_connect_temp.data

    PDM_write_meshb(filename,
                    <int *> n_elt_table.data,
                    tag_buffer,
                    vtx_connect_buffer,
                    <double *> vtx_coords.data)


def write_solb(char *filename,
               int n_vtx,
               int n_field,
               NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_write_gamma_sol(filename,
                      n_vtx,
                      n_field,
           <double *> field.data)


def write_matsym_solb(char *filename,
                      int n_vtx,
                      NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_write_gamma_matsym(filename,
                         n_vtx,
              <double *> field.data)


def read_solb(char *filename,
              int n_vtx,
              int n_field,
              NPY.ndarray[NPY.double_t, mode='c', ndim=1] field):
  """
  """
  PDM_read_gamma_sol(filename,
                 n_vtx,
                 n_field,
      <double *> field.data)


def meshb_to_dmesh_nodal(char *filename, MPI.Comm comm, int fix_orientation_2d, int fix_orientation_3d):
  """
  """
  cdef PDM_dmesh_nodal_t* dmn
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

  dmn = PDM_reader_gamma_dmesh_nodal(pdm_comm, filename, fix_orientation_2d, fix_orientation_3d)

  py_caps = PyCapsule_New(dmn, NULL, NULL);

  return DistributedMeshNodalCapsule(py_caps) # The free is inside the class


def read_sol_at_vertices(char *filename):
  """
  read_sol_at_vertices(filename)

  Parameters:
    filename (str) : Solution file name

  Returns:
    Dictionary
      - ``"n_vtx"``        (int)                               : Number of vertices
      - ``"field_stride"`` (`np.ndarray[np.int32_t]`)          : Field strides
      - ``"field_values"`` (`list of np.ndarray[np.double_t]`) : Field values
  """
  cdef int      n_field
  cdef int     *field_stride
  cdef double **field_values
  cdef int      n_vtx

  n_vtx = PDM_read_gamma_sol_at_vertices(filename,
                                         &n_field,
                                         &field_stride,
                                         &field_values)

  return {
    "n_vtx"        : n_vtx,
    "field_stride" : create_numpy_i(field_stride, n_vtx),
    "field_values" : [create_numpy_d(field_values[i], field_stride[i]*n_vtx) for i in range(n_field)]
  }
