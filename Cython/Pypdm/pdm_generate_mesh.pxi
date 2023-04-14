cdef extern from "pdm_generate_mesh.h":
  void PDM_generate_mesh_rectangle_simplified(PDM_MPI_Comm   comm,
                                              PDM_g_num_t    n_vtx_seg,
                                              int           *n_vtx,
                                              int           *n_elt,
                                              double       **coords,
                                              int          **elt_vtx_idx,
                                              int          **elt_vtx);

# ------------------------------------------------------------------

cdef extern from "pdm_mpi.h":
    PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm (void *mpi_comm)

# ------------------------------------------------------------------

def generate_mesh_rectangle_simplified(MPI.Comm      comm,
                                       PDM_g_num_t   n_vtx_seg):

  """
  Get a simple rectangular mesh.

  Parameters:
    n_vtx_seg     (int) : Number of vertices on the rectangle sides

  Returns:
    n_vtx       (int)                       : Number of vertices
    n_elt       (int)                       : Number of elements
    coords      (NPY.ndarray[NPY.double_t]) : Coordinates
    elt_vtx_idx (NPY.ndarray[NPY.int32_t])  : Index of element-vertex connectivity
    elt_vtx     (NPY.ndarray[NPY.int32_t])  : Element-vertex connectivity
  """

  # MPI communicator
  cdef MPI.MPI_Comm c_comm = comm.ob_mpi

  # Get
  cdef int n_vtx = -1
  cdef int n_elt = -1
  cdef double *coords
  cdef int    *elt_vtx_idx
  cdef int    *elt_vtx

  PDM_generate_mesh_rectangle_simplified(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                         n_vtx_seg,
                                         &n_vtx,
                                         &n_elt,
                                         &coords,
                                         &elt_vtx_idx,
                                         &elt_vtx)

  return {
          'n_vtx'       : n_vtx,
          'n_elt'       : n_elt,
          'coords'      : create_numpy_d(coords, 3*n_vtx),
          'elt_vtx_idx' : create_numpy_i(elt_vtx_idx, n_elt+1),
          'elt_vtx'     : create_numpy_i(elt_vtx, elt_vtx_idx[n_elt])
         }
