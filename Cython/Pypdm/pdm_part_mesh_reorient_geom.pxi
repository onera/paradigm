cdef extern from "pdm_part_mesh_reorient_geom.h":

  int PDM_part_mesh_reorient_geom(int          mesh_dim,
                                  double      *mesh_inv,
                                  int          n_part,
                                  int         *n_cell,
                                  int         *n_face,
                                  int         *n_face_group,
                                  int        **cell_face_idx,
                                  int        **cell_face,
                                  int        **face_cell,
                                  int        **face_vtx_idx,
                                  int        **face_vtx,
                                  double     **vtx_coord,
                                  int        **face_group_idx,
                                  int        **face_group,
                                  int          comm)

def part_mesh_reorient_geom(int                                         mesh_dim,
                            NPY.ndarray[NPY.double_t, mode='c', ndim=1] mesh_inv,
                            list                                        n_cell not None,
                            list                                        n_face not None,
                            list                                        n_face_group,
                            list                                        cell_face_idx,
                            list                                        cell_face,
                            list                                        face_cell,
                            list                                        face_vtx_idx,
                            list                                        face_vtx not None,
                            list                                        vtx_coord not None,
                            list                                        face_group_idx,
                            list                                        face_group,
                            MPI.Comm                                    comm):

  """
  Fix mesh orientation geometrically

  Parameters:
    mesh_dim       (int)                  Mesh dimension (2 or 3)
    mesh_inv       (numpy.double_t array) Mesh invariance direction (2D)
    n_cell         (list)                 Number of cells
    n_face         (list)                 Number of faces
    n_face_group   (list)                 Number of face groups or NULL
    cell_face_idx  (list)                 Cell face connectivity index or NULL
    cell_face      (list)                 Cell face connectivity or NULL
    face_cell      (list)                 Face cell connectivity or NULL
    face_vtx_idx   (list)                 Face to vertex connectivity index
    face_vtx       (list)                 Face to vertex connectivity
    vtx_coord      (list)                 Vertex coordinates
    face_group_idx (list)                 Index of faces list of each group or NULL
    face_group     (list)                 Faces list of each group or NULL
    comm           (MPI.Comm)             MPI communicator

  Returns:
    Reorientation indicator
  """

  cdef int n_part = len(n_cell)

  # Convert mpi4py -> PDM_MPI
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm PDM_comm = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  cdef int* _n_cell       = list_to_int_pointer(n_cell)
  cdef int* _n_face       = list_to_int_pointer(n_face)
  cdef int* _n_face_group = list_to_int_pointer(n_face_group) if n_face_group is not None else NULL

  cdef int** _cell_face_idx  = np_list_to_int_pointers(cell_face_idx) if cell_face_idx is not None else NULL
  cdef int** _cell_face      = np_list_to_int_pointers(cell_face) if cell_face is not None else NULL
  cdef int** _face_cell      = np_list_to_int_pointers(face_cell) if face_cell is not None else NULL
  cdef int** _face_vtx_idx   = np_list_to_int_pointers(face_vtx_idx) if face_vtx_idx is not None else NULL
  cdef int** _face_vtx       = np_list_to_int_pointers(face_vtx)
  cdef int** _face_group_idx = np_list_to_int_pointers(face_group_idx) if face_group_idx is not None else NULL
  cdef int** _face_group     = np_list_to_int_pointers(face_group) if face_group is not None else NULL

  cdef double** _vtx_coord = np_list_to_double_pointers(vtx_coord)

  cdef int fixed = PDM_part_mesh_reorient_geom(mesh_dim,
                                    <double *> mesh_inv.data,
                                               n_part,
                                               _n_cell,
                                               _n_face,
                                               _n_face_group,
                                               _cell_face_idx,
                                               _cell_face,
                                               _face_cell,
                                               _face_vtx_idx,
                                               _face_vtx,
                                               _vtx_coord,
                                               _face_group_idx,
                                               _face_group,
                                               PDM_comm)

  free(_n_cell      )
  free(_n_face      )
  free(_n_face_group)

  free(_cell_face_idx )
  free(_cell_face     )
  free(_face_cell     )
  free(_face_vtx_idx  )
  free(_face_vtx      )
  free(_face_group_idx)
  free(_face_group    )

  free(_vtx_coord)

  return fixed
