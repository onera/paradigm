cdef extern from "pdm_laplacian_smoothing.h":

  void PDM_laplacian_smoothing_idw_weights_compute(int       n_part,
                                                   double  **p_vtx_coord,
                                                   int      *p_n_edge,
                                                   int     **p_edge_vtx,
                                                   int       exponent,
                                                   double ***out_p_edge_weight);

  void PDM_laplacian_smoothing_cotangent_weights_compute(int                      n_part,
                                                         double                 **p_vtx_coord,
                                                         int                     *p_n_face,
                                                         int                    **p_face_edge,
                                                         int                     *p_n_edge,
                                                         int                    **p_edge_vtx,
                                                         PDM_part_comm_graph_t   *pcg_edge,
                                                         double                ***out_p_edge_weight);

  void PDM_laplacian_smoothing_fields(const PDM_MPI_Comm            comm,
                                            int                     n_part,
                                            int                    *p_n_vtx,
                                            PDM_part_comm_graph_t  *pcg_vtx,
                                            int                    *p_n_vtx_frozen,
                                            int                   **p_vtx_frozen,
                                            int                    *p_n_edge,
                                            int                   **p_edge_vtx,
                                            double                **p_edge_weight,
                                            PDM_part_comm_graph_t  *pcg_edge,
                                            double                  damping,
                                            int                     n_iter,
                                            double                  tol,
                                            int                     stride,
                                            double                **p_vtx_field);

  double PDM_laplacian_smoothing_fields_one_iteration(const PDM_MPI_Comm            comm,
                                                            int                     n_part,
                                                            int                    *p_n_vtx,
                                                            PDM_part_comm_graph_t  *pcg_vtx,
                                                            int                    *p_n_vtx_frozen,
                                                            int                   **p_vtx_frozen,
                                                            int                    *p_n_edge,
                                                            int                   **p_edge_vtx,
                                                            double                **p_edge_weight,
                                                            PDM_part_comm_graph_t  *pcg_edge,
                                                            double                  damping,
                                                            double                  tol,
                                                            int                     stride,
                                                            double                **p_vtx_field_prev,
                                                            double                **p_vtx_field_current);

cdef extern from "pdm_part_comm_graph.h":
  ctypedef struct PDM_part_comm_graph_t:
    pass

def compute_idw_weights(list p_vtx_coord,
                        list p_edge_vtx,
                        int  exponent):
  """
  Compute Inverse Distance Weighting-like edge weights.

  Parameters:
    p_vtx_coord   (list) : Vertex coordinates (size=n_part)
    p_edge_vtx    (list) : Edge→vertex connectivity (size=n_part)
    exponent      (int ) : Edge length exponent

  Returns:
    p_edge_weight (list) : Edge weight (size=n_part)
  """

  # Get sizes
  cdef int  n_part = len(p_vtx_coord)
  cdef int *n_edge = <int *> malloc(sizeof(int) * n_part)
  for i_part in range(n_part):
    n_edge[i_part] = len(p_edge_vtx[i_part])//2

  # Convert
  cdef double **c_p_vtx_coord = np_list_to_double_pointers(p_vtx_coord)
  cdef int    **c_p_edge_vtx  = np_list_to_int_pointers(p_edge_vtx)

  cdef double **c_p_edge_weight
  PDM_laplacian_smoothing_idw_weights_compute(n_part,
                                              c_p_vtx_coord,
                                              n_edge,
                                              c_p_edge_vtx,
                                              exponent,
                                              &c_p_edge_weight)

  p_edge_weight = list()
  for i_part in range(n_part):
    p_edge_weight.append(create_numpy_or_none_d(c_p_edge_weight[i_part], n_edge[i_part]))

  # Free
  free(n_edge         )
  free(c_p_edge_vtx   )
  free(c_p_vtx_coord  )
  free(c_p_edge_weight)

  return p_edge_weight

def compute_cotangent_weights(list p_vtx_coord,
                              list p_face_edge,
                              list p_edge_vtx,
                              PyPartCommGraph pypcg_edge):
  """
  Compute cotangent edge weights for triangular meshes.

  Parameters:
    p_vtx_coord   (list           ) : Vertex coordinates (size=n_part)
    p_face_edge   (list           ) : Face→edge connectivity (size=n_part)
    p_edge_vtx    (list           ) : Edge→vertex connectivity (size=n_part)
    pypcg_edge    (PyPartCommGraph) : Edge part comm graph or None

  Returns:
    p_edge_weight (list           ) : Edge weight (size=n_part)
  """

  # Get sizes
  cdef int  n_part = len(p_vtx_coord)
  cdef int *n_face = <int *> malloc(sizeof(int) * n_part)
  cdef int *n_edge = <int *> malloc(sizeof(int) * n_part)
  for i_part in range(n_part):
    n_face[i_part] = len(p_face_edge[i_part])//3
    n_edge[i_part] = len(p_edge_vtx [i_part])//2

  # Convert
  cdef double **c_p_vtx_coord = np_list_to_double_pointers(p_vtx_coord)
  cdef int    **c_p_face_edge = np_list_to_int_pointers(p_face_edge)
  cdef int    **c_p_edge_vtx  = np_list_to_int_pointers(p_edge_vtx)

  cdef double **c_p_edge_weight
  PDM_laplacian_smoothing_cotangent_weights_compute(n_part,
                                                    c_p_vtx_coord,
                                                    n_face,
                                                    c_p_face_edge,
                                                    n_edge,
                                                    c_p_edge_vtx,
                                                    pypcg_edge.pcg,
                                                    &c_p_edge_weight)

  p_edge_weight = list()
  for i_part in range(n_part):
    p_edge_weight.append(create_numpy_or_none_d(c_p_edge_weight[i_part], n_edge[i_part]))

  # Free
  free(n_face         )
  free(n_edge         )
  free(c_p_face_edge  )
  free(c_p_edge_vtx   )
  free(c_p_vtx_coord  )
  free(c_p_edge_weight)

  return p_edge_weight

def laplacian_smoothing_fields(MPI.Comm        comm,
                               list            p_vtx_frozen,
                               PyPartCommGraph pypcg_vtx,
                               list            p_edge_vtx,
                               list            p_edge_weight,
                               PyPartCommGraph pypcg_edge,
                               float           damping,
                               int             n_iter,
                               float           tol,
                               int             stride,
                               list            p_vtx_field):
  """
  Apply laplacian smoothing to strided fields (interlaced).

  Parameters:
    comm          (MPI.Comm       ) : MPI communicator
    p_vtx_frozen  (list           ) : Frozen vertex list (size=n_part)
    pypcg_vtx     (PyPartCommGraph) : Vertex part comm graph
    p_edge_vtx    (list           ) : Edge→vertex connectivity (size=n_part)
    p_edge_weight (list           ) : Edge weights (size=n_part) or None
    pypcg_edge    (PyPartCommGraph) : Edge part comm graph or None
    damping       (float          ) : Damping constant (between 0. and 1.)
    n_iter        (int            ) : Number of smoothing iterations
    tol           (float          ) : Relative tolerance for convergence (ignored if negative)
    stride        (int            ) : Field stride (interlaced values)
    p_vtx_field   (list           ) : Fields (size=n_part)
  """

  # Convert mpi4py -> PDM_MPI
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm PDM_comm = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  # Get sizes
  cdef int  n_part = len(p_vtx_field)
  cdef int *n_vtx  = <int *> malloc(sizeof(int) * n_part)
  cdef int *n_edge = <int *> malloc(sizeof(int) * n_part)
  for i_part in range(n_part):
    n_vtx [i_part] = len(p_vtx_field[i_part])//stride
    n_edge[i_part] = len(p_edge_vtx[i_part])//2

  # Convert
  cdef int    **c_p_edge_vtx  = np_list_to_int_pointers(p_edge_vtx)
  cdef double **c_p_vtx_field = np_list_to_double_pointers(p_vtx_field)

  cdef int  *n_vtx_frozen   = NULL
  cdef int **c_p_vtx_frozen = NULL
  if p_vtx_frozen is not None:
    n_vtx_frozen = <int *> malloc(sizeof(int) * n_part)
    for i_part in range(n_part):
      n_vtx_frozen[i_part] = len(p_vtx_frozen[i_part])
    c_p_vtx_frozen = np_list_to_int_pointers(p_vtx_frozen)

  cdef double **c_p_edge_weight = NULL
  if p_edge_weight is not None:
    c_p_edge_weight = np_list_to_double_pointers(p_edge_weight)

  cdef PDM_part_comm_graph_t *pcg_edge = NULL
  if pypcg_edge is not None:
    pcg_edge = pypcg_edge.pcg

  PDM_laplacian_smoothing_fields(PDM_comm,
                                 n_part,
                                 n_vtx,
                                 pypcg_vtx.pcg,
                                 n_vtx_frozen,
                                 c_p_vtx_frozen,
                                 n_edge,
                                 c_p_edge_vtx,
                                 c_p_edge_weight,
                                 pcg_edge,
                                 damping,
                                 n_iter,
                                 tol,
                                 stride,
                                 c_p_vtx_field)

  # Free
  free(n_vtx          )
  free(n_edge         )
  free(n_vtx_frozen   )
  free(c_p_edge_vtx   )
  free(c_p_edge_weight)
  free(c_p_vtx_frozen )
  free(c_p_vtx_field  )

def laplacian_smoothing_fields_one_iteration(MPI.Comm        comm,
                                             list            p_vtx_frozen,
                                             PyPartCommGraph pypcg_vtx,
                                             list            p_edge_vtx,
                                             list            p_edge_weight,
                                             PyPartCommGraph pypcg_edge,
                                             float           damping,
                                             int             n_iter,
                                             float           tol,
                                             int             stride,
                                             list            p_vtx_field_prev,
                                             list            p_vtx_field_current):
  """
  Apply one iteration of Laplacian smoothing to strided fields (interlaced).

  Parameters:
    comm                (MPI.Comm       ) : MPI communicator
    p_vtx_frozen        (list           ) : Frozen vertex list (size=n_part)
    pypcg_vtx           (PyPartCommGraph) : Vertex part comm graph
    p_edge_vtx          (list           ) : Edge→vertex connectivity (size=n_part)
    p_edge_weight       (list           ) : Edge weights (size=n_part) or None
    pypcg_edge          (PyPartCommGraph) : Edge part comm graph or None
    damping             (float          ) : Damping constant (between 0. and 1.)
    tol                 (float          ) : Relative tolerance for convergence (ignored if negative)
    stride              (int            ) : Field stride (interlaced values)
    p_vtx_field_prev    (list           ) : Previous fields (size=n_part)
    p_vtx_field_current (list           ) : Current fields (size=n_part)

  Returns:
    eps                 (float          ) : Maximal absolute relative variation of fields (if tol is positive)
  """

  # Convert mpi4py -> PDM_MPI
  cdef MPI.MPI_Comm c_comm   = comm.ob_mpi
  cdef PDM_MPI_Comm PDM_comm = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)

  # Get sizes
  cdef int  n_part = len(p_vtx_field_prev)
  cdef int *n_vtx  = <int *> malloc(sizeof(int) * n_part)
  cdef int *n_edge = <int *> malloc(sizeof(int) * n_part)
  for i_part in range(n_part):
    n_vtx [i_part] = len(p_vtx_field_prev[i_part])//stride
    n_edge[i_part] = len(p_edge_vtx[i_part])//2

  # Convert
  cdef int    **c_p_edge_vtx          = np_list_to_int_pointers(p_edge_vtx)
  cdef double **c_p_vtx_field_prev    = np_list_to_double_pointers(p_vtx_field_prev)
  cdef double **c_p_vtx_field_current = np_list_to_double_pointers(p_vtx_field_current)
  cdef double **c_p_edge_weight       = np_list_to_double_pointers(p_edge_weight)

  cdef int  *n_vtx_frozen   = NULL
  cdef int **c_p_vtx_frozen = NULL
  if p_vtx_frozen is not None:
    n_vtx_frozen = <int *> malloc(sizeof(int) * n_part)
    for i_part in range(n_part):
      n_vtx_frozen[i_part] = len(p_vtx_frozen[i_part])
    c_p_vtx_frozen = np_list_to_int_pointers(p_vtx_frozen)

  cdef double eps
  eps = PDM_laplacian_smoothing_fields_one_iteration(PDM_comm,
                                                     n_part,
                                                     n_vtx,
                                                     pypcg_vtx.pcg,
                                                     n_vtx_frozen,
                                                     c_p_vtx_frozen,
                                                     n_edge,
                                                     c_p_edge_vtx,
                                                     c_p_edge_weight,
                                                     pypcg_edge.pcg,
                                                     damping,
                                                     tol,
                                                     stride,
                                                     c_p_vtx_field_prev,
                                                     c_p_vtx_field_current)

  # Free
  free(n_vtx          )
  free(n_edge         )
  free(n_vtx_frozen   )
  free(c_p_edge_vtx   )
  free(c_p_edge_weight)
  free(c_p_vtx_frozen )
  free(c_p_vtx_field_prev   )
  free(c_p_vtx_field_current)

  return eps
