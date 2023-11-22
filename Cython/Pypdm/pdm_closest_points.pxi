
cdef extern from "pdm_closest_points.h":
  # > Wrapping of Ppart Structures
  ctypedef struct PDM_closest_point_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  PDM_closest_point_t* PDM_closest_points_create(PDM_MPI_Comm    comm,
                                int             n_closest,
                                PDM_ownership_t owner)

  void PDM_closest_points_n_part_cloud_set(PDM_closest_point_t* cls,
                                           int                  n_part_cloud_src,
                                           int                  n_part_cloud_tgt)

  void PDM_closest_points_tgt_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum)

  void PDM_closest_points_src_cloud_set(PDM_closest_point_t *cls,
                                        int                  i_part,
                                        int                  n_points,
                                        double              *coords,
                                        PDM_g_num_t         *gnum)

  void PDM_closest_points_compute(PDM_closest_point_t  *cls)

  void PDM_closest_points_get(PDM_closest_point_t  *cls,
                              int                   i_part_tgt,
                              PDM_g_num_t         **closest_src_gnum,
                              double              **closest_src_distance)

  void PDM_closest_points_tgt_in_src_get(PDM_closest_point_t  *cls,
                                         int                   i_part_src,
                                         int                 **tgt_in_src_idx,
                                         PDM_g_num_t         **tgt_in_src)

  void PDM_closest_points_tgt_in_src_dist_get(PDM_closest_point_t  *cls,
                                              int                   i_part_src,
                                              int                 **tgt_in_src_idx,
                                              double              **tgt_in_src_dist)

  void PDM_closest_points_free(PDM_closest_point_t *cls)

  void PDM_closest_points_dump_times(PDM_closest_point_t *cls)

  void PDM_closest_points_part_to_part_get(PDM_closest_point_t  *cls,
                                           PDM_part_to_part_t  **ptp,
                                           PDM_ownership_t       ownership);


# ------------------------------------------------------------------
cdef class ClosestPoints:
  """
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_closest_point_t *_cls
  cdef int _size
  cdef int _rank
  cdef int* src_n_points
  cdef int* tgt_n_points
  cdef int n_closest
  cdef MPI.Comm py_comm
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, MPI.Comm comm,
                     int      n_closest):
    """
    """
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.n_closest = n_closest
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.py_comm = comm
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.src_n_points = NULL
    self.tgt_n_points = NULL
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._cls = PDM_closest_points_create(PDMC, n_closest, PDM_OWNERSHIP_USER) # Python take ownership
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int n_part_cloud_src,
                             int n_part_cloud_tgt):
    """
    """
    assert(self.src_n_points  == NULL)
    assert(self.tgt_n_points  == NULL)
    self.src_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_src )
    self.tgt_n_points = <int *> malloc(sizeof(int *) * n_part_cloud_tgt )
    PDM_closest_points_n_part_cloud_set(self._cls,
                                        n_part_cloud_src,
                                        n_part_cloud_tgt)

  # ------------------------------------------------------------------------
  def tgt_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    """
    self.tgt_n_points[i_part] = n_points
    PDM_closest_points_tgt_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)


  # ------------------------------------------------------------------------
  def src_cloud_set(self, int                                           i_part,
                          int                                           n_points,
                          NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                          NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum):
    """
    """
    self.src_n_points[i_part] = n_points
    PDM_closest_points_src_cloud_set(self._cls,
                                     i_part,
                                     n_points,
                           <double*> coords.data,
                      <PDM_g_num_t*> gnum.data)

  # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    PDM_closest_points_compute(self._cls)

  # ------------------------------------------------------------------------
  def points_get(self, int i_part_tgt):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef PDM_g_num_t *closest_src_gnum
    cdef double      *closest_src_distance
    # ************************************************************************

    # > Get
    PDM_closest_points_get(self._cls,
                           i_part_tgt,
                           &closest_src_gnum,
                           &closest_src_distance)

    dim = self.tgt_n_points[i_part_tgt] * self.n_closest
    
    return {'closest_src_gnum'     : create_numpy_g(closest_src_gnum, dim),
            'closest_src_distance' : create_numpy_d(closest_src_distance, dim)}

  # ------------------------------------------------------------------------
  def tgt_in_src_get(self, int i_part_src):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef int *tgt_in_src_idx
    cdef PDM_g_num_t *tgt_in_src
    cdef double *tgt_in_src_dist
    # ************************************************************************

    # > Get
    PDM_closest_points_tgt_in_src_get(self._cls,
                                      i_part_src,
                                      &tgt_in_src_idx,
                                      &tgt_in_src)

    PDM_closest_points_tgt_in_src_dist_get(self._cls,
                                           i_part_src,
                                           &tgt_in_src_idx,
                                           &tgt_in_src_dist)

    src_n_pts = self.src_n_points[i_part_src]
    return {
            'tgt_in_src_idx'  : create_numpy_i(tgt_in_src_idx, src_n_pts+1),
            'tgt_in_src'      : create_numpy_g(tgt_in_src, tgt_in_src_idx[src_n_pts]),
            'tgt_in_src_dist2': create_numpy_d(tgt_in_src_dist, tgt_in_src_idx[src_n_pts])
            }

  # ------------------------------------------------------------------------
  def part_to_part_get(self):
    """
    """
    cdef PDM_part_to_part_t *ptpc
    PDM_closest_points_part_to_part_get(self._cls,
                                        &ptpc,
                                        PDM_OWNERSHIP_USER)

    py_caps = PyCapsule_New(ptpc, NULL, NULL)
    return PartToPartCapsule(py_caps, self.py_comm) # The free is inside the class

  # ------------------------------------------------------------------------
  def dump_times(self):
    """
    """
    PDM_closest_points_dump_times(self._cls)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
    """
    free(self.tgt_n_points)
    free(self.src_n_points)
    PDM_closest_points_free(self._cls)
