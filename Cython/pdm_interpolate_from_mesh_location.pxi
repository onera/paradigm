
cdef extern from "pdm_interpolate_from_mesh_location.h":

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_interpolate_from_mesh_location_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  PDM_interpolate_from_mesh_location_t* PDM_interpolate_from_mesh_location_create(int                    n_part_src,
                                                                                  int                    n_cloud_target,
                                                                                  PDM_interpolate_kind_t interp_kind,
                                                                                  PDM_MPI_Comm           comm);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_compute(PDM_interpolate_from_mesh_location_t  *interp_from_ml)
  void PDM_interpolate_from_mesh_location_exch(PDM_interpolate_from_mesh_location_t  *interp_from_ml,
                                               int i_point_cloud,
                                               size_t s_data,
                                               double **part_data_in,
                                               double ***cloud_data_out)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_n_part_cloud_set(PDM_interpolate_from_mesh_location_t *interp_from_ml,
                                                           int                                   i_point_cloud,
                                                           int                                   n_part);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_cloud_set(PDM_interpolate_from_mesh_location_t *interp_from_ml,
                                                    int                                   i_point_cloud,
                                                    int                                   i_part,
                                                    int                                   n_points,
                                                    double                               *coords,
                                                    PDM_g_num_t                          *gnum)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_part_set(PDM_interpolate_from_mesh_location_t *interp_from_ml,
                                                   int                                   i_part,
                                                   int                                   n_cell,
                                                   int                                  *cell_face_idx,
                                                   int                                  *cell_face,
                                                   PDM_g_num_t                          *cell_ln_to_gn,
                                                   int                                   n_face,
                                                   int                                  *face_vtx_idx,
                                                   int                                  *face_vtx,
                                                   PDM_g_num_t                          *face_ln_to_gn,
                                                   int                                   n_vtx,
                                                   double                               *coords,
                                                   PDM_g_num_t                          *vtx_ln_to_gn);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_points_in_elt_set(PDM_interpolate_from_mesh_location_t *interp_from_ml,
                                                            int                                   i_part,
                                                            int                                   i_point_cloud,
                                                            int                                   *elt_pts_inside_idx,
                                                            PDM_g_num_t                           *points_gnum,
                                                            double                                *points_coords,
                                                            double                                *points_uvw,
                                                            int                                   *points_weights_idx,
                                                            double                                *points_weights,
                                                            double                                *points_dist2,
                                                            double                                *points_projected_coords);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_interpolate_from_mesh_location_free(PDM_interpolate_from_mesh_location_t *interp_from_ml);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class InterpolateFromMeshLocation:
  """
     PointsMerge: Interface to build connection between multiple cloud in parallel
     Useful for get connection between partiton from faces coordinates
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_interpolate_from_mesh_location_t* _interp_from_ml
  cdef int _size
  cdef int _rank
  cdef int _n_part_src
  cdef int _n_point_cloud
  # ************************************************************************

  # ------------------------------------------------------------------------
  def __init__(self, PDM_mesh_nature_t mesh_nature,
                     int               n_part_src,
                     int               n_point_cloud,
                     MPI.Comm          comm):
    """
    """
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._n_part_src = n_part_src
    self._n_point_cloud = n_point_cloud
    self.n_part_cloud = dict()
    self.n_point_cloud = [[dict()]*n_point_cloud]
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._rank = comm.Get_rank()
    self._size = comm.Get_size()
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(<void *> &c_comm)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self._interp_from_ml = PDM_interpolate_from_mesh_location_create(n_part_src, n_point_cloud, PDM_INTERPOLATE_KIND_FROM_CENTER, PDMC)
    # ::::::::::::::::::::::::::::::::::::::::::::::::::

  # ------------------------------------------------------------------------
  def n_part_cloud_set(self, int i_point_cloud,
                             int n_part):
    """
    """
    self.n_part_cloud[i_point_cloud] = n_part
    PDM_interpolate_from_mesh_location_n_part_cloud_set(self._interp_from_ml,
                                                        i_point_cloud,
                                                        n_part);

  # # ------------------------------------------------------------------------
  def cloud_set(self, int i_point_cloud,
                      int i_part,
                      int n_points,
                      NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                      NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] gnum,
                      ):
    """
    """
    self.n_point_cloud[i_point_cloud][i_part] = n_points
    PDM_interpolate_from_mesh_location_cloud_set(self._interp_from_ml,
                                                 i_point_cloud,
                                                 i_part,
                                                 n_points,
                                  <double*>      coords.data,
                                  <PDM_g_num_t*> gnum.data);

  # ------------------------------------------------------------------------
  def part_set(self, int i_part,
                     int n_cell,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                     int n_face,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx,
                     NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                     int n_vtx,
                     NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords,
                     NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn):
    """
    """
    PDM_interpolate_from_mesh_location_part_set(self._interp_from_ml,
                                                i_part,
                                                n_cell,
                                 <int*>         cell_face_idx.data,
                                 <int*>         cell_face.data,
                                 <PDM_g_num_t*> cell_ln_to_gn.data,
                                                n_face,
                                 <int*>         face_vtx_idx.data,
                                 <int*>         face_vtx.data,
                                 <PDM_g_num_t*> face_ln_to_gn.data,
                                                n_vtx,
                                 <double*>      coords.data,
                                 <PDM_g_num_t*> vtx_ln_to_gn.data);

  # # ------------------------------------------------------------------------
  def compute(self):
    """
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_interpolate_from_mesh_location_compute(self._interp_from_ml)

  # ------------------------------------------------------------------------
  def exch(self,
           int i_point_cloud,
           list_part_data_in):
    """
    """
    # ************************************************************************
    # > Declaration
    cdef double** pdata_data_in
    cdef double** cloud_data_out
    cdef NPY.npy_intp dim
    # ************************************************************************

    assert(list_part_data_in == self._n_part_src)

    pdata_data_in = <double **> malloc(self._n_part_src * sizeof(double **))

    PDM_interpolate_from_mesh_location_exch(self._interp_from_ml,
                                            i_point_cloud,
                                            sizeof(double),
                                            pdata_data_in,
                                            &cloud_data_out)

    n_part_out = self.n_part_cloud[i_point_cloud]
    list_cloud_data_out = list()
    # > Il faut le nombre de point par cloud et par partition ...
    for i_part in range(n_part_out):
      n_point_cloud = self.n_point_cloud[i_point_cloud][i_part]
      dim = <NPY.npy_intp> n_point_cloud
      np_cloud_data_out = NPY.PyArray_SimpleNewFromData(1,
                                                        &dim,
                                                        NPY.NPY_DOUBLE,
                                               <void *> cloud_data_out[i_part])
      PyArray_ENABLEFLAGS(np_cloud_data_out, NPY.NPY_OWNDATA);

    free(cloud_data_out) # Free the part indirection other is ok with numpy

    return list_cloud_data_out


  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    print('PDM_interpolate_from_mesh_location')
    PDM_interpolate_from_mesh_location_free(self._interp_from_ml)
