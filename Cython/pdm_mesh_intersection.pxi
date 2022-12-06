
cdef extern from "pdm_mesh_intersection.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure
    ctypedef struct PDM_mesh_intersection_t:
        pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_mesh_intersection_t* PDM_mesh_intersection_create(PDM_mesh_intersection_kind_t intersection_kind,
                                                          int                          dim_mesh_a,
                                                          int                          dim_mesh_b,
                                                          int                          n_part_mesh_a,
                                                          int                          n_part_mesh_b,
                                                          double                       project_coeff,
                                                          PDM_MPI_Comm                 comm);
    void PDM_mesh_intersection_part_set(PDM_mesh_intersection_t  *mi,
                                        PDM_ol_mesh_t             i_mesh,
                                        int                       i_part,
                                        int                       n_cell,
                                        int                       n_face,
                                        int                       n_edge,
                                        int                       n_vtx,
                                        int                      *cell_face_idx,
                                        int                      *cell_face,
                                        int                      *face_edge_idx,
                                        int                      *face_edge,
                                        int                      *edge_vtx,
                                        int                      *face_vtx_idx,
                                        int                      *face_vtx,
                                        PDM_g_num_t              *cell_ln_to_gn,
                                        PDM_g_num_t              *face_ln_to_gn,
                                        PDM_g_num_t              *edge_ln_to_gn,
                                        PDM_g_num_t              *vtx_ln_to_gn,
                                        double                   *vtx_coord);
    void PDM_mesh_intersection_compute(PDM_mesh_intersection_t  *mi);
    void PDM_mesh_intersection_free(PDM_mesh_intersection_t* mi);

# ------------------------------------------------------------------
cdef class MeshIntersection:
    """
    Define a method to compute the distance from a cloud to a surface
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_mesh_intersection_t* _mi
    cdef MPI.Comm py_comm
    cdef dict ptp_objects
    cdef list keep_alive
    # ************************************************************************
    # ------------------------------------------------------------------
    def __cinit__(self,
                  MPI.Comm                     comm,
                  PDM_mesh_intersection_kind_t intersection_kind,
                  int                          dim_mesh_a,
                  int                          dim_mesh_b,
                  int                          n_part_mesh_a,
                  int                          n_part_mesh_b,
                  double                       project_coeff=1.e-6):
        """
        Compute the distance from point clouds to a surface
        """
        self.keep_alive  = list()

        # Convert mpi4py -> PDM_MPI
        self.py_comm = comm
        cdef MPI.MPI_Comm c_comm   = self.py_comm.ob_mpi
        cdef PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)


        self._mi = PDM_mesh_intersection_create(intersection_kind,
                                                dim_mesh_a,
                                                dim_mesh_b,
                                                n_part_mesh_a,
                                                n_part_mesh_b,
                                                project_coeff,
                                                pdm_comm);
    # ------------------------------------------------------------------
    def part_set(self,
                 int                                           i_mesh,
                 int                                           i_part,
                 int                                           n_cell,
                 int                                           n_face,
                 int                                           n_edge,
                 int                                           n_vtx,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face_idx,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] cell_face    ,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge_idx,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_edge    ,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] edge_vtx     ,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx_idx ,
                 NPY.ndarray[NPY.int32_t   , mode='c', ndim=1] face_vtx     ,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] cell_ln_to_gn,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] face_ln_to_gn,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] edge_ln_to_gn,
                 NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1] vtx_ln_to_gn ,
                 NPY.ndarray[NPY.double_t  , mode='c', ndim=1] coords):
      """
      """
      self.keep_alive.append(cell_face_idx)
      self.keep_alive.append(cell_face)
      self.keep_alive.append(face_edge_idx)
      self.keep_alive.append(face_edge)
      self.keep_alive.append(edge_vtx)
      self.keep_alive.append(face_vtx_idx)
      self.keep_alive.append(face_vtx)
      self.keep_alive.append(cell_ln_to_gn)
      self.keep_alive.append(face_ln_to_gn)
      self.keep_alive.append(edge_ln_to_gn)
      self.keep_alive.append(vtx_ln_to_gn)
      self.keep_alive.append(coords)

      PDM_mesh_intersection_part_set(self._mi,
                    <PDM_ol_mesh_t>  i_mesh,
                                     i_part,
                                     n_cell,
                                     n_face,
                                     n_edge,
                                     n_vtx,
                    <int         *>  cell_face_idx.data,
                    <int         *>  cell_face    .data,
                    <int         *>  face_edge_idx.data,
                    <int         *>  face_edge    .data,
                    <int         *>  edge_vtx     .data,
                    <int         *>  face_vtx_idx .data,
                    <int         *>  face_vtx     .data,
                    <PDM_g_num_t *>  cell_ln_to_gn.data,
                    <PDM_g_num_t *>  face_ln_to_gn.data,
                    <PDM_g_num_t *>  edge_ln_to_gn.data,
                    <PDM_g_num_t *>  vtx_ln_to_gn .data,
                    <double      *>  coords       .data)


    # ------------------------------------------------------------------
    def compute(self):
        """
        """
        PDM_mesh_intersection_compute(self._mi)

    # ------------------------------------------------------------------
    def __dealloc__(self):
        """
        Free a mesh_intersection structure
        """
        PDM_mesh_intersection_free(self._mi)
