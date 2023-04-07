cdef extern from "pdm_point_cloud_gen.h":
   void PDM_point_cloud_gen_random(PDM_MPI_Comm        comm,
                                   const int           seed,
                                   const int           geometric_g_num,
                                   const PDM_g_num_t   gn_pts,
                                   const double        x_min,
                                   const double        y_min,
                                   const double        z_min,
                                   const double        x_max,
                                   const double        y_max,
                                   const double        z_max,
                                   int                *pn_pts,
                                   double            **coord,
                                   PDM_g_num_t       **g_num);

   void PDM_point_cloud_gen_cartesian(PDM_MPI_Comm        comm,
                                      const int           nx,
                                      const int           ny,
                                      const int           nz,
                                      const double        x_min,
                                      const double        y_min,
                                      const double        z_min,
                                      const double        x_max,
                                      const double        y_max,
                                      const double        z_max,
                                      int                *n_pts,
                                      double            **pts_coord,
                                      PDM_g_num_t       **pts_ln_to_gn);

   void PDM_dpoint_cloud_gen_random(PDM_MPI_Comm        comm,
                                    const int           seed,
                                    const PDM_g_num_t   gn_pts,
                                    const double        x_min,
                                    const double        y_min,
                                    const double        z_min,
                                    const double        x_max,
                                    const double        y_max,
                                    const double        z_max,
                                    double            **dpts_coord,
                                    PDM_g_num_t       **distrib_pts);

   void PDM_point_cloud_gen_cartesian(PDM_MPI_Comm        comm,
                                      const int           nx,
                                      const int           ny,
                                      const int           nz,
                                      const double        x_min,
                                      const double        y_min,
                                      const double        z_min,
                                      const double        x_max,
                                      const double        y_max,
                                      const double        z_max,
                                      int                *n_pts,
                                      double            **pts_coord,
                                      PDM_g_num_t       **pts_ln_to_gn);

   void PDM_dpoint_cloud_gen_cartesian(PDM_MPI_Comm        comm,
                                       const int           nx,
                                       const int           ny,
                                       const int           nz,
                                       const double        x_min,
                                       const double        y_min,
                                       const double        z_min,
                                       const double        x_max,
                                       const double        y_max,
                                       const double        z_max,
                                       double            **dpts_coord,
                                       PDM_g_num_t       **distrib_pts);

# ------------------------------------------------------------------------
def part_point_cloud_gen_random(MPI.Comm       comm,
                                int            seed,
                                int            geometric_g_num,
                                npy_pdm_gnum_t gn_pts,
                                NPY.double_t   x_min,
                                NPY.double_t   y_min,
                                NPY.double_t   z_min,
                                NPY.double_t   x_max,
                                NPY.double_t   y_max,
                                NPY.double_t   z_max):
   """
   Generate an partition cloud point
   """
   cdef int          pn_pts       = 0
   cdef double      *pts_coord    = NULL
   cdef PDM_g_num_t *pts_ln_to_gn = NULL

   cdef MPI.MPI_Comm c_comm = comm.ob_mpi
   PDM_point_cloud_gen_random(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                              seed,
                              geometric_g_num,
                              gn_pts,
                              x_min,
                              y_min,
                              z_min,
                              x_max,
                              y_max,
                              z_max,
                              &pn_pts,
                              &pts_coord,
                              &pts_ln_to_gn);

   dim = <NPY.npy_intp> (pn_pts * 3)
   np_pts_coord = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, pts_coord)
   PyArray_ENABLEFLAGS(np_pts_coord, NPY.NPY_OWNDATA)

   dim = <NPY.npy_intp> (pn_pts )
   np_pts_ln_to_gn = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, pts_ln_to_gn)
   PyArray_ENABLEFLAGS(np_pts_ln_to_gn, NPY.NPY_OWNDATA)

   return {"np_pts_ln_to_gn" : np_pts_ln_to_gn,
   "np_pts_coord"    : np_pts_coord}

# ------------------------------------------------------------------------
def part_point_cloud_gen_cartesian(MPI.Comm      comm,
                                   int           nx,
                                   int           ny,
                                   int           nz,
                                   NPY.double_t  x_min,
                                   NPY.double_t  y_min,
                                   NPY.double_t  z_min,
                                   NPY.double_t  x_max,
                                   NPY.double_t  y_max,
                                   NPY.double_t  z_max):
   """
   Generate an partition cloud point
   """
   cdef int          pn_pts       = 0
   cdef double      *pts_coord    = NULL
   cdef PDM_g_num_t *pts_ln_to_gn = NULL

   cdef MPI.MPI_Comm c_comm = comm.ob_mpi
   PDM_point_cloud_gen_cartesian(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                 nx,
                                 ny,
                                 nz,
                                 x_min,
                                 y_min,
                                 z_min,
                                 x_max,
                                 y_max,
                                 z_max,
                                 &pn_pts,
                                 &pts_coord,
                                 &pts_ln_to_gn);

   dim = <NPY.npy_intp> (pn_pts * 3)
   np_pts_coord = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, pts_coord)
   PyArray_ENABLEFLAGS(np_pts_coord, NPY.NPY_OWNDATA)

   dim = <NPY.npy_intp> (pn_pts )
   np_pts_ln_to_gn = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, pts_ln_to_gn)
   PyArray_ENABLEFLAGS(np_pts_ln_to_gn, NPY.NPY_OWNDATA)

   return {"np_pts_ln_to_gn" : np_pts_ln_to_gn,
   "np_pts_coord"    : np_pts_coord}



# ------------------------------------------------------------------------
def dpoint_cloud_gen_random(MPI.Comm       comm,
                            int            seed,
                            npy_pdm_gnum_t gn_pts,
                            NPY.double_t   x_min,
                            NPY.double_t   y_min,
                            NPY.double_t   z_min,
                            NPY.double_t   x_max,
                            NPY.double_t   y_max,
                            NPY.double_t   z_max):
   """
   Generate an partition cloud point
   """
   i_rank = comm.rank
   n_rank = comm.size
   cdef double      *dpts_coord  = NULL
   cdef PDM_g_num_t *distrib_pts = NULL

   cdef MPI.MPI_Comm c_comm = comm.ob_mpi
   PDM_dpoint_cloud_gen_random(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                               seed,
                               gn_pts,
                               x_min,
                               y_min,
                               z_min,
                               x_max,
                               y_max,
                               z_max,
                               &dpts_coord,
                               &distrib_pts);

   dn_pts  = distrib_pts [i_rank+1] - distrib_pts [i_rank]
   dim = <NPY.npy_intp> (dn_pts * 3)
   np_dpts_coord = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, dpts_coord)
   PyArray_ENABLEFLAGS(np_dpts_coord, NPY.NPY_OWNDATA)

   dim = <NPY.npy_intp> (n_rank + 1 )
   np_distrib_pts = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, distrib_pts)
   PyArray_ENABLEFLAGS(np_distrib_pts, NPY.NPY_OWNDATA)

   return {"np_distrib_pts" : np_distrib_pts,
           "np_dpts_coord"    : np_dpts_coord}

# ------------------------------------------------------------------------
def dpoint_cloud_gen_cartesian(MPI.Comm      comm,
                               int           nx,
                               int           ny,
                               int           nz,
                               NPY.double_t  x_min,
                               NPY.double_t  y_min,
                               NPY.double_t  z_min,
                               NPY.double_t  x_max,
                               NPY.double_t  y_max,
                               NPY.double_t  z_max):
   """
   Generate an partition cloud point
   """
   i_rank = comm.rank
   n_rank = comm.size
   cdef double      *dpts_coord  = NULL
   cdef PDM_g_num_t *distrib_pts = NULL

   cdef MPI.MPI_Comm c_comm = comm.ob_mpi
   PDM_dpoint_cloud_gen_cartesian(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                  nx,
                                  ny,
                                  nz,
                                  x_min,
                                  y_min,
                                  z_min,
                                  x_max,
                                  y_max,
                                  z_max,
                                  &dpts_coord,
                                  &distrib_pts);

   dn_pts  = distrib_pts [i_rank+1] - distrib_pts [i_rank]
   dim = <NPY.npy_intp> (dn_pts * 3)
   np_dpts_coord = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, dpts_coord)
   PyArray_ENABLEFLAGS(np_dpts_coord, NPY.NPY_OWNDATA)

   dim = <NPY.npy_intp> (n_rank + 1 )
   np_distrib_pts = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, distrib_pts)
   PyArray_ENABLEFLAGS(np_distrib_pts, NPY.NPY_OWNDATA)

   return {"np_distrib_pts" : np_distrib_pts,
           "np_dpts_coord"  : np_dpts_coord}
