cdef extern from "pdm_sphere_surf_gen.h":
    void PDM_sphere_surf_icosphere_gen(const PDM_MPI_Comm   comm,
                                       const PDM_g_num_t    n,
                                       const double         x_center,
                                       const double         y_center,
                                       const double         z_center,
                                       const double         radius,
                                             double       **dvtx_coord,
                                             int          **dface_vtx_idx,
                                             PDM_g_num_t  **dface_vtx,
                                             PDM_g_num_t  **distrib_vtx,
                                             PDM_g_num_t  **distrib_face)

    void PDM_sphere_surf_icosphere_gen_part(const PDM_MPI_Comm        comm,
                                            const PDM_g_num_t         n,
                                            const double              x_center,
                                            const double              y_center,
                                            const double              z_center,
                                            const double              radius,
                                            const int                 n_part,
                                            const PDM_split_dual_t    part_method,
                                                  int               **pn_vtx,
                                                  double           ***pvtx_coord,
                                                  PDM_g_num_t      ***pvtx_ln_to_gn,
                                                  int               **pn_face,
                                                  int              ***pface_vtx_idx,
                                                  int              ***pface_vtx,
                                                  PDM_g_num_t      ***pface_ln_to_gn)


# ------------------------------------------------------------------------

def sphere_surf_icosphere_gen(MPI.Comm       comm,
                              npy_pdm_gnum_t n,
                              NPY.double_t   x_center,
                              NPY.double_t   y_center,
                              NPY.double_t   z_center,
                              NPY.double_t   radius):

  i_rank = comm.rank
  n_rank = comm.size

  cdef double      *dvtx_coord    = NULL
  cdef int         *dface_vtx_idx = NULL
  cdef PDM_g_num_t *dface_vtx     = NULL
  cdef PDM_g_num_t *distrib_vtx   = NULL
  cdef PDM_g_num_t *distrib_face  = NULL

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  PDM_sphere_surf_icosphere_gen(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                n,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &dvtx_coord,
                                &dface_vtx_idx,
                                &dface_vtx,
                                &distrib_vtx,
                                &distrib_face)

  dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank]
  dn_face = distrib_face[i_rank+1] - distrib_face[i_rank]

  dim = <NPY.npy_intp> (dn_vtx * 3)
  np_dvtx_coord = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, dvtx_coord)
  PyArray_ENABLEFLAGS(np_dvtx_coord, NPY.NPY_OWNDATA)

  dim = <NPY.npy_intp> (dn_face+1)
  np_dface_vtx_idx = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT32, dface_vtx_idx)
  PyArray_ENABLEFLAGS(np_dface_vtx_idx, NPY.NPY_OWNDATA)

  dim = <NPY.npy_intp> dface_vtx_idx[dn_face]
  np_dface_vtx = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT32, dface_vtx)
  PyArray_ENABLEFLAGS(np_dface_vtx, NPY.NPY_OWNDATA)

  dim = <NPY.npy_intp> (n_rank+1)
  np_distrib_vtx = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, distrib_vtx)
  PyArray_ENABLEFLAGS(np_distrib_vtx, NPY.NPY_OWNDATA)

  dim = <NPY.npy_intp> (n_rank+1)
  np_distrib_face = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, distrib_face)
  PyArray_ENABLEFLAGS(np_distrib_face, NPY.NPY_OWNDATA)

  return {
  "dvtx_coord"    : np_dvtx_coord,
  "dface_vtx_idx" : np_dface_vtx_idx,
  "dface_vtx"     : np_dface_vtx,
  "distrib_vtx"   : np_distrib_vtx,
  "distrib_face"  : np_distrib_face
  }

# ------------------------------------------------------------------------

def sphere_surf_icosphere_gen_part(MPI.Comm         comm,
                                   npy_pdm_gnum_t   n,
                                   NPY.double_t     x_center,
                                   NPY.double_t     y_center,
                                   NPY.double_t     z_center,
                                   NPY.double_t     radius,
                                   int              n_part,
                                   PDM_split_dual_t part_method):

  cdef int          *_pn_vtx         = NULL
  cdef double      **_pvtx_coord     = NULL
  cdef PDM_g_num_t **_pvtx_ln_to_gn  = NULL
  cdef int          *_pn_face        = NULL
  cdef int         **_pface_vtx_idx  = NULL
  cdef int         **_pface_vtx      = NULL
  cdef PDM_g_num_t **_pface_ln_to_gn = NULL

  cdef MPI.MPI_Comm c_comm = comm.ob_mpi
  PDM_sphere_surf_icosphere_gen_part(PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                     n,
                                     x_center,
                                     y_center,
                                     z_center,
                                     radius,
                                     n_part,
                                     part_method,
                                     &_pn_vtx,
                                     &_pvtx_coord,
                                     &_pvtx_ln_to_gn,
                                     &_pn_face,
                                     &_pface_vtx_idx,
                                     &_pface_vtx,
                                     &_pface_ln_to_gn)

  pn_vtx  = create_numpy_i(_pn_vtx,  n_part)
  pn_face = create_numpy_i(_pn_face, n_part)

  pvtx_coord     = []
  pvtx_ln_to_gn  = []
  pface_vtx_idx  = []
  pface_vtx      = []
  pface_ln_to_gn = []
  for i_part in range(n_part):
    dim = <NPY.npy_intp> (pn_vtx[i_part] * 3)
    data = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_DOUBLE, _pvtx_coord[i_part])
    PyArray_ENABLEFLAGS(data, NPY.NPY_OWNDATA);
    pvtx_coord.append(data)

    dim = <NPY.npy_intp> pn_vtx[i_part]
    data = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, _pvtx_ln_to_gn[i_part])
    PyArray_ENABLEFLAGS(data, NPY.NPY_OWNDATA);
    pvtx_ln_to_gn.append(data)


    dim = <NPY.npy_intp> (pn_face[i_part]+1)
    data = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT32, _pface_vtx_idx[i_part])
    PyArray_ENABLEFLAGS(data, NPY.NPY_OWNDATA);
    pface_vtx_idx.append(data)

    dim = <NPY.npy_intp> pface_vtx_idx[i_part][pn_face[i_part]]
    data = NPY.PyArray_SimpleNewFromData(1, &dim, NPY.NPY_INT32, _pface_vtx[i_part])
    PyArray_ENABLEFLAGS(data, NPY.NPY_OWNDATA);
    pface_vtx.append(data)

    dim = <NPY.npy_intp> pn_face[i_part]
    data = NPY.PyArray_SimpleNewFromData(1, &dim, PDM_G_NUM_NPY_INT, _pface_ln_to_gn[i_part])
    PyArray_ENABLEFLAGS(data, NPY.NPY_OWNDATA);
    pface_ln_to_gn.append(data)

  return {
  "pn_vtx"         : pn_vtx,
  "pvtx_coord"     : pvtx_coord,
  "pvtx_ln_to_gn"  : pvtx_ln_to_gn,
  "pn_face"        : pn_face,
  "pface_vtx_idx"  : pface_vtx_idx,
  "pface_vtx"      : pface_vtx,
  "pface_ln_to_gn" : pface_ln_to_gn
  }

