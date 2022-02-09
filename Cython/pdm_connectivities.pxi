cdef extern from "pdm_dconnectivity_transform.h":
    void PDM_dconnectivity_transpose(const PDM_MPI_Comm     comm,
                                     const PDM_g_num_t     *entity1_distrib,
                                           PDM_g_num_t     *entity2_distrib,
                                     const int             *dentity1_entity2_idx,
                                     const PDM_g_num_t     *dentity1_entity2,
                                           int              is_signed,
                                           int            **dentity2_entity1_idx,
                                           PDM_g_num_t    **dentity2_entity1)

def dconnectivity_transpose(MPI.Comm comm,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    entity1_distrib,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    entity2_distrib,
                            NPY.ndarray[int           , mode='c', ndim=1]    dentity1_entity2_idx,
                            NPY.ndarray[npy_pdm_gnum_t, mode='c', ndim=1]    dentity1_entity2,
                            bint                                             is_signed):
                                  
    ## entity2_distrib can be recomputed if allocated but entity2_distrib[0] = -1
    # TODO : manage case entity1_distrib recomputed
    # > Convert mpi4py -> PDM_MPI
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    cdef PDM_MPI_Comm PDMC   = PDM_MPI_mpi_2_pdm_mpi_comm(&c_comm)

    cdef PDM_g_num_t * _entity1_distrib
    cdef PDM_g_num_t * _entity2_distrib
    _entity1_distrib = <PDM_g_num_t *> entity1_distrib.data if entity1_distrib is not None else NULL
    _entity2_distrib = <PDM_g_num_t *> entity2_distrib.data if entity2_distrib is not None else NULL

    cdef int *_dentity1_entity2_idx
    _dentity1_entity2_idx = <int *> dentity1_entity2_idx.data if dentity1_entity2_idx is not None else NULL
    cdef PDM_g_num_t *_dentity1_entity2
    _dentity1_entity2 = <PDM_g_num_t *> dentity1_entity2.data if dentity1_entity2 is not None else NULL


    cdef        int * _dentity2_entity1_idx
    cdef PDM_g_num_t* _dentity2_entity1

    PDM_dconnectivity_transpose(PDMC,
                                _entity1_distrib,
                                _entity2_distrib,
                                _dentity1_entity2_idx,
                                _dentity1_entity2,
                         <int>  is_signed,
                              &_dentity2_entity1_idx,
                              &_dentity2_entity1)
    
    dn_entity2 = entity2_distrib[comm.Get_rank()+1] - entity2_distrib[comm.Get_rank()]

    np_dentity2_entity1_idx = create_numpy_i(_dentity2_entity1_idx, dn_entity2 + 1)
    PyArray_ENABLEFLAGS(np_dentity2_entity1_idx, NPY.NPY_OWNDATA)

    np_dentity2_entity1 = create_numpy_pdm_gnum(_dentity2_entity1, np_dentity2_entity1_idx[dn_entity2])
    PyArray_ENABLEFLAGS(np_dentity2_entity1, NPY.NPY_OWNDATA)

    return np_dentity2_entity1_idx, np_dentity2_entity1
