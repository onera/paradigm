
cdef extern from "pdm_dmesh_nodal.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of Ppart Structure 
    ctypedef struct PDM_DMesh_nodal_t:
      pass
      
    ctypedef enum PDM_Mesh_nodal_elt_t: 
      PDM_MESH_NODAL_POINT    = 1   
      PDM_MESH_NODAL_BAR2     = 2  
      PDM_MESH_NODAL_TRIA3    = 3   
      PDM_MESH_NODAL_QUAD4    = 4   
      PDM_MESH_NODAL_POLY_2D  = 5     
      PDM_MESH_NODAL_TETRA4   = 6    
      PDM_MESH_NODAL_PYRAMID5 = 7      
      PDM_MESH_NODAL_PRISM6   = 8    
      PDM_MESH_NODAL_HEXA8    = 9   
      PDM_MESH_NODAL_POLY_3D  = 10 
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function 
    int PDM_DMesh_nodal_create(PDM_MPI_Comm comm)
    void PDM_DMesh_nodal_free(int handle)
    
    void PDM_DMesh_nodal_coord_set(int handle, int n_vtx, double* coords)

    PDM_g_num_t *PDM_DMesh_nodal_distrib_vtx_get(int handle)
    PDM_g_num_t *PDM_DMesh_nodal_distrib_section_get(int handle, int id_section)
    
    int                  PDM_DMesh_nodal_n_vtx_get(int handle)
    int                  PDM_DMesh_nodal_n_sections_get(int handle)
    int*                 PDM_DMesh_nodal_sections_id_get(int handle)
    PDM_Mesh_nodal_elt_t PDM_Mesh_nodal_section_type_get(int handle, int id_section)
    double*              PDM_DMesh_nodal_vtx_get(int handle)
    
    int                  PDM_Mesh_nodal_section_add(int handle, PDM_Mesh_nodal_elt_t t_elt)
    void                 PDM_Mesh_nodal_section_std_set(int          handle, 
                                                        int          id_section,
                                                        int          n_elmts, 
                                                        PDM_g_num_t* connec)
    
    PDM_g_num_t* PDM_DMesh_nodal_section_std_get(int          handle, int id_section)
    int PDM_DMesh_nodal_section_n_elt_get(int          handle, int id_section)
    
    void PDM_DMesh_nodal_section_poly2d_set(int handle, int id_section, PDM_l_num_t n_elt, 
                                            PDM_l_num_t* connec_idx, 
                                            PDM_g_num_t   *connec)
    
    PDM_g_num_t PDM_DMesh_nodal_total_n_cell_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_face_get(int handle)
    PDM_g_num_t PDM_DMesh_nodal_total_n_vtx_get(int handle)
    void PDM_DMesh_nodal_cell_face_compute(int handle)
    int PDM_DMesh_nodal_cell_face_get(int handle, int** cell_face_idx, PDM_g_num_t **cell_face)
    int PDM_DMesh_nodal_face_vtx_get(int handle, int** dface_vtx_idx, PDM_g_num_t **dface_vtx)
    
    PDM_g_num_t* PDM_DMesh_nodal_distrib_cell_get(int handle)
    PDM_g_num_t* PDM_DMesh_nodal_distrib_face_get(int handle)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# ------------------------------------------------------------------
cdef class DistributedMeshNodal:
    """
       DistributedMeshNodal: Interface to build face from Element->Vtx connectivity
    """
    # ************************************************************************
    # > Class attributes
    cdef PDM_part_to_block_t *PTB
    cdef int                  partN
    cdef int                  Size
    cdef int                  Rank

    cdef int                 *NbElmts
    cdef PDM_g_num_t        **LNToGN

    cdef PDM_part_to_block_distrib_t t_distrib
    cdef PDM_part_to_block_post_t    t_post
    cdef PDM_writer_part_stride_t    t_stride
    # ************************************************************************
    # ------------------------------------------------------------------------
    def __cinit__(self, MPI.Comm comm, list pLNToGN, int partN,
                        PDM_part_to_block_distrib_t t_distrib = <PDM_part_to_block_distrib_t> (0),
                        PDM_part_to_block_post_t    t_post    = <PDM_part_to_block_post_t   > (0),
                        PDM_writer_part_stride_t    t_stride  = <PDM_writer_part_stride_t   > (0),
                        float partActiveNode = 1.):
        """
        TODOUX
        """
        # ************************************************************************
        # > Declaration
        cdef int      nElts
        cdef int      idx
        # > Numpy array
        cdef NPY.ndarray[npy_pdm_gnum_t, ndim=1, mode='fortran'] partLNToGN
        # ************************************************************************
        
        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # ::::::::::::::::::::::::::::::::::::::::::::::::::

    # ------------------------------------------------------------------------
    def DistributedMeshNodal_Compute(self):
        """
        """
        # ************************************************************************
        # > Declaration
        # ************************************************************************

        # ::::::::::::::::::::::::::::::::::::::::::::::::::
        # ::::::::::::::::::::::::::::::::::::::::::::::::::


    # ------------------------------------------------------------------------
    def getDistributionCopy(self):
      """
         Return a copy of the distrisbution array compute in library
         Copy because remove of PTB object can made a core ...
      """
      # ************************************************************************
      # > Declaration
      cdef PDM_g_num_t* Distrib
      # ************************************************************************

      # ::::::::::::::::::::::::::::::::::::::::::::::::::
      # > Get
      # Distrib = PDM_part_to_block_distrib_index_get(self.PTB)
      # ::::::::::::::::::::::::::::::::::::::::::::::::::


    # ------------------------------------------------------------------------
    def __dealloc__(self):
      """
         Use the free method of PDM Lib
      """
      # ************************************************************************
      # > Declaration
      # cdef PDM_part_to_block_t *a
      # ************************************************************************

      # > Free Ppart Structure
      # a = PDM_part_to_block_free(self.PTB)

      # > Free allocated array
      # free(self.LNToGN)
      # free(self.NbElmts)    

