
cdef extern from "pdm_dmesh_nodal_to_dmesh.h":
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of structure
    ctypedef struct PDM_dmesh_nodal_to_dmesh_t:
      pass
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # > Wrapping of function
    PDM_dmesh_nodal_to_dmesh_t* PDM_dmesh_nodal_to_dmesh_create(int             n_mesh,
                                                                PDM_MPI_Comm    comm,
                                                                PDM_ownership_t owner)

    void PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(PDM_dmesh_nodal_to_dmesh_t *dmn_to_dm,
                                                  int                         i_mesh,
                                                  PDM_DMesh_nodal_t          *dmn)

    # void PDM_dmesh_nodal_to_dmesh_compute(PDM_dmesh_nodal_to_dmesh_t*          dmn_to_dm,
    #                                       PDM_dmesh_nodal_to_dmesh_transform_t transform_kind);

    void PDM_dmesh_nodal_to_dmesh_free(PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm)
    # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# ------------------------------------------------------------------
cdef class DMeshNodalToDMesh:
  """
     DistributedMesh: Distributed mesh structure
  """
  # ************************************************************************
  # > Class attributes
  cdef PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm
  # ************************************************************************
  # ------------------------------------------------------------------------
  def __cinit__(self, n_mesh,
                     MPI.Comm comm):
    """
    TODOUX
    """
    # ************************************************************************
    # > Declaration
    cdef MPI.MPI_Comm c_comm = comm.ob_mpi
    # ************************************************************************

    # ::::::::::::::::::::::::::::::::::::::::::::::::::
    self.dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_mesh,
                                                     PDM_MPI_mpi_2_pdm_mpi_comm (<void *> &c_comm),
                                                     PDM_OWNERSHIP_USER) # Python take ownership);
    # ::::::::::::::::::::::::::::::::::::::::::::::::::
  # ------------------------------------------------------------------------
  def add_dmesh_nodal(self, int i_mesh, DistributedMeshNodal dmn):
    """
    """
    PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(self.dmn_to_dm,
                                             i_mesh,
                                             dmn.dmn)

  # ------------------------------------------------------------------------
  def __dealloc__(self):
    """
       Use the free method of PDM Lib
    """
    # ************************************************************************
    # > Declaration
    # ************************************************************************
    PDM_dmesh_nodal_to_dmesh_free(self.dmn_to_dm)

