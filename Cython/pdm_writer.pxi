
cdef extern from "pdm_writer.h":
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping of Ppart Structure
  ctypedef struct PDM_writer_t:
    pass
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C
  ctypedef enum PDM_writer_elt_geom_t:
    PDM_WRITER_POINT    = PDM_MESH_NODAL_POINT
    PDM_WRITER_BAR2     = PDM_MESH_NODAL_BAR2
    PDM_WRITER_TRIA3    = PDM_MESH_NODAL_TRIA3
    PDM_WRITER_QUAD4    = PDM_MESH_NODAL_QUAD4
    PDM_WRITER_POLY_2D  = PDM_MESH_NODAL_POLY_2D
    PDM_WRITER_TETRA4   = PDM_MESH_NODAL_TETRA4
    PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5
    PDM_WRITER_PRISM6   = PDM_MESH_NODAL_PRISM6
    PDM_WRITER_HEXA8    = PDM_MESH_NODAL_HEXA8
    PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D


  ctypedef enum PDM_writer_fmt_fic_t: 
    PDM_WRITER_FMT_BIN   = 0
    PDM_WRITER_FMT_ASCII = 1


  ctypedef enum PDM_writer_var_dim_t:
    PDM_WRITER_VAR_CST           = 0
    PDM_WRITER_VAR_SCALAR        = 1
    PDM_WRITER_VAR_VECTOR        = 3
    PDM_WRITER_VAR_TENSOR_SYM    = 6
    PDM_WRITER_VAR_TENSOR_ASYM   = 9

  ctypedef enum PDM_writer_var_loc_t:
    PDM_WRITER_VAR_VERTICES     = 0
    PDM_WRITER_VAR_ELEMENTS     = 1
    PDM_WRITER_VAR_PARTICLES    = 2


  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  PDM_writer_t * PDM_writer_create(char                   *fmt,
                                   PDM_writer_fmt_fic_t    fmt_fic,
                                   PDM_writer_topology_t  topologie,
                                   PDM_writer_status_t     st_reprise,
                                   char                   *rep_sortie,
                                   char                   *nom_sortie,
                                   PDM_MPI_Comm            pdm_mpi_comm,
                                   PDM_io_acces_t          acces,
                                   double                  prop_noeuds_actifs,
                                   char                   *options)
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  void PDM_writer_free(PDM_writer_t *cs)

  void PDM_writer_step_beg(PDM_writer_t  *cs, double   physical_time)

  void PDM_writer_step_end(PDM_writer_t  *cs)

  int PDM_writer_geom_create(PDM_writer_t               *cs,
                             char                 *nom_geom,
                             PDM_writer_status_t   st_decoup_poly2d,
                             PDM_writer_status_t   st_decoup_poly3d,
                             int                   n_part)

  int PDM_writer_geom_create_from_mesh_nodal(PDM_writer_t        *cs,
                                            char                *nom_geom,
                                            PDM_writer_status_t  st_decoup_poly2d,
                                            PDM_writer_status_t  st_decoup_poly3d,
                                            PDM_Mesh_nodal_t    *mesh)


  void PDM_writer_geom_coord_set(PDM_writer_t  *cs,
                                 int            id_geom,
                                 int            id_part,
                                 int            n_som,
                                 PDM_real_t    *coords,
                                 PDM_g_num_t   *numabs)


  void PDM_writer_geom_coord_from_parent_set(PDM_writer_t *cs,
                                             int           id_geom,
                                             int           id_part,
                                             int           n_som,
                                             int           n_som_parent,
                                             PDM_g_num_t  *numabs,
                                             int          *num_parent,
                                             PDM_real_t   *coords_parent,
                                             PDM_g_num_t  *numabs_parent)


  int PDM_writer_geom_bloc_add(PDM_writer_t          *cs,
                               int                    id_geom,
                               PDM_writer_status_t    st_free_data,
                               PDM_writer_elt_geom_t  t_elt)
  

  void PDM_writer_geom_bloc_std_set(PDM_writer_t  *cs,
                                    int      id_geom,
                                    int      id_bloc,
                                    int      id_part,
                                    int      n_elt,
                                    PDM_l_num_t   *connec,
                                    PDM_g_num_t   *numabs)
  

  void PDM_writer_geom_bloc_poly2d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t   *connec_idx,
                                       PDM_l_num_t   *connec,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_bloc_poly3d_set(PDM_writer_t  *cs,
                                       int            id_geom,
                                       int            id_bloc,
                                       int            id_part,
                                       PDM_l_num_t    n_elt,
                                       PDM_l_num_t    n_face,
                                       PDM_l_num_t   *facsom_idx,
                                       PDM_l_num_t   *facsom,
                                       PDM_l_num_t   *cellfac_idx,
                                       PDM_l_num_t   *cellfac,
                                       PDM_g_num_t   *numabs)
  
  void PDM_writer_geom_cell3d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)
  

  void PDM_writer_geom_cell2d_cellface_add(PDM_writer_t *cs,
                                           int           id_geom,
                                           int           id_part,
                                           int           n_cell,
                                           int           n_face,
                                           PDM_l_num_t  *face_som_idx,
                                           PDM_l_num_t  *face_som_nb,
                                           PDM_l_num_t  *face_som,
                                           PDM_l_num_t  *cell_face_idx,
                                           PDM_l_num_t  *cell_face_nb,
                                           PDM_l_num_t  *cell_face,
                                           PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_faces_facesom_add(PDM_writer_t *cs,
                                         int     id_geom,
                                         int     id_part,
                                         int     n_face,
                                         PDM_l_num_t  *face_som_idx,
                                         PDM_l_num_t  *face_som_nb,
                                         PDM_l_num_t  *face_som,
                                         PDM_g_num_t  *numabs)

  
  void PDM_writer_geom_write(PDM_writer_t *cs, int id_geom)

  
  void PDM_writer_geom_data_free(PDM_writer_t *cs, int id_geom)
  

  void PDM_writer_geom_free(PDM_writer_t *cs, int id_geom)

  
  int
  PDM_writer_var_create
  (
   PDM_writer_t               *cs,
    PDM_writer_status_t   st_dep_tps,
    PDM_writer_var_dim_t  dim,
    PDM_writer_var_loc_t  loc,
    char                 *nom_var
  )
  
  
  
  void
  PDM_writer_name_map_add
  (
   PDM_writer_t *cs,
    char   *public_name,
    char   *private_name
  )
  
  
  
  void
  PDM_writer_var_write
  (
   PDM_writer_t *cs,
    int     id_var
  )
  
  
  void
  PDM_writer_var_set
  (
   PDM_writer_t     *cs,
    int         id_var,
    int         id_geom,
    int         id_part,
    PDM_real_t *val
  )
  
  
  void
  PDM_writer_var_data_free
  (
   PDM_writer_t *cs,
    int     id_var
  )
  
  
  void
  PDM_writer_var_free
  (
   PDM_writer_t *cs,
    int     id_var
  )
  
   */
  
  void
  PDM_writer_fmt_add
  (
    char                  *name,            
    PDM_writer_fct_t       create_fct,      
    PDM_writer_fct_t       free_fct,        
    PDM_writer_fct_t       beg_step_fct,    
    PDM_writer_fct_t       end_step_fct,    
    PDM_writer_geom_fct_t  geom_create_fct, 
    PDM_writer_geom_fct_t  geom_write_fct,  
    PDM_writer_geom_fct_t  geom_free_fct,   
    PDM_writer_var_fct_t   var_create_fct,  
    PDM_writer_var_fct_t   var_write_fct,   
    PDM_writer_var_fct_t   var_free_fct     
  )
  
  
  void
  PDM_writer_fmt_free
  (
   void
  )
  
  
  void
  PDM_writer_geom_data_reset
  (
   PDM_writer_t *cs,
    int     id_geom
   )
  