
cdef extern from "pdm_mesh_location.h":

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # > Wrapping for C
  ctypedef enum PDM_mesh_location_method_t:
    PDM_MESH_LOCATION_OCTREE  = 1
    PDM_MESH_LOCATION_DBBTREE = 2
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  int PDM_mesh_location_create(PDM_mesh_nature_t mesh_nature,
                               int               n_point_cloud,
                               PDM_MPI_Comm      comm);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_n_part_cloud_set(int          id,
                                          int          i_point_cloud,
                                          int          n_part);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_cloud_set(int          id,
                                   int          i_point_cloud,
                                   int          i_part,
                                   int          n_points,
                                   double      *coords,
                                   PDM_g_num_t *gnum);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_shared_nodal_mesh_set(int  id, int  mesh_nodal_id);
  void PDM_mesh_location_mesh_global_data_set (int  id, int  n_part);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set(int          id,
                                  int          i_part,
                                  int          n_cell,
                                  int         *cell_face_idx,
                                  int         *cell_face,
                                  PDM_g_num_t *cell_ln_to_gn,
                                  int          n_face,
                                  int         *face_vtx_idx,
                                  int         *face_vtx,
                                  PDM_g_num_t *face_ln_to_gn,
                                  int          n_vtx,
                                  double      *coords,
                                  PDM_g_num_t *vtx_ln_to_gn);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_part_set_2d(int          id,
                                     int          i_part,
                                     int          n_cell,
                                     int         *cell_edge_idx,
                                     int         *cell_edge,
                                     PDM_g_num_t *cell_ln_to_gn,
                                     int          n_edge,
                                     int         *edge_vtx_idx,
                                     int         *edge_vtx,
                                     PDM_g_num_t *edge_ln_to_gn,
                                     int          n_vtx,
                                     double      *coords,
                                     PDM_g_num_t *vtx_ln_to_gn);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_tolerance_set(int    id, double tol);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_method_set(int id, PDM_mesh_location_method_t method);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_compute(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_get(int           id,
                             int           i_point_cloud,
                             int           i_part,
                             int          *n_points,
                             double      **coord,
                             PDM_g_num_t **g_num,
                             PDM_g_num_t **location,
                             int         **weights_idx,
                             double      **weights);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_free(int id, int partial);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_dump_times(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  void PDM_mesh_location_mesh_nodal_id_get(int id);
  # :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
