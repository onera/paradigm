#include "doctest/extensions/doctest_mpi.h"
#include <array>
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_dcube_nodal_gen.h"

MPI_TEST_CASE("[1p] dcube_nodal_gen HEXA ",1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  std::vector<std::vector<PDM_g_num_t>> connec_expected = {
    {1, 2, 11, 10, 4, 5, 14, 13, 2, 3, 12, 11, 5, 6, 15, 14, 4, 5, 14, 13, 7, 8, 17, 16, 5, 6, 15, 14, 8, 9, 18, 17, 10, 11, 20, 19, 13, 14, 23, 22, 11, 12, 21, 20, 14, 15, 24, 23, 13, 14, 23, 22, 16, 17, 26, 25, 14, 15, 24, 23, 17, 18, 27, 26}, // HEXA
    {1, 10, 13, 4, 4, 13, 16, 7, 10, 19, 22, 13, 13, 22, 25, 16}, // QUAD IMIN
    {3, 12, 15, 6, 6, 15, 18, 9, 12, 21, 24, 15, 15, 24, 27, 18}, // QUAD IMAX
    {1, 2, 11, 10, 2, 3, 12, 11, 10, 11, 20, 19, 11, 12, 21, 20}, // QUAD JMIN
    {7, 8, 17, 16, 8, 9, 18, 17, 16, 17, 26, 25, 17, 18, 27, 26}, // QUAD JMAX
    {1, 4, 5, 2, 2, 5, 6, 3, 4, 7, 8, 5, 5, 8, 9, 6}, // QUAD KMIN
    {19, 22, 23, 20, 20, 23, 24, 21, 22, 25, 26, 23, 23, 26, 27, 24}  // QUAD KMAX
  };

  std::vector<std::vector<PDM_g_num_t>> distrib_expected = {
    {0, 8}, // HEXA
    {0, 4}, // QUAD IMIN
    {0, 4}, // QUAD IMAX
    {0, 4}, // QUAD JMIN
    {0, 4}, // QUAD JMAX
    {0, 4}, // QUAD KMIN
    {0, 4}  // QUAD KMAX
  };

  int n_vtx_seg = 3;

  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_init(pdm_comm,
                                                      n_vtx_seg,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      PDM_MESH_NODAL_HEXA8,
                                                      PDM_OWNERSHIP_KEEP);

  PDM_dmesh_nodal_t* dmesh_nodal = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube); /* It will be free by PDM_dcube_nodal_gen_free because PDM_OWNERSHIP_KEEP */

  /* Verif global */
  PDM_g_num_t n_cell_abs = -100;
  PDM_g_num_t n_face_abs = -100;
  PDM_g_num_t n_edge_abs = -100;
  PDM_g_num_t n_vtx_abs  = -100;
  PDM_DMesh_nodal_section_g_dims_get(dmesh_nodal, &n_cell_abs, &n_face_abs, &n_edge_abs, &n_vtx_abs);

  CHECK( n_cell_abs ==  8);
  CHECK( n_face_abs == -1);
  CHECK( n_edge_abs == -1);
  CHECK( n_vtx_abs  == 27);

  int n_section = PDM_DMesh_nodal_n_section_get(dmesh_nodal);

  CHECK( n_section == 7); // HEXA + 6 * QUAD

  int* sections_id = PDM_DMesh_nodal_sections_id_get(dmesh_nodal);

  /* Verification HEXA */
  int id_hexa = sections_id[0];
  PDM_Mesh_nodal_elt_t t_elmt_s1 = PDM_DMesh_nodal_section_type_get(dmesh_nodal, id_hexa);

  CHECK(t_elmt_s1 == PDM_MESH_NODAL_HEXA8);

  // PDM_g_num_t* section_distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, id_hexa);
  // PDM_g_num_t* connect         = PDM_DMesh_nodal_section_std_get(dmesh_nodal, id_hexa);

  // PDM_log_trace_array_long(section_distrib, n_rank , "section_distrib:: ");

  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = sections_id[i_section];
    PDM_Mesh_nodal_elt_t t_elmt = PDM_DMesh_nodal_section_type_get(dmesh_nodal, id_section);

    PDM_g_num_t* distrib = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal, id_section);
    PDM_g_num_t* connect = PDM_DMesh_nodal_section_std_get(dmesh_nodal, id_section);

    int dn_elmt = distrib[i_rank+1] - distrib[i_rank];

    int n_vtx_per_elmt = PDM_Mesh_nodal_n_vertices_element(t_elmt, 1);

    /* Keep it to make easyly ref */
    // PDM_log_trace_array_long(distrib, n_rank+1                 , "distrib:: ");
    // PDM_log_trace_array_long(connect, n_vtx_per_elmt * dn_elmt , "connect:: ");

    CHECK_EQ_C_ARRAY( distrib, distrib_expected[i_section], n_rank+1                );
    CHECK_EQ_C_ARRAY( connect, connec_expected [i_section], n_vtx_per_elmt * dn_elmt);


  }


  PDM_dcube_nodal_gen_free(dcube);

}
