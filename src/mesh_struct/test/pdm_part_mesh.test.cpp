
#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_part_mesh.h"
#include "pdm_mpi.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_multipart.h"
#include "pdm_priv.h"
#include "pdm_logging.h"


static
PDM_part_mesh_t*
_generate_mesh
(
  PDM_MPI_Comm pdm_comm
)
{
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        4,
                                                        4,
                                                        0,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        PDM_MESH_NODAL_QUAD4,
                                                        1,
                                                        PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build(dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);
  // free
  PDM_dcube_nodal_gen_free(dcube);

  int n_part = 1;
  PDM_multipart_t* mpart = PDM_multipart_create(1,
                                                &n_part,
                                                PDM_FALSE,
                                                PDM_SPLIT_DUAL_WITH_HILBERT,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);
  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);
  PDM_DMesh_nodal_free(dmn);

  PDM_part_mesh_t* pm = NULL;
  PDM_multipart_get_part_mesh(mpart, 0, &pm, PDM_OWNERSHIP_USER);

  PDM_multipart_free(mpart);
  return pm;
}


MPI_TEST_CASE("[pdm_part_mesh] Constructor",1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_t* pm = _generate_mesh(pdm_comm);


  PDM_part_mesh_free(pm);
}

MPI_TEST_CASE("[pdm_part_mesh] - PDM_part_mesh_part_comm_graph_compute_from_gnum", 2) {


  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_t* pm = _generate_mesh(pdm_comm);

  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX);
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_FACE);

  PDM_part_comm_graph_t *pcg_vtx = NULL;
  PDM_part_mesh_part_comm_graph_get(pm, PDM_MESH_ENTITY_VTX, &pcg_vtx, PDM_OWNERSHIP_KEEP);

  /* Check vertices */
  int *pvtx_bound = NULL;
  int n_vtx_part_bound = PDM_part_comm_graph_entity_graph_get(pcg_vtx,
                                                              0,
                                                              &pvtx_bound,
                                                              PDM_OWNERSHIP_USER);

  if(1 == 0) {
    PDM_log_trace_array_int(pvtx_bound, 4 * n_vtx_part_bound, "pvtx_bound ::");
  }

  int p0_expected_pvtx_bound[20] = {3, 1, 1, 1, 6, 1, 1, 3, 7, 1, 1, 5, 8, 1, 1, 6, 9, 1, 1, 7};
  int p1_expected_pvtx_bound[20] = {1, 0, 1, 3, 3, 0, 1, 6, 5, 0, 1, 7, 6, 0, 1, 8, 7, 0, 1, 9};

  CHECK(n_vtx_part_bound == 5);

  MPI_CHECK_EQ_C_ARRAY(0, pvtx_bound, p0_expected_pvtx_bound, 20);
  MPI_CHECK_EQ_C_ARRAY(1, pvtx_bound, p1_expected_pvtx_bound, 20);


  /* Check edge */
  PDM_part_comm_graph_t *pcg_edge = NULL;
  PDM_part_mesh_part_comm_graph_get(pm, PDM_MESH_ENTITY_EDGE, &pcg_edge, PDM_OWNERSHIP_USER);

  int *pedge_bound = NULL;
  int n_edge_part_bound = PDM_part_comm_graph_entity_graph_get(pcg_edge,
                                                               0,
                                                               &pedge_bound,
                                                               PDM_OWNERSHIP_KEEP);


  if(0 == 1) {
    PDM_log_trace_array_int(pedge_bound, 4 * n_edge_part_bound, "pedge_bound ::");
  }

  int p0_expected_pedge_bound[16] = {5, 1, 1, 2, 10, 1, 1, 5, 11, 1, 1, 6, 12, 1, 1, 8};
  int p1_expected_pedge_bound[16] = {2, 0, 1, 5, 5, 0, 1, 10, 6, 0, 1, 11, 8, 0, 1, 12};

  CHECK(n_edge_part_bound == 4);

  MPI_CHECK_EQ_C_ARRAY(0, pedge_bound, p0_expected_pedge_bound, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pedge_bound, p1_expected_pedge_bound, 16);


  /* Check face */
  PDM_part_comm_graph_t *pcg_face = NULL;
  PDM_part_mesh_part_comm_graph_get(pm, PDM_MESH_ENTITY_FACE, &pcg_face, PDM_OWNERSHIP_KEEP);

  int *pface_bound = NULL;
  int n_face_part_bound = PDM_part_comm_graph_entity_graph_get(pcg_face,
                                                               0,
                                                               &pface_bound,
                                                               PDM_OWNERSHIP_KEEP);

  CHECK(n_face_part_bound == 0);

  PDM_free(pvtx_bound);
  PDM_part_comm_graph_free(pcg_edge);

  PDM_part_mesh_free(pm);
}


MPI_TEST_CASE("[pdm_part_mesh] - PDM_part_mesh_gnum_compute_from_part_comm_graph", 2) {


  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_t* pm = _generate_mesh(pdm_comm);

  PDM_g_num_t *pvtx_ln_to_gn  = NULL;
  PDM_g_num_t *pedge_ln_to_gn = NULL;
  PDM_g_num_t *pface_ln_to_gn = NULL;

  int n_vtx  = PDM_part_mesh_n_entity_get(pm, 0, PDM_MESH_ENTITY_VTX );
  int n_edge = PDM_part_mesh_n_entity_get(pm, 0, PDM_MESH_ENTITY_EDGE);
  int n_face = PDM_part_mesh_n_entity_get(pm, 0, PDM_MESH_ENTITY_FACE);

  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_VTX , &pvtx_ln_to_gn , PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_EDGE, &pedge_ln_to_gn, PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_FACE, &pface_ln_to_gn, PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_log_trace_array_long(pvtx_ln_to_gn , n_vtx , "pvtx_ln_to_gn  ::");
    PDM_log_trace_array_long(pedge_ln_to_gn, n_edge, "pedge_ln_to_gn ::");
    PDM_log_trace_array_long(pface_ln_to_gn, n_face, "pface_ln_to_gn ::");
  }

  int p0_vtx_ln_to_gn_before [9 ] = {1, 2, 3, 5, 6, 7, 9, 10, 11};
  int p0_edge_ln_to_gn_before[12] = {1, 2, 3, 5, 6, 7, 9, 10, 12, 13, 14, 16};
  int p0_face_ln_to_gn_before[4 ] = {1, 2, 4, 5};

  int p1_vtx_ln_to_gn_before [12] = {3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
  int p1_edge_ln_to_gn_before[16] = {4, 6, 8, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  int p1_face_ln_to_gn_before[5 ] = {3, 6, 7, 8, 9};


  MPI_CHECK_EQ_C_ARRAY(0, pvtx_ln_to_gn , p0_vtx_ln_to_gn_before , 9 );
  MPI_CHECK_EQ_C_ARRAY(0, pedge_ln_to_gn, p0_edge_ln_to_gn_before, 12);
  MPI_CHECK_EQ_C_ARRAY(0, pface_ln_to_gn, p0_face_ln_to_gn_before, 4 );

  MPI_CHECK_EQ_C_ARRAY(1, pvtx_ln_to_gn , p1_vtx_ln_to_gn_before , 12);
  MPI_CHECK_EQ_C_ARRAY(1, pedge_ln_to_gn, p1_edge_ln_to_gn_before, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pface_ln_to_gn, p1_face_ln_to_gn_before, 5 );


  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_FACE);


  PDM_part_mesh_gnum_compute_from_part_comm_graph(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_gnum_compute_from_part_comm_graph(pm, PDM_MESH_ENTITY_EDGE);
  PDM_part_mesh_gnum_compute_from_part_comm_graph(pm, PDM_MESH_ENTITY_FACE);

  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_VTX , &pvtx_ln_to_gn , PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_EDGE, &pedge_ln_to_gn, PDM_OWNERSHIP_KEEP);
  PDM_part_mesh_entity_ln_to_gn_get(pm, 0, PDM_MESH_ENTITY_FACE, &pface_ln_to_gn, PDM_OWNERSHIP_KEEP);

  if(0 == 1) {
    PDM_log_trace_array_long(pvtx_ln_to_gn , n_vtx , "pvtx_ln_to_gn  ::");
    PDM_log_trace_array_long(pedge_ln_to_gn, n_edge, "pedge_ln_to_gn ::");
    PDM_log_trace_array_long(pface_ln_to_gn, n_face, "pface_ln_to_gn ::");
  }

  int p0_vtx_ln_to_gn_after [9 ] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  int p0_edge_ln_to_gn_after[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  int p0_face_ln_to_gn_after[4 ] = {1, 2, 3, 4};

  int p1_vtx_ln_to_gn_after [12] = {3, 10, 6, 11, 7, 8, 9, 12, 13, 14, 15, 16};
  int p1_edge_ln_to_gn_after[16] = {13, 5, 14, 15, 10, 11, 16, 12, 17, 18, 19, 20, 21, 22, 23, 24};
  int p1_face_ln_to_gn_after[5 ] = {5, 6, 7, 8, 9};


  MPI_CHECK_EQ_C_ARRAY(0, pvtx_ln_to_gn , p0_vtx_ln_to_gn_after , 9 );
  MPI_CHECK_EQ_C_ARRAY(0, pedge_ln_to_gn, p0_edge_ln_to_gn_after, 12);
  MPI_CHECK_EQ_C_ARRAY(0, pface_ln_to_gn, p0_face_ln_to_gn_after, 4 );

  MPI_CHECK_EQ_C_ARRAY(1, pvtx_ln_to_gn , p1_vtx_ln_to_gn_after , 12);
  MPI_CHECK_EQ_C_ARRAY(1, pedge_ln_to_gn, p1_edge_ln_to_gn_after, 16);
  MPI_CHECK_EQ_C_ARRAY(1, pface_ln_to_gn, p1_face_ln_to_gn_after, 5 );

  PDM_part_mesh_free(pm);

}
