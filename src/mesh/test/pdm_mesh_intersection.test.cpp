#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_generate_mesh.h"
#include "pdm_logging.h"
#include "pdm_mesh_intersection.h"
#include "pdm_priv.h"


MPI_TEST_CASE("PDM_mesh_intersection - surface / surface - full overlap", 2) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  /* Generate 1st mesh */
  PDM_part_mesh_nodal_t *mesh_a = PDM_generate_mesh_rectangle(comm,
                                                              PDM_MESH_NODAL_QUAD4,
                                                              1,
                                                              NULL,
                                                              0.,
                                                              0.,
                                                              0.,
                                                              1.,
                                                              1.,
                                                              8,
                                                              14,
                                                              1,
                                                              PDM_SPLIT_DUAL_WITH_HILBERT);

  /* Generate 2nd mesh */
  int          *pn_vtx         = NULL;
  int          *pn_edge        = NULL;
  int          *pn_face        = NULL;
  double      **pvtx_coord     = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  int         **pface_vtx      = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_generate_mesh_rectangle_ngon(comm,
                                   PDM_MESH_NODAL_TRIA3,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   11,
                                   11,
                                   1,
                                   PDM_SPLIT_DUAL_WITH_IMPLICIT,
                                   0.,
                                   &pn_vtx,
                                   &pn_edge,
                                   &pn_face,
                                   &pvtx_coord,
                                   &pedge_vtx,
                                   &pface_edge_idx,
                                   &pface_edge,
                                   &pface_vtx,
                                   &pvtx_ln_to_gn,
                                   &pedge_ln_to_gn,
                                   &pface_ln_to_gn);



  /* Compute intersection */
  PDM_mesh_intersection_t *mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_WEIGHT,
                                                             2,
                                                             2,
                                                             0.,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  PDM_mesh_intersection_mesh_nodal_set(mi, 0, mesh_a);

  PDM_mesh_intersection_n_part_set(mi, 1, 1);

  SUBCASE("With edges") {
    PDM_mesh_intersection_part_set(mi,
                                   1,
                                   0,
                                   0,
                                   pn_face[0],
                                   pn_edge[0],
                                   pn_vtx [0],
                                   NULL,
                                   NULL,
                                   pface_edge_idx[0],
                                   pface_edge    [0],
                                   pedge_vtx     [0],
                                   NULL,
                                   NULL,
                                   NULL,
                                   pface_ln_to_gn[0],
                                   pedge_ln_to_gn[0],
                                   pvtx_ln_to_gn [0],
                                   pvtx_coord    [0]);
  }
  SUBCASE("Without edges") {
    PDM_mesh_intersection_part_set(mi,
                                   1,
                                   0,
                                   0,
                                   pn_face[0],
                                   0,
                                   pn_vtx [0],
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   NULL,
                                   pface_edge_idx[0],
                                   pface_vtx     [0],
                                   NULL,
                                   pface_ln_to_gn[0],
                                   NULL,
                                   pvtx_ln_to_gn [0],
                                   pvtx_coord    [0]);
  }


  PDM_mesh_intersection_compute(mi);

  /* Check result */
  int         *a_b_idx;
  PDM_g_num_t *a_b;
  double      *a_b_area;
  PDM_mesh_intersection_result_from_a_get(mi,
                                          0,
                                          &a_b_idx,
                                          &a_b,
                                          &a_b_area);

  int n_elt_a = PDM_part_mesh_nodal_n_elmts_get(mesh_a,
                                                PDM_GEOMETRY_KIND_SURFACIC,
                                                0);

  const double a_area = (1./7) * (1./13);

  for (int i_elt_a = 0; i_elt_a < n_elt_a; i_elt_a++) {
    double sum_area = 0;
    for (int i = a_b_idx[i_elt_a]; i < a_b_idx[i_elt_a+1]; i++) {
      sum_area += a_b_area[i];
    }
    CHECK(abs(sum_area - a_area) < 1e-12);
  }


  PDM_part_to_part_t *ptp = NULL;
  PDM_mesh_intersection_part_to_part_get(mi,
                                         &ptp,
                                         PDM_OWNERSHIP_KEEP);

  int  n_ref_b = 0;
  int *ref_b   = NULL;
  PDM_part_to_part_ref_lnum2_single_part_get(ptp,
                                             0,
                                             &n_ref_b,
                                             &ref_b);

  int         *b_a_idx = NULL;
  PDM_g_num_t *b_a     = NULL;
  PDM_part_to_part_gnum1_come_from_single_part_get(ptp,
                                                   0,
                                                   &b_a_idx,
                                                   &b_a);

  double *b_a_area;
  PDM_mesh_intersection_result_from_b_get(mi,
                                          0,
                                          &b_a_area);

  CHECK(n_ref_b == pn_face[0]);

  const double b_area = 0.5 * (1./10) * (1./10);

  for (int i_ref_b = 0; i_ref_b < n_ref_b; i_ref_b++) {
    double sum_area = 0;
    for (int i = b_a_idx[i_ref_b]; i < b_a_idx[i_ref_b+1]; i++) {
      sum_area += b_a_area[i];
    }
    CHECK(abs(sum_area - b_area) < 1e-12);
  }


  /* Free memory */
  PDM_mesh_intersection_free(mi);
  PDM_part_mesh_nodal_free(mesh_a);
  PDM_free(pvtx_coord    [0]);
  PDM_free(pedge_vtx     [0]);
  PDM_free(pface_edge_idx[0]);
  PDM_free(pface_edge    [0]);
  PDM_free(pface_vtx     [0]);
  PDM_free(pvtx_ln_to_gn [0]);
  PDM_free(pedge_ln_to_gn[0]);
  PDM_free(pface_ln_to_gn[0]);
  PDM_free(pn_vtx        );
  PDM_free(pn_edge       );
  PDM_free(pn_face       );
  PDM_free(pvtx_coord    );
  PDM_free(pedge_vtx     );
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge    );
  PDM_free(pface_vtx     );
  PDM_free(pvtx_ln_to_gn );
  PDM_free(pedge_ln_to_gn);
  PDM_free(pface_ln_to_gn);
}



MPI_TEST_CASE("PDM_mesh_intersection - surface / surface - partial overlap", 2) {

  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_Mesh_nodal_elt_t elt_type_a = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t elt_type_b = PDM_MESH_NODAL_TRIA3;

  SUBCASE("TRIA3 - TRIA3") {
    elt_type_a = PDM_MESH_NODAL_TRIA3;
    elt_type_b = PDM_MESH_NODAL_TRIA3;
  }
  SUBCASE("TRIA3 - QUAD4") {
    elt_type_a = PDM_MESH_NODAL_TRIA3;
    elt_type_b = PDM_MESH_NODAL_QUAD4;
  }
  SUBCASE("TRIA3 - POLY_2D") {
    elt_type_a = PDM_MESH_NODAL_TRIA3;
    elt_type_b = PDM_MESH_NODAL_POLY_2D;
  }
  SUBCASE("QUAD4 - TRIA3") {
    elt_type_a = PDM_MESH_NODAL_QUAD4;
    elt_type_b = PDM_MESH_NODAL_TRIA3;
  }
  SUBCASE("QUAD4 - QUAD4") {
    elt_type_a = PDM_MESH_NODAL_QUAD4;
    elt_type_b = PDM_MESH_NODAL_QUAD4;
  }
  SUBCASE("QUAD4 - POLY_2D") {
    elt_type_a = PDM_MESH_NODAL_QUAD4;
    elt_type_b = PDM_MESH_NODAL_POLY_2D;
  }
  SUBCASE("POLY_2D - TRIA3") {
    elt_type_a = PDM_MESH_NODAL_POLY_2D;
    elt_type_b = PDM_MESH_NODAL_TRIA3;
  }
  SUBCASE("POLY_2D - QUAD4") {
    elt_type_a = PDM_MESH_NODAL_POLY_2D;
    elt_type_b = PDM_MESH_NODAL_QUAD4;
  }
  SUBCASE("POLY_2D - POLY_2D") {
    elt_type_a = PDM_MESH_NODAL_POLY_2D;
    elt_type_b = PDM_MESH_NODAL_POLY_2D;
  }

  /* Generate 1st mesh */
  int          *pn_vtx         = NULL;
  int          *pn_edge        = NULL;
  int          *pn_face        = NULL;
  double      **pvtx_coord     = NULL;
  int         **pedge_vtx      = NULL;
  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;
  int         **pface_vtx      = NULL;
  PDM_g_num_t **pvtx_ln_to_gn  = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;
  PDM_generate_mesh_rectangle_ngon(comm,
                                   elt_type_a,
                                   0.,
                                   0.,
                                   0.,
                                   1.,
                                   1.,
                                   13,
                                   13,
                                   1,
                                   PDM_SPLIT_DUAL_WITH_IMPLICIT,
                                   0.,
                                   &pn_vtx,
                                   &pn_edge,
                                   &pn_face,
                                   &pvtx_coord,
                                   &pedge_vtx,
                                   &pface_edge_idx,
                                   &pface_edge,
                                   &pface_vtx,
                                   &pvtx_ln_to_gn,
                                   &pedge_ln_to_gn,
                                   &pface_ln_to_gn);

  /* Generate 2nd mesh */
  PDM_part_mesh_nodal_t *mesh_b = PDM_generate_mesh_rectangle(comm,
                                                              elt_type_b,
                                                              1,
                                                              NULL,
                                                              0.,
                                                              0.,
                                                              0.,
                                                              1.,
                                                              1.,
                                                              17,
                                                              17,
                                                              1,
                                                              PDM_SPLIT_DUAL_WITH_HILBERT);

  int     n_vtx_b     = PDM_part_mesh_nodal_n_vtx_get(mesh_b, 0)                                                              ;
  double *vtx_coord_b = PDM_part_mesh_nodal_vtx_coord_get(mesh_b, 0, PDM_OWNERSHIP_BAD_VALUE);

  const double angle = atan2(1., 2.);
  const double c     = cos(angle);
  const double s     = sin(angle);

  for (int i_vtx_b = 0; i_vtx_b < n_vtx_b; i_vtx_b++) {
    // Rotate
    double x = vtx_coord_b[3*i_vtx_b  ];
    double y = vtx_coord_b[3*i_vtx_b+1];
    vtx_coord_b[3*i_vtx_b  ] = c*x - s*y;
    vtx_coord_b[3*i_vtx_b+1] = s*x + c*y;

    // Offset
    vtx_coord_b[3*i_vtx_b  ] += 0.5;
    vtx_coord_b[3*i_vtx_b+1] += 0.5;
  }


  /* Compute intersection */
  PDM_mesh_intersection_t *mi = PDM_mesh_intersection_create(PDM_MESH_INTERSECTION_KIND_WEIGHT,
                                                             2,
                                                             2,
                                                             0.,
                                                             comm,
                                                             PDM_OWNERSHIP_KEEP);

  PDM_mesh_intersection_n_part_set(mi, 0, 1);

  PDM_mesh_intersection_part_set(mi,
                                 0,
                                 0,
                                 0,
                                 pn_face[0],
                                 0,
                                 pn_vtx [0],
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 NULL,
                                 pface_edge_idx[0],
                                 pface_vtx     [0],
                                 NULL,
                                 pface_ln_to_gn[0],
                                 NULL,
                                 pvtx_ln_to_gn [0],
                                 pvtx_coord    [0]);

  PDM_mesh_intersection_mesh_nodal_set(mi, 1, mesh_b);


  PDM_mesh_intersection_compute(mi);

  /* Check result */
  int         *a_b_idx;
  PDM_g_num_t *a_b;
  double      *a_b_area;
  PDM_mesh_intersection_result_from_a_get(mi,
                                          0,
                                          &a_b_idx,
                                          &a_b,
                                          &a_b_area);

  double l_total = 0.;
  for (int i = a_b_idx[0]; i < a_b_idx[pn_face[0]]; i++) {
    l_total += a_b_area[i];
  }

  double g_total;
  PDM_MPI_Reduce(&l_total, &g_total, 1, PDM_MPI_DOUBLE, PDM_MPI_SUM, 0, comm);

  if (test_rank == 0) {
    CHECK(abs(g_total - 0.25) < 1e-14);
  }


  /* Free memory */
  PDM_mesh_intersection_free(mi);
  PDM_free(pvtx_coord    [0]);
  PDM_free(pedge_vtx     [0]);
  PDM_free(pface_edge_idx[0]);
  PDM_free(pface_edge    [0]);
  PDM_free(pface_vtx     [0]);
  PDM_free(pvtx_ln_to_gn [0]);
  PDM_free(pedge_ln_to_gn[0]);
  PDM_free(pface_ln_to_gn[0]);
  PDM_free(pn_vtx        );
  PDM_free(pn_edge       );
  PDM_free(pn_face       );
  PDM_free(pvtx_coord    );
  PDM_free(pedge_vtx     );
  PDM_free(pface_edge_idx);
  PDM_free(pface_edge    );
  PDM_free(pface_vtx     );
  PDM_free(pvtx_ln_to_gn );
  PDM_free(pedge_ln_to_gn);
  PDM_free(pface_ln_to_gn);
  PDM_part_mesh_nodal_free(mesh_b);
}
