#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_part_connectivity_transform.h"

static PDM_part_mesh_nodal_t *
_generate_mesh
(
 PDM_MPI_Comm         comm,
 PDM_Mesh_nodal_elt_t elt_type
 )
{
  int              n_part      = 1;
  PDM_g_num_t      n_vtx_seg   = 4;
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;

  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(elt_type);

  if (mesh_dimension == 2) {
    return PDM_generate_mesh_rectangle(comm,
                                       elt_type,
                                       1,
                                       NULL,
                                       0.,
                                       0.,
                                       0.,
                                       1.,
                                       1.,
                                       n_vtx_seg,
                                       n_vtx_seg,
                                       n_part,
                                       part_method);
  }
  else {
    return PDM_generate_mesh_parallelepiped(comm,
                                            elt_type,
                                            1,
                                            NULL,
                                            0.,
                                            0.,
                                            0.,
                                            1.,
                                            1.,
                                            1.,
                                            n_vtx_seg,
                                            n_vtx_seg,
                                            n_vtx_seg,
                                            n_part,
                                            part_method);
  }
}

MPI_TEST_CASE("pdm_part_mesh_nodal_to_part_mesh", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_Mesh_nodal_elt_t elt_type                 = PDM_MESH_NODAL_N_ELEMENT_TYPES;
  PDM_bool_t           keep_link_elmt_to_entity = PDM_FALSE;

  PDM_bool_t enable_connectivity[PDM_CONNECTIVITY_TYPE_MAX];
  PDM_bool_t enable_gnum        [PDM_MESH_ENTITY_MAX];
  PDM_bool_t enable_groups      [PDM_BOUND_TYPE_MAX];

  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    enable_connectivity[i] = PDM_FALSE;
  }

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    enable_gnum[i] = PDM_FALSE;
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    enable_groups[i] = PDM_FALSE;
  }


  /* Subcases */
  SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("TRIA3, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("TRIA3, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("TRIA3, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TRIA3, FACE_VTX, VTX, keep_link = 0") {
      printf("TRIA3, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_VTX, VTX, keep_link = 1") {
      printf("TRIA3, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TRIA3, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_TRIA3;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("TRIA3, FACE_VTX, keep_link = 0") {
      printf("TRIA3, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TRIA3, FACE_VTX, keep_link = 1") {
      printf("TRIA3, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("QUAD4, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("QUAD4, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("QUAD4, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("QUAD4, FACE_VTX, VTX, keep_link = 0") {
      printf("QUAD4, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_VTX, VTX, keep_link = 1") {
      printf("QUAD4, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_QUAD4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("QUAD4, FACE_VTX, keep_link = 0") {
      printf("QUAD4, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("QUAD4, FACE_VTX, keep_link = 1") {
      printf("QUAD4, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("POLY_2D, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_VTX, VTX, keep_link = 0") {
      printf("POLY_2D, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_VTX, VTX, keep_link = 1") {
      printf("POLY_2D, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_POLY_2D;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("POLY_2D, FACE_VTX, keep_link = 0") {
      printf("POLY_2D, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("POLY_2D, FACE_VTX, keep_link = 1") {
      printf("POLY_2D, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4, CELL_FACE, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_TETRA4;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, keep_link = 0") {
      printf("TETRA4, CELL_FACE, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("TETRA4, CELL_FACE, FACE_VTX, keep_link = 1") {
      printf("TETRA4, CELL_FACE, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, keep_link = 0") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PYRAMID5, CELL_FACE, FACE_VTX, keep_link = 1") {
      printf("PYRAMID5, CELL_FACE, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6, CELL_FACE, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_PRISM6;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, keep_link = 0") {
      printf("PRISM6, CELL_FACE, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("PRISM6, CELL_FACE, FACE_VTX, keep_link = 1") {
      printf("PRISM6, CELL_FACE, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_EDGE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, EDGE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_EDGE, EDGE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_VTX, FACE, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_FACE]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_VTX, FACE, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_VTX, CELL, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_CELL]            = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_VTX, CELL, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_VTX, VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    enable_gnum        [PDM_MESH_ENTITY_VTX ]            = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_VTX, VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_VTX, VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8, CELL_FACE, FACE_VTX") {
    elt_type = PDM_MESH_NODAL_HEXA8;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_CELL_FACE] = PDM_TRUE;
    enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_VTX ] = PDM_TRUE;
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, keep_link = 0") {
      printf("HEXA8, CELL_FACE, FACE_VTX, keep_link = 0\n");
      keep_link_elmt_to_entity = PDM_FALSE;
    };
    SUBCASE("HEXA8, CELL_FACE, FACE_VTX, keep_link = 1") {
      printf("HEXA8, CELL_FACE, FACE_VTX, keep_link = 1\n");
      keep_link_elmt_to_entity = PDM_TRUE;
    };
  }



  /* Test */
  PDM_part_mesh_nodal_t *pmesh_nodal = _generate_mesh(pdm_comm, elt_type);

  PDM_part_mesh_nodal_to_part_mesh_t *pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmesh_nodal,
                                                                                          keep_link_elmt_to_entity,
                                                                                          PDM_OWNERSHIP_USER);

  for (int i = 0; i < PDM_CONNECTIVITY_TYPE_MAX; i++) {
    if (enable_connectivity[i]) {
      PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                           (PDM_connectivity_type_t) i);
    }
  }

  for (int i = 0; i < PDM_MESH_ENTITY_MAX; i++) {
    if (enable_gnum[i]) {
      PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm,
                                                     (PDM_mesh_entities_t) i);
    }
  }

  for (int i = 0; i < PDM_BOUND_TYPE_MAX; i++) {
    if (enable_groups[i]) {
      PDM_part_mesh_nodal_to_part_mesh_groups_enable(pmn_to_pm,
                                                     (PDM_bound_type_t) i);
    }
  }

  PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

  PDM_part_mesh_t *pmesh = NULL;
  PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                 &pmesh,
                                                 PDM_OWNERSHIP_KEEP);

  if (enable_connectivity[PDM_CONNECTIVITY_TYPE_FACE_EDGE] == PDM_TRUE &&
      enable_connectivity[PDM_CONNECTIVITY_TYPE_EDGE_VTX ] == PDM_TRUE) {
    // Check face->edge + edge->vtx => face->vtx
    int n_face = PDM_part_mesh_n_entity_get(pmesh,
                                            0,
                                            PDM_MESH_ENTITY_FACE);
    int *face_edge_idx = NULL;
    int *face_edge     = NULL;
    PDM_part_mesh_connectivity_get(pmesh,
                                   0,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &face_edge,
                                   &face_edge_idx,
                                   PDM_OWNERSHIP_KEEP);

    int *edge_vtx_idx = NULL;
    int *edge_vtx     = NULL;
    PDM_part_mesh_connectivity_get(pmesh,
                                   0,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &edge_vtx,
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_KEEP);


    int *face_vtx = NULL;
    PDM_compute_face_vtx_from_face_and_edge(n_face,
                                            face_edge_idx,
                                            face_edge,
                                            edge_vtx,
                                            &face_vtx);
    PDM_free(face_vtx);
  }

  PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);

  PDM_part_mesh_nodal_free(pmesh_nodal);
}


// The subcases were generated using the following Python code :

/*
count = 0

dim_elt_types = {
  2: [
  "TRIA3",
  "QUAD4",
  "POLY_2D"
  ],
  3: [
  "TETRA4",
  "PYRAMID5",
  "PRISM6",
  "HEXA8"
  ]
}

dim_entities = {
  0: "VTX",
  1: "EDGE",
  2: "FACE",
  3: "CELL"
}

dim_connectivities = {
  3: [
    ["CELL_FACE", "FACE_EDGE", "EDGE_VTX"],
    ["CELL_FACE", "FACE_VTX"]
  ],
  2: [
    ["FACE_EDGE", "EDGE_VTX"],
    ["FACE_VTX"],
  ]
}


def partitions(a, current=[], p=[]):
  if len(a) == 0:
    p.append(current)
  else:
    p = partitions(a[1:], current + [a[0]], p)
    p = partitions(a[1:], current,          p)
  return p

pdm_bool = ["PDM_FALSE", "PDM_TRUE"]


def gen_subcase(subcase, name=None):
  global count
  count += 1

  if name is None:
    name = str(count)

  s = ''
  s += '  SUBCASE("%s") {\n' % name

  s += '    elt_type = PDM_MESH_NODAL_%s;\n' % subcase["elt_type"]
  for c in subcase["connectivity"]:
    s += '    enable_connectivity[PDM_CONNECTIVITY_TYPE_%s%s] = PDM_TRUE;\n' % (c, (9-len(c))*" ")
  for e in subcase["gnum"]:
    s += '    enable_gnum        [PDM_MESH_ENTITY_%s%s]            = PDM_TRUE;\n' % (e, (4-len(e))*" ")

  for keep_link in range(2):
    _name = '%s, keep_link = %d' % (name, keep_link)
    s += '    SUBCASE("%s") {\n' % _name
    s += '      printf("%s\\n");\n' % _name
    s += '      keep_link_elmt_to_entity = %s;\n' % pdm_bool[keep_link]
    s += '    };\n'

    count += 1

  s += '  }\n'

  return s




text = ""

n_vtx_seg = 3

for dim in dim_elt_types.keys():

  for elt_type in dim_elt_types[dim]:

    for connectivities in dim_connectivities[dim]:

      # get list of available entities in current scenario (from built connectivities)
      entities = list(set([e for c in connectivities for e in c.split("_")]))

      # generate partitions of the set of available entities
      p = partitions(entities, [], [])

      for gnum in p:
        # eliminate subcases that require gnums for some entities but not for vertices
        if len(gnum) > 0 and "VTX" not in gnum:
          continue

        subcase = dict()
        subcase["elt_type"]     = elt_type
        subcase["connectivity"] = connectivities
        subcase["gnum"]         = gnum

        name = ", ".join([elt_type] + connectivities + gnum)

        text += gen_subcase(subcase, name)
        text += "\n"


with open("pmn_to_pm_unit_test_subcases.cpp", "w") as out:
  out.write(text)

print(f"generated {count} subcases")
*/
