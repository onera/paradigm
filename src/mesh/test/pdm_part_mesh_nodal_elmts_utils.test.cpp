#include <vector>
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_doctest.h"
#include "pdm_part_mesh_nodal_elmts_utils.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_logging.h"
#include "pdm_error.h"

static PDM_part_mesh_nodal_t *
_generate_mesh
(
 PDM_MPI_Comm         comm,
 PDM_Mesh_nodal_elt_t elt_type
 )
{
  int              n_part      = 1;
  PDM_g_num_t      n_vtx_seg   = 5;
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


MPI_TEST_CASE("[pdm_part_mesh_nodal_elmts_utils] - part_mesh_nodal_elmts_compute_child_entities", 1) {

  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);

  PDM_Mesh_nodal_elt_t elt_type        = PDM_MESH_NODAL_N_ELEMENT_TYPES;
  PDM_bool_t           only_child_link = PDM_FALSE;


  SUBCASE("TRIA3") {
    printf("TRIA3\n");
    elt_type = PDM_MESH_NODAL_TRIA3;
    SUBCASE("TRIA3, only child link") {
      printf("TRIA3, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("QUAD4") {
    printf("QUAD4\n");
    elt_type = PDM_MESH_NODAL_QUAD4;
    SUBCASE("QUAD4, only child link") {
      printf("QUAD4, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("POLY_2D") {
    printf("POLY_2D\n");
    elt_type = PDM_MESH_NODAL_POLY_2D;
    SUBCASE("POLY_2D, only child link") {
      printf("POLY_2D, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("TETRA4") {
    printf("TETRA4\n");
    elt_type = PDM_MESH_NODAL_TETRA4;
    SUBCASE("TETRA4, only child link") {
      printf("TETRA4, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("PYRAMID5") {
    printf("PYRAMID5\n");
    elt_type = PDM_MESH_NODAL_PYRAMID5;
    SUBCASE("PYRAMID5, only child link") {
      printf("PYRAMID5, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("PRISM6") {
    printf("PRISM6\n");
    elt_type = PDM_MESH_NODAL_PRISM6;
    SUBCASE("PRISM6, only child link") {
      printf("PRISM6, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }

  SUBCASE("HEXA8") {
    printf("HEXA8\n");
    elt_type = PDM_MESH_NODAL_HEXA8;
    SUBCASE("HEXA8, only child link") {
      printf("HEXA8, only child link\n");
      only_child_link = PDM_TRUE;
    };
  }




  /* Generate mesh */
  PDM_part_mesh_nodal_t *pmn = _generate_mesh(pdm_comm, elt_type);


  /* Compute link between child and parent entities */
  PDM_geometry_kind_t geom_kind_parent = PDM_GEOMETRY_KIND_MAX;
  int mesh_dimension = PDM_Mesh_nodal_elt_dim_get(elt_type);
  if (mesh_dimension == 3) {
    geom_kind_parent = PDM_GEOMETRY_KIND_VOLUMIC;
  }
  else if (mesh_dimension == 2) {
    geom_kind_parent = PDM_GEOMETRY_KIND_SURFACIC;
  }
  else {
    PDM_error(__FILE__, __LINE__, 0, "1D not yet available\n");
  }

  PDM_part_mesh_nodal_elmts_t *pmne_parent = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                           geom_kind_parent);

  for (int geom_kind = (int) geom_kind_parent + 1; geom_kind < (int) PDM_GEOMETRY_KIND_CORNER; geom_kind++) {

    log_trace("geom_kind %d\n", geom_kind);
    PDM_part_mesh_nodal_elmts_t *pmne_child = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                            (PDM_geometry_kind_t) geom_kind);

    int **child_to_parent_idx  = NULL;
    int **child_to_parent      = NULL;
    int  *n_entity             = NULL;
    int **entity_to_vtx_idx    = NULL;
    int **entity_to_vtx        = NULL;
    int **parent_to_entity_idx = NULL;
    int **parent_to_entity     = NULL;

    PDM_part_mesh_nodal_elmts_compute_child_entities(pmne_parent,
                                                     pmne_child,
                                                     only_child_link,
                                                     &child_to_parent_idx,
                                                     &child_to_parent,
                                                     &n_entity,
                                                     &entity_to_vtx_idx,
                                                     &entity_to_vtx,
                                                     &parent_to_entity_idx,
                                                     &parent_to_entity);

    /* Check */
    // TODOUX

    /* Free memory */
    PDM_free(child_to_parent_idx[0]);
    PDM_free(child_to_parent    [0]);
    if (!only_child_link) {
      PDM_free(entity_to_vtx_idx   [0]);
      PDM_free(entity_to_vtx       [0]);
      PDM_free(parent_to_entity_idx[0]);
      PDM_free(parent_to_entity    [0]);
    }

    PDM_free(child_to_parent_idx );
    PDM_free(child_to_parent     );
    PDM_free(n_entity            );
    PDM_free(entity_to_vtx_idx   );
    PDM_free(entity_to_vtx       );
    PDM_free(parent_to_entity_idx);
    PDM_free(parent_to_entity    );
  }

  PDM_part_mesh_nodal_free(pmn);
}

