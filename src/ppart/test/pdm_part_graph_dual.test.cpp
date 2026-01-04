#include <stddef.h>
#include <vector>

#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_doctest.h"
#include "pdm_logging.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_graph_dual.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_part_mesh_nodal_to_part_mesh.h"
#include "pdm_vtk.h"


static
PDM_part_mesh_nodal_t*
_generate_mesh
(
  PDM_MPI_Comm    pdm_comm,
  int             n_vtx_seg
)
{
  int              n_part       = 1;
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;

  /* Warmup */
  PDM_dcube_nodal_t* dcube = PDM_dcube_nodal_gen_create(pdm_comm,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        n_vtx_seg,
                                                        2.,
                                                        -1,
                                                        -1,
                                                        -1,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        PDM_OWNERSHIP_KEEP);
  PDM_dcube_nodal_gen_build (dcube);
  PDM_dmesh_nodal_t* dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  int n_domain = 1;
  int n_part_domains = n_part;
  PDM_multipart_t *mpart = PDM_multipart_create(n_domain,
                                                &n_part_domains,
                                                PDM_FALSE,
                                                split_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                pdm_comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);

  if(0 == 1) {
    PDM_dmesh_nodal_dump_vtk(dmn, PDM_GEOMETRY_KIND_SURFACIC, "out_surfacic");
  }
  PDM_multipart_compute(mpart);

  PDM_part_mesh_nodal_t* pmesh_nodal = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart, 0, &pmesh_nodal, PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_complete_part_comm_graph(pmesh_nodal);

  PDM_multipart_free(mpart);
  PDM_dcube_nodal_gen_free(dcube);
  return pmesh_nodal;
}

static
PDM_part_mesh_t*
_part_mesh_nodal_to_part_mesh
(
  int                    dim,
  PDM_part_mesh_nodal_t *pmn
)
{
  PDM_part_mesh_nodal_to_part_mesh_t* pmn_to_pm = PDM_part_mesh_nodal_to_part_mesh_create(pmn,
                                                                                          PDM_FALSE,
                                                                                          PDM_OWNERSHIP_USER);

  PDM_part_mesh_t *pm = NULL;
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_nodal_to_part_mesh_g_nums_enable(pmn_to_pm, PDM_MESH_ENTITY_EDGE);
  if(dim == 3) {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX);
  } else {
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE);
    PDM_part_mesh_nodal_to_part_mesh_connectivity_enable(pmn_to_pm,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX);
  }

  PDM_part_mesh_nodal_to_part_mesh_compute(pmn_to_pm);

  PDM_part_mesh_nodal_to_part_mesh_part_mesh_get(pmn_to_pm,
                                                 &pm,
                                                 PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_to_part_mesh_free(pmn_to_pm);

  // Finalize part_mesh
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_VTX );
  PDM_part_mesh_part_comm_graph_compute_from_gnum(pm, PDM_MESH_ENTITY_EDGE);

  return pm;
}


MPI_TEST_CASE("[PDM_part_assembly_dual_graph] Full ", 2) {

  int n_vtx_a = 5;
  PDM_part_mesh_nodal_t* pmn = _generate_mesh(test_comm, n_vtx_a);

  PDM_part_mesh_t* pm = _part_mesh_nodal_to_part_mesh(2, pmn);




  PDM_part_mesh_nodal_free(pmn);
  PDM_part_mesh_free(pm);

}


MPI_TEST_CASE("[PDM_part_assembly_dual_graph] select ", 2) {


}

// Same but with pcg_node == NULL (cell centered)

MPI_TEST_CASE("[PDM_transfer_entity1_part_id_to_entity2_part_id] ", 2) {


}

