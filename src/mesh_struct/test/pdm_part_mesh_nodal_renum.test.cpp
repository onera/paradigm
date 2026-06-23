
#include <vector>
#include <stddef.h>
#include "pdm_doctest.h"
#include "doctest/doctest.h"
#include "doctest/extensions/doctest_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_renum.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_part_mesh_nodal_algorithm.h"
#include "pdm_generate_mesh.h"
#include "pdm_priv.h"
#include "pdm_logging.h"



MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_renum", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  // PDM_plugin_load();

  // Pour le part_mesh_ndal :
  // - Si la renumératation change
  int n_vtx_seg = 40;
  int n_part    = 1;
  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_rectangle(pdm_comm,
                                                           PDM_MESH_NODAL_TRIA3,
                                                           1,
                                                           NULL,
                                                           -1,
                                                           -1,
                                                           0,
                                                           2,
                                                           2,
                                                           n_vtx_seg,
                                                           n_vtx_seg,
                                                           n_part,
                                                           PDM_SPLIT_DUAL_WITH_PTSCOTCH);

  PDM_part_mesh_nodal_dump_vtk(pmn,
                               PDM_GEOMETRY_KIND_SURFACIC,
                               "pmn_in");


  // 2 chemins possibles :
  //   - On fait le tri sur un part_mesh classique et on applique le tri aux part_mesh_nodal
  //   - Tri "dédié" aux nodal :
  //      + On calcul le cell_vtx + vtx_cell et on agence pour les élement sont proches (bfs par dimension)
  //      + On reordonne les sommets pour soi également proche les uns des autres
  //        => On parcours la connectivité cell_vtx, on renumérote les sommets par ordre d'apparition

  const int renum_properties_face[9] = {100, 0, 1, 64, 2, 1, 1, 1, 1};
  // Ou on stocke le old_to_new ?
  const char *renum_method[4] = {"PDM_PART_RENUM_VTX_NONE",
                                 "PDM_PART_RENUM_EDGE_NONE",
                                 "PDM_PART_RENUM_FACE_HPC",
                                 "PDM_PART_RENUM_CELL_NONE"
                                 };
  const int *renum_properties[4] = {NULL,
                                    NULL,
                                    renum_properties_face,
                                    NULL};
  PDM_part_mesh_nodal_renum(pmn,
                            renum_method,
                            renum_properties);

  PDM_part_mesh_nodal_dump_vtk(pmn,
                               PDM_GEOMETRY_KIND_SURFACIC,
                               "pmn_out");

  PDM_part_mesh_nodal_free(pmn);
}



MPI_TEST_CASE("[pdm_part_mesh_nodal] - PDM_part_mesh_nodal_renum - mixed ", 1) {

  int i_rank = -1;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  PDM_MPI_Comm_rank(pdm_comm, &i_rank);

  PDM_part_mesh_nodal_t* pmn = PDM_generate_mesh_nodal_from_file(pdm_comm,
                                                                 1,
                                                                 PDM_SPLIT_DUAL_WITH_HILBERT,
                                                                 PDM_MESH_DIR"mixed_elements_2d.mesh");

  PDM_part_mesh_nodal_dump_vtk(pmn,
                               PDM_GEOMETRY_KIND_SURFACIC,
                               "pmn_in");


  // 2 chemins possibles :
  //   - On fait le tri sur un part_mesh classique et on applique le tri aux part_mesh_nodal
  //   - Tri "dédié" aux nodal :
  //      + On calcul le cell_vtx + vtx_cell et on agence pour les élement sont proches (bfs par dimension)
  //      + On reordonne les sommets pour soi également proche les uns des autres
  //        => On parcours la connectivité cell_vtx, on renumérote les sommets par ordre d'apparition

  const int renum_properties_face[9] = {100, 0, 1, 64, 2, 1, 1, 1, 1};
  // Ou on stocke le old_to_new ?
  const char *renum_method[4] = {"PDM_PART_RENUM_VTX_NONE",
                                 "PDM_PART_RENUM_EDGE_NONE",
                                 "PDM_PART_RENUM_FACE_HPC",
                                 "PDM_PART_RENUM_CELL_NONE"
                                 };
  const int *renum_properties[4] = {NULL,
                                    NULL,
                                    renum_properties_face,
                                    NULL};
  PDM_part_mesh_nodal_renum(pmn,
                            renum_method,
                            renum_properties);

  PDM_part_mesh_nodal_dump_vtk(pmn,
                               PDM_GEOMETRY_KIND_SURFACIC,
                               "pmn_out");


  PDM_part_mesh_nodal_free(pmn);

}
