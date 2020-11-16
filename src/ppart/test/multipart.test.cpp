#include "doctest/extensions/doctest_mpi.h"

#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "std_e/utils/concatenate.hpp"
#include "std_e/future/span.hpp"
#include "std_e/algorithm/algorithm.hpp"
#include "pdm_dmesh.h"

#include <vector>

struct simple_mesh2_for_test {
/*


              ________________
             /14              /\15
            / |              /|  \
           /_______________ /      \
          /2  |           3/\ |      \
         /                /|  \        \
        /__|__|_________7/ |  | \        \
       6|\             /|\ |      \        \
        |  |  | _ _ _ _ |_ \_ |_ _ _\ _ _ _ _\ 11
        |  \   13     / |    \ 12     \      /
        |  |  /         |  | / \        \   /
        |   _\_ _ _ _/_ | _|_ _ _\_ _ _ _ \/
        |   0           |   1      \      /8
        | /   \     /   | /          \   /
        |______________5| _____________\/9
       4 \     \   /    /             /
          \            /           /
           \    \ /   /         / 
            \        /      /
             \   /  /   /
              \    / /
               \ |//  
                 10
*/
  std::vector<PDM_g_num_t> hex = {
     0, 1, 3, 2, 4, 5, 7, 6,
     13,12,15,14, 0, 1, 3, 2
  };
  std::vector<PDM_g_num_t> prism = {
      1, 8, 3,5,9,7,
     12,11,15,1,8,3
  };
  std::vector<PDM_g_num_t> pyra = {
     4,5,7,6,10
  };
  std::vector<PDM_g_num_t> tet = {
     5,9,7,10
  };

};


MPI_TEST_CASE("part by elt", 1) {
  // create multipart
  int n_zone = 1;
  std::vector<int> n_part_per_zone = {1};
  PDM_bool_t merge_block = PDM_FALSE; 
  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  PDM_part_size_t part_weight_method = PDM_PART_SIZE_HOMOGENEOUS;
  double* part_weight = nullptr;
  PDM_ownership_t ownership = PDM_OWNERSHIP_USER;
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int multi_part_id = PDM_multipart_create(n_zone, n_part_per_zone.data(), merge_block, split_method, part_weight_method, part_weight, pdm_comm, ownership);

  // create dmesh
  // 0. mesh construction {
  simple_mesh2_for_test m;
  std_e::offset(m.hex,1);
  std_e::offset(m.prism,1);
  std_e::offset(m.pyra,1);
  std_e::offset(m.tet,1);

  // 1. mesh distribution {
  int n_vtx = 16;
  int n_cell = 6;
  // 1. mesh distribution }

  // 2. dmesh_nodal creation {
  PDM_dmesh_nodal_t* dmesh_nodal = PDM_DMesh_nodal_create(pdm_comm,3,n_vtx,n_cell,0,0); // 3: dim, 0: n_faces, 0: n_edge
  //PDM_DMesh_nodal_coord_set(dmesh_nodal, dn_vtx, dvtx_coord.data());
  int id_tet_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_TETRA4  );
  int id_pyra_section  = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PYRAMID5);
  int id_prism_section = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PRISM6  );
  int id_hex_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_HEXA8   );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_tet_section  ,1,m.tet.data()  );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_pyra_section ,1,m.pyra.data() );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_prism_section,2,m.prism.data());
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_hex_section  ,2,m.hex.data()  );
  // 2. dmesh_nodal creation }

  // populate multipart with dmesh_nodal
  int zone_id = 0;
  PDM_multipart_register_dmesh_nodal(multi_part_id, zone_id, dmesh_nodal);

  //int* renum_properties_cell_data = nullptr;
  //const char* renum_cell_method = "PDM_PART_RENUM_CELL_NONE";
  //const char* renum_face_method = "PDM_PART_RENUM_FACE_NONE";
  //PDM_multipart_set_reordering_options(multi_part_id,
  //                                     -1, // -1: all zones
  //                                     renum_cell_method,
  //                                     renum_properties_cell_data,
  //                                     renum_face_method);
  PDM_multipart_run_ppart(multi_part_id);

  PDM_DMesh_nodal_free(dmesh_nodal,0);
}


