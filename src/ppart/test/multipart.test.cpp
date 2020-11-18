#include "std_e/unit_test/doctest.hpp"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm_para_graph_dual.h"
#include "pdm_dmesh_nodal.h"
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
                              // 0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10,11,12,13,14,15
  std::vector<double> coord_X = {0.,1.,0.,1.,0.,1.,0.,1.,2.,2.,1., 2., 1., 0., 0., 1.};
  std::vector<double> coord_Y = {0.,0.,1.,1.,0.,0.,1.,1.,0.,0.,0., 0., 0., 0., 1., 1.};
  std::vector<double> coord_Z = {0.,0.,0.,0.,1.,1.,1.,1.,0.,1.,2.,-1.,-1.,-1.,-1.,-1.};

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
  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
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
  auto vtx_coord = std_e::concatenate(m.coord_X,m.coord_Y,m.coord_Z);
  // 1. mesh distribution }

  // 2. dmesh_nodal creation {
  PDM_dmesh_nodal_t* dmesh_nodal = PDM_DMesh_nodal_create(pdm_comm,3,n_vtx,n_cell,0,0); // 3: dim, 0: n_faces, 0: n_edge
  PDM_DMesh_nodal_coord_set(dmesh_nodal, n_vtx, vtx_coord.data());
  int id_tet_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_TETRA4  );
  int id_pyra_section  = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PYRAMID5);
  int id_prism_section = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PRISM6  );
  int id_hex_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_HEXA8   );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_tet_section  ,1,m.tet.data()  );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_pyra_section ,1,m.pyra.data() );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_prism_section,2,m.prism.data());
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_hex_section  ,2,m.hex.data()  );
  // 2. dmesh_nodal creation }

  //SUBCASE("dual graph") {
  //  PDM_g_num_t* cell_cell_idx;
  //  PDM_g_num_t* cell_cell;
  //  PDM_g_num_t* cell_dist;
  //  PDM_dmesh_nodal_dual_graph(dmesh_nodal,&cell_cell_idx,&cell_cell,3,pdm_comm,&cell_dist);

  //  auto cc_idx = std_e::make_span(cell_cell_idx,n_cell+1);
  //  auto cc = std_e::make_span(cell_cell,cc_idx.back());

  //  std::vector<int> cc_idx_expected     = {0      , 3      , 6            ,11      ,14            ,19      , 22};
  //  std::vector<PDM_g_num_t> cc_expected = {5, 2, 3, 1, 5, 3, 1, 2, 6, 4, 5, 3, 6, 5, 1, 2, 3, 4, 6, 3, 4, 5,};

  //  CHECK( cc_idx == cc_idx_expected );
  //  CHECK( cc == cc_expected );

  //  free(cell_cell_idx);
  //  free(cell_cell);
  //  free(cell_dist);
  //}

  SUBCASE("run ppart") {
    _run_ppart_zone_nodal(dmesh_nodal,PDM_SPLIT_DUAL_WITH_PTSCOTCH,pdm_comm);
  }

  PDM_DMesh_nodal_free(dmesh_nodal,0);
}


