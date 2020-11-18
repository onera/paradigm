#include "std_e/unit_test/doctest.hpp"
#include "doctest/extensions/doctest_mpi.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "std_e/utils/concatenate.hpp"
#include "std_e/future/span.hpp"
#include "std_e/algorithm/algorithm.hpp"
#include "pdm_dmesh.h"
#include "pdm_distrib.h"

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

  std::vector<PDM_g_num_t> tet = {
     5,9,7,10
  };
  std::vector<PDM_g_num_t> pyra = {
     4,5,7,6,10
  };
  std::vector<PDM_g_num_t> prism = {
      1, 8, 3,5,9,7,
     12,11,15,1,8,3
  };
  std::vector<PDM_g_num_t> hex = {
     0, 1, 3, 2, 4, 5, 7, 6,
     13,12,15,14, 0, 1, 3, 2
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
  int n_elt = 6;
  auto vtx_coord = std_e::concatenate(m.coord_X,m.coord_Y,m.coord_Z);
  // 1. mesh distribution }

  // 2. dmesh_nodal creation {
  PDM_dmesh_nodal_t* dmesh_nodal = PDM_DMesh_nodal_create(pdm_comm,3,n_vtx,n_elt,0,0); // 3: dim, 0: n_faces, 0: n_edge
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

  //SUBCASE("concat_elt_sections") {
  //  int* section_idx;
  //  int* delt_vtx_idx;
  //  PDM_g_num_t* delt_vtx;
  //  int n_section = PDM_concat_elt_sections(dmesh_nodal,&section_idx,&delt_vtx_idx,&delt_vtx);

  //  auto _section_idx = std_e::make_span(section_idx,n_section+1);
  //  auto _delt_vtx_idx = std_e::make_span(delt_vtx_idx,_section_idx.back()+1);
  //  auto _delt_vtx = std_e::make_span(delt_vtx,_delt_vtx_idx.back());
  //  std_e::offset(_delt_vtx,-1);

  //  CHECK(n_section == 4);
  //  CHECK( _section_idx   == std::vector{0       ,  1         ,  2                           ,  4                                    ,6 } );
  //  CHECK( _delt_vtx_idx == std::vector{0       ,  4         ,  9          ,  15            ,  21             ,  29                 ,37} );
  //  CHECK( _delt_vtx ==     std::vector{5,9,7,10,  4,5,7,6,10,  1,8,3,5,9,7,  12,11,15,1,8,3,  0,1,3,2,4,5,7,6,  13,12,15,14,0,1,3,2   } );

  //  free(section_idx);
  //  free(delt_vtx_idx);
  //  free(delt_vtx);
  //}

  //SUBCASE("dual graph") {
  //  int n_rank = test_nb_procs;

  //  int* section_idx;
  //  int* delt_vtx_idx;
  //  PDM_g_num_t* delt_vtx;
  //  int n_section = PDM_concat_elt_sections(dmesh_nodal,&section_idx,&delt_vtx_idx,&delt_vtx);


  //  PDM_g_num_t* mesh_vtx_dist = PDM_dmesh_nodal_vtx_distrib_get(dmesh_nodal);
  //  PDM_g_num_t*  vtx_dist = (PDM_g_num_t*) malloc(n_rank+1 * sizeof(PDM_g_num_t));
  //  for (int i=0; i<n_rank+1; ++i) {
  //    vtx_dist[i] = mesh_vtx_dist[i]+1;
  //  }

  //  int dn_elt = section_idx[n_section];
  //  PDM_g_num_t* elt_dist = PDM_compute_entity_distribution(pdm_comm, dn_elt);

  //  PDM_g_num_t* elt_elt_idx;
  //  PDM_g_num_t* elt_elt;
  //  PDM_dmesh_nodal_dual_graph(vtx_dist,elt_dist,delt_vtx_idx,delt_vtx,&elt_elt_idx,&elt_elt,pdm_comm);

  //  auto ee_idx = std_e::make_span(elt_elt_idx,n_elt+1);
  //  auto ee = std_e::make_span(elt_elt,ee_idx.back());

  //  std::vector<int> ee_idx_expected     = {0      , 3      , 6            ,11      ,14            ,19      , 22};
  //  std::vector<PDM_g_num_t> ee_expected = {5, 2, 3, 1, 5, 3, 1, 2, 6, 4, 5, 3, 6, 5, 1, 2, 3, 4, 6, 3, 4, 5,};

  //  CHECK( ee_idx == ee_idx_expected );
  //  CHECK( ee == ee_expected );

  //  free(section_idx);
  //  free(delt_vtx_idx);
  //  free(delt_vtx);

  //  free(elt_dist);
  //  free(vtx_dist);
  //  free(elt_elt_idx);
  //  free(elt_elt);
  //}

  SUBCASE("run ppart") {
    int dn_part = 1;
    _run_ppart_zone_nodal(dmesh_nodal,PDM_SPLIT_DUAL_WITH_PTSCOTCH,dn_part,pdm_comm);
  }

  PDM_DMesh_nodal_free(dmesh_nodal,0);
}


