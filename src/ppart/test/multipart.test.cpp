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
#include "pdm_partitioning_algorithm.h"

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
  PDM_MPI_Comm comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
  int n_rank = test_nb_procs;
  // create dmesh
  // 0. mesh construction {
  simple_mesh2_for_test m;
  std_e::offset(m.hex,1);
  std_e::offset(m.prism,1);
  std_e::offset(m.pyra,1);
  std_e::offset(m.tet,1);

  // 1. mesh distribution {
  int n_vtx = 16;
  int n_tet   = m.tet.size()/4;
  int n_pyra  = m.pyra.size()/5;
  int n_prism = m.prism.size()/6;
  int n_hex   = m.hex.size()/8;
  int n_elt = n_tet + n_pyra + n_prism + n_hex;
  auto vtx_coord = std_e::concatenate(m.coord_X,m.coord_Y,m.coord_Z);
  // 1. mesh distribution }

  // 2. dmesh_nodal creation {
  PDM_dmesh_nodal_t* dmesh_nodal = PDM_DMesh_nodal_create(comm,3,n_vtx,n_elt,0,0); // 3: dim, 0: n_faces, 0: n_edge
  PDM_DMesh_nodal_coord_set(dmesh_nodal, n_vtx, vtx_coord.data());
  int id_tet_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_TETRA4  );
  int id_pyra_section  = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PYRAMID5);
  int id_prism_section = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_PRISM6  );
  int id_hex_section   = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_HEXA8   );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_tet_section  ,n_tet  ,m.tet.data()  );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_pyra_section ,n_pyra ,m.pyra.data() );
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_prism_section,n_prism,m.prism.data());
  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_hex_section  ,n_hex  ,m.hex.data()  );
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
  //  PDM_g_num_t* elt_dist = PDM_compute_entity_distribution(comm, dn_elt);

  //  PDM_g_num_t* elt_elt_idx;
  //  PDM_g_num_t* elt_elt;
  //  PDM_dmesh_nodal_dual_graph(vtx_dist,elt_dist,delt_vtx_idx,delt_vtx,&elt_elt_idx,&elt_elt,comm);

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

  //// 0. concat sections
  //SUBCASE("PDM_part_multi_dconnectivity_to_pconnectivity_sort") {
  //  auto split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
  //  int dn_part = 1;

  //  int* section_idx;
  //  int* delt_vtx_idx;
  //  PDM_g_num_t* delt_vtx;
  //  int n_section = PDM_concat_elt_sections(dmesh_nodal,&section_idx,&delt_vtx_idx,&delt_vtx);

  //  // 1. distributions
  //  PDM_g_num_t* mesh_vtx_dist = PDM_dmesh_nodal_vtx_distrib_get(dmesh_nodal);
  //  PDM_g_num_t* vtx_dist = (PDM_g_num_t*) malloc(n_rank+1 * sizeof(PDM_g_num_t));
  //  for (int i=0; i<n_rank+1; ++i) {
  //    vtx_dist[i] = mesh_vtx_dist[i]+1;
  //  }

  //  int dn_elt = section_idx[n_section];
  //  PDM_g_num_t* elt_dist = PDM_compute_entity_distribution(comm, dn_elt);

  //  // 2. elt_elt graph
  //  PDM_g_num_t* elt_elt_idx;
  //  PDM_g_num_t* elt_elt;
  //  PDM_dmesh_nodal_dual_graph(vtx_dist,elt_dist,delt_vtx_idx,delt_vtx,&elt_elt_idx,&elt_elt,comm);

  //  // 3. partitioning
  //  double* part_fractions = NULL; // TODO gen
  //  int* elt_part = (int*) malloc(dn_elt * sizeof(int));
  //  for (int i=0; i<elt_elt_idx[dn_elt]; ++i) {
  //    elt_elt[i]--;
  //  }
  //  PDM_para_graph_split(split_method,
  //                       elt_dist,
  //                       elt_elt_idx,
  //                       elt_elt,
  //                       NULL, NULL,
  //                       dn_part,
  //                       part_fractions,
  //                       elt_part,
  //                       comm);

  //  // 4. reconstruct elts on partitions
  //  PDM_g_num_t* part_distri = PDM_compute_entity_distribution(comm, dn_part );
  //  //int n_part = part_distri[n_rank]-1;
  //  int** pn_elt_section = (int**)malloc(n_section * sizeof(int*));
  //  PDM_g_num_t*** pelt_section_ln_to_gn = (PDM_g_num_t***)malloc(n_section * sizeof(PDM_g_num_t**));
  //  PDM_g_num_t** elt_section_distri = (PDM_g_num_t**)malloc(n_section * sizeof(PDM_g_num_t*));
  //  for (int i_section=0; i_section<n_section; ++i_section) {
  //    elt_section_distri[i_section] = PDM_DMesh_nodal_section_distri_std_get(dmesh_nodal,i_section);
  //  }

  //  for (int i_section=0; i_section<n_section; ++i_section) {
  //    int* elt_section_part = elt_part + section_idx[i_section];

  //    PDM_part_assemble_partitions(comm,
  //                                 part_distri,
  //                                 elt_section_distri[i_section],
  //                                 elt_section_part,
  //                                &pn_elt_section[i_section],
  //                                &pelt_section_ln_to_gn[i_section]);
  //    //printf("\npelt_section_ln_to_gn:");
  //    //for (int i = 0; i < pn_elt_section[i_section][0]; i++)
  //    //  printf(" %d ", pelt_section_ln_to_gn[i_section][0][i]);
  //    //printf("\n");
  //  }

  //  int* pn_vtx;
  //  int** pvtx_ln_to_gn;
  //  int*** pelt_vtx_idx;
  //  int*** pelt_vtx;
  //  PDM_part_multi_dconnectivity_to_pconnectivity_sort(comm,
  //                                                     dn_part,
  //                                                     n_section,
  //                                                     section_idx,
  //                                                     elt_section_distri,
  //                                                     delt_vtx_idx,
  //                                                     delt_vtx,
  //                                                     pn_elt_section,
  //                                                     pelt_section_ln_to_gn,
  //                                                    &pn_vtx,
  //                                                    &pvtx_ln_to_gn,
  //                                                    &pelt_vtx_idx,
  //                                                    &pelt_vtx);
  //  // the [0] are there because only one partition
  //  CHECK( pn_vtx[0] == 16 );

  //  auto _pvtx_ln_to_gn = std_e::make_span(pvtx_ln_to_gn[0],pn_vtx[0]);
  //  CHECK( _pvtx_ln_to_gn == std::vector{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} );

  //  auto _ptet_vtx_idx = std_e::make_span(pelt_vtx_idx[0][0],n_tet+1);
  //  auto _ptet_vtx = std_e::make_span(pelt_vtx[0][0],_ptet_vtx_idx.back());
  //  CHECK( _ptet_vtx_idx == std::vector{0,4} );
  //  CHECK( _ptet_vtx == std::vector{6,10,8,11} );

  //  auto _phex_vtx_idx = std_e::make_span(pelt_vtx_idx[3][0],n_hex+1);
  //  auto _phex_vtx = std_e::make_span(pelt_vtx[3][0],_phex_vtx_idx.back());
  //  CHECK( _phex_vtx_idx == std::vector{0,8,16} );
  //  CHECK( _phex_vtx == std::vector{1, 2, 4, 3, 5, 6, 8, 7,  14,13,16,15, 1, 2, 4, 3} );
  //}
     
  SUBCASE("run ppart") {
    int dn_part = 1;
    _run_ppart_zone_nodal(dmesh_nodal,PDM_SPLIT_DUAL_WITH_PTSCOTCH,dn_part,comm);
  }

  PDM_DMesh_nodal_free(dmesh_nodal,0);
}


