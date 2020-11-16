#include "doctest/extensions/doctest_mpi.h"

#include "pdm_para_graph_dual.h"
#include "pdm_multipart.h"
#include "std_e/utils/concatenate.hpp"
#include "std_e/future/span.hpp"
#include "std_e/algorithm/algorithm.hpp"
#include "pdm_dmesh.h"

#include <vector>

template<class Contiguous_range, class Distribution> auto
select_block(Contiguous_range& r, const Distribution& dist, int rank) {
  return std_e::make_span(r,dist[rank],dist[rank+1]);
}

struct simple_mesh_for_test {
  // SEE meshes in simple_meshes.hpp
  std::vector<double> coord_X = 
    { 0.,1.,2.,3.,
      0.,1.,2.,3.,
      0.,1.,2.,3.,
      0.,1.,2.,3.,
      0.,1.,2.,3.,
      0.,1.,2.,3. };
  std::vector<double> coord_Y =
    { 0.,0.,0.,0.,
      1.,1.,1.,1.,
      2.,2.,2.,2.,
      0.,0.,0.,0.,
      1.,1.,1.,1.,
      2.,2.,2.,2. };
  std::vector<double> coord_Z =
    { 0.,0.,0.,0.,
      0.,0.,0.,0.,
      0.,0.,0.,0.,
      1.,1.,1.,1.,
      1.,1.,1.,1.,
      1.,1.,1.,1. };

  std::vector<PDM_g_num_t> cell_vtx = {
     0, 1, 5, 4,12,13,17,16,
     1, 3, 6, 5,13,14,18,17,
     2, 3, 7, 6,14,15,19,18,
     4, 5, 9, 8,16,17,21,20,
     5, 6,10, 9,17,18,22,21,
     6, 7,11,10,18,19,23,22
  };

  std::vector<PDM_g_num_t> face_vtx_i = {
     0, 4,16,12,   1, 5,17,13,   2, 6,18,14,   3, 7,19,15,
     4, 8,20,16,   5, 9,21,17,   6,10,22,18,   7,11,23,19
  };
  std::vector<PDM_g_num_t> face_vtx_j = {
     0,12,13, 1,   1,13,14, 2,   2,14,15, 3,
     4,16,17, 5,   5,17,18, 6,   6,18,19, 7,
     8,20,21, 9,   9,21,22,10,  10,22,23,11
  };
  std::vector<PDM_g_num_t> face_vtx_k = {
     0, 1, 5, 4,   1, 2, 6, 5,   2, 3, 7, 6,
     4, 5, 9, 8,   5, 6,10, 9,   6, 7,11,10,
    12,13,17,16,  13,14,18,17,  14,15,19,18,
    16,17,21,20,  17,18,22,21,  18,19,23,22
  };
  std::vector<PDM_g_num_t> face_cell_i = { 
    -1, 0      ,   0, 1      ,   1, 2      ,   2,-1      ,
    -1, 3      ,   3, 4      ,   4, 5      ,   5,-1
  };
  std::vector<PDM_g_num_t> face_cell_j = { 
    -1, 0      ,  -1, 1      ,  -1, 2      ,
     0, 3      ,   1, 4      ,   2, 5      ,
     3,-1      ,   4,-1      ,   5,-1
  };
  std::vector<PDM_g_num_t> face_cell_k = { 
    -1, 0      ,  -1, 1      ,  -1, 2      ,  
    -1, 3      ,  -1, 4      ,  -1, 5      ,  
     0,-1      ,   1,-1      ,   2,-1      ,
     3,-1      ,   4,-1      ,   5,-1      
  };
  std::vector<PDM_g_num_t> face_vtx  = std_e::concatenate(face_vtx_i ,face_vtx_j ,face_vtx_k );
  //auto face_cell = std_e::concatenate(face_cell_i,face_cell_j,face_cell_k);
  std::vector<PDM_g_num_t> face_cell_l = {
    -1         ,   0         ,   1         ,   2         ,
    -1         ,   3         ,   4         ,   5         ,

    -1         ,  -1         ,  -1         ,
     0         ,   1         ,   2         ,
     3         ,   4         ,   5         ,

    -1         ,  -1         ,  -1         ,  
    -1         ,  -1         ,  -1         ,  
     0         ,   1         ,   2         ,
     3         ,   4         ,   5         
  };
  std::vector<PDM_g_num_t> face_cell_r = {
        0      ,      1      ,      2      ,     -1      ,
        3      ,      4      ,      5      ,     -1      ,

        0      ,      1      ,      2      ,
        3      ,      4      ,      5      ,
       -1      ,     -1      ,     -1      ,

        0      ,      1      ,      2      ,  
        3      ,      4      ,      5      ,  
       -1      ,     -1      ,     -1      ,
       -1      ,     -1      ,     -1      
  };
  std::vector<PDM_g_num_t> face_cell = std_e::concatenate(face_cell_l,face_cell_r);
};


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
        |  \   13     / |    \ 12      \      /
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
  //std::vector<double> coord_X = 
  //  { 0.,1.,2.,3.,
  //    0.,1.,2.,3. };
  //std::vector<double> coord_Y =
  //  { 0.,0.,0.,0.,
  //    2.,2.,2.,2. };
  //std::vector<double> coord_Z =
  //  { 0.,0.,0.,0.,
  //    1.,1.,1.,1. };

  std::vector<PDM_g_num_t> hex = {
     0, 1, 3, 2, 4, 5, 7, 6,
     13,12,15,14, 0, 1, 3, 2
  };
  std::vector<PDM_g_num_t> prism = {
      1, 8, 3,5,9,7,
     12,11,15,1,8,3
  };
  std::vector<PDM_g_num_t> pyra = {
     4,5,5,6,10
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
  int n_vtx = 11;
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
  //int zone_id = 0;
  //PDM_multipart_register_dmesh_nodal(multi_part_id, zone_id, dmesh_nodal);

  ////int* renum_properties_cell_data = nullptr;
  ////const char* renum_cell_method = "PDM_PART_RENUM_CELL_NONE";
  ////const char* renum_face_method = "PDM_PART_RENUM_FACE_NONE";
  ////PDM_multipart_set_reordering_options(multi_part_id,
  ////                                     -1, // -1: all zones
  ////                                     renum_cell_method,
  ////                                     renum_properties_cell_data,
  ////                                     renum_face_method);
  //PDM_multipart_run_ppart(multi_part_id);

  //PDM_DMesh_nodal_free(dmesh_nodal,0);
}


//MPI_TEST_CASE("part by elt hex", 1) {
//  // create multipart
//  int n_zone = 1;
//  std::vector<int> n_part_per_zone = {1};
//  PDM_bool_t merge_block = PDM_FALSE; 
//  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
//  PDM_part_size_t part_weight_method = PDM_PART_SIZE_HOMOGENEOUS;
//  double* part_weight = nullptr;
//  PDM_ownership_t ownership = PDM_OWNERSHIP_USER;
//  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
//  int multi_part_id = PDM_multipart_create(n_zone, n_part_per_zone.data(), merge_block, split_method, part_weight_method, part_weight, pdm_comm, ownership);
//
//  // create dmesh
//  // 0. mesh construction {
//  simple_mesh_for_test m;
//  std_e::offset(m.cell_vtx,1);
//
//  // 1. mesh distribution {
//  int n_vtx = 4*3*2;
//  int n_cell = 3*2*1;
//  int n_face = 2*4 + 3*3 + 6*2;
//
//  //std::vector<int> vtx_dist = {0,n_vtx/2,n_vtx};
//  //std::vector<int> face_dist = {0,n_face/2,n_face};
//  //std::vector<int> cell_dist = {0,n_cell/2,n_cell};
//  std::vector<int> vtx_dist = {0,n_vtx};
//  std::vector<int> face_dist = {0,n_face};
//  std::vector<int> cell_dist = {0,n_cell};
//
//  std::vector<int> face_vtx_dist = face_dist;
//  std_e::scale(face_vtx_dist,4);
//  std::vector<int> face_cell_dist = face_dist;
//  std_e::scale(face_cell_dist,2);
//
//  int dn_vtx  = vtx_dist[test_rank+1] - vtx_dist[test_rank];
//  int dn_face = face_dist[test_rank+1] - face_dist[test_rank];
//  int dn_cell = cell_dist[test_rank+1] - cell_dist[test_rank];
//
//  int dn_face_vtx = face_vtx_dist[test_rank+1] - face_vtx_dist[test_rank];
//
//  auto dcoord_X = select_block(m.coord_X,vtx_dist,test_rank);
//  auto dcoord_Y = select_block(m.coord_Y,vtx_dist,test_rank);
//  auto dcoord_Z = select_block(m.coord_Z,vtx_dist,test_rank);
//  auto dvtx_coord = std_e::concatenate(dcoord_X,dcoord_Y,dcoord_Z);
//
//  auto dface_vtx = select_block(m.face_vtx,face_vtx_dist,test_rank);
//  auto dface_cell = select_block(m.face_cell,face_cell_dist,test_rank);
//  // 1. mesh distribution }
//
//  // 2. dmesh_nodal creation {
//  int dn_bnd = 0;
//  int dn_join = 0;
//  PDM_dmesh_nodal_t* dmesh_nodal = PDM_DMesh_nodal_create(pdm_comm,3,n_vtx,n_cell,0,0); // 3: dim, 0: n_faces, 0: n_edge
//  PDM_DMesh_nodal_coord_set(dmesh_nodal, dn_vtx, dvtx_coord.data());
//  int id_hex_section = PDM_DMesh_nodal_section_add(dmesh_nodal,PDM_MESH_NODAL_HEXA8);
//  PDM_DMesh_nodal_section_std_set(dmesh_nodal,id_hex_section,6,m.cell_vtx.data());
//  // 2. dmesh_nodal creation }
//
//  // populate multipart with dmesh_nodal
//  //int zone_id = 0;
//  //PDM_multipart_register_dmesh_nodal(multi_part_id, zone_id, dmesh_nodal);
//
//  ////int* renum_properties_cell_data = nullptr;
//  ////const char* renum_cell_method = "PDM_PART_RENUM_CELL_NONE";
//  ////const char* renum_face_method = "PDM_PART_RENUM_FACE_NONE";
//  ////PDM_multipart_set_reordering_options(multi_part_id,
//  ////                                     -1, // -1: all zones
//  ////                                     renum_cell_method,
//  ////                                     renum_properties_cell_data,
//  ////                                     renum_face_method);
//  //PDM_multipart_run_ppart(multi_part_id);
//
//  //PDM_DMesh_nodal_free(dmesh_nodal,0);
//}

//MPI_TEST_CASE("part by faces", 1) {
//  // create multipart
//  int n_zone = 1;
//  std::vector<int> n_part_per_zone = {1};
//  PDM_bool_t merge_block = PDM_FALSE; 
//  PDM_split_dual_t split_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
//  PDM_part_size_t part_weight_method = PDM_PART_SIZE_HOMOGENEOUS;
//  double* part_weight = nullptr;
//  PDM_ownership_t ownership = PDM_OWNERSHIP_USER;
//  PDM_MPI_Comm pdm_comm = PDM_MPI_mpi_2_pdm_mpi_comm(&test_comm);
//  int multi_part_id = PDM_multipart_create(n_zone, n_part_per_zone.data(), merge_block, split_method, part_weight_method, part_weight, pdm_comm, ownership);
//
//  // create dmesh
//  // 0. mesh construction {
//  simple_mesh_for_test m;
//  std_e::offset(m.face_vtx,1);
//  std_e::offset(m.face_cell,1);
//
//  // 1. mesh distribution {
//  int n_vtx = 4*3*2;
//  int n_cell = 3*2*1;
//  int n_face = 2*4 + 3*3 + 6*2;
//
//  //std::vector<int> vtx_dist = {0,n_vtx/2,n_vtx};
//  //std::vector<int> face_dist = {0,n_face/2,n_face};
//  //std::vector<int> cell_dist = {0,n_cell/2,n_cell};
//  std::vector<int> vtx_dist = {0,n_vtx};
//  std::vector<int> face_dist = {0,n_face};
//  std::vector<int> cell_dist = {0,n_cell};
//
//  std::vector<int> face_vtx_dist = face_dist;
//  std_e::scale(face_vtx_dist,4);
//  std::vector<int> face_cell_dist = face_dist;
//  std_e::scale(face_cell_dist,2);
//
//  int dn_vtx  = vtx_dist[test_rank+1] - vtx_dist[test_rank];
//  int dn_face = face_dist[test_rank+1] - face_dist[test_rank];
//  int dn_cell = cell_dist[test_rank+1] - cell_dist[test_rank];
//
//  int dn_face_vtx = face_vtx_dist[test_rank+1] - face_vtx_dist[test_rank];
//
//  auto dcoord_X = select_block(m.coord_X,vtx_dist,test_rank);
//  auto dcoord_Y = select_block(m.coord_Y,vtx_dist,test_rank);
//  auto dcoord_Z = select_block(m.coord_Z,vtx_dist,test_rank);
//  auto dvtx_coord = std_e::concatenate(dcoord_X,dcoord_Y,dcoord_Z);
//
//  auto dface_vtx = select_block(m.face_vtx,face_vtx_dist,test_rank);
//  auto dface_cell = select_block(m.face_cell,face_cell_dist,test_rank);
//  // 1. mesh distribution }
//
//  // 2. dmesh creation {
//  int dn_bnd = 0;
//  int dn_join = 0;
//  PDM_dmesh_t* dmesh = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, dn_bnd, dn_join);
//  int* dface_bound_idx = nullptr;
//  PDM_g_num_t* dface_bound = nullptr;
//  int* joins_ids = nullptr;
//  int* dface_join_idx = nullptr;
//  int* dface_join = nullptr;
//  std::vector<int> dface_vtx_idx(dn_face_vtx,4);
//  PDM_dmesh_set(
//    dmesh,
//    dvtx_coord.data(),
//    dface_vtx_idx.data(), dface_vtx.data(), 
//    dface_cell.data(),
//    dface_bound_idx, dface_bound,
//    joins_ids,
//    dface_join_idx, dface_join
//  );
//  // 2. dmesh creation }
//
//  // populate multipart with dmesh
//  int zone_id = 0;
//  PDM_multipart_register_block(multi_part_id, zone_id, dmesh);
//
//  int* renum_properties_cell_data = nullptr;
//  const char* renum_cell_method = "PDM_PART_RENUM_CELL_NONE";
//  const char* renum_face_method = "PDM_PART_RENUM_FACE_NONE";
//  PDM_multipart_set_reordering_options(multi_part_id,
//                                       -1, // -1: all zones
//                                       renum_cell_method,
//                                       renum_properties_cell_data,
//                                       renum_face_method);
//  PDM_multipart_run_ppart(multi_part_id);
//  PDM_dmesh_free(dmesh);
//}
