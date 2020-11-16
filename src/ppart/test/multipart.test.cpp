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

MPI_TEST_CASE("part", 1) {
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

  std::vector<int> face_vtx_i = {
     0, 4,16,12,   1, 5,17,13,   2, 6,18,14,   3, 7,19,15,
     4, 8,20,16,   5, 9,21,17,   6,10,22,18,   7,11,23,19
  };
  std::vector<int> face_vtx_j = {
     0,12,13, 1,   1,13,14, 2,   2,14,15, 3,
     4,16,17, 5,   5,17,18, 6,   6,18,19, 7,
     8,20,21, 9,   9,21,22,10,  10,22,23,11
  };
  std::vector<int> face_vtx_k = {
     0, 1, 5, 4,   1, 2, 6, 5,   2, 3, 7, 6,
     4, 5, 9, 8,   5, 6,10, 9,   6, 7,11,10,
    12,13,17,16,  13,14,18,17,  14,15,19,18,
    16,17,21,20,  17,18,22,21,  18,19,23,22
  };
  std::vector<int> face_cell_i = { 
    -1, 0      ,   0, 1      ,   1, 2      ,   2,-1      ,
    -1, 3      ,   3, 4      ,   4, 5      ,   5,-1
  };
  std::vector<int> face_cell_j = { 
    -1, 0      ,  -1, 1      ,  -1, 2      ,
     0, 3      ,   1, 4      ,   2, 5      ,
     3,-1      ,   4,-1      ,   5,-1
  };
  std::vector<int> face_cell_k = { 
    -1, 0      ,  -1, 1      ,  -1, 2      ,  
    -1, 3      ,  -1, 4      ,  -1, 5      ,  
     0,-1      ,   1,-1      ,   2,-1      ,
     3,-1      ,   4,-1      ,   5,-1      
  };
  auto face_vtx  = std_e::concatenate(face_vtx_i ,face_vtx_j ,face_vtx_k );
  auto face_cell = std_e::concatenate(face_cell_i,face_cell_j,face_cell_k);
 
  // offseting to 1
  std_e::offset(face_vtx,1);
  std_e::offset(face_cell,1);
  // mesh construction }

  // 1. mesh distribution {
  int n_vtx = 4*3*2;
  int n_cell = 3*2*1;
  int n_face = 2*4 + 3*3 + 6*2;

  std::vector<int> vtx_dist = {0,n_vtx/2,n_vtx};
  std::vector<int> face_dist = {0,n_face/2,n_face};
  std::vector<int> cell_dist = {0,n_cell/2,n_cell};

  std::vector<int> face_vtx_dist = face_dist;
  std_e::scale(face_vtx_dist,4);
  std::vector<int> face_cell_dist = face_dist;
  std_e::scale(face_cell_dist,2);

  int dn_vtx  = vtx_dist[test_rank+1] - vtx_dist[test_rank];
  int dn_face = face_dist[test_rank+1] - face_dist[test_rank];
  int dn_cell = cell_dist[test_rank+1] - cell_dist[test_rank];

  int dn_face_vtx = face_vtx_dist[test_rank+1] - face_vtx_dist[test_rank];

  auto dcoord_X = select_block(coord_X,vtx_dist,test_rank);
  auto dcoord_Y = select_block(coord_Y,vtx_dist,test_rank);
  auto dcoord_Z = select_block(coord_Z,vtx_dist,test_rank);
  auto dvtx_coord = std_e::concatenate(dcoord_X,dcoord_Y,dcoord_Z);

  auto dface_vtx = select_block(face_vtx,face_vtx_dist,test_rank);
  auto dface_cell = select_block(face_cell,face_cell_dist,test_rank);
  // 1. mesh distribution }

  // 2. dmesh creation {
  int dn_bnd = 0;
  int dn_join = 0;
  int dmesh_id = PDM_dmesh_create(dn_cell, dn_face, dn_vtx, dn_bnd, dn_join);
  int* dface_bound_idx = nullptr;
  int* dface_bound = nullptr;
  int* joins_ids = nullptr;
  int* dface_join_idx = nullptr;
  int* dface_join = nullptr;
  std::vector<int> dface_vtx_idx(dn_face_vtx,4);
  PDM_dmesh_set(
    dmesh_id,
    dvtx_coord.data(),
    dface_vtx_idx.data(), dface_vtx.data(), 
    dface_cell.data(),
    dface_bound_idx, dface_bound,
    joins_ids,
    dface_join_idx, dface_join
  );
  // 2. dmesh creation }

  // populate multipart with dmesh
  int zone_id = 0;
  PDM_multipart_register_block(multi_part_id, zone_id, dmesh_id);
}
