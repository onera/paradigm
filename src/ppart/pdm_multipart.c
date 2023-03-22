/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2017       ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*============================================================================
 * TODO : write module description here
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_distrib.h"
#include "pdm_timer.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_para_graph_dual.h"
#include "pdm_part_geom.h"
#include "pdm_part_renum.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_mesh_nodal.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_binary_search.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_multi_block_to_part.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_para_graph_dual.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_unique.h"
#include "pdm_partitioning_nodal_algorithm.h"
#include "pdm_dmesh_nodal_to_dmesh.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_multipart_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/


/*=============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Translate in _part_t structure (only mapping)
 */
static
_part_t**
_map_part_t_with_part_mesh
(
  PDM_part_mesh_t* pm
)
{
  int n_part = pm->n_part;
  _part_t **pdm_part = (_part_t **) malloc(n_part * sizeof(_part_t *));
  for(int i_part = 0; i_part < n_part; ++i_part) {

    pdm_part[i_part] = _part_create();

    pdm_part[i_part]->n_cell                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_CELL  );
    pdm_part[i_part]->n_face                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_FACE  );
    pdm_part[i_part]->n_edge                   = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_EDGE  );
    pdm_part[i_part]->n_vtx                    = PDM_part_mesh_n_entity_get(pm, i_part, PDM_MESH_ENTITY_VERTEX);
    pdm_part[i_part]->n_section                = 0;
    pdm_part[i_part]->n_elt                    = NULL;

    pdm_part[i_part]->n_face_group             = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_FACE);
    pdm_part[i_part]->n_edge_group             = PDM_part_mesh_n_bound_get(pm, PDM_BOUND_TYPE_EDGE);

    pdm_part[i_part]->n_face_part_bound        = 0;
    pdm_part[i_part]->n_vtx_part_bound         = 0;

    PDM_part_mesh_vtx_coord_get(pm, i_part, &pdm_part[i_part]->vtx, PDM_OWNERSHIP_BAD_VALUE);
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   &pdm_part[i_part]->face_vtx,
                                   &pdm_part[i_part]->face_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                   &pdm_part[i_part]->cell_face,
                                   &pdm_part[i_part]->cell_face_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                   &pdm_part[i_part]->face_edge,
                                   &pdm_part[i_part]->face_edge_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                   &pdm_part[i_part]->edge_face,
                                   &pdm_part[i_part]->edge_face_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    int* face_cell_idx = NULL;
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                   &pdm_part[i_part]->face_cell,
                                   &face_cell_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    int* edge_vtx_idx = NULL;
    PDM_part_mesh_connectivity_get(pm,
                                   i_part,
                                   PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                   &pdm_part[i_part]->edge_vtx,
                                   &edge_vtx_idx,
                                   PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &pdm_part[i_part]->cell_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pdm_part[i_part]->face_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &pdm_part[i_part]->edge_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &pdm_part[i_part]->vtx_ln_to_gn,
                                      PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_bound_concat_get(pm,
                                   i_part,
                                   PDM_BOUND_TYPE_FACE,
                                   &pdm_part[i_part]->face_bound_idx,
                                   &pdm_part[i_part]->face_bound,
                                   &pdm_part[i_part]->face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_KEEP);

    PDM_part_mesh_bound_concat_get(pm,
                                   i_part,
                                   PDM_BOUND_TYPE_EDGE,
                                   &pdm_part[i_part]->edge_bound_idx,
                                   &pdm_part[i_part]->edge_bound,
                                   &pdm_part[i_part]->edge_bound_ln_to_gn,
                                   PDM_OWNERSHIP_KEEP);

    pdm_part[i_part]->face_group_idx           = NULL;
    pdm_part[i_part]->face_group               = NULL;
    pdm_part[i_part]->face_group_ln_to_gn      = NULL;

    pdm_part[i_part]->face_part_bound_proc_idx = NULL;
    pdm_part[i_part]->face_part_bound_part_idx = NULL;
    pdm_part[i_part]->face_part_bound          = NULL;

    pdm_part[i_part]->edge_part_bound_proc_idx = NULL;
    pdm_part[i_part]->edge_part_bound_part_idx = NULL;
    pdm_part[i_part]->edge_part_bound          = NULL;

    pdm_part[i_part]->vtx_part_bound_proc_idx  = NULL;
    pdm_part[i_part]->vtx_part_bound_part_idx  = NULL;
    pdm_part[i_part]->vtx_part_bound           = NULL;

    pdm_part[i_part]->face_join_ln_to_gn       = NULL;
    pdm_part[i_part]->face_join_idx            = NULL;
    pdm_part[i_part]->face_join                = NULL;
    pdm_part[i_part]->edge_join_idx            = NULL;
    pdm_part[i_part]->edge_join                = NULL;

    pdm_part[i_part]->elt_section_ln_to_gn     = NULL;

    pdm_part[i_part]->cell_tag                 = NULL;
    pdm_part[i_part]->face_tag                 = NULL;
    pdm_part[i_part]->edge_tag                 = NULL;
    pdm_part[i_part]->vtx_tag                  = NULL;

    pdm_part[i_part]->cell_weight              = NULL;
    pdm_part[i_part]->face_weight              = NULL;

    pdm_part[i_part]->cell_color               = NULL;
    pdm_part[i_part]->face_color               = NULL;
    pdm_part[i_part]->face_hp_color            = NULL;
    pdm_part[i_part]->edge_color               = NULL;
    pdm_part[i_part]->vtx_color                = NULL;
    pdm_part[i_part]->thread_color             = NULL;
    pdm_part[i_part]->hyperplane_color         = NULL;

    pdm_part[i_part]->vtx_ghost_information    = NULL;

    pdm_part[i_part]->new_to_old_order_cell    = NULL;
    pdm_part[i_part]->new_to_old_order_face    = NULL;
    pdm_part[i_part]->new_to_old_order_edge    = NULL;
    pdm_part[i_part]->new_to_old_order_vtx     = NULL;

    pdm_part[i_part]->subpartlayout            = NULL;
  }

  return pdm_part;
}


// /**
//  *
//  * \brief Map each pair of (join, opposite join) to a same global id and count
//  *        the total number of faces in this unified join. Return a distribution.
//  *        Arrays are allocated in this function.
//  *
//  * \param [in]   _multipart          multipart object
//  * \param [out]  join_to_ref_join    Unique join id associated to each join
//  *                                     (size = n_total_join)
//  * \param [out]  face_in_join_distri Distribution of join faces over the ref
//  *                                   join ids (size = n_unique_joins+1)
//  */
// static void
// _build_join_uface_distribution
// (
//  PDM_multipart_t   *multipart,
//  int              **join_to_ref_join,
//  int              **face_in_join_distri
// )
// {

//   int i_rank;
//   int n_rank;
//   PDM_MPI_Comm_rank(multipart->comm, &i_rank);
//   PDM_MPI_Comm_size(multipart->comm, &n_rank);

//   // PDM_printf("pdm::_build_join_uface_distribution\n");
//   int n_total_joins  = multipart->n_total_joins;
//   int n_unique_joins = n_total_joins/2;
//   *join_to_ref_join    = (int *) malloc(n_total_joins  * sizeof(int));
//   *face_in_join_distri = (int *) malloc((n_unique_joins+1) * sizeof(int));
//   int* _face_in_join_distri = *face_in_join_distri;
//   int* _join_to_ref_join    = *join_to_ref_join;

//   //Build join_to_ref_join : we want the join and opposite join to have the same shift index,
//   // so we take the smaller join global id as the reference
//   int ref_join_gid = 0;
//   for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
//   {
//     int opp_join = multipart->join_to_opposite[ijoin];
//     if (ijoin < opp_join)
//     {
//       _join_to_ref_join[ijoin] = ref_join_gid;
//       _join_to_ref_join[opp_join] = ref_join_gid;
//       ref_join_gid ++;
//     }
//   }
//   /*
//   PDM_printf("Join to reference join :");
//   for (int ijoin = 0; ijoin < n_total_joins; ijoin++)
//    PDM_printf(" %d ", _join_to_ref_join[ijoin]);
//   PDM_printf("\n");
//   */

//   //Count faces in joins
//   int *nb_face_in_joins = PDM_array_zeros_int(n_unique_joins);

//   for (int i_zone = 0; i_zone < multipart->n_zone; i_zone++){
//     for (int i_part = 0; i_part < multipart->n_part[i_zone]; i_part++){
//       int *pface_join_idx = multipart->pmeshes[i_zone].parts[i_part]->face_join_idx;
//       for (int ijoin=0; ijoin < multipart->pmeshes[i_zone].n_joins; ijoin ++){
//         int join_gid     = multipart->pmeshes[i_zone].joins_ids[ijoin];
//         int join_opp_gid = multipart->join_to_opposite[join_gid];
//         //Paired joins must be counted only once
//         if (join_gid < join_opp_gid)
//           nb_face_in_joins[_join_to_ref_join[join_gid]] += pface_join_idx[ijoin+1] - pface_join_idx[ijoin];
//       }
//     }
//   }
//   /*
//   PDM_printf("[%d] nb_face_joins : ", i_rank);
//   for (int i = 0; i < n_unique_joins ; i++)
//     PDM_printf(" %d ", nb_face_in_joins[i]);
//   PDM_printf("\n");
//   */

//   //Sum faces and build distribution
//   PDM_MPI_Allreduce(nb_face_in_joins, &_face_in_join_distri[1], n_unique_joins,
//                     PDM_MPI_INT, PDM_MPI_SUM, multipart->comm);

//   _face_in_join_distri[0] = 0;
//   PDM_array_accumulate_int(_face_in_join_distri, n_unique_joins+1);

//   /*
//   PDM_printf("[%d] _face_in_join_distri : ", i_rank);
//   for (int i = 0; i < n_unique_joins + 1; i++)
//     PDM_printf(" %d ", _face_in_join_distri[i]);
//   PDM_printf("\n");
//   */

//   free(nb_face_in_joins);
//   // PDM_printf("pdm::_build_join_uface_distribution end \n");
// }

// /**
//  *
//  * \brief Complete join data, which originally only contains face local id,
//  *        with the connecting data opp proc, opp part, opp face local id.
//  *
//  * \param [inout]   _multipart          multipart object
//  */
// static void
// _search_matching_joins
// (
//  PDM_multipart_t *multipart
// )
// {
//   int i_rank;
//   int n_rank;
//   PDM_MPI_Comm_rank(multipart->comm, &i_rank);
//   PDM_MPI_Comm_size(multipart->comm, &n_rank);

//   //Construction of (unique) join distribution
//   int *join_to_ref_join;
//   int *face_in_join_distri;
//   _build_join_uface_distribution(multipart, &join_to_ref_join, &face_in_join_distri);

//   //Count total nb of join_faces
//   int nb_of_joins = 0;
//   for (int izone = 0 ; izone < multipart->n_zone; izone ++) {
//     nb_of_joins += multipart->n_part[izone] * multipart->pmeshes[izone].n_joins;
//   }

//   // Prepare lntogn numbering and partitioned data
//   PDM_g_num_t **shifted_lntogn = (PDM_g_num_t **) malloc(nb_of_joins * sizeof(PDM_g_num_t*));
//   int              **part_data = (int **)         malloc(nb_of_joins * sizeof(int *));
//   int        *nb_face_per_join = (int *)          malloc(nb_of_joins * sizeof(int));

//   int ijoin_pos  = 0;
//   for (int izone = 0 ; izone < multipart->n_zone; izone ++) {
//     int n_join            = multipart->pmeshes[izone].n_joins;
//     _part_mesh_t _pmeshes = multipart->pmeshes[izone];
//     for (int i_part = 0; i_part < multipart->n_part[izone]; i_part++) {
//       int         *face_join_idx    = _pmeshes.parts[i_part]->face_join_idx;
//       int         *face_join        = _pmeshes.parts[i_part]->face_join;
//       PDM_g_num_t *face_join_lntogn = _pmeshes.parts[i_part]->face_join_ln_to_gn;
//       for (int ijoin = 0; ijoin < n_join; ijoin++) {
//         int join_size = face_join_idx[ijoin + 1] - face_join_idx[ijoin];
//         nb_face_per_join[ijoin_pos] = join_size;
//         PDM_g_num_t *shifted_lntogn_loc = (PDM_g_num_t *) malloc(join_size * sizeof(PDM_g_num_t));
//         int         *part_data_loc      = (int *)         malloc(3 * join_size * sizeof(int));
//         //Get shift value from join unique distribution
//         int join_gid    = multipart->pmeshes[izone].joins_ids[ijoin];
//         int shift_value = face_in_join_distri[join_to_ref_join[join_gid]];
//         int j = 0;
//         //Prepare partitioned data : (PL, i_rank, i_part)
//         for (int iface = face_join_idx[ijoin]; iface < face_join_idx[ijoin + 1]; iface ++) {
//           shifted_lntogn_loc[j] = (PDM_g_num_t) shift_value + face_join_lntogn[iface];
//           part_data_loc[3*j]    = face_join[4*iface];
//           part_data_loc[3*j+1]  = i_rank;
//           part_data_loc[3*j+2]  = i_part;
//           j++;
//         }
//         shifted_lntogn[ijoin_pos] = shifted_lntogn_loc;
//         part_data[ijoin_pos]      = part_data_loc;
//         ijoin_pos += 1;
//       }
//     }
//   }
//   /*
//   PDM_printf("[%d] nb_face_per_join : ", i_rank);
//   for (int i = 0; i < nb_of_joins; i++)
//     PDM_printf(" %d ", nb_face_per_join[i]);
//   PDM_printf("\n");
//   */

//   //Now exchange join information using part_to_block / block_to_part
//   PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
//                                                        PDM_PART_TO_BLOCK_POST_NOTHING,
//                                                        1.,
//                                                        shifted_lntogn,
//                                                        NULL,
//                                                        nb_face_per_join,
//                                                        nb_of_joins,
//                                                        multipart->comm);

//   PDM_g_num_t *distrib_index = PDM_part_to_block_distrib_index_get(ptb);

//   /*
//   PDM_printf("[%d] PTB distri : ", i_rank);
//   for (int i=0; i < n_rank + 1; i++)
//     PDM_printf(" %d ", distrib_index[i]);
//   PDM_printf("\n");
//   */

//   int         *block_data;
//   PDM_part_to_block_exch(ptb,
//                          sizeof(int),
//                          PDM_STRIDE_CST_INTERLACED,
//                          3,
//                          NULL,
//                          (void **) part_data,
//                          NULL,
//                          (void **) &block_data);

//   /*
//   int n_elt_block = PDM_part_to_block_n_elt_block_get(ptb);
//   PDM_printf("[%d] PTB nb_elem : %d\n", i_rank, n_elt_block);
//   if (i_rank == 1)
//   {
//     PDM_g_num_t *glob_num = PDM_part_to_block_block_gnum_get(ptb);
//     PDM_printf("[%d] PTB globnum : ", i_rank);
//     for (int i = 0; i < n_elt_block; i++)
//       printf(" %d ", glob_num[i]);
//     PDM_printf("\n");
//     PDM_printf("[%d] PTB data : ", i_rank);
//     for (int i = 0; i < n_elt_block; i++)
//       printf(" (%d %d %d) ", block_data[3*i],
//                              block_data[3*i+1],
//                              block_data[3*i+2]);
//     PDM_printf("\n");
//   }
//   */

//   // Don't free ptb now since we need the distribution and the block_data
//   PDM_block_to_part_t *btp = PDM_block_to_part_create(distrib_index,
//                                (const PDM_g_num_t **) shifted_lntogn,
//                                                       nb_face_per_join,
//                                                       nb_of_joins,
//                                                       multipart->comm);

//   int **new_part_data = (int **) malloc(nb_of_joins * sizeof(int *));
//   for (int ijoin = 0; ijoin < nb_of_joins; ijoin ++){
//     new_part_data[ijoin] = (int *) malloc(6*nb_face_per_join[ijoin]*sizeof(int));
//   }
//   int cst_stride = 6;

//   PDM_block_to_part_exch_in_place(btp,
//                          sizeof(int),
//                          PDM_STRIDE_CST_INTERLACED,
//                          &cst_stride,
//                          (void *) block_data,
//                          NULL,
//                          (void **) new_part_data);

//   free(block_data);
//   PDM_part_to_block_free(ptb);
//   PDM_block_to_part_free(btp);

//   /*
//   if (i_rank == 0)
//   {
//     PDM_printf("[%d] BTP data : \n",  i_rank);
//     for (int ijoin = 0; ijoin < nb_of_joins; ijoin++)
//     {
//       PDM_printf("  ijoin %d(%d) :", ijoin, nb_face_per_join[ijoin]);
//       for (int iface = 0; iface < nb_face_per_join[ijoin]; iface++)
//         PDM_printf(" (%d %d %d %d %d %d) ", new_part_data[ijoin][6*iface],
//                                             new_part_data[ijoin][6*iface+1],
//                                             new_part_data[ijoin][6*iface+2],
//                                             new_part_data[ijoin][6*iface+3],
//                                             new_part_data[ijoin][6*iface+4],
//                                             new_part_data[ijoin][6*iface+5]);
//       PDM_printf("\n");
//     }
//   }
//   */


//   //Process received data
//   ijoin_pos = 0;
//   for (int izone = 0 ; izone < multipart->n_zone; izone ++) {
//     int n_join = multipart->pmeshes[izone].n_joins;
//     _part_mesh_t _pmeshes = multipart->pmeshes[izone];
//     for (int i_part = 0; i_part < multipart->n_part[izone]; i_part++) {
//       int *face_join_idx = _pmeshes.parts[i_part]->face_join_idx;
//       int *face_join     = _pmeshes.parts[i_part]->face_join;
//       for (int ijoin = 0; ijoin < n_join; ijoin++) {
//         int join_size = face_join_idx[ijoin + 1] - face_join_idx[ijoin];
//         int *part_data_loc = new_part_data[ijoin_pos];
//         for (int i = 0; i < join_size; i++) {
//           int opp_proc = -1;
//           int opp_part = -1;
//           int opp_pl   = -1;
//           if (part_data_loc[6*i + 1] != i_rank)
//           {
//             opp_proc = part_data_loc[6*i + 1];
//             opp_part = part_data_loc[6*i + 2];
//             opp_pl   = part_data_loc[6*i + 0];
//           }
//           else if (part_data_loc[6*i + 4] != i_rank)
//           {
//             opp_proc = part_data_loc[6*i + 4];
//             opp_part = part_data_loc[6*i + 5];
//             opp_pl   = part_data_loc[6*i + 3];
//           }
//           // The two joins are on the same proc, look at the parts
//           else
//           {
//             opp_proc = i_rank;
//             if (part_data_loc[6*i + 2] != i_part)
//             {
//               opp_part = part_data_loc[6*i + 2];
//               opp_pl   = part_data_loc[6*i + 0];
//             }
//             else if (part_data_loc[6*i + 5] != i_part)
//             {
//               opp_part = part_data_loc[6*i + 5];
//               opp_pl   = part_data_loc[6*i + 3];
//             }
//             // The two joins have the same proc id / part id, we need to check original pl
//             else
//             {
//               opp_part = i_part;
//               int original_pl = face_join[4*(face_join_idx[ijoin] + i)];
//               if (part_data_loc[6*i] != original_pl)
//                 opp_pl = part_data_loc[6*i];
//               else
//                 opp_pl = part_data_loc[6*i+3];
//             }
//           }
//           //Fill values opp_proc, opp_part, opp_plvalue
//           face_join[4*(face_join_idx[ijoin] + i) + 1] = opp_proc;
//           face_join[4*(face_join_idx[ijoin] + i) + 2] = opp_part;
//           face_join[4*(face_join_idx[ijoin] + i) + 3] = opp_pl;
//         }
//         ijoin_pos += 1;
//       }
//     }
//   }

//   //Deallocate
//   for (int i = 0; i < nb_of_joins; i++) {
//     free(shifted_lntogn[i]);
//     free(part_data[i]);
//     free(new_part_data[i]);
//   }
//   free(shifted_lntogn);
//   free(part_data);
//   free(new_part_data);
//   free(join_to_ref_join);
//   free(face_in_join_distri);
//   free(nb_face_per_join);
// }

/**
 *
 * \brief Free the memory occuped by a partition structure
 *
 * \param [inout]   part          _part_t object
 */
static void
_part_free
(
 _part_t         *part
)
{
  /* Following is not results but internal array */
  if (part->new_to_old_order_cell != NULL)
    free(part->new_to_old_order_cell);
  part->new_to_old_order_cell = NULL;

  if (part->new_to_old_order_face != NULL)
    free(part->new_to_old_order_face);
  part->new_to_old_order_face = NULL;

  if (part->new_to_old_order_edge != NULL)
    free(part->new_to_old_order_edge);
  part->new_to_old_order_edge = NULL;


  if (part->new_to_old_order_vtx != NULL)
    free(part->new_to_old_order_vtx);
  part->new_to_old_order_vtx = NULL;

  if(part->subpartlayout != NULL){
    if(part->subpartlayout->cell_tile_idx!= NULL)
      free(part->subpartlayout->cell_tile_idx);
    if(part->subpartlayout->face_tile_idx!= NULL)
      free(part->subpartlayout->face_tile_idx);
    if(part->subpartlayout->face_bnd_tile_idx!= NULL)
      free(part->subpartlayout->face_bnd_tile_idx);
    if(part->subpartlayout->mask_tile_idx!= NULL)
      free(part->subpartlayout->mask_tile_idx);
    if(part->subpartlayout->cell_vect_tile_idx!= NULL)
      free(part->subpartlayout->cell_vect_tile_idx);
    if(part->subpartlayout->mask_tile_n!= NULL)
      free(part->subpartlayout->mask_tile_n);
    if(part->subpartlayout->cell_vect_tile_n!= NULL)
      free(part->subpartlayout->cell_vect_tile_n);
    if(part->subpartlayout->mask_tile!= NULL)
      free(part->subpartlayout->mask_tile);
    free(part->subpartlayout);
  }


  free(part);
}


static void
_setup_ghost_information
(
const int           i_rank,
const int           n_part,
const int          *pn_vtx,
      int         **pinternal_vtx_priority,
      PDM_g_num_t  *distrib_partition
)
{
  /* 0 : Interior / 1 : owner join (at least one) / 2 : not owner */
  // pinternal_vtx_priority contains value between i_rank + i_part
  for (int ipart = 0; ipart < n_part; ipart++) {
    for(int ivtx = 0; ivtx < pn_vtx[ipart]; ++ivtx) {

      int g_part = pinternal_vtx_priority[ipart][ivtx];
      int t_part = -1;
      if( g_part >= distrib_partition[i_rank] && g_part < distrib_partition[i_rank+1]) {
        t_part = g_part - distrib_partition[i_rank];
      }

      // if(pinternal_vtx_priority[ipart][ivtx] == i_rank){
      if(t_part == ipart){
        pinternal_vtx_priority[ipart][ivtx] = 1;
      } else if(pinternal_vtx_priority[ipart][ivtx] == -1) {
        pinternal_vtx_priority[ipart][ivtx] = 0;
      } else { /* Pas owner / pas interieur donc sur 1 autre proc */
         pinternal_vtx_priority[ipart][ivtx] = 2;
      }
    }
  }
}

static
void
_create_dparent_num_corner
(
  PDM_dmesh_nodal_elmts_t  *dmn_elts,
  PDM_g_num_t             **dparent_gnum_corner,
  PDM_g_num_t             **distrib_corner
)
{
  // Just repart all section in one block
  int n_rank = -1;
  PDM_MPI_Comm_size (dmn_elts->comm, &n_rank);

  int n_section = dmn_elts->n_section;

  int          *dn_corner       = malloc(n_section * sizeof(int          ));
  PDM_g_num_t **corner_ln_to_gn = malloc(n_section * sizeof(PDM_g_num_t *));
  PDM_g_num_t **corner_vtx      = malloc(n_section * sizeof(PDM_g_num_t *));
  int         **corner_vtx_n    = malloc(n_section * sizeof(int         *));

  for (int i_section = 0; i_section < n_section; i_section++) {

    int id_section = dmn_elts->sections_id[i_section];
    const PDM_g_num_t* distrib = PDM_DMesh_nodal_elmts_distrib_section_get(dmn_elts, id_section);

    PDM_g_num_t beg_elmt_gnum = distrib[dmn_elts->i_rank] + dmn_elts->section_distribution[i_section];

    PDM_Mesh_nodal_elt_t t_elt = PDM_DMesh_nodal_elmts_section_type_get(dmn_elts, id_section);

    assert(t_elt == PDM_MESH_NODAL_POINT);

    int n_elt           = PDM_DMesh_nodal_elmts_section_n_elt_get(dmn_elts, id_section);
    PDM_g_num_t* connec = PDM_DMesh_nodal_elmts_section_std_get(dmn_elts, id_section);

    dn_corner[i_section] = n_elt;
    corner_ln_to_gn[i_section] = malloc( n_elt * sizeof(PDM_g_num_t));
    corner_vtx_n   [i_section] = malloc( n_elt * sizeof(int        ));
    corner_vtx     [i_section] = connec;
    for(int i = 0; i < n_elt; ++i) {
      corner_vtx_n   [i_section][i] = 1;
      corner_ln_to_gn[i_section][i] = beg_elmt_gnum + i + 1;
    }

    // PDM_log_trace_array_long(connec, n_elt, "corner_vtx ::");
    // PDM_log_trace_array_long(corner_ln_to_gn[i_section], n_elt, "corner_ln_to_gn ::");

  }

  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      corner_ln_to_gn,
                                                      NULL,
                                                      dn_corner,
                                                      n_section,
                                                      dmn_elts->comm);


  PDM_g_num_t* distrib_ptb = PDM_part_to_block_distrib_index_get(ptb);

  PDM_g_num_t* _distrib_corner = malloc((n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    _distrib_corner[i] = distrib_ptb[i];
  }

  int         *blk_child_n    = NULL;
  PDM_g_num_t *blk_child_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         -1,
                         corner_vtx_n,
             (void **)   corner_vtx,
                        &blk_child_n,
             (void **)  &blk_child_gnum);

  *dparent_gnum_corner = blk_child_gnum;
  free(blk_child_n);

  PDM_part_to_block_free(ptb);
  for (int i_section = 0; i_section < n_section; i_section++) {
    free(corner_ln_to_gn[i_section]);
    free(corner_vtx_n[i_section]);
  }
  free(corner_ln_to_gn);
  free(corner_vtx_n);
  free(corner_vtx);
  free(dn_corner);

  *distrib_corner = _distrib_corner;
}

static PDM_part_mesh_nodal_t*
_compute_part_mesh_nodal_3d
(
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part,
 PDM_ownership_t    ownership
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  int           *pn_cell        = (int          * ) malloc( n_part * sizeof(int          ));
  int           *pn_face        = (int          * ) malloc( n_part * sizeof(int          ));
  int           *pn_vtx         = (int          * ) malloc( n_part * sizeof(int          ));
  PDM_g_num_t  **pcell_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  double       **pvtx_coord     = (double      ** ) malloc( n_part * sizeof(double      *));
  for(int i_part = 0; i_part < n_part; ++i_part){

    pn_cell[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_CELL  );
    pn_face[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_FACE  );
    pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_VERTEX);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_CELL,
                                      &pcell_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pface_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &pvtx_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_vtx_coord_get(pm->pmesh,
                                i_part,
                                &pvtx_coord[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_vol = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->volumic,
                                                                                        n_part,
                                                                                        pn_vtx,
                                                                                        pvtx_ln_to_gn,
                                                                                        pn_cell,
                                                                                        pcell_ln_to_gn,
                                                                                        NULL);

  int          *pn_surf             = NULL;
  PDM_g_num_t **psurf_gnum          = NULL;
  PDM_g_num_t **psurf_to_face_g_num = NULL;
  PDM_reverse_dparent_gnum(dmn->surfacic->dparent_gnum,
                           NULL, // dparent_sign
                           dmn->surfacic->delmt_child_distrib,
                           n_part,
                           pn_face,
                           pface_ln_to_gn,
                          &pn_surf,
                          &psurf_gnum,
                          &psurf_to_face_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  if(0 == 1) {
    PDM_log_trace_array_long(dmn->surfacic->dparent_gnum, dmn->surfacic->delmt_child_distrib[i_rank+1] - dmn->surfacic->delmt_child_distrib[i_rank], "dmn->surfacic->dparent_gnum : ");
    for(int i_part = 0; i_part < n_part; ++i_part){
      PDM_log_trace_array_long(pface_ln_to_gn     [i_part], pn_face[i_part], "pface_ln_to_gn      : ");
      PDM_log_trace_array_long(psurf_to_face_g_num[i_part], pn_surf[i_part], "psurf_to_face_g_num : ");
    }
  }

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_surf,
                                                                                         psurf_gnum,
                                                                                         psurf_to_face_g_num);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(psurf_gnum         [i_part]);
    free(psurf_to_face_g_num[i_part]);
  }
  free(pn_surf);
  free(psurf_gnum);
  free(psurf_to_face_g_num);
  free(pface_ln_to_gn);
  free(pn_face);

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = NULL;

  if(dmn->ridge != NULL) {

    PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
    int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_edge[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_EDGE);
      PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_EDGE,
                                        &pedge_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);
    }
    int          *pn_ridge             = NULL;
    PDM_g_num_t **pridge_gnum          = NULL;
    PDM_g_num_t **pridge_to_edge_g_num = NULL;
    PDM_reverse_dparent_gnum(dmn->ridge->dparent_gnum,
                             NULL, // dparent_sign
                             dmn->ridge->delmt_child_distrib,
                             n_part,
                             pn_edge,
                             pedge_ln_to_gn,
                            &pn_ridge,
                            &pridge_gnum,
                            &pridge_to_edge_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                               n_part,
                                                               pn_vtx,
                                                               pvtx_ln_to_gn,
                                                               pn_ridge,
                                                               pridge_gnum,
                                                               pridge_to_edge_g_num);

    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pridge_gnum[i_part]);
      free(pridge_to_edge_g_num[i_part]);
    }
    free(pn_ridge);
    free(pridge_gnum);
    free(pridge_to_edge_g_num);
    free(pedge_ln_to_gn);
    free(pn_edge);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_corner = NULL;

  if(dmn->corner != NULL) {

    PDM_g_num_t *dparent_gnum_corner = NULL;
    PDM_g_num_t *distrib_corner      = NULL;
    _create_dparent_num_corner(dmn->corner, &dparent_gnum_corner, &distrib_corner);

    int          *pn_corner             = NULL;
    PDM_g_num_t **pcorner_gnum          = NULL;
    PDM_g_num_t **pcorner_to_vtx_g_num  = NULL;
    PDM_reverse_dparent_gnum(dparent_gnum_corner,
                             NULL, // dparent_sign
                             distrib_corner,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                            &pn_corner,
                            &pcorner_gnum,
                            &pcorner_to_vtx_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_corner = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->corner,
                                                                n_part,
                                                                pn_vtx,
                                                                pvtx_ln_to_gn,
                                                                pn_corner,
                                                                pcorner_gnum,
                                                                pcorner_to_vtx_g_num);
    free(dparent_gnum_corner);
    free(distrib_corner);
    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pcorner_gnum[i_part]);
      free(pcorner_to_vtx_g_num[i_part]);
    }
    free(pn_corner);
    free(pcorner_gnum);
    free(pcorner_to_vtx_g_num);
  }

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_vol , ownership);
  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf, ownership);
  if(pmn_ridge != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge, ownership);
  }
  if(pmn_corner != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_corner, ownership);
  }
  for(int i_part = 0; i_part < n_part; ++i_part) {

    // Copy coordinates because ownership between part_mesh and part_mesh_nodal is complicated
    double      *lvtx_coords   = (double      *) malloc(3 * pn_vtx[i_part] * sizeof(double     ));
    PDM_g_num_t *lvtx_ln_to_gn = (PDM_g_num_t *) malloc(3 * pn_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < 3 * pn_vtx[i_part]; ++i_vtx) {
      lvtx_coords[i_vtx] = pvtx_coord[i_part][i_vtx];
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      lvtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_part][i_vtx];
    }
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  pn_vtx[i_part],
                                  lvtx_coords,
                                  lvtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
  }

  free(pcell_ln_to_gn);
  free(pn_cell);
  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_vtx);

  return pmn;
}

static PDM_part_mesh_nodal_t*
_compute_part_mesh_nodal_2d
(
 PDM_dmesh_nodal_t *dmn,
 _part_mesh_t      *pm,
 int                n_part,
 PDM_ownership_t    ownership
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  /*
   * Rebuild the volumic part from cell
   */
  PDM_g_num_t  **pface_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_face        = (int *  )         malloc( n_part * sizeof(int          ));

  PDM_g_num_t  **pvtx_ln_to_gn  = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_vtx         = (int *  )         malloc( n_part * sizeof(int          ));
  double       **pvtx_coord     = (double      ** ) malloc( n_part * sizeof(double      *));
  for(int i_part = 0; i_part < n_part; ++i_part){
    pn_face[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_FACE  );
    pn_vtx [i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_VERTEX);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_FACE,
                                      &pface_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      &pvtx_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);

    PDM_part_mesh_vtx_coord_get(pm->pmesh,
                                i_part,
                                &pvtx_coord[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }

  PDM_part_mesh_nodal_elmts_t* pmn_surf = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->surfacic,
                                                                                         n_part,
                                                                                         pn_vtx,
                                                                                         pvtx_ln_to_gn,
                                                                                         pn_face,
                                                                                         pface_ln_to_gn,
                                                                                         NULL);

  PDM_g_num_t  **pedge_ln_to_gn = (PDM_g_num_t ** ) malloc( n_part * sizeof(PDM_g_num_t *));
  int           *pn_edge        = (int *  )         malloc( n_part * sizeof(int          ));
  for(int i_part = 0; i_part < n_part; ++i_part) {
    pn_edge[i_part] = PDM_part_mesh_n_entity_get(pm->pmesh, i_part, PDM_MESH_ENTITY_EDGE);
    PDM_part_mesh_entity_ln_to_gn_get(pm->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      &pedge_ln_to_gn[i_part], PDM_OWNERSHIP_BAD_VALUE);
  }
  int          *pn_ridge;
  PDM_g_num_t **pridge_gnum;
  PDM_g_num_t **pridge_to_edge_g_num;
  PDM_reverse_dparent_gnum(dmn->ridge->dparent_gnum,
                           NULL, // dparent_sign
                           dmn->ridge->delmt_child_distrib,
                           n_part,
                           pn_edge,
                           pedge_ln_to_gn,
                          &pn_ridge,
                          &pridge_gnum,
                          &pridge_to_edge_g_num,
                           NULL, // pchild_parent_sign
                           dmn->comm);

  PDM_part_mesh_nodal_elmts_t* pmn_ridge = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->ridge,
                                                                                          n_part,
                                                                                          pn_vtx,
                                                                                          pvtx_ln_to_gn,
                                                                                          pn_ridge,
                                                                                          pridge_gnum,
                                                                                          pridge_to_edge_g_num);

  for(int i_part = 0; i_part < n_part; ++i_part){
    free(pridge_gnum[i_part]);
    free(pridge_to_edge_g_num[i_part]);
  }
  free(pn_ridge);
  free(pridge_gnum);
  free(pridge_to_edge_g_num);

  PDM_part_mesh_nodal_elmts_t* pmn_corner = NULL;

  if(dmn->corner != NULL) {

    PDM_g_num_t *dparent_gnum_corner = NULL;
    PDM_g_num_t *distrib_corner      = NULL;
    _create_dparent_num_corner(dmn->corner, &dparent_gnum_corner, &distrib_corner);

    int          *pn_corner             = NULL;
    PDM_g_num_t **pcorner_gnum          = NULL;
    PDM_g_num_t **pcorner_to_vtx_g_num  = NULL;
    PDM_reverse_dparent_gnum(dparent_gnum_corner,
                             NULL, // dparent_sign
                             distrib_corner,
                             n_part,
                             pn_vtx,
                             pvtx_ln_to_gn,
                            &pn_corner,
                            &pcorner_gnum,
                            &pcorner_to_vtx_g_num,
                             NULL, // pchild_parent_sign
                             dmn->comm);


    pmn_corner = PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmn->corner,
                                                                n_part,
                                                                pn_vtx,
                                                                pvtx_ln_to_gn,
                                                                pn_corner,
                                                                pcorner_gnum,
                                                                pcorner_to_vtx_g_num);
    free(dparent_gnum_corner);
    free(distrib_corner);
    for(int i_part = 0; i_part < n_part; ++i_part){
      free(pcorner_gnum[i_part]);
      free(pcorner_to_vtx_g_num[i_part]);
    }
    free(pn_corner);
    free(pcorner_gnum);
    free(pcorner_to_vtx_g_num);
  }

  /* Create top structure */
  PDM_part_mesh_nodal_t* pmn = PDM_part_mesh_nodal_create(dmn->mesh_dimension,
                                                          n_part,
                                                          dmn->comm);

  PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_surf , ownership);
  if(pmn_ridge != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_ridge, ownership);
  }
  if(pmn_corner != NULL) {
    PDM_part_mesh_nodal_add_part_mesh_nodal_elmts(pmn, pmn_corner, ownership);
  }
  // TO DO : corners?

  for(int i_part = 0; i_part < n_part; ++i_part) {
    // Copy coordinates because ownership between part_mesh and part_mesh_nodal is complicated
    double      *lvtx_coords   = (double      *) malloc(3 * pn_vtx[i_part] * sizeof(double     ));
    PDM_g_num_t *lvtx_ln_to_gn = (PDM_g_num_t *) malloc(    pn_vtx[i_part] * sizeof(PDM_g_num_t));
    for(int i_vtx = 0; i_vtx < 3 * pn_vtx[i_part]; ++i_vtx) {
      lvtx_coords[i_vtx] = pvtx_coord[i_part][i_vtx];
    }
    for(int i_vtx = 0; i_vtx < pn_vtx[i_part]; ++i_vtx) {
      lvtx_ln_to_gn[i_vtx] = pvtx_ln_to_gn[i_part][i_vtx];
    }
    PDM_part_mesh_nodal_coord_set(pmn,
                                  i_part,
                                  pn_vtx[i_part],
                                  lvtx_coords,
                                  lvtx_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);
    // free(lvtx_ln_to_gn);
  }

  free(pedge_ln_to_gn);
  free(pn_edge);
  free(pface_ln_to_gn);
  free(pn_face);

  free(pvtx_ln_to_gn);
  free(pvtx_coord);
  free(pn_vtx);
  return pmn;
}

static
PDM_g_num_t*
_split_graph
(
      PDM_MPI_Comm       comm,
      PDM_dmesh_t       *dmesh,
      _part_mesh_t      *pmeshes,
      int                n_part,
      PDM_split_dual_t   split_method,
      PDM_part_size_t    part_size_method,
const double            *part_fraction,
      PDM_g_num_t       *distrib_node,
      int              **node_part
)
{
  int verbose = 0;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Split only graph by the most high entity : (cell for 3D / face for 2D)
   *  The most high level entitvy is named elmt
   */
  int  dn_cell = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_CELL  );
  int  dn_face = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_FACE  );
  int  dn_edge = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_EDGE  );
  int  dn_vtx  = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VERTEX);

  if(verbose && i_rank == 0) {
    printf(" dn_cell = %i \n", dn_cell);
    printf(" dn_face = %i \n", dn_face);
    printf(" dn_edge = %i \n", dn_edge);
    printf(" dn_vtx  = %i \n", dn_vtx );
  }

  int         *darc_to_elmt_idx = NULL; // Donc face_cell OU edge_face
  PDM_g_num_t *darc_to_elmt_tmp = NULL;
  PDM_g_num_t *darc_to_elmt     = NULL;
  int         *delmt_to_arc_idx = NULL; // Donc cell_face OU face_edge
  PDM_g_num_t *delmt_to_arc     = NULL;
  int dn_node = 0;
  int dn_arc  = 0;

  if(dmesh->dn_cell == 0) { // Donc 2D
    dn_node = dmesh->dn_face;
    dn_arc  = dmesh->dn_edge;

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                               &darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                               &delmt_to_arc,
                               &delmt_to_arc_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
  } else {
    dn_node = dmesh->dn_cell;
    dn_arc  = dmesh->dn_face;
    assert(dmesh->dn_cell > 0);
    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_CELL,
                               &darc_to_elmt_tmp,
                               &darc_to_elmt_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               &delmt_to_arc,
                               &delmt_to_arc_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
  }
  assert(darc_to_elmt_idx == NULL);

  PDM_setup_connectivity_idx(dn_arc,
                             2,
                             darc_to_elmt_tmp,
                             &darc_to_elmt_idx,
                             &darc_to_elmt);
  PDM_g_num_t *distrib_arc  = PDM_compute_entity_distribution(comm, dn_arc );

  if(delmt_to_arc == NULL) {
    PDM_dconnectivity_transpose(comm,
                                distrib_arc,
                                distrib_node,
                                darc_to_elmt_idx,
                                darc_to_elmt,
                                0,
                                &delmt_to_arc_idx,
                                &delmt_to_arc);
    if(dmesh->dn_cell == 0) { // Donc 2D
      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 delmt_to_arc,
                                 delmt_to_arc_idx,
                                 PDM_OWNERSHIP_KEEP);
    } else {
      PDM_dmesh_connectivity_set(dmesh, PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 delmt_to_arc,
                                 delmt_to_arc_idx,
                                 PDM_OWNERSHIP_KEEP);
    }
  }

  PDM_g_num_t *dual_graph_idx = NULL;
  PDM_g_num_t *dual_graph     = NULL;
  PDM_deduce_combine_connectivity_dual(comm,
                                       distrib_node,
                                       distrib_arc,
                                       delmt_to_arc_idx,
                                       delmt_to_arc,
                                       darc_to_elmt_idx,
                                       darc_to_elmt,
                                       1, // is signed
                                       &dual_graph_idx,
                                       &dual_graph);
  free(darc_to_elmt_idx);
  free(darc_to_elmt);
  free(distrib_arc);

  /* Shift to 0 dual */
  for(int i = 0; i < dual_graph_idx[dn_node]; ++i) {
    dual_graph[i] = dual_graph[i] - 1;
  }

  // Compute total number of partitions for this zone
  int tn_part;
  PDM_MPI_Allreduce(&n_part, &tn_part, 1, PDM_MPI_INT, PDM_MPI_SUM, comm);
  pmeshes->tn_part = tn_part;

  PDM_g_num_t *distrib_partition = PDM_compute_entity_distribution(comm, n_part );
  double *part_fractions = NULL;
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS){
    int *n_part_per_rank = (int    *) malloc( n_rank * sizeof(int   ));
    int *displ           = (int    *) malloc( n_rank * sizeof(int   ));
    part_fractions       = (double *) malloc(tn_part * sizeof(double));
    for (int i =0; i < n_rank; i++){
      n_part_per_rank[i] = distrib_partition[i+1] - distrib_partition[i];
      displ[i] = distrib_partition[i];
    }

    PDM_MPI_Allgatherv((void*) part_fraction,
                       n_part,
                       PDM_MPI_DOUBLE,
                       part_fractions,
                       n_part_per_rank,
                       displ,
                       PDM_MPI_DOUBLE,
                       comm);
    free(n_part_per_rank);
    free(displ);
  }
  int *_node_part = malloc(dn_node * sizeof(int));
  if (split_method == PDM_SPLIT_DUAL_WITH_HILBERT) {

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);

    if(dmesh->dn_cell == 0) { // Donc 2D

      const int          *dedge_vtx_idx = NULL;
      const PDM_g_num_t  *dedge_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                (PDM_g_num_t **) &dedge_vtx,
                (int         **) &dedge_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

      PDM_g_num_t *edge_distri = PDM_compute_entity_distribution(comm, dn_edge);
      PDM_g_num_t *vtx_distri  = PDM_compute_entity_distribution(comm, dn_vtx );
      PDM_part_geom (PDM_PART_GEOM_HILBERT,
                     n_part,
                     comm,
                     dn_node,
                     delmt_to_arc_idx,
                     delmt_to_arc,
                     NULL, //cell_weight
                     dedge_vtx_idx,
                     dedge_vtx,
                     edge_distri,
                     dvtx_coord,
                     vtx_distri,
                     _node_part);

      free(edge_distri);
      free(vtx_distri);

    } else {
      int         *dface_vtx_idx = NULL;
      PDM_g_num_t *dface_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh, PDM_CONNECTIVITY_TYPE_FACE_VTX,
                (PDM_g_num_t **) &dface_vtx,
                (int         **) &dface_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

      PDM_g_num_t *face_distri = PDM_compute_entity_distribution(comm, dn_face);
      PDM_g_num_t *vtx_distri  = PDM_compute_entity_distribution(comm, dn_vtx );
      PDM_part_geom (PDM_PART_GEOM_HILBERT,
                     n_part,
                     comm,
                     dn_node,
                     delmt_to_arc_idx,
                     delmt_to_arc,
                     NULL, //cell_weight
                     dface_vtx_idx,
                     dface_vtx,
                     face_distri,
                     dvtx_coord,
                     vtx_distri,
                     _node_part);

      free(face_distri);
      free(vtx_distri);
    }

  } else {
    PDM_para_graph_split (split_method,
                          distrib_node,
                          dual_graph_idx,
                          dual_graph,
                          NULL,
                          NULL,
                          tn_part,
                          part_fractions,
                          _node_part,
                          comm);
  }

  // PDM_log_trace_array_int (_node_part, dn_node, "_node_part :: ");

  free(dual_graph_idx);
  free(dual_graph);
  if (part_size_method == PDM_PART_SIZE_HETEROGENEOUS) {
    free(part_fractions);
  }

  *node_part = _node_part;

  return distrib_partition;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Build a multipart structure
 *
 * \param [in]   n_zone       Number of zones in the original mesh
 * \param [in]   n_part       Number of partition per proc in each zone
 * \param [in]   merge_blocks Merge or not the zones before splitting
 * \param [in]   split_method Choice of library used to split the mesh
 * \param [in]   part_size_method Choice of homogeneous or heterogeneous partitions
 * \param [in]   part_weight  Weight (in %) of each partition in heterogeneous case
 * \param [in]   comm         PDM_MPI communicator
 *
 * \return     Pointer to a new \ref PDM_multipart_t object
 */

PDM_multipart_t *
PDM_multipart_create
(
 const int              n_zone,
 const int             *n_part,
 const PDM_bool_t       merge_blocks,
 const PDM_split_dual_t split_method,
 const PDM_part_size_t  part_size_method,
 const double          *part_fraction,
 const PDM_MPI_Comm     comm,
 const PDM_ownership_t  owner
)
{
  PDM_multipart_t *multipart = (PDM_multipart_t *) malloc(sizeof(PDM_multipart_t));

  multipart->n_zone           = n_zone;
  multipart->n_part           = (int * ) malloc( multipart->n_zone * sizeof(int));

  for (int i = 0; i < multipart->n_zone; ++i) {
    multipart->n_part[i] = n_part[i];
  }

  multipart->merge_blocks     = merge_blocks;
  multipart->split_method     = split_method;
  multipart->part_size_method = part_size_method;
  multipart->part_fraction    = part_fraction;
  multipart->comm             = comm;
  multipart->owner            = owner;

  multipart->n_total_joins    = 0;
  multipart->join_to_opposite = NULL;

  // multipart->dmeshes_ids = (int *) malloc(multipart->n_zone * sizeof(int));

  multipart->dmeshes       = (PDM_dmesh_t                **) malloc(multipart->n_zone * sizeof(PDM_dmesh_t                *));
  multipart->dmeshes_nodal = (PDM_dmesh_nodal_t          **) malloc(multipart->n_zone * sizeof(PDM_dmesh_nodal_t          *));
  multipart->dmn_to_dm     = (PDM_dmesh_nodal_to_dmesh_t **) malloc(multipart->n_zone * sizeof(PDM_dmesh_nodal_to_dmesh_t *));
  for (int izone = 0; izone < multipart->n_zone; ++izone) {
    multipart->dmeshes_nodal[izone] = NULL;
    multipart->dmeshes      [izone] = NULL;
    multipart->dmn_to_dm    [izone] = NULL;
  }

  multipart->pmeshes       = (_part_mesh_t *) malloc(multipart->n_zone * sizeof(_part_mesh_t));

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get("PDM_PART_RENUM_CELL_NONE");
  int _renum_face_method = PDM_part_renum_method_face_idx_get("PDM_PART_RENUM_FACE_NONE");
  int _renum_edge_method = PDM_part_renum_method_edge_idx_get("PDM_PART_RENUM_EDGE_NONE");
  int _renum_vtx_method  = PDM_part_renum_method_vtx_idx_get ("PDM_PART_RENUM_VTX_NONE" );
  for (int izone = 0; izone < multipart->n_zone; izone++) {
    multipart->pmeshes[izone].renum_cell_method = _renum_cell_method;
    multipart->pmeshes[izone].renum_face_method = _renum_face_method;
    multipart->pmeshes[izone].renum_edge_method = _renum_edge_method;
    multipart->pmeshes[izone].renum_vtx_method  = _renum_vtx_method;
    multipart->pmeshes[izone].renum_cell_properties = NULL;
    multipart->pmeshes[izone].joins_ids = NULL;
    multipart->pmeshes[izone].pmesh     = PDM_part_mesh_create(n_part[izone], comm);
    multipart->pmeshes[izone].vtx_ghost_information = malloc(n_part[izone] * sizeof(int *));
    multipart->pmeshes[izone].hyperplane_color      = malloc(n_part[izone] * sizeof(int *));
    multipart->pmeshes[izone].thread_color          = malloc(n_part[izone] * sizeof(int *));
    for(int i_part = 0; i_part < n_part[izone]; ++i_part) {
      multipart->pmeshes[izone].vtx_ghost_information[i_part] = NULL;
      multipart->pmeshes[izone].hyperplane_color     [i_part] = NULL;
      multipart->pmeshes[izone].thread_color         [i_part] = NULL;
    }
    multipart->pmeshes[izone].is_owner_vtx_ghost_information = PDM_TRUE;
    multipart->pmeshes[izone].is_owner_hyperplane_color      = PDM_TRUE;
    multipart->pmeshes[izone].is_owner_thread_color          = PDM_TRUE;
  }

  return (PDM_multipart_t *) multipart;
}

/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */
void PDM_multipart_register_block
(
 PDM_multipart_t   *multipart,
 const int          zone_id,
       PDM_dmesh_t *dmesh
)
{
  assert(zone_id < multipart->n_zone);
  multipart->dmeshes[zone_id] = dmesh;
}

/**
 *
 * \brief Set distributed mesh data for the input zone
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   zone_id        Global zone id
 * \param [in]   dmesh_id       Id of the distributed mesh structure to use
 */
void PDM_multipart_register_dmesh_nodal
(
 PDM_multipart_t         *multipart,
 const int                zone_id,
       PDM_dmesh_nodal_t *dmesh_nodal
)
{
  assert(zone_id < multipart->n_zone);
  assert(multipart->dmeshes_nodal[zone_id] == NULL);
  multipart->dmeshes_nodal[zone_id] = dmesh_nodal;
}

/**
 * \brief Set number of element in the block entity
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  entity_type           Type of entity (can be cell/face/edge/vtx)
 * \param [in]  dn_entity             Distributed number of entity in current process
 *
 */
void
PDM_multipart_dn_entity_set
(
       PDM_multipart_t     *multipart,
 const int                  i_zone,
       PDM_mesh_entities_t  entity_type,
       int                  dn_entity
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(i_zone);
  PDM_UNUSED(entity_type);
  PDM_UNUSED(dn_entity);
  abort();
}

/**
 * \brief Set number connectivity for current block
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  connectivity_type     Type of connectivity
 * \param [in]  connect               connectivity (size = connect_idx[dn_entity] )
 * \param [in]  connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity )
 *
 */
void
PDM_multipart_dconnectivity_set
(
       PDM_multipart_t         *multipart,
 const int                      i_zone,
       PDM_connectivity_type_t  connectivity_type,
       PDM_g_num_t             *dconnect,
       int                     *dconnect_idx
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(i_zone);
  PDM_UNUSED(connectivity_type);
  PDM_UNUSED(dconnect);
  PDM_UNUSED(dconnect_idx);
  abort();
}

/**
 * \brief Set group connectivity by kind
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  bound_type            Type of bound
 * \param [in]  connect               connectivity (size = connect_idx[dn_entity] )
 * \param [in]  connect_idx           Index of connectivity or NULL if face_cell for example  (size = dn_entity )
 *
 */
void
PDM_multipart_dgroup_set
(
       PDM_multipart_t          *multipart,
 const int                       i_zone,
       PDM_bound_type_t          bound_type,
       PDM_g_num_t              *dconnect,
       int                      *dconnect_idx
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(i_zone);
  PDM_UNUSED(bound_type);
  PDM_UNUSED(dconnect);
  PDM_UNUSED(dconnect_idx);
  abort();
}

/**
 * \brief Set group connectivity by kind
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [in]  dvtx_coord            Mesh coordinates (size = 3 * dn_vtx)
 */
void
PDM_multipart_dvtx_coord_set
(
       PDM_multipart_t *multipart,
 const int              i_zone,
 const double          *dvtx_coord
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(i_zone);
  PDM_UNUSED(dvtx_coord);
  abort();

}

void
PDM_multipart_domain_interface_shared_set
(
  PDM_multipart_t        *multipart,
  PDM_domain_interface_t *ditrf
)
{
  PDM_UNUSED(multipart);
  PDM_UNUSED(ditrf);
  abort();
}


/**
 *
 * \brief Set connecting data between all the zones
 *
 * \param [in]   multipart         Pointer to \ref PDM_multipart_t object
 * \param [in]   n_total_joins     Total number of interfaces
 * \param [in]   join_to_opposite  For each global join id, give the global id
 *                                   of the opposite join (size = n_total_joins)
 *
 * \note Join global id numbering must start at 0 and be continuous.
 */
void PDM_multipart_register_joins
(
 PDM_multipart_t *multipart,
 const int        n_total_joins,
 const int       *join_to_opposite
)
{

  multipart->n_total_joins    = n_total_joins;
  multipart->join_to_opposite = join_to_opposite;
}

/**
 *
 * \brief Set the reordering methods to be used after partitioning
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone                Id of zone which parameters apply (or -1 for all zones)
 * \param [in]   renum_cell_method     Choice of renumbering method for cells
 * \param [in]   renum_cell_properties Parameters used by cacheblocking method :
 *                                     [n_cell_per_cache_wanted, is_asynchrone, is_vectorisation,
                                        n_vect_face, split_method]
 * \param [in]   renum_face_method     Choice of renumbering method for faces
 *
 */
void PDM_multipart_set_reordering_options
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const char      *renum_cell_method,
 const int       *renum_cell_properties,
 const char      *renum_face_method
)
{

  int _renum_cell_method = PDM_part_renum_method_cell_idx_get(renum_cell_method);
  int _renum_face_method = PDM_part_renum_method_face_idx_get(renum_face_method);
  if (_renum_cell_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering cell method\n", renum_cell_method);
  }
  if (_renum_face_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering face method\n", renum_face_method);
  }

  if(i_zone < 0) {
    for (int izone = 0; izone < multipart->n_zone; izone++) {
      multipart->pmeshes[izone].renum_cell_method = _renum_cell_method;
      multipart->pmeshes[izone].renum_face_method = _renum_face_method;
      multipart->pmeshes[izone].renum_cell_properties = renum_cell_properties;
    }
  }
  else {
    assert(i_zone < multipart->n_zone);
    multipart->pmeshes[i_zone].renum_cell_method = _renum_cell_method;
    multipart->pmeshes[i_zone].renum_face_method = _renum_face_method;
    multipart->pmeshes[i_zone].renum_cell_properties = renum_cell_properties;
  }
}
void PDM_multipart_set_reordering_options_vtx
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 const char      *renum_vtx_method
)
{

  int _renum_vtx_method = PDM_part_renum_method_vtx_idx_get(renum_vtx_method);
  if (_renum_vtx_method == -1) {
    PDM_error (__FILE__, __LINE__, 0, "'%s' is an unknown renumbering vtx method\n", renum_vtx_method);
  }

  if(i_zone < 0) {
    for (int izone = 0; izone < multipart->n_zone; izone++) {
      multipart->pmeshes[izone].renum_vtx_method = _renum_vtx_method;
    }
  }
  else {
    assert(i_zone < multipart->n_zone);
    multipart->pmeshes[i_zone].renum_vtx_method = _renum_vtx_method;
  }
}


static
void
_run_ppart_zone2
(
PDM_dmesh_t       *dmesh,
PDM_dmesh_nodal_t *dmesh_nodal,
_part_mesh_t      *pmeshes,
int                n_part,
PDM_split_dual_t   split_method,
PDM_part_size_t    part_size_method,
const double*      part_fraction,
PDM_MPI_Comm       comm
)
{

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // int  dn_cell = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_CELL  );
  int  dn_face = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_FACE  );
  int  dn_edge = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_EDGE  );
  int  dn_vtx  = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VERTEX);

  int dn_node = 0;
  if(dmesh->dn_cell == 0) { // Donc 2D
    dn_node = dmesh->dn_face;
  } else {
    assert(dmesh->dn_cell > 0);
    dn_node = dmesh->dn_cell;
  }
  PDM_g_num_t *distrib_node = PDM_compute_entity_distribution(comm, dn_node);

  /*
   *  Split graph (manage 3D/2D automaticaly)
   */
  int *node_part = NULL;
  PDM_g_num_t* distrib_partition = _split_graph(comm,
                                                dmesh,
                                                pmeshes,
                                                n_part,
                                                split_method,
                                                part_size_method,
                                                part_fraction,
                                                distrib_node,
                                                &node_part);

  /*
   * Deduce node_ln_to_gn
   */
  int  *pn_node;
  PDM_g_num_t **pnode_ln_to_gn;
  PDM_part_assemble_partitions(comm,
                               distrib_partition,
                               distrib_node,
                               node_part,
                               NULL,
                               NULL,
                              &pn_node,
                              &pnode_ln_to_gn,
                               NULL);


  free(node_part);

  if(0 == 1) {
    for(int i_part = 0; i_part < n_part; ++i_part) {
      PDM_log_trace_array_long(pnode_ln_to_gn[i_part], pn_node[i_part], "pnode_ln_to_gn :: ");
    }
  }

  /*
   *  Deduce all required connectivity by descending connectivity
   */
  int          *pn_cell        = NULL;
  PDM_g_num_t **pcell_ln_to_gn = NULL;

  int          *pn_face        = NULL;
  PDM_g_num_t **pface_ln_to_gn = NULL;

  int         **pcell_face_idx = NULL;
  int         **pcell_face     = NULL;

  PDM_g_num_t *face_distrib    = NULL;

  if(dmesh->dn_cell != 0) { // Donc 3D
    PDM_g_num_t *cell_distri    = distrib_node;
    PDM_g_num_t *dcell_face     = NULL;
    int         *dcell_face_idx = NULL;
    PDM_dmesh_connectivity_get(dmesh,
                               PDM_CONNECTIVITY_TYPE_CELL_FACE,
                               &dcell_face,
                               &dcell_face_idx,
                               PDM_OWNERSHIP_BAD_VALUE);

    pn_cell        = pn_node;
    pcell_ln_to_gn = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;

    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 cell_distri,
                                                 dcell_face_idx,
                                                 dcell_face,
                                                 n_part,
                                                 pn_cell,
                          (const PDM_g_num_t **) pcell_ln_to_gn,
                                                &pn_face,
                                                &pface_ln_to_gn,
                                                &pcell_face_idx,
                                                &pcell_face);

  } else {
    face_distrib   = distrib_node;
    pn_face        = pn_node;
    pface_ln_to_gn = pnode_ln_to_gn;
    pn_node        = NULL;
    pnode_ln_to_gn = NULL;
  }


  // Fill _part_t structures with temporary arrays
  for (int i_part = 0; i_part < n_part; i_part++) {
    if(pn_cell != NULL) {

      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_CELL, pn_cell[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                     pcell_face    [i_part],
                                     pcell_face_idx[i_part],
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        pcell_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
      // PDM_log_trace_array_long(pcell_ln_to_gn[i_part], pn_cell       [i_part]             , "(in part ) cell_ln_to_gn ::");
      // PDM_log_trace_array_int (pcell_face_idx[i_part], pn_cell[i_part]+1             , "(in part ) cell_face_idx ::");
      // PDM_log_trace_array_int (pcell_face    [i_part], pcell_face_idx[i_part][pn_cell[i_part]], "(in part ) cell_face ::");
    }
  }


  int from_face_edge = 0;
  int from_face_vtx  = 0;
  // face vtx
  PDM_g_num_t *dface_vtx     = NULL;
  int         *dface_vtx_idx = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_VTX,
                             &dface_vtx,
                             &dface_vtx_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_vtx_idx != NULL){
    from_face_vtx = 1;
  }

  // face edge
  PDM_g_num_t *dface_edge     = NULL;
  int         *dface_edge_idx = NULL;
  PDM_dmesh_connectivity_get(dmesh,
                             PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                             &dface_edge,
                             &dface_edge_idx,
                             PDM_OWNERSHIP_BAD_VALUE);

  if(dface_edge_idx != NULL) {
    from_face_edge = 1;
  }

  // PDM_log_trace_connectivity_long(dface_edge_idx, dface_edge, dmesh->dn_face, "dface_edge ::");

  PDM_g_num_t *edge_distrib = NULL;
  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &edge_distrib);

  int own_edge_distrib = 0;
  if(edge_distrib == NULL) {
    own_edge_distrib = 1;
    edge_distrib = PDM_compute_entity_distribution(comm, dn_edge);
  }

  PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &face_distrib);

  int own_face_distrib = 0;
  if(face_distrib == NULL) {
    own_face_distrib = 1;
    face_distrib = PDM_compute_entity_distribution(comm, dn_face);
  }

  int          *pn_edge        = NULL;
  PDM_g_num_t **pedge_ln_to_gn = NULL;
  int         **pedge_vtx_idx  = NULL;
  int         **pedge_vtx      = NULL;

  int         **pface_edge_idx = NULL;
  int         **pface_edge     = NULL;

  int         **pface_vtx_idx = NULL;
  int         **pface_vtx     = NULL;

  int          *pn_vtx        = NULL;
  PDM_g_num_t **pvtx_ln_to_gn = NULL;
  if(from_face_edge == 1) {

    PDM_log_trace_array_long(face_distrib,  n_rank+1, "face_distribb ::");
    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 face_distrib,
                                                 dface_edge_idx,
                                                 dface_edge,
                                                 n_part,
                                                 pn_face,
                          (const PDM_g_num_t **) pface_ln_to_gn,
                                                 &pn_edge,
                                                 &pedge_ln_to_gn,
                                                 &pface_edge_idx,
                                                 &pface_edge);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_FACE, pn_face[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                     pface_edge    [i_part],
                                     pface_edge_idx[i_part],
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        pface_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);

      // PDM_log_trace_array_long(pface_ln_to_gn[i_part], pn_face       [i_part]             , "(in part ) face_ln_to_gn ::");
      // PDM_log_trace_array_int (pface_edge_idx[i_part], pn_face[i_part]+1             , "(in part ) face_edge_idx ::");
      // PDM_log_trace_array_int (pface_edge    [i_part], pface_edge_idx[i_part][pn_face[i_part]], "(in part ) face_edge ::");

      // PDM_log_trace_part_connectivity_gnum(pmeshes->parts[i_part]->face_edge_idx,
      //                                      pmeshes->parts[i_part]->face_edge,
      //                                      pface_ln_to_gn[i_part],
      //                                      pedge_ln_to_gn[i_part],
      //                                      pn_face[i_part],
      //                                      "face_edge_gnum");

    }

    // edge_vtx
    PDM_g_num_t *dedge_vtx     = NULL;
    int         *dedge_vtx_idx = NULL;
    PDM_dmesh_connectivity_get(dmesh,
                               PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                               &dedge_vtx,
                               &dedge_vtx_idx,
                               PDM_OWNERSHIP_BAD_VALUE);
    int *_dedge_vtx_idx = NULL;
    if(dedge_vtx_idx == NULL)  {
      _dedge_vtx_idx = malloc( (dn_edge+1) * sizeof(int));
      for(int i_edge = 0; i_edge < dn_edge+1; ++i_edge) {
        _dedge_vtx_idx[i_edge] = 2*i_edge;
      }
    } else {
      _dedge_vtx_idx = dedge_vtx_idx;
    }

    PDM_log_trace_connectivity_long(dedge_vtx_idx, dedge_vtx, dmesh->dn_edge, "dedge_vtx ::");
    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 edge_distrib,
                                                 _dedge_vtx_idx,
                                                 dedge_vtx,
                                                 n_part,
                                                 pn_edge,
                           (const PDM_g_num_t **) pedge_ln_to_gn,
                                                 &pn_vtx,
                                                 &pvtx_ln_to_gn,
                                                 &pedge_vtx_idx,
                                                 &pedge_vtx);
    if(dedge_vtx_idx == NULL)  {
      free(_dedge_vtx_idx);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_EDGE, pn_edge[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                     pedge_vtx     [i_part],
                                     NULL,
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_EDGE,
                                        pedge_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
      free(pedge_vtx_idx [i_part]);
    }
  } else if(from_face_vtx == 1) {
    // PDM_log_trace_connectivity_long(dface_vtx_idx, dface_vtx, dmesh->dn_face, "dface_vtx ::");
    PDM_part_dconnectivity_to_pconnectivity_sort(comm,
                                                 face_distrib,
                                                 dface_vtx_idx,
                                                 dface_vtx,
                                                 n_part,
                                                 pn_face,
                          (const PDM_g_num_t **) pface_ln_to_gn,
                                                 &pn_vtx,
                                                 &pvtx_ln_to_gn,
                                                 &pface_vtx_idx,
                                                 &pface_vtx);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_FACE, pn_face[i_part]);
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                     pface_vtx     [i_part],
                                     pface_vtx_idx [i_part],
                                     PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        pface_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
    }
  }

  // Vertex
  PDM_g_num_t *vtx_distrib = PDM_compute_entity_distribution(comm, dn_vtx);

  double *dvtx_coord = NULL;
  PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_BAD_VALUE);
  double      **pvtx_coord = NULL;
  PDM_part_dcoordinates_to_pcoordinates(comm,
                                        n_part,
                                        vtx_distrib,
                                        dvtx_coord,
                                        pn_vtx,
                 (const PDM_g_num_t **) pvtx_ln_to_gn,
                                       &pvtx_coord);


  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_VERTEX, pn_vtx[i_part]);
    PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_MESH_ENTITY_VERTEX,
                                      pvtx_ln_to_gn[i_part],
                                      PDM_OWNERSHIP_KEEP);
    PDM_part_mesh_vtx_coord_set(pmeshes->pmesh,
                                i_part,
                                pvtx_coord   [i_part],
                                PDM_OWNERSHIP_KEEP);

    // PDM_log_trace_array_long(pvtx_ln_to_gn[i_part], pn_vtx[i_part], "pvtx_ln_to_gn ::");
  }

  /*
   * Toutes les informations sont transferer aux partitions
   * Maintenant on souahite faire des post_traitement (genre des transpose poour fitter aux codes)
   * Egalement tout renumerotés !!
   */
  if(pn_cell != NULL){
    int **pface_cell = NULL;
    PDM_part_reverse_pcellface(n_part,
                               pn_cell,
                               pn_face,
                (const int **) pcell_face_idx,
                (const int **) pcell_face,
                              &pface_cell);
    // PDM_part_reorient_bound_faces(n_part,
    //                               pn_face,
    //                               pface_cell,
    //                (const int **) pcell_face_idx,
    //                               pcell_face,
    //                               NULL, // pface_vtx_idx
    //                               NULL, // pface_vtx
    //                               pface_edge_idx,  //
    //                               pface_edge); //

    for (int i_part = 0; i_part < n_part; i_part++) {
      // pmeshes->parts[i_part]->face_cell    = pface_cell   [i_part];
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                     pface_cell     [i_part],
                                     NULL,
                                     PDM_OWNERSHIP_KEEP);
    }

    free(pface_cell);
  } else {
    int **pedge_face = NULL;
    PDM_part_reverse_pcellface(n_part,
                               pn_face,
                               pn_edge,
                (const int **) pface_edge_idx,
                (const int **) pface_edge,
                              &pedge_face);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_connectivity_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_CONNECTIVITY_TYPE_EDGE_FACE,
                                     pedge_face     [i_part],
                                     NULL,
                                     PDM_OWNERSHIP_KEEP);
    }
    free(pedge_face);
  }

  /*
   * Force // ordering of vertex (needed by other ordering method)
   */
  int         **pinternal_vtx_bound_proc_idx  = NULL;
  int         **pinternal_vtx_bound_part_idx  = NULL;
  int         **pinternal_vtx_bound           = NULL;
  int         **pinternal_vtx_priority        = NULL;
  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      vtx_distrib,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                     &pinternal_vtx_priority);

  _setup_ghost_information(i_rank,
                           n_part,
                           pn_vtx,
                           pinternal_vtx_priority,
                           distrib_partition);

  /* Free in order to be correclty */
  for (int ipart = 0; ipart < n_part; ipart++) {
    free(pinternal_vtx_bound_proc_idx[ipart]);
    free(pinternal_vtx_bound_part_idx[ipart]);
    free(pinternal_vtx_bound[ipart]);
    // free(pinternal_vtx_priority[ipart]);
  }
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);

  /*
   * Group can be necessary for renumumbering
   */
  if(pn_cell != NULL){
    PDM_g_num_t *dface_bound = NULL;
    int         *dface_bound_idx = NULL;
    int n_face_group = PDM_dmesh_bound_get(dmesh,
                                           PDM_BOUND_TYPE_FACE,
                                           &dface_bound,
                                           &dface_bound_idx,
                                           PDM_OWNERSHIP_BAD_VALUE);

    int         **pface_bound_idx               = NULL;
    int         **pface_bound                   = NULL;
    PDM_g_num_t **pface_bound_ln_to_gn          = NULL;
    PDM_part_distgroup_to_partgroup(comm,
                                    face_distrib,
                                    n_face_group,
                                    dface_bound_idx,
                                    dface_bound,
                                    n_part,
                                    pn_face,
             (const PDM_g_num_t **) pface_ln_to_gn,
                                   &pface_bound_idx,
                                   &pface_bound,
                                   &pface_bound_ln_to_gn);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_bound_concat_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_BOUND_TYPE_FACE,
                                     n_face_group,
                                     pface_bound_idx     [i_part],
                                     pface_bound         [i_part],
                                     pface_bound_ln_to_gn[i_part],
                                     PDM_OWNERSHIP_KEEP);
    }

    free(pface_bound_idx     );
    free(pface_bound         );
    free(pface_bound_ln_to_gn);
  } else {

    PDM_g_num_t *dedge_bound = NULL;
    int         *dedge_bound_idx = NULL;
    int n_edge_group = PDM_dmesh_bound_get(dmesh,
                                           PDM_BOUND_TYPE_EDGE,
                                           &dedge_bound,
                                           &dedge_bound_idx,
                                           PDM_OWNERSHIP_BAD_VALUE);


    int         **pedge_bound_idx               = NULL;
    int         **pedge_bound                   = NULL;
    PDM_g_num_t **pedge_bound_ln_to_gn          = NULL;
    PDM_part_distgroup_to_partgroup(comm,
                                    edge_distrib,
                                    n_edge_group,
                                    dedge_bound_idx,
                                    dedge_bound,
                                    n_part,
                                    pn_edge,
             (const PDM_g_num_t **) pedge_ln_to_gn,
                                   &pedge_bound_idx,
                                   &pedge_bound,
                                   &pedge_bound_ln_to_gn);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_bound_concat_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_BOUND_TYPE_EDGE,
                                     n_edge_group,
                                     pedge_bound_idx     [i_part],
                                     pedge_bound         [i_part],
                                     pedge_bound_ln_to_gn[i_part],
                                     PDM_OWNERSHIP_KEEP);
    }

    free(pedge_bound_idx               );
    free(pedge_bound                   );
    free(pedge_bound_ln_to_gn          );

  }


  /* MAP on _part_t to reuse ordering */
  _part_t** parts = _map_part_t_with_part_mesh(pmeshes->pmesh);
  for (int i_part = 0; i_part < n_part; i_part++) {
    parts[i_part]->vtx_ghost_information = pinternal_vtx_priority[i_part];
    pmeshes->vtx_ghost_information[i_part] = pinternal_vtx_priority[i_part];
  }

  /*
   * Real re-numebering
   */
  PDM_part_renum_cell(parts, n_part, pmeshes->renum_cell_method, (void *) pmeshes->renum_cell_properties);
  PDM_part_renum_face(parts, n_part, pmeshes->renum_face_method, NULL);
  PDM_part_renum_edge(parts, n_part, pmeshes->renum_edge_method, NULL);
  PDM_part_renum_vtx (parts, n_part, pmeshes->renum_vtx_method , (void *) pinternal_vtx_priority);
  free(pinternal_vtx_priority);

  for (int i_part = 0; i_part < n_part; i_part++) {

    if(parts[i_part]->cell_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_CELL,
                                     parts[i_part]->cell_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->face_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_FACE,
                                     parts[i_part]->face_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->edge_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_EDGE,
                                     parts[i_part]->edge_color,
                                     PDM_OWNERSHIP_KEEP);
    }
    if (parts[i_part]->vtx_color != NULL) {
      PDM_part_mesh_entity_color_set(pmeshes->pmesh,
                                     i_part,
                                     PDM_MESH_ENTITY_VERTEX,
                                     parts[i_part]->vtx_color,
                                     PDM_OWNERSHIP_KEEP);
    }

    pmeshes->hyperplane_color[i_part] = parts[i_part]->hyperplane_color;
    pmeshes->thread_color    [i_part] = parts[i_part]->thread_color;

    _part_free(parts[i_part]);
  }
  free(parts);

  /*
   * All entities are reorder - In case of HO mesh we need to append all ho vtx in vtx_ln_to_gn AND pvtx_coord
   */
  int have_ho = 0;
  if(dmesh_nodal != NULL) {
    have_ho =  PDM_dmesh_nodal_have_ho(dmesh_nodal);
  }

  if(have_ho == 1) {

    /* Deduce vtx from the connectivity by elmt */
    int          *pn_vtx_all       = NULL;
    PDM_g_num_t **vtx_all_ln_to_gn = NULL;
    PDM_generate_ho_vtx_ln_to_gn(dmesh_nodal,
                                 n_part,
                                 pn_cell,
                                 pcell_ln_to_gn,
                                 pn_face,
                                 pface_ln_to_gn,
                                 pn_edge,
                                 pedge_ln_to_gn,
                                 pn_vtx,
                                 pvtx_ln_to_gn,
                                 &pn_vtx_all,
                                 &vtx_all_ln_to_gn);
    if(0 == 1) {
      for(int i_part = 0; i_part < n_part; ++i_part) {
        log_trace("pn_vtx_all[%i] = %i \n", i_part, pn_vtx_all[i_part]);
      }
    }

    for (int i_part = 0; i_part < n_part; i_part++) {

      PDM_g_num_t *tmp_vtx_ln_to_gn = NULL;
      double      *tmp_vtx_coord    = NULL;
      PDM_part_mesh_entity_ln_to_gn_get(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VERTEX,
                                        &tmp_vtx_ln_to_gn,
                                        PDM_OWNERSHIP_USER);
      PDM_part_mesh_vtx_coord_get(pmeshes->pmesh,
                                  i_part,
                                  &tmp_vtx_coord,
                                  PDM_OWNERSHIP_USER);
      free(tmp_vtx_ln_to_gn);
      free(tmp_vtx_coord);

      pvtx_ln_to_gn[i_part] = vtx_all_ln_to_gn[i_part];
    }
    // free(pvtx_ln_to_gn);
    free(pvtx_coord);
    PDM_part_dcoordinates_to_pcoordinates(comm,
                                          n_part,
                                          vtx_distrib,
                                          dvtx_coord,
                                          pn_vtx_all,
                   (const PDM_g_num_t **) vtx_all_ln_to_gn,
                                         &pvtx_coord);


    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_n_entity_set(pmeshes->pmesh, i_part, PDM_MESH_ENTITY_VERTEX, pn_vtx_all[i_part]);
      PDM_part_mesh_entity_ln_to_gn_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_MESH_ENTITY_VERTEX,
                                        pvtx_ln_to_gn[i_part],
                                        PDM_OWNERSHIP_KEEP);
      PDM_part_mesh_vtx_coord_set(pmeshes->pmesh,
                                  i_part,
                                  pvtx_coord   [i_part],
                                  PDM_OWNERSHIP_KEEP);
    }

    free(vtx_all_ln_to_gn);
    free(pn_vtx_all);
  }



  /*
   *  All data has been reorder, we can now and only now setup desired comm graph
   */
  if(pn_cell != NULL){
    int         **pinternal_face_bound_proc_idx = NULL;
    int         **pinternal_face_bound_part_idx = NULL;
    int         **pinternal_face_bound          = NULL;
    PDM_part_generate_entity_graph_comm(comm,
                                        distrib_partition,
                                        face_distrib,
                                        n_part,
                                        pn_face,
                 (const PDM_g_num_t **) pface_ln_to_gn,
                                        NULL,
                                       &pinternal_face_bound_proc_idx,
                                       &pinternal_face_bound_part_idx,
                                       &pinternal_face_bound,
                                        NULL);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_BOUND_TYPE_FACE,
                                        pinternal_face_bound_proc_idx[i_part],
                                        pinternal_face_bound_part_idx[i_part],
                                        pinternal_face_bound[i_part],
                                        PDM_OWNERSHIP_KEEP);

    }
    free(pinternal_face_bound_proc_idx );
    free(pinternal_face_bound_part_idx );
    free(pinternal_face_bound          );

  } else if(pn_edge  != NULL) {
    int **pinternal_edge_bound_proc_idx = NULL;
    int **pinternal_edge_bound_part_idx = NULL;
    int **pinternal_edge_bound          = NULL;
    PDM_part_generate_entity_graph_comm(comm,
                                        distrib_partition,
                                        edge_distrib,
                                        n_part,
                                        pn_edge,
                 (const PDM_g_num_t **) pedge_ln_to_gn,
                                        NULL,
                                       &pinternal_edge_bound_proc_idx,
                                       &pinternal_edge_bound_part_idx,
                                       &pinternal_edge_bound,
                                        NULL);

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                        i_part,
                                        PDM_BOUND_TYPE_EDGE,
                                        pinternal_edge_bound_proc_idx[i_part],
                                        pinternal_edge_bound_part_idx[i_part],
                                        pinternal_edge_bound[i_part],
                                        PDM_OWNERSHIP_KEEP);
    }
    free(pinternal_edge_bound_proc_idx );
    free(pinternal_edge_bound_part_idx );
    free(pinternal_edge_bound          );
  }

  PDM_part_generate_entity_graph_comm(comm,
                                      distrib_partition,
                                      vtx_distrib,
                                      n_part,
                                      pn_vtx,
               (const PDM_g_num_t **) pvtx_ln_to_gn,
                                      NULL,
                                     &pinternal_vtx_bound_proc_idx,
                                     &pinternal_vtx_bound_part_idx,
                                     &pinternal_vtx_bound,
                                      NULL);

  // Finally complete parts structure with internal join data and bounds
  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_part_mesh_part_graph_comm_set(pmeshes->pmesh,
                                      i_part,
                                      PDM_BOUND_TYPE_VTX,
                                      pinternal_vtx_bound_proc_idx[i_part],
                                      pinternal_vtx_bound_part_idx[i_part],
                                      pinternal_vtx_bound[i_part],
                                      PDM_OWNERSHIP_KEEP);

  }
  free(pinternal_vtx_bound_proc_idx);
  free(pinternal_vtx_bound_part_idx);
  free(pinternal_vtx_bound);

  if(own_edge_distrib == 1) {
    free(edge_distrib);
  }
  if(own_face_distrib == 1) {
    free(face_distrib);
  }

  if(pn_cell != NULL) {
    free(pn_cell);
    free(pcell_ln_to_gn);
    free(pcell_face_idx);
    free(pcell_face);
  }
  if(pn_face != NULL) {
    free(pn_face);
    free(pface_ln_to_gn);
    free(pface_edge_idx);
    free(pface_edge);
    if(pface_vtx_idx  !=  NULL) {
      free(pface_vtx_idx);
      free(pface_vtx);
    }
  }
  if(pn_edge != NULL) {
    free(pn_edge);
    free(pedge_ln_to_gn);
    free(pedge_vtx_idx);
    free(pedge_vtx);
  }
  if(pn_vtx != NULL) {
    free(pn_vtx);
    free(pvtx_ln_to_gn);
    free(pvtx_coord);
  }
  free(distrib_node);
  free(vtx_distrib);
  free(distrib_partition);
}

/**
 * \brief Construct the partitioned meshes on every zones
 *
 * \param [in]   multipart             Pointer to \ref PDM_multipart_t object
 */
void
PDM_multipart_run_ppart
(
 PDM_multipart_t *multipart
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(multipart->comm, &i_rank);
  PDM_MPI_Comm_size(multipart->comm, &n_rank);


  /*
   * Step 1 : Split the graph (If we have a dmesh_nodal and prepare the dmesh before the treatment for faces and elements are the same)
   * Step 2 : Rebuild all connectivity in coherent manner
   * Step 3 : Apply all reordering
   * Step 4 : Deduce all mesh_nodal connectivity
   */

  if (multipart->merge_blocks)
  {
    // 1. Generate global numerotation using all blocks
    // 2. Call the partitionner once on the global numbering
  } else {
    PDM_timer_t *timer = PDM_timer_create();
    double cum_elapsed_time = 0;
    int *starting_part_idx =  PDM_array_new_idx_from_sizes_int(multipart->n_part, multipart->n_zone);

    // int is_by_elt = 0;
    for (int i_zone = 0; i_zone < multipart->n_zone; ++i_zone) {
      PDM_dmesh_nodal_t* dmesh_nodal = multipart->dmeshes_nodal[i_zone];
      if (dmesh_nodal != NULL) { // element representation
        // is_by_elt = 1;
        // PDM_printf("Partitionning elt zone %d/%d \n", i_zone+1, multipart->n_zone);
        PDM_MPI_Comm comm = multipart->comm;
        PDM_split_dual_t split_method = multipart->split_method;
        int n_part = multipart->n_part[i_zone];
        _part_mesh_t* pmesh = &(multipart->pmeshes[i_zone]);

        // _run_ppart_zone_nodal(dmesh_nodal,pmesh,split_method,n_part,comm);

        const double* part_fraction      = &multipart->part_fraction[starting_part_idx[i_zone]];
        PDM_part_size_t part_size_method = multipart->part_size_method;

        //Convert dmesh nodal to dmesh
        PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(1, comm, PDM_OWNERSHIP_KEEP);
        PDM_dmesh_nodal_generate_distribution(dmesh_nodal);
        PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, 0, dmesh_nodal);

        // printf("dmesh_nodal->n_cell_abs = "PDM_FMT_G_NUM" \n", dmesh_nodal->n_cell_abs );
        if(dmesh_nodal->n_cell_abs != 0) {
          PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);
        } else {
          PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_EDGE,
                                           PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_EDGE);
        }

        PDM_dmesh_t  *_dmesh = NULL;
        PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &_dmesh);
        _run_ppart_zone2(_dmesh, dmesh_nodal, pmesh, n_part, split_method, part_size_method, part_fraction, comm);
        multipart->dmeshes  [i_zone] = _dmesh;
        multipart->dmn_to_dm[i_zone] = dmn_to_dm; /* Store it - We need it for PDM_multipart_get_part_mesh_nodal */
        // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);

      } else { // face representation
        // PDM_printf("Partitionning face zone %d/%d \n", i_zone+1, multipart->n_zone);

        PDM_MPI_Comm comm = multipart->comm;

        PDM_split_dual_t split_method    = multipart->split_method;
        PDM_part_size_t part_size_method = multipart->part_size_method;

        const double* part_fraction = &multipart->part_fraction[starting_part_idx[i_zone]];

        PDM_dmesh_t  *_dmeshes =   multipart->dmeshes[i_zone];
        _part_mesh_t *_pmeshes = &(multipart->pmeshes[i_zone]);

        int n_part = multipart->n_part[i_zone];


        if (0 && i_rank == 0)
          PDM_printf("Running partitioning for block %i...\n", i_zone+1);
        PDM_timer_resume(timer);
        // _run_ppart_zone(_dmeshes, _pmeshes, n_part, split_method, part_size_method, part_fraction, comm);
        _run_ppart_zone2(_dmeshes, NULL, _pmeshes, n_part, split_method, part_size_method, part_fraction, comm);
        PDM_timer_hang_on(timer);
        if (0 && i_rank == 0)
          PDM_printf("...completed (elapsed time : %f)\n", PDM_timer_elapsed(timer) - cum_elapsed_time);
        cum_elapsed_time = PDM_timer_elapsed(timer);
      }
    }
    PDM_timer_free(timer);

    free(starting_part_idx);

    // Now rebuild joins over the zones
    // On commente temporairement
    // if (!is_by_elt) { // WARNING TODO also implement if element representation
    //   _search_matching_joins(multipart);
    // }



  }
}

/**
 * \brief Retreive the partitionned nodal mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [out] pmesh_nodal           Nodal partitionned mesh
 * \param [in]  ownership             Who is responsible to free retreived data ?
 *
 */

void
PDM_multipart_get_part_mesh_nodal
(
PDM_multipart_t        *multipart,
const int               i_zone,
PDM_part_mesh_nodal_t **pmesh_nodal,
PDM_ownership_t         ownership
)
{
  assert(i_zone < multipart->n_zone);

  _part_mesh_t      *pmesh       = &(multipart->pmeshes    [i_zone]);
  PDM_dmesh_nodal_t *dmesh_nodal = multipart->dmeshes_nodal[i_zone];
  if (dmesh_nodal == NULL) {
    *pmesh_nodal = NULL;
  }
  else {
    int n_part = multipart->n_part[i_zone];
    if(dmesh_nodal->mesh_dimension == 3){
      *pmesh_nodal = _compute_part_mesh_nodal_3d(dmesh_nodal, pmesh, n_part, ownership);
    } else if(dmesh_nodal->mesh_dimension == 2){
      *pmesh_nodal = _compute_part_mesh_nodal_2d(dmesh_nodal, pmesh, n_part, ownership);
    } else {
      PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_compute_part_mesh_nodal error : Bad dmesh_nodal dimension \n");
    }
  }
}

/**
 * \brief Retreive the partitionned mesh
 *
 * \param [in]  multipart             Pointer to \ref PDM_multipart_t object
 * \param [in]  i_zone                Id of zone
 * \param [out] pmesh                 Partitionned mesh
 *
 */


void
PDM_multipart_get_part_mesh
(
       PDM_multipart_t  *multipart,
 const int               i_zone,
       PDM_part_mesh_t **pmesh
)
{
  _pdm_multipart_t *_multipart = (_pdm_multipart_t *) multipart;
  assert(i_zone < _multipart->n_zone);

  *pmesh = &(_multipart->pmeshes    [i_zone]);
}

/**
 *
 * \brief Returns the dimensions of a given partition
 *
 * \param [in]
 *
 */
void
PDM_multipart_part_dim_get
(
PDM_multipart_t *multipart,
const int        i_zone,
const int        i_part,
      int       *n_cell,
      int       *n_face,
      int       *n_face_part_bound,
      int       *n_vtx,
      int       *n_proc,
      int       *n_total_part,
      int       *s_cell_face,
      int       *s_face_vtx,
      int       *s_face_bound,
      int       *n_bound_groups
)
{

  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  *n_cell = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
  *n_face = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE  );
  *n_vtx  = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VERTEX);

  PDM_MPI_Comm_size(multipart->comm, n_proc);
  *n_total_part = _pmeshes.tn_part;

  *s_cell_face = 1;

  int *cell_face     = NULL;
  int *cell_face_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 &cell_face,
                                 &cell_face_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  int *face_vtx     = NULL;
  int *face_vtx_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 &face_vtx,
                                 &face_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  if(*n_cell > 0) {
    *s_cell_face = cell_face_idx[*n_cell];
  }
  if(face_vtx_idx != NULL) {
    *s_face_vtx  = face_vtx_idx[*n_face];
  } else {
    *s_face_vtx  = 0;
  }

  int                     *face_part_bound_proc_idx;
  int                     *face_part_bound_part_idx;
  int                     *face_part_bound;
  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_BOUND_TYPE_FACE,
                                    &face_part_bound_proc_idx,
                                    &face_part_bound_part_idx,
                                    &face_part_bound,
                                    PDM_OWNERSHIP_BAD_VALUE);


  *n_face_part_bound = 0;
  if(face_part_bound_part_idx != NULL) {
    *n_face_part_bound = face_part_bound_part_idx[*n_total_part];
  }
  *n_bound_groups = PDM_part_mesh_n_bound_get(_pmeshes.pmesh, PDM_BOUND_TYPE_FACE);

  int         *face_bound_idx;
  int         *face_bound;
  PDM_g_num_t *face_bound_ln_to_gn;
  PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_BOUND_TYPE_FACE,
                                 &face_bound_idx,
                                 &face_bound,
                                 &face_bound_ln_to_gn,
                                 PDM_OWNERSHIP_BAD_VALUE);

  *s_face_bound = 0;
  if(face_bound_idx !=NULL) {
    *s_face_bound   = face_bound_idx[*n_bound_groups];
  }
}


void
PDM_multipart_part_graph_comm_get
(
 PDM_multipart_t    *multipart,
 const int           i_zone,
 const int           i_part,
 PDM_bound_type_t    bound_type,
 int               **ppart_bound_proc_idx,
 int               **ppart_bound_part_idx,
 int               **ppart_bound,
 PDM_ownership_t     ownership
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    bound_type,
                                    ppart_bound_proc_idx,
                                    ppart_bound_part_idx,
                                    ppart_bound,
                                    ownership);
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
void
PDM_multipart_part_val_get
(
PDM_multipart_t     *multipart,
const int            i_zone,
const int            i_part,
      int         ***elt_vtx_idx,
      int         ***elt_vtx,
      PDM_g_num_t ***elt_section_ln_to_gn,
      int          **cell_tag,
      int          **cell_face_idx,
      int          **cell_face,
      PDM_g_num_t  **cell_ln_to_gn,
      int          **face_tag,
      int          **face_cell,
      int          **face_vtx_idx,
      int          **face_vtx,
      PDM_g_num_t  **face_ln_to_gn,
      int          **face_part_bound_proc_idx,
      int          **face_part_bound_part_idx,
      int          **face_part_bound,
      int          **vtx_tag,
      double       **vtx,
      PDM_g_num_t  **vtx_ln_to_gn,
      int          **face_bound_idx,
      int          **face_bound,
      PDM_g_num_t  **face_bound_ln_to_gn,
      int          **face_join_idx,
      int          **face_join,
      PDM_g_num_t  **face_join_ln_to_gn
)
{

  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  *cell_tag = NULL;
  *face_tag = NULL;
  *vtx_tag  = NULL;

  // *cell_ln_to_gn = _pmeshes.parts[i_part]->cell_ln_to_gn;
  // *face_ln_to_gn = _pmeshes.parts[i_part]->face_ln_to_gn;
  // *vtx_ln_to_gn  = _pmeshes.parts[i_part]->vtx_ln_to_gn;


  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_CELL,
                                    cell_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_MESH_ENTITY_VERTEX,
                                    vtx_ln_to_gn,
                                    PDM_OWNERSHIP_BAD_VALUE);

  // *cell_face_idx = _pmeshes.parts[i_part]->cell_face_idx;
  // *cell_face     = _pmeshes.parts[i_part]->cell_face;
  // *face_cell     = _pmeshes.parts[i_part]->face_cell;
  // *face_vtx_idx  = _pmeshes.parts[i_part]->face_vtx_idx;
  // *face_vtx      = _pmeshes.parts[i_part]->face_vtx;

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 cell_face,
                                 cell_face_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 face_vtx,
                                 face_vtx_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);

  int *face_cell_idx = NULL;
  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 PDM_CONNECTIVITY_TYPE_FACE_CELL,
                                 face_cell,
                                 &face_cell_idx,
                                 PDM_OWNERSHIP_BAD_VALUE);
  assert(face_cell_idx == NULL);

  PDM_part_mesh_vtx_coord_get(_pmeshes.pmesh, i_part, vtx, PDM_OWNERSHIP_BAD_VALUE);

  // *face_part_bound_proc_idx = _pmeshes.parts[i_part]->face_part_bound_proc_idx;
  // *face_part_bound_part_idx = _pmeshes.parts[i_part]->face_part_bound_part_idx;
  // *face_part_bound          = _pmeshes.parts[i_part]->face_part_bound;

  PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                    i_part,
                                    PDM_BOUND_TYPE_FACE,
                                    face_part_bound_proc_idx,
                                    face_part_bound_part_idx,
                                    face_part_bound,
                                    PDM_OWNERSHIP_BAD_VALUE);

  int n_cell = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
  if (n_cell > 0) {
    // *face_bound_idx       = _pmeshes.parts[i_part]->face_bound_idx;
    // *face_bound           = _pmeshes.parts[i_part]->face_bound;
    // *face_bound_ln_to_gn  = _pmeshes.parts[i_part]->face_bound_ln_to_gn;
    PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                   i_part,
                                   PDM_BOUND_TYPE_FACE,
                                   face_bound_idx,
                                   face_bound,
                                   face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_BAD_VALUE);
  } else {
    // *face_bound_idx       = _pmeshes.parts[i_part]->edge_bound_idx;
    // *face_bound           = _pmeshes.parts[i_part]->edge_bound;
    // *face_bound_ln_to_gn  = _pmeshes.parts[i_part]->edge_bound_ln_to_gn;
    PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                   i_part,
                                   PDM_BOUND_TYPE_EDGE,
                                   face_bound_idx,
                                   face_bound,
                                   face_bound_ln_to_gn,
                                   PDM_OWNERSHIP_BAD_VALUE);
  }
  *face_join_idx        = NULL; // _pmeshes.parts[i_part]->face_join_idx;
  *face_join            = NULL; // _pmeshes.parts[i_part]->face_join;
  *face_join_ln_to_gn   = NULL; // _pmeshes.parts[i_part]->face_join_ln_to_gn;

  *elt_vtx_idx          = NULL; // _pmeshes.parts[i_part]->elt_vtx_idx;
  *elt_vtx              = NULL; // _pmeshes.parts[i_part]->elt_vtx;
  *elt_section_ln_to_gn = NULL; // _pmeshes.parts[i_part]->elt_section_ln_to_gn;
}


int
PDM_multipart_part_tn_part_get
(
PDM_multipart_t *multipart,
const int        i_zone
)
{
  assert(i_zone < multipart->n_zone);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];
  return PDM_part_mesh_tn_part_get(_pmeshes.pmesh);
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
int
PDM_multipart_part_connectivity_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      PDM_connectivity_type_t   connectivity_type,
      int                     **connect,
      int                     **connect_idx,
      PDM_ownership_t           ownership
)
{
  PDM_UNUSED(ownership);
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];
  int pn_entity = -1;

  if( connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_ELMT ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_CELL ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_FACE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_EDGE ||
      connectivity_type == PDM_CONNECTIVITY_TYPE_CELL_VTX)
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_FACE_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_EDGE_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_EDGE  );
  } else if( connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_ELMT ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_CELL ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_FACE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_EDGE ||
             connectivity_type == PDM_CONNECTIVITY_TYPE_VTX_VTX )
  {
    pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VERTEX);
  } else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_part_connectivity_get error : Wrong connectivity_type \n");
  }

  PDM_part_mesh_connectivity_get(_pmeshes.pmesh,
                                 i_part,
                                 connectivity_type,
                                 connect,
                                 connect_idx,
                                 ownership);

  return pn_entity;
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
int
PDM_multipart_part_n_entity_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  int pn_entity = 0;
  switch (entity_type) {
    case PDM_MESH_ENTITY_CELL:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );
      break;
    case PDM_MESH_ENTITY_FACE:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_FACE  );
      break;
    case PDM_MESH_ENTITY_EDGE:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_EDGE  );
      break;
    case PDM_MESH_ENTITY_VERTEX:
       pn_entity = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VERTEX);
      break;
    default:
      PDM_error(__FILE__, __LINE__, 0, "PDM_multipart_part_ln_to_gn_get error : Wrong entity_type \n");
      break;
  }

  return pn_entity;
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
int
PDM_multipart_part_ln_to_gn_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      PDM_g_num_t         **entity_ln_to_gn,
      PDM_ownership_t       ownership
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  int pn_entity = PDM_multipart_part_n_entity_get(multipart, i_zone, i_part, entity_type);

  PDM_part_mesh_entity_ln_to_gn_get(_pmeshes.pmesh,
                                    i_part,
                                    entity_type,
                                    entity_ln_to_gn,
                                    ownership);

  return pn_entity;
}

/**
 *
 * \brief Returns the data arrays of a given partition
 */
int
PDM_multipart_partition_color_get
(
PDM_multipart_t            *multipart,
const int                   i_zone,
const int                   i_part,
      PDM_mesh_entities_t   entity_type,
      int                 **entity_color,
      PDM_ownership_t       ownership
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  int pn_entity = PDM_multipart_part_n_entity_get(multipart, i_zone, i_part, entity_type);
  PDM_part_mesh_entity_color_get(_pmeshes.pmesh,
                                 i_part,
                                 entity_type,
                                 entity_color,
                                 ownership);

  return pn_entity;
}

void
PDM_multipart_part_hyperplane_color_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **hyperplane_color,
      PDM_ownership_t   ownership
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_zone]);

  *hyperplane_color = _pmeshes->hyperplane_color[i_part];
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_zone].is_owner_hyperplane_color = PDM_FALSE;
  } else {
    multipart->pmeshes[i_zone].is_owner_hyperplane_color = PDM_TRUE;
  }
}

void
PDM_multipart_part_thread_color_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **thread_color,
      PDM_ownership_t   ownership
)
{
  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_zone]);

  *thread_color = _pmeshes->thread_color[i_part];
    if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_zone].is_owner_thread_color = PDM_FALSE;
  } else {
    multipart->pmeshes[i_zone].is_owner_thread_color = PDM_TRUE;
  }
}

void
PDM_multipart_part_ghost_infomation_get
(
PDM_multipart_t        *multipart,
const int               i_zone,
const int               i_part,
      int             **vtx_ghost_information,
      PDM_ownership_t   ownership
)
{

  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t* _pmeshes = (&multipart->pmeshes[i_zone]);

  *vtx_ghost_information = _pmeshes->vtx_ghost_information[i_part];
  if(ownership == PDM_OWNERSHIP_USER || ownership == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE) {
    multipart->pmeshes[i_zone].is_owner_vtx_ghost_information = PDM_FALSE;
  } else {
    multipart->pmeshes[i_zone].is_owner_vtx_ghost_information = PDM_TRUE;
  }
}


/**
 *
 * \brief Return times for a given zone
 * (NOT IMPLEMENTED)
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 * \param [in]   i_zone         Id of current zone
 * \param [out]  elapsed        Elapsed time
 * \param [out]  cpu            CPU time
 * \param [out]  cpu_user       User CPU time
 * \param [out]  cpu_sys        System CPU time
 *
 */

void
PDM_multipart_time_get
(
 PDM_multipart_t *multipart,
 const int        i_zone,
 double         **elapsed,
 double         **cpu,
 double         **cpu_user,
 double         **cpu_sys
)
{
  assert(i_zone < multipart->n_zone);

  // PDM_printf("PDM_multipart_time_get: Not implemented\n");
  *elapsed  = NULL;
  *cpu      = NULL;
  *cpu_user = NULL;
  *cpu_sys  = NULL;

}

/**
 *
 * \brief Free the structure
 *
 * \param [in]   multipart      Pointer to \ref PDM_multipart_t object
 */

void
PDM_multipart_free
(
 PDM_multipart_t *multipart
)
{
  // free(multipart->dmeshes_ids);
  for (int i_zone = 0; i_zone < multipart->n_zone; i_zone++) {
    if (multipart->pmeshes[i_zone].joins_ids != NULL) {
      free(multipart->pmeshes[i_zone].joins_ids);
    }

    for (int i_part = 0; i_part < multipart->n_part[i_zone]; i_part++) {
      if(multipart->pmeshes[i_zone].vtx_ghost_information[i_part] != NULL) {
        if(multipart->pmeshes[i_zone].is_owner_vtx_ghost_information == PDM_TRUE) {
          free(multipart->pmeshes[i_zone].vtx_ghost_information[i_part]);
        }
      }

      if(multipart->pmeshes[i_zone].hyperplane_color[i_part] != NULL) {
        if(multipart->pmeshes[i_zone].is_owner_hyperplane_color == PDM_TRUE) {
          free(multipart->pmeshes[i_zone].hyperplane_color[i_part]);
        }
      }

      if(multipart->pmeshes[i_zone].thread_color[i_part] != NULL) {
        if(multipart->pmeshes[i_zone].is_owner_thread_color == PDM_TRUE) {
          free(multipart->pmeshes[i_zone].thread_color[i_part]);
        }
      }
    }
    free(multipart->pmeshes[i_zone].vtx_ghost_information);
    free(multipart->pmeshes[i_zone].hyperplane_color);
    free(multipart->pmeshes[i_zone].thread_color);

    // if (multipart->pmeshes[i_zone].parts != NULL) {
    //   for (int ipart = 0; ipart < multipart->n_part[i_zone]; ipart++) {
    //     _part_free(multipart->pmeshes[i_zone].parts[ipart], multipart->owner);
    //   }
    //   free(multipart->pmeshes[i_zone].parts);
    // }
    PDM_part_mesh_free(multipart->pmeshes[i_zone].pmesh);

    if(multipart->dmn_to_dm[i_zone] != NULL) {
      PDM_dmesh_nodal_to_dmesh_free(multipart->dmn_to_dm[i_zone]);
      multipart->dmn_to_dm[i_zone] = NULL;
    }
  }
  free(multipart->pmeshes);
  free(multipart->dmeshes);
  free(multipart->dmeshes_nodal);
  free(multipart->dmn_to_dm);
  free(multipart->n_part);

  //PDM_part_renum_method_purge();
  free (multipart);
  multipart = NULL;

  // PDM_printf("Cleaned from PDM_multipart_free\n");
}


int
PDM_multipart_part_vtx_coord_get
(
PDM_multipart_t                *multipart,
const int                       i_zone,
const int                       i_part,
      double                  **vtx_coord,
      PDM_ownership_t           ownership
)
{
  PDM_UNUSED(ownership);

  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);

  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  // *vtx_coord = _pmeshes.parts[i_part]->vtx;
  PDM_part_mesh_vtx_coord_get(_pmeshes.pmesh,
                              i_part,
                              vtx_coord,
                              ownership);

  return PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_VERTEX);
}


void PDM_multipart_bound_get
(
 PDM_multipart_t   *multipart,
 const int          i_zone,
 const int          i_part,
 PDM_bound_type_t   bound_type,
 int               *n_bound,
 int              **bound_idx,
 int              **bound,
 PDM_g_num_t      **bound_ln_to_gn,
 PDM_ownership_t    ownership
)
{

  assert(i_zone < multipart->n_zone && i_part < multipart->n_part[i_zone]);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  *n_bound = PDM_part_mesh_n_bound_get(_pmeshes.pmesh, bound_type);

  PDM_part_mesh_bound_concat_get(_pmeshes.pmesh,
                                 i_part,
                                 bound_type,
                                 bound_idx,
                                 bound,
                                 bound_ln_to_gn,
                                 ownership);
}

void
PDM_multipart_stat_get
(
 PDM_multipart_t  *multipart,
 int               i_zone,
 int              *cells_average,
 int              *cells_median,
 double           *cells_std_deviation,
 int              *cells_min,
 int              *cells_max,
 int              *bound_part_faces_average,
 int              *bound_part_faces_median,
 double           *bound_part_faces_std_deviation,
 int              *bound_part_faces_min,
 int              *bound_part_faces_max,
 int              *bound_part_faces_sum
)
{
  int n_rank;
  PDM_MPI_Comm_size(multipart->comm, &n_rank);

  assert(i_zone < multipart->n_zone);
  _part_mesh_t _pmeshes = multipart->pmeshes[i_zone];

  int* dpart_proc = (int *) malloc((n_rank + 1) * sizeof(int));
  PDM_MPI_Allgather((void *) &multipart->n_part[i_zone],
                    1,
                    PDM_MPI_INT,
           (void *) (&dpart_proc[1]),
                    1,
                    PDM_MPI_INT,
                    multipart->comm);

  dpart_proc[0] = 0;
  for (int i = 1; i < n_rank+1; i++) {
    dpart_proc[i] = dpart_proc[i] + dpart_proc[i-1];
  }

  int *n_loc = (int *) malloc(multipart->n_part[i_zone]  * sizeof(int));
  int *n_tot = (int *) malloc(dpart_proc[n_rank]          * sizeof(int));

  int *s_loc = (int *) malloc(multipart->n_part[i_zone]  * sizeof(int));
  int *s_tot = (int *) malloc(dpart_proc[n_rank]          * sizeof(int));

  for (int i = 0; i < multipart->n_part[i_zone]; i++) {
    n_loc[i] = 0;
    s_loc[i] = 0;
  }

  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    n_tot[i] = 0;
    s_tot[i] = 0;
  }


  int tn_part = dpart_proc[n_rank];
  for (int i_part = 0; i_part < multipart->n_part[i_zone]; i_part++) {
    n_loc[i_part] = PDM_part_mesh_n_entity_get(_pmeshes.pmesh, i_part, PDM_MESH_ENTITY_CELL  );

    int                     *face_part_bound_proc_idx;
    int                     *face_part_bound_part_idx;
    int                     *face_part_bound;
    PDM_part_mesh_part_graph_comm_get(_pmeshes.pmesh,
                                      i_part,
                                      PDM_BOUND_TYPE_FACE,
                                      &face_part_bound_proc_idx,
                                      &face_part_bound_part_idx,
                                      &face_part_bound,
                                      PDM_OWNERSHIP_BAD_VALUE);
    if(face_part_bound_part_idx != NULL) {
      s_loc[i_part] = face_part_bound_part_idx[tn_part];
    }
  }

  int *n_part_proc = (int *) malloc((n_rank) * sizeof(int));

  for (int i = 0; i < n_rank; i++) {
    n_part_proc[i] = dpart_proc[i+1] - dpart_proc[i];
  }

  PDM_MPI_Allgatherv(n_loc,
                     multipart->n_part[i_zone],
                     PDM_MPI_INT,
                     n_tot,
                     n_part_proc,
                     dpart_proc,
                     PDM_MPI_INT,
                     multipart->comm);

  PDM_MPI_Allgatherv(s_loc,
                     multipart->n_part[i_zone],
                     PDM_MPI_INT,
                     s_tot,
                     n_part_proc,
                     dpart_proc,
                     PDM_MPI_INT,
                     multipart->comm);

  PDM_quick_sort_int(s_tot, 0, dpart_proc[n_rank]-1);
  PDM_quick_sort_int(n_tot, 0, dpart_proc[n_rank]-1);

  double   _cells_average;
  double   _bound_part_faces_average;

  *bound_part_faces_min = -1;
  *bound_part_faces_max = -1;
  *cells_min = -1;
  *cells_max = -1;
  _cells_average = 0;
  _bound_part_faces_average = 0;

  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    if (*bound_part_faces_min < 0)
      *bound_part_faces_min = s_tot[i];
    else
      *bound_part_faces_min = PDM_MIN(*bound_part_faces_min, s_tot[i]);
    if (*bound_part_faces_max < 0)
      *bound_part_faces_max = s_tot[i];
    else
      *bound_part_faces_max = PDM_MAX(*bound_part_faces_max, s_tot[i]);
    if (*cells_min < 0)
      *cells_min = n_tot[i];
    else
      *cells_min = PDM_MIN(*cells_min, n_tot[i]);
    if (*cells_max < 0)
      *cells_max = n_tot[i];
    else
      *cells_max = PDM_MAX(*cells_max, n_tot[i]);

    _cells_average += n_tot[i];
    _bound_part_faces_average += s_tot[i];
  }

  _cells_average = (_cells_average/((double) dpart_proc[n_rank]));
  *bound_part_faces_sum = (int) _bound_part_faces_average;
  _bound_part_faces_average =
    _bound_part_faces_average/((double) dpart_proc[n_rank]);

  *cells_average = (int) round(_cells_average);
  *bound_part_faces_average = (int) round(_bound_part_faces_average);

  *cells_std_deviation = 0.;
  *bound_part_faces_std_deviation = 0.;
  for (int i = 0; i < dpart_proc[n_rank]; i++) {
    *cells_std_deviation += (n_tot[i] - _cells_average) * (n_tot[i] - _cells_average);
    *bound_part_faces_std_deviation += (s_tot[i] - _bound_part_faces_average) *
                                      (s_tot[i] - _bound_part_faces_average);
  }

  *cells_std_deviation = sqrt(*cells_std_deviation/dpart_proc[n_rank]);
  *bound_part_faces_std_deviation =
    sqrt(*bound_part_faces_std_deviation/dpart_proc[n_rank]);

  int mid = dpart_proc[n_rank]/2;
  if (dpart_proc[n_rank] % 2 == 1) {
    *cells_median = n_tot[mid];
    *bound_part_faces_median = s_tot[mid];
  }

  else {
    *cells_median =(int) round((n_tot[mid-1] + n_tot[mid])/2.);
    *bound_part_faces_median = (int) ((s_tot[mid-1] + s_tot[mid])/2.);
  }

  free(n_part_proc);
  free(n_tot);
  free(s_tot);
  free(n_loc);
  free(s_loc);
  free(dpart_proc);
}


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
