/*
 * \file
 */

#ifndef __PDM_PART_EXTENSION_ALGORITHM_H__
#define __PDM_PART_EXTENSION_ALGORITHM_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_to_part.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

void
PDM_part_extension_entity1_to_entity2
(
  PDM_g_num_t                   shift_by_domain_entity2,
  int                           n_part,
  int                          *pn_entity1,
  PDM_g_num_t                 **pentity1_ln_to_gn,
  int                         **pentity1_to_pentity1_idx,
  int                         **pentity1_to_pentity1_triplet,
  int                         **pentity1_to_pentity1_interface,
  int                          *pn_entity2,
  PDM_g_num_t                 **pentity2_ln_to_gn,
  int                         **pentity2_alrdy_sent,
  int                         **pentity2_entity1_idx,
  int                         **pentity2_entity1,
  int                           prev_dentity2_itrf_n_blk,
  PDM_g_num_t                  *prev_dentity2_itrf_blk_gnum,
  PDM_g_num_t                  *prev_dentity2_itrf_blk_ancstr,
  int                          *prev_dentity2_itrf_blk_path_itrf_strd,
  int                          *prev_dentity2_itrf_blk_path_itrf,
  int                          *prev_dentity2_itrf_gnum_and_itrf_strid,
  PDM_g_num_t                  *prev_dentity2_itrf_gnum_and_itrf_data,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_alrdy_sent_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  int                          *next_dentity2_itrf_n_blk_out,
  PDM_g_num_t                 **next_dentity2_itrf_blk_gnum_out,
  PDM_g_num_t                 **next_dentity2_itrf_blk_ancstr_out,
  int                         **next_dentity2_itrf_blk_path_itrf_strd_out,
  int                         **next_dentity2_itrf_blk_path_itrf_out,
  int                         **next_dentity2_itrf_gnum_and_itrf_strid_out,
  PDM_g_num_t                 **next_dentity2_itrf_gnum_and_itrf_data_out,
  PDM_MPI_Comm                  comm
);

void
PDM_part_extension_interface_by_entity1_to_interface_by_entity2
(
  PDM_part_domain_interface_t  *pdi,
  PDM_bound_type_t              entity1_bound,
  int                           n_domain,
  PDM_g_num_t                  *shift_by_domain_entity2,
  int                          *n_part,
  int                         **pn_entity1_in,
  PDM_g_num_t                ***pentity1_ln_to_gn_in,
  int                        ***pentity1_hint_in,
  int                         **pn_entity2_in,
  PDM_g_num_t                ***pentity2_ln_to_gn_in,
  int                        ***pentity2_entity1_idx_in,
  int                        ***pentity2_entity1_in,
  int                         **pn_entity2_extented_out,
  PDM_g_num_t                ***pentity2_extented_ln_to_gn_out,
  int                        ***pentity2_extented_to_pentity2_idx_out,
  int                        ***pentity2_extented_to_pentity2_triplet_out,
  int                        ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                  comm
);



void
PDM_part_extension_build_entity1_graph
(
  PDM_part_domain_interface_t   *pdi,
  PDM_bound_type_t               entity1_bound,
  int                            n_domain,
  int                           *n_part,
  int                          **pn_entity1_in,
  PDM_g_num_t                 ***pentity1_ln_to_gn_in,
  int                         ***pentity1_hint_in,
  int                         ***pentity1_extented_to_pentity1_idx_out,
  int                         ***pentity1_extented_to_pentity1_triplet_out,
  int                         ***pentity1_extented_to_pentity1_interface_out,
  PDM_MPI_Comm                   comm
);

void
PDM_part_extension_pentity1_entity2_to_extended_pentity1_entity2
(
 int                    n_part,
 int                    n_interface,
 PDM_g_num_t            shift_by_domain_entity2,
 int                    prev_dentity2_itrf_n_blk,
 PDM_g_num_t           *prev_dentity2_itrf_blk_gnum,
 PDM_g_num_t           *prev_dentity2_itrf_blk_ancstr,
 int                   *prev_dentity2_itrf_blk_path_itrf_strd,
 int                   *prev_dentity2_itrf_blk_path_itrf,
 int                   *prev_dentity2_itrf_gnum_and_itrf_strid,
 PDM_g_num_t           *prev_dentity2_itrf_gnum_and_itrf_data,
 int                   *prev_dentity2_itrf_gnum_and_itrf_sens,
 int                   *pn_entity1,
 PDM_g_num_t          **pentity1_ln_to_gn,
 int                   *pn_entity2,
 PDM_g_num_t          **pentity2_ln_to_gn,
 int                  **pentity1_entity2_idx,
 int                  **pentity1_entity2,
 int                   *pn_entity1_extended,
 PDM_g_num_t          **pentity1_extended_ln_to_gn,
 int                  **pentity1_extended_to_pentity1_idx,
 int                  **pentity1_extended_to_pentity1_triplet,
 int                  **pentity1_extended_to_pentity1_interface,
 int                  **pn_entity2_extended_out,
 PDM_g_num_t         ***pentity2_extended_ln_to_gn_out,
 int                 ***pextended_entity1_entity2_idx_out,
 int                 ***pextended_entity1_entity2_out,
 int                 ***pentity2_extended_to_pentity2_idx_out,
 int                 ***pentity2_extended_to_pentity2_triplet_out,
 int                 ***pentity2_extended_to_pentity2_interface_out,
 int                   *next_dentity2_itrf_n_blk_out,
 PDM_g_num_t          **next_dentity2_itrf_blk_gnum_out,
 PDM_g_num_t          **next_dentity2_itrf_blk_ancstr_out,
 int                  **next_dentity2_itrf_blk_path_itrf_strd_out,
 int                  **next_dentity2_itrf_blk_path_itrf_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_strid_out,
 PDM_g_num_t          **next_dentity2_itrf_gnum_and_itrf_data_out,
 int                  **next_dentity2_itrf_gnum_and_itrf_sens_out,
 PDM_MPI_Comm           comm
);

void
PDM_part_extension_pconnectivity_to_extented_pconnectivity
(
  PDM_part_domain_interface_t    *pdi,
  PDM_bound_type_t                entity2_bound,
  int                             n_domain,
  PDM_g_num_t                    *shift_by_domain_entity2,
  int                            *n_part,
  int                           **pn_entity1_in,
  PDM_g_num_t                  ***pentity1_ln_to_gn_in,
  int                           **pn_entity2_in,
  PDM_g_num_t                  ***pentity2_ln_to_gn_in,
  int                          ***pentity1_entity2_idx_in,
  int                          ***pentity1_entity2_in,
  int                            *pn_entity1_extented,
  PDM_g_num_t                   **pentity1_extented_ln_to_gn,
  int                           **pentity1_extented_to_pentity1_idx,
  int                           **pentity1_extented_to_pentity1_triplet,
  int                           **pentity1_extented_to_pentity1_interface,
  int                           **pn_entity2_extented_out,
  PDM_g_num_t                  ***pentity2_extented_ln_to_gn_out,
  int                          ***pextented_entity1_entity2_idx_out,
  int                          ***pextented_entity1_entity2_out,
  int                          ***pentity2_extented_to_pentity2_idx_out,
  int                          ***pentity2_extented_to_pentity2_triplet_out,
  int                          ***pentity2_extented_to_pentity2_interface_out,
  PDM_MPI_Comm                    comm
);



// Prevoir un concatenate




#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_EXTENSION_ALGORITHM_H__ */
