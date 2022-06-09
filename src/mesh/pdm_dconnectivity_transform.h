#ifndef PDM_DCONNECTIVITY_TRANSFORM_H_
#define PDM_DCONNECTIVITY_TRANSFORM_H_

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Types definition
 *============================================================================*/

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

void
PDM_deduce_combine_connectivity
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       int            **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
);


/**
 *
 * \brief Compute the dual connectivty of entity1
 *
 * \param [in]   comm                  PDM_MPI communicator
 * \param [in]   entity1_distrib       distribution of entity1 over the procs (size=n_rank+1)
 * \param [in]   entity2_distrib       distribution of entity2 over the procs (size=n_rank+1)
 * \param [in]   dentity1_entity2_idx
 * \param [in]   dentity1_entity2
 * \param [in]   dentity1_entity2      is array is signed
 * \param [out]   dentity2_entity1_idx
 * \param [out]   dentity2_entity1
 */
void
PDM_dconnectivity_transpose
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
       PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
       int              is_signed,
       int            **dentity2_entity1_idx,
       PDM_g_num_t    **dentity2_entity1
);

void PDM_dfacecell_to_dcellface
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const PDM_g_num_t* dface_cell,
  int**              dcell_face_idx,
  PDM_g_num_t**      dcell_face,
  PDM_MPI_Comm       comm
);
void PDM_dcellface_to_dfacecell
(
  const PDM_g_num_t* face_distri,
  const PDM_g_num_t* cell_distri,
  const int*         dcell_face_idx,
  const PDM_g_num_t* dcell_face,
  PDM_g_num_t**      dface_cell,
  PDM_MPI_Comm       comm
);

void
PDM_deduce_combine_connectivity_dual
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity1_distrib,
 const PDM_g_num_t     *entity2_distrib,
 const int             *dentity1_entity2_idx,
 const PDM_g_num_t     *dentity1_entity2,
 const int             *dentity2_entity3_idx,
 const PDM_g_num_t     *dentity2_entity3,
 const int              is_signed,
       PDM_g_num_t    **dentity1_entity3_idx,
       PDM_g_num_t    **dentity1_entity3
);


void
PDM_dorder_reverse
(
 const PDM_MPI_Comm     comm,
 const PDM_g_num_t     *entity_distrib,
 const PDM_g_num_t     *dentity1_entity2,
       PDM_g_num_t    **dentity2_entity1
);


void
PDM_dgroup_entity_transpose
(
 int            n_group,
 int           *dgroup_entity_idx,
 PDM_g_num_t   *dgroup_entity,
 PDM_g_num_t   *distrib_entity,
 int          **dentity_group_idx,
 int          **dentity_group,
 PDM_MPI_Comm   comm
);

void
PDM_dentity_group_transpose
(
 int            n_group,
 int           *dentity_group_idx,
 int           *dentity_group,
 PDM_g_num_t   *distrib_entity,
 int          **dgroup_entity_idx,
 PDM_g_num_t  **dgroup_entity,
 PDM_MPI_Comm   comm
);

void
PDM_dconnectivity_to_extract_dconnectivity
(
 const PDM_MPI_Comm    comm,
       int             n_selected_entity1,
       PDM_g_num_t    *select_entity1,
       PDM_g_num_t    *entity1_distribution,
       int            *dentity1_entity2_idx,
       PDM_g_num_t    *dentity1_entity2,
       PDM_g_num_t   **extract_entity1_distribution,
       PDM_g_num_t   **extract_entity2_distribution,
       int           **dextract_entity1_entity2_idx,
       PDM_g_num_t   **dextract_entity1_entity2,
       PDM_g_num_t   **dparent_entity1_g_num,
       PDM_g_num_t   **dparent_entity2_g_num,
       PDM_g_num_t   **entity1_old_to_new
);


#ifdef  __cplusplus
}
#endif

#endif  /* PDM_DCONNECTIVITY_TRANSFORM_H_ */