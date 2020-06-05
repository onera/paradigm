#ifndef __PDM_PARTITIONING_ALGORITHM_H__
#define __PDM_PARTITIONING_ALGORITHM_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

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

/**
 *  \brief Gather the entities splitted by the partitioner
 *   (usually cells) to their attributed partition, using the array mapping
 *   entities id to their assigned partition number.
*/
int
PDM_part_assemble_partitions
(
 const PDM_MPI_Comm    comm,
 const PDM_g_num_t    *part_distribution,
 const PDM_g_num_t    *entity_distribution,
 const int            *dentity_to_part,
       int           **pn_entity,
       PDM_g_num_t  ***pentity_ln_to_gn
);

/**
 *  \brief Construct the face->cell connectivity from the cell->face connectivity
 */
void
PDM_part_reverse_pcellface
(
  const int         n_part,
  const int        *n_cell,
  const int        *n_face,
  const int       **pcell_face_idx,
  const int       **pcell_face,
        int      ***pface_cell
);

/**
 *  \brief Setup cell_ln_to_gn
 */
void
PDM_generate_part_face_group_ln_to_gn
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *face_distribution,
 int                  *dface_group_idx,
 PDM_g_num_t          *dface_group,
 int                   n_part,
 int                   n_face_group,
 int                  *n_faces,
 PDM_g_num_t         **pface_ln_to_gn,
 PDM_g_num_t        ***pface_group_ln_to_gn,
 int                ***pface_group,
 int                ***pface_group_idx
);

/**
 *  \brief Setup cell_ln_to_gn
 */
void
PDM_generate_part_entity_ln_to_gn_sort
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *cell_distribution,
 int                  *dcell_face_idx,
 PDM_g_num_t          *dcell_face,
 int                   n_part,
 int                  *n_elmts,
 PDM_g_num_t         **pcell_ln_to_gn,
 int                 **n_faces,
 PDM_g_num_t        ***pface_ln_to_gn,
 int                ***pcell_face_idx,
 int                ***pcell_face
);

/**
 *  \brief Setup cell_ln_to_gn
 */
void
PDM_generate_part_entity_ln_to_gn_hash
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *cell_distribution,
 int                  *dcell_face_idx,
 PDM_g_num_t          *dcell_face,
 int                   n_part,
 int                  *n_elmts,
 PDM_g_num_t         **pcell_ln_to_gn,
 int                 **n_faces,
 PDM_g_num_t        ***pface_ln_to_gn,
 int                ***pcell_face_idx,
 int                ***pcell_face
);


/**
 *  \brief
 */
void
PDM_generate_entity_graph_comm
(
 const PDM_MPI_Comm    comm,
 PDM_g_num_t          *part_distribution,
 PDM_g_num_t          *entity_distribution,
 int                   n_part,
 int                  *n_entities,
 PDM_g_num_t         **pentity_ln_to_gn,
 int                ***pproc_bound_idx,
 int                ***ppart_bound_idx,
 int                ***pentity_bound_idx
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_PARTITIONING_ALGORITHM_H__ */
