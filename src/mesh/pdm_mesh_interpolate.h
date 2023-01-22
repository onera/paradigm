#ifndef __PDM_MESH_INTERPOLATE_H__
#define __PDM_MESH_INTERPOLATE_H__


/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_part_domain_interface.h"
#include "pdm_part_mesh_nodal.h"

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

typedef struct _pdm_mesh_interpolate_t PDM_mesh_interpolate_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

PDM_mesh_interpolate_t*
PDM_mesh_interpolate_create
(
 const int            n_domain,
 const int           *n_part,
 const int           *n_group,
 const int            interp_kind,
 const PDM_MPI_Comm   comm
);


/**
 *
 * \brief Compute a part extension structure
 *
 * \param [in]   part_ext          PDM_part_extension_t
 *
 */
void
PDM_mesh_interpolate_compute
(
  PDM_mesh_interpolate_t *part_ext
);



void
PDM_mesh_interpolate_part_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  int                       n_cell,
  int                       n_face,
  int                       n_edge,
  int                       n_vtx,
  int                      *cell_face_idx,
  int                      *cell_face,
  int                      *face_edge_idx,
  int                      *face_edge,
  int                      *edge_vtx,
  int                      *face_vtx_idx,
  int                      *face_vtx,
  PDM_g_num_t              *cell_ln_to_gn,
  PDM_g_num_t              *face_ln_to_gn,
  PDM_g_num_t              *edge_ln_to_gn,
  PDM_g_num_t              *vtx_ln_to_gn,
  double                   *vtx_coord
);


void
PDM_mesh_interpolate_graph_comm_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  PDM_mesh_entities_t       mesh_entity,
  int                      *entity_part_bound_proc_idx,
  int                      *entity_part_bound_part_idx,
  int                      *entity_part_bound
);

void
PDM_mesh_interpolate_part_group_set
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_domain,
  int                       i_part,
  int                       i_group,
  PDM_bound_type_t          bound_type,
  int                       n_group_entity,
  int                      *group_entity
);

void
PDM_mesh_interpolate_part_domain_interface_shared_set
(
  PDM_mesh_interpolate_t      *mi,
  PDM_part_domain_interface_t *pdi
);

void
PDM_mesh_interpolate_set_field
(
  PDM_mesh_interpolate_t   *mi,
  int                       i_field
);

void
PDM_mesh_interpolate_part_mesh_nodal_set
(
  PDM_mesh_interpolate_t     *mi,
  PDM_part_mesh_nodal_t      *pmn
);


void
PDM_mesh_interpolate_exch
(
  PDM_mesh_interpolate_t      *mi,
  const int                    stride,
        double              ***local_field,
        double             ****bound_field,
        double             ****result_field
);




void
PDM_mesh_interpolate_free
(
 PDM_mesh_interpolate_t *mi
);



/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_INTERPOLATE_H__ */

