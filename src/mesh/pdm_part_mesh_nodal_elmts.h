#ifndef __PDM_PART_MESH_NODAL_ELMTS_H__
#define __PDM_PART_MESH_NODAL_ELMTS_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Macro definition
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct _pdm_part_mesh_nodal_elmts_t PDM_part_mesh_nodal_elmts_t;

/*=============================================================================
 * Global variables
 *============================================================================*/

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */
PDM_part_mesh_nodal_elmts_t*
PDM_part_mesh_nodal_elmts_create
(
const int          mesh_dimension,
const int          n_part,
const PDM_MPI_Comm comm
);

int
PDM_part_mesh_nodal_elmts_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt
);

int
PDM_part_mesh_nodal_elmts_ho_add
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const PDM_Mesh_nodal_elt_t         t_elt,
const int                          order,
const char                        *ho_ordering
);

void
PDM_part_mesh_nodal_elmts_free
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

void
PDM_part_mesh_nodal_elmts_std_set
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part,
const int                          n_elt,
const int                         *connec,
const PDM_g_num_t                 *numabs,
const int                         *parent_num,
const PDM_g_num_t                 *parent_entity_g_num,
      PDM_ownership_t              owner
);

void
PDM_part_mesh_nodal_elmts_block_std_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num
);

void
PDM_part_mesh_nodal_elmts_block_std_ho_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      int                         **connec,
      PDM_g_num_t                 **numabs,
      int                         **parent_num,
      PDM_g_num_t                 **parent_entity_g_num,
      int                          *order,
const char                        **ho_ordering
);

/**
 * \brief Return a polygon block description
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly2d_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_block,
 const int                           id_part,
       int                         **connec_idx,
       int                         **connec
);

/**
 * \brief Get the cell-vertex connectivity of a polyhedra block
 *
 * \param [in]  pmne           Pointer to \ref PDM_part_mesh_nodal_elmts object
 * \param [in]  id_block       Block identifier
 * \param [in]  id_part        Partition identifier
 * \param [out] cell_vtx_idx   Index of cell vertex connectivity
 * \param [out] cell_vtx       Cell vertex connectivity
 *
 */

void
PDM_part_mesh_nodal_elmts_block_poly3d_cell_vtx_connect_get
(
       PDM_part_mesh_nodal_elmts_t  *pmne,
 const int                           id_block,
 const int                           id_part,
       int                         **cell_vtx_idx,
       int                         **cell_vtx
);


void
PDM_part_mesh_nodal_elmts_block_poly3d_get
(
      PDM_part_mesh_nodal_elmts_t  *pmne,
const int                           id_block,
const int                           id_part,
      PDM_l_num_t                  *n_face,
      PDM_l_num_t                 **face_vtx_idx,
      PDM_l_num_t                 **face_vtx,
      PDM_l_num_t                 **cell_face_idx,
      PDM_l_num_t                 **cell_face
);

int
PDM_part_mesh_nodal_elmts_block_n_elt_get
(
      PDM_part_mesh_nodal_elmts_t *pmne,
const int                          id_block,
const int                          id_part
);

int
PDM_part_mesh_nodal_elmts_n_section_get
(
 PDM_part_mesh_nodal_elmts_t *pmne
);

int *
PDM_part_mesh_nodal_elmts_sections_id_get
(
  PDM_part_mesh_nodal_elmts_t *pmne
);

PDM_Mesh_nodal_elt_t
PDM_part_mesh_nodal_elmts_block_type_get
(
      PDM_part_mesh_nodal_elmts_t *mesh,
const int                          id_block
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_H__ */
