#ifndef __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__
#define __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_mesh_nodal_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \struct  PDM_Mesh_nodal_geom_prepa_blocks_t
 *
 * \brief   Used to build blocks from cell to face face to edge connectivity
 *
 */

typedef struct _pdm_part_mesh_nodal_elmts_t _pdm_part_mesh_nodal_elmts_t;
struct _pdm_part_mesh_nodal_elmts_t {

  PDM_MPI_Comm                         pdm_mpi_comm;              /*!< MPI Communicator            */
  PDM_Mesh_nodal_prepa_blocks_t       *prepa_blocks;              /*!< Blocks preparation          */

  int                                  n_part;
  PDM_l_num_t                         *n_elmts;                   /*!< Nombre de blocs d'elements  */

  int                                  n_blocks;                  /*!< Total number of blocks      */
  int                                 *blocks_id;                 /*!< Blocks identifier           */
  PDM_Mesh_nodal_block_std_t         **blocks_std;                /*!< Standard blocks             */
  PDM_Mesh_nodal_block_poly2d_t      **blocks_poly2d;             /*!< Polygon blocks              */
  PDM_Mesh_nodal_block_poly3d_t      **blocks_poly3d;             /*!< Polyhedron blocks           */

  PDM_l_num_t                        **num_elmt_parent_to_local;  /*!< Initial local numbering to local numbering
                                                                   *   imposed by blocks */

  int                      is_vtx_def_from_parent;                /*<! Are the points defined from parents */
  PDM_g_num_t                          **numabs;                  /*<! Global numbering per elmts per partition */
};

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PART_MESH_NODAL_ELMTS_PRIV_H__ */
