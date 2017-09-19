#ifndef __PDM_DMESH_NODAL_PRIV_H__
#define __PDM_DMESH_NODAL_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_handles.h"

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
 * \struct PDM_Mesh_nodal_som_t
 * \brief  Vertices of a mesh partition 
 *
 */

typedef struct PDM_DMesh_nodal_vtx_t PDM_DMesh_nodal_vtx_t;

struct PDM_Mesh_nodal_vtx_t {
  int          n_proc;  /*!< Number of proceses */ 
  PDM_l_num_t  n_vtx;   /*!< Number of vertices */
  double      *_coords;  /*!< Coordinates 
                          * (Memory mapping) (size = 3 * \ref n_vtx) */
  PDM_g_num_t *distrib;  /*!< Distribution on the processes 
                          * (size = \ref n_proc + 1) */  
};

/**
 * \struct PDM_Mesh_nodal_block_std_t
 * \brief  Standard geometric block 
 *
 */

typedef struct PDM_DMesh_nodal_block_std_t {

  int                     n_proc;  /*!< Number of proceses */ 
  PDM_Mesh_nodal_elt_t    t_elt;   /*!< Element type */
  PDM_l_num_t             n_elt;   /*!< Number elements */
  PDM_g_num_t            *_connec; /*!< Connectivity (Memory mapping) 
                                    *   (size = Number of vertices per element * \ref n_elt)  */
  PDM_g_num_t            *distrib; /*!< Distribution on the processes (size = \ref n_proc + 1) */  
  
} PDM_DMesh_nodal_block_std_t;


/**
 * \struct PDM_Mesh_nodal_block_poly2d_t
 * \brief  Polygon geometric block 
 *
 */

typedef struct PDM_DMesh_nodal_block_poly2d_t {

  int          n_proc;        /*!< Number of proceses */ 
  PDM_l_num_t  n_elt;         /*!< Number of elements of each partition */
  PDM_l_num_t  *_connec_idx;  /*!< Index of elements connectivity 
                               *  (Memory mapping) (size = \ref n_elt + 1) */
  PDM_g_num_t  *_connec;      /*!< Elements connectivity 
                               * (Memory mapping) (size = \ref connec_idx[\ref n_elt]) */
  PDM_g_num_t  *distrib;      /*!< Distribution on the processes (size = \ref n_proc + 1) */  

} PDM_DMesh_nodal_block_poly2d_t;


/**
 * \struct PDM_Mesh_nodal_block_poly3d_t
 * \brief  Polyhedron geometric block 
 *
 */

typedef struct PDM_DMesh_nodal_block_poly3d_t{

  int          n_proc;       /*!< Number of proceses */ 
  PDM_l_num_t  n_elt;        /*!< Number of elements */
  PDM_l_num_t  n_face;       /*!< Number of faces */
  PDM_l_num_t  *_facvtx_idx; /*!< Index of faces connectivity 
                              * (Memory mapping) (Size = \ref n_face + 1) */

  PDM_g_num_t  *_facvtx;      /*!< Faces connectivity 
                               * (Memory mapping) (Size = \ref _facvtx_idx[\ref n_face])*/             
  PDM_l_num_t  *_cellfac_idx; /*!< Index of cell->face connectivity 
                               * (Memory mapping) (Size = \ref n_cell + 1) */

  PDM_g_num_t  *_cellfac;     /*!< cell->face connectivity 
                               * (Memory mapping) (Size = \ref _cellfac[\ref n_cell]) */
  PDM_g_num_t  *distrib;      /*!< Distribution on the processes (size = \ref n_proc + 1) */  

} PDM_DMesh_nodal_block_poly3d_t;


/**
 * \struct  PDM_Mesh_nodal_geom_prepa_blocks_t
 *
 * \brief   Used to build blocks from cell to face face to edge connectivity
 *
 */

struct _PDM_Mesh_nodal_t {

  PDM_g_num_t                         n_som_abs;                /*!< Global number of vertices */
  PDM_g_num_t                         n_elt_abs;                /*!< Global number of elements */
  PDM_DMesh_nodal_vtx_t               *vtx;                    /*!< Description des sommmets de chaque partition */
  PDM_Handles_t                       *blocks_std;              /*!< Standard blocks */
  PDM_Handles_t                       *blocks_poly2d;           /*!< Polygon blocks */
  PDM_Handles_t                       *blocks_poly3d;           /*!< Polyhedron blocks */
  PDM_MPI_Comm                         pdm_mpi_comm;            /*!< MPI Communicator */
  int                                 *blocks_id;               /*!< Blocks identifier */
  int                                  n_blocks;                /*!< Total number of blocks */
} ;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_PRIV_H__ */
