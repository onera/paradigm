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
#include "pdm_dmesh_nodal_elmts_priv.h"

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

struct PDM_DMesh_nodal_vtx_t {
  PDM_l_num_t       n_vtx;          /*!< Number of vertices */
  PDM_real_t       *_coords;        /*!< Coordinates
                                       * (Memory mapping) (size = 3 * \ref n_vtx) */
  PDM_g_num_t      *distrib;        /*!< Distribution on the processes
                               * (size = \ref n_rank + 1) */
  PDM_ownership_t   owner;
};

/**
 * \struct  PDM_Mesh_nodal_geom_prepa_sections_t
 *
 * \brief   Used to build sections from cell to face face to edge connectivity
 *
 */
struct _pdm_dmesh_nodal_t {

  PDM_MPI_Comm           comm;                     /*!< MPI Communicator */
  int                    n_rank;                   /*!< Number of processes */
  int                    i_rank;                   /*!< Number of processes */
  int mesh_dimension;                              /*! Principal dimension of meshes */

  PDM_g_num_t            n_cell_abs;               /*!< Global number of elements */
  PDM_g_num_t            n_face_abs;               /*!< Global number of faces    */
  PDM_g_num_t            n_edge_abs;               /*!< Global number of edges    */
  PDM_g_num_t            n_vtx_abs;                /*!< Global number of vertices */

  PDM_DMesh_nodal_vtx_t *vtx;                      /*!< Description des sommmets de chaque partition */

  _pdm_dmesh_nodal_elts_t* volumic;
  _pdm_dmesh_nodal_elts_t* surfacic;
  _pdm_dmesh_nodal_elts_t* ridge;
  _pdm_dmesh_nodal_elts_t* corner;

  // To move in pdm_dmesh
  PDM_l_num_t            dn_cell;                  /*!< Local number of cells in the local block */
  PDM_l_num_t           *dcell_face_idx;           /*!< Index of the cell to face connectivity
                                                    * (size = \ref dn_cell) */
  PDM_g_num_t           *dcell_face;               /*!< Cell to face connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */
  PDM_g_num_t           *_dface_cell;              /*!< Face to cell connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */
  PDM_g_num_t           *cell_distrib;             /*!< Distribution of cells (size = number of processes + 1) */
  PDM_l_num_t            dn_face;                  /*!< Local number of faces in the local block */
  PDM_l_num_t           *_dface_vtx_idx;           /*!< Index of the cell to face connectivity
                                                    * (size = \ref dn_cell) */
  PDM_g_num_t           *_dface_vtx;               /*!< Cell to face connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */
  PDM_g_num_t           *face_distrib;             /*!< Distribution of faces (size = number of processes + 1) */

};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_PRIV_H__ */
