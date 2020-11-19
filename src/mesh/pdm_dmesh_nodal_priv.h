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

//-->>
typedef enum {

  PDM_SECTION_TYPE_STD3D  = 0,
  PDM_SECTION_TYPE_STD2D  = 1,
  PDM_SECTION_TYPE_STD1D  = 2,
  PDM_SECTION_TYPE_POLY2D = 3,
  PDM_SECTION_TYPE_POLY3D = 4

} PDM_section_type_t;
//<<--

/**
 * \struct PDM_Mesh_nodal_som_t
 * \brief  Vertices of a mesh partition
 *
 */

typedef struct PDM_DMesh_nodal_vtx_t PDM_DMesh_nodal_vtx_t;

struct PDM_DMesh_nodal_vtx_t {
  PDM_l_num_t       n_vtx;          /*!< Number of vertices */
  const PDM_real_t *_coords;        /*!< Coordinates
                                       * (Memory mapping) (size = 3 * \ref n_vtx) */
  PDM_g_num_t      *distrib;        /*!< Distribution on the processes
                                     * (size = \ref n_rank + 1) */
};

/**
 * \struct PDM_Mesh_nodal_section_std_t
 * \brief  Standard geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_std_t {

  PDM_Mesh_nodal_elt_t    t_elt;   /*!< Element type */
  PDM_l_num_t             n_elt;   /*!< Number elements */
  PDM_g_num_t            *_connec; /*!< Connectivity (Memory mapping)
                                    *   (size = Number of vertices per element * \ref n_elt)  */
  PDM_g_num_t            *distrib; /*!< Distribution on the processes (size = \ref n_rank + 1) */

} PDM_DMesh_nodal_section_std_t;


/**
 * \struct PDM_Mesh_nodal_section_poly2d_t
 * \brief  Polygon geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_poly2d_t {

  PDM_l_num_t  n_elt;         /*!< Number of elements of each partition */
  PDM_l_num_t  *_connec_idx;  /*!< Index of elements connectivity
                               *  (Memory mapping) (size = \ref n_elt + 1) */
  PDM_g_num_t  *_connec;      /*!< Elements connectivity
                               * (Memory mapping) (size = \ref connec_idx[\ref n_elt]) */
  PDM_g_num_t  *distrib;      /*!< Distribution on the processes (size = \ref n_rank + 1) */

} PDM_DMesh_nodal_section_poly2d_t;


/**
 * \struct PDM_Mesh_nodal_section_poly3d_t
 * \brief  Polyhedron geometric section
 *
 */

typedef struct PDM_DMesh_nodal_section_poly3d_t{

  PDM_l_num_t  n_elt;          /*!< Number of elements */
  PDM_l_num_t  n_face;         /*!< Number of faces */
  PDM_l_num_t  *_face_vtx_idx; /*!< Index of faces connectivity
                                * (Memory mapping) (Size = \ref n_face + 1) */

  PDM_g_num_t  *_face_vtx;      /*!< Faces connectivity
                                 * (Memory mapping) (Size = \ref _face_vtx_idx[\ref n_face])*/
  PDM_l_num_t  *_cell_face_idx; /*!< Index of cell->face connectivity
                                 * (Memory mapping) (Size = \ref n_cell + 1) */

  PDM_g_num_t  *_cell_face;     /*!< cell->face connectivity
                                 * (Memory mapping) (Size = \ref _cell_face_idx[\ref n_cell]) */
  PDM_g_num_t  *distrib;        /*!< Distribution on the processes (size = \ref n_rank + 1) */

} PDM_DMesh_nodal_section_poly3d_t;


/**
 * \struct  PDM_Mesh_nodal_geom_prepa_sections_t
 *
 * \brief   Used to build sections from cell to face face to edge connectivity
 *
 */
struct _pdm_dmesh_nodal_t {

  int mesh_dimension;                               /*! Principal dimension of meshes */

  PDM_g_num_t            n_cell_abs;               /*!< Global number of elements */
  PDM_g_num_t            n_face_abs;               /*!< Global number of faces    */
  PDM_g_num_t            n_edge_abs;               /*!< Global number of edges    */
  PDM_g_num_t            n_vtx_abs;                /*!< Global number of vertices */

  PDM_g_num_t            n_elmt_tot;
  PDM_g_num_t            n_elmt;
  PDM_g_num_t            n_elmt_l1;
  PDM_g_num_t            n_elmt_l2;

  int          n_group_elmt;
  int         *dgroup_elmt_idx;
  PDM_g_num_t *dgroup_elmt;

  // A gere dans la function chapeau - Lien a faire au dessus - Optionel
  // int         *delmt_join_idx;
  // PDM_g_num_t *delmt_join;


  PDM_DMesh_nodal_vtx_t *vtx;                   /*!< Description des sommmets de chaque partition */

  PDM_section_type_t    *section_type;
  int                   *section_idx;
  int                    n_section_tot;                       /*!< Total number of sections */

  int                    n_section;                           /*!< Total number of sections           */

  int                    n_section_std;                       /*!< Total number of standard sections  */
  int                    n_section_poly3d;                    /*!< Total number of olyhedron sections */

  PDM_DMesh_nodal_section_std_t    **sections_std;             /*!< Standard sections                  */
  PDM_DMesh_nodal_section_poly3d_t **sections_poly3d;          /*!< Polyhedron sections                */
  PDM_g_num_t                       *section_distribution;     /*!< Element distribution               */

  int                                n_section_l1;             /*!< Total number of sections           */
  int                                n_section_std_l1;         /*!< Total number of standard sections  */
  int                                n_section_poly2d_l1;      /*!< Total number of polygon sections   */
  PDM_g_num_t                       *section_distribution_l1;  /*!< Element distribution               */

  PDM_DMesh_nodal_section_std_t    **sections_std_l1;          /*!< Standard sections                  */
  PDM_DMesh_nodal_section_poly2d_t **sections_poly2d_l1;       /*!< Polygon sections                   */

  int                                n_section_l2;             /*!< Total number of sections           */
  int                                n_section_std_l2;         /*!< Total number of standard sections  */
  PDM_g_num_t                       *section_distribution_l2;  /*!< Element distribution               */

  PDM_DMesh_nodal_section_std_t    **sections_std_l2;          /*!< Standard sections                  */

  PDM_MPI_Comm           pdm_mpi_comm;             /*!< MPI Communicator */
  int                    n_rank;                   /*!< Number of processes */
  int                    i_rank;                   /*!< Number of processes */

  // To move in pdm_dmesh
  PDM_l_num_t            dn_elmt;                  /*!< Local number of cells in the local block */
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


  PDM_g_num_t           *edge_distrib;             /*!< Distribution of cells (size = number of processes + 1) */
  PDM_l_num_t            dn_edge;                  /*!< Local number of faces in the local block */
  PDM_l_num_t           *_dedge_vtx_idx;           /*!< Index of the cell to face connectivity
                                                    * (size = \ref dn_cell) */
  PDM_g_num_t           *_dedge_vtx;               /*!< Cell to face connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */
  PDM_l_num_t           *dface_edge_idx;           /*!< Index of the cell to face connectivity
                                                    * (size = \ref dn_cell) */
  PDM_g_num_t           *dface_edge;               /*!< Cell to face connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */
  PDM_g_num_t           *_dedge_face;              /*!< Cell to face connectivity
                                                    * (size = \ref dcell_face_idx[\ref dn_cell] */





};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_NODAL_PRIV_H__ */
