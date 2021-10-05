#ifndef __PDM_DCUBE_NODAL_GEN2_PRIV_H__
#define __PDM_DCUBE_NODAL_GEN2_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _pdm_dcube_t
 * \brief  Distributed cube
 *
 * _dcube_t define a distributed mesh of a cube
 *
 */

struct _pdm_dcube_nodal2_t {
  PDM_MPI_Comm          comm;                   /*!< MPI communicator                          */
  PDM_ownership_t       owner;                  /*!< Which have the responsabilities of results*/
  PDM_ownership_t       owner_for_dmesh_nodal;  /*!< Which have the responsabilities of results*/
  PDM_g_num_t           nx;                     /*!< Number of elements in segments along x    */
  PDM_g_num_t           ny;                     /*!< Number of elements in segments along y    */
  PDM_g_num_t           nz;                     /*!< Number of elements in segments along z    */
  double                length;                 /*!< Segment length                            */
  double                zero_x;                 /*!< Coordinates of the origin                 */
  double                zero_y;                 /*!< Coordinates of the origin                 */
  double                zero_z;                 /*!< Coordinates of the origin                 */

  int                   order;                  /*!< Order of elements                         */

  int                   dn_cell;                /*!< Number of cells stored in this process    */
  int                   dn_face;                /*!< Number of faces stored in this process    */
  int                   dn_edge;                /*!< Number of edges stored in this process    */
  int                   dn_vtx;                 /*!< Number of vertices stored in this process */

  int                   n_face_group;           /*!< Number of faces groups                    */
  //int                   dn_cell;                /*!< Number of cells stored in this process    */
  //int                   dn_vtx;                 /*!< Number of vertices stored in this process */

  PDM_g_num_t          *distrib_hexa;
  PDM_g_num_t          *distrib_quad;
  PDM_g_num_t          *distrib_bar;

  int                   dn_hexa;
  int                   dn_quad;
  int                   dn_bar;



  PDM_g_num_t          *delmt_vtx;              /*!< Faces from vertices connectivity          */
  PDM_g_num_t          *delmt_lim_vtx;          /*!< Faces from vertices connectivity          */
  double               *dvtx_coord;             /*!< Vertices coordinates                      */
  int                  *dface_group_idx;        /*!< Faces groups index                        */
  PDM_g_num_t          *dface_group;            /*!< Faces groups                              */
  PDM_Mesh_nodal_elt_t  t_elt;                  /*!< Type of elements to generate              */

  PDM_dmesh_nodal_t    *dmesh_nodal;            /*!< Results                                   */

} ;


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DCUBE_NODAL_GEN2_PRIV_H__ */
