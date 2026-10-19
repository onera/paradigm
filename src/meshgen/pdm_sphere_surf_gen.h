#ifndef __PDM_SPHERE_SURF_GEN_H__
#define __PDM_SPHERE_SURF_GEN_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/
#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a sphere surface mesh (UV sphere)
 *
 * \note The output mesh is *block-distributed*.
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dvtx_coord      Vertex coordinates
 * \param[out] dface_vtx_idx   Index for face -> vertex connectivity
 * \param[out] dface_vtx       Face -> vertex connectivity (global ids)
 * \param[out] distrib_vtx     Distribution of vertices (size : *n_rank* + 1)
 * \param[out] distrib_face    Distribution of faces (size : *n_rank* + 1)
 *
 */

void
PDM_sphere_surf_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
);


/**
 *
 * \brief Create a sphere surface mesh (UV sphere)
 *
 * \note The output mesh is *block-distributed*,
 * in the form of a \ref PDM_dmesh_nodal_t object
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  nu              Number of points in longitude
 * \param[in]  nv              Number of points in latitude
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dmn             Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         nu,
 const PDM_g_num_t         nv,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **dmn
);


/**
 *
 * \brief Create a sphere surface mesh (icosphere)
 *
 * \note The output mesh is *block-distributed*.
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dvtx_coord      Vertex coordinates
 * \param[out] dface_vtx_idx   Index for face -> vertex connectivity
 * \param[out] dface_vtx       Face -> vertex connectivity (global ids)
 * \param[out] distrib_vtx     Vertex distribution (size : *n_rank* + 1)
 * \param[out] distrib_face    Face distribution (size : *n_rank* + 1)
 *
 */

void
PDM_sphere_surf_icosphere_gen
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       double            **dvtx_coord,
       int               **dface_vtx_idx,
       PDM_g_num_t       **dface_vtx,
       PDM_g_num_t       **distrib_vtx,
       PDM_g_num_t       **distrib_face
);


/**
 *
 * \brief Create a sphere surface mesh (icosphere)
 *
 * \note The output mesh is *block-distributed*,
 * in the form of a \ref PDM_dmesh_nodal_t object
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[out] dmn             Sphere mesh in the form of a distributed nodal mesh
 *
 */

void
PDM_sphere_surf_icosphere_gen_nodal
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
       PDM_dmesh_nodal_t **dmn
);


/**
 *
 * \brief Create a sphere surface mesh (icosphere)
 *
 * \note The output mesh is *partitioned*.
 *
 * \param[in]  comm            MPI communicator
 * \param[in]  n               Icosphere subdivision level
 * \param[in]  x_center        x-coordinate of the sphere center
 * \param[in]  y_center        y-coordinate of the sphere center
 * \param[in]  z_center        z-coordinate of the sphere center
 * \param[in]  radius          Sphere radius
 * \param[in]  n_part          Number of partitions
 * \param[in]  part_method     Partitioning method
 * \param[out] pn_vtx          Number of vertices (size : \p n_part)
 * \param[out] pvtx_coord      Vertex coordinates (for each part, size : 3 * \p pn_vtx)
 * \param[out] pvtx_ln_to_gn   Vertex global ids (for each part, size : \p pn_vtx)
 * \param[out] pn_face         Number of faces (size : \p n_part)
 * \param[out] pface_vtx_idx   Index for face->vertex connectivity (for each part, size : \p pn_face + 1)
 * \param[out] pface_vtx       Face->vertex connectivity (for each part, size : \p pface_vtx_idx[\p pn_face])
 * \param[out] pface_ln_to_gn  Face global ids (for each part, size : \p pn_face)
 *
 */

void
PDM_sphere_surf_icosphere_gen_part
(
 const PDM_MPI_Comm        comm,
 const PDM_g_num_t         n,
 const double              x_center,
 const double              y_center,
 const double              z_center,
 const double              radius,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
       int               **pn_vtx,
       double           ***pvtx_coord,
       PDM_g_num_t      ***pvtx_ln_to_gn,
       int               **pn_face,
       int              ***pface_vtx_idx,
       int              ***pface_vtx,
       PDM_g_num_t      ***pface_ln_to_gn
);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* __PDM_SPHERE_SURF_GEN_H__ */
