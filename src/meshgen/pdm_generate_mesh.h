#ifndef __PDM_GENERATE_MESH_H__
#define __PDM_GENERATE_MESH_H__

#include "pdm.h"
#include "pdm_part_mesh.h"
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

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

// TO DO: deformation methods on partitionned data(warning for randomization, same seed)

/**
 *
 * \brief Create a simple partitionned sphere mesh (2D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   order       Mesh element order
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_sphere_simplified
(
 const PDM_MPI_Comm   comm,
 const int            order,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitionned rectangle mesh (2D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   order       Mesh element order
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_rectangle_simplified
(
 const PDM_MPI_Comm   comm,
 const int            order,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitionned ball mesh (3D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   order       Mesh element order
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_ball_simplified
(
 const PDM_MPI_Comm   comm,
 const int            order,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a simple partitionned parallelepiped mesh (3D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   order       Mesh element order
 * \param [out]  n_vtx       Number of vertices
 * \param [out]  n_elt       Number of elements
 * \param [out]  coords      Array of vertex coordinates
 * \param [out]  elt_vtx_idx Index array of the element vertex connectivity
 * \param [out]  elt_vtx     Array of the element vertex connectivity
 *
 */

void
PDM_generate_mesh_parallelepiped_simplified
(
 const PDM_MPI_Comm   comm,
 const int            order,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
);

/**
 *
 * \brief Create a partitionned sphere mesh (2D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering ?
 * \param [in]  radius      Radius of the sphere
 * \param [in]  center_x    x-coordinate of the sphere center
 * \param [in]  center_y    y-coordinate of the sphere center
 * \param [in]  center_z    z-coordinate of the sphere center
 * \param [in]  n_u         Number of points in longitude
 * \param [in]  n_v         Number of points in latitude
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitionning method
 *
 * \return PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_sphere
(
 const PDM_MPI_Comm           comm,
 const PDM_Mesh_nodal_elt_t   elt_type,
 const int                    order,
 const char                 **ho_ordering,
 const double                 radius,
 const double                 center_x,
 const double                 center_y,
 const double                 center_z,
 const PDM_g_num_t            n_u,
 const PDM_g_num_t            n_v,
 const int                    n_part,
 const PDM_split_dual_t       part_method
);

/**
 *
 * \brief Create a partitionned rectangle mesh (2D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering ?
 * \param [in]  xmin        x-coordinate of the rctangle minimum corner
 * \param [in]  ymin        y-coordinate of the rctangle minimum corner
 * \param [in]  zmin        z-coordinate of the rctangle minimum corner
 * \param [in]  lengthx     Length of the rectangle in the x-direction
 * \param [in]  lengthy     Length of the rectangle in the y-direction
 * \param [in]  n_x         Number of points in the x-direction
 * \param [in]  n_y         Number of points in the y-direction
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitionning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_rectangle
(
  const PDM_MPI_Comm     comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char            **ho_ordering,
 double                  xmin,
 double                  ymin,
 double                  zmin,
 double                  lengthx,
 double                  lengthy,
 PDM_g_num_t             n_x,
 PDM_g_num_t             n_y,
 const int               n_part,
 const PDM_split_dual_t  part_method
);

/**
 *
 * \brief Create a partitionned ball mesh (3D).
 *
 * \param [in]  comm            MPI communicator
 * \param [in]  elt_type        Mesh element type
 * \param [in]  order           Mesh element order
 * \param [in]  ho_ordering     ?
 * \param [in]  radius          Radius of the ball
 * \param [in]  hole_radius     Radius of the hole of the ball
 * \param [in]  center_x        x-coordinate of the ball center
 * \param [in]  center_y        y-coordinate of the ball center
 * \param [in]  center_z        z-coordinate of the ball center
 * \param [in]  n_x             Number of vertices on segments in x-direction
 * \param [in]  n_y             Number of vertices on segments in y-direction
 * \param [in]  n_z             Number of vertices on segments in z-direction
 * \param [in]  n_layer         Number of extrusion layers
 * \param [in]  geometric_ratio Geometric ratio for layer thickness
 * \param [in]  n_part          Number of mesh partitions
 * \param [in]  part_method     Mesh partitionning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_ball
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char            **ho_ordering,
 const double            radius,
 const double            hole_radius,
 const double            center_x,
 const double            center_y,
 const double            center_z,
 const PDM_g_num_t       n_x,
 const PDM_g_num_t       n_y,
 const PDM_g_num_t       n_z,
 const PDM_g_num_t       n_layer,
 const double            geometric_ratio,
 const int               n_part,
 const PDM_split_dual_t  part_method
);

/**
 *
 * \brief Create a partitionned parallelepiped mesh (3D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering ?
 * \param [in]  xmin        x-coordinate of the rctangle minimum corner
 * \param [in]  ymin        y-coordinate of the rctangle minimum corner
 * \param [in]  zmin        z-coordinate of the rctangle minimum corner
 * \param [in]  lengthx     Length of the rectangle in the x-direction
 * \param [in]  lengthy     Length of the rectangle in the y-direction
 * \param [in]  lengthz     Length of the rectangle in the z-direction
 * \param [in]  n_x         Number of points in the x-direction
 * \param [in]  n_y         Number of points in the y-direction
 * \param [in]  n_z         Number of points in the z-direction
 * \param [in]  n_part      Number of mesh partitions
 * \param [in]  part_method Mesh partitionning method
 *
 * \return PDM_part_mesh_t or PDM_part_mesh_nodal_t
 *
 */

PDM_part_mesh_nodal_t *
PDM_generate_mesh_parallelepiped
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char            **ho_ordering,
 double                  xmin,
 double                  ymin,
 double                  zmin,
 double                  lengthx,
 double                  lengthy,
 double                  lengthz,
 PDM_g_num_t             n_x,
 PDM_g_num_t             n_y,
 PDM_g_num_t             n_z,
 const int               n_part,
 const PDM_split_dual_t  part_method
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_GENERATE_MESH__ */
