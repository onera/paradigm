/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_generate_mesh.h"

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_multipart.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_sphere_vol_gen.h"
#include "pdm_dcube_nodal_gen.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static PDM_part_mesh_nodal_t*
_dmn_to_pmn
(
  const PDM_MPI_Comm       comm,
  const PDM_split_dual_t   part_method,
  const int                n_part,
        PDM_dmesh_nodal_t *dmn
)
{

 int n_zone = 1;
 PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                               &n_part,
                                               PDM_FALSE,
                                               part_method,
                                               PDM_PART_SIZE_HOMOGENEOUS,
                                               NULL,
                                               comm,
                                               PDM_OWNERSHIP_KEEP);

 PDM_multipart_set_reordering_options(mpart,
                                      -1,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      NULL,
                                      "PDM_PART_RENUM_FACE_NONE");

 PDM_multipart_register_dmesh_nodal(mpart, 0, dmn);

 PDM_multipart_run_ppart(mpart);

 PDM_part_mesh_nodal_t *pmn = NULL;
 PDM_multipart_get_part_mesh_nodal(mpart,
                                   0,
                                   &pmn,
                                   PDM_OWNERSHIP_USER);

  // free
  // PDM_multipart_free(mpart);

 return pmn;

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a partitionned sphere mesh (2D).
 *
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
)
{
  if (elt_type == PDM_MESH_NODAL_TRIA3) {

    assert(order == 1);

    // generate distributed sphere mesh
    PDM_dmesh_nodal_t *dmn = NULL;
    if (n_u != n_v) {

      PDM_sphere_surf_gen_nodal(comm,
                                n_u,
                                n_v,
                                center_x,
                                center_y,
                                center_z,
                                radius,
                                &dmn);

    } else {

      PDM_sphere_surf_icosphere_gen_nodal(comm,
                                          n_u,
                                          center_x,
                                          center_y,
                                          center_z,
                                          radius,
                                          &dmn);

    }

    // generate partionned sphere mesh
    PDM_part_mesh_nodal_t *pmn = _dmn_to_pmn(comm,
                                             part_method,
                                             n_part,
                                             dmn);

    // frees
    PDM_DMesh_nodal_free(dmn);

    return pmn;

  } else {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
  }

}

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
)
{
  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_sphere(comm,
                                                        PDM_MESH_NODAL_TRIA3,
                                                        1,
                                                        NULL,
                                                        1.,
                                                        0.,
                                                        0.,
                                                        0.,
                                                        100,
                                                        100,
                                                        1,
                                                        PDM_SPLIT_DUAL_WITH_HILBERT);

  *n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                         0);

  *coords = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                              0);

  PDM_g_num_t *pelt_ln_to_gn       = NULL;
  int         *parent_num          = NULL;
  PDM_g_num_t *parent_entity_g_num = NULL;
  PDM_part_mesh_nodal_section_std_get(pmn,
                                      0,
                                      0,
                                      elt_vtx,
                                      &pelt_ln_to_gn,
                                      &parent_num,
                                      &parent_entity_g_num);

  *n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                 0,
                                                 0);

  *elt_vtx_idx = malloc(sizeof(int) * ((*n_elt)+1));
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }

}

/**
 *
 * \brief Create a partitionned ball mesh (3D).
 *
 * \param [in]  comm        MPI communicator
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
)
{
  PDM_dmesh_nodal_t *dmn = NULL;

  // ball without a hole
  if (hole_radius == 0) {

    // generate distributed ball mesh
    if (n_x != n_y && n_x != n_z) {

       PDM_sphere_vol_gen_nodal(comm,
                                n_x,
                                n_y,
                                n_z,
                                radius,
                                center_x,
                                center_y,
                                center_z,
                                elt_type,
                                order,
                                &dmn);
    } else {

      if (elt_type == PDM_MESH_NODAL_TETRA4) {

        assert(order == 1);

        PDM_sphere_vol_icosphere_gen_nodal(comm,
                                           n_x,
                                           center_x,
                                           center_y,
                                           center_z,
                                           radius,
                                           &dmn);

      } else {
        PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
      }

    }

  // ball with a hole
  } else {

    if (elt_type == PDM_MESH_NODAL_PRISM6) {

      assert(order == 1);
      assert(n_x == n_y && n_x == n_z);

      PDM_sphere_vol_hollow_gen_nodal(comm,
                                      n_x,
                                      n_layer,
                                      center_x,
                                      center_y,
                                      center_z,
                                      hole_radius,
                                      radius,
                                      geometric_ratio,
                                      &dmn);

    } else {
      PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
    }
  }

  // generate partionned ball mesh
  PDM_part_mesh_nodal_t *pmn = _dmn_to_pmn(comm,
                                           part_method,
                                           n_part,
                                           dmn);

  // frees
  PDM_DMesh_nodal_free(dmn);

  return pmn;

}

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
)
{
  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_ball(comm,
                                                      PDM_MESH_NODAL_TETRA4,
                                                      1,
                                                      NULL,
                                                      1.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      0.,
                                                      100,
                                                      100,
                                                      100,
                                                      0,
                                                      0.,
                                                      1,
                                                      PDM_SPLIT_DUAL_WITH_HILBERT);

  *n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                         0);

  *coords = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                              0);

  PDM_g_num_t *pelt_ln_to_gn       = NULL;
  int         *parent_num          = NULL;
  PDM_g_num_t *parent_entity_g_num = NULL;
  PDM_part_mesh_nodal_section_std_get(pmn,
                                      0,
                                      0,
                                      elt_vtx,
                                      &pelt_ln_to_gn,
                                      &parent_num,
                                      &parent_entity_g_num);

  *n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                 0,
                                                 0);

  *elt_vtx_idx = malloc(sizeof(int) * ((*n_elt)+1));
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 4; // because PDM_MESH_NODAL_TETRA4
  }
}

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
 const PDM_MPI_Comm      comm,
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
)
{
  // generate distributed rectangle mesh
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create (comm,
                                                         n_x,
                                                         n_y,
                                                         0.,
                                                         lengthx,
                                                         xmin,
                                                         ymin,
                                                         zmin,
                                                         elt_type,
                                                         order,
                                                         PDM_OWNERSHIP_USER);
  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // scale to rectangle is necessary
  if (lengthx != lengthy) {

    int dn_vtx = PDM_DMesh_nodal_n_vtx_get(dmn);

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn);

    double y = 0.;
    for (int i = 0; i < dn_vtx; i++) {
      y = (dvtx_coord[3*i+1] - ymin) / lengthx;
      dvtx_coord[3*i+1] = y * lengthy + ymin; // change y coordinate
    }

  }

  // generate partionned ball mesh
  PDM_part_mesh_nodal_t *pmn = _dmn_to_pmn(comm,
                                           part_method,
                                           n_part,
                                           dmn);

  // frees
  PDM_DMesh_nodal_free(dmn);

  return pmn;
}

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
)
{
  PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_rectangle(comm,
                                                           PDM_MESH_NODAL_TRIA3,
                                                           1,
                                                           NULL,
                                                           0.,
                                                           0.,
                                                           0.,
                                                           10.,
                                                           5.,
                                                           100,
                                                           100,
                                                           1,
                                                           PDM_SPLIT_DUAL_WITH_HILBERT);

  *n_vtx = PDM_part_mesh_nodal_n_vtx_get(pmn,
                                         0);

  *coords = PDM_part_mesh_nodal_vtx_coord_get(pmn,
                                              0);

  PDM_g_num_t *pelt_ln_to_gn       = NULL;
  int         *parent_num          = NULL;
  PDM_g_num_t *parent_entity_g_num = NULL;
  PDM_part_mesh_nodal_section_std_get(pmn,
                                      0,
                                      0,
                                      elt_vtx,
                                      &pelt_ln_to_gn,
                                      &parent_num,
                                      &parent_entity_g_num);

  *n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                 0,
                                                 0);

  *elt_vtx_idx = malloc(sizeof(int) * ((*n_elt)+1));
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
