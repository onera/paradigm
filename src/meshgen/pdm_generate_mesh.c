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
#include "pdm_part_connectivity_transform.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static void
_dmn_to_multipart
(
  const PDM_MPI_Comm            comm,
  const PDM_split_dual_t        part_method,
  const int                     n_part,
        PDM_dmesh_nodal_t      *dmn,
        PDM_multipart_t       **mpart
)
{

 int n_zone = 1;
 *mpart = PDM_multipart_create(n_zone,
                               &n_part,
                               PDM_FALSE,
                               part_method,
                               PDM_PART_SIZE_HOMOGENEOUS,
                               NULL,
                               comm,
                               PDM_OWNERSHIP_KEEP);

 PDM_multipart_set_reordering_options(*mpart,
                                      -1,
                                      "PDM_PART_RENUM_CELL_NONE",
                                      NULL,
                                      "PDM_PART_RENUM_FACE_NONE");

 PDM_multipart_register_dmesh_nodal(*mpart, 0, dmn);

 PDM_multipart_run_ppart(*mpart);

}

// sphere mesh (2D)

static void
_generate_mesh_sphere
(
 const PDM_MPI_Comm            comm,
 const PDM_Mesh_nodal_elt_t    elt_type,
 const int                     order,
 const char                   *ho_ordering,
 const double                  radius,
 const double                  center_x,
 const double                  center_y,
 const double                  center_z,
 const PDM_g_num_t             n_u,
 const PDM_g_num_t             n_v,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
       PDM_dmesh_nodal_t     **dmn,
       PDM_multipart_t       **mpart
)
{
  PDM_UNUSED(ho_ordering);

  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
  assert(dim == 2);

  if (elt_type == PDM_MESH_NODAL_TRIA3) {

    assert(order == 1);

    // generate distributed sphere mesh
    if (n_u != n_v) {

      PDM_sphere_surf_gen_nodal(comm,
                                n_u,
                                n_v,
                                center_x,
                                center_y,
                                center_z,
                                radius,
                                dmn);

    } else {

      PDM_sphere_surf_icosphere_gen_nodal(comm,
                                          n_u,
                                          center_x,
                                          center_y,
                                          center_z,
                                          radius,
                                          dmn);

    }

    // generate partionned sphere mesh
    _dmn_to_multipart(comm,
                      part_method,
                      n_part,
                      *dmn,
                      mpart);


  } else {
    PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
  }
}

// ball mesh (3D)

static void
_generate_mesh_ball
(
 const PDM_MPI_Comm        comm,
 PDM_Mesh_nodal_elt_t      elt_type,
 int                       order,
 const char               *ho_ordering,
 const double              radius,
 const double              hole_radius,
 const double              center_x,
 const double              center_y,
 const double              center_z,
 const PDM_g_num_t         n_x,
 const PDM_g_num_t         n_y,
 const PDM_g_num_t         n_z,
 const PDM_g_num_t         n_layer,
 const double              geometric_ratio,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
       PDM_dmesh_nodal_t **dmn,
       PDM_multipart_t   **mpart
)
{
  PDM_UNUSED(ho_ordering);

  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
  assert(dim == 3);

  // ball without a hole
  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (hole_radius == 0) {
  PDM_GCC_SUPPRESS_WARNING_POP

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
                                dmn);
    } else {

      if (elt_type == PDM_MESH_NODAL_TETRA4) {

        assert(order == 1);

        PDM_sphere_vol_icosphere_gen_nodal(comm,
                                           n_x,
                                           center_x,
                                           center_y,
                                           center_z,
                                           radius,
                                           dmn);

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
                                      dmn);

    } else {
      PDM_error(__FILE__, __LINE__, 0, "Not implemented yet for element type %d\n", (int) elt_type);
    }
  }

  // generate partionned ball mesh
  _dmn_to_multipart(comm,
                    part_method,
                    n_part,
                    *dmn,
                    mpart);

}

// rectangle mesh (2D)
static void
_generate_mesh_rectangle
(
 const PDM_MPI_Comm      comm,
 PDM_Mesh_nodal_elt_t    elt_type,
 int                     order,
 const char             *ho_ordering,
 double                  xmin,
 double                  ymin,
 double                  zmin,
 double                  lengthx,
 double                  lengthy,
 PDM_g_num_t             n_x,
 PDM_g_num_t             n_y,
 const int               n_part,
 const PDM_split_dual_t  part_method,
       PDM_dmesh_nodal_t **dmn,
       PDM_multipart_t   **mpart
)
{

  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
  assert(dim == 2);

  // generate distributed rectangle mesh
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
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

  if (order > 1) {
    PDM_dcube_nodal_gen_ordering_set(dcube,
                                     ho_ordering);
  }

  PDM_dcube_nodal_gen_build (dcube);

  *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(*dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // scale to rectangle is necessary
  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (lengthx != lengthy) {
  PDM_GCC_SUPPRESS_WARNING_POP

    int dn_vtx = PDM_DMesh_nodal_n_vtx_get(*dmn);

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn);

    double y = 0.;
    for (int i = 0; i < dn_vtx; i++) {
      y = (dvtx_coord[3*i+1] - ymin) / lengthx;
      dvtx_coord[3*i+1] = y * lengthy + ymin; // change y coordinate
    }

  }

  // generate partionned rectangle mesh
  _dmn_to_multipart(comm,
                    part_method,
                    n_part,
                    *dmn,
                    mpart);
}

// parallelepiped mesh (3D)
static void
_generate_mesh_parallelepiped
(
 const PDM_MPI_Comm        comm,
 PDM_Mesh_nodal_elt_t      elt_type,
 int                       order,
 const char               *ho_ordering,
 double                    xmin,
 double                    ymin,
 double                    zmin,
 double                    lengthx,
 double                    lengthy,
 double                    lengthz,
 PDM_g_num_t               n_x,
 PDM_g_num_t               n_y,
 PDM_g_num_t               n_z,
 const int                 n_part,
 const PDM_split_dual_t    part_method,
       PDM_dmesh_nodal_t **dmn,
       PDM_multipart_t   **mpart
)
{

  int dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
  assert(dim == 3);

  // generate distributed parallelepiped mesh
  PDM_dcube_nodal_t *dcube = PDM_dcube_nodal_gen_create(comm,
                                                        n_x,
                                                        n_y,
                                                        n_z,
                                                        lengthx,
                                                        xmin,
                                                        ymin,
                                                        zmin,
                                                        elt_type,
                                                        order,
                                                        PDM_OWNERSHIP_USER);

  if (order > 1) {
    PDM_dcube_nodal_gen_ordering_set(dcube,
                                     ho_ordering);
  }

  PDM_dcube_nodal_gen_build (dcube);

  *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(*dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // scale to parallelepiped is necessary
  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (lengthx != lengthy) {
  PDM_GCC_SUPPRESS_WARNING_POP

    int dn_vtx = PDM_DMesh_nodal_n_vtx_get(*dmn);

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn);

    double y = 0.;
    for (int i = 0; i < dn_vtx; i++) {
      y = (dvtx_coord[3*i+1] - ymin) / lengthx;
      dvtx_coord[3*i+1] = y * lengthy + ymin; // change y coordinate
    }

  }

  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (lengthx != lengthz) {
  PDM_GCC_SUPPRESS_WARNING_POP

    int dn_vtx = PDM_DMesh_nodal_n_vtx_get(*dmn);

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn);

    double z = 0.;
    for (int i = 0; i < dn_vtx; i++) {
      z = (dvtx_coord[3*i+2] - zmin) / lengthx;
      dvtx_coord[3*i+2] = z * lengthz + zmin; // change z coordinate
    }

  }

  // generate partionned parallelepiped mesh
  _dmn_to_multipart(comm,
                    part_method,
                    n_part,
                    *dmn,
                    mpart);

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
 * \param [in]  ho_ordering High order nodes ordering type
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
 const char                  *ho_ordering,
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

  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_sphere(comm,
                        elt_type,
                        order,
                        ho_ordering,
                        radius,
                        center_x,
                        center_y,
                        center_z,
                        n_u,
                        n_v,
                        n_part,
                        part_method,
                        &dmn,
                        &mpart);

  // get partionned sphere mesh
  PDM_part_mesh_nodal_t *pmn   = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;

}

/**
 *
 * \brief Create a simple partitionned sphere mesh (2D).
 *
 * \param [in]   comm        MPI communicator
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
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
)
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_sphere(comm,
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
                        PDM_SPLIT_DUAL_WITH_HILBERT,
                        &dmn,
                        &mpart);

  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);

  // get elt-vtx connectivity
  int  *face_edge     = NULL;
  int  *face_edge_idx = NULL;
  *n_elt = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               0,
                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                               &face_edge,
                                               &face_edge_idx,
                                               PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_compute_face_vtx_from_face_and_edge(*n_elt,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          elt_vtx);

  *elt_vtx_idx = malloc(sizeof(int) * ((*n_elt)+1));
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

}

/**
 *
 * \brief Create a partitionned ball mesh (3D).
 *
 * \param [in]  comm            MPI communicator
 * \param [in]  elt_type        Mesh element type
 * \param [in]  order           Mesh element order
 * \param [in]  ho_ordering     High order nodes ordering type
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
 const char             *ho_ordering,
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
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_ball(comm,
                      elt_type,
                      order,
                      ho_ordering,
                      radius,
                      hole_radius,
                      center_x,
                      center_y,
                      center_z,
                      n_x,
                      n_y,
                      n_z,
                      n_layer,
                      geometric_ratio,
                      n_part,
                      part_method,
                      &dmn,
                      &mpart);

  // get partionned ball mesh
  PDM_part_mesh_nodal_t *pmn   = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;

}

/**
 *
 * \brief Create a simple partitionned ball mesh (3D).
 *
 * \param [in]   comm        MPI communicator
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
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
)
{

  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_ball(comm,
                      PDM_MESH_NODAL_TETRA4,
                      1,
                      NULL,
                      1.,
                      0.,
                      0.,
                      0.,
                      0.,
                      10,
                      10,
                      10,
                      0,
                      0.,
                      1,
                      PDM_SPLIT_DUAL_WITH_HILBERT,
                      &dmn,
                      &mpart);

  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);

  // get elt-vtx connectivity
  int  *cell_face     = NULL;
  int  *cell_face_idx = NULL;
  *n_elt = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               0,
                                               PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                               &cell_face,
                                               &cell_face_idx,
                                               PDM_OWNERSHIP_KEEP);

  int  *face_edge     = NULL;
  int  *face_edge_idx = NULL;
  int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                   0,
                                                   0,
                                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                   &face_edge,
                                                   &face_edge_idx,
                                                   PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

  int  *face_vtx  = NULL;
  PDM_compute_face_vtx_from_face_and_edge(n_face,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          &face_vtx);

  int *face_vtx_idx = malloc(sizeof(int) * (n_face + 1));
  face_vtx_idx[0] = 0;
  for (int i = 0; i < n_face; i++) {
    face_vtx_idx[i+1] = face_vtx_idx[i] + 3; // face of a thetrahedron is a tirangle
  }

  PDM_combine_connectivity(*n_elt,
                           cell_face_idx,
                           cell_face,
                           face_vtx_idx,
                           face_vtx,
                           elt_vtx_idx,
                           elt_vtx);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
  free(face_vtx_idx);
  free(face_vtx);

}

/**
 *
 * \brief Create a partitionned rectangle mesh (2D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering High order nodes ordering type
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
 const char             *ho_ordering,
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
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_rectangle(comm,
                           elt_type,
                           order,
                           ho_ordering,
                           xmin,
                           ymin,
                           zmin,
                           lengthx,
                           lengthy,
                           n_x,
                           n_y,
                           n_part,
                           part_method,
                           &dmn,
                           &mpart);

  // get partionned rectangle mesh
  PDM_part_mesh_nodal_t *pmn   = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;
}

/**
 *
 * \brief Create a simple partitionned rectangle mesh (2D).
 *
 * \param [in]   comm        MPI communicator
 * \param [in]   n_vtx_seg   Number of vertices along each side of the rectangle
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
 const PDM_g_num_t    n_vtx_seg,
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
)
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_rectangle(comm,
                           PDM_MESH_NODAL_TRIA3,
                           1,
                           NULL,
                           0.,
                           0.,
                           0.,
                           10.,
                           5.,
                           n_vtx_seg,
                           n_vtx_seg,
                           1,
                           PDM_SPLIT_DUAL_WITH_HILBERT,
                           &dmn,
                           &mpart);

  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);

  // get elt-vtx connectivity
  int  *face_edge     = NULL;
  int  *face_edge_idx = NULL;
  *n_elt = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               0,
                                               PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                               &face_edge,
                                               &face_edge_idx,
                                               PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_compute_face_vtx_from_face_and_edge(*n_elt,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          elt_vtx);

  *elt_vtx_idx = malloc(sizeof(int) * ((*n_elt)+1));
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

}

/**
 *
 * \brief Create a partitionned parallelepiped mesh (3D).
 *
 * \param [in]  comm        MPI communicator
 * \param [in]  elt_type    Mesh element type
 * \param [in]  order       Mesh element order
 * \param [in]  ho_ordering High order nodes ordering type
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
 const char             *ho_ordering,
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
)
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_parallelepiped(comm,
                                elt_type,
                                order,
                                ho_ordering,
                                xmin,
                                ymin,
                                zmin,
                                lengthx,
                                lengthy,
                                lengthz,
                                n_x,
                                n_y,
                                n_z,
                                n_part,
                                part_method,
                                &dmn,
                                &mpart);

  // get partionned rectangle mesh
  PDM_part_mesh_nodal_t *pmn   = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;
}

/**
 *
 * \brief Create a simple partitionned parallelepiped mesh (3D).
 *
 * \param [in]   comm        MPI communicator
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
 int                 *n_vtx,
 int                 *n_elt,
 double             **coords,
 int                **elt_vtx_idx,
 int                **elt_vtx
)
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;
  _generate_mesh_parallelepiped(comm,
                                PDM_MESH_NODAL_TETRA4,
                                1,
                                NULL,
                                0.,
                                0.,
                                0.,
                                10.,
                                10.,
                                10.,
                                100,
                                100,
                                100,
                                1,
                                PDM_SPLIT_DUAL_WITH_HILBERT,
                                &dmn,
                                &mpart);

  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);

  // get elt-vtx connectivity
  int  *cell_face     = NULL;
  int  *cell_face_idx = NULL;
  *n_elt = PDM_multipart_part_connectivity_get(mpart,
                                               0,
                                               0,
                                               PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                               &cell_face,
                                               &cell_face_idx,
                                               PDM_OWNERSHIP_KEEP);

  int  *face_edge     = NULL;
  int  *face_edge_idx = NULL;
  int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                   0,
                                                   0,
                                                   PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                   &face_edge,
                                                   &face_edge_idx,
                                                   PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx,
                                      &edge_vtx_idx,
                                      PDM_OWNERSHIP_KEEP);

  int  *face_vtx  = NULL;
  PDM_compute_face_vtx_from_face_and_edge(n_face,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          &face_vtx);

  int *face_vtx_idx = malloc(sizeof(int) * (n_face + 1));
  face_vtx_idx[0] = 0;
  for (int i = 0; i < n_face; i++) {
    face_vtx_idx[i+1] = face_vtx_idx[i] + 3; // face of a thetrahedron is a tirangle
  }

  // /!\ Ca marche pas vraiment ça
  PDM_combine_connectivity(*n_elt,
                           cell_face_idx,
                           cell_face,
                           face_vtx_idx,
                           face_vtx,
                           elt_vtx_idx,
                           elt_vtx);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
  free(face_vtx_idx);
  free(face_vtx);
}



void
PDM_generate_mesh_rectangle_ngon
(
 const PDM_MPI_Comm            comm,
 const PDM_Mesh_nodal_elt_t    elt_type,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  lengthx,
 const double                  lengthy,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn
)
{
  PDM_dmesh_nodal_t *dmn   = NULL;
  PDM_multipart_t   *mpart = NULL;
  _generate_mesh_rectangle(comm,
                           elt_type,
                           1,
                           NULL,
                           xmin,
                           ymin,
                           zmin,
                           lengthx,
                           lengthy,
                           n_x,
                           n_y,
                           n_part,
                           part_method,
                           &dmn,
                           &mpart);
  PDM_DMesh_nodal_free(dmn);

  *pn_vtx         = malloc(sizeof(int          ) * n_part);
  *pn_edge        = malloc(sizeof(int          ) * n_part);
  *pn_face        = malloc(sizeof(int          ) * n_part);
  *pvtx_coord     = malloc(sizeof(double      *) * n_part);
  *pedge_vtx      = malloc(sizeof(int         *) * n_part);
  *pface_edge_idx = malloc(sizeof(int         *) * n_part);
  *pface_edge     = malloc(sizeof(int         *) * n_part);
  *pvtx_ln_to_gn  = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pedge_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);
  *pface_ln_to_gn = malloc(sizeof(PDM_g_num_t *) * n_part);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VERTEX,
                                                       &(*pvtx_ln_to_gn)[ipart],
                                                       PDM_OWNERSHIP_USER);

    PDM_multipart_part_vtx_coord_get(mpart,
                                     0,
                                     ipart,
                                     &(*pvtx_coord)[ipart],
                                     PDM_OWNERSHIP_USER);

    (*pn_edge)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_MESH_ENTITY_EDGE,
                                                        &(*pedge_ln_to_gn)[ipart],
                                                        PDM_OWNERSHIP_USER);

    (*pn_face)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_MESH_ENTITY_FACE,
                                                        &(*pface_ln_to_gn)[ipart],
                                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &(*pface_edge)    [ipart],
                                        &(*pface_edge_idx)[ipart],
                                        PDM_OWNERSHIP_USER);

    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &(*pedge_vtx)[ipart],
                                        &edge_vtx_idx,
                                        PDM_OWNERSHIP_USER);
  }

  PDM_multipart_free(mpart);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
