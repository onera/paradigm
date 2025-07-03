/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <assert.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
#include "pdm_generate_mesh.h"
#include "pdm.h"
#include "pdm_array.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_error.h"
#include "pdm_mem_tool.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_connectivity_transform.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_priv.h"
#include "pdm_reader_gamma.h"
#include "pdm_reader_stl.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_sphere_vol_gen.h"
#include "pdm_vtk.h"

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

 int n_domain = 1;
 *mpart = PDM_multipart_create(n_domain,
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

 PDM_multipart_dmesh_nodal_set(*mpart, 0, dmn);

 PDM_multipart_compute(*mpart);

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
 const double            random_factor,
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
                                                        0,
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

  PDM_dcube_nodal_gen_random_factor_set(dcube, random_factor);

  PDM_dcube_nodal_gen_build(dcube);

  *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(*dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // scale to rectangle is necessary
  PDM_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wfloat-equal")
  if (lengthx != lengthy) {
  PDM_GCC_SUPPRESS_WARNING_POP

    int dn_vtx = PDM_DMesh_nodal_n_vtx_get(*dmn);

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn, PDM_OWNERSHIP_BAD_VALUE);

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

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn, PDM_OWNERSHIP_BAD_VALUE);

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

    double* dvtx_coord = PDM_DMesh_nodal_vtx_get(*dmn, PDM_OWNERSHIP_BAD_VALUE);

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


/**
 *
 * \brief  Get file extension from its name
 *
 */
// https://stackoverflow.com/questions/5309471/getting-file-extension-in-c
static const char *
_get_file_extension
(
  const char *filename
)
{
  const char *dot = strrchr(filename, '.');
  if (!dot || dot == filename) {
    return "";
  }
  else {
    return dot + 1;
  }
}


static void
_read_mesh_file
(
  const PDM_MPI_Comm        comm,
  const int                 n_part,
  const PDM_split_dual_t    part_method,
  const char               *filename,
        PDM_dmesh_nodal_t **out_dmn,
        PDM_multipart_t   **out_mpart
)
{
  // Get file extension
  const char *file_extension = _get_file_extension(filename);

  // Use appropriate reader
  if (strcmp(file_extension, "stl") == 0) {
    // STL
    *out_dmn = PDM_reader_stl_dmesh_nodal(comm,
                                          filename);
  }

  else if (strcmp(file_extension, "mesh") == 0) {
    // GAMMA
    *out_dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                            filename,
                                            0,
                                            0);
  }

  else if (strcmp(file_extension, "vtk") == 0) {
    // VTK
    int          n_vtx_field      = 0;
    char       **vtx_field_name   = NULL;
    PDM_data_t  *vtx_field_type   = NULL;
    int         *vtx_field_stride = NULL;
    void       **vtx_field_value  = NULL;
    int          n_elt_field      = 0;
    char       **elt_field_name   = NULL;
    PDM_data_t  *elt_field_type   = NULL;
    int         *elt_field_stride = NULL;
    void       **elt_field_value  = NULL;
    *out_dmn = PDM_vtk_read_to_dmesh_nodal(comm,
                                            filename,
                                            &n_vtx_field,
                                            &vtx_field_name,
                                            &vtx_field_type,
                                            &vtx_field_stride,
                                            &vtx_field_value,
                                            &n_elt_field,
                                            &elt_field_name,
                                            &elt_field_type,
                                            &elt_field_stride,
                                            &elt_field_value);
    for (int i_field = 0; i_field < n_vtx_field; i_field++) {
      PDM_free(vtx_field_name [i_field]);
      PDM_free(vtx_field_value[i_field]);
    }
    for (int i_field = 0; i_field < n_elt_field; i_field++) {
      PDM_free(elt_field_name [i_field]);
      PDM_free(elt_field_value[i_field]);
    }
    PDM_free(vtx_field_name  );
    PDM_free(vtx_field_type  );
    PDM_free(vtx_field_stride);
    PDM_free(vtx_field_value );
    PDM_free(elt_field_name  );
    PDM_free(elt_field_type  );
    PDM_free(elt_field_stride);
    PDM_free(elt_field_value );
  }

  else {
    PDM_error(__FILE__, __LINE__, 0, "PDM_generate_mesh_from_file: Unknown mesh format %s\n", file_extension);
  }

  // Partition mesh
  _dmn_to_multipart(comm,
                    part_method,
                    n_part,
                    *out_dmn,
                    out_mpart);
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/

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
  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;

}


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
                        20,
                        20,
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
                                               &face_edge_idx,
                                               &face_edge,
                                               PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx_idx,
                                      &edge_vtx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_compute_face_vtx_from_face_and_edge(*n_elt,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          elt_vtx);

  PDM_malloc(*elt_vtx_idx, *n_elt + 1, int);
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

}


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

  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_KEEP);

  *n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, 0);

  int         *connec              = NULL;
  PDM_g_num_t *numabs              = NULL;
  int         *parent_num          = NULL;
  PDM_g_num_t *parent_entity_g_num = NULL;
  PDM_part_mesh_nodal_section_std_get(pmn,
                                      0,
                                      0,
                                      &connec,
                                      &numabs,
                                      &parent_num,
                                      &parent_entity_g_num,
                                      PDM_OWNERSHIP_KEEP);

  *elt_vtx_idx = PDM_array_new_idx_from_const_stride_int(4, *n_elt);
  PDM_malloc(*elt_vtx, (*elt_vtx_idx)[*n_elt], int);
  memcpy(*elt_vtx, connec, sizeof(int) * (*elt_vtx_idx)[*n_elt]);



  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);
  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
  PDM_part_mesh_nodal_free(pmn);
}


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
                           0.,
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
                           0.,
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
                                               &face_edge_idx,
                                               &face_edge,
                                               PDM_OWNERSHIP_KEEP);

  int  *edge_vtx      = NULL;
  int  *edge_vtx_idx  = NULL;
  PDM_multipart_part_connectivity_get(mpart,
                                      0,
                                      0,
                                      PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                      &edge_vtx_idx,
                                      &edge_vtx,
                                      PDM_OWNERSHIP_KEEP);

  PDM_compute_face_vtx_from_face_and_edge(*n_elt,
                                          face_edge_idx,
                                          face_edge,
                                          edge_vtx,
                                          elt_vtx);

  PDM_malloc(*elt_vtx_idx, ((*n_elt)+1), int);
  (*elt_vtx_idx)[0] = 0;
  for (int i = 0; i < (*n_elt); i++) {
    (*elt_vtx_idx)[i+1] = (*elt_vtx_idx)[i] + 3; // because PDM_MESH_NODAL_TRIA3
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

}


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


void
PDM_generate_mesh_parallelepiped_simplified
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
                                n_vtx_seg,
                                n_vtx_seg,
                                n_vtx_seg,
                                1,
                                PDM_SPLIT_DUAL_WITH_HILBERT,
                                &dmn,
                                &mpart);

  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_KEEP);

  *n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, 0);

  int         *connec              = NULL;
  PDM_g_num_t *numabs              = NULL;
  int         *parent_num          = NULL;
  PDM_g_num_t *parent_entity_g_num = NULL;
  PDM_part_mesh_nodal_section_std_get(pmn,
                                      0,
                                      0,
                                      &connec,
                                      &numabs,
                                      &parent_num,
                                      &parent_entity_g_num,
                                      PDM_OWNERSHIP_KEEP);

  *elt_vtx_idx = PDM_array_new_idx_from_const_stride_int(4, *n_elt);
  PDM_malloc(*elt_vtx, (*elt_vtx_idx)[*n_elt], int);
  memcpy(*elt_vtx, connec, sizeof(int) * (*elt_vtx_idx)[*n_elt]);

  // get coordinates
  *n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                            0,
                                            0,
                                            coords,
                                            PDM_OWNERSHIP_USER);

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);
  PDM_part_mesh_nodal_free(pmn);
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
 const double                  random_factor,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
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
                           random_factor,
                           &dmn,
                           &mpart);
  PDM_DMesh_nodal_free(dmn);

  PDM_malloc(*pn_vtx        , n_part, int          );
  PDM_malloc(*pn_edge       , n_part, int          );
  PDM_malloc(*pn_face       , n_part, int          );
  PDM_malloc(*pvtx_coord    , n_part, double      *);
  PDM_malloc(*pedge_vtx     , n_part, int         *);
  PDM_malloc(*pface_edge_idx, n_part, int         *);
  PDM_malloc(*pface_edge    , n_part, int         *);
  PDM_malloc(*pface_vtx     , n_part, int         *);
  PDM_malloc(*pvtx_ln_to_gn , n_part, PDM_g_num_t *);
  PDM_malloc(*pedge_ln_to_gn, n_part, PDM_g_num_t *);
  PDM_malloc(*pface_ln_to_gn, n_part, PDM_g_num_t *);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
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
                                        &(*pface_edge_idx)[ipart],
                                        &(*pface_edge)    [ipart],
                                        PDM_OWNERSHIP_USER);

    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &(*pedge_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
    if (edge_vtx_idx != NULL) {
      PDM_free(edge_vtx_idx);
    }

    int *face_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                        &face_vtx_idx,
                                        &(*pface_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
    PDM_free(face_vtx_idx);
  }

  PDM_multipart_free(mpart);
}

void
PDM_generate_mesh_sphere_ngon
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
 const PDM_split_dual_t       part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn
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
  PDM_DMesh_nodal_free(dmn);

  PDM_malloc(*pn_vtx        , n_part, int          );
  PDM_malloc(*pn_edge       , n_part, int          );
  PDM_malloc(*pn_face       , n_part, int          );
  PDM_malloc(*pvtx_coord    , n_part, double      *);
  PDM_malloc(*pedge_vtx     , n_part, int         *);
  PDM_malloc(*pface_edge_idx, n_part, int         *);
  PDM_malloc(*pface_edge    , n_part, int         *);
  PDM_malloc(*pface_vtx     , n_part, int         *);
  PDM_malloc(*pvtx_ln_to_gn , n_part, PDM_g_num_t *);
  PDM_malloc(*pedge_ln_to_gn, n_part, PDM_g_num_t *);
  PDM_malloc(*pface_ln_to_gn, n_part, PDM_g_num_t *);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
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
                                        &(*pface_edge_idx)[ipart],
                                        &(*pface_edge)    [ipart],
                                        PDM_OWNERSHIP_USER);

    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &(*pedge_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
    if (edge_vtx_idx != NULL) {
      PDM_free(edge_vtx_idx);
    }

    PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                            (*pface_edge_idx)[ipart],
                                            (*pface_edge)[ipart],
                                            (*pedge_vtx)[ipart],
                                            &(*pface_vtx)[ipart]);

  }

  PDM_multipart_free(mpart);
}

void
PDM_generate_mesh_ball_ngon
(
 const PDM_MPI_Comm            comm,
 PDM_Mesh_nodal_elt_t          elt_type,
 int                           order,
 const char                   *ho_ordering,
 const double                  radius,
 const double                  hole_radius,
 const double                  center_x,
 const double                  center_y,
 const double                  center_z,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const PDM_g_num_t             n_z,
 const PDM_g_num_t             n_layer,
 const double                  geometric_ratio,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 int                         **pn_cell,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pcell_ln_to_gn,
 int                         **pn_surface,
 int                        ***psurface_face_idx,
 int                        ***psurface_face,
 PDM_g_num_t                ***psurface_face_ln_to_gn
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
  PDM_DMesh_nodal_free(dmn);

  PDM_malloc(*pn_vtx                , n_part, int          );
  PDM_malloc(*pn_edge               , n_part, int          );
  PDM_malloc(*pn_face               , n_part, int          );
  PDM_malloc(*pn_cell               , n_part, int          );
  PDM_malloc(*pvtx_coord            , n_part, double      *);
  PDM_malloc(*pedge_vtx             , n_part, int         *);
  PDM_malloc(*pface_edge_idx        , n_part, int         *);
  PDM_malloc(*pface_edge            , n_part, int         *);
  PDM_malloc(*pface_vtx             , n_part, int         *);
  PDM_malloc(*pcell_face_idx        , n_part, int         *);
  PDM_malloc(*pcell_face            , n_part, int         *);
  PDM_malloc(*pvtx_ln_to_gn         , n_part, PDM_g_num_t *);
  PDM_malloc(*pedge_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pface_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pcell_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pn_surface            , n_part, int          );
  PDM_malloc(*psurface_face_idx     , n_part, int         *);
  PDM_malloc(*psurface_face         , n_part, int         *);
  PDM_malloc(*psurface_face_ln_to_gn, n_part, PDM_g_num_t *);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
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

    (*pn_cell)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_MESH_ENTITY_CELL,
                                                        &(*pcell_ln_to_gn)[ipart],
                                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &(*pcell_face_idx)[ipart],
                                        &(*pcell_face)    [ipart],
                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &(*pface_edge_idx)[ipart],
                                        &(*pface_edge)    [ipart],
                                        PDM_OWNERSHIP_USER);

    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &(*pedge_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
    if (edge_vtx_idx != NULL) {
      PDM_free(edge_vtx_idx);
    }

    PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                            (*pface_edge_idx)[ipart],
                                            (*pface_edge)[ipart],
                                            (*pedge_vtx)[ipart],
                                            &(*pface_vtx)[ipart]);

    PDM_multipart_group_get(mpart,
                            0,
                            ipart,
                            PDM_MESH_ENTITY_FACE,
                            &(*pn_surface)[ipart],
                            &(*psurface_face_idx)[ipart],
                            &(*psurface_face)[ipart],
                            &(*psurface_face_ln_to_gn)[ipart],
                            PDM_OWNERSHIP_USER);
  }

  PDM_multipart_free(mpart);
}



void
PDM_generate_mesh_parallelepiped_ngon
(
 const PDM_MPI_Comm            comm,
 PDM_Mesh_nodal_elt_t          elt_type,
 int                           order,
 const char                   *ho_ordering,
 const double                  xmin,
 const double                  ymin,
 const double                  zmin,
 const double                  lengthx,
 const double                  lengthy,
 const double                  lengthz,
 const PDM_g_num_t             n_x,
 const PDM_g_num_t             n_y,
 const PDM_g_num_t             n_z,
 const int                     n_part,
 const PDM_split_dual_t        part_method,
 int                         **pn_vtx,
 int                         **pn_edge,
 int                         **pn_face,
 int                         **pn_cell,
 double                     ***pvtx_coord,
 int                        ***pedge_vtx,
 int                        ***pface_edge_idx,
 int                        ***pface_edge,
 int                        ***pface_vtx,
 int                        ***pcell_face_idx,
 int                        ***pcell_face,
 PDM_g_num_t                ***pvtx_ln_to_gn,
 PDM_g_num_t                ***pedge_ln_to_gn,
 PDM_g_num_t                ***pface_ln_to_gn,
 PDM_g_num_t                ***pcell_ln_to_gn,
 int                         **pn_surface,
 int                        ***psurface_face_idx,
 int                        ***psurface_face,
 PDM_g_num_t                ***psurface_face_ln_to_gn,
 int                         **pn_ridge,
 int                        ***pridge_edge_idx,
 int                        ***pridge_edge,
 PDM_g_num_t                ***pridge_edge_ln_to_gn
 )
{
  PDM_dmesh_nodal_t *dmn = NULL;
  PDM_multipart_t *mpart = NULL;

  _generate_mesh_parallelepiped (comm,
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

  PDM_DMesh_nodal_free(dmn);

  PDM_malloc(*pn_vtx                , n_part, int          );
  PDM_malloc(*pn_edge               , n_part, int          );
  PDM_malloc(*pn_face               , n_part, int          );
  PDM_malloc(*pn_cell               , n_part, int          );
  PDM_malloc(*pvtx_coord            , n_part, double      *);
  PDM_malloc(*pedge_vtx             , n_part, int         *);
  PDM_malloc(*pface_edge_idx        , n_part, int         *);
  PDM_malloc(*pface_edge            , n_part, int         *);
  PDM_malloc(*pface_vtx             , n_part, int         *);
  PDM_malloc(*pcell_face_idx        , n_part, int         *);
  PDM_malloc(*pcell_face            , n_part, int         *);
  PDM_malloc(*pvtx_ln_to_gn         , n_part, PDM_g_num_t *);
  PDM_malloc(*pedge_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pface_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pcell_ln_to_gn        , n_part, PDM_g_num_t *);
  PDM_malloc(*pn_surface            , n_part, int          );
  PDM_malloc(*psurface_face_idx     , n_part, int         *);
  PDM_malloc(*psurface_face         , n_part, int         *);
  PDM_malloc(*psurface_face_ln_to_gn, n_part, PDM_g_num_t *);
  PDM_malloc(*pn_ridge              , n_part, int          );
  PDM_malloc(*pridge_edge_idx       , n_part, int         *);
  PDM_malloc(*pridge_edge           , n_part, int         *);
  PDM_malloc(*pridge_edge_ln_to_gn  , n_part, PDM_g_num_t *);

  for (int ipart = 0; ipart < n_part; ipart++) {
    (*pn_vtx)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                       0,
                                                       ipart,
                                                       PDM_MESH_ENTITY_VTX,
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

    (*pn_cell)[ipart] = PDM_multipart_part_ln_to_gn_get(mpart,
                                                        0,
                                                        ipart,
                                                        PDM_MESH_ENTITY_CELL,
                                                        &(*pcell_ln_to_gn)[ipart],
                                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &(*pcell_face_idx)[ipart],
                                        &(*pcell_face)    [ipart],
                                        PDM_OWNERSHIP_USER);

    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &(*pface_edge_idx)[ipart],
                                        &(*pface_edge)    [ipart],
                                        PDM_OWNERSHIP_USER);


    int *edge_vtx_idx = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        0,
                                        ipart,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx_idx,
                                        &(*pedge_vtx)[ipart],
                                        PDM_OWNERSHIP_USER);
    if (edge_vtx_idx != NULL) {
      PDM_free(edge_vtx_idx);
    }

    PDM_compute_face_vtx_from_face_and_edge((*pn_face)[ipart],
                                            (*pface_edge_idx)[ipart],
                                            (*pface_edge)[ipart],
                                            (*pedge_vtx)[ipart],
                                            &(*pface_vtx)[ipart]);

    PDM_multipart_group_get(mpart,
                            0,
                            ipart,
                            PDM_MESH_ENTITY_FACE,
                            &(*pn_surface)[ipart],
                            &(*psurface_face_idx)[ipart],
                            &(*psurface_face)[ipart],
                            &(*psurface_face_ln_to_gn)[ipart],
                            PDM_OWNERSHIP_USER);

    PDM_multipart_group_get(mpart,
                            0,
                            ipart,
                            PDM_MESH_ENTITY_EDGE,
                            &(*pn_ridge)[ipart],
                            &(*pridge_edge_idx)[ipart],
                            &(*pridge_edge)[ipart],
                            &(*pridge_edge_ln_to_gn)[ipart],
                            PDM_OWNERSHIP_USER);


  }

  PDM_multipart_free(mpart);

}



PDM_part_mesh_nodal_t *
PDM_generate_mesh_nodal_from_file
(
  const PDM_MPI_Comm      comm,
  const int               n_part,
  const PDM_split_dual_t  part_method,
  const char             *filename
)
{
  // Read and partition mesh
  PDM_dmesh_nodal_t *dmn   = NULL;
  PDM_multipart_t   *mpart = NULL;
  _read_mesh_file(comm,
                  n_part,
                  part_method,
                  filename,
                  &dmn,
                  &mpart);

  // Retrieve partitioned mesh
  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  // Free memory
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmn;
}



PDM_part_mesh_t *
PDM_generate_mesh_from_file
(
  const PDM_MPI_Comm      comm,
  const int               n_part,
  const PDM_split_dual_t  part_method,
  const char             *filename
)
{
  // Read and partition mesh
  PDM_dmesh_nodal_t *dmn   = NULL;
  PDM_multipart_t   *mpart = NULL;
  _read_mesh_file(comm,
                  n_part,
                  part_method,
                  filename,
                  &dmn,
                  &mpart);

  // Retrieve partitioned mesh
  PDM_part_mesh_t *pmesh = NULL;
  PDM_multipart_get_part_mesh(mpart,
                              0,
                              &pmesh,
                              PDM_OWNERSHIP_USER);

  // Free memory
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  return pmesh;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
