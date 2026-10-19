#ifndef __PDM_ISOSURFACE_TEST_UTILS_H__
#define __PDM_ISOSURFACE_TEST_UTILS_H__

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_dmesh.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_isosurface.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_mesh.h"
#include "pdm_part_mesh_nodal.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct PDM_isosurface_test_utils_params_t {

  char                 *mesh_name;
  int                   n_part_in;
  int                   n_part_out;
  int                   n_isovalues;
  double               *isovalues;
  PDM_Mesh_nodal_elt_t  elt_type;
  int                   randomize;
  PDM_g_num_t           n_vtx_seg;
  int                   use_part_mesh;
  int                   generate_edges;
  int                   local;
  int                   use_groups;
  int                   visu;

} PDM_isosurface_test_utils_params_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Read arguments from the command line
 *
 */

void
PDM_isosurface_test_utils_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_part_in,
  int                   *n_part_out,
  char                 **mesh_name,
  int                   *visu,
  int                   *n_isovalues,
  double               **isovalues,
  PDM_Mesh_nodal_elt_t  *elt_type,
  int                   *randomize,
  PDM_g_num_t           *n_vtx_seg,
  int                   *use_part_mesh,
  int                   *generate_edges,
  int                   *local,
  int                   *use_groups
);


/**
 *
 * \brief Mesh generation for ngon cases
 *
 */

void
PDM_isosurface_test_utils_analytic_field_function
(
 const double  x,
 const double  y,
 const double  z,
 double       *value
);


/**
 *
 * \brief Mesh generation for ngon cases
 *
 */

void
PDM_isosurface_test_utils_gen_mesh
(
  PDM_MPI_Comm          comm,
  const char           *filename,
  int                   n_part,
  PDM_g_num_t           n_vtx_seg,
  int                   randomize,
  PDM_Mesh_nodal_elt_t  elt_type,
  int                   generate_edges,
  int                   use_groups,
  PDM_multipart_t     **mpart,
  PDM_part_mesh_t      *pmesh,
  PDM_dmesh_t         **out_dmesh
);


/**
 *
 * \brief Mesh generation for nodal cases
 *
 */

void
PDM_isosurface_test_utils_gen_mesh_nodal
(
  PDM_MPI_Comm            comm,
  const char             *filename,
  int                     n_part,
  PDM_g_num_t             n_vtx_seg,
  int                     randomize,
  PDM_Mesh_nodal_elt_t    elt_type,
  int                     use_groups,
  PDM_part_mesh_nodal_t **out_pmn,
  PDM_dmesh_nodal_t     **out_dmn
);


/**
 *
 * \brief Compute iso field from vertex coordinates
 *
 */

void
PDM_isosurface_test_utils_compute_iso_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
);


/**
 *
 * \brief Compute interpolate field from vertex coordinates
 *
 */

void
PDM_isosurface_test_utils_compute_itp_field
(
  int     n_vtx,
  double *vtx_coord,
  double *vtx_field
);


/**
 *
 * \brief Interpolation source field onto dist iso mesh
 *
 */

void
PDM_isosurface_test_utils_dist_interpolation
(
  PDM_isosurface_t *isos,
  int               i_iso,
  double           *itp_dfield_vtx,
  double           *itp_dfield_face,
  double           *itp_dfield_cell,
  double          **iso_itp_dfield_vtx,
  double          **iso_itp_dfield_edge,
  double          **iso_itp_dfield_face
);


/**
 *
 * \brief Interpolation source field onto part iso mesh
 *
 */

void
PDM_isosurface_test_utils_part_interpolation
(
  PDM_isosurface_t *isos,
  int               i_iso,
  int               n_part,
  int               local,
  double          **itp_field_vtx,
  double          **itp_field_face,
  double          **itp_field_cell,
  double         ***iso_itp_field_vtx,
  double         ***iso_itp_field_edge,
  double         ***iso_itp_field_face
);



/**
 *
 * \brief Write vtk output from isosurface dist result
 *
 */

void
PDM_isosurface_test_utils_dist_vtk
(
  PDM_isosurface_t *isos,
  int               i_iso,
  const double     *iso_vtx_fld,
  const double     *iso_edge_fld,
  const double     *iso_face_fld,
  PDM_MPI_Comm      comm
);


/**
 *
 * \brief Write vtk output from isosurface part result
 *
 */

void
PDM_isosurface_test_utils_part_vtk
(
  PDM_isosurface_t  *isos,
  int                id_iso,
  int                n_part,
  double           **iso_vtx_fld,
  double           **iso_edge_fld,
  double           **iso_face_fld,
  PDM_MPI_Comm       comm
);


/**
 *
 * \brief Get isosurface size (part or dist (n_part <= 0))
 *
 */

void
PDM_isosurface_test_utils_isosurface_size_get
(
  PDM_isosurface_t   *isos,
  int                 id_iso,
  int                 n_part,
  PDM_g_num_t        *gn_iso_vtx,
  PDM_g_num_t        *gn_iso_edge,
  PDM_g_num_t        *gn_iso_face,
  PDM_MPI_Comm        comm
);


/**
 * \brief Dump test parameters and equivalent command line
 */
void
PDM_isosurface_test_utils_isosurface_params_dump
(
  PDM_MPI_Comm                        comm,
  const char                         *test_name,
  PDM_isosurface_test_utils_params_t  params
);


#ifdef  __cplusplus
}
#endif

#endif /* __PDM_ISOSURFACE_TEST_UTILS_H__ */
