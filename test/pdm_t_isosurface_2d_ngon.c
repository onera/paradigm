#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"

#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_vtk.h"

#include "pdm_array.h"
#include "pdm_reader_gamma.h"
#include "pdm_dcube_nodal_gen.h"

#include "pdm_isosurface.h"
#include "pdm_multipart.h"
#include "pdm_multipart_priv.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_dconnectivity_transform.h"

#include "pdm_isosurface_test_utils.c"

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int main
(
  int   argc,
  char *argv[]
)
{
  /*
   *  Init MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Init(&argc, &argv);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  /*
   *  Read args
   */
  char                 *mesh_name      = NULL;
  char                 *sol_name       = NULL;
  int                   n_part         = 1;
  int                   visu           = 0;
  int                   n_isovalues    = 1;
  double               *isovalues      = NULL;
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_QUAD4;
  PDM_g_num_t           n_vtx_seg      = 10;
  int                   randomize      = 0;
  int                   use_part_mesh  = 0;
  int                   generate_edges = 0;
  int                   local          = 0;

  PDM_isosurface_test_utils_read_args(argc,
                                      argv,
                                     &n_part,
                                     &mesh_name,
                                     &sol_name,
                                     &visu,
                                     &n_isovalues,
                                     &isovalues,
                                     &elt_type,
                                     &randomize,
                                     &n_vtx_seg,
                                     &use_part_mesh,
                                     &generate_edges,
                                     &local);

  if (isovalues == NULL) {
    n_isovalues = 1;
    PDM_malloc(isovalues, n_isovalues, double);
    isovalues[0] = 0.;
  }

  // Generate/read mesh
  PDM_multipart_t *mpart = NULL;
  PDM_part_mesh_t *pmesh = NULL;
  PDM_dmesh_t     *dmesh = NULL;
  if (n_part > 0) {
    if (use_part_mesh) {
      pmesh = PDM_part_mesh_create(n_part, comm);
    }
  }

  PDM_isosurface_test_utils_gen_mesh(comm,
                                     mesh_name,
                                     n_part,
                                     n_vtx_seg,
                                     randomize,
                                     elt_type,
                                     generate_edges,
                                    &mpart,
                                     pmesh,
                                    &dmesh);


  /**
   * Compute scalar field
   */
  double  *iso_dfield      = NULL;
  double  *itp_dfield_vtx  = NULL;
  double  *itp_dfield_face = NULL;
  double **iso_field       = NULL;
  double **itp_field_vtx   = NULL;
  double **itp_field_face  = NULL;
  if (n_part > 0) {
    // Partitioned
    PDM_malloc(iso_field     , n_part, double *);
    PDM_malloc(itp_field_vtx , n_part, double *);
    PDM_malloc(itp_field_face, n_part, double *);
    for (int i_part = 0; i_part < n_part; i_part++) {

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_malloc(iso_field     [i_part], n_vtx, double);
      PDM_malloc(itp_field_vtx [i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field    [i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field_vtx[i_part]);

      PDM_g_num_t *face_gnum = NULL;
      int n_face = PDM_multipart_part_ln_to_gn_get(mpart, 0, i_part, PDM_MESH_ENTITY_FACE,
                                                  &face_gnum, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_face[i_part], n_face, double);
      for (int i_face=0; i_face<n_face; ++i_face) {
        itp_field_face[i_part][i_face] = (double) face_gnum[i_face];
      }
    }
  }
  else {
    // Block-distributed
    int dn_vtx = PDM_dmesh_dn_entity_get(dmesh, PDM_MESH_ENTITY_VTX);

    double *dvtx_coord = NULL;
    PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);

    PDM_malloc(iso_dfield    , dn_vtx, double);
    PDM_malloc(itp_dfield_vtx, dn_vtx, double);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, iso_dfield);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, itp_dfield_vtx);


    PDM_g_num_t *face_distri = NULL;
    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &face_distri);
    int dn_face = face_distri[i_rank+1]-face_distri[i_rank];
    PDM_malloc(itp_dfield_face, dn_face, double);
    for (int i_face=0; i_face<dn_face; ++i_face) {
      itp_dfield_face[i_face] = (double) (face_distri[i_rank]+i_face);
    }
  }


  /**
   * Create Isosurface instance
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 2);

  if (n_part > 0) {
    PDM_extract_part_kind_t extract_kind = PDM_EXTRACT_PART_KIND_REEQUILIBRATE;
    if (local) {
      extract_kind = PDM_EXTRACT_PART_KIND_LOCAL;
    }
    PDM_isosurface_redistribution_set(isos,
                                      extract_kind,
                                      PDM_SPLIT_DUAL_WITH_HILBERT);
  }


  /* Set mesh */
  int         *dface_vtx_idx = NULL;
  PDM_g_num_t *dface_vtx     = NULL;
  int own_dface_vtx = 0;
  if (n_part > 0) {
    // Partitioned
    if (use_part_mesh) {
      PDM_isosurface_part_mesh_set(isos, pmesh);
    }
    else {
      PDM_isosurface_n_part_set(isos, n_part);

      for (int i_part = 0; i_part < n_part; i_part++) {

        // Connectivities
        if (generate_edges) {
          int *face_edge_idx = NULL;
          int *face_edge     = NULL;
          int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                           0,
                                                           i_part,
                                                           PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                                           &face_edge_idx,
                                                           &face_edge,
                                                           PDM_OWNERSHIP_KEEP);
          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                          n_face,
                                          face_edge_idx,
                                          face_edge);

          int *edge_vtx_idx = NULL;
          int *edge_vtx     = NULL;
          int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                           0,
                                                           i_part,
                                                           PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                           &edge_vtx_idx,
                                                           &edge_vtx,
                                                           PDM_OWNERSHIP_KEEP);
          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          n_edge,
                                          NULL,
                                          edge_vtx);
        }
        else {
          int *face_vtx_idx = NULL;
          int *face_vtx     = NULL;
          int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                           0,
                                                           i_part,
                                                           PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                           &face_vtx_idx,
                                                           &face_vtx,
                                                           PDM_OWNERSHIP_KEEP);

          PDM_isosurface_connectivity_set(isos,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                          n_face,
                                          face_vtx_idx,
                                          face_vtx);
        }



        // Coordinates
        double *vtx_coord = NULL;
        int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
        PDM_isosurface_vtx_coord_set(isos,
                                     i_part,
                                     n_vtx,
                                     vtx_coord);

        // Global IDs
        PDM_g_num_t *face_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *vtx_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_VTX,
                                        &vtx_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    vtx_ln_to_gn);

        if (generate_edges) {
          PDM_g_num_t *edge_ln_to_gn = NULL;
          PDM_multipart_part_ln_to_gn_get(mpart,
                                          0,
                                          i_part,
                                          PDM_MESH_ENTITY_EDGE,
                                          &edge_ln_to_gn,
                                          PDM_OWNERSHIP_KEEP);
          PDM_isosurface_ln_to_gn_set(isos,
                                      i_part,
                                      PDM_MESH_ENTITY_EDGE,
                                      edge_ln_to_gn);
        }
      }
    }
  }
  else {
    // Block-distributed
    if (use_part_mesh) {
      PDM_isosurface_dmesh_set(isos, dmesh);
    }
    else {
      for (int i_entity = PDM_MESH_ENTITY_FACE; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *distrib = NULL;
        PDM_dmesh_distrib_get     (dmesh, i_entity, &distrib);
        PDM_isosurface_distrib_set(isos,  i_entity,  distrib);
      }

      int         *dface_edge_idx = NULL;
      PDM_g_num_t *dface_edge     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 &dface_edge,
                                 &dface_edge_idx,
                                 PDM_OWNERSHIP_KEEP);
      int         *dedge_vtx_idx = NULL;
      PDM_g_num_t *dedge_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                 &dedge_vtx,
                                 &dedge_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);
      if (generate_edges) {
        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         dface_edge_idx,
                                         dface_edge);

        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         NULL,
                                         dedge_vtx);
      }
      else {
        PDM_dmesh_connectivity_get(dmesh,
                                   PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                   &dface_vtx,
                                   &dface_vtx_idx,
                                   PDM_OWNERSHIP_KEEP);

        if (dface_vtx == NULL) {
          own_dface_vtx = 1;
          PDM_g_num_t *distrib_face = NULL;
          PDM_g_num_t *distrib_edge = NULL;
          PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_FACE, &distrib_face);
          PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_EDGE, &distrib_edge);
          PDM_dconnectivity_dface_vtx_from_face_and_edge(comm,
                                                         distrib_face,
                                                         distrib_edge,
                                                         dface_edge_idx,
                                                         dface_edge,
                                                         dedge_vtx,
                                                         &dface_vtx);
          dface_vtx_idx = dface_edge_idx;
        }

        PDM_isosurface_dconnectivity_set(isos,
                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         dface_vtx_idx,
                                         dface_vtx);
      }

      double *dvtx_coord = NULL;
      PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dvtx_coord_set(isos, dvtx_coord);
    }
  }


  // Plane slice
  double plane_equation [4] = {1.,0.,0.,0.5};
  double plane_isovalues[3] = {-0.30,0.,0.30};
  int iso1 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_PLANE,
                                3,
                                plane_isovalues);

  PDM_isosurface_equation_set(isos,
                              iso1,
                              plane_equation);

  // Scalar field isosurface
  int iso2 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_FIELD,
                                n_isovalues,
                                isovalues);

  if (n_part > 0) {
    // Partitioned
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_isosurface_field_set(isos,
                               iso2,
                               i_part,
                               iso_field[i_part]);
    }
  }
  else {
    // Block-distributed
    PDM_isosurface_dfield_set(isos,
                              iso2,
                              iso_dfield);
  }

  int n_iso = iso2 + 1;
  for (int i_iso = 0; i_iso < n_iso; i_iso++) {
    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_VTX,
                                       0);

    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_EDGE,
                                       0);
  }

  /*
   *  Compute isosurface
   */
  // PDM_isosurface_compute(isos, iso1);
  // PDM_isosurface_reset  (isos, iso1);
  // PDM_isosurface_compute(isos, iso1);
  // PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, -1);
  PDM_isosurface_reset  (isos, -1);
  PDM_isosurface_compute(isos, -1);



  /*
   *  Interpolate field
   */
  double  **iso_itp_dfield_vtx  = NULL;
  double  **iso_itp_dfield_edge = NULL;
  double ***iso_itp_field_vtx   = NULL;
  double ***iso_itp_field_edge  = NULL;

  if (n_part > 0) {
    PDM_malloc(iso_itp_field_vtx  , n_iso, double **);
    PDM_malloc(iso_itp_field_edge , n_iso, double **);
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
      PDM_isosurface_test_utils_part_interpolation(isos, i_iso, n_part, local,
                                                   itp_field_vtx ,
                                                   itp_field_face,
                                                   NULL,
                                                  &iso_itp_field_vtx [i_iso],
                                                  &iso_itp_field_edge[i_iso],
                                                   NULL);
    }
  }
  else {
    PDM_malloc(iso_itp_dfield_vtx , n_iso, double  *);
    PDM_malloc(iso_itp_dfield_edge, n_iso, double  *);
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
      PDM_isosurface_test_utils_dist_interpolation(isos, i_iso,
                                                   itp_dfield_vtx ,
                                                   itp_dfield_face,
                                                   NULL,
                                                  &iso_itp_dfield_vtx [i_iso],
                                                  &iso_itp_dfield_edge[i_iso],
                                                   NULL);
    }
  }



  /*
   *  Visu isosurfaces
   */
  if (visu) {

    for (int i_iso = 0; i_iso < n_iso; i_iso++) {

      if (n_part > 0) {
        PDM_isosurface_test_utils_part_vtk(isos, i_iso, n_part,
                                           iso_itp_field_vtx [i_iso],
                                           iso_itp_field_edge[i_iso],
                                           NULL,
                                           comm);
      }
      else {
        PDM_isosurface_test_utils_dist_vtk(isos, i_iso, 
                                           iso_itp_dfield_vtx [i_iso],
                                           iso_itp_dfield_edge[i_iso],
                                           NULL,
                                           comm);
      }

    }
  }



  /* Free memory */
  PDM_isosurface_free(isos);
  PDM_free(isovalues);

  if (n_part > 0) {
    PDM_multipart_free(mpart);

    if (use_part_mesh) {
      PDM_part_mesh_free(pmesh);
    }

    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_free(iso_field     [i_part]);
      PDM_free(itp_field_vtx [i_part]);
      PDM_free(itp_field_face[i_part]);
    }
    PDM_free(iso_field);
    PDM_free(itp_field_vtx);
    PDM_free(itp_field_face);
    
    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      for (int i_part = 0; i_part < n_part; i_part++) {
        PDM_free(iso_itp_field_vtx [i_iso][i_part]);
        PDM_free(iso_itp_field_edge[i_iso][i_part]);
      }
      PDM_free(iso_itp_field_vtx [i_iso]);
      PDM_free(iso_itp_field_edge[i_iso]);
    }
    PDM_free(iso_itp_field_vtx);
    PDM_free(iso_itp_field_edge);
  }
  else {
    if (own_dface_vtx) {
      PDM_free(dface_vtx);
    }
    PDM_free(iso_dfield);
    PDM_free(itp_dfield_vtx);
    PDM_free(itp_dfield_face);
    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      PDM_free(iso_itp_dfield_vtx [i_iso]);
      // PDM_free(iso_itp_dfield_edge[i_iso]);
    }
    PDM_free(iso_itp_dfield_vtx);
    PDM_free(iso_itp_dfield_edge);
    PDM_dmesh_free(dmesh);
  }


  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }



  PDM_MPI_Finalize();


  return EXIT_SUCCESS;
}
