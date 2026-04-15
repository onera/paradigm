#include <stdio.h>
#include <stdlib.h>

#include "pdm_dmesh.h"
#include "pdm_error.h"
#include "pdm_isosurface_test_utils.h"
#include "pdm_isosurface.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_part_mesh.h"
#include "pdm.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

static void _run_test_3d_ngon(PDM_isosurface_test_utils_params_t params);

/*============================================================================
 * Public function definitions
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
  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);


  double default_isovalues[] = {0, -0.25, 0.25};

  PDM_isosurface_test_utils_params_t params;
  int n_tests = 0;

  if (argc <= 1) {
    /* Run matrix of tests */

    // Set fixed parameters
    params.mesh_name = NULL;
    params.n_vtx_seg = 10;
    params.randomize = 0;
    params.elt_type  = PDM_MESH_NODAL_HEXA8;
    params.visu      = 0;

    // Test combinations of variable parameters
    for (int n_part_in = 0; n_part_in <= 2; n_part_in++) {

      params.n_part_in = n_part_in;

      for (int n_isovalues = 1; n_isovalues <= 3; n_isovalues += 2) {

        params.n_isovalues = n_isovalues;
        params.isovalues   = default_isovalues;

        for (int use_part_mesh = 0; use_part_mesh <= 1; use_part_mesh++) {

          params.use_part_mesh = use_part_mesh;

          for (int generate_edges = 0; generate_edges <= 1; generate_edges++) {

            params.generate_edges = generate_edges;

            for (int local = 0; local <= (n_part_in > 0); local++) {

              params.local = local;

              int n_part_out_min = n_part_in;
              int n_part_out_max = n_part_in;
              if (n_part_in > 0 && !local) {
                n_part_out_min = 1;
              }

              for (int n_part_out = n_part_out_min; n_part_out <= n_part_out_max; n_part_out++) {

                params.n_part_out = n_part_out;

                for (int use_groups = 0; use_groups <= 1; use_groups++) {

                  params.use_groups = use_groups;

                  _run_test_3d_ngon(params);
                  n_tests++;

                } // End loop on use_group

              } // End loop on n_part_out

            } // End loop on local/redistribute

          } // End loop on generate_edges

        } // End loop on use_part_mesh

      } // End loop on n_isovalues

    } // End loop on n_part_in
  }

  else {
    /* Run test with user-provided parameters */

    // Initialize parameters to their default values
    params.mesh_name      = NULL;
    params.n_part_in      = 1;
    params.n_part_out     = 1;
    params.visu           = 0;
    params.n_isovalues    = 1;
    params.isovalues      = default_isovalues;
    params.elt_type       = PDM_MESH_NODAL_HEXA8;
    params.n_vtx_seg      = 10;
    params.randomize      = 0;
    params.use_part_mesh  = 0;
    params.generate_edges = 0;
    params.local          = 0;
    params.use_groups     = 0;

    // Overwrite with user parameter values
    PDM_isosurface_test_utils_read_args(argc,
                                        argv,
                                        &params.n_part_in,
                                        &params.n_part_out,
                                        &params.mesh_name,
                                        &params.visu,
                                        &params.n_isovalues,
                                        &params.isovalues,
                                        &params.elt_type,
                                        &params.randomize,
                                        &params.n_vtx_seg,
                                        &params.use_part_mesh,
                                        &params.generate_edges,
                                        &params.local,
                                        &params.use_groups);

    // Run test
    _run_test_3d_ngon(params);
    n_tests++;

    if (params.isovalues != default_isovalues) {
      PDM_free(params.isovalues);
    }
  }


  if (i_rank == 0) {
    printf("\nAll %d tests passed :D\n", n_tests);
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return EXIT_SUCCESS;
}


/**
 * Run test driven by custom parameters
 */
static void
_run_test_3d_ngon
(
  PDM_isosurface_test_utils_params_t params
)
{
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  PDM_isosurface_test_utils_isosurface_params_dump(comm, "pdm_t_isosurface_3d_ngon", params);
  PDM_MPI_Barrier(comm);


  char                 *mesh_name      = params.mesh_name;
  int                   n_part_in      = params.n_part_in;
  int                   n_part_out     = params.n_part_out;
  int                   visu           = params.visu;
  int                   n_isovalues    = params.n_isovalues;
  double               *isovalues      = params.isovalues;
  PDM_Mesh_nodal_elt_t  elt_type       = params.elt_type;
  PDM_g_num_t           n_vtx_seg      = params.n_vtx_seg;
  int                   randomize      = params.randomize;
  int                   use_part_mesh  = params.use_part_mesh;
  int                   generate_edges = params.generate_edges;
  int                   local          = params.local;
  int                   use_groups     = params.use_groups;


  /*
   *  Generate mesh
   */
  PDM_multipart_t *mpart = NULL;
  PDM_part_mesh_t *pmesh = NULL;
  PDM_dmesh_t     *dmesh = NULL;
  if (n_part_in > 0) {
    if (use_part_mesh) {
      pmesh = PDM_part_mesh_create(n_part_in, comm);
    }
  }

  PDM_isosurface_test_utils_gen_mesh(comm,
                                     mesh_name,
                                     n_part_in,
                                     n_vtx_seg,
                                     randomize,
                                     elt_type,
                                     generate_edges,
                                     use_groups,
                                    &mpart,
                                     pmesh,
                                    &dmesh);


  /**
   * Compute scalar field
   */
  double  *iso_dfield      = NULL;
  double  *itp_dfield_vtx  = NULL;
  double  *itp_dfield_face = NULL;
  double  *itp_dfield_cell = NULL;
  double **iso_field       = NULL;
  double **itp_field_vtx   = NULL;
  double **itp_field_face  = NULL;
  double **itp_field_cell  = NULL;
  if (n_part_in > 0) {
    // Partitioned
    PDM_malloc(iso_field     , n_part_in, double *);
    PDM_malloc(itp_field_vtx , n_part_in, double *);
    PDM_malloc(itp_field_face, n_part_in, double *);
    PDM_malloc(itp_field_cell, n_part_in, double *);
    for (int i_part = 0; i_part < n_part_in; i_part++) {

      double *vtx_coord = NULL;
      int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                   0,
                                                   i_part,
                                                   &vtx_coord,
                                                   PDM_OWNERSHIP_KEEP);

      PDM_malloc(iso_field    [i_part], n_vtx, double);
      PDM_malloc(itp_field_vtx[i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field    [i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field_vtx[i_part]);

      PDM_g_num_t *face_gnum = NULL;
      int n_face = PDM_multipart_part_ln_to_gn_get(mpart, 0, i_part, PDM_MESH_ENTITY_FACE,
                                                  &face_gnum, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_face[i_part], n_face, double);
      for (int i_face = 0; i_face < n_face; ++i_face) {
        itp_field_face[i_part][i_face] = (double) face_gnum[i_face];
      }

      PDM_g_num_t *cell_gnum = NULL;
      int n_cell = PDM_multipart_part_ln_to_gn_get(mpart, 0, i_part, PDM_MESH_ENTITY_CELL,
                                                  &cell_gnum, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_cell[i_part], n_cell, double);
      for (int i_cell = 0; i_cell < n_cell; ++i_cell) {
        itp_field_cell[i_part][i_cell] = (double) cell_gnum[i_cell];
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
    for (int i_face = 0; i_face < dn_face; ++i_face) {
      itp_dfield_face[i_face] = (double) (face_distri[i_rank]+i_face);
    }

    PDM_g_num_t *cell_distri = NULL;
    PDM_dmesh_distrib_get(dmesh, PDM_MESH_ENTITY_CELL, &cell_distri);
    int dn_cell = cell_distri[i_rank+1]-cell_distri[i_rank];
    PDM_malloc(itp_dfield_cell, dn_cell, double);
    for (int i_cell = 0; i_cell < dn_cell; ++i_cell) {
      itp_dfield_cell[i_cell] = (double) (cell_distri[i_rank]+i_cell);
    }
  }


  /**
   * Create Isosurface instance
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 3);

  if (n_part_in > 0 && local == 0) {
    PDM_isosurface_redistribution_set(isos, PDM_EXTRACT_PART_KIND_REEQUILIBRATE, PDM_SPLIT_DUAL_WITH_HILBERT);
    PDM_isosurface_n_part_out_set(isos, n_part_out);
  }


  /* Set mesh */
  if (n_part_in > 0) {
    // Partitioned
    if (use_part_mesh) {
      PDM_isosurface_part_mesh_set(isos, pmesh);
    }
    else {
      PDM_isosurface_n_part_set(isos, n_part_in);

      for (int i_part = 0; i_part < n_part_in; i_part++) {

        // Connectivities
        int *cell_face_idx = NULL;
        int *cell_face     = NULL;
        int n_cell = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                                         &cell_face_idx,
                                                         &cell_face,
                                                         PDM_OWNERSHIP_KEEP);

        int *face_vtx_idx = NULL;
        int *face_vtx     = NULL;
        int n_face = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                                         &face_vtx_idx,
                                                         &face_vtx,
                                                         PDM_OWNERSHIP_KEEP);

        int *face_edge_idx = NULL;
        int *face_edge     = NULL;
        PDM_multipart_part_connectivity_get(mpart,
                                            0,
                                            i_part,
                                            PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                            &face_edge_idx,
                                            &face_edge,
                                            PDM_OWNERSHIP_KEEP);

        int *edge_vtx_idx = NULL;
        int *edge_vtx     = NULL;
        int n_edge = PDM_multipart_part_connectivity_get(mpart,
                                                         0,
                                                         i_part,
                                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                                         &edge_vtx_idx,
                                                         &edge_vtx,
                                                         PDM_OWNERSHIP_KEEP);

        PDM_isosurface_pconnectivity_set(isos,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                         n_cell,
                                         cell_face_idx,
                                         cell_face);

        PDM_isosurface_pconnectivity_set(isos,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                         n_face,
                                         face_vtx_idx,
                                         face_vtx);

        PDM_isosurface_pconnectivity_set(isos,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                         n_face,
                                         face_edge_idx,
                                         face_edge);

        PDM_isosurface_pconnectivity_set(isos,
                                         i_part,
                                         PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                         n_edge,
                                         NULL,
                                         edge_vtx);

        // Coordinates
        double *vtx_coord = NULL;
        int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                     0,
                                                     i_part,
                                                     &vtx_coord,
                                                     PDM_OWNERSHIP_KEEP);
        PDM_isosurface_pvtx_coord_set(isos,
                                      i_part,
                                      n_vtx,
                                      vtx_coord);

        // Global IDs
        PDM_g_num_t *cell_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_CELL,
                                        &cell_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *face_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_FACE,
                                        &face_ln_to_gn,
                                        PDM_OWNERSHIP_KEEP);

        PDM_g_num_t *edge_ln_to_gn = NULL;
        PDM_multipart_part_ln_to_gn_get(mpart,
                                        0,
                                        i_part,
                                        PDM_MESH_ENTITY_EDGE,
                                        &edge_ln_to_gn,
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
                                    PDM_MESH_ENTITY_CELL,
                                    cell_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    face_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_EDGE,
                                    edge_ln_to_gn);

        PDM_isosurface_ln_to_gn_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_VTX,
                                    vtx_ln_to_gn);

        if (use_groups) {
          // Groups
          int          n_surface             = 0;
          int         *surface_face_idx      = NULL;
          int         *surface_face          = NULL;
          PDM_g_num_t *surface_face_ln_to_gn = NULL;
          PDM_multipart_group_get(mpart,
                                  0,
                                  i_part,
                                  PDM_MESH_ENTITY_FACE,
                                  &n_surface,
                                  &surface_face_idx,
                                  &surface_face,
                                  &surface_face_ln_to_gn,
                                  PDM_OWNERSHIP_KEEP);

          PDM_isosurface_n_group_set(isos,
                                     PDM_MESH_ENTITY_FACE,
                                     n_surface);

          PDM_isosurface_pgroup_set(isos,
                                    i_part,
                                    PDM_MESH_ENTITY_FACE,
                                    surface_face_idx,
                                    surface_face,
                                    surface_face_ln_to_gn);
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
      for (int i_entity = PDM_MESH_ENTITY_CELL; i_entity < PDM_MESH_ENTITY_MAX; i_entity++) {
        PDM_g_num_t *distrib = NULL;
        PDM_dmesh_distrib_get     (dmesh, (PDM_mesh_entities_t) i_entity, &distrib);
        PDM_isosurface_distrib_set(isos,  (PDM_mesh_entities_t) i_entity,  distrib);
      }

      int         *dcell_face_idx = NULL;
      PDM_g_num_t *dcell_face     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                 &dcell_face,
                                 &dcell_face_idx,
                                 PDM_OWNERSHIP_KEEP);

      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                       dcell_face_idx,
                                       dcell_face);

      int         *dface_edge_idx = NULL;
      PDM_g_num_t *dface_edge     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                 &dface_edge,
                                 &dface_edge_idx,
                                 PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                       dface_edge_idx,
                                       dface_edge);

      int         *dface_vtx_idx = NULL;
      PDM_g_num_t *dface_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                 &dface_vtx,
                                 &dface_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_FACE_VTX,
                                       dface_vtx_idx,
                                       dface_vtx);

      int         *dedge_vtx_idx = NULL;
      PDM_g_num_t *dedge_vtx     = NULL;
      PDM_dmesh_connectivity_get(dmesh,
                                 PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                 &dedge_vtx,
                                 &dedge_vtx_idx,
                                 PDM_OWNERSHIP_KEEP);

      PDM_isosurface_dconnectivity_set(isos,
                                       PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                       NULL,
                                       dedge_vtx);

      double *dvtx_coord = NULL;
      PDM_dmesh_vtx_coord_get(dmesh, &dvtx_coord, PDM_OWNERSHIP_KEEP);
      PDM_isosurface_dvtx_coord_set(isos, dvtx_coord);

      int         *dsurface_face_idx = NULL;
      PDM_g_num_t *dsurface_face     = NULL;
      int n_surface = PDM_dmesh_bound_get(dmesh,
                                          PDM_BOUND_TYPE_FACE,
                                          &dsurface_face,
                                          &dsurface_face_idx,
                                          PDM_OWNERSHIP_KEEP);

      if (use_groups) {
        PDM_isosurface_n_group_set(isos,
                                   PDM_MESH_ENTITY_FACE,
                                   n_surface);

        PDM_isosurface_dgroup_set(isos,
                                  PDM_MESH_ENTITY_FACE,
                                  dsurface_face_idx,
                                  dsurface_face);
      }
    }
  }



  /*
   *  Add isosurface parameters
   */

  // Plane slice
  double plane_equation [4] = {1.,-1.,0.,0.};
  double plane_isovalues[3] = {-0.30,0.,1.};
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
  if (n_part_in > 0) { // Partitioned
    for (int i_part = 0; i_part < n_part_in; i_part++) {
      PDM_isosurface_pfield_set(isos, iso2, i_part, iso_field[i_part]);
    }
  }
  else { // Block-distributed
    PDM_isosurface_dfield_set(isos, iso2, iso_dfield);
  }

  // Analytic field isosurface
  double iso3_isovalue = 0.3;
  int iso3 = PDM_isosurface_add(isos,
                                PDM_ISO_SURFACE_KIND_FUNCTION,
                                1,
                                &iso3_isovalue);

  PDM_isosurface_field_function_set(isos,
                                    iso3,
                                    &PDM_isosurface_test_utils_analytic_field_function);



  /*
   *  Compute isosurface
   */
  int n_iso = iso3 + 1;

  for (int i_iso = 0; i_iso < n_iso; i_iso++) {
    PDM_isosurface_part_to_part_enable(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_VTX,
                                       0);

    PDM_isosurface_part_to_part_enable(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_EDGE,
                                       0);

    PDM_isosurface_part_to_part_enable(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_FACE,
                                       0);
  }

  PDM_MPI_Barrier(comm);

  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset  (isos, iso1);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, -1);
  PDM_isosurface_reset  (isos, -1);
  PDM_isosurface_compute(isos, -1);

  PDM_isosurface_dump_times(isos);

  for (int i_iso = 0; i_iso < n_iso; i_iso++) {
    PDM_g_num_t gn_iso_vtx, gn_iso_edge, gn_iso_face;
    PDM_isosurface_test_utils_isosurface_size_get(isos,
                                                  i_iso,
                                                  n_part_out,
                                                  &gn_iso_vtx,
                                                  &gn_iso_edge,
                                                  &gn_iso_face,
                                                  comm);
    if (i_rank == 0) {
      printf("isosurface %d : "PDM_FMT_G_NUM" vtx, "PDM_FMT_G_NUM" edges, "PDM_FMT_G_NUM" faces\n",
             i_iso, gn_iso_vtx, gn_iso_edge, gn_iso_face);
      fflush(stdout);
    }
  }


  /*
   *  Interpolate field
   */
  double  **iso_itp_dfield_vtx  = NULL;
  double  **iso_itp_dfield_edge = NULL;
  double  **iso_itp_dfield_face = NULL;
  double ***iso_itp_field_vtx   = NULL;
  double ***iso_itp_field_edge  = NULL;
  double ***iso_itp_field_face  = NULL;

  if (n_part_out > 0) {
    PDM_malloc(iso_itp_field_vtx , n_iso, double **);
    PDM_malloc(iso_itp_field_edge, n_iso, double **);
    PDM_malloc(iso_itp_field_face, n_iso, double **);
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
      PDM_isosurface_test_utils_part_interpolation(isos, i_iso, n_part_out, local,
                                                   itp_field_vtx,
                                                   itp_field_face,
                                                   itp_field_cell,
                                                  &iso_itp_field_vtx [i_iso],
                                                  &iso_itp_field_edge[i_iso],
                                                  &iso_itp_field_face[i_iso]);
    }
  }
  else {
    PDM_malloc(iso_itp_dfield_vtx , n_iso, double *);
    PDM_malloc(iso_itp_dfield_edge, n_iso, double *);
    PDM_malloc(iso_itp_dfield_face, n_iso, double *);
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
      PDM_isosurface_test_utils_dist_interpolation(isos, i_iso,
                                                   itp_dfield_vtx ,
                                                   itp_dfield_face,
                                                   itp_dfield_cell,
                                                  &iso_itp_dfield_vtx [i_iso],
                                                  &iso_itp_dfield_edge[i_iso],
                                                  &iso_itp_dfield_face[i_iso]);
    }
  }


  /*
   *  Visu isosurfaces
   */
  if (visu) {

    for (int i_iso = 0; i_iso < n_iso; i_iso++) {

      if (n_part_out > 0) {
        PDM_isosurface_test_utils_part_vtk(isos, i_iso, n_part_out,
                                           iso_itp_field_vtx [i_iso],
                                           iso_itp_field_edge[i_iso],
                                           iso_itp_field_face[i_iso],
                                           comm);
      }
      else {
        PDM_isosurface_test_utils_dist_vtk(isos, i_iso,
                                           iso_itp_dfield_vtx [i_iso],
                                           iso_itp_dfield_edge[i_iso],
                                           iso_itp_dfield_face[i_iso],
                                           comm);
      }

    }
  }


  /*
   *  Free memory
   */
  PDM_isosurface_free(isos);

  if (n_part_in > 0) {
    PDM_multipart_free(mpart);

    if (use_part_mesh) {
      PDM_part_mesh_free(pmesh);
    }

    for (int i_part = 0; i_part < n_part_in; i_part++) {
      PDM_free(iso_field     [i_part]);
      PDM_free(itp_field_vtx [i_part]);
      PDM_free(itp_field_face[i_part]);
      PDM_free(itp_field_cell[i_part]);
    }
    PDM_free(iso_field);
    PDM_free(itp_field_vtx);
    PDM_free(itp_field_face);
    PDM_free(itp_field_cell);

    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      for (int i_part = 0; i_part < n_part_out; i_part++) {
        PDM_free(iso_itp_field_vtx [i_iso][i_part]);
        PDM_free(iso_itp_field_edge[i_iso][i_part]);
        PDM_free(iso_itp_field_face[i_iso][i_part]);
      }
      PDM_free(iso_itp_field_vtx [i_iso]);
      PDM_free(iso_itp_field_edge[i_iso]);
      PDM_free(iso_itp_field_face[i_iso]);
    }
    PDM_free(iso_itp_field_vtx);
    PDM_free(iso_itp_field_edge);
    PDM_free(iso_itp_field_face);
  }
  else {
    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      PDM_free(iso_itp_dfield_vtx [i_iso]);
      PDM_free(iso_itp_dfield_edge[i_iso]);
      PDM_free(iso_itp_dfield_face[i_iso]);
    }
    PDM_free(iso_itp_dfield_vtx);
    PDM_free(iso_itp_dfield_edge);
    PDM_free(iso_itp_dfield_face);

    PDM_free(iso_dfield);
    PDM_free(itp_dfield_vtx);
    PDM_free(itp_dfield_face);
    PDM_free(itp_dfield_cell);
    PDM_dmesh_free(dmesh);
  }

  if (i_rank == 0) {
    printf("\nOK! :)\n");
    printf("-------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }
}