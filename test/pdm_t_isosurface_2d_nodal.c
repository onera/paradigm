
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "pdm_distrib.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_error.h"
#include "pdm_isosurface_test_utils.h"
#include "pdm_isosurface.h"
#include "pdm_mem_tool.h"
#include "pdm_mesh_nodal.h"
#include "pdm_mpi.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

 static void _run_test_2d_nodal(PDM_isosurface_test_utils_params_t params);

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
     params.mesh_name      = NULL;
     params.sol_name       = NULL;
     params.n_vtx_seg      = 10;
     params.randomize      = 0;
     params.use_part_mesh  = 0; // irrelevant here
     params.generate_edges = 0; // irrelevant here
     params.visu           = 0;

     // Test combinations of variable parameters
     for (int n_part_in = 0; n_part_in <= 2; n_part_in++) {

       params.n_part_in = n_part_in;

       for (int n_isovalues = 1; n_isovalues <= 3; n_isovalues += 2) {

         params.n_isovalues = n_isovalues;
         params.isovalues   = default_isovalues;

         for (int simplex = 0; simplex <= 1; simplex++) {

           if (simplex) {
             params.elt_type = PDM_MESH_NODAL_TRIA3;
           }
           else {
             params.elt_type = PDM_MESH_NODAL_POLY_2D;
           }

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

                 _run_test_2d_nodal(params);
                 n_tests++;

               } // End loop on use_group

             } // End loop on n_part_out

           } // End loop on local/redistribute

         } // End loop on element type

       } // End loop on n_isovalues

     } // End loop on n_part_in
   }

   else {
     /* Run test with user-provided parameters */

     // Initialize parameters to their default values
     params.mesh_name      = NULL;
     params.sol_name       = NULL;
     params.n_part_in      = 1;
     params.n_part_out     = 1;
     params.visu           = 0;
     params.n_isovalues    = 1;
     params.isovalues      = default_isovalues;
     params.elt_type       = PDM_MESH_NODAL_TRIA3;
     params.n_vtx_seg      = 10;
     params.randomize      = 0;
     params.use_part_mesh  = 0; // irrelevant here
     params.generate_edges = 0; // irrelevant here
     params.local          = 0;
     params.use_groups     = 0;

     // Overwrite with user parameter values
     PDM_isosurface_test_utils_read_args(argc,
                                         argv,
                                         &params.n_part_in,
                                         &params.n_part_out,
                                         &params.mesh_name,
                                         &params.sol_name,
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
     _run_test_2d_nodal(params);
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
 _run_test_2d_nodal
 (
   PDM_isosurface_test_utils_params_t params
 )
 {
   PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
   int i_rank;
   int n_rank;
   PDM_MPI_Comm_rank(comm, &i_rank);
   PDM_MPI_Comm_size(comm, &n_rank);

   PDM_isosurface_test_utils_isosurface_params_dump(comm, "pdm_t_isosurface_2d_nodal", params);
   PDM_MPI_Barrier(comm);


   char                 *mesh_name   = params.mesh_name;
   // char                 *sol_name    = params.sol_name;
   int                   n_part_in   = params.n_part_in;
   int                   n_part_out  = params.n_part_out;
   int                   n_isovalues = params.n_isovalues;
   double               *isovalues   = params.isovalues;
   PDM_Mesh_nodal_elt_t  elt_type    = params.elt_type;
   PDM_g_num_t           n_vtx_seg   = params.n_vtx_seg;
   int                   randomize   = params.randomize;
   int                   local       = params.local;
   int                   use_groups  = params.use_groups;
   int                   visu        = params.visu;

   int elt_dim = PDM_Mesh_nodal_elt_dim_get(elt_type);
   if (elt_dim != 2) {
     PDM_error(__FILE__, __LINE__, 0,
               "Invalid element dimension for element type %d (expected 2, got %d)",
               elt_type, elt_dim);
   }


  /*
   *  Generate mesh
   */
  PDM_dmesh_nodal_t     *dmn = NULL;
  PDM_part_mesh_nodal_t *pmn = NULL;

  PDM_isosurface_test_utils_gen_mesh_nodal(comm,
                                           mesh_name,
                                           n_part_in,
                                           n_vtx_seg,
                                           randomize,
                                           elt_type,
                                           use_groups,
                                           &pmn,
                                           &dmn);


  /*
   *  Create isosurface object
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 2);
  if (n_part_in == 0) {
    PDM_isosurface_dmesh_nodal_set(isos, dmn);
  }
  else {
    PDM_isosurface_part_mesh_nodal_set(isos, pmn);
    if (local == 0) {
      PDM_isosurface_redistribution_set(isos, PDM_EXTRACT_PART_KIND_REEQUILIBRATE, PDM_SPLIT_DUAL_WITH_HILBERT); // TODO: Test various partitioning ?
      PDM_isosurface_n_part_out_set(isos, n_part_out);
    }
  }


  /*
   *  Compute isosurface and interpolation field
   */

  double  *iso_dfield      = NULL;
  double  *itp_dfield_vtx  = NULL;
  double  *itp_dfield_face = NULL;
  double **iso_field       = NULL;
  double **itp_field_vtx   = NULL;
  double **itp_field_face  = NULL;

  if (n_part_in == 0) {
    // Block-distributed
    int     dn_vtx     = PDM_DMesh_nodal_n_vtx_get(dmn);
    double *dvtx_coord = PDM_DMesh_nodal_vtx_get(dmn, PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(iso_dfield    , dn_vtx, double);
    PDM_malloc(itp_dfield_vtx, dn_vtx, double);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, iso_dfield    );
    PDM_isosurface_test_utils_compute_itp_field(dn_vtx, dvtx_coord, itp_dfield_vtx);

    int dn_face = 0;
    int  n_section   = PDM_DMesh_nodal_n_section_get  (dmn, PDM_GEOMETRY_KIND_SURFACIC);
    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmn, PDM_GEOMETRY_KIND_SURFACIC);
    for (int i_section = 0; i_section < n_section; i_section++) {
      int n_elt = PDM_DMesh_nodal_section_n_elt_get(dmn, PDM_GEOMETRY_KIND_SURFACIC, sections_id[i_section]);
      dn_face += n_elt;
    }

    PDM_g_num_t *face_distri = PDM_compute_entity_distribution(comm, dn_face);
    PDM_malloc(itp_dfield_face, dn_face, double);
    for (int i_face = 0; i_face < dn_face; ++i_face) {
      itp_dfield_face[i_face] = (double) (face_distri[i_rank]+i_face+1);
    }
    PDM_free(face_distri);

  }
  else {
    // Partitioned
    PDM_malloc(iso_field     , n_part_in, double *);
    PDM_malloc(itp_field_vtx , n_part_in, double *);
    PDM_malloc(itp_field_face, n_part_in, double *);

    for (int i_part = 0; i_part < n_part_in; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part, PDM_OWNERSHIP_BAD_VALUE);

      PDM_malloc(iso_field    [i_part], n_vtx, double);
      PDM_malloc(itp_field_vtx[i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field    [i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field_vtx[i_part]);

      int n_face = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part);
      PDM_g_num_t *face_gnum = PDM_part_mesh_nodal_g_num_get_from_part(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_face[i_part], n_face, double);
      for (int i_face = 0; i_face < n_face; ++i_face) {
        itp_field_face[i_part][i_face] = (double) face_gnum[i_face];
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
   *  Compute isosurfaces
   */
  int n_iso = iso3+1;
  for (int i_iso = 0; i_iso < n_iso; i_iso++) {
    PDM_isosurface_part_to_part_enable(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_VTX,
                                       0);

    PDM_isosurface_part_to_part_enable(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_EDGE,
                                       0);
  }

  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset  (isos, iso1);
  PDM_isosurface_reset  (isos, iso2);

  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);
  PDM_isosurface_compute(isos, iso3);

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
      printf("isosurface %d : "PDM_FMT_G_NUM" vtx, "PDM_FMT_G_NUM" edges\n",
             i_iso, gn_iso_vtx, gn_iso_edge);
      fflush(stdout);
    }
  }


  /*
   *  Interpolate field
   */
  double  **iso_itp_dfield_vtx  = NULL;
  double  **iso_itp_dfield_edge = NULL;
  double ***iso_itp_field_vtx   = NULL;
  double ***iso_itp_field_edge  = NULL;

  if (n_part_out > 0) {
    PDM_malloc(iso_itp_field_vtx  , n_iso, double **);
    PDM_malloc(iso_itp_field_edge , n_iso, double **);
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
      PDM_isosurface_test_utils_part_interpolation(isos, i_iso, n_part_out, local,
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
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
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

      if (n_part_out > 0) {
        PDM_isosurface_test_utils_part_vtk(isos, i_iso, n_part_out,
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

  /*
   *  Free objects
   */
  PDM_isosurface_free(isos);

  if (n_part_in == 0) {
    PDM_DMesh_nodal_free(dmn);
  } else {
    PDM_part_mesh_nodal_free(pmn);
  }

  PDM_free(iso_dfield);
  PDM_free(itp_dfield_vtx);
  PDM_free(itp_dfield_face);

  if (iso_field != NULL) {
    for (int i_part = 0; i_part < n_part_in; ++i_part) {
      PDM_free(iso_field[i_part]);
    }
    PDM_free(iso_field);
  }
  if (itp_field_vtx != NULL) {
    for (int i_part = 0; i_part < n_part_in; ++i_part) {
      PDM_free(itp_field_vtx [i_part]);
      PDM_free(itp_field_face[i_part]);
    }
    PDM_free(itp_field_vtx );
    PDM_free(itp_field_face);
  }

  if (iso_itp_dfield_vtx != NULL) {
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
      PDM_free(iso_itp_dfield_vtx [i_iso]);
      PDM_free(iso_itp_dfield_edge[i_iso]);
    }
  }
  PDM_free(iso_itp_dfield_vtx );
  PDM_free(iso_itp_dfield_edge);

  if (iso_itp_field_vtx != NULL) {
    for (int i_iso = 0; i_iso < n_iso; ++i_iso) {
      for (int i_part = 0; i_part < n_part_out; ++i_part) {
        PDM_free(iso_itp_field_vtx [i_iso][i_part]);
        PDM_free(iso_itp_field_edge[i_iso][i_part]);
      }
      PDM_free(iso_itp_field_vtx [i_iso]);
      PDM_free(iso_itp_field_edge[i_iso]);
    }
  }
  PDM_free(iso_itp_field_vtx );
  PDM_free(iso_itp_field_edge);

  if (i_rank == 0) {
    printf("\nOK! :)\n");
    printf("-------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }
}