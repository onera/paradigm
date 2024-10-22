#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_config.h"
#include "pdm_mpi.h"

#include "pdm_logging.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#include "pdm_isosurface.h"
#include "pdm_generate_mesh.h"

#include "pdm_isosurface_test_utils.h"

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
int main(int argc, char *argv[])
{

  /*
   *  Init MPI
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;
  PDM_MPI_Init(&argc, &argv);
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
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_TRIA3;
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


  /*
   *  Generate mesh
   */
  PDM_dmesh_nodal_t     *dmn = NULL;
  PDM_part_mesh_nodal_t *pmn = NULL;

  PDM_isosurface_test_utils_gen_mesh_nodal(comm,
                                           mesh_name,
                                           n_part,
                                           n_vtx_seg,
                                           randomize,
                                           elt_type,
                                           &pmn,
                                           &dmn);

  if (n_part==0) {
    // Block-distributed
  }
  else {
    // Partitioned
  }


  /*
   *  Creating isosurface object
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 2);
  if (n_part==0) {
    PDM_isosurface_dmesh_nodal_set(isos, dmn);
  } else {
    PDM_isosurface_mesh_nodal_set(isos, pmn);
    if (local==0) {
      PDM_isosurface_redistribution_set(isos, PDM_EXTRACT_PART_KIND_REEQUILIBRATE, PDM_SPLIT_DUAL_WITH_HILBERT); // TODO: Test various partitioning ?
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

  if (n_part==0) {
    int     dn_vtx     = PDM_DMesh_nodal_n_vtx_get(dmn);
    double *dvtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(iso_dfield     , dn_vtx, double);
    PDM_malloc(itp_dfield_vtx , dn_vtx, double);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, iso_dfield     );
    PDM_isosurface_test_utils_compute_itp_field(dn_vtx, dvtx_coord, itp_dfield_vtx );

    PDM_g_num_t *face_distri = PDM_DMesh_nodal_section_distri_std_get(dmn, PDM_GEOMETRY_KIND_SURFACIC, 0);
    int dn_face = face_distri[i_rank+1]-face_distri[i_rank];
    PDM_malloc(itp_dfield_face, dn_face, double);
    for (int i_face=0; i_face<dn_face; ++i_face) {
      itp_dfield_face[i_face] = (double) (face_distri[i_rank]+i_face);
    }

  } else {
    PDM_malloc(iso_field     , n_part, double *);
    PDM_malloc(itp_field_vtx , n_part, double *);
    PDM_malloc(itp_field_face, n_part, double *);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);
      
      PDM_malloc(iso_field     [i_part], n_vtx, double);
      PDM_malloc(itp_field_vtx [i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field     [i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field_vtx [i_part]);

      int n_face = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part);
      PDM_g_num_t *face_gnum = PDM_part_mesh_nodal_g_num_get(pmn, 0, i_part, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_face[i_part], n_face, double);
      for (int i_face=0; i_face<n_face; ++i_face) {
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
  if (n_part > 0) { // Partitioned
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_isosurface_field_set(isos, iso2, i_part, iso_field[i_part]);
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
    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_VTX,
                                       0);

    PDM_isosurface_enable_part_to_part(isos,
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
  if (visu==1) {
    
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


  /*
   * TODO: 
   *   - if default test conf, put assert ?
   */


  /*
   *  Free objects
   */
  PDM_isosurface_free(isos);
  PDM_free(isovalues);

  if (n_part==0) {
    PDM_DMesh_nodal_free(dmn);
  } else {
    PDM_part_mesh_nodal_free(pmn);
  }

  PDM_free(iso_dfield);
  PDM_free(itp_dfield_vtx);
  PDM_free(itp_dfield_face);

  if (iso_field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (iso_field[i_part]!=NULL) {
        PDM_free(iso_field[i_part]);
      }
    }
    PDM_free(iso_field);
  }
  if (itp_field_vtx!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      PDM_free(itp_field_vtx [i_part]);
      PDM_free(itp_field_face[i_part]);
    }
    PDM_free(itp_field_vtx );
    PDM_free(itp_field_face);
  }

  if (iso_itp_dfield_vtx!=NULL) {
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
      PDM_free(iso_itp_dfield_vtx [i_iso]);
      PDM_free(iso_itp_dfield_edge[i_iso]);
    }
  }
  PDM_free(iso_itp_dfield_vtx );
  PDM_free(iso_itp_dfield_edge);

  if (iso_itp_field_vtx!=NULL) {
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
      for (int i_part=0; i_part<n_part; ++i_part) {
        PDM_free(iso_itp_field_vtx [i_iso][i_part]);
        PDM_free(iso_itp_field_edge[i_iso][i_part]);
      }
      PDM_free(iso_itp_field_vtx [i_iso]);
      PDM_free(iso_itp_field_edge[i_iso]);
    }
  }
  PDM_free(iso_itp_field_vtx );
  PDM_free(iso_itp_field_edge);

  PDM_MPI_Finalize();

  return 0;
}
