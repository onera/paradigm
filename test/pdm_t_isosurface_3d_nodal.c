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
#include "pdm_vtk.h"

#include "pdm_array.h"

#include "pdm_isosurface.h"
#include "pdm_dcube_nodal_gen.h"
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
  PDM_Mesh_nodal_elt_t  elt_type       = PDM_MESH_NODAL_TETRA4;
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


  /*
   *  Generate mesh
   */
  PDM_dcube_nodal_t     *dcube_nodal = NULL;
  PDM_dmesh_nodal_t     *dmn         = NULL;
  PDM_part_mesh_nodal_t *pmn         = NULL;

  double  *iso_dfield      = NULL;
  double  *itp_dfield_vtx  = NULL;
  double  *itp_dfield_face = NULL;
  double  *itp_dfield_cell = NULL;
  double **iso_field       = NULL;
  double **itp_field_vtx   = NULL;
  double **itp_field_face  = NULL;
  double **itp_field_cell  = NULL;

  if (n_part==0) {
    // Block-distributed
    dcube_nodal = PDM_dcube_nodal_gen_create(comm,
                                             n_vtx_seg,
                                             n_vtx_seg,
                                             n_vtx_seg,
                                             1.,
                                             0.,
                                             0.,
                                             0.,
                                             elt_type,
                                             1,
                                             PDM_OWNERSHIP_KEEP);
    PDM_dcube_nodal_gen_build(dcube_nodal);

    dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube_nodal);

    int     dn_vtx     = PDM_DMesh_nodal_n_vtx_get(dmn);
    double *dvtx_coord = PDM_DMesh_nodal_coord_get(dmn, PDM_OWNERSHIP_BAD_VALUE);

    PDM_malloc(iso_dfield     , dn_vtx, double);
    PDM_malloc(itp_dfield_vtx , dn_vtx, double);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, iso_dfield    );
    PDM_isosurface_test_utils_compute_itp_field(dn_vtx, dvtx_coord, itp_dfield_vtx);

    PDM_g_num_t *face_distri = PDM_DMesh_nodal_section_distri_std_get(dmn, PDM_GEOMETRY_KIND_SURFACIC, 0);
    int dn_face = face_distri[i_rank+1]-face_distri[i_rank];
    PDM_malloc(itp_dfield_face, dn_face, double);
    for (int i_face=0; i_face<dn_face; ++i_face) {
      itp_dfield_face[i_face] = (double) (face_distri[i_rank]+i_face);
    }

    PDM_g_num_t *cell_distri = PDM_DMesh_nodal_section_distri_std_get(dmn, PDM_GEOMETRY_KIND_VOLUMIC, 0);
    int dn_cell = cell_distri[i_rank+1]-cell_distri[i_rank];
    PDM_malloc(itp_dfield_cell, dn_cell, double);
    for (int i_cell=0; i_cell<dn_cell; ++i_cell) {
      itp_dfield_cell[i_cell] = (double) (cell_distri[i_rank]+i_cell);
    }

  }
  else {
    // Partitioned
    pmn = PDM_generate_mesh_parallelepiped(comm,
                                           elt_type,
                                           1,
                                           NULL,
                                           -0.5,
                                           -0.5,
                                           -0.5,
                                           1.,
                                           1.,
                                           1.,
                                           n_vtx_seg,
                                           n_vtx_seg,
                                           n_vtx_seg,
                                           n_part,
                                           PDM_SPLIT_DUAL_WITH_PARMETIS); // TODO: Allow various partitioning ?
    
    // > Fields initialisation
    PDM_malloc(iso_field     , n_part, double *);
    PDM_malloc(itp_field_vtx , n_part, double *);
    PDM_malloc(itp_field_face, n_part, double *);
    PDM_malloc(itp_field_cell, n_part, double *);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);

      PDM_malloc(iso_field     [i_part], n_vtx, double);
      PDM_malloc(itp_field_vtx [i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field    [i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field_vtx[i_part]);

      int n_face = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_SURFACIC, i_part);
      PDM_g_num_t *face_gnum = PDM_part_mesh_nodal_g_num_get(pmn, 0, i_part, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_face[i_part], n_face, double);
      for (int i_face=0; i_face<n_face; ++i_face) {
        itp_field_face[i_part][i_face] = (double) face_gnum[i_face];
      }

      int n_cell = PDM_part_mesh_nodal_n_elmts_get(pmn, PDM_GEOMETRY_KIND_VOLUMIC, i_part);
      PDM_g_num_t *cell_gnum = PDM_part_mesh_nodal_g_num_get(pmn, 0, i_part, PDM_OWNERSHIP_KEEP);
      PDM_malloc(itp_field_cell[i_part], n_cell, double);
      for (int i_cell=0; i_cell<n_cell; ++i_cell) {
        itp_field_cell[i_part][i_cell] = (double) cell_gnum[i_cell];
      }

    }
  }


  /*
   *  TODO:
   *    - extract bc
   *    - plusieurs isovalues
   *    - reequilibrate/local
   *    - test reset (et chp variable ?)
   */



  /*
   *  Creating isosurface object
   */
  PDM_isosurface_t *isos = PDM_isosurface_create(comm, 3);
  if (n_part==0) {
    PDM_isosurface_dmesh_nodal_set(isos, dmn);
  } else {
    PDM_isosurface_mesh_nodal_set(isos, pmn);
    if (local==0) {
      PDM_isosurface_redistribution_set(isos, PDM_EXTRACT_PART_KIND_REEQUILIBRATE, PDM_SPLIT_DUAL_WITH_HILBERT); // TODO: Test various partitioning ?
    }
  }


  /*
   *  Add isosurface parameters
   */

  // > Plane isosurface
  double plane_equation [4] = {1.,0.,0.,0.5};
  double plane_isovalues[3] = {-0.50,0.,0.50};
  int iso1 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_PLANE,
                                3,
                                plane_isovalues);
  PDM_isosurface_equation_set(isos,
                              iso1,
                              plane_equation);

  // > User field isosurface

  double field_isovalues[1] = {0.1};
  int iso2 = PDM_isosurface_add(isos, 
                                PDM_ISO_SURFACE_KIND_FIELD,
                                1,
                                field_isovalues);
  if (n_part==0) {
    PDM_isosurface_dfield_set(isos,
                              iso2,
                              iso_dfield);
  }
  else {
    for (int i_part=0; i_part<n_part; ++i_part) {
      PDM_isosurface_field_set(isos,
                               iso2,
                               i_part,
                               iso_field[i_part]);
    }
  }

  /*
   *  Compute isosurface
   */
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

    PDM_isosurface_enable_part_to_part(isos,
                                       i_iso,
                                       PDM_MESH_ENTITY_FACE,
                                       0);
  }

  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_reset(isos, iso1);
  double plane_isovalues2[2] = {-0.30,0.30};
  PDM_isosurface_set_isovalues(isos, iso1, 2, plane_isovalues2);
  PDM_isosurface_compute(isos, iso1);
  PDM_isosurface_compute(isos, iso2);


  /*
   *  Interpolate field
   */
  double  **iso_itp_dfield_vtx  = NULL;
  double  **iso_itp_dfield_edge = NULL;
  double  **iso_itp_dfield_face = NULL;
  double ***iso_itp_field_vtx   = NULL;
  double ***iso_itp_field_edge  = NULL;
  double ***iso_itp_field_face  = NULL;

  if (n_part > 0) {
    PDM_malloc(iso_itp_field_vtx  , n_iso, double **);
    PDM_malloc(iso_itp_field_edge , n_iso, double **);
    PDM_malloc(iso_itp_field_face , n_iso, double **);
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
      PDM_isosurface_test_utils_part_interpolation(isos, i_iso, n_part, local,
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
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {
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

      if (n_part > 0) {
        PDM_isosurface_test_utils_part_vtk(isos, i_iso, n_part,
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
  PDM_free(isovalues);

  if (n_part > 0) {
    PDM_part_mesh_nodal_free(pmn);

    for (int i_part = 0; i_part < n_part; i_part++) {
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
      for (int i_part = 0; i_part < n_part; i_part++) {
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
    PDM_dcube_nodal_gen_free(dcube_nodal);

    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      PDM_free(iso_itp_dfield_vtx [i_iso]);
      // PDM_free(iso_itp_dfield_edge[i_iso]);
      // PDM_free(iso_itp_dfield_face[i_iso]);
    }
    PDM_free(iso_itp_dfield_vtx);
    PDM_free(iso_itp_dfield_edge);
    PDM_free(iso_itp_dfield_face);

    PDM_free(iso_dfield);
    PDM_free(itp_dfield_vtx);
    PDM_free(itp_dfield_face);
    PDM_free(itp_dfield_cell);
  }


  if (i_rank == 0) {
    printf("End :D\n");
    fflush(stdout);
  }


  PDM_MPI_Finalize();


  return EXIT_SUCCESS;
}
