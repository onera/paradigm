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

  double  *iso_dfield = NULL;
  double  *itp_dfield = NULL;
  double **iso_field  = NULL;
  double **itp_field  = NULL;

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

    PDM_malloc(iso_dfield, dn_vtx, double);
    PDM_malloc(itp_dfield, dn_vtx, double);
    PDM_isosurface_test_utils_compute_iso_field(dn_vtx, dvtx_coord, iso_dfield);
    PDM_isosurface_test_utils_compute_itp_field(dn_vtx, dvtx_coord, itp_dfield);

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
    PDM_malloc(iso_field, n_part, double *);
    PDM_malloc(itp_field, n_part, double *);

    for (int i_part=0; i_part<n_part; ++i_part) {
      int     n_vtx     = PDM_part_mesh_nodal_n_vtx_get(pmn, i_part);
      double *vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, i_part);

      PDM_malloc(iso_field[i_part], n_vtx, double);
      PDM_malloc(itp_field[i_part], n_vtx, double);
      PDM_isosurface_test_utils_compute_iso_field(n_vtx, vtx_coord, iso_field[i_part]);
      PDM_isosurface_test_utils_compute_itp_field(n_vtx, vtx_coord, itp_field[i_part]);
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

  int n_iso = iso2+1;
  if (n_part>0) {
    for (int i_iso = 0; i_iso < n_iso; i_iso++) {
      PDM_isosurface_enable_part_to_part(isos,
                                         i_iso,
                                         PDM_MESH_ENTITY_VTX,
                                         0);
    }
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
  double  **iso_itp_dfield = NULL;
  double ***iso_itp_field  = NULL;
  PDM_malloc(iso_itp_dfield, n_iso, double  *);
  PDM_malloc(iso_itp_field , n_iso, double **);
  for (int i_iso=0; i_iso<n_iso; ++i_iso) {
    iso_itp_dfield[i_iso] = NULL;
    iso_itp_field [i_iso] = NULL;
  } 
  if (n_part==0) {
    for (int i_iso=0; i_iso<n_iso; ++i_iso) {

      PDM_isosurface_test_utils_dist_interpolation(isos, i_iso,
                                                       itp_dfield,
                                                   NULL,
                                                   NULL,
                                                  &iso_itp_dfield[i_iso],
                                                   NULL,
                                                   NULL);
    }
  }
  else {
    if (local==1) {
      // Local
      for (int i_iso=0; i_iso<n_iso; ++i_iso) {
        PDM_malloc(iso_itp_field[i_iso], n_part, double *);
        for (int i_part=0; i_part<n_part; ++i_part) {
          
          int    *vtx_parent_idx  = NULL;
          int    *vtx_parent_lnum = NULL;
          double *vtx_parent_wght = NULL;
          int iso_n_vtx = PDM_isosurface_local_parent_get(isos, i_iso, i_part, PDM_MESH_ENTITY_VTX, 
                                                         &vtx_parent_idx, &vtx_parent_lnum,
                                                          PDM_OWNERSHIP_KEEP);
          PDM_isosurface_vtx_parent_weight_get(isos, i_iso, i_part, &vtx_parent_idx, &vtx_parent_wght, PDM_OWNERSHIP_KEEP);
          
          iso_itp_field[i_iso][i_part] = PDM_array_zeros_double(iso_n_vtx);
          for (int i_iso_vtx=0; i_iso_vtx<iso_n_vtx; ++i_iso_vtx) {
            int i_beg_parent = vtx_parent_idx[i_iso_vtx  ];
            int i_end_parent = vtx_parent_idx[i_iso_vtx+1];
            for (int i_parent=i_beg_parent; i_parent<i_end_parent; ++i_parent) {
              int    parent_lnum = vtx_parent_lnum[i_parent];
              double parent_wght = vtx_parent_wght[i_parent];
              iso_itp_field[i_iso][i_part][i_iso_vtx] += itp_field[i_part][parent_lnum-1]*parent_wght;
            }
          }

          printf("DO for edge and faces\n");
        }
      }
    }
    else {
      // Reequilibrate
      for (int i_iso=0; i_iso<n_iso; ++i_iso) {
        PDM_malloc(iso_itp_field[i_iso], n_part, double *);

        PDM_part_to_part_t *ptp_vtx = NULL;
        PDM_isosurface_part_to_part_get(isos,
                                        i_iso,
                                        PDM_MESH_ENTITY_VTX,
                                        &ptp_vtx,
                                        PDM_OWNERSHIP_KEEP);

        double **recv_vtx_field = NULL;
        int request_vtx = -1;
        PDM_part_to_part_reverse_iexch(ptp_vtx,
                                       PDM_MPI_COMM_KIND_P2P,
                                       PDM_STRIDE_CST_INTERLACED,
                                       PDM_PART_TO_PART_DATA_DEF_ORDER_PART2,
                                       1,
                                       sizeof(double),
                                       NULL,
                      (const void  **) itp_field,
                                       NULL,
                      (      void ***) &recv_vtx_field,
                                       &request_vtx);

        PDM_part_to_part_reverse_iexch_wait(ptp_vtx, request_vtx);

        for (int i_part = 0; i_part < n_part; i_part++) {
          int    *iso_vtx_parent_idx;
          double *iso_vtx_parent_weight;
          int iso_n_vtx = PDM_isosurface_vtx_parent_weight_get(isos,
                                                               i_iso,
                                                               i_part,
                                                               &iso_vtx_parent_idx,
                                                               &iso_vtx_parent_weight,
                                                               PDM_OWNERSHIP_KEEP);

          PDM_calloc(iso_itp_field[i_iso][i_part], iso_n_vtx, double);
          for (int i_vtx = 0; i_vtx < iso_n_vtx; i_vtx++) {
            for (int i = iso_vtx_parent_idx[i_vtx]; i < iso_vtx_parent_idx[i_vtx+1]; i++) {
              iso_itp_field[i_iso][i_part][i_vtx] += iso_vtx_parent_weight[i] * recv_vtx_field[i_part][i];
            }
          }

          PDM_free(recv_vtx_field[i_part]);
        } // End loop on parts
        PDM_free(recv_vtx_field);
      }

      // PDM_error(__FILE__, __LINE__, 0, "Part entry with local=0 not implemented\n");
    }
  }


  /*
   *  Visu isosurfaces
   */
  if (visu==1) {
    if (n_part==0) {
      for (int i_iso = 0; i_iso < n_iso; i_iso++) {
        PDM_isosurface_test_utils_dist_vtk(isos, i_iso, comm);
      }
    }
    else {
      // > iso line output
      for (int i_iso = 0; i_iso < n_iso; i_iso++) {
        PDM_isosurface_test_utils_part_vtk(isos, i_iso, n_part, iso_itp_field, comm);
      }
    }
  }


  /*
   *  Free objects
   */
  PDM_isosurface_free(isos);
  if (n_part==0) {
    PDM_dcube_nodal_gen_free(dcube_nodal);
    PDM_free(iso_dfield);
    PDM_free(itp_dfield);
  } else {
    PDM_part_mesh_nodal_free(pmn);
  }

  if (iso_field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (iso_field[i_part]!=NULL) {
        PDM_free(iso_field[i_part]);
      }
    }
    PDM_free(iso_field);
  }
  if (itp_field!=NULL) {
    for (int i_part=0; i_part<n_part; ++i_part) {
      if (itp_field[i_part]!=NULL) {
        PDM_free(itp_field[i_part]);
      }
    }
    PDM_free(itp_field);
  }
  for (int id_iso=0; id_iso<n_iso; ++id_iso) {
    PDM_free(iso_itp_dfield[id_iso]);
    if (iso_itp_field[id_iso]!=NULL) {
      for (int i_part=0; i_part<n_part; ++i_part) {
        if (iso_itp_field[id_iso][i_part]!=NULL) {
          PDM_free(iso_itp_field[id_iso][i_part]);
        }
      }
      PDM_free(iso_itp_field[id_iso]);
    }
  }
  PDM_free(iso_itp_field);
  PDM_free(iso_itp_dfield);


  PDM_MPI_Finalize();

  return 0;
}
