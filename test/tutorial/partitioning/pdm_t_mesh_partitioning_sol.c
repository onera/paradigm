#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_dcube_nodal_gen.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"
#include "pdm_vtk.h"
#include "pdm_printf.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -vtk   Mesh vtk output.\n\n"
     "  -fe    Get FE mesh structure instead of FV.\n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int    argc,
           char **argv,
           int   *vtk,
           int   *fe)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-vtk") == 0) {
      *vtk = 1;
    }

    else if (strcmp(argv[i], "-fe") == 0) {
      *fe = 1;
    }

    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

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
  int vtk = 0;
  int fe  = 0;

  _read_args(argc,
             argv,
             &vtk,
             &fe);

  // Initialize MPI environment
  int          i_rank = -1;
  int          n_rank = -1;
  PDM_MPI_Comm comm   = PDM_MPI_COMM_WORLD;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Generate block-distributed parallelepided mesh
  int                  n_x      = 10;
  int                  n_y      = 10;
  int                  n_z      = 10;
  double               lengthx  = 1.;
  double               xmin     = 0.;
  double               ymin     = 0.;
  double               zmin     = 0.;
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_TETRA4;
  int                  order    = 1; // call PDM_dcube_nodal_gen_ordering_set if order > 1
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

  PDM_dcube_nodal_gen_build (dcube);

  PDM_dmesh_nodal_t *dmn = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube);

  PDM_dmesh_nodal_generate_distribution(dmn);

  // free
  PDM_dcube_nodal_gen_free(dcube);

  // Create partitioning object
  int              n_zone      = 1; // fixed
  int              n_part      = 1; // fixed
  PDM_split_dual_t part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
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

  int i_section = 0; // fixed
  int i_part    = 0; // fixed

  // Get mesh arrrays in FE structure
  if (fe) {
    PDM_part_mesh_nodal_t *pmn  = NULL;
    PDM_multipart_get_part_mesh_nodal(mpart,
                                      0,
                                      &pmn,
                                      PDM_OWNERSHIP_USER);

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn,
                                                      i_section,
                                                      i_part);

    int         *connec              = NULL;
    PDM_g_num_t *numabs              = NULL;
    int         *parent_num          = NULL;
    PDM_g_num_t *parent_entity_g_num = NULL;
    PDM_part_mesh_nodal_section_std_get(pmn,
                                        i_section,
                                        i_part,
                                        &connec,
                                        &numabs,
                                        &parent_num,
                                        &parent_entity_g_num,
                                        PDM_OWNERSHIP_KEEP);

    int *elt_vtx_idx = PDM_array_new_idx_from_const_stride_int(4, n_elt);
    int *elt_vtx = malloc(sizeof(int) * elt_vtx_idx[n_elt]);
    memcpy(elt_vtx, connec, sizeof(int) * elt_vtx_idx[n_elt]);

    PDM_g_num_t *elt_ln_to_gn = PDM_part_mesh_nodal_g_num_get(pmn,
                                                              i_section,
                                                              i_part,
                                                              PDM_OWNERSHIP_USER);

    double *coords = NULL;
    int n_vtx = PDM_multipart_part_vtx_coord_get(mpart,
                                                 i_section,
                                                 i_part,
                                                 &coords,
                                                 PDM_OWNERSHIP_USER);

    PDM_g_num_t *vtx_ln_to_gn = PDM_part_mesh_nodal_vtx_g_num_get(pmn,
                                                                  i_part);

    // Visualisation cell->vtx
    if (vtk) {
      char filename[999];
      sprintf(filename, "partitioning_FE_%3.3d.vtk", i_rank);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 coords,
                                 vtx_ln_to_gn,
                                 elt_type,
                                 n_elt,
                                 elt_vtx,
                                 elt_ln_to_gn,
                                 0,
                                 NULL,
                                 NULL);
    }

    // free
    free(elt_ln_to_gn);
    free(elt_vtx_idx);
    free(elt_vtx);
    free(coords);
    PDM_part_mesh_nodal_free(pmn);
  }

  // Get mesh arrrays in FV structure
  else {
    PDM_g_num_t *vtx_ln_to_gn = NULL;
    int n_vtx = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_section, // i_zone == i_section ?
                                                i_part,
                                                PDM_MESH_ENTITY_VERTEX,
                                                &vtx_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    double *coords = NULL;
    PDM_multipart_part_vtx_coord_get(mpart,
                                     i_section,
                                     i_part,
                                     &coords,
                                     PDM_OWNERSHIP_USER);

    PDM_g_num_t *edge_ln_to_gn = NULL;
    int n_edge = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_section,
                                                i_part,
                                                PDM_MESH_ENTITY_EDGE,
                                                &edge_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    int *edge_vtx_idx = NULL;
    int *edge_vtx     = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_section,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                        &edge_vtx,
                                        &edge_vtx_idx,
                                        PDM_OWNERSHIP_USER);

    if (edge_vtx_idx != NULL) free (edge_vtx_idx);

    PDM_g_num_t *face_ln_to_gn = NULL;
    int n_face = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_section,
                                                i_part,
                                                PDM_MESH_ENTITY_FACE,
                                                &face_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    int *face_edge_idx = NULL;
    int *face_edge     = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_section,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_FACE_EDGE,
                                        &face_edge,
                                        &face_edge_idx,
                                        PDM_OWNERSHIP_USER);

    PDM_g_num_t *cell_ln_to_gn = NULL;
    int n_cell = PDM_multipart_part_ln_to_gn_get(mpart,
                                                i_section,
                                                i_part,
                                                PDM_MESH_ENTITY_CELL,
                                                &cell_ln_to_gn,
                                                PDM_OWNERSHIP_USER);

    int *cell_face_idx = NULL;
    int *cell_face     = NULL;
    PDM_multipart_part_connectivity_get(mpart,
                                        i_section,
                                        i_part,
                                        PDM_CONNECTIVITY_TYPE_CELL_FACE,
                                        &cell_face,
                                        &cell_face_idx,
                                        PDM_OWNERSHIP_USER);

    // Visualisation edge->vtx
    if (vtk) {
      char filename[999];
      sprintf(filename, "partitioning_FV_%3.3d.vtk", i_rank);
      PDM_vtk_write_std_elements(filename,
                                 n_vtx,
                                 coords,
                                 vtx_ln_to_gn,
                                 PDM_MESH_NODAL_BAR2,
                                 n_edge,
                                 edge_vtx,
                                 edge_ln_to_gn,
                                 0,
                                 NULL,
                                 NULL);
    }


    // Use PDM_compute_face_vtx_from_face_and_edge if you need face->vtx connectivity

    // free
    free(vtx_ln_to_gn);
    free(coords);
    free(edge_ln_to_gn);
    free(edge_vtx_idx);
    free(edge_vtx);
    free(face_ln_to_gn);
    free(face_edge_idx);
    free(face_edge);
    free(cell_ln_to_gn);
    free(cell_face_idx);
    free(cell_face);
  }

  // free
  PDM_DMesh_nodal_free(dmn);
  PDM_multipart_free(mpart);

  // Finalize MPI environment
  PDM_MPI_Finalize();

  return 0;
}
